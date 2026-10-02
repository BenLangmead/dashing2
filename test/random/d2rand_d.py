"""Helpers for suite D (input processing and k-mer generation modes), Python standard library only.

This module reimplements, in Python, the parts of dashing2's k-mer pipeline
that decide which items enter a sketch, so that exact modes can be compared
item by item:

* the encodings of DNA and of the four protein alphabets (bonsai alphabet.h,
  encoder.h): 2 bits per base for DNA, first character most significant;
  base-20, base-14 and base-6 integers for --protein, --protein14 and
  --protein6; 3 bits per residue for --protein8;
* the mask applied before k-mers are stored (src/enums.h maskfn): XOR with a
  seed-dependent mask (0 for the default seed) and then the Wang 64-bit hash,
  applied to each 64-bit half for 128-bit k-mers;
* the window order used by -w (bonsai qmap.h, encoder.h): a window holds the
  last w - c + 1 valid k-mers of the record (c is the k-mer span; invalid
  k-mers are skipped, not counted), one item is emitted per full window, a
  record with fewer valid k-mers than the window emits its best one, and the
  best item has the smallest (score, value) pair; the score is the FRev64 (64
  bit) or CEHasher (128 bit) permutation of the k-mer, or for --entmin that
  hash divided by the Shannon entropy of the k-mer's characters;
* the --downsample decision (src/d2.h downsample_hash).

For k beyond the exact-encoding limit dashing2 uses a polynomial rolling hash
whose constants are not reproduced here; those paths are checked through
cardinalities and through structural relations between dashing2's own outputs.
"""
import bz2
import gzip
import lzma
import math
import os
import shutil
import struct
import subprocess

M64 = (1 << 64) - 1
M128 = (1 << 128) - 1

# ---------------------------------------------------------------------------
# Hashes


def wang(key):
    key = ((~key) + (key << 21)) & M64
    key ^= key >> 24
    key = (key + (key << 3) + (key << 8)) & M64
    key ^= key >> 14
    key = (key + (key << 2) + (key << 4)) & M64
    key ^= key >> 28
    key = (key + (key << 31)) & M64
    return key


def wang_inv(key):
    tmp = (key - (key << 31)) & M64
    key = (key - (tmp << 31)) & M64
    tmp = key ^ key >> 28
    key = key ^ tmp >> 28
    key = (key * 14933078535860113213) & M64
    tmp = key ^ key >> 14
    tmp = key ^ tmp >> 14
    tmp = key ^ tmp >> 14
    key = key ^ tmp >> 14
    key = (key * 15244667743933553977) & M64
    tmp = key ^ key >> 24
    key = key ^ tmp >> 24
    tmp = (~key) & M64
    tmp = (~(key - (tmp << 21))) & M64
    tmp = (~(key - (tmp << 21))) & M64
    key = (~(key - (tmp << 21))) & M64
    return key


class Mask:
    """dashing2's maskfn and its inverse for a given --seed (src/enums.cpp seed_mask)."""

    def __init__(self, seed=0):
        seed = int(seed) & M64
        if seed == 0:
            self.x64 = 0
            self.x128 = 0
        else:
            self.x64 = wang(seed)
            self.x128 = self.x64 | (wang(self.x64) << 64)

    def m64(self, x):
        return wang(x ^ self.x64)

    def m128(self, x):
        x ^= self.x128
        return wang(x & M64) | (wang(x >> 64) << 64)

    def mask(self, x, long):
        return self.m128(x) if long else self.m64(x)

    def unmask(self, x, long):
        if long:
            return (wang_inv(x & M64) | (wang_inv(x >> 64) << 64)) ^ self.x128
        return wang_inv(x) ^ self.x64


_C1, _C2, _C3 = 0x533f8c2151b20f97, 0x9a98567ed20c127d | 1, 0x691a9d706391077a
_C1H, _C2H, _C3H = wang(0x533f8c2151b20f97), wang(0x9a98567ed20c127d) | 1, wang(0x691a9d706391077a)


def lex64(x):
    """FRev64: xor, multiply, rotate left 31, xor (bonsai encoder.h lex_score for 64-bit k-mers)."""
    h = ((x ^ _C1) * _C2) & M64
    h = ((h << 31) | (h >> 33)) & M64
    return h ^ _C3


def lex128(x):
    """CEHasher on a 128-bit k-mer: xor/multiply/xor on each 64-bit lane, the high lane with rehashed constants."""
    lo, hi = x & M64, x >> 64
    lo = (((lo ^ _C1) * _C2) & M64) ^ _C3
    hi = (((hi ^ _C1H) * _C2H) & M64) ^ _C3H
    return lo | (hi << 64)


def entropy(codes):
    """Shannon entropy in nats of a k-mer's character codes (bonsai entropy.h CircusEnt::value)."""
    n = len(codes)
    cnt = {}
    for c in codes:
        cnt[c] = cnt.get(c, 0) + 1
    return -sum(v / n * math.log(v / n) for v in cnt.values())


def ent_score(x, ent, long):
    """--entmin score: the top 48 bits of the lexical score divided by (entropy + 1e-4); lower wins."""
    h = (lex128(x) >> 80) if long else (lex64(x) >> 16)
    return int(h / (max(ent, 0.0) + 1e-4))


def downsample_hash(x, long):
    """Salted hash deciding whether --downsample keeps a (masked) k-mer (src/d2.h)."""
    def h(v):
        v ^= 0x9e3779b97f4a7c15
        v = ((v ^ (v >> 33)) * 0xff51afd7ed558ccd) & M64
        v = ((v ^ (v >> 33)) * 0xc4ceb9fe1a85ec53) & M64
        return v ^ (v >> 33)
    if long:
        return h((x & M64) ^ h(x >> 64))
    return h(x)


def downsample_threshold(f):
    # ceil((2^64 - 1) * f) in long double; on arm64 long double is double, and
    # double(2^64 - 1) is 2^64. Values within one part in 2^53 of the boundary
    # are vanishingly rare, so double arithmetic decides the same way.
    return int(math.ceil(float(M64) * f))


# ---------------------------------------------------------------------------
# Alphabets


class Alphabet:
    def __init__(self, name, flags, groups, mul, bitpacked, cap64, cap128, aliases):
        self.name, self.flags = name, list(flags)
        self.mul, self.bitpacked = mul, bitpacked
        self.cap64, self.cap128 = cap64, cap128
        self.code = {}
        for i, g in enumerate(groups.split(",")):
            for ch in g:
                self.code[ch] = i
                self.code[ch.lower()] = i
        for a, b in aliases:
            if a not in self.code:
                self.code[a] = self.code[b]
                self.code[a.lower()] = self.code[b]
        self.letters = sorted({c for c in self.code if c.isupper()})

    def cap(self, long):
        return self.cap128 if long else self.cap64

    def encode(self, codes):
        v = 0
        if self.bitpacked:
            sh = {4: 2, 8: 3}[self.mul]
            for c in codes:
                v = (v << sh) | c
        else:
            for c in codes:
                v = v * self.mul + c
        return v


DNA = Alphabet("dna", [], "A,C,G,T", 4, True, 32, 64, [("U", "T")])
AMINO20 = Alphabet("protein20", ["--protein"], "A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y", 20, False, 14, 29,
                   [("O", "K"), ("U", "C")])
SEB14 = Alphabet("protein14", ["--protein14"], "A,C,D,EQ,FY,G,H,IV,KR,LM,N,P,ST,W", 14, False, 16, 33,
                 [("O", "K"), ("U", "C")])
SEB8 = Alphabet("protein8", ["--protein8"], "AST,C,DHN,EKQR,FWY,G,ILMV,P", 8, True, 21, 42, [("O", "K"), ("U", "C")])
SEB6 = Alphabet("protein6", ["--protein6"], "AST,CP,DHNEKQR,FWY,G,ILMV", 6, False, 24, 49, [("O", "K"), ("U", "C")])
PROTEIN_ALPHABETS = [AMINO20, SEB14, SEB8, SEB6]
# Every spelling of the protein flags that dashing2 accepts (src/options.h).
PROTEIN_ALIASES = {"protein20": ["--protein", "--protein20", "--enable-protein"], "protein14": ["--protein14"],
                   "protein8": ["--protein8"], "protein6": ["--protein6"]}
_RC = str.maketrans("ACGTacgtUu", "TGCAtgcaAa")


def dna_rc_value(v, k):
    r = 0
    for _ in range(k):
        r = (r << 2) | (3 - (v & 3))
        v >>= 2
    return r


def offsets_from_spacing(k, spacing):
    """Positions of the kept characters for a --spacing list of k - 1 gaps (0 = adjacent)."""
    offs = [0]
    for g in (spacing if spacing else [0] * (k - 1)):
        offs.append(offs[-1] + g + 1)
    return offs


def kmer_items(rec, k, alph=DNA, canon=True, spacing=None):
    """Valid k-mers of one record in position order as (value, entropy-codes) pairs.

    Characters outside the alphabet break k-mers. canon applies only to
    unspaced DNA (spaced seeds and protein are never canonicalized).
    """
    offs = offsets_from_spacing(k, spacing)
    span = offs[-1] + 1
    codes = [alph.code.get(ch, -1) for ch in rec]
    out = []
    use_canon = canon and alph is DNA and not spacing
    for i in range(len(rec) - span + 1):
        cs = [codes[i + o] for o in offs]
        if min(cs) < 0:
            continue
        v = alph.encode(cs)
        if use_canon:
            r = dna_rc_value(v, k)
            v = min(v, r)
        out.append((v, cs))
    return out


def select_window(items, m, scorefn):
    """dashing2's window order: one pick per full window of m valid items, best (score, value) wins.

    items are (value, codes); scorefn(value, codes) gives the score. With m <= 1
    every item is emitted. A record with 0 < len(items) < m emits its best item.
    Returns the list of picked values.
    """
    if m <= 1:
        return [v for v, _ in items]
    keyed = [(scorefn(v, cs), v) for v, cs in items]
    out = []
    # Monotone deque of indices for the sliding minimum.
    from collections import deque
    dq = deque()
    for j, key in enumerate(keyed):
        while dq and keyed[dq[-1]] >= key:
            dq.pop()
        dq.append(j)
        if dq[0] <= j - m:
            dq.popleft()
        if j >= m - 1:
            out.append(keyed[dq[0]][1])
    if 0 < len(keyed) < m:
        out.append(min(keyed)[1])
    return out


def window_picks(recs, k, w, alph=DNA, canon=True, spacing=None, long=False, entmin=False):
    """Per record, the raw (unmasked) items dashing2 emits for these parsing options (exact-encoding path)."""
    offs = offsets_from_spacing(k, spacing)
    span = offs[-1] + 1
    m = (max(w, span) - span + 1) if w and w > 0 else 1
    if entmin:
        def sc(v, cs):
            return ent_score(v, entropy(cs), long)
    else:
        lx = lex128 if long else lex64

        def sc(v, cs):
            return lx(v)
    return [select_window(kmer_items(r, k, alph, canon, spacing), m, sc) for r in recs]


def distinct_kmer_strings(recs, k, alph=DNA, canon=True):
    """Distinct k-mers as strings (any k); used for the rolling-hash paths."""
    out = set()
    for r in recs:
        cur = []
        for ch in r + "\0":
            c = alph.code.get(ch, -1)
            if c < 0:
                seg = cur
                cur = []
                n = len(seg)
                if n < k:
                    continue
                t = tuple(seg)
                if canon and alph is DNA:
                    rc = tuple(3 - x for x in reversed(seg))
                    for i in range(n - k + 1):
                        f = t[i:i + k]
                        b = rc[n - i - k:n - i]
                        out.add(min(f, b))
                else:
                    for i in range(n - k + 1):
                        out.add(t[i:i + k])
            else:
                cur.append(c)
    return out


def n_valid_kmers(rec, k, alph=DNA):
    n = run = 0
    for ch in rec:
        if alph.code.get(ch, -1) < 0:
            run = 0
        else:
            run += 1
            if run >= k:
                n += 1
    return n


# ---------------------------------------------------------------------------
# Files


def write_records(path, recs, names=None, fmt="fa", width=None, comment=True):
    """Writes FASTA or FASTQ, compressed by extension (.gz, .bz2, .xz, .zst)."""
    names = names or ["r%d" % i for i in range(len(recs))]
    lines = []
    for nm, r in zip(names, recs):
        if fmt == "fq":
            lines.append("@%s%s\n%s\n+\n%s\n" % (nm, " c" if comment else "", r, "I" * len(r)))
        else:
            lines.append(">%s%s\n" % (nm, " some description" if comment else ""))
            w = width or max(1, len(r))
            for j in range(0, len(r), w):
                lines.append(r[j:j + w] + "\n")
            if not r:
                pass
    data = "".join(lines).encode()
    if path.endswith(".gz"):
        with gzip.open(path, "wb") as f:
            f.write(data)
    elif path.endswith(".bz2"):
        with bz2.open(path, "wb") as f:
            f.write(data)
    elif path.endswith(".xz"):
        with lzma.open(path, "wb") as f:
            f.write(data)
    elif path.endswith(".zst"):
        raw = path[:-4]
        with open(raw, "wb") as f:
            f.write(data)
        subprocess.run(["zstd", "-q", "-f", "--rm", raw, "-o", path], check=True)
    else:
        with open(path, "wb") as f:
            f.write(data)


def have_zstd():
    return shutil.which("zstd") is not None


def read_kmerset(path, long):
    with open(path, "rb") as f:
        b = f.read()
    card = struct.unpack("<d", b[:8])[0]
    body = b[8:]
    if long:
        vals = [int.from_bytes(body[i:i + 16], "little") for i in range(0, len(body), 16)]
    else:
        vals = list(struct.unpack("<%dQ" % (len(body) // 8), body))
    return card, vals


def read_stream(path, long):
    with open(path, "rb") as f:
        b = f.read()
    if long:
        return [int.from_bytes(b[i:i + 16], "little") for i in range(0, len(b), 16)]
    return list(struct.unpack("<%dQ" % (len(b) // 8), b))


def read_byseq_stream_file(path, long):
    """Parses `sketch -G --parse-by-seq -o` output.

    Layout: uint64 n, uint32 k, uint32 w, uint32 (alphabet | canonical << 8),
    n doubles (stream lengths in 64-bit words), then the streams back to back.
    """
    with open(path, "rb") as f:
        b = f.read()
    n, k, w, dt = struct.unpack("<QIII", b[:20])
    lens = struct.unpack("<%dd" % n, b[20:20 + 8 * n])
    pos = 20 + 8 * n
    streams = []
    for ln in lens:
        ln = int(ln)
        chunk = b[pos:pos + 8 * ln]
        pos += 8 * ln
        if long:
            streams.append([int.from_bytes(chunk[i:i + 16], "little") for i in range(0, len(chunk), 16)])
        else:
            streams.append(list(struct.unpack("<%dQ" % (len(chunk) // 8), chunk)))
    return dict(n=n, k=k, w=w, dtype=dt, lens=lens, streams=streams, trailing=len(b) - pos)


def levenshtein(a, b):
    if len(a) < len(b):
        a, b = b, a
    prev = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        cur = [i] + [0] * len(b)
        for j, y in enumerate(b, 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (x != y))
        prev = cur
    return prev[-1]


def hp_collapse(xs):
    out = []
    for x in xs:
        if not out or out[-1] != x:
            out.append(x)
    return out


def files_written(cwd, before):
    return sorted(set(os.listdir(cwd)) - before)
