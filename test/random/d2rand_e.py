"""Helpers for suite E (workflow and output modes), Python standard library only.

Suite E checks the parts of dashing2 that sit around the comparison core:
nearest-neighbour and clustering outputs, binary and text output layouts,
sketch files and their round trips, `contain`, `wsketch`, `printmin` and BED
input. This module holds what those checks share: k-mer identifiers exactly as
dashing2 stores them, decoders for every file layout the suite reads, a
generator of clustered random genomes, and a small runner that keeps binary
stdout intact.

It builds on d2rand.py (k-mer counting, FASTA writing, trial machinery) and
does not change it.
"""
import math
import os
import random
import re
import struct
import subprocess
from collections import Counter

from d2rand import (add_repeats, decorate, kmer_counts, mutate, random_dna, read_fastx, write_fastx)

M64 = (1 << 64) - 1
GOLDEN = 0x9e3779b97f4a7c15
DEFAULT_XORMASK = 0x724526e320f9967d
FLT_MAX = 3.4028234663852886e+38


# ---------------------------------------------------------------------------
# K-mer identifiers as dashing2 stores them
#
# dashing2 encodes DNA k-mers in 2 bits per base (A=0, C=1, G=2, T=3, first
# base most significant), takes the numerically smaller of the k-mer and its
# reverse complement when canonicalizing (the same choice as the
# lexicographically smaller string), and stores maskfn(kmer) =
# WangHash(kmer ^ XORMASK), where XORMASK is 0 for --seed 0 (the default) and
# WangHash(seed) otherwise (src/enums.h, src/enums.cpp). With --long-kmers the
# 128-bit k-mer is masked with XORMASK2, each half is hashed, and the result
# is folded to 64 bits by fold64 (src/enums.h).

def wang64(key):
    key = ((~key) + (key << 21)) & M64
    key ^= key >> 24
    key = (key + (key << 3) + (key << 8)) & M64
    key ^= key >> 14
    key = (key + (key << 2) + (key << 4)) & M64
    key ^= key >> 28
    key = (key + (key << 31)) & M64
    return key


def xormasks(seed):
    """(XORMASK, XORMASK2) that dashing2 derives from --seed."""
    if seed == 0:
        return 0, 0
    x = wang64(seed)
    return x, x | (wang64(x) << 64)


_ENC = {"A": 0, "C": 1, "G": 2, "T": 3}


def encode(s):
    v = 0
    for c in s:
        v = (v << 2) | _ENC[c]
    return v


def kmer_id(s, seed=0, use128=False):
    """The 64-bit identifier dashing2 stores for k-mer string s (already canonical if wanted)."""
    x = encode(s)
    m64, m128 = xormasks(seed)
    if not use128:
        return wang64(x ^ m64)
    x ^= m128
    lo = wang64(x & M64)
    hi = wang64(x >> 64)
    return lo ^ wang64(hi ^ GOLDEN)


def kmer_id128_raw(s, seed=0):
    """The unfolded 128-bit value dashing2 stores in .kmerset128 files for --long-kmers (low half first)."""
    x = encode(s) ^ xormasks(seed)[1]
    return wang64(x & M64) | (wang64(x >> 64) << 64)


def id_counts(counts, seed=0, use128=False):
    """Maps a Counter of k-mer strings to a Counter of dashing2 identifiers."""
    out = Counter()
    for s, n in counts.items():
        out[kmer_id(s, seed, use128)] += n
    return out


def u128s(path, off=0):
    with open(path, "rb") as f:
        b = f.read()
    v = struct.unpack_from("<%dQ" % ((len(b) - off) // 8), b, off)
    return [v[i] | (v[i + 1] << 64) for i in range(0, len(v) - 1, 2)]


def kmer_sequence(recs, k, canon=True):
    """Every k-mer of every record in order (A/C/G/T only, any other character breaks), canonicalized if asked.

    Returns one list per record, as dashing2 --parse-by-seq sees them.
    """
    from d2rand import revcomp
    out = []
    for r in recs:
        seq = []
        for seg in re.split(r"[^ACGT]+", r.upper()):
            n = len(seg)
            for i in range(n - k + 1):
                f = seg[i:i + k]
                if canon:
                    b = revcomp(f)
                    f = f if f <= b else b
                seq.append(f)
        out.append(seq)
    return out


# ---------------------------------------------------------------------------
# float32 helpers (dashing2 computes and prints results as float32)

def f32(x):
    """Rounds a Python float to the nearest float32."""
    if math.isinf(x) or math.isnan(x):
        return x
    if abs(x) > FLT_MAX:
        return math.copysign(float("inf"), x)
    return struct.unpack("<f", struct.pack("<f", x))[0]


def same32(a, b):
    """Equality of two values that should be the same float32 (text printed with full precision)."""
    if a is None or b is None:
        return False
    if (math.isinf(a) or a >= FLT_MAX) and (math.isinf(b) or b >= FLT_MAX):
        return True
    return f32(a) == f32(b)


def near32(a, b, rel=1.2e-7):
    """Equality up to printing with 8 significant digits ({:0.8g}), which cannot always round-trip a float32."""
    if a is None or b is None:
        return False
    if (math.isinf(a) or a >= FLT_MAX) and (math.isinf(b) or b >= FLT_MAX):
        return True
    return abs(a - b) <= rel * max(abs(a), abs(b)) + 1e-30


def between32(lo, hi):
    """A float32-representable value strictly between two distinct float32 values lo < hi, or None."""
    t = f32((lo + hi) / 2.0)
    if lo < t < hi:
        return t
    return None


# ---------------------------------------------------------------------------
# Running dashing2 with binary output intact

class Run:
    def __init__(self, rc, out, err, cmd):
        self.rc, self.out, self.err, self.cmd = rc, out, err, cmd

    @property
    def text(self):
        return self.out.decode("utf-8", "replace")

    def brief(self):
        tail = self.err.strip().splitlines()[-3:]
        return "exit %d (%s): %s" % (self.rc, " ".join(self.cmd[1:])[:300], " | ".join(tail)[-400:])


def run(binary, args, cwd, timeout=600, env=None):
    """Runs dashing2 with args in cwd; stdout is kept as bytes."""
    cmd = [binary] + [str(a) for a in args]
    e = None
    if env:
        e = dict(os.environ)
        e.update(env)
    try:
        p = subprocess.run(cmd, cwd=cwd, capture_output=True, timeout=timeout, env=e)
    except subprocess.TimeoutExpired:
        return Run(-999, b"", "timeout after %ds" % timeout, cmd)
    return Run(p.returncode, p.stdout, p.stderr.decode("utf-8", "replace"), cmd)


def cmdline(r, cwd=None):
    """A shell command reproducing a call (for messages)."""
    import shlex
    s = " ".join(shlex.quote(c) for c in ["dashing2"] + r.cmd[1:])
    return ("(cd %s && %s)" % (shlex.quote(cwd), s)) if cwd else s


def snapshot(wd):
    return set(os.listdir(wd))


def remove_new(wd, before):
    """Deletes files that appeared in wd since snapshot `before` (dashing2 writes k-mer files beside inputs)."""
    import shutil
    for fn in os.listdir(wd):
        if fn in before:
            continue
        p = os.path.join(wd, fn)
        if os.path.isdir(p):
            shutil.rmtree(p, ignore_errors=True)
        else:
            os.remove(p)


# ---------------------------------------------------------------------------
# Decoders for dashing2 outputs

def parse_text_matrix(text):
    """Parses the human-readable full, square and panel outputs.

    Returns (kind, rows, cols, vals) with kind 'sym', 'square' or 'panel' and
    vals mapping (row index, column index) to float.
    """
    from d2rand import parse_matrix
    M = parse_matrix(text)
    return M.kind, M.rows, M.cols, M.vals


def parse_phylip(text):
    """Parses --phylip output: a count line, then for row i the values for columns i+1..n-1."""
    lines = [l for l in text.splitlines() if l.strip() and not l.startswith("#")]
    n = int(lines[0].strip())
    rows, vals = [], {}
    for r, line in enumerate(lines[1:]):
        name = line[:9].strip() if "\t" not in line[:10] else line.split("\t")[0].strip()
        rest = line.split("\t")[1:]
        rows.append(name)
        for t, v in enumerate(rest):
            vals[(r, r + 1 + t)] = float(v)
    return n, rows, vals


def floats(b, off=0, n=None):
    if n is None:
        n = (len(b) - off) // 4
    return list(struct.unpack_from("<%df" % n, b, off))


def parse_condensed(b, n):
    """Binary upper-triangular output: rows 0..n-1, columns i+1..n-1, float32, row-major."""
    v = floats(b)
    out, t = {}, 0
    for i in range(n):
        for j in range(i + 1, n):
            out[(i, j)] = v[t] if t < len(v) else None
            t += 1
    return out, len(v)


def parse_dense(b, nr, nc):
    v = floats(b)
    return {(i, j): (v[i * nc + j] if i * nc + j < len(v) else None) for i in range(nr) for j in range(nc)}, len(v)


def parse_nn_text(text):
    """Top-k and threshold text output: one line per input, 'name<TAB>neighbor:value...'."""
    rows = []
    for line in text.splitlines():
        if line.startswith("#") or not line.strip():
            continue
        t = line.split("\t")
        nb = []
        for x in t[1:]:
            name, _, v = x.rpartition(":")
            nb.append((name, float(v)))
        rows.append((t[0], nb))
    return rows


def parse_csr(b, idsize=4, dsize=4):
    """CSR layout of --topk/--similarity-threshold binary output (python/parse.py parse_knn)."""
    if len(b) < 16:
        raise ValueError("CSR output has %d bytes" % len(b))
    nids, nnz = struct.unpack_from("<QQ", b, 0)
    off = 16
    indptr = list(struct.unpack_from("<%dQ" % (nids + 1), b, off))
    off += 8 * (nids + 1)
    idx = list(struct.unpack_from("<%d%s" % (nnz, "I" if idsize == 4 else "Q"), b, off))
    off += idsize * nnz
    data = list(struct.unpack_from("<%df" % nnz, b, off))
    off += 4 * nnz
    rows = [[(idx[p], data[p]) for p in range(indptr[i], indptr[i + 1])] for i in range(nids)]
    return rows, off, len(b)


def parse_greedy_text(text):
    """--greedy text output: 'Cluster-<i><TAB>rep:id<TAB>member:id...' after a '#Clustering' header."""
    clusters = []
    header = None
    for line in text.splitlines():
        if line.startswith("#"):
            header = line
            continue
        if not line.strip():
            continue
        t = line.split("\t")
        members = []
        for x in t[1:]:
            name, _, i = x.rpartition(":")
            members.append((name, int(i)))
        clusters.append((t[0], members))
    return header, clusters


def parse_greedy_bin(b):
    """--greedy binary output: (nclusters, nitems) u64, indptr u64[nclusters + 1], ids u32[nitems] (rep first)."""
    nc, nnz = struct.unpack_from("<QQ", b, 0)
    indptr = list(struct.unpack_from("<%dQ" % (nc + 1), b, 16))
    off = 16 + 8 * (nc + 1)
    ids = list(struct.unpack_from("<%dI" % nnz, b, off))
    off += 4 * nnz
    return [ids[indptr[i]:indptr[i + 1]] for i in range(nc)], off, len(b)


def parse_stacked(path):
    """Stacked sketches written by -o: (n, sketch size) u64, n cardinalities f64, then the registers."""
    with open(path, "rb") as f:
        b = f.read()
    n, s = struct.unpack_from("<QQ", b, 0)
    cards = list(struct.unpack_from("<%dd" % n, b, 16))
    return n, s, cards, b[16 + 8 * n:]


def parse_names(path):
    """<out>.names.txt: '#Name<TAB>Cardinality' header, then name, cardinality and optionally a count file."""
    out = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            t = line.rstrip("\n").split("\t")
            out.append((t[0], float(t[1]) if len(t) > 1 and t[1] else None, t[2] if len(t) > 2 else None))
    return out


def parse_kmerdb(path):
    """<out>.kmer64 written by sketch --save-kmers -o: a 24-byte header, then sketch-size ids per input.

    Header: u32 alphabet | canonical << 8 | long-kmers << 9, u32 sketch size,
    u32 k, u32 w, u64 seed.
    """
    with open(path, "rb") as f:
        b = f.read()
    dtype, s, k, w, seed = struct.unpack_from("<4IQ", b, 0)
    n = (len(b) - 24) // 8
    v = list(struct.unpack_from("<%dQ" % n, b, 24))
    rows = [v[i * s:(i + 1) * s] for i in range(n // s)] if s else []
    return dict(alphabet=dtype & 0xff, canon=bool(dtype >> 8 & 1), long=bool(dtype >> 9 & 1), sketchsize=s, k=k,
                w=w, seed=seed, rows=rows, extra=(len(b) - 24) % (8 * s) if s else 0)


def parse_contain_text(text):
    """contain text output: '##References:' line, then per query 'name<TAB>cov%:depth...'."""
    refs, rows = None, []
    for line in text.splitlines():
        if line.startswith("##References:"):
            refs = line.split("\t")[1:]
            continue
        if line.startswith("#") or not line.strip():
            continue
        t = line.split("\t")
        vals = []
        for x in t[1:]:
            c, _, d = x.partition("%:")
            vals.append((float(c), float(d)))
        rows.append((t[0], vals))
    return refs, rows


def parse_contain_bin(b):
    """contain -b output: (nref, nquery) u64, coverage f32[nquery][nref], depth f32[nquery][nref]."""
    nref, nq = struct.unpack_from("<QQ", b, 0)
    v = floats(b, 16)
    cov = [v[q * nref:(q + 1) * nref] for q in range(nq)]
    dep = [v[nref * nq + q * nref: nref * nq + (q + 1) * nref] for q in range(nq)]
    return nref, nq, cov, dep, len(v)


def parse_mmerseq_stacked(path):
    """Minimizer sequences from sketch -G --parse-by-seq -o: n u64, k, w, dtype u32, n lengths f64, u64 items."""
    with open(path, "rb") as f:
        b = f.read()
    n = struct.unpack_from("<Q", b, 0)[0]
    k, w, dtype = struct.unpack_from("<3I", b, 8)
    lens = list(struct.unpack_from("<%dd" % n, b, 20))
    off = 20 + 8 * n
    tot = int(sum(lens))
    items = list(struct.unpack_from("<%dQ" % min(tot, (len(b) - off) // 8), b, off))
    seqs, p = [], 0
    for L in lens:
        seqs.append(items[p:p + int(L)])
        p += int(L)
    return dict(n=n, k=k, w=w, dtype=dtype, lens=lens, seqs=seqs, size=len(b), expected=off + 8 * tot)


def u64s(path, off=0):
    with open(path, "rb") as f:
        b = f.read()
    return list(struct.unpack_from("<%dQ" % ((len(b) - off) // 8), b, off))


def f64s(path, off=0):
    with open(path, "rb") as f:
        b = f.read()
    return list(struct.unpack_from("<%dd" % ((len(b) - off) // 8), b, off))


# ---------------------------------------------------------------------------
# Clustered random genomes

SNP_LEVELS = [0.0, 0.002, 0.005, 0.01, 0.02, 0.04, 0.08, 0.15]


def clustered_inputs(rng, wd, n, k, minlen=400, maxlen=4000, features=True, prefix="g", dup_prob=0.08,
                     slice_prob=0.2, empty_prob=0.0, fmt_mix=True):
    """Writes n FASTA/FASTQ files drawn from a few random root genomes and returns their names and records.

    Members of a cluster are mutated copies of its root (SNP rate drawn from
    SNP_LEVELS, small indels), sometimes a slice of it (so containment is
    asymmetric), sometimes with repeats (k-mer counts above 1). A file is
    occasionally an exact copy of an earlier one, which creates ties.
    """
    ncl = rng.randint(1, max(1, min(5, n // 2 + 1)))
    roots = [random_dna(rng, rng.randint(minlen, maxlen)) for _ in range(ncl)]
    names, recs_all, meta = [], [], []
    for i in range(n):
        if i and rng.random() < dup_prob:
            j = rng.randrange(i)
            recs = list(recs_all[j])
            m = dict(meta[j]); m["copy_of"] = j
        else:
            c = rng.randrange(ncl)
            snp = rng.choice(SNP_LEVELS)
            s = mutate(rng, roots[c], snp, rng.choice([0, 0, 0.001, 0.003]))
            sl = False
            if rng.random() < slice_prob and len(s) > 4 * k:
                a = rng.randrange(len(s) // 2)
                s = s[a:a + rng.randint(len(s) // 5, len(s))]
                sl = True
            if features and rng.random() < 0.4:
                s = add_repeats(rng, s, rng.randint(1, 4))
            p = {}
            if features:
                p = dict(nrec=rng.randint(1, 4), short_recs=rng.randint(0, 1), rc_prob=rng.choice([0, 0.5]),
                         n_runs=rng.randint(0, 2), lower_runs=rng.randint(0, 1))
            recs = decorate(rng, s, k, p) if p else [s]
            if rng.random() < empty_prob:
                recs = [s[:max(1, k - 1)]]
            m = dict(cluster=c, snp=snp, slice=sl)
        fmt = "fq" if fmt_mix and rng.random() < 0.15 else "fa"
        gz = fmt_mix and rng.random() < 0.1
        name = "%s%02d.%s%s" % (prefix, i, fmt, ".gz" if gz else "")
        write_fastx(os.path.join(wd, name), recs, rng, fmt, rng.choice([None, 60]))
        names.append(name)
        recs_all.append(recs)
        meta.append(m)
    return names, recs_all, meta


def exact_counts(wd, names, k, canon=True):
    return [kmer_counts(read_fastx(os.path.join(wd, f)), k, canon) for f in names]


def write_list(wd, name, items):
    with open(os.path.join(wd, name), "w") as f:
        f.write("\n".join(items) + "\n")
    return name


# ---------------------------------------------------------------------------
# Small numeric helpers

def jaccard(a, b):
    """Set Jaccard of two collections (1 for two empty sets)."""
    a, b = set(a), set(b)
    u = len(a | b)
    return len(a & b) / u if u else 1.0


def weighted_jaccard(wa, wb):
    keys = set(wa) | set(wb)
    num = sum(min(wa.get(x, 0), wb.get(x, 0)) for x in keys)
    den = sum(max(wa.get(x, 0), wb.get(x, 0)) for x in keys)
    return num / den if den else 1.0


def prob_jaccard_w(wa, wb):
    """Probability Jaccard of two weight maps (Moulton and Jiang), quadratic but fine for small maps."""
    sa, sb = sum(wa.values()), sum(wb.values())
    if not sa or not sb:
        return float(sa == sb)
    keys = set(wa) | set(wb)
    tot = 0.0
    for x in set(wa) & set(wb):
        px, qx = wa[x] / sa, wb[x] / sb
        s = 0.0
        for y in keys:
            s += max(wa.get(y, 0) / sa / px, wb.get(y, 0) / sb / qx)
        tot += 1.0 / s
    return tot


def eq_z(neq, m, J):
    """z-score of an equal-register count against probability J, with a one-register resolution floor."""
    se = math.sqrt(max(J * (1 - J), 0.0) / m + 1.0 / m ** 2)
    return (neq / m - J) / se
