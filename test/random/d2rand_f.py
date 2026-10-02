"""Helpers for suite F (documentation audit), Python standard library only.

Suite F checks that what dashing2's usage messages and README.md say is what
the binary does. This module holds the pieces that are specific to that
suite: running a subcommand and capturing everything it wrote, parsers for the
non-matrix outputs (neighbour lists, greedy clusters, contain tables, stacked
sketch files), exact oracles for alphabets and seeds that d2rand does not
cover (protein k-mers, spaced seeds, minimizer sequences), the parser's own
option tables read from src/options.h, and the registry of documentation
bugs (DOC_BUGS) that the suite reports as expected failures.
"""
import os
import re
import struct
import subprocess

from d2rand import D2Error, parse_matrix, revcomp

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "..", "src")

AMINO20 = "ACDEFGHIKLMNPQRSTVWY"


# ---------------------------------------------------------------------------
# Running dashing2

class Run:
    """Result of one dashing2 call: exit status, stdout, stderr (all text)."""

    def __init__(self, cmd, p):
        self.cmd = cmd
        self.rc = p.returncode
        self.out = p.stdout
        self.err = p.stderr

    @property
    def crashed(self):
        # Signals (segfault, abort) show up as negative return codes.
        return self.rc < 0

    def brief(self):
        tail = self.err.strip().splitlines()[-2:] if self.err.strip() else []
        return "exit %d%s" % (self.rc, (": " + " | ".join(tail))[:300] if tail else "")


def run(d2, args, cwd, stdin=None, timeout=300, binary=False):
    """Runs dashing2 with args in cwd. With binary=True stdout is kept as bytes."""
    cmd = [d2.binary] + [str(a) for a in args]
    p = subprocess.run(cmd, cwd=cwd, input=stdin, capture_output=True, text=not binary, timeout=timeout)
    if binary:
        p.stderr = p.stderr.decode("utf-8", "replace")
    return Run(cmd, p)


def matrix(r):
    """Parses the matrix printed by a successful run; raises D2Error otherwise."""
    if r.rc != 0:
        raise D2Error(r.brief())
    return parse_matrix(r.out)


def same_matrix(m1, m2, rel=0.0):
    """True if two parsed matrices have the same rows, columns and values (rel tolerance)."""
    if m1.rows != m2.rows or m1.cols != m2.cols or set(m1.vals) != set(m2.vals):
        return False
    for key, v in m1.vals.items():
        w = m2.vals[key]
        if v != w and abs(v - w) > rel * max(abs(v), abs(w)):
            return False
    return True


def header_options(text):
    """The `#Dashing2Options:` header as a dict (k, sketchsize, sketchtype, ...)."""
    m = re.search(r"^#Dashing2Options: Dashing2Options;(.*)$", text, re.M)
    if not m:
        return {}
    d = {}
    for item in m.group(1).split(";"):
        if ":" in item:
            a, b = item.split(":", 1)
            d[a] = b
        elif item:
            d[item] = True
    return d


def cmd_line(d2, args):
    return " ".join([os.path.basename(d2.binary)] + [str(a) for a in args])


# ---------------------------------------------------------------------------
# Output parsers

def parse_lists(text):
    """Neighbour lists from --topk / --similarity-threshold: {row: [(name, value), ...]}."""
    out = {}
    for line in text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        t = line.split("\t")
        lst = []
        for item in t[1:]:
            if not item:
                continue
            name, val = item.rsplit(":", 1)
            lst.append((name, float(val)))
        out[t[0].strip()] = lst
    return out


def parse_greedy(text):
    """Clusters from --greedy text output: list of lists of member names."""
    cl = []
    for line in text.splitlines():
        if line.startswith("Cluster-"):
            cl.append([x.rsplit(":", 1)[0] for x in line.split("\t")[1:] if x])
    return cl


def parse_contain(text):
    """`dashing2 contain` table: (reference names, {input: [(coverage percent, depth), ...]})."""
    refs, rows = None, {}
    for line in text.splitlines():
        if line.startswith("##References:"):
            refs = line.split("\t")[1:]
        elif line.startswith("#") or not line.strip():
            continue
        else:
            t = line.split("\t")
            vals = []
            for item in t[1:]:
                cov, dep = item.split(":")
                vals.append((float(cov.rstrip("%")), float(dep)))
            rows[t[0]] = vals
    return refs, rows


def read_stacked(path, regbytes=8):
    """A stacked sketch file (sketch -o): (n, sketch size, cardinalities, register bytes)."""
    with open(path, "rb") as f:
        data = f.read()
    n, S = struct.unpack_from("<QQ", data, 0)
    cards = list(struct.unpack_from("<%dd" % n, data, 16))
    return n, S, cards, data[16 + 8 * n:]


def condensed(m, n):
    """Values of a symmetric matrix in condensed (upper triangle, row-major) order."""
    return [m.get(i, j) for i in range(n) for j in range(i + 1, n)]


def floats(data, offset=0, count=None):
    count = (len(data) - offset) // 4 if count is None else count
    return list(struct.unpack_from("<%df" % count, data, offset))


# ---------------------------------------------------------------------------
# Exact oracles

def kmer_set(recs, k, canon=True):
    """Distinct DNA k-mers (canonical: lexicographic minimum with the reverse complement)."""
    from d2rand import kmer_counts
    return set(kmer_counts(recs, k, canon))


def protein_kmers(recs, k):
    """Distinct k-mers of records over the 20 standard amino acids (records hold only those)."""
    s = set()
    for r in recs:
        for i in range(len(r) - k + 1):
            s.add(r[i:i + k])
    return s


def spacing_offsets(spec):
    """Positions kept by a --spacing spec: gaps between consecutive kept positions, with AxN run-length terms."""
    gaps = []
    for t in spec.split(","):
        if "x" in t:
            a, b = t.split("x")
            gaps += [int(a)] * int(b)
        else:
            gaps.append(int(t))
    offs = [0]
    for g in gaps:
        offs.append(offs[-1] + 1 + g)
    return offs


def spaced_kmers(recs, offs):
    """Distinct strand-specific spaced k-mers (A/C/G/T only; other characters break)."""
    span = offs[-1] + 1
    s = set()
    for r in recs:
        for seg in re.split(r"[^ACGT]+", r.upper()):
            for i in range(len(seg) - span + 1):
                s.add("".join(seg[i + o] for o in offs))
    return s


def canonical_kmer_sequence(seq, k, canon=True):
    """The k-mers of one record in order (each canonicalized), as `printmin` prints them without -w."""
    out = []
    for i in range(len(seq) - k + 1):
        f = seq[i:i + k]
        out.append(min(f, revcomp(f)) if canon else f)
    return out


def dedupe_runs(xs):
    out = []
    for x in xs:
        if not out or out[-1] != x:
            out.append(x)
    return out


def random_protein(rng, n):
    return "".join(rng.choice(AMINO20) for _ in range(n))


def mutate_protein(rng, s, p):
    return "".join(c if rng.random() >= p else rng.choice(AMINO20) for c in s)


# ---------------------------------------------------------------------------
# The option parser, read from the source tree

def parser_long_options(src=SRC):
    """Long options accepted by sketch/cmp, from SHARED_OPTS in src/options.h (None if absent).

    Entries that are commented out (/* ... */) are left out, as the compiler
    does.
    """
    path = os.path.join(src, "options.h")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        text = f.read()
    m = re.search(r"#define SHARED_OPTS(.*?)\n\n", text, re.S)
    if not m:
        return None
    body = re.sub(r"/\*.*?\*/", "", m.group(1), flags=re.S)
    names = set(re.findall(r'LO_(?:ARG|NO|FLAG)\("([^"]+)"', body))
    names |= set(re.findall(r'\{"([^"]+)",\s*(?:no|required)_argument', body))
    return names


def help_long_options(text):
    """Long options named anywhere in a usage text (tokens of the form --name)."""
    return set(re.findall(r"--([A-Za-z0-9][A-Za-z0-9-]*)", text))


# ---------------------------------------------------------------------------
# Documentation bugs found by this suite, by id. Each entry: what the text says,
# what the binary does (with the numbers seen on the fixed build) and the
# likely cause. Checks that fail on one of these report XFAIL[id]; --strict
# turns them into failures.

DOC_BUGS = {
    "help-k-limits": "help: 'If k is greater than this limit (31 for DNA, ... 22 for --protein8 ...)'; the direct "
                     "limits and defaults are 32 for DNA and 21 for --protein8 (bonsai rhtraits.h nper64), and "
                     "--protein14 (16) is not listed (options.h:484). Previously noted in bughunt inputs/findings.md.",
    "default-not-phylip": "help: '1. Upper Triangular PHYLIP (default)'; the default output is the tab-separated "
                          "'#Dashing2 Symmetric pairwise' matrix; PHYLIP comes only from the undocumented --phylip "
                          "(options.h:588, sketch_main.cpp:41, emitrect.cpp:139).",
    "header-sketchtype-label": "the long forms --multiset/--bagminhash/--bmh/--prob/--pminhash/--probminhash are "
                               "LO_FLAGs that set only sketch_space, so cmp's #Dashing2Options header says "
                               "sketchtype:onepermsetsketch while the values are BagMinHash/ProbMinHash ones (-B/-P "
                               "print the right label); -J prints sketchtype:fullyunknown (options.h:88-114, "
                               "d2.cpp:20-25).",
    "topk-large-k": "help: '--topk <arg> ... If <arg> is greater than N - 1, pairwise distances are instead "
                    "emitted.'; a neighbour list with all N-1 neighbours is printed instead. Previously noted in "
                    "bughunt workflow/findings.md.",
    "sparse-binary-format": "help: 'For top-k filtered, this emits a matrix of min(k, |N|) x |N| of IDs and "
                            "distances' and CSR '64-bit indptr, 32-bit indices, 32-bit floats'; both top-k and "
                            "threshold output are CSR preceded by an undocumented 16-byte [nrows, nnz] header. "
                            "Previously noted in bughunt workflow/findings.md.",
    "greedy-binary-ids": "help: '--greedy ... machine-readable ... followed by nsets 64-bit integers'; member ids "
                         "are 32-bit. Previously noted in bughunt workflow/findings.md.",
    "presketched-mash-k": "cmp --presketched loads a stacked file that does not record k, and --mash-distance then "
                          "uses the default k (32) unless -k is repeated; the help ('use this flag and pass in a "
                          "single positional argument') does not say so: k=15 sketches give 0.0024958584 instead of "
                          "0.005324498 (cmp_main.cpp:265, 386).",
    "save-kmers-names-path": "help: 'names will be written to <arg>.kmer.names.txt'; they go to "
                             "<arg>.kmer64.names.txt. Previously noted in bughunt workflow/findings.md.",
    "save-kmers-header-size": "help: '-s/--save-kmers ... This has a 16-byte header'; the header is 24 bytes (alphabet, "
                              "sketch size, k, w, then a 64-bit seed), as contain_main.cpp:192 reads it.",
    "save-kmers-no-file": "help: '-s/--save-kmers: ... puts the k-mers saved into .kmer files' and '-N ... into "
                          ".kmercounts.f64 files'; for the default one-permutation sketch and --full no per-input "
                          "file is written (with or without --cache) unless -o is given, while -B/--prob always "
                          "write them (fastxsketch.cpp:626-629 keep ids only for the stacked file).",
    "seq-header": "help (-G/--seq): 'header: [uint64_t nitems, uint32_t k, uint32_t w]'; the file has a fourth "
                  "field, a uint32 alphabet/flags word, so the header is 20 bytes (printminmain.cpp:31-41).",
    "fastcmp-presets-need-full": "help: '--setsketch-ab ... only supported for the SetSketch' and calls the default "
                                 "'SetSketch (one-permutation)'; --fastcmp-bytes/-shorts/-words and --setsketch-ab "
                                 "abort ('Sketch compressed is only available for FullSetSketch') unless --full is "
                                 "also given (cmp_main.h:117).",
    "seqs-in-ram-ignored": "--seqs-in-ram is accepted but has no effect: fastxsketch.h:19 declares "
                           "`static bool seqs_in_memory` in a header, so options.h:123 sets the copy in "
                           "sketch_main.cpp/cmp_main.cpp while fastxsketchbyseq.cpp:173 reads its own copy (the "
                           "debug log still says it is swapping to RAM because the flag is off). The help also says "
                           "--parse-by-seq spills to $TMPDIR, but inputs under 2e9 bases are always kept in RAM.",
    "readme-use7-no-output": "README Use 7: 'dashing2 sketch ... --set --topk 25 -o input_sequence_set.topk.tsv' writes "
                             "the stacked sketches (binary) to the .tsv and no top-k table at all; --cmpout is needed.",
    "readme-canon-default": "README: 'Canonicalization is off by default.'; DNA k-mers are canonical by default "
                            "(help item 6, header ';canon'): a sequence and its reverse complement give similarity 1.",
    "wsketch-U-1d": "wsketch help: '-U: Read 32-bit data weights'; with one or two paths the weights are read as "
                    "float64 (total weight 0 instead of 4952 for 200 weights), because wmh_from_file has no branch "
                    "for -U (wsketch.cpp:202-208). Listed as a follow-up in the fix plan; still reproduces.",
    "wsketch-tw-garbage": "wsketch writes '<prefix>.sampled.tw.txt' ending in one garbage character instead of "
                          "';d;L': `';' + 'd' + ';' + 'L'` adds chars (wsketch.cpp:365).",
    "wsketch-example-order": "wsketch help examples 3-5 pass the data/weights file before the indices "
                             "('... g1.fastq.k31.data64 fq.fastq.k31.indices64 ...', '... - indices64 indptr64'); "
                             "the usage line and code take indices first, so the '-' example aborts (mmap of '-') "
                             "and the others silently sketch the weights as ids (wsketch.cpp:288-292).",
    "wsketch-example-kmerset-header": "wsketch help examples 1-2 sketch g1.fastq.k31.kmerset64 directly, but .kmerset64 "
                                      "files start with an 8-byte cardinality, which wsketch reads as one more id "
                                      "(2973 ids for 2972 k-mers; with counts the vectors differ in length).",
    "short-opts-differ": "-C (no-canon) is in sketch's getopt string but not cmp's, so `cmp -C` prints the usage "
                         "and exits (sketch_main.cpp:63 vs cmp_main.cpp:258); neither documents it.",
    "undocumented-options": "options accepted by the sketch/cmp parser (options.h SHARED_OPTS) that the usage text "
                            "never names.",
}

KNOWN = {k: v.split(";")[0][:110] for k, v in DOC_BUGS.items()}
