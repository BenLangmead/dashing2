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
README = os.path.join(HERE, "..", "..", "README.md")

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


def readme_text(path=README):
    """README.md of the source tree, or None if it is absent."""
    try:
        with open(path) as f:
            return f.read()
    except OSError:
        return None


def help_long_options(text):
    """Long options named anywhere in a usage text (tokens of the form --name)."""
    return set(re.findall(r"--([A-Za-z0-9][A-Za-z0-9-]*)", text))


# ---------------------------------------------------------------------------
# Documentation bugs found by this suite, by id. Each entry: what the text says,
# what the binary does and the likely cause. Checks that fail on one of these
# report XFAIL[id]; --strict turns them into failures. Every issue this suite
# found is fixed in the tested code (the ids are in d2rand.FIXED_IDS, and
# README.md lists them), so the table is empty.

DOC_BUGS = {}

KNOWN = {k: v.split(";")[0][:110] for k, v in DOC_BUGS.items()}
