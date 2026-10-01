"""Shared helpers for the randomized dashing2 test suites (Python standard library only).

The suites in this directory generate small random FASTA/FASTQ inputs, run
dashing2 on them and compare its output with exact values computed either by
KMC3 (k-mer counting and kmc_tools set operations) or by the pure Python
k-mer counter in this module. Every trial is derived from a master seed and a
trial index, so any failure can be replayed with --seed and --trial.
"""
import argparse
import gzip
import math
import os
import random
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
_COMP = str.maketrans("ACGTacgt", "TGCAtgca")
_SPLIT = re.compile(r"[^ACGT]+")
# Characters that dashing2 and KMC both treat as k-mer breaks. U is left out on
# purpose: dashing2 reads it as T (RNA alias) while KMC treats it as invalid.
BREAK_CHARS = "NNNNnRYKMSWBDHV-."


def revcomp(s):
    return s.translate(_COMP)[::-1]


# ---------------------------------------------------------------------------
# Random inputs

def random_dna(rng, n):
    return "".join(rng.choices("ACGT", k=n))


def mutate(rng, seq, snp=0.0, indel=0.0, max_indel=4):
    """Applies independent SNPs and small insertions/deletions to seq."""
    if snp <= 0 and indel <= 0:
        return seq
    out = []
    i, n = 0, len(seq)
    total = snp + indel
    while i < n:
        # Jump to the next event position with a geometric draw.
        if total <= 0:
            out.append(seq[i:]); break
        gap = int(math.log(1.0 - rng.random()) / math.log(1.0 - total)) if total < 1 else 0
        if i + gap >= n:
            out.append(seq[i:]); break
        out.append(seq[i:i + gap]); i += gap
        if rng.random() < snp / total:
            out.append(rng.choice([c for c in "ACGT" if c != seq[i].upper()])); i += 1
        elif rng.random() < 0.5:
            out.append(random_dna(rng, rng.randint(1, max_indel)))
        else:
            i += rng.randint(1, max_indel)
    return "".join(out)


def add_repeats(rng, seq, nrep, unit=(5, 60), copies=(2, 8), rc_prob=0.3):
    """Inserts tandem repeats and dispersed copies so that k-mer counts exceed 1."""
    for _ in range(nrep):
        if not seq:
            break
        if rng.random() < 0.5:
            u = random_dna(rng, rng.randint(*unit))
            ins = u * rng.randint(*copies)
        else:
            a = rng.randrange(len(seq))
            ins = seq[a:a + rng.randint(*unit) * 3]
            ins = ins * rng.randint(1, 3)
        if rng.random() < rc_prob:
            ins = revcomp(ins)
        p = rng.randrange(len(seq) + 1)
        seq = seq[:p] + ins + seq[p:]
    return seq


def add_homopolymers(rng, seq, n, k):
    """Inserts single-base runs around length k (all-A and all-T k-mers encode as 0 and all ones)."""
    for _ in range(n):
        ln = max(1, k + rng.choice([-1, 0, 0, 1, 2, k]))
        p = rng.randrange(len(seq) + 1)
        seq = seq[:p] + rng.choice("ACGT") * ln + seq[p:]
    return seq


def decorate(rng, seq, k, p):
    """Splits seq into records and applies the features chosen in p.

    p keys (all optional): nrec, short_recs, rc_prob, n_runs, lower_runs.
    Returns the list of record strings.
    """
    nrec = max(1, p.get("nrec", 1))
    cuts = sorted(rng.sample(range(1, len(seq)), min(nrec - 1, max(0, len(seq) - 1)))) if len(seq) > 1 else []
    recs, last = [], 0
    for c in cuts + [len(seq)]:
        if c > last:
            recs.append(seq[last:c])
        last = c
    for _ in range(p.get("short_recs", 0)):
        if k > 1 and seq:
            ln = rng.randint(1, k - 1)
            a = rng.randrange(max(1, len(seq) - ln))
            recs.insert(rng.randrange(len(recs) + 1), seq[a:a + ln])
    out = []
    for r in recs:
        if rng.random() < p.get("rc_prob", 0.0):
            r = revcomp(r)
        r = list(r)
        for _ in range(p.get("n_runs", 0)):
            if not r:
                break
            a = rng.randrange(len(r)); ln = rng.randint(1, 5)
            ch = rng.choice(BREAK_CHARS)
            for j in range(a, min(len(r), a + ln)):
                r[j] = ch
        for _ in range(p.get("lower_runs", 0)):
            if not r:
                break
            a = rng.randrange(len(r)); ln = rng.randint(1, 3 * k + 10)
            for j in range(a, min(len(r), a + ln)):
                r[j] = r[j].lower()
        out.append("".join(r))
    return [r for r in out if r]


def write_fastx(path, recs, rng, fmt="fa", width=None):
    """Writes records as FASTA (wrapped at width, or one line) or FASTQ, gzipped if path ends in .gz."""
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "wt") as f:
        for i, r in enumerate(recs):
            if fmt == "fq":
                f.write("@r%d\n%s\n+\n%s\n" % (i, r, "I" * len(r)))
            else:
                f.write(">r%d some description\n" % i)
                w = width or len(r)
                for j in range(0, len(r), w):
                    f.write(r[j:j + w] + "\n")


def read_fastx(path):
    """Reads back the sequences written by write_fastx."""
    op = gzip.open if path.endswith(".gz") else open
    recs, cur = [], []
    with op(path, "rt") as f:
        lines = f.read().split("\n")
    if lines and lines[0].startswith("@"):
        for i in range(1, len(lines), 4):
            recs.append(lines[i])
        return recs
    for ln in lines:
        if ln.startswith(">"):
            if cur:
                recs.append("".join(cur))
            cur = []
        elif ln:
            cur.append(ln.strip())
    if cur:
        recs.append("".join(cur))
    return recs


# ---------------------------------------------------------------------------
# Exact oracle

def kmer_counts(recs, k, canon=True):
    """Counts k-mers over A/C/G/T (either case); any other character breaks k-mers.

    Canonical k-mers are the lexicographically smaller of the k-mer and its
    reverse complement, which is also the form KMC prints.
    """
    c = Counter()
    for r in recs:
        for seg in _SPLIT.split(r.upper()):
            n = len(seg)
            if n < k:
                continue
            if canon:
                rc = revcomp(seg)
                for i in range(n - k + 1):
                    f = seg[i:i + k]
                    b = rc[n - i - k:n - i]
                    c[f if f <= b else b] += 1
            else:
                for i in range(n - k + 1):
                    c[seg[i:i + k]] += 1
    return c


_P61 = (1 << 61) - 1


def kmer_counts_hashed(recs, k, canon=True, seed=12345):
    """Like kmer_counts, but keys are 61-bit polynomial hashes of the canonical k-mer.

    Used for very large k, where slicing every k-mer would be quadratic. Two
    distinct k-mers collide with probability about k / 2^61.
    """
    r = random.Random(seed)
    B = r.randrange(1 << 40, _P61)
    Binv = pow(B, _P61 - 2, _P61)
    top = pow(B, k - 1, _P61)
    val = {"A": 1, "C": 2, "G": 3, "T": 4}
    cval = {"A": 4, "C": 3, "G": 2, "T": 1}
    c = Counter()
    for rec in recs:
        for seg in _SPLIT.split(rec.upper()):
            n = len(seg)
            if n < k:
                continue
            # hf hashes the window read forward, hr its reverse complement:
            # hf = sum_j s[j] B^(k-1-j) and hr = sum_j comp(s[j]) B^j.
            hf = hr = 0
            pw = 1
            for i in range(k):
                hf = (hf * B + val[seg[i]]) % _P61
                hr = (hr + cval[seg[i]] * pw) % _P61
                pw = pw * B % _P61
            for i in range(n - k + 1):
                if i:
                    o, x = seg[i - 1], seg[i + k - 1]
                    hf = ((hf - val[o] * top) * B + val[x]) % _P61
                    hr = ((hr - cval[o]) * Binv + cval[x] * top) % _P61
                c[min(hf, hr) if canon else hf] += 1
    return c


def threshold(cnt, m):
    return cnt if m <= 1 else Counter({x: v for x, v in cnt.items() if v >= m})


def pair_truth(ca, cb):
    """Exact set and weighted quantities for two k-mer count tables."""
    if len(ca) > len(cb):
        small, big = cb, ca
    else:
        small, big = ca, cb
    I = 0; wI = 0
    for x, v in small.items():
        w = big.get(x)
        if w is not None:
            I += 1; wI += min(v, w)
    A, B = len(ca), len(cb)
    wA, wB = sum(ca.values()), sum(cb.values())
    return dict(A=A, B=B, I=I, U=A + B - I, wA=wA, wB=wB, wI=wI, wU=wA + wB - wI)


def prob_jaccard(ca, cb):
    """Probability Jaccard index of the k-mer frequency distributions (Moulton and Jiang)."""
    sa, sb = sum(ca.values()), sum(cb.values())
    if not sa or not sb:
        return float(sa == sb)
    shared = [x for x in ca if x in cb]
    if not shared:
        return 0.0
    # J_P = sum over shared x of 1 / sum_y max(p_y / p_x, q_y / q_x).
    # Group keys by (p, q) to keep this quadratic only in distinct ratio pairs.
    pairs = Counter()
    for x in set(ca) | set(cb):
        pairs[(ca.get(x, 0), cb.get(x, 0))] += 1
    items = list(pairs.items())
    tot = 0.0
    cache = {}
    for x in shared:
        key = (ca[x], cb[x])
        if key not in cache:
            px, qx = key[0] / sa, key[1] / sb
            s = 0.0
            for (p, q), n in items:
                s += n * max((p / sa) / px, (q / sb) / qx)
            cache[key] = 1.0 / s
        tot += cache[key]
    return tot


def measure_value(meas, t, weighted=False):
    """Exact value of a dashing2 measure for row operand a and column operand b."""
    if weighted:
        A, B, I, U = t["wA"], t["wB"], t["wI"], t["wU"]
    else:
        A, B, I, U = t["A"], t["B"], t["I"], t["U"]
    if meas == "intersection":
        return I
    if meas == "union":
        return U
    if A == 0 or B == 0:
        # dashing2 defines ratios with an empty operand as 1 for two empty inputs, else 0.
        sim = float(A == B)
        if meas == "mash":
            return 0.0 if sim else float("inf")
        return sim
    if meas == "similarity":
        return I / U
    if meas == "containment":
        return I / A
    if meas == "symcontain":
        return I / min(A, B)
    if meas == "mash":
        return mash_from_j(I / U, t["k"])
    raise ValueError(meas)


def mash_from_j(j, k):
    if j <= 0:
        return float("inf")
    return -math.log(2 * j / (1 + j)) / k


def j_from_mash(d, k):
    if math.isinf(d) or d > 1e30:
        return 0.0
    return 1.0 / (2.0 * math.exp(k * d) - 1.0)


MEAS_FLAGS = {
    "similarity": [],
    "intersection": ["--intersection"],
    "union": ["--union-size"],
    "containment": ["--containment"],
    "symcontain": ["--symmetric-containment"],
    "mash": ["--mash-distance"],
}


# ---------------------------------------------------------------------------
# dashing2

class D2Error(Exception):
    pass


class D2RowsError(D2Error):
    """dashing2 exited normally but its matrix lacks rows or has them out of order."""
    pass


ROW_LOSS = ("emit-row-loss", "dashing2 cmp occasionally exits 0 with a matrix row missing (about 1 call in "
            "several thousand; see README.md)")


def run_error_known(e, default=None):
    """Known-issue id for a failed dashing2 call: the row-loss symptom, else default."""
    return ROW_LOSS[0] if isinstance(e, D2RowsError) else default


class Matrix:
    """A parsed dashing2 text result. vals maps (row, col) to float; rows/cols are input names."""

    def __init__(self, rows, cols, vals, kind):
        self.rows, self.cols, self.vals, self.kind = rows, cols, vals, kind
        self.stderr = ""

    def compression_base(self):
        """The base b dashing2 chose for log-compressed registers (--fastcmp < 8), if it reported one."""
        m = re.search(r"Truncated via setsketch, a = ([-0-9.eE+inf]+) and b = ([-0-9.eE+inf]+)", self.stderr)
        return float(m.group(2)) if m else None

    def get(self, i, j):
        """Value for row i and column j, reading the upper triangle for symmetric output."""
        if (i, j) in self.vals:
            return self.vals[(i, j)]
        if self.kind == "sym" and (j, i) in self.vals:
            return self.vals[(j, i)]
        return None


def parse_matrix(text):
    rows, vals, sources, kind = [], {}, None, None
    for line in text.splitlines():
        if line.startswith("#Dashing2 Symmetric"):
            kind = "sym"
        elif line.startswith("#Dashing2 Asymmetric"):
            kind = "square"
        elif line.startswith("#Dashing2 Panel"):
            kind = "panel"
        if line.startswith("#Sources"):
            sources = line.split("\t")[1:]
            continue
        if line.startswith("#") or not line.strip():
            continue
        t = line.split("\t")
        r = len(rows)
        rows.append(t[0].strip())
        for j, v in enumerate(t[1:]):
            v = v.strip()
            if v in ("-", ""):
                continue
            vals[(r, j)] = float(v)
    if kind is None or sources is None:
        raise D2Error("could not parse dashing2 output:\n" + text[:2000])
    cols = sources[len(rows):] if kind == "panel" else sources
    return Matrix(rows, cols, vals, kind)


class Dashing2:
    def __init__(self, binary):
        self.binary = os.path.abspath(binary)
        if not os.access(self.binary, os.X_OK):
            raise SystemExit("dashing2 binary not found or not executable: %s" % binary)

    def run(self, args, cwd, timeout=600):
        cmd = [self.binary] + [str(a) for a in args]
        p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout)
        return p

    def cmp(self, args, files, cwd, qfile=None, timeout=600, rows=None):
        """Runs `dashing2 cmp` and returns a Matrix; raises D2Error on failure.

        rows, if given, are the row names dashing2 must print, in order
        (defaults to files); a mismatch raises D2RowsError.
        """
        a = ["cmp"] + list(args)
        if qfile:
            a += ["-Q", qfile]
        a += list(files)
        before = set(os.listdir(cwd))
        try:
            p = self.run(a, cwd, timeout)
        finally:
            remove_new_files(cwd, before)
        if p.returncode != 0:
            raise D2Error("exit %d: %s" % (p.returncode, p.stderr[-800:]))
        M = parse_matrix(p.stdout)
        M.stderr = p.stderr
        want = list(rows if rows is not None else files)
        if want and M.rows != want:
            raise D2RowsError("exit 0 but printed rows %s, expected %s" % (M.rows, want))
        return M

    def cardinalities(self, args, files, cwd, timeout=600):
        """Runs `dashing2 sketch -o` and returns the per-input cardinalities from <out>.names.txt."""
        out = "card_out"
        before = set(os.listdir(cwd))
        res = None
        try:
            p = self.run(["sketch"] + list(args) + ["-o", out] + list(files), cwd, timeout)
            if p.returncode == 0:
                res = []
                with open(os.path.join(cwd, out + ".names.txt")) as f:
                    for ln in f:
                        if ln.startswith("#") or not ln.strip():
                            continue
                        res.append(float(ln.rstrip("\n").split("\t")[1]))
        finally:
            remove_new_files(cwd, before)
        if res is None:
            raise D2Error("exit %d: %s" % (p.returncode, p.stderr[-800:]))
        return res


def remove_new_files(cwd, before):
    """Removes files dashing2 wrote into cwd (k-mer sets, counts, sketches); it does so even without --cache."""
    for fn in os.listdir(cwd):
        if fn in before:
            continue
        p = os.path.join(cwd, fn)
        if os.path.isfile(p):
            os.remove(p)


# ---------------------------------------------------------------------------
# KMC3

class KMC:
    """Thin wrapper around kmc, kmc_dump and kmc_tools.

    KMC is run with -ci1 (its default minimum count is 2), a 32-bit counter
    (-cs4294967295) and -b for non-canonical counting.
    """

    def __init__(self, bindir=None):
        self.kmc = self.tools = self.dump_bin = None
        cands = [bindir] if bindir else []
        cands += os.environ.get("PATH", "").split(os.pathsep)
        for d in cands:
            if not d:
                continue
            if os.path.isfile(d) and os.path.basename(d) == "kmc":
                d = os.path.dirname(d)
            k = os.path.join(d, "kmc")
            if os.access(k, os.X_OK) and os.access(os.path.join(d, "kmc_tools"), os.X_OK):
                self.kmc = k
                self.tools = os.path.join(d, "kmc_tools")
                self.dump_bin = os.path.join(d, "kmc_dump")
                break
        self.available = self.kmc is not None

    def _run(self, cmd, cwd):
        p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
        if p.returncode != 0:
            raise RuntimeError("KMC command failed (%d): %s\n%s" % (p.returncode, " ".join(cmd), p.stderr[-800:]))
        return p

    def _count_cmd(self, path, k, canon, cwd, name):
        fmt = "q" if ".fq" in path else "m"  # -fq FASTQ, -fm multi-line FASTA
        tmp = "tmp_" + name
        os.makedirs(os.path.join(cwd, "kmc", tmp), exist_ok=True)
        cmd = [self.kmc, "-k%d" % k, "-ci1", "-cs4294967295", "-f" + fmt, "-t2", "-hp"]
        if not canon:
            cmd.append("-b")
        return cmd + [os.path.abspath(os.path.join(cwd, path)), name, tmp]

    def count(self, path, k, canon, cwd, name):
        """Counts k-mers of cwd/path into the database cwd/kmc/name and returns name."""
        self._run(self._count_cmd(path, k, canon, cwd, name), os.path.join(cwd, "kmc"))
        return name

    def count_many(self, jobs, cwd, par=8):
        """Runs several counts concurrently; jobs are (path, k, canon, name) tuples.

        KMC spends most of a sub-second run idle in its second stage, so running
        the counts side by side saves wall time.
        """
        kd = os.path.join(cwd, "kmc")
        os.makedirs(kd, exist_ok=True)
        pending = list(jobs)
        while pending:
            batch, pending = pending[:par], pending[par:]
            procs = [(j, subprocess.Popen(self._count_cmd(j[0], j[1], j[2], cwd, j[3]),
                                          cwd=kd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)) for j in batch]
            for j, p in procs:
                out, err = p.communicate()
                if p.returncode != 0:
                    raise RuntimeError("KMC count failed (%d) for %s: %s" % (p.returncode, j, err[-800:]))
        return [j[3] for j in jobs]

    def dump(self, db, cwd, ci=1):
        """Returns the k-mer counts of database cwd/kmc/db with count >= ci."""
        kd = os.path.join(cwd, "kmc")
        out = db + ".dump.txt"
        p = subprocess.run([self.dump_bin, "-ci%d" % ci, db, out], cwd=kd, capture_output=True, text=True)
        # kmc_dump 3.2.4 sometimes exits with status 1 and no message on an
        # empty database; accept that case when the (empty) dump was written.
        if p.returncode != 0 and (p.stderr.strip() or not os.path.exists(os.path.join(kd, out))
                                  or os.path.getsize(os.path.join(kd, out)) != 0):
            raise RuntimeError("kmc_dump failed (%d) on %s: %s" % (p.returncode, db, p.stderr[-800:]))
        c = Counter()
        with open(os.path.join(kd, out)) as f:
            for ln in f:
                s, v = ln.split("\t")
                c[s] = int(v)
        os.remove(os.path.join(kd, out))
        return c

    def pair(self, db1, db2, cwd, ci=1):
        """Uses kmc_tools simple to compute |A n B|, |A u B|, the sum of min counts and the sum of max counts.

        Each output gets an explicit -ci1: when the inputs are filtered with
        -ci > 1 and the output has no -ci of its own, kmc_tools 3.2.4 writes an
        empty output database.
        """
        kd = os.path.join(cwd, "kmc")
        i, u = db1 + "_x_" + db2 + "_I", db1 + "_x_" + db2 + "_U"
        self._run([self.tools, "-t2", "-hp", "simple", db1, "-ci%d" % ci, db2, "-ci%d" % ci,
                   "intersect", i, "-ocmin", "-ci1", "-cs4294967295",
                   "union", u, "-ocmax", "-ci1", "-cs4294967295"], kd)
        ci_ = self.dump(i, cwd)
        cu = self.dump(u, cwd)
        for d in (i, u):
            for ext in (".kmc_pre", ".kmc_suf"):
                try:
                    os.remove(os.path.join(kd, d + ext))
                except OSError:
                    pass
        return dict(I=len(ci_), U=len(cu), wI=sum(ci_.values()), wU=sum(cu.values()))


# ---------------------------------------------------------------------------
# Command line, reporting, statistics

def common_args(desc, presets):
    ap = argparse.ArgumentParser(description=desc)
    ap.add_argument("--dashing2", default=os.environ.get("DASHING2", os.path.join(HERE, "..", "..", "dashing2")),
                    help="dashing2 binary (env DASHING2; default ../../dashing2)")
    ap.add_argument("--kmc-bin", default=os.environ.get("KMC_BIN"),
                    help="directory holding kmc, kmc_tools and kmc_dump (env KMC_BIN; default: search PATH)")
    ap.add_argument("--no-kmc", action="store_true", help="skip KMC comparisons even if KMC is available")
    ap.add_argument("--seed", type=int, default=1, help="master seed (default 1)")
    ap.add_argument("--preset", choices=sorted(presets), default="quick")
    ap.add_argument("--trial", type=int, action="append", help="run only this trial index (repeatable)")
    ap.add_argument("--jobs", "-j", type=int, default=4, help="trials run in parallel (default 4)")
    ap.add_argument("--keep", action="store_true", help="keep the work directories of failing trials")
    ap.add_argument("--workdir", default=None, help="parent directory for temporary trial directories")
    ap.add_argument("--strict", action="store_true", help="count failures matching known issues as failures")
    ap.add_argument("-v", "--verbose", action="store_true")
    return ap


def trial_rng(suite, master, idx):
    # Seeding with a string is deterministic across Python processes and versions >= 3.2.
    return random.Random("%s:%d:%d" % (suite, master, idx))


def repro_cmd(script, args, idx):
    """One-line command that reruns trial idx alone, verbosely, keeping its inputs."""
    parts = ["DASHING2=" + shlex.quote(os.path.abspath(args.dashing2))]
    if getattr(args, "kmc_bin", None):
        parts.append("KMC_BIN=" + shlex.quote(args.kmc_bin))
    parts += ["python3", os.path.relpath(os.path.join(HERE, script)), "--seed", str(args.seed),
              "--preset", args.preset, "--trial", str(idx), "-v", "--keep"]
    if getattr(args, "no_kmc", False):
        parts.append("--no-kmc")
    parts += [shlex.quote(x) for x in getattr(args, "repro_extra", [])]
    return " ".join(parts)


STRICT = [False]  # set by --strict: known issues count as failures


class TrialResult:
    def __init__(self, idx, desc):
        self.idx = idx
        self.desc = desc
        self.checks = 0
        self.fails = []
        self.skips = []
        self.data = []  # suite-specific records (for aggregate statistics)
        self.tags = {}  # feature -> value, used to break failures down by feature
        self.xfails = []  # (known issue id, message) for expected failures
        self.workdir = None
        self.elapsed = 0.0

    def check(self, ok, name, detail="", known=None):
        """Records one check. A failure matching a known issue id (known) is an expected failure."""
        self.checks += 1
        if not ok:
            if known and not STRICT[0]:
                self.xfails.append((known, "%s: %s" % (name, detail)))
            else:
                self.fails.append("%s: %s%s" % (name, detail, " (known issue %s)" % known if known else ""))
        return ok

    def fail(self, name, detail):
        self.checks += 1
        self.fails.append("%s: %s" % (name, detail))


def close(a, b, rel=2e-6, abs_=1e-6):
    """Equality up to the float32 rounding dashing2 applies when printing."""
    if a is None or b is None:
        return False
    if math.isinf(b):
        return a > 1e30 or math.isinf(a)
    return abs(a - b) <= max(abs_, rel * abs(b))


def run_trials(suite, script, args, ntrials, trial_fn):
    """Runs trial_fn(idx, rng, workdir) for each trial (optionally in parallel) and prints results.

    Returns the list of TrialResult objects.
    """
    STRICT[0] = getattr(args, "strict", False)
    if STRICT[0]:
        args.repro_extra = getattr(args, "repro_extra", []) + ["--strict"]
    idxs = args.trial if args.trial else list(range(ntrials))
    results = []
    t0 = time.time()

    def one(idx):
        rng = trial_rng(suite, args.seed, idx)
        wd = tempfile.mkdtemp(prefix="d2rand_%s_%d_" % (suite, idx), dir=args.workdir)
        st = time.time()
        try:
            res = trial_fn(idx, rng, wd)
        except Exception as e:  # harness errors are failures too
            import traceback
            res = TrialResult(idx, "exception")
            res.fail("harness", "%s\n%s" % (e, traceback.format_exc()))
        res.elapsed = time.time() - st
        if res.fails and args.keep:
            res.workdir = wd
        else:
            shutil.rmtree(wd, ignore_errors=True)
        return res

    if args.jobs > 1 and len(idxs) > 1:
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=args.jobs) as ex:
            for res in ex.map(one, idxs):
                results.append(res)
                _print_trial(res, args, script)
    else:
        for idx in idxs:
            res = one(idx)
            results.append(res)
            _print_trial(res, args, script)
    results.sort(key=lambda r: r.idx)
    return results, time.time() - t0


def _print_trial(res, args, script):
    status = "FAIL" if res.fails else ("XFAIL" if res.xfails else "ok")
    if res.fails or res.xfails or args.verbose:
        print("[trial %d] %s (%d checks, %.1fs) %s" % (res.idx, status, res.checks, res.elapsed, res.desc))
        for s in res.skips:
            print("    SKIP " + s)
        for kid, f in res.xfails[:5 if not args.verbose else None]:
            print("    XFAIL[%s] %s" % (kid, f))
        for f in res.fails[:20 if not args.verbose else None]:
            print("    FAIL " + f)
        if len(res.fails) > 20 and not args.verbose:
            print("    ... %d more failures" % (len(res.fails) - 20))
        if res.fails or res.xfails:
            print("    repro: " + repro_cmd(script, args, res.idx))
            if res.workdir:
                print("    inputs kept in " + res.workdir)
        sys.stdout.flush()


def k_band(k, long=False):
    """Coarse k ranges matching dashing2's encoders (direct 64-bit, direct 128-bit, rolling hash)."""
    lim = 64 if long else 32
    if k <= lim:
        return "k<=%d" % lim
    if k <= 256:
        return "%d<k<=256" % lim
    return "k>256"


def failure_breakdown(results):
    """Prints, for each tagged feature, how many trials with each value failed."""
    if not any(r.fails for r in results):
        return
    print("failure breakdown (failed trials / trials):")
    keys = sorted({k for r in results for k in r.tags})
    for key in keys:
        tot, bad = Counter(), Counter()
        for r in results:
            if key in r.tags:
                v = r.tags[key]
                tot[v] += 1
                bad[v] += bool(r.fails)
        print("  %-10s " % key + "  ".join("%s:%d/%d" % (v, bad[v], tot[v]) for v in sorted(tot, key=str)))
    def cname(f):
        return re.sub(r"\[.*", "", f.split(":")[0])
    names = Counter()
    for r in results:
        for f in r.fails:
            names[cname(f)] += 1
    print("  failed checks by name: " + ", ".join("%s=%d" % kv for kv in names.most_common(25)))
    print("  most frequent failed checks per feature value:")
    for key in keys:
        for v in sorted({r.tags[key] for r in results if key in r.tags}, key=str):
            c = Counter(cname(f) for r in results if r.tags.get(key) == v for f in r.fails)
            if c:
                print("    %s=%s: %s" % (key, v, ", ".join("%s=%d" % kv for kv in c.most_common(4))))


def summarize(name, results, elapsed, extra_fail=0, known=None):
    failure_breakdown(results)
    xf = Counter(kid for r in results for kid, _ in r.xfails)
    for kid, nchk in sorted(xf.items()):
        ntr = sum(1 for r in results if any(k == kid for k, _ in r.xfails))
        print("XFAIL %s: %d checks in %d trials match known issue: %s" % (kid, nchk, ntr, (known or {}).get(kid, "")))
    n = len(results)
    nf = sum(1 for r in results if r.fails)
    nc = sum(r.checks for r in results)
    nfc = sum(len(r.fails) for r in results)
    skips = sorted({s for r in results for s in r.skips})
    for s in skips:
        print("SKIP " + s)
    print("%s: %d trials, %d checks, %d failed checks in %d trials, %d aggregate failures, %.1fs" %
          (name, n, nc, nfc, nf, extra_fail, elapsed))
    ok = nf == 0 and extra_fail == 0
    print("%s: %s" % (name, "PASS" if ok else "FAIL"))
    return 0 if ok else 1


def mean_sd(xs):
    n = len(xs)
    if n == 0:
        return float("nan"), float("nan")
    m = sum(xs) / n
    if n < 2:
        return m, float("nan")
    return m, math.sqrt(sum((x - m) ** 2 for x in xs) / (n - 1))


# ---------------------------------------------------------------------------
# Standard errors of sketch estimates
#
# Every sketch estimate dashing2 prints is a function of three noisy
# quantities: the register-agreement estimate of the (set, weighted or
# probability) Jaccard index J, and the two cardinality estimates A and B:
#   similarity = J,  intersection = J (A + B) / (1 + J),  union = (A + B) / (1 + J),
#   containment = intersection / A,  symmetric containment = intersection / min(A, B),
#   Mash distance = -ln(2J / (1 + J)) / k.
# We model the covariance of (J, A, B) for m registers as for m independent
# min-hash registers (exponential minima over A \ B, B \ A and A n B):
#   Var J = p (1 - p) / (m (1 - c)^2), p = J + (1 - J) c,
#   Var A = kappa^2 (A^2 / m + 1), Var B likewise, Corr(A, B) = J,
#   Cov(J, A) = kappa J A |B \ A| / (|A u B| m), and symmetrically for B,
# where c is the probability that two registers holding different minima
# compare equal after compression (0 for full registers, 2^-bits for b-bit
# signatures, about ln(b) / 4 for registers log-compressed with base b) and
# kappa is the relative standard error constant of the cardinality estimator
# (0 when the cardinality is an exact total, as for -B and --prob).
# Standard errors of the derived measures follow by the delta method.

def collision_prob(compress, bits=None, base=None):
    if compress == "bbit":
        return 2.0 ** -bits
    if compress == "log":
        return math.log(base) / 4.0 if base and base > 1 else 0.0
    return 0.0


def jcov(J, nA, nB, nI, m, kappa, c):
    """Covariance matrix of (J_hat, A_hat, B_hat) under the model above."""
    nU = nA + nB - nI
    p = J + (1 - J) * c
    vJ = p * (1 - p) / (m * (1 - c) ** 2)
    # The "+ 1" is one element of resolution: for sets much smaller than m the
    # estimators count occupied registers, and their error is a whole number
    # of register collisions rather than a continuous quantity.
    vA = kappa ** 2 * (nA * nA / m + 1)
    vB = kappa ** 2 * (nB * nB / m + 1)
    cAB = kappa ** 2 * J * nA * nB / m
    cJA = kappa * J * nA * (nB - nI) / (nU * m) if nU else 0.0
    cJB = kappa * J * nB * (nA - nI) / (nU * m) if nU else 0.0
    return [[vJ, cJA, cJB], [cJA, vA, cAB], [cJB, cAB, vB]]


def derived(meas, J, A, B):
    S = A + B
    I = J * S / (1 + J)
    if meas == "similarity":
        return J
    if meas == "intersection":
        return I
    if meas == "union":
        return S - I
    if meas == "containment":
        return I / A
    if meas == "symcontain":
        return I / min(A, B)
    raise ValueError(meas)


def delta_se(meas, J, nA, nB, nI, m, kappa, c):
    """Delta-method standard error of a derived measure at the true point."""
    cov = jcov(J, nA, nB, nI, m, kappa, c)
    x = [J, float(nA), float(nB)]
    grad = []
    for i in range(3):
        h = 1e-6 * max(abs(x[i]), 1e-3)
        lo, hi = list(x), list(x)
        lo[i] -= h; hi[i] += h
        if i == 0:
            lo[0] = max(lo[0], 0.0); hi[0] = min(hi[0], 1.0)
        grad.append((derived(meas, *hi) - derived(meas, *lo)) / (hi[i] - lo[i]))
    v = sum(grad[i] * cov[i][j] * grad[j] for i in range(3) for j in range(3))
    return math.sqrt(max(v, 0.0))


class SketchMode:
    """A dashing2 sketch configuration and the error model that goes with it.

    kind is "set" (OPH, FullSetSketch: set Jaccard, estimated cardinalities),
    "weighted" (BagMinHash: weighted Jaccard, exact total counts) or "prob"
    (ProbMinHash: probability Jaccard, exact total counts). comp is None,
    "bbit" (b-bit signatures with `bits` bits) or "log" (registers
    log-compressed with base `base`; dashing2 picks the base itself unless the
    flags fix it, and reports it on stderr).
    """

    def __init__(self, name, flags, kind, comp=None, bits=None, base=None):
        self.name, self.flags, self.kind = name, list(flags), kind
        self.comp, self.bits, self.base = comp, bits, base

    @property
    def kappa(self):
        return 1.0 if self.kind == "set" else 0.0

    def measures(self):
        """Measures with a defined exact counterpart (ProbMinHash has only the similarity)."""
        if self.kind == "prob":
            return ["similarity", "mash"]
        return list(MEAS_FLAGS)


def sketch_modes():
    L = []
    L.append(SketchMode("oph", [], "set"))
    L.append(SketchMode("full", ["--full"], "set"))
    L.append(SketchMode("bmh", ["-B"], "weighted"))
    L.append(SketchMode("pmh", ["--prob"], "prob"))
    for base, bf, kind in (("oph", [], "set"), ("full", ["--full"], "set"), ("bmh", ["-B"], "weighted"),
                           ("pmh", ["--prob"], "prob")):
        for nb in (1, 2, 4):
            if kind == "prob" and nb == 4:
                continue
            L.append(SketchMode("%s-fc%d" % (base, nb), bf + ["--fastcmp", str(nb)], kind, "log"))
            if kind != "prob":
                L.append(SketchMode("%s-bb%d" % (base, nb), bf + ["--bbit-sigs", "--fastcmp", str(nb)], kind,
                                    "bbit", bits=8 * nb))
    # Directly sketched log-compressed SetSketches with preset (a, b).
    L.append(SketchMode("full-ab-bytes", ["--full", "--fastcmp-bytes"], "set", "log", base=1.2))
    L.append(SketchMode("full-ab-shorts", ["--full", "--fastcmp-shorts"], "set", "log", base=1.0005))
    L.append(SketchMode("full-ab-words", ["--full", "--fastcmp-words"], "set", "log", base=1.0000000109723500835))
    return L


SIGMA_FLOOR_REL = 2e-6  # dashing2 prints float32 values


def sketch_truth(mode, ca, cb, k):
    """Exact (J, nA, nB, nI) for the quantity mode estimates; ca and cb are k-mer count tables."""
    t = pair_truth(ca, cb)
    if mode.kind == "set":
        nA, nB, nI = t["A"], t["B"], t["I"]
    else:
        nA, nB, nI = t["wA"], t["wB"], t["wI"]
    nU = nA + nB - nI
    if mode.kind == "prob":
        J = prob_jaccard(ca, cb)
    else:
        J = nI / nU if nU else 1.0
    return J, nA, nB, nI


def sketch_z(mode, meas, est, J, nA, nB, nI, m, k, base=None, diag=False, shift=0.0):
    """Returns (z, regular) for one printed estimate, or None if the pair has no defined exact value.

    diag selects the cardinality check (the diagonal of --union-size --square).
    regular is False for pairs too close to J = 0 or 1 for the normal
    approximation (fewer than 5 expected equal or unequal registers); those
    still get the per-trial bound but are left out of aggregate statistics.
    """
    if est is None:
        return None
    if diag:
        if nA == 0:
            return (0.0 if abs(est) < 1e-9 else float("inf")), False
        if mode.kind != "set":
            # Exact total counts: z only reflects printing precision, so keep it out of aggregates.
            return (est - nA) / (SIGMA_FLOOR_REL * max(1.0, nA)), False
        se = mode.kappa * math.sqrt(nA * nA / m + 1) + SIGMA_FLOOR_REL * nA
        return (est - nA) / se, True
    if nA == 0 or nB == 0:
        return None
    if meas == "mash":
        # Compare on the Jaccard scale: invert d = -ln(2J / (1 + J)) / k.
        est = j_from_mash(est, k)
        meas = "similarity"
    c = collision_prob(mode.comp, mode.bits, base if base is not None else mode.base)
    J = J + shift  # shift moves the expectation for a documented bias (shift is 0 otherwise)
    if mode.kind == "prob":
        truth = J
        se1 = math.sqrt((J + (1 - J) * c) * (1 - J - (1 - J) * c) / m) / (1 - c)
        sef = 1.0 / m
    else:
        truth = derived(meas, J, nA, nB)
        se1 = delta_se(meas, J, nA, nB, nI, m, mode.kappa, c)
        # One register's worth of resolution on the J scale, a continuity
        # correction for the discrete count of equal registers near J = 0 or 1.
        sef = abs(derived(meas, min(1.0, J + 1.0 / m), nA, nB) - truth)
    se = math.sqrt(se1 ** 2 + sef ** 2 + (SIGMA_FLOOR_REL * abs(truth)) ** 2)
    regular = m * J >= 5 and m * (1 - J) >= 5
    if se == 0:
        return (0.0 if est == truth else float("inf")), False
    return (est - truth) / se, regular


def aggregate_check(zs, label, bias_allow=0.1, sd_range=(0.5, 1.6), min_n=8, n_eff=None, sd_min_units=30):
    """Bias and spread checks on a list of z-scores; returns (ok, message).

    n_eff is the number of independent units (trials) behind zs. The bias
    bound is 3 / sqrt(n_eff) + bias_allow; this stays conservative when the
    z-scores of one trial are correlated, because the mean z of a trial has
    variance at most 1. The spread check needs at least sd_min_units units:
    with fewer, a single trial's correlated z-scores dominate the SD.
    """
    n = len(zs)
    if n < min_n:
        return True, "%s: n=%d (too few for aggregate checks)" % (label, n)
    mu, sd = mean_sd(zs)
    ne = n_eff or n
    bound = 3.0 / math.sqrt(ne) + bias_allow
    ok = abs(mu) <= bound
    if ne >= sd_min_units:
        ok = ok and sd_range[0] <= sd <= sd_range[1]
        sdmsg = "SD z %.3f (range %.1f-%.1f)" % (sd, sd_range[0], sd_range[1])
    else:
        sdmsg = "SD z %.3f (not checked, fewer than %d trials)" % (sd, sd_min_units)
    return ok, "%s: n=%d (trials %d) mean z %+.3f (bound %.3f) %s" % (label, n, ne, mu, bound, sdmsg)
