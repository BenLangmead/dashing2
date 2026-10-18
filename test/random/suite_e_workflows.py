#!/usr/bin/env python3
"""Suite E: dashing2 workflow and output modes checked for internal consistency and against exact oracles.

Each trial belongs to one family (round robin by trial index):

  nn        --topk K and --similarity-threshold T against the full matrix of the same mode
  greedy    --greedy T and --greedy TE against the clustering rule and an exact simulation
  formats   text, binary, --cmpout, --phylip, square, panel and condensed layouts agree
  roundtrip sketch -o / --cache / --outprefix files reloaded with --presketched
  kmers     --set/-J k-mer files and --save-kmers/--save-kmercounts files decoded exactly
  contain   dashing2 contain coverage and depth against an exact recomputation
  wsketch   dashing2 wsketch (CSR and 1-D) sampled ids, weights and registers
  seq       -G minimizer sequences and printmin against the exact k-mer sequence
  bed       BED interval sketches against exact interval sets
  measures  -Q panels, --square and measure identities across all measures

Usage: DASHING2=/path/to/dashing2 python3 test/random/suite_e_workflows.py [--preset quick|thorough] [--seed N]
"""
import math
import os
import struct
import sys
import time
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from d2rand import (common_args, kmer_counts, read_fastx, revcomp, run_trials, summarize, TrialResult,
                    write_fastx, random_dna, mutate)
import d2rand_e as E

SUITE = "e"
SCRIPT = "suite_e_workflows.py"
FAMILIES = ["nn", "greedy", "formats", "roundtrip", "kmers", "contain", "wsketch", "seq", "bed", "measures"]
PRESETS = {
    # trials (round robin over FAMILIES), largest number of inputs, longest root genome
    "quick": dict(trials=150, maxn=14, maxlen=3000),
    "thorough": dict(trials=800, maxn=32, maxlen=8000),
}

# Known issues in the tested code. A failing check tagged with one of these ids is
# an expected failure (XFAIL) unless --strict is given. See README.md, suite E. The
# issues this suite found are fixed; their ids are in d2rand.FIXED_IDS.
KNOWN = {}

# Sketch modes for the comparison families: name, flags, kind.
MODES = [
    ("oph", [], "set"),
    ("full", ["--full"], "set"),
    ("bmh", ["-B"], "weighted"),
    ("pmh", ["--prob"], "prob"),
    ("set", ["--set"], "exact"),
    ("cd", ["-J"], "exact"),
    ("full-fc2", ["--full", "--fastcmp", "2"], "set"),
    ("full-fc1", ["--full", "--fastcmp", "1"], "set"),
    ("oph-bb2", ["--bbit-sigs", "--fastcmp", "2"], "set"),
    ("full-shorts", ["--full", "--fastcmp-shorts"], "set"),
]
MODE = {m[0]: m for m in MODES}
MEASURE_FLAGS = {
    "similarity": [],
    "mash": ["--mash-distance"],
    "containment": ["--containment"],
    "symcontain": ["--symmetric-containment"],
    "intersection": ["--intersection"],
    "union": ["--union-size"],
}
DISTANCES = {"mash"}  # documented distances; every other measure ranks larger values first


def is_better(meas, a, b):
    return a < b if meas in DISTANCES else a > b


def is_null(meas, v):
    """True for a value that means 'nothing shared' (similarity 0 or distance at the float maximum)."""
    if meas in DISTANCES:
        return v >= E.FLT_MAX or math.isinf(v)
    return v == 0


class Ctx:
    """Per-trial state: binary, work directory, result object and a log of commands for messages."""

    def __init__(self, args, wd, res):
        self.bin = os.path.abspath(args.dashing2)
        self.wd = wd
        self.res = res
        self.args = args

    def run(self, a, env=None, timeout=600):
        return E.run(self.bin, a, self.wd, timeout=timeout, env=env)

    def ok(self, r, name):
        """Records a failed call as a failure and returns False; True if it exited 0."""
        if r.rc != 0:
            self.res.fail(name, r.brief())
            return False
        return True

    def path(self, *p):
        return os.path.join(self.wd, *p)


def base_flags(p):
    a = ["-k", p["k"], "-S", p["S"]]
    if p.get("seed") is not None:
        a += ["--seed", p["seed"]]
    if not p.get("canon", True):
        a.append("--no-canon")
    if p.get("long"):
        a.append("-2")
    return a + list(MODE[p["mode"]][1])


def square_matrix(ctx, flags, files, name="square"):
    """Runs --square and returns {(i, j): value} over input indices, or None after recording a failure."""
    r = ctx.run(["cmp", "--square"] + flags + files)
    if not ctx.ok(r, name):
        return None
    kind, rows, cols, vals = E.parse_text_matrix(r.text)
    ctx.res.check(rows == files and cols == files, name + "-labels", "rows %s cols %s expected %s" % (rows, cols, files))
    return vals


# ---------------------------------------------------------------------------
# Family nn: --topk and --similarity-threshold

def fam_nn(ctx, rng, preset):
    res = ctx.res
    n = rng.randint(4, preset["maxn"])
    p = dict(k=rng.randint(12, 28), S=rng.choice([64, 100, 128, 256, 512, 1000, 1024]), seed=rng.choice([None, rng.randint(1, 10 ** 6)]),
             canon=rng.random() < 0.85, mode=rng.choice([m[0] for m in MODES]))
    meas = rng.choice(["similarity", "similarity", "mash", "containment", "symcontain", "intersection"])
    if MODE[p["mode"]][2] == "prob":
        meas = rng.choice(["similarity", "mash"])
    K = rng.choice([1, 1, 2, 3, n - 1, rng.randint(1, n + 1)])
    pt = rng.randint(2, 8)
    # Sets with many identical inputs fill the candidate lists with ties.
    dups = rng.random() < 0.3
    files, _, _ = E.clustered_inputs(rng, ctx.wd, n, p["k"], maxlen=preset["maxlen"], dup_prob=0.6 if dups else 0.08)
    res.desc = ("dup-heavy " if dups else "") + "nn n=%d mode=%s measure=%s K=%d k=%d S=%d seed=%s%s -p%d" % (
        n, p["mode"], meas, K, p["k"], p["S"], p["seed"], "" if p["canon"] else " --no-canon", pt)
    res.tags = {"family": "nn", "mode": p["mode"], "measure": meas}
    flags = base_flags(p) + MEASURE_FLAGS[meas]
    known = "symcontain-as-distance" if meas == "symcontain" else None
    M = square_matrix(ctx, flags, files)
    if M is None:
        return
    Msim = M if meas == "similarity" else square_matrix(ctx, base_flags(p), files, "square-sim")
    if Msim is None:
        return
    strict_mode = p["mode"] in STRICT_LSH
    # ntoquery = min(n - 1, 3.5 K) candidates per item (src/index_build.cpp); with n - 1 <= 3.5 K every
    # item sharing an LSH bucket is compared, so the lists are exact apart from bucket misses.
    full_cand = n - 1 <= int(K * 3.5)

    # Top-k: text with -p 1, binary with -p 1, text with -p N.
    rt = ctx.run(["cmp", "--topk", K, "-p", 1] + flags + files)
    rb = ctx.run(["cmp", "--topk", K, "-p", 1, "--binary-output"] + flags + files)
    rp = ctx.run(["cmp", "--topk", K, "-p", pt] + flags + files)
    if ctx.ok(rt, "topk-run") and ctx.ok(rb, "topk-bin-run") and ctx.ok(rp, "topk-pN-run"):
        T = E.parse_nn_text(rt.text)
        check_nn_rows(res, "topk", T, files, M, meas, known)
        try:
            B, used, size = E.parse_csr(rb.out)
            res.check(used == size, "topk-bin-size", "CSR uses %d of %d bytes" % (used, size))
            same = len(B) == len(T) and all(
                len(B[i]) == len(T[i][1]) and all(files[j] == nm and E.near32(v, tv) for (j, v), (nm, tv) in zip(B[i], T[i][1]))
                for i in range(len(B)))
            res.check(same, "topk-bin-vs-text", "binary rows %s text rows %s" % (B[:3], [t[1] for t in T[:3]]))
            for i, row in enumerate(B):
                for j, v in row:
                    if not E.same32(v, M.get((i, j))):
                        res.check(False, "topk-bin-value[%d,%d]" % (i, j), "binary %r matrix %r" % (v, M.get((i, j))), known)
        except (struct.error, ValueError) as e:
            res.fail("topk-bin-parse", str(e))
        m1 = check_topk(res, "topk", T, files, M, Msim, meas, K, strict_mode, full_cand, known)
        mN = check_topk(res, "topkN", E.parse_nn_text(rp.text), files, M, Msim, meas, K, strict_mode, full_cand, known)
        # Which candidate a capped query drops depends on the order of ids in the buckets, which the
        # threads change; when neither run misses a true neighbor the lists must not depend on -p.
        res.check(rp.text == rt.text or not full_cand, "topk-threads", "-p %d output differs from -p 1" % pt,
                  known or ("lsh-self-candidate" if m1 + mN else None))

    # Similarity threshold: a value strictly between two observed off-diagonal values.
    offd = sorted({v for (i, j), v in M.items() if i != j and not is_null(meas, v)})
    thr = None
    if len(offd) >= 2:
        a = rng.randrange(len(offd) - 1)
        thr = E.between32(offd[a], offd[a + 1])
    if thr is not None and thr > 0:
        rt = ctx.run(["cmp", "--similarity-threshold", repr(thr), "-p", 1] + flags + files)
        rb = ctx.run(["cmp", "--similarity-threshold", repr(thr), "-p", 1, "--binary-output"] + flags + files)
        rp = ctx.run(["cmp", "--similarity-threshold", repr(thr), "-p", pt] + flags + files)
        if ctx.ok(rt, "thresh-run") and ctx.ok(rb, "thresh-bin-run") and ctx.ok(rp, "thresh-pN-run"):
            T = E.parse_nn_text(rt.text)
            check_nn_rows(res, "thresh", T, files, M, meas, known)
            try:
                B, used, size = E.parse_csr(rb.out)
                res.check(used == size, "thresh-bin-size", "CSR uses %d of %d bytes" % (used, size))
                tb = [[(files[j], v) for j, v in row] for row in B]
                tt = [r[1] for r in T]
                same = len(tb) == len(tt) and all(len(x) == len(y) and all(a[0] == b[0] and E.near32(a[1], b[1]) for a, b in zip(x, y))
                                                for x, y in zip(tb, tt))
                res.check(same, "thresh-bin-vs-text", "binary %s text %s" % (tb[:3], tt[:3]))
            except (struct.error, ValueError) as e:
                res.fail("thresh-bin-parse", str(e))
            m1 = check_threshold(res, "thresh", T, files, M, Msim, meas, thr, strict_mode, known)
            mN = check_threshold(res, "threshN", E.parse_nn_text(rp.text), files, M, Msim, meas, thr, strict_mode, known)
            res.check(rp.text == rt.text, "thresh-threads", "-p %d output differs from -p 1" % pt,
                      known or ("lsh-self-candidate" if m1 + mN else None))
    res.data.append(("nn", p["mode"], meas, full_cand))


def check_nn_rows(res, tag, T, files, M, meas, known):
    """Structure shared by top-k and threshold lists: one row per input in order, values from the matrix."""
    res.check([t[0] for t in T] == files, tag + "-rows", "rows %s expected %s" % ([t[0] for t in T], files))
    idx = {f: i for i, f in enumerate(files)}
    for i, (name, nb) in enumerate(T[:len(files)]):
        names = [x[0] for x in nb]
        res.check(name not in names, "%s-self[%d]" % (tag, i), "row lists itself: %s" % nb)
        res.check(len(set(names)) == len(names), "%s-dup[%d]" % (tag, i), "neighbor listed twice: %s" % nb)
        for nm, v in nb:
            j = idx.get(nm)
            if j is None:
                res.fail("%s-name[%d]" % (tag, i), "unknown neighbor %r" % nm)
                continue
            if not E.near32(v, M.get((i, j))):
                res.check(False, "%s-value[%d,%d]" % (tag, i, j), "listed %r matrix %r" % (v, M.get((i, j))), known)
        vals = [v for _, v in nb]
        order_ok = all(not is_better(meas, vals[t + 1], vals[t]) for t in range(len(vals) - 1))
        res.check(order_ok, "%s-order[%d]" % (tag, i), "%s not best first for %s: %s" % (tag, meas, vals), known)


def miss_is_strong(strict, Msim, i, j):
    """Whether a missing neighbor counts as a failure rather than an LSH approximation.

    In modes whose LSH keys are the registers that the similarity counts
    (strict), every pair with a positive value shares a bucket and must be
    found. Elsewhere (exact modes index bottom-k hashes, --fastcmp and
    --bbit-sigs index the full registers while comparing compressed ones) a
    pair is only expected to be found when its similarity is at least
    STRONG_J, where sharing no bucket has probability below (1 - 0.25)^64.
    """
    return strict or (Msim.get((i, j)) or 0) >= STRONG_J


STRONG_J = 0.25
# Modes in which the LSH index is keyed by exactly the registers the similarity compares.
STRICT_LSH = {"oph", "full", "bmh", "pmh", "full-shorts"}


def check_topk(res, tag, T, files, M, Msim, meas, K, strict, full_cand, known):
    """Each row must hold the K best neighbors by the matrix (ties at the K-th value kept), nothing null."""
    n = len(files)
    idx = {f: i for i, f in enumerate(files)}
    nwant = nmiss = nstrong = 0
    for i, (name, nb) in enumerate(T[:n]):
        cand = [(M[(i, j)], j) for j in range(n) if j != i and (i, j) in M and not is_null(meas, M[(i, j)])]
        cand.sort(key=lambda x: (-x[0] if meas not in DISTANCES else x[0]))
        if len(cand) > K:
            kth = cand[K - 1][0]
            want = {j for v, j in cand if not is_better(meas, kth, v)}
        else:
            want = {j for v, j in cand}
        # Distance lists may be padded with neighbors at the maximum distance; similarity lists drop zeros.
        listed = [(idx[nm], v) for nm, v in nb if nm in idx]
        got = {j for j, v in listed if not is_null(meas, v)}
        if meas not in DISTANCES:
            res.check(all(not is_null(meas, v) for _, v in listed), "%s-null[%d]" % (tag, i), "zero listed: %s" % nb, known)
        vals = [v for _, v in nb]
        if len(vals) > K:
            res.check(all(E.near32(v, vals[K - 1]) for v in vals[K:]), "%s-len[%d]" % (tag, i),
                      "%d listed for K=%d beyond ties: %s" % (len(vals), K, vals), known)
        miss = want - got
        extra = got - want
        # A listed item outside the true top K is only acceptable when a missed true neighbor made room for it.
        if extra and not miss:
            res.check(False, "%s-extra[%d]" % (tag, i), "listed %s not in true top %d %s (values %s)" % (
                sorted(extra), K, sorted(want), {j: M[(i, j)] for j in extra | want}), known)
        nwant += len(want)
        nmiss += len(miss)
        # With fewer candidates than inputs (n - 1 > 3.5 K) the index keeps the first candidates it meets
        # rather than the best ones, so misses there only enter the recall statistics.
        strong = sorted(j for j in miss if full_cand and miss_is_strong(strict, Msim, i, j))
        if strong:
            nstrong += len(strong)
            res.check(False, "%s-miss[%d]" % (tag, i), "true top-%d neighbors %s missing (values %s, similarity %s; listed %s)" % (
                K, strong, {j: M[(i, j)] for j in strong}, {j: Msim.get((i, j)) for j in strong}, nb), known or "lsh-self-candidate")
    if tag == "topk":
        res.data.append(("recall", "topk-all" if full_cand else "topk-capped", nwant, nmiss))
    return nmiss


def check_threshold(res, tag, T, files, M, Msim, meas, thr, strict, known):
    """Each row must list exactly the j with M[i, j] >= thr (similarities) or < thr (distances)."""
    n = len(files)
    idx = {f: i for i, f in enumerate(files)}
    nwant = nmiss = nstrong = 0
    for i, (name, nb) in enumerate(T[:n]):
        if meas in DISTANCES:
            want = {j for j in range(n) if j != i and M.get((i, j), E.FLT_MAX) < thr}
        else:
            want = {j for j in range(n) if j != i and M.get((i, j), 0) >= thr}
        got = {idx[nm] for nm, _ in nb if nm in idx}
        extra = got - want
        res.check(not extra, "%s-extra[%d]" % (tag, i), "T=%r listed %s not qualifying (values %s)" % (
            thr, sorted(extra), {j: M[(i, j)] for j in extra}), known)
        miss = want - got
        nwant += len(want)
        nmiss += len(miss)
        # Refinement stops after 20 consecutive candidates below the threshold (src/refine.cpp), so with
        # more than 21 inputs only well-matched pairs must be found.
        strong = sorted(j for j in miss if miss_is_strong(strict and n <= 21, Msim, i, j))
        if strong:
            nstrong += len(strong)
            res.check(False, "%s-miss[%d]" % (tag, i), "T=%r qualifying %s missing (values %s, similarity %s)" % (
                thr, strong, {j: M[(i, j)] for j in strong}, {j: Msim.get((i, j)) for j in strong}), known or "lsh-self-candidate")
    if tag == "thresh":
        res.data.append(("recall", "thresh", nwant, nmiss))
    return nmiss


# ---------------------------------------------------------------------------
# Family greedy: --greedy T (LSH candidates) and --greedy TE (every cluster)

def simulate_exhaustive(M, n, meas, thr):
    """The --greedy TE rule of src/dedup_core.cpp: inputs in order join the best existing cluster within thr.

    Scores are value * mult (mult = 1 for distances, -1 for similarities) in
    float32; the smallest (score, cluster index) wins, and an input starts a
    new cluster when no cluster exists or the best score exceeds mult * thr.
    """
    mult = 1.0 if meas in DISTANCES else -1.0
    t = E.f32(mult * E.f32(thr))
    reps, members = [], []
    for i in range(n):
        best = None
        for c, r in enumerate(reps):
            v = M.get((i, r))
            v = E.FLT_MAX if v is None or math.isinf(v) else v
            key = (E.f32(v * mult), c)
            if best is None or key < best:
                best = key
        if best is None or best[0] > t:
            reps.append(i)
            members.append([])
        else:
            members[best[1]].append(i)
    return [[r] + m for r, m in zip(reps, members)]


def joins(meas, v, thr):
    """Whether value v passes the --greedy threshold (similarity >= thr, distance <= thr)."""
    if v is None:
        return False
    return v <= E.f32(thr) if meas in DISTANCES else v >= E.f32(thr)


def fam_greedy(ctx, rng, preset):
    res = ctx.res
    n = rng.randint(4, preset["maxn"])
    p = dict(k=rng.randint(12, 28), S=rng.choice([64, 128, 256, 512, 1000]), seed=rng.choice([None, rng.randint(1, 10 ** 6)]),
             canon=rng.random() < 0.85, mode=rng.choice([m[0] for m in MODES]))
    meas = rng.choice(["similarity", "similarity", "mash", "containment", "symcontain"])
    if MODE[p["mode"]][2] == "prob":
        meas = rng.choice(["similarity", "mash"])
    pt = rng.randint(2, 8)
    files, _, _ = E.clustered_inputs(rng, ctx.wd, n, p["k"], maxlen=preset["maxlen"])
    res.tags = {"family": "greedy", "mode": p["mode"], "measure": meas}
    flags = base_flags(p) + MEASURE_FLAGS[meas]
    known = "symcontain-as-distance" if meas == "symcontain" else None
    M = square_matrix(ctx, flags, files)
    if M is None:
        return
    offd = sorted({v for (i, j), v in M.items() if i != j and not is_null(meas, v) and 0 < v <= 1})
    thr = None
    if len(offd) >= 2:
        a = rng.randrange(len(offd) - 1)
        thr = E.between32(offd[a], offd[a + 1])
    if thr is None:
        thr = 0.5 if meas not in DISTANCES else 0.05
    res.desc = "greedy n=%d mode=%s measure=%s T=%r k=%d S=%d seed=%s%s -p%d" % (
        n, p["mode"], meas, thr, p["k"], p["S"], p["seed"], "" if p["canon"] else " --no-canon", pt)
    idx = {f: i for i, f in enumerate(files)}

    def clusters_of(r, tag):
        header, cl = E.parse_greedy_text(r.text)
        out = []
        for name, members in cl:
            ids = [i for _, i in members]
            res.check(all(files[i] == nm if i < n else False for nm, i in members), tag + "-labels",
                      "%s: names do not match ids %s" % (name, members))
            out.append(ids)
        return out

    for variant in ("E", ""):
        tag = "greedyE" if variant else "greedy"
        arg = repr(thr) + variant
        r1 = ctx.run(["cmp", "--greedy", arg, "-p", 1] + flags + files)
        rN = ctx.run(["cmp", "--greedy", arg, "-p", pt] + flags + files)
        rB = ctx.run(["cmp", "--greedy", arg, "-p", 1, "--binary-output"] + flags + files)
        if not (ctx.ok(r1, tag + "-run") and ctx.ok(rN, tag + "-pN-run") and ctx.ok(rB, tag + "-bin-run")):
            continue
        C = clusters_of(r1, tag)
        res.check(clusters_of(rN, tag) == C, tag + "-threads", "-p %d %s vs -p 1 %s" % (pt, clusters_of(rN, tag), C))
        try:
            Bc, used, size = E.parse_greedy_bin(rB.out)
            res.check(used == size, tag + "-bin-size", "binary clustering uses %d of %d bytes" % (used, size))
            res.check(Bc == C, tag + "-bin-vs-text", "binary %s text %s" % (Bc, C))
        except struct.error as e:
            res.fail(tag + "-bin-parse", str(e))
        flat = sorted(i for c in C for i in c)
        res.check(flat == list(range(n)), tag + "-partition", "clusters %s do not partition 0..%d" % (C, n - 1))
        for c in C:
            for m in c[1:]:
                if not joins(meas, M.get((m, c[0])), thr):
                    res.check(False, "%s-member[%d]" % (tag, m), "member %d of cluster with representative %d has %s %r, "
                              "threshold %r" % (m, c[0], meas, M.get((m, c[0])), thr), known)
        if variant:
            want = simulate_exhaustive(M, n, meas, thr)
            res.check(C == want, tag + "-rule", "clusters %s, the rule gives %s" % (C, want), known)
        else:
            # LSH greedy compares each input with at most a few candidate clusters, so a representative
            # may have an earlier representative within the threshold; count these for the statistics.
            reps = [c[0] for c in C]
            missed = sum(1 for a in range(len(reps)) for b in range(a) if joins(meas, M.get((reps[a], reps[b])), thr))
            res.data.append(("greedy-missed-merges", len(reps), missed))


# ---------------------------------------------------------------------------
# Family formats: one comparison in every output layout

SKETCH_MODES = [m[0] for m in MODES]


def draw_cmp_params(rng, preset, nmin=3, nmax=9, modes=None):
    p = dict(k=rng.randint(12, 28), S=rng.choice([64, 100, 128, 256, 512, 1024]),
             seed=rng.choice([None, rng.randint(1, 10 ** 6)]), canon=rng.random() < 0.85,
             mode=rng.choice(modes or SKETCH_MODES), n=rng.randint(nmin, min(nmax, preset["maxn"])))
    meas = rng.choice(list(MEASURE_FLAGS))
    if MODE[p["mode"]][2] == "prob":
        meas = rng.choice(["similarity", "mash"])
    return p, meas


def compare_vals(res, tag, got, want, exact=True, known=None, limit=6):
    """Compares two {(i, j): value} maps on the keys of want; records at most `limit` mismatches."""
    bad = []
    for key, v in sorted(want.items()):
        g = got.get(key)
        ok = E.same32(g, v) if exact else E.near32(g, v)
        if not ok:
            bad.append((key, g, v))
    res.check(not bad, tag, "%d of %d values differ, e.g. %s" % (len(bad), len(want), bad[:limit]), known)
    return not bad


def fam_formats(ctx, rng, preset):
    res = ctx.res
    if rng.random() < 0.3:
        return formats_large_batch(ctx, rng)
    p, meas = draw_cmp_params(rng, preset)
    n = p["n"]
    pt = rng.randint(1, 8)
    files, _, _ = E.clustered_inputs(rng, ctx.wd, n, p["k"], maxlen=preset["maxlen"])
    flags = base_flags(p) + MEASURE_FLAGS[meas] + ["-p", pt]
    nq = rng.randint(1, n)
    qidx = rng.sample(range(n), nq)
    qfile = E.write_list(ctx.wd, "queries.txt", [files[q] for q in qidx])
    use_F = rng.random() < 0.5
    pos = []
    if use_F:
        flags += ["-F", E.write_list(ctx.wd, "refs.txt", files)]
    else:
        pos = list(files)
    res.desc = "formats n=%d mode=%s measure=%s k=%d S=%d seed=%s -p%d%s nq=%d" % (
        n, p["mode"], meas, p["k"], p["S"], p["seed"], pt, " -F" if use_F else "", nq)
    res.tags = {"family": "formats", "mode": p["mode"], "measure": meas}
    sq = ctx.run(["cmp", "--square"] + flags + pos)
    if not ctx.ok(sq, "square-run"):
        return
    kind, rows, cols, M = E.parse_text_matrix(sq.text)
    res.check(kind == "square" and rows == files and cols == files, "square-header", "kind %s rows %s" % (kind, rows))
    layouts = {
        "sym": ([], {(i, j): M[(i, j)] for i in range(n) for j in range(i + 1, n)}),
        "square": (["--square"], M),
        "panel": (["-Q", qfile], {(r, c): M[(r, q)] for r in range(n) for c, q in enumerate(qidx)}),
        "phylip": (["--phylip"], {(i, j): M[(i, j)] for i in range(n) for j in range(i + 1, n)}),
    }
    for name, (extra, want) in layouts.items():
        rt = ctx.run(["cmp"] + extra + flags + pos)
        rb = ctx.run(["cmp", "--binary-output"] + extra + flags + pos)
        rf = ctx.run(["cmp", "--cmpout", "out_%s.txt" % name] + extra + flags + pos)
        rfb = ctx.run(["cmp", "--binary-output", "--cmpout", "out_%s.bin" % name] + extra + flags + pos)
        if not all(ctx.ok(r, "%s-run" % name) for r in (rt, rb, rf, rfb)):
            continue
        if name == "phylip":
            cnt, prow, got = E.parse_phylip(rt.text)
            res.check(cnt == n and prow == [f[:9].strip() for f in files], "phylip-header",
                      "count %d rows %s for %d inputs" % (cnt, prow, n))
        else:
            k2, prow, pcol, got = E.parse_text_matrix(rt.text)
            res.check(prow == files and pcol == ([files[q] for q in qidx] if name == "panel" else files),
                      "%s-labels" % name, "rows %s columns %s" % (prow, pcol))
        compare_vals(res, "%s-text-values" % name, got, want)
        with open(ctx.path("out_%s.txt" % name), "rb") as fh:
            res.check(fh.read() == rt.out, "%s-cmpout-text" % name, "--cmpout file differs from stdout")
        with open(ctx.path("out_%s.bin" % name), "rb") as fh:
            res.check(fh.read() == rb.out, "%s-cmpout-bin" % name, "--cmpout binary file differs from stdout")
        if name in ("sym", "phylip"):
            gb, nv = E.parse_condensed(rb.out, n)
            nexp = n * (n - 1) // 2
        elif name == "square":
            gb, nv = E.parse_dense(rb.out, n, n)
            nexp = n * n
        else:
            gb, nv = E.parse_dense(rb.out, n, nq)
            nexp = n * nq
        res.check(nv == nexp, "%s-bin-size" % name, "%d floats, expected %d" % (nv, nexp), "binary-tail-32k" if nv < nexp else None)
        compare_vals(res, "%s-bin-values" % name, gb, want)
    # sketch --cmpout computes the same comparison as cmp.
    rs = ctx.run(["sketch", "--square", "--cmpout", "sk_cmp.txt"] + flags + pos)
    if ctx.ok(rs, "sketch-cmpout-run"):
        with open(ctx.path("sk_cmp.txt")) as fh:
            k2, prow, pcol, got = E.parse_text_matrix(fh.read())
        res.check(prow == files, "sketch-cmpout-labels", "rows %s" % prow)
        compare_vals(res, "sketch-cmpout-values", got, M)


def formats_large_batch(ctx, rng):
    """Binary output with blocks of exactly 32768 floats (n rows times a batch of rows per block)."""
    res = ctx.res
    n, pt = rng.choice([(256, 128), (512, 64), (256, 100), (128, 256)])
    names = []
    for i in range(n):
        nm = "t%03d.fa" % i
        with open(ctx.path(nm), "w") as fh:
            fh.write(">r\n%s\n" % random_dna(rng, rng.randint(30, 80)))
        names.append(nm)
    lst = E.write_list(ctx.wd, "all.txt", names)
    flags = ["-k", 15, "-S", 16, "--full", "-F", lst, "-p", pt, "--batch-size", pt]
    res.desc = "formats large-batch n=%d -p %d --batch-size %d --square (blocks of %d floats)" % (n, pt, pt, n * min(pt, n))
    res.tags = {"family": "formats", "mode": "large-batch", "measure": "similarity"}
    rt = ctx.run(["cmp", "--square"] + flags)
    rb = ctx.run(["cmp", "--square", "--binary-output", "--cmpout", "sq.bin"] + flags)
    if not (ctx.ok(rt, "square-run") and ctx.ok(rb, "square-bin-run")):
        return
    _, rows, cols, M = E.parse_text_matrix(rt.text)
    res.check(rows == names, "square-labels", "rows differ from inputs")
    with open(ctx.path("sq.bin"), "rb") as fh:
        gb, nv = E.parse_dense(fh.read(), n, n)
    known = "binary-tail-32k" if (n * min(pt, n)) % 32768 == 0 else None
    res.check(nv == n * n, "square-bin-size", "%d floats, expected %d (repro: %s)" % (nv, n * n, E.cmdline(rb, ctx.wd)), known)
    compare_vals(res, "square-bin-values", gb, M, known=known)


# ---------------------------------------------------------------------------
# Family roundtrip: sketch files written and read back

STACK_EXT = {"oph": ".opss", "full": ".ss", "bmh": ".bmh", "pmh": ".pmh"}


def fam_roundtrip(ctx, rng, preset):
    res = ctx.res
    p, meas = draw_cmp_params(rng, preset, 2, 6, modes=list(STACK_EXT) + ["set"])
    n = p["n"]
    files, _, _ = E.clustered_inputs(rng, ctx.wd, n, p["k"], maxlen=preset["maxlen"])
    inputs = set(files)
    flags = base_flags(p) + MEASURE_FLAGS[meas]
    res.desc = "roundtrip n=%d mode=%s measure=%s k=%d S=%d seed=%s" % (n, p["mode"], meas, p["k"], p["S"], p["seed"])
    res.tags = {"family": "roundtrip", "mode": p["mode"], "measure": meas}
    before = E.snapshot(ctx.wd)
    M0 = square_matrix(ctx, flags, files, "direct")
    E.remove_new(ctx.wd, before)
    if M0 is None:
        return
    # Sketch files do not record k, which --mash-distance needs, so it is passed again.
    mf = ["-k", p["k"]] + MEASURE_FLAGS[meas]
    if p["mode"] in STACK_EXT:
        stack = "stack" + STACK_EXT[p["mode"]]
        rs = ctx.run(["sketch", "-o", stack] + flags + files)
        rc = ctx.run(["cmp", "--square", "-o", "cmp" + STACK_EXT[p["mode"]]] + flags + files)
        if ctx.ok(rs, "sketch-o") and ctx.ok(rc, "cmp-o"):
            nn_, s, cards, regs = E.parse_stacked(ctx.path(stack))
            res.check(nn_ == n and s == p["S"] and len(regs) == 8 * n * s, "stack-header",
                      "n=%d S=%d %d register bytes (expected %d inputs, S=%d)" % (nn_, s, len(regs), n, p["S"]))
            names = E.parse_names(ctx.path(stack + ".names.txt"))
            res.check([x[0] for x in names] == files and [x[1] for x in names] == cards, "stack-names",
                      "names.txt %s, header cardinalities %s" % (names, cards))
            if p["mode"] != "oph":
                # cmp densifies one-permutation sketches in place, so only the other types are byte-identical.
                with open(ctx.path("cmp" + STACK_EXT[p["mode"]]), "rb") as a, open(ctx.path(stack), "rb") as b:
                    res.check(a.read() == b.read(), "cmp-o-vs-sketch-o", "cmp -o and sketch -o stacked files differ")
            for st in (stack, "cmp" + STACK_EXT[p["mode"]]):
                rp = ctx.run(["cmp", "--presketched", "--square"] + mf + [st])
                if ctx.ok(rp, "presketched-stack"):
                    _, rows, _, got = E.parse_text_matrix(rp.text)
                    res.check(rows == files, "presketched-stack-labels", "rows %s" % rows)
                    compare_vals(res, "presketched-stack-values", got, M0)
        E.remove_new(ctx.wd, before)
    # Per-input files from --cache, compared with --presketched.
    rs = ctx.run(["sketch", "--cache"] + flags + files)
    if ctx.ok(rs, "sketch-cache"):
        made = sorted(set(os.listdir(ctx.wd)) - before)
        per = []
        for f in files:
            c = [m for m in made if m.startswith(f + ".") and (m.endswith(STACK_EXT.get(p["mode"], ".kmerset64")))]
            per.append(c[0] if len(c) == 1 else None)
        res.check(all(per), "cache-files", "one cache file per input expected, found %s" % made)
        if all(per):
            rp = ctx.run(["cmp", "--presketched", "--square"] + mf + per)
            if ctx.ok(rp, "presketched-files"):
                _, rows, _, got = E.parse_text_matrix(rp.text)
                res.check(rows == per, "presketched-files-labels", "rows %s" % rows)
                compare_vals(res, "presketched-files-values", got, M0)
            # A second run with --cache reads these files and must give the same answer.
            stamp = [os.stat(ctx.path(x)).st_mtime_ns for x in per]
            M1 = square_matrix(ctx, flags + ["--cache"], files, "cache-reuse")
            if M1 is not None:
                compare_vals(res, "cache-reuse-values", M1, M0)
                res.check([os.stat(ctx.path(x)).st_mtime_ns for x in per] == stamp, "cache-reused",
                          "cache files were rewritten instead of read")
    E.remove_new(ctx.wd, before)
    # --outprefix: cache files go into the directory, none beside the inputs.
    os.mkdir(ctx.path("op"))
    before2 = E.snapshot(ctx.wd)
    for rep in range(2):
        M2 = square_matrix(ctx, flags + ["--cache", "--outprefix", "op"], files, "outprefix-run%d" % rep)
        if M2 is not None:
            compare_vals(res, "outprefix-values%d" % rep, M2, M0)
    inop = os.listdir(ctx.path("op"))
    beside = sorted(set(os.listdir(ctx.wd)) - before2)
    res.check(len(inop) >= n, "outprefix-files", "%d files in --outprefix directory for %d inputs: %s" % (len(inop), n, inop))
    res.check(not beside, "outprefix-beside", "files written beside the inputs despite --outprefix: %s" % beside)
    E.remove_new(ctx.wd, before)


# ---------------------------------------------------------------------------
# Family measures: -Q panels, symmetric output and identities between measures

def fam_measures(ctx, rng, preset):
    res = ctx.res
    p, _ = draw_cmp_params(rng, preset, 3, 7)
    n = p["n"]
    files, _, _ = E.clustered_inputs(rng, ctx.wd, n, p["k"], maxlen=preset["maxlen"])
    extra_q = []
    if rng.random() < 0.5:
        extra_q, _, _ = E.clustered_inputs(rng, ctx.wd, rng.randint(1, 3), p["k"], maxlen=preset["maxlen"], prefix="q")
    qidx = rng.sample(range(n), rng.randint(1, n))
    queries = [files[q] for q in qidx] + extra_q
    qfile = E.write_list(ctx.wd, "queries.txt", queries)
    pt = rng.randint(1, 6)
    flags = base_flags(p) + ["-p", pt]
    meas_list = ["similarity", "mash"] if MODE[p["mode"]][2] == "prob" else list(MEASURE_FLAGS)
    res.desc = "measures n=%d queries=%d mode=%s k=%d S=%d seed=%s -p%d" % (n, len(queries), p["mode"], p["k"], p["S"], p["seed"], pt)
    res.tags = {"family": "measures", "mode": p["mode"]}
    allf = files + extra_q
    sq = {}
    for meas in meas_list:
        mf = MEASURE_FLAGS[meas]
        M = square_matrix(ctx, flags + mf, allf, "square-" + meas)
        if M is None:
            continue
        sq[meas] = M
        rp = ctx.run(["cmp", "-Q", qfile] + flags + mf + files)
        if ctx.ok(rp, "panel-run-" + meas):
            kind, rows, cols, P = E.parse_text_matrix(rp.text)
            res.check(kind == "panel" and rows == files and cols == queries, "panel-labels-" + meas,
                      "kind %s rows %s cols %s" % (kind, rows, cols))
            want = {(r, c): M[(r, allf.index(q))] for r in range(n) for c, q in enumerate(queries)}
            compare_vals(res, "panel-vs-square-" + meas, P, want)
        rs = ctx.run(["cmp"] + flags + mf + files)
        # --fastcmp fits its compression to the registers of all inputs of a run, so a run over a
        # different set of inputs is only comparable without compression.
        fitted = "--fastcmp" in MODE[p["mode"]][1] and extra_q
        if not fitted and ctx.ok(rs, "sym-run-" + meas):
            _, rows, _, S_ = E.parse_text_matrix(rs.text)
            compare_vals(res, "sym-vs-square-" + meas, S_, {(i, j): M[(i, j)] for i in range(n) for j in range(i + 1, n)})
        if meas != "containment":
            asym = [(i, j, M[(i, j)], M[(j, i)]) for i in range(len(allf)) for j in range(i) if not E.near32(M[(i, j)], M[(j, i)], 2e-6)]
            res.check(not asym, "symmetric-" + meas, "M[i,j] != M[j,i] for %s" % asym[:4])
    N = len(allf)
    if "similarity" in sq and "mash" in sq:
        bad = []
        for (i, j), s in sq["similarity"].items():
            d = sq["mash"][(i, j)]
            want = E.FLT_MAX if s <= 0 else -math.log(2 * s / (1 + s)) / p["k"]
            if not (abs(d - want) <= 1e-5 * max(abs(want), 1e-3) or (want >= E.FLT_MAX and d >= E.FLT_MAX) or abs(d - want) < 2e-7):
                bad.append(((i, j), s, d, want))
        res.check(not bad, "mash-vs-similarity", "%s" % bad[:4])
    if all(m in sq for m in ("containment", "intersection", "union")):
        card = [sq["union"][(i, i)] for i in range(N)]
        bad = []
        for i in range(N):
            for j in range(N):
                c, I = sq["containment"][(i, j)], sq["intersection"][(i, j)]
                if card[i] > 0 and abs(c * card[i] - I) > 1e-4 * max(1.0, I):
                    bad.append(((i, j), c, I, card[i]))
        res.check(not bad, "containment-identity", "containment * |row| != intersection: %s" % bad[:4])
        if "symcontain" in sq:
            bad = [((i, j), sq["symcontain"][(i, j)], sq["intersection"][(i, j)], card[i], card[j]) for i in range(N) for j in range(N)
                   if min(card[i], card[j]) > 0 and abs(sq["symcontain"][(i, j)] * min(card[i], card[j]) - sq["intersection"][(i, j)])
                   > 1e-4 * max(1.0, sq["intersection"][(i, j)])]
            res.check(not bad, "symcontain-identity", "symmetric containment * min card != intersection: %s" % bad[:4])


# ---------------------------------------------------------------------------
# Family kmers: exact k-mer files and --save-kmers/--save-kmercounts files decoded

def fam_kmers(ctx, rng, preset):
    res = ctx.res
    long_ = rng.random() < 0.3
    k = rng.randint(8, 60 if long_ else 32)
    seed = rng.choice([0, 0, rng.randint(1, 10 ** 6)])
    canon = rng.random() < 0.8
    m = rng.choice([1, 1, 1, 2, 3])
    n = rng.randint(1, 4)
    S = rng.choice([16, 64, 100, 256])
    files, _, _ = E.clustered_inputs(rng, ctx.wd, n, k, maxlen=preset["maxlen"])
    counts = [Counter({x: c for x, c in cc.items() if c >= m}) for cc in E.exact_counts(ctx.wd, files, k, canon)]
    common = ["-k", k, "--seed", seed] + ([] if canon else ["--no-canon"]) + (["-2"] if long_ else [])
    res.desc = "kmers n=%d k=%d seed=%d%s%s -m %d S=%d" % (n, k, seed, "" if canon else " --no-canon", " -2" if long_ else "", m, S)
    res.tags = {"family": "kmers", "-2": long_, "canon": canon}
    ext = ".kmerset128" if long_ else ".kmerset64"
    # --set and --countdict files: [cardinality f64][sorted ids], counts in a parallel f64 file.
    for mode in ("--set", "-J"):
        before = E.snapshot(ctx.wd)
        r = ctx.run(["sketch", mode, "-m", m] + common + files)
        if not ctx.ok(r, "sketch%s" % mode):
            continue
        made = sorted(set(os.listdir(ctx.wd)) - before)
        for i, f in enumerate(files):
            ks = [x for x in made if x.startswith(f + ".") and x.endswith(ext)]
            if not res.check(len(ks) == 1, "%s-file[%d]" % (mode, i), "k-mer file for %s not found among %s" % (f, made)):
                continue
            card = E.f64s(ctx.path(ks[0]))[0] if os.path.getsize(ctx.path(ks[0])) >= 8 else None
            if long_:
                got = E.u128s(ctx.path(ks[0]), 8)
                want = {E.kmer_id128_raw(s, seed): c for s, c in counts[i].items()}
            else:
                got = E.u64s(ctx.path(ks[0]), 8)
                want = {E.kmer_id(s, seed): c for s, c in counts[i].items()}
            res.check(sorted(got) == got, "%s-sorted[%d]" % (mode, i), "ids not sorted")
            res.check(set(got) == set(want) and len(got) == len(want), "%s-ids[%d]" % (mode, i),
                      "%d ids, %d exact k-mers, %d in common" % (len(got), len(want), len(set(got) & set(want))))
            wcard = len(want) if mode == "--set" else sum(want.values())
            res.check(card == wcard, "%s-card[%d]" % (mode, i), "header %s expected %d" % (card, wcard))
            if mode == "-J":
                cf = [x for x in made if x.startswith(f + ".") and x.endswith(".kmercounts.f64")]
                if res.check(len(cf) == 1, "-J-countfile[%d]" % i, "count file not found among %s" % made):
                    cv = E.f64s(ctx.path(cf[0]))
                    res.check(len(cv) == len(got) and all(want.get(x) == c for x, c in zip(got, cv)), "-J-counts[%d]" % i,
                              "%d counts for %d ids; mismatches %s" % (len(cv), len(got),
                                                                     [(x, c, want.get(x)) for x, c in zip(got, cv) if want.get(x) != c][:4]))
        E.remove_new(ctx.wd, before)
    # Sketch-mode k-mer databases: every sampled id is a k-mer of its input, its count the exact count.
    idc = [E.id_counts(c, seed, long_) for c in counts]
    for mname in rng.sample(["oph", "full", "bmh", "pmh"], 2):
        before = E.snapshot(ctx.wd)
        r = ctx.run(["sketch", "-N", "-S", S, "-m", m, "-o", "db"] + MODE[mname][1] + common + files)
        if not ctx.ok(r, "sketch-N-" + mname):
            continue
        db = E.parse_kmerdb(ctx.path("db.kmer64"))
        res.check((db["k"], db["sketchsize"], db["seed"], db["canon"], db["long"], db["alphabet"], db["w"]) ==
                  (k, S, seed, canon, long_, 0, k) and len(db["rows"]) == n and db["extra"] == 0, "kmerdb-header-" + mname,
                  "header %s, %d rows" % ({x: db[x] for x in ("k", "sketchsize", "seed", "canon", "long", "alphabet", "w")}, len(db["rows"])))
        names = open(ctx.path("db.kmer64.names.txt")).read().split()
        res.check(names == files, "kmerdb-names-" + mname, "names %s" % names)
        size = os.path.getsize(ctx.path("db.kmercounts.f64"))
        res.check(size == 8 * n * S, "kmercounts-f64-size-" + mname, "<out>.kmercounts.f64 holds %d bytes for %d x %d counts "
                  "(float32 rather than float64)" % (size, n, S), "stacked-kmercounts-f32" if size == 4 * n * S else None)
        with open(ctx.path("db.kmercounts.f64"), "rb") as fh:
            b = fh.read()
        cv = list(struct.unpack("<%df" % (len(b) // 4), b)) if size == 4 * n * S else list(struct.unpack("<%dd" % (len(b) // 8), b))
        for i, row in enumerate(db["rows"]):
            bad = [x for x in row if x not in idc[i]]
            empty = set(bad)
            # One-permutation sketches keep one placeholder id in empty buckets. With n distinct k-mers and S buckets
            # the expected number of empty buckets is about S * exp(-n / S), which is negligible only above about 20 * S.
            okbad = not bad or (mname == "oph" and len(empty) == 1 and len(idc[i]) < 20 * S) or (not idc[i] and len(empty) == 1)
            res.check(okbad, "kmerdb-ids-%s[%d]" % (mname, i), "%d of %d sampled ids are not k-mers of the input (%s)" % (
                len(bad), S, sorted(empty)[:3]), "bmh-stale-ids" if mname == "bmh" and not idc[i] else None)
            cbad = [(x, cv[i * S + j], idc[i][x]) for j, x in enumerate(row) if x in idc[i] and i * S + j < len(cv)
                    and cv[i * S + j] != idc[i][x]]
            res.check(not cbad, "kmerdb-counts-%s[%d]" % (mname, i), "saved counts differ from exact counts (id, saved, exact): %s" % cbad[:4],
                      "full-mincount-counts" if mname == "full" and m > 1 else None)
        E.remove_new(ctx.wd, before)


# ---------------------------------------------------------------------------
# Family contain: coverage and depth of sampled reference k-mers in query streams

def fam_contain(ctx, rng, preset):
    res = ctx.res
    win = rng.random() < 0.3
    long_ = not win and rng.random() < 0.25
    k = rng.randint(10, 60 if long_ else 32)
    w = k + rng.randint(1, 20) if win else 0
    seed = rng.choice([0, 0, rng.randint(1, 10 ** 6)])
    canon = rng.random() < 0.8
    S = rng.choice([32, 64, 128, 256])
    mname = rng.choice(["oph", "full", "full", "bmh", "pmh"])
    nr, nq = rng.randint(1, 4), rng.randint(1, 4)
    pt = rng.randint(2, 6)
    refs, _, _ = E.clustered_inputs(rng, ctx.wd, nr + nq, k, maxlen=preset["maxlen"])
    queries = refs[nr:] + rng.sample(refs[:nr], rng.randint(0, nr))
    refs = refs[:nr]
    if rng.random() < 0.25:
        # A reference without any k-mer (only records shorter than k).
        write_fastx(ctx.path("short.fa"), [random_dna(rng, k - 1)], rng, "fa", None)
        refs.insert(rng.randint(0, nr), "short.fa")
        nr += 1
    common = ["-k", k, "--seed", seed] + ([] if canon else ["--no-canon"]) + (["-2"] if long_ else []) + (["-w", w] if win else [])
    res.desc = "contain refs=%d queries=%d mode=%s k=%d%s S=%d seed=%d%s%s -p%d" % (
        nr, len(queries), mname, k, " -w %d" % w if win else "", S, seed, "" if canon else " --no-canon", " -2" if long_ else "", pt)
    res.tags = {"family": "contain", "mode": mname, "-2": long_, "window": win}
    r = ctx.run(["sketch", "-s", "-S", S, "-o", "db"] + MODE[mname][1] + common + refs)
    if not ctx.ok(r, "sketch-db"):
        return
    db = E.parse_kmerdb(ctx.path("db.kmer64"))
    if win:
        # With -w, k-mers are minimizers. Their stream is taken from sketch -G with the same options (a
        # per-input file of masked ids, checked exactly against the k-mers by the seq family); contain
        # counts a minimizer once per run of consecutive windows that select it.
        streams = {}
        for f in sorted(set(refs + queries)):
            before = E.snapshot(ctx.wd)
            rg = ctx.run(["sketch", "-G", "-o", "stream.tmp"] + common + [f])
            made = [x for x in set(os.listdir(ctx.wd)) - before if x.endswith(".mmerseq64") and x.startswith(f + ".")]
            if not (ctx.ok(rg, "minimizer-stream") and res.check(len(made) == 1, "minimizer-file", "found %s" % made)):
                return
            st = E.u64s(ctx.path(made[0]))
            E.remove_new(ctx.wd, before)
            streams[f] = Counter(x for t, x in enumerate(st) if t == 0 or st[t - 1] != x)
        rid = [streams[f] for f in refs]
        qid = [streams[f] for f in queries]
        pt = 1
    else:
        rid = [E.id_counts(c, seed, long_) for c in E.exact_counts(ctx.wd, refs, k, canon)]
        qid = [E.id_counts(c, seed, long_) for c in E.exact_counts(ctx.wd, queries, k, canon)]
    want_cov, want_dep = [], []
    for q in range(len(queries)):
        cr, dr = [], []
        for i, row in enumerate(db["rows"]):
            sampled = [x for x in row if x in rid[i]]  # placeholders of empty buckets hold no k-mer
            hit = [qid[q][x] for x in sampled if x in qid[q]]
            cr.append(len(hit) / len(sampled) if sampled and hit else 0.0)
            dr.append(E.f32(float(sum(hit)) / len(hit)) if hit else 0.0)
        want_cov.append(cr)
        want_dep.append(dr)

    def stale(i):
        # BagMinHash keeps the sampled ids of the previous input when an input has no k-mers.
        return "bmh-stale-ids" if mname == "bmh" and not rid[i] else None

    rt = ctx.run(["contain", "db.kmer64"] + queries)
    rb = ctx.run(["contain", "-b", "-p", pt, "-o", "c.bin", "db.kmer64"] + queries)
    rp = ctx.run(["contain", "-p", pt, "db.kmer64"] + queries)
    if ctx.ok(rt, "contain-run"):
        refnames, rows = E.parse_contain_text(rt.text)
        res.check(refnames == refs and [x[0] for x in rows] == queries, "contain-labels",
                  "references %s queries %s" % (refnames, [x[0] for x in rows]))
        for q, (_, vals) in enumerate(rows):
            for i, (c, d) in enumerate(vals):
                ok = abs(c - 100 * want_cov[q][i]) <= 1e-5 * 100 and abs(d - want_dep[q][i]) <= 1e-5 * max(1, want_dep[q][i])
                res.check(ok, "contain-text[%d,%d]" % (q, i), "printed %g%%:%g, exact %g%%:%g" % (
                    c, d, 100 * want_cov[q][i], want_dep[q][i]), stale(i))
    if ctx.ok(rp, "contain-pN-run"):
        res.check(rp.out == rt.out, "contain-threads", "-p %d output differs from -p 1" % pt)
    if ctx.ok(rb, "contain-bin-run"):
        with open(ctx.path("c.bin"), "rb") as fh:
            nref, nqq, cov, dep, nv = E.parse_contain_bin(fh.read())
        res.check(nref == nr and nqq == len(queries) and nv == 2 * nr * nqq, "contain-bin-header",
                  "nref %d nquery %d floats %d" % (nref, nqq, nv))
        for q in range(min(nqq, len(queries))):
            for i in range(min(nref, nr)):
                ok = E.same32(cov[q][i], want_cov[q][i]) and E.same32(dep[q][i], want_dep[q][i])
                res.check(ok, "contain-bin[%d,%d]" % (q, i), "binary %r:%r exact %r:%r" % (cov[q][i], dep[q][i], want_cov[q][i], want_dep[q][i]),
                          stale(i))


# ---------------------------------------------------------------------------
# Family seq: -G minimizer sequences and printmin

def fam_seq(ctx, rng, preset):
    res = ctx.res
    k = rng.randint(3, 31)
    canon = rng.random() < 0.7
    hpc = rng.random() < 0.5
    seed = rng.choice([0, rng.randint(1, 10 ** 6)])
    base = random_dna(rng, rng.randint(50, 600))
    recs = []
    for _ in range(rng.randint(1, 6)):
        s = mutate(rng, base, 0.02, 0.0)[:rng.randint(1, 400)]
        if rng.random() < 0.3:
            s = s[:len(s) // 2] + "A" * rng.randint(k, 2 * k + 3) + s[len(s) // 2:]
        if rng.random() < 0.3 and len(s) > 4:
            a = rng.randrange(len(s) - 2)
            s = s[:a] + "NN" + s[a + 2:]
        recs.append(s)
    write_fastx(ctx.path("s.fa"), recs, rng, "fa", None)
    flags = ["-k", k, "-G", "--seed", seed] + ([] if canon else ["--no-canon"]) + (["--hp-compress"] if hpc else [])
    res.desc = "seq records=%d k=%d seed=%d%s%s" % (len(recs), k, seed, "" if canon else " --no-canon", " --hp-compress" if hpc else "")
    res.tags = {"family": "seq", "hpc": hpc, "canon": canon}
    seqs = E.kmer_sequence(recs, k, canon)

    def compress(xs):
        return [x for t, x in enumerate(xs) if not hpc or t == 0 or xs[t - 1] != x]

    # Per record (--parse-by-seq): the stacked file holds raw k-mer encodings.
    r = ctx.run(["sketch", "--parse-by-seq", "-o", "byseq.mm"] + flags + ["s.fa"])
    if ctx.ok(r, "byseq-run"):
        st = E.parse_mmerseq_stacked(ctx.path("byseq.mm"))
        res.check(st["n"] == len(recs) and st["k"] == k and st["size"] == st["expected"] and (st["dtype"] >> 8 & 1) == int(canon),
                  "byseq-header", "n=%d k=%d dtype=%#x size %d expected %d" % (st["n"], st["k"], st["dtype"], st["size"], st["expected"]))
        for i, sq in enumerate(seqs[:st["n"]]):
            want = compress([E.encode(x) for x in sq])
            got = st["seqs"][i]
            res.check(got == want, "byseq-seq[%d]" % i, "%d items, expected %d; first difference at %s" % (
                len(got), len(want), next((t for t, (a, b) in enumerate(zip(got, want)) if a != b), min(len(got), len(want)))))
        rp = ctx.run(["printmin", "byseq.mm"])
        if ctx.ok(rp, "printmin-run"):
            lines = [l for l in rp.text.splitlines() if l.startswith("MinimizerSequence")]
            want = ["MinimizerSequence%d%s" % (i, "".join(" " + x for x in compress(sq))) for i, sq in enumerate(seqs)]
            res.check(lines == want, "printmin-text", "printmin %s expected %s" % ([l[:80] for l in lines[:3]], [w[:80] for w in want[:3]]))
        rf = ctx.run(["printmin", "-f", "byseq.mm"])
        if ctx.ok(rf, "printmin-fasta-run"):
            want = []
            for i, sq in enumerate(seqs):
                for t, x in enumerate(compress(sq)):
                    want += [">MinimizerSequence%d-Minimizer#%d" % (i, t), x]
            res.check([l for l in rf.text.splitlines() if l and not l.startswith("#")] == want, "printmin-fasta", "FASTA output differs")
    # Whole file: a per-input .mmerseq64 file of masked ids; -o is also required and must be readable by printmin.
    before = E.snapshot(ctx.wd)
    r = ctx.run(["sketch", "-o", "file.mm"] + flags + ["s.fa"])
    if ctx.ok(r, "file-run"):
        made = [x for x in set(os.listdir(ctx.wd)) - before if x.endswith(".mmerseq64")]
        if res.check(len(made) == 1, "file-mmerseq", "per-input sequence file not found: %s" % made):
            got = E.u64s(ctx.path(made[0]))
            want = compress([E.kmer_id(x, seed) for sq in seqs for x in sq])
            res.check(got == want, "file-seq", "%d ids, expected %d" % (len(got), len(want)))
        st = E.parse_mmerseq_stacked(ctx.path("file.mm"))
        tot = sum(len(compress([x for sq in seqs for x in sq])) for _ in [0])
        ok = st["size"] == st["expected"] and st["n"] == 1 and st["lens"] == [float(tot)]
        res.check(ok, "file-o-layout", "sketch -G -o in file mode wrote %d bytes, header expects %d (n=%d lengths %s)" % (
            st["size"], st["expected"], st["n"], st["lens"]), "seq-o-file-mode")
    E.remove_new(ctx.wd, before)


# ---------------------------------------------------------------------------
# Family wsketch: weighted sketches of id/weight vectors (CSR and 1-D)

WS_TYPES = {"pmh": ([], ".pmh", ["--prob"]), "bmh": (["-B"], ".bmh", ["-B"]), "fss": (["-q"], ".ss", ["--full"])}


def write_bin(path, fmt, vals):
    with open(path, "wb") as fh:
        fh.write(struct.pack("<%d%s" % (len(vals), fmt), *vals))


def fam_wsketch(ctx, rng, preset):
    res = ctx.res
    tname = rng.choice(list(WS_TYPES))
    targs, ext, cmpflag = WS_TYPES[tname]
    S = rng.choice([64, 128, 256, 512])
    nr = rng.randint(2, 6)
    big = rng.random() < 0.5
    universe = rng.sample(range(1, (1 << 40) if big else 100000), 3000)
    base = rng.sample(universe, rng.randint(20, 400))
    rows = []
    for _ in range(nr):
        keep = [x for x in base if rng.random() < rng.choice([0.3, 0.7, 0.95, 1.0])]
        extra = rng.sample(universe, rng.randint(0, 200))
        row = sorted(set(keep + extra)) or [universe[0]]
        rows.append(row)
    weights = [[float(rng.randint(1, 6)) for _ in r] for r in rows]
    ids = [x for r in rows for x in r]
    w = [x for r in weights for x in r]
    ip = [0]
    for r in rows:
        ip.append(ip[-1] + len(r))
    write_bin(ctx.path("ids.u64"), "Q", ids)
    write_bin(ctx.path("w.f64"), "d", w)
    write_bin(ctx.path("ip.u64"), "Q", ip)
    # The same data in the narrower encodings selected by -u, -f, -H and -P.
    var = []
    if not big and rng.random() < 0.5:
        write_bin(ctx.path("ids.u32"), "I", ids); var.append("-u")
    if rng.random() < 0.5:
        write_bin(ctx.path("w.f32"), "f", w); var.append("-f")
    elif rng.random() < 0.3:
        write_bin(ctx.path("w.u16"), "H", [int(x) for x in w]); var.append("-H")
    if rng.random() < 0.5:
        write_bin(ctx.path("ip.u32"), "I", ip); var.append("-P")
    res.desc = "wsketch type=%s S=%d rows=%d ids<2^%d variant=%s" % (tname, S, nr, 40 if big else 17, var)
    res.tags = {"family": "wsketch", "mode": tname}
    r = ctx.run(["wsketch", "-S", S] + targs + ["-o", "ref", "ids.u64", "w.f64", "ip.u64"])
    if not ctx.ok(r, "wsketch-csr"):
        return
    pre = "ref.sampled.%%s.stacked.%d.%d" % (nr, S)
    regpath = ctx.path((pre % "regs") + ".f64")
    n_, s_, cards, regs = E.parse_stacked(regpath)
    res.check(n_ == nr and s_ == S and len(regs) == 8 * nr * S, "ws-header", "n=%d S=%d %d register bytes" % (n_, s_, len(regs)))
    wantcard = [sum(x) for x in weights] if tname != "fss" else [float(len(x)) for x in rows]
    res.check(cards == wantcard, "ws-weights", "total weights %s expected %s" % (cards, wantcard))
    info = [float(x) for x in open(ctx.path("ref.sampled.info.txt")).read().split()]
    res.check(info == cards, "ws-info", "info.txt %s, header %s" % (info, cards))
    sid = E.u64s(ctx.path((pre % "indices") + ".i64"))
    res.check(len(sid) == nr * S, "ws-ids-size", "%d sampled ids for %d x %d" % (len(sid), nr, S))
    for i in range(nr):
        rs = set(rows[i])
        bad = [x for x in sid[i * S:(i + 1) * S] if x not in rs]
        res.check(not bad, "ws-ids[%d]" % i, "%d sampled ids are not ids of row %d: %s" % (len(bad), i, bad[:3]))
    R = [regs[i * 8 * S:(i + 1) * 8 * S] for i in range(nr)]
    eqfrac = {}
    for i in range(nr):
        for j in range(i + 1, nr):
            neq = sum(1 for t in range(S) if R[i][8 * t:8 * t + 8] == R[j][8 * t:8 * t + 8])
            wa = dict(zip(rows[i], weights[i]))
            wb = dict(zip(rows[j], weights[j]))
            J = {"pmh": E.prob_jaccard_w, "bmh": E.weighted_jaccard}.get(tname, lambda a, b: E.jaccard(a, b))(wa, wb)
            z = E.eq_z(neq, S, J)
            eqfrac[(i, j)] = neq / S
            res.check(abs(z) <= 6, "ws-similarity[%d,%d]" % (i, j), "%d of %d registers equal, exact %s %.4f, z=%.2f" % (
                neq, S, tname, J, z))
    # The stacked registers load as a sketch file of the matching type.
    os.link(regpath, ctx.path("ws" + ext))
    rp = ctx.run(["cmp", "--presketched", "--square"] + cmpflag + ["ws" + ext])
    if ctx.ok(rp, "ws-presketched"):
        _, rows_, _, P = E.parse_text_matrix(rp.text)
        bad = [(key, P.get(key), v) for key, v in eqfrac.items() if not E.near32(P.get(key), v, 2e-6)]
        res.check(not bad, "ws-presketched-values", "cmp --presketched gives %s, equal-register fractions %s" % (bad[:3], [x[2] for x in bad[:3]]))
    if var:
        a = ["wsketch", "-S", S] + targs + var + ["-o", "var", "ids.u32" if "-u" in var else "ids.u64",
                                                  "w.f32" if "-f" in var else "w.u16" if "-H" in var else "w.f64",
                                                  "ip.u32" if "-P" in var else "ip.u64"]
        rv = ctx.run(a)
        if ctx.ok(rv, "wsketch-variant"):
            for kind, sfx in (("regs", ".f64"), ("indices", ".i64"), ("hashes", ".i64")):
                with open(ctx.path(("var.sampled.%s.stacked.%d.%d" % (kind, nr, S)) + sfx), "rb") as f1, \
                        open(ctx.path((pre % kind) + sfx), "rb") as f2:
                    res.check(f1.read() == f2.read(), "ws-variant-" + kind, "%s output with %s differs from 64-bit inputs" % (kind, var))
    # 1-D input: one row as ids and weights, which must sketch exactly like that row of the matrix.
    i = rng.randrange(nr)
    write_bin(ctx.path("row.u64"), "Q", rows[i])
    write_bin(ctx.path("roww.f64"), "d", weights[i])
    r1 = ctx.run(["wsketch", "-S", S] + targs + ["-o", "one", "row.u64", "roww.f64"])
    if ctx.ok(r1, "wsketch-1d"):
        res.check(E.u64s(ctx.path("one.sampled.ids.u64")) == sid[i * S:(i + 1) * S], "ws-1d-ids",
                  "1-D sampled ids differ from CSR row %d" % i)
        with open(ctx.path("one.sampled.hashes.f64"), "rb") as fh:
            b = fh.read()
        tw = struct.unpack_from("<d", b, 0)[0]
        res.check(tw == wantcard[i], "ws-1d-weight", "total weight %r expected %r" % (tw, wantcard[i]))
        res.check(b[8:] == R[i], "ws-1d-regs", "1-D registers differ from CSR row %d" % i)
        tw_line = open(ctx.path("one.sampled.tw.txt")).read().strip()
        res.check(tw_line.endswith(";row.u64;roww.f64;d;L"), "ws-1d-twtxt", "tw.txt reads %r" % tw_line, "wsketch-tw-suffix")


# ---------------------------------------------------------------------------
# Family bed: interval sets

def fam_bed(ctx, rng, preset):
    res = ctx.res
    mname = rng.choice(["full", "oph", "bmh"])
    S = rng.choice([256, 512, 1000, 1024])
    nf = rng.randint(2, 5)
    chroms = {"chr%d" % c: rng.randint(3000, 20000) for c in range(1, rng.randint(2, 4))}
    pool = []
    for _ in range(60):
        c = rng.choice(list(chroms))
        a = rng.randrange(chroms[c])
        pool.append((c, a, min(chroms[c], a + rng.randint(1, 600))))
    files, sets, multis = [], [], []
    for f in range(nf):
        iv = [x for x in pool if rng.random() < rng.choice([0.2, 0.5, 0.9])]
        for _ in range(rng.randint(0, 30)):
            c = rng.choice(list(chroms))
            a = rng.randrange(chroms[c])
            iv.append((c, a, min(chroms[c], a + rng.randint(1, 600))))
        if rng.random() < 0.5 and iv:
            iv += rng.sample(iv, rng.randint(1, max(1, len(iv) // 4)))  # duplicates raise multiplicities
        iv = [x for x in iv if x[2] > x[1]] or [("chr1", 10, 20)]
        rng.shuffle(iv)
        name = "b%d.bed" % f
        with open(ctx.path(name), "w") as fh:
            if rng.random() < 0.3:
                fh.write("# header comment\n")
            for c, a, b in iv:
                fh.write("%s\t%d\t%d%s\n" % (c, a, b, "\tname\t0\t+" if rng.random() < 0.3 else ""))
        cnt = Counter()
        for c, a, b in iv:
            for x in range(a, b):
                cnt[(c, x)] += 1
        files.append(name)
        sets.append(set(cnt))
        multis.append(cnt)
    flags = ["--bed", "-S", S] + MODE[mname][1]
    res.desc = "bed files=%d mode=%s S=%d sizes=%s" % (nf, mname, S, [len(s) for s in sets])
    res.tags = {"family": "bed", "mode": mname}
    M = square_matrix(ctx, flags, files)
    if M is None:
        return
    for i in range(nf):
        for j in range(i + 1, nf):
            J = E.weighted_jaccard(multis[i], multis[j]) if mname == "bmh" else E.jaccard(sets[i], sets[j])
            v = M[(i, j)]
            z = (v - J) / math.sqrt(J * (1 - J) / S + 1.0 / S ** 2)
            res.check(abs(z) <= 6, "bed-similarity[%d,%d]" % (i, j), "%s %.5f exact %.5f z=%.2f" % (mname, v, J, z))
    rs = ctx.run(["sketch", "-o", "bed.stack"] + flags + files)
    if ctx.ok(rs, "bed-sketch-o"):
        n_, s_, cards, regs = E.parse_stacked(ctx.path("bed.stack"))
        res.check(n_ == nf and s_ == S and len(regs) == 8 * nf * S, "bed-stack-header", "n=%d S=%d bytes=%d" % (n_, s_, len(regs)))
        names = [x[0] for x in E.parse_names(ctx.path("bed.stack.names.txt"))]
        res.check(names == files, "bed-stack-names", "names %s" % names)
        for i in range(nf):
            if mname == "bmh":
                res.check(cards[i] == sum(multis[i].values()), "bed-card[%d]" % i, "%r expected %d" % (cards[i], sum(multis[i].values())))
            else:
                z = (cards[i] - len(sets[i])) / (len(sets[i]) / math.sqrt(S) + 1)
                res.check(abs(z) <= 6, "bed-card[%d]" % i, "%r estimate of %d, z=%.2f" % (cards[i], len(sets[i]), z))
        ext = {"full": ".ss", "oph": ".opss", "bmh": ".bmh"}[mname]
        os.link(ctx.path("bed.stack"), ctx.path("bedstack" + ext))
        rp = ctx.run(["cmp", "--presketched", "--square"] + MODE[mname][1] + ["bedstack" + ext])
        if ctx.ok(rp, "bed-presketched"):
            _, rows, _, P = E.parse_text_matrix(rp.text)
            compare_vals(res, "bed-presketched-values", P, M)


# ---------------------------------------------------------------------------

FAMILY_FN = {"nn": fam_nn, "greedy": fam_greedy, "formats": fam_formats, "roundtrip": fam_roundtrip,
             "measures": fam_measures, "kmers": fam_kmers, "contain": fam_contain, "seq": fam_seq,
             "wsketch": fam_wsketch, "bed": fam_bed}


def run_one(args, preset):
    def trial(idx, rng, wd):
        fam = FAMILIES[idx % len(FAMILIES)]
        res = TrialResult(idx, fam)
        res.tags = {"family": fam}
        ctx = Ctx(args, wd, res)
        fn = FAMILY_FN.get(fam)
        if fn is None:
            res.skips.append("family %s not implemented" % fam)
            return res
        fn(ctx, rng, preset)
        return res
    return trial


def aggregate(results):
    """Prints recall of the LSH-assisted modes across trials (informational); returns 0.

    Misses that the per-trial checks treat as failures are already counted
    there; what remains here is the documented approximation: candidates
    capped at 3.5 K for --topk, bottom-k keys in exact modes, full-register
    keys for --fastcmp and --bbit-sigs, and few candidate clusters for --greedy.
    """
    rec = defaultdict(lambda: [0, 0])
    gm = [0, 0]
    for r in results:
        for d in r.data:
            if d[0] == "recall":
                rec[d[1]][0] += d[2]
                rec[d[1]][1] += d[3]
            elif d[0] == "greedy-missed-merges":
                gm[0] += d[1]
                gm[1] += d[2]
    for kind in sorted(rec):
        want, miss = rec[kind]
        if want:
            print("recall %-12s %d of %d true neighbors listed (%.4f)" % (kind, want - miss, want, (want - miss) / want))
    if gm[0]:
        print("greedy (LSH): %d representatives, %d pairs of representatives within the threshold of each other" % tuple(gm))
    return 0


def main():
    ap = common_args(__doc__.split("\n")[0], PRESETS)
    ap.add_argument("--family", action="append", choices=FAMILIES, help="run only trials of this family (repeatable)")
    args = ap.parse_args()
    preset = PRESETS[args.preset]
    if args.family and not args.trial:
        args.trial = [i for i in range(preset["trials"]) if FAMILIES[i % len(FAMILIES)] in args.family]
    if not os.access(args.dashing2, os.X_OK):
        raise SystemExit("dashing2 binary not found or not executable: %s" % args.dashing2)
    results, el = run_trials(SUITE, SCRIPT, args, preset["trials"], run_one(args, preset))
    extra = aggregate(results)
    return summarize("suite E (workflows and outputs)", results, el, extra_fail=extra, known=KNOWN)


if __name__ == "__main__":
    sys.exit(main())
