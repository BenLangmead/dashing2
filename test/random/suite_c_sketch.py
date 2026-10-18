#!/usr/bin/env python3
"""Suite C: dashing2's sketch estimates must agree with the exact values within their standard errors.

Each trial picks one sketch mode (round robin over OPH, FullSetSketch,
BagMinHash, ProbMinHash and their --fastcmp, --bbit-sigs and preset
SetSketch-(a,b) compressed variants, optionally with -2), random inputs, k, -S,
--seed and -p, and compares every printed estimate (similarity, intersection,
union, containment, symmetric containment, Mash distance, and the cardinality
from the diagonal of --union-size --square and from `sketch -o`) with the exact
value computed by the Python oracle. See README.md for the error model.

Checks: per estimate |z| <= --zmax (default 6); deterministic identities between
the measures dashing2 prints; and, over all trials, the mean and SD of z per
mode and measure (bias and spread checks).

Usage: DASHING2=/path/to/dashing2 python3 test/random/suite_c_sketch.py [--preset quick|thorough|calibrate] [--seed N]
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from d2rand import (MEAS_FLAGS, add_repeats, aggregate_check, close, common_args, D2Error,
                    Dashing2, decorate, j_from_mash, k_band, kmer_counts, mean_sd, mutate,
                    random_dna, read_fastx, run_error_known, run_trials, sketch_modes, sketch_truth,
                    sketch_z, summarize, TrialResult, write_fastx)

SUITE = "c"
SCRIPT = "suite_c_sketch.py"
PRESETS = {
    "quick": dict(per_mode=10, maxlen=30000),
    "thorough": dict(per_mode=80, maxlen=100000),
    "calibrate": dict(per_mode=200, maxlen=40000),
}
SIZES = [128, 200, 256, 300, 500, 512, 1000, 1001, 1024, 1500, 2048, 4096]
# Known issues in the code under test. Failures that match one are reported as
# XFAIL (or as failures with --strict). README.md has the details.
KNOWN = {}
# Allowed |mean z| on top of 3 / sqrt(trials). Cardinality estimators carry
# an O(1/m) upward bias (they invert a sum of m minima), which reaches the
# union and intersection, and symmetric containment divides by the smaller of
# two noisy cardinalities (Jensen). The calibration in README.md motivates the values.
BIAS_ALLOW = {"card": 0.15, "union": 0.15, "intersection": 0.15, "symcontain": 0.25}
DEFAULT_ALLOW = 0.1


def draw_k(rng):
    r = rng.random()
    if r < 0.6:
        return rng.randint(12, 32)
    if r < 0.85:
        return rng.randint(33, 64)
    return rng.randint(65, 160)


def make_inputs(rng, p, wd, mode):
    L = p["L"]
    g = random_dna(rng, L)
    if mode.kind != "set" or rng.random() < 0.3:
        g = add_repeats(rng, g, rng.randint(1, 12), unit=(5, 200))
    b = mutate(rng, g, rng.choice([0, 0.002, 0.01, 0.03, 0.08]), rng.choice([0, 0.001, 0.003]))
    if rng.random() < 0.4:
        a = rng.randrange(len(b) // 2)
        b = b[a:a + rng.randint(len(b) // 4, len(b))]
    if rng.random() < 0.3:
        b += random_dna(rng, rng.randint(1, L))
    r = rng.random()
    if r < 0.4:
        a = rng.randrange(len(g) // 2)
        c = g[a:a + rng.randint(len(g) // 10 + 1, len(g))]
    elif r < 0.7:
        c = random_dna(rng, rng.randint(100, L))
    else:
        c = mutate(rng, g, 0.05, 0.005)
    files = []
    for i, s in enumerate((g, b, c)):
        feats = dict(nrec=rng.randint(1, 4), rc_prob=rng.choice([0, 0.5]), n_runs=rng.randint(0, 3),
                     lower_runs=rng.randint(0, 2), short_recs=rng.randint(0, 2))
        recs = decorate(rng, s, p["k"], feats)
        fn = "in_%d.fa" % i
        write_fastx(os.path.join(wd, fn), recs, rng, "fa", rng.choice([None, 80]))
        files.append(fn)
    return files


def run_one(args, d2, preset, modes, zmax):
    def trial(idx, rng, wd):
        mode = modes[idx % len(modes)]
        p = dict(k=draw_k(rng), S=rng.choice(SIZES), seed=rng.randint(1, 2 ** 31 - 1),
                 threads=1 if rng.random() < 0.5 else rng.randint(2, 8), long=rng.random() < 0.25,
                 canon=rng.random() < 0.85,
                 L=int(math.exp(rng.uniform(math.log(300), math.log(preset["maxlen"])))))
        desc = "mode=%s flags=%s k=%d -S %d --seed %d -p%d%s%s L=%d" % (
            mode.name, " ".join(mode.flags), p["k"], p["S"], p["seed"], p["threads"], " -2" if p["long"] else "",
            "" if p["canon"] else " --no-canon", p["L"])
        res = TrialResult(idx, desc)
        res.tags = {"mode": mode.name, "k": k_band(p["k"], p["long"]), "-2": p["long"],
                    "threads": "p1" if p["threads"] == 1 else "p>1", "canon": p["canon"]}
        files = make_inputs(rng, p, wd, mode)
        k, m = p["k"], p["S"]
        counts = [kmer_counts(read_fastx(os.path.join(wd, f)), k, p["canon"]) for f in files]
        base = ["-k", k, "-S", m, "--seed", p["seed"], "-p", p["threads"]] + mode.flags
        if p["long"]:
            base.append("-2")
        if not p["canon"]:
            base.append("--no-canon")
        n = len(files)
        truth = {(i, j): sketch_truth(mode, counts[i], counts[j], k) for i in range(n) for j in range(n)}
        meas_list = list(mode.measures())  # for set and weighted modes the union diagonal gives |A|
        printed = {}
        for meas in meas_list:
            try:
                M = d2.cmp(base + ["--square"] + MEAS_FLAGS[meas], files, wd)
            except D2Error as e:
                res.check(False, "run %s" % meas, str(e), known=run_error_known(e))
                continue
            cbase = M.compression_base()
            for i in range(n):
                for j in range(n):
                    est = M.get(i, j)
                    printed[(meas, i, j)] = est
                    J, nA, nB, nI = truth[(i, j)]
                    if i == j:
                        if meas != "union" or mode.kind == "prob":
                            continue
                        r = sketch_z(mode, "card", est, J, nA, nB, nI, m, k, diag=True)
                        key = "card"
                    else:
                        r = sketch_z(mode, meas, est, J, nA, nB, nI, m, k, base=cbase)
                        key = meas
                    if r is None:
                        continue
                    z, regular = r
                    ok = res.check(abs(z) <= zmax, "%s[%d,%d]" % (key, i, j),
                                   "z=%.2f estimate %s exact %s (J=%.4g |A|=%g |B|=%g |AnB|=%g m=%d%s)" % (
                                       z, est, nA if i == j else (J if meas in ("similarity", "mash") else None),
                                       J, nA, nB, nI, m, " b=%g" % cbase if cbase else ""))
                    if not (i > j and meas in ("similarity", "intersection", "union", "mash", "symcontain")):
                        # Symmetric measures are recorded once per unordered pair.
                        res.data.append((mode.name, key, z, regular, k, p["long"]))
        # Deterministic identities between the printed measures.
        if mode.kind != "prob":
            for i in range(n):
                for j in range(n):
                    if i == j:
                        continue
                    g = lambda ms: printed.get((ms, i, j))
                    J, I, U = g("similarity"), g("intersection"), g("union")
                    cA = printed.get(("union", i, i))
                    cB = printed.get(("union", j, j))
                    if None in (J, I, U, cA, cB) or truth[(i, i)][1] == 0 or truth[(j, j)][1] == 0:
                        # With an empty operand the ratio measures are defined as 0 while the
                        # compressed intersection may hold a stray register match.
                        continue
                    if U > 0:
                        res.check(close(I / U, J, rel=2e-4, abs_=2e-6), "identity sim=I/U[%d,%d]" % (i, j),
                                  "similarity %s but intersection/union %s/%s" % (J, I, U))
                    res.check(close(I + U, cA + cB, rel=2e-4), "identity I+U=|A|+|B|[%d,%d]" % (i, j),
                              "I+U=%s but |A|+|B|=%s" % (I + U, cA + cB))
                    C = g("containment")
                    if C is not None and cA > 0:
                        res.check(close(C, I / cA, rel=2e-4, abs_=2e-6), "identity C=I/|A|[%d,%d]" % (i, j),
                                  "containment %s but I/|row| = %s" % (C, I / cA))
                    SC = g("symcontain")
                    if SC is not None and min(cA, cB) > 0:
                        res.check(close(SC, I / min(cA, cB), rel=2e-4, abs_=2e-6),
                                  "identity SC=I/min[%d,%d]" % (i, j), "symmetric containment %s but I/min = %s" % (
                                      SC, I / min(cA, cB)))
                    d = g("mash")
                    if d is not None and J > 0:
                        res.check(close(j_from_mash(d, k), J, rel=2e-4, abs_=2e-6), "identity mash(J)[%d,%d]" % (i, j),
                                  "mash %s gives J=%s, similarity %s" % (d, j_from_mash(d, k), J))
        # Thread independence: with a fixed --seed, -p N must print what -p 1 prints.
        if p["threads"] > 1:
            meas = mode.measures()[0]
            pb = [x for x in base]
            pb[pb.index("-p") + 1] = 1
            try:
                M1 = d2.cmp(pb + ["--square"] + MEAS_FLAGS[meas], files, wd)
                for i in range(n):
                    for j in range(n):
                        a, b = printed.get((meas, i, j)), M1.get(i, j)
                        res.check(a == b, "-p%d=-p1 %s[%d,%d]" % (p["threads"], meas, i, j),
                                  "-p %d printed %s, -p 1 printed %s" % (p["threads"], a, b))
            except D2Error as e:
                res.check(False, "run -p 1", str(e), known=run_error_known(e))
        # Cardinalities written by `dashing2 sketch -o` (names.txt).
        try:
            cards = d2.cardinalities(base, files, wd)
            for i in range(n):
                J, nA, nB, nI = truth[(i, i)]
                diag = printed.get(("union", i, i))
                if mode.kind == "set" and diag is not None:
                    res.check(close(cards[i], diag, rel=2e-5), "names.txt=diag[%d]" % i,
                              "names.txt %s, --union-size diagonal %s" % (cards[i], diag))
                r = sketch_z(mode, "card", cards[i], J, nA, nB, nI, m, k, diag=True)
                res.check(abs(r[0]) <= zmax, "names.txt card[%d]" % i, "z=%.2f estimate %s exact %s" % (
                    r[0], cards[i], nA))
        except D2Error as e:
            res.check(False, "sketch -o", str(e))
        return res
    return trial


def aggregates(results, args, verbose):
    """Bias and spread checks over all trials; returns the number of failing groups."""
    groups = {}
    for r in results:
        for (mode, key, z, regular, k, lng) in r.data:
            if regular and math.isfinite(z):
                groups.setdefault((mode, key), {}).setdefault(r.idx, []).append(z)
    nfail = 0
    lines = []
    bykey = {}
    for (mode, key), per in sorted(groups.items()):
        zs = [z for v in per.values() for z in v]
        bykey.setdefault(key, {}).update({(mode, t): v for t, v in per.items()})
        # One bias check per mode and measure (well over 100 groups per run) needs a stricter
        # bound than the pooled check below to keep chance failures rare.
        ok, msg = aggregate_check(zs, "%-15s %-12s" % (mode, key), BIAS_ALLOW.get(key, DEFAULT_ALLOW),
                                  n_eff=len(per), z_crit=4.0)
        if not ok:
            nfail += 1
        lines.append((ok, msg))
    print("aggregate checks (regular pairs only):")
    for ok, msg in lines:
        if not ok or verbose:
            print("  %s %s" % ("ok  " if ok else "FAIL", msg))
    for key, per in sorted(bykey.items()):
        zs = [z for v in per.values() for z in v]
        ok, msg = aggregate_check(zs, "%-15s %-12s" % ("all modes", key), BIAS_ALLOW.get(key, DEFAULT_ALLOW),
                                  n_eff=len(per))
        nfail += not ok
        print("  %s %s" % ("ok  " if ok else "FAIL", msg))
    return nfail, groups


def calibration_table(results):
    """Per mode and measure: trials, z count, mean z, SD z, max |z| (all pairs and regular pairs)."""
    groups = {}
    for r in results:
        for (mode, key, z, regular, k, lng) in r.data:
            g = groups.setdefault((mode, key), dict(all=[], reg=[], trials=set()))
            g["all"].append(z)
            if regular:
                g["reg"].append(z)
                g["trials"].add(r.idx)
    print("calibration: mode measure | trials n mean_z sd_z max|z| "
          "(regular pairs) | n max|z| (all pairs)")
    for (mode, key), g in sorted(groups.items()):
        mu, sd = mean_sd(g["reg"])
        print("  %-15s %-12s | %4d %5d %+.3f %.3f %5.2f | %5d %5.2f" % (
            mode, key, len(g["trials"]), len(g["reg"]), mu, sd, max([abs(z) for z in g["reg"]] or [0]),
            len(g["all"]), max(abs(z) for z in g["all"])))


def main():
    ap = common_args(__doc__.split("\n")[0], PRESETS)
    ap.add_argument("--zmax", type=float, default=6.0, help="per-estimate |z| bound (default 6)")
    ap.add_argument("--modes", default=None, help="comma-separated subset of sketch modes")
    args = ap.parse_args()
    preset = PRESETS[args.preset]
    d2 = Dashing2(args.dashing2)
    modes = sketch_modes()
    if args.modes:
        want = args.modes.split(",")
        modes = [m for m in modes if m.name in want]
        args.repro_extra = ["--modes", args.modes]
    if args.zmax != 6.0:
        args.repro_extra = getattr(args, "repro_extra", []) + ["--zmax", str(args.zmax)]
    ntr = preset["per_mode"] * len(modes)
    results, el = run_trials(SUITE, SCRIPT, args, ntr, run_one(args, d2, preset, modes, args.zmax))
    if args.preset == "calibrate":
        calibration_table(results)
    nagg, _ = aggregates(results, args, args.verbose)
    if nagg:
        print("aggregate failures are reproduced by rerunning the whole preset: " +
              "DASHING2=%s python3 %s --seed %d --preset %s" % (os.path.abspath(args.dashing2),
                                                                 os.path.relpath(os.path.join(os.path.dirname(
                                                                     os.path.abspath(__file__)), SCRIPT)),
                                                                 args.seed, args.preset))
    return summarize("suite C (sketch estimates vs exact)", results, el, nagg, KNOWN)


if __name__ == "__main__":
    sys.exit(main())
