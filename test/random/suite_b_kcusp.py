#!/usr/bin/env python3
"""Suite B: no cusp in k. Exact modes are exact at every k, and sketch error does not depend on k.

dashing2 encodes DNA k-mers directly in 64 bits for k <= 32 (128 bits with -2
for k <= 64) and switches to a rolling polynomial hash above that. Each trial
takes one k from a list that is dense around 31/32/33, 63/64/65, 127/128/129,
255/256/257 and 511/512/513 and also reaches far beyond KMC's limit of 256,
and then:
  * runs --set and -J, canonical and --no-canon, with and without -2, and
    checks similarity, intersection, union, containment and the cardinality
    against KMC3 (k <= 256) or the exact Python oracle (k > 256);
  * runs sketch modes (OPH, FullSetSketch, BagMinHash, optionally -2) with
    several seeds and records z-scores of the similarity and cardinality, which
    are then compared on both sides of every boundary (bias shift and spread).
The inputs contain homopolymer runs near length k (all-A and all-T k-mers are
the extreme 64- and 128-bit encodings), records of length k - 1, k and k + 1,
N runs, lowercase, reverse-complemented records and repeats.

Usage: DASHING2=/path/to/dashing2 python3 test/random/suite_b_kcusp.py [--preset quick|thorough] [--seed N]
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from d2rand import (KMC, MEAS_FLAGS, add_homopolymers, add_repeats, aggregate_check,
                    close, common_args, D2Error, Dashing2, decorate, k_band, kmer_counts,
                    kmer_counts_hashed, mean_sd, measure_value, mutate, pair_truth, random_dna,
                    read_fastx, run_error_known, run_trials, sketch_truth, sketch_z, SketchMode,
                    summarize, TrialResult, write_fastx)

SUITE = "b"
SCRIPT = "suite_b_kcusp.py"
BOUNDARIES = [32, 64, 128, 256, 512]


def k_list(preset):
    ks = set()
    for b in BOUNDARIES:
        ks.update(range(b - 2, b + 3))
    ks.update([1, 2, 3, 4, 7, 8, 15, 16, 17, 21, 47, 96, 200, 300, 1000, 4097, 100000])
    if preset == "thorough":
        ks.update(range(1, 80))
        ks.update(range(120, 136))
        ks.update(range(248, 264))
        ks.update([160, 192, 384, 640, 768, 1023, 1024, 1025, 2048, 10000, 1000000])
    return sorted(ks)


PRESETS = {
    # sketch replicate seeds per mode and k
    "quick": dict(reps=3),
    "thorough": dict(reps=8),
}
EXACT_MEAS = ["similarity", "intersection", "union", "containment"]
SKETCH_MODES = [SketchMode("oph", [], "set"), SketchMode("full", ["--full"], "set"),
                SketchMode("bmh", ["-B"], "weighted")]


def make_inputs(rng, k, wd):
    L = max(int(2.5 * k) + 200, rng.randint(400, 2500))
    g = random_dna(rng, L)
    g = add_repeats(rng, g, rng.randint(0, 4))
    g = add_homopolymers(rng, g, rng.randint(1, 3), k)
    seqs = [g]
    b = mutate(rng, g, rng.choice([0.002, 0.01, 0.03]), rng.choice([0, 0.002]))
    seqs.append(add_homopolymers(rng, b, 1, k))
    a = rng.randrange(len(g) // 2)
    seqs.append(g[a:a + rng.randint(len(g) // 3, len(g))])
    files = []
    for i, s in enumerate(seqs):
        feats = dict(nrec=rng.randint(1, 3), rc_prob=0.4, n_runs=rng.randint(0, 3), lower_runs=rng.randint(0, 2),
                     short_recs=rng.randint(0, 2))
        recs = decorate(rng, s, k, feats)
        # Records right at the length boundary: k - 1 (no k-mer), k (one) and k + 1 (two).
        for d in (-1, 0, 1):
            if k + d >= 1 and rng.random() < 0.5:
                st = rng.randrange(max(1, len(s) - k - 1))
                recs.append(s[st:st + k + d])
        rng.shuffle(recs)
        fn = "in_%d.fa" % i
        write_fastx(os.path.join(wd, fn), [r for r in recs if r], rng, "fa", rng.choice([None, 70]))
        files.append(fn)
    return files


def run_one(args, d2, kmc, preset, ks):
    def trial(idx, rng, wd):
        k = ks[idx]
        res = TrialResult(idx, "k=%d" % k)
        res.tags = {"k": k, "band": k_band(k)}
        files = make_inputs(rng, k, wd)
        recs = [read_fastx(os.path.join(wd, f)) for f in files]
        n = len(files)
        use_kmc = kmc is not None and k <= 256
        if not use_kmc:
            res.skips.append("k > 256 (beyond KMC): exact values from the Python oracle" if k > 256 else
                             "KMC not found: exact values from the Python oracle")
        counts = {}
        for canon in (True, False):
            if k > 2000:
                oracle = [kmer_counts_hashed(r, k, canon) for r in recs]
            else:
                oracle = [kmer_counts(r, k, canon) for r in recs]
            if use_kmc:
                tag = "c" if canon else "n"
                dbs = kmc.count_many([(f, k, canon, "%s%d" % (tag, i)) for i, f in enumerate(files)], wd)
                kc = [kmc.dump(db, wd) for db in dbs]
                for i in range(n):
                    res.check(kc[i] == oracle[i], "oracle-vs-kmc", "canon=%s file %d: KMC %d k-mers, Python %d" % (
                        canon, i, len(kc[i]), len(oracle[i])))
                for i in range(n):
                    for j in range(i + 1, n):
                        kt = kmc.pair(dbs[i], dbs[j], wd)
                        t = pair_truth(kc[i], kc[j])
                        res.check(all(kt[x] == t[x] for x in kt), "kmc_tools-vs-dump",
                                  "pair %d,%d: %s vs %s" % (i, j, kt, {x: t[x] for x in kt}))
                counts[canon] = kc
            else:
                counts[canon] = oracle
        res.desc = "k=%d |A|=%d |B|=%d |C|=%d (canonical)" % (k, len(counts[True][0]), len(counts[True][1]),
                                                              len(counts[True][2]))

        # Exact modes.
        for canon in (True, False):
            truth = {}
            for i in range(n):
                for j in range(n):
                    t = pair_truth(counts[canon][i], counts[canon][j])
                    t["k"] = k
                    truth[(i, j)] = t
            for lng in (False, True):
                for mode in ("--set", "-J"):
                    thr = 1 if rng.random() < 0.5 else rng.randint(2, 8)
                    base = ["-k", k, mode, "-p", thr, "--seed", rng.randint(1, 10 ** 6)]
                    if not canon:
                        base.append("--no-canon")
                    if lng:
                        base.append("-2")
                    cfg = "%s %s%s" % (mode, "canon" if canon else "--no-canon", " -2" if lng else "")
                    for meas in EXACT_MEAS:
                        try:
                            M = d2.cmp(base + ["--square"] + MEAS_FLAGS[meas], files, wd)
                        except D2Error as e:
                            res.check(False, "%s %s" % (cfg, meas), str(e), known=run_error_known(e))
                            continue
                        for i in range(n):
                            for j in range(n):
                                got = M.get(i, j)
                                exp = measure_value(meas, truth[(i, j)], mode == "-J")
                                res.check(close(got, exp, rel=1e-5), "%s %s[%d,%d]" % (cfg, meas, i, j),
                                          "dashing2 %s expected %s (-p %d)" % (got, exp, thr))
                    try:
                        cards = d2.cardinalities(base, files, wd)
                        for i in range(n):
                            exp = truth[(i, i)]["wA" if mode == "-J" else "A"]
                            res.check(close(cards[i], exp), "%s cardinality[%d]" % (cfg, i),
                                      "dashing2 %s expected %d (-p %d)" % (cards[i], exp, thr))
                    except D2Error as e:
                        res.fail("%s sketch -o" % cfg, str(e))

        # Sketch modes: z-scores of similarity and cardinality for the boundary statistics.
        for sm in SKETCH_MODES:
            for lng in (False, True):
                for rep in range(preset["reps"]):
                    S = rng.choice([256, 512, 1024])
                    seed = rng.randint(1, 2 ** 31 - 1)
                    base = ["-k", k, "-S", S, "--seed", seed, "-p", rng.choice([1, 4])] + sm.flags + (
                        ["-2"] if lng else [])
                    cfg = "%s%s -S %d --seed %d" % (sm.name, " -2" if lng else "", S, seed)
                    for meas in ("similarity", "union"):
                        try:
                            M = d2.cmp(base + ["--square"] + MEAS_FLAGS[meas], files, wd)
                        except D2Error as e:
                            res.check(False, "%s %s" % (cfg, meas), str(e), known=run_error_known(e))
                            continue
                        for i in range(n):
                            for j in range(n):
                                # Off-diagonal similarities, and the --union-size diagonal as the cardinality.
                                if (meas == "similarity" and i >= j) or (meas == "union" and i != j):
                                    continue
                                J, nA, nB, nI = sketch_truth(sm, counts[True][i], counts[True][j], k)
                                r = sketch_z(sm, "card" if i == j else meas, M.get(i, j), J, nA, nB, nI, S, k,
                                             diag=(i == j))
                                if r is None:
                                    continue
                                z, regular = r
                                res.check(abs(z) <= args.zmax, "%s %s[%d,%d]" % (cfg, "card" if i == j else meas, i, j),
                                          "z=%.2f estimate %s (J=%.4g |A|=%d |B|=%d)" % (z, M.get(i, j), J, nA, nB))
                                if regular:
                                    res.data.append((sm.name, lng, "card" if i == j else "similarity", k, rep, z))
        return res
    return trial


def boundary_checks(results, verbose):
    """Compares sketch z-scores for k in [b-4, b] with those for k in (b, b+4] at every boundary b.

    Groups are each sketch mode and also all modes pooled (per -2 setting and
    measure). The mean shift must stay within 3 standard errors; the spread
    must stay in 0.5-1.6 on each side when a side has at least 20 units
    (distinct k, replicate and mode).
    """
    nfail = 0
    by = {}
    for r in results:
        for (mode, lng, meas, k, rep, z) in r.data:
            for key in ((mode, lng, meas), ("all", lng, meas)):
                by.setdefault(key, []).append((k, rep, mode, z))
    print("boundary checks (sketch z-scores below vs above each boundary):")
    nchk = 0
    for (mode, lng, meas), rows in sorted(by.items()):
        label = "%-4s%-3s %-10s" % (mode, " -2" if lng else "", meas)
        ok, msg = aggregate_check([z for *_, z in rows], label + " all k", 0.15 if meas == "card" else 0.1,
                                  n_eff=len({(k, rep, md) for k, rep, md, _ in rows}))
        nfail += not ok
        nchk += 1
        if not ok or verbose:
            print("  %s %s" % ("ok  " if ok else "FAIL", msg))
        for b in BOUNDARIES:
            lo = [x for x in rows if b - 4 <= x[0] <= b]
            hi = [x for x in rows if b < x[0] <= b + 4]
            if len(lo) < 6 or len(hi) < 6:
                continue
            nlo = len({x[:3] for x in lo})
            nhi = len({x[:3] for x in hi})
            mlo, slo = mean_sd([x[3] for x in lo])
            mhi, shi = mean_sd([x[3] for x in hi])
            bound = 3.0 * math.sqrt(1.0 / nlo + 1.0 / nhi)
            ok = abs(mlo - mhi) <= bound
            for sd, nn in ((slo, nlo), (shi, nhi)):
                if nn >= 20:
                    ok = ok and 0.5 <= sd <= 1.6
            nfail += not ok
            nchk += 1
            if not ok or verbose:
                print("  %s %s k=%d: below mean %+.3f SD %.3f (n=%d), above mean %+.3f SD %.3f (n=%d), "
                      "|diff| bound %.3f" % ("ok  " if ok else "FAIL", label, b, mlo, slo, len(lo), mhi, shi,
                                              len(hi), bound))
    if not nfail:
        print("  all %d checks pass (use -v for the table)" % nchk)
    return nfail


def main():
    ap = common_args(__doc__.split("\n")[0], PRESETS)
    ap.add_argument("--zmax", type=float, default=6.0, help="per-estimate |z| bound for sketch modes (default 6)")
    args = ap.parse_args()
    preset = PRESETS[args.preset]
    d2 = Dashing2(args.dashing2)
    kmc = None if args.no_kmc else KMC(args.kmc_bin)
    if kmc is not None and not kmc.available:
        print("SKIP KMC comparisons: kmc/kmc_tools/kmc_dump not found (set KMC_BIN); using the Python oracle")
        kmc = None
    ks = k_list(args.preset)
    results, el = run_trials(SUITE, SCRIPT, args, len(ks), run_one(args, d2, kmc, preset, ks))
    nb = boundary_checks(results, args.verbose)
    if nb:
        print("boundary failures are reproduced by rerunning the whole preset: DASHING2=%s python3 %s --seed %d "
              "--preset %s" % (os.path.abspath(args.dashing2), os.path.relpath(os.path.join(
                  os.path.dirname(os.path.abspath(__file__)), SCRIPT)), args.seed, args.preset))
    return summarize("suite B (no cusp in k)", results, el, nb, known={})


if __name__ == "__main__":
    sys.exit(main())
