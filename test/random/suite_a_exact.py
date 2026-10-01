#!/usr/bin/env python3
"""Suite A: dashing2's exact modes (--set and -J/--countdict) must match KMC3 exactly.

Each trial draws random inputs and options, counts k-mers with KMC (kmc,
kmc_dump, kmc_tools simple intersect/union) and checks every measure dashing2
reports (similarity, intersection, union, containment, symmetric containment,
Mash distance and the per-input cardinality) against the values derived from
the KMC databases. The Python k-mer counter is checked against KMC as well, and
replaces KMC when KMC is not installed.

Usage: DASHING2=/path/to/dashing2 python3 test/random/suite_a_exact.py [--preset quick|thorough] [--seed N]
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from d2rand import (KMC, MEAS_FLAGS, ROW_LOSS, add_repeats, close, common_args, D2Error, Dashing2,
                    decorate, k_band, kmer_counts, measure_value, mutate, pair_truth, random_dna,
                    read_fastx, run_error_known, run_trials, summarize, threshold, TrialResult,
                    write_fastx)

SUITE = "a"
SCRIPT = "suite_a_exact.py"
PRESETS = {
    # trials, max base genome length
    "quick": dict(trials=48, maxlen=4000),
    "thorough": dict(trials=400, maxlen=20000),
}


def draw_k(rng):
    r = rng.random()
    if r < 0.12:
        return rng.randint(1, 12)
    if r < 0.45:
        return rng.randint(13, 32)
    if r < 0.70:
        return rng.randint(33, 64)
    if r < 0.85:
        return rng.randint(65, 128)
    return rng.randint(129, 256)


def draw_params(rng, preset):
    p = {}
    p["k"] = draw_k(rng)
    p["canon"] = rng.random() < 0.7
    p["mode"] = rng.choice(["--set", "-J"])
    p["long"] = rng.random() < 0.3
    p["nfiles"] = rng.randint(2, 4)
    p["layout"] = rng.choices(["sym", "square", "panel"], [0.4, 0.35, 0.25])[0]
    p["threads"] = 1 if rng.random() < 0.4 else rng.randint(2, 8)
    p["seed"] = rng.randint(1, 10 ** 6) if rng.random() < 0.6 else None
    p["m"] = rng.choice([2, 3]) if rng.random() < 0.15 else 1
    p["ffile"] = rng.random() < 0.2
    p["base_len"] = rng.randint(300, preset["maxlen"])
    p["files"] = []
    for i in range(p["nfiles"]):
        f = dict(
            snp=rng.choice([0, 0.001, 0.01, 0.05]),
            indel=rng.choice([0, 0.001, 0.005]),
            slice=rng.random() < 0.4,
            nrep=rng.randint(0, 5),
            nrec=rng.randint(1, 6),
            short_recs=rng.randint(0, 3),
            rc_prob=rng.choice([0, 0.3, 0.5, 1.0]),
            n_runs=rng.randint(0, 5),
            lower_runs=rng.randint(0, 3),
            fmt="fq" if rng.random() < 0.2 else "fa",
            gz=rng.random() < 0.2,
            width=rng.choice([None, 60, 13]),
            empty=rng.random() < 0.04,
        )
        p["files"].append(f)
    return p


def describe(p):
    fl = ";".join("%s%s snp=%g indel=%g rep=%d rec=%d short=%d rc=%g N=%d low=%d%s" % (
        f["fmt"], ".gz" if f["gz"] else "", f["snp"], f["indel"], f["nrep"], f["nrec"], f["short_recs"],
        f["rc_prob"], f["n_runs"], f["lower_runs"], " empty" if f["empty"] else "") for f in p["files"])
    return "k=%d %s %s%s layout=%s -p%d seed=%s m=%d%s base=%d files=[%s]" % (
        p["k"], p["mode"], "canon" if p["canon"] else "--no-canon", " -2" if p["long"] else "", p["layout"],
        p["threads"], p["seed"], p["m"], " -F" if p["ffile"] else "", p["base_len"], fl)


def make_inputs(rng, p, wd):
    base = random_dna(rng, p["base_len"])
    names = []
    for i, f in enumerate(p["files"]):
        s = mutate(rng, base, f["snp"], f["indel"])
        if f["slice"] and len(s) > 50:
            a = rng.randrange(len(s) // 2)
            s = s[a:a + rng.randint(len(s) // 4, len(s))]
        s = add_repeats(rng, s, f["nrep"])
        recs = decorate(rng, s, p["k"], f)
        if f["empty"]:
            recs = [r[:max(1, p["k"] - 1)] for r in recs[:2]] if p["k"] > 1 else ["N"]
        name = "in_%d.%s%s" % (i, f["fmt"], ".gz" if f["gz"] else "")
        write_fastx(os.path.join(wd, name), recs, rng, f["fmt"], f["width"])
        names.append(name)
    return names


def run_one(args, d2, kmc, preset):
    def trial(idx, rng, wd):
        p = draw_params(rng, preset)
        res = TrialResult(idx, describe(p))
        res.tags = {"k": k_band(p["k"], p["long"]), "mode": p["mode"], "canon": p["canon"], "-2": p["long"],
                    "threads": "p1" if p["threads"] == 1 else "p>1", "layout": p["layout"], "m": p["m"]}
        files = make_inputs(rng, p, wd)
        k, canon, m = p["k"], p["canon"], p["m"]
        weighted = p["mode"] == "-J"

        # Exact counts: KMC if available, checked against the Python oracle.
        oracle = [threshold(kmer_counts(read_fastx(os.path.join(wd, f)), k, canon), m) for f in files]
        if kmc is not None:
            dbs = kmc.count_many([(f, k, canon, "db%d" % i) for i, f in enumerate(files)], wd)
            counts = [kmc.dump(db, wd, ci=m) for db in dbs]
            for i in range(len(files)):
                res.check(counts[i] == oracle[i], "oracle-vs-kmc",
                          "file %d: KMC %d k-mers, Python %d" % (i, len(counts[i]), len(oracle[i])))
        else:
            res.skips.append("KMC not found: exact values come from the Python oracle only")
            counts = oracle
        n = len(files)
        truth = {}
        for i in range(n):
            for j in range(n):
                t = pair_truth(counts[i], counts[j])
                t["k"] = k
                truth[(i, j)] = t
        if kmc is not None:
            for i in range(n):
                for j in range(i + 1, n):
                    kt = kmc.pair(dbs[i], dbs[j], wd, ci=m)
                    t = truth[(i, j)]
                    res.check(all(kt[x] == t[x] for x in ("I", "U", "wI", "wU")), "kmc_tools-vs-dump",
                              "pair %d,%d kmc_tools %s dumps %s" % (i, j, kt, {x: t[x] for x in kt}))

        base = ["-k", k, p["mode"], "-p", p["threads"]]
        if not canon:
            base.append("--no-canon")
        if p["long"]:
            base.append("-2")
        if p["seed"] is not None:
            base += ["--seed", p["seed"]]
        if m > 1:
            base += ["-m", m]
        pos = files
        if p["ffile"]:
            with open(os.path.join(wd, "list_F.txt"), "w") as fh:
                fh.write("\n".join(files) + "\n")
            pos = []
            base += ["-F", "list_F.txt"]

        # Layout: which (row, column) pairs exist and how indices map to files.
        qfile = None
        if p["layout"] == "panel":
            nq = rng.randint(1, n)
            qidx = rng.sample(range(n), nq)
            ridx = list(range(n))
            with open(os.path.join(wd, "list_Q.txt"), "w") as fh:
                fh.write("\n".join(files[q] for q in qidx) + "\n")
            qfile = "list_Q.txt"
            cells = [(r, c, ridx[r], qidx[c]) for r in range(n) for c in range(nq)]
        elif p["layout"] == "square":
            base.append("--square")
            cells = [(i, j, i, j) for i in range(n) for j in range(n)]
        else:
            cells = [(i, j, i, j) for i in range(n) for j in range(i + 1, n)]

        for meas, flags in MEAS_FLAGS.items():
            try:
                M = d2.cmp(base + flags, pos, wd, qfile=qfile, rows=files)
            except D2Error as e:
                res.check(False, "run %s" % meas, str(e), known=run_error_known(e))
                continue
            for (r, c, fi, fj) in cells:
                got = M.get(r, c)
                exp = measure_value(meas, truth[(fi, fj)], weighted)
                res.check(close(got, exp, rel=1e-5), "%s[%d,%d]" % (meas, fi, fj),
                          "dashing2 %s expected %s (A=%d B=%d I=%d U=%d%s)" % (
                              got, exp, truth[(fi, fj)]["A"], truth[(fi, fj)]["B"], truth[(fi, fj)]["I"],
                              truth[(fi, fj)]["U"], " wA=%d wB=%d wI=%d" % (
                                  truth[(fi, fj)]["wA"], truth[(fi, fj)]["wB"], truth[(fi, fj)]["wI"]) if weighted else ""))
        # Cardinality from `dashing2 sketch -o`: distinct k-mers for --set, total count for -J.
        try:
            sk = [a for a in base if a != "--square"]
            cards = d2.cardinalities(sk, pos, wd)
            for i in range(n):
                exp = truth[(i, i)]["wA" if weighted else "A"]
                res.check(len(cards) == n and close(cards[i], exp), "cardinality[%d]" % i,
                          "dashing2 %s expected %d" % (cards[i] if i < len(cards) else None, exp))
        except D2Error as e:
            res.fail("sketch -o", str(e))
        return res
    return trial


def main():
    ap = common_args(__doc__.split("\n")[0], PRESETS)
    args = ap.parse_args()
    preset = PRESETS[args.preset]
    d2 = Dashing2(args.dashing2)
    kmc = None if args.no_kmc else KMC(args.kmc_bin)
    if kmc is not None and not kmc.available:
        print("SKIP KMC comparisons: kmc/kmc_tools/kmc_dump not found (set KMC_BIN); using the Python oracle")
        kmc = None
    results, el = run_trials(SUITE, SCRIPT, args, preset["trials"], run_one(args, d2, kmc, preset))
    return summarize("suite A (exact modes vs %s)" % ("KMC3" if kmc else "Python oracle"), results, el,
                     known=dict([ROW_LOSS]))


if __name__ == "__main__":
    sys.exit(main())
