#!/usr/bin/env python3
"""Suite F: dashing2 does what its usage messages and README.md say.

Each trial runs one documentation check (round robin over CHECKS) on freshly
drawn random inputs and parameters: every documented option is accepted and
has its documented effect (compared with an exact oracle, with the output
without the option, or with an equivalent spelling), defaults are as
documented, documented output formats have the documented layout, and the
examples in the help and README run. Mismatches that are documentation bugs
are registered in d2rand_f.DOC_BUGS and reported as XFAIL[id].

Usage: DASHING2=/path/to/dashing2 python3 test/random/suite_f_docs.py [--preset quick|thorough] [--seed N]
"""
import math
import os
import re
import struct
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from d2rand import (D2Error, Dashing2, close, common_args, kmer_counts, measure_value, mutate, pair_truth,
                    random_dna, read_fastx, revcomp, run_trials, summarize, TrialResult, write_fastx)
from d2rand_f import (KNOWN, canonical_kmer_sequence, condensed, dedupe_runs, floats, header_options,
                      help_long_options, kmer_set, matrix, mutate_protein, parse_contain, parse_greedy,
                      parse_lists, parser_long_options, protein_kmers, random_protein, read_stacked, run,
                      same_matrix, spaced_kmers, spacing_offsets)

SUITE = "f"
SCRIPT = "suite_f_docs.py"
PRESETS = {
    # rounds: how many times each check runs (with fresh random inputs each time)
    "quick": dict(rounds=2),
    "thorough": dict(rounds=8),
}


# ---------------------------------------------------------------------------
# Inputs

def write_fa(wd, name, recs, rng=None):
    with open(os.path.join(wd, name), "w") as f:
        for i, r in enumerate(recs):
            f.write(">%s_%d\n%s\n" % (os.path.splitext(name)[0], i, r))
    return name


def related_dna(rng, wd, n, lo=2000, hi=5000, rates=(0.0, 0.08), prefix="g", nrec=(1, 3)):
    """n FASTA files derived from one random genome (SNPs/indels at random rates), split into records."""
    base = random_dna(rng, rng.randint(lo, hi))
    names = []
    for i in range(n):
        s = mutate(rng, base, rng.uniform(*rates), rng.uniform(0, rates[1] / 4))
        if rng.random() < 0.3 and len(s) > 400:
            a = rng.randrange(len(s) // 3)
            s = s[a:a + rng.randint(len(s) // 2, len(s))]
        k = rng.randint(*nrec)
        cuts = sorted(rng.sample(range(1, len(s)), k - 1)) if k > 1 else []
        recs = [s[a:b] for a, b in zip([0] + cuts, cuts + [len(s)]) if b > a]
        names.append(write_fa(wd, "%s%d.fa" % (prefix, i), recs))
    return names


def families(rng, wd, nfam, per, lo=2500, hi=4000, rates=(0.002, 0.02)):
    """nfam unrelated groups of `per` closely related genomes; returns (files, family index per file)."""
    files, fam = [], []
    for f in range(nfam):
        base = random_dna(rng, rng.randint(lo, hi))
        for j in range(per):
            s = mutate(rng, base, rng.uniform(*rates))
            files.append(write_fa(wd, "f%d_%d.fa" % (f, j), [s]))
            fam.append(f)
    return files, fam


def recs_of(wd, files):
    return [read_fastx(os.path.join(wd, f)) for f in files]


def write_list(wd, name, lines):
    with open(os.path.join(wd, name), "w") as f:
        f.write("\n".join(lines) + "\n")
    return name


# ---------------------------------------------------------------------------
# Check context

class Ctx:
    def __init__(self, d2, rng, wd, res):
        self.d2, self.rng, self.wd, self.res = d2, rng, wd, res

    def run(self, args, stdin=None, binary=False, timeout=300):
        return run(self.d2, args, self.wd, stdin=stdin, binary=binary, timeout=timeout)

    def ok(self, cond, name, detail="", known=None):
        return self.res.check(bool(cond), name, detail, known=known)

    def cmd(self, args):
        return "dashing2 " + " ".join(str(a) for a in args)

    def mat(self, args, name):
        """Runs dashing2 and parses a matrix; records a failure and returns None if that fails."""
        r = self.run(args)
        try:
            return matrix(r)
        except D2Error as e:
            self.res.fail(name, "%s: %s" % (self.cmd(args), e))
            return None

    def path(self, name):
        return os.path.join(self.wd, name)


def ret_desc(res, s):
    res.desc = res.desc + " " + s


# ---------------------------------------------------------------------------
# Checks. Each takes a Ctx and records checks on ctx.res.

DOCUMENTED_SHORT = ["k", "w", "2", "m", "F", "Q", "S", "L", "s", "N", "o", "B", "H", "J", "G", "p"]


def chk_help(c):
    """Usage messages exist for every subcommand; documented options exist; undocumented options."""
    top = c.run([])
    for sub in ("sketch", "cmp", "wsketch", "contain", "printmin"):
        c.ok(sub in top.err, "top-level-lists-%s" % sub, "dashing2 (no args) does not mention %s" % sub)
    expect = {"sketch": "dashing2 sketch <opts>", "cmp": "dashing2 cmp <opts>", "dist": "dashing2 cmp <opts>",
              "contain": "Usage: dashing2 contain", "wsketch": "Usage: dashing2 wsketch", "printmin": "Usage:"}
    texts = {}
    for sub, want in expect.items():
        for h in ("-h", "--help"):
            r = c.run([sub, h])
            c.ok(want in r.err and not r.crashed, "help[%s %s]" % (sub, h), r.brief())
            texts[(sub, h)] = r.err
    strip = lambda t: [l for l in t.splitlines()[2:] if not l.startswith("--presketched")]
    c.ok(strip(texts[("sketch", "-h")]) == strip(texts[("cmp", "-h")]), "sketch-and-cmp-share-help")
    c.ok(texts[("cmp", "-h")].splitlines()[1:] == texts[("dist", "-h")].splitlines()[1:], "dist-is-cmp")
    help_text = texts[("sketch", "-h")]
    # Every documented short option is accepted by sketch and cmp.
    files = related_dna(c.rng, c.wd, 2, 1500, 2500)
    write_list(c.wd, "L.txt", files)
    arg = {"k": ["-k", "15"], "w": ["-k", "15", "-w", "20"], "2": ["-2"], "m": ["-m", "1"], "F": ["-F", "L.txt"],
           "Q": ["-Q", "L.txt"], "S": ["-S", "64"], "L": ["-L", "6"], "s": ["-s", "-o", "o_s"],
           "N": ["-N", "-o", "o_n"], "o": ["-o", "o_o"], "B": ["-B"], "H": ["-H"], "J": ["-J"],
           "G": ["-G", "-o", "o_g"], "p": ["-p", "2"]}
    for opt in DOCUMENTED_SHORT:
        c.ok(re.search(r"(^|[\s(])-%s[/:\s]" % re.escape(opt), help_text, re.M), "short-documented[-%s]" % opt)
        for sub in ("sketch", "cmp"):
            a = [sub] + arg[opt] + (["--cmpout", "-"] if sub == "sketch" else [])
            pos = [] if opt == "F" else files
            if sub == "cmp" and opt == "G":
                continue  # cmp -G compares minimizer sequences, which needs -o; covered by chk_seq
            r = c.run(a + pos)
            c.ok(r.rc == 0, "short-accepted[%s -%s]" % (sub, opt), "%s: %s" % (c.cmd(a + pos), r.brief()))
    # Long options: named in the help but unknown to the parser, and the reverse.
    parser = parser_long_options()
    if parser is None:
        c.res.skips.append("src/options.h not found: parser option table not compared with the help")
    else:
        named = help_long_options(help_text) | help_long_options(texts[("cmp", "-h")])
        unknown = sorted(named - parser - {"presketched"})
        c.ok(not unknown, "help-options-exist", "help names options the parser lacks: %s" % unknown)
        undoc = sorted(parser - named - {"help"})
        c.ok(not undoc, "undocumented-options", "parser accepts but help never names: %s" % " ".join(
            "--" + x for x in undoc), known="undocumented-options")
    # Short options behave the same in sketch and cmp.
    rs = c.run(["sketch", "-C", "-k", "15", "--cmpout", "-"] + files)
    rc = c.run(["cmp", "-C", "-k", "15"] + files)
    c.ok((rs.rc == 0) == (rc.rc == 0), "short-opts-differ",
         "sketch -C exit %d, cmp -C exit %d" % (rs.rc, rc.rc), known="short-opts-differ")


def chk_defaults(c):
    files = related_dna(c.rng, c.wd, 3)
    r = c.run(["cmp"] + files)
    h = header_options(r.out)
    c.ok(h.get("k") == "32", "default-k", "header k=%s, help says 32 for DNA" % h.get("k"))
    c.ok(h.get("sketchsize") == "1024", "default-S", "header sketchsize=%s, help says 1024" % h.get("sketchsize"))
    c.ok(h.get("sketchtype") == "onepermsetsketch", "default-sketch-type", str(h.get("sketchtype")))
    c.ok(h.get("canon") is True, "default-canonical", "header %s" % h)
    m0 = c.mat(["cmp"] + files, "default-run")
    m1 = c.mat(["cmp", "-k", "32", "-S", "1024"] + files, "explicit-run")
    if m0 and m1:
        c.ok(same_matrix(m0, m1), "defaults-equal-explicit", "cmp vs cmp -k 32 -S 1024")
        c.ok(m0.kind == "sym", "default-output-symmetric", m0.kind)
    h2 = header_options(c.run(["cmp", "-2"] + files).out)
    c.ok(h2.get("k") == "64", "default-k-long", "-2 header k=%s" % h2.get("k"))
    # "If k is greater than this limit (31 for DNA, 14 for --protein, 22 for --protein8, 24 for --protein6)":
    # the default k is that limit (item 2: "Maximum expressible in uint64_t").
    p = [write_fa(c.wd, "p%d.fa" % i, [random_protein(c.rng, 300)]) for i in range(2)]
    for flag, lim, fl in (("", 31, files), ("--protein", 14, p), ("--protein8", 22, p), ("--protein6", 24, p)):
        hh = header_options(c.run(["cmp"] + ([flag] if flag else []) + fl).out)
        c.ok(hh.get("k") == str(lim), "help-k-limit[%s]" % (flag or "DNA"),
             "help says %d, default k is %s" % (lim, hh.get("k")), known="help-k-limits")
    # "1. Upper Triangular PHYLIP (default)": a PHYLIP matrix starts with the number of taxa.
    first = r.out.splitlines()[0] if r.out else ""
    c.ok(first.strip() == str(len(files)), "default-is-phylip",
         "first line %r; %s" % (first[:40], c.cmd(["cmp"] + files)), known="default-not-phylip")


def chk_sketch_types(c):
    files = related_dna(c.rng, c.wd, 3)
    k = c.rng.randint(12, 31)
    S = c.rng.choice([128, 256, 512, 1000])
    base = ["cmp", "-k", k, "-S", S]
    groups = [("full", ["--full", "--full-setsketch"], "fullsetsketch"),
              ("bmh", ["-B", "--multiset", "--bagminhash"], "bagminhash"),
              ("pmh", ["-P", "--prob", "--pminhash"], "probminhash"),
              ("set", ["--set", "-H"], "mmerset64"),
              ("countdict", ["-J", "--countdict"], "mmerset64,kmercountsf64")]
    ret_desc(c.res, "k=%d S=%d" % (k, S))
    counts = [kmer_counts(r, k) for r in recs_of(c.wd, files)]
    for g, flags, label in groups:
        ms = []
        for fl in flags:
            r = c.run(base + [fl] + files)
            try:
                ms.append((fl, matrix(r)))
            except D2Error as e:
                c.res.fail("run[%s]" % fl, str(e)); continue
            got = header_options(r.out).get("sketchtype")
            c.ok(got == label, "header-label[%s]" % fl, "sketchtype:%s, expected %s (%s)" % (
                got, label, c.cmd(base + [fl] + files)), known="header-sketchtype-label")
        for fl, m in ms[1:]:
            c.ok(same_matrix(ms[0][1], m), "alias-equal[%s=%s]" % (ms[0][0], fl))
        if g in ("set", "countdict") and ms:
            m = ms[0][1]
            for i in range(3):
                for j in range(i + 1, 3):
                    exp = measure_value("similarity", pair_truth(counts[i], counts[j]), g == "countdict")
                    c.ok(close(m.get(i, j), exp, 1e-5), "exact[%s][%d,%d]" % (g, i, j), "%s vs %s" % (m.get(i, j), exp))


MEASURE_ALIASES = {"similarity": [[]], "intersection": [["--intersection"], ["--intersection-size"]],
                   "union": [["--union-size"]], "containment": [["--containment"]],
                   "symcontain": [["--symmetric-containment"]],
                   "mash": [["--mash-distance"], ["--distance"], ["--poisson-distance"]]}


def chk_measures(c):
    files = related_dna(c.rng, c.wd, 3, 800, 3000, rates=(0.0, 0.1))
    k = c.rng.randint(9, 32)
    mode = c.rng.choice(["--set", "-J"])
    ret_desc(c.res, "k=%d %s" % (k, mode))
    counts = [kmer_counts(r, k) for r in recs_of(c.wd, files)]
    for meas, alias in MEASURE_ALIASES.items():
        prev = None
        for fl in alias:
            a = ["cmp", "-k", k, mode, "--square"] + fl + files
            m = c.mat(a, "run[%s]" % fl)
            if m is None:
                continue
            if prev is not None:
                c.ok(same_matrix(prev, m), "alias-equal%s" % fl)
            prev = m
            for i in range(3):
                for j in range(3):
                    t = pair_truth(counts[i], counts[j]); t["k"] = k
                    exp = measure_value(meas, t, mode == "-J")
                    c.ok(close(m.get(i, j), exp, 1e-5), "%s%s[%d,%d]" % (meas, fl, i, j),
                         "%s vs exact %s (%s)" % (m.get(i, j), exp, c.cmd(a)))


def chk_layouts(c):
    n = c.rng.randint(3, 5)
    files = related_dna(c.rng, c.wd, n)
    k = c.rng.randint(14, 31)
    base = ["cmp", "-k", k, "--containment"]
    sq = [c.mat(base + [fl] + files, "run" + fl) for fl in ("--square", "--asymmetric-all-pairs", "--asymmetric")]
    if all(sq):
        c.ok(same_matrix(sq[0], sq[1]) and same_matrix(sq[0], sq[2]), "square-aliases-equal")
        c.ok(sq[0].kind == "square" and len(sq[0].rows) == n and all(sq[0].get(i, j) is not None
             for i in range(n) for j in range(n)), "square-shape")
    write_list(c.wd, "F.txt", files)
    # "Performing --asymmetric-all-pairs with the same input for -F and -Q should yield equivalent results."
    pq = c.mat(base + ["-F", "F.txt", "-Q", "F.txt"], "run-FQ")
    if pq and sq[0]:
        c.ok(pq.kind == "panel" and all(close(pq.get(i, j), sq[0].get(i, j), 1e-6) for i in range(n) for j in range(n)),
             "F-Q-equals-square")
    # "The output shape then has |F| rows and |Q| columns"
    qs = c.rng.sample(files, c.rng.randint(1, n))
    write_list(c.wd, "Q.txt", qs)
    pm = c.mat(base + ["-F", "F.txt", "-Q", "Q.txt"], "run-panel")
    if pm:
        c.ok(len(pm.rows) == n and pm.cols == qs and all(pm.get(i, j) is not None for i in range(n)
             for j in range(len(qs))), "panel-shape", "rows %d cols %s" % (len(pm.rows), pm.cols))
    # "-F/--ffile: read paths from file in addition to positional arguments."
    write_list(c.wd, "F2.txt", files[1:])
    mf = c.mat(["cmp", "-k", k, files[0], "-F", "F2.txt"], "run-F-plus-positional")
    if mf:
        c.ok(mf.rows == files, "F-in-addition-to-positional", "rows %s" % mf.rows)
    # "you can put multiple files separated by spaces into a single line to place them all into a single sketch"
    write_list(c.wd, "J.txt", [files[0] + " " + files[1]] + files[2:])
    r = c.run(["sketch", "-k", k, "--set", "-F", "J.txt", "-o", "jo"])
    rc = recs_of(c.wd, files[:2])
    exp = len(kmer_set(rc[0] + rc[1], k))
    try:
        with open(c.path("jo.names.txt")) as f:
            got = [float(l.split("\t")[1]) for l in f if not l.startswith("#")][0]
    except (OSError, IndexError, ValueError):
        got = None
    c.ok(got == exp, "joint-line-single-sketch", "cardinality %s, |A u B| = %d" % (got, exp))


def chk_sparse(c):
    nfam = c.rng.randint(2, 3)
    per = c.rng.randint(3, 4)
    files, fam = families(c.rng, c.wd, nfam, per, rates=(0.002, 0.03))
    N = len(files)
    k = c.rng.randint(15, 25)
    base = ["cmp", "-k", k, "--full", "-S", c.rng.choice([512, 1024])]
    ret_desc(c.res, "k=%d N=%d" % (k, N))
    M = c.mat(base + files, "matrix")
    if not M:
        return
    sim = lambda i, j: M.get(min(i, j), max(i, j))
    K = per - 1
    for fl in ("--topk", "--top-k"):
        r = c.run(base + [fl, K] + files)
        L = parse_lists(r.out)
        c.ok(r.rc == 0 and len(L) == N, "topk-rows[%s]" % fl, r.brief())
        for i, f in enumerate(files):
            lst = L.get(f, [])
            mates = {files[j] for j in range(N) if fam[j] == fam[i] and j != i}
            c.ok({x for x, _ in lst} == mates, "topk-members[%s %d]" % (fl, i),
                 "%s vs family %s (%s)" % ([x for x, _ in lst], sorted(mates), c.cmd(base + [fl, K] + files)))
            vals = [v for _, v in lst]
            c.ok(vals == sorted(vals, reverse=True), "topk-sorted[%d]" % i)
            for x, v in lst:
                j = files.index(x)
                c.ok(close(v, sim(i, j), 1e-6), "topk-value[%d,%d]" % (i, j), "%s vs matrix %s" % (v, sim(i, j)))
    # "If <arg> is greater than N - 1, pairwise distances are instead emitted."
    r = c.run(base + ["--topk", N + c.rng.randint(0, 3)] + files)
    c.ok(r.out.startswith("#Dashing2 Symmetric") or r.out.startswith("#Dashing2 Asym"), "topk-large-k-is-matrix",
         "first line %r (%s)" % (r.out.splitlines()[0][:60] if r.out else "", c.cmd(base + ["--topk", N] + files)),
         known="topk-large-k")
    # "--similarity-threshold <arg>: only pairwise similarities over <arg> will be emitted"
    within = sorted({round(sim(i, j), 7) for i in range(N) for j in range(i + 1, N) if fam[i] == fam[j]})
    if len(within) >= 2:
        a = c.rng.randrange(len(within) - 1)
        T = (within[a] + within[a + 1]) / 2
    else:
        T = within[0] / 2
    r = c.run(base + ["--similarity-threshold", "%.6f" % T] + files)
    L = parse_lists(r.out)
    for i, f in enumerate(files):
        exp = {files[j] for j in range(N) if j != i and sim(i, j) > T}
        got = {x for x, _ in L.get(f, [])}
        c.ok(got == exp, "threshold[%d]" % i, "T=%.6f listed %s expected %s" % (T, sorted(got), sorted(exp)))


def chk_greedy(c):
    nfam = c.rng.randint(2, 4)
    files, fam = families(c.rng, c.wd, nfam, c.rng.randint(1, 3), rates=(0.0002, 0.002))
    t = c.rng.choice(["0.5", "0.6", "0.5E", "0.7E"])
    k = c.rng.randint(15, 25)
    ret_desc(c.res, "--greedy %s k=%d" % (t, k))
    # The expected clusters are the families only if every within-family similarity clears the
    # threshold and every cross-family one stays far below it; check that on the matrix first.
    M = c.mat(["cmp", "-k", k, "--full"] + files, "matrix")
    if not M:
        return
    th = float(t.rstrip("E"))
    n = len(files)
    sims = [(M.get(i, j), fam[i] == fam[j]) for i in range(n) for j in range(i + 1, n)]
    if any(s < th + 0.05 for s, same in sims if same) or any(s > th - 0.2 for s, same in sims if not same):
        c.res.skips.append("greedy: families not separated by the threshold in this draw")
        return
    a = ["cmp", "-k", k, "--full", "--greedy", t] + files
    r = c.run(a)
    cl = parse_greedy(r.out)
    exp = sorted(sorted(files[j] for j in range(len(files)) if fam[j] == f) for f in range(nfam))
    c.ok(sorted(sorted(x) for x in cl) == exp, "greedy-clusters", "%s vs families %s (%s)" % (cl, exp, c.cmd(a)))
    # "You may want to emit fasta-formatted output ... This is only allowed for --parse-by-seq."
    r = c.run(["cmp", "-k", k, "--greedy", "0.5F"] + files)
    c.ok(r.rc != 0, "greedy-F-requires-parse-by-seq", r.brief())
    with open(c.path("all.fa"), "w") as out:
        for f in files:
            out.write(open(c.path(f)).read())
    r = c.run(["cmp", "-k", k, "--full", "--greedy", "0.5F", "--parse-by-seq", "all.fa"])
    seqs = {s for f in files for s in read_fastx(c.path(f))}
    got = [l for l in r.out.splitlines() if l and not l.startswith(">") and not l.startswith("#")]
    c.ok(r.rc == 0 and got and all(s in seqs for s in got) and len(got) == nfam, "greedy-F-fasta",
         "%d sequences, %d families: %s" % (len(got), nfam, r.brief()))


def chk_binary(c):
    n = c.rng.randint(3, 6)
    files = related_dna(c.rng, c.wd, n)
    k = c.rng.randint(14, 31)
    base = ["cmp", "-k", k]
    m = c.mat(base + files, "text")
    sq = c.mat(base + ["--square"] + files, "text-square")
    if not (m and sq):
        return
    f32 = lambda x: struct.unpack("<f", struct.pack("<f", x))[0]
    outs = {}
    for fl in ("--binary-output", "--emit-binary", "--binary"):
        r = c.run(base + [fl] + files, binary=True)
        outs[fl] = r.out
    c.ok(len(set(outs.values())) == 1, "binary-aliases-equal")
    b = outs["--binary-output"]
    exp = condensed(m, n)
    c.ok(len(b) == 4 * len(exp) and all(close(x, f32(y), 1e-6) for x, y in zip(floats(b), exp)),
         "binary-condensed", "%d bytes for %d pairs" % (len(b), len(exp)))
    b = c.run(base + ["--square", "--binary-output"] + files, binary=True).out
    exp = [sq.get(i, j) for i in range(n) for j in range(n)]
    c.ok(len(b) == 4 * n * n and all(close(x, f32(y), 1e-6) for x, y in zip(floats(b), exp)), "binary-square")
    qs = files[:c.rng.randint(1, n)]
    write_list(c.wd, "Q.txt", qs)
    b = c.run(base + ["-Q", "Q.txt", "--binary-output"] + files, binary=True).out
    c.ok(len(b) == 4 * n * len(qs), "binary-panel", "%d bytes for %dx%d" % (len(b), n, len(qs)))
    # --cmpout aliases, and '-' for stdout (README).
    got = []
    for fl in ("--cmpout", "--distout", "--cmp-outfile"):
        c.run(["sketch", "-k", k, fl, "o" + fl] + files)
        try:
            got.append(open(c.path("o" + fl)).read())
        except OSError:
            got.append(None)
    c.ok(got[0] is not None and len(set(got)) == 1, "cmpout-aliases-equal")
    r = c.run(["sketch", "-k", k, "--cmpout", "-"] + files)
    c.ok(got[0] is not None and r.out == got[0], "cmpout-dash-is-stdout")
    # Top-k: "this emits a matrix of min(k, |N|) x |N| of IDs and distances"
    K = c.rng.randint(1, n - 1)
    b = c.run(base + ["--topk", K, "--binary-output"] + files, binary=True).out
    c.ok(len(b) == min(K, n) * n * 8, "topk-binary-size", "%d bytes, help implies %d (%s)" % (
        len(b), min(K, n) * n * 8, c.cmd(base + ["--topk", K, "--binary-output"] + files)), known="sparse-binary-format")
    # Greedy: 2 u64, (nclusters + 1) u64, nsets u64.
    r = c.run(base + ["--greedy", "0.5", "--binary-output"] + files, binary=True)
    b = r.out
    if len(b) >= 16:
        nc, ns = struct.unpack_from("<QQ", b)
        c.ok(len(b) == 16 + 8 * (nc + 1) + 8 * ns, "greedy-binary-size", "%d bytes for %d clusters %d sets, help "
             "implies %d" % (len(b), nc, ns, 16 + 8 * (nc + 1) + 8 * ns), known="greedy-binary-ids")
    else:
        c.res.fail("greedy-binary", r.brief())


def chk_stacked(c):
    files = related_dna(c.rng, c.wd, 3)
    k = c.rng.choice([c.rng.randint(12, 31), c.rng.randint(33, 60)])  # 32 is the default k
    S = c.rng.choice([64, 100, 256, 512])
    mode = c.rng.choice([[], ["--full"], ["-B"], ["--prob"]])
    ret_desc(c.res, "k=%d S=%d %s" % (k, S, mode))
    a = ["sketch", "-k", k, "-S", S] + mode + ["-o", "st"] + files
    r = c.run(a)
    try:
        n, s, cards, regs = read_stacked(c.path("st"))
        names = [l.rstrip("\n").split("\t") for l in open(c.path("st.names.txt")) if not l.startswith("#")]
    except (OSError, struct.error) as e:
        c.res.fail("stacked-file", "%s: %s %s" % (c.cmd(a), r.brief(), e)); return
    c.ok(n == 3 and s == S and len(regs) == n * S * 8, "stacked-layout", "n=%d S=%d regs=%d bytes" % (n, s, len(regs)))
    c.ok([x[0] for x in names] == files and all(close(float(x[1]), y, 1e-9) for x, y in zip(names, cards)),
         "stacked-names")
    direct = c.mat(["cmp", "-k", k, "-S", S] + mode + files, "direct")
    pre = c.mat(["cmp"] + mode + ["--presketched", "st"], "presketched")
    if direct and pre:
        c.ok(all(close(pre.get(i, j), direct.get(i, j), 1e-6) for i in range(3) for j in range(i + 1, 3)),
             "presketched-equals-direct")
    if k != 32:
        # The Mash distance depends on k, which the stacked file does not record.
        c.run(["sketch", "-k", k, "-S", S, "-o", "stm"] + files)
        d = c.mat(["cmp", "-k", k, "-S", S, "--mash-distance"] + files, "direct-mash")
        p = c.mat(["cmp", "--mash-distance", "--presketched", "stm"], "presketched-mash")
        if d and p:
            c.ok(close(p.get(0, 1), d.get(0, 1), 1e-5), "presketched-mash-distance",
                 "%s vs direct %s (k=%d; sketch -k %d -S %d -o stm, cmp --mash-distance --presketched stm)" % (
                     p.get(0, 1), d.get(0, 1), k, k, S), known="presketched-mash-k")


def chk_save_kmers(c):
    files = related_dna(c.rng, c.wd, c.rng.randint(2, 4))
    k = c.rng.randint(12, 31)
    S = c.rng.choice([16, 64, 100, 256])
    n = len(files)
    ret_desc(c.res, "k=%d S=%d" % (k, S))
    a = ["sketch", "-k", k, "-S", S, "-s", "-o", "db"] + files
    r = c.run(a)
    c.ok(r.rc == 0 and os.path.exists(c.path("db.kmer64")), "save-kmers-file", r.brief())
    if os.path.exists(c.path("db.kmer64")):
        sz = os.path.getsize(c.path("db.kmer64"))
        c.ok(sz == 16 + n * S * 8, "save-kmers-header", "%d bytes = %d + %d x %d x 8; help says a 16-byte header (%s)"
             % (sz, sz - n * S * 8, n, S, c.cmd(a)), known="save-kmers-header-size")
    c.ok(os.path.exists(c.path("db.kmer.names.txt")), "save-kmers-names",
         "files: %s" % sorted(x for x in os.listdir(c.wd) if x.startswith("db")), known="save-kmers-names-path")
    r = c.run(["sketch", "-k", k, "-S", S, "-N", "-o", "dc"] + files)
    c.ok(os.path.exists(c.path("dc.kmercounts.f64")), "save-kmercounts-file", r.brief())
    mode = c.rng.choice([[], ["--full"]])
    before = set(os.listdir(c.wd))
    a = ["sketch", "-k", k, "-S", S, "-s"] + mode + files
    c.run(a)
    new = sorted(set(os.listdir(c.wd)) - before)
    c.ok(any(".kmer" in x for x in new), "save-kmers-without-o", "%s wrote %s" % (c.cmd(a), new),
         known="save-kmers-no-file")


def chk_contain(c):
    files = related_dna(c.rng, c.wd, c.rng.randint(2, 4))
    k = c.rng.randint(13, 31)
    S = c.rng.choice([32, 64, 128])
    ret_desc(c.res, "k=%d S=%d" % (k, S))
    r = c.run(["sketch", "-k", k, "-S", S, "-s", "-o", "db"] + files)
    if r.rc:
        c.res.fail("contain-db", r.brief()); return
    a = ["contain", "db.kmer64"] + files
    r = c.run(a)
    refs, rows = parse_contain(r.out)
    c.ok(refs == files and list(rows) == files, "contain-shape", r.brief())
    for i, f in enumerate(files):
        v = rows.get(f, [])
        c.ok(len(v) == len(files) and close(v[i][0], 100.0, 1e-6), "contain-self[%d]" % i, "%s (%s)" % (v, c.cmd(a)))
    r2 = c.run(["contain", "-o", "cout.txt", "db.kmer64"] + files)
    c.ok(os.path.exists(c.path("cout.txt")) and open(c.path("cout.txt")).read() == r.out, "contain-o", r2.brief())
    write_list(c.wd, "L.txt", files)
    c.ok(c.run(["contain", "-F", "L.txt", "db.kmer64"]).out == r.out, "contain-F")
    c.ok(c.run(["contain", "-p", "3", "db.kmer64"] + files).out == r.out, "contain-threads")
    rs = c.run(["contain", "db.kmer64"], stdin=open(c.path(files[0])).read())
    c.ok(parse_contain(rs.out)[1].get("/dev/stdin") == rows.get(files[0]), "contain-stdin-default")
    b = c.run(["contain", "-b", "db.kmer64"] + files, binary=True).out
    nr, nq = struct.unpack_from("<QQ", b) if len(b) >= 16 else (0, 0)
    ok = nr == len(files) and nq == len(files) and len(b) == 16 + 8 * nr * nq
    if ok:
        cov = floats(b, 16, nr * nq)
        ok = all(close(100 * cov[q * nr + j], rows[files[q]][j][0], 1e-5) for q in range(nq) for j in range(nr))
    c.ok(ok, "contain-binary", "%d bytes" % len(b))


def chk_seq(c):
    k = c.rng.randint(3, 15)
    canon = c.rng.random() < 0.7
    recs = [random_dna(c.rng, c.rng.randint(k, 200)) for _ in range(c.rng.randint(1, 4))]
    if c.rng.random() < 0.5:
        recs.append("A" * c.rng.randint(k, 40) + "C" * c.rng.randint(k, 40))
    write_fa(c.wd, "s.fa", recs)
    flags = ["-k", k] + ([] if canon else ["--no-canon"])
    ret_desc(c.res, "k=%d canon=%s nrec=%d" % (k, canon, len(recs)))
    a = ["sketch"] + flags + ["--seq", "--parse-by-seq", "-o", "mm", "s.fa"]
    r = c.run(a)
    t = c.run(["printmin", "mm"])
    rows = [l.split()[1:] for l in t.out.splitlines()]
    exp = [canonical_kmer_sequence(s, k, canon) for s in recs]
    c.ok(rows == exp, "printmin-equals-kmers", "%s / printmin mm: %s" % (c.cmd(a), t.brief()))
    fa = c.run(["printmin", "-f", "mm"]).out.splitlines()
    c.ok([x for x in fa if not x.startswith(">")] == [x for row in exp for x in row], "printmin-fasta")
    c.run(["printmin", "-o", "pm.txt", "mm"])
    c.ok(os.path.exists(c.path("pm.txt")) and open(c.path("pm.txt")).read() == t.out, "printmin-o")
    c.run(["sketch"] + flags + ["--seq", "--hp-compress", "--parse-by-seq", "-o", "mh", "s.fa"])
    rh = [l.split()[1:] for l in c.run(["printmin", "mh"]).out.splitlines()]
    c.ok(rh == [dedupe_runs(x) for x in exp], "hp-compress")
    # "header: [uint64_t nitems, uint32_t k, uint32_t w], followed by `nitems` [double]"
    if os.path.exists(c.path("mm")):
        sz = os.path.getsize(c.path("mm"))
        body = 8 * len(recs) + 8 * sum(len(x) for x in exp)
        c.ok(sz == 16 + body, "seq-header-size", "%d bytes = %d-byte header + lengths + minimizers (%s)" % (
            sz, sz - body, c.cmd(a)), known="seq-header")


def chk_filterset(c):
    files = related_dna(c.rng, c.wd, 3, rates=(0.0, 0.05))
    k = c.rng.randint(9, 32)
    long = c.rng.random() < 0.3
    fl = ["-k", k, "--set"] + (["-2"] if long else [])
    ret_desc(c.res, "k=%d%s" % (k, " -2" if long else ""))
    R = recs_of(c.wd, files)
    S = [kmer_set(r, k) for r in R]
    a = ["sketch"] + fl + ["--filterset", files[2], "-o", "fo"] + files[:2]
    c.run(a)
    try:
        got = [float(l.split("\t")[1]) for l in open(c.path("fo.names.txt")) if not l.startswith("#")]
    except OSError:
        got = []
    exp = [len(S[0] - S[2]), len(S[1] - S[2])]
    c.ok(got == exp, "filterset-cardinality", "%s vs %s (%s)" % (got, exp, c.cmd(a)))
    m = c.mat(["cmp"] + fl + ["--intersection", "--filterset", files[2]] + files[:2], "filterset-cmp")
    if m:
        e = len((S[0] & S[1]) - S[2])
        c.ok(close(m.get(0, 1), e), "filterset-intersection", "%s vs %d" % (m.get(0, 1), e))


def card_of(c, args, name):
    a = ["sketch"] + args + ["-o", name]
    c.run(a)
    try:
        return [float(l.split("\t")[1]) for l in open(c.path(name + ".names.txt")) if not l.startswith("#")]
    except OSError:
        return None


def chk_window_downsample(c):
    files = related_dna(c.rng, c.wd, 2, 4000, 8000)
    k = c.rng.randint(11, 31)
    full = card_of(c, ["-k", k, "--set"] + files, "w0")
    exp = [len(kmer_set(r, k)) for r in recs_of(c.wd, files)]
    c.ok(full == exp, "set-cardinality", "%s vs %s" % (full, exp))
    # "-w/--window-size ... [Default: k-mer length]"
    c.ok(card_of(c, ["-k", k, "-w", k, "--set"] + files, "w1") == full, "window-default-is-k")
    W = k + c.rng.randint(5, 30)
    wc = card_of(c, ["-k", k, "-w", W, "--set"] + files, "w2")
    dens = 2.0 / (W - k + 2)
    c.ok(wc and all(0.5 * dens * e < x < 1.6 * dens * e for x, e in zip(wc, exp)), "window-density",
         "-w %d: %s for %s k-mers (random-minimizer density %.3f)" % (W, wc, exp, dens))
    f = round(c.rng.uniform(0.1, 0.9), 3)
    dc = card_of(c, ["-k", k, "--downsample", f, "--set"] + files, "w3")
    ok = dc is not None and all(abs(x - f * e) <= 6 * math.sqrt(e * f * (1 - f)) + 1 for x, e in zip(dc, exp))
    c.ok(ok, "downsample-fraction", "--downsample %s: %s of %s" % (f, dc, exp))


def chk_spacing(c):
    nk = c.rng.randint(4, 12)
    gaps = [c.rng.choice([0, 0, 1, 2, 3]) for _ in range(nk - 1)]
    spec = ",".join(map(str, gaps))
    # Run-length form: equal neighbours written as AxN.
    rle, i = [], 0
    while i < len(gaps):
        j = i
        while j < len(gaps) and gaps[j] == gaps[i]:
            j += 1
        rle.append("%dx%d" % (gaps[i], j - i) if j - i > 1 else str(gaps[i]))
        i = j
    rle = ",".join(rle)
    files = related_dna(c.rng, c.wd, 2, 1000, 3000)
    ret_desc(c.res, "k=%d --spacing %s (%s)" % (nk, spec, rle))
    offs = spacing_offsets(spec)
    c.ok(spacing_offsets(rle) == offs, "oracle-rle")
    R = recs_of(c.wd, files)
    S = [spaced_kmers(r, offs) for r in R]
    got = card_of(c, ["-k", nk, "--spacing", spec, "--set"] + files, "sp")
    c.ok(got == [len(x) for x in S], "spaced-cardinality", "%s vs %s" % (got, [len(x) for x in S]))
    m = c.mat(["cmp", "-k", nk, "--spacing", spec, "--set", "--intersection"] + files, "spaced-cmp")
    if m:
        c.ok(close(m.get(0, 1), len(S[0] & S[1])), "spaced-intersection", "%s vs %d" % (m.get(0, 1), len(S[0] & S[1])))
    c.ok(card_of(c, ["-k", nk, "--spacing", rle, "--set"] + files, "sp2") == got, "spacing-rle-equivalent")
    c.ok(card_of(c, ["-k", nk, "--spacing", spec, "--no-canon", "--set"] + files, "sp3") == got,
         "spacing-disables-canon")


def chk_protein(c):
    base = random_protein(c.rng, c.rng.randint(300, 1500))
    files = [write_fa(c.wd, "p%d.fa" % i, [mutate_protein(c.rng, base, c.rng.uniform(0, 0.1))]) for i in range(2)]
    k = c.rng.choice([c.rng.randint(3, 14), c.rng.randint(15, 24)])
    ret_desc(c.res, "k=%d" % k)
    R = recs_of(c.wd, files)
    S = [protein_kmers(r, k) for r in R]
    prev = None
    for fl in ("--protein", "--protein20", "--enable-protein"):
        got = card_of(c, [fl, "-k", k, "--set"] + files, "pr")
        c.ok(got == [len(x) for x in S], "protein-cardinality[%s]" % fl, "%s vs %s" % (got, [len(x) for x in S]))
        m = c.mat(["cmp", fl, "-k", k, "--set", "--intersection"] + files, "protein-cmp")
        if m:
            c.ok(close(m.get(0, 1), len(S[0] & S[1])), "protein-intersection[%s]" % fl,
                 "%s vs %d" % (m.get(0, 1), len(S[0] & S[1])))
            if prev:
                c.ok(same_matrix(prev, m), "protein-aliases-equal")
            prev = m


def chk_countmin(c):
    files = related_dna(c.rng, c.wd, 3)
    k = c.rng.randint(12, 31)
    cs = c.rng.choice([16, 256, 4096])
    ret_desc(c.res, "k=%d -c %d" % (k, cs))
    base = ["cmp", "-k", k]
    m0 = c.mat(base + files, "oph")
    r = c.run(base + ["--countmin-size", cs] + files)
    lines = r.out.splitlines()
    pre = lines[:next((i for i, l in enumerate(lines) if l.startswith("#Sources")), len(lines))]
    c.ok(all(l.startswith("#") for l in pre), "countmin-header-lines", "lines before #Sources: %s" % pre,
         known="countsketch-header-newline")
    # Remove the stray header line to read the values.
    clean = "\n".join(l for l in lines if not l.startswith(";"))
    try:
        from d2rand import parse_matrix
        m1 = parse_matrix(clean)
    except D2Error:
        m1 = None
    c.ok(m0 and m1 and all(close(m0.get(0, j), m1.get(0, j), 0) for j in (1, 2)), "countmin-no-effect-on-oph")
    # "This is only relevant to WeightedSetSketch and DiscreteProbabilitySetSketch."
    mode = c.rng.choice(["--set", "-J"])
    e = c.mat(base + [mode] + files, "exact")
    r = c.run(base + [mode, "-c", cs] + files)
    try:
        x = parse_matrix("\n".join(l for l in r.out.splitlines() if not l.startswith(";")))
    except D2Error:
        x = None
    c.ok(e and x and all(close(e.get(0, j), x.get(0, j), 1e-6) for j in (1, 2)), "countmin-exact-unchanged",
         "%s: %s vs %s without -c" % (c.cmd(base + [mode, "-c", cs] + files), x and [x.get(0, 1), x.get(0, 2)],
                                      e and [e.get(0, 1), e.get(0, 2)]), known="countmin-exact-modes")
    b0 = c.mat(base + ["-B"] + files, "bmh")
    try:
        b1 = parse_matrix("\n".join(l for l in c.run(base + ["-B", "-c", 16] + files).out.splitlines()
                                    if not l.startswith(";")))
    except D2Error:
        b1 = None
    c.ok(b0 and b1 and not same_matrix(b0, b1), "countmin-affects-bagminhash")


def chk_cache(c):
    files = related_dna(c.rng, c.wd, 2)
    k = c.rng.randint(12, 31)
    mode = c.rng.choice([[], ["--full"], ["-B"]])
    opt = c.rng.choice(["--outprefix", "--prefix"])
    cflag = c.rng.choice(["--cache", "--cache-sketches"])
    os.mkdir(c.path("cdir"))
    before = set(os.listdir(c.wd))
    a = ["cmp", "-k", k] + mode + [cflag, opt, "cdir"] + files
    m1 = c.mat(a, "cache-1")
    new_here = sorted(set(os.listdir(c.wd)) - before)
    cached = os.listdir(c.path("cdir"))
    c.ok(len(cached) == len(files) and not new_here, "cache-in-outprefix", "cdir %s, next to inputs %s (%s)" % (
        cached, new_here, c.cmd(a)))
    m2 = c.mat(a, "cache-2")
    m3 = c.mat(["cmp", "-k", k] + mode + files, "nocache")
    c.ok(m1 and m2 and m3 and same_matrix(m1, m2) and same_matrix(m1, m3), "cache-reuse-equal")


def chk_sizes(c):
    files = related_dna(c.rng, c.wd, 2)
    k = c.rng.randint(12, 31)
    s = c.rng.choice([c.rng.randint(16, 3000), 1 << c.rng.randint(4, 11)])
    h = header_options(c.run(["cmp", "-k", k, "-S", s] + files).out)
    c.ok(h.get("sketchsize") == str(s), "S-sets-sketchsize", "-S %d header %s" % (s, h.get("sketchsize")))
    l2 = c.rng.randint(4, 11)
    mL = c.mat(["cmp", "-k", k, "-L", l2] + files, "L")
    mS = c.mat(["cmp", "-k", k, "-S", 1 << l2] + files, "S")
    c.ok(mL and mS and same_matrix(mL, mS), "L-equals-S-power")
    for bad in (0, 64):
        c.ok(c.run(["cmp", "-k", k, "-L", bad] + files).rc != 0, "L-bounds[%d]" % bad)
    nb = c.rng.choice([1, 2, 4])
    ms = [c.mat(["cmp", "-k", k, fl, nb] + files, "fastcmp" + fl) for fl in ("--fastcmp", "--regsize", "--regbytes")]
    c.ok(all(ms) and same_matrix(ms[0], ms[1]) and same_matrix(ms[0], ms[2]), "fastcmp-aliases")
    c.ok(c.run(["cmp", "-k", k, "--fastcmp", 3] + files).rc != 0, "fastcmp-rejects-3")
    presets = {"--fastcmp-bytes": ("20,1.2", 1), "--fastcmp-shorts": (".06,1.0005", 2),
               "--fastcmp-words": ("19.77,1.0000000109723500835", 4)}
    p, (ab, n) = c.rng.choice(sorted(presets.items()))
    m1 = c.mat(["cmp", "-k", k, "--full", p] + files, "preset")
    m2 = c.mat(["cmp", "-k", k, "--full", "--setsketch-ab", ab, "--fastcmp", n] + files, "preset-ab")
    c.ok(m1 and m2 and same_matrix(m1, m2), "preset-equals-ab[%s]" % p)
    r = c.run(["cmp", "-k", k, p] + files)
    c.ok(r.rc == 0, "preset-with-default-sketch[%s]" % p, "%s: %s" % (c.cmd(["cmp", "-k", k, p] + files), r.brief()),
         known="fastcmp-presets-need-full")


def chk_seqs_in_ram(c):
    files = related_dna(c.rng, c.wd, 1, nrec=(2, 4))
    a = ["sketch", "-v", "-v", "-v", "--seqs-in-ram", "--parse-by-seq", "--cmpout", "t.tbl", "-k", 15] + files
    r = c.run(a)
    c.ok(r.rc == 0 and "Swapping to keep sequences in RAM" not in r.err, "seqs-in-ram-has-effect",
         "%s logs the !seqs_in_memory branch" % c.cmd(a), known="seqs-in-ram-ignored")


def chk_readme(c):
    files = related_dna(c.rng, c.wd, 3)
    k = c.rng.randint(15, 31)
    write_list(c.wd, "F.txt", files)
    # Use 1
    r1 = c.run(["sketch", "-k", k, "--cmpout", "u1.txt"] + files)
    r1b = c.run(["sketch", "-k", k, "--cmpout", "u1b.txt", "-F", "F.txt"])
    m = c.mat(["cmp", "-k", k] + files, "cmp")
    try:
        from d2rand import parse_matrix
        u1 = parse_matrix(open(c.path("u1.txt")).read())
        u1b = parse_matrix(open(c.path("u1b.txt")).read())
    except (OSError, D2Error):
        u1 = u1b = None
    c.ok(u1 and u1b and m and same_matrix(u1, m) and same_matrix(u1b, m), "readme-use1", r1.brief() + r1b.brief())
    # Use 2 and 3
    for fl, v in (("--topk", 2), ("--similarity-threshold", 0.1)):
        c.run(["sketch", "-k", k, "--cmpout", "u2.txt", fl, v, "-F", "F.txt"])
        try:
            L = parse_lists(open(c.path("u2.txt")).read())
        except (OSError, ValueError):
            L = {}
        c.ok(list(L) == files, "readme-%s" % fl)
    # Use 7: '-o input_sequence_set.topk.tsv' is expected to hold the top-k table.
    a = ["sketch", "--cache-sketches", "-p8", "-F", "F.txt", "-k31", "-S1024", c.rng.choice(["--set", "--countdict"]),
         "--topk", "25", "-o", "input_sequence_set.topk.tsv"]
    r = c.run(a)
    try:
        txt = open(c.path("input_sequence_set.topk.tsv"), errors="replace").read()
    except OSError:
        txt = ""
    c.ok(r.out.strip() or txt.startswith("#Collection"), "readme-use7",
         "%s: stdout empty and the .tsv holds %d bytes of stacked sketches" % (c.cmd(a), len(txt)),
         known="readme-use7-no-output")
    # Use 4: protein, by sequence.
    base = random_protein(c.rng, 600)
    write_fa(c.wd, "uniref50.fa", [mutate_protein(c.rng, base, 0.05) for _ in range(3)])
    fl = c.rng.choice(["--protein", "--protein14", "--protein8", "--protein6"])
    r = c.run(["sketch", "-S256", "--cmpout", "prot.k5.table", "--parse-by-seq", "-k5", fl, "uniref50.fa"])
    try:
        pm = parse_matrix(open(c.path("prot.k5.table")).read())
        ok = pm.rows == ["uniref50_%d" % i for i in range(3)]
    except (OSError, D2Error):
        ok = False
    c.ok(ok, "readme-use4[%s]" % fl, r.brief())
    # "Canonicalization is off by default."
    recs = read_fastx(c.path(files[0]))
    write_fa(c.wd, "rc.fa", [revcomp(x) for x in recs])
    m = c.mat(["cmp", "-k", k, "--set", files[0], "rc.fa"], "rc")
    c.ok(m and m.get(0, 1) < 0.5, "readme-canon-off-by-default",
         "similarity of a file and its reverse complement is %s" % (m and m.get(0, 1)), known="readme-canon-default")


def chk_edit_distance(c):
    files = related_dna(c.rng, c.wd, 3, 1000, 3000, prefix="e")
    with open(c.path("input.fasta"), "w") as out:
        for f in files:
            out.write(open(c.path(f)).read())
    a = ["sketch", "-p8", "--cmpout", "knn.edit-distance.tbl", "-k7", "--parse-by-seq", "--edit-distance",
         "--compute-edit-distance", "input.fasta"]
    r = c.run(a)
    c.ok(r.rc == 0, "readme-use6", "%s: %s" % (c.cmd(a), r.brief()), known="readme-edit-distance-crash")
    outs = []
    for _ in range(3):
        rr = c.run(["sketch", "-k7", "--parse-by-seq", "--edit-distance", "--cmpout", "-", "input.fasta"])
        outs.append(rr.out)
    c.ok(len(set(outs)) == 1, "edit-distance-deterministic", "3 runs of sketch -k7 --parse-by-seq --edit-distance "
         "--cmpout - input.fasta gave %d different outputs" % len(set(outs)), known="edit-distance-nondeterministic")


def chk_wsketch(c):
    n = c.rng.randint(20, 300)
    ids = c.rng.sample(range(1, 1 << 40), n)
    w = [c.rng.randint(1, 50) for _ in range(n)]
    nrow = c.rng.randint(1, 4)
    cuts = sorted(c.rng.sample(range(1, n), nrow - 1)) if nrow > 1 else []
    ip = [0] + cuts + [n]
    open(c.path("ids.u64"), "wb").write(struct.pack("<%dQ" % n, *ids))
    open(c.path("ids.u32"), "wb").write(struct.pack("<%dI" % n, *[x & 0xffffffff for x in ids]))
    open(c.path("ids32.u64"), "wb").write(struct.pack("<%dQ" % n, *[x & 0xffffffff for x in ids]))
    open(c.path("w.f64"), "wb").write(struct.pack("<%dd" % n, *w))
    open(c.path("w.f32"), "wb").write(struct.pack("<%df" % n, *w))
    open(c.path("w.u16"), "wb").write(struct.pack("<%dH" % n, *w))
    open(c.path("w.u32"), "wb").write(struct.pack("<%dI" % n, *w))
    open(c.path("ip.u64"), "wb").write(struct.pack("<%dQ" % len(ip), *ip))
    open(c.path("ip.u32"), "wb").write(struct.pack("<%dI" % len(ip), *ip))
    S = c.rng.choice([16, 64, 100])
    ret_desc(c.res, "n=%d rows=%d S=%d" % (n, nrow, S))

    def tw(args):
        r = c.run(["wsketch", "-S", S, "-o", "o"] + args)
        try:
            txt = open(c.path("o.sampled.tw.txt"), "rb").read().decode("latin-1")
            os.remove(c.path("o.sampled.tw.txt"))
        except OSError:
            return None, r
        m = re.match(r"Total weight: ([0-9.eE+-]+);", txt)
        return (float(m.group(1)) if m else None, txt), r

    for flags, wf in (([], "w.f64"), (["-f"], "w.f32"), (["-H"], "w.u16"), (["-U"], "w.u32")):
        got, r = tw(flags + ["ids.u64", wf])
        val = got[0] if got else None
        c.ok(val is not None and close(val, sum(w)), "wsketch-1d-weight%s" % flags,
             "total weight %s, sum %d (dashing2 wsketch %s ids.u64 %s)" % (val, sum(w), " ".join(flags), wf),
             known="wsketch-U-1d" if flags == ["-U"] else None)
        if flags == [] and got:
            c.ok(re.search(r";[fdH];[WL]\n$", got[1]) is not None, "wsketch-tw-format",
                 "tw.txt is %r" % got[1][-30:], known="wsketch-tw-garbage")
    got, _ = tw(["ids.u64"])
    c.ok(got and close(got[0], n), "wsketch-1d-unweighted", str(got and got[0]))
    for flags, wf, ipf in (([], "w.f64", "ip.u64"), (["-f"], "w.f32", "ip.u64"), (["-H"], "w.u16", "ip.u64"),
                           (["-U"], "w.u32", "ip.u64"), (["-P"], "w.f64", "ip.u32"), ([], "-", "ip.u64")):
        for x in os.listdir(c.wd):
            if x.startswith("o.sampled"):
                os.remove(c.path(x))
        r = c.run(["wsketch", "-S", S, "-o", "o"] + flags + ["ids.u64", wf, ipf])
        try:
            rows = [float(x) for x in open(c.path("o.sampled.info.txt")).read().split()]
        except (OSError, ValueError):
            rows = None
        exp = [float(sum(w[a:b]) if wf != "-" else b - a) for a, b in zip(ip, ip[1:])]
        c.ok(rows is not None and len(rows) == nrow and all(close(x, y) for x, y in zip(rows, exp)),
             "wsketch-csr-row-weights%s[%s]" % (flags, wf), "%s vs %s: %s" % (rows, exp, r.brief()))
    # -u reads 32-bit ids: the same ids stored in 64 bits give the same sketch.
    c.run(["wsketch", "-S", S, "-o", "u", "-u", "ids.u32", "w.f64"])
    c.run(["wsketch", "-S", S, "-o", "v", "ids32.u64", "w.f64"])
    try:
        same = open(c.path("u.sampled.hashes.f64"), "rb").read() == open(c.path("v.sampled.hashes.f64"), "rb").read()
    except OSError:
        same = False
    c.ok(same, "wsketch-u32-ids")
    # Help example 5: "dashing2 wsketch -S 64 -o g1.k31.k64.mat - fq.fastq.k31.indices64 fq.fastq.k31.indptr64"
    r = c.run(["wsketch", "-S", 64, "-o", "ex5", "-", "ids.u64", "ip.u64"])
    c.ok(r.rc == 0 and os.path.exists(c.path("ex5.sampled.info.txt")), "wsketch-example5",
         "dashing2 wsketch -S 64 -o ex5 - ids.u64 ip.u64: %s" % r.brief(), known="wsketch-example-order")
    # Help example 2: "dashing2 wsketch -S 64 -o g1.k31.k64 g1.fastq.k31.kmerset64 # sketches k-mers only"
    files = related_dna(c.rng, c.wd, 1, 500, 1500)
    k = c.rng.randint(12, 31)
    c.run(["sketch", "-k", k, "--set", "--cache"] + files)
    ks = [x for x in os.listdir(c.wd) if x.endswith(".kmerset64")]
    if ks:
        nk = len(kmer_set(read_fastx(c.path(files[0])), k))
        got, r = tw([ks[0]])
        c.ok(got and close(got[0], nk), "wsketch-example2", "total weight %s for %d k-mers (dashing2 wsketch -S %d -o o "
             "%s)" % (got and got[0], nk, S, ks[0]), known="wsketch-example-kmerset-header")
    else:
        c.res.fail("kmerset", "no .kmerset64 written")


CHECKS = [chk_help, chk_defaults, chk_sketch_types, chk_measures, chk_layouts, chk_sparse, chk_greedy, chk_binary,
          chk_stacked, chk_save_kmers, chk_contain, chk_seq, chk_filterset, chk_window_downsample, chk_spacing,
          chk_protein, chk_countmin, chk_cache, chk_sizes, chk_seqs_in_ram, chk_readme, chk_edit_distance,
          chk_wsketch]


def main():
    ap = common_args(__doc__.split("\n")[0], PRESETS)
    ap.add_argument("--check", action="append", help="run only this check (name without chk_; repeatable)")
    args = ap.parse_args()
    d2 = Dashing2(args.dashing2)
    checks = CHECKS
    if args.check:
        checks = [f for f in CHECKS if f.__name__[4:] in args.check]
        args.repro_extra = sum((["--check", x] for x in args.check), [])
    ntr = PRESETS[args.preset]["rounds"] * len(checks)

    def trial(idx, rng, wd):
        fn = checks[idx % len(checks)]
        res = TrialResult(idx, fn.__name__[4:])
        res.tags = {"check": fn.__name__[4:]}
        fn(Ctx(d2, rng, wd, res))
        return res

    results, el = run_trials(SUITE, SCRIPT, args, ntr, trial)
    return summarize("suite F (documentation)", results, el, known=KNOWN)


if __name__ == "__main__":
    sys.exit(main())
