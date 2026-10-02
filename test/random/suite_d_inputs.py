#!/usr/bin/env python3
"""Suite D: dashing2's input-processing and k-mer-generation modes against exact oracles.

Each trial belongs to one family (chosen round robin by trial index) and draws
its own inputs and options:

  byseq      --parse-by-seq rows equal one file per record, for every sketch type and measure
  minimizer  -w minimizer sets and -G streams (direct and rolling paths, 64/128 bit, canonical or not)
  entmin     --entmin picks: exact oracle plus entropy, strand and fallback properties
  spacing    --spacing exact spaced-k-mer sets and streams, run-length syntax, rejections
  protein    --protein, --protein14, --protein8, --protein6 (and aliases) exact sets
  seqmode    -G streams, by-seq stream files, Hamming and --compute-edit-distance comparisons
  filter     --filterset (sequence and binary forms) on direct and rolling paths
  downsample --downsample: exact hash oracle, kept fraction and Jaccard preservation
  hpseed     --hp-compress and --seed invariance
  inputs     FASTA/FASTQ, gzip/bzip2/xz/zstd, -F/-Q lists with joint entries, checked against KMC3

Usage: DASHING2=/path/to/dashing2 python3 test/random/suite_d_inputs.py [--preset quick|thorough] [--seed N]
"""
import math
import os
import random
import shutil
import struct
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from d2rand import (KMC, MEAS_FLAGS, D2Error, Dashing2, TrialResult, add_homopolymers, add_repeats, close,
                    common_args, decorate, kmer_counts, measure_value, mutate, pair_truth, parse_matrix,
                    random_dna, revcomp, run_trials, summarize)
from d2rand_d import (AMINO20, DNA, M64, M128, PROTEIN_ALIASES, PROTEIN_ALPHABETS, Mask, distinct_kmer_strings,
                      downsample_hash, downsample_threshold, ent_score, entropy, have_zstd, hp_collapse,
                      kmer_items, levenshtein, lex64, lex128, n_valid_kmers, offsets_from_spacing,
                      read_byseq_stream_file, read_kmerset, read_stream, select_window, window_picks,
                      write_records)

SUITE = "d"
SCRIPT = "suite_d_inputs.py"
FAMILIES = ["byseq", "minimizer", "entmin", "spacing", "protein", "seqmode", "filter", "downsample", "hpseed",
            "inputs"]
PRESETS = {
    # trials per family, maximum base length
    "quick": dict(per_family=20, maxlen=2500),
    "thorough": dict(per_family=100, maxlen=6000),
}

KNOWN = {
    "uncanon-window-polyT": "with --no-canon (or --spacing) and k filling the k-mer word (k = 32, or 64 with -2), "
                            "a window whose minimizer is poly-T emits nothing: the all-ones encoding doubles as the "
                            "'window not full' marker (bonsai encoder.h for_each_uncanon_unspaced_windowed, "
                            "for_each_uncanon_spaced; qmap.h next_value)",
    "spaced-window-invalid": "--spacing with -w: k-mers that overlap an invalid character enter the window as the "
                             "all-ones marker with a real score, so windows count positions rather than valid k-mers "
                             "and a window won by the marker emits nothing; records shorter than the window emit "
                             "nothing (bonsai encoder.h next_minimizer and for_each_uncanon_spaced have no validity "
                             "check and no partial-window tail)",
    "byseq-exact-card": "--parse-by-seq with OPH or --full replaces the estimated cardinality by an exact count when "
                        "it is below 10 x S (src/fastxsketchbyseq.cpp, 'exact counting fall-back'), so intersection, "
                        "union, containment and symmetric containment differ from the same record sketched as a file",
    "filterset-windowed": "with -w, --filterset removes only the filter file's own minimizers rather than every k-mer "
                          "in it, so input minimizers that occur in the filter file as non-minimizers are kept "
                          "(src/d2.cpp Dashing2Options::filterset builds the set with the windowed encoder)",
    "spacing-long-filename": "--spacing with irregular gaps and large k: exact modes (--set, -J, -G) and --cache "
                             "abort with 'Failed to open' because the output file name spells out the whole seed "
                             "and exceeds the 255-byte file name limit (src/fastxmerge.cpp:94-95 makedest appends "
                             "the spacing string to the name)",
    "byseq-seq-stale-window": "sketch -G --parse-by-seq on the rolling path without canonicalization: a record "
                              "shorter than k gets the previous record's last window minimum (RollingHasher::"
                              "for_each_uncanon returns before resetting its window, and the by-seq fallback for "
                              "short records reads that stale window; src/fastxsketchbyseq.cpp, bonsai encoder.h)",
}


# ---------------------------------------------------------------------------
# Running dashing2 in a trial directory


class Runner:
    """Runs dashing2 in a directory and collects (then removes) every file it writes there, recursively."""

    def __init__(self, d2, wd):
        self.d2, self.wd = d2, wd
        self.calls = 0

    def _snapshot(self):
        out = set()
        for root, _, files in os.walk(self.wd):
            for f in files:
                out.add(os.path.relpath(os.path.join(root, f), self.wd))
        return out

    def run(self, args, keep_files=False):
        """Returns (proc, {relative path: bytes}) for the files dashing2 created."""
        before = self._snapshot()
        self.calls += 1
        # Decode leniently: names read from undecoded compressed input can be arbitrary bytes.
        p = subprocess.run([self.d2.binary] + [str(a) for a in args], cwd=self.wd, capture_output=True,
                           text=True, errors="replace", timeout=300)
        new = {}
        for rel in sorted(self._snapshot() - before):
            fp = os.path.join(self.wd, rel)
            with open(fp, "rb") as f:
                new[rel] = f.read()
            if not keep_files:
                os.remove(fp)
        return p, new

    def ok(self, args):
        p, new = self.run(args)
        if p.returncode != 0:
            raise D2Error("exit %d: dashing2 %s\n%s" % (p.returncode, " ".join(map(str, args)), p.stderr[-600:]))
        return p, new


def pick_file(new, key):
    c = [k for k in new if key in k]
    if len(c) != 1:
        raise D2Error("expected one file containing %r, dashing2 wrote %s" % (key, sorted(new)))
    return new[c[0]]


def parse_kmerset_bytes(b, long):
    card = struct.unpack("<d", b[:8])[0]
    body = b[8:]
    if long:
        return card, [int.from_bytes(body[i:i + 16], "little") for i in range(0, len(body), 16)]
    return card, list(struct.unpack("<%dQ" % (len(body) // 8), body))


def parse_stream_bytes(b, long):
    if long:
        return [int.from_bytes(b[i:i + 16], "little") for i in range(0, len(b), 16)]
    return list(struct.unpack("<%dQ" % (len(b) // 8), b))


def kset(R, flags, path):
    """dashing2 sketch --set: (cardinality, sorted stored values)."""
    long = "-2" in flags
    p, new = R.ok(["sketch", "--set"] + flags + [path])
    return parse_kmerset_bytes(pick_file(new, "kmerset"), long)


def kcounts(R, flags, path):
    """dashing2 sketch -J: {stored value: count}."""
    long = "-2" in flags
    p, new = R.ok(["sketch", "-J"] + flags + [path])
    _, vals = parse_kmerset_bytes(pick_file(new, "kmerset"), long)
    cb = pick_file(new, "kmercounts.f64")
    cnt = struct.unpack("<%dd" % (len(cb) // 8), cb)
    return dict(zip(vals, cnt)), len(vals) == len(cnt)


def gstream(R, flags, path):
    """dashing2 sketch -G -o: the emitted stream of one input."""
    long = "-2" in flags
    p, new = R.ok(["sketch", "-G", "-o", "gout"] + flags + [path])
    return parse_stream_bytes(pick_file(new, "mmerseq"), long)


def cmp_matrix(R, args, rows=None):
    p, _ = R.run(["cmp"] + args)
    if p.returncode != 0:
        raise D2Error("exit %d: dashing2 cmp %s\n%s" % (p.returncode, " ".join(map(str, args)), p.stderr[-600:]))
    M = parse_matrix(p.stdout)
    M.stderr = p.stderr
    if rows is not None and M.rows != list(rows):
        raise D2Error("printed rows %s, expected %s" % (M.rows, list(rows)))
    return M


def sketch_names(R, args, path):
    """dashing2 sketch -o: [(name, cardinality)] from <out>.names.txt."""
    p, new = R.ok(["sketch"] + args + ["-o", "sk_out", path])
    txt = new.get("sk_out.names.txt", b"").decode(errors="replace")
    out = []
    for ln in txt.splitlines():
        if ln.startswith("#") or not ln.strip():
            continue
        t = ln.split("\t")
        out.append((t[0], float(t[1]) if len(t) > 1 else None))
    return out, new


def drop_option(args, opt):
    """args without opt and the value that follows it."""
    out, skip = [], False
    for x in args:
        if skip:
            skip = False
        elif x == opt:
            skip = True
        else:
            out.append(x)
    return out


def crashed(p):
    """A segmentation fault or bus error. Aborts after an uncaught exception (exit 134, signal 6) count as
    rejections: dashing2 reports most invalid options that way (a known, separately listed issue)."""
    return p.returncode in (-11, -10, 138, 139)


# ---------------------------------------------------------------------------
# Inputs


def dna_records(rng, k, maxlen, nrec=None, breaks=True, homopolymers=True, rna=False):
    """Random DNA records with repeats, breaks, lowercase, homopolymers near k, short and empty records."""
    n = rng.randint(80, maxlen)
    base = random_dna(rng, n)
    base = add_repeats(rng, base, rng.randint(0, 4))
    if homopolymers and rng.random() < 0.6:
        base = add_homopolymers(rng, base, rng.randint(1, 3), max(1, k))
    p = dict(nrec=nrec or rng.randint(1, 5), short_recs=rng.randint(0, 2), rc_prob=rng.choice([0, 0.3]),
             n_runs=rng.randint(0, 4) if breaks else 0, lower_runs=rng.randint(0, 2))
    recs = decorate(rng, base, max(k, 2), p)
    if rng.random() < 0.3:
        recs.insert(rng.randrange(len(recs) + 1), "")
    if rna and rng.random() < 0.5:
        recs = [r.replace("T", "U") if rng.random() < 0.5 else r for r in recs]
    return base, recs


AA = "ACDEFGHIKLMNPQRSTVWY"


def protein_records(rng, maxlen):
    recs = []
    for _ in range(rng.randint(1, 4)):
        n = rng.randint(1, maxlen // 3)
        r = list("".join(rng.choices(AA, k=n)))
        for _ in range(rng.randint(0, 3)):
            if r:
                r[rng.randrange(len(r))] = rng.choice("OUouBZXJ*")
        if rng.random() < 0.3:
            a = rng.randrange(len(r))
            for j in range(a, min(len(r), a + rng.randint(1, 40))):
                r[j] = r[j].lower()
        recs.append("".join(r))
    if rng.random() < 0.2:
        recs.append("")
    return recs


def flags_for(k, w=None, canon=True, long=False, entmin=False, spacing=None, alph=DNA, seed=None, extra=()):
    fl = ["-k", str(k)]
    if w:
        fl += ["-w", str(w)]
    if not canon:
        fl.append("--no-canon")
    if long:
        fl.append("-2")
    if entmin:
        fl.append("--entmin")
    if spacing is not None:
        fl += ["--spacing", spacing if isinstance(spacing, str) else ",".join(map(str, spacing))]
    fl += alph.flags
    if seed is not None:
        fl += ["--seed", str(seed)]
    return fl + list(extra)


def spacing_cache_name(gaps):
    """The spacing part of dashing2's output file names: runs of equal offsets as 'offset x count'."""
    out, i = [], 0
    while i < len(gaps):
        j = i
        while j < len(gaps) and gaps[j] == gaps[i]:
            j += 1
        out.append("%dx%d" % (gaps[i] + 1, j - i))
        i = j
    return ",".join(out)


def spacing_string(rng, gaps):
    """Writes a gap list with random use of the run-length 'gap x count' syntax."""
    out, i = [], 0
    while i < len(gaps):
        j = i
        while j < len(gaps) and gaps[j] == gaps[i]:
            j += 1
        n = j - i
        if n > 1 and rng.random() < 0.6:
            c = rng.randint(1, n)
            out.append("%dx%d" % (gaps[i], c) if c > 1 or rng.random() < 0.5 else str(gaps[i]))
            i += c
        else:
            out.append(str(gaps[i]))
            i += 1
    return ",".join(out)


# ---------------------------------------------------------------------------
# Oracle helpers with emulations of registered issues


def all_ones(long):
    return M128 if long else M64


def exact_picks(recs, k, w, alph, canon, spacing, long, entmin):
    return window_picks(recs, k, w, alph, canon, spacing, long, entmin and not spacing)


def emulated_picks(recs, k, w, alph, canon, spacing, long, entmin):
    """What dashing2 emits given the registered issues uncanon-window-polyT and spaced-window-invalid."""
    offs = offsets_from_spacing(k, spacing)
    span = offs[-1] + 1
    m = (max(w or 0, span) - span + 1) if w else 1
    ones = all_ones(long)
    lx = lex128 if long else lex64
    if spacing:
        out = []
        for r in recs:
            codes = [alph.code.get(ch, -1) for ch in r]
            keyed = []
            for i in range(len(r) - span + 1):
                cs = [codes[i + o] for o in offs]
                v = ones if min(cs) < 0 else alph.encode(cs)
                keyed.append((lx(v), v))
            picks = []
            for j in range(m - 1, len(keyed)):
                b = min(keyed[j - m + 1:j + 1])
                if b[1] != ones:
                    picks.append(b[1])
            out.append(picks)
        return out
    picks = window_picks(recs, k, w, alph, canon, None, long, entmin)
    if m > 1 and not (canon and alph is DNA):
        # Full-window picks equal to all ones are dropped; the tail of a short record is not.
        res = []
        for r, p in zip(recs, picks):
            nv = len(kmer_items(r, k, alph, canon, None))
            if nv >= m:
                p = [v for v in p if v != ones]
            res.append(p)
        return res
    return picks


def flat(picks):
    return [v for p in picks for v in p]


def compare_known(res, name, got, exp, emul, kid, fmt=lambda x: "%d items" % len(x)):
    """Checks got == exp; a mismatch that equals the emulation of a registered issue is an expected failure."""
    if got == exp:
        res.check(True, name)
        return True
    detail = "dashing2 %s, oracle %s" % (fmt(got), fmt(exp))
    if emul is not None and got == emul and emul != exp:
        res.check(False, name, detail + " (matches the emulation of %s)" % kid, known=kid)
    else:
        res.check(False, name, detail)
    return False


def set_detail(got, exp):
    g, e = set(got), set(exp)
    return "dashing2 %d, oracle %d, %d only in dashing2, %d only in oracle" % (len(g), len(e), len(g - e), len(e - g))


def check_set(res, name, got_card, got_vals, exp_vals, emul_vals=None, kid=None):
    got_vals = sorted(got_vals)
    exp = sorted(set(exp_vals))
    emul = sorted(set(emul_vals)) if emul_vals is not None else None
    ok = compare_known(res, name, got_vals, exp, emul, kid, fmt=lambda x: "%d values" % len(x))
    if not ok and (emul is None or got_vals != emul):
        res.fails[-1] = "%s: %s" % (name, set_detail(got_vals, exp))
    res.check(got_card == len(got_vals), name + "-card", "header %s, %d values" % (got_card, len(got_vals)))


# ---------------------------------------------------------------------------
# Families


def fam_minimizer(rng, wd, R, res, preset):
    long = rng.random() < 0.35
    cap = 64 if long else 32
    rolling = rng.random() < 0.3
    if rolling:
        k = rng.randint(cap + 1, cap + 60)
    else:
        k = rng.choice([rng.randint(1, cap), cap, cap, rng.randint(max(1, cap - 4), cap)])
    canon = rng.random() < 0.6
    w = rng.choice([None, k - rng.randint(0, 3), k + 1, k + rng.randint(2, 12), k + rng.randint(13, 60)])
    if w is not None and w < 1:
        w = None
    seed = rng.choice([None, None, rng.randint(1, 10 ** 9)])
    base, recs = dna_records(rng, k, preset["maxlen"], rna=True)
    if k == cap and rng.random() < 0.5:
        # Poly-T is the all-ones encoding at full width.
        recs.append(random_dna(rng, 20) + "T" * (k + rng.randint(0, 40)) + random_dna(rng, 20))
    path = "in.fa" if rng.random() < 0.7 else "in.fq.gz"
    write_records(os.path.join(wd, path), recs, fmt="fq" if ".fq" in path else "fa", width=rng.choice([None, 60]))
    res.desc = "minimizer k=%d w=%s %s%s seed=%s %s recs=%d lens=%s" % (
        k, w, "canon" if canon else "--no-canon", " -2" if long else "", seed, path, len(recs),
        [len(r) for r in recs])
    res.tags = {"family": "minimizer", "path": "rolling" if rolling else "direct", "-2": long, "canon": canon,
                "w": "none" if not w or w <= k else "w>k"}
    mk = Mask(seed or 0)
    fl = flags_for(k, w, canon, long, seed=seed)
    if not rolling:
        exp = exact_picks(recs, k, w, DNA, canon, None, long, False)
        emu = emulated_picks(recs, k, w, DNA, canon, None, long, False)
        card, vals = kset(R, fl, path)
        check_set(res, "set", card, vals, [mk.mask(v, long) for v in flat(exp)],
                  [mk.mask(v, long) for v in flat(emu)], "uncanon-window-polyT")
        st = gstream(R, fl, path)
        compare_known(res, "stream", st, [mk.mask(v, long) for v in flat(exp)],
                      [mk.mask(v, long) for v in flat(emu)], "uncanon-window-polyT")
        # -J counts: how often each item was emitted.
        cnt, aligned = kcounts(R, fl, path)
        from collections import Counter
        ec = Counter(mk.mask(v, long) for v in flat(exp))
        emc = Counter(mk.mask(v, long) for v in flat(emu))
        got = {x: int(c) for x, c in cnt.items()}
        res.check(aligned, "countdict-files", "k-mer and count files differ in length")
        compare_known(res, "countdict", got, dict(ec), dict(emc), "uncanon-window-polyT",
                      fmt=lambda d: "%d k-mers, total %d" % (len(d), sum(d.values())))
        # Density sanity for long random records: the oracle must not be degenerate.
        m = (w - k + 1) if w and w > k else 1
        nk = sum(len(kmer_items(r, k, DNA, canon)) for r in recs)
        if m > 1 and nk > 20 * m and k >= 12:
            dens = len(set(flat(exp))) / nk
            res.check(0.3 * 2 / (m + 1) <= dens <= 3.0 * 2 / (m + 1) + 0.02, "density",
                      "distinct minimizers per k-mer %.4f, expected about %.4f" % (dens, 2 / (m + 1)))
    else:
        # Rolling hash: unwindowed set size is the number of distinct k-mers; windowed output is the sliding
        # minimum of the unwindowed stream under the same score.
        card, vals = kset(R, flags_for(k, None, canon, long, seed=seed), path)
        nd = len(distinct_kmer_strings(recs, k, DNA, canon))
        res.check(card == nd and len(vals) == nd, "rolling-card", "dashing2 %s (%d values), oracle %d" % (
            card, len(vals), nd))
        full = gstream(R, flags_for(k, None, canon, long, seed=seed), path)
        nv = [n_valid_kmers(r, k) for r in recs]
        res.check(len(full) == sum(nv), "rolling-stream-length", "dashing2 %d, valid k-mers %d" % (len(full), sum(nv)))
        if w and len(full) == sum(nv):
            raw = [mk.unmask(x, long) for x in full]
            sc = (lambda v, cs: lex128(v) & M64) if long else (lambda v, cs: lex64(v))
            m = max(w - k + 1, 1)
            exp, pos = [], 0
            for n in nv:
                exp += select_window([(v, None) for v in raw[pos:pos + n]], m, sc)
                pos += n
            got = gstream(R, fl, path)
            compare_known(res, "rolling-window-stream", got, [mk.mask(v, long) for v in exp], None, None)
            card2, vals2 = kset(R, fl, path)
            res.check(set(vals2) == {mk.mask(v, long) for v in exp}, "rolling-window-set",
                      set_detail(vals2, [mk.mask(v, long) for v in exp]))
            res.check(set(vals2) <= set(vals), "rolling-window-subset", "windowed set not within unwindowed set")


def fam_entmin(rng, wd, R, res, preset):
    alph = DNA if rng.random() < 0.7 else rng.choice(PROTEIN_ALPHABETS)
    long = rng.random() < 0.35
    cap = alph.cap(long)
    k = rng.randint(max(2, cap // 2), cap)
    w = k + rng.randint(1, 30)
    canon = rng.random() < 0.6
    if alph is DNA:
        base, recs = dna_records(rng, k, preset["maxlen"])
        # Low-entropy stretches make entropy weighting matter.
        for _ in range(rng.randint(1, 4)):
            unit = random_dna(rng, rng.randint(1, 3))
            p = rng.randrange(len(base) + 1)
            base = base[:p] + unit * rng.randint(k // 2, 2 * k) + base[p:]
        recs = recs + [base]
    else:
        recs = protein_records(rng, preset["maxlen"])
        recs.append("".join(rng.choices("AAAAG", k=rng.randint(k, 3 * k))) + "".join(rng.choices(AA, k=200)))
    path = "in.fa"
    write_records(os.path.join(wd, path), recs)
    res.desc = "entmin %s k=%d w=%d %s%s recs=%d" % (alph.name, k, w, "canon" if canon else "--no-canon",
                                                    " -2" if long else "", len(recs))
    res.tags = {"family": "entmin", "alph": alph.name, "-2": long, "canon": canon}
    fl = flags_for(k, w, canon, long, entmin=True, alph=alph)
    exp = exact_picks(recs, k, w, alph, canon, None, long, True)
    emu = emulated_picks(recs, k, w, alph, canon, None, long, True)
    mk = Mask(0)
    st = gstream(R, fl, path)
    expm = [mk.mask(v, long) for v in flat(exp)]
    if st != expm and len(st) == len(expm):
        # The score divides by a floating-point entropy summed in hash-map order; allow a pick whose score is
        # within rounding of the best one, but nothing else.
        exact = window_items_scores(recs, k, w, alph, canon, long)
        bad = 0
        for got, want, (cands, _) in zip(st, expm, exact):
            if got == want:
                continue
            gv = mk.unmask(got, long)
            sc = dict(cands)
            best = min(s for s, _ in cands)
            if gv not in sc or sc[gv] > best * (1 + 1e-12) + 2:
                bad += 1
        res.check(bad == 0, "entmin-stream", "%d windows pick an item that is not entropy-minimal" % bad)
    else:
        compare_known(res, "entmin-stream", st, expm, [mk.mask(v, long) for v in flat(emu)], "uncanon-window-polyT")
    # Property: picks are k-mers of the input.
    allk = {mk.mask(v, long) for r in recs for v, _ in kmer_items(r, k, alph, canon)}
    res.check(set(st) <= allk, "entmin-valid-kmers", "%d picks are not k-mers of the input" % len(set(st) - allk))
    # Property: entropy weighting prefers higher-entropy k-mers than plain minimizers.
    lexp = flat(exact_picks(recs, k, w, alph, canon, None, long, False))
    ent_of = {}
    for r in recs:
        for v, cs in kmer_items(r, k, alph, canon):
            ent_of[v] = entropy(cs)
    if lexp and st:
        e_ent = sum(ent_of[mk.unmask(x, long)] for x in st if mk.unmask(x, long) in ent_of) / len(st)
        e_lex = sum(ent_of[v] for v in lexp) / len(lexp)
        res.check(e_ent >= e_lex - 1e-9, "entmin-prefers-entropy",
                  "mean entropy of --entmin picks %.4f below plain minimizers %.4f" % (e_ent, e_lex))
    # Property: reverse-complement symmetry (canonical DNA).
    if alph is DNA and canon:
        write_records(os.path.join(wd, "rc.fa"), [revcomp(r) for r in recs])
        c1, v1 = kset(R, fl, path)
        c2, v2 = kset(R, fl, "rc.fa")
        res.check(v1 == v2, "entmin-revcomp", set_detail(v2, v1))
    # Property: without -w, --entmin changes nothing; above the exact limit it falls back to plain minimizers.
    if rng.random() < 0.5:
        a = kset(R, flags_for(k, None, canon, long, entmin=True, alph=alph), path)
        b = kset(R, flags_for(k, None, canon, long, alph=alph), path)
        res.check(a == b, "entmin-unwindowed", "--entmin without -w changed the k-mer set")
    else:
        kk = cap + rng.randint(1, 10)
        a = gstream(R, flags_for(kk, kk + 10, canon, long, entmin=True, alph=alph), path)
        b = gstream(R, flags_for(kk, kk + 10, canon, long, alph=alph), path)
        res.check(a == b, "entmin-fallback", "--entmin above the exact limit differs from plain minimizers "
                  "(%d vs %d items)" % (len(a), len(b)))


def window_items_scores(recs, k, w, alph, canon, long):
    """Per full window (in emission order): the candidate (value, score) pairs, for tolerance checks."""
    out = []
    m = w - k + 1
    for r in recs:
        items = [(v, ent_score(v, entropy(cs), long)) for v, cs in kmer_items(r, k, alph, canon)]
        for j in range(m - 1, len(items)):
            out.append((items[j - m + 1:j + 1], None))
        if 0 < len(items) < m:
            out.append((items, None))
    return out


def fam_spacing(rng, wd, R, res, preset):
    long = rng.random() < 0.3
    alph = DNA if rng.random() < 0.8 else rng.choice(PROTEIN_ALPHABETS)
    cap = alph.cap(long)
    k = rng.choice([rng.randint(2, min(cap, 20)), rng.randint(2, cap), cap])
    gaps = [rng.choice([0, 0, 1, 2, 3, 5]) for _ in range(k - 1)]
    if rng.random() < 0.3:
        g = rng.randint(0, 3)
        gaps = [g] * (k - 1)
    sstr = spacing_string(rng, gaps)
    span = sum(gaps) + k
    w = rng.choice([None, None, span + rng.randint(1, 20)])
    canon = rng.random() < 0.5
    if alph is DNA:
        base, recs = dna_records(rng, span, preset["maxlen"])
    else:
        recs = protein_records(rng, preset["maxlen"])
    path = "in.fa"
    write_records(os.path.join(wd, path), recs)
    res.desc = "spacing %s k=%d --spacing %s (span %d) w=%s %s%s" % (alph.name, k, sstr, span, w,
                                                                   "canon" if canon else "--no-canon",
                                                                   " -2" if long else "")
    res.tags = {"family": "spacing", "alph": alph.name, "-2": long, "w": bool(w)}
    fl = flags_for(k, w, canon, long, spacing=sstr, alph=alph)
    mk = Mask(0)
    exp = flat(exact_picks(recs, k, w, alph, canon, gaps, long, False))
    emu = flat(emulated_picks(recs, k, w, alph, canon, gaps, long, False))
    # Exact modes write <input><spacing>.k<k>... files whose name spells out every gap.
    namelen = len(path) + len(spacing_cache_name(gaps)) + 40
    try:
        card, vals = kset(R, fl, path)
    except D2Error as e:
        long_name = namelen > 255 and "Failed to open" in str(e)
        res.check(False, "spaced-set-run", str(e)[-300:], known="spacing-long-filename" if long_name else None)
        return
    check_set(res, "spaced-set", card, vals, [mk.mask(v, long) for v in exp], [mk.mask(v, long) for v in emu],
              "spaced-window-invalid" if w else "uncanon-window-polyT")
    st = gstream(R, fl, path)
    compare_known(res, "spaced-stream", st, [mk.mask(v, long) for v in exp], [mk.mask(v, long) for v in emu],
                  "spaced-window-invalid" if w else "uncanon-window-polyT")
    # The expanded list and the run-length form are the same seed.
    if sstr != ",".join(map(str, gaps)):
        c2, v2 = kset(R, flags_for(k, w, canon, long, spacing=gaps, alph=alph), path)
        res.check(v2 == vals, "spacing-syntax", "run-length %s and expanded forms differ" % sstr)
    # All-zero spacing is the contiguous seed but is never canonicalized.
    if rng.random() < 0.3 and alph is DNA:
        c3, v3 = kset(R, flags_for(k, None, True, long, spacing=[0] * (k - 1)), path)
        c4, v4 = kset(R, flags_for(k, None, False, long), path)
        res.check(v3 == v4, "spacing-zero-is-uncanonical", set_detail(v3, v4))
    # Rejections: k above the exact limit, wrong number of gaps.
    kk = cap + rng.randint(1, 8)
    p, new = R.run(["sketch", "--set"] + flags_for(kk, None, True, long, spacing="1x%d" % (kk - 1), alph=alph) +
                   [path])
    res.check(p.returncode != 0 and not crashed(p) and not any("kmerset" in f for f in new), "spacing-reject-capacity",
              "k=%d above capacity %d: exit %d, files %s" % (kk, cap, p.returncode, sorted(new)))
    bad = ",".join(["0"] * (k + rng.randint(0, 3)))
    p, new = R.run(["sketch", "--set"] + flags_for(k, None, True, long, spacing=bad, alph=alph) + [path])
    res.check(p.returncode != 0 and not crashed(p) and not any("kmerset" in f for f in new), "spacing-reject-length",
              "%d gaps for k=%d: exit %d, files %s" % (len(bad.split(",")), k, p.returncode, sorted(new)))


def fam_protein(rng, wd, R, res, preset):
    alph = rng.choice(PROTEIN_ALPHABETS)
    flag = rng.choice(PROTEIN_ALIASES[alph.name])
    long = rng.random() < 0.4
    cap = alph.cap(long)
    rolling = rng.random() < 0.25
    k = rng.randint(cap + 1, cap + 25) if rolling else rng.choice([rng.randint(1, 6), rng.randint(1, cap), cap])
    w = None if rng.random() < 0.6 or rolling else k + rng.randint(1, 15)
    canon = rng.random() < 0.5
    recs = protein_records(rng, preset["maxlen"])
    path = rng.choice(["in.fa", "in.fa.gz", "in.fq"])
    write_records(os.path.join(wd, path), recs, fmt="fq" if ".fq" in path else "fa", width=rng.choice([None, 50]))
    res.desc = "protein %s k=%d w=%s %s%s %s" % (flag, k, w, "canon" if canon else "--no-canon",
                                                " -2" if long else "", path)
    res.tags = {"family": "protein", "alph": alph.name, "-2": long, "path": "rolling" if rolling else "direct"}
    fl = ["-k", str(k), flag] + (["-w", str(w)] if w else []) + ([] if canon else ["--no-canon"]) + \
         (["-2"] if long else [])
    mk = Mask(0)
    card, vals = kset(R, fl, path)
    if rolling:
        nd = len(distinct_kmer_strings(recs, k, alph, False))
        res.check(card == nd and len(vals) == nd, "protein-rolling-card", "dashing2 %s (%d values), oracle %d" % (
            card, len(vals), nd))
    else:
        exp = flat(exact_picks(recs, k, w, alph, False, None, long, False))
        check_set(res, "protein-set", card, vals, [mk.mask(v, long) for v in exp])
    # Protein k-mers are never canonicalized: --no-canon must not matter.
    fl2 = [x for x in fl if x != "--no-canon"] + ([] if not canon else ["--no-canon"])
    c2, v2 = kset(R, fl2, path)
    res.check(v2 == vals, "protein-canon-irrelevant", set_detail(v2, vals))
    # Aliases: O reads as K and U as C; writing them out must give the same k-mers.
    swapped = [r.replace("K", "O").replace("C", "U").replace("k", "o").replace("c", "u") for r in recs]
    write_records(os.path.join(wd, "alias.fa"), swapped)
    c3, v3 = kset(R, fl, "alias.fa")
    res.check(v3 == vals, "protein-aliases", set_detail(v3, vals))


def fam_byseq(rng, wd, R, res, preset):
    alph = DNA if rng.random() < 0.8 else rng.choice(PROTEIN_ALPHABETS)
    long = rng.random() < 0.25
    cap = alph.cap(long)
    k = rng.choice([rng.randint(5, min(cap, 21)), rng.randint(5, cap), cap + rng.randint(1, 30)])
    modes = [("oph", []), ("full", ["--full"]), ("bmh", ["-B"]), ("pmh", ["--prob"]),
             ("full-fc2", ["--full", "--fastcmp", "2"]), ("full-bb1", ["--full", "--bbit-sigs", "--fastcmp", "1"]),
             ("full-ab-shorts", ["--full", "--fastcmp-shorts"])]
    mname, mflags = rng.choice(modes)
    extra = []
    feature = rng.choice(["plain", "plain", "window", "entmin", "spacing", "filter", "downsample", "canon"])
    w = None
    if feature in ("window", "entmin"):
        w = k + rng.randint(1, 20)
        extra += ["-w", str(w)]
        if feature == "entmin":
            extra.append("--entmin")
    if feature == "spacing" and k <= cap:
        gaps = [rng.choice([0, 1, 2]) for _ in range(k - 1)]
        extra += ["--spacing", spacing_string(rng, gaps)]
    if feature == "canon" or rng.random() < 0.2:
        extra.append("--no-canon")
    if long:
        extra.append("-2")
    if alph is not DNA:
        extra += alph.flags
    S = rng.choice([64, 128, 256, 512])
    seed = rng.randint(1, 10 ** 6)
    if alph is DNA:
        base, recs = dna_records(rng, k, preset["maxlen"], nrec=rng.randint(2, 6))
        recs = [r for r in recs]
        # Records shorter than k, empty records and repeated names.
        if rng.random() < 0.5:
            recs.insert(rng.randrange(len(recs) + 1), base[:max(1, k - rng.randint(1, 3))])
    else:
        recs = protein_records(rng, preset["maxlen"])
        if len(recs) < 2:
            recs.append(recs[0][::-1] if recs[0] else "ACD")
    names = ["s%d" % i for i in range(len(recs))]
    if rng.random() < 0.3 and len(names) > 1:
        names[rng.randrange(1, len(names))] = names[0]
    if feature == "filter":
        frecs = [random_dna(rng, 300)] if alph is not DNA else [recs[0][len(recs[0]) // 3:]]
        write_records(os.path.join(wd, "flt.fa"), frecs)
        extra += ["--filterset", "flt.fa"]
    if feature == "downsample":
        extra += ["--downsample", str(rng.choice([0.25, 0.5, 0.75]))]
    ext = rng.choice(["fa", "fa", "fa", "fq", "fa.gz", "fq.gz"] + (["fa.bz2", "fa.xz"] if rng.random() < 0.3 else []))
    path = "recs." + ext
    fmt = "fq" if ext.startswith("fq") else "fa"
    write_records(os.path.join(wd, path), recs, names, fmt=fmt)
    files = []
    for i, r in enumerate(recs):
        fn = "rec%d.%s" % (i, "fq" if fmt == "fq" else "fa")
        write_records(os.path.join(wd, fn), [r], [names[i]], fmt=fmt)
        files.append(fn)
    common = ["-k", str(k), "-S", str(S), "--seed", str(seed)] + mflags + extra
    res.desc = "byseq %s mode=%s feature=%s k=%d S=%d %s recs=%d lens=%s names=%s" % (
        alph.name, mname, feature, k, S, " ".join(extra), len(recs), [len(r) for r in recs], names)
    res.tags = {"family": "byseq", "mode": mname, "feature": feature, "ext": ext.split(".")[-1]}
    compressed_unsupported = ext.endswith(("bz2", "xz", "zst"))
    # Exact modes are rejected cleanly.
    p, _ = R.run(["cmp", "--parse-by-seq", rng.choice(["--set", "-J"]), "-k", str(k), path])
    res.check(p.returncode != 0 and not crashed(p), "byseq-exact-rejected", "exit %d" % p.returncode)
    measures = rng.sample(list(MEAS_FLAGS), 3)
    seq_ok = True
    for meas in measures:
        mf = MEAS_FLAGS[meas]
        try:
            A = cmp_matrix(R, ["--parse-by-seq", "--square"] + common + mf + [path])
        except D2Error as e:
            res.check(False, "byseq-run-%s" % meas, str(e))
            seq_ok = False
            continue
        B = cmp_matrix(R, ["--square"] + common + mf + files, rows=files)
        if A.rows != names:
            res.check(False, "byseq-rows", "printed rows %s, expected record names %s" % (A.rows[:8], names),
                      known="byseq-compressed" if compressed_unsupported else None)
            seq_ok = False
            break
        res.check(True, "byseq-rows")
        n = len(recs)
        diffs = []
        for i in range(n):
            for j in range(n):
                a, b = A.get(i, j), B.get(i, j)
                if not close(a, b, rel=2e-6, abs_=1e-6):
                    diffs.append((i, j, a, b))
        if not diffs:
            res.check(True, "byseq-%s" % meas)
            continue
        kid = None
        if "--downsample" in extra:
            # Emulation: the by-seq result equals the same run without --downsample.
            nod = drop_option(common, "--downsample")
            C = cmp_matrix(R, ["--parse-by-seq", "--square"] + nod + mf + [path])
            if all(close(A.get(i, j), C.get(i, j)) for i in range(n) for j in range(n)):
                kid = "byseq-downsample"
        if kid is None and mname in ("oph", "full", "full-fc2", "full-bb1", "full-ab-shorts") and meas in ("intersection", "union", "containment",
                                                                           "symcontain"):
            if byseq_card_emulation(R, common, path, files, A, B, meas, k, S):
                kid = "byseq-exact-card"
        i, j, a, b = diffs[0]
        res.check(False, "byseq-%s" % meas, "%d of %d cells differ, e.g. [%s,%s] by-seq %s, files %s" % (
            len(diffs), n * n, names[i], names[j], a, b), known=kid)
    # By-seq -G -o stream file equals each record's own stream (the per-record streams are checked against the
    # oracle in the seqmode family).
    if seq_ok and "--filterset" not in extra and not compressed_unsupported and rng.random() < 0.6:
        gl = ["-k", str(k)] + drop_option(extra, "--downsample")
        if w is None and rng.random() < 0.5:
            gl += ["-w", str(k + rng.randint(1, 15))]
        p, new = R.ok(["sketch", "-G", "--parse-by-seq", "-o", "bs_out"] + gl + [path])
        stream_file = os.path.join(wd, "bs_out")
        with open(stream_file, "wb") as f:
            f.write(new["bs_out"])
        info = read_byseq_stream_file(stream_file, long)
        per = []
        bmask = Mask(0)
        for fn in files:
            per.append([bmask.unmask(x, long) for x in gstream(R, gl, fn)])
        got = info["streams"]
        res.check(info["n"] == len(recs) and info["trailing"] == 0, "byseq-G-layout",
                  "n=%d records=%d trailing bytes %d" % (info["n"], len(recs), info["trailing"]))
        kid = None
        if got != per:
            rolling = k > cap
            noncanon = "--no-canon" in gl or alph is not DNA or "--spacing" in gl
            if rolling and noncanon and "-w" in gl:
                # A record shorter than k (but not empty) gets the minimum of the last window of the most
                # recent record that reset the window, written through invmaskfn a second time.
                emu = [list(s) for s in per]
                last = None
                for i in range(len(recs)):
                    if len(recs[i]) >= k:
                        # Records of length >= k reset the window; it then holds their last window.
                        last = per[i][-1] if per[i] else None
                    elif recs[i] and last is not None and not per[i]:
                        emu[i] = [bmask.unmask(last, long)]
                if emu == got:
                    kid = "byseq-seq-stale-window"
        bad = [i for i in range(min(len(got), len(per))) if got[i] != per[i]]
        res.check(got == per, "byseq-G-streams", "records %s differ, e.g. %s: by-seq %d items, file %d items" % (
            bad[:5], bad[:1], len(got[bad[0]]) if bad else -1, len(per[bad[0]]) if bad else -1), known=kid)


def byseq_card_emulation(R, common, path, files, A, B, meas, k, S):
    """True if the by-seq matrix equals the by-file similarity combined with exact cardinalities."""
    try:
        J = cmp_matrix(R, ["--square"] + common + files, rows=files)
        cards, _ = sketch_names(R, ["--parse-by-seq"] + common, path)
    except D2Error:
        return False
    c = [x[1] for x in cards]
    n = len(files)
    for i in range(n):
        for j in range(n):
            j_ = J.get(i, j)
            a, b = c[i], c[j]
            if a == 0 or b == 0:
                continue
            I = j_ * (a + b) / (1 + j_)
            e = {"intersection": I, "union": a + b - I, "containment": I / a, "symcontain": I / min(a, b)}[meas]
            if not close(A.get(i, j), e, rel=1e-4, abs_=1e-4):
                return False
    return True


def fam_seqmode(rng, wd, R, res, preset):
    long = rng.random() < 0.3
    cap = 64 if long else 32
    k = rng.randint(5, cap)
    w = rng.choice([None, k + rng.randint(1, 20)])
    canon = rng.random() < 0.7
    hp = rng.random() < 0.4
    base, recs1 = dna_records(rng, k, preset["maxlen"])
    recs2 = [mutate(rng, r, 0.01, 0.003) for r in recs1]
    if rng.random() < 0.5:
        recs2 = recs2[1:] + recs2[:1]
    write_records(os.path.join(wd, "a.fa"), recs1)
    write_records(os.path.join(wd, "b.fa"), recs2)
    res.desc = "seqmode k=%d w=%s %s%s%s recs=%d" % (k, w, "canon" if canon else "--no-canon", " -2" if long else "",
                                                     " --hp-compress" if hp else "", len(recs1))
    res.tags = {"family": "seqmode", "-2": long, "hp": hp, "w": bool(w)}
    mk = Mask(0)
    fl = flags_for(k, w, canon, long) + (["--hp-compress"] if hp else [])
    streams, emus = [], []
    for fn, recs in (("a.fa", recs1), ("b.fa", recs2)):
        exp = [mk.mask(v, long) for v in flat(exact_picks(recs, k, w, DNA, canon, None, long, False))]
        emu = [mk.mask(v, long) for v in flat(emulated_picks(recs, k, w, DNA, canon, None, long, False))]
        if hp:
            exp, emu = hp_collapse(exp), hp_collapse(emu)
        got = gstream(R, fl, fn)
        compare_known(res, "G-stream-%s" % fn, got, exp, emu, "uncanon-window-polyT")
        streams.append(exp)
        emus.append(got)
    # cmp -G: exact_kmer_dist is forced on for -G (src/cmp_main.h), so every measure prints the longer stream
    # length minus the item-level edit distance, and --compute-edit-distance prints the edit distance itself.
    # (The Hamming comparison in src/wcompare.cpp is not reachable from the command line.)
    sa, sb = emus
    if len(sa) * len(sb) <= 4e6:
        ed = levenshtein(sa, sb)
        sim = max(len(sa), len(sb)) - ed
        meas = rng.choice(["", "--intersection", "--mash-distance"])
        M = cmp_matrix(R, ["-G", "-o", "cmp_out", "--square"] + fl + ([meas] if meas else []) + ["a.fa", "b.fa"],
                       rows=["a.fa", "b.fa"])
        res.check(close(M.get(0, 1), sim) and close(M.get(1, 0), sim), "G-similarity",
                  "dashing2 %s/%s, max length minus edit distance %d" % (M.get(0, 1), M.get(1, 0), sim))
        res.check(close(M.get(0, 0), len(sa)) and close(M.get(1, 1), len(sb)), "G-similarity-diag",
                  "diagonal %s/%s, lengths %d/%d" % (M.get(0, 0), M.get(1, 1), len(sa), len(sb)))
        M = cmp_matrix(R, ["-G", "-o", "cmp_out", "--compute-edit-distance"] + fl + ["a.fa", "b.fa"],
                       rows=["a.fa", "b.fa"])
        res.check(close(M.get(0, 1), ed), "G-edit-distance", "dashing2 %s, Levenshtein %d" % (M.get(0, 1), ed))
    # The cardinality of a -G sketch is its stream length.
    cards, _ = sketch_names(R, ["-G"] + fl, "a.fa")
    res.check(len(cards) == 1 and cards[0][1] == len(sa), "G-cardinality", "names.txt %s, stream %d" % (cards, len(sa)))
    # By-seq stream file: header and per-record streams against the oracle.
    names = ["q%d" % i for i in range(len(recs1))]
    write_records(os.path.join(wd, "bs.fa"), recs1, names)
    p, new = R.ok(["sketch", "-G", "--parse-by-seq", "-o", "bs_out"] + fl + ["bs.fa"])
    with open(os.path.join(wd, "bs_out"), "wb") as f:
        f.write(new["bs_out"])
    info = read_byseq_stream_file(os.path.join(wd, "bs_out"), long)
    picks = exact_picks(recs1, k, w, DNA, canon, None, long, False)
    epicks = emulated_picks(recs1, k, w, DNA, canon, None, long, False)
    # By-seq stream files hold the k-mers themselves (src/fastxsketchbyseq.cpp applies invmaskfn before
    # writing), whereas file-mode .mmerseq files hold masked values.
    exp = [list(pk) for pk in picks]
    emu = [list(pk) for pk in epicks]
    if hp:
        exp = [hp_collapse(x) for x in exp]
        emu = [hp_collapse(x) for x in emu]
    res.check(info["k"] == k and info["w"] in ((w or k), 0xffffffff) and info["n"] == len(recs1) and info["trailing"] == 0,
              "byseq-G-header", "n=%d k=%d w=%d trailing=%d" % (info["n"], info["k"], info["w"], info["trailing"]))
    res.check((info["dtype"] >> 8) == int(canon) and (info["dtype"] & 0xff) == 0, "byseq-G-dtype",
              "dtype field 0x%x for canon=%s" % (info["dtype"], canon))
    compare_known(res, "byseq-G-streams", info["streams"], exp, emu, "uncanon-window-polyT",
                  fmt=lambda x: "%d streams, %d items" % (len(x), sum(len(s) for s in x)))
    lens_ok = all(int(l) == len(s) * (2 if long else 1) for l, s in zip(info["lens"], info["streams"]))
    res.check(lens_ok, "byseq-G-lengths", "lengths %s are not stream sizes in 64-bit words" % (info["lens"][:6],))
    nm = new.get("bs_out.names.txt", b"").decode().splitlines()[1:]
    res.check([x.split("\t")[0] for x in nm] == names, "byseq-G-names", "names %s" % nm[:6])


def fam_filter(rng, wd, R, res, preset):
    long = rng.random() < 0.35
    cap = 64 if long else 32
    rolling = rng.random() < 0.25
    k = rng.randint(cap + 1, cap + 40) if rolling else rng.randint(8, cap)
    w = None if rolling or rng.random() < 0.6 else k + rng.randint(1, 20)
    canon = rng.random() < 0.7
    base, recs = dna_records(rng, k, preset["maxlen"], homopolymers=False)
    # The filter: fragments of the base (some long, some short) and unrelated sequence.
    frecs = []
    for _ in range(rng.randint(1, 8)):
        a = rng.randrange(len(base))
        frecs.append(base[a:a + rng.choice([k + rng.randint(0, 10), rng.randint(k, 600)])])
    frecs.append(random_dna(rng, rng.randint(0, 200)))
    write_records(os.path.join(wd, "in.fa"), recs)
    write_records(os.path.join(wd, "flt.fa"), frecs)
    form = "seq" if rolling or rng.random() < 0.5 else rng.choice(["bin", "bin.gz"])
    res.desc = "filter k=%d w=%s %s%s form=%s nfilter=%d" % (k, w, "canon" if canon else "--no-canon",
                                                           " -2" if long else "", form, len(frecs))
    res.tags = {"family": "filter", "-2": long, "path": "rolling" if rolling else "direct", "form": form,
                "w": bool(w)}
    mk = Mask(0)
    fl = flags_for(k, w, canon, long)
    if rolling:
        card, vals = kset(R, fl + ["--filterset", "flt.fa"], "in.fa")
        exp = distinct_kmer_strings(recs, k, DNA, canon) - distinct_kmer_strings(frecs, k, DNA, canon)
        res.check(card == len(exp), "filter-rolling-card", "dashing2 %s, oracle %d" % (card, len(exp)))
        return
    allf = {mk.mask(v, long) for r in frecs for v, _ in kmer_items(r, k, DNA, canon)}
    if form == "seq":
        farg = "flt.fa"
    else:
        vals = sorted(allf)
        data = b"".join(v.to_bytes(16 if long else 8, "little") for v in vals)
        import gzip as _gz
        fn = "flt.bin" + (".gz" if form == "bin.gz" else "")
        with (_gz.open if form == "bin.gz" else open)(os.path.join(wd, fn), "wb") as f:
            f.write(data)
        farg = fn + ":b"
    picks = flat(exact_picks(recs, k, w, DNA, canon, None, long, False))
    epicks = flat(emulated_picks(recs, k, w, DNA, canon, None, long, False))
    exp = {mk.mask(v, long) for v in picks} - allf
    card, vals = kset(R, fl + ["--filterset", farg], "in.fa")
    emul = None
    kid = "uncanon-window-polyT"
    if w and form == "seq":
        minf = {mk.mask(v, long) for v in flat(exact_picks(frecs, k, w, DNA, canon, None, long, False))}
        emul = sorted({mk.mask(v, long) for v in epicks} - minf)
        kid = "filterset-windowed"
    elif epicks != picks:
        emul = sorted({mk.mask(v, long) for v in epicks} - allf)
    check_set(res, "filter-set", card, vals, exp, emul, kid)
    # The filtered -G stream drops filtered items in place.
    if rng.random() < 0.5:
        st = gstream(R, fl + ["--filterset", farg], "in.fa")
        filt = allf if not (w and form == "seq") else allf
        exps = [x for x in (mk.mask(v, long) for v in picks) if x not in filt]
        emus = None
        if w and form == "seq":
            emus = [x for x in (mk.mask(v, long) for v in epicks) if x not in minf]
        compare_known(res, "filter-stream", st, exps, emus, "filterset-windowed" if emus is not None else None)


def fam_downsample(rng, wd, R, res, preset):
    long = rng.random() < 0.3
    cap = 64 if long else 32
    k = rng.randint(11, cap)
    f = rng.choice([0.05, 0.1, 0.25, 0.5, 0.7, 0.9, rng.random()])
    seed = rng.choice([None, rng.randint(1, 10 ** 6)])
    canon = rng.random() < 0.7
    n = rng.randint(2000, max(2500, preset["maxlen"] * 2))
    base = random_dna(rng, n)
    other = mutate(rng, base, rng.choice([0.002, 0.01, 0.03]), 0.001)
    write_records(os.path.join(wd, "a.fa"), [base])
    write_records(os.path.join(wd, "b.fa"), [other])
    res.desc = "downsample f=%.4f k=%d %s%s seed=%s n=%d" % (f, k, "canon" if canon else "--no-canon",
                                                           " -2" if long else "", seed, n)
    res.tags = {"family": "downsample", "-2": long, "f": "<0.3" if f < 0.3 else ">=0.3"}
    mk = Mask(seed or 0)
    fl = flags_for(k, None, canon, long, seed=seed)
    thr = downsample_threshold(f)
    sets = []
    for fn, rec in (("a.fa", base), ("b.fa", other)):
        full = {mk.mask(v, long) for v, _ in kmer_items(rec, k, DNA, canon)}
        exp = {x for x in full if downsample_hash(x, long) < thr}
        card, vals = kset(R, fl + ["--downsample", repr(f)], fn)
        check_set(res, "downsample-set-%s" % fn, card, vals, exp)
        nf = len(full)
        sd = math.sqrt(nf * f * (1 - f))
        res.check(abs(len(vals) - nf * f) <= 5 * sd + 2, "downsample-fraction-%s" % fn,
                  "kept %d of %d (f=%.3f, expected %.1f +- %.1f)" % (len(vals), nf, f, nf * f, sd))
        sets.append((full, set(vals)))
    (fa, da), (fb, db) = sets
    J = len(fa & fb) / len(fa | fb)
    U = len(da | db)
    Jd = len(da & db) / U if U else 1.0
    se = math.sqrt(max(J * (1 - J), 1e-12) / max(U, 1)) + 1.0 / max(U, 1)
    res.check(abs(Jd - J) <= 5 * se, "downsample-jaccard", "J %.4f, downsampled J %.4f (SE %.4f, U=%d)" % (
        J, Jd, se, U))
    # The printed similarity of the downsampled exact sets.
    M = cmp_matrix(R, ["--set"] + fl + ["--downsample", repr(f), "a.fa", "b.fa"], rows=["a.fa", "b.fa"])
    res.check(close(M.get(0, 1), Jd, rel=1e-5), "downsample-cmp", "dashing2 %s, downsampled J %.6f" % (M.get(0, 1), Jd))
    # Edge values: 1 keeps everything, 0 keeps nothing, values outside [0, 1] are rejected.
    c1, v1 = kset(R, fl + ["--downsample", "1"], "a.fa")
    res.check(set(v1) == fa, "downsample-one", set_detail(v1, fa))
    c0, v0 = kset(R, fl + ["--downsample", "0"], "a.fa")
    res.check(c0 == 0 and not v0, "downsample-zero", "kept %d" % len(v0))
    p, _ = R.run(["sketch", "--set"] + fl + ["--downsample", rng.choice(["1.5", "-0.2"]), "a.fa"])
    res.check(p.returncode != 0 and not crashed(p), "downsample-reject", "exit %d" % p.returncode)
    # The -G stream keeps the sampled items in order.
    if rng.random() < 0.5:
        st = gstream(R, fl + ["--downsample", repr(f)], "a.fa")
        exps = [x for x in (mk.mask(v, long) for v, _ in kmer_items(base, k, DNA, canon)) if downsample_hash(x, long) < thr]
        compare_known(res, "downsample-stream", st, exps, None, None)


def fam_hpseed(rng, wd, R, res, preset):
    long = rng.random() < 0.3
    cap = 64 if long else 32
    k = rng.randint(3, cap)
    w = rng.choice([None, k + rng.randint(1, 30)])
    canon = rng.random() < 0.7
    base, recs = dna_records(rng, k, preset["maxlen"])
    # Low-complexity stretches give runs of equal minimizers across windows.
    recs.append(random_dna(rng, 3) * rng.randint(10, 60) + random_dna(rng, 50))
    write_records(os.path.join(wd, "a.fa"), recs)
    b = [mutate(rng, r, 0.02, 0.0) for r in recs]
    write_records(os.path.join(wd, "b.fa"), b)
    seeds = [0] + [rng.choice([1, 42, rng.randint(2, 10 ** 9), 2 ** 64 - 1, 2 ** 63])
                   for _ in range(2)]
    res.desc = "hpseed k=%d w=%s %s%s seeds=%s recs=%d" % (k, w, "canon" if canon else "--no-canon",
                                                         " -2" if long else "", seeds, len(recs))
    res.tags = {"family": "hpseed", "-2": long, "w": bool(w)}
    picks = flat(exact_picks(recs, k, w, DNA, canon, None, long, False))
    epicks = flat(emulated_picks(recs, k, w, DNA, canon, None, long, False))
    outs = []
    for s in seeds:
        mk = Mask(s)
        fl = flags_for(k, w, canon, long, seed=s if s else None)
        # --hp-compress: the -G stream with consecutive repeats collapsed, across record boundaries.
        st = gstream(R, fl + ["--hp-compress"], "a.fa")
        compare_known(res, "hp-stream-seed%d" % s, st, hp_collapse([mk.mask(v, long) for v in picks]),
                      hp_collapse([mk.mask(v, long) for v in epicks]), "uncanon-window-polyT")
        # --hp-compress does not change exact sets; the stored values follow the seed's mask.
        c1, v1 = kset(R, fl + ["--hp-compress"], "a.fa")
        check_set(res, "seed-set-%d" % s, c1, v1, [mk.mask(v, long) for v in picks],
                  [mk.mask(v, long) for v in epicks], "uncanon-window-polyT")
        mode = rng.choice(["--set", "-J"])
        meas = rng.choice(list(MEAS_FLAGS))
        M = cmp_matrix(R, [mode] + fl + MEAS_FLAGS[meas] + ["--square", "a.fa", "b.fa"], rows=["a.fa", "b.fa"])
        outs.append((mode, meas, [M.get(i, j) for i in range(2) for j in range(2)]))
        # Sketches with a seed are valid estimates (similarity in [0, 1], 1 on the diagonal).
        S = rng.choice([128, 512])
        M2 = cmp_matrix(R, [rng.choice(["--full", "-B"]), "-S", str(S)] + fl + ["--square", "a.fa", "b.fa"],
                        rows=["a.fa", "b.fa"])
        vals = [M2.get(i, j) for i in range(2) for j in range(2)]
        res.check(all(v is not None and -1e-6 <= v <= 1 + 1e-6 for v in vals) and close(vals[0], 1) and
                  close(vals[3], 1), "seed-sketch-valid", "similarities %s" % vals)
    # Exact-mode measures do not depend on the seed (same mode and measure for every seed).
    mode, meas = outs[0][0], outs[0][1]
    ref = None
    for s in seeds:
        fl = flags_for(k, w, canon, long, seed=s if s else None)
        M = cmp_matrix(R, [mode] + fl + MEAS_FLAGS[meas] + ["--square", "a.fa", "b.fa"], rows=["a.fa", "b.fa"])
        v = [M.get(i, j) for i in range(2) for j in range(2)]
        if ref is None:
            ref = v
        res.check(all(close(x, y) for x, y in zip(v, ref)), "seed-invariant-%s-%s" % (mode, meas),
                  "seed %d gives %s, seed 0 gives %s" % (s, v, ref))
    # --hp-compress with --parse-by-seq collapses within each record.
    names = ["h%d" % i for i in range(len(recs))]
    write_records(os.path.join(wd, "bs.fa"), recs, names)
    mk = Mask(0)
    fl = flags_for(k, w, canon, long)
    p, new = R.ok(["sketch", "-G", "--parse-by-seq", "--hp-compress", "-o", "bs_out"] + fl + ["bs.fa"])
    with open(os.path.join(wd, "bs_out"), "wb") as f:
        f.write(new["bs_out"])
    info = read_byseq_stream_file(os.path.join(wd, "bs_out"), long)
    exp = [hp_collapse(list(pk)) for pk in exact_picks(recs, k, w, DNA, canon, None, long, False)]
    emu = [hp_collapse(list(pk)) for pk in emulated_picks(recs, k, w, DNA, canon, None, long, False)]
    compare_known(res, "byseq-hp-streams", info["streams"], exp, emu, "uncanon-window-polyT",
                  fmt=lambda x: "%d streams, %d items" % (len(x), sum(len(s) for s in x)))


def fam_inputs(rng, wd, R, res, preset, kmc):
    k = rng.randint(5, 31)
    canon = rng.random() < 0.7
    n = rng.randint(2, 4)
    base = random_dna(rng, rng.randint(300, preset["maxlen"]))
    exts = ["fa", "fq", "fa.gz", "fq.gz", "fa.bz2", "fa.xz"] + (["fa.zst"] if have_zstd() else [])
    files, recsets = [], []
    for i in range(n):
        s = mutate(rng, base, rng.choice([0, 0.005, 0.02]), 0.002)
        s = add_repeats(rng, s, rng.randint(0, 3))
        recs = decorate(rng, s, k, dict(nrec=rng.randint(1, 4), short_recs=rng.randint(0, 2), n_runs=rng.randint(0, 3),
                                        lower_runs=rng.randint(0, 2)))
        ext = rng.choice(exts)
        sub = rng.random() < 0.3
        fn = ("sub/" if sub else "") + "in%d.%s" % (i, ext)
        if sub:
            os.makedirs(os.path.join(wd, "sub"), exist_ok=True)
        write_records(os.path.join(wd, fn), recs, fmt="fq" if ext.startswith("fq") else "fa",
                      width=rng.choice([None, 60, 7]), comment=rng.random() < 0.5)
        files.append(fn)
        recsets.append(recs)
    # Entries: single files and sometimes one joint entry (two files on one line, sketched together).
    entries = [[f] for f in files]
    if n >= 3 and rng.random() < 0.6:
        entries = [[files[0], files[1]]] + [[f] for f in files[2:]]
    names = [" ".join(e) for e in entries]
    erecs = [sum((recsets[files.index(f)] for f in e), []) for e in entries]
    layout = rng.choice(["F", "FQ", "pos"])
    res.desc = "inputs k=%d %s files=%s entries=%s layout=%s" % (k, "canon" if canon else "--no-canon", files,
                                                                 names, layout)
    res.tags = {"family": "inputs", "layout": layout, "joint": len(entries) < len(files)}
    # Counts: KMC on an uncompressed FASTA copy of each entry, checked against the Python oracle.
    counts = []
    for i, recs in enumerate(erecs):
        c = kmer_counts(recs, k, canon)
        if kmc is not None:
            plain = "kmc_in%d.fa" % i
            write_records(os.path.join(wd, plain), recs)
            db = kmc.count(plain, k, canon, wd, "db%d" % i)
            kc = kmc.dump(db, wd)
            res.check(kc == c, "oracle-vs-kmc", "entry %d: KMC %d k-mers, Python %d" % (i, len(kc), len(c)))
            c = kc
            os.remove(os.path.join(wd, plain))
        counts.append(c)
    # Lists with blank lines, extra spaces and CRLF line ends.
    def listfile(name, ents):
        """Writes a list file and returns its entries as dashing2 names them (the trimmed line)."""
        lines, shown = [], []
        for e in ents:
            sep = " " if rng.random() < 0.7 else "  "
            ln = sep.join(e)
            if rng.random() < 0.2:
                ln = " " + ln + " "
            lines.append(ln)
            shown.append(ln.strip())
            if rng.random() < 0.2:
                lines.append("")
        eol = "\r\n" if rng.random() < 0.3 else "\n"
        with open(os.path.join(wd, name), "w", newline="") as f:
            f.write(eol.join(lines) + eol)
        return shown
    mode = rng.choice(["--set", "-J"])
    weighted = mode == "-J"
    fl = ["-k", str(k), mode] + ([] if canon else ["--no-canon"])
    meas = rng.choice(list(MEAS_FLAGS))
    if layout == "pos":
        if len(entries) < len(files):
            entries = [[f] for f in files]
            names = files
            erecs = recsets
            counts = [kmer_counts(r, k, canon) for r in erecs]
        args = fl + MEAS_FLAGS[meas] + ["--square"] + names
        cols = list(range(len(names)))
    elif layout == "F":
        names = listfile("list.txt", entries)
        args = fl + MEAS_FLAGS[meas] + ["--square", "-F", "list.txt"]
        cols = list(range(len(names)))
    else:
        names = listfile("list.txt", entries)
        q = rng.sample(range(len(entries)), rng.randint(1, len(entries)))
        listfile("q.txt", [entries[i] for i in q])
        args = fl + MEAS_FLAGS[meas] + ["-F", "list.txt", "-Q", "q.txt"]
        cols = q
    try:
        M = cmp_matrix(R, args, rows=names)
    except D2Error as e:
        res.check(False, "inputs-run", str(e))
        return
    for r in range(len(entries)):
        for ci, c in enumerate(cols):
            t = pair_truth(counts[r], counts[c])
            t["k"] = k
            exp = measure_value(meas, t, weighted)
            got = M.get(r, ci)
            res.check(close(got, exp, rel=1e-5), "inputs-%s[%d,%d]" % (meas, r, c),
                      "dashing2 %s expected %s (A=%d B=%d I=%d)" % (got, exp, t["A"], t["B"], t["I"]))
    # Exact sets and counts of each entry against the counts.
    mk = Mask(0)
    for i, e in enumerate(entries):
        if len(e) > 1:
            continue
        cnt, aligned = kcounts(R, ["-k", str(k)] + ([] if canon else ["--no-canon"]), e[0])
        exp = {mk.m64(DNA.encode([DNA.code[ch] for ch in s])): v for s, v in counts[i].items()}
        got = {x: int(v) for x, v in cnt.items()}
        res.check(got == exp, "inputs-countdict-%d" % i, "dashing2 %d k-mers (total %d), truth %d (total %d)" % (
            len(got), sum(got.values()), len(exp), sum(exp.values())))


# ---------------------------------------------------------------------------


def run_one(args, d2, kmc, preset):
    def trial(idx, rng, wd):
        fam = FAMILIES[idx % len(FAMILIES)]
        res = TrialResult(idx, fam)
        R = Runner(d2, wd)
        try:
            if fam == "inputs":
                fam_inputs(rng, wd, R, res, preset, kmc)
            else:
                globals()["fam_" + fam](rng, wd, R, res, preset)
        except D2Error as e:
            res.fail("%s-run" % fam, str(e))
        if not res.tags:
            res.tags = {"family": fam}
        res.desc += " [%d dashing2 calls]" % R.calls
        return res
    return trial


def main():
    ap = common_args(__doc__.split("\n")[0], PRESETS)
    ap.add_argument("--family", action="append", choices=FAMILIES, help="run only trials of this family")
    args = ap.parse_args()
    preset = PRESETS[args.preset]
    d2 = Dashing2(args.dashing2)
    kmc = None if args.no_kmc else KMC(args.kmc_bin)
    if kmc is not None and not kmc.available:
        print("SKIP KMC comparisons: kmc/kmc_tools/kmc_dump not found (set KMC_BIN); using the Python oracle")
        kmc = None
    ntrials = preset["per_family"] * len(FAMILIES)
    if args.family and not args.trial:
        args.trial = [i for i in range(ntrials) if FAMILIES[i % len(FAMILIES)] in args.family]
    if not have_zstd():
        print("SKIP .zst inputs: zstd not found")
    results, el = run_trials(SUITE, SCRIPT, args, ntrials, run_one(args, d2, kmc, preset))
    return summarize("suite D (input processing and k-mer generation)", results, el, known=KNOWN)


if __name__ == "__main__":
    sys.exit(main())
