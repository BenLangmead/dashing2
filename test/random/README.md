# Randomized test suites for dashing2

Three randomized, seed-reproducible suites that cross-validate dashing2 against
KMC3 and against exact Python oracles. They need only Python 3 (standard
library) and, for the KMC comparisons, KMC 3.x (`kmc`, `kmc_tools`,
`kmc_dump`).

| Script | What it checks |
|---|---|
| `suite_a_exact.py` | The exact modes (`--set`, `-J/--countdict`) agree exactly with KMC3. |
| `suite_b_kcusp.py` | No cusp in k: exact modes are exact at every k, and sketch errors do not shift across the encoder boundaries. |
| `suite_c_sketch.py` | Sketch estimates agree with the exact values within their standard errors. |
| `d2rand.py` | Shared module: input generators, exact oracles, KMC and dashing2 wrappers, error model, reporting. |

## Running

```sh
export DASHING2=/path/to/dashing2          # or pass --dashing2
export KMC_BIN=/path/to/kmc/bin            # or put kmc on PATH; optional
python3 test/random/suite_a_exact.py                     # quick preset, master seed 1
python3 test/random/suite_b_kcusp.py --preset thorough --seed 7
python3 test/random/suite_c_sketch.py -j 8               # trials run 8 at a time
```

Common options: `--preset quick|thorough` (suite C also has `calibrate`),
`--seed N` (master seed), `--trial I` (run one trial; repeatable), `-j N`
(parallel trials, default 4), `--keep` (keep the inputs of failing trials),
`--workdir DIR`, `--no-kmc`, `--strict` (see known issues), `-v`.

Each suite exits nonzero on any failure. Every failing trial prints its
parameters and a one-line reproduction command, for example

```
[trial 27] FAIL (11 checks, 2.1s) k=21 --set --no-canon layout=sym -p5 seed=670378 m=1 base=2012 files=[...]
    FAIL union[0,1]: dashing2 0.0 expected 4046 (A=1915 B=2131 I=0 U=4046)
    repro: DASHING2=/.../dashing2 KMC_BIN=... python3 test/random/suite_a_exact.py --seed 1 --preset quick --trial 27 -v --keep
```

A trial is a pure function of (suite, master seed, trial index): its random
generator is seeded with the string `"<suite>:<seed>:<index>"`, so rerunning
with the same `--seed`, `--preset` and `--trial` regenerates the same inputs
and options, and `--keep` leaves them on disk. Aggregate failures (suites B
and C) depend on all trials and are reproduced by rerunning the whole preset
with the same seed. At the end each suite prints, for every randomized
feature (k range, mode, `-2`, threads, ...), how many trials with each value
failed and which checks failed most, which usually points at the cause.

Runtimes on an Apple M-series laptop (12 cores, `-j 4`): quick presets about
30 s (A), 60 s (B) and 45 s (C); thorough presets are listed at the end of
this file.

dashing2 writes `*.kmerset64/128`, `*.kmercounts.f64` and other files next to
its inputs even without `--cache`; every trial runs in its own temporary
directory and the suites delete everything dashing2 wrote after each call.

## Inputs

Every trial draws its own inputs: a random base genome (length randomized),
derived genomes with SNPs and small indels at random rates, slices, tandem
repeats and dispersed (sometimes reverse-complemented) copies so that k-mer
counts exceed 1, homopolymer runs near length k (suite B), split into 1 to 6
records, plus records shorter than k, runs of non-ACGT characters (N, n,
IUPAC codes, `-`, `.`), lowercase stretches and reverse-complemented records.
Files are FASTA (one line or wrapped) or FASTQ, optionally gzipped.

U/u is not used: dashing2 reads it as T (RNA alias) while KMC treats it as an
invalid character, a documented semantic difference rather than a bug.

## Exact values

`d2rand.kmer_counts` counts k-mers over A/C/G/T in either case, treating
every other character as a break, and canonicalizes to the lexicographically
smaller of a k-mer and its reverse complement (KMC's form). For very large k
`kmer_counts_hashed` keys k-mers by a 61-bit rolling polynomial hash instead
(linear time in k). From count tables the oracle derives |A|, |B|, |A n B|,
|A u B|, the weighted analogues (sum of counts, sum of min counts, sum of max
counts) and the probability Jaccard index (Moulton and Jiang), and from those
every dashing2 measure: similarity I/U, `--intersection` I, `--union-size` U,
`--containment` I/|row|, `--symmetric-containment` I/min(|A|, |B|) and
`--mash-distance` -ln(2J/(1+J))/k. For `-J` the same formulas use the
weighted quantities (similarity is the weighted Jaccard index, sum of min
counts over sum of max counts). Ratios with an empty operand follow dashing2's
definition: 1 for two empty inputs, 0 otherwise.

When KMC is available (suites A and B, k <= 256), the counts come from KMC:
`kmc -k<k> -ci1 -cs4294967295 -fm|-fq [-b]`, read back with `kmc_dump`, and
pairwise |A n B|, |A u B| and the sums of min and max counts come from
`kmc_tools simple A B intersect -ocmin union -ocmax`. Every trial also checks
that KMC's dump equals the Python oracle and that kmc_tools agrees with the
dumps, so the oracle is validated against KMC continuously. When KMC is
missing a `SKIP` line says so and the Python oracle is used.

KMC behaviours the wrappers rely on, each verified on small examples:
the default minimum count is 2, so `-ci1` is required; `-fm` reads multi-line
FASTA (`-fa` expects one line per record); `-b` turns off canonical k-mers;
lowercase acgt are valid and every other character breaks k-mers; `kmc_dump`
3.2.4 may exit with status 1 and no message on an empty database (accepted
when the dump file is empty); when the inputs of `kmc_tools simple` are
filtered with `-ci2` or more, an output without its own `-ci1` comes out empty,
so every output gets an explicit `-ci1`.

For k above the direct-encoding limit (32, or 64 with `-2`) dashing2 stores
64-bit (128-bit with `-2`) polynomial hashes of the k-mers rather than the
k-mers themselves, so its exact modes are exact only up to hash collisions.
With at most about 10^5 distinct k-mers per trial the chance of any collision
is about n^2 / 2^62, below 10^-8, so the suites require exact equality.
Values are compared with a relative tolerance of 1e-5 because dashing2 prints
float32; a single differing k-mer changes every measure by far more than that
at these input sizes.

## Suite A: exact modes match KMC3

Each trial draws k (1 to 256, weighted toward 13 to 64), canonical or
`--no-canon`, `--set` or `-J`, with or without `-2`, 2 to 4 files, the output
layout (upper-triangular default, `--square`, or a `-Q` panel with a random
query subset), `-p 1` or `-p 2..8`, an optional `--seed`, an optional
`-m 2|3` count threshold (KMC side: `kmc_dump -ci`, `kmc_tools -ci`), and
positional inputs or `-F`. It checks similarity, intersection, union,
containment, symmetric containment and Mash distance for every cell of the
output against the KMC-derived values, and the cardinality written by
`dashing2 sketch -o` (`<out>.names.txt`: distinct k-mers for `--set`, total
count for `-J`). Quick: 48 trials; thorough: 400 trials with larger genomes.

## Suite B: no cusp in k

dashing2 encodes DNA k-mers directly in 64 bits for k <= 32 and in 128 bits
with `-2` for k <= 64, and switches to a rolling polynomial hash (one 61-bit
lane, or two with `-2`) above that; those are the only k-dependent code paths
for DNA (src/d2.cpp, src/fastxsketch.cpp, bonsai encoder.h). Each trial takes
one k from a list that is dense around 30-34, 62-66, 126-130, 254-258 and
510-514 and also includes 1-4, 7, 8, 15-17, 21, 47, 96, 200, 300, 1000, 4097
and 100000 (thorough: every k up to 79, 120-135, 248-263, more large values up
to 10^6). Per k it runs `--set` and `-J`, canonical and `--no-canon`, with and
without `-2` (random `-p` and `--seed`), checks similarity, intersection,
union, containment (`--square`) and the cardinality exactly against KMC
(k <= 256) or the Python oracle (k > 256, hashed oracle for k > 2000), and
then sketches with OPH, `--full` and `-B` (with and without `-2`, 3 seeds per
k quick, 8 thorough). The sketch z-scores (similarity and cardinality, error
model below) feed boundary checks: for every boundary b in 32, 64, 128, 256,
512 the z-scores for k in [b-4, b] and in (b, b+4] must have means that differ
by less than 3 standard errors and SDs in 0.5-1.6 (when a side has 20 or more
distinct (k, replicate, mode) units), per mode and pooled over modes, and the pooled z over all k must pass
the suite C aggregate checks. Inputs include homopolymer runs (all-A and all-T
k-mers are the extreme encodings, 0 and all ones, at k = 32 and 64) and
records of length k-1, k and k+1.

### Largest k

dashing2 parses `-k` with `atoi` into an `int`. Every k from 1 to
2147483647 (INT_MAX) is accepted, but the spaced-seed machinery allocates
about 8 bytes per unit of k before reading any input (k = 10^8 uses 0.8 GB,
k = 2^31-1 uses 9.6 GB of resident memory), so the practical limit is memory.
Exact modes stay exact far beyond KMC's 256 (tested to 10^6 in the thorough
preset). Values beyond INT_MAX are not rejected: `atoi` overflows, so
`-k 2147483648` (and any negative k) silently falls back to the default k
(32, or 64 with `-2`), `-k 4294967296` becomes 0 and aborts with a spacer
error after allocating about 10 GB, and `-k 4294967297` silently runs with
k = 1. `-k 0` is rejected. Behaviour is the same with `-2`. (The help text
also says the direct limit is 31 for DNA; the code uses 32.)

## Suite C: sketch estimates agree with the exact values

Each trial picks a sketch mode round robin from 27 configurations: OPH
(default), `--full`, `-B` (BagMinHash, weighted Jaccard), `--prob`
(ProbMinHash, probability Jaccard), each with `--fastcmp 1|2|4` (registers
log-compressed after sketching; `--prob` uses 1 and 2, and dashing2 turns its
registers, which hold hashes of the selected items, into b-bit signatures
instead) and `--bbit-sigs --fastcmp 1|2|4` (b-bit signatures; not for
`--prob`), and the directly sketched compressed
SetSketches `--full --fastcmp-bytes|-shorts|-words`. It draws k (12 to 160),
`-S` (128 to 4096, including sizes that are not powers of two), `--seed`,
`-p`, `-2` (25%) and `--no-canon` (15%), three related inputs, and runs every
measure with `--square` plus `dashing2 sketch -o`. Checks:

* per estimate: |z| <= 6, z = (estimate - exact) / SE;
* the cardinality from the diagonal of `--union-size --square` and from
  `<out>.names.txt` (equal to each other for set modes; exact totals for `-B`
  and `--prob`, and exact counts for sets of at most 10 x S distinct k-mers);
* deterministic identities between printed measures: similarity = I/U,
  I + U = |A| + |B|, containment = I/|row|, symmetric containment =
  I/min(|A|, |B|), Mash distance = f(similarity);
* thread independence: with a fixed `--seed`, `-p N` prints exactly what
  `-p 1` prints;
* aggregates over all trials, per mode and measure and pooled over modes:
  |mean z| <= 3/sqrt(trials) + allowance, and SD of z in 0.5-1.6 for groups
  of at least 30 trials (with fewer, one trial's correlated z-scores can
  dominate the SD; in the quick preset only the pooled groups qualify).

Only "regular" pairs (at least 5 expected equal and 5 expected unequal
registers, m J >= 5 and m (1 - J) >= 5) enter the aggregates; pairs at J = 0
or 1 still get the per-estimate bound. For `--prob` only similarity and Mash
distance have an exact counterpart; the other measures are not checked.

### Error model

Every estimate dashing2 prints is a function of three noisy quantities: the
register-agreement estimate of J (set, weighted or probability Jaccard) and
the cardinality estimates A and B (src/cmp_core.cpp):
similarity J, intersection J (A + B)/(1 + J), union (A + B)/(1 + J),
containment I/A, symmetric containment I/min(A, B), Mash -ln(2J/(1+J))/k.
The suites model (J, A, B) for m registers as m independent min-hash
registers whose minima come from A \ B, B \ A and A n B (exponential minima;
the correlations below follow from that model):

* Var J = p (1 - p) / (m (1 - c)^2) with p = J + (1 - J) c, where c is the
  probability that two registers holding different minima compare equal after
  compression: 0 for full registers, 2^-bits for b-bit signatures, ln(b)/4
  for registers log-compressed with base b (two independent Gumbel log-minima
  differ by a logistic variable; the chance that they fall in the same bucket
  of width ln b is about ln(b)/4 for equal cardinalities and smaller
  otherwise). dashing2 prints the base it fits for `--fastcmp N`; the presets
  fix b = 1.2, 1.0005 and 1.0000000109724. `--prob --fastcmp N` uses b-bit
  signatures with 8N bits.
* Var A = kappa^2 (A^2 / m + 1), Var B likewise, Corr(A, B) = J,
  Cov(J, A) = kappa J A |B \ A| / (|A u B| m), with kappa = 1 for OPH and
  FullSetSketch and kappa = 0 for `-B` and `--prob`, whose cardinalities are
  exact total counts. The "+ 1" is one element of resolution: for sets much
  smaller than m, OPH effectively counts occupied registers and errs by whole
  register collisions (2 k-mers in 256 registers collide 1 time in 256 and
  read as 1), which a continuous 1/sqrt(m) model would score as z = -8.
* Small sets are counted exactly: OPH, `--full` and their compressed forms
  report the number of distinct sketch ids as the cardinality when it is at
  most 10 x m (src/fastxsketch.h `SmallSetCounter`, in file mode and with
  `--parse-by-seq`), so kappa is 0 for such a set (`card_kappa` and
  `EXACT_CARD_FACTOR` in d2rand.py). Its cardinality (the `--union-size`
  diagonal and names.txt) must then equal the exact count up to float32
  printing and stays out of the aggregates, and its A or B terms drop out of
  the SE of the derived measures, which keeps the spread check meaningful
  for them. Larger sets keep the calibrated model below.
* Derived measures get their SE by the delta method. Two floors are added in
  quadrature: one register's worth of resolution (the change of the measure
  when J moves by 1/m), a continuity correction for the discrete count of
  equal registers near J = 0 or 1, and 2e-6 relative for float32 printing.
* Mash distances are inverted to J (J = 1/(2 e^{kd} - 1)) and compared on
  the J scale.

Suite C uses m >= 128 because cardinality estimators of this kind behave like
(m - 1)/Gamma(m) and have a heavy right tail: at m = 64 a 6-sigma normal bound
would be exceeded with probability about 2e-5 per check, at m >= 128 well
below 1e-6.

### Calibration

The formulas were not trusted blindly. Before the suite was written, 300
trials per mode were collected and the z distributions inspected; that showed
the need for the collision term c (1-byte modes had SD of z up to 3.7 at J
near 0 without it) and for the one-register floor (z up to 10 when a single
register matched by chance at J = 0 with 2-byte registers). The final model
was then checked with the `calibrate` preset (200 trials per mode, master
seed 101, 5400 trials, 14 minutes with `-j 6`):

```sh
python3 test/random/suite_c_sketch.py --preset calibrate --seed 101 -j 6
```

Mean z / SD z over regular pairs (trials matching a known issue, below,
excluded); "exact" means the cardinality is an exact total. This calibration
was run on `test/all-fixes`, before the fixes listed under "Fixed issues"
below; in particular `pmh-fc1` and `pmh-fc2` still log-compressed the
`--prob` registers, and trials of the `full-ab-*` modes with -p > 1 or a
sketch size that does not fill whole 64-bit words were excluded.

| mode | trials | similarity | intersection | union | containment | symcontain | card | max abs z (all pairs) |
|---|---|---|---|---|---|---|---|---|
| oph | 200 | +0.09 / 0.96 | +0.12 / 0.90 | +0.03 / 0.81 | +0.09 / 0.95 | +0.13 / 0.95 | +0.06 / 0.85 | 3.3 |
| oph-fc1 | 199 | -0.19 / 0.93 | -0.10 / 0.91 | +0.09 / 0.91 | -0.16 / 0.94 | -0.12 / 1.01 | +0.04 / 0.87 | 3.0 |
| oph-fc2 | 200 | -0.07 / 0.99 | +0.02 / 0.95 | +0.14 / 0.93 | -0.04 / 0.96 | +0.01 / 0.98 | +0.06 / 0.85 | 3.1 |
| oph-fc4 | 200 | -0.00 / 0.91 | -0.04 / 0.88 | -0.03 / 0.87 | -0.02 / 0.91 | +0.02 / 0.92 | -0.01 / 0.85 | 2.8 |
| oph-bb1 | 200 | -0.02 / 0.96 | -0.03 / 0.91 | +0.04 / 0.82 | -0.03 / 0.96 | +0.01 / 0.91 | +0.02 / 0.83 | 3.4 |
| oph-bb2 | 200 | +0.01 / 0.91 | -0.03 / 0.87 | -0.02 / 0.94 | +0.00 / 0.92 | +0.02 / 0.94 | -0.00 / 0.87 | 3.6 |
| oph-bb4 | 200 | +0.09 / 0.95 | +0.05 / 0.95 | -0.04 / 0.85 | +0.09 / 0.96 | +0.13 / 0.97 | +0.00 / 0.87 | 3.2 |
| full | 200 | -0.02 / 1.01 | +0.00 / 0.93 | +0.10 / 0.93 | -0.02 / 0.99 | +0.02 / 0.96 | +0.07 / 0.95 | 3.9 |
| full-fc1 | 200 | +0.04 / 1.02 | +0.05 / 1.04 | +0.05 / 1.06 | +0.03 / 1.00 | +0.10 / 0.96 | +0.07 / 1.02 | 3.3 |
| full-fc2 | 199 | -0.07 / 0.97 | -0.06 / 0.99 | +0.01 / 1.05 | -0.06 / 0.98 | +0.00 / 0.97 | -0.01 / 1.03 | 3.2 |
| full-fc4 | 200 | +0.10 / 0.97 | +0.04 / 0.93 | -0.08 / 0.96 | +0.08 / 0.98 | +0.11 / 0.97 | -0.01 / 0.95 | 3.9 |
| full-bb1 | 200 | -0.09 / 0.95 | -0.08 / 0.90 | +0.04 / 1.04 | -0.07 / 0.95 | +0.01 / 0.94 | +0.02 / 0.95 | 3.9 |
| full-bb2 | 200 | -0.14 / 1.01 | -0.18 / 1.02 | -0.01 / 1.01 | -0.15 / 1.01 | -0.15 / 0.98 | -0.01 / 1.00 | 4.4 |
| full-bb4 | 200 | +0.03 / 1.02 | +0.05 / 1.01 | +0.10 / 1.02 | +0.03 / 1.02 | +0.07 / 1.02 | +0.06 / 1.01 | 3.4 |
| full-ab-bytes | 73 | -0.19 / 0.98 | -0.13 / 1.14 | +0.14 / 1.11 | -0.14 / 1.00 | -0.04 / 1.02 | -0.01 / 1.13 | 3.5 |
| full-ab-shorts | 90 | +0.15 / 1.10 | +0.13 / 1.02 | -0.01 / 1.07 | +0.15 / 1.05 | +0.24 / 1.00 | -0.04 / 1.02 | 8.8* |
| full-ab-words | 104 | +0.15 / 1.08 | +0.16 / 1.12 | +0.08 / 1.05 | +0.11 / 1.05 | +0.14 / 1.02 | +0.09 / 0.99 | 3.2 |
| bmh | 197 | +0.07 / 0.97 | +0.06 / 0.97 | -0.06 / 0.97 | +0.06 / 0.97 | +0.06 / 0.97 | exact | 3.1 |
| bmh-fc1 | 191 | +0.02 / 0.97 | +0.00 / 0.97 | -0.00 / 0.97 | +0.00 / 0.97 | +0.00 / 0.97 | exact | 3.4 |
| bmh-fc2 | 193 | +0.03 / 1.07 | +0.02 / 1.07 | -0.02 / 1.07 | +0.02 / 1.07 | +0.02 / 1.07 | exact | 3.0 |
| bmh-fc4 | 186 | -0.01 / 1.07 | -0.02 / 1.06 | +0.02 / 1.06 | -0.02 / 1.06 | -0.02 / 1.06 | exact | 3.8 |
| bmh-bb1 | 186 | +0.01 / 0.94 | +0.00 / 0.95 | -0.00 / 0.95 | +0.00 / 0.94 | +0.00 / 0.95 | exact | 3.2 |
| bmh-bb2 | 192 | +0.01 / 0.97 | +0.00 / 0.97 | -0.00 / 0.97 | +0.00 / 0.97 | +0.00 / 0.97 | exact | 3.8 |
| bmh-bb4 | 189 | -0.07 / 0.99 | -0.09 / 0.99 | +0.09 / 0.99 | -0.09 / 0.99 | -0.09 / 0.99 | exact | 3.0 |
| pmh | 191 | +0.05 / 0.95 | n/a | n/a | n/a | n/a | exact | 3.5 |
| pmh-fc1 | 192 | -0.03 / 1.03 | n/a | n/a | n/a | n/a | exact | 5.4 |
| pmh-fc2 | 189 | -0.02 / 0.96 | n/a | n/a | n/a | n/a | exact | 2.8 |

Pooled over modes: similarity -0.003 / 0.986 (10274 z, 4746 trials),
intersection -0.005 / 0.973, union +0.024 / 0.970, containment -0.008 / 0.983,
symmetric containment +0.022 / 0.980, cardinality +0.025 / 0.935.
(*) The 8.8 for `full-ab-shorts` is a tiny-input pair affected by the
`ab-card-floor` issue (fixed since). The OPH cardinality SD of 0.85 is below 1 because
the OPH estimator is more precise than 1/sqrt(m) for sets smaller than m.

Thresholds chosen from this:

* per estimate |z| <= 6 (`--zmax`): the largest |z| in 5400 calibration
  trials (about 377000 checks) was 5.4, and under a normal model a 6-sigma
  excursion has probability 2e-9 per check;
* bias: |mean z| <= 3/sqrt(trials) + allowance for the pooled groups (4/sqrt(trials) for the per-mode groups, of which a run has well over a hundred, so that chance failures stay rare), allowance 0.1, raised to
  0.15 for cardinality, union and intersection (cardinality estimators that
  invert a sum of m minima carry an O(1/m) upward bias) and to 0.25 for
  symmetric containment (it divides by the smaller of two noisy
  cardinalities, a Jensen bias); the calibration means stay within +-0.2;
* spread: SD of z in 0.5-1.6 (groups of at least 30 trials); calibration
  SDs are 0.81-1.14.

## Known issues in the tested code (expected failures)

None at present. The `XFAIL[id]` mechanism (and `--strict`, which turns expected failures into failures) remains available for registering a newly found discrepancy without loosening tolerances.

## Fixed issues

The suites found these discrepancies in `test/all-fixes` and registered them
as expected failures; each is fixed on its own branch (merged into
`test/all-fixes-2`), and the suites now check the affected cases normally.

* `emit-row-loss` (`fix/emit-queue-lock`): rarely, `dashing2 cmp` exited 0 but printed a matrix with a row missing (seen twice in about 30 suite B runs, never reproduced on demand). The writer thread in src/emitrect.cpp read the output deque without holding its lock while the comparison threads appended to it. It now takes items off the deque under the lock. ThreadSanitizer reported races in that function in 8 of 8 runs before the change and none after. The suites still check every printed row, and a missing row now counts as a failure.

* `ab-pad` (`fix/ab-sketchsize-padding`): the directly sketched compressed
  SetSketch modes (`--fastcmp-bytes`, `--fastcmp-shorts`, `--fastcmp-words`
  and `--setsketch-ab A,B --fastcmp N`) gave near-zero similarity whenever
  the sketch size times the register width was not a multiple of 8 bytes,
  because the `Dashing2DistOptions` constructor (src/cmp_main.h) padded the
  caller's options instead of its own sketch size. The sketch size is now
  padded to whole 64-bit words.
* `ab-race` (`fix/ab-thread-safety`): those modes gave results that changed
  from run to run with `-p > 1` (and occasionally segfaulted with
  `--fastcmp-shorts`). `SetSketch::harmean()` shared a static table of powers
  between threads, and the per-thread sketches kept their skip cutoff and
  (with `-m`) their pending k-mer counts from one file to the next. The table
  is gone and `clear()` resets all of that state.
* `ab-card-floor` (`fix/ab-cardinality-floor`): those modes could not report
  a cardinality below about 1/a (an empty input read as 0.046, 0.051 or 16.7,
  and 10 k-mers as 22 with `--fastcmp-shorts`). `SetSketch::cardinality()`
  now treats registers still at 0 as censored, so an empty input gives 0 and
  small inputs are unbiased; estimates without zero registers are unchanged.
* `fc-all-empty` (`fix/fastcmp-all-empty`): with `--fastcmp N` and no k-mer
  in any input, the fitted a and b were NaN and every measure printed inf.
  `make_compressed()` (src/cmp_core.cpp) now uses the preset a and b when
  there is nothing to fit, so the output equals the uncompressed output.
* `prob-fc1-bias` (`fix/prob-fastcmp-collisions`): `--prob --fastcmp 1`
  did not remove chance register agreements (similarity about 0.007 at
  J_P = 0 for small inputs), because `--prob` registers hold hashes of the
  selected items, to which the SetSketch collision correction of log
  compression does not apply. `--prob --fastcmp N` now compresses them to
  b-bit signatures, whose chance agreements are removed exactly.

Minor observation not treated as a failure: with an empty operand, compressed
modes may report a nonzero intersection (one chance register match, e.g.
0.88 for a 880-k-mer set at S = 1000) while the ratio measures are defined as
0; the identity checks skip pairs with an empty operand.


The wider suites D, E and F found further issues; these are fixed in the tested code (`test/all-fixes` plus the branches named here). The suites still detect each by its id, and `FIXED_IDS` in d2rand.py makes a recurrence a failure reported as a regression rather than an expected failure.

* `symcontain-as-distance` (`fix/symcontain-similarity`): `--topk`, `--similarity-threshold` and `--greedy` ranked and thresholded `--symmetric-containment` as a distance.
* `binary-tail-32k` (`fix/binary-tail-32k`): the final flush of queued `cmp` output dropped a block whose size was an exact multiple of the 128 KiB write chunk, in both binary and text output, with exit status 0.
* `lsh-self-candidate` (`fix/lsh-self-candidate`): the LSH candidate cap counted the query itself, so a qualifying pair could be missing from `--topk` and `--similarity-threshold` output.
* `countmin-exact-modes` (`fix/countmin-exact-modes`): `-c/--countmin-size` made `--set` and `-J` inexact; exact modes now ignore it with a note on stderr.
* `countsketch-header-newline` (`fix/countsketch-header-newline`): with `-c` the `#Dashing2Options` header line was split in two.
* `readme-edit-distance-crash` (`fix/sketch-compute-edit-distance`): `dashing2 sketch --parse-by-seq --edit-distance --compute-edit-distance` (README Use 6) crashed because the sequences were freed before comparison.
* `byseq-downsample` (`fix/byseq-downsample`): `--parse-by-seq` ignored `--downsample`.
* `byseq-compressed` (`fix/byseq-compressed`): `--parse-by-seq` read `.bz2`, `.xz` and `.zst` inputs as raw bytes.

* `edit-distance-nondeterministic` (`fix/edit-distance-nondeterministic`, with bonsai and sketch branches `fix/omh-finalize-index-type`): `--edit-distance` (OrderMinHash) gave different values for the same pair from run to run, because `OMHasher::finalize` in the sketch library stored sequence positions in the element type (one byte for `const char*`), so positions of 128 or more wrapped and were read from before the record.

The remaining issues of suites D, E and F are fixed by the branches merged into `test/all-fixes-7` (dashing2, with bonsai and sketch branches of the same name where the fix is in those libraries). Details of each issue are in the suite sections below.

* Suite D: `uncanon-window-polyT` (bonsai `fix/uncanon-window-polyT`), `spaced-window-invalid` (bonsai `fix/spaced-window-invalid`), `byseq-seq-stale-window` (bonsai `fix/rolling-uncanon-window-reset`), `filterset-windowed` (`fix/filterset-windowed`), `spacing-long-filename` (`fix/spacing-long-filename`), `byseq-exact-card` (`fix/byseq-exact-card`: file mode now counts small sets exactly too, see the error model of suite C).
* Suite E: `stacked-kmercounts-f32` (`fix/stacked-kmercounts-f64`), `full-mincount-counts` (`fix/full-mincount-counts`, with bonsai and sketch branches), `bmh-stale-ids` (`fix/bmh-reset-ids`, with bonsai and sketch branches), `seq-o-file-mode` (`fix/seq-o-file-mode`), `wsketch-tw-suffix` (`fix/wsketch-tw-suffix`).
* Suite F, behaviour: `wsketch-U-1d` (`fix/wsketch-1d-u32-weights`), `wsketch-tw-garbage` (`fix/wsketch-tw-suffix`), `header-sketchtype-label` (`fix/header-sketchtype-label`), `short-opts-differ` (`fix/cmp-no-canon-short`), `fastcmp-presets-need-full` (`fix/fastcmp-presets-default-sketch`: the presets and `--setsketch-ab` select `--full`), `seqs-in-ram-ignored` (`fix/seqs-in-ram`), `presketched-mash-k` (`fix/presketched-k-warning`: cmp warns that `--presketched --mash-distance` uses the default k unless `-k` is given, and the help says so).
* Suite F, usage text (`docs/help-text`): `help-k-limits`, `default-not-phylip`, `topk-large-k`, `sparse-binary-format`, `greedy-binary-ids`, `save-kmers-names-path`, `save-kmers-header-size`, `save-kmers-no-file`, `seq-header`, `wsketch-example-order`, `wsketch-example-kmerset-header`, `undocumented-options` (every option is documented except the internal `--by-chrom`, `--exact-kmer-dist`, `--pairlist` and `--sig-ram-limit`).
* Suite F, README (`docs/readme`): `readme-use7-no-output`, `readme-canon-default`, and the README side of `default-not-phylip`.

## Runtimes

Measured on the fixed build with `-j 4` on a 12-core Apple Silicon laptop
(KMC 3.2.4 present):

| preset | suite A | suite B | suite C |
|---|---|---|---|
| quick | 48 trials, about 30 s | 42 k values, about 60 s | 270 trials, about 45 s |
| thorough | 400 trials, 216 s | 133 k values, 358 s | 2160 trials, 388 s |
| calibrate | | | 5400 trials, 861 s (`-j 6`) |

KMC's own run time (about 0.85 s per `kmc` call regardless of input size,
mostly idle in its second stage) dominates suites A and B; the wrappers run
the counts of one trial side by side.

## Suite D: input processing and k-mer generation

`suite_d_inputs.py` (helpers in `d2rand_d.py`) checks the modes that decide
which items reach a sketch. Each trial belongs to one of ten families, chosen
round robin by trial index (`--family NAME` runs one family), and draws its
own inputs (random records with repeats, N and IUPAC runs, lowercase, RNA U,
homopolymers, records shorter than k, empty records) and options. Quick: 20
trials per family (200 trials, about 30 s with `-j 6`); thorough: 100 per
family (1000 trials, about 85 s with `-j 8`).

`d2rand_d.py` reimplements the selection exactly, so most checks compare
dashing2's stored values item by item rather than counts: the DNA and protein
encodings (2 bits per base; base-20, base-14 and base-6 integers; 3 bits for
`--protein8`; first character most significant; aliases U as T, O as K, U as
C), the stored form (XOR with a mask derived from `--seed`, 0 by default,
then the Wang hash, per 64-bit half for `-2`), the window order (a window is
the last w - c + 1 valid k-mers of the record, c the k-mer span; invalid
k-mers are skipped; one item per full window; a record with fewer valid
k-mers emits its best one; best = smallest (score, k-mer), score = FRev64 or
the 128-bit CEHasher, for `--entmin` the top 48 bits of that hash divided by
the k-mer's Shannon entropy plus 1e-4) and the `--downsample` hash. The
oracle reproduces `sketch --set`, `-J` and `-G` output exactly on every
direct-encoding path. For k beyond the exact limit (rolling polynomial hash)
the checks are the number of distinct k-mers and the windowed stream as the
sliding minimum of dashing2's own unwindowed stream (score: FRev64 of the
hash, or the low 64 bits of CEHasher with `-2`).

| Family | What it checks |
|---|---|
| byseq | `cmp --parse-by-seq --square` equals the same records as one file each, cell by cell, for OPH, `--full`, `-B`, `--prob`, `--fastcmp 2`, `--bbit-sigs`, `--fastcmp-shorts` and three random measures, with windows, `--entmin`, `--spacing`, `--filterset`, `--downsample`, `--no-canon`, `-2`, protein, FASTQ/gzip; row names are the record names (duplicates kept); `--set`/`-J` are rejected cleanly; the `sketch -G --parse-by-seq -o` file holds each record's stream |
| minimizer | `--set`, `-J` (emission counts) and `-G` against the exact window oracle, k 1-64 with and without `-2`, canonical or not, w below, at and above k, `--seed`; density sanity; rolling path as above |
| entmin | exact `--entmin` streams (a pick within float rounding of the best score is accepted); picks are input k-mers; mean entropy of picks is at least that of plain minimizers; reverse-complement symmetry; no effect without `-w`; fallback to plain minimizers above the exact limit |
| spacing | spaced sets and streams (never canonicalized), run-length syntax equals the expanded list, all-zero spacing equals contiguous `--no-canon`, rejection above capacity and of a wrong number of gaps |
| protein | `--protein`/`--protein20`/`--enable-protein`, `--protein14`, `--protein8`, `--protein6` exact sets up to the 64/128-bit capacity (14/29, 16/33, 21/42, 24/49), distinct counts above it; `--no-canon` irrelevant; O and U aliases; B, Z, X, J and `*` break k-mers |
| seqmode | file-mode `.mmerseq` streams (with `--hp-compress`), `cmp -G` (every measure prints the longer stream length minus the item edit distance; `--compute-edit-distance` prints the edit distance), the `-G` cardinality, and the by-seq stream file layout |
| filter | `--filterset` from FASTA, from a binary file of stored values (`path:b`, also gzipped), with `-2`, windows and the rolling path; filtered `-G` streams |
| downsample | kept set equals the hash oracle, kept fraction within 5 binomial SD, Jaccard of the downsampled sets within 5 SE of the full Jaccard, printed `--set` similarity, `-G` stream, values 0, 1 and out of range |
| hpseed | `--hp-compress` collapses repeats across record boundaries in file mode and within records with `--parse-by-seq`, leaves `--set` alone; stored values follow each seed's mask; exact-mode measures are identical for seeds 0, 1, 42, 2^63, 2^64 - 1 and random; seeded sketches are valid |
| inputs | FASTA (wrapped or not, with or without header comments), FASTQ, gzip, bzip2, xz, zstd, files in subdirectories, `-F` and `-Q` lists with blank lines, extra spaces and CRLF, joint entries (two files on one line, row name = the trimmed line); every measure and the `-J` counts against KMC3 (which is checked against the Python oracle) |

Format notes (behaviour, not failures): the `--parse-by-seq -G -o` file has a
20-byte header (n, k, w, then alphabet | canonical << 8, as the help
describes; w is 2^32 - 1 without `-w`) and stores the k-mers
themselves, while file-mode `.mmerseq` files store masked values; the
Hamming comparison in `wcompare.cpp` is unreachable because `-G` always
compares by edit distance.

### Issues found by suite D (fixed)

Each was registered with an emulation: a mismatch counted as an expected
failure only when dashing2's output equaled what the described defect
predicts, so any other discrepancy in the same check still failed. All are
fixed in the tested code (branches under "Fixed issues"); the emulations stay,
and a recurrence is reported as a regression of its id.

* `uncanon-window-polyT`: with `--no-canon` (or `--spacing`) and k filling
  the k-mer word (32, or 64 with `-2`), a window whose minimizer is poly-T
  emits nothing, because the all-ones encoding is also the "window not full"
  marker (bonsai encoder.h `for_each_uncanon_unspaced_windowed`,
  `for_each_uncanon_spaced`; qmap.h `next_value`). Repro: a 40-bp T run,
  `sketch -G --no-canon -k 32 -w 40` emits 354 items, the oracle 355.
* `spaced-window-invalid`: `--spacing` with `-w` pushes k-mers overlapping
  an invalid character into the window as the all-ones marker with a real
  score, so windows count positions, and a window won by the marker emits
  nothing; records shorter than the window emit nothing (encoder.h
  `next_minimizer`, `for_each_uncanon_spaced`). Repro: `-k 8 --spacing
  0,1,0,2,0,0,1 -w 20`, 5 records with N runs: 479 stream items and 95
  distinct, oracle 463 and 92.
* `spacing-long-filename`: exact modes and caches name their files with the
  whole seed (`2x1,3x1,...`, src/fastxmerge.cpp:94-95), so an irregular seed
  with large k exceeds the 255-byte name limit and `sketch --set -k 64 -2
  --spacing 1,2,1,2,...` aborts with "Failed to open" (exit 134).
* `byseq-exact-card`: `--parse-by-seq` with OPH or `--full` (also log or
  b-bit compressed, `--fastcmp-shorts`) replaces cardinalities below 10 x S by
  exact counts (src/fastxsketchbyseq.cpp, "exact counting fall-back"); file
  mode keeps the estimate, so intersection, union and the containments of a
  record differ from the same record as a file (e.g. 1186 vs 1136.09). The
  similarities agree exactly.
* `byseq-downsample`: `--parse-by-seq` ignores `--downsample`: the per-record
  callbacks never call `downsample_pass`. Repro: `cmp --parse-by-seq --full
  --downsample 0.5 -k 15` gives intersection 1186 for a record whose full set
  has 1186 k-mers (file mode 608.7).
* `byseq-compressed`: `--parse-by-seq` opens `.bz2`, `.xz` and `.zst`
  inputs with `gzopen`, which passes compressed bytes through: one garbage
  record with a binary name, exit 0 (file mode decompresses them).
* `filterset-windowed`: with `-w`, `--filterset FASTA` removes only the
  filter file's minimizers (src/d2.cpp builds the set with the windowed
  encoder), so input minimizers present in the filter file as non-minimizers
  survive. Repro: filter = 25-bp fragments of the input, `-k 15 -w 30`: 185
  of 223 minimizers kept, 163 expected.
* `byseq-seq-stale-window`: `sketch -G --parse-by-seq` on the rolling path
  without canonicalization (k > 32, or > 64 with `-2`; `--no-canon` or
  protein) gives a record shorter than k the previous record's last window
  minimum: `RollingHasher::for_each_uncanon` returns before resetting its
  window, and the by-seq fallback for short records reads it. Repro: a 200-bp
  record then a 10-bp record, `-k 40 -w 50 --no-canon`: the second stream has
  1 item (0 expected).

On stock v2.1.20 the quick preset (seed 1) fails 163 of 200 trials (391
checks): every family is hit, mainly by the fixed `--entmin`, `--downsample`,
`--spacing` capacity, protein encoding, `--parse-by-seq --set` and rolling
hash bugs.

## Suite E: workflow and output modes

`suite_e_workflows.py` (helpers in `d2rand_e.py`) checks the parts of
dashing2 around the comparison core: nearest-neighbour and clustering
outputs, every output layout, sketch and k-mer files and their round trips,
`contain`, `wsketch`, `printmin` and BED input. It uses the shared options
(`--preset quick|thorough`, `--seed`, `--trial`, `-j`, `--keep`, `--strict`,
`-v`) and adds `--family NAME` (repeatable) to run only the trials of some
families. It does not need KMC.

```sh
DASHING2=/path/to/dashing2 python3 test/random/suite_e_workflows.py               # quick, 150 trials
DASHING2=/path/to/dashing2 python3 test/random/suite_e_workflows.py --family nn --preset thorough -j 8
```

Trial i belongs to family `FAMILIES[i % 10]`. Inputs are clustered random
genomes (`d2rand_e.clustered_inputs`): a few random roots, members with SNP
rates from 0 to 15% and small indels, some members sliced (so containment is
asymmetric), some with repeats, some exact copies of earlier members (ties),
plus the record features of the other suites (several records, short
records, N runs, lowercase, FASTQ, gzip).

### Exact identifiers

dashing2 stores a DNA k-mer as maskfn(kmer) = WangHash(kmer ^ XORMASK), with
2 bits per base (A=0, C=1, G=2, T=3), the numerically smaller of the k-mer
and its reverse complement when canonical, and XORMASK = 0 for `--seed 0`
and WangHash(seed) otherwise; with `-2` each 64-bit half is hashed and the
pair is folded to 64 bits for sketches (`fold64`) or stored unfolded in
`.kmerset128`. `d2rand_e.kmer_id` reproduces this (verified against `--set`,
`-J`, `--save-kmers`, with seeds and `-2`), which turns every saved k-mer file
into an exact oracle for k <= 32 (k <= 64 with `-2`).

### Families and checks

* `nn`: `--topk K` and `--similarity-threshold T` against the `--square`
  matrix of the same mode and measure (OPH, `--full`, `-B`, `--prob`, `--set`,
  `-J`, `--fastcmp 1|2`, `--bbit-sigs --fastcmp 2`, `--fastcmp-shorts`;
  similarity, Mash distance, containment, symmetric containment,
  intersection). Every listed value must equal the matrix entry (binary CSR
  exactly, text to 8 digits), no self or duplicate entries, best first, ties at
  the K-th value only, binary equal to text, `-p N` equal to `-p 1`. Thresholds
  are chosen strictly between two observed values, so a listed pair must
  qualify. Misses are classified: when every input is a candidate
  (n - 1 <= 3.5 K for top-k; always for thresholds with up to 21 inputs) and
  the LSH keys are the compared registers (OPH, `--full`, `-B`, `--prob`,
  `--fastcmp-shorts`), any missed positive pair fails; in the other modes a
  missed pair fails only if its similarity is at least 0.25. All other misses
  (candidates capped at 3.5 K, bottom-k keys of exact modes, full-register keys
  of `--fastcmp`/`--bbit-sigs`) enter recall statistics printed at the end.
  A third of the trials use many identical inputs, which fills candidate lists
  with ties.
* `greedy`: `--greedy T` and `--greedy TE`. The clusters must partition the
  inputs, every member must be within T of its representative (similarity
  >= T, distance <= T), binary output must equal text and `-p N` must equal
  `-p 1`. `TE` is simulated exactly from the matrix (inputs in input order
  join the best existing cluster, ties to the lowest cluster index, float32
  scores). For the LSH variant, pairs of representatives within T of each
  other are counted (approximation, reported).
* `formats`: the default upper triangle, `--square`, `-Q` panel and
  `--phylip`, each as text, `--binary-output` on stdout, and `--cmpout` text
  and binary files: identical values (binary decoded per python/parse.py),
  correct labels and sizes, `--cmpout` files equal to stdout, panel equal to
  the square's submatrix, `sketch --cmpout` equal to `cmp`. A third of the
  trials compare 128 to 512 tiny inputs with `-p`/`--batch-size` chosen so
  that a queued block holds exactly 32768 floats.
* `roundtrip`: `sketch -o` stacked files (header, cardinalities equal to
  `.names.txt`, `cmp -o` equal to `sketch -o` except OPH, which `cmp`
  densifies in place), `cmp --presketched` on the stacked file and on per-input
  `--cache` files equal to the direct comparison (k is passed again: sketch
  files do not record it), a second `--cache` run reads the files without
  rewriting them, and `--outprefix` keeps every file in its directory.
* `kmers`: `--set`/`-J` files (sorted ids equal to the exact k-mer set,
  cardinality header, counts aligned with ids), with `-m`, `--no-canon`,
  `--seed` and `-2`; `sketch -N -o` k-mer databases (`.kmer64` header fields,
  names file, every sampled id a k-mer of its input, saved counts equal to the
  exact counts).
* `contain`: coverage and depth for OPH, `--full`, `-B` and `--prob`
  databases, canonical and `--no-canon`, seeds, `-2`, against a recomputation
  from the decoded database and exact query counts (binary output exactly;
  text, `-b -o` and `-p N`). With `-w`, minimizer streams come from `sketch -G`
  with the same options (itself checked exactly by `seq`). A quarter of the
  trials add a reference without k-mers.
* `wsketch`: CSR input with `-B`, `-q` and the default ProbMinHash: header,
  total weights, `info.txt`, every sampled id an id of its row, equal-register
  fractions against the exact weighted, probability or set Jaccard (|z| <= 6),
  `cmp --presketched` on the stacked registers, `-u`/`-f`/`-H`/`-P` inputs
  byte-identical to 64-bit inputs, and 1-D input equal to the same CSR row.
* `seq`: `sketch -G --parse-by-seq -o` (raw k-mer codes per record, with and
  without `--hp-compress`), `printmin` text and `-f` output, and the per-input
  `.mmerseq64` file of file mode (masked ids), all exactly.
* `bed`: BED files with comments, extra columns and overlapping or repeated
  intervals: similarity against the exact Jaccard of covered positions (or
  weighted Jaccard with multiplicities for `-B`), labels, `sketch --bed -o`
  stacked file, cardinalities and its `--presketched` reload.
* `measures`: for every measure, `--square` (with extra query-only inputs),
  `-Q` panels (including `--mash-distance -Q`), the default symmetric output,
  symmetry of the symmetric measures, Mash distance as a function of
  similarity, containment x |row| = intersection and symmetric containment x
  min(|A|, |B|) = intersection. `--fastcmp` fits its compression to all
  inputs of a run, so runs over different input sets are not compared there.

### Issues found by suite E (fixed)

All are fixed in the tested code (branches under "Fixed issues"); a
recurrence is reported as a regression of its id.

* `symcontain-as-distance` (high for its users): `--symmetric-containment`
  is ranked as a distance by `--topk`, `--similarity-threshold` and
  `--greedy`, because `distance()` (src/cmp_main.h:44) returns true for every
  measure other than similarity, containment, intersection and union. Top-k
  lists the least contained neighbours first, the threshold keeps pairs below
  T, and greedy joins unrelated inputs (symmetric containment 0) into one
  cluster. Repro on four genomes in two clusters:
  `dashing2 cmp -k 17 --full --symmetric-containment --similarity-threshold 0.5 ...`
  lists c0_1 with c0_3 (0.483) but not c0_0 (0.890); `--greedy 0.5E` puts
  all 8 unrelated genomes in c0_0's cluster.
* `lsh-self-candidate` (medium): `--topk` and `--similarity-threshold`
  request n - 1 (or 3.5 K) candidates per input, but the input's own id
  counts toward that cap (`query_candidates`, src/ssi.h:431 and :453; caller
  src/index_build.cpp:61), so one real candidate is lost per truncated query;
  a pair disappears when both of its queries drop each other, and which pair
  is dropped depends on `-p`. Repro: six inputs where four are identical
  (g00, g01, g02, g04) and g03, g05 differ, `--full -k 20 -S 512`: the square
  matrix gives g03 vs g04 = 0.371, but `--similarity-threshold 0.3 -p 1` and
  `--topk 5 -p 1` omit that pair from both rows (every other pair is listed);
  with `-p 2` a different pair is dropped.
* `binary-tail-32k` (medium): binary matrix output loses a whole block when a
  queued block holds a multiple of 32768 floats and is still queued when the
  computation ends: the final write loop writes `(nwritten & 32767)` floats
  for the last 128 KiB chunk, which is 0 (src/emitrect.cpp:395). Repro: 256
  inputs, `cmp --square --binary-output -p 128 --batch-size 128` writes
  131072 of 262144 bytes (exit 0); `-p 64` is complete. Realistic triggers:
  2048 genomes with `-p 16`, or 4096 with `-p 8`. The text path has the same
  pattern at src/emitrect.cpp:369 for formatted buffers that are an exact
  multiple of 128 KiB (not reproduced; it needs that exact size).
* `stacked-kmercounts-f32` (low): `sketch -N -o out` writes
  `out.kmercounts.f64` as float32 (`kmercounts_` is `std::vector<float>`,
  src/fastxsketch.h:45; written at src/sketch_core.cpp:200), while the
  per-input `.kmercounts.f64` files hold float64; a reader trusting the name
  gets garbage (the values are correct as float32).
* `full-mincount-counts` (low): `--full -N -m M` with M > 1 saves count 1 for
  every sampled k-mer (OPH, `-B` and `--prob` save the true counts):
  `CountFilteredCSetSketch::update` (bonsai hll/include/sketch/setsketch.h:1080
  and :1097) sets the count to 1 when an item reaches M and returns early for
  later occurrences.
* `bmh-stale-ids` (low): with `-B --save-kmers`, an input without k-mers
  (all records shorter than k, or nothing reaching `-m`) gets the sampled ids
  and counts of the previous input sketched by the same thread, because
  `bmh_t::reset()` (bonsai hll/include/sketch/bmh.h:404) does not clear
  `track_ids_`. `contain` then reports coverage for it: a database of a 3 kbp
  genome a and a 6 bp record b, queried with a, prints `100%:1` for b.
* `seq-o-file-mode` (low): `sketch -G -o out` without `--parse-by-seq`
  (where `-o` is required) writes the minimizer-sequence header over a file
  sized for sketch registers (src/fastxsketch.cpp:271, src/sketch_core.cpp:145):
  the body is zeros and `printmin` rejects it (8216 bytes, header implies
  372). The per-input `.mmerseq64` files are correct. Also in v2.1.20.
* `wsketch-tw-suffix` (cosmetic): 1-D `wsketch` ends `out.sampled.tw.txt`
  with a garbage character, because `';' + 'd' + ';' + 'L'` is integer
  arithmetic appended as one char (src/wsketch.cpp:365).

### Observations not treated as failures

* `--fastcmp N` fits its log-compression parameters to the registers of all
  inputs in a run, so the value for one pair changes with the other inputs.
* Sketch files do not record k, so `cmp --presketched --mash-distance` uses
  the default k unless `-k` is given again; cmp warns about it and the
  `--presketched` help says to pass the sketching `-k` (`presketched-mash-k`).
* `cmp -o` densifies one-permutation sketches in the stacked file, so it is
  not byte-identical to `sketch -o` (both reload to the same values).
* With fewer candidates than inputs (n - 1 > 3.5 K) top-k keeps the first
  candidates the index meets rather than the best (known caveat of the exact
  LSH change); recall of true top-K neighbours was about 0.65 to 0.75 in that
  regime in thorough runs, and 0.95 when every input is a candidate (which
  includes exact and compressed modes, whose low-similarity pairs may share
  no LSH key).
* `--phylip` prints the upper triangle with 9-character padded names, and
  `inf` for unrelated pairs.
* python/parse.py: `parse_binary_clustering` reads `fpath` (undefined) and an
  indptr of nclusters rather than nclusters + 1 entries; `parse_binary_kmers`
  indexes the function `alphabetcvt` with brackets. The suite decodes the
  layouts in its own code.

### Not covered

BigWig input (no BigWig writer in the Python standard library); LeafCutter
input; protein alphabets in these workflows; windowed minimizers are not
simulated (the `-G` stream is the reference for `contain -w`); exact oracles
for saved k-mers stop at k = 32 (64 with `-2`), beyond which dashing2 stores
rolling hashes.

### Runtime

Quick preset: 150 trials, 3 to 12 s with `-j 4`; thorough: 800 trials, about
45 s with `-j 8` (Apple Silicon, 12 cores). On stock v2.1.20 the quick preset
(seed 1) fails in 87 of 150 trials: wsketch row offsets hashed instead of ids,
`contain` with OPH and `-2` databases, BED labels and `sketch --bed -o`,
`--set`/`-J`/`--prob` crashes with `-m` or `-s`, `--presketched` labels,
containment direction, `--fastcmp-shorts` crashes, inverted `--greedy` for
distances and duplicate top-k entries under `-p`.

## Suite F: documentation audit

`suite_f_docs.py` (helpers in `d2rand_f.py`) checks that dashing2 does what
its usage messages (`dashing2`, `dashing2 <subcommand> -h/--help`) and
README.md say. Each trial runs one of 23 checks round robin on freshly drawn
inputs and parameters; quick runs every check twice (46 trials, about 2 s
with `-j 6`), thorough eight times (184 trials, about 8 s). `--check NAME`
(repeatable) runs only the named checks. KMC is not used.

| check | what it verifies |
|---|---|
| help | every subcommand prints its usage for `-h` and `--help`; `dist` is `cmp`; each documented short option (`-k -w -2 -m -F -Q -S -L -s -N -o -B -H -J -G -p -C -W -Z -P -c -f -v`) is documented and accepted by `sketch` and `cmp`; every long option named in the help exists in the parser (`src/options.h`, when present); the parser accepts no option the help never names, except the internal `--by-chrom`, `--exact-kmer-dist`, `--pairlist` and `--sig-ram-limit`; `-C` behaves alike in `sketch` and `cmp` |
| defaults | k 32 (64 with `-2`), S 1024, one-permutation sketch, canonical, symmetric output, equal to the explicit flags; the default k of every alphabet, with and without `-2`, equals the limit the help gives; the default output has the format comparison item 1 of the help names, and `--phylip` prints PHYLIP |
| sketch_types | aliases `--full/--full-setsketch`, `-B/--multiset/--bagminhash`, `-P/--prob/--pminhash`, `--set/-H`, `-J/--countdict` give identical values and the right header label; `--set` and `-J` are exact |
| measures | every documented measure and alias (`--intersection(-size)`, `--union-size`, `--containment` as I/|row|, `--symmetric-containment`, `--mash-distance/--distance/--poisson-distance`) against the exact oracle |
| layouts | `--square/--asymmetric-all-pairs/--asymmetric`; `-F X -Q X` equals `--square`; `-Q` gives |F| x |Q|; `-F` adds to positional arguments; a `-F` line with several files is one sketch |
| sparse | `--topk/--top-k` lists the true top neighbours with the matrix values, sorted; `--topk` of N-1 or more gives what the `--topk` help line says (every other item that shares a k-mer); `--similarity-threshold` lists exactly the entries above the threshold (inputs are families of related genomes so that LSH recall is not at issue) |
| greedy | `--greedy t` and `tE` recover well-separated families; `F` needs `--parse-by-seq` and then prints input records |
| binary | `--binary-output/--emit-binary/--binary` layouts for symmetric, square and panel output equal the text values; `--cmpout/--distout/--cmp-outfile`, `--cmpout -`; binary top-k, threshold (CSR) and greedy output in the layout the help describes, holding the text output's values |
| stacked | `sketch -o` file layout and names file; `cmp --presketched` equals the direct run (OPH, `--full`, `-B`, `--prob`); for Mash distances the `--presketched` help says to pass the sketching `-k`, cmp warns without it, and with it the value equals the direct run |
| save_kmers | `-s/-N` output files; the names path, header size and header fields (sketch size, k, seed) the help gives; which sketch types write per-input k-mer files without `-o`, as the help says |
| contain | self-coverage 100%, `-o`, `-b` layout, `-F`, `-p`, stdin default |
| seq | `--seq --parse-by-seq -o` plus `printmin` (tabular, `-f`, `-o`) equals the k-mer sequence oracle; `--hp-compress`; the header fields the help lists (including the canonical bit) |
| filterset, window_downsample, spacing, protein | exact oracles for `--filterset`, `-w` (default k, minimizer density), `--downsample` (binomial), `--spacing` (gap and run-length forms; disables canonicalization), `--protein/--protein20/--enable-protein` |
| countmin, cache, sizes, seqs_in_ram | `-c/--countmin-size` scope and header; `--cache/--cache-sketches` with `--outprefix/--prefix`; `-S`, `-L` bounds, `--fastcmp/--regsize/--regbytes`, `--fastcmp-bytes/-shorts/-words` equal their `--setsketch-ab` values and also select `--full` without it; `--seqs-in-ram` |
| readme, edit_distance | README Uses 1, 2, 3, 4 and 6; the Use 7 commands, read from README.md and run, write the neighbour table to the file they name; the README's canonicalization claim |
| wsketch | total weights for every weight type (`-f`, `-H`, `-U`), CSR row sums (`-P`, `-` for uniform weights), `-u` ids, `tw.txt` format; the help's examples, read from `wsketch -h` and run with files of the roles their names suggest (k-mer files prepared as the help says), give the right weights and the same sketch as the usage-line order |

Checks of documented values read the text from the binary's usage messages
and from README.md and src/options.h of the source tree, so they fail when
text and behaviour disagree, whichever of the two is wrong. Each failing check
prints the dashing2 command that shows it. With this suite, the binary and
source tree of `test/all-fixes` (before the documentation fixes) fail every
id below that `test/all-fixes-7` fixes, and `test/all-fixes-7` passes.

The issues this suite found, each described as it was before its fix. All are
fixed (branches under "Fixed issues"); a recurrence is reported as a
regression of its id, and `d2rand_f.DOC_BUGS` remains for registering new
ones:

* `help-k-limits`: help says the direct k limit is 31 for DNA and 22 for `--protein8`; it is 32 and 21 (and `--protein14`, 16, is not listed).
* `default-not-phylip`: help says the default is "Upper Triangular PHYLIP"; the default is the tab-separated symmetric matrix, PHYLIP needs the undocumented `--phylip`.
* `header-sketchtype-label`: `cmp --multiset/--bagminhash/--prob/--pminhash` print `sketchtype:onepermsetsketch` in the header (values are right; `-B`/`-P` print the right label); `-J` prints `fullyunknown`.
* `topk-large-k`: `--topk` above N-1 prints neighbour lists, not the pairwise matrix the help promises.
* `sparse-binary-format`, `greedy-binary-ids`: binary top-k/threshold output is CSR with a 16-byte header, not the documented matrix; greedy member ids are 32-bit, not 64-bit.
* `presketched-mash-k`: `cmp --presketched --mash-distance` uses the default k (32) unless `-k` is repeated, because stacked files do not store k.
* `save-kmers-names-path`, `save-kmers-header-size`, `save-kmers-no-file`: `-s -o X` writes names to `X.kmer64.names.txt` (help: `X.kmer.names.txt`) with a 24-byte header (help: 16); without `-o`, `-s/-N` write nothing for the default sketch and `--full`.
* `seq-header`: `--seq` files have a 20-byte header (an extra alphabet word), not the documented 16.
* `countmin-exact-modes`, `countsketch-header-newline`: `--countmin-size` also turns `--set`/`-J` into count-sketch bucket comparisons (help: only BagMinHash/ProbMinHash); it also puts a newline into the `#Dashing2Options` header.
* `fastcmp-presets-need-full`: `--fastcmp-bytes/-shorts/-words` and `--setsketch-ab` abort unless `--full` is given.
* `seqs-in-ram-ignored`: `--seqs-in-ram` has no effect (`static bool` in a header).
* `readme-use7-no-output`, `readme-edit-distance-crash`, `edit-distance-nondeterministic`, `readme-canon-default`: README Use 7 writes binary sketches to its `.tsv` and no table; Use 6 segfaults; `--edit-distance` results change between runs; README says canonicalization is off by default (it is on).
* `wsketch-U-1d`, `wsketch-tw-garbage`, `wsketch-example-order`, `wsketch-example-kmerset-header`: `-U` weights read as float64 with one or two paths; `tw.txt` ends in a garbage character; help examples put the weights before the ids (the `-` example aborts); `.kmerset64` files used in the examples carry an 8-byte header that is read as an id.
* `short-opts-differ`, `undocumented-options`: `-C` works in `sketch` but not `cmp`; 26 long options accepted by the parser are never named in the help.

On stock v2.1.20 the suite additionally fails exact `--union-size` (55 vs 2270), `--protein --set` (cardinality 0), `cmp --presketched` (aborts), `--full --fastcmp-bytes` (aborts) and `--filterset` with `-2` (no effect), all fixed in the tested build.
