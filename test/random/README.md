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
  and `--prob`);
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
* bias: |mean z| <= 3/sqrt(trials) + allowance, allowance 0.1, raised to
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

### Known issues (expected failures)

Run with `--strict` to count these as failures.

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
  the default k unless `-k` is given again (part of the known issue that sketch
  files carry no parameters).
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
