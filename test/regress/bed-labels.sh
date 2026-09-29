#!/usr/bin/env bash
# BED rows must be labelled with the file whose values they hold. Exact
# interval Jaccard: a~b = 500/1500 (overlap 500 bp of a 1500 bp union); c is on
# another chromosome, so a~c = b~c = 0. c.bed is the largest file, so it is
# sketched first; a permutation of labels shows up as c~b = 1/3.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
printf 'chr1\t0\t1000\n' > a.bed
printf 'chr1\t500\t1500\n' > b.bed
printf 'chr2\t0\t1000\nchr2\t2000\t2500\nchr2\t4000\t4100\n' > c.bed
"$D" cmp --bed --full --square a.bed b.bed c.bed 2>/dev/null > out.tsv
python3 - <<'PY'
import sys
rows = [l.rstrip('\n').split('\t') for l in open('out.tsv') if l.strip()]
cols = next(r for r in rows if r[0] == '#Sources')[1:]
m = {r[0].strip(): dict(zip(cols, r[1:])) for r in rows if not r[0].startswith('#')}
fail = 0
for x, y, want in [('a.bed', 'b.bed', 500 / 1500), ('a.bed', 'c.bed', 0.0), ('b.bed', 'c.bed', 0.0)]:
    got = m.get(x, {}).get(y)
    ok = got is not None and abs(float(got) - want) < 0.03
    print(('PASS' if ok else 'FAIL') + f' {x}~{y}: got {got} expected {want:.4f}')
    fail |= not ok
sys.exit(fail)
PY
