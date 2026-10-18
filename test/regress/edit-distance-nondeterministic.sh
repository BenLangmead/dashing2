#!/usr/bin/env bash
# `sketch --parse-by-seq --edit-distance` (OrderMinHash) must give the same
# similarities on every run and at every thread count, and the similarity of
# two records must not depend on the order of the records in the file. Before
# the fix, sketches of records longer than 128 characters read memory before
# the sequence, so the values changed from run to run.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
python3 - <<'PY'
import random
r = random.Random(5)
base = ''.join(r.choice('ACGT') for _ in range(2000))
def mutate(s, p):
    out = []
    for c in s:
        x = r.random()
        if x < p / 3: continue
        elif x < 2 * p / 3: out.append(r.choice('ACGT'))
        elif x < p: out += [c, r.choice('ACGT')]
        else: out.append(c)
    return ''.join(out)
recs = [mutate(base, p) for p in (0.02, 0.05, 0.2)]
with open('input.fasta', 'w') as f:
    for i in (0, 1, 2): f.write(f'>r{i}\n{recs[i]}\n')
with open('reversed.fasta', 'w') as f:
    for i in (2, 1, 0): f.write(f'>r{i}\n{recs[i]}\n')
PY
pairs() { # file threads: prints "ri rj similarity" for each pair, sorted by name
    "$D2" sketch -p"$2" -k7 --parse-by-seq --edit-distance --cmpout - "$1" 2>/dev/null |
        awk -F'\t' '/^#Sources/ {for (i = 2; i <= NF; i++) n[i] = $i; next}
                    /^#/ {next}
                    {sub(/ +$/, "", $1); for (i = 2; i <= NF; i++) if ($i != "-") {
                        a = $1; b = n[i]; if (b < a) {t = a; a = b; b = t}; print a, b, $i}}' | sort
}
ref=$(pairs input.fasta 1)
for run in "input.fasta 1" "input.fasta 1" "input.fasta 2" "input.fasta 4" "input.fasta 8" \
           "reversed.fasta 1" "reversed.fasta 4"; do
    got=$(pairs $run)
    if [ -n "$ref" ] && [ "$got" = "$ref" ]; then echo "PASS $run threads: got $(echo $got) expected the first -p1 run"
    else echo "FAIL $run threads: got $(echo $got) expected $(echo $ref)"; fail=1; fi
done
exit $fail
