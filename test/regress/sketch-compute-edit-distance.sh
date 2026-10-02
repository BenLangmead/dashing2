#!/usr/bin/env bash
# README Use 6: `dashing2 sketch --parse-by-seq --edit-distance
# --compute-edit-distance --cmpout` must exit 0 and report the exact edit
# distance between every pair of records. The expected values come from a
# plain Levenshtein dynamic program in Python.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
python3 - <<'PY'
import random
r = random.Random(11)
base = ''.join(r.choice('ACGT') for _ in range(1200))
def mutate(s, p):
    out = []
    for c in s:
        x = r.random()
        if x < p / 3: continue                                  # deletion
        elif x < 2 * p / 3: out.append(r.choice('ACGT'))        # substitution
        elif x < p: out += [c, r.choice('ACGT')]                # insertion
        else: out.append(c)
    return ''.join(out)
recs = [mutate(base, p) for p in (0.02, 0.05, 0.2)]
with open('input.fasta', 'w') as f:
    for i, s in enumerate(recs):
        f.write(f'>r{i}\n{s}\n')
def lev(a, b):
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb)))
        prev = cur
    return prev[-1]
with open('expected.txt', 'w') as f:
    for i in range(3):
        for j in range(i + 1, 3):
            f.write(f'r{i} r{j} {lev(recs[i], recs[j])}\n')
PY
check() { # description outfile args...
    local name=$1 out=$2; shift 2
    "$D2" sketch "$@" input.fasta >stdout.txt 2>err.txt; local rc=$?
    [ "$out" = - ] && out=stdout.txt
    local got
    got=$(awk -F'\t' '/^#Sources/ {for (i = 2; i <= NF; i++) n[i] = $i; next}
                      /^#/ {next}
                      {sub(/ +$/, "", $1); for (i = 2; i <= NF; i++) if ($i != "-") print $1, n[i], $i + 0}' "$out" 2>/dev/null)
    if [ $rc -eq 0 ] && [ "$got" = "$(cat expected.txt)" ]; then
        echo "PASS $name: got $(echo $got) expected $(echo $(cat expected.txt))"
    else
        echo "FAIL $name: status $rc, got '$(echo $got)' expected $(echo $(cat expected.txt))"; fail=1
    fi
}
check "README Use 6 (-p8, --cmpout file)" knn.edit-distance.tbl -p8 --cmpout knn.edit-distance.tbl -k7 --parse-by-seq --edit-distance --compute-edit-distance
check "-p1, --cmpout -" - -p1 --cmpout - -k7 --parse-by-seq --edit-distance --compute-edit-distance
check "--seqs-in-ram" - -p4 --cmpout - -k7 --parse-by-seq --edit-distance --compute-edit-distance --seqs-in-ram
exit $fail
