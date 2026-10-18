#!/usr/bin/env bash
# --countdict must pair each file's k-mers with that file's own counts. Files
# are sketched in decreasing size order, so the inputs below (listed smallest
# first) exercise the mapping between processing order and input order.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
check() { # name got expected
    if python3 -c "import sys; sys.exit(abs(float('$2') - float('$3')) > 1e-4)"; then
        echo "PASS $1: got $2 expected $3"
    else
        echo "FAIL $1: got $2 expected $3"; fail=1
    fi
}
# a: sequence s once; b: s three times; c: unrelated u six times; d: s[:3000] twice.
python3 - <<'PY'
import random
r = random.Random(1)
s = ''.join(r.choice('ACGT') for _ in range(5000))
u = ''.join(r.choice('ACGT') for _ in range(5000))
for name, recs in [('a', [s]), ('b', [s] * 3), ('c', [u] * 6), ('d', [s[:3000]] * 2)]:
    with open(name + '.fa', 'w') as f:
        f.writelines(f'>{name}{i}\n{x}\n' for i, x in enumerate(recs))
PY
# Exact weighted Jaccard (sum of min counts / sum of max counts, canonical 21-mers).
read -r exp_ad exp_bd <<<"$(python3 - <<'PY'
from collections import Counter
comp = str.maketrans('ACGT', 'TGCA')
def counts(fn):
    c = Counter()
    for line in open(fn):
        if line.startswith('>'): continue
        x = line.strip()
        for i in range(len(x) - 20):
            k = x[i:i + 21]; c[min(k, k.translate(comp)[::-1])] += 1
    return c
def wj(a, b):
    ks = set(a) | set(b)
    return sum(min(a[k], b[k]) for k in ks) / sum(max(a[k], b[k]) for k in ks)
a, b, d = counts('a.fa'), counts('b.fa'), counts('d.fa')
print(wj(a, d), wj(b, d))
PY
)"
for p in 1 4; do
    "$D2" cmp -p $p -k 21 --countdict --square a.fa b.fa c.fa d.fa >out.txt 2>/dev/null
    rm -f ./*.kmer*
    # Rows follow input order: a, b, c, d; column 1 is the name.
    get() { grep -v '^#' out.txt | awk -v r="$1" -v c="$2" 'NR == r {print $(c + 1)}'; }
    check "p$p self(c)" "$(get 3 3)" 1
    check "p$p self(d)" "$(get 4 4)" 1
    check "p$p a~d" "$(get 1 4)" "$exp_ad"
    check "p$p b~d" "$(get 2 4)" "$exp_bd"
done
# The count-file column of <out>.names.txt must name each row's own file.
"$D2" sketch -k 21 --countdict -o sk.out a.fa b.fa c.fa d.fa >/dev/null 2>&1
rm -f ./*.kmer*
bad=$(awk -F'\t' '!/^#/ {n = $1; sub(/\.fa$/, "", n); if (index($3, n ".fa") == 0) c++} END {print c + 0}' sk.out.names.txt)
check "names.txt count files match rows (mismatches)" "$bad" 0
exit $fail
