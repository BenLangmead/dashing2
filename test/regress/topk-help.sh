#!/usr/bin/env bash
# The --topk help line must say that a value of N - 1 or more lists every other
# item except those sharing no sampled k-mers (no LSH key), since such items are
# never LSH candidates. Fixtures: a random 3 kbp genome a, a copy b with 1%
# substitutions, and an unrelated random genome c. With --topk 2 (N - 1), a
# lists b, and c, which shares no k-mers with either, lists nobody.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
r = random.Random(5)
a = ''.join(r.choice('ACGT') for _ in range(3000))
b = list(a)
for p in r.sample(range(len(b)), 30):
    b[p] = r.choice('ACGT'.replace(b[p], ''))
c = ''.join(r.choice('ACGT') for _ in range(3000))
for n, s in (('a', a), ('b', ''.join(b)), ('c', c)):
    open(n + '.fa', 'w').write('>' + n + '\n' + s + '\n')
PY
fail=0
check() {  # check <name> <got> <expected>
    if [ "$2" = "$3" ]; then echo "PASS $1: got $2 expected $3"
    else echo "FAIL $1: got $2 expected $3"; fail=1; fi
}
line=$("$D" cmp --help 2>&1 | grep '^--topk/--top-k')
check "help line names the exception" "$(grep -c 'share no sampled k-mers' <<<"$line")" 1
"$D" cmp -k 21 --full -S 512 --topk 2 a.fa b.fa c.fa > tk.txt 2>/dev/null
# Each row: the item's name, then one tab-separated name:similarity per listed neighbor
listed() { awk -F'\t' -v q="$1" '!/^#/ && $1 ~ q { n = 0; for (i = 2; i <= NF; ++i) if ($i != "") { sub(/:[^:]*$/, "", $i); printf "%s%s", (n++ ? "," : ""), $i } }' tk.txt; }
check "a lists its relative" "$(listed '^a.fa')" b.fa
check "unrelated c lists nobody" "$(listed '^c.fa')" ""
exit $fail
