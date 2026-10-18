#!/usr/bin/env bash
# Stacked and per-input sketch files do not record k, and the Mash distance depends
# on it: cmp --presketched --mash-distance uses the default k (32) unless -k is
# given again. With -k the result must equal the direct run, and without it cmp
# must warn on stderr.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
r = random.Random(1)
a = ''.join(r.choice('ACGT') for _ in range(5000))
b = list(a)
for i in range(0, len(b), 97):
    b[i] = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[b[i]]
open('a.fa', 'w').write('>a\n' + a + '\n')
open('b.fa', 'w').write('>b\n' + ''.join(b) + '\n')
PY
fail=0
check() {  # check <name> <got> <expected>
    if [ "$2" = "$3" ]; then echo "PASS $1: got $2 expected $3"
    else echo "FAIL $1: got ${2:-nothing} expected $3"; fail=1; fi
}
cell() { awk 'NR==4{print $3}'; }
"$D" sketch --full -k 15 -o stacked a.fa b.fa >/dev/null 2>&1
direct=$("$D" cmp --full -k 15 --mash-distance a.fa b.fa 2>/dev/null | cell)
check "presketched -k 15 equals direct k=15" "$("$D" cmp --presketched -k 15 --mash-distance stacked 2>/dev/null | cell)" "$direct"
"$D" cmp --presketched --mash-distance stacked >/dev/null 2>err.txt
check "warning without -k" "$(grep -c 'do not record k' err.txt)" 1
"$D" cmp --presketched -k 15 --mash-distance stacked >/dev/null 2>err.txt
check "no warning with -k" "$(grep -c 'do not record k' err.txt)" 0
"$D" cmp --presketched stacked >/dev/null 2>err.txt
check "no warning for similarity" "$(grep -c 'do not record k' err.txt)" 0
exit $fail
