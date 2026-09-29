#!/usr/bin/env bash
# --filterset F removes F's k-mers before sketching. With -2 (128-bit k-mers)
# it must remove exactly as many as with 64-bit k-mers: the --set union of a
# file with itself is the number of its distinct k-mers not in F.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
# x: 50 kbp random; half: its first 25 kbp, used as the filter.
python3 - <<'PY'
import random
r = random.Random(6)
x = ''.join(r.choice('ACGT') for _ in range(50000))
open('x.fa', 'w').write(f'>x\n{x}\n')
open('half.fa', 'w').write(f'>half\n{x[:25000]}\n')
PY
# Exact count of distinct canonical k-mers of x that are absent from half.
truth() {
    python3 - "$1" <<'PY'
import sys
k = int(sys.argv[1]); comp = str.maketrans('ACGT', 'TGCA')
def kset(fn):
    s = ''.join(l.strip() for l in open(fn) if not l.startswith('>'))
    return {min(s[i:i + k], s[i:i + k].translate(comp)[::-1]) for i in range(len(s) - k + 1)}
print(len(kset('x.fa') - kset('half.fa')))
PY
}
# -k 21: 64-bit control; -2 -k 50: 128-bit encoded k-mers; -2 -k 70: 128-bit rolling hashes.
for opts in "-k 21" "-2 -k 50" "-2 -k 70"; do
    k=${opts##* }
    got=$("$D2" cmp $opts --set --union-size --filterset half.fa x.fa x.fa 2>/dev/null | grep -v '^#' | awk 'NR == 1 {print $3}')
    rm -f ./*.kmer*
    exp=$(truth "$k")
    if [ "$got" = "$exp" ]; then echo "PASS $opts --filterset union: got $got expected $exp"
    else echo "FAIL $opts --filterset union: got $got expected $exp"; fail=1; fi
done
exit $fail
