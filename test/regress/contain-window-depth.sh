#!/usr/bin/env bash
# contain reports <coverage%:mean depth>. Depth is the mean number of times the
# query contains each matched reference k-mer. Fixtures: x.fa (random 20 kbp),
# x2.fa (two copies of x as two records), x25.fa (two copies plus the first
# half of x). Expected depths: 1, 2, and 2.5 (half of x's sampled k-mers occur
# three times, the rest twice; accept 2.4..2.6). Checked without windows and
# with minimizer windows -w 30 and -w 50, where each window's minimizer must
# not be counted once per window.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
r = random.Random(1)
x = ''.join(r.choice('ACGT') for _ in range(20000))
y = ''.join(r.choice('ACGT') for _ in range(20000))
open('x.fa', 'w').write('>x\n' + x + '\n')
open('y.fa', 'w').write('>y\n' + y + '\n')
open('x2.fa', 'w').write('>a\n' + x + '\n>b\n' + x + '\n')
open('x25.fa', 'w').write('>a\n' + x + '\n>b\n' + x + '\n>c\n' + x[:10000] + '\n')
PY
fail=0
check() { # name got lo hi
    if awk -v g="$2" -v l="$3" -v h="$4" 'BEGIN{exit !(g != "" && g >= l && g <= h)}'; then
        echo "PASS $1: got depth $2 expected $3..$4"
    else
        echo "FAIL $1: got depth $2 expected $3..$4"; fail=1
    fi
}
for w in "" "-w 30" "-w 50"; do
    "$D" sketch --full -k 21 $w --save-kmers -o db x.fa y.fa >/dev/null 2>&1
    for p in 1 2; do
        out=$("$D" contain -p $p db.kmer64 x.fa x2.fa x25.fa 2>/dev/null)
        depth() { echo "$out" | awk -v q="$1" '$1==q{split($2,a,":"); print a[2]}'; }
        check "k=21 $w -p $p x.fa" "$(depth x.fa)" 1 1
        check "k=21 $w -p $p x2.fa" "$(depth x2.fa)" 2 2
        check "k=21 $w -p $p x25.fa" "$(depth x25.fa)" 2.4 2.6
    done
    rm -f db db.*
done
exit $fail
