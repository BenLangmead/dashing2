#!/usr/bin/env bash
# contain with a database built by the default one-permutation sketch (OPH).
# s.fa is a 2 kbp random sequence (about 1,980 k-mers), so with 1024 or 4096
# buckets many buckets stay empty; empty buckets hold no k-mer and must not
# count against coverage. Invariants: s.fa against itself is 100%; h.fa (the
# first half of s.fa) covers about half of the sampled k-mers (accept 40..60%).
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
r = random.Random(3)
s = ''.join(r.choice('ACGT') for _ in range(2000))
y = ''.join(r.choice('ACGT') for _ in range(20000))
open('s.fa', 'w').write('>s\n' + s + '\n')
open('h.fa', 'w').write('>h\n' + s[:1000] + '\n')
open('y.fa', 'w').write('>y\n' + y + '\n')
PY
fail=0
check() { # name got lo hi
    if awk -v g="$2" -v l="$3" -v h="$4" 'BEGIN{exit !(g != "" && g >= l && g <= h)}'; then
        echo "PASS $1: got $2% expected $3..$4%"
    else
        echo "FAIL $1: got $2% expected $3..$4%"; fail=1
    fi
}
for S in 1024 4096; do
    "$D" sketch -k 21 -S $S --save-kmers -o db s.fa y.fa >/dev/null 2>&1
    out=$("$D" contain db.kmer64 s.fa h.fa 2>/dev/null)
    cov() { echo "$out" | awk -v q="$1" '$1==q{split($2,a,"%"); print a[1]}'; }
    check "OPH -S $S s.fa in s" "$(cov s.fa)" 100 100
    check "OPH -S $S h.fa in s" "$(cov h.fa)" 40 60
    rm -f db db.*
done
exit $fail
