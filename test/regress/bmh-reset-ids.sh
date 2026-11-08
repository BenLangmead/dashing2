#!/usr/bin/env bash
# With -B --save-kmers, an input without k-mers must not get the sampled ids
# and counts of the input sketched before it by the same thread. Inputs (one
# thread; sketched largest first): a repeats a unit twice, c holds one copy of
# another unit (no k-mer reaches -m 2) and b is a 6 bp record (no 15-mer). The
# saved ids and counts of b and c must all be zero, and a's must not be.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
python3 - <<'PY'
import random
r = random.Random(3)
u, v = (''.join(r.choice('ACGT') for _ in range(n)) for n in (3000, 2000))
open('a.fa', 'w').write('>a1\n%s\n>a2\n%s\n' % (u, u))
open('c.fa', 'w').write('>c\n%s\n' % v)
open('b.fa', 'w').write('>b\nACGTAC\n')
PY
"$D2" sketch -p 1 -k 15 -S 16 -m 2 -B -N --cache a.fa c.fa b.fa >/dev/null 2>&1
for f in a c b; do
    got=$(python3 - $f <<'PY'
import glob, struct, sys
f = sys.argv[1]
ids = glob.glob(f + '.fa*.kmer.u64'); cts = glob.glob(f + '.fa*.kmercounts.f64')
if not ids or not cts:
    print('missing'); sys.exit()
i = open(ids[0], 'rb').read(); c = open(cts[0], 'rb').read()
i = struct.unpack('<%dQ' % (len(i) // 8), i); c = struct.unpack('<%dd' % (len(c) // 8), c)
print('%d nonzero ids, %d nonzero counts' % (sum(x != 0 for x in i), sum(x != 0 for x in c)))
PY
)
    if [ $f = a ]; then expected="16 nonzero ids, 16 nonzero counts"; else expected="0 nonzero ids, 0 nonzero counts"; fi
    if [ "$got" = "$expected" ]; then echo "PASS $f.fa: got $got expected $expected"
    else echo "FAIL $f.fa: got $got expected $expected"; fail=1; fi
done
rm -f ./*.bmh ./*.kmer* ./*.kmercounts*
exit $fail
