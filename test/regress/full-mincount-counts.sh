#!/usr/bin/env bash
# `sketch --full -N -m M -o out` saves the count of every sampled k-mer, as the
# default sketch, -B and --prob do. Each input holds one random unit repeated
# as separate records (a: 3 copies, b: 4) and one unit present once, so with
# -m 2 every sampled k-mer of a has count 3 and every one of b count 4. Before
# the fix --full saved count 1 for every sampled k-mer when M > 1.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
python3 - <<'PY'
import random
r = random.Random(9)
for name, reps in (('a', 3), ('b', 4)):
    u, once = (''.join(r.choice('ACGT') for _ in range(1500)) for _ in range(2))
    with open(name + '.fa', 'w') as f:
        f.writelines(f'>{name}{i}\n{u}\n' for i in range(reps))
        f.write(f'>{name}once\n{once}\n')
PY
S=32
for mode in --full default; do
    flag=$mode; [ "$mode" = default ] && flag=
    rm -f out out.*
    "$D2" sketch -k 21 -S $S -m 2 -N $flag -o out a.fa b.fa >/dev/null 2>&1
    # The stacked counts are one value per register; accept float32 or float64 files.
    got=$(python3 - $S <<'PY'
import struct, sys
S = int(sys.argv[1])
try:
    d = open('out.kmercounts.f64', 'rb').read()
except OSError:
    print('missing'); sys.exit()
fmt = 'd' if len(d) == 2 * S * 8 else 'f' if len(d) == 2 * S * 4 else None
if fmt is None:
    print('%d bytes' % len(d)); sys.exit()
v = struct.unpack('<%d%s' % (2 * S, fmt), d)
print('a=%s b=%s' % (sorted(set(v[:S])), sorted(set(v[S:]))))
PY
)
    expected="a=[3.0] b=[4.0]"
    if [ "$got" = "$expected" ]; then echo "PASS $mode -m 2 saved counts: got $got expected $expected"
    else echo "FAIL $mode -m 2 saved counts: got $got expected $expected"; fail=1; fi
done
rm -f ./*.kmer* ./*.kmercounts*
exit $fail
