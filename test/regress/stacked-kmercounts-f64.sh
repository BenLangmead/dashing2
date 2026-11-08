#!/usr/bin/env bash
# `sketch -N -o out` writes the counts of the sampled k-mers to
# out.kmercounts.f64, one float64 per sketch register, like the per-input
# .kmercounts.f64 files. Each input repeats one random 1.5 kbp unit as
# separate records (a: 3 copies, b: 5), so every k-mer of a has count 3 and
# every k-mer of b count 5. Before the fix the stacked file held float32 values
# (half the expected size), which read as float64 are meaningless.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
python3 - <<'PY'
import random
r = random.Random(4)
for name, reps in (('a', 3), ('b', 5)):
    u = ''.join(r.choice('ACGT') for _ in range(1500))
    with open(name + '.fa', 'w') as f:
        f.writelines(f'>{name}{i}\n{u}\n' for i in range(reps))
PY
S=32
for mode in default --full; do
    flag=$mode; [ "$mode" = default ] && flag=
    rm -f out*
    "$D2" sketch -k 21 -S $S -N $flag -o out a.fa b.fa >/dev/null 2>&1
    got=$(python3 - $S <<'PY'
import os, struct, sys
S = int(sys.argv[1])
try:
    d = open('out.kmercounts.f64', 'rb').read()
except OSError:
    print('missing'); sys.exit()
if len(d) != 2 * S * 8:
    print('%d bytes' % len(d)); sys.exit()
v = struct.unpack('<%dd' % (2 * S), d)
print('a=%s b=%s' % (sorted(set(v[:S])), sorted(set(v[S:]))))
PY
)
    expected="a=[3.0] b=[5.0]"
    if [ "$got" = "$expected" ]; then echo "PASS $mode stacked counts: got $got expected $expected"
    else echo "FAIL $mode stacked counts: got $got expected $expected ($((2 * S * 8)) bytes)"; fail=1; fi
done
rm -f ./*.kmer* ./*.kmercounts*
exit $fail
