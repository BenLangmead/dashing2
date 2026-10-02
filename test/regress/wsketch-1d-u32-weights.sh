#!/usr/bin/env bash
# `wsketch -U` reads 32-bit integer weights (see wsketch --help). A two-path
# (one-row) input with -U must give the same total weight and the same sketch
# as the same weights stored as float64 (the total is the sum of the weights,
# or the number of ids for -q, which ignores weights). Before the fix the
# one-row path ignored -U and read the 32-bit weights as float64.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
python3 - <<'PY'
import random, struct
r = random.Random(11)
n = 200
ids = r.sample(range(1, 1 << 40), n)
w = [r.randint(1, 50) for _ in range(n)]
open('ids.u64', 'wb').write(struct.pack('<%dQ' % n, *ids))
open('w.f64', 'wb').write(struct.pack('<%dd' % n, *w))
open('w.u32', 'wb').write(struct.pack('<%dI' % n, *w))
open('sum.txt', 'w').write('%d\n' % sum(w))
PY
for type in "" -B -q; do
    "$D2" wsketch -S 64 $type -o f64 ids.u64 w.f64 >/dev/null 2>&1
    "$D2" wsketch -S 64 $type -U -o u32 ids.u64 w.u32 >/dev/null 2>&1
    want=$(cat sum.txt); [ "$type" = -q ] && want=200
    got=$(sed -n 's/^Total weight: \([0-9]*\)\..*/\1/p' u32.sampled.tw.txt 2>/dev/null)
    if [ "$got" = "$want" ]; then echo "PASS -U ${type:-default} total weight: got $got expected $want"
    else echo "FAIL -U ${type:-default} total weight: got $got expected $want"; fail=1; fi
    same=yes
    for f in hashes.f64 ids.u64; do cmp -s f64.sampled.$f u32.sampled.$f || same=no; done
    if [ $same = yes ]; then echo "PASS -U ${type:-default} sketch: identical to float64 weights"
    else echo "FAIL -U ${type:-default} sketch: differs from float64 weights"; fail=1; fi
    rm -f f64.sampled.* u32.sampled.*
done
exit $fail
