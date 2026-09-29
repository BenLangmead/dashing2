#!/usr/bin/env bash
# wsketch must sketch the ids it is given, not their offsets within a row or
# file. Rows A = {0..1999} and B = {1000..2999} with unit weights have Jaccard
# 1000/3000 = 1/3 (also the probability Jaccard for ProbMinHash and the weighted
# Jaccard for BagMinHash), so the fraction of equal registers must be near 1/3
# (1024 registers, standard error about 0.015; accept 0.25..0.42). Before the
# fix both rows hashed offsets 0..1999 and all registers were equal. Checked
# for CSR input (-q, -B, default ProbMinHash) and for two 1-D inputs.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import struct
a = list(range(0, 2000)); b = list(range(1000, 3000))
open('ind.u64', 'wb').write(struct.pack('<4000Q', *(a + b)))
open('indptr.u64', 'wb').write(struct.pack('<3Q', 0, 2000, 4000))
open('w.f64', 'wb').write(struct.pack('<4000d', *([1.0] * 4000)))
open('a.u64', 'wb').write(struct.pack('<2000Q', *a))
open('b.u64', 'wb').write(struct.pack('<2000Q', *b))
PY
fail=0
check() { # name got
    if awk -v g="$2" 'BEGIN{exit !(g != "" && g >= 0.25 && g <= 0.42)}'; then
        echo "PASS $1: got $2 expected 0.25..0.42"
    else
        echo "FAIL $1: got $2 expected 0.25..0.42"; fail=1
    fi
}
# CSR: the regs file has a 16-byte header, n cardinalities, then n x m registers.
csrfrac() { python3 -c "
import struct, glob
d = open(glob.glob('o.sampled.regs.stacked.*')[0], 'rb').read()
n, m = struct.unpack('<2Q', d[:16]); v = struct.unpack('<%dd' % (n * m), d[16 + 8 * n:16 + 8 * n + 8 * n * m])
print('%.3f' % (sum(x == y for x, y in zip(v[:m], v[m:2 * m])) / m))"; }
for f in -q -B ""; do
    "$D" wsketch $f -S 1024 -o o ind.u64 w.f64 indptr.u64 >/dev/null 2>&1
    check "CSR wsketch '$f'" "$(csrfrac)"
    rm -f o.*
done
# 1-D: each hashes file holds the total weight, then m registers.
for f in -q -B ""; do
    "$D" wsketch $f -S 1024 -o oa a.u64 >/dev/null 2>&1
    "$D" wsketch $f -S 1024 -o ob b.u64 >/dev/null 2>&1
    got=$(python3 -c "
import struct
x = open('oa.sampled.hashes.f64', 'rb').read()[8:]; y = open('ob.sampled.hashes.f64', 'rb').read()[8:]
m = len(x) // 8; x = struct.unpack('<%dd' % m, x); y = struct.unpack('<%dd' % m, y)
print('%.3f' % (sum(p == q for p, q in zip(x, y)) / m))")
    check "1-D wsketch '$f'" "$got"
    rm -f oa.* ob.*
done
exit $fail
