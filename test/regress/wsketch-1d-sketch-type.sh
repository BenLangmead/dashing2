#!/usr/bin/env bash
# wsketch's 1-D path (ids [weights]) and CSR path (ids weights|- indptr) must use
# the same sketch for the same flag: ProbMinHash by default, BagMinHash with -B,
# SetSketch with -q. A CSR matrix with a single row [0, n) holds the same data
# as the 1-D input, so both paths must emit identical registers.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
python3 - <<'PY'
import random, struct
r = random.Random(4)
ids = sorted(r.sample(range(1, 10**9), 300))
w = [r.random() * 5 + 0.1 for _ in ids]
open('ids.u64', 'wb').write(struct.pack(f'<{len(ids)}Q', *ids))
open('w.f64', 'wb').write(struct.pack(f'<{len(w)}d', *w))
open('ip.u64', 'wb').write(struct.pack('<2Q', 0, len(ids)))
PY
# 1-D output: [double total weight][S registers]; CSR output: [u64 n][u64 S][n doubles][n x S registers].
same_regs() {
    python3 - "$1" "$2" <<'PY'
import sys
S = 64
one = open(sys.argv[1], 'rb').read()[8:8 + 8 * S]
csr = open(sys.argv[2], 'rb').read()[24:24 + 8 * S]
print('identical' if one == csr else 'different')
PY
}
for f in "" "-B" "-q"; do
    for wt in w.f64 -; do
        if [ "$wt" = - ]; then one_args="ids.u64"; else one_args="ids.u64 $wt"; fi
        "$D2" wsketch $f -S 64 -o one $one_args 2>/dev/null
        "$D2" wsketch $f -S 64 -o csr ids.u64 $wt ip.u64 2>/dev/null
        got=$(same_regs one.sampled.hashes.f64 csr.sampled.regs.stacked.1.64.f64)
        name="wsketch flag '${f:-default}' weights '$wt' 1-D vs CSR registers"
        if [ "$got" = identical ]; then echo "PASS $name: got $got expected identical"
        else echo "FAIL $name: got $got expected identical"; fail=1; fi
        rm -f one.sampled.* csr.sampled.*
    done
done
exit $fail
