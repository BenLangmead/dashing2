#!/bin/bash
# A byte in 0x80-0xff (for example from UTF-8 text or a corrupt file) is not a
# base. A sequence with such a byte at offset 150 must have exactly the k-mers
# of the same sequence with N there, so the exact Jaccard (cmp --set) must be 1.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
DASHING2=${DASHING2:-$HERE/../../dashing2}
D=$(mktemp -d)
trap 'rm -rf "$D"' EXIT
cd "$D"

BYTES="80 c1 e5 ff"
python3 - $BYTES <<'EOF'
import random, sys
random.seed(13)
s = "".join(random.choice("ACGT") for _ in range(301)).encode()
open("N.fa", "wb").write(b">s\n" + s[:150] + b"N" + s[151:] + b"\n")
for h in sys.argv[1:]:
    open(h + ".fa", "wb").write(b">s\n" + s[:150] + bytes([int(h, 16)]) + s[151:] + b"\n")
EOF

fail=0
check() { # name expected flags... file1 file2
    local name=$1 exp=$2; shift 2
    local got
    got=$("$DASHING2" cmp --set "$@" 2>/dev/null | awk -F'\t' '!/^#/ {print $3; exit}' | tr -d ' ')
    rm -f ./*.kmerset*
    if [ "$got" = "$exp" ]; then echo "PASS $name: got $got expected $exp"
    else echo "FAIL $name: got $got expected $exp"; fail=1; fi
}
for b in $BYTES; do
    check "byte_0x${b}_direct_k21"         1 -k 21 N.fa $b.fa
    check "byte_0x${b}_direct128_k40"      1 -2 -k 40 N.fa $b.fa
    check "byte_0x${b}_rolling_k40"        1 -k 40 N.fa $b.fa
    check "byte_0x${b}_rolling_k40_nocanon" 1 --no-canon -k 40 N.fa $b.fa
done
exit $fail
