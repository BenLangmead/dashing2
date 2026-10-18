#!/bin/bash
# A record with fewer k-mers than the window (40 bp, k = 21, w = 100) must still
# yield one minimizer, as it does with --no-canon and on the rolling path. The
# canonical direct path used to yield none, leaving an empty sketch.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
DASHING2=${DASHING2:-$HERE/../../dashing2}
D=$(mktemp -d)
trap 'rm -rf "$D"' EXIT
cd "$D"

python3 - <<'EOF'
import random
random.seed(15)
open("s.fa", "w").write(">s\n" + "".join(random.choice("ACGT") for _ in range(40)) + "\n")
EOF

fail=0
check() { # name expected flags...
    local name=$1 exp=$2; shift 2
    rm -f ./*.kmerset*
    "$DASHING2" sketch --set "$@" s.fa 2>/dev/null
    local f got
    f=$(ls s.fa*.kmerset* 2>/dev/null | head -1)
    # File layout: 8-byte cardinality, then one 8-byte (16 with -2) value per minimizer.
    if [ -n "$f" ]; then
        got=$(python3 -c "import os,sys; f=sys.argv[1]; print((os.path.getsize(f) - 8) // (16 if f.endswith('128') else 8))" "$f")
    else got=missing; fi
    rm -f ./*.kmerset*
    if [ "$got" = "$exp" ]; then echo "PASS $name: got $got expected $exp"
    else echo "FAIL $name: got $got expected $exp"; fail=1; fi
}
check short_record_nocanon     1 --no-canon -k 21 -w 100
check short_record_canon       1 -k 21 -w 100
check short_record_canon_128   1 -2 -k 21 -w 100
exit $fail
