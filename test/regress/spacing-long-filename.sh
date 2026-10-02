#!/bin/bash
# Exact modes and --cache name their files after the spaced seed. A long
# irregular seed spelled out in full made the name longer than the 255-byte
# file name limit, so sketch --set aborted with "Failed to open". Such seeds
# are named by a digest; short seeds keep their spelled-out names, so
# existing caches still match. Expected --set size: the number of distinct
# spaced k-mers of x (random, so every placement is distinct), counted below.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
DASHING2=${DASHING2:-$HERE/../../dashing2}
D=$(mktemp -d)
trap 'rm -rf "$D"' EXIT
cd "$D"

python3 - <<'EOF'
import random
random.seed(5)
open("x.fa", "w").write(">x\n" + "".join(random.choice("ACGT") for _ in range(3000)) + "\n")
EOF
LONG=$(python3 -c "print(','.join(str(1 + i % 2) for i in range(63)))")   # k = 64, gaps 1,2,1,2,...

fail=0
report() { # name got expected
    if [ "$2" = "$3" ]; then echo "PASS $1: got $2 expected $3"; else echo "FAIL $1: got $2 expected $3"; fail=1; fi
}
nspaced() { # spacing -> number of distinct spaced k-mers of x
    python3 - "$1" <<'EOF'
import sys
gaps = [int(g) for g in sys.argv[1].split(",")]
offs = [0]
for g in gaps: offs.append(offs[-1] + g + 1)
s = open("x.fa").read().split("\n")[1]
print(len({"".join(s[i + o] for o in offs) for i in range(len(s) - offs[-1])}))
EOF
}
setsize() { # kmerset file -> number of stored items
    python3 -c "import os, sys; f = sys.argv[1]; print((os.path.getsize(f) - 8) // (16 if f.endswith('128') else 8))" "$1"
}

("$DASHING2" sketch --set -k 64 -2 --spacing "$LONG" x.fa) >/dev/null 2>&1
report long_seed_set_exit $? 0
f=$(ls x.fa*.kmerset* 2>/dev/null | head -1)
report long_seed_set_size "$( [ -n "$f" ] && setsize "$f")" "$(nspaced "$LONG")"
report long_seed_name_length "$( [ -n "$f" ] && [ ${#f} -le 255 ] && echo ok)" ok
rm -f ./*.kmerset*

("$DASHING2" sketch --cache -k 64 -2 --spacing "$LONG" x.fa) >/dev/null 2>&1
report long_seed_cache_write $? 0
("$DASHING2" sketch --cache -k 64 -2 --spacing "$LONG" x.fa) >/dev/null 2>&1
report long_seed_cache_reload $? 0
rm -f x.fa?*

("$DASHING2" sketch --set -k 8 --spacing 0,1,0,2,0,0,1 x.fa) >/dev/null 2>&1
report short_seed_name "$(ls x.fa*.kmerset* 2>/dev/null)" "x.fa1x1,2x1,1x1,3x1,1x2,2x1.k8.FullMmerSet.DNA.kmerset64"
rm -f ./*.kmerset*
exit $fail
