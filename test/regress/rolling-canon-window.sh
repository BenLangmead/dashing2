#!/bin/bash
# Canonical windowed sketching with k above the exact-encoding limit (rolling
# hashes) must pick one canonical hash per window. Two checks per setting:
#   - every windowed item is also in the unwindowed canonical set;
#   - "sketch -G" emits exactly one item per full window, i.e.
#     (L - k + 1) - (w - k) items for a record of length L.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
DASHING2=${DASHING2:-$HERE/../../dashing2}
D=$(mktemp -d)
trap 'rm -rf "$D"' EXIT
cd "$D"

python3 -c '
import random
random.seed(4)
open("a.fa", "w").write(">a\n%s\n" % "".join(random.choice("ACGT") for _ in range(3000)))'

cat > sets.py <<'EOF'
import sys, struct
def read(path):
    b = open(path, "rb").read()[8:]
    n = 16 if path.endswith("128") else 8
    return {b[i:i + n] for i in range(0, len(b), n)}
full, win = read(sys.argv[1]), read(sys.argv[2])
print(len(win - full))
EOF

fail=0
report() { # name got expected
    if [ "$2" = "$3" ]; then echo "PASS $1: got $2 expected $3"
    else echo "FAIL $1: got $2 expected $3"; fail=1; fi
}
run() { # name k w extra-flags
    local name=$1 k=$2 w=$3 extra=$4
    "$DASHING2" sketch --set -k "$k" $extra a.fa >/dev/null 2>&1; mv a.fa*.kmerset* full.set
    "$DASHING2" sketch --set -k "$k" -w "$w" $extra a.fa >/dev/null 2>&1; mv a.fa*.kmerset* win.set
    report "${name}_windowed_items_not_in_full_set" "$(python3 sets.py full.set win.set)" 0
    "$DASHING2" sketch -G -k "$k" -w "$w" $extra -o out a.fa >/dev/null 2>&1
    local bytes; bytes=$(wc -c < a.fa*.mmerseq* | tr -d ' ')
    local width=8; [ -n "$extra" ] && width=16
    report "${name}_items_emitted" "$((bytes / width))" "$(( (3000 - k + 1) - (w - k) ))"
    rm -f full.set win.set a.fa*.mmerseq* out
}
run k40_w60 40 60 ""
run k100_w120_128bit 100 120 "-2"
exit $fail
