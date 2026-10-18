#!/bin/bash
# Canonical windowed sketches (-w) of a sequence with one N must contain only
# k-mers of the input. The k-mers that span the N used to enter the window as
# the all-A k-mer. With seed 0, sketch --set writes WangHash of each 2-bit
# canonical k-mer (A=0 C=1 G=2 T=3; with -2, WangHash of each 64-bit half), so
# the expected superset is computed exactly below.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
DASHING2=${DASHING2:-$HERE/../../dashing2}
D=$(mktemp -d)
trap 'rm -rf "$D"' EXIT
cd "$D"

python3 - <<'EOF'
import random
random.seed(14)
s = "".join(random.choice("ACGT") for _ in range(400))
open("s.fa", "w").write(">s\n" + s[:200] + "N" + s[201:] + "\n")
EOF

fail=0
check() { # name k flags...
    local name=$1 k=$2; shift 2
    rm -f ./*.kmerset*
    "$DASHING2" sketch --set -k "$k" "$@" s.fa 2>/dev/null
    local res
    res=$(K=$k python3 - <<'EOF'
import glob, os, struct
M = (1 << 64) - 1
def wang(key):
    key = ((~key) + (key << 21)) & M; key ^= key >> 24
    key = (key + (key << 3) + (key << 8)) & M; key ^= key >> 14
    key = (key + (key << 2) + (key << 4)) & M; key ^= key >> 28
    return (key + (key << 31)) & M
k = int(os.environ["K"])
s = open("s.fa").read().split("\n")[1]
exp = set()
for i in range(len(s) - k + 1):
    w = s[i:i + k]
    if "N" in w: continue
    r = w[::-1].translate(str.maketrans("ACGT", "TGCA"))
    v = 0
    for c in min(w, r): v = (v << 2) | "ACGT".index(c)
    exp.add(wang(v) if k <= 32 else wang(v & M) | (wang(v >> 64) << 64))
f = glob.glob("s.fa*.kmerset*")[0]
b = open(f, "rb").read()[8:]
w = 16 if f.endswith("128") else 8
got = {int.from_bytes(b[i:i + w], "little") for i in range(0, len(b), w)}
print(len(got - exp))
EOF
)
    rm -f ./*.kmerset*
    if [ "$res" = "0" ]; then echo "PASS $name: got $res minimizers outside the input's k-mers, expected 0"
    else echo "FAIL $name: got $res minimizers outside the input's k-mers, expected 0"; fail=1; fi
}
check canon_window_N_k21     21 -w 26
check canon_window_N_k32     32 -w 37
check canon_window_N_k40_128 40 -w 45 -2
exit $fail
