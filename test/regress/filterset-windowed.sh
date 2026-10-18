#!/bin/bash
# --filterset F removes every k-mer of F from the inputs, also when sketching
# minimizers (-w). The filter used to hold only F's own minimizers, so input
# minimizers that occur in F as non-minimizers were kept.
# Expected: the unfiltered windowed --set of x minus every k-mer of F. For
# k = 15 the k-mers of F are computed below (with seed 0, --set stores
# WangHash of each 2-bit k-mer, A=0 C=1 G=2 T=3, canonical unless --no-canon);
# for k = 40 (rolling hash) they are the unwindowed --set of F.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
DASHING2=${DASHING2:-$HERE/../../dashing2}
D=$(mktemp -d)
trap 'rm -rf "$D"' EXIT
cd "$D"

python3 - <<'EOF'
import random
random.seed(31)
s = "".join(random.choice("ACGT") for _ in range(3000))
open("x.fa", "w").write(">x\n" + s + "\n")
# The filter holds 50 bp pieces of x, one every 100 bp.
open("f.fa", "w").write("".join(">f%d\n%s\n" % (i, s[i:i + 50]) for i in range(0, 3000, 100)))
EOF

stored() { # sketch --set options... -> stored values of the result, one per line
    rm -f ./*.kmerset*
    "$DASHING2" sketch --set "$@" 2>/dev/null
    python3 -c "
import glob
b = open(glob.glob('*.kmerset*')[0], 'rb').read()[8:]
print('\n'.join(str(int.from_bytes(b[i:i + 8], 'little')) for i in range(0, len(b), 8)))"
    rm -f ./*.kmerset*
}
python_filter() { # canon(0/1) -> WangHash of every 15-mer of f.fa
    CANON=$1 python3 - <<'EOF'
import os
M = (1 << 64) - 1
def wang(key):
    key = ((~key) + (key << 21)) & M; key ^= key >> 24
    key = (key + (key << 3) + (key << 8)) & M; key ^= key >> 14
    key = (key + (key << 2) + (key << 4)) & M; key ^= key >> 28
    return (key + (key << 31)) & M
k = 15
for line in open("f.fa"):
    if line.startswith(">"): continue
    f = line.strip()
    for i in range(len(f) - k + 1):
        w = f[i:i + k]
        if os.environ["CANON"] == "1": w = min(w, w[::-1].translate(str.maketrans("ACGT", "TGCA")))
        v = 0
        for c in w: v = (v << 2) | "ACGT".index(c)
        print(wang(v))
EOF
}
fail=0
check() { # name filter-values-file sketch options...
    local name=$1 filt=$2; shift 2
    stored "$@" x.fa > all.txt
    stored "$@" --filterset f.fa x.fa > kept.txt
    local res
    res=$(python3 -c "
allm = set(open('all.txt').read().split()); kept = set(open('kept.txt').read().split())
filt = set(open('$filt').read().split())
print(len(kept), len(allm - filt))")
    set -- $res
    if [ "$1" = "$2" ] && [ "$1" -gt 0 ]; then echo "PASS $name: got $1 minimizers kept expected $2"
    else echo "FAIL $name: got $1 minimizers kept expected $2"; fail=1; fi
}
python_filter 1 > f15.txt
check filterset_w30_k15 f15.txt -k 15 -w 30
python_filter 0 > f15nc.txt
check filterset_w30_k15_nocanon f15nc.txt -k 15 -w 30 --no-canon
stored -k 40 f.fa > f40.txt
check filterset_w60_k40_rolling f40.txt -k 40 -w 60
exit $fail
