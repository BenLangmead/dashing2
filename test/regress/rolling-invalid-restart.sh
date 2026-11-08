#!/bin/bash
# The rolling hasher (k > 32, or k > 64 with -2) must resume right after an N.
# s.fa is 600 bp with N at offsets 150, 160, 300 and 520; runs.fa holds its
# N-free runs as separate records. Their k-mer sets must be identical (J = 1),
# and the cardinality must equal the exact number of distinct k-mers.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
DASHING2=${DASHING2:-$HERE/../../dashing2}
D=$(mktemp -d)
trap 'rm -rf "$D"' EXIT
cd "$D"

python3 - <<'EOF'
import random
random.seed(21)
s = list("".join(random.choice("ACGT") for _ in range(600)))
for p in (150, 160, 300, 520): s[p] = "N"
s = "".join(s)
open("s.fa", "w").write(">s\n" + s + "\n")
open("runs.fa", "w").write("".join(">r%d\n%s\n" % (i, r) for i, r in enumerate(s.split("N"))))
def rc(w): return w[::-1].translate(str.maketrans("ACGT", "TGCA"))
for k in (40, 100):
    for canon in (True, False):
        n = len({min(w, rc(w)) if canon else w for w in (s[i:i + k] for i in range(len(s) - k + 1)) if "N" not in w})
        open("exact_k%d_%s" % (k, "canon" if canon else "nocanon"), "w").write(str(n))
EOF

fail=0
report() { # name got expected
    if [ "$2" = "$3" ]; then echo "PASS $1: got $2 expected $3"
    else echo "FAIL $1: got $2 expected $3"; fail=1; fi
}
check() { # name k exactfile flags...
    local name=$1 k=$2 exact; exact=$(cat "$3"); shift 3
    local j card f
    j=$("$DASHING2" cmp --set -k "$k" "$@" s.fa runs.fa 2>/dev/null | awk -F'\t' '!/^#/ {print $3; exit}' | tr -d ' ')
    rm -f ./*.kmerset*
    "$DASHING2" sketch --set -k "$k" "$@" s.fa 2>/dev/null
    f=$(ls s.fa*.kmerset* | head -1)
    card=$(python3 -c "import struct,sys; print(int(struct.unpack('<d', open(sys.argv[1], 'rb').read(8))[0]))" "$f")
    rm -f ./*.kmerset*
    report "${name}_J_vs_runs" "$j" 1
    report "${name}_cardinality" "$card" "$exact"
}
check rolling_k40          40  exact_k40_canon
check rolling_k40_nocanon  40  exact_k40_nocanon --no-canon
check rolling128_k100      100 exact_k100_canon -2
check rolling128_k100_nocanon 100 exact_k100_nocanon -2 --no-canon
exit $fail
