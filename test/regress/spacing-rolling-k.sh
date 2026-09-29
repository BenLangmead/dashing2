#!/bin/bash
# --spacing needs exact k-mer encoding; longer k-mers are rolling-hashed, and
# the rolling hasher has no spaced seeds. Such a k must be rejected rather than
# silently sketching contiguous k-mers.
#   - DNA -k 40 --spacing 1x39 (64-bit limit 32): must fail with an error.
#   - protein -k 20 --spacing 0x19 --protein with a real gap (64-bit limit 14):
#     must fail with an error.
#   - DNA -2 -k 40 --spacing 1x39 (128-bit limit 64): must work, and the --set
#     size must equal the number of distinct spaced k-mers (offsets 0, 2, ...,
#     78; spaced k-mers are not canonicalized), computed in python.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
DASHING2=${DASHING2:-$HERE/../../dashing2}
D=$(mktemp -d)
trap 'rm -rf "$D"' EXIT
cd "$D"

python3 - > expected.txt <<'EOF'
import random
random.seed(11)
s = "".join(random.choice("AC") for _ in range(800))
open("dna.fa", "w").write(">d\n%s\n" % s)
p = "".join(random.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(800))
open("prot.fa", "w").write(">p\n%s\n" % p)
print(len({s[i:i + 79:2] for i in range(len(s) - 78)}))
EOF

fail=0
report() { # name got expected
    if [ "$2" = "$3" ]; then echo "PASS $1: got $2 expected $3"
    else echo "FAIL $1: got $2 expected $3"; fail=1; fi
}
rejects() { # name args...
    local name=$1; shift
    local rc=0
    "$DASHING2" sketch --set "$@" > out.txt 2> err.txt || rc=$?
    rm -f ./*.kmerset*
    if [ $rc -ne 0 ] && grep -q -- "--spacing requires" err.txt; then report "$name" "error" "error"
    else report "$name" "rc=$rc, no error" "error"; fi
}
rejects dna_k40_spacing_64bit -k 40 --spacing 1x39 dna.fa
rejects protein_k20_spacing_64bit --protein -k 20 --spacing 0,1,0x17 prot.fa
"$DASHING2" sketch --set -2 -k 40 --spacing 1x39 dna.fa > /dev/null 2>&1
n=$(( ($(wc -c < dna.fa*.kmerset128) - 8) / 16 ))
rm -f ./*.kmerset*
report dna_k40_spacing_128bit_count "$n" "$(cat expected.txt)"
exit $fail
