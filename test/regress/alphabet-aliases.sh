#!/bin/bash
# RNA U must read as T, and protein O and U as K and C (aliases declared in
# bonsai alphabet.h). Each check compares a sequence against a copy in which
# one letter is replaced everywhere by its alias; the exact Jaccard must be 1.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
DASHING2=${DASHING2:-$HERE/../../dashing2}
D=$(mktemp -d)
trap 'rm -rf "$D"' EXIT
cd "$D"

python3 - <<'EOF'
import random
random.seed(23)
dna = "".join(random.choice("ACGT") for _ in range(400))
aa = "".join(random.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(400))
def w(name, s): open(name, "w").write(">%s\n%s\n" % (name, s))
w("dna.fa", dna); w("rna.fa", dna.replace("T", "U"))
w("prot.fa", aa); w("prot_O.fa", aa.replace("K", "O")); w("prot_U.fa", aa.replace("C", "U"))
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
check rna_U_as_T_k21         1 -k 21 dna.fa rna.fa
check rna_U_as_T_k40_rolling 1 -k 40 dna.fa rna.fa
check protein_O_as_K         1 --protein -k 5 prot.fa prot_O.fa
check protein_U_as_C         1 --protein -k 5 prot.fa prot_U.fa
exit $fail
