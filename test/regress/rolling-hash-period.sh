#!/bin/bash
# Rolling hashes (k above the exact-encoding limit) must not collide
# structurally for k beyond the hash word size. Oracles:
#   - a sequence with period p has exactly p distinct k-mers once k >= p, and
#     two unrelated period-64 repeats share none, so exact J = 0;
#   - two 100-mers differing by swapping positions 0 and 64 are distinct
#     k-mers, so their exact J at k = 100 is 0.
# Checks use --no-canon; the canonical path gives the same results once the
# reverse-strand hash is correct.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
DASHING2=${DASHING2:-$HERE/../../dashing2}
D=$(mktemp -d)
trap 'rm -rf "$D"' EXIT
cd "$D"

python3 - <<'EOF'
import random
random.seed(8)
r = lambda n: "".join(random.choice("ACGT") for _ in range(n))
def w(name, s): open(name, "w").write(">%s\n%s\n" % (name, s))
u, v, x, y = r(64), r(64), r(128), r(128)
w("p64a.fa", u * 20); w("p64b.fa", v * 20)
w("p128a.fa", x * 10); w("p128b.fa", y * 10)
s = list(r(100))
if s[0] == s[64]: s[0] = "A" if s[64] != "A" else "C"
w("s1.fa", "".join(s)); s[0], s[64] = s[64], s[0]; w("s2.fa", "".join(s))
EOF

fail=0
report() { # name got expected
    if [ "$2" = "$3" ]; then echo "PASS $1: got $2 expected $3"
    else echo "FAIL $1: got $2 expected $3"; fail=1; fi
}
# Prints "card(a) card(b) J(a,b)" from exact --set comparisons.
cmpset() {
    local u j
    u=$("$DASHING2" cmp --set --no-canon --union-size --square "$@" 2>/dev/null | awk -F'\t' '!/^#/ {printf "%s ", $(++n + 1)}')
    j=$("$DASHING2" cmp --set --no-canon "$@" 2>/dev/null | awk -F'\t' '!/^#/ {print $3; exit}')
    rm -f ./*.kmerset*
    echo "$u$j" | tr -s ' '
}
report period64_k128 "$(cmpset -k 128 p64a.fa p64b.fa)" "64 64 0"
report period128_k256_128bit "$(cmpset -2 -k 256 p128a.fa p128b.fa)" "128 128 0"
report swap_0_64_k100 "$(cmpset -k 100 s1.fa s2.fa)" "1 1 0"
exit $fail
