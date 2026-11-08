#!/bin/bash
# End-to-end checks of the core bonsai fixes (canonical hashing for k > 32,
# 128-bit canonicalization, and protein k-mer encoding), each against an exact
# Python count. All of these fail when dashing2 is built against bonsai
# without those fixes.
# Usage: DASHING2=/path/to/dashing2 test/regress/bonsai-core.sh
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - > exp.txt <<'PY'
import random
r = random.Random(11)
rc = lambda s: s[::-1].translate(str.maketrans("ACGT", "TGCA"))
x = "".join(r.choice("ACGT") for _ in range(3000))
y = "".join(r.choice("ACGT") for _ in range(3000))
xm = "".join(r.choice("ACGT".replace(c, "")) if i % 60 == 0 else c for i, c in enumerate(x))
aa = "ACDEFGHIKLMNPQRSTVWY"
p = "".join(r.choice(aa) for _ in range(3000))
pm = "".join(r.choice(aa.replace(c, "")) if r.random() < 0.02 else c for c in p)
for name, s in [("x", x), ("x_rc", rc(x)), ("y", y), ("xm", xm), ("p", p), ("pm", pm)]:
    open(name + ".fa", "w").write(f">{name}\n{s}\n")
def canon(s, k):
    t = rc(s); n = len(s)
    return {min(s[i:i+k], t[n-i-k:n-i]) for i in range(n - k + 1)}
def plain(s, k): return {s[i:i+k] for i in range(len(s) - k + 1)}
J = lambda a, b: len(a & b) / len(a | b)
print("B1_rc_k40", 1.0)
print("B2_unrelated_-2_k36", J(canon(x, 36), canon(y, 36)))
print("B2_mutated_-2_k48", J(canon(x, 48), canon(xm, 48)))
print("B6_protein_k14", J(plain(p, 14), plain(pm, 14)))
PY
v() { rm -f ./*.kmerset*; "$D" cmp --set "$@" 2>/dev/null | awk -F'\t' '!/^#/{print $3; exit}'; }
fail=0
check() {
    exp=$(awk -v n="$1" '$1==n{print $2}' exp.txt)
    if python3 -c "import sys; sys.exit(0 if abs(float('$2') - $exp) < 1e-3 else 1)" 2>/dev/null; then
        echo "PASS $1: got $2 expected $exp"
    else
        echo "FAIL $1: got $2 expected $exp"; fail=1
    fi
}
check B1_rc_k40 "$(v -k 40 x.fa x_rc.fa)"
check B2_unrelated_-2_k36 "$(v -2 -k 36 x.fa y.fa)"
check B2_mutated_-2_k48 "$(v -2 -k 48 x.fa xm.fa)"
check B6_protein_k14 "$(v --protein -k 14 p.fa pm.fa)"
exit $fail
