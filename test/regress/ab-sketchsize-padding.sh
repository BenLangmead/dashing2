#!/usr/bin/env bash
# Directly sketched compressed SetSketches (--full with --fastcmp-bytes,
# --fastcmp-shorts, --fastcmp-words, or --setsketch-ab A,B --fastcmp N) must
# estimate Jaccard for any sketch size, including sizes whose registers do not
# fill whole 64-bit words (-S 500 and -S 1001 bytes, -S 1001 shorts).
# Fixtures: a = 20 kbp random DNA, b = a[:12000] + 5 kbp of new sequence (exact
# Jaccard of canonical 21-mers computed in Python, about 0.48), and a2 = an
# identical copy of a, which must give exactly 1.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
D2="${DASHING2:-$HERE/../../dashing2}"
T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
cd "$T"
J=$(python3 - <<'PY'
import random
R = random.Random(1)
rnd = lambda n: ''.join(R.choice('ACGT') for _ in range(n))
a = rnd(20000)
b = a[:12000] + rnd(5000)
for name, s in (('a', a), ('b', b), ('a2', a)):
    open(name + '.fa', 'w').write('>' + name + '\n' + s + '\n')
def kmers(x, k=21):
    rc = x.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    n = len(x)
    return {min(x[i:i+k], rc[n-i-k:n-i]) for i in range(n - k + 1)}
A, B = kmers(a), kmers(b)
print(len(A & B) / len(A | B))
PY
)
fail=0
check() { # name got expected tolerance
    if [ -n "$2" ] && python3 -c "import sys; sys.exit(abs($2 - $3) > $4)"; then
        echo "PASS $1: got $2 expected $3 (tolerance $4)"
    else
        echo "FAIL $1: got '$2' expected $3 (tolerance $4)"; fail=1
    fi
}
for S in 500 1001; do
    for opt in --fastcmp-bytes --fastcmp-shorts --fastcmp-words "--setsketch-ab 20,1.2 --fastcmp 1"; do
        out=$("$D2" cmp -k 21 -S $S -p 1 --full $opt a.fa b.fa a2.fa 2>/dev/null)
        name="S$S${opt// /_}"
        # The standard error of a Jaccard estimate from 500 registers is about 0.022.
        check "related_$name" "$(echo "$out" | awk '$1 == "a.fa" {print $3}')" "$J" 0.1
        check "identical_$name" "$(echo "$out" | awk '$1 == "a.fa" {print $4}')" 1 0
    done
done
rm -f ./*.kmerset* ./*.ss ./*.opss
exit $fail
