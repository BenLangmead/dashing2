#!/usr/bin/env bash
# --prob --fastcmp 1 compresses ProbMinHash registers to one byte, where two
# registers holding different items sometimes agree by chance, and the estimate
# must remove those chance agreements as --prob --bbit-sigs, -B --fastcmp 1 and
# --full --fastcmp 1 do. Fixtures (Python, random DNA): x and y are built from
# unrelated 300 bp units repeated 1, 1, 2 and 3 times (probability Jaccard 0);
# z = x's units repeated 2, 1, 1, 3 times plus one unit of y (exact probability
# Jaccard computed in Python with Ertl's J_P formula). Estimates are averaged
# over 10 seeds at -S 4096, where the standard error of one estimate is below 0.01.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
D2="${DASHING2:-$HERE/../../dashing2}"
T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT
cd "$T"
JP=$(python3 - <<'PY'
import random
from collections import Counter
R = random.Random(5)
rnd = lambda n: ''.join(R.choice('ACGT') for _ in range(n))
ux = [rnd(300) for _ in range(4)]
uy = [rnd(300) for _ in range(4)]
spec = {'x': [(ux[j], c) for j, c in enumerate((1, 1, 2, 3))],
        'y': [(uy[j], c) for j, c in enumerate((1, 1, 2, 3))],
        'z': [(ux[j], c) for j, c in enumerate((2, 1, 1, 3))] + [(uy[0], 1)]}
def canon_counts(recs, k=21):
    c = Counter()
    for s in recs:
        rc = s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
        for i in range(len(s) - k + 1):
            c[min(s[i:i+k], rc[len(s)-i-k:len(s)-i])] += 1
    return c
C = {}
for name, units in spec.items():
    recs = [u for u, n in units for _ in range(n)]
    open(name + '.fa', 'w').write(''.join('>%s_%d\n%s\n' % (name, r, s) for r, s in enumerate(recs)))
    C[name] = canon_counts(recs)
A, B = C['x'], C['z']
sa, sb = sum(A.values()), sum(B.values())
pairs = Counter((A.get(x, 0) / sa, B.get(x, 0) / sb) for x in set(A) | set(B))
print(sum(1 / sum(n * max(p / (A[x] / sa), q / (B[x] / sb)) for (p, q), n in pairs.items())
          for x in set(A) & set(B)))
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
# Mean over seeds 1..10 of the estimate for the pair (first file, second file).
mean() {
    for s in $(seq 10); do
        "$D2" cmp -k 21 -S 4096 --seed $s "$@" 2>/dev/null | awk '!/^#/ {print $3; exit}'
    done | awk '{s += $1; n++} END {if (n == 10) printf "%.5f\n", s / n}'
}
# Clamping at 0 leaves a small positive mean at probability Jaccard 0 (about 0.0005).
check unrelated_fastcmp1 "$(mean --prob --fastcmp 1 x.fa y.fa)" 0 0.003
check related_fastcmp1 "$(mean --prob --fastcmp 1 x.fa z.fa)" "$JP" 0.01
check unrelated_fastcmp2 "$(mean --prob --fastcmp 2 x.fa y.fa)" 0 0.003
rm -f ./*.kmerset* ./*.pmh
exit $fail
