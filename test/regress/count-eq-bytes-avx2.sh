#!/usr/bin/env bash
# With --bbit-sigs, --fastcmp 1 compares one-byte signatures and must estimate
# the Jaccard similarity as well as --fastcmp 2 and 4 do. AVX2 builds count the
# equal bytes in blocks of 32, so the default sketch size (1024 registers)
# reaches the vector code.
# Fixtures: a = 200 kbp random DNA, b = a[:100000] + 100 kbp of new sequence,
# so J is about 1/3. The exact J over canonical 21-mers is computed in Python.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
J=$(python3 - <<'PY'
import random
R = random.Random(9)
rnd = lambda n: ''.join(R.choice('ACGT') for _ in range(n))
a = rnd(200000)
b = a[:100000] + rnd(100000)
open('a.fa', 'w').write('>a\n' + a + '\n')
open('b.fa', 'w').write('>b\n' + b + '\n')
def kmers(x, k=21):
    rc = x.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    n = len(x)
    return {min(x[i:i+k], rc[n-i-k:n-i]) for i in range(n - k + 1)}
A, B = kmers(a), kmers(b)
print(len(A & B) / len(A | B))
PY
)
fail=0
for n in 1 2 4; do
    got=$("$D" cmp -k 21 --bbit-sigs --fastcmp $n a.fa b.fa 2>/dev/null | awk '$1 == "a.fa" {print $3}')
    # 1024 registers give a standard error of about 0.015 at J = 1/3; allow 0.05.
    if python3 -c "import sys; sys.exit(not abs(float('${got:-nan}') - $J) <= 0.05)"; then
        echo "PASS similarity_fastcmp$n: got $got expected $J"
    else
        echo "FAIL similarity_fastcmp$n: got ${got:-nothing} expected $J"; fail=1
    fi
done
rm -f ./*.kmerset* ./*.opss
exit $fail
