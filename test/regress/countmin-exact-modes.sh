#!/usr/bin/env bash
# -c/--countmin-size approximates weights for BagMinHash and ProbMinHash only.
# The exact k-mer modes (--set, -J/--countdict) must give the same, exact
# answers with or without it.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
check() { # name got expected tolerance
    if python3 -c "import sys; sys.exit(not abs(float('$2') - float('$3')) <= $4)" 2>/dev/null; then
        echo "PASS $1: got $2 expected $3"
    else
        echo "FAIL $1: got $2 expected $3"; fail=1
    fi
}
# x: 3 kbp random; y: x with every 20th base changed.
python3 - <<'PY'
import random
r = random.Random(7)
x = ''.join(r.choice('ACGT') for _ in range(3000))
y = ''.join('ACGT'[('ACGT'.index(c) + 1) % 4] if i % 20 == 10 else c for i, c in enumerate(x))
for name, s in [('x', x), ('y', y)]:
    open(name + '.fa', 'w').write(f'>{name}\n{s}\n')
PY
# Exact Jaccard of the canonical 15-mer sets.
exact=$(python3 - <<'PY'
def kmers(s, k=15):
    rc = s[::-1].translate(str.maketrans('ACGT', 'TGCA'))
    n = len(s)
    return {min(s[i:i + k], rc[n - i - k:n - i]) for i in range(n - k + 1)}
seq = lambda f: open(f).read().split('\n')[1]
a, b = kmers(seq('x.fa')), kmers(seq('y.fa'))
print(len(a & b) / len(a | b))
PY
)
entry() { awk '$1 == "x.fa" {print $3}'; }
run() { "$D2" cmp -k 15 "$@" x.fa y.fa 2>/dev/null; rm -f ./*.kmer*; }
check "--set exact" "$(run --set | entry)" "$exact" 1e-6
check "--set -c 16" "$(run --set -c 16 | entry)" "$exact" 1e-6
check "-J -c 256 equals -J" "$(run -J -c 256 | entry)" "$(run -J | entry)" 1e-9
exit $fail
