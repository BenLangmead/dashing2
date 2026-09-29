#!/usr/bin/env bash
# --downsample F must keep or drop each distinct k-mer as a whole, so that a
# file and an identical copy stay identical (J = 1), the sampled Jaccard stays
# close to the exact one, and the result does not depend on -p.
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
# x: 50 kbp random; xcopy: identical; xhalf: its first 25 kbp.
python3 - <<'PY'
import random
r = random.Random(4)
x = ''.join(r.choice('ACGT') for _ in range(50000))
for name, s in [('x', x), ('xcopy', x), ('xhalf', x[:25000])]:
    open(name + '.fa', 'w').write(f'>{name}\n{s}\n')
PY
# xhalf's 21-mers are a subset of x's, so the exact Jaccard is |xhalf| / |x|.
exact=$(python3 -c "print((25000 - 20) / (50000 - 20))")
run() { "$D2" cmp -k 21 -S 8192 --downsample 0.5 "$@" x.fa xcopy.fa xhalf.fa 2>/dev/null; rm -f ./*.kmer*; }
entry() { grep -v '^#' | awk -v c="$1" 'NR == 1 {print $(c + 1)}'; }
for mode in --set --full -B; do
    out=$(run $mode -p 1)
    check "$mode x~xcopy" "$(entry 2 <<<"$out")" 1 0
    check "$mode x~xhalf (+-0.03)" "$(entry 3 <<<"$out")" "$exact" 0.03
    same=0
    for rep in 1 2 3; do [ "$(run $mode -p 4)" = "$out" ] && same=$((same + 1)); done
    check "$mode -p 4 runs identical to -p 1 (of 3)" "$same" 3 0
done
exit $fail
