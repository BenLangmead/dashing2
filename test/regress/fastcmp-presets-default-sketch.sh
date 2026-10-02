#!/usr/bin/env bash
# --fastcmp-bytes/-shorts/-words and --setsketch-ab build compressed SetSketches,
# which exist only as full SetSketches. Without --full they aborted ("Sketch
# compressed is only available for FullSetSketch"); they now select --full, so
# each run must equal the same run with --full.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
r = random.Random(1)
a = ''.join(r.choice('ACGT') for _ in range(4000))
b = a[:2500] + ''.join(r.choice('ACGT') for _ in range(1500))
open('a.fa', 'w').write('>a\n' + a + '\n')
open('b.fa', 'w').write('>b\n' + b + '\n')
PY
fail=0
v() { "$D" "$@" a.fa b.fa 2>/dev/null | awk '/^a.fa/{print $3}'; }
check() {  # check <name> <got> <expected>
    if [ -n "$3" ] && [ "$2" = "$3" ]; then echo "PASS $1: got $2 expected $3"
    else echo "FAIL $1: got ${2:-nothing} expected ${3:-a value}"; fail=1; fi
}
for p in --fastcmp-bytes --fastcmp-shorts --fastcmp-words "--setsketch-ab 0.4,1.005 --fastcmp 2"; do
    check "cmp $p" "$(v cmp -k 21 $p)" "$(v cmp -k 21 --full $p)"
done
check "sketch --fastcmp-bytes" "$(v sketch -k 21 --fastcmp-bytes --cmpout -)" "$(v sketch -k 21 --full --fastcmp-bytes --cmpout -)"
exit $fail
