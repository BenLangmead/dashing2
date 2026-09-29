#!/usr/bin/env bash
# Options that change a sketch (--downsample, --entmin, --filterset) must be
# part of the --cache key. Invariant: after a run without the option has
# populated the cache, a cached run with the option equals a fresh run with it.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
r = random.Random(1)
a = ''.join(r.choice('ACGT') for _ in range(5000))
b = list(a)
for i in range(0, len(b), 97):
    b[i] = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[b[i]]
open('a.fa', 'w').write('>a\n' + a + '\n')
open('b.fa', 'w').write('>b\n' + ''.join(b) + '\n')
open('f.fa', 'w').write('>f\n' + a[:2500] + '\n')   # filter set: first half of a
PY
val() { "$D" cmp --full -k 21 "$@" a.fa b.fa 2>/dev/null | awk 'NR==4{print $3}'; }
fail=0
check() {  # check <name> <options shared by both runs> <option under test>
    fresh=$(val $2 $3)
    val $2 --cache >/dev/null
    cached=$(val $2 $3 --cache)
    rm -f a.fa.* b.fa.*
    if [ -n "$fresh" ] && [ "$cached" = "$fresh" ]; then
        echo "PASS $1: got $cached expected $fresh"
    else
        echo "FAIL $1: got $cached expected $fresh"; fail=1
    fi
}
check downsample "" "--downsample 0.3"
check entmin "-w 30" "--entmin"
check filterset "" "--filterset f.fa"
exit $fail
