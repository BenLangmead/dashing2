#!/usr/bin/env bash
# The first line of -Q (query by reference) output names the layout; it read
# "#Dashing2 Panel (Query/Refernce) Output".
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 -c "
import random
r = random.Random(1)
for n in 'ab':
    open(n + '.fa', 'w').write('>' + n + '\n' + ''.join(r.choice('ACGT') for _ in range(2000)) + '\n')
"
echo b.fa > q.txt
got=$("$D" cmp -k 15 -S 64 -Q q.txt a.fa 2>/dev/null | head -n 1)
want="#Dashing2 Panel (Query/Reference) Output"
if [ "$got" = "$want" ]; then echo "PASS panel header: got $got expected $want"; exit 0; fi
echo "FAIL panel header: got ${got:-nothing} expected $want"; exit 1
