#!/usr/bin/env bash
# The contain header is printed with fmt, which does not use printf escapes, so
# "%%" must not appear in the output.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 -c "
import random
r = random.Random(1)
open('x.fa', 'w').write('>x\n' + ''.join(r.choice('ACGT') for _ in range(5000)) + '\n')
"
"$D" sketch --full -k 21 --save-kmers -o db x.fa >/dev/null 2>&1
n=$("$D" contain db.kmer64 x.fa 2>/dev/null | grep '^#' | grep -c '%%')
if [ "$n" -eq 0 ]; then echo "PASS header: got $n lines with %% expected 0"; exit 0
else echo "FAIL header: got $n lines with %% expected 0"; exit 1; fi
