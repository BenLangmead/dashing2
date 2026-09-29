#!/usr/bin/env bash
# contain reads reference names from <db>.names.txt. Without that file it must
# still run, with one column per reference (2 here), labelled by index, and x.fa
# must cover its own entry (column 0) 100%.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
r = random.Random(1)
for name in 'xy':
    open(name + '.fa', 'w').write('>' + name + '\n' + ''.join(r.choice('ACGT') for _ in range(5000)) + '\n')
PY
"$D" sketch --full -k 21 --save-kmers -o db x.fa y.fa >/dev/null 2>&1
rm -f db.kmer64.names.txt
out=$("$D" contain db.kmer64 x.fa 2>/dev/null)
refs=$(echo "$out" | awk -F'\t' '/^##References/{$1=""; print substr($0, 2)}')
self=$(echo "$out" | awk '$1=="x.fa"{print $2}')
if [ "$refs" = "0 1" ] && [ "$self" = "100%:1" ]; then
    echo "PASS names fallback: got references '$refs', x.fa $self expected '0 1', 100%:1"; exit 0
else
    echo "FAIL names fallback: got references '$refs', x.fa '$self' expected '0 1', 100%:1"; exit 1
fi
