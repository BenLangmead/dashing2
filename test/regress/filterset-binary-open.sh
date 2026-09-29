#!/usr/bin/env bash
# `--filterset path:b` reads a binary file of k-mers (for example the
# .mmerseq64 file written by `sketch -G`) and removes them before sketching.
# The --set union of a file with itself is then the number of its distinct
# k-mers not in the filter, the same as when the filter is given as sequence.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
# x: 50 kbp random; half: its first 25 kbp, used as the filter.
python3 - <<'PY'
import random
r = random.Random(6)
x = ''.join(r.choice('ACGT') for _ in range(50000))
open('x.fa', 'w').write(f'>x\n{x}\n')
open('half.fa', 'w').write(f'>half\n{x[:25000]}\n')
PY
# Exact count of distinct canonical 21-mers of x that are absent from half.
exp=$(python3 - <<'PY'
k = 21; comp = str.maketrans('ACGT', 'TGCA')
def kset(fn):
    s = ''.join(l.strip() for l in open(fn) if not l.startswith('>'))
    return {min(s[i:i + k], s[i:i + k].translate(comp)[::-1]) for i in range(len(s) - k + 1)}
print(len(kset('x.fa') - kset('half.fa')))
PY
)
"$D2" sketch -G -k 21 -o half.stacked half.fa 2>/dev/null
mv half.fa.*.mmerseq64 half.u64
gzip -c half.u64 > half.u64.gz
for fs in half.fa half.u64:b half.u64.gz:b; do
    got=$("$D2" cmp -k 21 --set --union-size --filterset $fs x.fa x.fa 2>/dev/null | grep -v '^#' | awk 'NR == 1 {print $3}')
    rm -f ./*.kmer*
    if [ "$got" = "$exp" ]; then echo "PASS --filterset $fs union: got $got expected $exp"
    else echo "FAIL --filterset $fs union: got ${got:-nothing} expected $exp"; fail=1; fi
done
exit $fail
