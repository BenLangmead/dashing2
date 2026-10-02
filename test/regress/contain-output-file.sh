#!/usr/bin/env bash
# With -o, contain must write its whole table to the file and nothing to
# stdout. x86 builds format the columns in blocks of 8 (AVX2) or 16 (AVX-512)
# and the rest one at a time, so this uses 21 references to reach both.
# Expected file contents: the table contain prints to stdout without -o.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
r = random.Random(5)
refs = [''.join(r.choice('ACGT') for _ in range(2000)) for _ in range(21)]
for i, s in enumerate(refs):
    open('r%02d.fa' % i, 'w').write('>r%02d\n%s\n' % (i, s))
open('q.fa', 'w').write(''.join('>q%d\n%s\n' % (i, s[:300 + 80 * i]) for i, s in enumerate(refs)))
PY
"$D" sketch --full -k 21 --save-kmers -o db r*.fa >/dev/null 2>&1
"$D" contain db.kmer64 q.fa > expected.tsv 2>/dev/null
"$D" contain -o out.tsv db.kmer64 q.fa > stdout.txt 2>/dev/null
fail=0
n=$(wc -c < stdout.txt | tr -d ' ')
if [ "$n" = 0 ]; then echo "PASS stdout with -o: got $n bytes expected 0"
else echo "FAIL stdout with -o: got $n bytes expected 0"; fail=1; fi
g=$(awk -F'\t' '$1=="q.fa"{print NF-1}' out.tsv); e=$(awk -F'\t' '$1=="q.fa"{print NF-1}' expected.tsv)
if cmp -s out.tsv expected.tsv; then echo "PASS file contents: got $g columns, same as stdout expected $e columns, same as stdout"
else echo "FAIL file contents: got $g columns, differs from stdout expected $e columns, same as stdout"; fail=1; fi
rm -f ./*.kmerset* ./*.kmer64* db*
exit $fail
