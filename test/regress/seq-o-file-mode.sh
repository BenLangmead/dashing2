#!/usr/bin/env bash
# `sketch -G -o out` without --parse-by-seq writes one minimizer sequence per
# input file to the stacked file `out`, which `printmin` reads. With w = k the
# sequence of a file is every canonical k-mer of each of its records in order,
# computed here in Python. Before the fix the file had the size of a sketch
# file, its body was zeros, and printmin rejected it.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
python3 - <<'PY'
import random
r = random.Random(8)
seq = lambda n: ''.join(r.choice('ACGT') for _ in range(n))
open('a.fa', 'w').write('>a1\n%s\n>a2\n%s\n' % (seq(700), seq(400)))
open('b.fa', 'w').write('>b1\n%s\n>short\nACGTACG\n>b2\n%s\n' % (seq(900), seq(30)))
comp = str.maketrans('ACGT', 'TGCA')
with open('expected.txt', 'w') as f:
    for i, fn in enumerate(('a.fa', 'b.fa')):
        recs = [l.strip() for l in open(fn) if not l.startswith('>')]
        kmers = [min(s[j:j + 15], s[j:j + 15].translate(comp)[::-1]) for s in recs for j in range(len(s) - 14)]
        f.write('MinimizerSequence%d %s\n' % (i, ' '.join(kmers)))
PY
for run in plain "--seed 7" "--cache" "--cache (reloaded)"; do
    opts=${run% (reloaded)}; [ "$opts" = plain ] && opts=
    [ "$run" = "--cache" ] && rm -f ./*.mmerseq64
    rm -f out out.*
    "$D2" sketch -k 15 -G $opts -o out a.fa b.fa >/dev/null 2>&1
    if "$D2" printmin out >got.txt 2>/dev/null && cmp -s got.txt expected.txt; then
        echo "PASS -G -o $run: printmin output equals the k-mer oracle ($(wc -c <out | tr -d ' ') bytes)"
    else
        echo "FAIL -G -o $run: printmin output differs from the k-mer oracle ($(wc -c <out | tr -d ' ') bytes; $(head -c 80 got.txt))"; fail=1
    fi
done
rm -f ./*.mmerseq64 out out.*
exit $fail
