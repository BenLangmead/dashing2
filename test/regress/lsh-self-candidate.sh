#!/usr/bin/env bash
# --topk and --similarity-threshold ask the LSH index for up to n - 1 candidates
# per input; the input's own id must not use up one of those slots. Fixtures: six
# inputs, g00, g01, g02 and g04 identical, g03 the same sequence with 2.5%
# substitutions and g05 with 6%, so every pair shares k-mers. With n = 6 each
# query can return all five other inputs, so --similarity-threshold 0.05 and
# --topk 5 must list every pair, with the values of the --square matrix. The
# lists are compared as sets at -p 1: with more threads, v2.1.20 can also list a
# neighbour twice (a separate race in the neighbour-list update).
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
r = random.Random(3)
a = ''.join(r.choice('ACGT') for _ in range(3000))
def mut(s, x): return ''.join(c if r.random() > x else r.choice('ACGT'.replace(c, '')) for c in s)
b = mut(a, 0.025); c = mut(a, 0.06)
for i, s in enumerate([a, a, a, b, a, c]):
    open('g%02d.fa' % i, 'w').write('>g%02d\n%s\n' % (i, s))
PY
F="g00.fa g01.fa g02.fa g03.fa g04.fa g05.fa"
O="--full -k 20 -S 512 -p 1"
"$D" cmp $O --square $F > square.txt 2>/dev/null
fail=0
check() { # name, expected, got
    if [ "$2" = "$3" ]; then echo "PASS $1: got $3 expected $2"; else echo "FAIL $1: got $3 expected $2"; fail=1; fi
}
# Pairs "a-b" (a < b) with their values, from a neighbour-list file or the matrix.
pairs() { python3 -c '
import sys
p = set()
for l in open(sys.argv[1]):
    if l.startswith("#"): continue
    f = l.split()
    if sys.argv[2] == "square":
        p |= {tuple(sorted((f[0], n))) + (round(float(x), 4),) for n, x in zip(sys.argv[3:], f[1:]) if n != f[0] and float(x) >= 0.05}
    else:
        p |= {tuple(sorted((f[0], x.rsplit(":", 1)[0]))) + (round(float(x.rsplit(":", 1)[1]), 4),) for x in f[1:]}
print(" ".join("%s-%s:%g" % t for t in sorted(p)))' "$@"; }
exp=$(pairs square.txt square $F)
"$D" cmp $O --similarity-threshold 0.05 $F > thr.txt 2>/dev/null
check "threshold 0.05 pairs" "$exp" "$(pairs thr.txt list)"
check "threshold 0.05 list sizes" "5 5 5 5 5 5" "$(grep -v '^#' thr.txt | awk -F'\t' '{print NF - 1}' | tr '\n' ' ' | sed 's/ $//')"
"$D" cmp $O --topk 5 $F > topk.txt 2>/dev/null
check "topk 5 pairs" "$exp" "$(pairs topk.txt list)"
check "topk 5 list sizes" "5 5 5 5 5 5" "$(grep -v '^#' topk.txt | awk -F'\t' '{print NF - 1}' | tr '\n' ' ' | sed 's/ $//')"
rm -f ./*.kmerset* ./*.kmercounts*
exit $fail
