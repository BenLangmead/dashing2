#!/usr/bin/env bash
# -Q/--qfile asks for references (-F) against queries (-Q). --topk,
# --similarity-threshold and --greedy used to run over the union of both sets
# (reference-to-reference neighbors, clusters of references), or, with -Q given
# last, silently print the panel and ignore the flag. Acceptable behavior: an
# error naming -Q, or output in which no line relates two references.
# Fixtures: r1/r2 are near-identical references (10 substitutions in 5 kbp),
# q1 is a query close to r1, q2 an unrelated query.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
r = random.Random(5)
sub = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}
a = ''.join(r.choice('ACGT') for _ in range(5000))
def mut(s, step, off):
    s = list(s)
    for i in range(off, len(s), step):
        s[i] = sub[s[i]]
    return ''.join(s)
open('r1.fa', 'w').write('>r1\n' + a + '\n')
open('r2.fa', 'w').write('>r2\n' + mut(a, 500, 0) + '\n')
open('q1.fa', 'w').write('>q1\n' + mut(a, 500, 250) + '\n')
open('q2.fa', 'w').write('>q2\n' + ''.join(r.choice('ACGT') for _ in range(5000)) + '\n')
PY
printf 'r1.fa\nr2.fa\n' > F.txt; printf 'q1.fa\nq2.fa\n' > Q.txt
fail=0
for args in "-F F.txt -Q Q.txt --topk 1" "--topk 1 -F F.txt -Q Q.txt" "-F F.txt -Q Q.txt --similarity-threshold 0.5" "-F F.txt -Q Q.txt --greedy 0.5"; do
    out=$("$D" cmp --full -k 21 $args 2>err.txt); rc=$?
    if [ $rc -ne 0 ]; then
        if grep -q -- '-Q' err.txt; then r="rejected naming -Q"; ok=1; else r="exit $rc without a clear error"; ok=0; fi
    else
        n=$(echo "$out" | grep -v '^#' | grep 'r1.fa' | grep -c 'r2.fa')
        r="$n lines relating r1 and r2"; ok=$([ "$n" -eq 0 ] && echo 1 || echo 0)
        if echo "$out" | grep -q '^#Dashing2 Panel'; then r="a full panel (flag ignored)"; ok=0; fi
    fi
    if [ $ok -eq 1 ]; then echo "PASS $args: got $r expected rejection or neighbor/cluster output without reference-reference lines"
    else echo "FAIL $args: got $r expected rejection or neighbor/cluster output without reference-reference lines"; fail=1; fi
done
exit $fail
