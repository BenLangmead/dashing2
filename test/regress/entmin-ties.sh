#!/usr/bin/env bash
# Canonical --entmin windows must not depend on the strand: a set of records
# and their reverse complements must give the same minimizers (similarity 1).
# In low-complexity sequence many k-mers share an entropy score, so the window
# minimum depends on the tie-break by k-mer value. A tie-break that compares a
# k-mer with itself picks the most recent of the tied k-mers instead, which
# differs between the two strands. k 25 (64-bit encoder) is a control.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 -c "
import random
r = random.Random(7)
comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
rnd = lambda n: ''.join(r.choice('ACGT') for _ in range(n))
# Homopolymer runs and tandem repeats fill windows with distinct k-mers that
# have equal entropy scores (for example A^59C and CA^59).
recs = []
for _ in range(6):
    unit = rnd(r.randint(30, 70))
    run = r.choice('ACGT') * r.randint(80, 120)
    recs.append(rnd(200) + run + rnd(100) + unit * r.randint(4, 8) + rnd(200))
with open('f.fa', 'w') as f, open('rc.fa', 'w') as g:
    for i, s in enumerate(recs):
        f.write('>r%d\\n%s\\n' % (i, s))
        g.write('>r%d\\n%s\\n' % (i, ''.join(comp[c] for c in reversed(s))))
"
fail=0
check() {
    name=$1; shift
    got=$("$D" cmp --set --entmin "$@" f.fa rc.fa 2>/dev/null | awk '$1 == "f.fa" {print $NF}')
    if awk -v g="$got" 'BEGIN {exit !(g + 0 == 1)}'; then echo "PASS $name: got $got expected 1"
    else echo "FAIL $name: got $got expected 1"; fail=1; fi
}
check "k25 w40" -k 25 -w 40
check "-2 k45 w49" -2 -k 45 -w 49
check "-2 k60 w65" -2 -k 60 -w 65
exit $fail
