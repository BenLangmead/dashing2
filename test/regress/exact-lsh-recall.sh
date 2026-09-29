#!/usr/bin/env bash
# LSH-assisted modes (--topk, --similarity-threshold) with exact sketches
# (--set, --countdict) must find the neighbors that the exhaustive --square
# matrix shows. Fixtures: 30 random 5 kbp genomes in 5 clusters of 6 (member j
# of cluster c has 0.2% * j * (c + 1) substitutions from member 0). Checks: every
# --topk 1 row names one of its exact nearest neighbors (ties allowed), and
# --similarity-threshold 0.3 lists exactly the off-diagonal --square entries
# that are >= 0.3.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
r = random.Random(7)
sub = {'A': 'CGT', 'C': 'AGT', 'G': 'ACT', 'T': 'ACG'}
for c in range(5):
    base = ''.join(r.choice('ACGT') for _ in range(5000))
    for j in range(6):
        s = list(base)
        for p in r.sample(range(5000), int(5000 * 0.002 * j * (c + 1))):
            s[p] = r.choice(sub[s[p]])
        open('c%d_%d.fa' % (c, j), 'w').write('>c%d_%d\n' % (c, j) + ''.join(s) + '\n')
PY
cat > check.py <<'PY'
# Prints: rows whose top-1 neighbor is not an exact nearest neighbor,
# thresholded entries emitted, thresholded entries expected.
rows = [l.rstrip('\n').split('\t') for l in open('sq.txt') if not l.startswith('#')]
names = [r[0].strip() for r in rows]
best = {}
for i, r in enumerate(rows):
    cand = [(float(v), names[j]) for j, v in enumerate(r[1:]) if j != i]
    m = max(v for v, _ in cand)
    best[names[i]] = set(n for v, n in cand if v >= m - 1e-9)
bad, seen = 0, set()
for l in open('tk.txt'):
    if l.startswith('#'): continue
    f = l.rstrip('\n').split('\t'); q = f[0].strip(); seen.add(q)
    nb = [x.rsplit(':', 1)[0] for x in f[1:] if x]
    bad += not nb or nb[0] not in best[q]
bad += len(set(names) - seen)
want = sum(1 for i, r in enumerate(rows) for j, v in enumerate(r[1:]) if i != j and float(v) >= 0.3)
got = sum(len([x for x in l.rstrip('\n').split('\t')[1:] if x]) for l in open('th.txt') if not l.startswith('#'))
print(bad, got, want)
PY
fail=0
for mode in --set --countdict; do
    "$D" cmp $mode -k 21 --square c*.fa > sq.txt 2>/dev/null
    "$D" cmp $mode -k 21 --topk 1 c*.fa > tk.txt 2>/dev/null
    "$D" cmp $mode -k 21 --similarity-threshold 0.3 c*.fa > th.txt 2>/dev/null
    rm -f ./*.kmerset* ./*.kmercounts*
    read -r bad got want < <(python3 check.py)
    if [ "$bad" -eq 0 ]; then echo "PASS $mode --topk 1: got $bad rows missing their nearest neighbor expected 0"
    else echo "FAIL $mode --topk 1: got $bad rows missing their nearest neighbor expected 0"; fail=1; fi
    if [ "$got" -eq "$want" ]; then echo "PASS $mode --similarity-threshold 0.3: got $got entries expected $want"
    else echo "FAIL $mode --similarity-threshold 0.3: got $got entries expected $want"; fail=1; fi
done
exit $fail
