#!/usr/bin/env bash
# --downsample must apply in --parse-by-seq mode as in file mode. k-mers are
# sampled by a hash of the k-mer, so a record compared by sequence and the same
# record saved as its own file must give identical similarities, the union of
# two records must shrink to about half, and the -G minimizer stream must be
# about half as long.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
check() { # name got expected tolerance
    if python3 -c "import sys; sys.exit(not abs(float('$2') - float('$3')) <= $4)" 2>/dev/null; then
        echo "PASS $1: got $2 expected $3"
    else
        echo "FAIL $1: got $2 expected $3"; fail=1
    fi
}
# r1: 2 kbp random; r2: r1 with every 25th base changed; r3: r1's first half
# followed by 1 kbp of new sequence. r.fa holds all three records.
python3 - <<'PY'
import random
r = random.Random(5)
x = ''.join(r.choice('ACGT') for _ in range(2000))
y = ''.join('ACGT'[('ACGT'.index(c) + 1) % 4] if i % 25 == 12 else c for i, c in enumerate(x))
z = x[:1000] + ''.join(r.choice('ACGT') for _ in range(1000))
with open('r.fa', 'w') as f:
    for n, s in [('r1', x), ('r2', y), ('r3', z)]:
        f.write(f'>{n}\n{s}\n')
        open(n + '.fa', 'w').write(f'>{n}\n{s}\n')
PY
# Half the exact number of distinct canonical 21-mers in r1 and r3.
half_union=$(python3 - <<'PY'
def kmers(s, k=21):
    rc = s[::-1].translate(str.maketrans('ACGT', 'TGCA'))
    n = len(s)
    return {min(s[i:i + k], rc[n - i - k:n - i]) for i in range(n - k + 1)}
seq = lambda f: open(f).read().split('\n')[1]
print(len(kmers(seq('r1.fa')) | kmers(seq('r3.fa'))) / 2)
PY
)
vals() { grep -v '^#' | cut -f2- | tr '\n\t' '  '; }
d2() { "$D2" "$@" 2>/dev/null; }
byseq=$(d2 cmp -k 21 --full --downsample 0.5 --parse-by-seq r.fa | vals)
byfile=$(d2 cmp -k 21 --full --downsample 0.5 r1.fa r2.fa r3.fa | vals)
if [ "$byseq" = "$byfile" ]; then echo "PASS by-seq similarities equal file mode: $byseq"
else echo "FAIL by-seq similarities equal file mode: got $byseq expected $byfile"; fail=1; fi
u13=$(d2 cmp -k 21 --full --union-size --downsample 0.5 --parse-by-seq r.fa | awk '$1 == "r1" {print $4}')
check "by-seq union r1,r3 (+-10%)" "$u13" "$half_union" "$(python3 -c "print($half_union * 0.1)")"
# -G writes one 8-byte minimizer per entry after a header and per-record counts.
d2 sketch -k 21 -G --parse-by-seq r.fa -o full.seq
d2 sketch -k 21 -G --parse-by-seq --downsample 0.5 r.fa -o half.seq
size() { wc -c < "$1" | tr -d ' '; }
ratio=$(python3 -c "h = 20 + 8 * 3; print(($(size half.seq) - h) / ($(size full.seq) - h))")
check "-G by-seq kept fraction (+-0.05)" "$ratio" 0.5 0.05
exit $fail
