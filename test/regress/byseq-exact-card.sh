#!/bin/bash
# Sets of at most 10 x the sketch size get an exact cardinality instead of the
# sketch estimate, in file mode and with --parse-by-seq alike, so a record
# compared with --parse-by-seq gives the same intersection, union and
# containments as the same record in a file of its own (by-seq mode used to
# count small sets exactly while file mode kept the estimate). Fixture: four
# records; the exact numbers of distinct canonical 21-mers are computed below
# and must equal the cardinalities dashing2 sketch writes for small sets.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
DASHING2=${DASHING2:-$HERE/../../dashing2}
D=$(mktemp -d)
trap 'rm -rf "$D"' EXIT
cd "$D"

python3 - <<'EOF' > exact.txt
import random
R = random.Random(23)
rnd = lambda n: "".join(R.choice("ACGT") for _ in range(n))
a = rnd(1200)
recs = {"r1": a, "r2": a[:800] + rnd(600), "r3": rnd(100), "r4": rnd(15000)}
open("all.fa", "w").write("".join(">%s\n%s\n" % kv for kv in recs.items()))
for n, s in recs.items():
    open(n + ".fa", "w").write(">%s\n%s\n" % (n, s))
    rc = s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    L = len(s)
    print(n + ".fa", len({min(s[i:i + 21], rc[L - i - 21:L - i]) for i in range(L - 20)}))
EOF

fail=0
# Cardinalities of the small sets are exact in file mode (r4, above 10 x 1024 k-mers, is estimated).
"$DASHING2" sketch -k 21 -o out r1.fa r2.fa r3.fa 2>/dev/null
while read -r name exact; do
    [ "$name" = r4.fa ] && continue
    got=$(awk -v f="$name" '$1 == f {print $2}' out.names.txt)
    if [ "$got" = "$exact" ]; then echo "PASS card_$name: got $got expected $exact"
    else echo "FAIL card_$name: got $got expected $exact"; fail=1; fi
done < exact.txt
# Every cell of --parse-by-seq on all.fa equals the same cell for one file per record.
for sk in "" --full; do
    for m in --intersection --union-size --containment --symmetric-containment; do
        "$DASHING2" cmp -k 21 $sk $m --square --parse-by-seq all.fa 2>/dev/null | grep -v '^#' > byseq.txt
        "$DASHING2" cmp -k 21 $sk $m --square r1.fa r2.fa r3.fa r4.fa 2>/dev/null | grep -v '^#' > files.txt
        res=$(python3 - <<'EOF'
rows = lambda f: [l.split("\t")[1:] for l in open(f).read().splitlines()]
a, b = rows("byseq.txt"), rows("files.txt")
bad = [(x, y) for ra, rb in zip(a, b) for x, y in zip(ra, rb) if abs(float(x) - float(y)) > 1e-6 * max(1., abs(float(y)))]
print(len(bad), len(a) * len(a[0]) if a else 0, *(bad[0] if bad else ()))
EOF
)
        set -- $res
        name="byseq_vs_file${sk:+_${sk#--}}_${m#--}"
        if [ "$1" = 0 ] && [ "$2" = 16 ]; then echo "PASS $name: got 0 of $2 cells differing expected 0"
        else echo "FAIL $name: got $1 of $2 cells differing expected 0 (first: by-seq ${3:-} file ${4:-})"; fail=1; fi
    done
done
rm -f ./*.kmerset* ./*.opss ./out*
exit $fail
