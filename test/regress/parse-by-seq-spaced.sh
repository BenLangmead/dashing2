#!/bin/bash
# --parse-by-seq with --spacing must sketch the spaced k-mers of each record.
# With --parse-by-seq --full-setsketch, the diagonal of --union-size --square
# is the exact number of distinct k-mers in each record, so it must equal the
# count of distinct spaced k-mers computed below in python.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
DASHING2=${DASHING2:-$HERE/../../dashing2}
D=$(mktemp -d)
trap 'rm -rf "$D"' EXIT
cd "$D"

# Seed 0,1,0,2,0x7 with k = 12 keeps offsets 0,1,3,4,7,8,...,14. Records use
# only A and C so that k-mers repeat and the spaced count differs from the
# contiguous one.
python3 - > expected.txt <<'EOF'
import random
random.seed(5)
recs = ["".join(random.choice("AC") for _ in range(2000)) for _ in range(2)]
with open("recs.fa", "w") as f:
    for i, s in enumerate(recs): f.write(">r%d\n%s\n" % (i, s))
offs = [0]
for g in [0, 1, 0, 2] + [0] * 7: offs.append(offs[-1] + g + 1)
for s in recs:
    print(len({"".join(s[i + o] for o in offs) for i in range(len(s) - offs[-1])}))
EOF

"$DASHING2" cmp -k 12 --spacing 0,1,0,2,0x7 --parse-by-seq --full-setsketch --union-size --square recs.fa 2>/dev/null \
    | awk -F'\t' '!/^#/ {++n; print $(n + 1)}' | tr -d ' ' > got.txt
fail=0
for i in 1 2; do
    got=$(sed -n "${i}p" got.txt); exp=$(sed -n "${i}p" expected.txt)
    if [ -n "$exp" ] && [ "$got" = "$exp" ]; then echo "PASS spaced_cardinality_record$i: got $got expected $exp"
    else echo "FAIL spaced_cardinality_record$i: got $got expected $exp"; fail=1; fi
done
exit $fail
