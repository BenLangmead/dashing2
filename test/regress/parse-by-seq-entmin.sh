#!/bin/bash
# --parse-by-seq --entmin must select the same canonical minimizers as sketching
# one file per record. Checks: a record and its reverse complement have J = 1,
# and every by-record similarity equals the by-file similarity.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
DASHING2=${DASHING2:-$HERE/../../dashing2}
D=$(mktemp -d)
trap 'rm -rf "$D"' EXIT
cd "$D"

python3 - <<'EOF'
import random
random.seed(9)
x = "".join(random.choice("ACGT") for _ in range(3000))
y = x[:1500] + "".join(random.choice("ACGT") for _ in range(1500))
rc = x[::-1].translate(str.maketrans("ACGT", "TGCA"))
recs = [("x", x), ("xrc", rc), ("y", y)]
open("all.fa", "w").write("".join(">%s\n%s\n" % r for r in recs))
for name, s in recs: open(name + ".fa", "w").write(">%s\n%s\n" % (name, s))
EOF

flags=(--full-setsketch --entmin -k 21 -w 27)
"$DASHING2" cmp "${flags[@]}" --parse-by-seq all.fa 2>/dev/null | awk -F'\t' '!/^#/' > byseq.txt
"$DASHING2" cmp "${flags[@]}" x.fa xrc.fa y.fa 2>/dev/null | awk -F'\t' '!/^#/' > byfile.txt
cell() { awk -F'\t' -v r="$2" -v c="$3" 'NR == r {print $(c + 1)}' "$1" | tr -d ' '; }
fail=0
check() { # name got expected
    if [ -n "$3" ] && [ "$2" = "$3" ]; then echo "PASS $1: got $2 expected $3"
    else echo "FAIL $1: got $2 expected $3"; fail=1; fi
}
check byseq_x_vs_revcomp "$(cell byseq.txt 1 2)" 1
check byseq_equals_byfile_x_xrc "$(cell byseq.txt 1 2)" "$(cell byfile.txt 1 2)"
check byseq_equals_byfile_x_y "$(cell byseq.txt 1 3)" "$(cell byfile.txt 1 3)"
check byseq_equals_byfile_xrc_y "$(cell byseq.txt 2 3)" "$(cell byfile.txt 2 3)"
exit $fail
