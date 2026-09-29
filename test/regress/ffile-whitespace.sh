#!/bin/bash
# -F path lists: a line "a.fa b.fa" sketches the union of both files as one
# entry. Extra spaces, a trailing space, a CRLF line ending or blank lines must
# not change what is read. Each check sketches with --set -k 21 and compares
# the cardinality of the entry for a.fa with the exact number of distinct
# canonical 21-mers, computed below.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
DASHING2=${DASHING2:-$HERE/../../dashing2}
D=$(mktemp -d)
trap 'rm -rf "$D"' EXIT
cd "$D"

python3 - <<'EOF'
import random
random.seed(24)
a = "".join(random.choice("ACGT") for _ in range(500))
seqs = {"a": a, "b": a[:250] + "".join(random.choice("ACGT") for _ in range(250))}
for n, s in seqs.items():
    open(n + ".fa", "w").write(">%s\n%s\n" % (n, s))
def kmers(s, k=21):
    rc = lambda w: w[::-1].translate(str.maketrans("ACGT", "TGCA"))
    return {min(s[i:i + k], rc(s[i:i + k])) for i in range(len(s) - k + 1)}
open("exact_a", "w").write(str(len(kmers(seqs["a"]))))
open("exact_ab", "w").write(str(len(kmers(seqs["a"]) | kmers(seqs["b"]))))
EOF

fail=0
card() { # cardinality stored in the kmerset file of the entry named after $1
    local f
    f=$(ls "$1".*.kmerset64 2>/dev/null | head -1)
    if [ -n "$f" ]; then python3 -c "import struct,sys; print(int(struct.unpack('<d', open(sys.argv[1], 'rb').read(8))[0]))" "$f"
    else echo missing; fi
}
check() { # name list-contents expected-card-of-a.fa-entry
    local name=$1 exp=$3 got rc
    printf "$2" > list.txt
    rm -f ./*.kmerset64
    "$DASHING2" sketch --set -k 21 -F list.txt >/dev/null 2>&1
    rc=$?
    got="rc=$rc card=$(card a.fa)"
    rm -f ./*.kmerset64
    if [ "$got" = "rc=0 card=$exp" ]; then echo "PASS $name: got $got expected rc=0 card=$exp"
    else echo "FAIL $name: got $got expected rc=0 card=$exp"; fail=1; fi
}
AB=$(cat exact_ab); A=$(cat exact_a)
check joint_one_space       'a.fa b.fa\n'     "$AB"
check joint_two_spaces      'a.fa  b.fa\n'    "$AB"
check joint_trailing_space  'a.fa b.fa \n'    "$AB"
check joint_leading_space   '  a.fa b.fa\n'   "$AB"
check joint_crlf            'a.fa b.fa\r\n'   "$AB"
check separate_blank_line   'a.fa\n\nb.fa\n'  "$A"
check separate_crlf         'a.fa\r\nb.fa\r\n' "$A"

# cmp reads -F the same way: two entries, and J equal to that of positional arguments.
printf 'a.fa\r\n\nb.fa \r\n' > list.txt
exp=$("$DASHING2" cmp --set -k 21 a.fa b.fa 2>/dev/null | awk -F'\t' '!/^#/ {print $3; exit}')
got=$("$DASHING2" cmp --set -k 21 -F list.txt 2>/dev/null | awk -F'\t' '!/^#/ {print $3; exit}')
rm -f ./*.kmerset64
if [ -n "$exp" ] && [ "$got" = "$exp" ]; then echo "PASS cmp_blank_line_crlf: got J=$got expected J=$exp"
else echo "FAIL cmp_blank_line_crlf: got J=$got expected J=$exp"; fail=1; fi
exit $fail
