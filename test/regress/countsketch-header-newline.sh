#!/usr/bin/env bash
# With -c/--countmin-size the #Dashing2Options header must stay one line:
# "counting=countsketchN" is followed by the remaining fields (";canon"), not
# by a newline that leaves a stray non-comment line before the matrix.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
check() { # name got expected
    if [ "$2" = "$3" ]; then echo "PASS $1: got $2 expected $3"; else echo "FAIL $1: got $2 expected $3"; fail=1; fi
}
python3 - <<'PY'
import random
r = random.Random(3)
for name in 'xy':
    s = ''.join(r.choice('ACGT') for _ in range(2000))
    open(name + '.fa', 'w').write(f'>{name}\n{s}\n')
PY
out=$("$D2" cmp -k 15 -B -c 64 x.fa y.fa 2>/dev/null)
rm -f ./*.kmer*
opts=$(grep '^#Dashing2Options' <<<"$out")
check "header has countsketch64;canon" "$(grep -c 'counting=countsketch64;canon' <<<"$opts")" 1
# Two inputs give exactly two lines that do not start with '#': the matrix rows.
check "non-comment lines" "$(grep -vc '^#' <<<"$out")" 2
exit $fail
