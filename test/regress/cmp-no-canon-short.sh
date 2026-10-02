#!/usr/bin/env bash
# -C is the short form of --no-canon. sketch accepted it but cmp printed its usage
# and exited. A sequence and its reverse complement are identical when k-mers are
# canonicalized (similarity 1) and share almost no k-mers when they are not.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
r = random.Random(1)
a = ''.join(r.choice('ACGT') for _ in range(3000))
rc = a[::-1].translate(str.maketrans('ACGT', 'TGCA'))
open('a.fa', 'w').write('>a\n' + a + '\n')
open('rc.fa', 'w').write('>rc\n' + rc + '\n')
PY
fail=0
v() { "$D" cmp -k 21 -S 256 "$@" a.fa rc.fa 2>/dev/null | awk 'NR==4{print $3}'; }
check() {  # check <name> <got> <expected>
    if [ "$2" = "$3" ]; then echo "PASS $1: got $2 expected $3"
    else echo "FAIL $1: got ${2:-nothing} expected $3"; fail=1; fi
}
long=$(v --no-canon)
check "cmp default (canonical)" "$(v)" 1
check "cmp --no-canon is strand-specific" "$(awk -v x="$long" 'BEGIN{print (x != "" && x < 0.1) ? "yes" : "no"}')" yes
check "cmp -C equals cmp --no-canon" "$(v -C)" "$long"
"$D" cmp -k 21 -S 256 -C a.fa rc.fa >/dev/null 2>&1
check "cmp -C exit status" $? 0
exit $fail
