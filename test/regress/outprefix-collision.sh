#!/usr/bin/env bash
# --outprefix must keep inputs with the same file name in different
# directories apart. d1/a.fa and d1/z.fa are identical (expected 1); d2/a.fa
# is unrelated random sequence (expected 0 at k 21).
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
mkdir d1 d2 op
python3 - <<'PY'
import random
r = random.Random(1)
s1 = ''.join(r.choice('ACGT') for _ in range(5000))
s2 = ''.join(r.choice('ACGT') for _ in range(5000))
open('d1/a.fa', 'w').write('>x\n' + s1 + '\n')
open('d1/z.fa', 'w').write('>x\n' + s1 + '\n')
open('d2/a.fa', 'w').write('>y\n' + s2 + '\n')
PY
fail=0
for mode in --full --set; do
    out=$("$D" cmp $mode -k 21 --cache --outprefix op d1/a.fa d2/a.fa d1/z.fa 2>/dev/null)
    rm -f op/*
    # Row for d1/a.fa: columns are d1/a.fa, d2/a.fa, d1/z.fa
    read -r a_d2a a_z <<<"$(echo "$out" | awk 'NR==4{print $3, $4}')"
    for chk in "d1/a~d2/a ${a_d2a:-none} 0" "d1/a~d1/z ${a_z:-none} 1"; do
        set -- $chk
        if [ "$2" = "$3" ]; then echo "PASS $mode $1: got $2 expected $3"
        else echo "FAIL $mode $1: got $2 expected $3"; fail=1; fi
    done
done
exit $fail
