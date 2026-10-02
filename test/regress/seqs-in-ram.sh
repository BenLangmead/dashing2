#!/usr/bin/env bash
# --seqs-in-ram keeps --parse-by-seq sequences in memory. The flag was set in a
# per-file copy of a static variable that the by-sequence sketcher never read, so
# it had no effect: the debug log still reported the switch it makes when the flag
# is off ("Swapping to keep sequences in RAM"). Inputs under 2 billion bases are
# kept in memory either way, so the results must not change.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
r = random.Random(1)
with open('s.fa', 'w') as f:
    for i in range(4):
        f.write('>s%d\n%s\n' % (i, ''.join(r.choice('ACGT') for _ in range(400))))
PY
fail=0
check() {  # check <name> <got> <expected>
    if [ "$2" = "$3" ]; then echo "PASS $1: got $2 expected $3"
    else echo "FAIL $1: got ${2:-nothing} expected $3"; fail=1; fi
}
a=(sketch -v -v -v --parse-by-seq -k 15 -S 128)
"$D" "${a[@]}" --seqs-in-ram --cmpout with.tbl s.fa 2>with.log
check "exit status with --seqs-in-ram" $? 0
check "log of --seqs-in-ram run mentions the swap" "$(grep -c 'Swapping to keep sequences in RAM' with.log)" 0
"$D" "${a[@]}" --cmpout without.tbl s.fa 2>without.log
check "log without --seqs-in-ram mentions the swap" "$(grep -c 'Swapping to keep sequences in RAM' without.log)" 1
check "same table with and without --seqs-in-ram" "$(cmp -s with.tbl without.tbl && echo same || echo differs)" same
exit $fail
