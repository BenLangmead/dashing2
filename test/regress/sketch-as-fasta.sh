#!/usr/bin/env bash
# Sketch files passed to cmp without --presketched used to be parsed as FASTA,
# so two unrelated inputs' sketches compared as similarity 1. They must be
# rejected with exit status 1 and a message; a real FASTA file that happens to
# use such an extension must still be read (a file against its copy gives 1).
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
r = random.Random(1)
for name in 'ab':
    open(name + '.fa', 'w').write('>' + name + '\n' + ''.join(r.choice('ACGT') for _ in range(5000)) + '\n')
PY
fail=0
report() {
    if [ "$2" = "$3" ]; then echo "PASS $1: got $2 expected $3"
    else echo "FAIL $1: got $2 expected $3"; fail=1; fi
}
"$D" sketch --full -k 21 --cache a.fa b.fa >/dev/null 2>&1
out=$("$D" cmp --full -k 21 a.fa.*.ss b.fa.*.ss 2>&1); rc=$?
if [ $rc -gt 0 ] && [ $rc -lt 128 ] && echo "$out" | grep -q 'looks like a dashing2 sketch'; then got=rejected
else got="accepted, similarity $(echo "$out" | awk '$1 ~ /^a\.fa/ && NF == 3 {print $3}')"; fi
report "unrelated .ss files without --presketched" "$got" rejected
cp a.fa a-copy.ss
report "FASTA file named .ss" "$("$D" cmp --full -k 21 a.fa a-copy.ss 2>/dev/null | awk 'NR==4{print $3}')" 1
exit $fail
