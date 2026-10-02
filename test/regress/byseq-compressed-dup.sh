#!/usr/bin/env bash
# --parse-by-seq reads .xz, .bz2 and .zst inputs through a pipe. zlib must get
# its own copy of the pipe's descriptor: if gzclose and pclose both close the
# same number, the second close can hit a descriptor another thread has just
# opened. A joint entry of many compressed files read with -p 8 must give the
# same output as the plain files on every run, and (where strace is installed)
# no close() may fail with EBADF, the sign of a descriptor closed twice.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
check() { # name got expected
    if [ "$2" = "$3" ]; then echo "PASS $1: got $2 expected $3"; else echo "FAIL $1: got $2 expected $3"; fail=1; fi
}
# Twelve files of three related 1 kbp records each.
python3 - <<'PY'
import random
r = random.Random(7)
for f in range(12):
    x = ''.join(r.choice('ACGT') for _ in range(1000))
    with open('f%02d.fa' % f, 'w') as o:
        for j in range(3):
            y = ''.join(r.choice('ACGT') if r.random() < 0.02 * j else c for c in x)
            o.write('>f%02d_r%d\n%s\n' % (f, j, y))
PY
plain=""; comp=""; tools=(xz bzip2 zstd); sfx=(.xz .bz2 .zst)
for f in $(seq -w 0 11); do
    t=$(( 10#$f % 3 ))
    command -v ${tools[$t]} >/dev/null || { echo "SKIP: ${tools[$t]} not installed"; exit 0; }
    ${tools[$t]} -c f$f.fa > f$f.fa${sfx[$t]}
    plain="$plain f$f.fa"; comp="$comp f$f.fa${sfx[$t]}"
done
run() { "$D2" cmp -p 8 -k 21 --full --parse-by-seq "$1" 2>/dev/null | grep -v '^#Calling'; }
want=$(run "${plain# }" | cksum)
check "plain records" "$(run "${plain# }" | grep -c '^f')" 36
same=0
for i in 1 2 3 4 5; do [ "$(run "${comp# }" | cksum)" = "$want" ] && same=$((same + 1)); done
check "compressed runs equal to plain" "$same/5" "5/5"
if command -v strace >/dev/null; then
    strace -f -qq -e trace=close -o trace.txt "$D2" cmp -p 8 -k 21 --full --parse-by-seq "${comp# }" >/dev/null 2>&1
    check "close() calls failing with EBADF" "$(grep -c 'EBADF' trace.txt)" 0
else
    echo "SKIP EBADF check: strace not installed"
fi
exit $fail
