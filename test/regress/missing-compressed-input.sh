#!/usr/bin/env bash
# A missing input is an error whatever its suffix. .xz, .bz2 and .zst inputs
# are read through a decompression pipe; before the fix a missing one only made
# the decompressor print a message, and dashing2 exited 0 with an empty sketch
# (file mode) or no records (--parse-by-seq), while a missing plain or .gz
# input was an error. Existing compressed inputs must still work.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
python3 - <<'PY'
import random
r = random.Random(6)
open('a.fa', 'w').write('>a\n%s\n' % ''.join(r.choice('ACGT') for _ in range(2000)))
PY
cp a.fa c.fa
have=""
command -v xz >/dev/null && xz -kf c.fa && have="$have .xz"
command -v bzip2 >/dev/null && bzip2 -kf c.fa && have="$have .bz2"
command -v zstd >/dev/null && zstd -qf c.fa && have="$have .zst"
for sfx in "" .gz .xz .bz2 .zst; do
    for mode in file --parse-by-seq; do
        if [ $mode = file ]; then args="a.fa missing.fa$sfx"; else args="missing.fa$sfx"; fi
        rc=$( ("$D2" sketch -k 15 ${mode#file} --cmpout out.txt $args >/dev/null 2>&1; echo $?) 2>/dev/null )
        if [ $rc -ne 0 ]; then echo "PASS missing.fa$sfx ($mode): exit status $rc"
        else echo "FAIL missing.fa$sfx ($mode): exit status 0, expected an error"; fail=1; fi
    done
done
# Self-similarity of a compressed copy of a.fa must be 1 (file mode).
for sfx in $have; do
    got=$("$D2" cmp -k 15 --square a.fa c.fa$sfx 2>/dev/null | awk '!/^#/ {print $3; exit}')
    if [ "$got" = 1 ]; then echo "PASS a.fa vs c.fa$sfx: got $got expected 1"
    else echo "FAIL a.fa vs c.fa$sfx: got $got expected 1"; fail=1; fi
done
rm -f ./*.kmer* ./*.opss
exit $fail
