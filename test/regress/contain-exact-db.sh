#!/usr/bin/env bash
# sketch --set/--countdict --save-kmers -o db writes db.kmer64, but exact modes
# have no per-register sampled k-mers to put there, so contain found nothing
# (0% self coverage). Acceptable behavior: either sketch refuses with an error
# naming --save-kmers and a normal nonzero exit status (not a signal such as
# SIGABRT), or the database works (x.fa against itself is 100%).
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
r = random.Random(1)
for name in 'xy':
    s = ''.join(r.choice('ACGT') for _ in range(20000))
    open(name + '.fa', 'w').write('>' + name + '\n' + s + '\n')
PY
fail=0
for mode in --set --countdict; do
    err=$("$D" sketch $mode -k 21 --save-kmers -o db x.fa y.fa 2>&1 >/dev/null); rc=$?
    if [ $rc -ge 128 ]; then
        echo "FAIL $mode: got exit $rc (killed by a signal) expected rejection or 100% self coverage"; fail=1
    elif [ $rc -ne 0 ]; then
        if echo "$err" | grep -q -- '--save-kmers'; then
            echo "PASS $mode: got rejected with an error naming --save-kmers expected rejection or 100% self coverage"
        else
            echo "FAIL $mode: got exit $rc without a clear error expected rejection or 100% self coverage"; fail=1
        fi
    else
        got=$("$D" contain db.kmer64 x.fa 2>/dev/null | awk '$1=="x.fa"{print $2}')
        case $got in
            100%*) echo "PASS $mode: got self coverage $got expected rejection or 100% self coverage" ;;
            *) echo "FAIL $mode: got self coverage $got expected rejection or 100% self coverage"; fail=1 ;;
        esac
    fi
    rm -f db db.* x.fa.* y.fa.*
done
exit $fail
