#!/usr/bin/env bash
# --parse-by-seq has no per-record exact k-mer sets or count dictionaries, and
# per-record minimizer sequences (-G) cannot be compared. These combinations
# must exit with status 1 and an error message instead of crashing, while the
# supported ones keep working.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
python3 - <<'PY'
import random
r = random.Random(5)
with open('two.fa', 'w') as f:
    for i in range(2):
        f.write(f'>r{i}\n' + ''.join(r.choice('ACGT') for _ in range(300)) + '\n')
PY
expect() { # expected-status description args...
    local want=$1 name=$2; shift 2
    "$D2" "$@" -k 21 --parse-by-seq two.fa >/dev/null 2>err.txt; local rc=$?
    rm -f ./*.kmer* o.* c.txt
    local msg=""
    [ "$want" = 1 ] && ! grep -q '^Error: --parse-by-seq' err.txt && msg=" (no error message)"
    if [ $rc -eq "$want" ] && [ -z "$msg" ]; then echo "PASS $name: got status $rc expected $want"
    else echo "FAIL $name: got status $rc$msg expected $want"; fail=1; fi
}
expect 1 "cmp --set" cmp --set
expect 1 "cmp --countdict" cmp --countdict
expect 1 "cmp -G -o" cmp -G -o o.seq
expect 1 "sketch --set -o" sketch --set -o o.out
expect 1 "sketch --set --cmpout" sketch --set --cmpout c.txt
expect 1 "sketch -G -o --cmpout" sketch -G -o o.seq --cmpout c.txt
expect 0 "sketch -G -o (supported)" sketch -G -o o.seq
expect 0 "cmp --full (supported)" cmp --full
expect 0 "cmp default (supported)" cmp
exit $fail
