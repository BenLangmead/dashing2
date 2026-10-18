#!/usr/bin/env bash
# -k was parsed with atoi, so out-of-range values wrapped around (4294967327
# became k 31 and ran without complaint) and non-numbers became 0. Invariant:
# -k accepts integers from 1 to 16777216 and rejects anything else with exit
# status 1 and a message naming the option.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 -c "import random; r=random.Random(9); open('x.fa','w').write('>x\n' + ''.join(r.choice('ACGT') for _ in range(2000)) + '\n')"
fail=0
for k in 4294967327 4294967296 -5 0 abc 21x 16777217; do
    # A binary that wraps k to 0 sizes its spaced-seed table by k - 1 bases;
    # the address-space limit keeps that from allocating gigabytes.
    (ulimit -v 2097152; timeout 60 "$D" sketch -k "$k" x.fa >/dev/null 2>err.txt); rc=$?
    if [ $rc -eq 1 ] && grep -q "kmer-length" err.txt; then
        echo "PASS -k $k rejected"
    else
        echo "FAIL -k $k: exit $rc, $(head -c 100 err.txt | tr '\n' ' ')"; fail=1
    fi
done
for k in 1 21 1000000; do
    if timeout 60 "$D" sketch -k "$k" x.fa >/dev/null 2>err.txt; then
        echo "PASS -k $k accepted"
    else
        echo "FAIL -k $k rejected: $(head -c 100 err.txt | tr '\n' ' ')"; fail=1
    fi
done
exit $fail
