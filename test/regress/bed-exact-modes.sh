#!/usr/bin/env bash
# Exact modes (--set, --countdict) have no k-mer files to compare for BED input.
# cmp used to segfault; it must stop with exit status 1 and a message.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
printf 'chr1\t0\t1000\n' > a.bed
printf 'chr1\t500\t1500\n' > b.bed
fail=0
for mode in --set --countdict; do
    out=$("$D" cmp --bed $mode a.bed b.bed 2>&1); rc=$?
    if [ $rc -gt 0 ] && [ $rc -lt 128 ] && echo "$out" | grep -q 'support sketches only'; then got=rejected; else got="exit $rc"; fi
    if [ "$got" = rejected ]; then echo "PASS cmp --bed $mode: got $got expected rejected"
    else echo "FAIL cmp --bed $mode: got $got expected rejected"; fail=1; fi
done
exit $fail
