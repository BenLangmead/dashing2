#!/usr/bin/env bash
# sketch --bed -o must write the stacked sketch file (it used to abort with
# "Failed to open file ... for in-place modification"), and comparing that file
# with --presketched must give the same matrix as cmp on the BED files.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
printf 'chr1\t0\t1000\n' > a.bed
printf 'chr1\t500\t1500\n' > b.bed
printf 'chr2\t0\t1000\nchr2\t2000\t2500\n' > c.bed
fail=0
"$D" sketch --bed --full -o stacked.ss a.bed b.bed c.bed >/dev/null 2>&1; rc=$?
if [ $rc -eq 0 ] && [ -s stacked.ss ]; then echo "PASS sketch --bed -o: exit $rc, file written"
else echo "FAIL sketch --bed -o: exit $rc, file $( [ -s stacked.ss ] && echo written || echo missing)"; fail=1; fi
direct=$("$D" cmp --bed --full a.bed b.bed c.bed 2>/dev/null | grep -v '^#')
loaded=$("$D" cmp --presketched stacked.ss 2>/dev/null | grep -v '^#')
if [ -n "$direct" ] && [ "$loaded" = "$direct" ]; then echo "PASS presketched matrix equals direct matrix"
else printf 'FAIL presketched matrix differs\n--- direct\n%s\n--- presketched\n%s\n' "$direct" "$loaded"; fail=1; fi
exit $fail
