#!/usr/bin/env bash
# With --outprefix, BED inputs with the same file name in different directories
# must not share a cache file. d1/a.bed and d2/a.bed are on different
# chromosomes; c.bed is a copy of d2/a.bed. After d1/a.bed is cached under op/,
# d2/a.bed against c.bed must be 1 (identical intervals). When both a.bed files
# were cached as op/a.bed.ss, d2/a.bed loaded d1/a.bed's sketch and gave 0.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
mkdir d1 d2
printf 'chr1\t0\t1000\n' > d1/a.bed
printf 'chr2\t0\t1000\n' > d2/a.bed
cp d2/a.bed c.bed
for opts in "--full" "--multiset"; do
    rm -rf op; mkdir op
    "$D2" cmp --bed $opts --cache --outprefix op d1/a.bed >/dev/null 2>&1
    got=$("$D2" cmp --bed $opts --cache --outprefix op d2/a.bed c.bed 2>/dev/null | awk -F'\t' '!/^#/ && NF > 2 {print $3; exit}')
    if [ "$got" = 1 ]; then echo "PASS bed $opts --outprefix d2/a~c after caching d1/a: got $got expected 1"
    else echo "FAIL bed $opts --outprefix d2/a~c after caching d1/a: got $got expected 1"; fail=1; fi
done
exit $fail
