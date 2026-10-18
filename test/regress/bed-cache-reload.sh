#!/usr/bin/env bash
# A BED sketch reloaded from its --cache file must equal the freshly computed
# one, so a second identical `cmp --bed --cache` run must print the same matrix
# as the first. (Row labels are compared as printed; only the values matter.)
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
# Exact Jaccard: a~b = 500/1500, b~c = 500/6000, a~c = 0.
printf 'chr1\t0\t1000\n' > a.bed
printf 'chr1\t500\t1500\n' > b.bed
printf 'chr2\t0\t5000\nchr1\t1000\t1500\n' > c.bed
for opts in "--full" "--multiset" "--prob"; do
    rm -f a.bed.* b.bed.* c.bed.*
    first=$("$D2" cmp --bed $opts --cache a.bed b.bed c.bed 2>/dev/null | grep -v '^#' | tr -s ' \t\n' ' ')
    second=$("$D2" cmp --bed $opts --cache a.bed b.bed c.bed 2>/dev/null | grep -v '^#' | tr -s ' \t\n' ' ')
    if [ -n "$first" ] && [ "$first" = "$second" ]; then echo "PASS bed $opts cached rerun: got '$second' expected '$first'"
    else echo "FAIL bed $opts cached rerun: got '$second' expected '$first'"; fail=1; fi
done
rm -f a.bed.* b.bed.* c.bed.*
exit $fail
