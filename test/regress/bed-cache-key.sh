#!/usr/bin/env bash
# A BED sketch cache written under one set of options must not be reused under
# another. The cache name lacked the sketch size and the other options that
# change a BED sketch, so `--cache -S 512` after a `-S 2048` run read the first
# 512 registers of the larger sketch, `--cache -S 2048` after a `-S 512` run
# aborted, and --normalize-intervals reused the unnormalized sketch. Each cached
# run must print the same values as an uncached run with the same options.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
printf 'chr1\t0\t1000\n' > a.bed
printf 'chr1\t500\t1500\nchr1\t3000\t3100\n' > b.bed
# Uncached reference values are computed on copies in fresh/, so that no cache file they
# might leave behind can be read by the runs below.
mkdir fresh; cp a.bed b.bed fresh/
# Reading a mis-sized cache could corrupt memory and hang, so runs are time-limited where possible.
TO=$(command -v timeout >/dev/null && echo "timeout 60")
vals() { $TO "$D2" cmp --bed "$@" 2>/dev/null | grep -v '^#' | grep '\.bed' | cut -f2- | tr -s ' \t\n' ' '; }
# Each case: options of the run that writes the cache, then options of the run that must not reuse it.
while IFS='|' read -r first second; do
    rm -f a.bed.* b.bed.*
    want=$(vals $second fresh/a.bed fresh/b.bed)
    vals $first --cache a.bed b.bed >/dev/null
    got=$(vals $second --cache a.bed b.bed)
    if [ -n "$want" ] && [ "$got" = "$want" ]; then echo "PASS cache from '$first' then '$second': got '$got' expected '$want'"
    else echo "FAIL cache from '$first' then '$second': got '$got' expected '$want'"; fail=1; fi
done <<'CASES'
--full -S 2048|--full -S 512
--full -S 512|--full -S 2048
--multiset -S 512|--multiset -S 512 --normalize-intervals
--multiset -S 512|--multiset -S 512 --countsketch-size 64
CASES
rm -f a.bed.* b.bed.*
exit $fail
