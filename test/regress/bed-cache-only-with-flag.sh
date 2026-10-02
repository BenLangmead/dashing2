#!/usr/bin/env bash
# BED sketches must be written next to the inputs (or under --outprefix) only with
# --cache, as for sequence files. Without --cache, cmp and sketch wrote a.bed.ss,
# a.bed.bmh, ... on every run.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
check() { # name got expected
    if [ "$2" = "$3" ]; then echo "PASS $1: got $2 expected $3"; else echo "FAIL $1: got $2 expected $3"; fail=1; fi
}
printf 'chr1\t0\t1000\n' > a.bed
printf 'chr1\t500\t1500\n' > b.bed
mkdir op
# Files next to the inputs other than the inputs themselves, plus files under op/.
ncache() { { ls a.bed.* b.bed.* 2>/dev/null; ls -A op; } | wc -l | tr -d ' '; }
run() { ("$D2" "$@"; :) >/dev/null 2>&1; }
for opts in "--full" "--multiset" "--prob"; do
    rm -f a.bed.* b.bed.* op/*
    run cmp --bed $opts a.bed b.bed
    check "cmp --bed $opts without --cache writes no cache" "$(ncache)" 0
    run sketch --bed $opts -o stack a.bed b.bed
    run sketch --bed $opts --outprefix op -o stack a.bed b.bed
    check "sketch --bed $opts without --cache writes no cache" "$(ncache)" 0
    run cmp --bed $opts --cache a.bed b.bed
    check "cmp --bed $opts --cache writes one cache per input" "$(ncache)" 2
done
rm -f a.bed.* b.bed.* op/*
exit $fail
