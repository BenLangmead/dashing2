#!/usr/bin/env bash
# With --outprefix, the BED sketch cache file must carry the sketch suffix
# (a.bed.ss, a.bed.bmh, ...) like the cache written without it. Named after the
# input alone, it replaced a.bed with binary sketch data when --outprefix was the
# input's own directory, and with --cache the BED text itself was read back as a
# cached sketch.
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
cp a.bed a.orig; cp b.bed b.orig
vals() { grep -v '^#' | cut -f2- | tr -s ' \t\n' ' '; }
for opts in "--full" "--multiset"; do
    fresh=$("$D2" cmp --bed $opts a.bed b.bed 2>/dev/null | vals)
    rm -f a.bed.* b.bed.*
    # The input's own directory as the output prefix must leave the inputs alone.
    "$D2" sketch --bed $opts --outprefix . a.bed b.bed >/dev/null 2>&1
    check "bed $opts --outprefix . keeps a.bed" "$(cmp -s a.bed a.orig && echo unchanged || echo overwritten)" unchanged
    check "bed $opts --outprefix . keeps b.bed" "$(cmp -s b.bed b.orig && echo unchanged || echo overwritten)" unchanged
    cp a.orig a.bed; cp b.orig b.bed; rm -f a.bed.* b.bed.*
    # A first --cache run must sketch the BED rows, not read the BED file as a cache.
    got=$("$D2" cmp --bed $opts --cache --outprefix . a.bed b.bed 2>/dev/null | vals)
    check "bed $opts --cache --outprefix . equals uncached" "'$got'" "'$fresh'"
    cp a.orig a.bed; cp b.orig b.bed; rm -f a.bed.* b.bed.*
    # Another prefix directory must not receive a file named like the input.
    mkdir -p out; rm -f out/*
    "$D2" sketch --bed $opts --outprefix out a.bed b.bed >/dev/null 2>&1
    check "bed $opts --outprefix out has no out/a.bed" "$( [ -e out/a.bed ] && echo present || echo absent)" absent
done
exit $fail
