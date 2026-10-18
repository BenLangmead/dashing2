#!/usr/bin/env bash
# --seed changes k-mer hashing only. BED (and BigWig) positions are hashed from the
# chromosome name and coordinate with fixed hash functions, so --seed leaves their
# sketches unchanged. The help text must say so, and the BED distance matrix with
# --seed 7 must be identical to the one with the default seed.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
for sub in sketch cmp; do
    n=$("$D2" $sub --help 2>&1 | grep -c 'does not change --bed or --bigwig sketches')
    if [ "$n" -ge 1 ]; then echo "PASS $sub --help states that --seed does not apply to --bed/--bigwig: got $n line(s) expected >= 1"
    else echo "FAIL $sub --help states that --seed does not apply to --bed/--bigwig: got $n line(s) expected >= 1"; fail=1; fi
done
printf 'chr1\t0\t1000\nchr2\t0\t300\n' > a.bed
printf 'chr1\t500\t1500\n' > b.bed
for opts in "--full" "--multiset" "--prob"; do
    d0=$("$D2" cmp --bed $opts a.bed b.bed 2>/dev/null | grep -v '^#')
    d7=$("$D2" cmp --bed $opts --seed 7 a.bed b.bed 2>/dev/null | grep -v '^#')
    if [ -n "$d0" ] && [ "$d0" = "$d7" ]; then echo "PASS bed $opts --seed 7 matrix equals default"
    else echo "FAIL bed $opts --seed 7 matrix differs from default"; fail=1; fi
done
exit $fail
