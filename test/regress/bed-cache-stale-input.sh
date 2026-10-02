#!/usr/bin/env bash
# With --cache, a BED cache file must not be reused once its input has been
# overwritten. a.bed and b.bed are on different chromosomes, so a~b is 0, and
# that run caches both. b.bed is then overwritten with a copy of a.bed, and c.bed
# is another copy of a.bed that was never cached. b~c must be 1 (identical
# intervals), as without --cache; when the old cache of b.bed was reused it was 0.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
cell() { awk -F'\t' '!/^#/ && NF > 2 {print $3; exit}'; }
for opts in "--full" "--multiset" "--prob"; do
    rm -f ./*.bed*
    printf 'chr1\t0\t1000\n' > a.bed
    printf 'chr2\t0\t1000\n' > b.bed
    first=$("$D2" cmp --bed $opts --cache a.bed b.bed 2>/dev/null | cell)
    cp a.bed b.bed; cp a.bed c.bed
    got=$("$D2" cmp --bed $opts --cache b.bed c.bed 2>/dev/null | cell)
    if [ "$first" = 0 ] && [ "$got" = 1 ]; then echo "PASS bed $opts --cache b~c after overwriting b: got $got expected 1 (a~b before: $first)"
    else echo "FAIL bed $opts --cache b~c after overwriting b: got $got expected 1 (a~b before: $first, expected 0)"; fail=1; fi
done
exit $fail
