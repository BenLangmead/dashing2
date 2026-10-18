#!/usr/bin/env bash
# dashing2 2.1.20 hashed k-mers above the exact limit (k > 32 for DNA) with a
# cyclic rolling hash, while this version uses a polynomial hash modulo
# 2^61-1, and -2 now folds 128-bit k-mers instead of truncating them. Cache
# names used to omit the scheme, so a 2.1.20 cache was reused at k > 32 and a
# comparison against a freshly sketched file reported similarity 0 (0.371
# expected at k 40, 0.224 at k 70). Invariant: with --cache, a cache written
# under the old (untagged) name does not change the result; cmp --cache equals
# cmp without a cache.
#
# The stale cache comes from DASHING2_OLD (a dashing2 2.1.20 binary) when it is
# set. Otherwise it is simulated by the binary under test: the sketch of an
# unrelated sequence is placed under the name the old scheme would have used
# (this binary's cache name without the scheme suffix).
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
OLD=${DASHING2_OLD:-}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
r = random.Random(3)
g = "".join(r.choice("ACGT") for _ in range(100000))
h = "".join(r.choice("ACGT") if r.random() < 0.02 else c for c in g)
z = "".join(r.choice("ACGT") for _ in range(100000))
open("x.fa", "w").write(">x\n" + g + "\n")
open("y.fa", "w").write(">y\n" + h + "\n")
open("z.fa", "w").write(">z\n" + z + "\n")
PY
sim() { "$D" cmp "$@" 2>/dev/null | awk -F'\t' '!/^#/{print $3; exit}'; }
# Prints the cache file that `sketch --cache` writes for $1 with options $2...
cachefile() {
    local f=$1; shift
    rm -rf c; mkdir c; cp "$f" c/
    "$1" sketch "${@:2}" --cache "c/$f" >/dev/null 2>&1
    ls c | grep -v "^$f\$" | head -1
}
fail=0
for args in "-k 40" "-k 70" "-2 -k 40"; do
    fresh=$(sim $args x.fa y.fa)
    new=$(cachefile x.fa "$D" $args)
    old=${new/.ph61/}
    rm -rf s; mkdir s; cp x.fa y.fa s/
    if [ -n "$OLD" ]; then
        src=$(cachefile x.fa "$OLD" $args); how="2.1.20 cache"
    else
        src=$(cachefile z.fa "$D" $args); how="simulated old cache"
    fi
    if [ -z "$new" ] || [ -z "$src" ]; then
        echo "FAIL $args: no cache file written (new '$new', old '$src')"; fail=1; continue
    fi
    cp "c/$src" "s/$old"; touch "s/$old"
    reused=$(sim $args --cache s/x.fa s/y.fa)
    if [ -n "$fresh" ] && [ "$reused" = "$fresh" ]; then
        echo "PASS $args with a $how at $old: got $reused, fresh $fresh"
    else
        echo "FAIL $args with a $how at $old: got $reused, fresh $fresh"; fail=1
    fi
done
exit $fail
