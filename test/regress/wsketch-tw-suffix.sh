#!/usr/bin/env bash
# One- and two-path `wsketch` writes <prefix>.sampled.tw.txt as
# "Total weight: W;<ids>[;<weights>];<weight type>;<id width>": the weight type
# is f (-f, float32), d (float64), H (-H, 16-bit) or U (-U, 32-bit), the id width
# W (-u, 32-bit) or L (64-bit). Before the fix the two fields were summed as
# characters and written as one garbage character.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
python3 - <<'PY'
import struct
ids, w = [5, 9, 1234567, 42, 77], [1, 2, 3, 4, 5]
open('ids.u64', 'wb').write(struct.pack('<5Q', *ids))
open('w.f64', 'wb').write(struct.pack('<5d', *w))
open('w.f32', 'wb').write(struct.pack('<5f', *w))
open('w.u16', 'wb').write(struct.pack('<5H', *w))
PY
check() { # expected-suffix wsketch-args...
    local want=$1; shift
    rm -f o.sampled.*
    "$D2" wsketch -S 16 -o o "$@" >/dev/null 2>&1
    local got; got=$(LC_ALL=C sed -n '1s/^Total weight: [0-9.]*;//p' o.sampled.tw.txt 2>/dev/null | LC_ALL=C cat -v)
    if [ "$got" = "$want" ]; then echo "PASS wsketch $*: got $got expected $want"
    else echo "FAIL wsketch $*: got $got expected $want"; fail=1; fi
}
check "ids.u64;w.f64;d;L" ids.u64 w.f64
check "ids.u64;w.f32;f;L" -f ids.u64 w.f32
check "ids.u64;w.u16;H;L" -H ids.u64 w.u16
check "ids.u64;d;W" -u ids.u64
exit $fail
