#!/usr/bin/env bash
# BigWig sketch caches follow the BED rules: written only with --cache, named
# <input>.sketchsize<N>[.ct_threshold<M>]<sketch suffix> (under --outprefix,
# <name>.<hash of absolute path><suffix>), and not reused once the input was
# modified after the cache. The input is libBigWig's test file; b.bw is a copy
# whose chromosomes are renamed ("1" -> "2", "10" -> "20") by patching the
# uncompressed chromosome B+ tree, so a~b share no positions.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
BW=${TEST_BW:-$here/../../libBigWig/test/test.bw}
[ -s "$BW" ] || { echo "FAIL missing BigWig test input $BW"; exit 1; }
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
check() { if [ "$2" = "$3" ]; then echo "PASS $1: got $2 expected $3"; else echo "FAIL $1: got $2 expected $3"; fail=1; fi; }
cell() { awk -F'\t' '!/^#/ && NF > 2 {print $3; exit}'; }
rename_chroms() {
    python3 - "$1" "$2" <<'PY'
import struct, sys
d = bytearray(open(sys.argv[1], 'rb').read())
ct = struct.unpack_from('<Q', d, 8)[0]                  # chromosome tree offset
ks = struct.unpack_from('<I', d, ct + 8)[0]             # key size
leaf, _, cnt = struct.unpack_from('<BBH', d, ct + 32)
assert leaf == 1
for i in range(cnt):
    off = ct + 36 + i * (ks + 8)
    assert d[off:off + 1] == b'1'
    d[off] = ord('2')
open(sys.argv[2], 'wb').write(d)
PY
}
mkdir in; cp "$BW" in/a.bw
for opts in "--full" "--multiset" "--prob"; do
    case $opts in --full) sfx=.ss;; --multiset) sfx=.bmh;; --prob) sfx=.pmh;; esac
    rm -rf work; mkdir work; cd work
    cp ../in/a.bw a.bw
    nocache=$("$D2" cmp --bigwig $opts a.bw a.bw 2>/dev/null | cell)
    check "bigwig $opts files written without --cache" "$(ls | tr '\n' ' ')" "a.bw "
    first=$("$D2" cmp --bigwig $opts --cache a.bw a.bw 2>/dev/null | cell)
    check "bigwig $opts --cache file" "$(ls | tr '\n' ' ')" "a.bw a.bw.sketchsize1024$sfx "
    again=$("$D2" cmp --bigwig $opts --cache a.bw a.bw 2>/dev/null | cell)
    check "bigwig $opts a~a uncached, cached, reloaded" "$nocache $first $again" "1 1 1"
    # Inputs with the same name in different directories get their own caches under --outprefix.
    mkdir d1 d2 op; cp a.bw d1/x.bw; cp a.bw d2/x.bw
    "$D2" cmp --bigwig $opts --cache --outprefix op d1/x.bw d2/x.bw >/dev/null 2>&1
    check "bigwig $opts --outprefix caches" "$(ls op | sed 's/\.[0-9a-f]\{16\}\./.HASH./' | sort | uniq -c | awk '{print $1, $2}' | tr '\n' ' ')" "2 x.bw.HASH.sketchsize1024$sfx "
    # A cache older than its input is not reused.
    rename_chroms a.bw b.bw
    ab=$("$D2" cmp --bigwig $opts --cache a.bw b.bw 2>/dev/null | cell)
    cp a.bw b.bw; cp a.bw c.bw
    bc=$("$D2" cmp --bigwig $opts --cache b.bw c.bw 2>/dev/null | cell)
    check "bigwig $opts --cache a~b with renamed chromosomes, then b~c after overwriting b with a" "$ab $bc" "0 1"
    cd ..
done
exit $fail
