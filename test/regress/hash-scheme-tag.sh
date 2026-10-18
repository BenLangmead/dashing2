#!/usr/bin/env bash
# Values produced by the polynomial rolling hash (k above the alphabet's exact
# 64-bit limit) or by 128-bit k-mers (-2) differ from what dashing2 2.1.20
# stored under the same names and header bits. Invariant: exactly those runs
# get the .ph61 cache-name suffix and bit 10 in the --save-kmers -o database
# header; exact 64-bit runs keep their old names and header bits, since their
# values never changed.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
r = random.Random(5)
open("x.fa", "w").write(">x\n" + "".join(r.choice("ACGT") for _ in range(3000)) + "\n")
open("p.fa", "w").write(">p\n" + "".join(r.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(3000)) + "\n")
PY
# Bit 10 of the header word at byte offset $2
bit10() { python3 -c "import struct,sys; o=int(sys.argv[2]); print((struct.unpack_from('<I', open(sys.argv[1],'rb').read(), o)[0] >> 10) & 1)" "$@"; }
fail=0
check() {
    local want=$1 f=$2; shift 2
    rm -rf c; mkdir c; cp "$f" c/
    "$D" sketch "$@" --cache "c/$f" >/dev/null 2>&1
    local name; name=$(ls c | grep -v "^$f\$" | head -1)
    case "$name" in *.ph61*) got=1;; *) got=0;; esac
    if [ -n "$name" ] && [ "$got" = "$want" ]; then
        echo "PASS cache name for $*: $name"
    else
        echo "FAIL cache name for $* (tag expected: $want): '$name'"; fail=1
    fi
    if [ "$*" = "${*/--set/}" ]; then
        rm -f db db.kmer64*
        "$D" sketch "$@" --save-kmers -o db "$f" >/dev/null 2>&1
        got=$( [ -s db.kmer64 ] && bit10 db.kmer64 0 || echo none)
        if [ "$got" = "$want" ]; then
            echo "PASS k-mer database header bit 10 for $*: $got"
        else
            echo "FAIL k-mer database header bit 10 for $* (expected $want): $got"; fail=1
        fi
    fi
}
check 0 x.fa -k 21
check 0 x.fa -k 32
check 1 x.fa -k 33
check 1 x.fa -k 70
check 1 x.fa -2 -k 21
check 1 x.fa -2 -k 70
check 0 x.fa --set -k 31
check 1 x.fa --set -k 40
check 0 p.fa --protein -k 14
check 1 p.fa --protein -k 15
exit $fail
