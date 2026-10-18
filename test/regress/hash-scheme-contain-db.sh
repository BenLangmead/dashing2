#!/usr/bin/env bash
# A k-mer database written by dashing2 2.1.20 (sketch --save-kmers -o) lacks
# header bit 10, and above the exact limit or with -2 its k-mers were hashed
# differently, so contain would match none of them and report zero coverage
# without complaint. Invariant: contain rejects such a database with a message
# naming the version, and still accepts untagged exact 64-bit databases.
#
# The old database comes from DASHING2_OLD (a dashing2 2.1.20 binary) when it
# is set; otherwise a database from the binary under test with bit 10 cleared
# stands in for it.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
OLD=${DASHING2_OLD:-}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
r = random.Random(7)
g = "".join(r.choice("ACGT") for _ in range(20000))
open("x.fa", "w").write(">x\n" + g + "\n")
open("q.fa", "w").write(">q\n" + g[5000:15000] + "\n")
PY
clear10() { python3 -c "import struct,sys; b=bytearray(open(sys.argv[1],'rb').read()); w=struct.unpack_from('<I',b,0)[0] & ~(1<<10); struct.pack_into('<I',b,0,w); open(sys.argv[1],'wb').write(b)" "$1"; }
fail=0
# 2.1.20 wrote -2 databases at k <= 32 without bit 9 and byte-identical to
# 64-bit ones, so those hold valid exact 64-bit k-mers and are accepted. A
# database from this version with bit 10 cleared keeps bit 9 and is rejected.
if [ -n "$OLD" ]; then short128="-2 -k 21:ok"; else short128="-2 -k 21:reject"; fi
for args in "-k 21:ok" "-k 40:reject" "-k 70:reject" "$short128" "-2 -k 40:reject"; do
    want=${args#*:}; a=${args%:*}
    rm -f db db.kmer64*
    if [ -n "$OLD" ]; then
        "$OLD" sketch $a --save-kmers -o db x.fa >/dev/null 2>&1; how="2.1.20 database"
    else
        "$D" sketch $a --save-kmers -o db x.fa >/dev/null 2>&1; clear10 db.kmer64; how="database without bit 10"
    fi
    if [ ! -s db.kmer64 ]; then echo "FAIL $a: no database written"; fail=1; continue; fi
    "$D" contain db.kmer64 q.fa >out.txt 2>err.txt; rc=$?
    if [ "$want" = reject ]; then
        if [ $rc -ne 0 ] && grep -q "2.1.20" err.txt; then
            echo "PASS $a $how rejected: $(head -c 120 err.txt)"
        else
            echo "FAIL $a $how not rejected (exit $rc): $(head -c 120 out.txt | tr '\n' ' ')"; fail=1
        fi
    else
        if [ $rc -eq 0 ] && [ -s out.txt ]; then
            echo "PASS $a $how accepted"
        else
            echo "FAIL $a $how not accepted (exit $rc): $(head -c 120 err.txt)"; fail=1
        fi
    fi
done
exit $fail
