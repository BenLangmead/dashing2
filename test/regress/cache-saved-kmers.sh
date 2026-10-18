#!/usr/bin/env bash
# With -B or -P and -s/--save-kmers or -N/--save-kmercounts, a --cache run that
# loads its sketches must write the same stacked files (X, X.kmer64 and
# X.kmercounts.f64) as the run that computed and cached them. The per-input
# .kmer.u64 and .kmercounts.f64 cache files have no cardinality header, but the
# loader read their first 8 bytes as the input's cardinality and shifted the
# k-mers and counts by one. Fixtures: two random 6 kbp genomes (the second has
# 3% substitutions), each with a 300 bp block repeated four times so that k-mer
# counts differ.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
python3 - <<'PY'
import random
r = random.Random(11)
rep = ''.join(r.choice('ACGT') for _ in range(300))
a = ''.join(r.choice('ACGT') for _ in range(6000)) + rep * 4
b = list(a)
for p in r.sample(range(len(b)), len(b) * 3 // 100):
    b[p] = r.choice('ACGT'.replace(b[p], ''))
open('a.fa', 'w').write('>a\n' + a + '\n')
open('b.fa', 'w').write('>b\n' + ''.join(b) + '\n')
PY
cards() {  # cardinalities in a stacked file: uint64 count, uint64 sketch size, then one double per input
    python3 -c 'import struct,sys; d=open(sys.argv[1],"rb").read(); n=struct.unpack_from("<Q",d)[0]; print(",".join("%g" % c for c in struct.unpack_from("<%dd" % n, d, 16)))' "$1"
}
fail=0
for mode in "-B -s" "-B -N" "-P -s" "-P -N"; do
    name=${mode//[- ]/}
    mkdir "$name"; cp a.fa b.fa "$name"
    (cd "$name" && timeout 120 "$D" sketch -k 21 -S 256 $mode --cache -o X1 a.fa b.fa &&
                   timeout 120 "$D" sketch -k 21 -S 256 $mode --cache -o X a.fa b.fa) >/dev/null 2>&1
    for f in "" .kmer64 .kmercounts.f64; do
        [ -e "$name/X1$f" ] || continue
        if cmp -s "$name/X1$f" "$name/X$f"; then
            echo "PASS $name X$f: loaded from cache identical to computed (cardinalities got $(cards "$name/X") expected $(cards "$name/X1"))"
        else
            echo "FAIL $name X$f: loaded from cache differs from computed (cardinalities got $(cards "$name/X") expected $(cards "$name/X1"))"; fail=1
        fi
    done
done
exit $fail
