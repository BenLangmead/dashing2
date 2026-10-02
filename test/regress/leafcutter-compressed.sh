#!/usr/bin/env bash
# `sketch --leafcutter` reads .xz, .bz2 and .zst LeafCutter count files through
# a decompression pipe, as sequence inputs are, so a compressed copy of a file
# gives the same sketch file as the plain file. Before the fix zlib passed the
# compressed bytes through unchanged and dashing2 crashed parsing the header.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
python3 - <<'PY'
import random
r = random.Random(5)
samples = ['s1', 's2', 's3']
with open('a_perind.counts', 'w') as f:
    f.write('chrom ' + ' '.join(samples) + '\n')
    for i in range(300):
        s = 1000 + 37 * i
        f.write('chr1:%d:%d:clu_%d ' % (s, s + 200, i // 3) + ' '.join('%d/10' % r.randint(0, 3) for _ in samples) + '\n')
PY
"$D2" sketch --leafcutter -o plain.sk a_perind.counts >/dev/null 2>&1 || { echo "FAIL plain input: exit status $?"; exit 1; }
for tool in gzip bzip2 xz zstd; do
    if ! command -v $tool >/dev/null; then echo "SKIP $tool: not installed"; continue; fi
    case $tool in gzip) sfx=.gz;; bzip2) sfx=.bz2;; xz) sfx=.xz;; zstd) sfx=.zst;; esac
    $tool -kqf a_perind.counts
    "$D2" sketch --leafcutter -o c.sk a_perind.counts$sfx >/dev/null 2>&1; rc=$?
    if [ $rc -eq 0 ] && cmp -s plain.sk c.sk; then echo "PASS a_perind.counts$sfx: sketch equals the plain file's"
    else echo "FAIL a_perind.counts$sfx: exit status $rc, sketch $(cmp -s plain.sk c.sk && echo equals || echo differs from) the plain file's"; fail=1; fi
    rm -f c.sk c.sk.*
done
exit $fail
