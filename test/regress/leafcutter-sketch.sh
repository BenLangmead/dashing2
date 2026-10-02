#!/usr/bin/env bash
# `sketch --leafcutter` and `cmp --leafcutter` on a three-sample LeafCutter
# count file. The stacked -o file must hold one sketch per sample (header of
# entity count and sketch size, then one cardinality and one sketch per sample),
# the cardinalities must estimate each sample's number of junctions, and cmp
# must compare the samples with or without -o. LFResult::cardinalities() used
# to return the per-file sample counts, so the file size and the cardinalities
# were wrong and cmp aborted.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
check() { # name got expected
    if [ "$2" = "$3" ]; then echo "PASS $1: got $2 expected $3"; else echo "FAIL $1: got $2 expected $3"; fail=1; fi
}
# 600 junctions. s1 and s2 are nonzero on the same 300 junctions (Jaccard 1);
# s3 is nonzero on the last 150 of those and on 150 others (Jaccard 150/450).
python3 - <<'PY'
with open('a_perind.counts', 'w') as f:
    f.write('chrom s1 s2 s3\n')
    for i in range(600):
        s = 1000 + 37 * i
        c1 = 1 + i % 3 if i < 300 else 0
        c3 = 2 if 150 <= i < 450 else 0
        f.write('chr1:%d:%d:clu_%d %d/10 %d/10 %d/10\n' % (s, s + 200, i // 3, c1, c1, c3))
with open('one.fa', 'w') as f:
    f.write('>x\n' + 'ACGTTGCA' * 50 + '\n')
PY
S=1024
# Bytes per register, from a one-entity stacked file (16-byte header, one cardinality).
"$D2" sketch -S $S -o one.sk one.fa >/dev/null 2>&1
w=$(( ($(wc -c < one.sk) - 24) / S ))
"$D2" sketch --leafcutter -S $S -o lf.sk a_perind.counts >/dev/null 2>&1
check "sketch exit status" $? 0
check "stacked file size" "$(wc -c < lf.sk 2>/dev/null)" $(( 16 + 3 * 8 + 3 * S * w ))
check "header" "$(od -An -tu8 -N16 lf.sk 2>/dev/null | xargs)" "3 $S"
# Estimated cardinalities from the names file: 300, 300 and 300 junctions, within 10%.
card=$(python3 -c '
import sys
rows = [l.split("\t") for l in open("lf.sk.names.txt") if not l.startswith("#")]
print(" ".join("%s:%s" % (r[0], "ok" if abs(float(r[1]) - 300) <= 30 else r[1].strip()) for r in rows))' 2>/dev/null)
check "cardinalities" "$card" "s1:a:ok s2:a:ok s3:a:ok"
# Pairwise similarities with --full (the default one-permutation sketch on
# v2.1.20 misjudges sets this sparse): s1~s2 exactly 1, s1~s3 and s2~s3 within
# 0.05 of 1/3.
sims() {
    python3 -c '
import sys
v = {}
for l in sys.stdin:
    if l.startswith("#"): continue
    f = l.split()
    v[f[0]] = f[1:]
def ok(x, e): return "ok" if x not in ("", "-") and abs(float(x) - e) <= 0.05 else x
try: print("s1s2:%s s1s3:%s s2s3:%s" % (ok(v["s1:a"][1], 1), ok(v["s1:a"][2], 1 / 3), ok(v["s2:a"][2], 1 / 3)))
except Exception as e: print("unparsable")'
}
check "cmp -o similarities" "$("$D2" cmp --leafcutter --full -S $S -o lf2.sk a_perind.counts 2>/dev/null | sims)" "s1s2:ok s1s3:ok s2s3:ok"
check "cmp without -o similarities" "$("$D2" cmp --leafcutter --full -S $S a_perind.counts 2>/dev/null | sims)" "s1s2:ok s1s3:ok s2s3:ok"
# Without -o or --cmpout, sketch has nowhere to write and must say so.
"$D2" sketch --leafcutter -S $S a_perind.counts >/dev/null 2>err.txt; rc=$?
check "sketch without -o exit status" $rc 1
check "sketch without -o message" "$(grep -c 'pass -o' err.txt)" 1
exit $fail
