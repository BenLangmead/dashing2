#!/bin/bash
# --entmin must weight minimizer selection by k-mer entropy, as the help text
# says: in every window the selected k-mer must minimize
#   score(x) = (top 48 bits of lex hash(x)) / (entropy(x) + 1e-4)
# where lex hash is bonsai's FRev64 and entropy is in nats over the k-mer's
# bases. The python oracle recomputes this for every window of "sketch -G"
# output (one wang-hashed k-mer per full window), with and without --no-canon.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
DASHING2=${DASHING2:-$HERE/../../dashing2}
D=$(mktemp -d)
trap 'rm -rf "$D"' EXIT
cd "$D"

cat > oracle.py <<'EOF'
import math, random, struct, sys, glob
M64 = (1 << 64) - 1
def wang(key):
    key = ((~key) + (key << 21)) & M64; key ^= key >> 24
    key = (key + (key << 3) + (key << 8)) & M64; key ^= key >> 14
    key = (key + (key << 2) + (key << 4)) & M64; key ^= key >> 28
    return (key + (key << 31)) & M64
def lex(x):  # bonsai FRev64: xor, multiply, rotate left 31, xor
    h = ((x ^ 0x533f8c2151b20f97) * 0x9a98567ed20c127d) & M64
    return (((h << 31) | (h >> 33)) & M64) ^ 0x691a9d706391077a
def enc(s): return int(s.translate(str.maketrans("ACGT", "0123")), 4)
def canon(s): return min(s, s[::-1].translate(str.maketrans("ACGT", "TGCA")))
def score(s):
    ent = -sum(s.count(c) / len(s) * math.log(s.count(c) / len(s)) for c in "ACGT" if c in s)
    return (lex(enc(s)) >> 16) / (ent + 1e-4)
if sys.argv[1] == "gen":
    random.seed(3)
    s = "".join("".join(random.choice("ACGT") for _ in range(40)) + "ACGT"[i % 4] * (10 + i) + "C" * (i % 7 + 3) + "AT" for i in range(20))
    open("a.fa", "w").write(">a\n%s\n" % s)
    sys.exit(0)
k, w, can = int(sys.argv[2]), int(sys.argv[3]), sys.argv[4] == "canon"
s = open("a.fa").read().split("\n")[1]
kms = [canon(s[i:i + k]) if can else s[i:i + k] for i in range(len(s) - k + 1)]
got = list(struct.unpack("<%dQ" % (len(open(sys.argv[1], "rb").read()) // 8), open(sys.argv[1], "rb").read()))
m = w - k + 1
bad = 0
for j in range(len(kms) - m + 1):
    win = kms[j:j + m]
    best = min(score(x) for x in win)
    picked = [x for x in win if wang(enc(x)) == (got[j] if j < len(got) else -1)]
    if not picked or score(picked[0]) > best * (1 + 1e-12) + 1: bad += 1
print("%d/%d" % (bad, len(kms) - m + 1))
EOF
python3 oracle.py gen

fail=0
for mode in canon nocanon; do
    extra=; [ $mode = nocanon ] && extra=--no-canon
    rm -f a.fa*.mmerseq*
    "$DASHING2" sketch -G --entmin -k 21 -w 42 $extra -o out a.fa >/dev/null 2>&1
    got=$(python3 oracle.py a.fa*.mmerseq* 21 42 $mode)
    total=${got#*/}
    if [ "$got" = "0/$total" ]; then echo "PASS entmin_windows_$mode: windows with a non-minimal pick $got expected 0/$total"
    else echo "FAIL entmin_windows_$mode: windows with a non-minimal pick $got expected 0/$total"; fail=1; fi
done
rm -f a.fa*.mmerseq* out
exit $fail
