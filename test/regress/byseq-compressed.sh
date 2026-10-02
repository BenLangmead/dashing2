#!/usr/bin/env bash
# --parse-by-seq must read every input format that file mode reads: plain,
# .gz, .bz2, .xz and .zst (the last three are decompressed through a pipe).
# Each compressed copy of a three-record file must give the same record names
# and the same matrix as the plain file, and those values must equal file mode
# run on the three records saved as separate files.
set -u
here=$(cd "$(dirname "$0")" && pwd)
D2=${DASHING2:-$here/../../dashing2}
T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cd "$T"
fail=0
check() { # name got expected
    if [ "$2" = "$3" ]; then echo "PASS $1: got $2 expected $3"; else echo "FAIL $1: got $2 expected $3"; fail=1; fi
}
# r1: 2 kbp random; r2: r1 with every 25th base changed; r3: r1's first half
# followed by 1 kbp of new sequence.
python3 - <<'PY'
import random
r = random.Random(5)
x = ''.join(r.choice('ACGT') for _ in range(2000))
y = ''.join('ACGT'[('ACGT'.index(c) + 1) % 4] if i % 25 == 12 else c for i, c in enumerate(x))
z = x[:1000] + ''.join(r.choice('ACGT') for _ in range(1000))
with open('r.fa', 'w') as f:
    for n, s in [('r1', x), ('r2', y), ('r3', z)]:
        f.write(f'>{n} comment\n{s}\n')
        open(n + '.fa', 'w').write(f'>{n}\n{s}\n')
PY
# Matrix values only (row labels dropped), and the names from the #Sources line.
vals() { grep -v '^#' | cut -f2- | tr '\n\t' '  '; }
names() { grep '^#Sources' | cut -f2- | tr '\t' ' ' | LC_ALL=C tr -cd '[:print:]'; }
run() { "$D2" cmp -k 21 --full "$@" 2>/dev/null; }
plain=$(run --parse-by-seq r.fa)
check "plain names" "$(names <<<"$plain")" "r1 r2 r3"
check "plain values equal file mode" "$(vals <<<"$plain")" "$(run r1.fa r2.fa r3.fa | vals)"
for z in "gzip .gz" "bzip2 .bz2" "xz .xz" "zstd .zst"; do
    set -- $z
    if ! command -v "$1" >/dev/null; then echo "SKIP $2: $1 not installed"; continue; fi
    "$1" -c r.fa > "r.fa$2"
    out=$(run --parse-by-seq "r.fa$2")
    check "$2 names" "$(names <<<"$out")" "r1 r2 r3"
    check "$2 values" "$(vals <<<"$out")" "$(vals <<<"$plain")"
done
exit $fail
