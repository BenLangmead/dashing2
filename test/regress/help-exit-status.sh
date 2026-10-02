#!/usr/bin/env bash
# Asking for the usage message (-h or --help) is not an error, so it exits with
# status 0; an unknown option still prints the usage and exits with status 1.
set -u
D=${DASHING2:-$(cd "$(dirname "$0")" && pwd)/../../dashing2}
fail=0
check() {  # check <name> <got> <expected>
    if [ "$2" = "$3" ]; then echo "PASS $1: got $2 expected $3"
    else echo "FAIL $1: got $2 expected $3"; fail=1; fi
}
rc() { "$D" "$@" >/dev/null 2>&1; echo $?; }
check "dashing2 --help" "$(rc --help)" 0
check "dashing2 -h" "$(rc -h)" 0
for sub in sketch cmp dist; do
    check "$sub -h" "$(rc $sub -h)" 0
    check "$sub --help" "$(rc $sub --help)" 0
    check "$sub unknown option" "$(rc $sub -y)" 1
done
for sub in contain wsketch printmin; do
    check "$sub -h" "$(rc $sub -h)" 0
done
check "contain unknown option" "$(rc contain -y)" 1
check "wsketch unknown option" "$(rc wsketch -y)" 1
check "dashing2 without a subcommand" "$(rc)" 1
exit $fail
