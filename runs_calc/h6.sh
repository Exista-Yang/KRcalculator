#!/bin/bash
cd "$(dirname "$0")/.."
K=/tmp/claude-0/bin/krhom
ulimit -v 14000000
export KRHOM_THREADS=4
run() { local tag=$1; shift; echo "=== $tag: krhom $*"; s=$(date +%s); timeout 7200 $K "$@" --json runs_calc/$tag.json > runs_calc/$tag.out 2> runs_calc/$tag.err; echo "exit $? ($(( $(date +%s)-s ))s)" >> runs_calc/$tag.out; }
run h6_61 2,3,4 --braid "1 1 2 -1 -3 2 -3" --strands 4 --colors 2 --check-library
run h6_815 2,3,4 --braid "1 1 -2 1 3 2 2 2 3" --strands 4 --colors 2 --check-library
echo ALLDONE
