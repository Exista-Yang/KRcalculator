#!/bin/bash
# H3-H5 acceptance runs (CALC-1), sequential, KRHOM_THREADS=4
cd "$(dirname "$0")/.."
K=/tmp/claude-0/bin/krhom
ulimit -v 14000000
run() { local tag=$1; shift; echo "=== $tag: krhom $*" ; timeout 7200 $K "$@" --json runs_calc/$tag.json > runs_calc/$tag.out 2> runs_calc/$tag.err; echo "exit $?" >> runs_calc/$tag.out; }
export KRHOM_THREADS=4
run h3_r 2,3,4 --braid "1 -2 1 -2" --strands 3 --colors 2 --reduced
run h4_u 2,3,4 --braid "1 1 1 2" --strands 3 --colors 2 --no-simplify
run h4_r 2,3,4 --braid "1 1 1 2" --strands 3 --colors 2 --no-simplify --reduced
run h4_ref 2,3,4 --braid "1 1 1" --strands 2 --colors 2
run h4_refr 2,3,4 --braid "1 1 1" --strands 2 --colors 2 --reduced
run h5a_u 2,3,4 --braid "1 2 3" --strands 4 --colors 2 --no-simplify
run h5a_r 2,3,4 --braid "1 2 3" --strands 4 --colors 2 --no-simplify --reduced
run h5b_u 2,3,4 --braid "1 1 1 2 3" --strands 4 --colors 2 --no-simplify
run h5b_r 2,3,4 --braid "1 1 1 2 3" --strands 4 --colors 2 --no-simplify --reduced
echo ALLDONE
