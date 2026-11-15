#!/bin/bash
# CALC-2: re-run of the CALC-1 acceptance runs H3-H6 (same commands as runs_calc/hrun.sh, h6.sh and
# the H3 unreduced run), sequential, KRHOM_THREADS=4; stdout compared with runs_calc/ by cmp_h.sh
cd "$(dirname "$0")/.."
K=target/release/krhom
O=runs_calc2/h
mkdir -p $O
ulimit -v 14000000
export KRHOM_THREADS=4
run() { local tag=$1; shift; echo "=== $tag: krhom $*"; s=$(date +%s); timeout 7200 $K "$@" --json $O/$tag.json > $O/$tag.out 2> $O/$tag.err; echo "exit $? ($(( $(date +%s)-s ))s)" >> $O/$tag.out; }
run h3_u 2,3,4 --braid "1 -2 1 -2" --strands 3 --colors 2
run h3_r 2,3,4 --braid "1 -2 1 -2" --strands 3 --colors 2 --reduced
run h4_u 2,3,4 --braid "1 1 1 2" --strands 3 --colors 2 --no-simplify
run h4_r 2,3,4 --braid "1 1 1 2" --strands 3 --colors 2 --no-simplify --reduced
run h4_ref 2,3,4 --braid "1 1 1" --strands 2 --colors 2
run h4_refr 2,3,4 --braid "1 1 1" --strands 2 --colors 2 --reduced
run h5a_u 2,3,4 --braid "1 2 3" --strands 4 --colors 2 --no-simplify
run h5a_r 2,3,4 --braid "1 2 3" --strands 4 --colors 2 --no-simplify --reduced
run h5b_u 2,3,4 --braid "1 1 1 2 3" --strands 4 --colors 2 --no-simplify
run h5b_r 2,3,4 --braid "1 1 1 2 3" --strands 4 --colors 2 --no-simplify --reduced
run h6_61 2,3,4 --braid "1 1 2 -1 -3 2 -3" --strands 4 --colors 2 --check-library
run h6_815 2,3,4 --braid "1 1 -2 1 3 2 2 2 3" --strands 4 --colors 2 --check-library
# H8: nothing written into the libraries
for d in ${KRHOM_LIB_DIR:-$HOME/.cache/krhom/lib}/*/*; do
  echo "H8 $d: newer than LIBRARY.txt: $(find $d -newer $d/LIBRARY.txt -type f | wc -l); layer_a $(ls $d/layer_a | wc -l); tmpl $(ls $d/tmpl | wc -l)"
done
tools/fetch_library.sh --verify-only l222_n4 l2222_n4; echo "verify-only exit $?"
echo ALLDONE
