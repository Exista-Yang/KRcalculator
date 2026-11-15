#!/bin/bash
# CALC-2 demonstrations of items 1, 2, 3, 5 (and 15, 16); stdout/stderr/exit recorded per command
cd "$(dirname "$0")/.."
K=target/release/krhom
D=runs_calc2/demo
mkdir -p $D
export KRHOM_THREADS=4
ulimit -v 14000000
run() { local tag=$1; shift; { echo "\$ krhom $(printf '%q ' "$@")"; } > $D/$tag.cmd; timeout 600 $K "$@" > $D/$tag.out 2> $D/$tag.err; echo "exit $?" > $D/$tag.exit; }
# item 1: memory cap on the exact 2-strand route (exit 5), --mem-cap-mb honoured, shape (exit 1)
run i1_cap 30 --braid "1 1 1" --colors 3
run i1_flag 2,3 --braid "1 1 1" --colors 2 --mem-cap-mb 1
run i1_shape 130 --braid "1 1 1" --colors 2
# item 2: N below the color before any refusal
run i2_nolib 1,2 --braid aBaB --colors 3
run i2_c4 2,3 --braid "1 1 1" --colors 4
# item 3: batch names before early returns / on the error line; item 4 CSV quoting
printf 'L6a5{0,1}; aaa\nbad; 1 0\n"q"; aBaB\n' > $D/i3.txt
run i3_batch 2 --file $D/i3.txt --csv $D/i3.csv --json $D/i3.json
# item 5: a batch killed while its second line runs keeps the first record
printf 'first; aaa\nslow; %s\n' "$(printf 'a%.0s' $(seq 80))" > $D/i5.txt
echo "\$ timeout -s KILL 20 krhom 5 --file $D/i5.txt --colors 2 --json $D/i5.json --csv $D/i5.csv" > $D/i5_kill.cmd
timeout -s KILL 20 $K 5 --file $D/i5.txt --colors 2 --json $D/i5.json --csv $D/i5.csv > $D/i5_kill.out 2> $D/i5_kill.err; echo "exit $?" > $D/i5_kill.exit
# items 15, 16
run i16_trefoil 2 --braid aaa --convention knotjob --table
run i16_trefoil_red 2 --braid aaa --convention knotjob --reduced --table
run i15_table 2,3,4 --braid "1 1 1" --colors 2 --table
echo DONE
