#!/bin/bash
cd /home/user/KRcalculator
for M in cube cable; do
  ( ulimit -v 5000000; KRHOM_THREADS=1 timeout 5400 /tmp/claude-0/-home-user-KRcalculator/934ab925-996a-59b6-8fc2-b741a70e6e9a/scratchpad/bin/krhom 4 --braid "1 -2 1 -2" --strands 3 --colors 2 --colored-method $M > runs_c2/ref/4_1_$M.out 2> runs_c2/ref/4_1_$M.err; echo "exit $?" >> runs_c2/ref/4_1_$M.err )
done
