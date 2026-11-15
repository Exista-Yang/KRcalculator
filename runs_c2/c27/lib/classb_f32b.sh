#!/bin/bash
# class-B jobs certified for F31: F32b runs, then replays of both primes (when the 3 big 7_2 shards are done)
cd /home/user/KRcalculator
while [ $(pgrep -f "sbox box" | wc -l) -gt 2 ]; do sleep 20; done
export COSTS=/home/user/work/costs_2222_c27.txt VLIM=4000000 WALL=5400
PAR=2 runs_c2/c27/lib/runjobs.sh 32b runs_c2/c27/lib/jobs_classb_ok.txt
PAR=2 runs_c2/c27/lib/runjobs.sh 31 runs_c2/c27/lib/jobs_classb_ok.txt replay
PAR=2 runs_c2/c27/lib/runjobs.sh 32b runs_c2/c27/lib/jobs_classb_ok.txt replay
