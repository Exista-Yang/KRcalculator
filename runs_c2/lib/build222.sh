#!/bin/bash
# (2,2,2) library: level N (prune N) run + replay for both primes. usage: build222.sh N
cd /home/user/KRcalculator
S=/tmp/claude-0/-home-user-KRcalculator/934ab925-996a-59b6-8fc2-b741a70e6e9a/scratchpad
N=$1
ulimit -v 14000000
for F in 31 32b; do
  KRHOM_COSTS_ANY_HOST=1 timeout 10000 $S/bin/libbuild run --lam 2,2,2 --n $N --field $F --cache $S/lib222${CSUF}/f$F --costs tests/data/homprov/costs_cloud.txt --jobs 4 --wall 900 --rss 3000 --json runs_c2/lib/l222_n${N}_f${F}_run.jsonl > runs_c2/lib/l222_n${N}_f${F}_run.out 2>&1
  echo "run $F $N exit $?"
  KRHOM_COSTS_ANY_HOST=1 timeout 10000 $S/bin/libbuild replay --lam 2,2,2 --n $N --field $F --cache $S/lib222${CSUF}/f$F --costs tests/data/homprov/costs_cloud.txt --jobs 4 --rss 4000 --json runs_c2/lib/l222_n${N}_f${F}_replay.jsonl > runs_c2/lib/l222_n${N}_f${F}_replay.out 2>&1
  echo "replay $F $N exit $?"
done
