#!/bin/bash
cd /home/user/KRcalculator
ulimit -v 6000000
for N in 4 8; do
  T0=$(date +%s.%N)
  timeout 7200 /tmp/claude-0/-home-user-KRcalculator/934ab925-996a-59b6-8fc2-b741a70e6e9a/scratchpad/bin/libbuild plan --lam 2,2,2,2 --n $N --json runs_c2/census/plan_2222_n$N.json > runs_c2/census/plan_2222_n$N.out 2>&1
  echo "n=$N exit $? secs $(echo "$(date +%s.%N) - $T0" | bc)" >> runs_c2/census/plan_2222_n$N.out
done
