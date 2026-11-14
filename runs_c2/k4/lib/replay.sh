#!/bin/bash
# replay (provider refusing every call) of level + corner, calibrated table + route stamp. usage: replay.sh FIELD
cd /home/user/KRcalculator
S=/tmp/claude-0/-home-user-KRcalculator/7ca6fc3a-b081-5be7-b7f0-2501a115a78c/scratchpad
F=$1; D=runs_c2/k4/lib
ulimit -v 14000000
C="--lam 2,2,2,2 --n 4 --field $F --cache $S/lib2222/f$F --costs runs_c2/k4/calib/costs_2222.txt ${SR---stamp-route}"
for SET in level corner; do
  SF=""; [ $SET = corner ] && SF="--corner markov-ext"
  T=l2222_n4${SET/level/}_f$F
  KRHOM_COSTS_ANY_HOST=1 timeout 20000 $S/bin/libbuild replay $C $SF --jobs ${JOBS:-4} --rss 3000 --json $D/${T}_replay.jsonl >> $D/${T}_replay.out 2>&1
  echo "$(date +%T) replay $F $SET exit $?"
done
