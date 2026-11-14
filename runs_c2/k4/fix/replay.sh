#!/bin/bash
# C2-5 replay (provider refusing every call, M4c-lf stamp, route stamp) of level + corner from the merged
# layer A (lib2222m). usage: replay.sh FIELD
cd /home/user/KRcalculator
S=/tmp/claude-0/-home-user-KRcalculator/4d3e75e2-4a4f-5c6d-a935-4d9b769b978c/scratchpad
F=$1; D=runs_c2/k4/fix/replay
mkdir -p $D
ulimit -v 14000000
C="--lam 2,2,2,2 --n 4 --field $F --cache $S/lib2222m/f$F --costs runs_c2/k4/fix/costs_2222_c25.txt --stamp-route --provider m4c-lf"
for SET in level corner; do
  SF=""; [ $SET = corner ] && SF="--corner markov-ext"
  T=l2222_n4${SET/level/}_f$F
  KRHOM_COSTS_ANY_HOST=1 timeout 20000 $S/bin/libbuild replay $C $SF --jobs ${JOBS:-4} --rss 3000 --json $D/${T}_replay.jsonl >> $D/${T}_replay.out 2>&1
  echo "$(date +%T) replay $F $SET exit $?"
done
