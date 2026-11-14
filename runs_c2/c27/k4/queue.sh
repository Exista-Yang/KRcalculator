#!/bin/bash
# C2-7 knot queue (sequential; each run uses 4 shards + a 4-thread finish). Lines: NAME;FIELD;MODES;WORD
cd /home/user/KRcalculator
while pgrep -f "box.sh 6_1 32b" > /dev/null; do sleep 20; done
while IFS=';' read NAME F M W; do
  [ -z "$NAME" ] && continue
  MODES=$M runs_c2/c27/k4/box.sh $NAME $F "$W" 4
done < ${1:-runs_c2/c27/k4/queue.txt}
