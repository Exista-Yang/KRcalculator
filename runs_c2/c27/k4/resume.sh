#!/bin/bash
# C2-7 after a container restart: resume the reduced queue (3 shards; box runs resume from checkpoints)
cd /home/user/KRcalculator
tail -n +2 runs_c2/c27/k4/queue2.txt | while IFS=';' read NAME F M WW; do
  [ -z "$NAME" ] && continue
  MODES=$M runs_c2/c27/k4/box.sh $NAME $F "$WW" 3
done
