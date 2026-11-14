#!/bin/bash
# runs after heavy.sh: rerun 12a_920 F31 (OOM-killed after 7 h) and 12a_1218 F31 (stopped to protect memory), alone.
cd /home/user/KRcalculator
while pgrep -f s0/heavy.sh >/dev/null || pgrep -x sbox >/dev/null; do sleep 60; done
for t in "12a_1218 31" "12a_920 31"; do
  set -- $t; N=$1; F=$2
  W=$(grep "^$N;" runs_c2/k3/slices/slice0.txt | cut -d';' -f2)
  if [ -f runs_c2/k3/s0/out/${N}_f$F.err ] && grep -q "^exit 137" runs_c2/k3/s0/out/${N}_f$F.err; then rm -f runs_c2/k3/s0/out/${N}_f$F.*; fi
  runs_c2/k3/s0/one.sh $N $F "$W"
done
