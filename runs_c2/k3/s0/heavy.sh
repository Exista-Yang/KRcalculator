#!/bin/bash
# After queue2 finishes: rerun the OOM-killed knots one job at a time (they peak at > 6-8 GB).
cd /home/user/KRcalculator
while pgrep -f s0/queue2.sh >/dev/null || pgrep -x sbox >/dev/null; do sleep 60; done
for t in "12a_920 32b" "12a_1215 31" "12a_1215 32b"; do
  set -- $t; N=$1; F=$2
  W=$(grep "^$N;" runs_c2/k3/slices/slice0.txt | cut -d';' -f2)
  if [ -f runs_c2/k3/s0/out/${N}_f$F.err ] && grep -q "^exit 137" runs_c2/k3/s0/out/${N}_f$F.err; then rm -f runs_c2/k3/s0/out/${N}_f$F.*; fi
  runs_c2/k3/s0/one.sh $N $F "$W"
done
