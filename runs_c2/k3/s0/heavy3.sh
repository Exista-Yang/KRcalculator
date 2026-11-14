#!/bin/bash
# after heavy2.sh: rerun every knot-run that was OOM-killed (exit 137), alone, in slice order.
cd /home/user/KRcalculator
while pgrep -f "s0/queue2.sh|s0/heavy.sh|s0/heavy2.sh" >/dev/null || pgrep -x sbox >/dev/null; do sleep 60; done
for e in runs_c2/k3/s0/out/*.err; do
  grep -q "^exit 137" $e || continue
  T=$(basename $e .err); N=${T%_f*}; F=${T##*_f}
  W=$(grep "^$N;" runs_c2/k3/slices/slice0.txt | cut -d';' -f2)
  rm -f runs_c2/k3/s0/out/$T.*
  runs_c2/k3/s0/one.sh $N $F "$W"
done
