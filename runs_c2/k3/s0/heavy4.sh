#!/bin/bash
# after heavy3.sh: knot-runs refused by the staged-rank memory ESTIMATE (exit 2, FATAL cap): retry alone with the cap
# raised to 200000 MB; ulimit -v 13 GB stays the real guard.
cd /home/user/KRcalculator
while pgrep -f "s0/queue2.sh|s0/heavy.sh|s0/heavy2.sh|s0/heavy3.sh" >/dev/null || pgrep -x sbox >/dev/null; do sleep 60; done
for e in runs_c2/k3/s0/out/*.err; do
  grep -q "FATAL: box: staged rank stage" $e || continue
  T=$(basename $e .err); N=${T%_f*}; F=${T##*_f}
  W=$(grep "^$N;" runs_c2/k3/slices/slice0.txt | cut -d';' -f2)
  rm -f runs_c2/k3/s0/out/$T.*
  MEMCAP=200000 MEMMB=200000 runs_c2/k3/s0/one.sh $N $F "$W"
done
