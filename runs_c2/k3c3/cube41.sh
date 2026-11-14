#!/bin/bash
cd /home/user/KRcalculator/runs_c2/k3c3
ulimit -v 10500000
date +%s > cube_41.start
timeout 14400 ../../target/release/krhom 4 --braid "1 -2 1 -2" --strands 3 --colors 2 --colored-method cube > cube_41.out 2> cube_41.err &
P=$!; HW=0
while kill -0 $P 2>/dev/null; do
  v=$(grep VmHWM /proc/$P/status 2>/dev/null | awk '{print $2}'); [ -n "$v" ] && [ "$v" -gt "$HW" ] && HW=$v
  echo "$(( $(date +%s)-$(cat cube_41.start) )) $HW" > cube_41.mon; sleep 10
done
wait $P; echo "exit $? elapsed $(( $(date +%s)-$(cat cube_41.start) ))s peakRSS_kB $HW" > cube_41.done
