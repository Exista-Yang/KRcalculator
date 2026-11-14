#!/bin/bash
# engine2 uncolored gl_3 of the reversed braid (N' = 3 duality reference) for the 4-strand knots
cd /home/user/KRcalculator
S=/tmp/claude-0/-home-user-KRcalculator/7ca6fc3a-b081-5be7-b7f0-2501a115a78c/scratchpad
ulimit -v 4000000
while IFS=';' read NAME LAM W NN; do
  f=runs_c2/k4/e2ref/$(echo $W | tr ' ' '_')_gl3.txt
  [ -s $f ] && continue
  t0=$(date +%s.%N); timeout 3600 nice -19 $S/bin/e2ref "$W" 3 reverse > $f 2>/dev/null; echo "$(echo "$(date +%s.%N) - $t0" | bc) s" > $f.time
  echo "$NAME $(tail -1 $f.time)"
done < ${1:-runs_c2/census/knots4_2222.txt}
