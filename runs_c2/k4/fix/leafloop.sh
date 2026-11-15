#!/bin/bash
# Leaves calibration loop for one pair (C2-5): run the job with m4c-lf; while it is refused on a
# Lemma A step call "> budget 100000 and no uncensored exact-shape cost record" of <= 400 000
# unknowns, measure that call (calib --windows, limit 125 s) and retry. usage: leafloop.sh FIELD JOB TABLE MAXIT
cd /home/user/KRcalculator
F=$1; J="$2"; T=$3; MAX=${4:-8}; D=runs_c2/k4/fix/loop_$(echo "$J" | tr ' ' '_')
mkdir -p $D
for i in $(seq 1 $MAX); do
  LIB=lib2222m JOBS=1 runs_c2/k4/fix/trace.sh $F $T $D/it$i m4c-lf -- "$J"
  L=$(ls $D/it$i/*.log)
  if grep -q '^RESULT {"ok":true' $L; then echo "$J certified at iteration $i"; exit 0; fi
  W=$(grep -E '^\[leafpairs\].*no uncensored exact-shape cost record' $L | tail -1 | python3 -c '
import re,sys
l=sys.stdin.read()
m=re.search(r"step call \[([\d, ]+)\] \[([^\]]*)\] -> \[([^\]]*)\] d=(-?\d+): (\d+) unknowns",l)
if m and int(m.group(5))<=400000: print(m.group(1).replace(" ",""),m.group(2),m.group(3),m.group(4),"unknowns="+m.group(5))')
  if [ -z "$W" ]; then echo "$J: stops at iteration $i: $(grep -E '^\[leafpairs\]' $L | tail -1 | cut -c1-300)"; exit 1; fi
  echo "$W" > $D/win$i.txt; echo "iteration $i: calibrate $W"
  cp $T $D/base$i.txt
  runs_c2/k4/fix/calib.sh $D/win$i.txt $D/base$i.txt 125 4000 > $D/calib$i.log 2>&1
  grep "^\[calib\] \[" $D/calib$i.log | cut -c1-200
  T=runs_c2/k4/fix/costs_2222_c25.txt
done
