#!/bin/bash
# third pass: memory-heavy knots, one job at a time (box ulimit 10 GB, finish cap 200000 MB with ulimit 9 GB)
cd /home/user/KRcalculator
export S=/tmp/claude-0/-home-user-KRcalculator/3a18306a-1d5f-50c6-9aa2-7362a60f64a8/scratchpad
for spec in "12n_767;1 -2 1 1 1 -2 -2 -1 -1 -1 -1 -2" "12a_1051;1 -2 1 1 1 1 1 -2 -2 1 1 -2" "12n_684;-1 -2 -2 -2 -1 2 -1 -1 -1 2 -1 2"; do
  N=${spec%%;*}; W=${spec#*;}
  for F in 31 32b; do runs_c2/k3/one_s1_heavy.sh $N $F "$W"; done
done
