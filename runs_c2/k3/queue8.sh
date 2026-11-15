#!/bin/bash
# C2-6: remaining <=10-crossing knots one at a time, one prime at a time, cap raised to 500 GB (estimate), real guard ulimit 14 GB
cd /home/user/KRcalculator
export MEMCAP=500000 TMO=40000
for spec in "10_91;-1 -1 -1 2 -1 2 2 -1 2 2" "10_112;-1 -1 -1 2 -1 2 -1 2 -1 2" "10_116;-1 -1 2 -1 -1 2 -1 2 -1 2"; do
  N=${spec%%;*}; W=${spec#*;}
  for F in 31 32b; do ULIM=14500000 MEMMB=14000 THREADS=3 runs_c2/k3/one2.sh $N $F "$W"; done
done
