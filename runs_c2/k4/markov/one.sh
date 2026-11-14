#!/bin/bash
# one 4-strand (2,2,2,2) word, one-word in-memory box loop (box N = 4, closures N' = 2,3,4, u and r),
# library from layer A only (--replay, calibrated table, route stamp). usage: one.sh TAG FIELD "WORD"
cd /home/user/KRcalculator
S=/tmp/claude-0/-home-user-KRcalculator/7ca6fc3a-b081-5be7-b7f0-2501a115a78c/scratchpad
TAG=$1; F=$2; W=$3; D=${OUTD:-runs_c2/k4/markov}
[ -f $D/${TAG}_f$F.json ] && exit 0
L=$(echo $W | wc -w)
ulimit -v ${VLIM:-14000000}
t0=$(date +%s)
KRHOM_MEM_CAP_MB=${MEMCAP:-4000} KRHOM_COSTS_ANY_HOST=1 KRHOM_THREADS=${THREADS:-4} timeout ${TMO:-14400} $S/bin/sbox --lam 2,2,2,2 --word "$W" --n 4 --ns 2,3,4 --field $F --seed 201 --prune 4 --mem-mb ${MEMMB:-3000} --replay --cache $S/lib2222/f$F --costs runs_c2/k4/calib/costs_2222.txt --stamp-route --stages $L --json $D/${TAG}_f$F.json > $D/${TAG}_f$F.out 2> $D/${TAG}_f$F.err
echo "exit $?" >> $D/${TAG}_f$F.err
echo "$(date +%T) ${TAG}_f$F $(tail -1 $D/${TAG}_f$F.err) $(( $(date +%s) - t0 ))s"
