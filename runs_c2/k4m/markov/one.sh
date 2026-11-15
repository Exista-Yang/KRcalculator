#!/bin/bash
# validation (b): a 4-strand (2,2,2,2) word with the C2-4 binary (M5i engine code), library from layer A
# only (--replay) of the C2-5 cache (M4c-lf stamp) plus C2-4's M5i additions; compared with C2-2's
# runs_c2/k4/markov results. usage: one.sh TAG FIELD "WORD"
cd /home/user/KRcalculator
TAG=$1; F=$2; W=$3; D=${OUTD:-runs_c2/k4m/markov}
[ -f $D/${TAG}_f$F.json ] && exit 0
L=$(echo $W | wc -w)
ulimit -v ${VLIM:-14000000}
t0=$(date +%s)
KRHOM_MEM_CAP_MB=${MEMCAP:-200000} KRHOM_COSTS_ANY_HOST=1 KRHOM_THREADS=${THREADS:-2} timeout ${TMO:-14400} ./target/release/sbox --lam 2,2,2,2 --word "$W" --n 4 --ns 2,3,4 --field $F --seed 201 --prune 4 --mem-mb ${MEMMB:-3000} --replay --cache ${CACHE:-/home/user/lib2222c/f$F} --costs /home/user/lib2222c/costs_2222_c25.txt --provider m4c-lf --stamp-route --stages $L --json $D/${TAG}_f$F.json > $D/${TAG}_f$F.out 2> $D/${TAG}_f$F.err
echo "exit $?" >> $D/${TAG}_f$F.err
echo "$(date +%T) ${TAG}_f$F $(tail -1 $D/${TAG}_f$F.err) $(( $(date +%s) - t0 ))s"
