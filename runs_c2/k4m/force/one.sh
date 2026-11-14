#!/bin/bash
# M5i engine validation (a): a (2,2,2) braid word with forced non-canonical *decomposable* ladders
# (KRHOM_FORCE_LADDER) for some types, so the engine runs through e_p and reduced coordinates;
# compared with the canonical-ladder results (C2-3, runs_c2/k3/out). Provider allowed (new keys).
# usage: one.sh TAG FIELD "WORD" "FORCE"
cd /home/user/KRcalculator
S=/tmp/claude-0/-home-user-KRcalculator/699a7aa7-71b5-5656-a90b-8ab1bea65187/scratchpad
TAG=$1; F=$2; W=$3; FORCE=$4; D=runs_c2/k4m/force
[ -f $D/${TAG}_f$F.json ] && exit 0
L=$(echo $W | wc -w)
ulimit -v 6000000
KRHOM_FORCE_LADDER="$FORCE" KRHOM_MEM_CAP_MB=20000 KRHOM_COSTS_ANY_HOST=1 KRHOM_THREADS=${THREADS:-2} timeout ${TMO:-7200} ./target/release/sbox --lam 2,2,2 --word "$W" --n 4 --ns 2,3,4 --field $F --seed 201 --prune 4 --mem-mb 3200 --cache $S/lib222f_f$F --costs tests/data/homprov/costs_cloud.txt --stamp-route --stages $L --json $D/${TAG}_f$F.json > $D/${TAG}_f$F.out 2> $D/${TAG}_f$F.err
echo "exit $?" >> $D/${TAG}_f$F.err
echo "${TAG}_f$F $(tail -1 $D/${TAG}_f$F.err)"
