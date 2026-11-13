#!/bin/bash
# one (knot, field) of the (2,2,2) gl_4 queue: prune 4, box N = 4, closures N' = 2,3,4, both modes,
# final stage; library from layer A only (--replay). Resumable: skips finished (knot, field).
# usage: one.sh NAME FIELD WORD
cd /home/user/KRcalculator
S=/tmp/claude-0/-home-user-KRcalculator/934ab925-996a-59b6-8fc2-b741a70e6e9a/scratchpad
NAME=$1; F=$2; W=$3
D=${K3DIR:-runs_c2/k3/out}; mkdir -p $D
TAG=${NAME}_f${F}
[ -f $D/$TAG.json ] && exit 0
[ -f $D/$TAG.err ] && grep -q "^exit" $D/$TAG.err && exit 0
L=$(echo $W | wc -w)
ulimit -v 6000000
KRHOM_MEM_CAP_MB=${MEMCAP:-3000} KRHOM_COSTS_ANY_HOST=1 KRHOM_THREADS=${THREADS:-2} timeout ${TMO:-7200} $S/bin/sbox --lam 2,2,2 --word "$W" --n 4 --ns 2,3,4 --field $F --seed ${SEED:-201} --prune 4 --mem-mb ${MEMMB:-3200} --replay --cache ${LIB:-$S/lib222}/f$F --costs tests/data/homprov/costs_cloud.txt --stages $L --json $D/$TAG.json > $D/$TAG.out 2> $D/$TAG.err
echo "exit $?" >> $D/$TAG.err
echo "$TAG $(tail -1 $D/$TAG.err) $(grep -o 'scan [0-9.]*s' $D/$TAG.err)"
