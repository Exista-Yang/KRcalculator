#!/bin/bash
# 4-strand (2,2,2,2) gl_4 scan of one (knot, field, modes) with the checkpointed sharded box driver:
# K shards in parallel (1 thread each), checkpoint every 600 s, resumable (rerun = resume), then finish.
# usage: box.sh NAME FIELD MODES(u,r|r) "WORD" [K]
cd /home/user/KRcalculator
S=/tmp/claude-0/-home-user-KRcalculator/7ca6fc3a-b081-5be7-b7f0-2501a115a78c/scratchpad
NAME=$1; F=$2; M=$3; W=$4; K=${5:-4}
TAG=${NAME}_f${F}_$(echo $M | tr -d ,)
D=$S/box/$TAG; O=runs_c2/k4/scan/out; mkdir -p $D $O
[ -f $O/$TAG.json ] && { echo "$TAG done"; exit 0; }
L=$(echo $W | wc -w)
COMMON="--lam 2,2,2,2 --word \"$W\" --n 4 --ns 2,3,4 --field $F --seed 201 --prune 4 --replay --cache ${LIB:-$S/lib2222}/f$F --costs ${COSTS:-tests/data/homprov/costs_cloud.txt} ${SR:---stamp-route} --stages $L --modes $M --mem-mb ${MEMMB:-3000}"
for k in $(seq 0 $((K-1))); do
  ( ulimit -v ${VLIM:-$((14000000 / K))}; eval KRHOM_MEM_CAP_MB=${MEMCAP:-3000} KRHOM_COSTS_ANY_HOST=1 KRHOM_THREADS=1 timeout ${TMO:-86400} $S/bin/sbox box $COMMON --dir $D/s$k --shard $k/$K --ckpt-secs 600 >> $O/$TAG.s$k.out 2>> $O/$TAG.s$k.err; echo "exit $?" >> $O/$TAG.s$k.err ) &
done
wait
DIRS=$(for k in $(seq 0 $((K-1))); do echo -n "$D/s$k,"; done | sed 's/,$//')
( ulimit -v 14000000; eval KRHOM_MEM_CAP_MB=${MEMCAP:-3000} KRHOM_COSTS_ANY_HOST=1 KRHOM_THREADS=4 timeout 20000 $S/bin/sbox finish $COMMON --dirs $DIRS --json $O/$TAG.json > $O/$TAG.finish.out 2> $O/$TAG.finish.err; echo "exit $?" >> $O/$TAG.finish.err )
echo "$(date +%T) $TAG $(tail -1 $O/$TAG.finish.err)"
