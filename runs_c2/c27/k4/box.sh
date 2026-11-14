#!/bin/bash
# C2-7: one (knot, field) of (2,2,2,2) gl_4 with the sharded box driver: prune 4, box N = 4,
# closures N' = 2,3,4, modes u,r, final stage; library from layer A only (--replay, provider refusing
# every call, M4c-lf stamp, route stamp). K shards in parallel (resumable: rerun the same command),
# then finish (χ = decat_chi and the D14 scaling check are mandatory inside sbox finish).
# usage: box.sh NAME FIELD "WORD" [K]
cd /home/user/KRcalculator
NAME=$1; F=$2; W=$3; K=${4:-4}
D=${OUTD:-runs_c2/c27/k4/out}; mkdir -p $D
TAG=${NAME}_f$F
[ -f $D/$TAG.json ] && exit 0
L=$(echo $W | wc -w)
ulimit -v ${VLIM:-14000000}
C="--lam 2,2,2,2 --word \"$W\" --n 4 --ns 2,3,4 --field $F --seed ${SEED:-201} --prune 4 --mem-mb ${MEMMB:-3000} --replay --cache ${CACHE:-/home/user/work/cache/f$F} --costs /home/user/work/costs_2222_c25.txt --provider m4c-lf --stamp-route --stages $L"
export KRHOM_MEM_CAP_MB=${MEMCAP:-200000} KRHOM_COSTS_ANY_HOST=1 KRHOM_THREADS=${THREADS:-1}
t0=$(date +%s)
dirs=""
for k in $(seq 0 $((K-1))); do
  SD=$D/${TAG}_s$k; mkdir -p $SD; dirs="$dirs,$SD"
  eval timeout ${TMO:-36000} ./target/release/sbox box $C --dir $SD --shard $k/$K --ckpt-secs ${CKPT:-600} > $SD/box.out 2> $SD/box.err &
done
wait
dirs=${dirs#,}
eval timeout 7200 ./target/release/sbox finish $C --dirs $dirs --json $D/$TAG.json > $D/$TAG.out 2> $D/$TAG.err
echo "exit $?" >> $D/$TAG.err
echo "$(date +%T) $TAG $(tail -1 $D/$TAG.err) $(( $(date +%s) - t0 ))s"
