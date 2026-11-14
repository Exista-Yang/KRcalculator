#!/bin/bash
# C2-S1 per (knot, field): sbox box (checkpointed, resumable, ulimit -v 6 GB, 1 thread) then sbox finish
# (closures + chi + scaling; serialized by a global lock, ulimit -v 13.5 GB, cap KRHOM_MEM_CAP_MB 40000)
# then conversion of the finish JSON to the single-run record format of gl4_222.json.
# usage: one_s1.sh NAME FIELD WORD
cd /home/user/KRcalculator
S=${S:-/tmp/claude-0/-home-user-KRcalculator/3a18306a-1d5f-50c6-9aa2-7362a60f64a8/scratchpad}
NAME=$1; F=$2; W=$3
D=${K3DIR:-runs_c2/k3/out_s1}; mkdir -p $D
TAG=${NAME}_f${F}
[ -f $D/$TAG.json ] && exit 0
rm -f "$D/$TAG.fail"
L=$(echo $W | wc -w)
COMMON="--lam 2,2,2 --word \"$W\" --n 4 --ns 2,3,4 --field $F --seed 201 --prune 4 --mem-mb 12000 --replay --cache $S/lib222/f$F --costs tests/data/homprov/costs_cloud.txt --stages $L"
export KRHOM_COSTS_ANY_HOST=1 KRHOM_THREADS=1
if ! [ -f $D/$TAG.box.done ]; then
  ( ulimit -v ${BOXULIM:-10000000}; eval timeout 172800 $S/bin/sbox box $COMMON --dir $D/$TAG.d --ckpt-secs 300 ) > $D/$TAG.box.out 2>> $D/$TAG.err
  rc=$?; echo "box exit $rc" >> $D/$TAG.err
  [ $rc = 0 ] && touch $D/$TAG.box.done || { echo "box exit $rc" > $D/$TAG.fail; exit 0; }
fi
(
  flock 9
  ( ulimit -v ${FINULIM:-9000000}; eval KRHOM_MEM_CAP_MB=200000 timeout 28800 $S/bin/sbox finish $COMMON --dirs $D/$TAG.d --json $D/$TAG.fin.json ) > $D/$TAG.fin.out 2>> $D/$TAG.err
  echo "finish exit $?" >> $D/$TAG.err
) 9> /tmp/s1_finish.lock
if [ -f $D/$TAG.fin.json ]; then
  python3 tools/c2/finconv.py $D/$TAG.fin.json $D/$TAG.err $D/$TAG.json && echo "$TAG ok"
else
  tail -n 2 $D/$TAG.err > $D/$TAG.fail; echo "$TAG FAILED $(tail -n 1 $D/$TAG.err)"
fi
