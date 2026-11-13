#!/bin/bash
# (3,3,3) knots of tools/ntrunc_t1/data/knots3_all.txt at prune 4, N = 4 box, closures N' = 3, 4,
# both modes, final stage; library from data-gl6/m5c (layer A only). Resumable per (knot, field).
cd /home/user/KRcalculator
S=/tmp/claude-0/-home-user-KRcalculator/d808422c-dfa4-592e-8d01-036c61ee01fe/scratchpad
ulimit -v 4000000
MAXLEN=${MAXLEN:-99}
awk -F';' '{n=split($3,a," "); print n";"$0}' tools/ntrunc_t1/data/knots3_all.txt | sort -t';' -k1,1n -s | while IFS=';' read L NAME LAM W NS; do
  [ "$L" -gt "$MAXLEN" ] && continue
  for F in 31 32b; do
    TAG=${NAME}_f${F}
    [ -f ${K3DIR:-runs_m6a1c/k3}/$TAG.json ] && continue
    [ -f ${K3DIR:-runs_m6a1c/k3}/$TAG.err ] && grep -q "^exit" ${K3DIR:-runs_m6a1c/k3}/$TAG.err && continue
    KRHOM_COSTS_ANY_HOST=1 KRHOM_THREADS=${THREADS:-1} timeout 3600 $S/sbox_c --lam 3,3,3 --word "$W" --n 4 --ns 3,4 --field $F --seed ${SEED:-201} --prune 4 --replay --stamp-route --cache $S/lib333/f$F --costs tests/data/homprov/costs_333_m5c.txt --scaling --stages $L --json ${K3DIR:-runs_m6a1c/k3}/$TAG.json > ${K3DIR:-runs_m6a1c/k3}/$TAG.out 2> ${K3DIR:-runs_m6a1c/k3}/$TAG.err
    echo "exit $?" >> ${K3DIR:-runs_m6a1c/k3}/$TAG.err
    echo "$TAG $(tail -1 ${K3DIR:-runs_m6a1c/k3}/$TAG.err)"
  done
done
