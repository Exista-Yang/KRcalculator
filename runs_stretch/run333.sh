#!/bin/bash
# (3,3,3) at prune 4, N = 4 box, closures N' = 3, 4, both modes; library replayed from data-gl6/m5b
cd /home/user/KRcalculator
for F in 31 32b; do for SEED in 201 202; do
 for W in "1 -2 1 -2" "1 2" "1 1 2" "1 1 1 2" "1 1 1 1 2" "1 1 1 1 1 2"; do
  TAG=$(echo "$W" | tr ' ' '_')_f${F}_seed${SEED}
  KRHOM_COSTS_ANY_HOST=1 KRHOM_THREADS=1 timeout 3000 /tmp/claude-0/-home-user-KRcalculator/578ee88c-5460-54e6-9ecb-80a93fee51bf/scratchpad/sbox_g2b --lam 3,3,3 --word "$W" --n 4 --ns 3,4 --field $F --seed $SEED --prune 4 --replay --cache /tmp/claude-0/-home-user-KRcalculator/578ee88c-5460-54e6-9ecb-80a93fee51bf/scratchpad/lib333/f$F --costs tests/data/homprov/costs_333.txt --scaling --json runs_stretch/$TAG.json > runs_stretch/$TAG.out 2> runs_stretch/$TAG.err
  echo "exit $?" >> runs_stretch/$TAG.err
 done
done; done
