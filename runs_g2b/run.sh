#!/bin/bash
# usage: run.sh FIELD SEED TAG  (σ^9 on (3,3), N = 6 box, closures N' = 3..6, both modes, scaling check)
cd /home/user/KRcalculator
ulimit -v 14000000
KRHOM_THREADS=${THREADS:-4} timeout 18000 ./target/release/sbox --lam 3,3 --word "1 1 1 1 1 1 1 1 1" --n 6 --ns 3,4,5,6 --field $1 --seed $2 --cache /tmp/claude-0/-home-user-KRcalculator/578ee88c-5460-54e6-9ecb-80a93fee51bf/scratchpad/cache33_f$1 --scaling --json runs_g2b/s9_f$1_seed$2.json > runs_g2b/s9_f$1_seed$2.out 2> runs_g2b/s9_f$1_seed$2.err
echo "exit $?" >> runs_g2b/s9_f$1_seed$2.err
