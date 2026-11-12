#!/bin/bash
# usage: runw.sh M FIELD SEED [extra sbox flags]   (σ^M on (3,3), N = 6 box, N' = 3..6, both modes, scaling check)
cd /home/user/KRcalculator
M=$1; FLD=$2; SEED=$3; shift 3
W=$(python3 -c "print(' '.join(['1']*$M))")
TAG=s${M}_f${FLD}_seed${SEED}
ulimit -v 14000000
KRHOM_THREADS=${THREADS:-4} timeout 14000 /tmp/claude-0/-home-user-KRcalculator/578ee88c-5460-54e6-9ecb-80a93fee51bf/scratchpad/sbox_g2b --lam 3,3 --word "$W" --n 6 --ns 3,4,5,6 --field $FLD --seed $SEED --cache /tmp/claude-0/-home-user-KRcalculator/578ee88c-5460-54e6-9ecb-80a93fee51bf/scratchpad/cache33_f$FLD --scaling --json runs_g2b/$TAG.json "$@" > runs_g2b/$TAG.out 2> runs_g2b/$TAG.err
echo "exit $?" >> runs_g2b/$TAG.err
