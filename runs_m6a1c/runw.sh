#!/bin/bash
# usage: runw.sh M FIELD SEED MODES [extra sbox flags]  (σ^M on (3,3), N = 6 box, N' = 3..6, scaling check)
cd /home/user/KRcalculator
S=/tmp/claude-0/-home-user-KRcalculator/d808422c-dfa4-592e-8d01-036c61ee01fe/scratchpad
M=$1; FLD=$2; SEED=$3; MODES=$4; shift 4
W=$(python3 -c "print(' '.join(['1']*$M))")
TAG=s${M}$( [ "$MODES" = r ] && echo r )_f${FLD}_seed${SEED}
[ -f runs_m6a1c/g2b/$TAG.json ] && exit 0
[ -f runs_m6a1c/g2b/$TAG.err ] && exit 0
ulimit -v 14000000
KRHOM_THREADS=${THREADS:-3} timeout 30000 $S/sbox_c --lam 3,3 --word "$W" --n 6 --ns 3,4,5,6 --field $FLD --seed $SEED --modes $MODES --mem-mb ${MEMMB:-3500} --cache $S/cache33_f$FLD --scaling --json runs_m6a1c/g2b/$TAG.json "$@" > runs_m6a1c/g2b/$TAG.out 2> runs_m6a1c/g2b/$TAG.err
echo "exit $?" >> runs_m6a1c/g2b/$TAG.err
