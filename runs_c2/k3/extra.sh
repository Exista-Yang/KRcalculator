#!/bin/bash
# extra checks: unpruned library (level 6) on small knots; Markov words sigma_1^m sigma_2^{+-1}
cd /home/user/KRcalculator
S=/tmp/claude-0/-home-user-KRcalculator/934ab925-996a-59b6-8fc2-b741a70e6e9a/scratchpad
ulimit -v 3400000
run() { # DIR TAG LIB PRUNE FIELD WORD
  [ -f $1/$2.json ] && return
  L=$(echo $6 | wc -w)
  KRHOM_COSTS_ANY_HOST=1 KRHOM_THREADS=1 timeout 7200 $S/bin/sbox --lam 2,2,2 --word "$6" --n 4 --ns 2,3,4 --field $5 --seed 201 --prune $4 --replay --cache $3/f$5 --costs tests/data/homprov/costs_cloud.txt --stages $L --json $1/$2.json > $1/$2.out 2> $1/$2.err
  echo "exit $?" >> $1/$2.err; echo "$2 $(tail -n1 $1/$2.err)"
}
for F in 31 32b; do
  run runs_c2/k3/markov s3p_f$F $S/lib222 4 $F "1 1 1 2"
  run runs_c2/k3/markov s3n_f$F $S/lib222 4 $F "1 1 1 -2"
  run runs_c2/k3/markov s5p_f$F $S/lib222 4 $F "1 1 1 1 1 2"
  run runs_c2/k3/markov s5n_f$F $S/lib222 4 $F "1 1 1 1 1 -2"
  run runs_c2/k3/markov s2p_f$F $S/lib222 4 $F "1 1 2"
  run runs_c2/k3/markov s4p_f$F $S/lib222 4 $F "1 1 1 1 2"
done
run runs_c2/k3/unpruned 4_1_f32b $S/lib222_n6 6 32b "1 -2 1 -2"
for K in "5_2;1 1 1 2 -1 2" "6_2;1 1 1 -2 1 -2" "6_3;1 1 -2 1 -2 -2" "7_3;1 1 1 1 1 2 -1 2"; do
  run runs_c2/k3/unpruned ${K%%;*}_f31 $S/lib222_n6 6 31 "${K#*;}"
done
