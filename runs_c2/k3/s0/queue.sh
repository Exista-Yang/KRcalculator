#!/bin/bash
# C2-S0: slice0 queue. usage: S=scratch P=4 runs_c2/k3/s0/queue.sh
cd /home/user/KRcalculator
while IFS=';' read NAME W PK; do
  for F in 31 32b; do printf '%s\t%s\t%s\n' "$NAME" "$F" "$W"; done
done < runs_c2/k3/slices/slice0.txt | xargs -d '\n' -P ${P:-4} -n 1 bash -c 'IFS=$'"'"'\t'"'"' read -r a b c <<< "$0"; runs_c2/k3/s0/one.sh "$a" "$b" "$c"'
