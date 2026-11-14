#!/bin/bash
# C2-S1: slice1 queue, 4 jobs x 1 thread, both primes per knot (adjacent in the queue).
cd /home/user/KRcalculator
while IFS=';' read NAME W PK; do
  for F in 31 32b; do printf '%s\t%s\t%s\n' "$NAME" "$F" "$W"; done
done < runs_c2/k3/slices/slice1.txt | xargs -d '\n' -P ${P:-4} -n 1 bash -c 'IFS=$'"'"'\t'"'"' read -r a b c <<< "$0"; runs_c2/k3/one_s1.sh "$a" "$b" "$c"'
