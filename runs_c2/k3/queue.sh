#!/bin/bash
# all knots of tools/ntrunc_t1/data/knots3_all.txt in order of the truncated envelope's peak
# Hom dimension at N = 4 (runs_c2/k3/order.txt, a cost proxy), both primes, P parallel workers
cd /home/user/KRcalculator
while IFS=';' read NAME W PK; do
  for F in 31 32b; do printf '%s\t%s\t%s\n' "$NAME" "$F" "$W"; done
done < runs_c2/k3/order.txt | xargs -d '\n' -P ${P:-2} -n 1 bash -c 'IFS=$'"'"'\t'"'"' read -r a b c <<< "$0"; runs_c2/k3/one.sh "$a" "$b" "$c"'
