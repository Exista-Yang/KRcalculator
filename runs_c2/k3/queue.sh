#!/bin/bash
# all knots of tools/ntrunc_t1/data/knots3_all.txt sorted by braid length, both primes, P parallel workers
cd /home/user/KRcalculator
awk -F';' '{n=split($3,a," "); print n";"$1";"$3}' tools/ntrunc_t1/data/knots3_all.txt | sort -t';' -k1,1n -s | while IFS=';' read L NAME W; do
  [ "$L" -gt "${MAXLEN:-99}" ] && continue
  for F in 31 32b; do printf '%s\t%s\t%s\n' "$NAME" "$F" "$W"; done
done | xargs -d '\n' -P ${P:-2} -n 1 bash -c 'IFS=$'"'"'\t'"'"' read -r a b c <<< "$0"; runs_c2/k3/one.sh "$a" "$b" "$c"'
