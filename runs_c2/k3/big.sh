#!/bin/bash
# C2-6: the 5 refused knots, raised memory cap, 2 jobs at a time
cd /home/user/KRcalculator
export ULIM=12000000 MEMCAP=${MEMCAP:-40000} MEMMB=12000 THREADS=2 TMO=${TMO:-28000}
while IFS=';' read NAME W PK; do
  for F in 31 32b; do printf '%s\t%s\t%s\n' "$NAME" "$F" "$W"; done
done < runs_c2/k3/refused5.txt | xargs -d '\n' -P 2 -n 1 bash -c 'IFS=$'"'"'\t'"'"' read -r a b c <<< "$0"; runs_c2/k3/one2.sh "$a" "$b" "$c"'
