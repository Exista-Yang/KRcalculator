#!/bin/bash
# C2-6: REVERSED remaining queue (second worker pair) (order_le10 first, then order.txt), raised memory cap, one2.sh skips finished runs
cd /home/user/KRcalculator
export ULIM=${ULIM:-8000000} MEMCAP=${MEMCAP:-40000} MEMMB=${MEMMB:-8000} THREADS=${THREADS:-1} TMO=${TMO:-28000}
( cat runs_c2/k3/order_le10.txt; cat runs_c2/k3/order.txt ) | awk -F';' '!s[$1]++' | tac | grep -v -E '^(8_18|10_85|10_100|10_148|10_155);' | while IFS=';' read NAME W PK; do
  for F in 31 32b; do printf '%s\t%s\t%s\n' "$NAME" "$F" "$W"; done
done | xargs -d '\n' -P ${P:-2} -n 1 bash -c 'IFS=$'"'"'\t'"'"' read -r a b c <<< "$0"; runs_c2/k3/one2.sh "$a" "$b" "$c"'
