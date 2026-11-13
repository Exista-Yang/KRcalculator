#!/bin/bash
# C2-3: like queue.sh but reads ORDER (default order_le10.txt); one.sh skips finished (knot, field).
cd /home/user/KRcalculator
while IFS=';' read NAME W PK; do
  for F in 31 32b; do printf '%s\t%s\t%s\n' "$NAME" "$F" "$W"; done
done < ${ORDER:-runs_c2/k3/order_le10.txt} | xargs -d '\n' -P ${P:-3} -n 1 bash -c 'IFS=$'"'"'\t'"'"' read -r a b c <<< "$0"; runs_c2/k3/one.sh "$a" "$b" "$c"'
