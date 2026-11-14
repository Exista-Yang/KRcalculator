#!/bin/bash
# C2-3: like queue.sh but reads ORDER (default slices/slice2.txt); one_s2.sh skips finished (knot, field).
cd /home/user/KRcalculator
while IFS=';' read NAME W PK; do
  for F in 31 32b; do printf '%s\t%s\t%s\n' "$NAME" "$F" "$W"; done
done < ${ORDER:-runs_c2/k3/slices/slice2.txt} | xargs -d '\n' -P ${P:-3} -n 1 bash -c 'IFS=$'"'"'\t'"'"' read -r a b c <<< "$0"; runs_c2/k3/one_s2.sh "$a" "$b" "$c"'
