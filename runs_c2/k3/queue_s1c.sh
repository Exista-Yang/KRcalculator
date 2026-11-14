#!/bin/bash
# third pass: everything not yet done, P=2, memory-heavy settings; wait until <= 1 old job is running
cd /home/user/KRcalculator
export S=/tmp/claude-0/-home-user-KRcalculator/3a18306a-1d5f-50c6-9aa2-7362a60f64a8/scratchpad
while [ $(pgrep -fc "one_s1(_heavy)?.sh [^ ]* [^ ]* ") -gt 2 ]; do sleep 30; done
while IFS=';' read NAME W PK; do
  for F in 31 32b; do printf '%s\t%s\t%s\n' "$NAME" "$F" "$W"; done
done < runs_c2/k3/slices/slice1.txt | xargs -d '\n' -P 2 -n 1 bash -c 'IFS=$'"'"'\t'"'"' read -r a b c <<< "$0"; runs_c2/k3/one_s1_heavy.sh "$a" "$b" "$c"'
