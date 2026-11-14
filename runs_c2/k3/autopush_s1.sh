#!/bin/bash
# commit + push run outputs to colored-2c4-k3-s1 every 30 min; also log RSS every minute
cd /home/user/KRcalculator
i=0
while true; do
  sleep 60; i=$((i+1))
  echo "$(date +%T) $(ps -C sbox -o rss= | awk '{s+=$1; if($1>m)m=$1} END{print NR" jobs total_rss_MB="int(s/1024)" max_MB="int(m/1024)}')" >> runs_c2/k3/rss_s1.log
  if [ $((i%30)) = 0 ]; then
    git add runs_c2/k3 2>/dev/null
    git commit -q -m "C2-S1: slice1 results checkpoint ($(ls runs_c2/k3/out_s1/*.json 2>/dev/null | wc -l) runs)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" && git push -q origin colored-2c4-k3-s1
  fi
done
