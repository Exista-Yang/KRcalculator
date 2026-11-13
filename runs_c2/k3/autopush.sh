#!/bin/bash
# commit + push run outputs to colored-2c4-k3 every 30 min (resumable)
cd /home/user/KRcalculator
while true; do
  sleep 1800
  git add runs_c2/k3 2>/dev/null
  git commit -q -m "C2-3: 3-strand results checkpoint ($(ls runs_c2/k3/out/*.json | wc -l) runs)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" && git push -q origin colored-2c4-k3
done
