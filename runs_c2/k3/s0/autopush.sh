#!/bin/bash
# commit + push slice-0 outputs to colored-2c4-k3-s0 every 30 min
cd /home/user/KRcalculator
while true; do
  sleep 1800
  git add runs_c2/k3/s0 runs_c2/k3/e2ref tools/c2/k3_c2_compare_s0.py 2>/dev/null
  git commit -q -m "C2-S0: slice0 results checkpoint ($(ls runs_c2/k3/s0/out/*.json 2>/dev/null | wc -l) runs)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" && git push -q origin colored-2c4-k3-s0
done
