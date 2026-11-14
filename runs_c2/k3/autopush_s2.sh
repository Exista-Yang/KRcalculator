#!/bin/bash
cd /home/user/KRcalculator
while true; do
  sleep 1800
  git add runs_c2/k3 docs/colored_design/c2 2>/dev/null
  git commit -q -m "C2-S2: slice2 results checkpoint ($(ls runs_c2/k3/out_s2/*.json 2>/dev/null | wc -l) runs)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" && git push -q origin colored-2c4-k3-s2
done
