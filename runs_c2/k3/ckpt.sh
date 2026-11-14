#!/bin/bash
cd /home/user/KRcalculator
python3 tools/c2/k3_c2_compare.py > runs_c2/k3/compare_last.txt 2>&1
git add runs_c2/k3 docs/colored_design/c2 2>/dev/null
git commit -qm "C2-6: checkpoint ($(ls runs_c2/k3/out/*.json | wc -l) runs)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" && git push -q origin colored-2c4-k3
