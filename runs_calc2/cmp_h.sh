#!/bin/bash
# compare runs_calc2/h/*.out with runs_calc/*.out, ignoring the "time" line and the timing in the exit line
cd "$(dirname "$0")/.."
norm() { grep -v '^time ' "$1" | sed 's/^exit \([0-9]*\).*/exit \1/'; }
for f in runs_calc2/h/*.out; do
  b=$(basename $f)
  if [ ! -f runs_calc/$b ]; then echo "$b: no CALC-1 reference"; continue; fi
  if diff <(norm runs_calc/$b) <(norm $f) > /dev/null; then echo "$b: identical"; else echo "$b: DIFFERS"; diff <(norm runs_calc/$b) <(norm $f) | head -10; fi
done
