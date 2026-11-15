#!/bin/bash
# CALC-2 item 14: replay (provider refusing every call) of the decomposable-source reroute jobs of
# data/libstatus on a scratch copy of the installed l2222_n4 library, both primes, with the
# libbuild of this branch (reroute_family returns None for a decomposable prefix type).
cd "$(dirname "$0")/.."
SCR=${SCR:-/tmp/claude-0/calc2_replay}
LIB=${KRHOM_LIB_DIR:-$HOME/.cache/krhom/lib}/l2222_n4
OUT=runs_calc2/replay14.jsonl
: > $OUT
for F in 31 32b; do
  rm -rf $SCR/f$F; mkdir -p $SCR; cp -a $LIB/f$F $SCR/f$F
  for J in 109 1789 1812 2052 2086 2275 2298; do
    t0=$(date +%s)
    R=$( (ulimit -v 14000000; KRHOM_COSTS_ANY_HOST=1 KRHOM_THREADS=4 timeout 7200 ./target/release/libbuild one --lam 2,2,2,2 --n 4 --field $F --cache $SCR/f$F --costs $SCR/f$F/costs.txt --provider m4c-lf --stamp-route --replay -- pair $J 2>>runs_calc2/replay14_f$F.stderr) | grep '^RESULT' | sed 's/^RESULT //')
    [ -z "$R" ] && R='{"ok":false,"errors":[{"what":"child","kind":"Killed","text":"no RESULT"}]}'
    echo "{\"field\":\"$F\",\"job\":\"pair $J\",\"wall\":$(( $(date +%s) - t0 )),\"result\":$R}" >> $OUT
    echo "$(date +%T) F$F pair $J $(echo $R | grep -o '"ok":[a-z]*' | head -1) $(( $(date +%s) - t0 ))s"
  done
done
echo ALLDONE
