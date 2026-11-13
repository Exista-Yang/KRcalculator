#!/bin/bash
# rerun the failed non-M5i jobs of a libbuild run file with the calibrated table (route stamp), one
# child per job (libbuild one), 4 in parallel; RESULT lines -> runs_c2/k4/lib/l2222_n4_fF_fix.jsonl
# usage: fix.sh FIELD RUNFILE...
cd /home/user/KRcalculator
S=/tmp/claude-0/-home-user-KRcalculator/7ca6fc3a-b081-5be7-b7f0-2501a115a78c/scratchpad
F=$1; shift
O=runs_c2/k4/lib/l2222_n4_f${F}_fix.jsonl
python3 - "$@" > $S/fixjobs_$F.txt <<'PY'
import json, sys
best = {}
for f in sys.argv[1:]:
    for l in open(f):
        r = json.loads(l)
        if 'attempt' in r:
            j = r['job']; o = r['outcome']
            if best.get(j) != 'ok': best[j] = o
for j, o in best.items():
    if o not in ('ok', 'NeedsIdempotent'): print(j)
PY
echo "$(wc -l < $S/fixjobs_$F.txt) jobs"
one() {
  J="$1"; F=$2; O=$3; S=$4
  CORNER=""; case "$J" in cross*) CORNER="--corner markov-ext";; esac
  t0=$(date +%s)
  R=$( (ulimit -v 6000000; KRHOM_COSTS_ANY_HOST=1 KRHOM_THREADS=1 timeout 1800 $S/bin/libbuild one --lam 2,2,2,2 --n 4 --field $F --cache $S/lib2222/f$F --costs runs_c2/k4/calib/costs_2222.txt --stamp-route $CORNER -- $J 2>/dev/null) | grep '^RESULT' | cut -c8-)
  echo "{\"job\":\"$J\",\"exit_wall\":$(( $(date +%s) - t0 )),\"result\":${R:-null}}" >> $O
  echo "$J -> $(echo "$R" | grep -o '"ok":[a-z]*' | head -1)"
}
export -f one
cat $S/fixjobs_$F.txt | xargs -P 4 -I{} bash -c "one '{}' $F $O $S"
