#!/bin/bash
# M5i library jobs of (2,2,2,2) prune 4 (the 43 NeedsIdempotent jobs of C2-2), one `libbuild one`
# child per job, on the restored C2-2 layer A (data-gl6 c2/l2222_n4) + the calibrated cost table.
# Resumable: a job with an ok RESULT line in OUT is skipped. usage: runjobs.sh FIELD [JOBSFILE] [MODE]
# MODE = run (provider on) | replay (provider refusing every call)
cd /home/user/KRcalculator
F=$1; JF=${2:-runs_c2/k4m/lib/jobs_m5i.txt}; MODE=${3:-run}
D=runs_c2/k4m/lib; OUT=$D/${TAGP:-m5i}_f${F}_${MODE}.jsonl
CACHE=${CACHE:-/home/user/lib2222c/f$F}
RF=""; [ $MODE = replay ] && RF="--replay"
one() {
  J="$1"
  grep -F "\"job\":\"$J\"" $OUT 2>/dev/null | grep -q '"ok":true' && return
  t0=$(date +%s)
  R=$( (ulimit -v ${VLIM:-7000000}; KRHOM_COSTS_ANY_HOST=1 KRHOM_THREADS=1 timeout ${WALL:-3600} ./target/release/libbuild one --lam 2,2,2,2 --n 4 --field $F --cache $CACHE --costs /home/user/lib2222c/costs_2222_c25.txt --provider m4c-lf --stamp-route $RF -- $J 2>>$D/${TAGP:-m5i}_f${F}_${MODE}.stderr) | grep '^RESULT' | sed 's/^RESULT //')
  [ -z "$R" ] && R="{\"ok\":false,\"errors\":[{\"what\":\"child\",\"kind\":\"Killed\",\"text\":\"no RESULT (timeout or memory)\"}]}"
  echo "{\"job\":\"$J\",\"exit_wall\":$(( $(date +%s) - t0 )),\"result\":$R}" >> $OUT
  echo "$(date +%T) $J $(echo $R | grep -o '"ok":[a-z]*') $(( $(date +%s) - t0 ))s"
}
export -f one; export F D OUT CACHE RF MODE TAGP
cat $JF | xargs -P ${PAR:-2} -I{} bash -c 'one "{}"'
