#!/bin/bash
# C2-7 library jobs of (2,2,2,2) prune 4, one `libbuild one` child per job, on the restored layer A
# (data-gl6 c2/l2222_n4 + c2/l2222_n4_c25 + c2/l2222_n4_c24) and costs_2222_c25.txt.
# KRHOM_LEAVES_C27=1 (joint junction words + Leaves pairs with e_q) unless C27=0.
# Resumable: a job with an ok RESULT line in OUT is skipped.
# usage: runjobs.sh FIELD JOBSFILE [MODE]   MODE = run | replay (provider refusing every call)
cd /home/user/KRcalculator
F=$1; JF=$2; MODE=${3:-run}
D=runs_c2/c27/lib; OUT=$D/${TAGP:-c27}_f${F}_${MODE}.jsonl
CACHE=${CACHE:-/home/user/work/cache/f$F}
RF=""; [ $MODE = replay ] && RF="--replay"
one() {
  J="$1"
  grep -F "\"job\":\"$J\"" $OUT 2>/dev/null | grep -q '"ok":true' && return
  t0=$(date +%s)
  R=$( (ulimit -v ${VLIM:-14000000}; KRHOM_LEAVES_C27=${C27:-1} KRHOM_COSTS_ANY_HOST=1 KRHOM_THREADS=${THR:-1} timeout ${WALL:-7200} ./target/release/libbuild one --lam 2,2,2,2 --n 4 --field $F --cache $CACHE --costs ${COSTS:-/home/user/work/costs_2222_c25.txt} --provider m4c-lf --stamp-route $RF -- $J 2>>$D/${TAGP:-c27}_f${F}_${MODE}.stderr) | grep '^RESULT' | sed 's/^RESULT //')
  [ -z "$R" ] && R="{\"ok\":false,\"errors\":[{\"what\":\"child\",\"kind\":\"Killed\",\"text\":\"no RESULT (timeout or memory)\"}]}"
  echo "{\"job\":\"$J\",\"exit_wall\":$(( $(date +%s) - t0 )),\"result\":$R}" >> $OUT
  echo "$(date +%T) $J $(echo $R | grep -o '"ok":[a-z]*' | head -1) $(( $(date +%s) - t0 ))s"
}
export -f one; export F D OUT CACHE RF MODE TAGP COSTS
cat $JF | xargs -P ${PAR:-1} -I{} bash -c 'one "{}"'
