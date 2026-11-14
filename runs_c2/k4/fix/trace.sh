#!/bin/bash
# trace one failing (2,2,2,2) job (KRHOM_TRACE) with a given cost table / provider.
# usage: trace.sh FIELD COSTS OUTDIR [PROVIDER] -- JOB...   (jobs like "pair 958" or "closure p")
cd /home/user/KRcalculator
S=/tmp/claude-0/-home-user-KRcalculator/4d3e75e2-4a4f-5c6d-a935-4d9b769b978c/scratchpad
F=$1; C=$2; O=$3; PV=${4:+--provider $4}; shift 5
mkdir -p $O
one() {
  J="$1"; F=$2; C=$3; O=$4; PV=$5; S=$6
  N=$(echo "$J" | tr ' ,/' '___')
  (ulimit -v 6500000; KRHOM_TRACE=1 KRHOM_COSTS_ANY_HOST=1 KRHOM_THREADS=1 timeout ${TMO:-1800} $S/bin/libbuild one --lam 2,2,2,2 --n 4 --field $F --cache $S/lib2222/f$F --costs $C --stamp-route $PV -- $J > $O/$N.log 2>&1)
  echo "$J -> $(grep '^RESULT' $O/$N.log | grep -o '"ok":[a-z]*' | head -1) $(grep -c 'refused' $O/$N.log)"
}
export -f one
printf "%s\n" "$@" | xargs -P ${JOBS:-4} -I{} bash -c "one '{}' $F $C $O '$PV' $S"
