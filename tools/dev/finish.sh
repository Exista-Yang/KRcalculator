#!/bin/bash
# Finish of a sharded (knot, field, modes) box run of tools/dev/box.sh (merge, completeness,
# close, chi, scaling, JSON), e.g. after the shards were run separately or the finish failed.
# usage: tools/dev/finish.sh NAME FIELD "WORD" MODES K
# Same environment as box.sh (CACHE required; COSTS, OUTD, SEED, MEMMB, MEMCAP, VLIM,
# THREADS_FIN, KRHOM_REPO optional).
set -u
NAME=$1; F=$2; W=$3; M=$4; K=$5
REPO=${KRHOM_REPO:-$(cd "$(dirname "$0")/../.." && pwd)}
SBOX=$REPO/target/release/sbox
[ -x "$SBOX" ] || { echo "finish.sh: $SBOX not found (cargo build --release --features dev)" >&2; exit 2; }
CACHE=${CACHE:-${KRHOM_LIB_DIR:+$KRHOM_LIB_DIR/l2222_n4/f$F}}
[ -n "${CACHE:-}" ] || { echo "finish.sh: set CACHE (library directory of this prime)" >&2; exit 2; }
COSTS=${COSTS:-$CACHE/costs.txt}
D=${OUTD:-out}; TAG=${NAME}_${M/,/}_f$F
L=$(echo $W | wc -w)
ulimit -v ${VLIM:-14000000}
dirs=$(for k in $(seq 0 $((K-1))); do printf "%s," "$D/${TAG}_s$k"; done); dirs=${dirs%,}
t0=$(date +%s)
KRHOM_MEM_CAP_MB=${MEMCAP:-200000} KRHOM_COSTS_ANY_HOST=1 KRHOM_THREADS=${THREADS_FIN:-4} timeout 7200 "$SBOX" finish --lam 2,2,2,2 --word "$W" --n 4 --ns 2,3,4 --field "$F" --seed "${SEED:-201}" --prune 4 --mem-mb "${MEMMB:-3000}" --replay --cache "$CACHE" --costs "$COSTS" --provider m4c-lf --stamp-route --stages "$L" --modes "$M" --dirs "$dirs" --json "$D/$TAG.json" > "$D/$TAG.out" 2> "$D/$TAG.err"
echo "exit $?" >> "$D/$TAG.err"
echo "$(date +%T) finish $TAG $(tail -1 "$D/$TAG.err") $(( $(date +%s) - t0 ))s"
