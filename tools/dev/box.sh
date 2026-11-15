#!/bin/bash
# Sharded 4-strand box driver (C2-7): one (knot, field) of (2,2,2,2) gl_4 with `sbox`: prune 4,
# box N = 4, closures N' = 2,3,4, modes u,r, final stage; library from layer A only (--replay,
# provider m4c-lf refusing every call, route stamp). K shards in parallel (resumable: rerun the
# same command), then finish (chi = decat_chi and the D14 scaling check are mandatory inside
# `sbox finish`). This is how 6_1, 7_6 and 7_7 in results/colored_gl4/gl4_2222_c27_compare.json
# were computed (logs on branch archive-pre-restructure, runs_c2/c27/k4/).
#
# usage: tools/dev/box.sh NAME FIELD "WORD" [K]          (FIELD = 31 | 32b, K shards, default 4)
#
# Environment:
#   CACHE   library cache directory of this prime (the installed l2222_n4 library, e.g.
#           ~/.cache/krhom/lib/l2222_n4/f31; default: $KRHOM_LIB_DIR/l2222_n4/f$FIELD if
#           KRHOM_LIB_DIR is set). Work on a copy if the installed library must stay
#           untouched (sbox may add cache files).
#   COSTS   cost table [$CACHE/costs.txt: the C2-5 calibration costs_2222_c25.txt that
#           tools/fetch_library.sh installs with the library; accepted on any host via
#           KRHOM_COSTS_ANY_HOST=1]
# Optional environment (defaults in brackets):
#   OUTD [out]  MODES [u,r]  SEED [201]  MEMMB [3000]  MEMCAP [200000]  VLIM [14000000 KB]
#   THREADS [1] per shard  THREADS_FIN [4]  TMO [36000 s] per shard  CKPT [600 s]
#   KRHOM_REPO [the repository containing this script]  (target/release/sbox is used;
#   build it with `cargo build --release --features dev`)
set -u
NAME=$1; F=$2; W=$3; K=${4:-4}
REPO=${KRHOM_REPO:-$(cd "$(dirname "$0")/../.." && pwd)}
SBOX=$REPO/target/release/sbox
[ -x "$SBOX" ] || { echo "box.sh: $SBOX not found (cargo build --release --features dev)" >&2; exit 2; }
CACHE=${CACHE:-${KRHOM_LIB_DIR:+$KRHOM_LIB_DIR/l2222_n4/f$F}}
[ -n "${CACHE:-}" ] || { echo "box.sh: set CACHE (library directory of this prime)" >&2; exit 2; }
COSTS=${COSTS:-$CACHE/costs.txt}
D=${OUTD:-out}; mkdir -p "$D"
M=${MODES:-u,r}
TAG=${NAME}_${M/,/}_f$F
[ -f "$D/$TAG.json" ] && exit 0
L=$(echo $W | wc -w)
ulimit -v ${VLIM:-14000000}
C=(--lam 2,2,2,2 --word "$W" --n 4 --ns 2,3,4 --field "$F" --seed "${SEED:-201}" --prune 4 --mem-mb "${MEMMB:-3000}" --replay --cache "$CACHE" --costs "$COSTS" --provider m4c-lf --stamp-route --stages "$L" --modes "$M")
export KRHOM_MEM_CAP_MB=${MEMCAP:-200000} KRHOM_COSTS_ANY_HOST=1
t0=$(date +%s)
dirs=""
for k in $(seq 0 $((K-1))); do
  SD=$D/${TAG}_s$k; mkdir -p "$SD"; dirs="$dirs,$SD"
  KRHOM_THREADS=${THREADS:-1} timeout ${TMO:-36000} "$SBOX" box "${C[@]}" --dir "$SD" --shard $k/$K --ckpt-secs ${CKPT:-600} > "$SD/box.out" 2> "$SD/box.err" &
done
wait
dirs=${dirs#,}
KRHOM_THREADS=${THREADS_FIN:-4} timeout 7200 "$SBOX" finish "${C[@]}" --dirs "$dirs" --json "$D/$TAG.json" > "$D/$TAG.out" 2> "$D/$TAG.err"
echo "exit $?" >> "$D/$TAG.err"
echo "$(date +%T) $TAG $(tail -1 "$D/$TAG.err") $(( $(date +%s) - t0 ))s"
