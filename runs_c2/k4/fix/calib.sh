#!/bin/bash
# C2-5 calibration round: `libbuild calib --windows` of the traced refused windows into
# runs_c2/k4/fix/costs_2222_c25.txt (base: the previous table). usage: calib.sh WINDOWS BASE [LIMIT] [RSS]
cd /home/user/KRcalculator
S=/tmp/claude-0/-home-user-KRcalculator/4d3e75e2-4a4f-5c6d-a935-4d9b769b978c/scratchpad
ulimit -v 14000000
KRHOM_COSTS_ANY_HOST=1 timeout ${TMO:-20000} $S/bin/libbuild calib --lam 2,2,2,2 --n 4 --field 31 --cache $S/calib_cache --base $2 --out runs_c2/k4/fix/costs_2222_c25.txt --windows $1 --limit ${3:-125} --budget ${BUDGET:-7200} --rss ${4:-3000} --jobs ${JOBS:-4} ${RETRY:+--retry-censored}
