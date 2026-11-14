#!/bin/bash
# calibration of the refused (2,2,2,2) calls (closure calls of the failing closure types, Search windows
# of refused pairs) into a new cost table based on costs_cloud.txt. usage: calib.sh WINDOWS_FILE [LIMIT] [BUDGET]
cd /home/user/KRcalculator
S=/tmp/claude-0/-home-user-KRcalculator/7ca6fc3a-b081-5be7-b7f0-2501a115a78c/scratchpad
ulimit -v 14000000
KRHOM_COSTS_ANY_HOST=1 timeout 20000 $S/bin/libbuild calib --lam 2,2,2,2 --n 4 --field 31 --cache $S/calib_cache --base tests/data/homprov/costs_cloud.txt --out runs_c2/k4/calib/costs_2222.txt --windows $1 --limit ${2:-125} --budget ${3:-3600} --rss 3000 --jobs ${JOBS:-4}
