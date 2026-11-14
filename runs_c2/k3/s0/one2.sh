#!/bin/bash
# one (knot, field) with admission control: start only if < MAXJ sbox running and MemAvailable > MINAVAIL_MB;
# skips if the same (word, field) is already running. Starts are serialized with flock.
cd /home/user/KRcalculator
NAME=$1; F=$2; W=$3
D=${K3DIR:-runs_c2/k3/s0/out}; TAG=${NAME}_f${F}
[ -f $D/$TAG.json ] && exit 0
[ -f $D/$TAG.err ] && grep -q "^exit" $D/$TAG.err && exit 0
pgrep -f -- "--word $W --n 4 --ns 2,3,4 --field $F " >/dev/null && exit 0
exec 9>/tmp/claude-0/s0_admit.lock
flock 9
while [ $(pgrep -c -x sbox) -ge ${MAXJ:-3} ] || [ $(awk '/MemAvailable/{print int($2/1024)}' /proc/meminfo) -lt ${MINAVAIL_MB:-6000} ]; do sleep 20; done
# launch in background, then hold the lock 120 s so memory growth of the new job becomes visible
runs_c2/k3/s0/one.sh "$NAME" "$F" "$W" 9>&- &
sleep 120
flock -u 9
wait
