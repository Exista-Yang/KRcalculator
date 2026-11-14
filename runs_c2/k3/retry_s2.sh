#!/bin/bash
# retry knot-primes refused by the memory-cap ESTIMATE (cap raised to 16000 MB; ulimit -v 13 GB stays) when a core is free
cd /home/user/KRcalculator
export S=/tmp/claude-0/-home-user-KRcalculator/f29f8b74-5f26-5b1e-b5d7-0d65e52edefb/scratchpad MEMCAP=16000 MEMMB=16000 THREADS=1 TMO=14400
while IFS=';' read NAME F W; do
  while [ $(pgrep -fc "^/tmp.*/bin/sbox --lam") -ge 4 ]; do sleep 20; done
  runs_c2/k3/one_s2.sh "$NAME" "$F" "$W" >> runs_c2/k3/queue_s2_retry.log 2>&1 &
  sleep 60
done < ${1:-runs_c2/k3/retry_s2.txt}
wait
