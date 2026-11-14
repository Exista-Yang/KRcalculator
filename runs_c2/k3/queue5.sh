#!/bin/bash
# C2-6: sequential knot-by-knot queue. cost < 700000: both primes in parallel (1 thread each, ulimit 6.8 GB);
# cost >= 700000: one prime at a time (2 threads, ulimit 12 GB) because real RSS reaches 7+ GB.
cd /home/user/KRcalculator
export MEMCAP=120000 TMO=28000
while pgrep -x sbox >/dev/null; do sleep 20; done
( cat runs_c2/k3/order_le10.txt; cat runs_c2/k3/order.txt ) | awk -F';' '!s[$1]++' | grep -v -E '^(8_18|10_85|10_100|10_148|10_155);' | while IFS=';' read NAME W PK; do
  [ -f runs_c2/k3/out/${NAME}_f31.err ] && continue
  if [ "$PK" -lt 700000 ]; then
    for F in 31 32b; do ULIM=6800000 MEMMB=6000 THREADS=2 runs_c2/k3/one2.sh $NAME $F "$W" & done; wait
  else
    for F in 31 32b; do ULIM=12000000 MEMMB=12000 THREADS=3 runs_c2/k3/one2.sh $NAME $F "$W"; done
  fi
done
