#!/bin/bash
# memory guard + retry runner for slice 2: kill the youngest sbox when total RSS > LIM_KB; retry killed/refused knot-primes later (one at a time, when RSS < LOW_KB and < 4 sbox)
cd /home/user/KRcalculator
LIM=${LIM_KB:-13300000}; LOW=${LOW_KB:-8500000}
export S=/tmp/claude-0/-home-user-KRcalculator/f29f8b74-5f26-5b1e-b5d7-0d65e52edefb/scratchpad MEMCAP=${CAP:-16000} MEMMB=${CAP:-16000} THREADS=1 TMO=14400
mkdir -p runs_c2/k3/out_s2/killed
declare -A started
while true; do
  tot=$(ps -eo rss,args | grep "[b]in/sbox --lam" | grep -v timeout | awk '{s+=$1} END{print s+0}')
  n=$(ps -eo args | grep "[b]in/sbox --lam" | grep -vc timeout)
  if [ "$tot" -gt "$LIM" ]; then
    pid=$(ps -eo pid,etimes,args | grep "[b]in/sbox --lam" | grep -v timeout | sort -k2 -n | head -1 | awk '{print $1}')
    tag=$(tr '\0' ' ' < /proc/$pid/cmdline | grep -o 'out_s2/[^ ]*\.json' | sed 's#out_s2/##; s#\.json##')
    echo "$(date +%T) guard: total ${tot} kB > limit, killing $tag (pid $pid)" >> runs_c2/k3/guard_s2.log
    kill $pid; sleep 6
    name=${tag%_f*}; f=${tag##*_f}
    w=$(grep "^$name;" runs_c2/k3/slices/slice2.txt | cut -d';' -f2)
    mv runs_c2/k3/out_s2/$tag.err runs_c2/k3/out_s2/killed/$tag.guard.$(date +%s).err 2>/dev/null; rm -f runs_c2/k3/out_s2/$tag.out
    echo "$name;$f;$w" >> runs_c2/k3/retry_s2.txt
    unset started[$tag]
  elif [ "$n" -lt ${MAXN:-3} ] && [ "$tot" -lt "$LOW" ]; then
    while IFS=';' read NAME F W; do
      tag=${NAME}_f$F
      [ -n "${started[$tag]}" ] && continue
      [ -f runs_c2/k3/out_s2/$tag.json ] && continue
      [ -f runs_c2/k3/out_s2/$tag.err ] && continue
      started[$tag]=1
      echo "$(date +%T) retry $tag" >> runs_c2/k3/guard_s2.log
      nohup runs_c2/k3/one_s2.sh "$NAME" "$F" "$W" >> runs_c2/k3/queue_s2_retry.log 2>&1 &
      sleep 90; break
    done < runs_c2/k3/retry_s2.txt
  fi
  sleep 10
done
