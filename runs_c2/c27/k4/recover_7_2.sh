#!/bin/bash
# 7_2 F31 reduced: shard 1 was OOM-killed (4 x 3.5 GB); wait for shards 0, 2, 3, rerun shard 1 alone
# (resumable box run), then finish; then the rest of queue2 with 3 shards per run (K3).
cd /home/user/KRcalculator
while pgrep -f "sbox box.*1 1 1 2 -1 2 3 -2 3" > /dev/null; do sleep 20; done
W="1 1 1 2 -1 2 3 -2 3"; D=runs_c2/c27/k4/out; TAG=7_2_r_f31
ulimit -v 14000000
export KRHOM_MEM_CAP_MB=200000 KRHOM_COSTS_ANY_HOST=1 KRHOM_THREADS=1
C="--lam 2,2,2,2 --word \"$W\" --n 4 --ns 2,3,4 --field 31 --seed 201 --prune 4 --mem-mb 3000 --replay --cache /home/user/work/cache/f31 --costs /home/user/work/costs_2222_c25.txt --provider m4c-lf --stamp-route --stages 9 --modes r"
eval timeout 36000 ./target/release/sbox box $C --dir $D/${TAG}_s1 --shard 1/4 --ckpt-secs 600 > $D/${TAG}_s1/box.out 2>> $D/${TAG}_s1/box.err
runs_c2/c27/k4/finish.sh 7_2 31 "$W" r 4
tail -n +2 runs_c2/c27/k4/queue2.txt | while IFS=';' read NAME F M WW; do
  [ -z "$NAME" ] && continue
  MODES=$M runs_c2/c27/k4/box.sh $NAME $F "$WW" 3
done
