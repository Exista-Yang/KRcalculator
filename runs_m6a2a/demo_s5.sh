#!/bin/bash
# M6a-2a demonstration: (3,3) sigma^5 at N = 6 (box N 6, N' = 3..6, u+r, F31 seed 101 = the M6a-1c
# G2B run) split into 2 shards; shard 0 is killed (SIGKILL) mid-run and resumed; finish merges.
cd /home/user/KRcalculator
S=${S:?scratch dir}
B=$S/sbox_m6a2a
cp target/release/sbox $B
ulimit -v 14000000
J="--lam 3,3 --word '1 1 1 1 1' --n 6 --ns 3,4,5,6 --field 31 --seed 101 --modes u,r --mem-mb 3500 --cache $S/cache33_f31"
O=runs_m6a2a/s5
mkdir -p $O
# library (built once into the cache; later invocations load layer A)
eval KRHOM_THREADS=4 timeout 7200 $B point --id 0 $J > $O/lib.out 2> $O/lib.err
# shard 1/2 in the background, shard 0/2 killed after 150 s, then resumed
( eval KRHOM_THREADS=2 timeout 7200 $B box --dir $O/shard1 --shard 1/2 --ckpt-secs 20 $J > $O/shard1.out 2> $O/shard1.err; echo "exit $?" >> $O/shard1.err ) &
eval KRHOM_THREADS=2 timeout -s KILL 150 $B box --dir $O/shard0 --shard 0/2 --ckpt-secs 20 $J > $O/shard0_killed.out 2> $O/shard0_killed.err; echo "exit $? (killed)" >> $O/shard0_killed.err
eval KRHOM_THREADS=2 timeout 7200 $B box --dir $O/shard0 --shard 0/2 --ckpt-secs 20 $J > $O/shard0_resumed.out 2> $O/shard0_resumed.err; echo "exit $?" >> $O/shard0_resumed.err
wait
eval KRHOM_THREADS=4 timeout 7200 $B finish --dirs $O/shard0,$O/shard1 --json $O/finish.json $J > $O/finish.out 2> $O/finish.err; echo "exit $?" >> $O/finish.err
echo DONE
