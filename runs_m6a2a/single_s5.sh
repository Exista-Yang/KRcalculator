#!/bin/bash
# uninterrupted single-shard run of the same job (digest comparison with the sharded+killed run)
cd /home/user/KRcalculator; S=${S:?}; B=$S/sbox_m6a2a; ulimit -v 14000000
J="--lam 3,3 --word '1 1 1 1 1' --n 6 --ns 3,4,5,6 --field 31 --seed 101 --modes u,r --mem-mb 3500 --cache $S/cache33_f31"
O=runs_m6a2a/s5
eval KRHOM_THREADS=4 timeout 7200 $B box --dir $S/s5_single --ckpt-secs 1000000 $J > $O/single.out 2> $O/single.err
eval KRHOM_THREADS=4 timeout 7200 $B finish --dirs $S/s5_single --json $O/finish_single.json $J > $O/finish_single.out 2> $O/finish_single.err; echo "exit $?" >> $O/finish_single.err
echo DONE
