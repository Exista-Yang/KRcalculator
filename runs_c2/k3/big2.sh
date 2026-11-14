#!/bin/bash
cd /home/user/KRcalculator
export ULIM=12000000 MEMCAP=120000 MEMMB=12000 THREADS=2 TMO=28000
printf '10_100\t31\t-1 -1 -1 2 -1 -1 2 -1 -1 2\n10_100\t32b\t-1 -1 -1 2 -1 -1 2 -1 -1 2\n' | xargs -d '\n' -P 2 -n 1 bash -c 'IFS=$'"'"'\t'"'"' read -r a b c <<< "$0"; runs_c2/k3/one2.sh "$a" "$b" "$c"'
