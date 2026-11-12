#!/bin/bash
cd /home/user/KRcalculator
while pgrep -f "sbox_c --lam 3,3 " > /dev/null; do sleep 20; done
R=runs_m6a1c/runw.sh
$R 6 31 101 u,r
$R 6 32b 102 u,r
$R 8 31 101 r
$R 8 32b 102 r
THREADS=2 MEMMB=7000 $R 9 31 101 r
$R 7 31 101 u,r
THREADS=2 MEMMB=7000 $R 9 32b 102 r
$R 7 32b 102 u,r
