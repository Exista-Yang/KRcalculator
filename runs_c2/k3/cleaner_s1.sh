#!/bin/bash
cd /home/user/KRcalculator/runs_c2/k3/out_s1
while true; do for j in *.json; do t=${j%.json}; case $t in *.fin) continue;; esac; [ -d $t.d ] && rm -rf $t.d; done; sleep 120; done
