#!/bin/bash
cd /home/user/KRcalculator
FA=$(cat runs_c2/k4m/force/force_A.txt)
FB="0,0,2/0,2,0/2,0,0=F0.2,F1.2,E1.2,E0.2;1,0,1/0,1,1/1,1,0=F0.1,F1.2,E1.2,E0.1"
echo "$FB" > runs_c2/k4m/force/force_B.txt
for F in 31 32b; do
  runs_c2/k4m/force/one.sh w41A $F "1 -2 1 -2" "$FA"
  runs_c2/k4m/force/one.sh w41B $F "1 -2 1 -2" "$FB"
  runs_c2/k4m/force/one.sh w52A $F "1 1 1 2 -1 2" "$FA"
  runs_c2/k4m/force/one.sh w52B $F "1 1 1 2 -1 2" "$FB"
  runs_c2/k4m/force/one.sh w31B $F "1 1 1 2" "$FB"
done
