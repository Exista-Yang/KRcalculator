#!/bin/bash
# second driver: wait until at most 1 old job is running, then P=3 (memory spikes of 4-5 GB seen in the box phase)
cd /home/user/KRcalculator
while [ $(pgrep -fc "one_s1.sh [^ ]* [^ ]* ") -gt 2 ]; do sleep 30; done
P=3 bash runs_c2/k3/queue_s1.sh
