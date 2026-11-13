#!/bin/bash
# Markov / stabilized 4-strand words, one field. usage: all.sh FIELD
cd /home/user/KRcalculator
while IFS=';' read NAME LAM W NN; do runs_c2/k4/markov/one.sh $NAME $1 "$W"; done < runs_c2/k4/markov_words.txt
