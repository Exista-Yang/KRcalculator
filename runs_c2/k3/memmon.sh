#!/bin/bash
# sample RSS (MB) of sbox processes every 20 s -> memlog.tsv  (time, pid, word-len, rssMB)
while true; do ps -C sbox -o pid=,rss=,args= | awk -v t=$(date +%s) '{w=$0; sub(/.*--word "/,"",w); sub(/".*/,"",w); n=split(w,a," "); print t"\t"$1"\t"n"\t"int($2/1024)"\t"$0}' | cut -c1-120 >> runs_c2/k3/memlog.tsv; sleep 20; done
