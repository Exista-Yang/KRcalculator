#!/bin/bash
# pack the merged (2,2,2,2) layer A (M4c-lf stamp) into data-gl6 c2/l2222_n4_c25/fF (<= 45 MB chunks + MANIFEST)
# usage: persist.sh FIELD COMMIT SETS_TEXT
set -euo pipefail
cd /home/user/KRcalculator
S=/tmp/claude-0/-home-user-KRcalculator/4d3e75e2-4a4f-5c6d-a935-4d9b769b978c/scratchpad
F=$1; COMMIT=$2; SETS=$3
OUT=$S/dgl6/c2/l2222_n4_c25/f$F
mkdir -p $OUT; rm -f $OUT/*.part-* $OUT/MANIFEST.txt
SNAP=$(mktemp -d -p $S); mkdir -p $SNAP/layer_a $SNAP/tmpl
cp -p $S/lib2222m/f$F/layer_a/*.txt $SNAP/layer_a/; cp -p $S/lib2222m/f$F/tmpl/*.txt $SNAP/tmpl/
STAMPS=$(head -1 "$(ls -d $SNAP/layer_a/* | head -1)" | sed 's/ key=.*//')
NA=$(ls $SNAP/layer_a | wc -l); NT=$(ls $SNAP/tmpl | wc -l)
tar -C $SNAP -czf - layer_a tmpl | split -b 45M -d -a 3 - $OUT/layerA_tmpl.tar.gz.part-
{
  echo "lambda: 2,2,2,2"
  echo "field: F$F"
  echo "label: l2222_n4_c25/f$F"
  echo "code commit: $COMMIT (branch colored-2c4-fix)"
  echo "layer-A stamps (first file): $STAMPS"
  echo "provider: m4c-lf (production chain + Stage 0 + Leaves closure (leaves-first) + Leaves pairs); all files carry the M4c-lf stamp (C2-2 files carried over with b4wstamp)"
  echo "files: layer_a $NA, tmpl $NT"
  echo "restore: cat layerA_tmpl.tar.gz.part-* | tar -C CACHE_DIR -xzf -   (then KRHOM_COSTS_ANY_HOST=1 libbuild replay / sbox --replay --cache CACHE_DIR --costs c2/costs_2222_c25.txt --stamp-route --provider m4c-lf)"
  echo "sets: $SETS"
  echo "sha256:"
  (cd $OUT && sha256sum layerA_tmpl.tar.gz.part-*)
} > $OUT/MANIFEST.txt
rm -rf $SNAP
cat $OUT/MANIFEST.txt
