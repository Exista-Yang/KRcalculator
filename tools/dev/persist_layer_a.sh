#!/bin/bash
# persist_layer_a.sh CACHE_DIR FIELD LABEL DEST_DIR COMMIT
# Packs CACHE_DIR/layer_a and CACHE_DIR/tmpl into tar.gz chunks of <= 45 MB under
# DEST_DIR/LABEL/ with a MANIFEST (lambda, sets, field, stamps, commit, file counts, sha256).
set -euo pipefail
CACHE=$1; FIELD=$2; LABEL=$3; DEST=$4; COMMIT=$5
OUT=$DEST/$LABEL
mkdir -p "$OUT"
rm -f "$OUT"/*.tar.gz.part-* "$OUT"/MANIFEST.txt
# snapshot first (a running libbuild may still be writing the cache)
SNAP=$(mktemp -d)
mkdir -p "$SNAP/layer_a" "$SNAP/tmpl"
cp -p "$CACHE"/layer_a/a-*.txt "$SNAP/layer_a/"
cp -p "$CACHE"/tmpl/*.txt "$SNAP/tmpl/"
CACHE=$SNAP
STAMPS=$(head -1 "$(ls -d $CACHE/layer_a/* | head -1)" | sed 's/ key=.*//')
NA=$(ls $CACHE/layer_a | wc -l); NT=$(ls $CACHE/tmpl | wc -l)
tar -C "$CACHE" -czf - layer_a tmpl | split -b 45M -d -a 3 - "$OUT/layerA_tmpl.tar.gz.part-"
{
  echo "lambda: 3,3,3"
  echo "field: F$FIELD"
  echo "label: $LABEL"
  echo "code commit: $COMMIT (branch colored-m5)"
  echo "layer-A stamps: $STAMPS"
  echo "files: layer_a $NA, tmpl $NT"
  echo "restore: cat layerA_tmpl.tar.gz.part-* | tar -C CACHE_DIR -xzf -   (then libbuild run/replay --cache CACHE_DIR --costs tests/data/homprov/costs_333.txt)"
  echo "sha256:"
  (cd "$OUT" && sha256sum layerA_tmpl.tar.gz.part-*)
} > "$OUT/MANIFEST.txt"
rm -rf "$SNAP"
cat "$OUT/MANIFEST.txt"
