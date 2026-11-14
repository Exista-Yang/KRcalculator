#!/bin/bash
# C2-7: pack the layer-A / template files written in this session into the restored (2,2,2,2) cache
# (C2-2 c2/l2222_n4 + C2-5 c2/l2222_n4_c25 + C2-4 c2/l2222_n4_c24) as a delta on top of them:
# data-gl6 c2/l2222_n4_c27/{f31,f32b}/ (≤ 45 MB chunks, MANIFEST, sha256). Files are selected by
# mtime (tar keeps archive mtimes on restore; the marker was written right after the restore).
set -e
cd /home/user/KRcalculator
DG=${DG:-/home/user/datagl6}
M=/home/user/work/costs_2222_c25.txt
for F in f31 f32b; do
  C=/home/user/work/cache/$F
  (cd $C && find layer_a tmpl -type f -newer $M | sort) > /tmp/new27_$F.lst
  O=$DG/c2/l2222_n4_c27/$F; mkdir -p $O; rm -f $O/*.part-*
  tar -C $C -czf - -T /tmp/new27_$F.lst | split -b 45m -d -a 3 - $O/layerA_tmpl_delta.tar.gz.part-
  {
    echo "lambda: 2,2,2,2"; echo "field: ${F^^}"; echo "label: l2222_n4_c27/$F (delta on top of l2222_n4 + l2222_n4_c25 + l2222_n4_c24, $F)"
    echo "code commit: $(git rev-parse --short HEAD) (branch colored-2c4-unblock)"
    echo "layer-A stamps: M4c-lf stamp (provider m4c-lf) plus the C2-7 route token c27=joint(sfx4,keep40,w12)+eqpairs (KRHOM_LEAVES_C27=1); sindec::route_header drops the token, so read with --stamp-route"
    echo "files: $(grep -c '^layer_a' /tmp/new27_$F.lst) layer_a, $(grep -c '^tmpl' /tmp/new27_$F.lst) tmpl (written in C2-7)"
    echo "restore: restore l2222_n4, l2222_n4_c25, l2222_n4_c24 ($F) in that order, then: cat layerA_tmpl_delta.tar.gz.part-* | tar -C CACHE_DIR -xzf -"
    echo "use: KRHOM_COSTS_ANY_HOST=1 libbuild one/replay | sbox --replay --cache CACHE_DIR --costs c2/costs_2222_c25.txt --provider m4c-lf --stamp-route"
    echo "jobs: the 6 class-A blockers (pairs 1235 1444 1470 1894 1769 1812) and later C2-7 jobs, see colored-2c4-unblock docs/colored_design/c2/C2_7_report.md"
    echo "sha256:"; (cd $O && sha256sum layerA_tmpl_delta.tar.gz.part-*)
  } > $O/MANIFEST.txt
  du -sh $O
done
