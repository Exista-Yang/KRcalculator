#!/bin/bash
# C2-4: pack the layer-A / template files added to the C2-5 (2,2,2,2) cache (data-gl6 c2/l2222_n4_c25) by the
# M5i jobs, as a delta on top of it: data-gl6 c2/l2222_n4_c24/{f31,f32b}/ (≤ 45 MB chunks, MANIFEST, sha256).
set -e
cd /home/user/KRcalculator
DG=/home/user/datagl6
for F in f31 f32b; do
  C=/home/user/lib2222c/$F
  git -C $DG show origin/data-gl6:c2/l2222_n4_c25/$F/layerA_tmpl.tar.gz.part-000 | tar -tzf - | sed 's|^\./||' | sort > /tmp/c25_$F.lst
  (cd $C && find layer_a tmpl -type f | sort) > /tmp/now_$F.lst
  comm -13 /tmp/c25_$F.lst /tmp/now_$F.lst > /tmp/new_$F.lst
  O=$DG/c2/l2222_n4_c24/$F; mkdir -p $O; rm -f $O/*.part-*
  tar -C $C -czf - -T /tmp/new_$F.lst | split -b 45m -d -a 3 - $O/layerA_tmpl_delta.tar.gz.part-
  {
    echo "lambda: 2,2,2,2"; echo "field: ${F^^}"; echo "label: l2222_n4_c24/$F (delta on top of l2222_n4_c25/$F)"
    echo "code commit: $(git rev-parse --short HEAD) (branch colored-2c4-m5i)"
    echo "layer-A stamps: as l2222_n4_c25 (provider m4c-lf, M4c-lf stamp); new keys kind=idem and eq=/ep=m5i-peel1 (M5i)"
    echo "files: $(grep -c '^layer_a' /tmp/new_$F.lst) layer_a, $(grep -c '^tmpl' /tmp/new_$F.lst) tmpl (new relative to c25)"
    echo "restore: restore l2222_n4_c25/$F first, then: cat layerA_tmpl_delta.tar.gz.part-* | tar -C CACHE_DIR -xzf -"
    echo "use: KRHOM_COSTS_ANY_HOST=1 libbuild one/replay | sbox --replay --cache CACHE_DIR --costs c2/costs_2222_c25.txt --provider m4c-lf --stamp-route"
    echo "jobs: the 43 C2-2 NeedsIdempotent jobs, see colored-2c4-m5i docs/colored_design/c2/C2_4_report.md"
    echo "sha256:"; (cd $O && sha256sum layerA_tmpl_delta.tar.gz.part-*)
  } > $O/MANIFEST.txt
  du -sh $O
done
