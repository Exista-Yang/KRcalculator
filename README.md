# data-gl6: data for the 3-colored gl_6 programme (SPEC_3C_GL6 rev1 D12)

Compressed chunks (<= 45 MB) with a MANIFEST per directory. No code here; the code is on colored-m5.

- m5b/f31, m5b/f32b: sindec layer A + template library of the (3,3,3) library (M5b), with the libbuild
  run logs (runs/*.jsonl, *.targets).
- m5c/f31, m5c/f32b: layer A + templates after M5c (level 5: 172/177 jobs both primes; level 6: 330/377 F31,
  329/377 F32b), run and replay logs (runs/; replay_n6_*_first.jsonl = replay before the resume pass), and
  m5c/costs_333_m5c.txt (= tests/data/homprov/costs_333_m5c.txt). Layer A mixes M5b- and M5c-table stamps: restore and
  load with `--stamp-route` (and KRHOM_COSTS_ANY_HOST=1 on another host class). Report: docs/colored_design/gl6/
  M5c_report.md on colored-m5.
