# c2: (2,2,2) decomposition library for 2-colored gl_4 (session C2-1, branch colored-2c4)

- l222_n4/{f31,f32b}: layer A + templates, prune 4 (level-n4 105/105 + corner-markov-ext-n4 21/21, replay 100 %).
- l222_n6_unpruned/{f31,f32b}: unpruned (level-n6 146/146 + corner 21/21, replay 100 %), for pruned-vs-unpruned checks.
- runs/: libbuild run/replay JSON lines and logs; costs_cloud.txt = the cost table used (copy of tests/data/homprov/costs_cloud.txt).
Restore: see each MANIFEST.txt (KRHOM_COSTS_ANY_HOST=1, --costs tests/data/homprov/costs_cloud.txt).

## C2-2: (2,2,2,2) library (session C2-2, branch colored-2c4)

- l2222_n4/{f31,f32b}: layer A + templates of the (2,2,2,2) prune-4 library (level-n4 + corner-markov-ext-n4).
  Level 1087/1157 per prime (70 failures, identical for both primes: 43 NeedsIdempotent = M5i types, 27 Refused /
  Exhausted), corner 144/144. See colored-2c4 docs/colored_design/c2/C2_2_report.md.
- costs_2222.txt: the calibrated cost table (costs_cloud.txt + 88 (2,2,2,2) closure calls; runs_c2/k4/calib).
  Layer-A files carry either table's stamp: read with `--costs costs_2222.txt --stamp-route`.
- runs2222/: libbuild run / fix / replay JSON lines and logs.
