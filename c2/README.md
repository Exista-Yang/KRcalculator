# c2: (2,2,2) decomposition library for 2-colored gl_4 (session C2-1, branch colored-2c4)

- l222_n4/{f31,f32b}: layer A + templates, prune 4 (level-n4 105/105 + corner-markov-ext-n4 21/21, replay 100 %).
- l222_n6_unpruned/{f31,f32b}: unpruned (level-n6 146/146 + corner 21/21, replay 100 %), for pruned-vs-unpruned checks.
- runs/: libbuild run/replay JSON lines and logs; costs_cloud.txt = the cost table used (copy of tests/data/homprov/costs_cloud.txt).
Restore: see each MANIFEST.txt (KRHOM_COSTS_ANY_HOST=1, --costs tests/data/homprov/costs_cloud.txt).
