#!/usr/bin/env python3
"""C2-7: append the 4-strand knot records (runs_c2/c27/k4/compare.json) to docs/colored_design/c2/gl4_2222.json
under key "c2_7" (replaced on rerun)."""
import json
p = 'docs/colored_design/c2/gl4_2222.json'
d = json.load(open(p))
recs = json.load(open('runs_c2/c27/k4/compare.json'))
d['c2_7'] = {
    'session': 'C2-7 (branch colored-2c4-unblock)',
    'library': 'data-gl6 c2/l2222_n4 + c2/l2222_n4_c25 + c2/l2222_n4_c24 + c2/l2222_n4_c27 (layer A only: sbox box/finish --replay --provider m4c-lf --stamp-route, costs_2222_c25.txt)',
    'box': 'prune 4, box N = 4 (1296 points; reduced-only runs use the 216 points with mu^0 = empty), closures N\' = 2, 3, 4, seed 201, sharded (4 shards + finish)',
    'checks': 'chi = decat_chi on every closure (fatal), D14 scaling, F31 = F32b, N\'=3 = uncolored gl_3 of the reversed braid (engine2, runs_c2/k4/e2ref), N\'=2 rank 1, generator determinant pass at every box point',
    'knots': recs,
}
json.dump(d, open(p, 'w'), indent=1)
print(len(recs), 'records')
