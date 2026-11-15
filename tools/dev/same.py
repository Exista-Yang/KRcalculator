#!/usr/bin/env python3
"""same.py A.json B.json: are all closures (stage, reduced, N) of two sbox/xscan JSONs equal?"""
import json, sys
def R(f):
    d = json.load(open(f))
    return {(r.get('stage'), r['reduced'], r['N']): {(t, q): c for t, q, c in r['P']['terms'] if c} for r in d['results']}
a, b = R(sys.argv[1]), R(sys.argv[2])
ka = {(k[1], k[2]) for k in a}; kb = {(k[1], k[2]) for k in b}
common = sorted(ka & kb)
A = {(k[1], k[2]): v for k, v in a.items()}; B = {(k[1], k[2]): v for k, v in b.items()}
bad = [k for k in common if A[k] != B[k]]
print(sys.argv[1], sys.argv[2], 'common', len(common), 'equal' if not bad else 'DIFF %s' % bad)
sys.exit(1 if bad else 0)
