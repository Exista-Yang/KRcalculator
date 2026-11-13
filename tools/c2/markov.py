#!/usr/bin/env python3
"""Markov check: sbox sigma_1^m sigma_2^{+-1} on (2,2,2) vs xscan sigma^m on (2,2) (closures N' = 2,3,4,
unreduced and reduced). Reports equality, or equality up to a (t, q) shift, per closure."""
import json, sys, os
def R(f):
    d = json.load(open(f))
    return {(r['reduced'], r['N']): {(t, q): c for t, q, c in r['P']['terms'] if c} for r in d['results']}
def shift(a, b):
    if sum(a.values()) != sum(b.values()) or not a: return None
    ka, kb = min(a), min(b); s = (ka[0]-kb[0], ka[1]-kb[1])
    return s if {(t+s[0], q+s[1]): c for (t, q), c in b.items()} == a else None
bad = 0
for m in (2, 3, 4, 5):
    for sg in 'pn':
        for F in ('31', '32b'):
            fs = 'runs_c2/k3/markov/s%d%s_f%s.json' % (m, sg, F); fx = 'runs_c2/xscan/s%d_f%s_p4.json' % (m, F)
            if not (os.path.exists(fs) and os.path.exists(fx)): continue
            a, b = R(fs), R(fx)
            res = []
            for k in sorted(set(a) & set(b)):
                if a[k] == b[k]: res.append('%s%d:eq' % ('r' if k[0] else 'u', k[1]))
                else:
                    s = shift(a[k], b[k]); res.append('%s%d:%s' % ('r' if k[0] else 'u', k[1], 'shift%s' % (s,) if s else 'DIFF'))
                    bad += s is None
            print('sigma1^%d sigma2%s F%s: %s' % (m, '+' if sg == 'p' else '-', F, ' '.join(res)))
sys.exit(1 if bad else 0)
