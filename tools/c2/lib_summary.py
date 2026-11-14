#!/usr/bin/env python3
"""C2-2: (2,2,2,2) library outcome per set and prime after the run (+ the fix pass with the calibrated
table): certified, failures by cause and job kind, max X, walls. usage: lib_summary.py FIELD"""
import json, sys, collections
F = sys.argv[1]
D = 'runs_c2/k4/lib'
def last(f):
    st, x, summ = {}, {}, None
    for l in open(f):
        r = json.loads(l)
        if 'attempt' in r:
            j = r['job']
            if st.get(j) != 'ok': st[j] = r['outcome']
            x[j] = r.get('x_flows', 0)
        elif 'summary' in r: summ = r['summary']
    return st, x, summ
fix = {}
try:
    for l in open('%s/l2222_n4_f%s_fix.jsonl' % (D, F)):
        r = json.loads(l); fix[r['job']] = bool(r['result'] and r['result'].get('ok'))
except FileNotFoundError: pass
for s in ('', 'corner'):
    f = '%s/l2222_n4%s_f%s_run.jsonl' % (D, s, F)
    st, x, summ = last(f)
    n = len(st); ok1 = sum(v == 'ok' for v in st.values())
    after = {j: ('ok' if v == 'ok' or fix.get(j) else v) for j, v in st.items()}
    ok2 = sum(v == 'ok' for v in after.values())
    c = collections.Counter((j.split()[0], v) for j, v in after.items() if v != 'ok')
    print('F%s %s: jobs %d, run %d (%.1f%%), after fix %d (%.1f%%); max X certified %d; failures %s' % (F, s or 'level', n, ok1, 100 * ok1 / n, ok2, 100 * ok2 / n,
          max(x[j] for j, v in after.items() if v == 'ok'), dict(c)))
    if summ: print('   summary:', {k: summ[k] for k in summ if k in ('wall', 'total', 'n_ok', 'jobs')})
