#!/usr/bin/env python3
"""Per X bucket: final outcome and wall of the jobs of libbuild run JSONL files (last attempt per job)."""
import json, sys
B = [0, 100, 500, 1000, 2000, 5000, 10000, 20000, 10**9]
for f in sys.argv[1:]:
    last = {}
    for line in open(f):
        try: r = json.loads(line)
        except Exception: continue
        if 'job' in r and 'outcome' in r: last[r['job']] = r  # attempts are in order
    print(f, len(last), 'jobs')
    for kind in ('pair', 'closure'):
        for i in range(len(B)-1):
            js = [r for k, r in last.items() if k.startswith(kind) and B[i] < r['x_flows'] <= B[i+1]]
            if not js: continue
            ok = [r for r in js if r['outcome'] == 'ok']
            w = sorted(r['wall'] for r in ok)
            print('  %-7s X(%d,%s]: %d jobs, ok %d (%.0f%%), ok wall mean %.1f s, max %.1f s' % (kind, B[i], B[i+1], len(js), len(ok), 100*len(ok)/len(js), sum(w)/max(1,len(w)), max(w, default=0)))
