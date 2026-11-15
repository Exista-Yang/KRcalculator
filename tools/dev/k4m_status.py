#!/usr/bin/env python3
"""C2-4: job status lines "JOB ok|fail" from libbuild JSONL (attempt lines) and runjobs JSONL
({"job":..,"result":{"ok":..}}). A job is ok if any given file certifies it.
usage: k4m_status.py OUT FILE1 [FILE2 ...]"""
import json, sys
st = {}
for f in sys.argv[2:]:
    for l in open(f):
        try:
            d = json.loads(l)
        except Exception:
            continue
        if 'attempt' in d:
            j = d.get('job'); ok = d.get('outcome') == 'ok'
        elif 'job' in d and 'result' in d:
            j = d['job']; ok = bool(d['result'].get('ok'))
        else:
            continue
        if j is None: continue
        st[j] = st.get(j, False) or ok
with open(sys.argv[1], 'w') as o:
    for j in sorted(st):
        o.write('%s %s\n' % (j, 'ok' if st[j] else 'fail'))
print(sum(st.values()), 'ok of', len(st))
