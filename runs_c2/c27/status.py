#!/usr/bin/env python3
"""C2-7: union status file ("JOB ok") from replay jsonl files (C2-5 full replay, C2-4 M5i replay,
C2-7 replay): a job is ok if any replay certified it. usage: status.py OUT IN1 IN2 ..."""
import json, sys
ok = set(); seen = set()
for f in sys.argv[2:]:
    for line in open(f):
        try:
            d = json.loads(line)
        except Exception:
            continue
        j = d.get('job')
        if not j:
            continue
        seen.add(j)
        r = d.get('result') or {}
        if r.get('ok') or d.get('outcome') == 'ok':
            ok.add(j)
with open(sys.argv[1], 'w') as o:
    for j in sorted(seen):
        o.write(f"{j} {'ok' if j in ok else 'fail'}\n")
print(f"{len(ok)} ok / {len(seen)} jobs")
