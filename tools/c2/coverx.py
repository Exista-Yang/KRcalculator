#!/usr/bin/env python3
"""X distribution of the library jobs each knot needs (libbuild cover --list output + plan JSON)
and an expected-failure projection from per-bucket success rates (default: (3,3,3) level 6, M5c F31)."""
import json, re, sys
plan = {j['job']: j for j in json.load(open(sys.argv[2]))['jobs']}
B = [0, 100, 500, 1000, 2000, 5000, 10000, 20000, 10**9]
RATE = {'pair': [1.0, 0.99, 0.92, 0.91, 0.73, 0.62, 0.75, 0.50], 'closure': [1.0, 0.87, 1.0, None, 0.0, 0.0, None, None]}
WALL = {'pair': [0.6, 7.0, 3.2, 4.8, 11.1, 7.7, 21.8, 16.1], 'closure': [0.3, 0.3, 0.8, None, None, None, None, None]}
name = None
for line in open(sys.argv[1]):
    m = re.match(r'\| (\S+) \| (\d+) \|', line)
    if m: name = m[1]; continue
    m = re.search(r'; jobs (.*)$', line)
    if not (m and name): continue
    jobs = [j.strip() for j in re.split(r', (?=pair |closure )', m[1])]
    out = []
    for kind in ('pair', 'closure'):
        js = [plan[j] for j in jobs if j.startswith(kind) and j in plan]
        if kind == 'pair': js = [j for j in js if j['want_targets'] > 0]
        h = [sum(1 for j in js if B[i] < j['x_flows'] <= B[i+1]) for i in range(len(B)-1)]
        exp_fail = sum(c * (1 - (RATE[kind][i] if RATE[kind][i] is not None else 0.5)) for i, c in enumerate(h))
        out.append('%s %d [%s] maxX %d expfail %.0f' % (kind, len(js), ' '.join(map(str, h)), max((j['x_flows'] for j in js), default=0), exp_fail))
    print(name, ' | '.join(out))
    name = None
