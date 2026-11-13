#!/usr/bin/env python3
"""X-flow distribution of libbuild plan JSONs (pair jobs with provider targets, closure jobs)."""
import json, sys
B = [0, 100, 500, 1000, 2000, 5000, 10000, 20000, 10**9]
for f in sys.argv[1:]:
    d = json.load(open(f))
    pairs = [j for j in d['jobs'] if j['job'].startswith('pair') and j['want_targets'] > 0]
    clos = [j for j in d['jobs'] if j['job'].startswith('closure')]
    print(f, d['plan'].split(':')[0], '| provider pairs', len(pairs), 'targets', sum(j['want_targets'] for j in pairs), '| closure', len(clos))
    for name, js in (('pairs', pairs), ('closure', clos)):
        h = [sum(1 for j in js if B[i] < j['x_flows'] <= B[i+1]) for i in range(len(B)-1)]
        print('  %-8s' % name, ' '.join('(%d,%s]:%d' % (B[i], B[i+1] if B[i+1] < 10**9 else 'inf', h[i]) for i in range(len(h))), ' max', max((j['x_flows'] for j in js), default=0))
