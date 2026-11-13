#!/usr/bin/env python3
"""Fit per-run scan seconds of the (2,2,2) queue against the truncated envelope (peak Hom dim at N = 4)
as a power law, and project the remaining 3-strand runs and 4-strand per-point costs."""
import json, glob, math, re, sys
def env(path):
    e = {}
    for l in open(path):
        f = l.rstrip().split('\t')
        if len(f) > 3 and f[2] == 'N=4': e[f[0]] = {k: int(v) for k, v in (x.split('=') for x in f[3:])}
    return e
e3 = env('runs_c2/census/env3_222.txt')
pts = []
for f in glob.glob('runs_c2/k3/out/*.json'):
    d = json.load(open(f)); k = f.split('/')[-1].rsplit('_f', 1)[0]
    pts.append((k, e3[k]['peakldim'], d['point_secs_mean']))
xs = [math.log(p[1]) for p in pts]; ys = [math.log(p[2]) for p in pts]
n = len(xs); mx = sum(xs)/n; my = sum(ys)/n
b = sum((x-mx)*(y-my) for x, y in zip(xs, ys))/sum((x-mx)**2 for x in xs); a = my - b*mx
pred = lambda pk: math.exp(a + b*math.log(pk))
print('fit: point_secs = exp(%.3f) * peakldim^%.3f over %d runs (1 thread, 4 jobs on 4 cores)' % (a, b, n))
res = sorted((math.exp(y)/pred(math.exp(x))) for x, y in zip(xs, ys)); print('residual ratio range %.2f .. %.2f' % (res[0], res[-1]))
done = {f.split('/')[-1][:-5] for f in glob.glob('runs_c2/k3/out/*.json')}
rem = [(k, 216*pred(v['peakldim'])) for k, v in e3.items() for F in ('31', '32b') if '%s_f%s' % (k, F) not in done]
print('remaining 3-strand runs %d, projected %.1f h on 4 cores' % (len(rem), sum(s for _, s in rem)/4/3600))
if len(sys.argv) > 1:
    e4 = env(sys.argv[1])
    for k, v in e4.items():
        s = pred(v['peakldim'])
        print('4-strand %-7s peakldim %9d  projected per point %8.1f s  unreduced 1296 pts %7.1f core-h  reduced-only 216 pts %6.1f core-h' % (k, v['peakldim'], s, 1296*s/3600, 216*s/3600))
