import json, os
from collections import defaultdict, Counter
HERE=os.path.dirname(os.path.abspath(__file__))
B2='/Users/existair/Documents/KRcalculator/.claude/worktrees/revrex/docs/colored_design/m4c/b2/'
recs=[json.loads(l) for f in ['gp0_n5b.jsonl','gp0_n6.jsonl'] for l in open(B2+f)]
steps=json.load(open(os.path.join(HERE,'steps.json')))
B=100000
def grp(o):
    if o['rexstep']: return 'E'
    if o['nk']==2: return 'A'
    j=o['joint'][0] if o['joint'] else None; ok=o['joint_ok'][0] if o['joint_ok'] else None
    if ok is not None and ok<=B: return 'B'
    if j is not None and j<=B: return 'C'
    if o['cf'] and o['cf'].startswith('Pu'): return 'D1'
    return 'D2'
stepgrp={}
for o in steps:
    for t in o['targets']:
        stepgrp[(t['n'],t['pair'],t['q'])]=(grp(o),o)
nonrex=[r for r in recs if not r['rex_class']]
inf=[r for r in nonrex if not r['feasible']]
print('non-rex targets', len(nonrex), 'not feasible', len(inf), 'distinct (n,pair,q) not feasible', len(set((r['n'],r['pair'],r['q']) for r in inf)))
G=defaultdict(lambda: dict(tg=0,pq=set(),pairs=set(),steps=set(),cen=[],j=[],loc=[],cf=set(),second=0))
for r in inf:
    k=(r['n'],r['pair'],r['q'])
    g,o=stepgrp.get(k,('?',None))
    e=G[g]; e['tg']+=1; e['pq'].add(k); e['pairs'].add((r['n'],r['pair'])); 
    if o:
        e['steps'].add((o['tp'],o['t'],o['rho'])); e['cen'].append(o['census'])
        jj=(o['joint_ok'][0] if o['joint_ok'] and g=='B' else (o['joint'][0] if o['joint'] else None)); e['j'].append(jj); e['loc'].append(o['loc']); e['cf'].add(o['cf'])
for g,e in sorted(G.items()):
    print(g, 'targets',e['tg'],'(n,pair,q)',len(e['pq']),'steps',len(e['steps']),'census',min(e['cen']),max(e['cen']),'3up',min(x for x in e['j'] if x is not None),max(x for x in e['j'] if x is not None),'loc4',min(e['loc']),max(e['loc']), sorted(e['pairs']))
# per pair verdict (worst-step level)
P=defaultdict(list)
for r in inf: P[(r['n'],r['pair'])].append(stepgrp[(r['n'],r['pair'],r['q'])][0])
for k,v in sorted(P.items()): print(k, Counter(v))
