import json, os, sys
HERE=os.path.dirname(os.path.abspath(__file__))
steps=json.load(open(os.path.join(HERE,'steps.json')))
B=100000
def cls(o):
    if o['rexstep']: return 'Id'
    j=o['joint'][0] if o['joint'] else None
    ok=o['joint_ok'][0] if o['joint_ok'] else None
    if ok is not None and ok<=B: return 'R'
    if j is not None and j<=B: return 'R?'
    if (ok or j) is not None and min(x for x in [ok,j] if x is not None)<=4*B: return 'M'
    return 'X'
T={}
for o in steps:
    c=cls(o)
    for t in o['targets']:
        key=(t['n'],t['pair'],t['q'])
        e=T.setdefault(key,dict(cls=c,final=o['final'],rep=t['rep'],nonfinal=t['nonfinal'],ss=[],step=f"{o['tp']}->{o['t']} {o['rho']}",ok=o['joint_ok'][0] if o['joint_ok'] else None,j=o['joint'][0] if o['joint'] else None,loc4=o['loc'],alt=t['alt']))
        e['ss'].append((t['s'],t['mult']))
for k,e in sorted(T.items()):
    second = (e['final'] and e['nonfinal']>B) or ((not e['final']) and e['alt'] and int(e['alt'].split('(')[-1].split()[0])>B)
    e['second']=second
    print(k, e['cls'], 'final' if e['final'] else 'nonfinal', e['rep'], 'nonfinal',e['nonfinal'], 'joint',e['j'],'ok',e['ok'],'loc4',e['loc4'],'SECOND' if second else '', e['ss'], e['alt'][:60] if e['alt'] else '')
from collections import Counter
print(Counter(e['cls'] for e in T.values()), sum(len(e['ss']) for e in T.values()))
print('with second known >budget:', sum(1 for e in T.values() if e['second']))
json.dump({f"{k[0]}|{k[1]}|{k[2]}":v for k,v in T.items()}, open(os.path.join(HERE,'targets.json'),'w'), indent=1)
