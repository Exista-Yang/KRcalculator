import json, os, sys, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kr
HERE=os.path.dirname(os.path.abspath(__file__)); LAM=[3,3,3]
steps=json.load(open(os.path.join(HERE,'steps.json')))
B2='/Users/existair/Documents/KRcalculator/.claude/worktrees/revrex/docs/colored_design/m4c/b2/'
recs=[json.loads(l) for f in ['gp0_n5b.jsonl','gp0_n6.jsonl'] for l in open(B2+f)]
lq={(r['n'],r['pair'],r['q']):r['lq'] for r in recs}
def u0(A,B):
    if A==B: return 0
    ws=kr.win_split(LAM,A,B)
    whole=max(kr.hom_u(LAM,A,B,0),kr.hom_u(LAM,B,A,0))
    if ws['at']>0 or ws['suffix']>0:
        return min(whole,max(kr.hom_u(ws['labels'],ws['a'],ws['b'],0),kr.hom_u(ws['labels'],ws['b'],ws['a'],0)))
    return whole
out=[]
for o in steps:
    pick=o['joint_ok'] if (o['joint_ok'] and o['joint_ok'][0]<=400000) else o['joint']
    if not pick: continue
    a=kr.parse_word(pick[1]); b=kr.parse_word(pick[2])
    rins=[]; routs=[]
    for wf,wt in o['words']:
        rins.append(u0(kr.parse_word(wf),a))
        if not o['final']: routs.append(u0(b,kr.parse_word(wt)))
    if o['final']:
        for t in o['targets']: routs.append(u0(b,kr.parse_word(lq[(t['n'],t['pair'],t['q'])])))
    o['rin']=max(rins); o['rout']=max(routs); o['pick']=pick
    out.append(o)
    print(f"{o['tp']}->{o['t']} {o['rho']} {'F' if o['final'] else 'N'} step={pick[0]} rin={o['rin']} rout={o['rout']}{' (exit to L_q)' if o['final'] else ''}")
json.dump(out,open(os.path.join(HERE,'steps2.json'),'w'),indent=1,default=str)
import collections
f=lambda v: '0' if v==0 else ('<=1e5' if v<=1e5 else ('<=4e5' if v<=4e5 else '>4e5'))
print('rin',collections.Counter(f(o['rin']) for o in out))
print('rout',collections.Counter(f(o['rout']) for o in out))
print('both<=1e5',sum(1 for o in out if o['rin']<=1e5 and o['rout']<=1e5),'of',len(out))
