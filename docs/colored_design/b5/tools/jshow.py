import sys, os, pickle, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kr
HERE = os.path.dirname(os.path.abspath(__file__)); LAM=[3,3,3]
res = pickle.load(open(os.path.join(HERE,'res.pkl'),'rb'))
JR = pickle.load(open(os.path.join(HERE,'joint.pkl'),'rb'))
def u0(A,B):
    if A==B: return 0
    ws=kr.win_split(LAM,A,B)
    whole=max(kr.hom_u(LAM,A,B,0),kr.hom_u(LAM,B,A,0))
    if ws['at']>0 or ws['suffix']>0:
        return min(whole, max(kr.hom_u(ws['labels'],ws['a'],ws['b'],0),kr.hom_u(ws['labels'],ws['b'],ws['a'],0)))
    return whole
rows={}
for x in res:
    k=(x['tp'],x['t'],x['rho']); r=x['rec']
    if k not in JR: continue
    e=rows.setdefault(k, dict(x=x, targets=[], qs=set()))
    e['targets'].append((r['n'],r['pair'],r['q'],r['s'],r['mult'])); e['qs'].add((r['n'],r['pair'],r['q'],r['lq'],r['wq']))
out=[]
for k,e in rows.items():
    x=e['x']; j=JR[k]
    if not j['top']:
        out.append(dict(k=k,x=x,targets=e['targets'],joint=None,nA=j['nA'],nB=j['nB'])); continue
    sc,fl,a,b,wl=j['top'][0]
    rin=u0(x['wfrom'],a)
    if x['final']:
        rout=max(u0(b, kr.parse_word(lq)) for (_,_,_,lq,_) in e['qs'])
        rout_kind='exit b->L_q'
    else:
        rout=u0(b,x['wto']); rout_kind='b->W_t'
    out.append(dict(k=k, x=x, targets=e['targets'], joint=sc, a=kr.wstr(a), b=kr.wstr(b), wl=wl, nA=j['nA'], nB=j['nB'], rin=rin, rout=rout, rout_kind=rout_kind, fa=kr.flows(LAM,a), fb=kr.flows(LAM,b)))
pickle.dump(out, open(os.path.join(HERE,'jshow.pkl'),'wb'))
for o in sorted(out, key=lambda o:-o['x']['score']):
    x=o['x']; tg=sorted(set((t[0],t[1]) for t in o['targets']))
    if o['joint'] is None:
        print(f"{kr.tfmt(x['tp'])}->{kr.tfmt(x['t'])} {kr.wstr([x['rho']])} final={x['final']} orig={x['score']} JOINT=none nA={o['nA']} nB={o['nB']} targets={tg} n={len(o['targets'])}"); continue
    tg=sorted(set((t[0],t[1]) for t in o['targets']))
    print(f"{kr.tfmt(x['tp'])}->{kr.tfmt(x['t'])} {kr.wstr([x['rho']])} final={x['final']} orig={x['score']} JOINT={o['joint']} win={o['wl']} a=[{o['a']}]({o['fa']}) b=[{o['b']}]({o['fb']}) nA={o['nA']} nB={o['nB']} rin={o['rin']} rout={o['rout']} ({o['rout_kind']}) targets={tg} n={len(o['targets'])}")
