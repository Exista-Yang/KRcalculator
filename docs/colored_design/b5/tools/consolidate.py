import sys, os, pickle, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kr
from gloc import char_type
HERE=os.path.dirname(os.path.abspath(__file__)); LAM=[3,3,3]
L=lambda n: pickle.load(open(os.path.join(HERE,n),'rb'))
res=L('res.pkl'); C2=L('chart2.pkl'); JR=L('joint.pkl'); FR=L('jfilter.pkl'); SR=L('seeded.pkl'); SO=L('seeded_ok.pkl'); HR=L('hard.pkl')
c2={}
for c in C2:
    k=c['key']; c2[(kr.parse_table(k[0]),kr.parse_table(k[1]),kr.parse_word(k[3])[0])]=c
def flagtype(lab,w):
    ch=kr.forward(lab,w); return all(sum(v.values())<=1 for v in ch.values())
def closed_form(a,b,rho,g):
    """B1 kind for the window step (a·rho -> b): Pull (X=a·rho flag-type) / Push (X=b flag-type) when
    G is one monomial with coefficient 1 at degree shift(Y)-shift(X)"""
    if a+[rho]==b: return 'Id'
    ws=kr.win_split(LAM,a+[rho],b); lab=ws['labels']; A=ws['a']; B=ws['b']
    if len(g)!=1 or list(g.values())[0]!=1:
        return 'none (|G|>1: Pull/Push+sprinkles [K])'
    d=list(g)[0]
    fa=flagtype(lab,A); fb=flagtype(lab,B)
    sa=kr.hu.shift(list(lab),A); sb=kr.hu.shift(list(lab),B)
    if fa and d==sb-sa: return 'Pull(Out)/Push(In)'
    if fb and d==sa-sb: return 'Push(Out)/Pull(In)'
    return f'none (flag {int(fa)}{int(fb)}, d={d}, shifts {sa},{sb})'
rows={}
for x in res:
    k=(x['tp'],x['t'],x['rho']); r=x['rec']
    e=rows.setdefault(k,dict(x=x,targets=[],words=set()))
    e['targets'].append(dict(n=r['n'],pair=r['pair'],q=r['q'],s=r['s'],mult=r['mult'],rep=r['max_step_unknowns'],nonfinal=r['max_nonfinal'],feasible=r['feasible'],alt=r['alt_final']))
    e['words'].add((kr.wstr(x['wfrom']),kr.wstr(x['wto'])))
out=[]
for k,e in rows.items():
    x=e['x']; c=c2.get(k)
    # census score (max over word variants of the same step)
    scores=[t['rep'] for t in e['targets']]
    j=None; jok=None
    if k in JR and JR[k]['top']:
        sc,fl,a,b,wl=JR[k]['top'][0]; j=(sc,kr.wstr(a),kr.wstr(b),wl)
        bo=FR[k]['best_ok']; jok=bo
    elif k in SR:
        j=SR[k][:4] if SR[k] else None; jok=SO.get(k); jok=jok[:3] if jok else None
    if k in HR and HR[k][0] and (j is None or HR[k][0][0]<j[0]): j=HR[k][0][:4]
    cf=None
    pick = jok if jok else (j[:3] if j else None)
    if pick:
        cf=closed_form(kr.parse_word(pick[1]),kr.parse_word(pick[2]),x['rho'],x['g'])
    if x['score']<=100000 and not any(t['rep']>100000 for t in e['targets']):
        continue
    out.append(dict(tp=kr.tfmt(x['tp']),t=kr.tfmt(x['t']),rho=kr.wstr([x['rho']]),lv=x['lv'],levels=f"{x['level_from']}->{x['level_to']}",
        final=x['final'],G=kr.lp_str(x['g']),census=x['score'],words=sorted(e['words']),
        wslab=x['ws']['labels'],whole=max(w for w,_ in x['det'].values()),
        nk=c['nk'] if c else None, o=c['o'] if c else None, llam=c.get('llam') if c else None, mirrored=c['mirrored'] if c else None,
        loc=c.get('loc') if c else None, chartZ=c.get('Z') if c else None, rin4=c.get('rin') if c else None, rout4=c.get('rout') if c else None,
        rexstep=x['rexstep'], joint=j, joint_ok=jok, cf=cf,
        targets=e['targets']))
json.dump(out,open(os.path.join(HERE,'steps.json'),'w'),indent=1,default=str)
for o in sorted(out,key=lambda o:-o['census']):
    tg=sorted(set((t['n'],t['pair']) for t in o['targets']))
    print(f"{o['tp']}->{o['t']} {o['rho']} lv{o['lv']} {'F' if o['final'] else 'N'} G={o['G']} census={o['census']} win={o['wslab']} nk={o['nk']} o={o['o']} llam={o['llam']} loc4={o['loc']} joint={o['joint'][0] if o['joint'] else None}{o['joint'][3] if o['joint'] else ''} ok={o['joint_ok'][0] if o['joint_ok'] else None} cf={o['cf']} rex={o['rexstep']} tg={tg} ntg={len(o['targets'])}")
print(len(out))
