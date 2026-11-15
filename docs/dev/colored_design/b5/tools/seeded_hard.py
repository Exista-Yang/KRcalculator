import sys, os, time, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kr
from seeded import suffixes
HERE=os.path.dirname(os.path.abspath(__file__)); LAM=[3,3,3]
res=pickle.load(open(os.path.join(HERE,'res.pkl'),'rb'))
C2=pickle.load(open(os.path.join(HERE,'chart2.pkl'),'rb'))
HP=os.path.join(HERE,'hard.pkl'); HR=pickle.load(open(HP,'rb')) if os.path.exists(HP) else {}
def mt(t): return tuple(tuple(reversed(r)) for r in reversed(t))
def mw(w): return [('F' if d=='E' else 'E',1-i,r) for d,i,r in w]
zof={}
for c in C2:
    k=c['key']; zof[(kr.parse_table(k[0]),kr.parse_table(k[1]),kr.parse_word(k[3])[0])]=(c.get('Z'),c['mirrored'])
hard=[("3,0,0/2,0,1/1,0,2","2,1,0/1,1,1/0,1,2","F0.3"),("1,1,1/0,2,1/0,2,1","2,0,1/1,1,1/0,2,1","E0.2"),
      ("1,2,0/1,1,1/1,1,1","1,2,0/1,0,2/1,1,1","F1.1"),("1,1,1/0,2,1/0,2,1","1,1,1/1,1,1/1,1,1","E0.2"),
      ("1,2,0/2,1,0/1,1,1","1,1,1/2,0,1/1,1,1","F1.2")]
T0=time.time()
for hs in hard:
    k=(kr.parse_table(hs[0]),kr.parse_table(hs[1]),kr.parse_word(hs[2])[0])
    if k in HR: continue
    if time.time()-T0>8: break
    x=[x for x in res if (x['tp'],x['t'],x['rho'])==k][0]
    Z,mir=zof[k]; Zw=kr.parse_word(Z)
    tp,t,rho=k
    if mir: tp2,t2,rho2=mt(tp),mt(t),mw([rho])[0]; WA,WB=mw(x['wfrom']),mw(x['wto'])
    else: tp2,t2,rho2=tp,t,rho; WA,WB=x['wfrom'],x['wto']
    Ps=[]
    for cand in [[r for r in Zw if r[1]!=2], Zw[:1], WA[:1], WA[:2], WB[:1], WB[:2], WA[:3]]:
        if cand and all(r[1]<2 for r in cand) and cand not in Ps: Ps.append(cand)
    g=kr.step_gf(tp2,rho2,t2); best=None
    for Pw in Ps:
        A=suffixes(Pw,tp2,4,15000); B=suffixes(Pw,t2,5,15000)
        for a in A:
            for b in B:
                if kr.top(LAM,a+[rho2])!=kr.top(LAM,b): continue
                det,ws=kr.step_unknowns(LAM,a,rho2,b,g)
                sc=max(min(w,v) if v is not None else w for w,v in det.values())
                if best is None or sc<best[0]:
                    aa,bb=(mw(a),mw(b)) if mir else (a,b)
                    best=(sc,kr.wstr(aa),kr.wstr(bb),ws['labels'],kr.wstr(mw(Pw) if mir else Pw))
    HR[k]=(best,[kr.wstr(p) for p in Ps]); pickle.dump(HR,open(HP,'wb'))
    print(hs, best, len(Ps), round(time.time()-T0,1), flush=True)
print('done', len(HR))
