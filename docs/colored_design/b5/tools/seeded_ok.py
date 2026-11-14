import sys, os, time, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kr
from seeded import suffixes
from jfilter_lib import gl_ok
HERE=os.path.dirname(os.path.abspath(__file__)); LAM=[3,3,3]
JR=pickle.load(open(os.path.join(HERE,'joint.pkl'),'rb'))
C2=pickle.load(open(os.path.join(HERE,'chart2.pkl'),'rb'))
SP=os.path.join(HERE,'seeded_ok.pkl'); SR=pickle.load(open(SP,'rb')) if os.path.exists(SP) else {}
def mt(t): return tuple(tuple(reversed(r)) for r in reversed(t))
def mw(w): return [('F' if d=='E' else 'E',1-i,r) for d,i,r in w]
zof={}
for c in C2:
    k=c['key']; zof[(kr.parse_table(k[0]),kr.parse_table(k[1]),kr.parse_word(k[3])[0])]=(c.get('Z'),c['mirrored'])
T0=time.time()
for k,v in JR.items():
    if v['top'] or k in SR: continue
    if time.time()-T0>11: break
    tp,t,rho=k; Z,mir=zof[k]
    Zw=kr.parse_word(Z)
    if mir: tp2,t2,rho2=mt(tp),mt(t),mw([rho])[0]
    else: tp2,t2,rho2=tp,t,rho
    P=[r for r in Zw if r[1]!=2]
    Ps=[]
    for cand in [P, P[:-1], Zw[:2]]:
        if all(r[1]<2 for r in cand) and cand not in Ps and cand: Ps.append(cand)
    g=kr.step_gf(tp2,rho2,t2); best=None
    for Pw in Ps:
        A=suffixes(Pw,tp2,4,15000); B=suffixes(Pw,t2,5,15000)
        for a in A:
            for b in B:
                if kr.top(LAM,a+[rho2])!=kr.top(LAM,b): continue
                det,ws=kr.step_unknowns(LAM,a,rho2,b,g)
                sc=max(min(w,x) if x is not None else w for w,x in det.values())
                if sc<=400000 and (best is None or sc<best[0]) and gl_ok(a,b,rho2,g):
                    aa,bb=(mw(a),mw(b)) if mir else (a,b)
                    best=(sc,kr.wstr(aa),kr.wstr(bb),ws['labels'],kr.wstr(mw(Pw) if mir else Pw),len(A),len(B))
    SR[k]=best; pickle.dump(SR,open(SP,'wb'))
    print(kr.tfmt(tp),kr.tfmt(t),kr.wstr([rho]),'mir' if mir else '', best, round(time.time()-T0,1), flush=True)
print('done', sum(1 for k,v in JR.items() if not v['top'] and k in SR), 'of', sum(1 for v in JR.values() if not v['top']))
