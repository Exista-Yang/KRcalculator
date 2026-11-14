# rescoring of the joint (cached) pairs with the G_local = G filter
import sys, os, pickle, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kr
from gloc import char_type
HERE=os.path.dirname(os.path.abspath(__file__)); LAM=[3,3,3]
res=pickle.load(open(os.path.join(HERE,'res.pkl'),'rb'))
WC=pickle.load(open(os.path.join(HERE,'wc.pkl'),'rb'))
JR=pickle.load(open(os.path.join(HERE,'joint.pkl'),'rb'))
FP=os.path.join(HERE,'jfilter.pkl'); FR=pickle.load(open(FP,'rb')) if os.path.exists(FP) else {}
_ct={}
def ct(lab,w):
    key=(tuple(lab),tuple(w))
    if key not in _ct: _ct[key]=char_type(lab,w)
    return _ct[key]
def gl_ok(a,b,rho,g):
    A=a+[rho]; ws=kr.win_split(LAM,A,b)
    if A==b: return True
    lab=ws['labels']; aw=ws['a']; bw=ws['b']; f=ws['first']
    if ws['suffix']>0 or not aw or aw[-1]!=(rho[0],rho[1]-f,rho[2]): return None
    mA=ct(lab,aw[:-1]); mB=ct(lab,bw)
    if mA is None or mB is None: return None
    return kr.step_gf(mA,aw[-1],mB)==g
T0=time.time()
for k,v in JR.items():
    if not v['top'] or k in FR: continue
    if time.time()-T0>12: break
    x=[y for y in res if (y['tp'],y['t'],y['rho'])==k][0]; g=x['g']; tp,t,rho=k
    A=WC[(tp,5)]; B=WC[(t,6)]
    cands=[]
    for a in A:
        for b in B:
            if kr.top(LAM,a+[rho])!=kr.top(LAM,b): continue
            if a+[rho]==b: cands.append((0,a,b)); continue
            det,ws=kr.step_unknowns(LAM,a,rho,b,g)
            cands.append((max(min(w,z) if z is not None else w for w,z in det.values()),a,b))
    cands.sort(key=lambda c:c[0])
    best_ok=None
    for c in cands:
        if c[0]>400000: break
        if gl_ok(c[1],c[2],rho,g):
            best_ok=(c[0],kr.wstr(c[1]),kr.wstr(c[2])); break
    FR[k]=dict(best=cands[0][0], best_ok=best_ok)
    pickle.dump(FR,open(FP,'wb'))
    print(kr.tfmt(tp),kr.tfmt(t),kr.wstr([rho]),'best',cands[0][0],'best G_loc=G',best_ok, round(time.time()-T0,1), flush=True)
print('done',len(FR),'of',sum(1 for v in JR.values() if v['top']))
