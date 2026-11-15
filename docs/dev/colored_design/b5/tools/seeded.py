# seeded joint search for big-type steps: words P·A (A <= 4 rungs, BFS cap) reduced for t' and P·B for t
import sys, os, time, pickle
from collections import deque
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kr
LAM=[3,3,3]
def suffixes(P, t, ml, cap):
    lam=list(LAM); mu=kr.cols(t); n=3; base=kr.top(lam,P)
    q=deque([[]]); seen=0; out=[]
    while q:
        w=q.popleft(); seen+=1
        if seen>cap: break
        tp=kr.top(base,w)
        if tp==mu and kr.reduced_for(lam,P+w,t): out.append(P+w)
        if len(w)>=ml: continue
        for i in range(n-1):
            for d in 'EF':
                g=tp[i+1] if d=='E' else tp[i]
                for r in range(1,g+1):
                    if w and w[-1][0]==d and w[-1][1]==i: continue
                    q.append(w+[(d,i,r)])
    return out
def run(tp_s,t_s,rho_s,Ps,ml=4,cap=25000):
    tp=kr.parse_table(tp_s); t=kr.parse_table(t_s); rho=kr.parse_word(rho_s)[0]; g=kr.step_gf(tp,rho,t)
    best=None
    for Ps_ in Ps:
        P=kr.parse_word(Ps_)
        A=suffixes(P,tp,ml,cap); B=suffixes(P,t,ml+1,cap)
        for a in A:
            for b in B:
                if kr.top(LAM,a+[rho])!=kr.top(LAM,b): continue
                det,ws=kr.step_unknowns(LAM,a,rho,b,g)
                sc=max(min(w,v) if v is not None else w for w,v in det.values())
                if best is None or sc<best[0]: best=(sc,kr.wstr(a),kr.wstr(b),ws['labels'],Ps_)
        print('  P',Ps_,len(A),len(B),best[:1] if best else None, flush=True)
    return best
if __name__=='__main__':
    t0=time.time()
    print(run("1,1,1/1,1,1/1,1,1","2,0,1/2,0,1/1,1,1","E0.2",["F0.1,F1.2","E0.1,F0.1","F0.1,F1.1"]), round(time.time()-t0,1))
