# joint search over 3-upright reduced words (W'_{t'}, W'_t) for a given step: min step score
import sys, os, time, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kr
LAM=[3,3,3]
def run(tp_s, t_s, rho_s, maxlen=5, cap=250000):
    t0=time.time()
    tp=kr.parse_table(tp_s); t=kr.parse_table(t_s); rho=kr.parse_word(rho_s)[0]
    g=kr.step_gf(tp,rho,t)
    A=kr.reduced_words(LAM,tp,maxlen,cap); B=kr.reduced_words(LAM,t,maxlen+1,cap)
    best=None
    for a in A:
        for b in B:
            if kr.top(LAM,a+[rho])!=kr.top(LAM,b): continue
            det,ws=kr.step_unknowns(LAM,a,rho,b,g)
            sc=max(min(w,x) if x is not None else w for w,x in det.values())
            if best is None or sc<best[0]: best=(sc,kr.wstr(a),kr.wstr(b),ws['labels'])
    return len(A),len(B),best,round(time.time()-t0,1)
cases=[("1,1,1/2,1,0/0,1,2","2,0,1/3,0,0/0,1,2","E0.2"),
       ("3,0,0/1,0,2/1,1,1","2,1,0/0,1,2/1,1,1","F0.2")]
i=int(sys.argv[1])
print(cases[i], run(*cases[i]))
