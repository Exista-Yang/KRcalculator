import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kr
from collections import deque
def bfs(lam, t, maxlen, cap):
    lam=list(lam); mu=kr.cols(t); n=len(lam)
    q=deque([[]]); seen=0; out=[]
    while q:
        w=q.popleft(); seen+=1
        if seen>cap: break
        tp=kr.top(lam,w)
        if tp==mu and kr.reduced_for(lam,w,t):
            out.append(w); 
            if len(out)>=3: break
            continue
        if len(w)>=maxlen: continue
        for i in range(n-1):
            for d in 'EF':
                g = tp[i+1] if d=='E' else tp[i]
                for r in range(1,g+1):
                    q.append(w+[(d,i,r)])
    return out, seen
t0=time.time()
zt=((1,2,0,0),(1,0,2,0),(1,0,0,2),(0,0,0,0))
w,seen=bfs([3,3,3,0],zt,6,400000)
print([kr.wstr(x) for x in w], seen, time.time()-t0)
