import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kr
t0=time.time()
zt = ((2,0,0,1),(0,3,0,0),(0,0,1,2),(0,0,0,0))
w = kr.reduced_words([3,3,3,0], zt, 4, 30000)
print(len(w), [kr.wstr(x) for x in w[:5]], time.time()-t0)
