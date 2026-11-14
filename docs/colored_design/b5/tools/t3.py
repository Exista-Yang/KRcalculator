import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kr
LAM=[3,3,3]
for ts, ml in [("1,1,1/2,1,0/0,1,2",5),("2,0,1/3,0,0/0,1,2",6)]:
    t0=time.time(); t=kr.parse_table(ts)
    w, seen = kr.reduced_words_pruned(LAM, t, ml)
    w2 = kr.reduced_words(LAM, t, ml, 250000)
    print(ts, len(w), seen, round(time.time()-t0,1), 'old', len(w2), set(map(tuple,w2)) <= set(map(tuple,w)))
