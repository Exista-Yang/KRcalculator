import pickle, sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kr
os.environ['ZB']='1e9'
HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE,'chart2.py')).read().split('def main():')[0]
exec(src)
t0=time.time()
need=[]
for x in res:
    if x['rho'][1] == 1: x, _ = mirror(x, x['rec'])
    zt, lam, kc, o = z_table(x['tp'], x['rho'], x['gr'])
    if (zt, tuple(lam)) not in ZC and (zt,lam) not in need: need.append((zt,lam))
print('missing', len(need))
for zt, lam in need:
    if time.time()-t0 > 12: break
    t1=time.time(); ZC[(zt,tuple(lam))] = bfs_words(lam, zt, 6, 300000, 3); pickle.dump(ZC, open(ZPK,'wb'))
    print(zt, [kr.wstr(w) for w in ZC[(zt,tuple(lam))]], round(time.time()-t1,1))
