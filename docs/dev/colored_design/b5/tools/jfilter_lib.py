# rescoring of the joint (cached) pairs with the G_local = G filter
import sys, os, pickle, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kr
from gloc import char_type
HERE=os.path.dirname(os.path.abspath(__file__)); LAM=[3,3,3]
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
