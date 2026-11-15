# G_local check for the joint/seeded pairs: is the window step (after the common prefix P) a step
# of local types m' -> m'' with G_local = G_global? (necessary for a priori completeness of the
# window family; cf. Lemma 3.9 / D1 §1.7 for the Grassmannian prefix)
import sys, os, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kr
HERE=os.path.dirname(os.path.abspath(__file__)); LAM=[3,3,3]
JR=pickle.load(open(os.path.join(HERE,'joint.pkl'),'rb'))
SR=pickle.load(open(os.path.join(HERE,'seeded.pkl'),'rb'))
HR=pickle.load(open(os.path.join(HERE,'hard.pkl'),'rb'))
res=pickle.load(open(os.path.join(HERE,'res.pkl'),'rb'))
def char_type(lab, w):
    ch=kr.forward(lab,w)
    m=max(ch,key=lambda u:kr.lmin(u))
    if ch[m]!={0:1}: return None
    pm=kr.min_perm(m)
    for u in ch:
        if u!=m and not kr.perm_bruhat_le(kr.min_perm(u),pm): return None
    return m
def lab_pad(lab):
    return list(lab)
def _main():
    out={}
    for k in JR:
        x=[y for y in res if (y['tp'],y['t'],y['rho'])==k][0]
        if JR[k]['top']:
            sc,fl,a,b,wl=JR[k]['top'][0]
        else:
            e=HR.get(k,(None,))[0] or SR.get(k)
            if k in HR and HR[k][0] and (SR.get(k) is None or HR[k][0][0]<SR[k][0]): e=HR[k][0]
            else: e=SR.get(k)
            sc=e[0]; a=kr.parse_word(e[1]); b=kr.parse_word(e[2])
        if k in HR and HR[k][0] and HR[k][0][0] < sc:
            e=HR[k][0]; sc=e[0]; a=kr.parse_word(e[1]); b=kr.parse_word(e[2])
        rho=k[2]
        A=a+[rho]; ws=kr.win_split(LAM,A,b)
        if sc==0:
            out[k]=dict(sc=0,ok='Id'); continue
        lab=ws['labels']; aw=ws['a']; bw=ws['b']; f=ws['first']
        if ws['suffix']>0 or not aw or aw[-1]!=(rho[0],rho[1]-f,rho[2]):
            out[k]=dict(sc=sc,ok='n/a (window does not end with rho)',lab=lab); continue
        # bottom table of the window: labels as rows (one block per upright)
        mA=char_type(lab,aw[:-1]); mB=char_type(lab,bw)
        if mA is None or mB is None:
            out[k]=dict(sc=sc,ok='local words not reduced',lab=lab); continue
        gl=kr.step_gf(mA,aw[-1],mB)
        out[k]=dict(sc=sc,ok=(gl==x['g']),gl=gl,g=x['g'],lab=lab,mA=kr.tfmt(mA),mB=kr.tfmt(mB))
    if __name__ == "__main__": pickle.dump(out,open(os.path.join(HERE,"gloc.pkl"),"wb"))
    from collections import Counter
    print(Counter(str(v['ok']) for v in out.values()))
    for k,v in out.items():
        if v['ok'] is not True: print(kr.tfmt(k[0]),kr.tfmt(k[1]),kr.wstr([k[2]]),v)
if __name__=="__main__":
    _main()
