# Chart words C_K(t) = Z(z)·T_K(m)(·U) of the Grassmannian reduction (D1 §1.7) for the worst
# steps, their flows and reducedness, and crude upper bounds for the rex detours W -> C.
import pickle, sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kr

HERE = os.path.dirname(os.path.abspath(__file__))
res = pickle.load(open(os.path.join(HERE, 'res.pkl'), 'rb'))
LAM = [3, 3, 3]

def pad_table(t):
    return tuple(tuple(list(r) + [0]) for r in t) + (tuple([0] * (len(t[0]) + 1)),)

def zword_cache():
    return {}
ZC = {}

def z_table(tp, rho, gr):
    k = rho[1]; oth = 2 if k == 0 else 0
    q = gr['q']  # merged table: columns (K, other) for k=0 or (other, K) for k=1
    qk = 0 if k == 0 else 1
    kcells = [(i, q[i][qk]) for i in range(3) if q[i][qk] > 0]
    otherrow = [tp[i][oth] for i in range(3)]
    o = sum(otherrow)
    nk = len(kcells)
    if nk == 3:
        ncol = 4; lam = [3, 3, 3, 0]
    else:
        ncol = nk + 1; lam = [3, 3, 3]
    rows = [[0] * ncol for _ in range(len(lam))]
    for a, (i, c) in enumerate(kcells):
        col = a if k == 0 else a + 1
        rows[i][col] = c
    ocol = nk if k == 0 else 0
    for i in range(3):
        rows[i][ocol] += otherrow[i]
    return tuple(tuple(r) for r in rows), lam, kcells, o, ocol

def shift(w, d):
    return [(a, i + d, r) for a, i, r in w]

def main():
    t0 = time.time()
    seen = {}
    out = []
    for x in res:
        r = x['rec']
        key = (kr.tfmt(x['tp']), kr.tfmt(x['t']), kr.wstr(x['wfrom']), kr.wstr([x['rho']]), kr.wstr(x['wto']))
        if key in seen:
            continue
        seen[key] = 1
        gr = x['gr']; rho = x['rho']; k = rho[1]
        zt, lam, kcells, o, ocol = z_table(x['tp'], rho, gr)
        zkey = (zt, tuple(lam))
        if zkey not in ZC:
            ZC[zkey] = kr.reduced_words(lam, zt, 4, 30000)[:3]
        zs = ZC[zkey]
        lb = x['local']['best']
        info = dict(key=key, zt=zt, lam=lam, o=o, nz=len(zs))
        if not zs or lb is None:
            out.append(info); continue
        Z = zs[0]
        ta = kr.parse_word(lb['a']); tb = kr.parse_word(lb['b'])
        nk = len(kcells)
        off = 0 if k == 0 else 1
        U = []
        if nk == 3 and k == 0:
            U = [('E', 2, o)] if o > 0 else []
        if nk == 2:
            # local words live on 2 uprights at the K position
            off = 0 if k == 0 else 1
        Ca = Z + shift(ta, off) + U
        Cb = Z + shift(tb, off) + U
        tpP = pad_table(x['tp']) if nk == 3 else x['tp']
        tP = pad_table(x['t']) if nk == 3 else x['t']
        reda = kr.reduced_for(lam, Ca, tpP)
        redb = kr.reduced_for(lam, Cb, tP)
        WA = x['wfrom']; WB = x['wto']
        fa = kr.flows(lam, Ca); fb = kr.flows(lam, Cb)
        # step window in chart form: C(t')·rho vs C(t) (rho commuted past U when K={0,1})
        rhoC = rho
        stepA = Z + shift(ta, off) + [rhoC] + U
        det, ws = kr.step_unknowns(lam, Z + shift(ta, off), rhoC, Z + shift(tb, off), x['g'])
        step_u = max(min(w, v) if v is not None else w for w, v in det.values())
        # rex detour bounds at d = 0: W_{t'}(·0) -> C(t') and C(t) -> W_t(·0)
        def rexb(A, B):
            d1, w1 = kr.step_unknowns(lam, A[:-1], A[-1], B, {0: 1}) if A else ({0: (kr.hom_u(lam, A, B, 0), None)}, None)
            w, v = d1[0]
            return (w, v)
        rin = rexb(WA, Ca); rout = rexb(Cb, WB)
        info.update(Z=kr.wstr(Z), Ca=kr.wstr(Ca), Cb=kr.wstr(Cb), reda=reda, redb=redb, fa=fa, fb=fb,
                    fwa=kr.flows(LAM, WA), fwb=kr.flows(LAM, WB), step_u=step_u, step_win=ws['labels'],
                    rin=rin, rout=rout, nk=nk, final=x['final'])
        out.append(info)
    print(f"{len(out)} distinct steps, {time.time()-t0:.1f}s", file=sys.stderr)
    pickle.dump(out, open(os.path.join(HERE, 'chart.pkl'), 'wb'))
    for i in out:
        print(i)

main()
