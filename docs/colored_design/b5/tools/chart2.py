# Chart words C_K(t) = Z(z)·T_K(m)(·U) (D1 §1.7) for every distinct worst step; choose Z and the
# local words to minimize the max of (rex detour in, local step, rex detour out); rex detours are
# bounded by one Exact Hom^0 call (min(whole, window), both directions) -- an upper bound for the
# rexplan chain. Final steps: the detour out is replaced by the exit chain C_K(q) -> L_q(·0).
import pickle, sys, os, time
from collections import deque
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kr

HERE = os.path.dirname(os.path.abspath(__file__))
res = pickle.load(open(os.path.join(HERE, 'res.pkl'), 'rb'))
LAM = [3, 3, 3]

def pad_table(t):
    return tuple(tuple(list(r) + [0]) for r in t) + (tuple([0] * (len(t[0]) + 1)),)

def bfs_words(lam, t, maxlen, cap, want=4):
    lam = list(lam); mu = kr.cols(t); n = len(lam)
    q = deque([[]]); seen = 0; out = []; best_len = None
    while q:
        w = q.popleft(); seen += 1
        if seen > cap or (best_len is not None and len(w) > best_len + 1):
            break
        tp = kr.top(lam, w)
        if tp == mu and kr.reduced_for(lam, w, t):
            out.append(w); best_len = best_len or len(w)
            if len(out) >= want:
                break
            continue
        if len(w) >= maxlen:
            continue
        for i in range(n - 1):
            for d in 'EF':
                g = tp[i + 1] if d == 'E' else tp[i]
                for r in range(1, g + 1):
                    q.append(w + [(d, i, r)])
    out.sort(key=lambda w: kr.flows(lam, w))
    return out

ZPK = os.path.join(HERE, 'zc.pkl')
ZC = pickle.load(open(ZPK, 'rb')) if os.path.exists(ZPK) else {}
ZBUDGET = [float(os.environ.get('ZB', '14'))]
T0 = time.time()
def z_words(zt, lam):
    key = (zt, tuple(lam))
    if key not in ZC:
        if time.time() - T0 > ZBUDGET[0]:
            raise SystemExit('ZBUDGET')
        ZC[key] = bfs_words(lam, zt, 6, 300000, 3)
        pickle.dump(ZC, open(ZPK, 'wb'))
    return ZC[key]

def z_table(tp, rho, gr):
    k = rho[1]; oth = 2 if k == 0 else 0
    q = gr['q']; qk = 0 if k == 0 else 1
    kcells = [(i, q[i][qk]) for i in range(3) if q[i][qk] > 0]
    otherrow = [tp[i][oth] for i in range(3)]
    o = sum(otherrow); nk = len(kcells)
    if nk == 3:
        ncol = 4; lam = [3, 3, 3, 0]
    else:
        ncol = 3; lam = [3, 3, 3]
    rows = [[0] * ncol for _ in range(len(lam))]
    for a, (i, c) in enumerate(kcells):
        rows[i][a if k == 0 else a + 1] += c
    ocol = (nk if k == 0 else 0) if nk == 3 else (2 if k == 0 else 0)
    for i in range(3):
        rows[i][ocol] += otherrow[i]
    return tuple(tuple(r) for r in rows), lam, kcells, o

def sh(w, d):
    return [(a, i + d, r) for a, i, r in w]

def u0(lam, A, B, degs):
    """max over degs of min(whole, window), both directions"""
    if A == B:
        return 0
    ws = kr.win_split(lam, A, B); mx = 0
    for d in degs:
        whole = max(kr.hom_u(lam, A, B, d), kr.hom_u(lam, B, A, d))
        u = whole
        if ws['at'] > 0 or ws['suffix'] > 0:
            u = min(whole, max(kr.hom_u(ws['labels'], ws['a'], ws['b'], d), kr.hom_u(ws['labels'], ws['b'], ws['a'], d)))
        mx = max(mx, u)
    return mx

def local_words(x, maxlen=5):
    gr = x['gr']; ll = list(gr['local_lam']); k = len(ll)
    mf, mt = gr['mK_from'], gr['mK_to']
    if k == 2:
        lam = ll; tf = tuple(tuple(r) for r in mf); tt = tuple(tuple(r) for r in mt)
    else:
        lam = ll; tf = tuple(tuple(r) + (0,) for r in mf); tt = tuple(tuple(r) + (0,) for r in mt)
    return lam, kr.reduced_words(lam, tf, maxlen, 60000)[:12], kr.reduced_words(lam, tt, maxlen, 60000)[:12]

def flag_type(lam, w):
    ch = kr.forward(lam, w)
    return all(sum(c.values()) <= 1 for c in ch.values()), ch

def mt(t):
    return tuple(tuple(reversed(r)) for r in reversed(t))
def mw(w):
    return [('F' if d == 'E' else 'E', 1 - i, r) for d, i, r in w]
def mirror(x, r):
    y = dict(x)
    y['tp'] = mt(x['tp']); y['t'] = mt(x['t']); y['rho'] = mw([x['rho']])[0]
    y['wfrom'] = mw(x['wfrom']); y['wto'] = mw(x['wto'])
    y['gr'] = kr.grass(y['tp'], y['t'], y['rho'])
    g2 = kr.step_gf(y['tp'], y['rho'], y['t'])
    assert g2 == x['g'], (g2, x['g'])
    return y, kr.wstr(mw(kr.parse_word(r['lq'])))

def main():
    t0 = time.time()
    seen = {}; out = []
    for x in res:
        r = x['rec']
        key = (kr.tfmt(x['tp']), kr.tfmt(x['t']), kr.wstr(x['wfrom']), kr.wstr([x['rho']]), kr.wstr(x['wto']))
        if key in seen:
            seen[key]['targets'].append((r['n'], r['pair'], r['q'], r['s'], r['mult']))
            continue
        mirrored = False
        lq_s = r['lq']
        if x['rho'][1] == 1:
            x, lq_s = mirror(x, r); mirrored = True
        rho = x['rho']; k = rho[1]
        zt, lam, kcells, o = z_table(x['tp'], rho, x['gr'])
        nk = len(kcells)
        zs = z_words(zt, lam)
        llam, wfs, wts = local_words(x)
        g = x['g']; degs = sorted(g)
        tpP = pad_table(x['tp']) if nk == 3 else x['tp']
        tP = pad_table(x['t']) if nk == 3 else x['t']
        qP = None
        lq = kr.parse_word(lq_s)
        WA = x['wfrom']; WB = x['wto']
        rl = (rho[0], 0, rho[2])
        best = None
        for Z in zs[:3]:
            off = 0 if k == 0 else 1
            U = [('E', 2, o)] if (nk == 3 and k == 0 and o > 0) else []
            for ta in wfs[:8]:
                Ca = Z + sh(ta, off) + U
                if not kr.reduced_for(lam, Ca, tpP):
                    continue
                rin = u0(lam, WA, Ca, [0])
                for tb in wts[:8]:
                    if kr.top(list(llam), ta + [rl]) != kr.top(list(llam), tb):
                        continue
                    loc = 0 if ta + [rl] == tb else u0(llam, ta + [rl], tb, degs)
                    Cb = Z + sh(tb, off) + U
                    if x['final']:
                        rout = u0(lam, Cb, lq, [0])  # exit chain bound (one Exact)
                    else:
                        rout = u0(lam, Cb, WB, [0])
                    cand = (max(rin, loc, rout), rin, loc, rout, Z, ta, tb, Ca, Cb)
                    if best is None or cand[0] < best[0]:
                        best = cand
        info = dict(key=key, targets=[(r['n'], r['pair'], r['q'], r['s'], r['mult'])], x=x, nk=nk, o=o, zt=zt, nz=len(zs), lam=lam,
                    final=x['final'], score=x['score'], g=g, mirrored=mirrored)
        if best:
            mx, rin, loc, rout, Z, ta, tb, Ca, Cb = best
            fta, _ = flag_type(llam, ta + [rl]); ftb, _ = flag_type(llam, tb)
            info.update(best=mx, rin=rin, loc=loc, rout=rout, Z=kr.wstr(Z), ta=kr.wstr(ta), tb=kr.wstr(tb), Ca=kr.wstr(Ca), Cb=kr.wstr(Cb),
                        fa=kr.flows(lam, Ca), fb=kr.flows(lam, Cb), flagA=fta, flagB=ftb,
                        shA=kr.hu.shift(list(llam), ta + [rl]), shB=kr.hu.shift(list(llam), tb), llam=llam)
        seen[key] = info; out.append(info)
    print(f"{len(out)} distinct steps, {time.time()-t0:.1f}s, {len(ZC)} z tables", file=sys.stderr)
    for i in out:
        i.pop('x', None)
    pickle.dump(out, open(os.path.join(HERE, 'chart2.pkl'), 'wb'))
    for i in out:
        print({k: v for k, v in i.items() if k not in ('targets',)}, len(i['targets']))

main()
