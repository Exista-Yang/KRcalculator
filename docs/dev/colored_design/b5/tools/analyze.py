import json, re, sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kr

B2 = '/Users/existair/Documents/KRcalculator/.claude/worktrees/revrex/docs/colored_design/m4c/b2/'
LAM = [3, 3, 3]
RX = re.compile(r"level (\d+) (\S+) -> (\S+) \[([^\]]*)\]·(\S+) -> \[([^\]]*)\]")

def load():
    recs = []
    for fn in ['gp0_n5b.jsonl', 'gp0_n6.jsonl']:
        for line in open(B2 + fn):
            recs.append(json.loads(line))
    return recs

_local_cache = {}
def local_window(gr, rho, g, maxlen=5):
    """unknowns of the rung-local K-window s in Hom(T'_K·rho_loc, T''_K), minimized over reduced
    local words of <= maxlen rungs (max over d in supp G, both directions; window trick)"""
    ll = list(gr['local_lam']); k = len(ll)
    mf, mt = gr['mK_from'], gr['mK_to']
    if k == 1:
        lam = [ll[0], 0]
        tf = (tuple(mf[0]), (0, 0)); tt = (tuple(mt[0]), (0, 0))
    elif k == 2:
        lam = ll; tf = tuple(tuple(r) for r in mf); tt = tuple(tuple(r) for r in mt)
    else:
        lam = ll; tf = tuple(tuple(r) + (0,) for r in mf); tt = tuple(tuple(r) + (0,) for r in mt)
    rl = (rho[0], 0, rho[2])
    key = (tuple(lam), tf, tt, rl, tuple(sorted(g.items())))
    if key in _local_cache:
        return _local_cache[key]
    gl = kr.step_gf(tf, rl, tt)
    wf = kr.reduced_words(lam, tf, maxlen, 60000)
    wt = kr.reduced_words(lam, tt, maxlen, 60000)
    best = None
    for a in wf[:40]:
        for b in wt[:40]:
            if kr.top(list(lam), a + [rl]) != kr.top(list(lam), b):
                continue
            if a + [rl] == b:
                sc = 0; det = {}
            else:
                det, ws = kr.step_unknowns(lam, a, rl, b, g)
                sc = max(min(w, x) if x is not None else w for w, x in det.values())
            cand = (sc, kr.flows(lam, a + [rl]) + kr.flows(lam, b), a, b)
            if best is None or cand[:2] < best[:2]:
                best = cand
    res = dict(lam=lam, tf=tf, tt=tt, g_local=gl, g_ok=(gl == g), nwf=len(wf), nwt=len(wt),
               best=None if best is None else dict(u=best[0], a=kr.wstr(best[2]), b=kr.wstr(best[3])))
    _local_cache[key] = res
    return res

def classify(gr, level_from):
    ll = gr['local_lam']; k = len(ll); nk = sum(ll)
    if k <= 2:
        return 'a'  # rung-local on 3 uprights (Z level = K row parts + other column = <= 3 blocks)
    if nk < 9:
        return 'a-pad'  # 3 row parts in K + nonempty other column: 4-block Z level, needs Pad
    return 'b'  # K column holds all 9 roots: local = global

def main():
    t0 = time.time()
    recs = load()
    out = []
    for r in recs:
        if r['rex_class']:
            continue
        m = RX.match(r['worst_step'])
        if not m:
            out.append(dict(rec=r, err='no worst step'))
            continue
        lv = int(m.group(1)); tp = kr.parse_table(m.group(2)); t = kr.parse_table(m.group(3))
        wfrom = kr.parse_word(m.group(4)); rho = kr.parse_word(m.group(5))[0]; wto = kr.parse_word(m.group(6))
        g = kr.step_gf(tp, rho, t)
        det, ws = kr.step_unknowns(LAM, wfrom, rho, wto, g)
        score = max(min(w, x) if x is not None else w for w, x in det.values())
        final = (kr.tfmt(t) == r['q'])
        # alt final word (gp0 re-chose W_q): recompute the step into wq
        alt = None
        if final and r['alt_final']:
            wq = kr.parse_word(r['wq'])
            det2, ws2 = kr.step_unknowns(LAM, wfrom, rho, wq, g)
            alt = dict(word=r['wq'], det=det2, ws=ws2, score=max(min(w, x) if x is not None else w for w, x in det2.values()))
        gr = kr.grass(tp, t, rho)
        cat = classify(gr, kr.cols(tp))
        # rex-type step: G = 1 and W_{t'}·rho reduced for t
        rexstep = (g == {0: 1}) and kr.reduced_for(LAM, wfrom + [rho], t)
        lw = local_window(gr, rho, g)
        out.append(dict(rec=r, lv=lv, tp=tp, t=t, wfrom=wfrom, rho=rho, wto=wto, g=g, det=det, ws=ws, score=score,
                        final=final, alt=alt, gr=gr, cat=cat, rexstep=rexstep, local=lw,
                        level_from=kr.cols(tp), level_to=kr.cols(t)))
    print(f"analyzed {len(out)} non-rex targets in {time.time()-t0:.1f}s", file=sys.stderr)
    return out

if __name__ == '__main__':
    import pickle
    res = main()
    pickle.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'res.pkl'), 'wb'))
    for x in res:
        r = x['rec']
        if 'err' in x:
            print(r['n'], r['pair'], r['q'], r['s'], 'ERR'); continue
        gr = x['gr']
        print(f"N{r['n']} p{r['pair']} q {r['q']} s {r['s']} | rep {r['max_step_unknowns']} nonfinal {r['max_nonfinal']} | lv {x['lv']} {kr.tfmt(x['tp'])}->{kr.tfmt(x['t'])} "
              f"[{kr.wstr(x['wfrom'])}]·{kr.wstr([x['rho']])}->[{kr.wstr(x['wto'])}] G={kr.lp_str(x['g'])} score={x['score']} final={x['final']} "
              f"win={x['ws']['labels']} [{kr.wstr(x['ws']['a'])}]->[{kr.wstr(x['ws']['b'])}] det={x['det']} | Q={kr.tfmt(gr['q'])} local={gr['local_lam']} "
              f"cat={x['cat']} rex={x['rexstep']} localwin={x['local']['best']} gok={x['local']['g_ok']}"
              + (f" | alt {x['alt']['word']} score {x['alt']['score']} win {x['alt']['ws']['labels']}" if x['alt'] else ''))
