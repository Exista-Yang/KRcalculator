# replica of homprov::hom_unknowns / fdeg_hist / word_shift and schur::emono_count, schur::flows
import sys, math
from functools import lru_cache
def parse(w):
    if w in ('', '-'): return []
    out = []
    for t in w.split(','):
        d = t[0]; i, r = t[1:].split('.')
        out.append((d, int(i), int(r)))
    return out
def ends(r):
    d, i, _ = r
    return (i, i+1) if d == 'E' else (i+1, i)   # (receiver, giver)
def apply(lab, r):
    lab = list(lab); rec, giv = ends(r)
    assert lab[giv] >= r[2] and r[2] > 0, (lab, r)
    lab[giv] -= r[2]; lab[rec] += r[2]; return lab
@lru_cache(None)
def gb(n, k):
    if k < 0 or k > n: return (0,)
    if k == 0 or k == n: return (1,)
    a = gb(n-1, k-1); b = gb(n-1, k)
    out = [0]*max(len(a), len(b)+k)
    for i, x in enumerate(a): out[i] += x
    for i, x in enumerate(b): out[i+k] += x
    return tuple(out)
def fdeg_hist(lab, rs):
    h = [1]
    for r in rs:
        rec, giv = ends(r); g = gb(lab[giv], r[2])
        out = [0]*(len(h)+len(g)-1)
        for i, x in enumerate(h):
            for j, y in enumerate(g): out[i+j] += x*y
        h = out; lab = apply(lab, r)
    return h
def shift(lab, rs):
    s = 0
    for r in rs:
        rec, giv = ends(r); s -= r[2]*(lab[giv]-r[2]); lab = apply(lab, r)
    return s
def emono(lam, w):
    ser = [0]*(w+1); ser[0] = 1
    for a in lam:
        for j in range(1, a+1):
            for d in range(j, w+1): ser[d] += ser[d-j]
    return ser[w]
def flows(lab, rs):
    f = 1
    for r in rs:
        rec, giv = ends(r); f *= math.comb(lab[giv], r[2]); lab = apply(lab, r)
    return f
def top(lab, rs):
    for r in rs: lab = apply(lab, r)
    return lab
def hu(lab, src, tgt, d):
    u = d - shift(list(lab), tgt) + shift(list(lab), src)
    hs = fdeg_hist(list(lab), src); ht = fdeg_hist(list(lab), tgt)
    tot = 0
    for i, a in enumerate(hs):
        if not a: continue
        for j, b in enumerate(ht):
            if not b: continue
            dd = u + 2*i - 2*j
            if dd < 0 or dd % 2: continue
            tot += a*b*emono(lab, dd//2)
    return tot
if __name__ == '__main__':
    lab = [int(x) for x in sys.argv[1].split(',')]
    src, tgt = parse(sys.argv[2]), parse(sys.argv[3])
    assert top(list(lab), src) == top(list(lab), tgt), (top(list(lab), src), top(list(lab), tgt))
    for d in [int(x) for x in sys.argv[4].split(',')]:
        print(f"lab={lab} src={sys.argv[2]}({flows(list(lab),src)} fl) tgt={sys.argv[3]}({flows(list(lab),tgt)} fl) d={d}: unknowns={hu(lab,src,tgt,d)}")
