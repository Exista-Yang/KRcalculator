"""Predicted colored gl_N homology of T(2,2m+1) (J. Wang, arXiv:2604.25174, Thm 1.2;
m=1 also arXiv:2211.08409 Prop 1.5) over Q, converted to krhom conventions:
ours(t,q) = Wang(1/t, 1/q).  H**(B_lambda) = Tor_{H*(BU(N))}(Q, H*(BJ_lambda)) computed
as Koszul homology mod p.  Usage: wang.py k N m  (prints Poincare polynomial, rank, chi)."""
import sys, itertools
from collections import defaultdict
P = 2147483647

def partitions(maxparts, maxpart):
    def rec(i, cap):
        if i == maxparts:
            yield ()
            return
        for v in range(cap, -1, -1):
            for rest in rec(i + 1, v):
                yield (v,) + rest
    for p in rec(0, maxpart):
        yield tuple(x for x in p if x > 0)

def pmul(a, b):
    out = defaultdict(int)
    for ma, ca in a.items():
        for mb, cb in b.items():
            m = tuple(x + y for x, y in zip(ma, mb))
            out[m] = (out[m] + ca * cb) % P
    return {m: c for m, c in out.items() if c}

def tor_poincare(groups, N):
    """groups: list of (size, multiplicity) meaning U(size) embedded `mult` times diagonally.
    Returns dict (i, w) -> dim Tor_i in internal weight w (q-degree 2w)."""
    # variables: for each group, e_1..e_size with weights 1..size
    wts = []
    gvars = []
    for size, mult in groups:
        idx = []
        for j in range(1, size + 1):
            idx.append(len(wts)); wts.append(j)
        gvars.append((idx, mult))
    nv = len(wts)
    zero = tuple([0] * nv)
    # c(t) = prod_g E_g(t)^mult ; coefficients as polynomials
    # represent polynomial in t: list of dict polys
    ct = [{zero: 1}]
    for idx, mult in gvars:
        E = [{zero: 1}]
        for v in idx:
            m = list(zero); m[v] = 1
            E.append({tuple(m): 1})
        for _ in range(mult):
            new = [dict() for _ in range(len(ct) + len(E) - 1)]
            for i, a in enumerate(ct):
                for j, b in enumerate(E):
                    pr = pmul(a, b)
                    acc = new[i + j]
                    for m, c in pr.items():
                        acc[m] = (acc.get(m, 0) + c) % P
            ct = [{m: c for m, c in d.items() if c} for d in new]
    assert len(ct) == N + 1, (len(ct), N)
    f = [dict(d) for d in ct]  # f[j] has weight j
    # eliminate f_j that are linear in an unused weight-j variable (regular, P/(f_j) polynomial)
    alive = [True] * nv
    seq = []
    for j in range(1, N + 1):
        x = None
        for v in range(nv):
            if alive[v] and wts[v] == j:
                mv = [0] * nv; mv[v] = 1; mv = tuple(mv)
                if f[j].get(mv, 0):
                    x = v; c = f[j][mv]; break
        if x is None:
            seq.append(j); continue
        inv = pow(c, P - 2, P)
        s = {m: (P - cf) * inv % P for m, cf in f[j].items() if m[x] == 0}  # x = s
        pows = [{zero: 1}]
        for l in range(j + 1, N + 1):
            newf = defaultdict(int)
            for m, cf in f[l].items():
                e = m[x]
                while len(pows) <= e: pows.append(pmul(pows[-1], s))
                base = list(m); base[x] = 0; base = tuple(base)
                for mp, cp in pows[e].items():
                    mm = tuple(a + b for a, b in zip(base, mp))
                    newf[mm] = (newf[mm] + cf * cp) % P
            f[l] = {m: cf for m, cf in newf.items() if cf}
        alive[x] = False
    dimR = N * N - sum(s * s for s, _ in groups)
    wmax = (dimR + len(seq)) // 2 + 1
    monos = defaultdict(list)
    def gen(i, w, cur):
        if i == nv:
            monos[w].append(tuple(cur)); return
        if not alive[i]:
            cur.append(0); gen(i + 1, w, cur); cur.pop(); return
        for e in range((wmax - w) // wts[i] + 1):
            cur.append(e); gen(i + 1, w + e * wts[i], cur); cur.pop()
    gen(0, 0, [])
    L = len(seq)
    subsets = {i: list(itertools.combinations(seq, i)) for i in range(L + 1)}
    def basis(i, w):
        out = []
        for S in subsets.get(i, []):
            ws = w - sum(S)
            if ws < 0: continue
            for m in monos.get(ws, []):
                out.append((m, S))
        return out
    def rank_d(i, w):
        if i == 0: return 0
        src = basis(i, w)
        if not src: return 0
        tgt = {b: n for n, b in enumerate(basis(i - 1, w))}
        rows = []
        for (m, S) in src:
            row = defaultdict(int)
            for l, s_ in enumerate(S):
                sign = 1 if l % 2 == 0 else P - 1
                T = S[:l] + S[l + 1:]
                for mf, cf in f[s_].items():
                    mm = tuple(x_ + y_ for x_, y_ in zip(m, mf))
                    key = tgt[(mm, T)]
                    row[key] = (row[key] + sign * cf) % P
            rows.append({k: v for k, v in row.items() if v})
        return rank_mod(rows)
    res = {}
    cache = {}
    for w in range(0, wmax + 1):
        for i in range(0, L + 1):
            dimK = len(basis(i, w))
            if dimK == 0: continue
            if (i, w) not in cache: cache[(i, w)] = rank_d(i, w)
            if (i + 1, w) not in cache: cache[(i + 1, w)] = rank_d(i + 1, w)
            ri = cache[(i, w)]; ri1 = cache[(i + 1, w)]
            h = dimK - ri - ri1
            if h: res[(i, w)] = h
    return res

def rank_mod(rows):
    piv = {}
    r = 0
    for row in rows:
        row = dict(row)
        while row:
            k = min(row)
            if k in piv:
                prow = piv[k]; c = row[k]
                for kk, vv in prow.items():
                    row[kk] = (row.get(kk, 0) - c * vv) % P
                    if row[kk] == 0: del row[kk]
            else:
                inv = pow(row[k], P - 2, P)
                piv[k] = {kk: vv * inv % P for kk, vv in row.items()}
                r += 1
                break
    return r

def predict(k, N, m):
    if 2 * k > N: k = N - k
    tot = defaultdict(int)  # (t, q) in Wang's convention
    for lam in partitions(k, m):
        r = len(lam)
        mults = defaultdict(int)
        for x in lam: mults[x] += 1
        groups = [(rr, 2) for _, rr in sorted(mults.items())]
        groups += [(k - r, 1), (N - k - r, 1)]
        groups = [g for g in groups if g[0] > 0]
        tor = tor_poincare(groups, N)
        sh_t = -2 * sum(lam)
        sh_q = 4 * sum((i + 1) * x for i, x in enumerate(lam)) + k * (N - k) * (2 * m - 1)
        for (i, w), d in tor.items():
            tot[(sh_t - i, sh_q + 2 * w)] += d
    ours = {(-t, -q): d for (t, q), d in tot.items()}
    return ours

if __name__ == "__main__":
    k, N, m = map(int, sys.argv[1:4])
    ours = predict(k, N, m)
    terms = sorted(ours.items())
    print(" + ".join(f"{d if d>1 else ''}t^{t}q^{q}" for (t, q), d in terms))
    chi = defaultdict(int)
    for (t, q), d in ours.items(): chi[q] += d * (-1) ** (t % 2)
    print("rank", sum(ours.values()))
    print("chi", " ".join(f"{c:+d}q^{q}" for q, c in sorted(chi.items()) if c))
