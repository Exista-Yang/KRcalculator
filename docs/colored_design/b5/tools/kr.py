# Python replica of the pieces of krhom needed for the B5 census analysis.
# dcoset.rs: step_terms / forward / reduced_for / min_perm / perm_bruhat_le / grass / reduced_words
# homprov.rs: hom_unknowns (via hu.py logic), win_split.  Read-only re-implementation; no cargo.
import math
from functools import lru_cache

HU = '/Users/existair/Documents/KRcalculator/.claude/worktrees/revrex/docs/colored_design/m4c/tools/c_tools/hu.py'
import importlib.util
_spec = importlib.util.spec_from_file_location('hu', HU)
hu = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(hu)

# ---------------------------------------------------------------- words
def parse_word(s):
    s = s.strip().strip('[]')
    if s in ('', '-'):
        return []
    out = []
    for t in s.split(','):
        d = t[0]; i, r = t[1:].split('.')
        out.append((d, int(i), int(r)))
    return out

def wstr(w):
    return ','.join(f"{d}{i}.{r}" for d, i, r in w)

def rg(r):
    d, i, _ = r
    return (i, i + 1) if d == 'E' else (i + 1, i)

def apply(lab, r):
    lab = list(lab); rec, giv = rg(r)
    if giv >= len(lab) or rec >= len(lab) or lab[giv] < r[2] or r[2] <= 0:
        return None
    lab[giv] -= r[2]; lab[rec] += r[2]
    return lab

def top(lab, w):
    for r in w:
        lab = apply(lab, r)
        if lab is None:
            return None
    return lab

def level(lab, w, k):
    return top(lab, w[:k])

def flows(lab, w):
    return hu.flows(list(lab), w)

def hom_u(lab, a, b, d):
    return hu.hu(list(lab), a, b, d)

# ---------------------------------------------------------------- tables (tuple of row tuples)
def parse_table(s):
    return tuple(tuple(int(x) for x in row.split(',')) for row in s.split('/'))

def tfmt(t):
    return '/'.join(','.join(str(x) for x in row) for row in t)

def rows(t):
    return [sum(r) for r in t]

def cols(t):
    return [sum(t[i][j] for i in range(len(t))) for j in range(len(t[0]))]

def c2(n):
    return n * (n - 1) // 2

def lmin(t):
    n, m = len(t), len(t[0]); s = 0
    for i in range(n):
        for j in range(m):
            a = t[i][j]
            if not a:
                continue
            for i2 in range(i):
                for j2 in range(j + 1, m):
                    s += a * t[i2][j2]
    return s

def lw(comp):
    return sum(c2(a) for a in comp)

def l_lr(t):
    return sum(c2(a) for row in t for a in row)

def lmax(t):
    return lmin(t) + lw(rows(t)) + lw(cols(t)) - l_lr(t)

def identity(lam):
    n = len(lam)
    return tuple(tuple(lam[i] if i == j else 0 for j in range(n)) for i in range(n))

@lru_cache(None)
def gauss(n, k):
    if k < 0 or k > n:
        return (0,)
    if k == 0 or k == n:
        return (1,)
    a = gauss(n - 1, k - 1); b = gauss(n - 1, k)
    out = [0] * (k * (n - k) + 1)
    for i, x in enumerate(a): out[i] += x
    for i, x in enumerate(b): out[i + k] += x
    return tuple(out)

# Laurent polys as dict exp -> coef
def lp_mul(a, b):
    o = {}
    for e1, c1 in a.items():
        for e2, c2_ in b.items():
            o[e1 + e2] = o.get(e1 + e2, 0) + c1 * c2_
    return {e: c for e, c in o.items() if c}

def lp_add(a, b):
    o = dict(a)
    for e, c in b.items():
        o[e] = o.get(e, 0) + c
    return {e: c for e, c in o.items() if c}

def lp_shift(a, s):
    return {e + s: c for e, c in a.items()}

def lp_str(a):
    return ' + '.join((f"{c}" if c != 1 else '') + f"v^{e}" for e, c in sorted(a.items()))

def step_terms(t, r):
    n, m = len(t), len(t[0])
    rec, giv = rg(r)
    if giv >= m or rec >= m:
        return []
    col = [t[i][giv] for i in range(n)]
    xs = []
    def rec_x(i, left, x):
        if i == n:
            if left == 0:
                xs.append(list(x))
            return
        for v in range(min(col[i], left) + 1):
            x.append(v); rec_x(i + 1, left - v, x); x.pop()
    rec_x(0, r[2], [])
    out = []
    for x in xs:
        if r[0] == 'E':
            first, second = x, [col[i] - x[i] for i in range(n)]
        else:
            first, second = [col[i] - x[i] for i in range(n)], x
        u = []
        for i in range(n):
            row = []
            for j in range(m):
                if j == giv:
                    row += [first[i], second[i]]
                else:
                    row.append(t[i][j])
            u.append(tuple(row))
        u = tuple(u)
        te = [list(row) for row in t]
        for i in range(n):
            te[i][giv] -= x[i]; te[i][rec] += x[i]
        to = tuple(tuple(row) for row in te)
        dd = lmax(t) - lmax(u)
        da = lmin(to) - lmin(u) - l_lr(to) + l_lr(u)
        poly = {0: 1}
        for i in range(n):
            a = to[i][rec]
            if x[i] > 0 and a > x[i]:
                g = gauss(a, x[i])
                poly = lp_mul(poly, {2 * j: c for j, c in enumerate(g) if c})
        out.append((to, lp_shift(poly, dd + da)))
    return out

def step_gf(t, r, t_to):
    for to, g in step_terms(t, r):
        if to == t_to:
            return g
    return None

def forward(lam, w):
    cur = {identity(lam): {0: 1}}
    for r in w:
        nxt = {}
        for t, c in cur.items():
            for u, g in step_terms(t, r):
                nxt[u] = lp_add(nxt.get(u, {}), lp_mul(c, g))
        cur = {k: v for k, v in nxt.items() if v}
    return cur

def min_perm(t):
    n, m = len(t), len(t[0])
    rs = rows(t); nxt = []; off = 0
    for i in range(n):
        nxt.append(off); off += rs[i]
    out = []
    for j in range(m):
        for i in range(n):
            for _ in range(t[i][j]):
                out.append(nxt[i]); nxt[i] += 1
    return out

def perm_bruhat_le(u, w):
    k = len(u)
    for j in range(k):
        cu = cw = 0
        for i in range(k):
            cu += u[i] >= j; cw += w[i] >= j
            if cu > cw:
                return False
    return True

def perm_comp(u, v):
    return [u[a] for a in v]

def perm_inv(u):
    o = [0] * len(u)
    for a, b in enumerate(u):
        o[b] = a
    return o

def table_of_perm(w, rws, cls):
    n, m = len(rws), len(cls)
    blk = []
    for i, a in enumerate(rws):
        blk += [i] * a
    e = [[0] * m for _ in range(n)]
    pos = 0
    for j, b in enumerate(cls):
        for _ in range(b):
            e[blk[w[pos]]][j] += 1; pos += 1
    return tuple(tuple(r) for r in e)

def reduced_for(lam, w, t):
    tp = top(list(lam), w)
    if tp is None or tp != cols(t):
        return False
    ch = forward(lam, w)
    if ch.get(t) != {0: 1}:
        return False
    mp = min_perm(t)
    return all(perm_bruhat_le(min_perm(u), mp) for u in ch)

def merge_columns(t, j):
    return tuple(tuple(list(row[:j]) + [row[j] + row[j + 1]] + list(row[j + 2:])) for row in t)

def rr(t):
    n, m = len(t), len(t[0])
    return [t[i][j] for j in range(m) for i in range(n) if t[i][j] > 0]

def grass(t_from, t_to, r):
    k = r[1]
    q = merge_columns(t_from, k)
    assert q == merge_columns(t_to, k)
    l = rr(q)
    qmin = min_perm(q); qinv = perm_inv(qmin)
    lam = rows(t_from)
    z = table_of_perm(qmin, lam, l)
    mf = table_of_perm(perm_comp(qinv, min_perm(t_from)), l, cols(t_from))
    mt = table_of_perm(perm_comp(qinv, min_perm(t_to)), l, cols(t_to))
    # cells of q in column-major order -> which ones are in the merged column k
    cellpos = [(i, j) for j in range(len(q[0])) for i in range(len(q)) if q[i][j] > 0]
    krows = [idx for idx, (i, j) in enumerate(cellpos) if j == k]
    local_lam = [q[i][k] for i in range(len(q)) if q[i][k] > 0]
    mK_from = tuple((mf[a][k], mf[a][k + 1]) for a in krows)
    mK_to = tuple((mt[a][k], mt[a][k + 1]) for a in krows)
    # sanity: rows of mf outside K stay in their own column
    return dict(q=q, l=l, z=z, mf=mf, mt=mt, local_lam=local_lam, mK_from=mK_from, mK_to=mK_to, krows=krows, cellpos=cellpos)

def reduced_words(lam, t, max_len, cap=200000):
    """all words of <= max_len rungs reduced for t (BFS), ascending flows"""
    from collections import deque
    lam = list(lam); mu = cols(t); n = len(lam)
    q = deque([[]]); seen = 0; out = []
    while q:
        w = q.popleft(); seen += 1
        if seen > cap:
            break
        tp = top(lam, w)
        if tp == mu and reduced_for(lam, w, t):
            out.append(w)
        if len(w) >= max_len:
            continue
        for i in range(n - 1):
            for d in 'EF':
                giver = tp[i + 1] if d == 'E' else tp[i]
                for rr_ in range(1, giver + 1):
                    q.append(w + [(d, i, rr_)])
    out.sort(key=lambda w: flows(lam, w))
    return out

# ---------------------------------------------------------------- win_split (homprov.rs:197)
def win_split(lam, a, b):
    pre = 0
    while pre < len(a) and pre < len(b) and a[pre] == b[pre]:
        pre += 1
    suf = 0
    while suf < len(a) - pre and suf < len(b) - pre and a[len(a) - 1 - suf] == b[len(b) - 1 - suf]:
        suf += 1
    wa, wb = a[pre:len(a) - suf], b[pre:len(b) - suf]
    lv = list(lam)
    for r in a[:pre]:
        lv = apply(lv, r)
    if not wa and not wb:
        return dict(at=pre, suffix=suf, first=0, labels=[], a=[], b=[])
    lo = min(r[1] for r in wa + wb); hi = max(r[1] + 1 for r in wa + wb)
    labels = lv[lo:hi + 1]
    sh = lambda rs: [(d, i - lo, r) for d, i, r in rs]
    return dict(at=pre, suffix=suf, first=lo, labels=labels, a=sh(wa), b=sh(wb))

def step_unknowns(lam, wfrom, rho, wto, g):
    """step_score per degree: dict d -> (whole, window or None)"""
    a = wfrom + [rho]
    ws = win_split(lam, a, wto)
    res = {}
    for d in sorted(g):
        whole = max(hom_u(lam, a, wto, d), hom_u(lam, wto, a, d))
        win = None
        if ws['at'] > 0 or ws['suffix'] > 0:
            win = max(hom_u(ws['labels'], ws['a'], ws['b'], d), hom_u(ws['labels'], ws['b'], ws['a'], d))
        res[d] = (whole, win)
    return res, ws
