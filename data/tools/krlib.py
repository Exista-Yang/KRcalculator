"""Small, dependency-free knot utilities used to normalise and cross-check the
KnotInfo / LinkInfo tables.

Conventions (same as krhom's src/pd.rs and tests/pd.rs):
  * PD X[i,j,k,l]: labels counterclockwise from the incoming under-strand (i -> k is the
    under-strand).  With i=S, j=E, k=N, l=W the crossing is positive iff the over-strand
    runs l -> j.  Orientation is found by walking along components; the direction of a
    component is fixed by its under-passages.
  * Kauffman bracket: A-smoothing of X[i,j,k,l] joins (i,j),(k,l); B-smoothing (i,l),(j,k).
    V = (-A^3)^{-w} <D> / delta at A = t^{-1/4};  V(right trefoil) = t + t^3 - t^4.
  * Braids: +g = sigma_g, strand at position g passes OVER strand g+1 (positive crossing);
    closure as in krhom::pd::braid_to_pd.
Polynomials are dicts {exponent(s): coefficient}.
"""
from fractions import Fraction
import re

# ---------------------------------------------------------------------------
# Laurent polynomials in one variable: dict int -> int
# ---------------------------------------------------------------------------

def padd(a, b, s=1):
    r = dict(a)
    for e, c in b.items():
        v = r.get(e, 0) + s * c
        if v:
            r[e] = v
        else:
            r.pop(e, None)
    return r


def pmul(a, b):
    r = {}
    for e1, c1 in a.items():
        for e2, c2 in b.items():
            e = e1 + e2
            r[e] = r.get(e, 0) + c1 * c2
    return {e: c for e, c in r.items() if c}


def pshift(a, k, s=1):
    return {e + k: s * c for e, c in a.items()}


def mirror1(a):
    return {-e: c for e, c in a.items()}


DELTA_A = {2: -1, -2: -1}  # -A^2 - A^-2


def div_delta(p):
    """exact division by delta = -A^2 - A^-2"""
    p = dict(p)
    if not p:
        return {}
    lo = min(p)
    q = {}
    while p:
        top = max(p)
        if top < lo + 2:
            raise ValueError("not divisible by delta")
        c = p.pop(top)
        q[top - 2] = -c
        v = p.get(top - 4, 0) - c
        if v:
            p[top - 4] = v
        else:
            p.pop(top - 4, None)
    return q


# ---------------------------------------------------------------------------
# PD codes
# ---------------------------------------------------------------------------

def ints(s):
    return [int(x) for x in re.findall(r'-?\d+', s)]


def parse_pd(s):
    v = ints(s)
    if len(v) % 4:
        raise ValueError("PD length not multiple of 4")
    return [v[i:i + 4] for i in range(0, len(v), 4)]


def orient(x):
    """Walk components. Returns (signs, comp_of_label, incoming) where incoming[(c,s)]
    is True if the edge in slot s enters crossing c.  Components are numbered in order
    of their smallest label."""
    occ = {}
    for c, xc in enumerate(x):
        for s in range(4):
            occ.setdefault(xc[s], []).append((c, s))
    for l, o in occ.items():
        if len(o) != 2:
            raise ValueError(f"label {l} occurs {len(o)} times")
    incoming = {}
    comp = {}
    nc = 0
    for l0 in sorted(occ):
        if l0 in comp:
            continue
        walk = []
        l, at = l0, occ[l0][0]
        while True:
            walk.append((l, at))
            out = (at[0], at[1] ^ 2)
            l2 = x[out[0]][out[1]]
            o = occ[l2]
            at = o[1] if o[0] == out else o[0]
            l = l2
            if l == l0 and at == occ[l0][0]:
                break
        fwd = any(w[1][1] == 0 for w in walk)
        bwd = any(w[1][1] == 2 for w in walk)
        if fwd and bwd:
            raise ValueError("inconsistent orientation")
        if fwd or bwd:
            rev = bwd
        else:
            m = len(walk)
            inc = sum(1 for k in range(m) if walk[(k + 1) % m][0] == walk[k][0] + 1)
            dec = sum(1 for k in range(m) if walk[k][0] == walk[(k + 1) % m][0] + 1)
            rev = dec > inc
        for (l, at) in walk:
            comp[l] = nc
            o = occ[l]
            other = o[1] if o[0] == at else o[0]
            incoming[at] = not rev
            incoming[other] = rev
        nc += 1
    signs = [1 if incoming[(c, 3)] else -1 for c in range(len(x))]
    for c in range(len(x)):
        assert incoming[(c, 0)] and not incoming[(c, 2)]
        assert incoming[(c, 1)] != incoming[(c, 3)]
    # renumber components by smallest label (walk started at sorted labels: already so)
    return signs, comp, incoming


def num_components(x):
    return len(set(orient(x)[1].values()))


def linking_matrix(x):
    signs, comp, _ = orient(x)
    n = len(set(comp.values()))
    m = [[0] * n for _ in range(n)]
    for c, xc in enumerate(x):
        a, b = comp[xc[0]], comp[xc[1]]
        if a != b:
            m[a][b] += signs[c]
            m[b][a] += signs[c]
    for i in range(n):
        for j in range(n):
            assert m[i][j] % 2 == 0
            m[i][j] //= 2
    return m


def _add_arc(m, u, v):
    if u == v:
        return True
    a, b = u, v
    w = m.pop(a, None)
    if w is not None:
        m.pop(w)
        if w == b:
            return True
        a = w
    w = m.pop(b, None)
    if w is not None:
        m.pop(w)
        if w == a:
            return True
        b = w
    m[a] = b
    m[b] = a
    return False


def bracket(x):
    """Kauffman bracket, every loop = delta (so unknot diagram gives delta)."""
    states = {(): {0: 1}}
    openc = {}
    done = [False] * len(x)
    for _ in range(len(x)):
        best = None
        for c in range(len(x)):
            if done[c]:
                continue
            k = sum(1 for l in x[c] if l in openc)
            if best is None or k > best[0]:
                best = (k, c)
        c = best[1]
        done[c] = True
        for l in x[c]:
            openc[l] = openc.get(l, 0) + 1
        for l in [l for l, v in openc.items() if v >= 2]:
            del openc[l]
        i, j, k, l = x[c]
        nxt = {}
        for key, coef in states.items():
            for arcs, aexp in (((i, j), (k, l)), 1), (((i, l), (j, k)), -1):
                m = {}
                for a, b in key:
                    m[a] = b
                    m[b] = a
                loops = 0
                for u, v in arcs:
                    if _add_arc(m, u, v):
                        loops += 1
                nk = tuple(sorted((a, b) for a, b in m.items() if a < b))
                t = pshift(coef, aexp)
                for _ in range(loops):
                    t = pmul(t, DELTA_A)
                nxt[nk] = padd(nxt.get(nk, {}), t)
        states = nxt
    assert len(states) == 1 and () in states
    return states[()]


def jones_pd(x):
    """Jones polynomial as dict {k: c} meaning c * t^(k/2)."""
    signs, _, _ = orient(x)
    w = sum(signs)
    b = div_delta(bracket(x))
    f = pshift(b, -3 * w, -1 if w % 2 else 1)
    out = {}
    for k, c in f.items():
        assert k % 2 == 0
        out[-k // 2] = c
    return out


# ---------------------------------------------------------------------------
# signature by Gordon-Litherland
# ---------------------------------------------------------------------------

def faces(x):
    occ = {}
    for c, xc in enumerate(x):
        for s in range(4):
            occ.setdefault(xc[s], []).append((c, s))

    def other(h):
        o = occ[x[h[0]][h[1]]]
        return o[1] if o[0] == h else o[0]

    face_of = {}
    nf = 0
    for c in range(len(x)):
        for s in range(4):
            if (c, s) in face_of:
                continue
            h = (c, s)
            while h not in face_of:
                face_of[h] = nf
                h2 = other(h)
                h = (h2[0], (h2[1] + 1) % 4)
            nf += 1
    return face_of, nf


def signature_symmetric(M):
    n = len(M)
    A = [[Fraction(v) for v in row] for row in M]
    idx = list(range(n))
    pos = neg = 0
    while idx:
        piv = next((i for i in idx if A[i][i] != 0), None)
        if piv is None:
            pair = next(((i, j) for i in idx for j in idx if i != j and A[i][j] != 0), None)
            if pair is None:
                break
            i, j = pair
            # row/col i += row/col j
            for k in range(n):
                A[i][k] += A[j][k]
            for k in range(n):
                A[k][i] += A[k][j]
            piv = i
            assert A[i][i] != 0
        p = A[piv][piv]
        if p > 0:
            pos += 1
        else:
            neg += 1
        idx.remove(piv)
        for a in idx:
            f = A[a][piv] / p
            if f:
                for b in idx:
                    A[a][b] -= f * A[piv][b]
        for a in idx:
            A[a][piv] = A[piv][a] = Fraction(0)
    return pos - neg


def signature_pd(x, color=0):
    """Signature with KnotInfo's convention (positive trefoil -> -2), by the
    Gordon-Litherland formula using the faces of colour `color` as white regions."""
    signs, _, _ = orient(x)
    face_of, nf = faces(x)
    if nf != len(x) + 2:
        raise ValueError(f"diagram not connected/planar: {nf} faces for {len(x)} crossings")
    # 2-colour the faces: at each crossing corners with departure slot 1,3 share a colour,
    # 0,2 share the other.
    parent = list(range(nf))
    par = [0] * nf

    def find(a):
        if parent[a] == a:
            return a, 0
        r, p = find(parent[a])
        parent[a] = r
        par[a] ^= p
        return r, par[a]

    def union(a, b, d):
        ra, pa = find(a)
        rb, pb = find(b)
        if ra == rb:
            if (pa ^ pb) != d:
                raise ValueError("not 2-colourable")
            return
        parent[ra] = rb
        par[ra] = pa ^ pb ^ d

    for c in range(len(x)):
        f = [face_of[(c, s)] for s in range(4)]
        union(f[1], f[3], 0)
        union(f[0], f[2], 0)
        union(f[0], f[1], 1)
    col = [find(a)[1] for a in range(nf)]
    white = [a for a in range(nf) if col[a] == color]
    wi = {a: i for i, a in enumerate(white)}
    G = [[0] * len(white) for _ in white]
    mu = 0
    for c in range(len(x)):
        # corners with departure slot d: A-regions are d in {0,2}
        white_is_A = col[face_of[(c, 0)]] == color
        eta = -1 if white_is_A else 1
        if white_is_A:
            fa, fb = face_of[(c, 0)], face_of[(c, 2)]
        else:
            fa, fb = face_of[(c, 1)], face_of[(c, 3)]
        if fa != fb:
            a, b = wi[fa], wi[fb]
            G[a][b] -= eta
            G[b][a] -= eta
            G[a][a] += eta
            G[b][b] += eta
        positive = signs[c] > 0
        # oriented smoothing = A-smoothing iff positive; it joins the A-regions.
        type2 = (white_is_A != positive)
        if type2:
            mu += eta
    Gr = [row[1:] for row in G[1:]]
    return signature_symmetric(Gr) - mu


# ---------------------------------------------------------------------------
# braids
# ---------------------------------------------------------------------------

def braid_to_pd(n, word):
    """Port of krhom::pd::braid_to_pd (closing arcs on the right)."""
    next_id = n
    cur = list(range(n))
    xs = []
    pos = []
    for g in word:
        i = abs(g) - 1
        assert 0 <= i and i + 1 < n, (g, n)
        l, r = cur[i], cur[i + 1]
        ne, nw = next_id, next_id + 1
        next_id += 2
        if g > 0:
            xs.append([r, ne, nw, l])
        else:
            xs.append([l, r, ne, nw])
        pos.append(g > 0)
        cur[i] = nw
        cur[i + 1] = ne
    rep = list(range(next_id))
    for p in range(n):
        assert cur[p] != p, "strand without crossings"
        rep[cur[p]] = p
    xs = [[rep[e] for e in x] for x in xs]
    nxt = [None] * next_id
    for c, x in enumerate(xs):
        nxt[x[0]] = x[2]
        if pos[c]:
            nxt[x[3]] = x[1]
        else:
            nxt[x[1]] = x[3]
    label = [0] * next_id
    lab = 0
    for p in range(n):
        if label[p]:
            continue
        e = p
        while True:
            lab += 1
            label[e] = lab
            e = nxt[e]
            if e == p:
                break
    return [[label[e] for e in x] for x in xs]


def braid_components(n, word):
    at = list(range(n))  # at[position] = strand (strand s starts at bottom position s)
    for g in word:
        i = abs(g) - 1
        at[i], at[i + 1] = at[i + 1], at[i]
    top = [0] * n
    for p, s in enumerate(at):
        top[s] = p
    seen = [False] * n
    nc = 0
    for p in range(n):
        if not seen[p]:
            nc += 1
            q = p
            while not seen[q]:
                seen[q] = True
                q = top[q]
    return nc


# ---------------------------------------------------------------------------
# polynomial expression parser (Mathematica-like / KnotInfo strings)
# ---------------------------------------------------------------------------

_TOK = re.compile(r'\s*(?:(\d+)|([A-Za-z]+)|(.))')


class MPoly:
    """multivariate Laurent polynomial with Fraction exponents; dict tuple->int"""

    def __init__(self, nv, terms=None):
        self.nv = nv
        self.t = terms or {}

    @staticmethod
    def const(nv, c):
        return MPoly(nv, {(Fraction(0),) * nv: c} if c else {})

    def __add__(self, o):
        r = dict(self.t)
        for e, c in o.t.items():
            v = r.get(e, 0) + c
            if v:
                r[e] = v
            else:
                r.pop(e, None)
        return MPoly(self.nv, r)

    def neg(self):
        return MPoly(self.nv, {e: -c for e, c in self.t.items()})

    def __mul__(self, o):
        r = {}
        for e1, c1 in self.t.items():
            for e2, c2 in o.t.items():
                e = tuple(a + b for a, b in zip(e1, e2))
                r[e] = r.get(e, 0) + c1 * c2
        return MPoly(self.nv, {e: c for e, c in r.items() if c})

    def pow(self, k):
        if len(self.t) == 1:
            (e, c), = self.t.items()
            if k.denominator != 1 and c != 1:
                raise ValueError("fractional power of non-monic monomial")
            cc = Fraction(c) ** k if k.denominator == 1 else 1
            if cc.denominator != 1:
                raise ValueError("non-integer coefficient")
            return MPoly(self.nv, {tuple(a * k for a in e): int(cc)})
        if k.denominator != 1 or k < 0:
            raise ValueError("power of polynomial")
        r = MPoly.const(self.nv, 1)
        for _ in range(int(k)):
            r = r * self
        return r

    def inv(self):
        if len(self.t) != 1:
            raise ValueError("division by non-monomial")
        (e, c), = self.t.items()
        if c not in (1, -1):
            raise ValueError("division by non-unit")
        return MPoly(self.nv, {tuple(-a for a in e): c})


def parse_poly(s, vars_):
    """Parse an expression in the variables vars_ (e.g. ['t'] or ['v','z']).
    Supports + - * / ^, parentheses, Sqrt[...], implicit multiplication."""
    toks = []
    for m in _TOK.finditer(s):
        if m.group(1):
            toks.append(('n', int(m.group(1))))
        elif m.group(2):
            toks.append(('v', m.group(2)))
        elif m.group(3) and not m.group(3).isspace():
            toks.append(('o', m.group(3)))
    nv = len(vars_)
    pos = [0]

    def peek():
        return toks[pos[0]] if pos[0] < len(toks) else (None, None)

    def take():
        t = toks[pos[0]]
        pos[0] += 1
        return t

    def expr():
        sign = 1
        r = MPoly(nv)
        first = True
        while True:
            k, v = peek()
            if k == 'o' and v in '+-':
                take()
                sign = 1 if v == '+' else -1
            elif not first:
                break
            t = term()
            r = r + (t if sign > 0 else t.neg())
            sign = 1
            first = False
            k, v = peek()
            if not (k == 'o' and v in '+-'):
                break
        return r

    def term():
        r = factor()
        while True:
            k, v = peek()
            if k == 'o' and v == '*':
                take()
                r = r * factor()
            elif k == 'o' and v == '/':
                take()
                r = r * factor().inv()
            elif k in ('n', 'v') or (k == 'o' and v == '('):
                r = r * factor()
            else:
                break
        return r

    def exponent():
        k, v = peek()
        if k == 'o' and v == '(':
            take()
            sgn = 1
            if peek() == ('o', '-'):
                take()
                sgn = -1
            elif peek() == ('o', '+'):
                take()
            num = take()[1]
            den = 1
            if peek() == ('o', '/'):
                take()
                den = take()[1]
            assert take() == ('o', ')')
            return Fraction(sgn * num, den)
        sgn = 1
        if k == 'o' and v == '-':
            take()
            sgn = -1
        k, v = take()
        assert k == 'n', s
        return Fraction(sgn * v)

    def factor():
        k, v = take()
        if k == 'n':
            b = MPoly.const(nv, v)
        elif k == 'v' and v == 'Sqrt':
            assert take() == ('o', '[')
            inner = expr()
            assert take() == ('o', ']')
            b = inner.pow(Fraction(1, 2))
        elif k == 'v':
            if v not in vars_:
                raise ValueError(f"unknown variable {v} in {s[:60]}")
            e = [Fraction(0)] * nv
            e[vars_.index(v)] = Fraction(1)
            b = MPoly(nv, {tuple(e): 1})
        elif k == 'o' and v == '(':
            b = expr()
            assert take() == ('o', ')')
        else:
            raise ValueError(f"unexpected token {v!r} in {s[:80]}")
        if peek() == ('o', '^'):
            take()
            b = b.pow(exponent())
        return b

    r = expr()
    if pos[0] != len(toks):
        raise ValueError(f"trailing tokens in {s[:80]}: {toks[pos[0]:pos[0]+5]}")
    return r


def split_top(s, sep=','):
    """split on separators at bracket depth 0"""
    out, depth, cur = [], 0, []
    for ch in s:
        if ch in '([{':
            depth += 1
        elif ch in ')]}':
            depth -= 1
        if ch == sep and depth == 0:
            out.append(''.join(cur))
            cur = []
        else:
            cur.append(ch)
    out.append(''.join(cur))
    return out


# ---------------------------------------------------------------------------
# HOMFLY-PT polynomial by the descending-diagram skein algorithm
# convention (KnotInfo/LinkInfo): v^-1 P(L+) - v P(L-) = z P(L0), P(unknot) = 1
# polynomials: dict (v_exp, z_exp) -> int
# ---------------------------------------------------------------------------

def _p2add(a, b):
    r = dict(a)
    for e, c in b.items():
        v = r.get(e, 0) + c
        if v:
            r[e] = v
        else:
            r.pop(e, None)
    return r


def _p2mul(a, b):
    r = {}
    for (a1, b1), c1 in a.items():
        for (a2, b2), c2 in b.items():
            e = (a1 + a2, b1 + b2)
            r[e] = r.get(e, 0) + c1 * c2
    return {e: c for e, c in r.items() if c}


_MU = {(-1, -1): 1, (1, -1): -1}  # (v^-1 - v)/z


def _mu_pow(k):
    r = {(0, 0): 1}
    for _ in range(k):
        r = _p2mul(r, _MU)
    return r


class Budget(Exception):
    pass


def _homfly_rec(cr, free, memo, budget=None):
    key = (cr, free)
    if key in memo:
        return memo[key]
    if budget is not None and len(memo) > budget:
        raise Budget()
    if not cr:
        res = _mu_pow(free - 1)
        memo[key] = res
        return res
    head = {}
    for idx, (ui, uo, oi, oo, s) in enumerate(cr):
        head[ui] = (idx, 0)
        head[oi] = (idx, 1)
    seen_l = set()
    seen_c = set()
    ncomp = 0
    bad = None
    for l0 in sorted(head):
        if l0 in seen_l:
            continue
        ncomp += 1
        e = l0
        while e not in seen_l:
            seen_l.add(e)
            idx, role = head[e]
            if idx not in seen_c:
                seen_c.add(idx)
                if role == 0:
                    bad = idx
                    break
            e = cr[idx][1] if role == 0 else cr[idx][3]
        if bad is not None:
            break
    if bad is None:
        res = _mu_pow(ncomp + free - 1)
        memo[key] = res
        return res
    ui, uo, oi, oo, s = cr[bad]
    sw = cr[:bad] + ((oi, oo, ui, uo, -s),) + cr[bad + 1:]
    ren = {}

    def r(a):
        while a in ren:
            a = ren[a]
        return a

    a, b = r(ui), r(oo)
    if a != b:
        ren[b] = a
    a, b = r(oi), r(uo)
    if a != b:
        ren[b] = a
    rest = tuple(tuple(r(e) for e in c[:4]) + (c[4],) for k, c in enumerate(cr) if k != bad)
    present = set()
    for c in rest:
        present.update(c[:4])
    loops = len({r(ui), r(oi)} - present)
    p_sw = _homfly_rec(sw, free, memo, budget)
    p_sm = _homfly_rec(rest, free + loops, memo, budget)
    if s > 0:   # D = L+ : P+ = v^2 P- + v z P0
        res = _p2add(_p2mul({(2, 0): 1}, p_sw), _p2mul({(1, 1): 1}, p_sm))
    else:       # D = L- : P- = v^-2 P+ - v^-1 z P0
        res = _p2add(_p2mul({(-2, 0): 1}, p_sw), _p2mul({(-1, 1): -1}, p_sm))
    memo[key] = res
    return res


def homfly_pd(x, budget=None):
    """HOMFLY-PT polynomial; with a budget (max memo size) returns None when exceeded"""
    signs, comp, incoming = orient(x)
    cr = []
    for c, X in enumerate(x):
        i, j, k, l = X
        if signs[c] > 0:
            cr.append((i, k, l, j, 1))
        else:
            cr.append((i, k, j, l, -1))
    try:
        return _homfly_rec(tuple(cr), 0, {}, budget)
    except Budget:
        return None


def mirror_homfly(p):
    return {(-a, b): c * (-1) ** (b % 2) for (a, b), c in p.items()}
