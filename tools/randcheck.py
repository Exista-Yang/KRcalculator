#!/usr/bin/env python3
"""randcheck.py: seeded random-example checker for the krhom calculator (CALC-1, spec section E
as amended by CALC_DECISIONS.md and CALC_SPEC_review.md).

usage: randcheck.py --seed S [--n 24 | --quick | --full] [--tiers 2,3,4] [--only ID]
                    [--krhom target/release/krhom] [--jobs 2] [--out runs_random] [--inject-fault]
                    [--no-anchor]

Draws random colored braids (color 2 on every strand) on 2, 3 and 4 strands, runs `krhom --json`
on them and on related words, and checks relations fixed in code before running (no tolerances):

  label C = certifying (a mismatch is FAIL), B = blocking (a mismatch is FAIL, a pass is not
  counted as confirmation), E = expected but not proven (a mismatch is WARN).

  O1  two primes agree (built into krhom: exit 3 = FAIL)                              C
  O2  seed independence: box runs use --seed 201 --seed-f32b S2 != 201                 C
  O3  chi = decat_chi recorded for every closure of every prime (box routes)           C
  O4  isotopy: rotation, conjugation, braid relation, cancelling pair: equal           C
  O5  Markov: a stabilization sigma_n^{+-1} of a core word: equal to the core          C
  O6  mirror: terms (t,q,k) -> (-t,-q,k)                                               C
  O7  reversal: equal                                                                  E
  O8  N'=3 duality: colored-2 N'=3 of w = uncolored N=3 of rev(w), same mode            B
  O9  N'=2 is exactly t^0q^0, both modes                                               C
  O10b Lewark-type: [N 2]_q P_red - P_unred = (1+t) Q with Q >= 0 (knots)              E
  O11 2-strand closed forms (Wang T(2,2m+1) via tools/wang.py, unknot, unlink)        C (B for m=4, m=3 N>=5)
  O12 unknot words: unreduced [N 2]_q, reduced t^0q^0                                  C
  O13a 2-strand cross-engine: --box-live = exact path (N' <= 4, both modes)           C

Tiers: T2 (2 strands, exact path), T3 (3 strands, library l222_n4), T4 (4 strands, library
l2222_n4, reduced only plus unreduced for length <= 3, coverage-filtered with the m5i_cover
example). Box-route calls always pass --no-simplify.

Report: OUT/SEED/examples.jsonl (one line per example), OUT/SEED/summary.md, and a final line
"seed=... n=... PASS a FAIL b WARN c SKIP d TIMEOUT e; checks: ...".
Exit 0 if no FAIL or ERROR, 1 on any FAIL or ERROR, 2 on a harness failure.
"""
import argparse
import concurrent.futures as cf
import json
import os
import random
import re
import resource
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# ------------------------------------------------------------------ polynomials (dict (t,q) -> k)

def poly_of(pj):
    return {(t, q): k for t, q, k in pj["terms"]}

def pstr(p):
    if not p:
        return "0"
    return " + ".join(("" if k == 1 else str(k)) + "t^%dq^%d" % (t, q) for (t, q), k in sorted(p.items()))

def parse_pstr(s):
    s = s.strip()
    if s == "0" or not s:
        return {}
    out = {}
    for term in s.split(" + "):
        m = re.fullmatch(r"(-?\d*)t\^(-?\d+)q\^(-?\d+)", term.strip())
        if not m:
            raise ValueError("bad term %r" % term)
        k = int(m.group(1)) if m.group(1) not in ("", "-") else (1 if m.group(1) == "" else -1)
        out[(int(m.group(2)), int(m.group(3)))] = k
    return out

def mirror_poly(p):
    return {(-t, -q): k for (t, q), k in p.items()}

def pdiff(a, b):
    keys = sorted(set(a) | set(b))
    d = [(k, a.get(k, 0), b.get(k, 0)) for k in keys if a.get(k, 0) != b.get(k, 0)]
    return "; ".join("t^%dq^%d: %d vs %d" % (k[0], k[1], x, y) for k, x, y in d[:8]) + (" ..." if len(d) > 8 else "")

def qbinom2(n):
    """[N 2]_q in krhom conventions (symmetric, steps of q^2), as {q: coeff}."""
    if n < 2:
        return {}
    # Gaussian binomial in v = q^2: coefficients of [n choose 2]_v
    c = [0] * (2 * (n - 2) + 1)
    for i in range(n):
        for j in range(i + 1, n):
            c[i + j - 1] += 1  # sum over 0<=i<j<n of v^(i+j-1)
    lo = -2 * (n - 2)
    return {lo + 2 * e: k for e, k in enumerate(c) if k}

def qpoly_t0(qp):
    return {(0, q): k for q, k in qp.items()}

def qmul(qp, p):
    out = {}
    for q1, c1 in qp.items():
        for (t, q), k in p.items():
            out[(t, q + q1)] = out.get((t, q + q1), 0) + c1 * k
    return {x: v for x, v in out.items() if v}

def lewark(n, pred, punred):
    """[N 2]_q P_red - P_unred = (1+t) Q, Q >= 0: returns None if it holds, else a detail."""
    d = qmul(qbinom2(n), pred)
    for k, v in punred.items():
        d[k] = d.get(k, 0) - v
    byq = {}
    for (t, q), k in d.items():
        if k:
            byq.setdefault(q, {})[t] = k
    for q, row in byq.items():
        ts = range(min(row), max(row) + 2)
        prev = 0
        for t in ts:
            cur = row.get(t, 0) - prev
            if cur < 0:
                return "negative Q coefficient %d at t^%d q^%d" % (cur, t, q)
            prev = cur
        if prev != 0:
            return "not divisible by (1+t) at q^%d" % q
    return None

# ------------------------------------------------------------------ braid words

def ncomp(word, n):
    at = list(range(n))
    for g in word:
        i = abs(g) - 1
        at[i], at[i + 1] = at[i + 1], at[i]
    perm = [0] * n
    for pos, s in enumerate(at):
        perm[s] = pos
    seen, c = [False] * n, 0
    for p in range(n):
        if not seen[p]:
            c += 1
            while not seen[p]:
                seen[p] = True
                p = perm[p]
    return c

def wstr(w):
    return " ".join(str(g) for g in w)

def cyc_adjacent_inverse(w):
    return any(w[i] == -w[(i + 1) % len(w)] for i in range(len(w))) if len(w) > 1 else False

def rotate(w, k):
    k %= max(1, len(w))
    return w[k:] + w[:k]

def braid_relation(w):
    """apply sigma_a sigma_b sigma_a = sigma_b sigma_a sigma_b (|a-b|=1, equal signs) at a cyclic
    position; None if the pattern does not occur"""
    n = len(w)
    for s in range(n):
        r = rotate(w, s)
        if n >= 3 and r[0] == r[2] and abs(abs(r[0]) - abs(r[1])) == 1 and (r[0] > 0) == (r[1] > 0):
            return [r[1], r[0], r[1]] + r[3:]
    return None

def free_word(rng, n, lo, hi):
    while True:
        L = rng.randint(lo, hi)
        w = [rng.choice([1, -1]) * rng.randint(1, n - 1) for _ in range(L)]
        if set(abs(g) for g in w) == set(range(1, n)) and not cyc_adjacent_inverse(w):
            return w

def isotopy_move(rng, w, n, kinds=("rotate", "conjugate", "relation", "cancel")):
    kind = rng.choice(kinds)
    if kind == "relation":
        r = braid_relation(w)
        if r is not None:
            return "relation", r
        kind = "cancel"
    if kind == "rotate":
        return "rotate", rotate(w, rng.randint(1, max(1, len(w) - 1)))
    g = rng.choice([1, -1]) * rng.randint(1, n - 1)
    if kind == "conjugate":
        return "conjugate", [g] + w + [-g]
    pos = rng.randint(0, len(w))
    return "cancel", w[:pos] + [g, -g] + w[pos:]

# ------------------------------------------------------------------ running krhom

class Ctx:
    def __init__(self, a):
        self.krhom = a.krhom
        self.jobs = a.jobs
        self.tmp = None
        self.inject = a.inject_fault
        self.injected = False

def limit_as():
    try:
        resource.setrlimit(resource.RLIMIT_AS, (8 << 30, 8 << 30))
    except Exception:
        pass

def run_krhom(ctx, ex, ns, word, strands, colors="2", reduced=False, extra=(), timeout=600, tag=""):
    """one krhom call; returns {N: poly} or raises RunFail"""
    jpath = os.path.join(ctx.tmp, "%s-%d.json" % (ex["id"], len(ex["runs"])))
    args = [ctx.krhom, ",".join(str(x) for x in ns), "--braid", wstr(word) if word else "", "--strands", str(strands), "--json", jpath]
    if colors is not None:
        args += ["--colors", colors]
    if reduced:
        args.append("--reduced")
    args += list(extra)
    env = dict(os.environ)
    env["KRHOM_THREADS"] = str(max(1, 4 // max(1, ctx.jobs)))
    t0 = time.time()
    rec = {"tag": tag, "word": wstr(word), "args": args[1:], "exit": None, "secs": None}
    ex["runs"].append(rec)
    try:
        p = subprocess.run(args, capture_output=True, text=True, timeout=timeout, env=env, preexec_fn=limit_as)
    except subprocess.TimeoutExpired:
        rec["exit"] = "timeout"
        rec["secs"] = round(time.time() - t0, 1)
        raise RunFail("TIMEOUT", "timeout after %ds: krhom %s" % (timeout, " ".join(args[1:])))
    rec["exit"] = p.returncode
    rec["secs"] = round(time.time() - t0, 1)
    j = None
    try:
        with open(jpath) as f:
            j = json.load(f)
    except Exception:
        pass
    if j is not None:
        rec["route"] = j.get("route")
        rec["results"] = {str(r["N"]): r["P"]["str"] for r in j.get("results", [])}
        rec["primes"] = [{k: pr.get(k) for k in ("field", "seed", "secs", "chi_checked", "closures", "points", "vmhwm_kb")} for pr in j.get("primes", [])]
    if p.returncode != 0:
        tail = p.stderr.strip().splitlines()[-3:]
        status = {3: "FAIL", 2: "FAIL", 4: "SKIP(library)", 5: "ERROR"}.get(p.returncode, "ERROR")
        raise RunFail(status, "exit %d: krhom %s: %s" % (p.returncode, " ".join(args[1:]), " | ".join(tail)), code=p.returncode)
    res = {r["N"]: poly_of(r["P"]) for r in j["results"]}
    # O3 (box routes): chi checked for every closure of every prime
    if j.get("route") in ("box", "box-live"):
        for pr in j["primes"]:
            ok = pr["chi_checked"] == pr["closures"] and pr["closures"] > 0
            add_check(ex, "O3", "C", "PASS" if ok else "FAIL", "" if ok else "%s: chi_checked %s of %s closures" % (pr["field"], pr["chi_checked"], pr["closures"]))
        if ctx.inject and not ctx.injected:
            ctx.injected = True
            n0 = sorted(res)[-1]
            p0 = dict(res[n0])
            k0 = sorted(p0)[0] if p0 else (0, 0)
            p0[k0] = p0.get(k0, 0) + 1
            res[n0] = p0
            rec["injected_fault"] = "N=%d term %s +1" % (n0, k0)
    return res

class RunFail(Exception):
    def __init__(self, status, msg, code=None):
        super().__init__(msg)
        self.status = status
        self.code = code

def add_check(ex, name, label, result, detail=""):
    ex["checks"].setdefault(name, []).append({"result": result, "label": label, "detail": detail})

def cmp_check(ex, name, label, got, want, ns, what):
    bad = []
    for n in ns:
        if n in got and n in want and got[n] != want[n]:
            bad.append("N=%d %s: %s" % (n, what, pdiff(got[n], want[n])))
    if bad:
        add_check(ex, name, label, "WARN" if label == "E" else "FAIL", " | ".join(bad))
    else:
        add_check(ex, name, label, "PASS", what)

def o9(ex, res, mode):
    if 2 in res:
        ok = res[2] == {(0, 0): 1}
        add_check(ex, "O9", "C", "PASS" if ok else "FAIL", "" if ok else "N'=2 %s: %s" % (mode, pstr(res[2])))

def box_seed_args(ex):
    s2 = 202 + ex["subseed"] % 97
    return ["--no-simplify", "--seed", "201", "--seed-f32b", str(s2)]

def run_box(ctx, ex, ns, word, strands, reduced, timeout, tag):
    """box route with O2 (different seeds per prime); on exit 3 rerun with equal seeds"""
    try:
        r = run_krhom(ctx, ex, ns, word, strands, reduced=reduced, extra=box_seed_args(ex), timeout=timeout, tag=tag)
        add_check(ex, "O2", "C", "PASS", "seeds 201/%d agree" % (202 + ex["subseed"] % 97))
        add_check(ex, "O1", "C", "PASS", tag)
        return r
    except RunFail as e:
        if e.code == 3:
            try:
                run_krhom(ctx, ex, ns, word, strands, reduced=reduced, extra=["--no-simplify"], timeout=timeout, tag=tag + "-sameseed")
                add_check(ex, "O2", "C", "FAIL", "primes disagree with different seeds but agree with equal seeds: " + str(e))
            except RunFail as e2:
                add_check(ex, "O1", "C", "FAIL", "primes disagree with equal seeds too: " + str(e2))
        raise

def uncolored_n3(ctx, ex, word, strands, reduced):
    return run_krhom(ctx, ex, [3], list(reversed(word)), strands, colors=None, reduced=reduced, timeout=600, tag="O8-uncolored-rev")

# ------------------------------------------------------------------ wang.py

_wang = {}
def wang(n, m):
    key = (n, m)
    if key not in _wang:
        out = subprocess.run([sys.executable, "-I", os.path.join(HERE, "wang.py"), "2", str(n), str(m)], capture_output=True, text=True, check=True).stdout
        _wang[key] = parse_pstr(out.splitlines()[0])
    return _wang[key]

def closed_form_2strand(w, n, reduced):
    """O11: (value, label) or None"""
    e = sum(1 if g > 0 else -1 for g in w)
    a = abs(e)
    if a == 1:
        return ({(0, 0): 1} if reduced else qpoly_t0(qbinom2(n))), "C"
    if a == 0:
        b = qbinom2(n)
        return (qpoly_t0(b) if reduced else qmul(b, qpoly_t0(b))), "C"
    if a % 2 == 1 and not reduced:
        m = (a - 1) // 2
        v = wang(n, m)
        if e < 0:
            v = mirror_poly(v)
        return v, ("B" if m >= 4 or (m == 3 and n >= 5) else "C")
    return None

def unknot_value(n, reduced):
    return {(0, 0): 1} if reduced else qpoly_t0(qbinom2(n))

# ------------------------------------------------------------------ tiers

def ex_t2(ctx, ex):
    rng = random.Random(ex["subseed"])
    L = rng.randint(1, 9)
    w = [rng.choice([1, -1]) for _ in range(L)]
    ex.update(word=wstr(w), lam=[2, 2], components=ncomp(w, 2), mode="T2")
    ns = [2, 3, 4, 5, 6]
    knot = ex["components"] == 1
    base = {}
    for red in (False, True):
        mode = "r" if red else "u"
        r = run_krhom(ctx, ex, ns, w, 2, reduced=red, timeout=120, tag="base-" + mode)
        add_check(ex, "O1", "C", "PASS", "base " + mode)
        base[red] = r
        o9(ex, r, mode)
        cf_ = [(n, closed_form_2strand(w, n, red)) for n in ns]
        for n, c in cf_:
            if c is not None:
                v, lab = c
                cmp_check(ex, "O11", lab, {n: r[n]}, {n: v}, [n], "closed form %s n=%d" % (mode, sum(1 if g > 0 else -1 for g in w)))
        # O6 mirror
        mw = [-g for g in w]
        rm = run_krhom(ctx, ex, ns, mw, 2, reduced=red, timeout=120, tag="mirror-" + mode)
        cmp_check(ex, "O6", "C", {n: mirror_poly(p) for n, p in rm.items()}, r, ns, "mirror " + mode)
        # O4 rotation (2 strands: one component through position 0 for knots; links: unreduced only)
        if L >= 2 and (knot or not red):
            rw = rotate(w, rng.randint(1, L - 1))
            rr = run_krhom(ctx, ex, ns, rw, 2, reduced=red, timeout=120, tag="rotate-" + mode)
            cmp_check(ex, "O4", "C", rr, r, ns, "rotation " + mode)
        # O8 N'=3 duality
        u3 = uncolored_n3(ctx, ex, w, 2, red)
        cmp_check(ex, "O8", "B", {3: r[3]}, u3, [3], "N'=3 duality " + mode)
        # O13a cross-engine
        lv = run_krhom(ctx, ex, [2, 3, 4], w, 2, reduced=red, extra=["--colored-method", "singular", "--box-live", "--no-simplify"], timeout=300, tag="box-live-" + mode)
        cmp_check(ex, "O13a", "C", lv, r, [2, 3, 4], "box-live " + mode)
    if knot:
        for n in (3, 4):
            d = lewark(n, base[True][n], base[False][n])
            add_check(ex, "O10b", "E", "PASS" if d is None else "WARN", "N=%d %s" % (n, d or "ok"))

def ex_t3(ctx, ex, mode_a, unknot_core):
    rng = random.Random(ex["subseed"])
    ns = [2, 3, 4]
    if mode_a:
        if unknot_core:
            core = [rng.choice([1, -1])]
        else:
            core = [rng.choice([1, -1]) for _ in range(rng.randint(1, 4))]
        pos = rng.randint(0, len(core))
        w = core[:pos] + [rng.choice([2, -2])] + core[pos:]
        moves = []
        for _ in range(rng.randint(1, 2)):
            k, w2 = isotopy_move(rng, w, 3)
            if len(w2) <= 7:
                moves.append(k)
                w = w2
        ex.update(word=wstr(w), lam=[2, 2, 2], components=ncomp(w, 3), mode="A", core={"word": wstr(core)}, partner={"kind": "markov+" + "+".join(moves), "word": wstr(w)})
        knot = ex["components"] == 1
        for red in ((False, True) if knot else (False,)):
            m = "r" if red else "u"
            cv = run_krhom(ctx, ex, ns, core, 2, reduced=red, timeout=120, tag="core-" + m)
            ex["core"]["value_" + m] = {str(n): pstr(p) for n, p in cv.items()}
            c = [(n, closed_form_2strand(core, n, red)) for n in ns]
            for n, cc in c:
                if cc is not None:
                    cmp_check(ex, "O11", cc[1], {n: cv[n]}, {n: cc[0]}, [n], "core closed form " + m)
            r = run_box(ctx, ex, ns, w, 3, red, 2400, "box-" + m)
            o9(ex, r, m)
            cmp_check(ex, "O5", "C", r, cv, ns, "Markov/isotopy vs core " + m)
            if unknot_core:
                cmp_check(ex, "O12", "C", r, {n: unknot_value(n, red) for n in ns}, ns, "unknot " + m)
            u3 = uncolored_n3(ctx, ex, w, 3, red)
            cmp_check(ex, "O8", "B", {3: r[3]}, u3, [3], "N'=3 duality " + m)
            ex.setdefault("_res", {})[red] = r
        if knot:
            for n in (3, 4):
                d = lewark(n, ex["_res"][True][n], ex["_res"][False][n])
                add_check(ex, "O10b", "E", "PASS" if d is None else "WARN", "N=%d %s" % (n, d or "ok"))
    else:
        w = free_word(rng, 3, 4, 6)
        kind = rng.choice(["isotopy", "mirror", "reverse"])
        if kind == "isotopy":
            k, pw = isotopy_move(rng, w, 3)
            if len(pw) > 6:
                k, pw = "rotate", rotate(w, rng.randint(1, len(w) - 1))
        elif kind == "mirror":
            k, pw = "mirror", [-g for g in w]
        else:
            k, pw = "reverse", list(reversed(w))
        ex.update(word=wstr(w), lam=[2, 2, 2], components=ncomp(w, 3), mode="B", partner={"kind": k, "word": wstr(pw)})
        knot = ex["components"] == 1
        res = {}
        for red in ((False, True) if knot else (False,)):
            m = "r" if red else "u"
            r = run_box(ctx, ex, ns, w, 3, red, 2400, "box-" + m)
            res[red] = r
            o9(ex, r, m)
            u3 = uncolored_n3(ctx, ex, w, 3, red)
            cmp_check(ex, "O8", "B", {3: r[3]}, u3, [3], "N'=3 duality " + m)
        pr = run_box(ctx, ex, ns, pw, 3, False, 2400, "partner-" + k)
        o9(ex, pr, "u")
        if k == "mirror":
            cmp_check(ex, "O6", "C", {n: mirror_poly(p) for n, p in pr.items()}, res[False], ns, "mirror u")
        elif k == "reverse":
            cmp_check(ex, "O7", "E", pr, res[False], [3, 4], "reversal u")
        else:
            cmp_check(ex, "O4", "C", pr, res[False], ns, k + " u")
        if knot:
            for n in (3, 4):
                d = lewark(n, res[True][n], res[False][n])
                add_check(ex, "O10b", "E", "PASS" if d is None else "WARN", "N=%d %s" % (n, d or "ok"))

def t4_candidates(rng, k=50):
    out = []
    for _ in range(k):
        e = rng.choice([1, -1, 3, -3])
        core = [1 if e > 0 else -1] * abs(e)
        p2 = rng.randint(0, len(core))
        w3 = core[:p2] + [rng.choice([2, -2])] + core[p2:]
        p3 = rng.randint(0, len(w3))
        w4 = w3[:p3] + [rng.choice([3, -3])] + w3[p3:]
        rot = rng.random() < 0.5
        if rot:
            w4 = rotate(w4, rng.randint(1, len(w4) - 1))
        out.append((core, w4, rot))
    return out

def ex_t4(ctx, ex, covered):
    rng = random.Random(ex["subseed"])
    cands = t4_candidates(rng)
    pick = [c for c in cands if wstr(c[1]) in covered]
    if not pick:
        ex.update(word="", lam=[2, 2, 2, 2], components=1, mode="A")
        ex["status"] = "SKIP(library)"
        return
    core, w, rot = pick[0]
    ex.update(word=wstr(w), lam=[2, 2, 2, 2], components=ncomp(w, 4), mode="A", core={"word": wstr(core)}, partner={"kind": "markov2" + ("+rotate" if rot else ""), "word": wstr(w)})
    ns = [2, 3, 4]
    modes = (True, False) if len(w) <= 3 else (True,)
    for red in modes:
        m = "r" if red else "u"
        cv = run_krhom(ctx, ex, ns, core, 2, reduced=red, timeout=120, tag="core-" + m)
        r = run_box(ctx, ex, ns, w, 4, red, 2400, "box-" + m)
        o9(ex, r, m)
        cmp_check(ex, "O5", "C", r, cv, ns, "Markov vs core " + m)
        if abs(sum(1 if g > 0 else -1 for g in core)) == 1:
            cmp_check(ex, "O12", "C", r, {n: unknot_value(n, red) for n in ns}, ns, "unknot " + m)
        u3 = uncolored_n3(ctx, ex, w, 4, red)
        cmp_check(ex, "O8", "B", {3: r[3]}, u3, [3], "N'=3 duality " + m)

# ------------------------------------------------------------------ harness

def finish_status(ex):
    if ex.get("status"):
        return ex["status"]
    res = [c["result"] for cs in ex["checks"].values() for c in cs]
    if "FAIL" in res:
        return "FAIL"
    if "WARN" in res:
        return "WARN"
    return "PASS"

def run_example(ctx, ex, covered):
    t0 = time.time()
    try:
        if ex["tier"] == 2:
            ex_t2(ctx, ex)
        elif ex["tier"] == 3:
            ex_t3(ctx, ex, ex["plan"] in ("A", "A-unknot"), ex["plan"] == "A-unknot")
        else:
            ex_t4(ctx, ex, covered)
        ex["status"] = finish_status(ex)
    except RunFail as e:
        ex["error"] = str(e)
        ex["status"] = e.status if not any(c["result"] == "FAIL" for cs in ex["checks"].values() for c in cs) else "FAIL"
    except Exception as e:  # harness bug: ERROR, never PASS
        ex["error"] = "harness: %r" % e
        ex["status"] = "ERROR"
    ex.pop("_res", None)
    ex["secs"] = round(time.time() - t0, 1)
    return ex

def self_tests(ctx, tiers, anchor):
    errs = []
    # 1. convention: exact path = wang.py (unreduced) for (2,2) sigma^3, N = 2,3,4
    ex = {"id": "selftest-convention", "runs": [], "checks": {}, "subseed": 0}
    r = run_krhom(ctx, ex, [2, 3, 4], [1, 1, 1], 2, timeout=120)
    for n in (2, 3, 4):
        if r[n] != wang(n, 1):
            errs.append("convention: krhom (2,2) sigma^3 N=%d %s != wang.py %s" % (n, pstr(r[n]), pstr(wang(n, 1))))
    # 3. libraries installed
    for t, lib in ((3, "l222_n4"), (4, "l2222_n4")):
        if t in tiers:
            p = subprocess.run(["bash", os.path.join(HERE, "fetch_library.sh"), "--verify-only", lib], capture_output=True, text=True)
            if p.returncode != 0:
                errs.append("tier %d needs library %s: run tools/fetch_library.sh %s (%s)" % (t, lib, lib, p.stderr.strip().splitlines()[-1:] if p.stderr else ""))
    # 2. anchor: 4_1 on 3 strands, both modes
    if anchor and 3 in tiers and not errs:
        want = {False: [1, 11, 46], True: [1, 5, 25]}
        n3 = {False: "t^-2q^6 + t^-2q^8 + t^-1q^0 + t^-1q^2 + t^0q^-2 + t^0q^0 + t^0q^2 + t^1q^-2 + t^1q^0 + t^2q^-8 + t^2q^-6", True: "t^-2q^6 + t^-1q^2 + t^0q^0 + t^1q^-2 + t^2q^-6"}
        ref = os.path.join(ROOT, "runs_c2/k3/out/4_1_f31.json")
        refj = json.load(open(ref)) if os.path.exists(ref) else None
        for red in (False, True):
            ex = {"id": "selftest-anchor-%s" % ("r" if red else "u"), "runs": [], "checks": {}, "subseed": 0}
            try:
                r = run_krhom(ctx, ex, [2, 3, 4], [1, -2, 1, -2], 3, reduced=red, extra=["--no-simplify"], timeout=2400)
            except RunFail as e:
                errs.append("anchor 4_1: %s" % e)
                continue
            if [sum(p.values()) for _, p in sorted(r.items())] != want[red] or pstr(r[3]) != n3[red]:
                errs.append("anchor 4_1 %s: %s" % ("reduced" if red else "unreduced", {n: pstr(p) for n, p in r.items()}))
            if refj:
                for x in refj["results"]:
                    if x["N"] == 4 and x["reduced"] == red and x["P"]["str"] != pstr(r[4]):
                        errs.append("anchor 4_1 N=4 %s differs from runs_c2/k3/out/4_1_f31.json" % ("reduced" if red else "unreduced"))
    return errs

def t4_coverage(ctx, cand_words):
    """words (strings) with 0 missing jobs in both primes' status files (m5i_cover example)"""
    exe = os.path.join(ROOT, "target/release/examples/m5i_cover")
    if not os.path.exists(exe):
        subprocess.run(["cargo", "build", "--release", "--example", "m5i_cover"], cwd=ROOT, check=True, capture_output=True)
    st = [open(os.path.join(ROOT, "data/libstatus/l2222_n4_%s.txt" % f)).read().splitlines() for f in ("f31", "f32b")]
    ok = set(l for l in st[0] if l.endswith(" ok")) & set(st[1])
    both = os.path.join(ctx.tmp, "both.txt")
    cands = os.path.join(ctx.tmp, "cands.txt")
    with open(both, "w") as f:
        f.write("\n".join(sorted(ok)) + "\n")
    with open(cands, "w") as f:
        for i, w in enumerate(cand_words):
            f.write("c%d;2,2,2,2;%s;4\n" % (i, w))
    out = subprocess.run([exe, "2,2,2,2", "4", cands, both], capture_output=True, text=True, check=True).stdout
    covered = set()
    for line in out.splitlines():
        m = re.match(r"c(\d+) \[.*\]: .*missing (\d+) ", line)
        if m and m.group(2) == "0":
            covered.add(cand_words[int(m.group(1))])
    return covered

CHECK_NAMES = [("O1", "primes"), ("O2", "seeds"), ("O3", "chi"), ("O4", "isotopy"), ("O5", "markov"), ("O6", "mirror"), ("O8", "duality-N3"), ("O7", "reversal"), ("O9", "N2"), ("O10b", "lewark"), ("O12", "unknot"), ("O11", "closed-form"), ("O13a", "cross-engine")]

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--n", type=int, default=24)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--full", action="store_true")
    ap.add_argument("--tiers", default="2,3,4")
    ap.add_argument("--only")
    ap.add_argument("--krhom", default=os.path.join(ROOT, "target/release/krhom"))
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--out", default=os.path.join(ROOT, "runs_random"))
    ap.add_argument("--inject-fault", action="store_true")
    ap.add_argument("--no-anchor", action="store_true", help="skip the 4_1 anchor self-test (about 3 min)")
    a = ap.parse_args()
    if not os.path.exists(a.krhom):
        print("randcheck: %s not found (cargo build --release --bin krhom)" % a.krhom, file=sys.stderr)
        return 2
    tiers = [int(x) for x in a.tiers.split(",") if x]
    if a.quick:
        sizes = {2: 4, 3: 6, 4: 0}
    elif a.full:
        sizes = {2: 10, 3: 20, 4: 10}
    else:
        n = a.n
        sizes = {2: n // 3, 3: n // 2, 4: n - n // 3 - n // 2}
    rng = random.Random(a.seed)
    exs = []
    for t in (2, 3, 4):
        for i in range(sizes[t]):
            ex = {"id": "T%d-%02d" % (t, i + 1), "subseed": rng.getrandbits(32), "tier": t, "runs": [], "checks": {}}
            if t == 3:
                nA = sizes[3] // 2
                ex["plan"] = ("A-unknot" if i < 2 else "A") if i < nA else "B"
            exs.append(ex)
    exs = [e for e in exs if e["tier"] in tiers]
    if a.only:
        exs = [e for e in exs if e["id"] == a.only]
        if not exs:
            print("randcheck: no example %s" % a.only, file=sys.stderr)
            return 2
    outdir = os.path.join(a.out, str(a.seed) + ("-fault" if a.inject_fault else "") + ("-" + a.only if a.only else ""))
    os.makedirs(outdir, exist_ok=True)
    ctx = Ctx(a)
    ctx.tmp = tempfile.mkdtemp(prefix="randcheck-", dir=outdir)
    t0 = time.time()
    ctx.inject = False
    errs = self_tests(ctx, tiers, not (a.no_anchor or a.only))
    ctx.inject = a.inject_fault
    if errs:
        for e in errs:
            print("randcheck: self-test FAILED: " + e, file=sys.stderr)
        return 2
    print("randcheck: self-tests passed (convention%s, libraries) in %.0fs" % ("" if (a.no_anchor or a.only or 3 not in tiers) else ", 4_1 anchor", time.time() - t0), file=sys.stderr)
    covered = set()
    t4 = [e for e in exs if e["tier"] == 4]
    if t4:
        words = []
        for e in t4:
            words += [wstr(c[1]) for c in t4_candidates(random.Random(e["subseed"]))]
        covered = t4_coverage(ctx, sorted(set(words)))
        print("randcheck: T4 coverage filter: %d of %d candidate words covered" % (len(covered), len(set(words))), file=sys.stderr)
    if a.inject_fault:
        ctx.jobs = 1
    # T4 and T3 first (longest), T2 last
    order = sorted(exs, key=lambda e: -e["tier"])
    with cf.ThreadPoolExecutor(max_workers=1 if a.inject_fault else a.jobs) as pool:
        futs = {pool.submit(run_example, ctx, e, covered): e for e in order}
        for f in cf.as_completed(futs):
            e = f.result()
            print("%s %-14s %-6s %s  %s (%.0fs)%s" % (e["id"], e["status"], e.get("mode", ""), e.get("word", ""), (e.get("partner") or {}).get("kind", ""), e["secs"], ("  " + e["error"]) if e.get("error") else ""), file=sys.stderr, flush=True)
    exs.sort(key=lambda e: e["id"])
    with open(os.path.join(outdir, "examples.jsonl"), "w") as f:
        for e in exs:
            e2 = {k: v for k, v in e.items() if k != "plan"}
            e2["checks"] = {k: {"result": ("FAIL" if any(c["result"] == "FAIL" for c in cs) else "WARN" if any(c["result"] == "WARN" for c in cs) else "PASS"), "label": cs[0]["label"], "n": len(cs), "detail": " | ".join(c["detail"] for c in cs if c["result"] != "PASS")} for k, cs in e["checks"].items()}
            f.write(json.dumps(e2) + "\n")
    cnt = {s: 0 for s in ("PASS", "FAIL", "WARN", "SKIP", "TIMEOUT", "ERROR")}
    for e in exs:
        cnt["SKIP" if e["status"].startswith("SKIP") else e["status"]] += 1
    chk = []
    for k, nm in CHECK_NAMES:
        p = sum(1 for e in exs for c in e["checks"].get(k, []) if c["result"] == "PASS")
        chk.append("%s %d" % (nm, p))
    final = "seed=%d n=%d PASS %d FAIL %d WARN %d SKIP %d TIMEOUT %d ERROR %d; checks: %s; %.0fs" % (a.seed, len(exs), cnt["PASS"], cnt["FAIL"], cnt["WARN"], cnt["SKIP"], cnt["TIMEOUT"], cnt["ERROR"], ", ".join(chk), time.time() - t0)
    with open(os.path.join(outdir, "summary.md"), "w") as f:
        f.write("# randcheck seed %d%s\n\nkrhom: `%s`; tiers %s; jobs %d\n\n" % (a.seed, " (--inject-fault)" if a.inject_fault else "", a.krhom, a.tiers, a.jobs))
        f.write("| id | tier | word | partner | comps | N'=4 ranks u/r | checks | status | secs |\n|---|---|---|---|---|---|---|---|---|\n")
        for e in exs:
            r4 = {}
            for r in e["runs"]:
                if r.get("tag", "").startswith(("base-", "box-")) and "results" in r and "4" in r["results"]:
                    r4[r["tag"][-1]] = sum(v for v in poly_of_str_ranks(r["results"]["4"]))
            checks = ", ".join("%s:%s" % (k, ("FAIL" if any(c["result"] == "FAIL" for c in cs) else "WARN" if any(c["result"] == "WARN" for c in cs) else "PASS")) for k, cs in sorted(e["checks"].items()))
            f.write("| %s | %d | `%s` | %s | %s | %s/%s | %s | %s | %s |\n" % (e["id"], e["tier"], e.get("word", ""), (e.get("partner") or {}).get("kind", ""), e.get("components", ""), r4.get("u", "-"), r4.get("r", "-"), checks, e["status"], e["secs"]))
        bad = [e for e in exs if e["status"] in ("FAIL", "ERROR", "WARN", "TIMEOUT")]
        if bad:
            f.write("\n## Details\n\n")
            for e in bad:
                f.write("- %s %s: %s\n" % (e["id"], e["status"], e.get("error", "")))
                for k, cs in e["checks"].items():
                    for c in cs:
                        if c["result"] != "PASS":
                            f.write("  - %s (%s) %s: %s\n" % (k, c["label"], c["result"], c["detail"]))
                for r in e["runs"]:
                    f.write("  - repro: `krhom %s`\n" % " ".join(("'%s'" % x if " " in x else x) for x in r["args"] if not x.endswith(".json") and x != "--json"))
        f.write("\n%s\n" % final)
    print(final)
    return 1 if cnt["FAIL"] or cnt["ERROR"] else 0

def poly_of_str_ranks(s):
    return parse_pstr(s).values()

if __name__ == "__main__":
    sys.exit(main())
