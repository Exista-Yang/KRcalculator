"""Build normalised knot/link reference tables from the KnotInfo / LinkInfo downloads.

usage: build.py <work_dir> <knottheory_dir> <out_dir>
  work_dir/knotinfo_sel.json, work_dir/linkinfo_sel.json : extracted by extract.py
  knottheory_dir : unpacked KnotTheory` package (Knot Atlas data files)
"""
import sys, os, json, re, time, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction
from krlib import (parse_pd, ints, jones_pd, signature_pd, linking_matrix, braid_to_pd,
                   braid_components, parse_poly, split_top, orient, mirror1, padd, pmul,
                   homfly_pd, mirror_homfly)

work, ktdir, outdir = sys.argv[1:4]
MAXC_KNOT = 12
BUDGET = 150000  # max memo size for HOMFLY of braid closures (long braids are skipped)
t0 = time.time()

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def intexp(f):
    assert f.denominator == 1, f
    return int(f)


def poly1(s, var, scale):
    """single-variable polynomial string -> {scale*exp: c} (exponents must become ints)"""
    p = parse_poly(s, [var])
    out = {}
    for (e,), c in p.t.items():
        out[intexp(e * scale)] = c
    return out


def poly2(s, vars_):
    p = parse_poly(s, vars_)
    return {tuple(intexp(e) for e in k): c for k, c in p.t.items()}


def braces_json(s):
    s = s.strip().replace('{', '[').replace('}', ']').replace(')', ']')
    return json.loads(s)


def terms1(d):
    return [[k, d[k]] for k in sorted(d)]


def relation(computed, listed):
    same = computed == listed
    mir = mirror1(computed) == listed
    if same and mir:
        return 'symmetric'
    if same:
        return 'same'
    if mir:
        return 'mirror'
    return 'mismatch'


def hrelation(computed, listed):
    if computed is None:
        return 'unknown'
    same = computed == listed
    mir = mirror_homfly(computed) == listed
    if same and mir:
        return 'symmetric'
    return 'same' if same else ('mirror' if mir else 'mismatch')


def combine(jrel, hrel, sig_c, sig_l, amphichiral=False, extra=()):
    """decide diagram vs listed invariants from Jones, HOMFLY, signature (+ linking numbers)"""
    cand = {'same', 'mirror'}
    for r in (jrel, hrel) + tuple(extra):
        if r == 'mismatch':
            return 'conflict'
        if r not in ('symmetric', 'unknown'):
            cand &= {r}
    if sig_c != sig_l:
        cand.discard('same')
    if sig_c != -sig_l:
        cand.discard('mirror')
    if cand == {'same'}:
        return 'same'
    if cand == {'mirror'}:
        return 'mirror'
    if not cand:
        return 'conflict'
    return 'same' if amphichiral else 'same_or_mirror'


def neg_matrix(m):
    return [[-v for v in row] for row in m]


def lk_relation(computed, listed, permute=False):
    import itertools
    n = len(listed)
    if len(computed) != n:
        return 'mismatch'
    perms = itertools.permutations(range(n)) if permute else [tuple(range(n))]
    same = mir = False
    for p in perms:
        pm = [[computed[p[i]][p[j]] for j in range(n)] for i in range(n)]
        same = same or pm == listed
        mir = mir or neg_matrix(pm) == listed
    if same and mir:
        return 'symmetric'
    return 'same' if same else ('mirror' if mir else 'mismatch')


def kh_mirror_symmetric(free, tors=()):
    return (sorted([-h, -q, r] for h, q, r in free) == sorted(free) and
            sorted([1 - h, -q, m, r] for h, q, m, r in tors) == sorted(tors))


def x_minus_xinv_pow(b):
    r = {0: 1}
    for _ in range(b):
        r = pmul(r, {1: 1, -1: -1})
    return r


def homfly_to_jones_check(hom, jones):
    """P(v=x^2, z=x-1/x) == V (x = t^(1/2)), compared after multiplying by (x-1/x)^B"""
    B = max(0, -min(b for (a, b) in hom))
    lhs = {}
    for (a, b), c in hom.items():
        lhs = padd(lhs, pmul({2 * a: c}, x_minus_xinv_pow(b + B)))
    rhs = pmul(jones, x_minus_xinv_pow(B))
    return lhs == rhs


def kh_euler(free):
    r = {}
    for h, q, rk in free:
        r = padd(r, {q: (-1) ** (h % 2) * rk})
    return r


def jones_at_minus_q(jones):
    return {k: c * (-1) ** (k % 2) for k, c in jones.items()}


def kh_vs_jones(free, jones, reduced):
    lhs = kh_euler(free)
    rhs = jones_at_minus_q(jones)
    if not reduced:
        rhs = pmul(rhs, {1: 1, -1: 1})
    return lhs == rhs


def s_vs_kh(free, s):
    have = {(h, q) for h, q, rk in free if rk > 0}
    return (0, s - 1) in have and (0, s + 1) in have


def kh_from_knotinfo_vector(vec):
    free = collections.Counter()
    tors = collections.Counter()
    for mod, rank, h, q in vec:
        if mod in (0, 1):
            free[(h, q)] += rank
        else:
            tors[(h, q, mod)] += rank
    return ([[h, q, r] for (h, q), r in sorted(free.items())],
            [[h, q, m, r] for (h, q, m), r in sorted(tors.items())])


def kh_from_knotinfo_string(s):
    p = poly2(s, ['t', 'q', 'T'])
    free = collections.Counter()
    tors = collections.Counter()
    for (h, q, T), c in p.items():
        if T == 0:
            free[(h, q)] += c
        else:
            tors[(h, q, T)] += c
    return ([[h, q, r] for (h, q), r in sorted(free.items())],
            [[h, q, m, r] for (h, q, m), r in sorted(tors.items())])


# ---------------------------------------------------------------------------
# Knot Atlas (KnotTheory` package data files)
# ---------------------------------------------------------------------------

def mlist(path, var):
    text = open(path).read()
    i = text.index(var + ' = {')
    j = text.index('\n\n', i)
    body = text[i + len(var) + 4: j].strip()
    assert body.endswith('}')
    return [s.strip() for s in split_top(body[:-1])]


def pd_file(path, kind):
    text = open(path).read()
    out = {}
    for chunk in re.split(r'\nPD\[' + kind + r'\[', text)[1:]:
        head, body = chunk.split(']]=', 1)
        body = body.split('\nEnd[]')[0]
        out[head.replace(' ', '')] = parse_pd(body) if 'Loop' not in body else []
    return out


ka_knots_order = []
for n, cnt in [(0, 1), (3, 1), (4, 1), (5, 2), (6, 3), (7, 7), (8, 21), (9, 49), (10, 165)]:
    ka_knots_order += [f"{n}_{k}" for k in range(1, cnt + 1)]
ka_jones = {}
for name, s in zip(ka_knots_order, mlist(os.path.join(ktdir, 'Jones4Knots.m'), 'Js')):
    ka_jones[name] = poly1(s, 'q', 2)
assert len(ka_jones) == 250
k11 = [f"11a_{k}" for k in range(1, 368)] + [f"11n_{k}" for k in range(1, 186)]
for name, s in zip(k11, mlist(os.path.join(ktdir, 'Jones4Knots11.m'), 'Js')):
    ka_jones[name] = poly1(s, 'q', 2)
assert len(ka_jones) == 250 + 552
ka_kh11 = {}
for name, s in zip(k11, mlist(os.path.join(ktdir, 'Kh4Knots11.m'), 'Is')):
    ka_kh11[name] = sorted([q, h, c] for (q, h), c in poly2(s, ['q', 't']).items())
ka_pd_knots = pd_file(os.path.join(ktdir, 'PD4Knots.m'), 'Knot')
ka_pd_k = {}
for head, pd in ka_pd_knots.items():
    n, k = head.split(',')
    ka_pd_k[f"{n}_{k}"] = pd
# validate: Jones of Knot Atlas PD == Knot Atlas Jones (checks ordering and our Jones code)
ka_selfcheck_fail = [nm for nm, pd in ka_pd_k.items() if pd and jones_pd(pd) != ka_jones[nm]]
print('Knot Atlas knots PD-vs-Jones self-check failures:', ka_selfcheck_fail, flush=True)

ka_pd_links = {}
for head, pd in pd_file(os.path.join(ktdir, 'PD4Links.m'), 'Link').items():
    n, t, k = head.split(',')
    ka_pd_links[f"L{n}{'a' if t == 'Alternating' else 'n'}{k}"] = pd
link_order = list(ka_pd_links.keys())
ka_link_jones = {}
for name, s in zip(link_order, mlist(os.path.join(ktdir, 'Jones4Links.m'), 'Js')):
    ka_link_jones[name] = poly1(s, 'q', 2)
assert len(ka_link_jones) == len(ka_pd_links) == 1424, (len(ka_link_jones), len(ka_pd_links))
ka_link_selfcheck_fail = [nm for nm, pd in ka_pd_links.items() if jones_pd(pd) != ka_link_jones[nm]]
print('Knot Atlas links PD-vs-Jones self-check failures:', ka_link_selfcheck_fail[:20], len(ka_link_selfcheck_fail), flush=True)
print('knot atlas loaded', round(time.time() - t0, 1), flush=True)

# ---------------------------------------------------------------------------
# knots
# ---------------------------------------------------------------------------
stats = collections.defaultdict(collections.Counter)
problems = []
knots_out = []
for d in json.load(open(os.path.join(work, 'knotinfo_sel.json'))):
    cn = int(d['crossing_number'])
    if cn == 0 or cn > MAXC_KNOT:
        continue
    name = d['name']
    pd = parse_pd(d['pd_notation'])
    assert len(pd) == cn, name
    signs, comp, _ = orient(pd)
    assert len(set(comp.values())) == 1, name
    bn = json.loads(d['braid_notation'].replace(' ', ''))
    words = bn if isinstance(bn[0], list) else [bn]
    bidx = int(d['braid_index'])
    # several braids are listed for a few knots: [minimal length, minimal strands];
    # keep the one with the fewest strands as 'braid', the other as 'braid_alt'
    words = sorted(words, key=lambda w: (max(abs(g) for g in w) + 1, len(w)))
    word = words[0]
    strands = max(abs(g) for g in word) + 1
    alt = None
    if len(words) > 1:
        assert len(words) == 2, name
        alt = {'strands': max(abs(g) for g in words[1]) + 1, 'word': words[1]}
    jones = poly1(d['jones_polynomial'], 't', 2)
    jv = ints(d['jones_polynomial_vector'])
    jvd = {2 * (jv[0] + i): c for i, c in enumerate(jv[2:]) if c}
    assert jv[1] - jv[0] + 1 == len(jv) - 2
    issues = []
    if jvd != jones:
        issues.append('jones_string_vs_vector')
    hom = poly2(d['homfly_polynomial'], ['v', 'z'])
    hv = json.loads(d['homfly_polynomial_vector'])
    hvd = {}
    for i, inner in enumerate(hv[2:]):
        for j, c in enumerate(inner[2:]):
            if c:
                hvd[(2 * (inner[0] + j), 2 * (hv[0] + i))] = c
    if hvd != hom:
        issues.append('homfly_string_vs_vector')
    kh, kht = kh_from_knotinfo_vector(json.loads(d['khovanov_unreduced_integral_vector']))
    if (kh, kht) != kh_from_knotinfo_string(d['khovanov_unreduced_integral_polynomial']):
        issues.append('kh_string_vs_vector')
    khr, khrt = kh_from_knotinfo_vector(json.loads(d['khovanov_reduced_integral_vector']))
    if (khr, khrt) != kh_from_knotinfo_string(d['khovanov_reduced_integral_polynomial']):
        issues.append('khred_string_vs_vector')
    khrq, _ = kh_from_knotinfo_vector(json.loads(d['khovanov_reduced_rational_vector']))
    if khrq != khr:
        issues.append('khred_Q_vs_Z_free_part')
    sig = int(d['signature'])
    s_inv = int(d['rasmussen_invariant']) if d['rasmussen_invariant'] not in ('', 'Not known') else None
    # internal consistency of the listed invariants
    if not homfly_to_jones_check(hom, jones):
        issues.append('homfly_vs_jones')
    if not kh_vs_jones(kh, jones, False):
        issues.append('kh_vs_jones')
    if not kh_vs_jones(khr, jones, True):
        issues.append('khred_vs_jones')
    if s_inv is not None and not s_vs_kh(kh, s_inv):
        issues.append('s_vs_kh')
    det = int(d['determinant'])
    vm1 = sum(c * (-1) ** (k // 2) for k, c in jones.items())  # V(-1) for knots
    if abs(vm1) != det:
        issues.append('determinant_vs_jones')
    # PD and braid against listed invariants
    jp = jones_pd(pd)
    sp = signature_pd(pd, 0)
    if signature_pd(pd, 1) != sp:
        issues.append('signature_colouring_disagree')
    amph = 'amphicheiral' in d['symmetry_type']
    pd_rel = combine(relation(jp, jones), hrelation(homfly_pd(pd), hom), sp, sig, amph)
    bpd = braid_to_pd(strands, word)
    if braid_components(strands, word) != 1:
        issues.append('braid_not_a_knot')
    jb = jones_pd(bpd)
    sb = signature_pd(bpd, 0)
    br_rel = combine(relation(jb, jones), hrelation(homfly_pd(bpd, BUDGET), hom), sb, sig, amph)
    alt_rel = None
    if alt is not None:
        apd = braid_to_pd(alt['strands'], alt['word'])
        if braid_components(alt['strands'], alt['word']) != 1:
            issues.append('braid_alt_not_a_knot')
        alt_rel = combine(relation(jones_pd(apd), jones), hrelation(homfly_pd(apd, BUDGET), hom),
                          signature_pd(apd, 0), sig, amph)
        if alt_rel != 'same':
            issues.append('braid_alt_' + alt_rel)
    for lab, rel in (('pd', pd_rel), ('braid', br_rel)):
        if rel == 'same_or_mirror':
            sym = kh_mirror_symmetric(kh, kht) and kh_mirror_symmetric(khr, khrt) and s_inv == 0
            issues.append(f"{lab}_chirality_undetermined" + ('' if sym else '_and_kh_or_s_not_mirror_symmetric'))
        elif rel != 'same':
            issues.append(f"{lab}_{rel}")
    if strands != bidx:
        issues.append('braid_strands_ne_braid_index')
    # Knot Atlas relation
    ka = None
    if name in ka_jones:
        r = relation(ka_jones[name], jones)
        if r == 'symmetric' and d['symmetry_type'].strip() in ('fully amphicheiral', 'negative amphicheiral', 'positive amphicheiral'):
            ka = 'amphichiral'
        elif r == 'symmetric':
            if name in ka_pd_k:
                s_ka = signature_pd(ka_pd_k[name], 0)
                ka = 'same' if s_ka == sig and sig != 0 else ('mirror' if s_ka == -sig and sig != 0 else 'undetermined')
            elif name in ka_kh11:
                kif = sorted([q, h, c] for h, q, c in kh)
                kim = sorted([-q, -h, c] for h, q, c in kh)
                ka = 'same' if ka_kh11[name] == kif and kif != kim else (
                    'mirror' if ka_kh11[name] == kim and kif != kim else 'undetermined')
            else:
                ka = 'undetermined'
        elif r in ('same', 'mirror'):
            ka = r
        else:
            ka = 'mismatch'
            issues.append('knot_atlas_jones_mismatch')
        stats['knot_atlas'][ka] += 1
    rec = collections.OrderedDict()
    rec['name'] = name
    rec['crossings'] = cn
    rec['components'] = 1
    rec['alternating'] = d['alternating'] == 'Y'
    rec['symmetry'] = d['symmetry_type'].strip()
    rec['pd'] = pd
    rec['writhe'] = sum(signs)
    rec['braid'] = {'strands': strands, 'word': word}
    if alt is not None:
        rec['braid_alt'] = alt
    rec['braid_index'] = bidx
    rec['jones'] = terms1(jones)
    rec['homfly'] = [[a, b, hom[(a, b)]] for (a, b) in sorted(hom, key=lambda k: (k[1], k[0]))]
    rec['kh'] = kh
    rec['kh_torsion'] = kht
    rec['kh_reduced'] = khr
    rec['kh_reduced_torsion'] = khrt
    rec['signature'] = sig
    rec['s'] = s_inv
    rec['determinant'] = det
    rec['knot_atlas'] = ka
    rec['pd_vs_invariants'] = pd_rel
    rec['braid_vs_invariants'] = br_rel
    rec['issues'] = issues
    knots_out.append(rec)
    stats['knot_pd'][pd_rel] += 1
    stats['knot_braid'][br_rel] += 1
    for i in issues:
        stats['knot_issue'][i] += 1
    if issues:
        problems.append((name, issues))
print('knots done', len(knots_out), round(time.time() - t0, 1), flush=True)

# ---------------------------------------------------------------------------
# links
# ---------------------------------------------------------------------------
links_out = []
braid_conflicts = []
for d in json.load(open(os.path.join(work, 'linkinfo_sel.json'))):
    name = d['name']
    cn = int(d['crossing_number'])
    ncomp = int(d['components'])
    pd = parse_pd(d['pd_notation_vector'])
    issues = []
    if parse_pd(d['pd_notation_math']) != pd:
        issues.append('pd_vector_vs_math')
    assert len(pd) == cn, name
    signs, comp, _ = orient(pd)
    if len(set(comp.values())) != ncomp:
        issues.append('pd_components')
    bl = ints(d['braid_notation'])
    strands, word = bl[0], bl[1:]
    if max(abs(g) for g in word) >= strands:
        issues.append('braid_generator_out_of_range')
    bidx = int(d['braid_index'])
    jones = poly1(d['jones_polynomial'], 'x', 1)
    jv = ints(d['jones_polynomial_vector'])
    jvd = {jv[0] + i: c for i, c in enumerate(jv[2:]) if c}
    if jvd != jones:
        issues.append('jones_string_vs_vector')
    hs = d['homflypt_polynomial'].strip()
    if hs.endswith('}'):
        hs = hs[:-1]
        issues.append('source_typo_homfly_string_trailing_brace')
    hom = poly2(hs, ['v', 'z'])
    try:
        hv = braces_json(d['homflypt_polynomial_vector'])
        ok = False
        for zs in (1, 2):
            hvd = {}
            for i, inner in enumerate(hv[2:]):
                for j, c in enumerate(inner[2:]):
                    if c:
                        hvd[(inner[0] + j, hv[0] + zs * i)] = c
            ok = ok or hvd == hom
        if not ok:
            issues.append('homfly_string_vs_vector')
    except Exception:
        issues.append('homfly_vector_unparsable')
    khp = poly2(d['khovanov_polynomial'], ['q', 't'])
    kh = sorted([h, q, c] for (q, h), c in khp.items())
    try:
        kv = braces_json(d['khovanov_polynomial_vector'])
        kvd = {}
        for i, inner in enumerate(kv[2:]):
            for j, c in enumerate(inner[2:]):
                if c:
                    kvd[(kv[0] + i, inner[0] + j)] = c
        if kvd != khp:
            issues.append('kh_string_vs_vector')
    except Exception:
        issues.append('kh_vector_unparsable')
    lm = ints(d['linking_matrix'])
    lmat = [lm[i * ncomp:(i + 1) * ncomp] for i in range(ncomp)]
    assert len(lm) == ncomp * ncomp
    sig = int(d['signature'])
    s_inv = int(d['rasmussen_s_invariant']) if d['rasmussen_s_invariant'] not in ('', 'Not known', 'undefined') else None
    nullity = int(d['nullity'])
    det = int(d['determinant'])
    if not homfly_to_jones_check(hom, jones):
        issues.append('homfly_vs_jones')
    if not kh_vs_jones(kh, jones, False):
        issues.append('kh_vs_jones')
    if s_inv is not None and not s_vs_kh(kh, s_inv):
        issues.append('s_vs_kh')
    vm1 = sum(c * (1j) ** k for k, c in jones.items())  # V at t^(1/2) = i, |V| = det
    if abs(round(abs(vm1)) - det) != 0:
        issues.append('determinant_vs_jones')
    # PD vs listed invariants
    jp = jones_pd(pd)
    sp = signature_pd(pd, 0)
    if signature_pd(pd, 1) != sp:
        issues.append('signature_colouring_disagree')
    lpd = linking_matrix(pd)
    lrel = lk_relation(lpd, lmat)
    if lrel == 'mismatch':
        issues.append('linking_matrix_vs_pd' + ('_up_to_component_order' if lk_relation(lpd, lmat, True) != 'mismatch' else ''))
    pd_rel = combine(relation(jp, jones), hrelation(homfly_pd(pd), hom), sp, sig, extra=(lrel,))
    # braid
    bpd = braid_to_pd(strands, word)
    if braid_components(strands, word) != ncomp:
        issues.append('braid_components')
    jb = jones_pd(bpd)
    sb = signature_pd(bpd, 0)
    lb = linking_matrix(bpd)
    blrel = lk_relation(lb, lmat, permute=True)
    hb = homfly_pd(bpd, BUDGET)
    br_rel = combine(relation(jb, jones), hrelation(hb, hom), sb, sig, extra=(blrel,))
    if br_rel == 'conflict':
        braid_conflicts.append((name, d['name_unoriented'], jb, hb, sb, lb))
    for lab, rel in (('pd', pd_rel), ('braid', br_rel)):
        if rel == 'same_or_mirror':
            sym = kh_mirror_symmetric(kh)
            issues.append(f"{lab}_chirality_undetermined" + ('' if sym else '_and_kh_not_mirror_symmetric'))
        elif rel != 'same':
            issues.append(f"{lab}_{rel}")
    if strands != bidx:
        issues.append('braid_strands_ne_braid_index')
    # Knot Atlas
    uname = d['name_unoriented']
    ka = None
    if uname in ka_pd_links:
        zero = all(v == 0 for v in ints(d['orientation']))
        if zero:
            ka = 'identical_pd' if ka_pd_links[uname] == pd else (
                'same_jones' if ka_link_jones[uname] == jones else 'differs')
        else:
            ka = 'reoriented'
        stats['link_knot_atlas'][ka] += 1
    rec = collections.OrderedDict()
    rec['name'] = name
    rec['name_unoriented'] = uname
    rec['orientation'] = ints(d['orientation'])
    rec['crossings'] = cn
    rec['components'] = ncomp
    rec['alternating'] = d['alternating'] == 'Y'
    rec['pd'] = pd
    rec['writhe'] = sum(signs)
    rec['braid'] = {'strands': strands, 'word': word}
    rec['braid_index'] = bidx
    rec['linking_matrix'] = lmat
    rec['jones'] = terms1(jones)
    rec['homfly'] = [[a, b, hom[(a, b)]] for (a, b) in sorted(hom, key=lambda k: (k[1], k[0]))]
    rec['kh'] = kh
    rec['signature'] = sig
    rec['nullity'] = nullity
    rec['s'] = s_inv
    rec['determinant'] = det
    rec['knot_atlas'] = ka
    rec['pd_vs_invariants'] = pd_rel
    rec['braid_vs_invariants'] = br_rel
    rec['issues'] = issues
    links_out.append(rec)
    stats['link_pd'][pd_rel] += 1
    stats['link_braid'][br_rel] += 1
    for i in issues:
        stats['link_issue'][i] += 1
    if issues:
        problems.append((name, issues))
# a braid whose closure is not the listed oriented link: look for the orientation it does give
by_name = {r['name']: r for r in links_out}
for name, uname, jb, hb, sb, lb in braid_conflicts:
    hits = []
    for r in links_out:
        if r['name_unoriented'] != uname or r['name'] == name:
            continue
        if (dict(map(tuple, r['jones'])) == jb and sb == r['signature'] and
                (hb is None or {(a, b): c for a, b, c in r['homfly']} == hb) and
                lk_relation(lb, r['linking_matrix'], permute=True) in ('same', 'symmetric')):
            hits.append(r['name'])
    rec = by_name[name]
    rec['braid_closure_matches'] = hits
    rec['issues'].append('braid_closure_has_invariants_of_' + ('/'.join(hits) if hits else 'no_listed_orientation'))
    stats['link_issue'][rec['issues'][-1]] += 1
    problems.append((name, rec['issues']))
print('links done', len(links_out), round(time.time() - t0, 1), flush=True)

os.makedirs(outdir, exist_ok=True)
with open(os.path.join(outdir, 'knots_upto12.jsonl'), 'w') as f:
    for r in knots_out:
        f.write(json.dumps(r, separators=(',', ':')) + '\n')
with open(os.path.join(outdir, 'links_upto11.jsonl'), 'w') as f:
    for r in links_out:
        f.write(json.dumps(r, separators=(',', ':')) + '\n')
summary = {
    'knots_per_crossing': dict(sorted(collections.Counter(r['crossings'] for r in knots_out).items())),
    'links_per_crossing_oriented': dict(sorted(collections.Counter(r['crossings'] for r in links_out).items())),
    'links_per_crossing_unoriented': dict(sorted(collections.Counter(r['crossings'] for r in links_out if all(v == 0 for v in r['orientation'])).items())),
    'links_per_components_oriented': dict(sorted(collections.Counter(r['components'] for r in links_out).items())),
    'stats': {k: dict(v) for k, v in stats.items()},
    'ka_selfcheck_fail': ka_selfcheck_fail,
    'ka_link_selfcheck_fail': ka_link_selfcheck_fail,
}
json.dump(summary, open(os.path.join(work, 'summary.json'), 'w'), indent=1)
json.dump(problems, open(os.path.join(work, 'problems.json'), 'w'), indent=0)
print(json.dumps(summary, indent=1))
print('total time', round(time.time() - t0, 1))
