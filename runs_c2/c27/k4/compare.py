#!/usr/bin/env python3
"""C2-7: 4-strand (2,2,2,2) gl_4 knots from `sbox finish` JSONs (runs_c2/c27/k4/out/NAME_MODES_fF.json).
Checks per knot: chi = decat_chi (fatal inside sbox; chi_ok on every record), D14 scaling pass,
F31 = F32b (closure by closure), N' = 3 = uncolored gl_3 of the reversed braid (engine2,
runs_c2/k4/e2ref; equal or equal up to a (t, q) shift), N' = 2 rank 1. Writes compare.json.
usage: compare.py [DIR]"""
import json, os, re, sys, glob
D = sys.argv[1] if len(sys.argv) > 1 else 'runs_c2/c27/k4/out'
def parse(s):
    d = {}
    if not s.strip(): return d
    for term in s.split(' + '):
        m = re.fullmatch(r'(\d*)t\^(-?\d+)q\^(-?\d+)', term.strip())
        k = (int(m[2]), int(m[3])); d[k] = d.get(k, 0) + (int(m[1]) if m[1] else 1)
    return d
def shift(a, b):
    if sum(a.values()) != sum(b.values()) or not a: return None
    ka, kb = min(a), min(b); s = (ka[0] - kb[0], ka[1] - kb[1])
    return s if {(t + s[0], q + s[1]): c for (t, q), c in b.items()} == a else None
def cmp(a, b):
    if a == b: return 'eq'
    s = shift(a, b)
    return 'shift%s' % (s,) if s else 'DIFF'
def e2(word):
    f = 'runs_c2/k4/e2ref/%s_gl3.txt' % '_'.join(map(str, word))
    res = {}
    if os.path.exists(f):
        for line in open(f):
            m = re.match(r'word .* gl_\d (unreduced|reduced): (.*)  \(rank', line)
            if m: res[m[1] == 'reduced'] = parse(m[2])
    return res
def fmt(p):
    return ' + '.join('%st^%dq^%d' % ('' if c == 1 else c, t, q) for (t, q), c in sorted(p.items()))
runs = {}
for f in sorted(glob.glob(D + '/*.json')):
    m = re.match(r'(.*)_(r|u|ur)_f(31|32b)\.json$', os.path.basename(f))
    if not m: continue
    d = json.load(open(f))
    runs.setdefault((m[1], m[2]), {})[m[3]] = d
out = []
for (name, modes), byf in sorted(runs.items()):
    F0 = sorted(byf)[0]; d0 = byf[F0]
    word = d0['results'][0]['word']
    b = {F: {(r['N'], r['reduced']): {(t, q): c for t, q, c in r['P']['terms'] if c} for r in d['results']} for F, d in byf.items()}
    rec = {'knot': name, 'modes': modes, 'word': word, 'primes': sorted(byf), 'library_hash': {F: d['library_hash'] for F, d in byf.items()},
           'det_pass': {F: d['det_pass'] for F, d in byf.items()}, 'commit': {F: d['commit'] for F, d in byf.items()}}
    rec['primes_identical'] = len(byf) == 2 and b['31'] == b['32b']
    rec['chi_ok'] = all(r['chi_ok'] for d in byf.values() for r in d['results'])
    rec['scaling_ok'] = all(d['scaling'] and all(s.get('pass') and s.get('checked', 0) > 0 for s in d['scaling']) for d in byf.values())
    ref = e2(word)
    bb = b[F0]
    reds = sorted({red for (_, red) in bb})
    rec['N3_gl3_reversed'] = {('r' if red else 'u'): (cmp(bb[(3, red)], ref[red]) if red in ref else 'no ref') for red in reds}
    rec['N2_rank1'] = all(sum(bb[(2, red)].values()) == 1 for red in reds)
    rec['ranks'] = {'%s%d' % ('r' if red else 'u', n): sum(bb[(n, red)].values()) for n in (2, 3, 4) for red in reds}
    rec['N4'] = {('r' if red else 'u'): fmt(bb[(4, red)]) for red in reds}
    rec['N3'] = {('r' if red else 'u'): fmt(bb[(3, red)]) for red in reds}
    ok = rec['chi_ok'] and rec['scaling_ok'] and rec['primes_identical'] and rec['N2_rank1'] and all(v in ('eq',) or v.startswith('shift') for v in rec['N3_gl3_reversed'].values())
    rec['status'] = 'PASS' if ok else ('PENDING (one prime)' if len(byf) < 2 else 'FAIL')
    out.append(rec)
    print('%-6s %-3s primes %-9s F31=F32b %-5s chi %s scaling %s N3~gl3rev %s N2r1 %s ranks %s  %s' % (name, modes, ','.join(rec['primes']), rec['primes_identical'], rec['chi_ok'], rec['scaling_ok'], rec['N3_gl3_reversed'], rec['N2_rank1'], rec['ranks'], rec['status']))
    for k, v in rec['N4'].items(): print('        N4 %s: %s' % (k, v))
json.dump(out, open(os.path.join(os.path.dirname(D.rstrip('/')) or '.', 'compare.json'), 'w'), indent=1)
