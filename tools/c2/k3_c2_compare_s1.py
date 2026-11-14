#!/usr/bin/env python3
"""C2-1: (2,2,2) knots at N' = 2, 3, 4 (runs_c2/k3/out/*.json, sbox prune 4, box N = 4) against
 - N' = 3: the uncolored gl_3 of the reversed braid (engine2, examples/e2ref; wedge^2 C^3 = (C^3)^* (x) det),
 - N' = 2: rank 1 (wedge^2 C^2 trivial),
 - chi = schur::decat_chi (chi_ok on every closure), D14 scaling record passed, the two primes identical.
 N' = 4 (wedge^2 C^4, self-dual) has no reference: it is the new value. Also reports mirror symmetry
 P(t,q) = P(1/t,1/q) of N' = 4 for amphichiral knots (informational).
Writes docs/colored_design/c2/gl4_222_s1.json and docs/colored_design/c2/gl4_222_s1_summary.md."""
import json, os, re, subprocess, sys
K = 'runs_c2/k3/out_s1'
E2 = os.environ.get('E2REF', 'target/release/examples/e2ref')
def parse(s):
    d = {}
    for term in s.split(' + '):
        m = re.fullmatch(r'(\d*)t\^(-?\d+)q\^(-?\d+)', term.strip())
        k = (int(m[2]), int(m[3])); d[k] = d.get(k, 0) + (int(m[1]) if m[1] else 1)
    return d
def T(P): return {(t, q): c for t, q, c in P['terms'] if c}
def e2(word, n):
    os.makedirs('runs_c2/k3/e2ref', exist_ok=True)
    f = 'runs_c2/k3/e2ref/%s_gl%d.txt' % ('_'.join(map(str, word)), n)
    if not os.path.exists(f):
        out = subprocess.run([E2, ' '.join(map(str, word)), str(n), 'reverse'], capture_output=True, text=True, timeout=3600).stdout
        open(f, 'w').write(out)
    res = {}
    for line in open(f):
        m = re.match(r'word .* gl_\d (unreduced|reduced): (.*)  \(rank', line)
        if m: res[m[1] == 'reduced'] = parse(m[2])
    return res
knots = [(a, None, b, c) for a, b, c in (l.strip().split(';') for l in open('runs_c2/k3/slices/slice1.txt') if l.strip())]
out, rows, fails, partial, ncomp = [], [], [], [], 0
for name, lam, w, _ in knots:
    word = list(map(int, w.split()))
    runs = {F: json.load(open('%s/%s_f%s.json' % (K, name, F))) for F in ('31', '32b') if os.path.exists('%s/%s_f%s.json' % (K, name, F))}
    if not runs:
        continue
    R = {F: {(r['N'], r['reduced']): r for r in d['results']} for F, d in runs.items()}
    base = R['31'] if '31' in R else next(iter(R.values()))
    primes_ok = len(runs) == 2 and all({k: T(v['P']) for k, v in r.items()} == {k: T(v['P']) for k, v in base.items()} for r in R.values())
    chi_ok = all(v['chi_ok'] for r in R.values() for v in r.values()) and all(len(r) == 6 for r in R.values())
    scal_ok = all(d['scaling'] and d['scaling'].get('pass') and d['scaling'].get('checked', 0) > 0 for d in runs.values())
    ref = e2(word, 3)
    rec = {'knot': name, 'word': word, 'length': len(word), 'primes': sorted(runs), 'primes_identical': primes_ok, 'chi_ok': chi_ok,
           'scaling_ok': scal_ok, 'secs_f31': runs.get('31', {}).get('total_secs'), 'points': next(iter(runs.values()))['points']}
    ok = primes_ok and chi_ok and scal_ok
    for red in (False, True):
        g4, g3, g2 = T(base[(4, red)]['P']), T(base[(3, red)]['P']), T(base[(2, red)]['P'])
        e3 = ref.get(red)
        eq3 = e3 is not None and e3 == g3
        n2 = g2 == {(0, 0): 1}
        mir = g4 == {(-t, -q): c for (t, q), c in g4.items()}
        rec['reduced' if red else 'unreduced'] = {'N4': base[(4, red)]['P']['str'], 'rank4': sum(g4.values()),
            'N3': base[(3, red)]['P']['str'], 'rank3': sum(g3.values()), 'gl3_reversed_equal': eq3, 'N2_rank1': n2, 'N4_mirror_symmetric': mir}
        ok = ok and eq3 and n2
    rec['red_rank4_eq_red_rank3_sq'] = rec['reduced']['rank4'] == rec['reduced']['rank3'] ** 2  # informational (exponential growth, thin knots)
    if len(runs) < 2 and ok is False and chi_ok and scal_ok and all(rec[m]['gl3_reversed_equal'] and rec[m]['N2_rank1'] for m in ('unreduced', 'reduced')):
        rec['status'] = 'PARTIAL (one prime)'
    else:
        rec['status'] = 'PASS' if ok else 'FAIL'
    if rec['status'] == 'FAIL': fails.append(name)
    if rec['status'] != 'PASS': partial.append(name)
    out.append(rec); ncomp += 1
    rows.append('| %s | %s | %d | %d | %d | %d | %s / %s | %s | %s | %s | %s |' % (name, w, rec['unreduced']['rank4'], rec['reduced']['rank4'], rec['unreduced']['rank3'], rec['reduced']['rank3'],
        'eq' if rec['unreduced']['gl3_reversed_equal'] else 'NE', 'eq' if rec['reduced']['gl3_reversed_equal'] else 'NE',
        'yes' if chi_ok else 'NO', 'yes' if scal_ok else 'NO', 'yes' if primes_ok else ('only ' + ','.join(runs) if len(runs) < 2 else 'NO'), rec['status']))
doc = {'lam': [2, 2, 2], 'prune': 4, 'box_N': 4, 'closures': [2, 3, 4], 'seed': 201,
       'library': 'data-gl6/c2/l222_n4 (layer A only, --replay)',
       'checks': "chi = decat_chi on every closure; D14 scaling (2 points); F31 = F32b; N'=3 = uncolored gl_3 of the reversed braid (engine2); N'=2 rank 1",
       'computed': ncomp, 'failed': fails, 'one_prime_only': [x for x in partial if x not in fails], 'knots': out}
os.makedirs('docs/colored_design/c2', exist_ok=True)
json.dump(doc, open('docs/colored_design/c2/gl4_222_s1.json', 'w'), indent=1)
hdr = ["| knot | word | N=4 unred | N=4 red | N=3 unred | N=3 red | N=3 vs gl_3 reversed (u / r) | χ | scaling | primes | status |", "|---|---|---|---|---|---|---|---|---|---|---|"]
open('docs/colored_design/c2/gl4_222_s1_summary.md', 'w').write('\n'.join(hdr + rows) + '\n\n%d knots computed, %d PASS, failures: %s\n' % (ncomp, ncomp - len(fails), fails))
import glob
nc = []
for name, _, w, c in knots:
    if not all(os.path.exists('%s/%s_f%s.json' % (K, name, F)) for F in ('31', '32b')):
        why = []
        for F in ('31', '32b'):
            if os.path.exists('%s/%s_f%s.json' % (K, name, F)): why.append('F%s done' % F); continue
            ff = '%s/%s_f%s.fail' % (K, name, F)
            why.append('F%s: %s' % (F, open(ff).read().strip().replace(chr(10), ' ')[:140] if os.path.exists(ff) else 'not finished / not started in this session'))
        nc.append('| %s | %s | %s |' % (name, w, '; '.join(why)))
with open('docs/colored_design/c2/gl4_222_s1_summary.md', 'a') as fh:
    fh.write('\nNot computed (%d of %d slice knots; NOT check failures - resource limits or session time):\n\n| knot | word | status |\n|---|---|---|\n' % (len(nc), len(knots)) + '\n'.join(nc) + '\n')
print('%d knots computed, %d PASS (both primes), failures: %s, one prime only: %s' % (ncomp, ncomp - len(partial), fails, [x for x in partial if x not in fails]))
print('reduced rank4 == (reduced rank3)^2 for', sum(r['red_rank4_eq_red_rank3_sq'] for r in out), 'of', ncomp, '; not for', [r['knot'] for r in out if not r['red_rank4_eq_red_rank3_sq']])
