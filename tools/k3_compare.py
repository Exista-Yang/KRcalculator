#!/usr/bin/env python3
"""M6a-1c: (3,3,3) knots at N = 3, 4 (runs_m6a1c/k3/*.json) against
 - the uncolored gl_4 of the reversed braid (engine2, examples/e2ref; "via duality (expected)"),
 - N = 3: rank 1 (wedge^3 C^3 trivial),
 - chi = schur::decat_chi (sbox's chi_ok), the two primes identical.
Writes docs/colored_design/gl6/gl4_333.json and prints a summary table (markdown)."""
import json, os, re, subprocess, sys
K = 'runs_m6a1c/k3'
E2 = 'target/release/examples/e2ref'
def parse(s):
    d = {}
    for term in s.split(' + '):
        m = re.fullmatch(r't\^(-?\d+)q\^(-?\d+)', term.strip())
        if m:
            k = (int(m[1]), int(m[2])); d[k] = d.get(k, 0) + 1
            continue
        m = re.fullmatch(r'(\d+)t\^(-?\d+)q\^(-?\d+)', term.strip())
        k = (int(m[2]), int(m[3])); d[k] = d.get(k, 0) + int(m[1])
    return d
def T(P): return {(t, q): c for t, q, c in P['terms']}
def e2(word):
    os.makedirs(K + '/e2ref', exist_ok=True)
    f = K + '/e2ref/' + '_'.join(map(str, word)) + '.txt'
    if not os.path.exists(f):
        out = subprocess.run([E2, ' '.join(map(str, word)), '4', 'reverse'], capture_output=True, text=True, timeout=3600).stdout
        open(f, 'w').write(out)
    res = {}
    for line in open(f):
        m = re.match(r'word .* gl_4 (unreduced|reduced): (.*)  \(rank', line)
        if m: res[m[1] == 'reduced'] = parse(m[2])
    return res
def shift(got, want):
    if not got or not want or sum(got.values()) != sum(want.values()): return None
    kg, kw = min(got), min(want); sh = (kg[0] - kw[0], kg[1] - kw[1])
    return sh if {(t + sh[0], q + sh[1]): c for (t, q), c in want.items()} == got else None
knots = [l.strip().split(';') for l in open('tools/ntrunc_t1/data/knots3_all.txt') if l.strip()]
out, rows, fails = [], [], []
for name, lam, w, _ in knots:
    word = list(map(int, w.split()))
    runs = {}
    for F in ('31', '32b'):
        f = '%s/%s_f%s.json' % (K, name, F)
        if os.path.exists(f):
            runs[F] = json.load(open(f))
    if not runs:
        errf = '%s/%s_f31.err' % (K, name)
        if os.path.exists(errf):
            last = [l for l in open(errf) if l.strip()][-3:]
            out.append({'knot': name, 'word': word, 'status': 'not computed', 'reason': ' | '.join(x.strip() for x in last)})
        continue
    R = {F: {(r['N'], r['reduced']): r for r in d['results']} for F, d in runs.items()}
    base = next(iter(R.values()))
    primes_ok = all({k: T(v['P']) for k, v in r.items()} == {k: T(v['P']) for k, v in base.items()} for r in R.values())
    chi_ok = all(v['chi_ok'] for r in R.values() for v in r.values())
    ref = e2(word)
    rec = {'knot': name, 'word': word, 'length': len(word), 'primes': sorted(runs), 'primes_identical': primes_ok, 'chi_ok': chi_ok,
           'secs_f31': runs.get('31', {}).get('total_secs'), 'points': next(iter(runs.values()))['points']}
    ok = primes_ok and chi_ok
    for red in (False, True):
        g4 = T(base[(4, red)]['P']); g3 = T(base[(3, red)]['P'])
        sh = shift(g4, ref[red])
        key = 'reduced' if red else 'unreduced'
        rec[key] = {'N4': base[(4, red)]['P']['str'], 'rank4': sum(g4.values()), 'N3': base[(3, red)]['P']['str'], 'rank3': sum(g3.values()),
                    'gl4_reversed_rank': sum(ref[red].values()), 'duality': 'equal' if sh == (0, 0) else ('shift t^%d q^%d' % sh if sh else 'MISMATCH'),
                    'n3_rank1': sum(g3.values()) == 1}
        ok = ok and sh is not None and sum(g3.values()) == 1
    rec['status'] = 'PASS' if ok else 'FAIL'
    if not ok: fails.append(name)
    out.append(rec)
    rows.append('| %s | %s | %d | %d | %s | %d | %s | %s | %s |' % (name, w, rec['unreduced']['rank4'], rec['reduced']['rank4'], rec['unreduced']['duality'] + ' / ' + rec['reduced']['duality'],
                rec['unreduced']['rank3'], 'yes' if chi_ok else 'NO', 'yes' if primes_ok else 'NO', rec['status']))
json.dump({'lam': [3, 3, 3], 'prune': 4, 'box_N': 4, 'closures': [3, 4], 'comparison': 'uncolored gl_4 of the reversed braid (engine2), via duality (expected)', 'knots': out}, open('docs/colored_design/gl6/gl4_333.json', 'w'), indent=1)
print('| knot | word | N=4 unred | N=4 red | vs gl_4 reversed (unred / red) | N=3 rank | χ | primes | status |')
print('|---|---|---|---|---|---|---|---|---|')
print('\n'.join(rows))
n = sum(1 for r in out if r['status'] in ('PASS', 'FAIL'))
print('\ncomputed %d, PASS %d, FAIL %d %s; not computed %d' % (n, n - len(fails), len(fails), fails, sum(1 for r in out if r['status'] == 'not computed')))
