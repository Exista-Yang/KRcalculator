#!/usr/bin/env python3
"""Stretch checks: (3,3,3) runs (runs_stretch/*.json) — primes x seeds identical; Markov words
sigma_1^m sigma_2 = 2-strand sigma^m (m2s) up to a global (t, q) shift."""
import json, glob, collections
def T(P): return {(t, q): c for t, q, c in P['terms']}
runs = collections.defaultdict(list)
for f in glob.glob('runs_stretch/*.json'):
    d = json.load(open(f)); runs[tuple(d['word'])].append((f, d))
ref = json.load(open('docs/colored_design/gl6/m2s/sigma9_33_F31.json'))
R = {(r['stage'], r['reduced'], r['N']): T(r['P']) for r in ref['results']}
bad = 0
for w, lst in sorted(runs.items()):
    # every record must carry chi_ok = true and a passed D14 scaling check (sengine review F1, F2)
    for f, d in lst:
        nchi = sum(1 for r in d['results'] if r.get('chi_ok') is not True)
        sc = d.get('scaling')
        scl = sc if isinstance(sc, list) else [sc]
        sc_ok = bool(scl) and all(isinstance(x, dict) and x.get('pass') is True for x in scl)
        if nchi or not sc_ok:
            bad += 1
            print(w, f, 'CERTIFICATION FAIL: %d results without chi_ok, scaling %s' % (nchi, sc))
    base = {(r['stage'], r['reduced'], r['N']): T(r['P']) for r in lst[0][1]['results']}
    same = all({(r['stage'], r['reduced'], r['N']): T(r['P']) for r in d['results']} == base for _, d in lst)
    bad += not same
    print(w, '%d runs (primes x seeds) identical: %s' % (len(lst), same))
    if w[-1] == 2 and all(x == 1 for x in w[:-1]):
        m = len(w) - 1
        st = len(w)
        for red in (False, True):
            for N in (3, 4):
                got = base[(st, red, N)]
                if m == 1:
                    want = {(0, 0): 1} if red or N == 3 else {(0, -3): 1, (0, -1): 1, (0, 1): 1, (0, 3): 1}
                else:
                    want = R.get((m, red, N))
                if want is None:
                    continue
                kg, kw = min(got), min(want)
                sh = (kg[0] - kw[0], kg[1] - kw[1])
                ok = {(t + sh[0], q + sh[1]): c for (t, q), c in want.items()} == got
                bad += not ok
                print('   [1^%d,2] %s N=%d vs 2-strand sigma^%d: %s (shift t^%d q^%d), rank %d' % (m, 'red' if red else 'unred', N, m, 'PASS' if ok else 'FAIL', sh[0], sh[1], sum(got.values())))
print('TOTAL failures', bad)
