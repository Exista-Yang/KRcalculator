#!/usr/bin/env python3
"""G2B: compare sbox JSON results (localized path) with the exact 2-strand references
(docs/colored_design/gl6/m2s/*.json) and the Wang blocks of tests/data/wang_ref.txt.
usage: g2b_compare.py REF.json SBOX.json [SBOX.json ...]   (exit 1 on any difference, any
chi_ok != true, or a missing / failed D14 scaling record)"""
import json, re, sys

def parse_poly(s):
    d = {}
    for m in re.finditer(r'([+-]?\s*\d*)t\^(-?\d+)q\^(-?\d+)', s.replace(' ', '')):
        c = m.group(1)
        c = 1 if c in ('', '+') else (-1 if c == '-' else int(c))
        d[(int(m.group(2)), int(m.group(3)))] = d.get((int(m.group(2)), int(m.group(3))), 0) + c
    return d

def wang(path='tests/data/wang_ref.txt'):
    out, key = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith('k N m ='):
            key = tuple(int(x) for x in line.split('=')[1].split())
        elif key and line.startswith('t^'):
            out[key] = parse_poly(line); key = None
    return out

def terms(P):
    return {(t, q): c for t, q, c in P['terms']}

ref = json.load(open(sys.argv[1]))
R = {(r['stage'], r['reduced'], r['N']): terms(r['P']) for r in ref['results']}
W = wang()
bad = 0
for f in sys.argv[2:]:
    d = json.load(open(f))
    for r in d['results']:
        k = (r['stage'], r['reduced'], r['N'])
        got = terms(r['P'])
        rank = sum(got.values())
        notes = []
        # chi = decat_chi must have passed (sengine review F2): a missing or false chi_ok fails
        if r.get('chi_ok') is not True:
            bad += 1
            notes.append('CHI FAIL (chi_ok %s)' % r.get('chi_ok'))
        if k in R:
            ok = R[k] == got
            notes.append('m2s ' + ('PASS' if ok else 'FAIL'))
            if not ok:
                bad += 1
                diff = {x: (got.get(x, 0), R[k].get(x, 0)) for x in set(got) | set(R[k]) if got.get(x, 0) != R[k].get(x, 0)}
                notes.append('diff (got, want) %s' % sorted(diff.items())[:12])
        else:
            notes.append('no m2s reference')
        st = r['stage']
        if not r['reduced'] and st % 2 == 1 and st >= 3 and d['lam'] == [3, 3]:
            wk = (3, r['N'], (st - 1) // 2)
            if wk in W:
                ok = W[wk] == got
                notes.append('Wang %s' % ('PASS' if ok else 'FAIL'))
                bad += not ok
        print('%s %s stage %d %s N=%d rank %d: %s' % (d['field'], f, st, 'red' if r['reduced'] else 'unred', r['N'], rank, '; '.join(notes)))
    sc = d.get('scaling')
    scl = sc if isinstance(sc, list) else [sc]
    sc_ok = bool(scl) and all(isinstance(x, dict) and x.get('pass') is True for x in scl)
    bad += not sc_ok
    print('%s %s scaling: %s%s' % (d['field'], f, sc, '' if sc_ok else '  SCALING MISSING OR FAILED'))
print('TOTAL differences:', bad)
sys.exit(1 if bad else 0)
