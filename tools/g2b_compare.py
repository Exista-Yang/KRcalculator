#!/usr/bin/env python3
"""G2B: compare sbox JSON results (localized path) with the exact 2-strand references
(docs/colored_design/gl6/m2s/*.json) and the Wang blocks of tests/data/wang_ref.txt.
usage: g2b_compare.py REF.json SBOX.json [SBOX.json ...]   (exit 1 on any difference)"""
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
    print('%s %s scaling: %s' % (d['field'], f, sc))
print('TOTAL differences:', bad)
sys.exit(1 if bad else 0)
