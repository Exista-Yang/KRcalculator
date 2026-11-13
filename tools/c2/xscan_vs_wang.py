#!/usr/bin/env python3
"""Compare xscan JSON (2-strand sigma^(2m+1), labels k,k) unreduced closures with tests/data/wang_ref.txt."""
import json, re, sys
def parse_poly(s):
    d = {}
    for m in re.finditer(r'([+-]?\s*\d*)t\^(-?\d+)q\^(-?\d+)', s.replace(' ', '')):
        c = m.group(1)
        c = 1 if c in ('', '+') else (-1 if c == '-' else int(c))
        k = (int(m.group(2)), int(m.group(3)))
        d[k] = d.get(k, 0) + c
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
W = wang()
def key_terms(P):
    return sorted((t, q, c) for t, q, c in P['terms'] if c)
bad = 0
for f in sys.argv[1:]:
    d = json.load(open(f))
    k = d['labels'][0]; L = len(d['word']); m = (L - 1) // 2
    for r in d['results']:
        if r['reduced']:
            continue
        ref = W.get((k, r['N'], m))
        if ref is None:
            print(f, r['N'], 'no ref'); continue
        mine = {(t, q): c for t, q, c in r['P']['terms'] if c}
        ok = mine == ref
        bad += not ok
        print(f"{f} k={k} N={r['N']} m={m} rank={r['P']['rank']} {'EQUAL' if ok else 'DIFF'}")
sys.exit(1 if bad else 0)
