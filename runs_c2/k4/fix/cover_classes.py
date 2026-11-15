#!/usr/bin/env python3
"""per-class union of blocking jobs from a `libbuild cover` markdown (C2-5). usage: cover_classes.py COVER.md"""
import re, sys
A = "6_1 7_2 7_4 7_6 7_7 8_4 8_6 8_8 8_11 8_13 8_14 9_10 9_11 9_13 9_17 9_18 9_20 9_23".split()
B = "8_15 9_22 11n_42".split()
C = ["Conway"]
rows = {}
for l in open(sys.argv[1]):
    c = l.split('|')
    if len(c) > 6 and c[1].strip() not in ('knot', '---'):
        rows[c[1].strip()] = re.findall(r"((?:p\d+|cl [\d,/]+)) \((\w+)\)", c[6])
for name, ks in (("A", A), ("B", B), ("C (Conway)", C)):
    u = {}
    for k in ks:
        for j, o in rows.get(k, []):
            u[j] = o
    m5i = sorted(j for j, o in u.items() if o == 'NeedsIdempotent')
    other = sorted((j for j, o in u.items() if o != 'NeedsIdempotent'), key=lambda x: (x[0], int(x[1:]) if x[1:].isdigit() else 0))
    print(f"| {name} | {len(ks)} | {len(u)} | {len(m5i)} | {len(other)}: {', '.join(other)} |")
    per = {k: [j for j, o in rows.get(k, []) if o != 'NeedsIdempotent'] for k in ks}
    print("   per knot non-M5i:", "; ".join(f"{k}: {' '.join(v) if v else '-'}" for k, v in per.items()))
