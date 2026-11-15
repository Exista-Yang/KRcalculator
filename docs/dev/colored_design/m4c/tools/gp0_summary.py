#!/usr/bin/env python3
"""Summarize gp0 census jsonl files (M4c B2 G-P0) per pair: markdown table on stdout.
usage: gp0_summary.py FILE.jsonl [FILE2.jsonl ...]"""
import json, sys
from collections import OrderedDict
rows = []
for f in sys.argv[1:]:
    rows += [json.loads(l) for l in open(f) if l.strip()]
by = OrderedDict()
for r in rows:
    by.setdefault((r["n"], r["pair"]), []).append(r)
def cls(r):
    if r["feasible"]:
        return "feasible"
    u = r["max_step_unknowns"]
    if r["chain"] != "none" and 100_000 < u <= 400_000 and "Lemma A" in r["reason"] and "chain" not in r["reason"]:
        return "measurable"  # step in the 1e5-4e5 band: admissible after a cost record
    return "infeasible"
print("| N | pair | X flows | targets (feasible/measurable/infeasible) | rex class | max Lemma A step (unknowns) | chain len | pad | max template flows | max Assoc flows (theorem moves) | max chain exact unknowns | max word flows | verdict |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
tot = {"feasible": 0, "measurable": 0, "infeasible": 0}
for (n, p), rs in by.items():
    c = {"feasible": 0, "measurable": 0, "infeasible": 0}
    for r in rs:
        c[cls(r)] += 1
    v = "feasible" if c["infeasible"] == 0 and c["measurable"] == 0 else ("measurable" if c["infeasible"] == 0 else "infeasible")
    tot[v] += 1
    rex = "yes" if any(r["rex_class"] for r in rs) else ""
    ms = max(r["max_step_unknowns"] for r in rs)
    cl = max(r["chain_len"] for r in rs)
    pad = "yes" if any(r["pad"] for r in rs) else ("no chain" if any(r["chain"] == "none" for r in rs) else "no")
    tf = max(r["max_template_flows"] for r in rs)
    af = max(r["max_assoc_flows"] for r in rs)
    th = sum(r["theorem_moves"] for r in rs)
    ex = max(r["max_exact_unknowns"] for r in rs)
    wf = max(r["max_word_flows"] for r in rs)
    reason = "; ".join(sorted(set(r["reason"] for r in rs if not r["feasible"])))[:120]
    print(f"| {n} | {p} | {rs[0]['x_flows']} | {len(rs)} ({c['feasible']}/{c['measurable']}/{c['infeasible']}) | {rex} | {ms:,} | {cl} | {pad} | {tf} | {af} ({th}) | {ex:,} | {wf:,} | {v}{' — ' + reason if reason else ''} |")
print()
print("pairs:", tot)
