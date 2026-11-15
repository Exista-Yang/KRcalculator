#!/usr/bin/env python3
"""Refused Search windows from KRHOM_TRACE logs -> `libbuild calib --windows` lines (M5c).

usage: trace_windows.py [--max-unknowns U] TRACE... > windows.txt
Each '[leaf_iso] refused LABELS [a] <-> [b] degree 0 (...unknowns U)' line gives Hom^0(a, b) and
Hom^0(b, a). Each '[choose_window] refused LABELS [rel] summand [lw]{s}: iota d A (...unknowns UA), pi d B
(unknowns UB)' line gives the two exact window calls Hom^A(lw, rel) and Hom^B(rel, lw); calls
above --max-unknowns are dropped (default: keep all). Output sorted by unknowns, distinct.
"""
import re, sys

iso = re.compile(r"\[leaf_iso\] refused \[([\d, ]+)\] \[([^\]]*)\] <-> \[([^\]]*)\] degree 0 \(\d+ x \d+, unknowns (\d+)\)")
pat = re.compile(r"\[choose_window\] refused \[([\d, ]+)\] \[([^\]]*)\] summand \[([^\]]*)\]\{(-?\d+)\}: iota d (-?\d+) \(\d+ x \d+, unknowns (\d+)\), pi d (-?\d+) \(unknowns (\d+)\)")


def word(w):
    w = w.strip()
    if not w:
        return "-"
    return ",".join(re.sub(r"([EF]\d+)\((\d+)\)", r"\1.\2", t) for t in w.split())


def main():
    args = sys.argv[1:]
    cap = None
    if args and args[0] == "--max-unknowns":
        cap = int(args[1])
        args = args[2:]
    out = {}
    for f in args:
        for l in open(f, errors="replace"):
            m = iso.search(l)
            if m:
                lab = ",".join(x.strip() for x in m.group(1).split(","))
                a, b, u = word(m.group(2)), word(m.group(3)), int(m.group(4))
                if cap is None or u <= cap:
                    out[(lab, a, b, 0)] = u
                    out[(lab, b, a, 0)] = u
                continue
            m = pat.search(l)
            if not m:
                continue
            lab = ",".join(x.strip() for x in m.group(1).split(","))
            rel, lw = word(m.group(2)), word(m.group(3))
            for (a, b, d, u) in ((lw, rel, int(m.group(5)), int(m.group(6))), (rel, lw, int(m.group(7)), int(m.group(8)))):
                if cap is None or u <= cap:
                    out[(lab, a, b, d)] = u
    for (lab, a, b, d), u in sorted(out.items(), key=lambda x: x[1]):
        print(f"{lab} {a} {b} {d} unknowns={u}")


main()
