#!/usr/bin/env python3
"""Summarize libbuild run files (JSON lines) for the M5b report.

usage: libbuild_report.py RUN.jsonl [RUN.jsonl ...]

Per file: the last header, jobs certified (last attempt per job), targets and Σ mult of the
certified pair jobs (plan counts, checked by libbuild), closure types, provider-call tier histogram
and seconds (mean, p90, max) from the RUN.jsonl.targets log (every fresh provider call, appended at
once, deduplicated by (p, q, s)), max X, peak RSS, failures and the M4c list.
"""
import json, re, sys


def tier_of(summary):
    m = re.search(r"tiers \[([0-9, ]+)\](?: exact (\d+))?", summary)
    if not m:
        return "?"
    hits = [int(x) for x in m.group(1).split(",")]
    ex = int(m.group(2)) if m.group(2) else 0
    names = ["S%d" % i for i, h in enumerate(hits) if h > 0]
    if not names:
        return "Exact"
    if ex > 0:
        names.append("Exact")
    return "+".join(names)


def pct(xs, p):
    if not xs:
        return 0.0
    xs = sorted(xs)
    k = min(len(xs) - 1, int(round(p * (len(xs) - 1))))
    return xs[k]


def main(paths):
    for path in paths:
        header, last, summ = None, {}, None
        peak = 0.0
        for l in open(path):
            if l.startswith('{"header"'):
                header = json.loads(l)["header"]
            elif l.startswith('{"attempt"'):
                a = json.loads(l)
                last[a["job"]] = a
                peak = max(peak, a["peak_rss_mb"])
            elif l.startswith('{"summary"'):
                summ = json.loads(l)["summary"]
        ok = [a for a in last.values() if a["outcome"] == "ok"]
        pairs_ok = [a for a in ok if a["job"].startswith("pair")]
        tg = sum(a["want_targets"] or 0 for a in pairs_ok)
        sm = sum(a["want_mult"] or 0 for a in pairs_ok)
        clos = sum(1 for a in ok if a["job"].startswith("closure"))
        cross = [a for a in ok if a["job"].startswith("cross")]
        cross_t = sum(a["result"]["stats"]["targets"] for a in cross)
        cross_m = sum(a["result"]["stats"]["sum_mult"] for a in cross)
        maxx = max([a["x_flows"] for a in last.values()] + [0])
        seen = {}
        try:
            for l in open(path + ".targets"):
                e = json.loads(l)
                seen[(e["p"], e["q"], e["s"])] = e
        except FileNotFoundError:
            pass
        hist = {}
        secs = []
        for e in seen.values():
            t = tier_of(e["summary"])
            hist[t] = hist.get(t, 0) + 1
            secs.append(e["secs"])
        jw = [a["wall"] for a in ok]
        fails = [(j, a["outcome"], a["scale"]) for j, a in last.items() if a["outcome"] != "ok"]
        print("== %s" % path)
        if header:
            print("set %s field F%s commit %s" % (header["set"], header["field"], header["commit"][:12]))
        print("jobs %d, certified %d (%.1f%%)" % (len(last), len(ok), 100.0 * len(ok) / max(1, len(last))))
        print("pair targets %d (Σ mult %d); cross jobs %d with %d targets (Σ mult %d); closure types %d" % (tg, sm, len(cross), cross_t, cross_m, clos))
        print("fresh provider calls %d; tiers %s" % (len(seen), dict(sorted(hist.items(), key=lambda x: -x[1]))))
        if secs:
            print("provider seconds per call: mean %.2f p90 %.2f max %.2f" % (sum(secs) / len(secs), pct(secs, 0.9), max(secs)))
        if jw:
            print("job wall (certified): mean %.1f p90 %.1f max %.1f, sum %.0f s" % (sum(jw) / len(jw), pct(jw, 0.9), max(jw), sum(jw)))
        print("max X %d flows, peak RSS %.0f MB" % (maxx, peak))
        for f in fails:
            print("FAIL %s: %s at x%s" % f)
        if summ:
            print("M4c: %s" % summ.get("m4c"))


if __name__ == "__main__":
    main(sys.argv[1:])
