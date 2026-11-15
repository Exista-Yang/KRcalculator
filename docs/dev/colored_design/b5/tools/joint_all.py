# Joint 3-upright search (A1 open problem 1): for every distinct over-budget worst step t'->t (rung rho),
# all reduced words of t' (<= 5 rungs) and of t (<= 6 rungs) found by BFS (cap 250k visits); the pair
# minimizing the step score (max over supp G of min(whole, window), both directions).
# Time-bounded batches with a persistent cache (each invocation < ~15 s).
import sys, os, time, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kr
HERE = os.path.dirname(os.path.abspath(__file__))
LAM = [3, 3, 3]
res = pickle.load(open(os.path.join(HERE, 'res.pkl'), 'rb'))
WC_P = os.path.join(HERE, 'wc.pkl'); JR_P = os.path.join(HERE, 'joint.pkl')
WC = pickle.load(open(WC_P, 'rb')) if os.path.exists(WC_P) else {}
JR = pickle.load(open(JR_P, 'rb')) if os.path.exists(JR_P) else {}
T0 = time.time(); BUDGET = 13.0
def words(t, ml):
    key = (t, ml)
    if key not in WC:
        if time.time() - T0 > BUDGET:
            raise TimeoutError
        WC[key] = kr.reduced_words(LAM, t, ml, 250000)
        pickle.dump(WC, open(WC_P, 'wb'))
    return WC[key]
steps = []
for x in res:
    if x['score'] <= 100000 and not x['alt']:
        continue
    k = (x['tp'], x['t'], x['rho'])
    if k not in [s[0] for s in steps]:
        steps.append((k, x))
done = 0
try:
    for k, x in steps:
        if k in JR:
            done += 1; continue
        tp, t, rho = k
        A = words(tp, 5); B = words(t, 6)
        g = x['g']; best = []
        for a in A:
            for b in B:
                if kr.top(LAM, a + [rho]) != kr.top(LAM, b):
                    continue
                if a + [rho] == b:
                    sc = 0; wl = []
                else:
                    det, ws = kr.step_unknowns(LAM, a, rho, b, g)
                    sc = max(min(w, v) if v is not None else w for w, v in det.values()); wl = ws['labels']
                best.append((sc, kr.flows(LAM, a) + kr.flows(LAM, b), a, b, wl))
        best.sort(key=lambda e: (e[0], e[1]))
        JR[k] = dict(nA=len(A), nB=len(B), top=best[:3])
        pickle.dump(JR, open(JR_P, 'wb')); done += 1
except TimeoutError:
    pass
print(f"done {done}/{len(steps)} in {time.time()-T0:.1f}s")
