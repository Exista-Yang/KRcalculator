# D_grass: window flows / max labels of the chosen 3-upright joint windows (pick) per step
import sys, json, time
sys.path.insert(0, sys.argv[1])
import kr
d = json.load(open(sys.argv[1] + '/steps2.json'))
def maxlab(lab, w):
    m = max(lab); l = list(lab)
    for r in w:
        l = kr.apply(l, r); m = max(m, max(l))
    return m
t0 = time.time()
for x in d:
    pairs = sorted(set(t['pair'] for t in x['targets']))
    for key in ('pick', 'joint'):
        p = x.get(key)
        if not p or p[0] == 0:
            continue
        lam = [3, 3, 3]
        a = kr.parse_word(p[1]) + kr.parse_word(x['rho'])
        b = kr.parse_word(p[2])
        ws = kr.win_split(lam, a, b)
        fa, fb = kr.flows(ws['labels'], ws['a']), kr.flows(ws['labels'], ws['b'])
        wa = kr.flows(lam, a); wb = kr.flows(lam, b)
        print(key, pairs, x['lv'], x['G'], 'u', p[0], 'lab', ws['labels'], 'winflows', fa, fb, 'maxlab', maxlab(ws['labels'], ws['a']), maxlab(ws['labels'], ws['b']), 'wordflows', wa, wb)
print('secs', round(time.time() - t0, 1))
