#!/usr/bin/env python3
"""C2-2: 4-strand (2,2,2,2) Markov / stabilized words (runs_c2/k4/markov/*.json, sbox prune 4, box N = 4)
against
 - sigma_1^m sigma_2^{+-1} sigma_3^{+-1}: the 2-strand (2,2) sigma^m xscan (runs_c2/xscan/s{m}_f{F}_p4.json);
   sigma_1 sigma_2 sigma_3: the colored unknot (N'=4: [4 choose 2]_q, rank 6; reduced 1);
 - (3-strand word) sigma_3: the C2-1 3-strand value (docs/colored_design/c2/gl4_222.json), N' = 3, 4;
 - N' = 3: uncolored gl_3 of the reversed braid (engine2, runs_c2/k4/e2ref); N' = 2: rank 1;
 - chi = decat_chi on every closure, D14 scaling record, F31 = F32b.
usage: k4_markov.py [DIR]  (prints a table, writes DIR/compare.txt and DIR/compare.json)"""
import json, os, re, sys
D = sys.argv[1] if len(sys.argv) > 1 else 'runs_c2/k4/markov'
def parse(s):
    d = {}
    if not s.strip(): return d
    for term in s.split(' + '):
        m = re.fullmatch(r'(\d*)t\^(-?\d+)q\^(-?\d+)', term.strip())
        k = (int(m[2]), int(m[3])); d[k] = d.get(k, 0) + (int(m[1]) if m[1] else 1)
    return d
def R(f):
    d = json.load(open(f))
    return d, {(r['N'], r['reduced']): {(t, q): c for t, q, c in r['P']['terms'] if c} for r in d['results']}
def shift(a, b):
    if sum(a.values()) != sum(b.values()) or not a: return None
    ka, kb = min(a), min(b); s = (ka[0] - kb[0], ka[1] - kb[1])
    return s if {(t + s[0], q + s[1]): c for (t, q), c in b.items()} == a else None
def cmp(a, b):
    if a == b: return 'eq'
    s = shift(a, b)
    return 'shift%s' % (s,) if s else 'DIFF'
def e2(word):
    f = 'runs_c2/k4/e2ref/%s_gl3.txt' % '_'.join(map(str, word))
    res = {}
    if os.path.exists(f):
        for line in open(f):
            m = re.match(r'word .* gl_\d (unreduced|reduced): (.*)  \(rank', line)
            if m: res[m[1] == 'reduced'] = parse(m[2])
    return res
k3 = {k['knot']: k for k in json.load(open('docs/colored_design/c2/gl4_222.json'))['knots']}
words = [l.strip().split(';') for l in open('runs_c2/k4/markov_words.txt') if l.strip()]
ref2 = {'u_s3': 2, 'm2': 2, 'm3': 3, 'm4': 4, 'm5': 5, 'm3n': 3, 'm3nn': 3, 'm3x': 3}
ref3 = {'w41s': '4_1', 'w52s': '5_2'}
unknot4 = {(0, -4): 1, (0, -2): 1, (0, 0): 2, (0, 2): 1, (0, 4): 1}  # [4 choose 2]_q
out, lines, bad = [], [], 0
for name, lam, w, _ in words:
    word = list(map(int, w.split()))
    runs = {F: R('%s/%s_f%s.json' % (D, name, F)) for F in ('31', '32b') if os.path.exists('%s/%s_f%s.json' % (D, name, F))}
    if not runs: continue
    F0 = sorted(runs)[0]; d0, b = runs[F0]
    rec = {'word_name': name, 'word': word, 'primes': sorted(runs), 'secs': {F: runs[F][0].get('total_secs') for F in runs},
           'points': d0.get('points'), 'point_secs_mean': d0.get('point_secs_mean'), 'vmhwm_mb': d0.get('vmhwm_kb', 0) // 1024}
    rec['primes_identical'] = len(runs) == 2 and runs['31'][1] == runs['32b'][1]
    rec['chi_ok'] = all(r['chi_ok'] for d, _ in runs.values() for r in d['results']) and all(len(d['results']) == 6 for d, _ in runs.values())
    rec['scaling_ok'] = all(d.get('scaling') and d['scaling'].get('pass') and d['scaling'].get('checked', 0) > 0 for d, _ in runs.values())
    ref = e2(word)
    rec['N3_gl3_reversed'] = {('r' if red else 'u'): (cmp(b[(3, red)], ref[red]) if red in ref else 'no ref') for red in (False, True)}
    rec['N2_rank1'] = all(b[(2, red)] == {(0, 0): 1} for red in (False, True))
    rec['ranks'] = {'%s%d' % ('r' if red else 'u', n): sum(b[(n, red)].values()) for n in (2, 3, 4) for red in (False, True)}
    rec['N4'] = {('r' if red else 'u'): ' + '.join('%st^%dq^%d' % ('' if c == 1 else c, t, q) for (t, q), c in sorted(b[(4, red)].items())) for red in (False, True)}
    chk = []
    if name in ref2:
        m = ref2[name]; fx = 'runs_c2/xscan/s%d_f%s_p4.json' % (m, F0)
        _, x = R(fx)
        chk = ['%s%d:%s' % ('r' if k[1] else 'u', k[0], cmp(b[k], x[k])) for k in sorted(set(b) & set(x))]
        rec['vs'] = '2-strand sigma^%d (xscan)' % m
    elif name == 'm1':
        chk = ['u4:%s' % cmp(b[(4, False)], unknot4), 'r4:%s' % cmp(b[(4, True)], {(0, 0): 1}), 'u3:%s' % cmp(b[(3, False)], {(0, -2): 1, (0, 0): 1, (0, 2): 1})]
        rec['vs'] = 'colored unknot'
    elif name in ref3:
        kk = k3[ref3[name]]
        for red in (False, True):
            md = 'reduced' if red else 'unreduced'
            for n in (3, 4):
                chk.append('%s%d:%s' % ('r' if red else 'u', n, cmp(b[(n, red)], parse(kk[md]['N%d' % n]))))
        rec['vs'] = '3-strand %s (C2-1)' % ref3[name]
    rec['markov'] = chk
    ok = rec['chi_ok'] and rec['scaling_ok'] and rec['N2_rank1'] and all(v.startswith(('eq', 'shift')) for v in rec['N3_gl3_reversed'].values()) and all('DIFF' not in c for c in chk)
    if len(runs) == 2: ok = ok and rec['primes_identical']
    rec['status'] = 'PASS' if ok else 'FAIL'
    bad += not ok
    out.append(rec)
    lines.append('%-5s %-22s %s  N3 vs gl3rev u/r %s/%s  N2 rank1 %s chi %s scal %s primes %s  vs %s: %s  [%s]' % (
        name, w, ' '.join('%s=%d' % kv for kv in rec['ranks'].items()), rec['N3_gl3_reversed']['u'], rec['N3_gl3_reversed']['r'],
        rec['N2_rank1'], rec['chi_ok'], rec['scaling_ok'], ','.join(rec['primes']) + (' =' if rec['primes_identical'] else ''), rec.get('vs', '-'), ' '.join(chk), rec['status']))
open(D + '/compare.txt', 'w').write('\n'.join(lines) + '\n')
json.dump(out, open(D + '/compare.json', 'w'), indent=1)
doc = {'lam': [2, 2, 2, 2], 'prune': 4, 'box_N': 4, 'box_points': 1296, 'closures': [2, 3, 4], 'seed': 201,
       'library': 'data-gl6 c2/l2222_n4 (layer A only, --replay, --costs costs_2222.txt --stamp-route)',
       'checks': "chi = decat_chi on every closure; D14 scaling; F31 = F32b; N'=3 = uncolored gl_3 of the reversed braid (engine2); N'=2 rank 1; Markov: = 2-strand sigma^m (xscan) / = 3-strand C2-1 value",
       'note': 'the 4-strand knots of knots4_cloud.txt need the M5i (idempotent) types and are not computed (C2_2_report.md section 2)',
       'computed': len(out), 'failed': [r['word_name'] for r in out if r['status'] != 'PASS'], 'words': out}
json.dump(doc, open('docs/colored_design/c2/gl4_2222.json', 'w'), indent=1)
print('\n'.join(lines)); print('%d words, %d not PASS' % (len(out), bad))
sys.exit(1 if bad else 0)
