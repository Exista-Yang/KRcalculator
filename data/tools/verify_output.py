"""Sanity checks of the generated tables against a few literature values and the schema.
usage: verify_output.py <dir with knots_upto12.jsonl and links_upto11.jsonl>"""
import sys, os, json, collections

d = sys.argv[1]
knots = [json.loads(l) for l in open(os.path.join(d, 'knots_upto12.jsonl'))]
links = [json.loads(l) for l in open(os.path.join(d, 'links_upto11.jsonl'))]
K = {r['name']: r for r in knots}
L = {r['name']: r for r in links}
fails = []


def check(cond, msg):
    if not cond:
        fails.append(msg)


# counts (known enumerations)
ck = collections.Counter(r['crossings'] for r in knots)
check([ck[n] for n in range(3, 13)] == [1, 1, 2, 3, 7, 21, 49, 165, 552, 2176], f'knot counts {ck}')
cu = collections.Counter(r['crossings'] for r in links if not any(r['orientation']))
check([cu[n] for n in (2, 4, 5, 6, 7, 8, 9, 10, 11)] == [1, 1, 1, 6, 9, 29, 83, 287, 1007], f'link counts {cu}')
for r in links:
    check(len(r['orientation']) == r['components'] - 1, r['name'] + ' orientation length')
per = collections.Counter(r['name_unoriented'] for r in links)
for r in links:
    check(per[r['name_unoriented']] == 2 ** (r['components'] - 1), r['name'] + ' number of orientations')

# orientation flags: linking matrix = {0,..,0} matrix with sign changes
base = {r['name_unoriented']: r for r in links if not any(r['orientation'])}
for r in links:
    b = base[r['name_unoriented']]
    eps = [1] + [-1 if o else 1 for o in r['orientation']]
    n = r['components']
    check(r['linking_matrix'] == [[eps[i] * eps[j] * b['linking_matrix'][i][j] for j in range(n)] for i in range(n)],
          r['name'] + ' orientation vs linking matrix')

# schema / PD sanity
for r in knots + links:
    labels = collections.Counter(l for x in r['pd'] for l in x)
    check(len(r['pd']) == r['crossings'], r['name'] + ' pd length')
    check(all(v == 2 for v in labels.values()) and sorted(labels) == list(range(1, 2 * r['crossings'] + 1)),
          r['name'] + ' pd labels')
    b = r['braid']
    check(all(0 < abs(g) < b['strands'] for g in b['word']), r['name'] + ' braid range')
    for k in ('jones', 'homfly', 'kh', 'signature', 'braid_index', 'pd_vs_invariants'):
        check(k in r, r['name'] + ' missing ' + k)

# literature values
t = K['3_1']
check(t['jones'] == [[2, 1], [6, 1], [8, -1]], '3_1 jones (right-handed: t+t^3-t^4)')
check(t['writhe'] == 3 and t['signature'] == -2 and t['s'] == 2, '3_1 writhe/signature/s')
check(t['kh'] == [[0, 1, 1], [0, 3, 1], [2, 5, 1], [3, 9, 1]] and t['kh_torsion'] == [[3, 7, 2, 1]], '3_1 Kh')
check(t['kh_reduced'] == [[0, 2, 1], [2, 6, 1], [3, 8, 1]], '3_1 reduced Kh')
check(sorted(map(tuple, t['homfly'])) == sorted([(2, 0, 2), (4, 0, -1), (2, 2, 1)]), '3_1 HOMFLY 2v^2-v^4+v^2z^2')
check(t['knot_atlas'] == 'mirror', '3_1 is the mirror of the Knot Atlas 3_1')
check(K['4_1']['jones'] == [[-4, 1], [-2, -1], [0, 1], [2, -1], [4, 1]], '4_1 jones')
check(K['5_1']['jones'] == [[4, 1], [8, 1], [10, -1], [12, 1], [14, -1]] and K['5_1']['signature'] == -4, '5_1')
p = K['12n_242']  # (-2,3,7) pretzel, a positive braid knot of genus 5
check(p['s'] == 10 and p['signature'] == -8 and all(g > 0 for g in p['braid']['word']), '12n_242')
c, kt = K['11n_34'], K['11n_42']  # Conway / Kinoshita-Terasaka mutants
check(c['jones'] == kt['jones'] and sorted(map(tuple, c['homfly'])) == sorted(map(tuple, kt['homfly'])), 'Conway/KT')
check(c['s'] == 0 and kt['s'] == 0, 'KT s=0; Conway s=0')
h = L['L2a1{1}']  # positive Hopf link
check(h['jones'] == [[1, -1], [5, -1]] and h['linking_matrix'] == [[0, 1], [1, 0]], 'positive Hopf')
check(h['kh'] == [[0, 0, 1], [0, 2, 1], [2, 4, 1], [2, 6, 1]] and h['signature'] == -1, 'positive Hopf Kh/sig')
check(L['L2a1{0}']['jones'] == [[-5, -1], [-1, -1]], 'negative Hopf')
bor = L['L6a4{0,0}']
check(bor['components'] == 3 and bor['linking_matrix'] == [[0, 0, 0]] * 3, 'Borromean rings')
check(L['L4a1{1}']['braid'] == {'strands': 2, 'word': [1, 1, 1, 1]}, 'T(2,4) braid')

print('knots', len(knots), 'links', len(links), 'unoriented links', len(per))
print('FAILURES:' if fails else 'all checks passed', fails)
