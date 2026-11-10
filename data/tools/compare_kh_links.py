import sys, os, json, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from krlib import parse_poly, split_top
kt = sys.argv[2]
def mlist(path, var):
    text = open(path).read(); i = text.index(var + ' = {'); j = text.index('\n\n', i)
    body = text[i + len(var) + 4: j].strip(); return [s.strip() for s in split_top(body[:-1])]
names = []
for chunk in re.split(r'\nPD\[Link\[', open(os.path.join(kt, 'PD4Links.m')).read())[1:]:
    n, t, k = chunk.split(']]=', 1)[0].replace(' ', '').split(',')
    names.append(f"L{n}{'a' if t == 'Alternating' else 'n'}{k}")
ka = {}
for nm, s in zip(names, mlist(os.path.join(kt, 'Kh4Links.m'), 'Is')):
    ka[nm] = sorted([int(h), int(q), c] for (q, h), c in parse_poly(s, ['q', 't']).t.items())
same = diff = 0; bad = []
for l in open(sys.argv[1]):
    r = json.loads(l)
    if any(r['orientation']): continue
    if ka[r['name_unoriented']] == r['kh']: same += 1
    else: diff += 1; bad.append(r['name'])
print('Kh LinkInfo{0..} vs Knot Atlas Kh4Links: same', same, 'diff', diff, bad[:10])
