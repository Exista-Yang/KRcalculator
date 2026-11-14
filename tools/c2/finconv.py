#!/usr/bin/env python3
"""Convert `sbox finish --json` output (single word) to the single-run record format of `sbox --word ... --json`
(as used by gl4_222.json / k3_c2_compare*.py). usage: finconv.py FIN.json ERRLOG OUT.json"""
import json, re, sys
fin, err, out = sys.argv[1:4]
d = json.load(open(fin))
log = open(err).read()
box = open(fin.replace('.fin.json', '.box.out')).read() if True else ''
d2 = dict(d)
d2['results'] = [{k: v for k, v in r.items() if k != 'word'} for r in d['results']]
sc = d['scaling']
d2['scaling'] = sc[0] if isinstance(sc, list) else sc
pts = sorted({x['points'] for x in d['digests']})
d2['points'] = max(pts) if pts else None
m = re.findall(r'(\d+) scanned \(mean ([\d.]+)s/point\).*?, ([\d.]+)s, VmHWM (\d+) MB', box)
secs = sum(float(x[2]) for x in m)
d2['scan_secs'] = secs
d2['box_vmhwm_mb'] = max([int(x[3]) for x in m] or [0])
m2 = re.search(r'close ([\d.]+)s, total ([\d.]+)s', log)
d2['close_secs'] = float(m2[1]) if m2 else None
d2['total_secs'] = secs + (float(m2[2]) if m2 else 0)
d2['via'] = 'sbox box + finish (checkpointed)'
json.dump(d2, open(out, 'w'))
