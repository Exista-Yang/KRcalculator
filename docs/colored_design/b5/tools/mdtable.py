import json, os
HERE=os.path.dirname(os.path.abspath(__file__))
steps=json.load(open(os.path.join(HERE,'steps2.json')))
B=100000
def grp(o):
    if o['rexstep']: return 'E'
    if o['nk']==2: return 'A'
    j=o['joint'][0] if o['joint'] else None; ok=o['joint_ok'][0] if o['joint_ok'] else None
    if ok is not None and ok<=B: return 'B'
    if j is not None and j<=B: return 'C'
    if o['cf'] and o['cf'].startswith('Pu'): return 'D1'
    return 'D2'
def fm(n):
    return '{:,}'.format(n).replace(',', ' ') if isinstance(n,int) else str(n)
rows=[]
for o in steps:
    g=grp(o)
    tg=sorted(set((t['n'],t['pair']) for t in o['targets']))
    tgs=', '.join(f"{p}{'(N5)' if n==5 else ''}" for n,p in tg)
    ntg=len(o['targets'])
    w=o['words'][0]
    j=o['joint']; ok=o['joint_ok']
    jtxt = f"{fm(j[0])} {j[3]}" if j else '-'
    oktxt = fm(ok[0]) if ok else 'none ≤ 4e5'
    if g=='E': jtxt='0 (Id)'; oktxt='0'
    K='{0,1}' if o['rho'][1]=='0' else '{1,2}'
    rows.append((g, -o['census'], f"| {g} | {tgs} ({ntg}) | {o['lv']}{'F' if o['final'] else ''} | {o['tp']} → {o['t']} | [{w[0]}]·{o['rho']} → [{w[1]}] | {o['G']} | {fm(o['whole'])} / {fm(o['census'])} {o['wslab']} | K={K}, cells {o['llam']}, other {o['o']}{' (mirror)' if o['mirrored'] else ''} | {fm(o['loc'])} | {jtxt} | {oktxt} | {fm(o['rin'])} / {fm(o['rout'])} | {o['cf']} |"))
rows.sort()
hdr="| grp | pairs (targets) | lv | t' → t | W_t'·ρ → W_t (census words) | G | census whole / step (window labels) | Q: K column cells, other col | 4-up chart local K-window | 3-up joint best (labels) | 3-up best with G_loc = G | detour Exact bound in / out | B1 closed form for the chosen window |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|"
open(os.path.join(HERE,'table.md'),'w').write(hdr+'\n'+'\n'.join(r[2] for r in rows)+'\n')
print(len(rows))
