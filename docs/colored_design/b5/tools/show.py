import pickle, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kr
res = pickle.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'res.pkl'), 'rb'))
seen = {}
for x in res:
    r = x['rec']
    key = (r['n'], r['pair'], r['q'])
    if key in seen: seen[key]['ss'].append(r['s']); continue
    x['ss'] = [r['s']]; seen[key] = x
for key, x in seen.items():
    r = x['rec']; gr = x['gr']
    ok = 'OK' if x['score'] == r['max_step_unknowns'] or x['alt'] else 'MISMATCH'
    print(f"N{r['n']} p{r['pair']} q {r['q']} s{x['ss']} mult? rep={r['max_step_unknowns']} nonfinal={r['max_nonfinal']} feas={r['feasible']}")
    print(f"   lv{x['lv']} {kr.tfmt(x['tp'])} -> {kr.tfmt(x['t'])}  levels {x['level_from']}->{x['level_to']}  final={x['final']} rex-step={x['rexstep']}")
    print(f"   [{kr.wstr(x['wfrom'])}]·{kr.wstr([x['rho']])} -> [{kr.wstr(x['wto'])}]  flows {kr.flows([3,3,3], x['wfrom']+[x['rho']])}x{kr.flows([3,3,3], x['wto'])}  G={kr.lp_str(x['g'])}")
    print(f"   per d (whole, window): {x['det']}  score={x['score']} {ok}")
    print(f"   window: at={x['ws']['at']} suf={x['ws']['suffix']} labels={x['ws']['labels']} [{kr.wstr(x['ws']['a'])}] -> [{kr.wstr(x['ws']['b'])}]")
    print(f"   Q={kr.tfmt(gr['q'])} local_lam={gr['local_lam']} mK {gr['mK_from']}->{gr['mK_to']} cat={x['cat']} local={x['local']['best']} g_ok={x['local']['g_ok']} nw={x['local']['nwf']},{x['local']['nwt']}")
    if x['alt']:
        a = x['alt']
        print(f"   ALT final word [{a['word']}] score={a['score']} det={a['det']} win labels={a['ws']['labels']} [{kr.wstr(a['ws']['a'])}]->[{kr.wstr(a['ws']['b'])}]")
