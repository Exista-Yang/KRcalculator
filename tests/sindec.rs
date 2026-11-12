//! M5a acceptance tests (SPEC_M4_M6a §5.3 D1, D3-D5, D7-D10; SPEC_3C_GL6 rev1 §2 M5a):
//! sindec with Foam_N pruning, the DecHom guards, layer A, freeze / restrict, and a mini
//! localized scan compared with the exact xengine.
//!
//! Default suite: (2,1), (2,2), (3,2), (1,1,1). The (3,3) comparisons and the (3,3,3) plan are
//! `#[ignore]` (cloud).
use krhom::crossing::Gen;
use krhom::field::{Field, F31};
use krhom::homprov::{Budget, HomProvider, Missing, Pairs, SearchStats};
use krhom::schur::{crossing_terms, Schur, Tid};
use krhom::sexpr::Expr;
use krhom::sindec::*;
use krhom::smap::{points, Arena, Eval, ExprId, MapId, MapNode, Point};
use krhom::web::Rung;
use krhom::xengine::{xscan, XOpts};
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;

type F = F31;

fn costs() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/homprov/costs_cloud.txt")
}
fn cache_dir(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("m5a_cache").join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}
fn builder(lam: &[u8], prune: Option<u8>, cache: Option<PathBuf>) -> SIndec<F> {
    let (prov, stamp) = provider_for_stamped::<F>(lam, Some(&costs()), 1.0);
    SIndec::with_stamp(lam, prov, stamp, SIndecOpts { cache_dir: cache, prune_n: prune, ..Default::default() })
}
fn tid(s: &mut Schur, t: &str) -> Tid {
    s.intern(krhom::schur::Table::parse(t))
}

// ------------------------------------------------------------------ mini localized scan

/// A complex of canonical ladders with Arena maps as components (test-only engine).
struct Cx {
    objs: BTreeMap<usize, (Tid, i32, i32)>,
    out: BTreeMap<usize, BTreeMap<usize, MapId>>,
    next: usize,
}
impl Cx {
    fn new() -> Cx {
        Cx { objs: BTreeMap::new(), out: BTreeMap::new(), next: 0 }
    }
    fn add_obj(&mut self, o: (Tid, i32, i32)) -> usize {
        self.objs.insert(self.next, o);
        self.out.insert(self.next, BTreeMap::new());
        self.next += 1;
        self.next - 1
    }
}
fn nonzero(ar: &Arena<F>, m: MapId, pts: &[Point<F>]) -> bool {
    pts.iter().any(|p| !Eval::new(ar, p.clone()).is_zero(m))
}
fn add_edge(sd: &mut SIndec<F>, cx: &mut Cx, a: usize, b: usize, m: MapId, pts: &[Point<F>]) {
    let m = match cx.out[&a].get(&b) {
        Some(&old) => sd.ar.add(MapNode::Lin(vec![(F::one(), old), (F::one(), m)]), "cx sum").unwrap(),
        None => m,
    };
    if nonzero(&sd.ar, m, pts) {
        cx.out.get_mut(&a).unwrap().insert(b, m);
    } else {
        cx.out.get_mut(&a).unwrap().remove(&b);
    }
}
fn kept_terms(c: u8, d: u8, positive: bool, i: usize, prune: Option<u8>) -> Vec<usize> {
    let nt = crossing_terms(c, d, positive, i).len();
    (0..nt).filter(|&k| prune.is_none_or(|n| c.max(d) as usize + k <= n as usize)).collect()
}
fn tensor_crossing(sd: &mut SIndec<F>, cx: &Cx, i: usize, positive: bool, pts: &[Point<F>]) -> Cx {
    let prune = sd.prune_n();
    let mut out = Cx::new();
    let mut place: BTreeMap<(usize, usize), Vec<usize>> = BTreeMap::new();
    let mut tops: BTreeMap<usize, (u8, u8)> = BTreeMap::new();
    for (&o, &(p, h, t)) in &cx.objs {
        let mu = sd.schur.table(p).cols();
        let (c, d) = (mu[i], mu[i + 1]);
        tops.insert(o, (c, d));
        let terms = crossing_terms(c, d, positive, i);
        for k in kept_terms(c, d, positive, i, prune) {
            let dec = sd.dec_word(p, &terms[k].rungs).unwrap();
            let ids = dec.summands.iter().map(|sm| out.add_obj((sm.q, h + terms[k].hdeg, t + terms[k].dq + sm.s))).collect();
            place.insert((o, k), ids);
        }
    }
    // d_C ⊗ id
    let edges: Vec<(usize, usize, MapId)> = cx.out.iter().flat_map(|(&a, m)| m.iter().map(move |(&b, &f)| (a, b, f))).collect();
    for (a, b, f) in edges {
        let (c, d) = tops[&a];
        let terms = crossing_terms(c, d, positive, i);
        let (pa, pb) = (cx.objs[&a].0, cx.objs[&b].0);
        for k in kept_terms(c, d, positive, i, prune) {
            let rs = terms[k].rungs.clone();
            let da = sd.dec_word(pa, &rs).unwrap();
            let db = sd.dec_word(pb, &rs).unwrap();
            let kf = if rs.is_empty() { f } else { sd.ar.add(MapNode::Kron(f, rs.clone()), "cx kron").unwrap() };
            for (si, sg) in da.summands.iter().enumerate() {
                let a1 = sd.ar.add(MapNode::Comp(kf, sg.iota), "cx").unwrap();
                for (ti, tu) in db.summands.iter().enumerate() {
                    let m = sd.ar.add(MapNode::Comp(tu.pi, a1), "cx").unwrap();
                    if nonzero(&sd.ar, m, pts) {
                        let (x, y) = (place[&(a, k)][si], place[&(b, k)][ti]);
                        add_edge(sd, &mut out, x, y, m, pts);
                    }
                }
            }
        }
    }
    // (−1)^h id ⊗ d
    for (&o, &(p, h, _)) in &cx.objs {
        let (c, d) = tops[&o];
        let kt = kept_terms(c, d, positive, i, prune);
        for j in 0..crossing_terms(c, d, positive, i).len().saturating_sub(1) {
            let (ks, ktt) = if positive { (j, j + 1) } else { (j + 1, j) };
            if !kt.contains(&ks) || !kt.contains(&ktt) {
                continue;
            }
            let comps = sd.rickard_comp_at(p, i, positive, j).unwrap();
            for &(si, ti, m) in comps.iter() {
                let m = if h.rem_euclid(2) == 1 { sd.ar.add(MapNode::Lin(vec![(F::zero().sub(&F::one()), m)]), "cx sign").unwrap() } else { m };
                let (x, y) = (place[&(o, ks)][si], place[&(o, ktt)][ti]);
                add_edge(sd, &mut out, x, y, m, pts);
            }
        }
    }
    out
}
/// Gaussian elimination of every isomorphism component (same type, degree 0, nonzero scalar)
fn simplify(sd: &mut SIndec<F>, cx: &mut Cx, pts: &[Point<F>]) {
    loop {
        let mut piv = None;
        'f: for (&a, outs) in &cx.out {
            let (qa, _, ta) = cx.objs[&a];
            for (&b, &m) in outs {
                let (qb, _, tb) = cx.objs[&b];
                if qa == qb && ta == tb {
                    let tab = sd.schur.table(qa).clone();
                    let cs = krhom::coset::CosetSpace::new(&sd.schur.lam.clone(), &tab.cols());
                    let c = cs.rank(&cs.rep(&tab));
                    let v = Eval::new(&sd.ar, pts[0].clone()).block(m, c).a[0].clone();
                    if !v.is_zero() {
                        piv = Some((a, b, v));
                        break 'f;
                    }
                }
            }
        }
        let Some((a, b, v)) = piv else { break };
        let cinv = v.inv();
        let ins_b: Vec<(usize, MapId)> = cx.out.iter().filter(|(x, o)| **x != a && o.contains_key(&b)).map(|(x, o)| (*x, o[&b])).collect();
        let outs_a: Vec<(usize, MapId)> = cx.out[&a].iter().filter(|(y, _)| **y != b).map(|(y, m)| (*y, *m)).collect();
        for &(x, dxb) in &ins_b {
            for &(y, day) in &outs_a {
                let comp = sd.ar.add(MapNode::Comp(day, dxb), "cx elim").unwrap();
                let m = sd.ar.add(MapNode::Lin(vec![(F::zero().sub(&cinv), comp)]), "cx elim").unwrap();
                add_edge(sd, cx, x, y, m, pts);
            }
        }
        for x in [a, b] {
            cx.objs.remove(&x);
            cx.out.remove(&x);
            for o in cx.out.values_mut() {
                o.remove(&x);
            }
        }
    }
}
/// final objects (table, h, t) of the mini scan of a 2- or 3-strand word (±(i+1))
fn mini_scan(sd: &mut SIndec<F>, word: &[i32]) -> Vec<(String, i32, i32)> {
    let pts = points::<F>(sd.ar.k(), 0x6d69, 2, 0);
    let lam = sd.schur.lam.clone();
    let p0 = sd.schur.intern(krhom::schur::Table::identity(&lam));
    let mut cx = Cx::new();
    cx.add_obj((p0, 0, 0));
    for &g in word {
        let i = (g.unsigned_abs() - 1) as usize;
        let mut n = tensor_crossing(sd, &cx, i, g > 0, &pts);
        simplify(sd, &mut n, &pts);
        cx = n;
    }
    let mut v: Vec<(String, i32, i32)> = cx.objs.values().map(|&(p, h, t)| (sd.schur.table(p).fmt(), h, t)).collect();
    v.sort();
    v
}
fn gens(w: &[i32]) -> Vec<Gen> {
    w.iter().map(|&g| Gen { i: (g.unsigned_abs() - 1) as usize, positive: g > 0 }).collect()
}
fn xengine_objects(lam: &[u8], word: &[i32], prune: Option<u8>) -> Vec<(String, i32, i32)> {
    let l: Vec<usize> = lam.iter().map(|&a| a as usize).collect();
    let s = xscan::<F>(&l, &gens(word), &XOpts { prune_n: prune.map(|n| n as usize), ..Default::default() }).unwrap();
    s.final_objects.clone()
}

// ------------------------------------------------------------------ D1 types

#[test]
fn d1_types_small() {
    for lam in [vec![2u8, 1], vec![1, 2], vec![2, 2], vec![3, 2], vec![1, 1, 1], vec![3, 1]] {
        let mut sd = builder(&lam, None, None);
        let mut n = 0;
        for mu in visited_blocks(&lam) {
            let order = sd.schur.canonical_basis(&mu).order.clone();
            for p in order {
                let nm = sd.schur.table(p).fmt();
                let st = sd.ty(p).unwrap_or_else(|e| panic!("{:?} {}: {}", lam, nm, e));
                assert!(st.e.is_none());
                n += 1;
            }
        }
        assert!(n > 0);
        assert!(decomposable_ladders(&mut sd.schur).is_empty(), "{:?} has a decomposable canonical ladder", lam);
    }
}

/// kept (3,3,3) types per visited block at N = 4, 5, 6 (R2 §1 table, NTRUNC: 40 of 55 at N = 6)
#[test]
fn d1_kill_counts_333() {
    let mut s = Schur::new(&[3, 3, 3]);
    let mut tot = [0usize; 10];
    let mut per6 = BTreeMap::new();
    let mut all = 0;
    for mu in visited_blocks(&[3, 3, 3]) {
        let ts = krhom::schur::tables(&[3, 3, 3], &mu);
        all += ts.len();
        let mut k6 = 0;
        for t in ts {
            let q = s.intern(t);
            for n in 3..10 {
                if !s.killed_at(q, n) {
                    tot[n] += 1;
                }
            }
            if !s.killed_at(q, 6) {
                k6 += 1;
            }
        }
        per6.insert(mu.clone(), k6);
    }
    assert_eq!(all, 219);
    assert_eq!((tot[4], tot[5], tot[6]), (12, 58, 138));
    assert_eq!(per6[&vec![3, 3, 3]], 40);
    assert_eq!(per6[&vec![4, 2, 3]], 30);
    assert_eq!(per6[&vec![5, 1, 3]], 15);
    assert_eq!(per6[&vec![6, 0, 3]], 4);
    assert_eq!(tot[9], 219);
}

// ------------------------------------------------------------------ D3, D4, D10

fn all_pairs(sd: &mut SIndec<F>, plan: &Plan) -> (usize, usize, usize) {
    let (mut ok, mut ak, mut kp) = (0, 0, 0);
    for pp in &plan.pairs {
        let r = sd.dec_rung(pp.p, pp.rung);
        if pp.p_killed {
            assert!(matches!(r, Err(DecErr::Killed { .. })), "killed source accepted");
            kp += 1;
            continue;
        }
        let d = r.unwrap_or_else(|e| panic!("{} ⊗ {:?}: {}", sd.schur.table(pp.p).fmt(), pp.rung, e));
        assert_eq!(d.all_killed, pp.all_killed);
        let mut got: BTreeMap<(Tid, i32), usize> = BTreeMap::new();
        for sm in &d.summands {
            *got.entry((sm.q, sm.s)).or_insert(0) += 1;
        }
        let want: BTreeMap<(Tid, i32), usize> = pp.kept.iter().map(|t| ((t.q, t.s), t.mult)).collect();
        assert_eq!(got, want);
        if d.all_killed {
            ak += 1;
        }
        ok += 1;
    }
    (ok, ak, kp)
}

#[test]
fn d3_d10_dec_rung_small() {
    for (lam, ns) in [(vec![2u8, 1], vec![None, Some(2u8), Some(3)]), (vec![2, 2], vec![None, Some(2), Some(3)]), (vec![3, 2], vec![None, Some(3), Some(4)]), (vec![1, 1, 1], vec![None, Some(1), Some(2)])] {
        for n in ns {
            let t0 = Instant::now();
            let mut sd = builder(&lam, n, None);
            let plan = pruned_plan(&mut sd.schur, n);
            let calls0 = sd.stats().provider_calls;
            let (ok, ak, kp) = all_pairs(&mut sd, &plan);
            let st = sd.stats();
            eprintln!("[d3] {:?} N={:?}: {} pairs ok, {} all-killed, {} killed sources; {} provider calls ({:.1}s), {}; {:.1}s", lam, n, ok, ak, kp, st.provider_calls - calls0, st.provider_secs, plan.summary(), t0.elapsed().as_secs_f64());
            assert_eq!(st.targets, plan.counts().4, "targets of the plan != provider targets");
            assert_eq!(st.sum_mult, plan.counts().5);
            if n == Some(1) && lam == vec![1, 1, 1] {
                assert!(ak > 0);
            }
        }
    }
}

/// D10: killed p, all-killed without provider calls, Rickard with a killed term
#[test]
fn d10_killed() {
    let lam = [2u8, 2];
    let mut sd = builder(&lam, Some(2), None);
    // the anti-diagonal (2,2) table (cell (4)) is killed at N = 2
    let anti = tid(&mut sd.schur, "0,2/2,0");
    assert!(matches!(sd.ty(anti), Err(DecErr::Killed { n: 2, .. })));
    assert!(matches!(sd.dec_rung(anti, Rung::e(0, 1)), Err(DecErr::Killed { .. })));
    // identity type: W_1, W_2 killed at N = 2 (c = d = 2: label 2 + k > 2)
    let p0 = tid(&mut sd.schur, "2,0/0,2");
    assert!(matches!(sd.rickard_comp_at(p0, 0, true, 0), Err(DecErr::Killed { .. })));
    // all-killed decompositions make no provider call
    let plan = pruned_plan(&mut sd.schur, Some(2));
    let before = sd.stats().provider_calls;
    let mut n_ak = 0;
    for pp in plan.pairs.iter().filter(|p| !p.p_killed && p.all_killed) {
        let d = sd.dec_rung(pp.p, pp.rung).unwrap();
        assert!(d.all_killed && d.summands.is_empty());
        n_ak += 1;
    }
    assert!(n_ak > 0);
    assert_eq!(sd.stats().provider_calls, before);
}

#[test]
fn d4_dec_word_crossing_terms() {
    for (lam, ns) in [(vec![2u8, 1], vec![None, Some(2)]), (vec![1, 2], vec![None]), (vec![2, 2], vec![None, Some(3)]), (vec![3, 1], vec![None, Some(3)]), (vec![3, 2], vec![None, Some(3), Some(4)]), (vec![1, 1, 1], vec![None, Some(2)])] {
        for n in ns {
            let mut sd = builder(&lam, n, None);
            let order = sd.schur.canonical_basis(&lam).order.clone();
            let mut cnt = 0;
            for p in order {
                if sd.killed(p) {
                    continue;
                }
                for i in 0..lam.len() - 1 {
                    for pos in [true, false] {
                        for t in crossing_terms(lam[i], lam[i + 1], pos, i) {
                            // the kept multiset assert runs inside dec_word
                            let d = sd.dec_word(p, &t.rungs).unwrap_or_else(|e| panic!("{:?} N={:?} {}: {}", lam, n, sd.schur.table(p).fmt(), e));
                            let mut got: BTreeMap<(Tid, i32), usize> = BTreeMap::new();
                            for sm in &d.summands {
                                *got.entry((sm.q, sm.s)).or_insert(0) += 1;
                            }
                            let mut want: BTreeMap<(Tid, i32), usize> = BTreeMap::new();
                            for (q, m) in sd.schur.term_mult(p, &t.rungs) {
                                if !sd.killed(q) {
                                    for (k, &c) in m.c.iter().enumerate() {
                                        if c != 0 {
                                            want.insert((q, m.lo + k as i32), c as usize);
                                        }
                                    }
                                }
                            }
                            assert_eq!(got, want);
                            cnt += 1;
                        }
                    }
                }
            }
            eprintln!("[d4] {:?} N={:?}: {} crossing-term decompositions, stats {:?}", lam, n, cnt, sd.stats());
        }
    }
}

// ------------------------------------------------------------------ D5 against xengine

fn d5_case(lam: &[u8], word: &[i32], prune: Option<u8>) {
    let t0 = Instant::now();
    let want = xengine_objects(lam, word, prune);
    let tx = t0.elapsed().as_secs_f64();
    let mut sd = builder(lam, prune, None);
    let got = mini_scan(&mut sd, word);
    eprintln!("[d5] {:?} {:?} N={:?}: {} final objects (xengine {:.1}s, sindec {:.1}s, provider {:.1}s)", lam, word, prune, got.len(), tx, t0.elapsed().as_secs_f64() - tx, sd.stats().provider_secs);
    assert_eq!(got, want, "{:?} {:?} N={:?}: sindec mini scan != xengine", lam, word, prune);
}

#[test]
fn d5_against_xengine_small() {
    d5_case(&[2, 1], &[1, 1], None);
    d5_case(&[2, 1], &[1, -1, 1, 1], None);
    d5_case(&[1, 1], &[1, 1, 1], None);
    d5_case(&[2, 2], &[1, 1], None);
    d5_case(&[2, 2], &[1, -1], None);
    d5_case(&[2, 2], &[1, 1, 1], Some(2));
    d5_case(&[2, 2], &[1, 1, 1], Some(3));
    d5_case(&[3, 2], &[1, 1], None);
    d5_case(&[3, 2], &[1, 1], Some(3));
    d5_case(&[3, 2], &[1, 1], Some(4));
}

/// 3 strands: the braid relation and Reidemeister II on the minimal complexes of the mini scan
#[test]
fn d5_three_strand_braid_relation() {
    for prune in [None, Some(2u8)] {
        let mut a = builder(&[1, 1, 1], prune, None);
        let x = mini_scan(&mut a, &[1, 2, 1]);
        let mut b = builder(&[1, 1, 1], prune, None);
        let y = mini_scan(&mut b, &[2, 1, 2]);
        assert_eq!(x, y, "N={:?}: σ1σ2σ1 != σ2σ1σ2", prune);
        let mut c = builder(&[1, 1, 1], prune, None);
        assert_eq!(mini_scan(&mut c, &[1, 2, -2, -1]).len(), 1);
        eprintln!("[d5] (1,1,1) N={:?}: [1,2,1] = [2,1,2], {} objects", prune, x.len());
    }
}

#[test]
#[ignore]
fn d5_against_xengine_33() {
    d5_case(&[3, 3], &[1], Some(3));
    d5_case(&[3, 3], &[1, 1], Some(3));
    d5_case(&[3, 3], &[1, 1], Some(4));
    d5_case(&[3, 3], &[1], None);
    d5_case(&[3, 3], &[1, 1], None);
    d5_case(&[3, 3], &[-1, -1], Some(4));
}

// ------------------------------------------------------------------ DecHom guards (negative)

fn target_of(sd: &mut SIndec<F>, p: Tid, r: Rung) -> (PairTarget, Vec<Rung>) {
    let mut xw = canonical_ladder(&mut sd.schur, p);
    xw.push(r);
    let lam = sd.schur.lam.clone();
    let x = sd.ar.expr(Expr::new(&lam, &xw).unwrap());
    let mult = sd.schur.term_mult(p, &[r]);
    let (q, m) = mult[0].clone();
    let k = m.c.iter().position(|&c| c != 0).unwrap();
    let lqw = canonical_ladder(&mut sd.schur, q);
    let lq = sd.ar.expr(Expr::new(&lam, &lqw).unwrap());
    (PairTarget { x, q, lq, eq: None, s: m.lo + k as i32, mult: m.c[k] as usize }, xw)
}

/// a provider that returns maps of the wrong degree (ι post-composed with a dot)
struct WrongDeg(Box<dyn HomProvider<F>>);
impl HomProvider<F> for WrongDeg {
    fn name(&self) -> &'static str {
        "WrongDeg"
    }
    fn full(&mut self, ar: &mut Arena<F>, s: &mut Schur, a: ExprId, b: ExprId, d: i32, bu: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        self.0.full(ar, s, a, b, d, bu)
    }
    fn pairs(&mut self, ar: &mut Arena<F>, s: &mut Schur, x: ExprId, q: Tid, lq: ExprId, eq: Option<MapId>, sdeg: i32, mult: usize, b: &Budget) -> Result<(Pairs<F>, SearchStats), Vec<Missing>> {
        let (mut p, st) = self.0.pairs(ar, s, x, q, lq, eq, sdeg, mult, b)?;
        let n = ar.expr_of(x).rungs.len();
        let dot = ar.add(MapNode::Dot { e: x, level: n, upright: 0, nu: vec![1] }, "dot").unwrap();
        for m in p.iotas.iter_mut() {
            *m = ar.add(MapNode::Comp(dot, *m), "wrong").unwrap();
        }
        Ok((p, st))
    }
    fn closure(&mut self, ar: &mut Arena<F>, s: &mut Schur, lp: ExprId, p: Tid, ep: Option<MapId>, outs: bool, b: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        self.0.closure(ar, s, lp, p, ep, outs, b)
    }
    fn end0(&mut self, ar: &mut Arena<F>, s: &mut Schur, l: ExprId, b: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        self.0.end0(ar, s, l, b)
    }
}

#[test]
fn dechom_negative() {
    let lam = [2u8, 2];
    let mut sd = builder(&lam, None, None);
    let mut prov = provider_for::<F>(&lam, Some(&costs()), 1.0);
    let b = Budget::default();
    let p0 = tid(&mut sd.schur, "2,0/0,2");
    let (t, _) = target_of(&mut sd, p0, Rung::e(0, 1));
    // the correct call passes
    let ok = DecHom::pairs(&mut *prov, &mut sd.ar, &mut sd.schur, &t, &b).unwrap();
    assert_eq!(ok.maps_in.len(), t.mult);
    // wrong mult
    let e = DecHom::pairs(&mut *prov, &mut sd.ar, &mut sd.schur, &PairTarget { mult: t.mult + 1, ..t.clone() }, &b).unwrap_err();
    assert!(e.guard && e.text.contains("wrong mult"), "{}", e.text);
    // wrong lq: the ladder of another type of the same block
    let mu = sd.schur.table(t.q).cols();
    let other = sd.schur.canonical_basis(&mu).order.iter().copied().find(|&r| r != t.q).unwrap();
    let ow = canonical_ladder(&mut sd.schur, other);
    let lo = sd.ar.expr(Expr::new(&lam, &ow).unwrap());
    let e = DecHom::pairs(&mut *prov, &mut sd.ar, &mut sd.schur, &PairTarget { lq: lo, ..t.clone() }, &b).unwrap_err();
    assert!(e.guard && e.text.contains("wrong lq"), "{}", e.text);
    // wrong degree: the provider's ι has degree s + 2
    let mut wd = WrongDeg(provider_for::<F>(&lam, Some(&costs()), 1.0));
    let e = DecHom::pairs(&mut wd, &mut sd.ar, &mut sd.schur, &t, &b).unwrap_err();
    assert!(e.guard && e.text.contains("wrong degree"), "{}", e.text);
    // degree scaling (D14): the true degree passes, a declared degree off by 1 or 2 fails
    let z = points::<F>(sd.ar.k(), 3, 1, 0).remove(0);
    let tab = sd.schur.table(t.q).clone();
    let cs = krhom::coset::CosetSpace::new(&lam, &tab.cols());
    let c = cs.rank(&cs.rep(&tab));
    for &m in &ok.maps_in {
        check_scaling(&sd.ar, m, c, t.s, &z).unwrap();
        assert!(check_scaling(&sd.ar, m, c, t.s + 2, &z).is_err());
        assert!(check_scaling(&sd.ar, m, c, t.s + 1, &z).is_err());
    }
}

/// a decomposable L_q (λ = (3,1,3), type 2,0,1/0,1,0/1,0,2) is refused without e_q
#[test]
fn dechom_decomposable_lq() {
    let lam = [3u8, 1, 3];
    let mut sd = builder(&lam, None, None);
    let q = tid(&mut sd.schur, "2,0,1/0,1,0/1,0,2");
    let w = canonical_ladder(&mut sd.schur, q);
    assert!(!is_indecomposable_of(&mut sd.schur, &w, q));
    assert!(matches!(sd.ty(q), Err(DecErr::NeedsIdempotent { .. })));
    let x = sd.ar.expr(Expr::new(&lam, &w).unwrap());
    let t = PairTarget { x, q, lq: x, eq: None, s: 0, mult: 1 };
    let mut prov = provider_for::<F>(&lam, Some(&costs()), 1.0);
    let e = DecHom::pairs(&mut *prov, &mut sd.ar, &mut sd.schur, &t, &Budget::default()).unwrap_err();
    assert!(e.guard && e.text.contains("decomposable"), "{}", e.text);
    let dl = decomposable_ladders(&mut sd.schur);
    assert!(dl.iter().any(|c| c.q == q));
}

// ------------------------------------------------------------------ D7 layer A, restrict

fn summand_lists(dd: &DecData<F>) -> BTreeMap<(String, String), Vec<(String, i32)>> {
    let mut out = BTreeMap::new();
    for (p, rs) in dd.dec_keys() {
        let d = dd.dec_word(p, &rs).unwrap();
        let mut v: Vec<(String, i32)> = d.summands.iter().map(|s| (format!("{}", s.q), s.s)).collect();
        v.sort();
        out.insert((format!("{}", p), krhom::smapio::rungs_str(&rs)), v);
    }
    out
}

fn build_terms(sd: &mut SIndec<F>) {
    let lam = sd.schur.lam.clone();
    let order = sd.schur.canonical_basis(&lam).order.clone();
    for p in order {
        if sd.killed(p) {
            continue;
        }
        for pos in [true, false] {
            for t in crossing_terms(lam[0], lam[1], pos, 0) {
                sd.dec_word(p, &t.rungs).unwrap();
            }
            for j in 0..crossing_terms(lam[0], lam[1], pos, 0).len() - 1 {
                let _ = sd.rickard_comp_at(p, 0, pos, j);
            }
        }
        if sd.schur.table(p).cols() == lam {
            sd.closure_gens(p).unwrap();
        }
    }
}

#[test]
fn d7_layer_a_round_trip_and_stale_stamp() {
    let lam = [3u8, 2];
    let dir = cache_dir("d7");
    let mut a = builder(&lam, None, Some(dir.clone()));
    build_terms(&mut a);
    let sa = a.stats();
    assert!(sa.layer_a_saved > 0 && sa.provider_calls > 0);
    eprintln!("[d7] {:?}: {} layer-A files, {} bytes", lam, sa.layer_a_saved, a.layer_a_bytes());
    let da = a.freeze();
    // reload: no provider call, same summand lists, same blocks at 2 points
    let mut b = builder(&lam, None, Some(dir.clone()));
    build_terms(&mut b);
    let sb = b.stats();
    assert_eq!(sb.provider_calls, 0, "layer A not used");
    assert_eq!(sb.layer_a_loaded, sa.layer_a_saved);
    let db = b.freeze();
    assert_eq!(summand_lists(&da), summand_lists(&db));
    let pts = points::<F>(da.ar().k(), 99, 2, 0);
    for (p, rs) in da.dec_keys() {
        let (x, y) = (da.dec_word(p, &rs).unwrap(), db.dec_word(p, &rs).unwrap());
        for pt in &pts {
            let (mut e1, mut e2) = (Eval::new(da.ar(), pt.clone()), Eval::new(db.ar(), pt.clone()));
            for (s1, s2) in x.summands.iter().zip(&y.summands) {
                for &c in &da.ar().flows(x.x).cosets {
                    assert_eq!(*e1.block(s1.iota, c), *e2.block(s2.iota, c));
                    assert_eq!(*e1.block(s1.pi, c), *e2.block(s2.pi, c));
                }
            }
        }
    }
    // a stale TMPL_ALG stamp is rejected (and the target recomputed)
    let la = dir.join("layer_a");
    for f in std::fs::read_dir(&la).unwrap() {
        let p = f.unwrap().path();
        let t = std::fs::read_to_string(&p).unwrap();
        std::fs::write(&p, t.replacen(&format!("tmpl_alg={}", krhom::tmpl::TMPL_ALG), "tmpl_alg=stale", 1)).unwrap();
    }
    let mut c = builder(&lam, None, Some(dir.clone()));
    build_terms(&mut c);
    let sc = c.stats();
    assert_eq!(sc.layer_a_loaded, 0);
    assert_eq!(sc.layer_a_rejected, sa.layer_a_saved);
    assert_eq!(sc.provider_calls, sa.provider_calls);
}

#[test]
fn d7_restrict() {
    let lam = [3u8, 2];
    let mut a = builder(&lam, Some(5), None);
    build_terms(&mut a);
    let da = a.freeze();
    assert!(matches!(da.restrict(6), Err(DecErr::Invalid(_))));
    for n in [3u8, 4] {
        let r = da.restrict(n).unwrap();
        assert_eq!(r.prune_n(), Some(n));
        let mut f = builder(&lam, Some(n), None);
        build_terms(&mut f);
        let df = f.freeze();
        let (lr, lf) = (summand_lists(&r), summand_lists(&df));
        // every decomposition of the fresh build is in the restricted data with the same summands
        for (k, v) in &lf {
            assert_eq!(lr.get(k), Some(v), "restrict({}) differs at {:?}", n, k);
        }
        // Rickard lists: same number of components
        for p in df.dec_keys().iter().map(|k| k.0).collect::<std::collections::BTreeSet<_>>() {
            for pos in [true, false] {
                for j in 0..2 {
                    if let Ok(x) = df.rickard_comp_at(p, 0, pos, j) {
                        assert_eq!(r.rickard_comp_at(p, 0, pos, j).unwrap().len(), x.len());
                    }
                }
            }
        }
    }
}

// ------------------------------------------------------------------ D8, D9

#[test]
fn d8_freeze_and_missing() {
    let lam = [2u8, 2];
    let mut sd = builder(&lam, None, None);
    let p0 = tid(&mut sd.schur, "2,0/0,2");
    let t = crossing_terms(2, 2, true, 0);
    let d = sd.dec_word(p0, &t[1].rungs).unwrap();
    let rc = sd.rickard_comp_at(p0, 0, true, 0).unwrap();
    let dd = sd.freeze();
    let d2 = dd.dec_word(p0, &t[1].rungs).unwrap();
    assert_eq!(d.summands, d2.summands);
    assert_eq!(*rc, *dd.rickard_comp_at(p0, 0, true, 1 - 1).unwrap());
    // a missing pair gives Err, never an empty result
    assert!(matches!(dd.dec_word(p0, &[Rung::f(0, 2)]), Err(DecErr::NotPrecomputed(_))));
    assert!(matches!(dd.rickard_comp_at(p0, 0, false, 1), Err(DecErr::NotPrecomputed(_))));
    assert!(matches!(dd.closure_gens(p0), Err(DecErr::NotPrecomputed(_))));
    assert!(matches!(dd.ty(9999), Err(DecErr::NotPrecomputed(_))));
}

fn _ss<T: Send + Sync>() {}
#[test]
fn d9_send_sync_threads() {
    _ss::<DecData<F>>();
    _ss::<Arc<DecData<F>>>();
    let lam = [2u8, 2];
    let mut sd = builder(&lam, None, None);
    let p0 = tid(&mut sd.schur, "2,0/0,2");
    let t = crossing_terms(2, 2, true, 0);
    let d = sd.dec_word(p0, &t[1].rungs).unwrap();
    let dd = sd.freeze();
    let pts = points::<F>(dd.ar().k(), 5, 2, 0);
    let seq: Vec<Vec<F>> = pts.iter().map(|pt| {
        let mut ev = Eval::new(dd.ar(), pt.clone());
        d.summands.iter().flat_map(|s| dd.ar().flows(d.x).cosets.iter().flat_map(|&c| ev.block(s.iota, c).a.clone()).collect::<Vec<_>>()).collect()
    }).collect();
    std::thread::scope(|sc| {
        let hs: Vec<_> = (0..2).map(|k| {
            let (dd, d, pt) = (dd.clone(), d.clone(), pts[k].clone());
            sc.spawn(move || {
                let mut ev = Eval::new(dd.ar(), pt);
                d.summands.iter().flat_map(|s| dd.ar().flows(d.x).cosets.iter().flat_map(|&c| ev.block(s.iota, c).a.clone()).collect::<Vec<_>>()).collect::<Vec<F>>()
            })
        }).collect();
        for (k, h) in hs.into_iter().enumerate() {
            assert_eq!(h.join().unwrap(), seq[k]);
        }
    });
}

// ------------------------------------------------------------------ (3,3) and (3,3,3) (cloud)

/// (3,3): all G1 pairs unpruned and pruned at N = 3, 4, 5; the DecData of every crossing term
/// with its Rickard components and closure generators (G2B needs the timing); layer-A size
#[test]
#[ignore]
fn d3_33_all_pairs_and_decdata_timing() {
    let lam = [3u8, 3];
    for n in [None, Some(5u8), Some(4), Some(3)] {
        let t0 = Instant::now();
        let mut sd = builder(&lam, n, None);
        let plan = pruned_plan(&mut sd.schur, n);
        let (ok, ak, kp) = all_pairs(&mut sd, &plan);
        let st = sd.stats();
        eprintln!("[33] N={:?}: {}; {} ok, {} all-killed, {} killed sources; provider {} calls {:.1}s; {:.1}s", n, plan.summary(), ok, ak, kp, st.provider_calls, st.provider_secs, t0.elapsed().as_secs_f64());
    }
    let dir = cache_dir("d33");
    let t0 = Instant::now();
    let mut sd = builder(&lam, None, Some(dir));
    build_terms(&mut sd);
    let st = sd.stats();
    let bytes = sd.layer_a_bytes();
    let ab = sd.ar.bytes().total();
    let dd = sd.freeze();
    eprintln!("[33] DecData (all crossing terms, Rickard, closure; unpruned): {:.1}s, provider {:.1}s in {} calls, sizes {:?}, layer A {} bytes, arena {} bytes, stats {:?}", t0.elapsed().as_secs_f64(), st.provider_secs, st.provider_calls, dd.sizes(), bytes, ab, st);
}

/// gate output: pruned_plan(3,3,3; N = 4, 5, 6) and idem_needed over the visited blocks
#[test]
#[ignore]
fn gate_plan_333() {
    let mut s = Schur::new(&[3, 3, 3]);
    for n in [Some(4u8), Some(5), Some(6), None] {
        let t0 = Instant::now();
        let p = pruned_plan(&mut s, n);
        eprintln!("[gate] {} ({:.1}s)", p.summary(), t0.elapsed().as_secs_f64());
        if n.is_none() {
            // R2 §1: 494 pairs (429 non-identity), 3 151 targets over the non-identity pairs
            let ni = p.pairs.iter().filter(|q| !q.identity).count();
            let tg: usize = p.pairs.iter().filter(|q| !q.identity).map(|q| q.kept.len()).sum();
            eprintln!("[gate] unpruned: {} non-identity pairs, {} targets", ni, tg);
            assert_eq!((p.pairs.len(), ni, tg), (494, 429, 3151));
        }
    }
    for c in decomposable_ladders(&mut s) {
        eprintln!("[gate] decomposable: block {:?} type {} ladder {} flows {} kill level {} char {:?}", c.mu, c.table, c.ladder, c.flows, c.kill_level, c.char_expansion);
    }
    for n in [4u8, 5, 6] {
        let v = idem_needed(&mut s, Some(n));
        eprintln!("[gate] idem_needed(3,3,3; N = {}) = {:?}", n, v.iter().map(|c| c.table.clone()).collect::<Vec<_>>());
    }
}

/// M5c route stamp (M5b open problem 4): a layer-A file written under another cost table is
/// rejected by default and accepted with `layer_a_route_stamp` (no provider call, maps
/// re-certified on load); a file of another route (chain token) is rejected in both modes.
#[test]
fn m5c_layer_a_route_stamp() {
    assert_eq!(route_header("v2 prov=homprov[hs-kl1] maxl=3 costs=hs-kl1:273recs:b8c4807fa14298b4 slowdown=1 chain=X key=k"), "v2 prov=homprov[hs-kl1] maxl=3 costs=hs-kl1 slowdown=1 chain=X key=k");
    assert_eq!(route_header("a costs=none b"), "a costs=none b");
    let lam = [3u8, 2];
    let dir = cache_dir("m5c_route");
    let mut a = builder(&lam, None, Some(dir.clone()));
    build_terms(&mut a);
    let sa = a.stats();
    assert!(sa.layer_a_saved > 0 && sa.provider_calls > 0);
    let la = dir.join("layer_a");
    let rewrite = |from: &str, to: &str| {
        for f in std::fs::read_dir(&la).unwrap() {
            let p = f.unwrap().path();
            let t = std::fs::read_to_string(&p).unwrap();
            let (h, body) = t.split_once('\n').unwrap();
            let h2: Vec<String> = h.split(' ').map(|w| if w.starts_with(from) { to.to_string() } else { w.to_string() }).collect();
            std::fs::write(&p, format!("{}\n{}", h2.join(" "), body)).unwrap();
        }
    };
    // another cost table (records and hash)
    rewrite("costs=", "costs=hs-kl1:999recs:00000000deadbeef");
    let mut b = builder(&lam, None, Some(dir.clone()));
    build_terms(&mut b);
    let sb = b.stats();
    assert_eq!((sb.layer_a_loaded, sb.layer_a_rejected), (0, sa.layer_a_saved), "default: exact stamp only");
    // b recomputed and rewrote every file with the current stamp; make them foreign again
    rewrite("costs=", "costs=hs-kl1:999recs:00000000deadbeef");
    let (prov, stamp) = provider_for_stamped::<F>(&lam, Some(&costs()), 1.0);
    let mut c = SIndec::<F>::with_stamp(&lam, prov, stamp, SIndecOpts { cache_dir: Some(dir.clone()), layer_a_route_stamp: true, ..Default::default() });
    build_terms(&mut c);
    let sc = c.stats();
    assert_eq!(sc.provider_calls, 0, "route stamp: layer A used");
    assert_eq!((sc.layer_a_loaded, sc.layer_a_route, sc.layer_a_rejected), (sa.layer_a_saved, sa.layer_a_saved, 0));
    // another route (chain) is rejected even with the route stamp
    rewrite("chain=", "chain=OTHER");
    let (prov, stamp) = provider_for_stamped::<F>(&lam, Some(&costs()), 1.0);
    let mut d = SIndec::<F>::with_stamp(&lam, prov, stamp, SIndecOpts { cache_dir: Some(dir.clone()), layer_a_route_stamp: true, ..Default::default() });
    build_terms(&mut d);
    let sd = d.stats();
    assert_eq!((sd.layer_a_loaded, sd.layer_a_rejected), (0, sa.layer_a_saved));
}

/// M5c: a leaf isomorphism whose only windows are refused by admission is reported as an
/// engineering (Refused) failure, not as Exhausted (level-5 pair 322 of (3,3,3): its leaf
/// [E0 E1(2) E0 F1 F0(2) F1] needs the (4,4,1) window iso with 5.4 M unknowns; M5b counted the
/// pair as Exhausted at 2x and sent it to the M4c list)
#[test]
fn m5c_leaf_iso_refusal_is_engineering() {
    let lam = [3u8, 3, 3];
    let mut s = Schur::new(&lam);
    let plan = pruned_plan(&mut s, Some(5));
    let pp = &plan.pairs[322];
    let (p, r) = (s.table(pp.p).fmt(), pp.rung);
    assert_eq!((p.as_str(), krhom::smapio::rung_str(&r).as_str()), ("1,2,0/1,1,1/1,1,1", "F1.1"));
    // the M5b provider (cost table v3; the table is stamped for the M5b host class)
    unsafe { std::env::set_var("KRHOM_COSTS_ANY_HOST", "1") };
    let costs = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/homprov/costs_333.txt");
    let (prov, stamp) = provider_for_stamped::<F>(&lam, Some(&costs), 1.0);
    assert!(stamp.contains("273recs"), "{}", stamp);
    let mut sd = SIndec::<F>::with_stamp(&lam, prov, stamp, SIndecOpts { prune_n: Some(5), ..Default::default() });
    let pt = tid(&mut sd.schur, &p);
    match sd.dec_rung(pt, r) {
        Err(DecErr::Provider { engineering, text }) => {
            assert!(engineering, "leaf iso admission refusal reported as exhausted: {}", text);
            assert!(text.contains("Refused"), "{}", text);
        }
        Err(e) => panic!("unexpected error {}", e),
        Ok(_) => panic!("pair 322 decomposed without a cost table"),
    }
}
