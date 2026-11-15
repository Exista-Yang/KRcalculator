//! M4c B1: closed-form generators (src/leafgen.rs) against exact `sbim::hom_space` windows,
//! `certify_added` on every template, switchback ψ∘φ = id, braid composites (CORE Lemma 3.12)
//! and far commutation (Lemma 3.13).
//!
//! Default suite: small blocks and cells ((1,1,1), (1,2,1), (2,1,2), (2,2,2)). Ignored (run once,
//! reported in docs/colored_design/m4c/B1_report.md): `certify_all_to_333` (every template up to
//! blocks (3,3,3)) and `braid_large`.
use krhom::coset::CosetSpace;
use krhom::dense::rank;
use krhom::field::{Field, F31};
use krhom::homprov::hom_unknowns;
use krhom::leafgen::*;
use krhom::schur::Table;
use krhom::sexpr::{Expr, FlowIndex};
use krhom::smap::*;
use krhom::tmpl::*;
use krhom::web::Rung;
use std::time::Instant;

type F = F31;

const EXACT_CAP: u128 = 100_000;

fn pt(k: usize, seed: u64) -> Vec<F> {
    points::<F>(k, seed, 1, 3)[0].x.clone()
}

fn kw(w: &Window) -> usize {
    w.labels.iter().map(|&a| a as usize).sum()
}

fn flat(w: &WinLoc<F>) -> Vec<F> {
    w.dense().into_iter().flatten().collect()
}

/// (dimension of the exact Hom^deg window, closed form in its span); None above EXACT_CAP unknowns
fn exact_check(lib: &mut TmplLib<F>, cf: &ClosedForm, t: TmplId) -> Option<(usize, bool)> {
    let u = hom_unknowns(&cf.src.labels, &cf.src.rungs, &cf.tgt.rungs, cf.deg);
    if u > EXACT_CAP {
        return None;
    }
    let ids = lib.try_hom_basis_cap(&cf.src, &cf.tgt, cf.deg, u128::MAX).expect("exact window");
    let ys = [pt(kw(&cf.src), 11), pt(kw(&cf.src), 12)];
    let vec_of = |lib: &TmplLib<F>, id: TmplId| -> Vec<F> { ys.iter().flat_map(|y| flat(&lib.localize(id, y))).collect() };
    let mut rows: Vec<Vec<F>> = ids.iter().map(|&i| vec_of(lib, i)).collect();
    let n = rows.first().map(|r| r.len()).unwrap_or(0);
    let r0 = rank(&rows, n);
    assert_eq!(r0, ids.len(), "{}: exact basis not independent after localization", cf.name);
    let c = vec_of(lib, t);
    let n = c.len();
    rows.push(c);
    let r1 = rank(&rows, n);
    Some((ids.len(), r1 == r0))
}

fn nonzero(lib: &TmplLib<F>, t: TmplId, k: usize) -> bool {
    lib.localize(t, &pt(k, 5)).vals.iter().any(|v| !v.is_zero())
}

/// build + certify + exact comparison; returns the template
fn check(lib: &mut TmplLib<F>, cf: &ClosedForm, seed: u64, want_deg: Option<i32>) -> TmplId {
    if let Some(d) = want_deg {
        assert_eq!(cf.deg, d, "{}: degree", cf.name);
    }
    let b = build(lib, cf, seed).unwrap_or_else(|e| panic!("{}", e));
    assert!(lib.is_certified(b.tmpl));
    assert!(nonzero(lib, b.tmpl, kw(&cf.src)), "{}: zero map", cf.name);
    if let Some((dim, ok)) = exact_check(lib, cf, b.tmpl) {
        assert!(ok, "{}: not in the exact Hom^{} window (dim {})", cf.name, cf.deg, dim);
        eprintln!("{}: deg {} flows {:?} exact dim {} in span", cf.name, cf.deg, b.flows, dim);
    } else {
        eprintln!("{}: deg {} flows {:?} (exact window above cap)", cf.name, cf.deg, b.flows);
    }
    b.tmpl
}

fn mul(a: &[Vec<F>], b: &[Vec<F>]) -> Vec<Vec<F>> {
    krhom::local::matmul(a, b)
}

fn is_identity(m: &[Vec<F>]) -> bool {
    m.iter().enumerate().all(|(i, r)| r.iter().enumerate().all(|(j, v)| if i == j { v.is_one() } else { v.is_zero() }))
}

// ------------------------------------------------------------------ R1: the hand derivations

/// R1 first test: the sideways crossing on blocks (1,1,1) and (2,1,2), both directions, against
/// the exact windows; degree ac; right∘left = multiplication by e(C, A∖T) (diagonal).
#[test]
fn r1_sideways_crossing_111_212() {
    let mut lib: TmplLib<F> = TmplLib::new(None);
    for &(a, b, c) in &[(1u8, 1u8, 1u8), (2, 1, 2), (1, 2, 1), (2, 2, 2)] {
        let (r, l) = sideways_abc(a, b, c).unwrap();
        assert_eq!(r.src.rungs, vec![Rung::e(0, c as usize), Rung::f(0, (b + c) as usize)]);
        assert_eq!(l.src.rungs, vec![Rung::f(0, b as usize)]);
        let ac = a as i32 * c as i32;
        let tr = check(&mut lib, &r, 100, Some(ac));
        let tl = check(&mut lib, &l, 101, Some(ac));
        // right ∘ left on Y = [F0.b]: diagonal e(C, A∖T)
        let y = pt((a + b + c) as usize, 7);
        let comp = mul(&lib.localize(tr, &y).dense(), &lib.localize(tl, &y).dense());
        let fi = FlowIndex::build(&l.src.expr());
        let cm: u64 = ((1u64 << c) - 1) << (a + b);
        for (i, row) in comp.iter().enumerate() {
            for (j, v) in row.iter().enumerate() {
                if i == j {
                    let t = fi.flow(i as u32)[0];
                    let rest = ((1u64 << (a + b)) - 1) & !t;
                    assert_eq!(*v, e_pair::<F>(cm, rest, &y), "right∘left diagonal at flow {}", i);
                } else {
                    assert!(v.is_zero());
                }
            }
        }
    }
}

/// R1 second test: the EKo (108) switchback E0.m F0.(N−b) ⇄ E0.b F0.a on (a, m): φ = Pull (all
/// ones), ψ = normalized push; ψ∘φ = id; both certified; both in the exact Hom^0 windows.
#[test]
fn r1_switchback() {
    let mut lib: TmplLib<F> = TmplLib::new(None);
    for &(a, m, b) in &[(1u8, 2u8, 1u8), (1, 3, 1), (1, 3, 2), (2, 2, 1), (2, 3, 1), (2, 3, 2), (2, 2, 2)] {
        let (phi, psi) = switchback_eko(a, m, b).unwrap();
        let n = (a + m) as usize;
        assert_eq!(phi.src.rungs, canonical_max(&[a, m], &[b, a + m - b]));
        let tp = check(&mut lib, &phi, 200, Some(0));
        let ts = check(&mut lib, &psi, 201, Some(0));
        let y = pt(n, 9);
        let c = mul(&lib.localize(ts, &y).dense(), &lib.localize(tp, &y).dense());
        assert!(is_identity(&c), "switchback ({},{},{}): psi∘phi != id", a, m, b);
        // φ∘ψ is an idempotent on Y, not the identity when Y has more flows than X
        let (ns, nt) = psi.flows();
        let e = mul(&lib.localize(tp, &y).dense(), &lib.localize(ts, &y).dense());
        assert_eq!(mul(&e, &e), e, "phi∘psi not idempotent");
        assert_eq!(is_identity(&e), ns == nt);
    }
}

// ------------------------------------------------------------------ the generator library

/// E2/E3 elementary maps are instances: Pull/Push reproduce Eps, Eta, EpsL, EtaL, Iota, Trace.
#[test]
fn elem_kinds_are_instances() {
    let y: Vec<F> = pt(5, 3);
    for &(a, b) in &[(2u8, 3u8), (3, 2), (1, 4)] {
        let l = [a, b];
        let w = Window::new(&l, &[Rung::e(0, b as usize), Rung::f(0, b as usize)]);
        let e = Window::new(&l, &[]);
        let pull = ClosedForm::new("eps", CfKind::Pull, w.clone(), e.clone()).unwrap();
        let push = ClosedForm::new("eta", CfKind::Push, e.clone(), w.clone()).unwrap();
        let (wf, ef) = (FlowIndex::build(&w.expr()), FlowIndex::build(&e.expr()));
        assert_eq!(pull.winloc(&y), elem_winloc(ElemKind::Eps, &wf, &ef, &y));
        assert_eq!(push.winloc(&y), elem_winloc(ElemKind::Eta, &ef, &wf, &y));
        assert_eq!(pull.deg, a as i32 * b as i32);
        assert_eq!(push.deg, a as i32 * b as i32);
        let wl = Window::new(&l, &[Rung::f(0, a as usize), Rung::e(0, a as usize)]);
        let wlf = FlowIndex::build(&wl.expr());
        let pl = ClosedForm::new("epsL", CfKind::Pull, wl.clone(), e.clone()).unwrap();
        let ql = ClosedForm::new("etaL", CfKind::Push, e.clone(), wl.clone()).unwrap();
        assert_eq!(pl.winloc(&y), elem_winloc(ElemKind::EpsL, &wlf, &ef, &y));
        assert_eq!(ql.winloc(&y), elem_winloc(ElemKind::EtaL, &ef, &wlf, &y));
    }
    for &(c, b) in &[(3u8, 1u8), (4, 2), (5, 3)] {
        let l = [c, 0];
        let d = Window::new(&l, &[Rung::f(0, b as usize), Rung::e(0, b as usize)]);
        let e = Window::new(&l, &[]);
        let (df, ef) = (FlowIndex::build(&d.expr()), FlowIndex::build(&e.expr()));
        let iota = ClosedForm::new("iota", CfKind::Pull, e.clone(), d.clone()).unwrap();
        let trace = ClosedForm::new("trace", CfKind::Push, d.clone(), e.clone()).unwrap();
        assert_eq!(iota.winloc(&y), elem_winloc(ElemKind::Iota, &ef, &df, &y));
        assert_eq!(trace.winloc(&y), elem_winloc(ElemKind::Trace, &df, &ef, &y));
        let ab = (b as i32) * (c as i32 - b as i32);
        assert_eq!((iota.deg, trace.deg), (-ab, -ab));
        // the general cap/cup on (c, 0) coincide with Trace / Iota∘(·e(S, ∅) = 1)
        let (cap, cup) = cup_cap(&l, Rung::f(0, b as usize)).unwrap();
        assert_eq!(cap.winloc(&y), trace.winloc(&y));
        assert_eq!(cup.winloc(&y), iota.winloc(&y));
    }
}

/// General cw/ccw cups and caps E^(r)F^(r) ⇄ 1 (r < block) at small labels.
#[test]
fn cups_and_caps_small() {
    let mut lib: TmplLib<F> = TmplLib::new(None);
    for &(a, b) in &[(1u8, 1u8), (1, 2), (2, 1), (2, 2), (1, 3), (3, 1), (2, 3)] {
        for r in 1..=b {
            let (cap, cup) = cup_cap(&[a, b], Rung::e(0, r as usize)).unwrap();
            let want = r as i32 * a as i32 - r as i32 * (b as i32 - r as i32);
            let tc = check(&mut lib, &cap, 300 + r as u64, Some(want));
            let tu = check(&mut lib, &cup, 310 + r as u64, Some(want));
            // cap ∘ cup on [] is multiplication by a polynomial of degree 2·want
            let y = pt((a + b) as usize, 4);
            let c = mul(&lib.localize(tc, &y).dense(), &lib.localize(tu, &y).dense());
            // (degree 2·want: a negative-degree endomorphism of R^λ vanishes)
            assert_eq!(c[0][0].is_zero(), want < 0, "cap∘cup at ({},{}) r = {}, degree {}", a, b, r, 2 * want);
        }
        for r in 1..=a {
            let (cap, cup) = cup_cap(&[a, b], Rung::f(0, r as usize)).unwrap();
            check(&mut lib, &cap, 320 + r as u64, None);
            check(&mut lib, &cup, 330 + r as u64, None);
        }
    }
}

/// Mirror (E-type) sideways crossings and switchbacks on other 2-upright words.
#[test]
fn mirror_crossings_and_switchbacks_small() {
    let mut lib: TmplLib<F> = TmplLib::new(None);
    for &(a, b, c) in &[(1u8, 1u8, 1u8), (1, 2, 1), (2, 1, 2), (2, 2, 1)] {
        let (r, l) = sideways(&[a, b + c], Rung::e(0, b as usize)).unwrap();
        check(&mut lib, &r, 400, Some(a as i32 * c as i32));
        check(&mut lib, &l, 401, Some(a as i32 * c as i32));
    }
    // mirrored EKo switchback F0.b E0.m' on (m, a): reduced for the maximum
    for &(m, a, b) in &[(2u8, 1u8, 1u8), (3, 2, 1), (2, 2, 1)] {
        let (phi, psi) = switchback(&[m, a], &[Rung::f(0, b as usize), Rung::e(0, a as usize)]).unwrap();
        check(&mut lib, &phi, 410, Some(0));
        check(&mut lib, &psi, 411, Some(0));
    }
    // a non-reduced word is refused by `switchback` (degree != 0)
    assert!(switchback(&[2, 2], &[Rung::e(0, 1), Rung::f(0, 1), Rung::e(0, 1), Rung::f(0, 1)]).is_err());
}

// ------------------------------------------------------------------ Assoc validator

#[test]
fn assoc_validator() {
    let lam = [2u8, 1, 1];
    let mut ar: Arena<F> = Arena::new(&lam, TmplLib::new(None));
    // up-up: E1.1 E0.2 ↔ E0.1 E0.1? no: merges in two orders through a zero upright
    let a = Expr::new(&lam, &[Rung::e(1, 1), Rung::e(0, 2)]).unwrap();
    let b = Expr::new(&lam, &[Rung::e(0, 1), Rung::e(1, 1), Rung::e(0, 1)]).unwrap();
    assert_eq!(multistep(&a.bottom, &a.rungs).unwrap(), multistep(&b.bottom, &b.rungs).unwrap());
    validate_assoc(&mut ar, &a, &b, 1).expect("up-up Assoc");
    assert_eq!(ar.assoc_certified(&a, &b), Some(F::one()));
    // down-down on (4,0,0): split as (1 | 3) then (1 | 2) vs (2 | 2) then (1 | 1)
    let l4 = [4u8, 0, 0];
    let c = Expr::new(&l4, &[Rung::f(0, 3), Rung::f(1, 2)]).unwrap();
    let d = Expr::new(&l4, &[Rung::f(0, 2), Rung::f(1, 2), Rung::f(0, 1)]).unwrap();
    assert_eq!(c.top(), d.top());
    let mut ar4: Arena<F> = Arena::new(&l4, TmplLib::new(None));
    validate_assoc(&mut ar4, &c, &d, 2).expect("down-down Assoc");
    // a genuine rex move (switchback) is not an Assoc
    let e1 = Expr::new(&[1, 2], &[Rung::e(0, 2), Rung::f(0, 2)]).unwrap();
    let e2 = Expr::new(&[1, 2], &[Rung::e(0, 1), Rung::f(0, 1)]).unwrap();
    let mut ar2: Arena<F> = Arena::new(&[1, 2], TmplLib::new(None));
    assert!(validate_assoc(&mut ar2, &e1, &e2, 3).is_err());
    assert_eq!(ar2.assoc_certified(&e1, &e2), None);
}

// ------------------------------------------------------------------ braid and far commutation

fn top_coset(lam: &[u8], top: &[u8], anti: bool) -> u64 {
    let n = lam.len();
    let cs = CosetSpace::new(lam, top);
    let mut e = vec![0u8; n * n];
    for i in 0..n {
        e[i * n + if anti { n - 1 - i } else { i }] = lam[i];
    }
    let t = Table { nr: n as u8, nc: n as u8, e };
    cs.rank(&cs.rep(&t))
}

/// The atomic braid composite A → B: certified (genuine) DAG, degree 0, 1⊗ ↦ 1⊗ at every coset,
/// 1×1 at the top coset with value 1; in the exact Hom^0(A, B) window where affordable; equal to
/// the direct Rex closed form on the whole window; B → A → B top scalar 1.
fn braid_case(x: u8, y: u8, z: u8, exact: bool, direct_max: usize) {
    let t0 = Instant::now();
    let bw = braid_words(x, y, z).unwrap();
    let lam = vec![x, y, z];
    let k = (x + y + z) as usize;
    let mut ar: Arena<F> = Arena::new(&lam, TmplLib::new(None));
    let (m, log) = braid_dag(&mut ar, &bw, 500).unwrap_or_else(|e| panic!("braid {:?}: {}", (x, y, z), e));
    for l in &log {
        eprintln!("braid {:?}: {}", (x, y, z), l);
    }
    ar.assert_genuine(m).expect("braid composite genuine");
    assert_eq!(ar.info(m).deg, 0);
    // the reverse braid B -> A (roles of the cells: a2 a1 a2 on (x,y,z) = a1 a2 a1 mirrored); built
    // as the direct Rex closed form when small, to test the 1×1 form at the top coset
    let p = points::<F>(k, 77, 1, 0).pop().unwrap();
    let (ea, eb) = (ar.info(m).src, ar.info(m).tgt);
    let mut ev = Eval::new(&ar, p.clone());
    let lay = ev.layout(ea, eb);
    for &c in &lay.cosets {
        let bl = ev.block(m, c);
        let ones = Block { rows: bl.cols, cols: 1, a: vec![F::one(); bl.cols as usize] };
        let img = bl.mul(&ones);
        assert!(img.a.iter().all(|v| v.is_one()), "braid: 1⊗ not mapped to 1⊗ at coset {}", c);
    }
    let top = top_coset(&lam, &[z, y, x], true);
    let bt = ev.block(m, top);
    assert_eq!((bt.rows, bt.cols), (1, 1), "top coset block not 1×1");
    assert!(bt.a[0].is_one());
    let nfa = FlowIndex::build(&Expr::new(&lam, &bw.a).unwrap()).nflows();
    eprintln!("braid {:?}: A {} flows, {} cosets, composite ok ({:.1} s)", (x, y, z), nfa, lay.cosets.len(), t0.elapsed().as_secs_f64());
    if nfa <= direct_max {
        // the direct closed form on the whole window agrees with the composite everywhere
        let wa = Window::new(&lam, &bw.a);
        let wb = Window::new(&lam, &bw.b);
        let direct = ClosedForm::new("braid direct", CfKind::Rex, wa.clone(), wb.clone()).unwrap();
        let mut lib: TmplLib<F> = TmplLib::new(None);
        let t = build(&mut lib, &direct, 501).unwrap_or_else(|e| panic!("{}", e)).tmpl;
        let wl = lib.localize(t, &p.x);
        let fa = FlowIndex::build(&wa.expr());
        let fb = FlowIndex::build(&wb.expr());
        let d = wl.dense();
        for &c in &lay.cosets {
            let bl = ev.block(m, c);
            let (ra, rb) = (fa.range(c), fb.range(c));
            for (ii, i) in rb.clone().enumerate() {
                for (jj, j) in ra.clone().enumerate() {
                    assert_eq!(bl.get(ii, jj), d[i as usize][j as usize], "composite != direct Rex at coset {}", c);
                }
            }
        }
        if exact {
            if let Some((dim, ok)) = exact_check(&mut lib, &direct, t) {
                assert!(ok, "braid {:?}: not in the exact Hom^0 window", (x, y, z));
                eprintln!("braid {:?}: exact Hom^0(A,B) dim {}, composite in span", (x, y, z), dim);
            }
        }
        // B -> A direct; B->A ∘ A->B at the top coset = 1 (rank of the 1×1 form)
        let back = ClosedForm::new("braid back", CfKind::Rex, wb, wa).unwrap();
        let tb = build(&mut lib, &back, 502).unwrap().tmpl;
        let prod = mul(&lib.localize(tb, &p.x).dense(), &d);
        let ia = fa.range(top).start as usize;
        assert!(prod[ia][ia].is_one(), "back∘braid top scalar != 1");
    }
}

#[test]
fn braid_small_cells() {
    braid_case(1, 1, 1, true, 400);
    braid_case(1, 2, 1, true, 400);
    braid_case(2, 1, 2, true, 400);
}

#[test]
fn braid_222() {
    // composite checks only (the direct 216-flow Rex comparison is in `braid_large`)
    braid_case(2, 2, 2, false, 0);
}

#[test]
#[ignore]
fn braid_large() {
    braid_case(2, 2, 2, false, 400);
    braid_case(3, 1, 3, false, 0);
    braid_case(2, 3, 2, false, 0);
}

const FAR_CELLS: &[(u8, u8, u8, u8)] = &[(1, 1, 1, 1), (1, 2, 1, 1), (2, 1, 1, 2)];

/// Far commutation at cells (3,3,2,1) (C_m4c C1: 192 484 exact unknowns): certify_added only.
#[test]
#[ignore]
fn far_commutation_3321() {
    let mut lib: TmplLib<F> = TmplLib::new(None);
    let t0 = Instant::now();
    let cf = far_commute(&[3, 3, 2, 1], &atomic(0, 3, 3), &atomic(2, 2, 1)).unwrap();
    let b = build(&mut lib, &cf, 601).unwrap_or_else(|e| panic!("{}", e));
    println!("CERT {} flows {:?} interp {:.1}s cert {:.1}s (pts {}, monos {}) total {:.1}s", cf.name, b.flows, b.interp_secs, b.cert.secs, b.cert.points, b.cert.max_monos, t0.elapsed().as_secs_f64());
}

/// Far commutation of atomic crossings on disjoint uprights: certified, degree 0, 1⊗ ↦ 1⊗, in
/// the exact window when affordable.
#[test]
fn far_commutation() {
    let mut lib: TmplLib<F> = TmplLib::new(None);
    for &(x, y, z, w) in FAR_CELLS {
        let l = [x, y, z, w];
        let a = atomic(0, x, y);
        let b = atomic(2, z, w);
        let cf = far_commute(&l, &a, &b).unwrap();
        let t = check(&mut lib, &cf, 600, Some(0));
        let wl = lib.localize(t, &pt((x + y + z + w) as usize, 2));
        assert!(wl.vals.iter().all(|v| v.is_one()));
    }
}

// ------------------------------------------------------------------ certify_added to (3,3,3)

/// Every template family up to blocks (3,3,3): sideways crossings (F- and E-type, all (a,b,c) ≤ 3),
/// EKo switchbacks on (a, m) with a + m ≤ 9, a ≤ 3, and cups/caps on labels ≤ (3, 6) / (6, 3).
#[test]
#[ignore]
fn certify_all_to_333() {
    let mut lib: TmplLib<F> = TmplLib::new(None);
    let mut n = 0;
    let t0 = Instant::now();
    let mut go = |lib: &mut TmplLib<F>, cf: ClosedForm| {
        let t1 = Instant::now();
        let b = build(lib, &cf, 900 + n).unwrap_or_else(|e| panic!("{}", e));
        n += 1;
        println!("CERT {} deg {} flows {:?} interp {:.2}s cert {:.2}s (pts {}, monos {}) total {:.2}s", cf.name, b.deg, b.flows, b.interp_secs, b.cert.secs, b.cert.points, b.cert.max_monos, t1.elapsed().as_secs_f64());
    };
    for a in 1..=3u8 {
        for b in 1..=3u8 {
            for c in 1..=3u8 {
                let (r, l) = sideways_abc(a, b, c).unwrap();
                go(&mut lib, r);
                go(&mut lib, l);
                let (r, l) = sideways(&[a, b + c], Rung::e(0, b as usize)).unwrap();
                go(&mut lib, r);
                go(&mut lib, l);
            }
        }
    }
    for a in 1..=3u8 {
        for m in 1..=(9 - a) {
            for b in 1..m {
                if let Ok((phi, psi)) = switchback_eko(a, m, b) {
                    go(&mut lib, phi);
                    go(&mut lib, psi);
                }
            }
        }
    }
    for &(a, b) in &[(3u8, 3u8), (3, 6), (6, 3)] {
        for r in 1..=b.min(3) {
            let (cap, cup) = cup_cap(&[a, b], Rung::e(0, r as usize)).unwrap();
            go(&mut lib, cap);
            go(&mut lib, cup);
        }
    }
    println!("CERT done in {:.1} s", t0.elapsed().as_secs_f64());
}

/// The switchback pieces of the braid route at large cells, certified through `certify_added`;
/// the middle up-up/down-down step M → M' is reported (Assoc validator: refused above 500 flows).
#[test]
#[ignore]
fn braid_pieces_large() {
    for &(x, y, z) in &[(3u8, 1u8, 3u8), (2, 3, 2), (3, 2, 3), (3, 3, 3)] {
        let bw = braid_words(x, y, z).unwrap();
        let lam = vec![x, y, z];
        let mut lib: TmplLib<F> = TmplLib::new(None);
        let (_, la, long_a, short_a) = bw.win_a.clone();
        let (_, lb, long_b, short_b) = bw.win_b.clone();
        let psi = ClosedForm::new("psi_A", CfKind::Rex, Window::new(&la, &long_a), Window::new(&la, &short_a)).unwrap();
        let phi = ClosedForm::new("phi_B", CfKind::Pull, Window::new(&lb, &short_b), Window::new(&lb, &long_b)).unwrap();
        for cf in [psi, phi] {
            let t0 = Instant::now();
            let b = build(&mut lib, &cf, 700).unwrap_or_else(|e| panic!("{:?} {}: {}", (x, y, z), cf.name, e));
            println!("CERT braid {:?} {} flows {:?} interp {:.1}s cert {:.1}s (monos {}) total {:.1}s", (x, y, z), cf.name, b.flows, b.interp_secs, b.cert.secs, b.cert.max_monos, t0.elapsed().as_secs_f64());
        }
        let (me, mpe) = (Expr::new(&lam, &bw.m).unwrap(), Expr::new(&lam, &bw.mp).unwrap());
        let mut ar: Arena<F> = Arena::new(&lam, TmplLib::new(None));
        let fa = FlowIndex::build(&Expr::new(&lam, &bw.a).unwrap()).nflows();
        let r = validate_assoc(&mut ar, &me, &mpe, 1);
        println!("CERT braid {:?}: A {} flows, M {} flows, M->M' Assoc: {:?}", (x, y, z), fa, FlowIndex::build(&me).nflows(), r.map(|c| c.secs));
    }
}

/// The two switchback pieces of the (3,3,3) braid (400 × 84 flows) through `certify_added`.
#[test]
#[ignore]
fn braid_pieces_333() {
    let bw = braid_words(3, 3, 3).unwrap();
    let mut lib: TmplLib<F> = TmplLib::new(None);
    let (_, la, long_a, short_a) = bw.win_a.clone();
    let (_, lb, long_b, short_b) = bw.win_b.clone();
    let psi = ClosedForm::new("psi_A", CfKind::Rex, Window::new(&la, &long_a), Window::new(&la, &short_a)).unwrap();
    let phi = ClosedForm::new("phi_B", CfKind::Pull, Window::new(&lb, &short_b), Window::new(&lb, &long_b)).unwrap();
    for cf in [psi, phi] {
        let t0 = Instant::now();
        let b = build(&mut lib, &cf, 710).unwrap_or_else(|e| panic!("{}: {}", cf.name, e));
        println!("CERT braid333 {} flows {:?} interp {:.1}s cert {:.1}s (monos {}) total {:.1}s", cf.name, b.flows, b.interp_secs, b.cert.secs, b.cert.max_monos, t0.elapsed().as_secs_f64());
    }
}

/// Switchback and cup/cap families of `certify_all_to_333` alone, smallest windows first
/// (separate process: the 9-root windows wait on one-time GrTable(r, 9) builds).
#[test]
#[ignore]
fn certify_switchbacks_and_cups() {
    let mut lib: TmplLib<F> = TmplLib::new(None);
    let mut list: Vec<ClosedForm> = vec![];
    for n in 3..=9u8 {
        for a in 1..=3u8.min(n - 1) {
            let m = n - a;
            for b in 1..m {
                if let Ok((phi, psi)) = switchback_eko(a, m, b) {
                    list.push(phi);
                    list.push(psi);
                }
            }
        }
        for &(a, b) in &[(3u8, 3u8), (3, 6), (6, 3)] {
            if a + b == n {
                for r in 1..=b.min(3) {
                    let (cap, cup) = cup_cap(&[a, b], Rung::e(0, r as usize)).unwrap();
                    list.push(cap);
                    list.push(cup);
                }
            }
        }
    }
    for (i, cf) in list.iter().enumerate() {
        let t1 = Instant::now();
        let b = build(&mut lib, cf, 2000 + i as u64).unwrap_or_else(|e| panic!("{}", e));
        println!("CERT {} deg {} flows {:?} interp {:.2}s cert {:.2}s (pts {}, monos {}) total {:.2}s", cf.name, b.deg, b.flows, b.interp_secs, b.cert.secs, b.cert.points, b.cert.max_monos, t1.elapsed().as_secs_f64());
    }
}

/// Review wave 4, finding 1: the interpolation / certificate memory gate is on monomials and
/// stored entries, computed from the words before any SBim is built. braid_dag at cells (3,2,3),
/// (3,3,2), (2,3,3) (middle step 560 flows, 39 148 monomials at the top entry degree, about 98 GB
/// of localization matrices) must return Err immediately instead of allocating.
#[test]
fn interp_gate_refuses_braid_fallback_323() {
    for &(x, y, z) in &[(3u8, 2u8, 3u8), (3, 3, 2), (2, 3, 3)] {
        let bw = braid_words(x, y, z).unwrap();
        let lam = vec![x, y, z];
        let (wm, wmp) = (Window::new(&lam, &bw.m), Window::new(&lam, &bw.mp));
        let fl = Expr::new(&lam, &bw.m).unwrap().nflows();
        let e = interp_gate(&wm, &wmp, 0).unwrap_err();
        eprintln!("gate {:?}: {} flows: {}", (x, y, z), fl, e);
        assert!(e.contains("monomials") || e.contains("entries"), "{}", e);
        let t0 = Instant::now();
        let mut ar: Arena<F> = Arena::new(&lam, TmplLib::new(None));
        let r = braid_dag(&mut ar, &bw, 500);
        assert!(r.is_err(), "braid_dag {:?} must refuse", (x, y, z));
        assert!(t0.elapsed().as_secs_f64() < 30.0, "braid_dag {:?} refusal took {:.1} s", (x, y, z), t0.elapsed().as_secs_f64());
        // interpolate and validate_assoc refuse through the same gate, before any SBim build
        let me = Expr::new(&lam, &bw.m).unwrap();
        let mpe = Expr::new(&lam, &bw.mp).unwrap();
        let mut f = |_: &[F]| -> WinLoc<F> { panic!("closed form must not be evaluated") };
        assert!(interpolate::<F>(&wm, &wmp, 0, &mut f, 1).is_err());
        assert!(validate_assoc(&mut ar, &me, &mpe, 1).is_err());
    }
    // the (3,2,3) count matches the review's generating-function count (39 148 at degree 42)
    assert_eq!(krhom::rexplan::monos_of_deg(&[3, 2, 3], 42), 39_148);
    // small windows pass the gate and still build (r1 shapes)
    let (a, b) = sideways_abc(1, 2, 1).unwrap();
    for cf in [&a, &b] {
        let c = interp_gate(&cf.src, &cf.tgt, cf.deg).unwrap();
        assert!(c.max_monos <= MAX_INTERP_MONOS && c.cells <= MAX_INTERP_CELLS);
        let mut lib: TmplLib<F> = TmplLib::new(None);
        build(&mut lib, cf, 9).unwrap();
    }
}
