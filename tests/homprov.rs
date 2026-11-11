//! M4b acceptance tests B1-B8 (SPEC_M4_M6a §4.6).
use krhom::field::{Field, F31};
use krhom::homprov::*;
use krhom::schur::{Schur, Tid};
use krhom::sexpr::Expr;
use krhom::smap::*;
use krhom::tmpl::*;
use krhom::web::{Dir, Rung};
use std::path::PathBuf;
use std::time::Instant;

type F = F31;

fn e(i: usize, r: usize) -> Rung {
    Rung::e(i, r)
}
fn f(i: usize, r: usize) -> Rung {
    Rung::f(i, r)
}
fn ex(b: &[u8], r: &[Rung]) -> Expr {
    Expr::new(b, r).unwrap_or_else(|x| panic!("{}", x))
}
/// a persistent per-test template cache (cleared at the start of each test)
fn cache_dir(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("m4b_cache").join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}
fn canonical(s: &mut Schur, mu: &[u8], p: Tid) -> Vec<Rung> {
    s.canonical_basis(mu).ladder_minflow[&p].clone()
}

// ------------------------------------------------------------------ B1

/// generator degrees of Hom(R^λ, L_p) and Hom(L_p, R^λ) from exact Nakayama (rlambda_bases)
fn nakayama_gens(lam: &[u8], w: &[Rung]) -> (Vec<i32>, Vec<i32>) {
    let x = ex(lam, w);
    let b = x.sbim::<F>();
    let basis = krhom::cclose::rlambda_bases(&b).unwrap();
    let sh = b.shift;
    let mut ins: Vec<i32> = basis.ins.iter().map(|g| g.0 + sh).collect();
    let mut outs: Vec<i32> = basis.outs.iter().map(|g| g.0 - sh).collect();
    ins.sort();
    outs.sort();
    (ins, outs)
}
fn gens_list(g: &std::collections::BTreeMap<i32, i128>) -> Vec<i32> {
    let mut v = vec![];
    for (&d, &c) in g {
        for _ in 0..c {
            v.push(d);
        }
    }
    v
}

#[test]
fn b1_hom_gens_equal_nakayama() {
    let mut n = 0;
    for (lam, maxf) in [(vec![2u8, 2], 1000u128), (vec![3, 3], 20), (vec![2, 2, 2], 18)] {
        let mut s = Schur::new(&lam);
        let order = s.canonical_basis(&lam).order.clone();
        for p in order {
            let w = canonical(&mut s, &lam, p);
            if krhom::schur::flows(&lam, &w) > maxf {
                continue;
            }
            let ch = s.ladder_char(&w);
            let unit = s.unit();
            let gi = gens_list(&hom_gens(&mut s, &unit, &ch));
            let go = gens_list(&hom_gens(&mut s, &ch, &unit));
            let (ni, no) = nakayama_gens(&lam, &w);
            assert_eq!(gi, ni, "B1 {:?} type {}: hom_gens(R, L_p) != Nakayama ins", lam, s.table(p).fmt());
            assert_eq!(go, no, "B1 {:?} type {}: hom_gens(L_p, R) != Nakayama outs", lam, s.table(p).fmt());
            n += 1;
        }
    }
    eprintln!("[B1] {} closed canonical ladders: hom_gens = Nakayama generator degrees (ins and outs)", n);
}

// ------------------------------------------------------------------ B2

fn parse_hom_rungs(s: &str) -> Vec<Rung> {
    s.split_whitespace()
        .map(|t| {
            let dir = if t.starts_with('E') { Dir::E } else { Dir::F };
            let (i, r) = t[1..].split_once(':').unwrap();
            Rung { i: i.parse().unwrap(), dir, r: r.parse().unwrap() }
        })
        .collect()
}

#[test]
fn b2_hom60_through_exact() {
    let text = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/schur/hom60.txt")).unwrap();
    let t0 = Instant::now();
    let (mut cases, mut spaces, mut maps) = (0, 0, 0);
    for line in text.lines().filter(|l| l.starts_with("HOM ")) {
        let fl: Vec<&str> = line[4..].split('|').collect();
        let lam: Vec<u8> = fl[0].split(',').map(|x| x.parse().unwrap()).collect();
        let (a, b) = (parse_hom_rungs(fl[1]), parse_hom_rungs(fl[2]));
        let dims: Vec<(i32, usize)> = fl[5].split_whitespace().map(|kv| {
            let (k, v) = kv.split_once(':').unwrap();
            (k.parse().unwrap(), v.parse().unwrap())
        }).collect();
        if lam.iter().all(|&x| x == 0) {
            continue;
        }
        let mut ar: Arena<F> = Arena::new(&lam, TmplLib::new(None));
        let mut s = Schur::new(&lam);
        let (sa, sb) = (ar.expr(ex(&lam, &a)), ar.expr(ex(&lam, &b)));
        // the hom60 fibres are small: admit labels up to 3 by the label rule for this oracle test
        let mut exact = Exact { max_label: 3, ..Exact::default() };
        for (d, want) in dims {
            let r = HomProvider::<F>::full(&mut exact, &mut ar, &mut s, sa, sb, d, &Budget::default());
            let (ms, _) = r.unwrap_or_else(|m| panic!("B2 {:?} {:?} -> {:?} degree {}: {:?}", lam, a, b, d, m));
            assert_eq!(ms.len(), want, "B2: certified dimension != Hom formula");
            spaces += 1;
            maps += ms.len();
        }
        cases += 1;
    }
    eprintln!("[B2] {} hom60 cases, {} graded pieces, {} certified maps through Exact in {:.1}s", cases, spaces, maps, t0.elapsed().as_secs_f64());
    assert_eq!(cases, 60);
}

/// all type pairs of the block (lambda, lambda), degrees 0..=dmax, through Chain::for_labels
fn chain_type_pairs(lam: &[u8], dmax: i32) -> (usize, usize, Vec<String>) {
    let mut s = Schur::new(lam);
    let order = s.canonical_basis(lam).order.clone();
    let mut ar: Arena<F> = Arena::new(lam, TmplLib::new(Some(cache_dir(&format!("b2_{:?}", lam)))));
    let mut prov: Chain<F> = Chain::for_labels(lam.iter().copied().max().unwrap(), Exact::default());
    let ids: Vec<(Tid, ExprId)> = order.iter().map(|&p| { let w = canonical(&mut s, lam, p); (p, ar.expr(ex(lam, &w))) }).collect();
    let (mut ok, mut tot) = (0, 0);
    let mut fails = vec![];
    for &(p, a) in &ids {
        for &(q, b) in &ids {
            for d in 0..=dmax {
                tot += 1;
                match prov.full(&mut ar, &mut s, a, b, d, &Budget::default()) {
                    Ok(_) => ok += 1,
                    Err(m) => fails.push(format!("{} -> {} degree {}: {:?}", s.table(p).fmt(), s.table(q).fmt(), d, m.iter().map(|x| format!("{:?} {}", x.kind, x.what)).collect::<Vec<_>>())),
                }
            }
        }
    }
    (ok, tot, fails)
}

#[test]
fn b2_type_pairs_22_212_through_chain() {
    for lam in [vec![2u8, 2], vec![2, 1, 2]] {
        let t0 = Instant::now();
        let (ok, tot, fails) = chain_type_pairs(&lam, 4);
        eprintln!("[B2] {:?}: {}/{} Hom^d (d = 0..4) certified through the Chain in {:.1}s", lam, ok, tot, t0.elapsed().as_secs_f64());
        for f in &fails {
            eprintln!("[B2]   missing: {}", f);
        }
        assert!(fails.is_empty(), "B2 {:?}: {} of {} not certified", lam, tot - ok, tot);
    }
}

#[test]
#[ignore = "B2 (3,3) and (2,2,2) type pairs through the Chain: long (reported once)"]
fn b2_type_pairs_33_222_through_chain() {
    for lam in [vec![3u8, 3], vec![2, 2, 2]] {
        let t0 = Instant::now();
        let (ok, tot, fails) = chain_type_pairs(&lam, 4);
        eprintln!("[B2] {:?}: {}/{} Hom^d (d = 0..4) certified through the Chain in {:.1}s", lam, ok, tot, t0.elapsed().as_secs_f64());
        for f in fails.iter().take(40) {
            eprintln!("[B2]   missing: {}", f);
        }
        assert!(fails.is_empty(), "B2 {:?}: {} of {} not certified", lam, tot - ok, tot);
    }
}

// ------------------------------------------------------------------ B3 / B4

struct FaithCase {
    lam: Vec<u8>,
    src: Vec<Rung>,
    tgt: Vec<Rung>,
    d: i32,
    dim: usize,
    m: usize,
    full1: usize,
    wit1: usize,
}

fn faith_cases() -> Vec<FaithCase> {
    let ef = |i: usize| vec![e(i, 1), f(i, 1)];
    let c = |lam: Vec<u8>, src: Vec<Rung>, tgt: Vec<Rung>, d, dim, m, full1, wit1| FaithCase { lam, src, tgt, d, dim, m, full1, wit1 };
    let bss = [ef(0), ef(0)].concat();
    let both = [ef(0), ef(1)].concat();
    vec![
        c(vec![1, 1], bss.clone(), bss.clone(), 0, 5, 2, 4, 4),
        c(vec![2, 2], ef(0), ef(0), 2, 5, 2, 4, 4),
        c(vec![2, 2], ef(0), ef(0), 4, 14, 5, 7, 5),
        c(vec![2, 2], vec![], ef(0), 3, 3, 2, 2, 2),
        c(vec![2, 2], ef(0), vec![e(0, 2), f(0, 2)], 3, 4, 2, 3, 3),
        c(vec![2, 2], bss.clone(), bss.clone(), 0, 48, 5, 36, 35),
        c(vec![1, 1, 1], both.clone(), both.clone(), 2, 5, 3, 3, 3),
        c(vec![2, 2, 2], ef(0), ef(0), 4, 21, 9, 7, 5),
        c(vec![2, 2, 2], both.clone(), both.clone(), 2, 9, 3, 7, 7),
    ]
}

/// B3 and B4 on one case: M from the oracle, points to full rank (full orbits and witnesses),
/// rank at one point
fn faith_check(c: &FaithCase) {
    let t0 = Instant::now();
    let mut ar: Arena<F> = Arena::new(&c.lam, TmplLib::new(None));
    let mut s = Schur::new(&c.lam);
    let (a, b) = (ar.expr(ex(&c.lam, &c.src)), ar.expr(ex(&c.lam, &c.tgt)));
    let mut exact = Exact { max_label: 3, max_fibre_product: u64::MAX, max_unknowns: u64::MAX, ..Exact::default() };
    let (ms, _) = HomProvider::<F>::full(&mut exact, &mut ar, &mut s, a, b, c.d, &Budget::default()).unwrap_or_else(|m| panic!("B3: {:?}", m));
    assert_eq!(ms.len(), c.dim, "B3: dim");
    let plan = detect_plan(&ar, &mut s, a, b, c.d);
    assert_eq!(plan.m, c.m, "B3: M(p,q,d) from the oracle for {:?} {:?} -> {:?} degree {}", c.lam, c.src, c.tgt, c.d);
    // full cosets: every coset of the joint support
    let lay = Layout::of(ar.flows(a), ar.flows(b));
    let full: Vec<u64> = lay.cosets.clone();
    let wit = detecting_cosets(&ar, &mut s, a, b);
    let pts = points::<F>(ar.k(), 31337, c.m + 3, 0);
    let first_full = (1..=pts.len()).find(|&k| stacked_rank(&ar, &ms, &pts[..k], &full) == c.dim).unwrap();
    let first_wit = (1..=pts.len()).find(|&k| stacked_rank(&ar, &ms, &pts[..k], &wit) == c.dim).unwrap();
    let r1f = stacked_rank(&ar, &ms, &pts[..1], &full);
    let r1w = stacked_rank(&ar, &ms, &pts[..1], &wit);
    eprintln!(
        "[B3/B4] {:?} {} -> {} degree {}: dim {}, M {}, points to full rank {} (full) / {} (witness), one point {}/{} (full) {}/{} (witness), detect {} with {} points; {:.1}s",
        c.lam, wstr(&c.src), wstr(&c.tgt), c.d, c.dim, plan.m, first_full, first_wit, r1f, c.dim, r1w, c.dim, plan.kind, plan.npoints, t0.elapsed().as_secs_f64()
    );
    assert_eq!(first_full, c.m, "B3: points to full rank (full cosets) != M");
    assert_eq!(first_wit, c.m, "B3: points to full rank (witness cosets) != M");
    assert_eq!(r1f, c.full1, "B4: one-point rank at the full cosets");
    assert_eq!(r1w, c.wit1, "B4: one-point rank at the witness cosets");
    certify_full(&ar, &mut s, a, b, c.d, &ms, 77).unwrap();
}

#[test]
fn b3_b4_faith_table() {
    for (k, c) in faith_cases().iter().enumerate() {
        if k == 5 {
            continue; // (2,2) E1F1E1F1 End^0: ignored test below
        }
        faith_check(c);
    }
}

#[test]
#[ignore = "B3/B4: (2,2) E1F1E1F1 End^0 (dim 48, about 20 s in debug)"]
fn b3_b4_faith_e1f1e1f1() {
    faith_check(&faith_cases()[5]);
}

#[test]
#[ignore = "B4b: (2,2) at degree 11 (hom_space of dim 69; long)"]
fn b4b_detecting_set_fallback() {
    // (2,2) types with Hom^11 of dim 69, M 20, E_w 3: the witness bound fails
    let lam = [2u8, 2];
    let mut s = Schur::new(&lam);
    let order = s.canonical_basis(&lam).order.clone();
    let mut done = false;
    for &p in &order {
        for &q in &order {
            let (wp, wq) = (canonical(&mut s, &lam, p), canonical(&mut s, &lam, q));
            let mut ar: Arena<F> = Arena::new(&lam, TmplLib::new(None));
            let (a, b) = (ar.expr(ex(&lam, &wp)), ar.expr(ex(&lam, &wq)));
            let plan = detect_plan(&ar, &mut s, a, b, 11);
            if plan.dim != 69 || plan.kind == "witness" {
                continue;
            }
            let mut exact = Exact { max_fibre_product: u64::MAX, max_unknowns: u64::MAX, ..Exact::default() };
            let t0 = Instant::now();
            let (ms, st) = HomProvider::<F>::full(&mut exact, &mut ar, &mut s, a, b, 11, &Budget::default()).unwrap();
            eprintln!("[B4b] {} -> {} degree 11: dim {}, M {}, E_w {}, detect {} ({} points), {:.1}s", s.table(p).fmt(), s.table(q).fmt(), ms.len(), plan.m, plan.rows_per_point, st.detect, plan.npoints, t0.elapsed().as_secs_f64());
            assert_ne!(st.detect, "witness");
            // forcing the witness set with M + 2 points fails the counting assertion
            let wit = detecting_cosets(&ar, &mut s, a, b);
            let rows: usize = Layout::of(ar.flows(a), ar.flows(b)).cosets.iter().filter(|c| wit.contains(c)).count();
            assert!(plan.rows_per_point * (plan.m + 2) < 69 || rows * (plan.m + 2) < 69);
            done = true;
        }
    }
    assert!(done, "B4b: no (2,2) pair with dim Hom^11 = 69 and a violated witness bound");
}

// ------------------------------------------------------------------ B5

#[test]
fn b5_genuineness() {
    let mut lib: TmplLib<F> = TmplLib::new(None);
    let (wins, ids) = lib.rickard(2, 1, true);
    let t = ids[0];
    let tm = lib.get(t).clone();
    // a copy of the map entered as a template is not a certified library leaf
    let copy = lib.add_template(Template::new(tm.src.clone(), tm.tgt.clone(), tm.deg, tm.phi.clone(), TmplSource::Rickard { c: 2, d: 1, positive: true, j: 0 }));
    let lab = [2u8, 1];
    let mut ar = Arena::new(&lab, lib);
    let (w0, w1) = (ar.expr(wins[0].expr()), ar.expr(wins[1].expr()));
    let genuine = ar.add(MapNode::Local { src: w0, tgt: w1, at: 0, first: 0, tmpl: t }, "rickard").unwrap();
    let fake = ar.add(MapNode::Local { src: w0, tgt: w1, at: 0, first: 0, tmpl: copy }, "copy").unwrap();
    assert_genuine(&ar, genuine).unwrap();
    let err = assert_genuine(&ar, fake).unwrap_err();
    assert!(err.contains("uncertified"), "{}", err);
    // the point-rescaled copy (A10 b): its polynomiality certificate fails
    let r = poly_certificate::<F>(&tm.src, &tm.tgt, tm.deg, &mut |y: &[F]| {
        let mut wl = ar.lib_ref().localize(t, y);
        let sc = wl.vals[0].inv();
        for v in wl.vals.iter_mut() {
            *v = v.mul(&sc);
        }
        wl
    }, 8);
    assert!(r.is_err(), "B5: point-rescaled map passed the certificate");
    // certify_pairs-style check: a stacked rank cannot tell them apart (rank is not genuineness)
    let pts = points::<F>(3, 5, 2, 0);
    let lay = Layout::of(ar.flows(w0), ar.flows(w1));
    assert_eq!(stacked_rank(&ar, &[genuine], &pts, &lay.cosets), stacked_rank(&ar, &[fake], &pts, &lay.cosets));
    eprintln!("[B5] rejected: {}; certificate of the rescaled copy: {}", err, r.unwrap_err());
}

// ------------------------------------------------------------------ B6

/// Search pairs on (2,2) and (2,1) decompositions: form ranks = term_mult, and Search ι, π lie
/// in the span of the exact Hom spaces (stacked rank of exact ∪ search = rank of exact)
#[test]
fn b6_search_pairs_match_exact() {
    let mut n = 0;
    for lam in [vec![2u8, 2], vec![2, 1], vec![1, 2]] {
        let mut s = Schur::new(&lam);
        let (pairs, _) = g1_targets(&mut s);
        let mut ar: Arena<F> = Arena::new(&lam, TmplLib::new(Some(cache_dir(&format!("b6_{:?}", lam)))));
        let mut search = Search::tiers(0, 6);
        let mut exact = Exact { max_fibre_product: u64::MAX, max_unknowns: u64::MAX, ..Exact::default() };
        for pr in pairs.iter().filter(|p| !p.identity) {
            let mut xw = canonical(&mut s, &pr.mu, pr.p);
            xw.push(pr.rung);
            let x = ar.expr(ex(&lam, &xw));
            for &(q, sd, mult) in &pr.targets {
                let mq = s.table(q).cols();
                let lqw = canonical(&mut s, &mq, q);
                let lq = ar.expr(ex(&lam, &lqw));
                let (p, st) = search.pairs(&mut ar, &mut s, x, q, lq, None, sd, mult, &Budget::default()).unwrap_or_else(|m| panic!("B6 {:?} {} {}: {:?}\n{:?}", lam, s.table(q).fmt(), sd, m, search.trace));
                assert_eq!(p.iotas.len(), mult);
                assert_eq!(krhom::dense::rank(&p.form, mult), mult, "B6: form rank != term_mult");
                // exact spaces
                let (ei, _) = HomProvider::<F>::full(&mut exact, &mut ar, &mut s, lq, x, sd, &Budget::default()).unwrap();
                let (ep, _) = HomProvider::<F>::full(&mut exact, &mut ar, &mut s, x, lq, -sd, &Budget::default()).unwrap();
                for (exm, sm, a, b, dd) in [(&ei, &p.iotas, lq, x, sd), (&ep, &p.pis, x, lq, -sd)] {
                    let plan = detect_plan(&ar, &mut s, a, b, dd);
                    let pts = points::<F>(ar.k(), 606, plan.npoints, 0);
                    let r0 = stacked_rank(&ar, exm, &pts, &plan.cosets);
                    let mut all = exm.clone();
                    all.extend(sm.iter().cloned());
                    assert_eq!(stacked_rank(&ar, &all, &pts, &plan.cosets), r0, "B6: a Search map is not in the span of the exact Hom space");
                }
                n += 1;
                let _ = st;
            }
        }
    }
    eprintln!("[B6] {} (q, s) targets of (2,2), (2,1), (1,2): Search pairs certified, in the exact span, form ranks = term_mult", n);
}

// ------------------------------------------------------------------ B7

#[test]
fn b7_determinism() {
    let run = || -> Vec<(Vec<MapId>, Vec<MapId>, usize, usize, [usize; 7])> {
        let lam = [2u8, 2, 2];
        let mut s = Schur::new(&lam);
        let (pairs, _) = g1_targets(&mut s);
        let mut ar: Arena<F> = Arena::new(&lam, TmplLib::new(None));
        let mut search = Search::tiers(0, 6);
        let mut out = vec![];
        for pr in pairs.iter().filter(|p| !p.identity).take(3) {
            let mut xw = canonical(&mut s, &pr.mu, pr.p);
            xw.push(pr.rung);
            let x = ar.expr(ex(&lam, &xw));
            for &(q, sd, mult) in &pr.targets {
                let mq = s.table(q).cols();
                let lqw = canonical(&mut s, &mq, q);
                let lq = ar.expr(ex(&lam, &lqw));
                let (p, st) = search.pairs(&mut ar, &mut s, x, q, lq, None, sd, mult, &Budget::default()).unwrap();
                out.push((p.iotas, p.pis, st.candidates, st.visited_words, st.tier_hits));
            }
        }
        out
    };
    let (a, b) = (run(), run());
    assert_eq!(a, b, "B7: two runs differ");
    eprintln!("[B7] {} targets: identical MapIds and counters in two runs", a.len());
}

// ------------------------------------------------------------------ B8

#[test]
fn b8_admission_order_and_scratch() {
    // (1) the (3,2) one-step request 144 -> 8 in degree 2 is refused without a call
    let lam = [3u8, 2];
    let mut s = Schur::new(&lam);
    let mut found = None;
    'outer: for mu in [vec![3u8, 2], vec![4, 1], vec![5, 0], vec![2, 3], vec![1, 4]] {
        let order = s.canonical_basis(&mu).order.clone();
        for &p in &order {
            let lp = canonical(&mut s, &mu, p);
            for t in krhom::schur::crossing_terms(mu[0], mu[1], true, 0) {
                let mut x = lp.clone();
                x.extend(t.rungs.iter().cloned());
                if Expr::new(&lam, &x).is_err() || krhom::schur::flows(&lam, &x) != 144 {
                    continue;
                }
                let top = Expr::new(&lam, &x).unwrap().top();
                for &q in &s.canonical_basis(&top).order.clone() {
                    let lq = canonical(&mut s, &top, q);
                    if krhom::schur::flows(&lam, &lq) == 8 {
                        found = Some((x.clone(), lq, q));
                        break 'outer;
                    }
                }
            }
        }
    }
    let (xw, lqw, q) = found.expect("B8: no (3,2) one-step 144 -> 8 request");
    let mut ar: Arena<F> = Arena::new(&lam, TmplLib::new(None));
    let (x, lq) = (ar.expr(ex(&lam, &xw)), ar.expr(ex(&lam, &lqw)));
    let mut exact = Exact::default();
    let t0 = Instant::now();
    let err = HomProvider::<F>::full(&mut exact, &mut ar, &mut s, x, lq, 2, &Budget::default()).unwrap_err();
    assert_eq!(err[0].kind, MissKind::Refused, "{:?}", err);
    assert!(t0.elapsed().as_secs_f64() < 1.0, "B8: the refusal took a call");
    assert_eq!(ar.lib_ref().stats().computed, 0, "B8: exact_calls != 0");
    let err = HomProvider::<F>::pairs(&mut exact, &mut ar, &mut s, x, q, lq, None, 0, 1, &Budget::default()).unwrap_err();
    assert_eq!(err[0].kind, MissKind::Refused);
    assert_eq!(ar.lib_ref().stats().computed, 0);
    eprintln!("[B8] (3,2) 144 -> 8 degree 2 refused without a call: {}", err[0].what);

    // (2) Chain::for_labels(3) on a (3,3) pair whose extremes certify: no Exact call
    let lam = [3u8, 3];
    let mut s = Schur::new(&lam);
    let (pairs, _) = g1_targets(&mut s);
    let pr = pairs.iter().find(|p| !p.identity).unwrap();
    let mut xw = canonical(&mut s, &pr.mu, pr.p);
    xw.push(pr.rung);
    let mut ar: Arena<F> = Arena::new(&lam, TmplLib::new(Some(cache_dir("b8_33"))));
    let x = ar.expr(ex(&lam, &xw));
    let mut chain: Chain<F> = Chain::for_labels(3, Exact::default());
    for &(q, sd, mult) in &pr.targets {
        let mq = s.table(q).cols();
        let lqw = canonical(&mut s, &mq, q);
        let lq = ar.expr(ex(&lam, &lqw));
        // (3) scratch region: the Arena grows by the certified sub-DAG only
        let n0 = ar.n_nodes();
        let (p, st) = chain.pairs(&mut ar, &mut s, x, q, lq, None, sd, mult, &Budget::default()).unwrap_or_else(|m| panic!("B8 (3,3): {:?}", m));
        assert_eq!(st.exact_calls, 0, "B8: Exact backend called for a (3,3) pair certified by S1/S2");
        assert!(st.tier_hits[1] + st.tier_hits[2] > 0, "B8: not certified by S1/S2 ({:?})", st.tier_hits);
        let mut reach = std::collections::BTreeSet::new();
        let mut stack: Vec<MapId> = p.iotas.iter().chain(&p.pis).cloned().collect();
        while let Some(m) = stack.pop() {
            if (m as usize) < n0 || !reach.insert(m) {
                continue;
            }
            match ar.node(m) {
                MapNode::Comp(g, h) => stack.extend([*g, *h]),
                MapNode::Lin(t) => stack.extend(t.iter().map(|z| z.1)),
                MapNode::Kron(g, _) => stack.push(*g),
                _ => {}
            }
        }
        assert_eq!(ar.n_nodes(), n0 + reach.len(), "B8: scratch nodes survived the call");
        eprintln!("[B8] (3,3) {} -> {} s {}: tiers {:?}, {} kept nodes", s.table(pr.p).fmt(), s.table(q).fmt(), sd, st.tier_hits, reach.len());
    }
}

// ------------------------------------------------------------------ closure (G1 shape, small)

#[test]
fn closure_22_through_chain() {
    for lam in [vec![2u8, 2], vec![2, 1, 2]] {
        let mut s = Schur::new(&lam);
        let order = s.canonical_basis(&lam).order.clone();
        let mut ar: Arena<F> = Arena::new(&lam, TmplLib::new(Some(cache_dir(&format!("cl_{:?}", lam)))));
        let mut chain: Chain<F> = Chain::for_labels(2, Exact::default());
        let mut search = Search::tiers(0, 6);
        for &p in &order {
            let w = canonical(&mut s, &lam, p);
            let lp = ar.expr(ex(&lam, &w));
            for outs in [false, true] {
                let (ms, _) = chain.closure(&mut ar, &mut s, lp, p, None, outs, &Budget::default()).unwrap_or_else(|m| panic!("closure {:?}: {:?}", lam, m));
                assert_eq!(ms.len() as i128, s.c_diag(p).at1());
                if lam.len() != 2 {
                    continue;
                }
                let (ms2, st) = search.closure(&mut ar, &mut s, lp, p, None, outs, &Budget::default()).unwrap_or_else(|m| panic!("Search closure {:?} {}: {:?}", lam, s.table(p).fmt(), m));
                assert_eq!(ms2.len(), ms.len());
                let _ = st;
            }
        }
        eprintln!("[closure] {:?}: {} types, ins and outs certified through the Chain{}", lam, order.len(), if lam.len() == 2 { " and through Search alone" } else { "" });
    }
}

#[test]
fn compile_time_send_sync() {
    fn _ss<T: Send + Sync>() {}
    _ss::<CostTable>();
    _ss::<Exact>();
    _ss::<Budget>();
}
