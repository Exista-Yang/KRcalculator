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

/// certified counts per provider / composition: "Exact" (no Search tier hit) or the Search
/// tiers that produced the basis ("S0", "S2", "S4", "S6", combinations joined by '+')
type Counts = std::collections::BTreeMap<String, usize>;
fn composition(st: &SearchStats, n: usize) -> String {
    if n == 0 {
        return "dim0".into();
    }
    let t: Vec<String> = st.tier_hits.iter().enumerate().filter(|(_, c)| **c > 0).map(|(i, _)| format!("S{}", i)).collect();
    if t.is_empty() { "Exact".into() } else { t.join("+") }
}

/// all type pairs of the block (lambda, lambda), degrees 0..=dmax, through Chain::for_labels
fn chain_type_pairs(lam: &[u8], dmax: i32) -> (usize, usize, Vec<String>) {
    let mut prov: Chain<F> = Chain::for_labels(lam.iter().copied().max().unwrap(), Exact::default());
    let (ok, tot, fails, counts) = type_pairs_with(lam, dmax, &mut prov, "b2");
    eprintln!("[B2] {:?} through the Chain, per composition: {:?}", lam, counts);
    assert_eq!(counts.values().sum::<usize>(), ok);
    (ok, tot, fails)
}

/// all type pairs of (lambda, lambda), degrees 0..=dmax, through one provider; counts per composition
fn type_pairs_with(lam: &[u8], dmax: i32, prov: &mut dyn HomProvider<F>, tag: &str) -> (usize, usize, Vec<String>, Counts) {
    let mut s = Schur::new(lam);
    let order = s.canonical_basis(lam).order.clone();
    let mut ar: Arena<F> = Arena::new(lam, TmplLib::new(Some(cache_dir(&format!("{}_{:?}", tag, lam)))));
    let ids: Vec<(Tid, ExprId)> = order.iter().map(|&p| { let w = canonical(&mut s, lam, p); (p, ar.expr(ex(lam, &w))) }).collect();
    let (mut ok, mut tot) = (0, 0);
    let mut fails = vec![];
    let mut counts = Counts::new();
    for &(p, a) in &ids {
        for &(q, b) in &ids {
            for d in 0..=dmax {
                tot += 1;
                match prov.full(&mut ar, &mut s, a, b, d, &Budget::default()) {
                    Ok((ms, st)) => {
                        ok += 1;
                        *counts.entry(composition(&st, ms.len())).or_insert(0) += 1;
                    }
                    Err(m) => fails.push(format!("{} -> {} degree {}: {:?}", s.table(p).fmt(), s.table(q).fmt(), d, m.iter().map(|x| format!("{:?} {}", x.kind, x.what)).collect::<Vec<_>>())),
                }
            }
        }
    }
    (ok, tot, fails, counts)
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

/// B2 through Search S0-S6 alone (review tests finding 2: the Chain run is now certified mostly
/// by Exact, so Search `full` completeness on these blocks is measured separately)
#[test]
fn b2_type_pairs_22_212_search_alone() {
    for lam in [vec![2u8, 2], vec![2, 1, 2]] {
        let t0 = Instant::now();
        let mut search = Search::tiers(0, 6);
        let (ok, tot, fails, counts) = type_pairs_with(&lam, 4, &mut search, "b2s");
        eprintln!("[B2 Search] {:?}: {}/{} Hom^d (d = 0..4) certified by Search S0-S6 alone in {:.1}s; per composition {:?}", lam, ok, tot, t0.elapsed().as_secs_f64(), counts);
        for f in fails.iter().take(40) {
            eprintln!("[B2 Search]   missing: {}", f);
        }
        // Measured 2026-11-11: (2,2) 45/45; (2,1,2) 581/605. The 24 misses are Exhausted (a Search
        // `full` completeness gap on endomorphism-like spaces of the types 1,0,1/0,1,0/1,0,1,
        // 0,0,2/1,0,0/1,1,0, ... in degrees 2-4: dots + S4 + S6 windows span too little), never
        // Time/Memory/Refused. The Chain covers them by Exact. This bound may only improve.
        let allowed = if lam.len() == 2 { 0 } else { 24 };
        assert!(fails.len() <= allowed, "B2 Search {:?}: {} of {} not certified (measured bound {})", lam, tot - ok, tot, allowed);
        assert!(fails.iter().all(|f| f.contains("\"Exhausted ")), "B2 Search {:?}: a miss that is not Exhausted", lam);
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
            // forcing the witness set with M + 2 points fails the counting assertion, before any
            // rank test (SPEC §4.6 B4b; review tests finding 6)
            let wit = detecting_cosets(&ar, &mut s, a, b);
            let ew: usize = {
                let (ca, cb) = (s.ladder_char(&wp), s.ladder_char(&wq));
                ca.t.iter().map(|(t, c)| (c.at1() * cb.coef(*t).at1()) as usize).sum()
            };
            assert!(ew * (plan.m + 2) < 69, "B4b: the witness bound is not violated");
            let forced = DetectPlan { cosets: wit, npoints: plan.m + 2, kind: "witness", dim: plan.dim, m: plan.m, rows_per_point: ew };
            let err = certify_full_at(&ar, &mut s, a, b, 11, &ms, 0x5eed_0b4b, forced).unwrap_err();
            assert!(err.contains("counting bound violated"), "B4b: forced witness set failed for another reason: {}", err);
            // the exact basis itself certifies through the fallback
            certify_full(&ar, &mut s, a, b, 11, &ms, 0x5eed_0b4c).unwrap();
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
        eprintln!("[B8] (3,3) target {} s {}: {} Search window calls ({:.2}s), 0 Exact-backend calls", s.table(q).fmt(), sd, st.window_calls, st.window_secs);
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

// ------------------------------------------------------------------ hardening (review 2026-11-11)

/// a (2,2) pair job and its first target, for the precondition tests
fn first_target(s: &mut Schur, ar: &mut Arena<F>, lam: &[u8]) -> (ExprId, Tid, ExprId, i32, usize, Vec<(Tid, i32, usize)>) {
    let (pairs, _) = g1_targets(s);
    let pr = pairs.iter().find(|p| !p.identity && p.targets.len() >= 2).unwrap();
    let (mu, pp, rung, targets) = (pr.mu.clone(), pr.p, pr.rung, pr.targets.clone());
    let mut xw = canonical(s, &mu, pp);
    xw.push(rung);
    let x = ar.expr(ex(lam, &xw));
    let (q, sd, mult) = targets[0];
    let mq = s.table(q).cols();
    let lqw = canonical(s, &mq, q);
    let lq = ar.expr(ex(lam, &lqw));
    (x, q, lq, sd, mult, targets)
}

/// review F1/F2: `pairs` refuses a wrong `mult`, an L_q of another type and a decomposable L_q
/// without e_q (both backends); with e_q = id on an indecomposable L_q it certifies; a non-
/// idempotent e_q is refused
#[test]
fn b9_pairs_preconditions() {
    let lam = [2u8, 2];
    let mut s = Schur::new(&lam);
    let mut ar: Arena<F> = Arena::new(&lam, TmplLib::new(Some(cache_dir("b9"))));
    let (x, q, lq, sd, mult, targets) = first_target(&mut s, &mut ar, &lam);
    let b = Budget::default();
    let mut exact = Exact::default();
    let mut search = Search::tiers(0, 6);
    let refused = |r: Result<(Pairs<F>, SearchStats), Vec<Missing>>, what: &str| {
        let m = match r {
            Ok(_) => panic!("B9 {}: certified", what),
            Err(m) => m,
        };
        assert_eq!(m[0].kind, MissKind::Refused, "B9 {}: {:?}", what, m);
        m[0].what.clone()
    };
    // (a) wrong multiplicity: mult + 1 and mult - 1 (never a partial or padded family)
    for bad in [mult + 1, mult.saturating_sub(1)] {
        if bad == mult {
            continue;
        }
        let w = refused(HomProvider::<F>::pairs(&mut exact, &mut ar, &mut s, x, q, lq, None, sd, bad, &b), "Exact wrong mult");
        assert!(w.contains("mult"), "{}", w);
        let w = refused(search.pairs(&mut ar, &mut s, x, q, lq, None, sd, bad, &b), "Search wrong mult");
        assert!(w.contains("mult"), "{}", w);
    }
    // (b) L_q of another type of the same top
    let other = targets.iter().map(|t| t.0).find(|&t| t != q && s.table(t).cols() == s.table(q).cols());
    if let Some(q2) = other {
        let w = refused(search.pairs(&mut ar, &mut s, x, q2, lq, None, sd, mult, &b), "Search L_q of another type");
        assert!(w.contains("e_q required"), "{}", w);
    }
    // (c) a decomposable "L_q" (the word X itself) without e_q
    let w = refused(HomProvider::<F>::pairs(&mut exact, &mut ar, &mut s, x, q, x, None, 0, 1, &b), "Exact decomposable L_q");
    assert!(w.contains("L_q decomposable: e_q required"), "{}", w);
    let w = refused(search.pairs(&mut ar, &mut s, x, q, x, None, 0, 1, &b), "Search decomposable L_q");
    assert!(w.contains("L_q decomposable: e_q required"), "{}", w);
    // (d) e_q = id_{L_q} passes the idempotent contract and certifies; 2·id is refused
    let id = ar.add(MapNode::Id(lq), "id").unwrap();
    let (p, _) = HomProvider::<F>::pairs(&mut exact, &mut ar, &mut s, x, q, lq, Some(id), sd, mult, &b).unwrap_or_else(|m| panic!("B9 e_q = id: {:?}", m));
    assert_eq!(p.iotas.len(), mult);
    let two = ar.add(MapNode::Lin(vec![(F::from_i64(2), id)]), "2 id").unwrap();
    let w = refused(HomProvider::<F>::pairs(&mut exact, &mut ar, &mut s, x, q, lq, Some(two), sd, mult, &b), "Exact e_q = 2 id");
    assert!(w.contains("e_q rejected"), "{}", w);
    let w = refused(search.pairs(&mut ar, &mut s, x, q, lq, Some(two), sd, mult, &b), "Search e_q = 2 id");
    assert!(w.contains("e_q rejected"), "{}", w);
    // the Chain reports the refusal of every provider
    let mut chain: Chain<F> = Chain::for_labels(2, Exact::default());
    let m = chain.pairs(&mut ar, &mut s, x, q, lq, None, sd, mult + 1, &b).unwrap_err();
    assert!(m.iter().all(|z| z.kind == MissKind::Refused), "{:?}", m);
    eprintln!("[B9] preconditions refused: wrong mult, other type, decomposable L_q, e_q = 2 id; e_q = id certifies");
}

/// review F3: certify_pairs / certify_closure check the ends and degrees of every map
#[test]
fn b9_certify_checks_ends_and_degrees() {
    let lam = [2u8, 2];
    let mut s = Schur::new(&lam);
    let mut ar: Arena<F> = Arena::new(&lam, TmplLib::new(Some(cache_dir("b9e"))));
    let (x, q, lq, sd, mult, _) = first_target(&mut s, &mut ar, &lam);
    let mut exact = Exact::default();
    let (p, _) = HomProvider::<F>::pairs(&mut exact, &mut ar, &mut s, x, q, lq, None, sd, mult, &Budget::default()).unwrap();
    let c = {
        let tab = s.table(q).clone();
        let cs = krhom::coset::CosetSpace::new(&lam, &tab.cols());
        cs.rank(&cs.rep(&tab))
    };
    certify_pairs(&ar, x, lq, sd, &p.iotas, &p.pis, c, mult, 7).unwrap();
    // a degree label off by 2 (the form cannot see it: π∘ι still has degree 0)
    let e = certify_pairs(&ar, x, lq, sd + 2, &p.iotas, &p.pis, c, mult, 7).unwrap_err();
    assert!(e.contains("degree"), "{}", e);
    // swapped ends
    let e = certify_pairs(&ar, lq, x, sd, &p.iotas, &p.pis, c, mult, 7).unwrap_err();
    assert!(e.contains("expected"), "{}", e);
    // closure: ins presented as outs
    let order = s.canonical_basis(&lam).order.clone();
    let pt = order[order.len() / 2];
    let lpw = canonical(&mut s, &lam, pt);
    let lp = ar.expr(ex(&lam, &lpw));
    let (ins, _) = HomProvider::<F>::closure(&mut exact, &mut ar, &mut s, lp, pt, None, false, &Budget::default()).unwrap();
    let r0 = ar.expr(ex(&lam, &[]));
    let mut want: Vec<i32> = ins.iter().map(|&m| ar.info(m).deg).collect();
    want.sort();
    certify_closure(&ar, &lam, r0, lp, None, &ins, &want, false, 9).unwrap();
    let e = certify_closure(&ar, &lam, r0, lp, None, &ins, &want, true, 9).unwrap_err();
    assert!(e.contains("expected"), "{}", e);
    eprintln!("[B9] certify_pairs / certify_closure reject wrong degree labels and ends");
}

/// review F3: Search's TmplId caches are bound to the template library of the Arena
#[test]
fn b9_search_caches_bound_to_library() {
    let lam = [2u8, 2];
    let mut s = Schur::new(&lam);
    let mut ar: Arena<F> = Arena::new(&lam, TmplLib::new(None));
    let (x, q, lq, sd, mult, _) = first_target(&mut s, &mut ar, &lam);
    let mut search = Search::tiers(0, 6);
    search.pairs(&mut ar, &mut s, x, q, lq, None, sd, mult, &Budget::default()).unwrap();
    let (a, c) = search.cached_windows();
    assert!(a + c > 0, "B9: nothing cached");
    // a second Arena with another library: the caches are dropped before any use
    let mut ar2: Arena<F> = Arena::new(&lam, TmplLib::new(None));
    let (x2, q2, lq2, sd2, mult2, _) = first_target(&mut s, &mut ar2, &lam);
    let _ = search.pairs(&mut ar2, &mut s, x2, q2, lq2, None, sd2, mult2 + 1, &Budget::default());
    assert_eq!(search.cached_windows(), (0, 0), "B9: TmplId caches survived a library change");
    let (p, _) = search.pairs(&mut ar2, &mut s, x2, q2, lq2, None, sd2, mult2, &Budget::default()).unwrap();
    assert_eq!(p.iotas.len(), mult2);
    eprintln!("[B9] Search caches ({} splits, {} comms) cleared on a library change", a, c);
}

/// review tests finding 3: censored records never admit (bucket +infinity), stamps, admission rules
#[test]
fn b9_cost_table_admission() {
    let rec = |labels: Vec<u8>, fs, ft, unknowns, deg, secs, censored| CostRec { labels, fs, ft, unknowns, deg, secs, rss_mb: 0.0, censored };
    let mut t = CostTable::new();
    // bucket A: completed 2 s; bucket B: censored 60 s (a lower bound)
    t.record(rec(vec![2, 2, 2], 30, 40, 9000, 0, 2.0, false));
    t.record(rec(vec![2, 2, 2], 40, 50, 12000, 1, 60.0, true));
    assert_eq!(t.predict(&[2, 2, 2], 1200, 0), Some(2.0));
    assert_eq!(t.predict(&[2, 2, 2], 2000, 1), Some(f64::INFINITY), "censored record must predict +infinity");
    assert_eq!(t.predict_exact(&[2, 2, 2], 40, 50, 12000, 1), Some(f64::INFINITY));
    // Exact::admit: cost table admits, a censored bucket vetoes the label rule, label rule alone
    let ex = Exact { costs: Some(t.clone()), max_label: 2, ..Exact::default() };
    assert_eq!(ex.admit(&[2, 2, 2], 30, 40, 9000, 0), Ok("cost-table"));
    assert!(ex.admit(&[2, 2, 2], 40, 40, 12000, 1).unwrap_err().contains("vetoed"));
    let ex0 = Exact::default();
    assert_eq!(ex0.admit(&[2, 2, 2], 30, 40, 9000, 0), Ok("label-rule"));
    assert!(ex0.admit(&[3, 3], 30, 40, 9000, 0).is_err());
    assert!(ex0.admit(&[2, 2, 2], 50, 50, 9000, 0).is_err(), "fibre product cap");
    assert!(ex0.admit(&[2, 2, 2], 30, 40, 40_001, 0).is_err(), "unknown cap");
    // Search window admission: a censored bucket does not admit a window above the unknown cap
    let lam = [2u8, 2];
    let ar: Arena<F> = Arena::new(&lam, TmplLib::new(None));
    let (a, b) = (vec![e(0, 1), f(0, 1), e(0, 1), f(0, 1)], vec![e(0, 2), f(0, 2)]);
    let (fa, fb) = (krhom::schur::flows(&lam, &a), krhom::schur::flows(&lam, &b));
    let unk = hom_unknowns(&lam, &a, &b, 0);
    let mut se = Search::tiers(0, 6);
    se.max_window_unknowns = unk - 1;
    assert!(!se.window_admitted(&ar, &lam, &a, &b, 0), "above the unknown cap");
    let mut tc = CostTable::new();
    tc.record(rec(lam.to_vec(), fa + 1, fb, unk + 5, 0, 60.0, true)); // same bucket, censored
    se.costs = Some(tc);
    assert!(!se.window_admitted(&ar, &lam, &a, &b, 0), "a censored bucket admitted a window");
    let mut tc = CostTable::new();
    tc.record(rec(lam.to_vec(), fa + 1, fb, unk + 5, 0, 1.0, false)); // same bucket, fast
    se.costs = Some(tc);
    assert!(se.window_admitted(&ar, &lam, &a, &b, 0), "a fast bucket may admit");
    let mut tc = CostTable::new();
    tc.record(rec(lam.to_vec(), fa, fb, unk, 0, 500.0, false)); // this exact shape, slow
    se.costs = Some(tc);
    se.max_window_unknowns = u128::MAX;
    assert!(!se.window_admitted(&ar, &lam, &a, &b, 0), "a measured slow shape must veto");
    // stamps: round trip; another host or solver or a v1 table is rejected
    let dir = cache_dir("b9_costs");
    let path = dir.join("c.txt");
    t.save(&path).unwrap();
    let back = CostTable::load(&path).unwrap();
    assert_eq!(back.len(), 2);
    assert_eq!((back.alg.as_str(), back.host.clone()), (COST_ALG, host_stamp()));
    let text = std::fs::read_to_string(&path).unwrap();
    std::fs::write(&path, text.replacen(&format!("host={}", host_stamp()), "host=other-machine", 1)).unwrap();
    if std::env::var("KRHOM_COSTS_ANY_HOST").is_err() {
        assert!(CostTable::load(&path).unwrap_err().to_string().contains("host stamp"));
    }
    std::fs::write(&path, text.replacen(&format!("alg={}", COST_ALG), "alg=old", 1)).unwrap();
    assert!(CostTable::load(&path).unwrap_err().to_string().contains("solver stamp"));
    std::fs::write(&path, "krhom-costs v1 fields=labels,fs,ft,unknowns,deg,secs,rss_mb,censored\n2,2 1 1 1 0 0.1 0.0 0\n").unwrap();
    assert!(CostTable::load(&path).unwrap_err().to_string().contains("v1"));
    eprintln!("[B9] cost table: censored = +inf, admission rules, stamps ({})", host_stamp());
}

/// review tests finding 4: a cache file with a stale header does not count as cached
#[test]
fn b9_has_hom_validates_header() {
    let dir = cache_dir("b9_has_hom");
    let lab = [2u8, 2];
    let (a, b) = (Window::new(&lab, &[e(0, 1), f(0, 1)]), Window::new(&lab, &[f(0, 1), e(0, 1)]));
    let lib: TmplLib<F> = TmplLib::new(Some(dir.clone()));
    let path = lib.hom_cache_path(&a, &b, 0).unwrap();
    std::fs::write(&path, "krhom-tmpl v1 field=F31 stale\nmaps 0\nend\n").unwrap();
    assert!(!lib.has_hom(&a, &b, 0), "a stale cache file counted as cached");
    let mut lib: TmplLib<F> = TmplLib::new(Some(dir.clone()));
    let ids = lib.hom_basis(&a, &b, 0);
    lib.save().unwrap();
    assert!(lib.has_hom(&a, &b, 0));
    let lib2: TmplLib<F> = TmplLib::new(Some(dir));
    assert!(lib2.has_hom(&a, &b, 0), "a valid cache file not recognized");
    eprintln!("[B9] has_hom: stale header rejected, valid file ({} maps) accepted", ids.len());
}

// ------------------------------------------------------------------ closure with e_p (review F4)

fn ptrim(mut p: Vec<F>) -> Vec<F> {
    while p.last().map_or(false, |c| c.is_zero()) {
        p.pop();
    }
    p
}
fn pmul(a: &[F], b: &[F]) -> Vec<F> {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let mut r = vec![F::zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            r[i + j].fma(x, y);
        }
    }
    ptrim(r)
}
fn psub(a: &[F], b: &[F]) -> Vec<F> {
    let n = a.len().max(b.len());
    let mut r = vec![F::zero(); n];
    for i in 0..n {
        if i < a.len() {
            r[i].add_assign(&a[i]);
        }
        if i < b.len() {
            r[i].sub_assign(&b[i]);
        }
    }
    ptrim(r)
}
fn pdivrem(a: &[F], b: &[F]) -> (Vec<F>, Vec<F>) {
    let b = ptrim(b.to_vec());
    let mut r = ptrim(a.to_vec());
    if r.len() < b.len() {
        return (vec![], r);
    }
    let mut q = vec![F::zero(); r.len() - b.len() + 1];
    let inv = b.last().unwrap().inv();
    while r.len() >= b.len() && !r.is_empty() {
        let k = r.len() - b.len();
        let c = r.last().unwrap().mul(&inv);
        q[k] = c.clone();
        for (i, y) in b.iter().enumerate() {
            let t = c.mul(y);
            r[k + i].sub_assign(&t);
        }
        r = ptrim(r);
    }
    (ptrim(q), r)
}
fn pxgcd(a: &[F], b: &[F]) -> (Vec<F>, Vec<F>, Vec<F>) {
    let (mut r0, mut r1) = (ptrim(a.to_vec()), ptrim(b.to_vec()));
    let (mut s0, mut s1) = (vec![F::one()], vec![]);
    let (mut t0, mut t1) = (vec![], vec![F::one()]);
    while !r1.is_empty() {
        let (q, r) = pdivrem(&r0, &r1);
        let s2 = psub(&s0, &pmul(&q, &s1));
        let t2 = psub(&t0, &pmul(&q, &t1));
        r0 = r1;
        r1 = r;
        s0 = s1;
        s1 = s2;
        t0 = t1;
        t1 = t2;
    }
    (r0, s0, t0)
}

/// e_p = P(a) from an exact End^0(L_p) basis, as an Arena map (SPEC §5.2 steps 2-5, the way M5
/// is to compute it): a = Σ c_i E_i with top scalar 1 at rep(p), minimal polynomial of a at a
/// point on the witness cosets of every table of char L_p, P = the idempotent polynomial of the
/// (t-1)-primary component, e = Σ P_k a^k
fn idempotent_from_end0(ar: &mut Arena<F>, s: &mut Schur, lam: &[u8], lp: ExprId, p: Tid, end0: &[MapId]) -> MapId {
    let w = ar.expr_of(lp).rungs.clone();
    let ch = s.ladder_char(&w);
    let rep = |s: &Schur, t: Tid| -> u64 {
        let tab = s.table(t).clone();
        let cs = krhom::coset::CosetSpace::new(lam, &tab.cols());
        cs.rank(&cs.rep(&tab))
    };
    let cos: Vec<u64> = ch.t.iter().map(|(t, _)| rep(s, *t)).collect();
    let ctop = rep(s, p);
    let pt = points::<F>(ar.k(), 0x1de3, 1, 0);
    let chi = |ar: &Arena<F>, m: MapId| -> F { Eval::new(ar, pt[0].clone()).block(m, ctop).a[0].clone() };
    let cs: Vec<F> = (0..end0.len()).map(|i| F::from_i64(1000 + 7919 * i as i64 % 100_003)).collect();
    let a0 = ar.add(MapNode::Lin(cs.iter().cloned().zip(end0.iter().cloned()).collect()), "a0").unwrap();
    let x = chi(ar, a0);
    assert!(!x.is_zero(), "random End^0 element has top scalar 0");
    let a = ar.add(MapNode::Lin(vec![(x.inv(), a0)]), "a").unwrap();
    let id = ar.add(MapNode::Id(lp), "id").unwrap();
    let mut pows: Vec<MapId> = vec![];
    let mut rows: Vec<Vec<F>> = vec![];
    let mut cur = id;
    let minpoly: Vec<F> = loop {
        let v = Eval::stacked(ar, &[cur], &pt, &cos).remove(0);
        let n = v.len();
        let mut r2 = rows.clone();
        r2.push(v.clone());
        if krhom::dense::rank(&r2, n) < r2.len() {
            let m = rows.len();
            let mat: Vec<Vec<F>> = (0..n).map(|r| (0..m).map(|i| rows[i][r].clone()).collect()).collect();
            let lam_ = krhom::dense::solve(&mat, m, &v).expect("krylov");
            let mut mp: Vec<F> = lam_.iter().map(|l| l.neg()).collect();
            mp.push(F::one());
            break mp;
        }
        rows.push(v);
        pows.push(cur);
        cur = ar.add(MapNode::Comp(a, cur), "a^k").unwrap();
    };
    let tm1 = vec![F::one().neg(), F::one()];
    let (mut q, mut k) = (minpoly, 0);
    loop {
        let (qq, r) = pdivrem(&q, &tm1);
        if r.is_empty() {
            q = qq;
            k += 1;
        } else {
            break;
        }
    }
    assert!(k > 0, "1 is not an eigenvalue of a");
    let mut pk = vec![F::one()];
    for _ in 0..k {
        pk = pmul(&pk, &tm1);
    }
    let (g, _, t) = pxgcd(&pk, &q);
    assert_eq!(g.len(), 1);
    let gi = g[0].inv();
    let poly: Vec<F> = pmul(&t, &q).into_iter().map(|c| c.mul(&gi)).collect();
    // powers beyond the Krylov ones
    while pows.len() < poly.len() {
        let last = *pows.last().unwrap();
        pows.push(ar.add(MapNode::Comp(a, last), "a^k").unwrap());
    }
    let terms: Vec<(F, MapId)> = poly.iter().cloned().zip(pows.iter().cloned()).filter(|(c, _)| !c.is_zero()).collect();
    ar.add(MapNode::Lin(terms), "e_p = P(a)").unwrap()
}

/// a type of (lambda, lambda) whose minimal-flow canonical ladder is decomposable
fn decomposable_type(s: &mut Schur, lam: &[u8]) -> Option<(Tid, Vec<Rung>)> {
    let order = s.canonical_basis(lam).order.clone();
    for p in order {
        let w = canonical(s, lam, p);
        let ch = s.ladder_char(&w);
        if s.expand(&ch).len() > 1 {
            return Some((p, w));
        }
    }
    None
}

/// review F4 / SPEC §4.3 item 7: closure through an idempotent e_p of a decomposable L_p. Without
/// e_p it is Refused; with e_p = P(a) from an exact End^0 basis, Exact and Search certify n maps
/// in the image of e_p, and their number and degrees are c_diag(b_p)
fn closure_with_idempotent(lam: &[u8]) {
    let mut s = Schur::new(lam);
    let Some((p, w)) = decomposable_type(&mut s, lam) else { panic!("no decomposable minimal-flow ladder in {:?}", lam) };
    let mut ar: Arena<F> = Arena::new(lam, TmplLib::new(Some(cache_dir(&format!("cl_ep_{:?}", lam)))));
    let lp = ar.expr(ex(lam, &w));
    let ch = s.ladder_char(&w);
    let parts: Vec<String> = s.expand(&ch).iter().map(|(t, c)| format!("{}{{{}}}", s.table(*t).fmt(), c.fmt_ec())).collect();
    let mut exact = Exact { max_label: 3, max_fibre_product: u64::MAX, max_unknowns: u64::MAX, ..Exact::default() };
    let b = Budget::default();
    for outs in [false, true] {
        let m = HomProvider::<F>::closure(&mut exact, &mut ar, &mut s, lp, p, None, outs, &b).unwrap_err();
        assert_eq!(m[0].kind, MissKind::Refused, "{:?}", m);
        assert!(m[0].what.contains("e_p required"), "{}", m[0].what);
    }
    let t0 = Instant::now();
    let (end0, _) = HomProvider::<F>::end0(&mut exact, &mut ar, &mut s, lp, &b).unwrap_or_else(|m| panic!("End^0: {:?}", m));
    let e = idempotent_from_end0(&mut ar, &mut s, lam, lp, p, &end0);
    check_idempotent(&ar, &mut s, lp, p, e, 77).unwrap_or_else(|m| panic!("e_p: {}", m));
    let n = s.c_diag(p).at1() as usize;
    let mut search = Search::tiers(0, 6);
    for outs in [false, true] {
        let (ms, _) = HomProvider::<F>::closure(&mut exact, &mut ar, &mut s, lp, p, Some(e), outs, &b).unwrap_or_else(|m| panic!("Exact closure with e_p ({}): {:?}", outs, m));
        assert_eq!(ms.len(), n);
        let (ms2, st) = search.closure(&mut ar, &mut s, lp, p, Some(e), outs, &b).unwrap_or_else(|m| panic!("Search closure with e_p ({}): {:?}\n{:?}", outs, m, search.trace));
        assert_eq!(ms2.len(), n);
        eprintln!("[closure e_p] {:?} {} ({}): {} maps, Search tiers {:?}", lam, s.table(p).fmt(), if outs { "outs" } else { "ins" }, n, st.tier_hits);
    }
    // a non-wrapped family of the right size is rejected by certify_closure
    let r0 = ar.expr(ex(lam, &[]));
    let (ms, _) = HomProvider::<F>::closure(&mut exact, &mut ar, &mut s, lp, p, Some(e), false, &b).unwrap();
    let unwrapped: Vec<MapId> = ms.iter().map(|&m| match ar.node(m) { MapNode::Comp(_, g) => *g, _ => m }).collect();
    let mut want: Vec<i32> = ms.iter().map(|&m| ar.info(m).deg).collect();
    want.sort();
    assert!(certify_closure(&ar, lam, r0, lp, Some(e), &unwrapped, &want, false, 5).is_err());
    eprintln!("[closure e_p] {:?} type {} with char L_p = {}: End^0 dim {}, e_p checked, closure certified with e_p (n = {}) in {:.1}s", lam, s.table(p).fmt(), parts.join(" + "), end0.len(), n, t0.elapsed().as_secs_f64());
}

#[test]
fn closure_with_idempotent_313() {
    // (3,1,3): type 2,0,1/0,1,0/1,0,2 has minimal-flow ladder character b_p + b_r (48 flows)
    closure_with_idempotent(&[3, 1, 3]);
}

/// M5c (M5b open problem 3): Exact cost buckets are keyed by call kind (one-sided closure calls
/// vs two-sided pair/window calls) and the exact-shape record is used first, so a censored
/// closure call no longer vetoes pair calls of the same fibre product and degree.
#[test]
fn m5c_exact_admission_by_call_kind() {
    let rec = |labels: Vec<u8>, fs, ft, unknowns, deg, secs, censored| CostRec { labels, fs, ft, unknowns, deg, secs, rss_mb: 0.0, censored };
    assert_eq!((call_kind(1, 4320), call_kind(4320, 1), call_kind(60, 72)), (0, 0, 1));
    let mut t = CostTable::new();
    // a censored whole-ladder closure call Hom^6(R, L_p), fibre product 4 320 (bucket 13)
    t.record(rec(vec![3, 3, 3], 1, 4320, 406_634, 6, 600.0, true));
    // a completed pair call in the same (labels, bucket, degree)
    t.record(rec(vec![3, 3, 3], 60, 72, 50_000, 6, 30.0, false));
    // pre-M5c key (both kinds): the censored closure record closes the bucket
    assert_eq!(t.predict(&[3, 3, 3], 4320, 6), Some(f64::INFINITY));
    assert_eq!(t.predict_kind(&[3, 3, 3], 1, 4320, 6), Some(f64::INFINITY));
    assert_eq!(t.predict_kind(&[3, 3, 3], 64, 70, 6), Some(30.0));
    let ex = Exact { costs: Some(t.clone()), max_label: 3, max_call_s: 600.0, ..Exact::default() };
    // a pair call of another shape in that bucket: admitted by its own kind's bucket
    assert_eq!(ex.admit(&[3, 3, 3], 64, 70, 60_000, 6), Ok("cost-table"));
    // the closure call itself stays refused (censored)
    assert!(ex.admit(&[3, 3, 3], 1, 4320, 406_634, 6).is_err());
    // the exact-shape record decides before the bucket: a fast measured shape in a censored bucket
    let mut t2 = t.clone();
    t2.record(rec(vec![3, 3, 3], 1, 4100, 300_000, 6, 120.0, false));
    let ex2 = Exact { costs: Some(t2), max_label: 3, max_call_s: 600.0, ..Exact::default() };
    assert_eq!(ex2.admit(&[3, 3, 3], 1, 4100, 300_000, 6), Ok("cost-table"));
    assert!(ex2.admit(&[3, 3, 3], 1, 4200, 310_000, 6).is_err(), "unmeasured closure shape: bucket of its kind (inf)");
    // a censored exact shape refuses even if the kind bucket were fast
    let mut t3 = CostTable::new();
    t3.record(rec(vec![3, 3, 3], 64, 64, 60_000, 6, 5.0, false));
    t3.record(rec(vec![3, 3, 3], 60, 70, 70_000, 6, 600.0, true));
    let ex3 = Exact { costs: Some(t3), max_label: 3, max_call_s: 600.0, ..Exact::default() };
    assert!(ex3.admit(&[3, 3, 3], 60, 70, 70_000, 6).is_err());
}
