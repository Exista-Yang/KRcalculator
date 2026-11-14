//! M6.0 acceptance (SYNTHESIS §3.1 and §9.2; R3 V1-V10): the exact small colored engine
//! `xengine.rs` on braids with at most 2 strands.
//!
//! * V1 colored unknots ([N choose a]) and unlinks; reduced unknot = q^0.
//! * V2 colors > N give 0; ∧^N components act as line bundles; dim Φ_N of every final type = MOY
//!   (0 when the canonical ladder has an internal label > N).
//! * V4 Hopf links (2,1), (1,2), (2,2) (and (1,1)), bigraded, N = 2..6, against the closed formula
//!   of R3 V4 (`hopf_ours.py`).
//! * V5 Wang: 2-colored T(2,3) at N = 3..8, 10 (ranks 7, 22, 46, 79, 121, 172, 301), T(2,5) at
//!   N = 3..6 (46 and 106 at N = 4, 5), T(2,7) at N = 3, 4 (78 at N = 4), bigraded against
//!   `wang.py` (tests/data/wang_ref.txt).
//! * V9 full cube: Hopf (2,1), (1,2), (2,2) and (2,2) kinks at N <= 4.
//! * V8 cabling: 2-colored T(2,n), n <= 5, N = 2..6, and (2,1)-colored T(2,2), T(2,4).
//! * V3 duality at N = a + 1: (2,2) at N = 3 = uncolored gl_3 (engine2) on the reversed word,
//!   labelled "via duality (expected)".
//! * χ = `schur::decat_chi` on every case; both primes F31 and F32b; mirror = negative braid;
//!   reduced: χ(unreduced) = [N choose a] χ(reduced) for knots, the Lewark-type relation
//!   [N choose a] P_red − P_unred = (1 + t) Q, Q ≥ 0,
//!   and rank-level reduced oracles: engine2 on (1,1), duality at N = 3 on (2,2), mirror,
//!   cancelling pairs, rotations, reduced unlinks (all asserted).
//! * 3 colors: (3,1) by default (formula, cabling); (3,2) and (3,3) `#[ignore]` (timings in the
//!   test output), against the Hopf formula, Wang and duality at N = a + 1.

use krhom::cclose::exact_type_dims_colored;
use krhom::colored::khr_colored_cable;
use krhom::crossing::{Gen, QConv};
use krhom::cube::{khr_full_cube, BraidLink};
use krhom::engine2::{khr_engine2_closure, Closure};
use krhom::field::{Field, F31, F32b};
use krhom::laurent::{LPoly, Poincare};
use krhom::mpoly::MonoCache;
use krhom::sbim::{hom_space, SBim};
use krhom::schur::{decat_chi, decat_scan, ScanOpts, Schur, Table};
use krhom::tower::EndPolicy;
use krhom::web::Ladder;
use krhom::xengine::*;
use std::any::Any;
use std::collections::{BTreeMap, HashMap};
use std::sync::{Arc, Mutex, OnceLock};

mod common;
use common::{hopf_formula, wang_ref};

// ------------------------------------------------------------------ helpers

fn gens(w: &[i32]) -> Vec<Gen> {
    w.iter().map(|&g| Gen { i: (g.unsigned_abs() - 1) as usize, positive: g > 0 }).collect()
}

type ScanKey = (String, Vec<usize>, Vec<i32>);
type ScanCell = Arc<OnceLock<Arc<dyn Any + Send + Sync>>>;

/// One scan per (field, labels, word), shared by all tests of this binary.
fn scan<F: Field>(labels: &[usize], word: &[i32]) -> Arc<XScan<F>> {
    static CACHE: OnceLock<Mutex<HashMap<ScanKey, ScanCell>>> = OnceLock::new();
    let cell = CACHE.get_or_init(Default::default).lock().unwrap().entry((F::name(), labels.to_vec(), word.to_vec())).or_default().clone();
    let any = cell
        .get_or_init(|| {
            let t0 = std::time::Instant::now();
            let s = xscan::<F>(labels, &gens(word), &XOpts::default()).unwrap_or_else(|e| panic!("xscan {:?} {:?} over {}: {}", labels, word, F::name(), e));
            eprintln!(
                "[scan] {:?} {:?} over {}: {} final objects, {:.2}s (hom_space calls {}, max {:.2}s)",
                labels,
                word,
                F::name(),
                s.data.objs.len(),
                t0.elapsed().as_secs_f64(),
                s.hom_calls,
                s.hom_max_secs
            );
            Arc::new(s) as Arc<dyn Any + Send + Sync>
        })
        .clone();
    any.downcast::<XScan<F>>().unwrap()
}

fn sing<F: Field>(labels: &[usize], word: &[i32], nn: usize, reduced: bool) -> Poincare {
    scan::<F>(labels, word).close(&[nn], reduced).remove(&nn).unwrap()
}

fn t0(p: &LPoly) -> Poincare {
    let mut out = Poincare::default();
    for (&e, &k) in &p.c {
        out.add_term(0, e, k);
    }
    out
}

fn one() -> Poincare {
    let mut p = Poincare::default();
    p.add_term(0, 0, 1);
    p
}

fn mirror(p: &Poincare) -> Poincare {
    let mut out = Poincare::default();
    for (&(t, q), &k) in &p.c {
        out.add_term(-t, -q, k);
    }
    out
}

fn shift(p: &Poincare, dt: i32, dq: i32) -> Poincare {
    let mut out = Poincare::default();
    for (&(t, q), &k) in &p.c {
        out.add_term(t + dt, q + dq, k);
    }
    out
}

fn qbinom(n: usize, k: usize) -> LPoly {
    LPoly::qbinom(n as i64, k as i64)
}

/// χ from the decategorified oracle (M3a).
fn chi_oracle(labels: &[usize], word: &[i32], nn: usize) -> LPoly {
    let lam: Vec<u8> = labels.iter().map(|&a| a as u8).collect();
    let mut s = Schur::new(&lam);
    let res = decat_scan(&mut s, word, &ScanOpts { peak: false, every_stage: false, marked_block: None });
    decat_chi(&mut s, &res, word, nn)
}

fn cube(labels: &[usize], word: &[i32], nn: usize) -> Poincare {
    khr_full_cube::<F31>(&BraidLink { n: nn, labels: labels.to_vec(), word: gens(word) }, QConv::A).poincare
}

/// Lewark-type check [N choose a] P_red − P_unred = (1 + t) Q with Q ≥ 0: returns Q, or an
/// explanation of the failure.
fn lewark_q(unred: &Poincare, red: &Poincare, nn: usize, a: usize) -> Result<Poincare, String> {
    let mut d = Poincare::default();
    for (&(t, q), &k) in &red.c {
        for (&e, &c) in &qbinom(nn, a).c {
            d.add_term(t, q + e, k * c);
        }
    }
    for (&(t, q), &k) in &unred.c {
        d.add_term(t, q, -k);
    }
    // divide by (1 + t) per q, lowest t first
    let mut by_q: BTreeMap<i32, BTreeMap<i32, i128>> = BTreeMap::new();
    for (&(t, q), &k) in &d.c {
        by_q.entry(q).or_default().insert(t, k);
    }
    let mut qq = Poincare::default();
    for (q, mut row) in by_q {
        while let Some((t, k)) = row.pop_first() {
            if k < 0 {
                return Err(format!("negative coefficient {} of Q at t^{} q^{}", k, t, q));
            }
            qq.add_term(t, q, k);
            let e = row.entry(t + 1).or_insert(0);
            *e -= k;
            if *e == 0 {
                row.remove(&(t + 1));
            }
        }
    }
    Ok(qq)
}

// ------------------------------------------------------------------ engine-level checks

#[test]
fn canonical_ladders_are_indecomposable() {
    // Every type used by the scans (intermediate stages and decomposition summands included, not
    // only the final types): End^0(L_p) = k and Hom^e(L_p, L_p) = 0 for every e < 0 (T-h2,
    // e_p = id), and dim End^e(L_p) = the Soergel Hom formula (b_p, b_p) in degree e
    // (`schur::hom_gd`, EKLP Thm 4.16) for e <= 2. Below e_min = min fdeg − max fdeg every entry
    // of a map would have negative degree, so the negative range is complete.
    let cases: Vec<(Vec<usize>, Vec<i32>)> = vec![
        (vec![1, 1], vec![1, 1, 1]),
        (vec![2, 2], vec![1, 1, 1]),
        (vec![2, 2], vec![-1, -1, -1]),
        (vec![2, 2], vec![1, 1, 1, 1, 1]),
        (vec![2, 1], vec![1, 1]),
        (vec![2, 1], vec![1, 1, 1, 1]),
        (vec![1, 2], vec![-1, -1]),
        (vec![3, 1], vec![1, 1]),
        (vec![1, 3], vec![-1, -1]),
    ];
    let mut checked = 0;
    for (labels, word) in cases {
        let sc = scan::<F31>(&labels, &word);
        assert_eq!(sc.scan_ladders.len(), sc.ladders.len(), "{:?} {:?}", labels, word);
        let k: usize = labels.iter().sum();
        let lam8: Vec<u8> = labels.iter().map(|&a| a as u8).collect();
        let mut s = Schur::new(&lam8);
        for (tab, rungs) in &sc.scan_ladders {
            let b = SBim::<F31>::identity(2 * k, &labels).add_rungs(rungs);
            let t = Table::parse(tab);
            let mu = t.cols();
            assert_eq!(b.top_labels, mu.iter().map(|&x| x as usize).collect::<Vec<_>>(), "type {} of {:?}", tab, labels);
            let tid = s.intern(t);
            let bp = s.canonical_basis(&mu).b[&tid].clone();
            let emin = b.fdeg.iter().min().unwrap() - b.fdeg.iter().max().unwrap();
            let ehi = 2;
            let gd = s.hom_gd(&bp, &bp, emin, ehi);
            let mut mc = MonoCache::new(b.ring.clone());
            for e in emin..=ehi {
                let d = hom_space(&b, &b, e, &mut mc).len();
                if e <= 0 {
                    assert_eq!(d, if e == 0 { 1 } else { 0 }, "type {} {:?} of {:?}: dim End^{} = {}", tab, rungs, labels, e, d);
                }
                assert_eq!(d as i128, gd[(e - emin) as usize], "type {} {:?} of {:?}: dim End^{} != Soergel Hom formula", tab, rungs, labels, e);
            }
            checked += 1;
        }
    }
    eprintln!("canonical ladders: {} (scan, type) pairs checked", checked);
}

#[test]
fn identity_decomposition_is_exercised() {
    // the shortcut ι = π = id (L_p ⊗ R literally the canonical ladder of its only summand) asserts
    // that both sides are the same SBim; make sure the default suite reaches it
    for (labels, word) in [(vec![2usize, 2], vec![1]), (vec![2, 2], vec![1, 1, 1]), (vec![2, 1], vec![1, 1])] {
        let sc = scan::<F31>(&labels, &word);
        assert!(sc.identity_decompositions > 0, "{:?} {:?}: no identity decomposition", labels, word);
    }
}

#[test]
fn stages_k0_and_lower_bound() {
    // the K_0 class is asserted inside xscan after every crossing; here: the final class equals
    // schur::decat_scan's, and every stage has at least LB objects (the 2-colored trefoil reaches
    // the bound 11 only up to cancelling pairs: 13 objects)
    for (labels, word) in [(vec![2usize, 2], vec![1, 1, 1]), (vec![2, 1], vec![1, 1]), (vec![2, 2], vec![-1, -1, -1]), (vec![2, 2], vec![1, 1, 1, 1, 1])] {
        let sc = scan::<F31>(&labels, &word);
        let lam: Vec<u8> = labels.iter().map(|&a| a as u8).collect();
        let mut s = Schur::new(&lam);
        let res = decat_scan(&mut s, &word, &ScanOpts { peak: false, every_stage: false, marked_block: None });
        let want: Vec<(String, i32, i128)> = res.k0.iter().map(|&(p, k, c)| (s.table(p).fmt(), k, c)).collect();
        let mut got = sc.k0.clone();
        got.sort();
        let mut want = want;
        want.sort();
        assert_eq!(got, want, "final K_0 class {:?} {:?}", labels, word);
        for st in &sc.stages {
            assert!(st.objects as u128 >= st.lower_bound);
        }
        let last = sc.stages.last().unwrap();
        eprintln!("{:?} {:?}: final objects {} (LB {}), stages {:?}", labels, word, last.objects, last.lower_bound, sc.stages.iter().map(|s| (s.unsimplified, s.objects)).collect::<Vec<_>>());
    }
    let sc = scan::<F31>(&[2, 2], &[1, 1, 1]);
    assert_eq!((sc.stages[2].objects, sc.stages[2].lower_bound), (13, 11));
}

#[test]
fn rungwise_and_direct_decompositions_agree() {
    for (labels, word) in [(vec![2usize, 1], vec![1, 1]), (vec![2, 2], vec![1, 1])] {
        let a = scan::<F31>(&labels, &word);
        let b = xscan::<F31>(&labels, &gens(&word), &XOpts { rungwise: false, ..Default::default() }).unwrap();
        assert_eq!(a.final_objects, b.final_objects, "{:?} {:?}", labels, word);
        for nn in 2..=5 {
            assert_eq!(a.close(&[nn], false), b.close(&[nn], false), "{:?} {:?} N={}", labels, word, nn);
            assert_eq!(a.close(&[nn], true), b.close(&[nn], true), "{:?} {:?} N={} reduced", labels, word, nn);
        }
    }
}

#[test]
fn rejects_unsupported_input() {
    let e = |l: &[usize], w: &[i32]| xscan::<F31>(l, &gens(w), &XOpts::default()).err().unwrap();
    assert!(e(&[1, 1, 1], &[1, 2]).contains("at most 2 strands"));
    assert!(e(&[2, 1], &[1]).contains("do not close up"));
    assert!(e(&[4, 1], &[1, 1]).contains("1..=3"));
    assert!(khr_singular_exact::<F31>(&[1, 1, 1], &gens(&[1, -2]), &[2], false).is_err());
}

// ------------------------------------------------------------------ V1, V2

#[test]
fn v1_colored_unknots_and_unlinks() {
    for a in 1..=3usize {
        for nn in 1..=7usize {
            let want = if a <= nn { t0(&qbinom(nn, a)) } else { Poincare::default() };
            assert_eq!(sing::<F31>(&[a], &[], nn, false), want, "unknot a={} N={}", a, nn);
            if a <= nn {
                assert_eq!(sing::<F31>(&[a], &[], nn, true), one(), "reduced unknot a={} N={}", a, nn);
            }
        }
    }
    // kinks σ_1^{±1} on (a, a): the colored unknot, reduced q^0
    for a in 1..=2usize {
        for w in [1, -1] {
            for nn in 1..=6usize {
                let want = if a <= nn { t0(&qbinom(nn, a)) } else { Poincare::default() };
                assert_eq!(sing::<F31>(&[a, a], &[w], nn, false), want, "kink a={} {} N={}", a, w, nn);
                assert_eq!(sing::<F32b>(&[a, a], &[w], nn, false), want, "kink a={} {} N={} (F32b)", a, w, nn);
                if a <= nn {
                    assert_eq!(sing::<F31>(&[a, a], &[w], nn, true), one(), "reduced kink a={} {} N={}", a, w, nn);
                }
            }
        }
    }
    // unlinks: two strands without crossings
    for (a, b) in [(1usize, 1usize), (2, 1), (2, 2), (3, 1)] {
        for nn in 1..=5usize {
            let want = if a <= nn && b <= nn { t0(&qbinom(nn, a).mul(&qbinom(nn, b))) } else { Poincare::default() };
            assert_eq!(sing::<F31>(&[a, b], &[], nn, false), want, "unlink ({},{}) N={}", a, b, nn);
        }
    }
}

#[test]
fn v2_labels_above_n_and_wedge_n() {
    // a color > N gives 0
    for (labels, word) in [(vec![2usize, 2], vec![1, 1, 1]), (vec![2, 1], vec![1, 1]), (vec![2, 2], vec![1, 1])] {
        assert_eq!(sing::<F31>(&labels, &word, 1, false), Poincare::default(), "{:?} {:?} N=1", labels, word);
    }
    // a component colored N is invisible up to q^{-2 lk a (N - a)} (R3 V2): Hopf (a, N)
    for (labels, nn) in [(vec![1usize, 2], 2usize), (vec![2, 1], 2), (vec![2, 2], 2), (vec![3, 1], 3), (vec![1, 3], 3)] {
        let a = labels.iter().copied().min().unwrap();
        let want = shift(&t0(&qbinom(nn, a)), 0, -2 * (a * (nn - a)) as i32);
        assert_eq!(sing::<F31>(&labels, &[1, 1], nn, false), want, "Hopf {:?} at N={}", labels, nn);
        assert_eq!(sing::<F31>(&labels, &[-1, -1], nn, false), mirror(&want), "negative Hopf {:?} at N={}", labels, nn);
    }
    // the 2-colored trefoil at N = 2 (∧^2 C^2 is a line): q^0
    assert_eq!(sing::<F31>(&[2, 2], &[1, 1, 1], 2, false), one());
    // internal labels > N: dim Φ_N of every final type = MOY of its canonical ladder (0 when the
    // ladder passes through a label > N, e.g. E(2)F(2) on (2,2) at N = 2, 3)
    let mut zero_seen = 0;
    for (labels, word) in [(vec![2usize, 2], vec![1, 1, 1]), (vec![2, 1], vec![1, 1]), (vec![3, 1], vec![1, 1])] {
        let sc = scan::<F31>(&labels, &word);
        for nn in *labels.iter().max().unwrap()..=5 {
            let dims = exact_type_dims_colored(&sc.data, nn, None);
            for (p, rungs) in &sc.type_ladders {
                let lad = Ladder { n: nn, bottom: labels.clone(), rungs: rungs.clone() };
                let moy = if lad.is_valid() { krhom::moy::closure_value(&lad).eval_at_one() } else { 0 };
                if moy == 0 {
                    zero_seen += 1;
                }
                assert_eq!(dims[p] as i128, moy, "type {} {:?} of {:?} at N={}", p, rungs, labels, nn);
            }
        }
    }
    assert!(zero_seen > 0, "no type with an internal label > N was exercised");
}

// ------------------------------------------------------------------ V4, V9

#[test]
fn v4_hopf_links_bigraded() {
    for labels in [[1usize, 1], [2, 1], [1, 2], [2, 2], [3, 1], [1, 3]] {
        for nn in 1..=6usize {
            let want = hopf_formula(nn, labels[0], labels[1]);
            let got = sing::<F31>(&labels, &[1, 1], nn, false);
            assert_eq!(got, want, "Hopf {:?} N={}", labels, nn);
            assert_eq!(sing::<F32b>(&labels, &[1, 1], nn, false), want, "Hopf {:?} N={} (F32b)", labels, nn);
            assert_eq!(sing::<F31>(&labels, &[-1, -1], nn, false), mirror(&want), "negative Hopf {:?} N={}", labels, nn);
            if nn >= 2 && labels == [2, 2] {
                eprintln!("Hopf {:?} N={}: {} (rank {})", labels, nn, got, got.total_rank());
            }
        }
    }
}

#[test]
fn v9_full_cube() {
    let t = std::time::Instant::now();
    for (labels, word, nmax) in [
        (vec![2usize, 1], vec![1, 1], 4usize),
        (vec![1, 2], vec![1, 1], 4),
        (vec![2, 2], vec![1, 1], 4),
        (vec![2, 2], vec![1], 4),
        (vec![2, 2], vec![-1], 4),
        (vec![2, 1], vec![-1, -1], 3),
        (vec![3, 1], vec![1, 1], 5),
        (vec![1, 3], vec![-1, -1], 5),
    ] {
        // (khr_full_cube needs every color <= N)
        for nn in *labels.iter().max().unwrap()..=nmax {
            assert_eq!(sing::<F31>(&labels, &word, nn, false), cube(&labels, &word, nn), "{:?} {:?} N={} vs full cube", labels, word, nn);
        }
    }
    // min(c, d) <= N < largest color: both give 0. khr_full_cube's closure_norm underflowed
    // (usize) here in debug builds before M6.0 made it signed; this guards that one-line change.
    // (N < min(c, d) is out of the cube's range: crossing::crossing_cx computes m(N − m) in usize.)
    for (labels, word) in [(vec![2usize, 1], vec![1, 1]), (vec![1, 2], vec![-1, -1]), (vec![3, 1], vec![1, 1]), (vec![1, 3], vec![-1, -1]), (vec![3, 1], vec![1, 1, 1, 1])] {
        for nn in 1..*labels.iter().max().unwrap() {
            assert_eq!(cube(&labels, &word, nn), Poincare::default(), "full cube {:?} {:?} N={} (color > N)", labels, word, nn);
            assert_eq!(sing::<F31>(&labels, &word, nn, false), Poincare::default(), "{:?} {:?} N={} (color > N)", labels, word, nn);
        }
    }
    eprintln!("v9 full cube comparisons: {:.1}s", t.elapsed().as_secs_f64());
}

// ------------------------------------------------------------------ V5 Wang

#[test]
fn v5_wang_two_colored_torus_knots() {
    let r = wang_ref();
    let ranks = [(3usize, 7i128), (4, 22), (5, 46), (6, 79), (7, 121), (8, 172), (10, 301)];
    for &(nn, rank) in &ranks {
        let p = sing::<F31>(&[2, 2], &[1, 1, 1], nn, false);
        assert_eq!(p.total_rank(), rank, "2-colored trefoil N={}", nn);
        assert_eq!(p, r[&(2, nn, 1)], "2-colored trefoil N={} bigraded vs wang.py", nn);
        eprintln!("2-colored trefoil N={}: rank {}", nn, rank);
    }
    for (nn, rank) in [(3usize, 11i128), (4, 46), (5, 106), (6, 191)] {
        let p = sing::<F31>(&[2, 2], &[1, 1, 1, 1, 1], nn, false);
        assert_eq!(p.total_rank(), rank, "2-colored T(2,5) N={}", nn);
        assert_eq!(p, r[&(2, nn, 2)], "2-colored T(2,5) N={} bigraded vs wang.py", nn);
    }
    for (nn, rank) in [(3usize, 15i128), (4, 78)] {
        let p = sing::<F31>(&[2, 2], &[1, 1, 1, 1, 1, 1, 1], nn, false);
        assert_eq!(p.total_rank(), rank, "2-colored T(2,7) N={}", nn);
        assert_eq!(p, r[&(2, nn, 3)], "2-colored T(2,7) N={} bigraded vs wang.py", nn);
    }
    // uncolored trefoil (k = 1)
    for nn in 2..=4 {
        assert_eq!(sing::<F31>(&[1, 1], &[1, 1, 1], nn, false), r[&(1, nn, 1)], "trefoil N={}", nn);
    }
    // the second prime on the trefoil and T(2,5)
    for nn in [3usize, 4, 5] {
        assert_eq!(sing::<F32b>(&[2, 2], &[1, 1, 1], nn, false), r[&(2, nn, 1)], "F32b trefoil N={}", nn);
        assert_eq!(sing::<F32b>(&[2, 2], &[1, 1, 1, 1, 1], nn, false), r[&(2, nn, 2)], "F32b T(2,5) N={}", nn);
    }
}

// ------------------------------------------------------------------ V8 cabling, V3 duality

#[test]
fn v8_cabling() {
    let t = std::time::Instant::now();
    let ns: Vec<usize> = (2..=6).collect();
    let mut cases: Vec<(Vec<usize>, Vec<i32>)> = (1..=5).map(|n| (vec![2usize, 2], vec![1; n])).collect();
    cases.push((vec![2, 1], vec![1, 1]));
    cases.push((vec![2, 1], vec![1, 1, 1, 1]));
    cases.push((vec![2, 2], vec![-1, -1, -1]));
    for (labels, word) in cases {
        let cab = khr_colored_cable::<F31>(&labels, &gens(&word), &ns, 7, Closure::Tower(EndPolicy::Right)).unwrap_or_else(|e| panic!("cable {:?} {:?}: {}", labels, word, e));
        for &nn in &ns {
            assert_eq!(sing::<F31>(&labels, &word, nn, false), cab.colored[&nn], "{:?} {:?} N={} vs cabling", labels, word, nn);
        }
    }
    eprintln!("v8 cabling comparisons: {:.1}s", t.elapsed().as_secs_f64());
}

#[test]
fn v3_duality_at_n_equals_a_plus_one() {
    // ∧^2 gl_3 of cl(β) = uncolored gl_3 of cl(rev β) (bigraded, no shift for knots, and for the
    // Hopf link with lk = 1: x(2,2) = x'(1,1) = 2 at N = 3). Labelled "via duality (expected)".
    for word in [vec![1, 1, 1], vec![1, 1, 1, 1, 1], vec![1, 1], vec![-1, -1, -1]] {
        let rev: Vec<i32> = word.iter().rev().copied().collect();
        let un = khr_engine2_closure::<F31>(2, &gens(&rev), &[3], &[], 7, false, Closure::Exact).poincare.remove(&3).unwrap();
        let col = sing::<F31>(&[2, 2], &word, 3, false);
        assert_eq!(col, un, "via duality (expected): (2,2) {:?} at N=3 vs uncolored gl_3 of the reversed word", word);
    }
}

// ------------------------------------------------------------------ χ, primes, mirror, reduced

fn standard_cases() -> Vec<(Vec<usize>, Vec<i32>)> {
    vec![
        (vec![2, 2], vec![1]),
        (vec![2, 2], vec![-1]),
        (vec![2, 1], vec![1, 1]),
        (vec![1, 2], vec![1, 1]),
        (vec![2, 1], vec![-1, -1]),
        (vec![2, 2], vec![1, 1]),
        (vec![2, 2], vec![1, 1, 1]),
        (vec![2, 2], vec![-1, -1, -1]),
        (vec![2, 2], vec![1, 1, 1, 1]),
        (vec![2, 2], vec![1, 1, 1, 1, 1]),
        (vec![2, 1], vec![1, 1, 1, 1]),
        (vec![2, 2], vec![1, -1, 1]),
        (vec![3, 1], vec![1, 1]),
        (vec![1, 3], vec![1, 1, 1, 1]),
    ]
}

#[test]
fn euler_characteristic_equals_decat_chi() {
    let mut cases = standard_cases();
    cases.push((vec![2, 2], vec![1; 7]));
    cases.push((vec![2, 2], vec![-1, -1, -1, -1, -1]));
    for (labels, word) in cases {
        for nn in 1..=6usize {
            let p = sing::<F31>(&labels, &word, nn, false);
            assert_eq!(p.euler(), chi_oracle(&labels, &word, nn), "χ {:?} {:?} N={}", labels, word, nn);
        }
    }
}

#[test]
fn two_primes_agree() {
    for (labels, word) in standard_cases() {
        for nn in 2..=6usize {
            for reduced in [false, true] {
                assert_eq!(sing::<F31>(&labels, &word, nn, reduced), sing::<F32b>(&labels, &word, nn, reduced), "{:?} {:?} N={} reduced {}", labels, word, nn, reduced);
            }
        }
        let (a, b) = (scan::<F31>(&labels, &word), scan::<F32b>(&labels, &word));
        assert_eq!(a.final_objects, b.final_objects, "final complexes over two primes {:?} {:?}", labels, word);
    }
}

#[test]
fn mirror_is_negative_braid() {
    for (labels, word) in [(vec![2usize, 2], vec![1, 1, 1]), (vec![2, 1], vec![1, 1]), (vec![2, 2], vec![1, 1]), (vec![2, 2], vec![1, 1, 1, 1, 1]), (vec![2, 1], vec![1, 1, 1, 1])] {
        let neg: Vec<i32> = word.iter().map(|g| -g).collect();
        for nn in 2..=6usize {
            assert_eq!(sing::<F31>(&labels, &neg, nn, false), mirror(&sing::<F31>(&labels, &word, nn, false)), "mirror {:?} {:?} N={}", labels, word, nn);
        }
    }
    // braid relation on 2 strands: σσ^{-1} insertion
    for nn in 2..=5usize {
        assert_eq!(sing::<F31>(&[2, 2], &[1, -1, 1], nn, false), sing::<F31>(&[2, 2], &[1], nn, false), "σσ^-1σ = σ N={}", nn);
    }
}

#[test]
fn reduced_theory() {
    // knots: χ(unreduced) = [N choose a] χ(reduced). Every case: the Lewark-type relation
    // [N choose a] P_red − P_unred = (1 + t) Q with Q ≥ 0 (a = color of the marked block 0).
    // This relation is asserted as a regression check: it rests on Wedrich Cor. 2.11 (freeness of
    // the closure over H*(Gr(a, N))) and a filtration spectral sequence, not on a theorem stated
    // for this closure, and it held on all 41 (case, N) pairs when M6.0 was written.
    let mut report = vec![];
    let mut failures = vec![];
    let mut pairs = 0;
    for (labels, word) in [(vec![2usize, 2], vec![1, 1, 1]), (vec![2, 2], vec![-1, -1, -1]), (vec![2, 2], vec![1, 1, 1, 1, 1]), (vec![1, 1], vec![1, 1, 1]), (vec![2, 2], vec![1]), (vec![2, 1], vec![1, 1]), (vec![1, 2], vec![1, 1]), (vec![2, 2], vec![1, 1])] {
        let knot = word.len() % 2 == 1;
        let a = labels[0];
        for nn in a.max(labels[1])..=6usize {
            let un = sing::<F31>(&labels, &word, nn, false);
            let red = sing::<F31>(&labels, &word, nn, true);
            if knot {
                assert_eq!(un.euler(), red.euler().mul(&qbinom(nn, a)), "χ(unred) = [N a] χ(red): {:?} {:?} N={}", labels, word, nn);
            }
            let lw = lewark_q(&un, &red, nn, a);
            pairs += 1;
            if let Err(e) = &lw {
                failures.push(format!("{:?} {:?} N={}: {}", labels, word, nn, e));
            }
            report.push(format!("{:?} {:?} N={}: reduced {} (rank {}); Lewark-type: {}", labels, word, nn, red, red.total_rank(), match &lw {
                Ok(q) => format!("holds, Q = {}", q),
                Err(e) => format!("FAILS ({})", e),
            }));
        }
    }
    for l in &report {
        eprintln!("{}", l);
    }
    assert_eq!(pairs, 41);
    assert!(failures.is_empty(), "Lewark-type relation [N a] P_red − P_unred = (1 + t) Q, Q ≥ 0, fails:\n{}", failures.join("\n"));
    // the reduced 2-colored trefoil at N = 3 equals the reduced uncolored gl_3 trefoil
    // (via duality (expected); more words in reduced_rank_level_oracles)
    let red_un = khr_engine2_closure::<F31>(2, &gens(&[1, 1, 1]), &[3], &[], 7, true, Closure::Exact).poincare.remove(&3).unwrap();
    assert_eq!(sing::<F31>(&[2, 2], &[1, 1, 1], 3, true), red_un, "via duality (expected): reduced (2,2) trefoil at N=3 vs reduced uncolored gl_3 trefoil");
}

#[test]
fn reduced_rank_level_oracles() {
    // Rank-level checks of the reduced theory (SYNTHESIS §4.2: χ agreement is not evidence about
    // differentials).
    let ns: Vec<usize> = (2..=6).collect();
    // (1,1), lengths 2..7, mixed signs: unreduced and reduced (strand at bottom position 0
    // marked in both engines) against engine2's exact closure
    let words11: Vec<Vec<i32>> = vec![
        vec![1, 1],
        vec![-1, -1],
        vec![1, -1],
        vec![1, 1, 1],
        vec![1, -1, 1],
        vec![-1, -1, -1],
        vec![1, 1, 1, 1],
        vec![1, 1, -1, 1],
        vec![-1, 1, -1, -1],
        vec![1, 1, 1, 1, 1],
        vec![1, -1, 1, 1, 1],
        vec![1, 1, 1, 1, 1, 1],
        vec![1, 1, -1, 1, 1, 1],
        vec![1, 1, 1, 1, 1, 1, 1],
        vec![1, -1, 1, 1, -1, 1, 1],
    ];
    for word in &words11 {
        for reduced in [false, true] {
            let e2 = khr_engine2_closure::<F31>(2, &gens(word), &ns, &[], 7, reduced, Closure::Exact).poincare;
            for &nn in &ns {
                assert_eq!(sing::<F31>(&[1, 1], word, nn, reduced), e2[&nn], "(1,1) {:?} N={} reduced {} vs engine2", word, nn, reduced);
            }
        }
    }
    // (2,2) at N = 3, reduced, against reduced uncolored gl_3 of the reversed word (via duality
    // (expected)): knots of both signs and mixed signs, and the Hopf link (lk = 1, no shift)
    for word in [vec![1, 1, 1], vec![-1, -1, -1], vec![1; 5], vec![-1; 5], vec![1; 7], vec![1, -1, 1], vec![1, 1, -1, 1, 1], vec![1, 1, 1, -1, 1], vec![-1, -1, 1, -1, -1], vec![1, 1]] {
        let rev: Vec<i32> = word.iter().rev().copied().collect();
        let un = khr_engine2_closure::<F31>(2, &gens(&rev), &[3], &[], 7, true, Closure::Exact).poincare.remove(&3).unwrap();
        assert_eq!(sing::<F31>(&[2, 2], &word, 3, true), un, "via duality (expected): reduced (2,2) {:?} at N=3 vs reduced uncolored gl_3 of the reversed word", word);
    }
    // mirror: the reduced homology of the negative braid is the mirror of the reduced homology
    for (labels, word) in [(vec![2usize, 2], vec![1, 1, 1]), (vec![2, 2], vec![1, 1, 1, 1, 1]), (vec![2, 2], vec![1, 1]), (vec![2, 1], vec![1, 1]), (vec![1, 2], vec![1, 1]), (vec![2, 1], vec![1, 1, 1, 1])] {
        let neg: Vec<i32> = word.iter().map(|g| -g).collect();
        for &nn in &ns {
            assert_eq!(sing::<F31>(&labels, &neg, nn, true), mirror(&sing::<F31>(&labels, &word, nn, true)), "reduced mirror {:?} {:?} N={}", labels, word, nn);
        }
    }
    // cancelling pairs and cyclic rotations (labels preserved, so the marked block is the same
    // component): σ^2σ^{-1}σ^2 and its rotations are the trefoil; σσ^{-1}σσ and σσσσ^{-1} are the
    // Hopf link
    for &nn in &ns {
        let tref = sing::<F31>(&[2, 2], &[1, 1, 1], nn, true);
        for w in [vec![1, 1, -1, 1, 1], vec![1, -1, 1, 1, 1], vec![-1, 1, 1, 1, 1], vec![1, 1, 1, 1, -1]] {
            assert_eq!(sing::<F31>(&[2, 2], &w, nn, true), tref, "reduced (2,2) {:?} = trefoil, N={}", w, nn);
        }
        let hopf = sing::<F31>(&[2, 1], &[1, 1], nn, true);
        for w in [vec![1, -1, 1, 1], vec![1, 1, 1, -1]] {
            assert_eq!(sing::<F31>(&[2, 1], &w, nn, true), hopf, "reduced (2,1) {:?} = Hopf, N={}", w, nn);
        }
    }
    // reduced unlinks (no crossings, or a cancelling pair): [N b] for the unmarked color b
    for (a, b) in [(1usize, 1usize), (2, 1), (1, 2), (2, 2), (3, 1)] {
        for nn in a.max(b)..=6usize {
            let want = t0(&qbinom(nn, b));
            assert_eq!(sing::<F31>(&[a, b], &[], nn, true), want, "reduced unlink ({},{}) N={}", a, b, nn);
            assert_eq!(sing::<F31>(&[a, b], &[1, -1], nn, true), want, "reduced (σσ^-1) ({},{}) N={}", a, b, nn);
        }
    }
}

// ------------------------------------------------------------------ CLI

#[test]
fn cli_singular_method_matches_cube_output() {
    let exe = env!("CARGO_BIN_EXE_krhom");
    let run = |method: &str| -> String {
        let out = std::process::Command::new(exe).args(["2,3", "--braid", "1 1", "--strands", "2", "--colors", "2,1", "--colored-method", method]).env("KRHOM_THREADS", "1").output().unwrap();
        assert!(out.status.success(), "krhom --colored-method {} failed: {}", method, String::from_utf8_lossy(&out.stderr));
        String::from_utf8_lossy(&out.stdout).lines().filter(|l| l.starts_with("N =") || l.trim_start().starts_with("Poincare") || l.trim_start().starts_with("total rank") || l.trim_start().starts_with("Euler")).collect::<Vec<_>>().join("\n")
    };
    let s = run("singular");
    assert_eq!(s, run("cube"));
    assert!(s.contains("Poincare polynomial: "), "{}", s);
    // the singular method computes over both primes and says so (SYNTHESIS §2.3, §4.2)
    let full = |extra: &[&str]| -> String {
        let mut args = vec!["2,3", "--braid", "1 1 1", "--strands", "2", "--colors", "2", "--colored-method", "singular"];
        args.extend_from_slice(extra);
        let out = std::process::Command::new(exe).args(&args).env("KRHOM_THREADS", "1").output().unwrap();
        assert!(out.status.success(), "krhom {:?} failed: {}", args, String::from_utf8_lossy(&out.stderr));
        String::from_utf8_lossy(&out.stdout).to_string()
    };
    let two = full(&[]);
    assert!(two.contains("coefficients: Q"), "{}", two);
    let red = full(&["--reduced"]);
    assert!(red.contains("coefficients: Q") && red.contains("t^0q^-4 + t^2q^-8 + t^3q^-12"), "{}", red);
    let one = full(&["--single-prime"]);
    assert!(one.contains("coefficients: F_2147483647 (developer option --single-prime: not certified over Q)"), "{}", one);
    // CALC-1: "1 2" on 3 strands destabilizes to the 1-strand unknot ([3 2]_q, rank 3); the
    // refusals of the 3-4 strand route are in tests/calc.rs
    let lib = std::env::temp_dir().join(format!("krhom_xengine_emptylib_{}", std::process::id()));
    std::fs::create_dir_all(&lib).unwrap();
    let out = std::process::Command::new(exe).args(["3", "--braid", "1 2", "--strands", "3", "--colors", "2", "--colored-method", "singular"]).env("KRHOM_THREADS", "1").env("KRHOM_LIB_DIR", &lib).output().unwrap();
    let so = String::from_utf8_lossy(&out.stdout);
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    assert!(so.contains("simplified braid:") && so.contains("total rank: 3"), "{}", so);
}

// ------------------------------------------------------------------ 3 colors

#[test]
fn three_colored_31() {
    // (3,1): Hopf against the formula (in v4) and cabling (4 cable strands); T(2,4)
    let ns: Vec<usize> = (2..=5).collect();
    for word in [vec![1, 1], vec![1, 1, 1, 1]] {
        let cab = khr_colored_cable::<F31>(&[3, 1], &gens(&word), &ns, 7, Closure::Tower(EndPolicy::Right)).unwrap();
        for &nn in &ns {
            let p = sing::<F31>(&[3, 1], &word, nn, false);
            assert_eq!(p, cab.colored[&nn], "(3,1) {:?} N={} vs cabling", word, nn);
            assert_eq!(p.euler(), chi_oracle(&[3, 1], &word, nn), "(3,1) {:?} N={} χ", word, nn);
        }
    }
}

#[test]
#[ignore]
fn three_colored_32_hopf() {
    // (3,2) Hopf: formula (Wang), duality at N = 4 (∧^3 ⊗ ∧^2 ↔ ∧^1 ⊗ ∧^2 with q^{2(x − x')},
    // x = 2·2, x' = 1·3), χ, two primes
    let t = std::time::Instant::now();
    let sc = scan::<F31>(&[3, 2], &[1, 1]);
    eprintln!("(3,2) Hopf scan {:.1}s (hom_space calls {}, max {:.1}s)", t.elapsed().as_secs_f64(), sc.hom_calls, sc.hom_max_secs);
    for nn in 2..=7usize {
        let p = sing::<F31>(&[3, 2], &[1, 1], nn, false);
        assert_eq!(p, hopf_formula(nn, 3, 2), "(3,2) Hopf N={}", nn);
        assert_eq!(p.euler(), chi_oracle(&[3, 2], &[1, 1], nn));
        eprintln!("(3,2) Hopf N={}: {} (rank {})", nn, p, p.total_rank());
    }
    assert_eq!(sing::<F31>(&[3, 2], &[1, 1], 3, false), cube(&[3, 2], &[1, 1], 3), "(3,2) Hopf N=3 vs full cube");
    let dual = shift(&sing::<F31>(&[1, 2], &[1, 1], 4, false), 0, -2);
    assert_eq!(sing::<F31>(&[3, 2], &[1, 1], 4, false), dual, "via duality (expected): (3,2) Hopf at N=4");
    let t = std::time::Instant::now();
    for nn in 3..=5usize {
        assert_eq!(sing::<F32b>(&[3, 2], &[1, 1], nn, false), hopf_formula(nn, 3, 2), "(3,2) Hopf N={} F32b", nn);
    }
    eprintln!("(3,2) Hopf over F32b: {:.1}s", t.elapsed().as_secs_f64());
}

#[test]
#[ignore]
fn three_colored_33_kink() {
    let t = std::time::Instant::now();
    for w in [1, -1] {
        for nn in 3..=6usize {
            assert_eq!(sing::<F31>(&[3, 3], &[w], nn, false), t0(&qbinom(nn, 3)), "(3,3) kink {} N={}", w, nn);
            assert_eq!(sing::<F31>(&[3, 3], &[w], nn, true), one(), "(3,3) kink {} N={} reduced", w, nn);
        }
    }
    eprintln!("(3,3) kinks: {:.1}s", t.elapsed().as_secs_f64());
}

#[test]
#[ignore]
fn three_colored_33_hopf_and_trefoil() {
    // (3,3) Hopf: formula; trefoil: Wang k = 3 (N = 4, 5, 6) and duality at N = 4 (uncolored gl_4
    // of the reversed word, via duality (expected)); χ
    let t = std::time::Instant::now();
    for nn in 3..=6usize {
        let p = sing::<F31>(&[3, 3], &[1, 1], nn, false);
        assert_eq!(p, hopf_formula(nn, 3, 3), "(3,3) Hopf N={}", nn);
        eprintln!("(3,3) Hopf N={}: rank {}", nn, p.total_rank());
    }
    eprintln!("(3,3) Hopf: {:.1}s", t.elapsed().as_secs_f64());
    let r = wang_ref();
    for nn in 3..=6usize {
        let p = sing::<F31>(&[3, 3], &[1, 1, 1], nn, false);
        assert_eq!(p, r[&(3, nn, 1)], "(3,3) trefoil N={} vs wang.py", nn);
        assert_eq!(p.euler(), chi_oracle(&[3, 3], &[1, 1, 1], nn));
        eprintln!("(3,3) trefoil N={}: {} (rank {})", nn, p, p.total_rank());
    }
    let un = khr_engine2_closure::<F31>(2, &gens(&[1, 1, 1]), &[4], &[], 7, false, Closure::Exact).poincare.remove(&4).unwrap();
    assert_eq!(sing::<F31>(&[3, 3], &[1, 1, 1], 4, false), un, "via duality (expected): (3,3) trefoil at N=4 vs uncolored gl_4");
    eprintln!("(3,3) Hopf + trefoil: {:.1}s", t.elapsed().as_secs_f64());
}
