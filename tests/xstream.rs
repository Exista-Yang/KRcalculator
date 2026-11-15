//! M2S′ (SPEC_3C_GL6_rev1 §2): the memory-first exact 2-strand engine.
//!
//! * streaming (`XOpts::stream`) = full build: same Poincaré polynomials at N = 3..6 and the same
//!   final object multisets on (2,1), (2,2), (3,1), (3,2) words; (3,3) σ^1, σ^2 pruned at N = 3, 4
//!   `#[ignore]`;
//! * the compact component store (`CMat`) round-trips to the same `PMat`, and compact = plain;
//! * rung-by-rung transport = direct transport;
//! * `close_stage(m)` equals a separate scan of σ^m for (2,2), m = 2, 3, and refuses N' > prune_n;
//! * gate G2X pieces on (3,3) `#[ignore]` (cloud): σ^5 with close_stages at N = 3..8 against
//!   Wang 3 N 1, 3 N 2 and the Hopf formula, both primes, χ, reduced base-point check.

use krhom::input::crossing::Gen;
use krhom::algebra::field::{Field, F31, F32b};
use krhom::algebra::laurent::{LPoly, Poincare};
use krhom::algebra::mpoly::Poly;
use krhom::algebra::pmat::PMat;
use krhom::combinat::schur::{decat_chi, decat_scan, ScanOpts, Schur};
use krhom::colored::xengine::*;
use rand::{rngs::StdRng, Rng, SeedableRng};
use std::collections::BTreeMap;

mod common;
use common::{hopf_formula, wang_ref};

fn gens(w: &[i32]) -> Vec<Gen> {
    w.iter().map(|&g| Gen { i: (g.unsigned_abs() - 1) as usize, positive: g > 0 }).collect()
}

fn opts(stream: bool, compact: bool, rung_transport: bool, prune_n: Option<usize>) -> XOpts {
    XOpts { stream, compact, rung_transport, prune_n, verbose: false, ..Default::default() }
}

fn polys<F: Field>(sc: &XScan<F>, ns: &[usize]) -> Vec<BTreeMap<usize, Poincare>> {
    vec![sc.close(ns, false), sc.close(ns, true)]
}

fn compare_builds<F: Field>(labels: &[usize], word: &[i32], ns: &[usize], prune: Option<usize>) {
    let w = gens(word);
    let full = xscan::<F>(labels, &w, &opts(false, false, false, prune)).unwrap();
    for (stream, compact, rt) in [(true, true, true), (true, false, false), (true, true, false), (false, true, false)] {
        let s = xscan::<F>(labels, &w, &opts(stream, compact, rt, prune)).unwrap();
        let tag = format!("{:?} {:?} prune {:?} stream {} compact {} rung transport {}", labels, word, prune, stream, compact, rt);
        assert_eq!(s.final_objects, full.final_objects, "final objects: {}", tag);
        assert_eq!(s.k0, full.k0, "K_0: {}", tag);
        assert_eq!(polys(&s, ns), polys(&full, ns), "Poincaré polynomials: {}", tag);
        for (a, b) in s.stages.iter().zip(&full.stages) {
            assert_eq!((a.unsimplified, a.objects), (b.unsimplified, b.objects), "stage {}: {}", a.step, tag);
            if stream {
                assert!(a.mem.peak_objects <= b.unsimplified, "stage {}: streaming held more objects than the full build: {}", a.step, tag);
            }
        }
    }
}

#[test]
fn streaming_equals_full_build() {
    let ns = [3, 4, 5, 6];
    let cases: Vec<(Vec<usize>, Vec<i32>)> = vec![
        (vec![2, 1], vec![1, 1]),
        (vec![2, 1], vec![1, -1, 1, 1]),
        (vec![2, 1], vec![1, 1, 1, 1]),
        (vec![2, 2], vec![1, 1, 1]),
        (vec![2, 2], vec![1, 1, 1, 1, 1]),
        (vec![2, 2], vec![1, -1, 1, 1, -1]),
        (vec![2, 2], vec![-1, -1, -1, -1]),
        (vec![3, 1], vec![1, 1, 1, 1]),
        (vec![3, 1], vec![-1, 1, -1, -1]),
        (vec![3, 2], vec![1, 1]),
    ];
    for (l, w) in &cases {
        compare_builds::<F31>(l, w, &ns, None);
    }
    compare_builds::<F32b>(&[2, 2], &[1, 1, 1, 1], &ns, None);
    // pruned scans (closed at N' <= prune_n)
    compare_builds::<F31>(&[2, 2], &[1, 1, 1, 1, 1], &[2, 3], Some(3));
    compare_builds::<F31>(&[3, 1], &[1, 1, 1, 1], &[3], Some(3));
}

#[test]
#[ignore]
fn streaming_equals_full_build_32() {
    compare_builds::<F31>(&[3, 2], &[1, 1, 1, 1], &[3, 4, 5, 6], None);
    compare_builds::<F31>(&[3, 2], &[1, -1, 1, 1], &[3, 4, 5], None);
}

#[test]
#[ignore]
fn streaming_equals_full_build_33_pruned() {
    for nn in [3, 4] {
        compare_builds::<F31>(&[3, 3], &[1], &[3, nn], Some(nn));
        compare_builds::<F31>(&[3, 3], &[1, 1], &[3, nn], Some(nn));
    }
}

fn random_pmat<F: Field>(rng: &mut StdRng, nr: usize, nc: usize, nvars: usize) -> PMat<F> {
    let mut m = PMat::zero(nr, nc);
    for j in 0..nc {
        let mut col = vec![];
        for i in 0..nr {
            if rng.gen_bool(0.4) {
                let mut t = vec![];
                let mut mono: u128 = 0;
                for _ in 0..rng.gen_range(1..6) {
                    mono += 1u128 << (8 * rng.gen_range(0..nvars)) as u32;
                    let c = F::from_i64(rng.gen_range(1..i64::from(u32::MAX)));
                    if !c.is_zero() {
                        t.push((mono, c));
                    }
                }
                if !t.is_empty() {
                    col.push((i as u32, Poly { t }));
                }
            }
        }
        m.cols[j] = col;
    }
    m
}

fn same_pmat<F: Field>(a: &PMat<F>, b: &PMat<F>) -> bool {
    a.nrows == b.nrows && a.ncols == b.ncols && a.cols == b.cols
}

#[test]
fn compact_store_round_trip() {
    fn run<F: Field>(seed: u64) {
        let mut rng = StdRng::seed_from_u64(seed);
        for _ in 0..50 {
            let (nr, nc) = (rng.gen_range(0..12), rng.gen_range(0..12));
            let m = random_pmat::<F>(&mut rng, nr, nc, 12);
            let c = CMat::new(m.clone(), true);
            assert!(matches!(c, CMat::Packed { .. }), "{}: not packed", F::name());
            assert!(same_pmat(&c.get(), &m), "{}: compact store does not round-trip", F::name());
            let (e, t, _) = c.size();
            let (e2, t2, _) = pmat_size(&m);
            assert_eq!((e, t), (e2, t2));
            let p = CMat::new(m.clone(), false);
            assert!(matches!(p, CMat::Plain(_)) && same_pmat(&p.get(), &m));
        }
        // extreme coefficients and monomials
        let big = F::from_i64(-1);
        let m = PMat { nrows: 2, ncols: 1, cols: vec![vec![(1, Poly { t: vec![(0, F::one()), (u128::MAX >> 8, big)] })]] };
        assert!(same_pmat(&CMat::new(m.clone(), true).get(), &m));
    }
    run::<F31>(1);
    run::<F32b>(2);
}

#[test]
fn close_stage_equals_separate_scan() {
    let ns = [2, 3, 4];
    let o = XOpts { close_stages: true, verbose: false, ..Default::default() };
    for word in [vec![1, 1, 1, 1], vec![-1, -1, -1, -1], vec![1, -1, 1, 1]] {
        let sc = xscan::<F31>(&[2, 2], &gens(&word), &o).unwrap();
        assert_eq!(sc.closed_stages(), vec![2, 3], "{:?}", word);
        for m in [2, 3] {
            let sep = xscan::<F31>(&[2, 2], &gens(&word[..m]), &XOpts { verbose: false, ..Default::default() }).unwrap();
            for reduced in [false, true] {
                assert_eq!(sc.close_stage(m, &ns, reduced).unwrap(), sep.close(&ns, reduced), "{:?} stage {} reduced {}", word, m, reduced);
            }
        }
        assert_eq!(sc.close_stage(4, &ns, false).unwrap(), sc.close(&ns, false));
        assert!(sc.close_stage(1, &ns, false).is_err());
    }
    // (2,1): odd stages do not close up and are not kept
    let sc = xscan::<F31>(&[2, 1], &gens(&[1, 1, 1, 1]), &o).unwrap();
    assert_eq!(sc.closed_stages(), vec![2]);
    let sep = xscan::<F31>(&[2, 1], &gens(&[1, 1]), &XOpts { verbose: false, ..Default::default() }).unwrap();
    assert_eq!(sc.close_stage(2, &[2, 3, 4], false).unwrap(), sep.close(&[2, 3, 4], false));
    // without close_stages no stage is kept
    let sc = xscan::<F31>(&[2, 2], &gens(&[1, 1, 1]), &XOpts { verbose: false, ..Default::default() }).unwrap();
    assert!(sc.closed_stages().is_empty() && sc.close_stage(2, &[3], false).is_err());
}

#[test]
fn close_stage_refuses_n_above_prune_n() {
    let o = XOpts { close_stages: true, prune_n: Some(3), verbose: false, ..Default::default() };
    let sc = xscan::<F31>(&[2, 2], &gens(&[1, 1, 1, 1]), &o).unwrap();
    let e = sc.close_stage(2, &[3, 4], false).unwrap_err();
    assert!(e.contains("pruned"), "{}", e);
    assert!(sc.close_stage(3, &[4], true).is_err());
    // N' <= prune_n is accepted and equals the unpruned scan
    let un = xscan::<F31>(&[2, 2], &gens(&[1, 1]), &XOpts { verbose: false, ..Default::default() }).unwrap();
    assert_eq!(sc.close_stage(2, &[2, 3], false).unwrap(), un.close(&[2, 3], false));
    // the final stage keeps close's guard (panics above prune_n)
    let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| sc.close_stage(4, &[4], false)));
    assert!(r.is_err() || r.unwrap().is_err());
}

// ------------------------------------------------------------------ gate G2X pieces (cloud)

fn chi_oracle(labels: &[usize], word: &[i32], nn: usize) -> LPoly {
    let lam: Vec<u8> = labels.iter().map(|&a| a as u8).collect();
    let mut s = Schur::new(&lam);
    let res = decat_scan(&mut s, word, &ScanOpts { peak: false, every_stage: false, marked_block: None });
    decat_chi(&mut s, &res, word, nn)
}

/// (3,3) σ^5 with close_stages: stage 2 = Hopf formula, stage 3 = Wang 3 N 1, stage 5 =
/// Wang 3 N 2, χ = decat_chi at every stage, at N = 3..6 (N = 7, 8 where the reference exists).
#[test]
#[ignore]
fn g2x_sigma5_33() {
    fn run<F: Field>() -> BTreeMap<(usize, bool), BTreeMap<usize, Poincare>> {
        let w = [1, 1, 1, 1, 1];
        let ns: Vec<usize> = (3..=8).collect();
        let sc = xscan::<F>(&[3, 3], &gens(&w), &XOpts { close_stages: true, ..Default::default() }).unwrap();
        let mut out = BTreeMap::new();
        for m in [2, 3, 4, 5] {
            for reduced in [false, true] {
                out.insert((m, reduced), sc.close_stage(m, &ns, reduced).unwrap());
            }
            for &nn in &ns {
                assert_eq!(out[&(m, false)][&nn].euler(), chi_oracle(&[3, 3], &w[..m], nn), "χ stage {} N={}", m, nn);
            }
        }
        for m in [3, 5] {
            assert_eq!(sc.close_stage_marked(m, &ns, Some(1)).unwrap(), out[&(m, true)], "base-point check stage {}", m);
            for &nn in &ns {
                // χ(unreduced) = [N choose 3] χ(reduced) for knots
                assert_eq!(out[&(m, false)][&nn].euler(), out[&(m, true)][&nn].euler().mul(&LPoly::qbinom(nn as i64, 3)), "reduced χ stage {} N={}", m, nn);
            }
        }
        for &nn in &ns {
            assert_eq!(out[&(2, false)][&nn], hopf_formula(nn, 3, 3), "Hopf N={}", nn);
            for (m, mm) in [(3, 1), (5, 2)] {
                if let Some(want) = wang_ref().get(&(3, nn, mm)) {
                    assert_eq!(&out[&(m, false)][&nn], want, "T(2,{}) N={}", 2 * mm + 1, nn);
                }
            }
        }
        out
    }
    let a = run::<F31>();
    let b = run::<F32b>();
    assert_eq!(a, b, "the two primes differ");
}
