//! CL1 acceptance (SYNTHESIS §9.1): the colored closure layer `cclose.rs`.
//!
//! 1. τ_N of the exact oracle = `sbim::Tau` on random elements of R^λ.
//! 2. Random small closed ladders (2–3 uprights, K ≤ 6): exact R^λ-bases by graded Nakayama
//!    (with the one-point determinant criterion of C2 M2 as an independent check), Gram by
//!    composition; the exact oracle's graded dimension = `StateSpace::compute` = MOY.
//!    Reduced mode on single types: unreduced = [N choose λ_m] · reduced.
//! 3. Exact oracle = (Δ G) tower adapter on the same cases (unreduced and reduced).
//! 4. On engine2 `ClosureData` (λ = 1^n) the colored oracle = `close_exact` (Conway knot).
//!
//! Extra: colored kinks (one crossing, labels (a, a)) assembled from exact Rickard
//! differentials and exact F matrices: exact oracle = adapter = `khr_full_cube`.

use krhom::colored::cclose::*;
use krhom::input::crossing::{crossing_cx, Gen, QConv};
use krhom::reference::cube::{khr_full_cube, BraidLink};
use krhom::uncolored::engine2::{close_with, prepare_closure_data, Closure};
use krhom::algebra::field::{Field, F31};
use krhom::algebra::laurent::{LPoly, Poincare};
use krhom::combinat::moy::closure_value;
use krhom::algebra::mpoly::MonoCache;
use krhom::soergel::sbim::{bottom_ring, hom_space, SBim, StateSpace, Tau};
use krhom::algebra::symfun::GrTable;
use krhom::uncolored::tower::{EndPolicy, TowerOpts};
use krhom::combinat::web::{random_closed_ladder, Ladder};
use rand::{rngs::StdRng, Rng, SeedableRng};
use std::collections::BTreeMap;

fn t0(p: &Poincare) -> LPoly {
    let mut r = LPoly::zero();
    for (&(t, q), &k) in &p.c {
        assert_eq!(t, 0, "single object closes in homological degree 0 only");
        r.add_term(q, k);
    }
    r
}

fn state_space_poly(lad: &Ladder, nn: usize) -> LPoly {
    let b = SBim::<F31>::identity(nn, &lad.bottom).add_rungs(&lad.rungs);
    let mut mc = MonoCache::new(b.ring.clone());
    let ss = StateSpace::compute(&b, &mut mc);
    let sh = b.shift - lad.bottom.iter().map(|&l| (l * (nn - l)) as i32).sum::<i32>();
    let mut got = LPoly::zero();
    for (d, p) in &ss.pieces {
        got.add_term(d + sh, p.reps.len() as i128);
    }
    got
}

fn rand_point(rng: &mut StdRng, k: usize) -> Vec<F31> {
    (0..k).map(|_| F31::from_i64(rng.gen_range(1..2_000_000_000))).collect()
}

// ------------------------------------------------------------------ ring

#[test]
fn ring_packing_and_weights() {
    let ring = ERing::new(&[2, 1, 3]);
    assert_eq!(ring.nvars, 6);
    assert_eq!(ring.var(2, 3), 5);
    let m = ring.pack(&[1, 0, 2, 0, 0, 1]);
    assert_eq!(ring.weight(m), 1 + 2 + 3);
    assert_eq!(ring.unpack(emono_mul(m, m)), vec![2, 0, 4, 0, 0, 2]);
    // weighted monomial counts of R^(2,2,2) by weight (SYNTHESIS A.4): 3, 9, 19, 39
    let r = ERing::new(&[2, 2, 2]);
    let counts: Vec<usize> = (1..=4).map(|w| r.monos_of_weight(w).len()).collect();
    assert_eq!(counts, vec![3, 9, 19, 39]);
    // 20 variables fit
    let big = ERing::new(&[5, 5, 5, 5]);
    let top = big.pack(&[63; 20]);
    assert_eq!(big.unpack(top), vec![63; 20]);
}

#[test]
#[should_panic(expected = "overflow")]
fn ring_exponent_overflow_is_caught() {
    let ring = ERing::new(&[1, 1]);
    let a = ring.pack(&[40, 1]);
    let _ = emono_mul(a, a);
}

#[test]
#[should_panic(expected = "e-variables")]
fn ring_too_many_variables() {
    let _ = ERing::new(&[5, 5, 5, 5, 1]);
}

// ------------------------------------------------------------------ 1. tau

#[test]
fn grassmannian_tables_equal_grtable_at_c0() {
    // the oracle's H*(Gr(a,N)) (Pieri in Schur coordinates) = symfun::GrTable(a,N) at c = 0
    for nn in 1..=7usize {
        for a in 1..=nn {
            let (basis, tab) = grassmannian_table::<F31>(a, nn);
            let g = GrTable::get(a, nn);
            let gb: Vec<Vec<u32>> = g.basis.iter().map(|al| al.iter().map(|&x| x as u32).collect()).collect();
            assert_eq!(basis, gb, "basis order Gr({},{})", a, nn);
            for j in 0..a {
                for i in 0..basis.len() {
                    let mut want: Vec<(u32, F31)> = g.e_mul[j][i]
                        .iter()
                        .filter_map(|(bi, cp)| {
                            let k: i64 = cp.iter().filter(|(e, _)| e.iter().all(|&z| z == 0)).map(|(_, k)| *k).sum();
                            let c = F31::from_i64(k);
                            (!c.is_zero()).then_some((*bi as u32, c))
                        })
                        .collect();
                    want.sort_by_key(|e| e.0);
                    assert_eq!(tab[j][i], want, "Gr({},{}) e_{} * e^{:?}", a, nn, j + 1, basis[i]);
                }
            }
        }
    }
}

#[test]
fn tau_matches_sbim_tau() {
    let mut rng = StdRng::seed_from_u64(11);
    // sbim::Tau builds symfun::GrTable(a, N), which is slow beyond N = 7 (Gr(3,8): 116 s)
    let cases: Vec<(Vec<u8>, usize)> = vec![(vec![1, 2], 6), (vec![2, 2], 6), (vec![3, 1, 2], 6), (vec![2, 1, 1], 6), (vec![4], 7), (vec![5, 2], 7), (vec![0, 3], 6), (vec![3, 3], 6)];
    let mut checked = 0;
    let mut nonzero = 0;
    for (lam, nmax) in cases {
        let ring = ERing::new(&lam);
        let labels: Vec<usize> = lam.iter().map(|&l| l as usize).collect();
        let (_, voff) = bottom_ring(&labels);
        let lo = *labels.iter().max().unwrap();
        for nn in lo..=nmax {
            let alg = GrAlg::<F31>::new(&lam, nn, None);
            let mut tau = Tau::<F31>::new(nn, &labels, &voff);
            assert_eq!(tau.top_degree(), 2 * alg.top_weight as i32);
            for it in 0..30 {
                // mostly top weight (where τ can be nonzero), some other weights
                let w = if it % 5 == 4 { rng.gen_range(0..=alg.top_weight + 2) } else { alg.top_weight };
                let ms = ring.monos_of_weight(w);
                if ms.is_empty() {
                    continue;
                }
                let mut terms = vec![];
                for _ in 0..rng.gen_range(1..=4) {
                    let m = ms[rng.gen_range(0..ms.len())];
                    terms.push((m, F31::from_i64(rng.gen_range(-50..50))));
                }
                let p = EPoly::from_terms(terms);
                let mine = alg.tau(&p);
                let theirs = tau.tau(&p.to_mpoly(ring.nvars));
                assert_eq!(mine, theirs, "lam {:?} N {} p {:?}", lam, nn, p);
                checked += 1;
                if !mine.is_zero() {
                    nonzero += 1;
                }
            }
        }
    }
    assert!(checked > 500 && nonzero > 200, "checked {} nonzero {}", checked, nonzero);
}

// ------------------------------------------------------------------ 2./3. closed ladders

fn random_ladders(count: usize, seed: u64) -> Vec<Ladder> {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut out: Vec<Ladder> = vec![];
    while out.len() < count {
        let m = rng.gen_range(2..=3usize);
        let bottom: Vec<usize> = (0..m).map(|_| rng.gen_range(1..=3usize)).collect();
        let k: usize = bottom.iter().sum();
        if k > 6 {
            continue;
        }
        let steps = rng.gen_range(1..=3usize);
        let l = random_closed_ladder(&mut rng, k, bottom, steps);
        if l.rungs.is_empty() || l.rungs.len() > 6 {
            continue;
        }
        if SBim::<F31>::identity(l.n, &l.bottom).add_rungs(&l.rungs).fsize() > 40 {
            continue;
        }
        if out.contains(&l) {
            continue;
        }
        out.push(l);
    }
    out
}

fn internal_max(l: &Ladder) -> usize {
    l.levels().iter().flat_map(|v| v.iter().copied()).max().unwrap()
}

/// N range per ladder: from the largest bottom label up to max(internal max, bottom max + 1), ≤ 6.
fn n_range(l: &Ladder) -> std::ops::RangeInclusive<usize> {
    let lo = *l.bottom.iter().max().unwrap();
    let hi = internal_max(l).max(lo + 1).min(6);
    lo.max(1)..=hi
}

/// `StateSpace::compute` (generic `hom_space` up to the top degree) is the expensive reference:
/// on (1,3,2) at N = 5 (top degree 32) it takes about 300 s. The default suite compares with
/// it when the top degree 2·Σ λ_i(N − λ_i) is at most this bound, and with MOY everywhere.
const SS_TOP_MAX: usize = 16;

#[test]
fn ladder_closure_equals_state_space_and_moy() {
    ladder_closure_check(SS_TOP_MAX);
}

#[test]
#[ignore]
fn ladder_closure_equals_state_space_all_n() {
    // heavy: StateSpace at every N of the range (minutes; a few hundred MB)
    ladder_closure_check(usize::MAX);
}

fn ladder_closure_check(ss_top_max: usize) {
    let mut rng = StdRng::seed_from_u64(2026);
    let mut compared_moy = 0;
    let mut compared_ss = 0;
    let mut compared_red = 0;
    let mut probes: Vec<String> = vec![];
    for l in random_ladders(30, 7) {
        let (data, basis) = closure_data_from_ladder::<F31>(&l).unwrap_or_else(|e| panic!("{:?}: {}", l, e));
        // the Nakayama counts equal the identity-coset flow count, and the generators pass the
        // one-point determinant criterion (C2 M2) at two generic points
        assert_eq!(basis.ins.len(), basis.rank, "{:?}", l);
        assert_eq!(basis.outs.len(), basis.rank, "{:?}", l);
        let b = SBim::<F31>::identity(l.n, &l.bottom).add_rungs(&l.rungs);
        let k: usize = l.bottom.iter().sum();
        for _ in 0..2 {
            let x = rand_point(&mut rng, k);
            assert_eq!(det_criterion(&b, &basis, &x), (true, true), "{:?}", l);
        }
        for nn in n_range(&l) {
            let got = t0(&close_exact_colored(&data, nn, None));
            let top: usize = l.bottom.iter().map(|&a| 2 * a * (nn - a)).sum();
            if top <= ss_top_max {
                assert_eq!(got, state_space_poly(&l, nn), "state space {:?} N={}", l, nn);
                compared_ss += 1;
            }
            if internal_max(&l) <= nn {
                let mut lm = l.clone();
                lm.n = nn;
                assert_eq!(got, closure_value(&lm), "MOY {:?} N={}", l, nn);
                compared_moy += 1;
            } else {
                // probe U4: internal labels > N; MOY of such a web is 0 at gl_N
                probes.push(format!("{:?} N={}: {}", l.rungs, nn, if got.is_zero() { "0".to_string() } else { format!("NONZERO {}", got) }));
            }
            // reduced mode on single types: unreduced = [N choose λ_m] · reduced
            for (m, &a) in l.bottom.iter().enumerate() {
                if a == 0 {
                    continue;
                }
                let red = t0(&close_exact_colored(&data, nn, Some(m)));
                assert_eq!(red.mul(&LPoly::qbinom(nn as i64, a as i64)), got, "reduced {:?} N={} marked {}", l, nn, m);
                compared_red += 1;
            }
        }
    }
    eprintln!("compared: {} state spaces, {} MOY values, {} reduced", compared_ss, compared_moy, compared_red);
    for p in &probes {
        eprintln!("internal label > N: {}", p);
    }
    assert!(compared_moy >= 30 && compared_ss >= 30, "moy {} ss {}", compared_moy, compared_ss);
}

#[test]
fn exact_oracle_equals_delta_tower_on_ladders() {
    let mut compared = 0;
    for l in random_ladders(30, 7) {
        let (data, _) = closure_data_from_ladder::<F31>(&l).unwrap();
        for nn in n_range(&l) {
            let mut marks = vec![None];
            marks.extend((0..l.bottom.len()).filter(|&m| l.bottom[m] > 0).map(Some));
            for marked in marks {
                if !l.bottom.iter().enumerate().all(|(i, &a)| Some(i) == marked || g4_certified(a, nn)) {
                    continue;
                }
                let exact = close_exact_colored(&data, nn, marked);
                for end in [EndPolicy::Right, EndPolicy::Left] {
                    let tower = close_delta_tower(&data, nn, marked, TowerOpts::new(end)).unwrap_or_else(|e| panic!("{:?} N={} {:?}: {}", l, nn, marked, e));
                    assert_eq!(tower, exact, "{:?} N={} marked {:?} end {:?}", l, nn, marked, end);
                }
                compared += 1;
            }
        }
    }
    eprintln!("adapter comparisons: {}", compared);
    assert!(compared >= 60, "only {} comparisons", compared);
}

#[test]
fn delta_tower_refuses_uncertified_pairs() {
    // a single 2-labelled upright with an E/F digon through label 3 (closed web on (2,1))
    let mut l = Ladder::new(3, vec![2, 1]);
    l.push(krhom::combinat::web::Rung::e(0, 1)).push(krhom::combinat::web::Rung::f(0, 1));
    let (data, _) = closure_data_from_ladder::<F31>(&l).unwrap();
    // N = 33 is outside every certified range (and beyond the tower's caps)
    assert!(!g4_certified(2, 33));
    let r = close_delta_tower(&data, 33, None, TowerOpts::new(EndPolicy::Right));
    assert!(r.unwrap_err().contains("not certified"));
}

// ------------------------------------------------------------------ 4. engine2 data (λ = 1^n)

const CONWAY: [i32; 11] = [1, 1, 2, -3, 2, 1, -3, -2, -2, -3, -3];

fn gens(w: &[i32]) -> Vec<Gen> {
    w.iter().map(|&k| Gen { i: (k.abs() - 1) as usize, positive: k > 0 }).collect()
}

fn engine2_case(n: usize, w: &[i32], seed: u64, ns: std::ops::RangeInclusive<usize>) {
    let data = prepare_closure_data::<F31>(n, &gens(w), &[], seed);
    let cdata = from_engine2(&data);
    cdata.validate().unwrap();
    for nn in ns {
        for reduced in [false, true] {
            let expect = close_with(&data, nn, reduced, Closure::Exact);
            let got = close_exact_colored(&cdata, nn, if reduced { Some(0) } else { None });
            assert_eq!(got, expect, "{:?} N={} reduced={}", w, nn, reduced);
            assert_eq!(cdata.global_q(nn, if reduced { Some(0) } else { None }), krhom::uncolored::engine2::closure_q_shift(n, &gens(w), nn, reduced));
        }
    }
}

#[test]
fn engine2_small_braids_match_close_exact() {
    engine2_case(2, &[1, 1, 1], 5, 2..=6);
    engine2_case(3, &[1, -2, 1, -2], 5, 2..=6);
    engine2_case(2, &[1, 1], 5, 2..=5);
    engine2_case(3, &[1, 1, -2, -2], 5, 2..=4);
}

#[test]
fn engine2_conway_matches_close_exact() {
    engine2_case(4, &CONWAY, 7, 2..=6);
}

#[test]
fn engine2_conway_type_dims_and_adapter() {
    // chain-level dim Φ_N(B_x) per type = engine2's, and the adapter (Δ = 1 for λ = 1^n, the
    // marked strand a dummy root with cap 0) = engine2's tower
    let data = prepare_closure_data::<F31>(4, &gens(&CONWAY), &[], 7);
    let cdata = from_engine2(&data);
    for nn in 2..=5usize {
        for reduced in [false, true] {
            let marked = if reduced { Some(0) } else { None };
            let mine: BTreeMap<usize, usize> = exact_type_dims_colored(&cdata, nn, marked).into_iter().map(|(x, d)| (x as usize, d)).collect();
            assert_eq!(mine, krhom::uncolored::engine2::exact_type_dims(&data, nn, reduced), "N={} reduced={}", nn, reduced);
            let tower = close_with(&data, nn, reduced, Closure::Tower(EndPolicy::Right));
            assert_eq!(close_delta_tower(&cdata, nn, marked, TowerOpts::new(EndPolicy::Right)).unwrap(), tower, "N={} reduced={}", nn, reduced);
        }
    }
}

/// Conway beyond N = 6 (measured on the laptop, KRHOM_THREADS=2: N = 2..8 in 42 s, 368 MB;
/// N = 8 alone: close_exact 12.4 s, colored oracle 15.0 s). N = 9, 10 (CONWAY_NMAX=10) are a
/// cloud job. Prints per-N timings of both closures.
#[test]
#[ignore]
fn engine2_conway_extended() {
    let t = std::time::Instant::now();
    let data = prepare_closure_data::<F31>(4, &gens(&CONWAY), &[], 7);
    eprintln!("Conway prepare_closure_data {:.2}s", t.elapsed().as_secs_f64());
    let cdata = from_engine2(&data);
    let nmax: usize = std::env::var("CONWAY_NMAX").ok().and_then(|v| v.parse().ok()).unwrap_or(8);
    for nn in 2..=nmax {
        for reduced in [false, true] {
            let t1 = std::time::Instant::now();
            let a = close_with(&data, nn, reduced, Closure::Exact);
            let ta = t1.elapsed().as_secs_f64();
            let t2 = std::time::Instant::now();
            let b = close_exact_colored(&cdata, nn, if reduced { Some(0) } else { None });
            let tb = t2.elapsed().as_secs_f64();
            assert_eq!(a, b, "Conway N={} reduced={}", nn, reduced);
            eprintln!("Conway N={} reduced={}: rank {}, close_exact {:.2}s, close_exact_colored {:.2}s", nn, reduced, a.total_rank(), ta, tb);
        }
    }
}

// ------------------------------------------------------------------ extra: colored kinks

/// Closure data of σ_1^{±1} on two uprights labelled (a, a): objects are the Rickard terms
/// (closed 2-rung ladders), differentials the unique degree-1 maps, F solved exactly.
fn kink_data(a: usize, positive: bool) -> ClosureDataC<F31> {
    let n = 2 * a; // no term is pruned: the data are N-independent
    let cx = crossing_cx(Gen { i: 0, positive }, a, a, n);
    let id = SBim::<F31>::identity(n, &[a, a]);
    let webs: Vec<SBim<F31>> = cx.terms.iter().map(|t| id.add_rungs(&t.rungs)).collect();
    let bases: Vec<LadderBasis<F31>> = webs.iter().map(|w| rlambda_bases(w).unwrap()).collect();
    let mut types = BTreeMap::new();
    let mut objs = vec![];
    for (k, (w, bs)) in webs.iter().zip(bases.iter()).enumerate() {
        let gram = gram_of(bs, w.ring.nvars);
        types.insert(k as u32, TypeBasisC { ins: bs.ins.iter().map(|g| g.0).collect(), outs: bs.outs.iter().map(|g| g.0).collect(), gram });
        let kk = cx.terms[k].k as i32;
        // N-independent part of the term shift: −k (positive), +k (negative); plus the web's MOY shift
        objs.push((k as u32, cx.terms[k].hdeg, w.shift + if positive { -kk } else { kk }));
    }
    let mut mc = MonoCache::new(id.ring.clone());
    let mut edges = vec![];
    for j in 0..webs.len() - 1 {
        let (s, t) = if positive { (j, j + 1) } else { (j + 1, j) };
        let u = 1 + webs[s].shift - webs[t].shift;
        let h = hom_space(&webs[s], &webs[t], u, &mut mc);
        assert_eq!(h.len(), 1);
        let f = solve_f(&h[0], u, &bases[s], &bases[t], &webs[t]).unwrap();
        edges.push((s, t, f));
    }
    let data = ClosureDataC { lam: vec![a as u8, a as u8], objs, types, edges, crossings: vec![(positive, a as u8)] };
    data.validate().unwrap();
    data
}

/// Copy of `cube::apply_local_poly` (private there): substitute the local bottom variables of
/// uprights (i, i+1) by the top operators of P.
fn apply_local_poly(p: &SBim<F31>, i: usize, c: usize, d: usize, poly: &krhom::algebra::mpoly::Poly<F31>, v: &krhom::algebra::pmat::PVec<F31>) -> krhom::algebra::pmat::PVec<F31> {
    let mut acc = krhom::algebra::pmat::PAcc::new();
    for (m, coef) in &poly.t {
        let mut w = v.clone();
        for j in 0..c {
            for _ in 0..krhom::algebra::mpoly::mono_exp(*m, j) {
                w = p.top[i][j].mul_vec(&w);
            }
        }
        for j in 0..d {
            for _ in 0..krhom::algebra::mpoly::mono_exp(*m, c + j) {
                w = p.top[i + 1][j].mul_vec(&w);
            }
        }
        acc.add_scaled_const(&w, coef);
    }
    acc.take()
}

/// The whole cube of resolutions of a colored braid as contract-v2 data: objects are the
/// vertex webs (N-independent: no pruning), edges id ⊗ d ⊗ id with the Koszul sign, F solved
/// exactly over R^λ. Mirrors `cube::khr_full_cube` up to the closure.
fn cube_data(labels: &[usize], word: &[Gen]) -> ClosureDataC<F31> {
    let n: usize = labels.iter().sum::<usize>() * 2;
    struct Loc {
        g: Gen,
        c: usize,
        d: usize,
        terms: Vec<(usize, i32, Vec<krhom::combinat::web::Rung>)>,
        maps: Vec<krhom::algebra::pmat::PMat<F31>>,
        fsizes: Vec<usize>,
    }
    let mut labs = labels.to_vec();
    let mut locs = vec![];
    let mut crossings = vec![];
    for g in word {
        let (c, d) = (labs[g.i], labs[g.i + 1]);
        let cx = crossing_cx(Gen { i: 0, positive: g.positive }, c, d, n);
        assert_eq!(cx.terms.len(), c.min(d) + 1, "no pruning");
        let id = SBim::<F31>::identity(n, &[c, d]);
        let webs: Vec<SBim<F31>> = cx.terms.iter().map(|t| id.add_rungs(&t.rungs)).collect();
        let mut mc = MonoCache::new(id.ring.clone());
        let mut maps = vec![];
        for j in 0..webs.len() - 1 {
            let (s, t) = if g.positive { (&webs[j], &webs[j + 1]) } else { (&webs[j + 1], &webs[j]) };
            let h = hom_space(s, t, 1 + s.shift - t.shift, &mut mc);
            assert_eq!(h.len(), 1);
            maps.push(h.into_iter().next().unwrap());
        }
        let terms = cx.terms.iter().map(|t| (t.k, t.hdeg, t.rungs.iter().map(|r| krhom::combinat::web::Rung { i: g.i, dir: r.dir, r: r.r }).collect())).collect();
        locs.push(Loc { g: *g, c, d, terms, maps, fsizes: webs.iter().map(|w| w.fsize()).collect() });
        crossings.push((g.positive, c.min(d) as u8));
        labs.swap(g.i, g.i + 1);
    }
    assert_eq!(labs, labels, "not closable");
    let nc = locs.len();
    // vertices and prefix bimodules
    let id = SBim::<F31>::identity(n, labels);
    let mut prefix: BTreeMap<Vec<usize>, SBim<F31>> = BTreeMap::from([(vec![], id.clone())]);
    let mut states: Vec<Vec<usize>> = vec![vec![]];
    for l in 0..nc {
        let mut next = vec![];
        for s in &states {
            let p = prefix[s].clone();
            for (ti, t) in locs[l].terms.iter().enumerate() {
                let mut s2 = s.clone();
                s2.push(ti);
                prefix.insert(s2.clone(), p.add_rungs(&t.2));
                next.push(s2);
            }
        }
        states = next;
    }
    let sidx: BTreeMap<Vec<usize>, usize> = states.iter().enumerate().map(|(i, s)| (s.clone(), i)).collect();
    let bases: Vec<LadderBasis<F31>> = states.iter().map(|s| rlambda_bases(&prefix[s]).unwrap()).collect();
    let mut types = BTreeMap::new();
    let mut objs = vec![];
    for (i, s) in states.iter().enumerate() {
        let b = &prefix[s];
        let gram = gram_of(&bases[i], b.ring.nvars);
        types.insert(i as u32, TypeBasisC { ins: bases[i].ins.iter().map(|g| g.0).collect(), outs: bases[i].outs.iter().map(|g| g.0).collect(), gram });
        let h: i32 = (0..nc).map(|l| locs[l].terms[s[l]].1).sum();
        let q: i32 = (0..nc).map(|l| { let k = locs[l].terms[s[l]].0 as i32; if locs[l].g.positive { -k } else { k } }).sum::<i32>() + b.shift;
        objs.push((i as u32, h, q));
    }
    let mut edges = vec![];
    for s in &states {
        for l in 0..nc {
            let loc = &locs[l];
            let kl = s[l];
            let (tk, mi) = if loc.g.positive {
                if kl + 1 >= loc.terms.len() { continue; }
                (kl + 1, kl)
            } else {
                if kl == 0 { continue; }
                (kl - 1, kl - 1)
            };
            let mut t = s.clone();
            t[l] = tk;
            let negate = (0..l).map(|l2| locs[l2].terms[s[l2]].1).sum::<i32>().rem_euclid(2) == 1;
            let p = &prefix[&s[..l].to_vec()];
            let fp = p.fsize();
            let (fs, ft) = (loc.fsizes[kl], loc.fsizes[tk]);
            let (bs, bt) = (&prefix[s], &prefix[&t]);
            let fa = bs.fsize() / (fp * fs);
            // ψ = id_P ⊗ L ⊗ id_above as a full polynomial matrix, column by column
            let mut cols = vec![];
            for beta in 0..bs.fsize() {
                let ba = beta % fa;
                let rest = beta / fa;
                let (bl, bp) = (rest % fs, rest / fs);
                let unit: krhom::algebra::pmat::PVec<F31> = vec![(bp as u32, krhom::algebra::mpoly::Poly::one())];
                let mut acc = krhom::algebra::pmat::PAcc::new();
                for (blp, lpoly) in &loc.maps[mi].cols[bl] {
                    let w = apply_local_poly(p, loc.g.i, loc.c, loc.d, lpoly, &unit);
                    let placed: krhom::algebra::pmat::PVec<F31> = w.into_iter().map(|(x, poly)| (((x as usize * ft + *blp as usize) * fa + ba) as u32, poly)).collect();
                    acc.add_scaled_const(&placed, &if negate { F31::one().neg() } else { F31::one() });
                }
                cols.push(acc.take());
            }
            let psi = krhom::algebra::pmat::PMat { nrows: bt.fsize(), ncols: bs.fsize(), cols };
            let (i, j) = (sidx[s], sidx[&t]);
            let deg = 1 + bs.shift - bt.shift;
            let f = solve_f(&psi, deg, &bases[i], &bases[j], bt).unwrap_or_else(|e| panic!("edge {:?} -> {:?}: {}", s, t, e));
            edges.push((i, j, f));
        }
    }
    let data = ClosureDataC { lam: labels.iter().map(|&l| l as u8).collect(), objs, types, edges, crossings };
    data.validate().unwrap();
    data
}

fn hopf_case(labels: &[usize], ns: std::ops::RangeInclusive<usize>) {
    let word = vec![Gen { i: 0, positive: true }, Gen { i: 0, positive: true }];
    let data = cube_data(labels, &word);
    for nn in ns {
        let cube = khr_full_cube::<F31>(&BraidLink { n: nn, labels: labels.to_vec(), word: word.clone() }, QConv::A).poincare;
        let exact = close_exact_colored(&data, nn, None);
        assert_eq!(exact, cube, "Hopf {:?} N={}", labels, nn);
        if labels.iter().all(|&a| g4_certified(a, nn)) {
            assert_eq!(close_delta_tower(&data, nn, None, TowerOpts::new(EndPolicy::Right)).unwrap(), exact, "adapter Hopf {:?} N={}", labels, nn);
            for m in 0..2 {
                assert_eq!(
                    close_delta_tower(&data, nn, Some(m), TowerOpts::new(EndPolicy::Left)).unwrap(),
                    close_exact_colored(&data, nn, Some(m)),
                    "adapter reduced Hopf {:?} N={} marked {}",
                    labels,
                    nn,
                    m
                );
            }
        }
        eprintln!("Hopf {:?} N={}: {}", labels, nn, exact);
    }
}

#[test]
fn colored_hopf_21_matches_full_cube() {
    hopf_case(&[2, 1], 2..=4);
}

#[test]
fn colored_hopf_22_matches_full_cube() {
    hopf_case(&[2, 2], 2..=4);
}

#[test]
fn colored_hopf_euler_characteristics() {
    // SYNTHESIS Appendix A.3 (the `chi` binary, cube of resolutions + MOY)
    let word = vec![Gen { i: 0, positive: true }, Gen { i: 0, positive: true }];
    let lp = |terms: &[(i32, i128)]| -> LPoly {
        let mut p = LPoly::zero();
        for &(e, k) in terms {
            p.add_term(e, k);
        }
        p
    };
    let d21 = cube_data(&[2, 1], &word);
    assert_eq!(close_exact_colored(&d21, 3, None).euler(), lp(&[(0, 1), (-2, 2), (-4, 2), (-6, 2), (-8, 1), (-10, 1)]));
    let d22 = cube_data(&[2, 2], &word);
    assert_eq!(
        close_exact_colored(&d22, 4, None).euler(),
        lp(&[(0, 1), (-2, 1), (-4, 3), (-6, 4), (-8, 6), (-10, 6), (-12, 6), (-14, 4), (-16, 3), (-18, 1), (-20, 1)])
    );
}

#[test]
fn colored_kinks_match_full_cube() {
    for a in 1..=2usize {
        for positive in [true, false] {
            let data = kink_data(a, positive);
            for nn in a..=4usize {
                let cube = khr_full_cube::<F31>(&BraidLink { n: nn, labels: vec![a, a], word: vec![Gen { i: 0, positive }] }, QConv::A).poincare;
                let exact = close_exact_colored(&data, nn, None);
                assert_eq!(exact, cube, "kink a={} positive={} N={}", a, positive, nn);
                // the colored unknot: [N choose a] at t^0
                let mut unknot = Poincare::default();
                for (&e, &k) in &LPoly::qbinom(nn as i64, a as i64).c {
                    unknot.add_term(0, e, k);
                }
                assert_eq!(exact, unknot, "kink a={} positive={} N={} is not the colored unknot", a, positive, nn);
                // reduced: q^0
                for m in 0..2 {
                    let red = close_exact_colored(&data, nn, Some(m));
                    let mut one = Poincare::default();
                    one.add_term(0, 0, 1);
                    assert_eq!(red, one, "reduced kink a={} positive={} N={} marked {}", a, positive, nn, m);
                }
                if g4_certified(a, nn) {
                    let tw = close_delta_tower(&data, nn, None, TowerOpts::new(EndPolicy::Right)).unwrap();
                    assert_eq!(tw, exact, "adapter kink a={} positive={} N={}", a, positive, nn);
                }
            }
        }
    }
}

// ------------------------------------------------------------------ review of CL1: gate G4, guards, errors

/// One type with a single in- and out-generator of degree 0 and G = (1) (the identity web on
/// λ), one object.
fn unit_data(lam: &[u8]) -> ClosureDataC<F31> {
    let mut types = BTreeMap::new();
    types.insert(0u32, TypeBasisC { ins: vec![0], outs: vec![0], gram: vec![vec![EPoly::one()]] });
    let d = ClosureDataC { lam: lam.to_vec(), objs: vec![(0, 0, 0)], types, edges: vec![], crossings: vec![] };
    d.validate().unwrap();
    d
}

#[test]
fn g4_certificate_holds_for_every_listed_pair() {
    // G4_CERTIFIED is hand-maintained: re-run the library certificate (A), (B) on every pair
    assert_eq!(G4_CERTIFIED.len(), 30);
    for &(a, nn) in G4_CERTIFIED {
        let r = g4::check(a as usize, nn as usize, false, false).unwrap();
        assert!(r.certified(), "G4_CERTIFIED lists ({}, {}) but the certificate fails: {:?}", a, nn, r);
        assert!(g4_certified(a as usize, nn as usize));
    }
    // C2's direct cross-checks (C), (D) on small pairs
    for (a, nn) in [(2usize, 4usize), (3, 5)] {
        let r = g4::check(a, nn, true, true).unwrap();
        assert_eq!(r.direct_dim, Some(g4::fact(nn) / g4::fact(nn - a)), "(C) for ({}, {})", a, nn);
        let want: BTreeMap<usize, usize> = g4::pi_a(a).into_iter().enumerate().filter(|&(_, c)| c > 0).map(|(k, c)| (k + a * (a - 1) / 2, c)).collect();
        assert_eq!(r.gens, Some(want), "(D) for ({}, {})", a, nn);
    }
    assert!(g4::check(9, 10, false, false).is_err());
    assert!(g4::check(3, 2, false, false).is_err());
    // N <= 20: the count N!/(N−a)! is computed in usize (for N >= 21 it wrapped around and
    // reported a false "NOT certified"); 20 itself works
    assert!(g4::check(2, 20, false, false).unwrap().certified());
    for (a, nn) in [(2usize, 21usize), (3, 21), (2, 22), (2, 64)] {
        let e = g4::check(a, nn, false, false).unwrap_err();
        assert!(e.contains("N <= 20"), "({}, {}): {}", a, nn, e);
    }
    // the front end reports it as a usage error (developer binary: cargo test --features dev)
    #[cfg(feature = "dev")]
    {
        let out = std::process::Command::new(env!("CARGO_BIN_EXE_g4_free")).arg("2:21").output().expect("cannot run g4_free");
        assert_eq!(out.status.code(), Some(2), "stdout: {}", String::from_utf8_lossy(&out.stdout));
        assert!(String::from_utf8_lossy(&out.stderr).contains("N <= 20"));
    }
}

#[test]
fn delta_tower_refuses_realistic_uncertified_pairs() {
    // a block of label 6 at N = 10: outside G4_CERTIFIED (a <= 5), inside the tower's N <= 32
    let d6 = unit_data(&[6]);
    assert!(!g4_certified(6, 10));
    let e = close_delta_tower(&d6, 10, None, TowerOpts::new(EndPolicy::Right)).unwrap_err();
    assert!(e.contains("not certified") && e.contains("(6, 10)"), "{}", e);
    // the marked block needs no certificate
    assert_eq!(close_delta_tower(&d6, 10, Some(0), TowerOpts::new(EndPolicy::Right)).unwrap(), close_exact_colored(&d6, 10, Some(0)));
    // label 2 at N = 11 (certified only up to N = 10): the closed (2,1) ladder with an E/F digon
    let mut l = Ladder::new(3, vec![2, 1]);
    l.push(krhom::combinat::web::Rung::e(0, 1)).push(krhom::combinat::web::Rung::f(0, 1));
    let (data, _) = closure_data_from_ladder::<F31>(&l).unwrap();
    assert!(!g4_certified(2, 11) && g4_certified(2, 10));
    let e = close_delta_tower(&data, 11, None, TowerOpts::new(EndPolicy::Right)).unwrap_err();
    assert!(e.contains("not certified") && e.contains("(2, 11)"), "{}", e);
    assert_eq!(close_delta_tower(&data, 10, None, TowerOpts::new(EndPolicy::Right)).unwrap(), close_exact_colored(&data, 10, None));
    // with the 2-labelled block marked only label 1 remains
    assert_eq!(close_delta_tower(&data, 11, Some(0), TowerOpts::new(EndPolicy::Left)).unwrap(), close_exact_colored(&data, 11, Some(0)));
}

/// The default cap in bytes: the guard tests use it explicitly, so they do not depend on
/// KRHOM_MEM_CAP_MB (see `mem_cap_reads_the_environment` for the variable itself).
fn default_cap() -> u128 {
    (DEFAULT_MEM_CAP_MB as u128) << 20
}

#[test]
fn delta_tower_resource_guard() {
    // the closed ladder [3,3] with rungs F2 E2 (random_ladders seed 99): at N = 7 the adapter
    // passed 1.2 GB RSS in the review of CL1, while the exact oracle takes 0.04 s
    let cap = default_cap();
    let mut l = Ladder::new(6, vec![3, 3]);
    l.push(krhom::combinat::web::Rung::f(0, 2)).push(krhom::combinat::web::Rung::e(0, 2));
    let (data, _) = closure_data_from_ladder::<F31>(&l).unwrap();
    let est = delta_tower_estimate(&data, 7, None);
    assert_eq!((est.roots, est.root_dim, est.max_ins), (6, 117_649, 3));
    assert!(est.bytes > cap, "estimate {} MB", est.bytes >> 20);
    let t = std::time::Instant::now();
    let e = close_delta_tower_capped(&data, 7, None, TowerOpts::new(EndPolicy::Right), cap).unwrap_err();
    assert!(e.contains("memory cap"), "{}", e);
    assert!(t.elapsed().as_secs_f64() < 1.0, "the guard must refuse before the expansion");
    let _ = try_close_exact_colored_capped(&data, 7, None, cap).unwrap(); // the oracle is cheap here
    // N = 5 (the top of its test range) is admitted and agrees with the oracle
    let est5 = delta_tower_estimate(&data, 5, None);
    assert!(est5.bytes <= cap, "estimate {} MB", est5.bytes >> 20);
    assert_eq!(close_delta_tower_capped(&data, 5, None, TowerOpts::new(EndPolicy::Right), cap).unwrap(), try_close_exact_colored_capped(&data, 5, None, cap).unwrap());
    // the cap is configurable
    let e = close_delta_tower_capped(&data, 5, None, TowerOpts::new(EndPolicy::Right), 1 << 20).unwrap_err();
    assert!(e.contains("memory cap"), "{}", e);
    // reduced mode drops the marked block's roots
    assert_eq!(delta_tower_estimate(&data, 7, Some(0)).roots, 3);
    // the review of the CL1 fixes: [1,1,3] with rungs E11 E01 F01 F11 E11 F11 (24 ins) took
    // 642 MB (N = 5) and 622 MB (N = 7, block 0 marked) above its data preparation, about 1.8x
    // the former estimate, inside the default cap. The heuristic bound is now at least twice
    // the measured peak (in fact 6-10x), so both calls are refused before any work.
    let mut l = Ladder::new(5, vec![1, 1, 3]);
    for (e, i) in [(true, 1), (true, 0), (false, 0), (false, 1), (true, 1), (false, 1)] {
        l.push(if e { krhom::combinat::web::Rung::e(i, 1) } else { krhom::combinat::web::Rung::f(i, 1) });
    }
    let (data, _) = closure_data_from_ladder::<F31>(&l).unwrap();
    for (nn, marked, measured_mb) in [(5usize, None, 642u128), (7, Some(0), 622)] {
        let est = delta_tower_estimate(&data, nn, marked);
        assert_eq!(est.max_ins_outs, 576);
        assert!(est.bytes >= 2 * (measured_mb << 20), "N = {} marked {:?}: estimate {} MB", nn, marked, est.bytes >> 20);
        let e = close_delta_tower_capped(&data, nn, marked, TowerOpts::new(EndPolicy::Right), cap).unwrap_err();
        assert!(e.contains("memory cap"), "{}", e);
    }
}

#[test]
fn delta_tower_refuses_degree_above_63() {
    // λ = (5,5) at N = 10 (certified): one type with G = e_{1,5}^5 e_{2,5}^5. Δ·G has total
    // degree 20 + 50 = 70 in the 10 roots with every exponent <= 9 = N − 1, so it survives the
    // truncation; tower.rs asserted ("Gram entry of degree 70 > 63") before the pre-check.
    let ring = ERing::new(&[5, 5]);
    let mut ex = vec![0u32; ring.nvars];
    ex[ring.var(0, 5)] = 5;
    ex[ring.var(1, 5)] = 5;
    let mut types = BTreeMap::new();
    types.insert(0u32, TypeBasisC { ins: vec![0], outs: vec![100], gram: vec![vec![EPoly::mono(ring.pack(&ex), F31::one())]] });
    let data = ClosureDataC { lam: vec![5, 5], objs: vec![(0, 0, 0)], types, edges: vec![], crossings: vec![] };
    data.validate().unwrap();
    assert!(g4_certified(5, 10));
    // at the default cap the size guard refuses first (N^roots = 10^10)
    let e = close_delta_tower_capped(&data, 10, None, TowerOpts::new(EndPolicy::Right), default_cap()).unwrap_err();
    assert!(e.contains("memory cap"), "{}", e);
    // without the size guard, the degree pre-check refuses before the tower
    let e = close_delta_tower_capped(&data, 10, None, TowerOpts::new(EndPolicy::Right), u128::MAX).unwrap_err();
    assert!(e.contains("total degree 70"), "{}", e);
}

#[test]
fn exact_oracle_returns_errors_instead_of_panicking() {
    // marked block of label 0 or out of range: Err (GrAlg::new asserted, close_exact_colored panicked)
    let d = unit_data(&[2, 0]);
    assert!(GrAlg::<F31>::try_new(&[2, 0], 3, Some(1)).err().unwrap().contains("marked block 1 invalid"));
    assert!(GrAlg::<F31>::try_new(&[2, 0], 3, Some(2)).is_err());
    assert!(GrAlg::<F31>::try_new(&[4, 1], 3, None).is_err());
    // H*(Gr(1, 300)) has top weight 299 > 255 (an assert in BlockAlg before)
    assert!(GrAlg::<F31>::try_new(&[1], 300, None).err().unwrap().contains("top weight"));
    assert!(try_close_exact_colored(&unit_data(&[1]), 300, None).unwrap_err().contains("top weight"));
    assert!(try_close_exact_colored(&d, 3, Some(1)).unwrap_err().contains("marked block 1 invalid"));
    assert!(try_exact_type_dims_colored(&d, 3, Some(5)).unwrap_err().contains("marked block 5 invalid"));
    assert!(close_delta_tower(&d, 3, Some(1), TowerOpts::new(EndPolicy::Right)).unwrap_err().contains("marked block 1 invalid"));
    let mut one = Poincare::default();
    one.add_term(0, 0, 1);
    assert_eq!(try_close_exact_colored(&d, 3, Some(0)).unwrap(), one);
    // exact_type_dims_colored validates the data (it skipped validate() before)
    let mut bad = unit_data(&[2, 1]);
    bad.types.get_mut(&0).unwrap().gram[0].push(EPoly::one());
    assert!(try_exact_type_dims_colored(&bad, 3, None).unwrap_err().contains("invalid closure data"));
    assert!(try_close_exact_colored(&bad, 3, None).unwrap_err().contains("invalid closure data"));
}

#[test]
#[should_panic(expected = "invalid closure data")]
fn exact_type_dims_colored_validates_its_data() {
    let mut bad = unit_data(&[2, 1]);
    bad.types.get_mut(&0).unwrap().gram[0].push(EPoly::one());
    let _ = exact_type_dims_colored(&bad, 3, None);
}

#[test]
fn exact_oracle_memory_cap() {
    // λ = (5,5,5) at N = 10: dim A_λ = 252^3 ≈ 1.6e7, by_weight alone about 0.5 GB
    let cap = default_cap();
    let d = unit_data(&[5, 5, 5]);
    let est = exact_mem_estimate(&d, 10, None);
    assert!(est > cap, "estimate {} MB", est >> 20);
    let t = std::time::Instant::now();
    assert!(try_close_exact_colored_capped(&d, 10, None, cap).unwrap_err().contains("memory cap"));
    assert!(try_exact_type_dims_colored_capped(&d, 10, None, cap).unwrap_err().contains("memory cap"));
    assert!(GrAlg::<F31>::try_new_capped(&[5, 5, 5], 10, None, cap).err().unwrap().contains("memory cap"));
    assert!(t.elapsed().as_secs_f64() < 1.0, "refused before allocating");
    // reduced mode (one block marked) fits: 252^2 box elements
    assert!(exact_mem_estimate(&d, 10, Some(0)) < cap);
    // the cap is configurable: refused under a tiny cap, accepted under a large one
    let d2 = unit_data(&[2, 2]);
    assert!(try_close_exact_colored_capped(&d2, 4, None, 1 << 10).unwrap_err().contains("memory cap"));
    assert!(GrAlg::<F31>::try_new_capped(&[2, 2], 4, None, 1 << 10).is_err());
    let p = try_close_exact_colored_capped(&d2, 4, None, u128::MAX).unwrap();
    assert_eq!(p, close_exact_colored(&d2, 4, None));
    assert_eq!(p.total_rank(), 36); // [4 choose 2]^2: two 2-colored unknots at gl_4
}

/// Helper of `mem_cap_reads_the_environment`: checks `mem_cap_bytes` only in the child
/// processes started by that test (KRHOM_TEST_EXPECT_CAP_MB set); a no-op otherwise.
#[test]
fn mem_cap_env_child() {
    if let Ok(want) = std::env::var("KRHOM_TEST_EXPECT_CAP_MB") {
        assert_eq!(mem_cap_bytes(), want.parse::<u128>().unwrap() << 20, "KRHOM_MEM_CAP_MB = {:?}", std::env::var("KRHOM_MEM_CAP_MB").ok());
    }
}

#[test]
fn mem_cap_reads_the_environment() {
    // KRHOM_MEM_CAP_MB, in child processes of this test binary (this process's environment is
    // left alone; the guard tests above use explicit caps); unparsable values give the default
    let exe = std::env::current_exe().expect("test binary path");
    for (val, want) in [(Some("100"), 100u64), (Some(" 4096 "), 4096), (Some("junk"), DEFAULT_MEM_CAP_MB), (None, DEFAULT_MEM_CAP_MB)] {
        let mut c = std::process::Command::new(&exe);
        c.args(["--exact", "mem_cap_env_child", "--test-threads", "1"]).env("KRHOM_TEST_EXPECT_CAP_MB", want.to_string());
        match val {
            Some(v) => c.env("KRHOM_MEM_CAP_MB", v),
            None => c.env_remove("KRHOM_MEM_CAP_MB"),
        };
        let out = c.output().expect("cannot run the test binary");
        let so = String::from_utf8_lossy(&out.stdout);
        assert!(out.status.success() && so.contains("1 passed"), "KRHOM_MEM_CAP_MB = {:?}:\n{}\n{}", val, so, String::from_utf8_lossy(&out.stderr));
    }
}
