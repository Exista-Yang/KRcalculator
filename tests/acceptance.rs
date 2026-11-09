//! Acceptance tests for the Soergel-bimodule + universal-construction model.

use krhom::crossing::{crossing_cx, Gen, QConv};
use krhom::cube::{khr_full_cube, BraidLink};
use krhom::field::{F31, Q};
use krhom::laurent::{LPoly, Poincare};
use krhom::moy::closure_value;
use krhom::mpoly::MonoCache;
use krhom::sbim::{hom_space, SBim, StateSpace};
use krhom::web::{random_closed_ladder, Dir, Ladder, Rung};
use rand::{rngs::StdRng, Rng, SeedableRng};

fn state_space_poly<F: krhom::field::Field>(lad: &Ladder) -> LPoly {
    let b = SBim::<F>::identity(lad.n, &lad.bottom).add_rungs(&lad.rungs);
    let mut mc = MonoCache::new(b.ring.clone());
    let ss = StateSpace::compute(&b, &mut mc);
    let sh = b.shift - lad.bottom.iter().map(|&l| (l * (lad.n - l)) as i32).sum::<i32>();
    let mut got = LPoly::zero();
    for (d, p) in &ss.pieces {
        got.add_term(d + sh, p.reps.len() as i128);
    }
    got
}

fn lad(n: usize, bottom: &[usize], rungs: &[(char, usize, usize)]) -> Ladder {
    let mut l = Ladder::new(n, bottom.to_vec());
    for &(d, i, r) in rungs {
        l.push(Rung { i, dir: if d == 'E' { Dir::E } else { Dir::F }, r });
    }
    l
}

#[test]
fn state_space_equals_moy_known_counterexamples() {
    // webs where truncated edge-ring models fail
    let cases = vec![
        lad(3, &[2, 1], &[('F', 0, 1), ('E', 0, 1)]),
        lad(3, &[1, 2, 0], &[('E', 0, 1), ('F', 0, 1), ('F', 1, 1), ('E', 0, 1), ('F', 0, 1), ('E', 1, 1)]),
        lad(3, &[3, 2, 1], &[('F', 1, 1), ('E', 1, 1), ('F', 1, 1), ('F', 0, 1), ('E', 0, 1), ('E', 1, 1)]),
    ];
    for l in cases {
        assert_eq!(state_space_poly::<F31>(&l), closure_value(&l), "{:?}", l);
    }
}

#[test]
fn state_space_equals_moy_random() {
    let mut rng = StdRng::seed_from_u64(42);
    for &(n, m, steps, count) in &[(2usize, 3usize, 4usize, 15usize), (3, 3, 3, 12), (4, 2, 4, 10), (4, 3, 2, 6)] {
        for _ in 0..count {
            let bottom: Vec<usize> = (0..m).map(|_| rng.gen_range(0..=n)).collect();
            let l = random_closed_ladder(&mut rng, n, bottom, steps);
            if l.rungs.len() > 6 || SBim::<F31>::identity(l.n, &l.bottom).add_rungs(&l.rungs).fsize() > 64 {
                continue;
            }
            assert_eq!(state_space_poly::<F31>(&l), closure_value(&l), "{:?}", l);
        }
    }
}

#[test]
fn state_space_over_q() {
    let l = lad(3, &[1, 2, 0], &[('E', 0, 1), ('F', 0, 1), ('F', 1, 1), ('E', 0, 1), ('F', 0, 1), ('E', 1, 1)]);
    assert_eq!(state_space_poly::<Q>(&l), closure_value(&l));
}

#[test]
fn rickard_differentials_unique() {
    for n in 2..=4usize {
        for c in 1..=n {
            for d in 1..=n {
                for pos in [true, false] {
                    let cx = crossing_cx(Gen { i: 0, positive: pos }, c, d, n);
                    let id = SBim::<F31>::identity(n, &[c, d]);
                    let webs: Vec<_> = cx.terms.iter().map(|t| id.add_rungs(&t.rungs)).collect();
                    let mut mc = MonoCache::new(id.ring.clone());
                    let mut maps = vec![];
                    for j in 0..webs.len().saturating_sub(1) {
                        let (s, t) = if pos { (&webs[j], &webs[j + 1]) } else { (&webs[j + 1], &webs[j]) };
                        let h = hom_space(s, t, 1 + s.shift - t.shift, &mut mc);
                        assert_eq!(h.len(), 1, "N={} c={} d={} pos={} j={}", n, c, d, pos, j);
                        maps.push(h.into_iter().next().unwrap());
                    }
                    for j in 0..maps.len().saturating_sub(1) {
                        let comp = if pos { maps[j + 1].matmul(&maps[j]) } else { maps[j].matmul(&maps[j + 1]) };
                        assert!(comp.is_zero());
                    }
                }
            }
        }
    }
}

fn gens(w: &[i32]) -> Vec<Gen> {
    w.iter().map(|&k| Gen { i: (k.abs() - 1) as usize, positive: k > 0 }).collect()
}

fn poincare(n: usize, labels: &[usize], w: &[i32]) -> Poincare {
    khr_full_cube::<F31>(&BraidLink { n, labels: labels.to_vec(), word: gens(w) }, QConv::A).poincare
}

#[test]
fn khovanov_trefoil_and_figure_eight() {
    let mut t = Poincare::default();
    for &(h, q) in &[(0, -1), (0, -3), (2, -5), (3, -9)] {
        t.add_term(h, q, 1);
    }
    assert_eq!(poincare(2, &[1, 1], &[1, 1, 1]), t);
    let mut f = Poincare::default();
    for &(h, q) in &[(-2, 5), (-1, 1), (0, -1), (0, 1), (1, -1), (2, -5)] {
        f.add_term(h, q, 1);
    }
    assert_eq!(poincare(2, &[1, 1, 1], &[1, -2, 1, -2]), f);
}

#[test]
fn trefoil_kr_formula_and_markov_invariance() {
    for n in 2..=4usize {
        // q^{2-2N}[N] + t^2 q^{-2N-1}[N-1] + t^3 q^{-4N-1}[N-1]
        let mut expect = Poincare::default();
        for (e, k) in LPoly::qint(n as i64).shift(2 - 2 * n as i32).c {
            expect.add_term(0, e, k);
        }
        for (e, k) in LPoly::qint(n as i64 - 1).shift(-2 * n as i32 - 1).c {
            expect.add_term(2, e, k);
        }
        for (e, k) in LPoly::qint(n as i64 - 1).shift(-4 * n as i32 - 1).c {
            expect.add_term(3, e, k);
        }
        assert_eq!(poincare(n, &[1, 1], &[1, 1, 1]), expect, "N={}", n);
        assert_eq!(poincare(n, &[1, 1, 1], &[1, 1, 1, 2]), expect, "N={} stabilized", n);
        assert_eq!(poincare(n, &[1, 1, 1], &[1, 1, 1, -2]), expect, "N={} neg stabilized", n);
    }
}

#[test]
fn colored_unknot_stabilization() {
    for n in 2..=4usize {
        for a in 1..n {
            let mut expect = Poincare::default();
            for (e, k) in LPoly::qbinom(n as i64, a as i64).c {
                expect.add_term(0, e, k);
            }
            assert_eq!(poincare(n, &[a, a], &[1]), expect);
            assert_eq!(poincare(n, &[a, a], &[-1]), expect);
        }
    }
}

#[test]
fn soergel_engine_matches_full_cube() {
    use krhom::engine2::khr_engine2;
    let cases: Vec<(usize, Vec<i32>)> = vec![(2, vec![1, 1, 1]), (3, vec![1, -2, 1, -2]), (3, vec![1, 1, 1, 2, -1, 2]), (4, vec![1, 1, 2, -1, -3, 2, -3])];
    for (n, w) in cases {
        let word = gens(&w);
        let r = khr_engine2::<F31>(n, &word, &[2, 3], &[], 5);
        for nn in [2usize, 3] {
            let cube = khr_full_cube::<F31>(&BraidLink { n: nn, labels: vec![1; n], word: word.clone() }, QConv::A).poincare;
            assert_eq!(r.poincare[&nn], cube, "braid {:?} N={}", w, nn);
        }
    }
}
