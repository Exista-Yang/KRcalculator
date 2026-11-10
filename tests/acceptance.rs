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
    use krhom::engine2::{close_with, khr_engine2, prepare_closure_data, Closure};
    let cases: Vec<(usize, Vec<i32>)> = vec![(2, vec![1, 1, 1]), (3, vec![1, -2, 1, -2]), (3, vec![1, 1, 1, 2, -1, 2]), (4, vec![1, 1, 2, -1, -3, 2, -3])];
    for (n, w) in cases {
        let word = gens(&w);
        let r = khr_engine2::<F31>(n, &word, &[2, 3], &[], 5);
        let data = prepare_closure_data::<F31>(n, &word, &[], 5);
        for nn in [2usize, 3] {
            let cube = khr_full_cube::<F31>(&BraidLink { n: nn, labels: vec![1; n], word: word.clone() }, QConv::A).poincare;
            assert_eq!(r.poincare[&nn], cube, "braid {:?} N={}", w, nn);
            // the restriction tower against the same independent full-cube computation
            for end in TOWER_ENDS {
                assert_eq!(close_with(&data, nn, false, Closure::Tower(end)), cube, "tower ({:?}) braid {:?} N={}", end, w, nn);
            }
        }
    }
}

// ---------------------------------------------------------------- M1b.1 restriction tower

use krhom::engine2::{close_tower_with, close_with, exact_type_dims, prepare_closure_data, tower_type_dims, Closure, ClosureData};
use krhom::tower::{EndPolicy, TowerOpts};

const TOWER_ENDS: [EndPolicy; 3] = [EndPolicy::Right, EndPolicy::Left, EndPolicy::Alternate];
const CONWAY: [i32; 11] = [1, 1, 2, -3, 2, 1, -3, -2, -2, -3, -3];

fn data_for(n: usize, w: &[i32], seed: u64) -> ClosureData<F31> {
    prepare_closure_data::<F31>(n, &gens(w), &[], seed)
}

/// [N]·P_reduced − P_unreduced is divisible by (1+t) with a nonnegative quotient
/// (Lewark-type spectral sequence check); returns the quotient.
fn lewark_quotient(nn: usize, red: &Poincare, unred: &Poincare) -> Result<Poincare, String> {
    use std::collections::BTreeMap;
    let mut d: BTreeMap<i32, BTreeMap<i32, i128>> = BTreeMap::new();
    for (&(h, q), &c) in &red.c {
        for k in 0..nn as i32 {
            *d.entry(q + nn as i32 - 1 - 2 * k).or_default().entry(h).or_default() += c;
        }
    }
    for (&(h, q), &c) in &unred.c {
        *d.entry(q).or_default().entry(h).or_default() -= c;
    }
    let mut quot = Poincare::default();
    for (q, row) in d {
        let (Some(&h0), Some(&h1)) = (row.keys().next(), row.keys().last()) else { continue };
        let mut prev = 0i128;
        for h in h0..=h1 {
            let e = row.get(&h).copied().unwrap_or(0) - prev;
            if e < 0 {
                return Err(format!("negative quotient coefficient at t^{} q^{}", h, q));
            }
            if h == h1 && e != 0 {
                return Err(format!("not divisible by (1+t) at q^{}", q));
            }
            quot.add_term(h, q, e);
            prev = e;
        }
    }
    Ok(quot)
}

#[test]
fn tower_matches_exact_on_small_braids_and_links() {
    let cases: Vec<(usize, Vec<i32>, &str)> = vec![
        (2, vec![1, 1, 1], "trefoil"),
        (3, vec![1, -2, 1, -2], "figure-eight"),
        (3, vec![1, 1, 1, 2, -1, 2], "3-strand knot"),
        (4, vec![1, 1, 2, -1, -3, 2, -3], "4-strand braid"),
        (2, vec![1, 1], "Hopf link"),
        (2, vec![-1, -1, -1, -1], "T(2,-4) link"),
        (3, vec![1, 1, -2, -2], "3-component chain"),
        (2, vec![1, -1], "2-component unlink"),
    ];
    for (n, w, name) in cases {
        let data = data_for(n, &w, 5);
        for nn in 2..=6usize {
            for reduced in [false, true] {
                let exact = close_with(&data, nn, reduced, Closure::Exact);
                for end in TOWER_ENDS {
                    let tower = close_with(&data, nn, reduced, Closure::Tower(end));
                    assert_eq!(tower, exact, "{} {:?} N={} reduced={} end={:?}", name, w, nn, reduced, end);
                }
            }
        }
    }
}

#[test]
fn tower_reduced_trefoil_formula() {
    // reduced trefoil: q^{-2(N-1)} + t^2 q^{-2(N+1)} + t^3 q^{-4N}
    let data = data_for(2, &[1, 1, 1], 5);
    for nn in 2..=10usize {
        let n = nn as i32;
        let mut expect = Poincare::default();
        expect.add_term(0, -2 * (n - 1), 1);
        expect.add_term(2, -2 * (n + 1), 1);
        expect.add_term(3, -4 * n, 1);
        for end in TOWER_ENDS {
            assert_eq!(close_with(&data, nn, true, Closure::Tower(end)), expect, "N={} end={:?}", nn, end);
        }
        if nn <= 6 {
            assert_eq!(close_with(&data, nn, true, Closure::Exact), expect, "exact N={}", nn);
        }
    }
}

#[test]
fn tower_conway_unreduced_matches_exact() {
    let data = data_for(4, &CONWAY, 7);
    let ranks = [34i128, 83, 132, 181, 230, 279, 328, 377, 426];
    for nn in 2..=10usize {
        let right = close_with(&data, nn, false, Closure::Tower(EndPolicy::Right));
        assert_eq!(right.total_rank(), ranks[nn - 2], "Conway rank N={}", nn);
        for end in [EndPolicy::Left, EndPolicy::Alternate] {
            assert_eq!(close_with(&data, nn, false, Closure::Tower(end)), right, "Conway N={} end={:?}", nn, end);
        }
        if nn <= 8 {
            assert_eq!(close_with(&data, nn, false, Closure::Exact), right, "Conway N={} tower vs exact", nn);
        }
    }
}

#[test]
fn tower_conway_reduced_and_lewark() {
    let data = data_for(4, &CONWAY, 7);
    for nn in 2..=10usize {
        let red = close_with(&data, nn, true, Closure::Tower(EndPolicy::Right));
        assert_eq!(red.total_rank(), if nn == 2 { 33 } else { 49 }, "Conway reduced rank N={}", nn);
        assert_eq!(close_with(&data, nn, true, Closure::Tower(EndPolicy::Left)), red, "Conway reduced N={} left end", nn);
        assert_eq!(close_with(&data, nn, true, Closure::Exact), red, "Conway reduced N={} tower vs exact", nn);
        if nn <= 6 {
            let unred = close_with(&data, nn, false, Closure::Tower(EndPolicy::Right));
            lewark_quotient(nn, &red, &unred).unwrap_or_else(|e| panic!("Conway N={}: {}", nn, e));
        }
    }
    // the same check on small knots
    for (n, w) in [(2usize, vec![1, 1, 1]), (3, vec![1, -2, 1, -2]), (3, vec![1, 1, 1, 2, -1, 2])] {
        let data = data_for(n, &w, 5);
        for nn in 2..=6usize {
            let red = close_with(&data, nn, true, Closure::Tower(EndPolicy::Right));
            let unred = close_with(&data, nn, false, Closure::Tower(EndPolicy::Right));
            lewark_quotient(nn, &red, &unred).unwrap_or_else(|e| panic!("{:?} N={}: {}", w, nn, e));
        }
    }
}

#[test]
fn tower_chain_dims() {
    // dim Φ_N(B_x) per type: tower decomposition = close_exact state spaces, and
    // dim_unreduced = N · dim_reduced (freeness over the marked strand)
    let data = data_for(4, &CONWAY, 7);
    for nn in 2..=5usize {
        for reduced in [false, true] {
            let exact = exact_type_dims(&data, nn, reduced);
            for end in TOWER_ENDS {
                for unit_pivots_only in [false, true] {
                    let opts = TowerOpts { policy: end, unit_pivots_only, eliminate_top: false };
                    assert_eq!(tower_type_dims(&data, nn, reduced, opts), exact, "N={} reduced={} {:?}", nn, reduced, opts);
                }
            }
        }
        let (u, r) = (exact_type_dims(&data, nn, false), exact_type_dims(&data, nn, true));
        for (x, du) in &u {
            assert_eq!(*du, nn * r[x], "type {} N={}: unreduced {} reduced {}", x, nn, du, r[x]);
        }
    }
}

#[test]
fn tower_correct_without_splitting() {
    // Stress mode: only constant pivots, so types that split by exact division stay
    // non-split blocks (non-cyclic objects, gdim-certified and Neumann-inverted
    // eliminations). Correctness must not depend on the splitting.
    let mut seen = [0usize; 4];
    let cases: Vec<(usize, Vec<i32>, Vec<usize>)> = vec![
        (2, vec![1, 1, 1], vec![2, 3, 4, 5]),
        (3, vec![1, -2, 1, -2], vec![2, 3, 4, 5]),
        (3, vec![1, 1, -2, -2], vec![2, 3, 4]),
        (4, vec![1, 1, 2, -1, -3, 2, -3], vec![2, 3, 4]),
        (4, CONWAY.to_vec(), vec![2, 3, 4]),
    ];
    for (n, w, nns) in cases {
        let data = data_for(n, &w, 5);
        for nn in nns {
            for reduced in [false, true] {
                let exact = close_with(&data, nn, reduced, Closure::Exact);
                for end in [EndPolicy::Right, EndPolicy::Left] {
                    let r = close_tower_with(&data, nn, reduced, TowerOpts { policy: end, unit_pivots_only: true, eliminate_top: false });
                    assert_eq!(r.poincare, exact, "{:?} N={} reduced={} end={:?} (unit pivots only)", w, nn, reduced, end);
                    for st in &r.stages {
                        seen[0] += st.nonsplit;
                        seen[1] += st.elims_noncyclic;
                        seen[2] += st.elims_neumann;
                        seen[3] += st.elims_gdim;
                    }
                }
            }
        }
    }
    // the non-split code path was actually exercised (real minimal complexes offer no
    // isomorphisms between non-cyclic objects; those are covered by the synthetic test)
    assert!(seen[0] > 0 && seen[3] > 0, "non-split / non-cyclic / Neumann / gdim counts {:?}", seen);
}

// Synthetic non-minimal complexes (reverse Gaussian elimination). Objects are
// Φ_N-images Im(G⊗A) of Gram-type matrices: X, Y of type A (Gram G), Z of type A'
// (Gram G·P1^{-1}) and Z' of type A'' (Gram G·P2^{-1}), with P1, P2 unipotent and
// non-constant. The complex
//     X ⊕ Z  --[[ (s_h + s_α s_β)·I , s_α P1^{-1} ], [ s_β P2 , P2 P1^{-1} ]]-->  Y ⊕ Z'
// is homotopy equivalent (cancel φ = P2 P1^{-1}: Z → Z') to X --s_h·I--> Y, so the answer
// depends on the elimination update d_YX − d_YZ φ^{-1} d_Z'X being exactly right; φ is a
// non-cyclic isomorphism between different generator matrices with a non-constant inverse.
mod synthetic {
    use krhom::engine::Obj;
    use krhom::engine2::{ClosureData, EPoly, TypeBasis};
    use krhom::field::{Field, F31};
    use std::collections::BTreeMap;

    type EP = EPoly<F31>;
    type M = Vec<Vec<EP>>;

    fn mono(n: usize, e: &[(usize, u32)], c: i64) -> EP {
        let mut v = vec![0u32; n];
        for &(i, k) in e {
            v[i] += k;
        }
        let mut p = EP::new();
        p.insert(v, F31::from_i64(c));
        p
    }
    fn add(a: &EP, b: &EP) -> EP {
        let mut r = a.clone();
        for (m, c) in b {
            r.entry(m.clone()).or_insert_with(F31::zero).add_assign(c);
        }
        r.retain(|_, c| !c.is_zero());
        r
    }
    fn mul(a: &EP, b: &EP) -> EP {
        let mut r = EP::new();
        for (m1, c1) in a {
            for (m2, c2) in b {
                let m: Vec<u32> = m1.iter().zip(m2).map(|(x, y)| x + y).collect();
                r.entry(m).or_insert_with(F31::zero).fma(c1, c2);
            }
        }
        r.retain(|_, c| !c.is_zero());
        r
    }
    fn sum(ps: &[EP]) -> EP {
        ps.iter().fold(EP::new(), |a, b| add(&a, b))
    }
    fn mmul(a: &M, b: &M) -> M {
        (0..a.len()).map(|i| (0..b[0].len()).map(|j| sum(&(0..b.len()).map(|k| mul(&a[i][k], &b[k][j])).collect::<Vec<_>>())).collect()).collect()
    }
    fn scal(s: &EP, a: &M) -> M {
        a.iter().map(|r| r.iter().map(|x| mul(s, x)).collect()).collect()
    }
    fn unipotent(n: usize, l: EP) -> (M, M) {
        let one = mono(n, &[], 1);
        let neg = mul(&mono(n, &[], -1), &l);
        (vec![vec![one.clone(), l], vec![EP::new(), one.clone()]], vec![vec![one.clone(), neg], vec![EP::new(), one]])
    }

    /// (big, small) closure data on n ≥ 2 strands (variables), see above.
    pub fn cone_pair(n: usize) -> (ClosureData<F31>, ClosureData<F31>) {
        let v = |i: usize, c: i64| mono(n, &[(i % n, 1)], c);
        let z = n - 1; // the last variable
        // G: coordinates of degrees (0, 2), relations of degrees (2, 2)
        let g: M = vec![
            vec![sum(&[v(0, 1), v(1, -1), v(z, 2)]), sum(&[mono(n, &[(0, 1), (1, 1)], 1), mono(n, &[(1, 2)], 1), mono(n, &[(z, 2)], -1)])],
            vec![sum(&[v(0, 1), v(1, 2)]), sum(&[mono(n, &[(0, 2)], 1), mono(n, &[(1, 1), (z, 1)], -3), mono(n, &[(1, 2)], 2)])],
        ];
        let (_p1, p1i) = unipotent(n, sum(&[v(0, 1), v(1, 2), v(z, -1)]));
        let (p2, p2i) = unipotent(n, sum(&[v(0, 3), v(1, -1)]));
        let s_h = sum(&[mono(n, &[(0, 2)], 1), mono(n, &[(1, 1), (z, 1)], 1)]);
        let s_a = v(0, 1);
        let s_b = sum(&[v(1, 1), v(0, -1), v(z, 1)]);
        let id: M = vec![vec![mono(n, &[], 1), EP::new()], vec![EP::new(), mono(n, &[], 1)]];
        let tb = |gram: M| TypeBasis { ins: vec![(0, 0), (1, 2)], outs: vec![(0, 2), (1, 2)], gram };
        let mut bases = BTreeMap::new();
        bases.insert(0usize, tb(g.clone()));
        bases.insert(1usize, tb(mmul(&g, &p1i)));
        bases.insert(2usize, tb(mmul(&g, &p2i)));
        let objs = vec![Obj { x: 0, q: 0, h: 0 }, Obj { x: 1, q: -2, h: 0 }, Obj { x: 0, q: -4, h: 1 }, Obj { x: 2, q: -2, h: 1 }];
        let big_edges = vec![(0usize, 2usize), (1, 2), (0, 3), (1, 3)];
        let big_f = vec![scal(&add(&s_h, &mul(&s_a, &s_b)), &id), scal(&s_a, &p1i), scal(&s_b, &p2), mmul(&p2, &p1i)];
        let mk = |objs: Vec<Obj>, edges: Vec<(usize, usize)>, fpolys: Vec<M>, bases: BTreeMap<usize, TypeBasis<F31>>| ClosureData {
            n,
            word: vec![],
            objs,
            edges,
            fpolys,
            bases,
            points: 0,
            max_objects: 0,
            timings: (0.0, 0.0, 0.0),
        };
        let mut small_bases = BTreeMap::new();
        small_bases.insert(0usize, tb(g.clone()));
        let small = mk(vec![Obj { x: 0, q: 0, h: 0 }, Obj { x: 0, q: -4, h: 1 }], vec![(0, 1)], vec![scal(&s_h, &id)], small_bases);
        (mk(objs, big_edges, big_f, bases), small)
    }
}

#[test]
fn tower_synthetic_noncyclic_elimination() {
    let mut seen = [0usize; 3];
    for n in [2usize, 3] {
        let (big, small) = synthetic::cone_pair(n);
        for nn in 2..=5usize {
            for reduced in [false, true] {
                let exact = close_with(&big, nn, reduced, Closure::Exact);
                assert_eq!(close_with(&small, nn, reduced, Closure::Exact), exact, "construction n={} N={}", n, nn);
                assert!(exact.total_rank() > 0);
                for policy in [EndPolicy::Right, EndPolicy::Left] {
                    for eliminate_top in [false, true] {
                        for unit_pivots_only in [false, true] {
                            let opts = TowerOpts { policy, unit_pivots_only, eliminate_top };
                            let r = close_tower_with(&big, nn, reduced, opts);
                            assert_eq!(r.poincare, exact, "n={} N={} reduced={} {:?}", n, nn, reduced, opts);
                            for st in &r.stages {
                                seen[0] += st.elims_noncyclic;
                                seen[1] += st.elims_neumann;
                                seen[2] += st.elims_gdim;
                            }
                        }
                    }
                }
            }
        }
    }
    assert!(seen.iter().all(|&c| c > 0), "non-cyclic / Neumann / gdim eliminations {:?}", seen);
}

#[test]
fn tower_second_prime() {
    // the same comparison over F_p, p = 4294967291 (two-prime agreement for Q-ranks)
    use krhom::field::F32b;
    let data = prepare_closure_data::<F32b>(4, &gens(&CONWAY), &[], 11);
    let ranks = [34i128, 83, 132, 181, 230];
    for nn in 2..=6usize {
        for reduced in [false, true] {
            let tower = close_with(&data, nn, reduced, Closure::Tower(EndPolicy::Right));
            if nn <= 5 {
                assert_eq!(tower, close_with(&data, nn, reduced, Closure::Exact), "F32b Conway N={} reduced={}", nn, reduced);
            }
            assert_eq!(tower.total_rank(), if reduced { if nn == 2 { 33 } else { 49 } } else { ranks.get(nn - 2).copied().unwrap_or(230) }, "N={}", nn);
        }
    }
}

#[test]
fn tower_six_strand_nonsplit_right_end() {
    // 6-strand braid whose right end does not split under 1x1 pivots (design critique):
    // the right-end tower keeps non-split blocks and must still agree with exact closure
    // and with the (fully splitting) left end.
    let data = data_for(6, &[-2, -3, -1, -5, 2, 4, -5, -2, 1, -2, 1, -2, 5], 7);
    let mut nonsplit = 0;
    for nn in 2..=5usize {
        let right = close_tower_with(&data, nn, false, TowerOpts::new(EndPolicy::Right));
        nonsplit += right.stages.iter().map(|s| s.nonsplit).sum::<usize>();
        assert_eq!(right.poincare.total_rank(), 3 * nn as i128 - 2, "N={}", nn);
        assert_eq!(close_with(&data, nn, false, Closure::Tower(EndPolicy::Left)), right.poincare, "left vs right N={}", nn);
        if nn <= 3 {
            assert_eq!(close_with(&data, nn, false, Closure::Exact), right.poincare, "exact vs tower N={}", nn);
        }
    }
    assert!(nonsplit > 0, "expected non-split blocks at the right end");
}

/// Closure data with one type of Gram matrix [[1]] (Φ_N(B) = A_n) and one map
/// X{0}[0] -> Y{-2k}[1] given by multiplication with the monomial x^e (|e| = k).
fn monomial_map_data(n: usize, e: &[u32]) -> ClosureData<F31> {
    use krhom::engine::Obj;
    use krhom::engine2::{EPoly, TypeBasis};
    use krhom::field::Field;
    use std::collections::BTreeMap;
    let k: u32 = e.iter().sum();
    let mono = |e: &[u32]| -> EPoly<F31> { [(e.to_vec(), F31::one())].into_iter().collect() };
    let mut bases = BTreeMap::new();
    bases.insert(0usize, TypeBasis { ins: vec![(0, 0)], outs: vec![(0, 0)], gram: vec![vec![mono(&vec![0; n])]] });
    ClosureData {
        n,
        word: vec![],
        objs: vec![Obj { x: 0, q: 0, h: 0 }, Obj { x: 0, q: -2 * k as i32, h: 1 }],
        edges: vec![(0, 1)],
        fpolys: vec![vec![vec![mono(e)]]],
        bases,
        points: 0,
        max_objects: 0,
        timings: (0.0, 0.0, 0.0),
    }
}

#[test]
fn exact_closure_large_n_no_packed_overflow() {
    // Regression: close_exact packed monomials with 5 bits per variable, so for N ≥ 17 the
    // sum of two exponents ≤ N-1 carried into the next field (bogus columns, index panic
    // in RankEchelon or silently wrong ranks). Multiplication by x^e on A_n has rank
    // Π_i max(N - e_i, 0) (unreduced) and the closure must have total rank 2(N^n - rank).
    for (n, e) in [(2usize, vec![16u32, 0]), (2, vec![0, 16]), (2, vec![9, 9]), (3, vec![0, 17, 0]), (3, vec![5, 0, 12]), (3, vec![1, 1, 1])] {
        let data = monomial_map_data(n, &e);
        for nn in [16usize, 17, 18, 20, 32] {
            for reduced in [false, true] {
                let caps = krhom::engine2::closure_caps(n, nn, reduced);
                let dim: i128 = caps.iter().map(|&c| c as i128 + 1).product();
                let rank: i128 = caps.iter().zip(&e).map(|(&c, &x)| (c as i128 + 1 - x as i128).max(0)).product();
                let exact = close_with(&data, nn, reduced, Closure::Exact);
                assert_eq!(exact.total_rank(), 2 * (dim - rank), "n={} x^{:?} N={} reduced={}", n, e, nn, reduced);
                assert_eq!(close_with(&data, nn, reduced, Closure::Tower(EndPolicy::Right)), exact, "tower n={} x^{:?} N={} reduced={}", n, e, nn, reduced);
            }
        }
    }
    // beyond the tower's range (N ≤ 32) the oracle still works
    let data = monomial_map_data(2, &[40, 3]);
    let exact = close_with(&data, 64, false, Closure::Exact);
    assert_eq!(exact.total_rank(), 2 * (64 * 64 - 24 * 61));
}

#[test]
fn exact_closure_large_n_braids() {
    // Regression for the same overflow on real closure data: the figure-eight knot and
    // "1 1 1 2" on 3 strands, where close_exact panicked at N = 32 (tower ranks 5N-4, 3N-2).
    for (w, a, b) in [(vec![1, -2, 1, -2], 5i128, -4i128), (vec![1, 1, 1, 2], 3, -2)] {
        let data = data_for(3, &w, 5);
        for nn in [17usize, 20, 32] {
            let exact = close_with(&data, nn, false, Closure::Exact);
            assert_eq!(exact.total_rank(), a * nn as i128 + b, "{:?} N={}", w, nn);
            assert_eq!(close_with(&data, nn, false, Closure::Tower(EndPolicy::Right)), exact, "{:?} N={}", w, nn);
            let red_exact = close_with(&data, nn, true, Closure::Exact);
            assert_eq!(close_with(&data, nn, true, Closure::Tower(EndPolicy::Left)), red_exact, "{:?} N={} reduced", w, nn);
        }
    }
}

#[test]
fn closure_data_save_load_roundtrip() {
    use krhom::engine2::{load_closure_data, save_closure_data};
    let data = data_for(3, &[1, -2, 1, -2], 5);
    let path = std::env::temp_dir().join(format!("krhom_roundtrip_{}.data", std::process::id()));
    let path = path.to_str().unwrap();
    save_closure_data(&data, path).unwrap();
    let back = load_closure_data::<F31>(path).unwrap();
    let wrong_field = load_closure_data::<krhom::field::F32b>(path);
    std::fs::remove_file(path).ok();
    assert!(wrong_field.is_err_and(|e| e.contains("field mismatch")));
    assert_eq!((back.n, &back.word, &back.objs, &back.edges, &back.fpolys), (data.n, &data.word, &data.objs, &data.edges, &data.fpolys));
    assert_eq!(back.bases.len(), data.bases.len());
    for (x, tb) in &data.bases {
        let tb2 = &back.bases[x];
        assert_eq!((&tb2.ins, &tb2.outs, &tb2.gram), (&tb.ins, &tb.outs, &tb.gram));
    }
    for nn in 2..=4usize {
        assert_eq!(close_with(&back, nn, false, Closure::Tower(EndPolicy::Right)), close_with(&data, nn, false, Closure::Exact));
    }
}

#[test]
fn exact_closure_degree_range_when_outs_below_ins() {
    // Regression: exact_state_space ran its degree loop up to topdeg - min(min ins, max outs),
    // which drops the top degrees when min(outs) < min(ins) (found by the re-verification).
    // One type with ins of degrees 0, 2 and outs of degrees 0, -2 and Gram = identity has
    // Φ_N = A_1 ⊕ A_1{2}, so the closure has total rank 2N (unreduced) and 2 (reduced).
    use krhom::engine::Obj;
    use krhom::engine2::{EPoly, TypeBasis};
    use krhom::field::Field;
    use std::collections::BTreeMap;
    let one: EPoly<F31> = [(vec![0u32], F31::one())].into_iter().collect();
    let zero: EPoly<F31> = EPoly::new();
    let mut bases = BTreeMap::new();
    bases.insert(0usize, TypeBasis { ins: vec![(0, 0), (1, 2)], outs: vec![(0, 0), (1, -2)], gram: vec![vec![one.clone(), zero.clone()], vec![zero, one]] });
    let data = ClosureData { n: 1, word: vec![], objs: vec![Obj { x: 0, q: 0, h: 0 }], edges: vec![], fpolys: vec![], bases, points: 0, max_objects: 0, timings: (0.0, 0.0, 0.0) };
    for nn in 2..=6usize {
        for reduced in [false, true] {
            let exact = close_with(&data, nn, reduced, Closure::Exact);
            let want = if reduced { 2 } else { 2 * nn as i128 };
            assert_eq!(exact.total_rank(), want, "exact N={} reduced={}", nn, reduced);
            for end in [EndPolicy::Right, EndPolicy::Left] {
                assert_eq!(close_with(&data, nn, reduced, Closure::Tower(end)), exact, "tower N={} reduced={}", nn, reduced);
            }
        }
    }
}

#[test]
fn tower_adds_parallel_edges() {
    // Regression: the tower kept only the last of several edges with the same source and
    // target, while close_exact adds them. X{0}[0] -> Y{-2}[1] with F = x0 and F = ±x0:
    // the sum 0 gives total rank 2N^2, the sum 2x0 behaves like x0 (rank 2N).
    use krhom::field::Field;
    for (sign, want) in [(-1i64, None), (1, Some(()))] {
        let mut data = monomial_map_data(2, &[1, 0]);
        let mut second = data.fpolys[0].clone();
        if sign < 0 {
            for row in second.iter_mut() {
                for p in row.iter_mut() {
                    for c in p.values_mut() {
                        *c = c.neg();
                    }
                }
            }
        }
        data.edges.push((0, 1));
        data.fpolys.push(second);
        for nn in 2..=5usize {
            let n2 = (nn * nn) as i128;
            let expected = if want.is_none() { 2 * n2 } else { 2 * nn as i128 };
            let exact = close_with(&data, nn, false, Closure::Exact);
            assert_eq!(exact.total_rank(), expected, "exact sign={} N={}", sign, nn);
            for end in [EndPolicy::Right, EndPolicy::Left, EndPolicy::Alternate] {
                assert_eq!(close_with(&data, nn, false, Closure::Tower(end)), exact, "tower sign={} N={}", sign, nn);
            }
        }
    }
}

// ---------------------------------------------------------------- colored homology by cabling

use krhom::colored::{colored_euler_moy, khr_colored_cable};

/// Colored homology by cabling (uncolored Soergel engine + tower closure), all N at once;
/// also checks the Euler characteristic against the MOY state sum.
fn cable_colored(labels: &[usize], w: &[i32], ns: &[usize], closure: Closure) -> std::collections::BTreeMap<usize, Poincare> {
    let word = gens(w);
    let r = khr_colored_cable::<F31>(labels, &word, ns, 7, closure).unwrap_or_else(|e| panic!("cable {:?} {:?}: {}", labels, w, e));
    for (&nn, p) in &r.colored {
        let chi = colored_euler_moy(labels, &word, nn, 100_000).unwrap().expect("too many resolutions");
        assert_eq!(p.euler(), chi, "Euler characteristic vs MOY: labels {:?} braid {:?} N={}", labels, w, nn);
    }
    r.colored
}

#[test]
fn cable_matches_cube_tiny_colored() {
    // (bottom labels, braid, N values): crossings (2,2) with m = 2; (2,1), (1,2) with m = 1, T = 1;
    // (3,1) with m = 1, T = 2; both signs; links with 2 and 3 components.
    let cases: Vec<(Vec<usize>, Vec<i32>, Vec<usize>)> = vec![
        (vec![2, 2], vec![1], vec![2, 3, 4]),
        (vec![2, 2], vec![-1], vec![2, 3, 4]),
        (vec![2, 1], vec![1, 1], vec![2, 3, 4]),
        (vec![2, 1], vec![-1, -1], vec![2, 3, 4]),
        (vec![1, 2], vec![1, 1], vec![2, 3, 4]),
        (vec![1, 2], vec![-1, -1], vec![2, 3, 4]),
        (vec![2, 1, 1], vec![1, 1, -2, -2], vec![2, 3, 4]),
        (vec![1, 2, 1], vec![-1, -1, 2, 2], vec![2, 3, 4]),
        (vec![3, 1], vec![1, 1], vec![2, 3, 4]),
        (vec![3, 1], vec![-1, -1], vec![3, 4]),
        (vec![2, 2], vec![1, 1], vec![3, 4]),
    ];
    for (labels, w, ns) in cases {
        let cab = cable_colored(&labels, &w, &ns, Closure::Tower(EndPolicy::Right));
        for &nn in &ns {
            let cube = poincare(nn, &labels, &w);
            assert_eq!(cab[&nn], cube, "cable vs cube: labels {:?} braid {:?} N={}", labels, w, nn);
        }
    }
    // the exact closure (oracle) on cabled data with a start object
    let ex = cable_colored(&[2, 1], &[1, 1], &[2, 3, 4], Closure::Exact);
    for nn in 2..=4usize {
        assert_eq!(ex[&nn], poincare(nn, &[2, 1], &[1, 1]), "exact closure N={}", nn);
    }
}

#[test]
fn cable_colored_unknot_trefoil_and_uncolored() {
    let ns: Vec<usize> = (2..=6).collect();
    // 2-colored unknot (one (2,2) crossing of either sign) = [N choose 2] in degree 0
    for w in [[1], [-1]] {
        let r = cable_colored(&[2, 2], &w, &ns, Closure::Tower(EndPolicy::Right));
        for &nn in &ns {
            let mut expect = Poincare::default();
            for (e, k) in LPoly::qbinom(nn as i64, 2).c {
                expect.add_term(0, e, k);
            }
            assert_eq!(r[&nn], expect, "2-colored unknot braid {:?} N={}", w, nn);
        }
    }
    // 2-colored trefoil, gl_3 (the full cube of the colored Rickard complexes, 42 min, gave this)
    // and gl_2 (∧^2 C^2 is the determinant: one-dimensional); the positive trefoil is the mirror.
    let mut neg3 = Poincare::default();
    for &(h, q) in &[(-3, 12), (-3, 14), (-2, 6), (-2, 8), (0, 2), (0, 4), (0, 6)] {
        neg3.add_term(h, q, 1);
    }
    let mut one = Poincare::default();
    one.add_term(0, 0, 1);
    let neg = cable_colored(&[2, 2], &[-1, -1, -1], &[2, 3], Closure::Tower(EndPolicy::Right));
    assert_eq!(neg[&3], neg3);
    assert_eq!(neg[&2], one);
    let pos = cable_colored(&[2, 2], &[1, 1, 1], &[3], Closure::Tower(EndPolicy::Right));
    let mut mirror = Poincare::default();
    for (&(h, q), &k) in &neg3.c {
        mirror.add_term(-h, -q, k);
    }
    assert_eq!(pos[&3], mirror);
    // all colors 1: the cable is the braid itself and the shifts vanish (trefoil KR formula)
    let un = cable_colored(&[1, 1], &[1, 1, 1], &[2, 3, 4], Closure::Tower(EndPolicy::Right));
    for nn in 2..=4usize {
        assert_eq!(un[&nn], poincare(nn, &[1, 1], &[1, 1, 1]), "uncolored N={}", nn);
    }
}

#[test]
fn cable_two_colored_t25_and_three_colored_t24_euler() {
    // larger cables, checked against the MOY Euler characteristic (inside cable_colored) and the
    // homological range of the colored complex (inside khr_colored_cable); colors > N give 0.
    let t25 = cable_colored(&[2, 2], &[1, 1, 1, 1, 1], &[2, 3, 4, 5, 6], Closure::Tower(EndPolicy::Right));
    let ranks: Vec<i128> = (2..=6).map(|nn| t25[&nn].total_rank()).collect();
    assert_eq!(ranks, vec![1, 11, 46, 106, 191]);
    // ∧^2 C^3 is the dual of C^3 and T(2,5) is invertible: 2-colored gl_3 = uncolored gl_3 (full cube)
    assert_eq!(t25[&3], poincare(3, &[1, 1], &[1, 1, 1, 1, 1]));
    let t24 = cable_colored(&[3, 1], &[1, 1, 1, 1], &[2, 3, 4, 5], Closure::Tower(EndPolicy::Right));
    assert!(t24[&2].c.is_empty());
}
