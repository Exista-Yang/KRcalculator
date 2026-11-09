//! Tests for PD codes -> braid words (src/pd.rs).
//!
//! The checks use an independent implementation (in this file) of
//!   * orientation / crossing signs of a PD code (walk along components), and
//!   * the Kauffman bracket by planar contraction, giving the Jones polynomial
//!     V = (-A^3)^{-w} <D> / δ at A = q^{-1/4}  (V(unknot) = 1, V(right trefoil) = q + q^3 - q^4),
//! and compare PD codes with the PD codes of the closures of the computed braids.

use krhom::crossing::Gen;
use krhom::pd::{braid_to_pd, parse_pd, pd_to_braid, pd_to_braid_vogel, BraidWord, Pd};
use std::collections::{BTreeMap, HashMap};

// PD codes and Jones polynomials from the Knot Atlas, https://katlas.org/wiki/<name>
// (KnotTheory`; retrieved 2026-11-09). Jones polynomials as published there (variable q,
// V(unknot) = 1), encoded as (2 * exponent of q, coefficient).
const DATA: &[(&str, &str, &[(i32, i64)])] = &[
    ("3_1", "X[1,4,2,5] X[3,6,4,1] X[5,2,6,3]", &[(-8, -1), (-6, 1), (-2, 1)]),
    ("4_1", "X[4,2,5,1] X[8,6,1,5] X[6,3,7,4] X[2,7,3,8]", &[(-4, 1), (-2, -1), (0, 1), (2, -1), (4, 1)]),
    (
        "5_1",
        "X[1,6,2,7] X[3,8,4,9] X[5,10,6,1] X[7,2,8,3] X[9,4,10,5]",
        &[(-14, -1), (-12, 1), (-10, -1), (-8, 1), (-4, 1)],
    ),
    (
        "5_2",
        "X[1,4,2,5] X[3,8,4,9] X[5,10,6,1] X[9,6,10,7] X[7,2,8,3]",
        &[(-12, -1), (-10, 1), (-8, -1), (-6, 2), (-4, -1), (-2, 1)],
    ),
    (
        "6_1",
        "X[1,4,2,5] X[7,10,8,11] X[3,9,4,8] X[9,3,10,2] X[5,12,6,1] X[11,6,12,7]",
        &[(-8, 1), (-6, -1), (-4, 1), (-2, -2), (0, 2), (2, -1), (4, 1)],
    ),
    (
        "6_2",
        "X[1,4,2,5] X[5,10,6,11] X[3,9,4,8] X[9,3,10,2] X[7,12,8,1] X[11,6,12,7]",
        &[(-10, 1), (-8, -2), (-6, 2), (-4, -2), (-2, 2), (0, -1), (2, 1)],
    ),
    (
        "6_3",
        "X[4,2,5,1] X[8,4,9,3] X[12,9,1,10] X[10,5,11,6] X[6,11,7,12] X[2,8,3,7]",
        &[(-6, -1), (-4, 2), (-2, -2), (0, 3), (2, -2), (4, 2), (6, -1)],
    ),
    (
        "7_4",
        "X[6,2,7,1] X[12,6,13,5] X[14,8,1,7] X[8,14,9,13] X[2,12,3,11] X[10,4,11,3] X[4,10,5,9]",
        &[(2, 1), (4, -2), (6, 3), (8, -2), (10, 3), (12, -2), (14, 1), (16, -1)],
    ),
    (
        "8_19",
        "X[4,2,5,1] X[8,4,9,3] X[9,15,10,14] X[5,13,6,12] X[13,7,14,6] X[11,1,12,16] X[15,11,16,10] X[2,8,3,7]",
        &[(6, 1), (10, 1), (16, -1)],
    ),
    (
        "8_20",
        "X[4,2,5,1] X[8,4,9,3] X[5,12,6,13] X[13,16,14,1] X[9,14,10,15] X[15,10,16,11] X[11,6,12,7] X[2,8,3,7]",
        &[(-10, -1), (-8, 1), (-6, -1), (-4, 2), (-2, -1), (0, 2), (2, -1)],
    ),
    (
        "9_42",
        "X[1,4,2,5] X[5,10,6,11] X[3,9,4,8] X[9,3,10,2] X[16,12,17,11] X[14,7,15,8] X[6,15,7,16] \
         X[18,14,1,13] X[12,18,13,17]",
        &[(-6, 1), (-4, -1), (-2, 1), (0, -1), (2, 1), (4, -1), (6, 1)],
    ),
    (
        "10_132",
        "X[4,2,5,1] X[8,4,9,3] X[5,12,6,13] X[15,18,16,19] X[9,16,10,17] X[17,10,18,11] X[13,20,14,1] \
         X[19,14,20,15] X[11,6,12,7] X[2,8,3,7]",
        &[(-14, -1), (-12, 1), (-10, -1), (-8, 1), (-4, 1)],
    ),
    (
        "K11n34",
        "X[4,2,5,1] X[8,4,9,3] X[12,5,13,6] X[2,8,3,7] X[9,17,10,16] X[11,18,12,19] X[6,13,7,14] \
         X[15,20,16,21] X[17,1,18,22] X[19,14,20,15] X[21,10,22,11]",
        &[(-12, 1), (-10, -2), (-8, 2), (-6, -2), (-4, 1), (2, 2), (4, -2), (6, 2), (8, -1)],
    ),
    (
        "K11n42",
        "X[4,2,5,1] X[8,4,9,3] X[12,5,13,6] X[2,8,3,7] X[9,18,10,19] X[11,21,12,20] X[6,13,7,14] \
         X[15,10,16,11] X[17,22,18,1] X[19,15,20,14] X[21,16,22,17]",
        &[(-12, 1), (-10, -2), (-8, 2), (-6, -2), (-4, 1), (2, 2), (4, -2), (6, 2), (8, -1)],
    ),
    ("L2a1", "X[4,1,3,2] X[2,3,1,4]", &[(-5, -1), (-1, -1)]),
    ("L4a1", "X[6,1,7,2] X[8,3,5,4] X[2,5,3,6] X[4,7,1,8]", &[(-9, -1), (-5, -1), (-3, 1), (-1, -1)]),
    (
        "L5a1",
        "X[6,1,7,2] X[10,7,5,8] X[4,5,1,6] X[2,10,3,9] X[8,4,9,3]",
        &[(-7, 1), (-5, -2), (-3, 1), (-1, -2), (1, 1), (3, -1)],
    ),
    (
        "L6a4",
        "X[6,1,7,2] X[12,8,9,7] X[4,12,1,11] X[10,5,11,6] X[8,4,5,3] X[2,9,3,10]",
        &[(-6, -1), (-4, 3), (-2, -2), (0, 4), (2, -2), (4, 3), (6, -1)],
    ),
];

/// Number of components of each test link (Knot Atlas).
fn expected_components(name: &str) -> usize {
    match name {
        "L2a1" | "L4a1" | "L5a1" => 2,
        "L6a4" => 3,
        _ => 1,
    }
}

// ---------------------------------------------------------------------------
// Laurent polynomials (exponent -> coefficient)
// ---------------------------------------------------------------------------

type Poly = BTreeMap<i32, i64>;

fn p_mono(e: i32, c: i64) -> Poly {
    let mut p = Poly::new();
    if c != 0 {
        p.insert(e, c);
    }
    p
}

fn p_add_to(a: &mut Poly, b: &Poly) {
    for (&e, &c) in b {
        let v = a.entry(e).or_insert(0);
        *v += c;
        if *v == 0 {
            a.remove(&e);
        }
    }
}

fn p_mul(a: &Poly, b: &Poly) -> Poly {
    let mut r = Poly::new();
    for (&ea, &ca) in a {
        for (&eb, &cb) in b {
            p_add_to(&mut r, &p_mono(ea + eb, ca * cb));
        }
    }
    r
}

fn p_pow(a: &Poly, k: usize) -> Poly {
    let mut r = p_mono(0, 1);
    for _ in 0..k {
        r = p_mul(&r, a);
    }
    r
}

/// δ = -A^2 - A^{-2}
fn delta_a() -> Poly {
    let mut d = p_mono(2, -1);
    p_add_to(&mut d, &p_mono(-2, -1));
    d
}

/// exact division by δ = -A^2 - A^{-2}
fn div_delta(p: &Poly) -> Poly {
    let mut r = p.clone();
    let mut q = Poly::new();
    while let Some((&top, &c)) = r.iter().next_back() {
        // c A^top = (-c A^{top-2}) * (-A^2) + ...
        let t = p_mono(top - 2, -c);
        p_add_to(&mut q, &t);
        let sub = p_mul(&t, &delta_a());
        let neg: Poly = sub.iter().map(|(&e, &v)| (e, -v)).collect();
        p_add_to(&mut r, &neg);
        assert!(r.keys().next_back().is_none_or(|&e| e < top), "division failed");
        if r.keys().next().is_some_and(|&e| e < -10_000) {
            panic!("not divisible by delta");
        }
    }
    q
}

// ---------------------------------------------------------------------------
// independent PD orientation and Kauffman bracket
// ---------------------------------------------------------------------------

/// Orientation by walking components; direction fixed by the under-strands (i -> k).
/// Returns (sign of each crossing (+1/-1), component of each label numbered by smallest label).
fn orient(x: &[[usize; 4]]) -> (Vec<i64>, BTreeMap<usize, usize>) {
    let mut occ: BTreeMap<usize, Vec<(usize, usize)>> = BTreeMap::new();
    for (c, xc) in x.iter().enumerate() {
        for s in 0..4 {
            occ.entry(xc[s]).or_default().push((c, s));
        }
    }
    assert!(occ.values().all(|o| o.len() == 2));
    // incoming[(c,s)] = true if the edge at slot s enters crossing c
    let mut incoming: HashMap<(usize, usize), bool> = HashMap::new();
    let mut comp: BTreeMap<usize, usize> = BTreeMap::new();
    let labels: Vec<usize> = occ.keys().copied().collect();
    for &l0 in &labels {
        if comp.contains_key(&l0) {
            continue;
        }
        let nc = comp.values().copied().max().map_or(0, |m| m + 1);
        // walk: (label, slot where it is entered)
        let mut walk: Vec<(usize, (usize, usize))> = vec![];
        let (mut l, mut at) = (l0, occ[&l0][0]);
        loop {
            walk.push((l, at));
            let out = (at.0, at.1 ^ 2);
            let l2 = x[out.0][out.1];
            let o = &occ[&l2];
            at = if o[0] == out { o[1] } else { o[0] };
            l = l2;
            if l == l0 && at == occ[&l0][0] {
                break;
            }
        }
        let fwd = walk.iter().any(|w| w.1 .1 == 0);
        let bwd = walk.iter().any(|w| w.1 .1 == 2);
        assert!(!(fwd && bwd));
        let rev = if fwd || bwd {
            bwd
        } else {
            let m = walk.len();
            let inc = (0..m).filter(|&k| walk[(k + 1) % m].0 == walk[k].0 + 1).count();
            let dec = (0..m).filter(|&k| walk[k].0 == walk[(k + 1) % m].0 + 1).count();
            dec > inc
        };
        for &(l, at) in &walk {
            comp.insert(l, nc);
            let o = &occ[&l];
            let other = if o[0] == at { o[1] } else { o[0] };
            incoming.insert(at, !rev);
            incoming.insert(other, rev);
        }
    }
    let signs = (0..x.len()).map(|c| if incoming[&(c, 3)] { 1 } else { -1 }).collect();
    (signs, comp)
}

/// KnotTheory's PositiveQ rule (valid when labels increase along components).
fn knottheory_positive(x: &[usize; 4]) -> bool {
    let (i, j, k, l) = (x[0], x[1], x[2], x[3]);
    i == j || k == l || j == l + 1 || l > j + 1
}

/// Kauffman bracket <D> in A, every closed loop counted as δ (so <unknot diagram> = δ).
/// A-smoothing of X[i,j,k,l] joins (i,j),(k,l); B-smoothing joins (i,l),(j,k).
fn bracket(x: &[[usize; 4]]) -> Poly {
    let mut states: HashMap<Vec<(usize, usize)>, Poly> = HashMap::new();
    states.insert(vec![], p_mono(0, 1));
    let mut open: HashMap<usize, usize> = HashMap::new(); // labels seen once -> count
    let mut done = vec![false; x.len()];
    for _ in 0..x.len() {
        // greedy: crossing with most labels already open
        let c = (0..x.len())
            .filter(|&c| !done[c])
            .max_by_key(|&c| (x[c].iter().filter(|l| open.contains_key(l)).count(), std::cmp::Reverse(c)))
            .unwrap();
        done[c] = true;
        for &l in &x[c] {
            *open.entry(l).or_insert(0) += 1;
        }
        open.retain(|_, v| *v < 2);
        let [i, j, k, l] = x[c];
        let mut next: HashMap<Vec<(usize, usize)>, Poly> = HashMap::new();
        for (key, coef) in &states {
            for (arcs, aexp) in [([(i, j), (k, l)], 1), ([(i, l), (j, k)], -1)] {
                let mut m: HashMap<usize, usize> = HashMap::new();
                for &(a, b) in key {
                    m.insert(a, b);
                    m.insert(b, a);
                }
                let mut loops = 0;
                for (u, v) in arcs {
                    if add_arc(&mut m, u, v) {
                        loops += 1;
                    }
                }
                let mut nk: Vec<(usize, usize)> = m.iter().filter(|(a, b)| a < b).map(|(&a, &b)| (a, b)).collect();
                nk.sort_unstable();
                let t = p_mul(&p_mul(coef, &p_mono(aexp, 1)), &p_pow(&delta_a(), loops));
                p_add_to(next.entry(nk).or_default(), &t);
            }
        }
        states = next;
    }
    assert_eq!(states.len(), 1);
    states.remove(&vec![]).unwrap()
}

/// Join an arc u–v (labels) into the matching m; returns true if a loop closed.
fn add_arc(m: &mut HashMap<usize, usize>, u: usize, v: usize) -> bool {
    if u == v {
        return true;
    }
    let (mut a, mut b) = (u, v);
    if let Some(w) = m.remove(&a) {
        m.remove(&w);
        if w == b {
            return true;
        }
        a = w;
    }
    if let Some(w) = m.remove(&b) {
        m.remove(&w);
        if w == a {
            return true;
        }
        b = w;
    }
    m.insert(a, b);
    m.insert(b, a);
    false
}

/// Jones polynomial (2 * exponent of q -> coefficient), V(unknot) = 1.
fn jones(x: &[[usize; 4]]) -> Poly {
    if x.is_empty() {
        return p_mono(0, 1);
    }
    let (signs, _) = orient(x);
    let w: i64 = signs.iter().sum();
    let b = div_delta(&bracket(x));
    // (-A^3)^{-w} = (-1)^w A^{-3w}
    let f = p_mul(&b, &p_mono((-3 * w) as i32, if w % 2 == 0 { 1 } else { -1 }));
    // A^k -> q^{-k/4}: doubled exponent -k/2
    f.iter()
        .map(|(&k, &c)| {
            assert!(k % 2 == 0, "odd power of A");
            (-k / 2, c)
        })
        .collect()
}

/// Jones polynomial of a braid closure (crossingless strands = split unknots).
fn jones_braid(b: &BraidWord) -> Poly {
    let n = b.strands;
    let mut touched = vec![false; n];
    for g in &b.word {
        touched[g.i] = true;
        touched[g.i + 1] = true;
    }
    let free = touched.iter().filter(|&&t| !t).count();
    let mut idx = vec![0; n];
    let mut k = 0;
    for p in 0..n {
        idx[p] = k;
        if touched[p] {
            k += 1;
        }
    }
    // V(L ⊔ O) = (-q^{1/2} - q^{-1/2}) V(L)
    let mut dq = p_mono(1, -1);
    p_add_to(&mut dq, &p_mono(-1, -1));
    if k == 0 {
        return p_pow(&dq, free - 1);
    }
    let word: Vec<Gen> = b.word.iter().map(|g| Gen { i: idx[g.i], positive: g.positive }).collect();
    let red = BraidWord { strands: k, word, components: vec![0; k] };
    p_mul(&jones(&braid_to_pd(&red).crossings), &p_pow(&dq, free))
}

/// The sub-braid formed by the strands whose component is in `keep`.
fn restrict_braid(b: &BraidWord, keep: &[usize]) -> BraidWord {
    let mut comps = b.components.clone();
    let mut word = vec![];
    for g in &b.word {
        let (ca, cb) = (comps[g.i], comps[g.i + 1]);
        if keep.contains(&ca) && keep.contains(&cb) {
            let below = comps[..g.i].iter().filter(|c| keep.contains(c)).count();
            word.push(Gen { i: below, positive: g.positive });
        }
        comps.swap(g.i, g.i + 1);
    }
    let components: Vec<usize> = b.components.iter().copied().filter(|c| keep.contains(c)).collect();
    BraidWord { strands: components.len(), word, components }
}

/// The sub-PD formed by the components in `keep` (other components deleted);
/// returns (crossings, number of kept components left without crossings).
fn restrict_pd(x: &[[usize; 4]], keep: &[usize]) -> (Vec<[usize; 4]>, usize) {
    let (_, comp) = orient(x);
    let mut parent: BTreeMap<usize, usize> = comp.keys().map(|&l| (l, l)).collect();
    fn find(p: &mut BTreeMap<usize, usize>, mut a: usize) -> usize {
        while p[&a] != a {
            a = p[&a];
        }
        a
    }
    let mut out = vec![];
    for xc in x {
        let under = keep.contains(&comp[&xc[0]]);
        let over = keep.contains(&comp[&xc[1]]);
        match (under, over) {
            (true, true) => out.push(*xc),
            (true, false) => {
                let (a, b) = (find(&mut parent, xc[0]), find(&mut parent, xc[2]));
                parent.insert(a, b);
            }
            (false, true) => {
                let (a, b) = (find(&mut parent, xc[1]), find(&mut parent, xc[3]));
                parent.insert(a, b);
            }
            (false, false) => {}
        }
    }
    let out: Vec<[usize; 4]> = out.iter().map(|xc| xc.map(|l| find(&mut parent, l))).collect();
    let mut with_crossings: Vec<usize> = out.iter().flatten().map(|l| comp[l]).collect();
    with_crossings.sort_unstable();
    with_crossings.dedup();
    let loops = keep.iter().filter(|c| !with_crossings.contains(c)).count();
    (out, loops)
}

fn jones_sub_pd(x: &[[usize; 4]], keep: &[usize]) -> Poly {
    let (sub, loops) = restrict_pd(x, keep);
    let mut dq = p_mono(1, -1);
    p_add_to(&mut dq, &p_mono(-1, -1));
    if sub.is_empty() {
        return p_pow(&dq, loops - 1);
    }
    p_mul(&jones(&sub), &p_pow(&dq, loops))
}

/// Linking numbers lk[a][b] (a < b) of a PD code, components numbered by smallest label.
fn linking_pd(x: &[[usize; 4]]) -> BTreeMap<(usize, usize), i64> {
    let (signs, comp) = orient(x);
    let mut lk = BTreeMap::new();
    for (c, xc) in x.iter().enumerate() {
        let (a, b) = (comp[&xc[0]], comp[&xc[1]]);
        if a != b {
            *lk.entry((a.min(b), a.max(b))).or_insert(0) += signs[c];
        }
    }
    lk.into_iter().filter(|(_, v)| *v != 0).map(|(k, v)| (k, v / 2)).collect()
}

/// Linking numbers of a braid closure using its component labels.
fn linking_braid(b: &BraidWord) -> BTreeMap<(usize, usize), i64> {
    let mut comps = b.components.clone();
    let mut lk = BTreeMap::new();
    for g in &b.word {
        let (a, c) = (comps[g.i], comps[g.i + 1]);
        if a != c {
            *lk.entry((a.min(c), a.max(c))).or_insert(0) += if g.positive { 1 } else { -1 };
        }
        comps.swap(g.i, g.i + 1);
    }
    lk.into_iter().filter(|(_, v)| *v != 0).map(|(k, v)| (k, v / 2)).collect()
}

/// Full comparison of a PD code with a braid claimed to represent it.
fn check_braid_against_pd(name: &str, x: &[[usize; 4]], b: &BraidWord) {
    let (_, comp) = orient(x);
    let ncomp = comp.values().copied().max().map_or(1, |m| m + 1);
    // closure consistency of the component labels
    assert_eq!(b.components.len(), b.strands, "{name}");
    assert_eq!(b.top_components(), b.components, "{name}: components not closure-consistent");
    let perm = b.permutation();
    for p in 0..b.strands {
        assert_eq!(b.components[perm[p]], b.components[p], "{name}");
    }
    // components of the braid = cycles of its permutation, labelled 0..ncomp
    assert_eq!(b.num_components(), ncomp, "{name}: number of components");
    let mut labels: Vec<usize> = b.components.clone();
    labels.sort_unstable();
    labels.dedup();
    assert_eq!(labels, (0..ncomp).collect::<Vec<_>>(), "{name}: component labels");
    // Jones polynomial of the closure
    let jb = jones_braid(b);
    assert_eq!(jb, jones(x), "{name}: Jones polynomial changed (braid {:?})", b.signed());
    // component identities: linking numbers and Jones polynomial of every component
    assert_eq!(linking_braid(b), linking_pd(x), "{name}: linking numbers");
    if ncomp > 1 {
        for c in 0..ncomp {
            assert_eq!(jones_braid(&restrict_braid(b, &[c])), jones_sub_pd(x, &[c]), "{name}: component {c}");
            let others: Vec<usize> = (0..ncomp).filter(|&d| d != c).collect();
            assert_eq!(jones_braid(&restrict_braid(b, &others)), jones_sub_pd(x, &others), "{name}: without {c}");
        }
    }
}

fn pd_of(s: &str) -> Pd {
    parse_pd(s).unwrap_or_else(|e| panic!("parse {s}: {e}"))
}

// ---------------------------------------------------------------------------
// tests
// ---------------------------------------------------------------------------

#[test]
fn bracket_basics() {
    // positive kink X[1,1,2,2]: <D> = -A^3 δ, V = 1
    assert_eq!(div_delta(&bracket(&[[1, 1, 2, 2]])), p_mono(3, -1));
    assert_eq!(jones(&[[1, 1, 2, 2]]), p_mono(0, 1));
    assert_eq!(orient(&[[1, 1, 2, 2]]).0, vec![1]);
    // right-handed trefoil = closure of σ_1^3 (left strand over): V = q + q^3 - q^4
    let b = BraidWord::from_signed(2, &[1, 1, 1]);
    let want: Poly = [(2, 1), (6, 1), (8, -1)].into_iter().collect();
    assert_eq!(jones(&braid_to_pd(&b).crossings), want);
    // positive Hopf link σ_1^2: V = -q^{1/2} - q^{5/2}
    let b = BraidWord::from_signed(2, &[1, 1]);
    let want: Poly = [(1, -1), (5, -1)].into_iter().collect();
    assert_eq!(jones(&braid_to_pd(&b).crossings), want);
}

#[test]
fn jones_matches_knot_atlas() {
    // validates the independent Jones implementation, the sign rule and chirality
    for (name, s, want) in DATA {
        let pd = pd_of(s);
        let want: Poly = want.iter().copied().collect();
        assert_eq!(jones(&pd.crossings), want, "{name}");
    }
}

#[test]
fn signs_and_components_from_pd() {
    for (name, s, _) in DATA {
        let pd = pd_of(s);
        let (signs, comp) = orient(&pd.crossings);
        let ours: Vec<i64> = pd.crossing_signs().iter().map(|&p| if p { 1 } else { -1 }).collect();
        assert_eq!(ours, signs, "{name}: signs");
        // KnotTheory's PositiveQ rule (Knot Atlas labels increase along components)
        let kt: Vec<i64> = pd.crossings.iter().map(|x| if knottheory_positive(x) { 1 } else { -1 }).collect();
        assert_eq!(kt, signs, "{name}: KnotTheory PositiveQ");
        assert_eq!(pd.label_components(), comp, "{name}: components");
        assert_eq!(pd.num_components(), expected_components(name), "{name}");
        assert_eq!(pd.writhe(), signs.iter().sum::<i64>());
    }
}

#[test]
fn vogel_braids_preserve_link() {
    for (name, s, _) in DATA {
        let pd = pd_of(s);
        let b = pd_to_braid_vogel(&pd);
        assert_eq!(b.strands, pd.seifert_circles(), "{name}: Vogel keeps the Seifert circles");
        assert_eq!(b.writhe(), pd.writhe(), "{name}: Vogel moves are Reidemeister II");
        check_braid_against_pd(name, &pd.crossings, &b);
    }
}

#[test]
fn simplified_braids_preserve_link() {
    for (name, s, _) in DATA {
        let pd = pd_of(s);
        let b = pd_to_braid(&pd);
        check_braid_against_pd(name, &pd.crossings, &b);
        // round trip through PD again
        let b2 = pd_to_braid(&braid_to_pd(&b));
        assert!(b2.strands <= b.strands, "{name}");
        assert_eq!(jones_braid(&b2), jones(&pd.crossings), "{name}");
    }
}

#[test]
fn report_strand_counts() {
    println!("{:<8} {:>9} {:>10} {:>9} {:>9}", "link", "crossings", "Seifert", "strands", "length");
    for (name, s, _) in DATA {
        let pd = pd_of(s);
        let v = pd_to_braid_vogel(&pd);
        let b = pd_to_braid(&pd);
        println!("{:<8} {:>9} {:>10} {:>9} {:>9}   {:?}  comps {:?}", name, pd.crossings.len(), v.strands, b.strands, b.word.len(), b.signed(), b.components);
        if *name == "K11n34" {
            println!("Conway knot 11n34: {} strands, word {:?}", b.strands, b.signed());
        }
    }
}

// Braid representatives from KnotInfo / LinkInfo (database_knotinfo csv, columns
// braid_index / braid_notation; retrieved 2026-11-09). KnotInfo's chirality convention
// differs from the Knot Atlas for some knots, so Jones polynomials are compared up to
// q <-> 1/q. For L4a1 the LinkInfo orientation L4a1{0} (braid index 3) is listed; the
// orientation L4a1{1} is the torus link T(2,4) with braid index 2.
const KNOTINFO_BRAIDS: &[(&str, usize, &[i64])] = &[
    ("3_1", 2, &[1, 1, 1]),
    ("4_1", 3, &[1, -2, 1, -2]),
    ("5_1", 2, &[1, 1, 1, 1, 1]),
    ("5_2", 3, &[1, 1, 1, 2, -1, 2]),
    ("6_1", 4, &[1, 1, 2, -1, -3, 2, -3]),
    ("6_2", 3, &[1, 1, 1, -2, 1, -2]),
    ("6_3", 3, &[1, 1, -2, 1, -2, -2]),
    ("7_4", 4, &[1, 1, 2, -1, 2, 2, 3, -2, 3]),
    ("8_19", 3, &[1, 1, 1, 2, 1, 1, 1, 2]),
    ("8_20", 3, &[1, 1, 1, -2, -1, -1, -1, -2]),
    ("9_42", 4, &[-1, -1, -1, 2, 1, 1, -3, 2, -3]),
    ("10_132", 4, &[1, 1, 1, -2, -1, -1, -2, -3, 2, -3, -3]),
    ("K11n34", 4, &[-1, -1, 2, -1, 2, -1, 3, -2, -2, 3, 3]),
    ("K11n42", 4, &[-1, -1, -1, -1, 2, 2, 1, -3, 2, 2, -3, 2, -3]),
    ("L2a1", 2, &[-1, -1]),
    ("L4a1", 3, &[-2, -2, -1, 2, -1]),
    ("L5a1", 3, &[-1, 2, -1, 2, -1]),
    ("L6a4", 3, &[1, -2, 1, -2, 1, -2]),
];

fn mirror(p: &Poly) -> Poly {
    p.iter().map(|(&e, &c)| (-e, c)).collect()
}

#[test]
fn knotinfo_braids_agree_up_to_mirror() {
    for (name, n, w) in KNOTINFO_BRAIDS {
        let (_, s, want) = DATA.iter().find(|d| d.0 == *name).unwrap();
        let want: Poly = want.iter().copied().collect();
        let got = jones_braid(&BraidWord::from_signed(*n, w));
        assert!(got == want || mirror(&got) == want, "{name}: {got:?} vs {want:?}");
        // our braid has no more strands than the published braid index
        let b = pd_to_braid(&pd_of(s));
        assert!(b.strands <= *n, "{name}: {} strands, braid index {n}", b.strands);
    }
}

/// Small deterministic RNG (xorshift64*).
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545F4914F6CDD1D)
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

#[test]
fn random_braids_roundtrip() {
    let mut rng = Rng(0x9E3779B97F4A7C15);
    let mut tested = 0;
    while tested < 150 {
        let n = 2 + rng.below(4) as usize;
        let len = 3 + rng.below(12) as usize;
        let w: Vec<i64> = (0..len)
            .map(|_| {
                let k = 1 + rng.below(n as u64 - 1) as i64;
                if rng.below(2) == 0 { k } else { -k }
            })
            .collect();
        let b0 = BraidWord::from_signed(n, &w);
        let mut touched = vec![false; n];
        for g in &b0.word {
            touched[g.i] = true;
            touched[g.i + 1] = true;
        }
        if touched.iter().any(|t| !t) {
            continue;
        }
        let pd = braid_to_pd(&b0);
        assert!(pd.validate().is_ok());
        // braid_to_pd numbers components by first bottom position, as from_signed does
        assert_eq!(linking_pd(&pd.crossings), linking_braid(&b0), "{w:?}");
        assert_eq!(jones(&pd.crossings), jones_braid(&b0));
        let name = format!("{w:?}");
        let b = pd_to_braid_vogel(&pd);
        assert_eq!(b.writhe(), b0.writhe());
        check_braid_against_pd(&name, &pd.crossings, &b);
        let b = pd_to_braid(&pd);
        check_braid_against_pd(&name, &pd.crossings, &b);
        tested += 1;
    }
}

#[test]
fn random_diagrams_need_vogel_moves() {
    // Mirror-reflect braid closures diagram-wise in a way that is not braided: take the
    // PD code of a braid closure and reverse the orientation of a subset of components
    // (relabel each such component backwards). The result is a different link diagram
    // with incoherent Seifert circles, so Vogel moves are required.
    let mut rng = Rng(12345);
    let mut moved = 0;
    let mut tested = 0;
    while tested < 120 {
        let n = 2 + rng.below(4) as usize;
        let len = 4 + rng.below(10) as usize;
        let w: Vec<i64> = (0..len)
            .map(|_| {
                let k = 1 + rng.below(n as u64 - 1) as i64;
                if rng.below(2) == 0 { k } else { -k }
            })
            .collect();
        let b0 = BraidWord::from_signed(n, &w);
        if (0..n).any(|p| !b0.word.iter().any(|g| g.i == p || g.i + 1 == p)) {
            continue;
        }
        let ncomp = b0.num_components();
        if ncomp < 2 {
            continue;
        }
        let pd = braid_to_pd(&b0);
        let (_, comp) = orient(&pd.crossings);
        let flip_mask = 1 + rng.below((1u64 << ncomp) - 1); // non-empty subset
        let flip = |c: usize| (flip_mask >> c) & 1 == 1 && (flip_mask != (1u64 << ncomp) - 1 || c > 0);
        // reverse a component: relabel its labels in reverse order, and for crossings where
        // it is the under-strand rotate the crossing by two slots (new incoming = old outgoing)
        let mut by_comp: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for (&l, &c) in &comp {
            by_comp.entry(c).or_default().push(l);
        }
        let mut relabel: BTreeMap<usize, usize> = comp.keys().map(|&l| (l, l)).collect();
        for (&c, ls) in &by_comp {
            if flip(c) {
                // labels min..=max along the component; reverse: l -> min + max - l, then
                // the edge previously min..min+1 etc. Reversal maps edge l (from node l to
                // node l+1) to an edge oriented the other way.
                let (lo, hi) = (*ls.iter().min().unwrap(), *ls.iter().max().unwrap());
                for &l in ls {
                    relabel.insert(l, lo + hi - l);
                }
            }
        }
        let x: Vec<[usize; 4]> = pd
            .crossings
            .iter()
            .map(|xc| {
                let y = xc.map(|l| relabel[&l]);
                if flip(comp[&xc[0]]) { [y[2], y[3], y[0], y[1]] } else { y }
            })
            .collect();
        let pd2 = Pd { crossings: x };
        pd2.validate().unwrap_or_else(|e| panic!("{w:?} mask {flip_mask}: {e}"));
        let s0 = pd2.seifert_circles();
        let name = format!("{w:?} reversed mask {flip_mask}");
        let b = pd_to_braid_vogel(&pd2);
        if b.word.len() > pd2.crossings.len() {
            moved += 1;
        }
        assert_eq!(b.strands, s0, "{name}");
        assert_eq!(b.writhe(), pd2.writhe(), "{name}");
        check_braid_against_pd(&name, &pd2.crossings, &b);
        check_braid_against_pd(&name, &pd2.crossings, &pd_to_braid(&pd2));
        tested += 1;
    }
    println!("{moved} of {tested} diagrams needed Vogel moves");
    assert!(moved > 50, "only {moved} diagrams needed Vogel moves");
}

#[test]
fn parser_formats_and_errors() {
    let a = parse_pd("[[1,4,2,5],[3,6,4,1],[5,2,6,3]]").unwrap();
    let b = parse_pd("X[1,4,2,5] X[3,6,4,1] X[5,2,6,3]").unwrap();
    let c = parse_pd("PD[X[1, 4, 2, 5], X[3, 6, 4, 1], X[5, 2, 6, 3]]").unwrap();
    assert_eq!(a, b);
    assert_eq!(a, c);
    assert_eq!(parse_pd(&a.to_string()).unwrap(), a);
    assert!(parse_pd("[[1,4,2,5],[3,6,4,1]]").is_err()); // labels not paired
    assert!(parse_pd("[[1,4,2,5],[3,6,4]]").is_err());
    assert!(parse_pd("[[1,4,2,y]]").is_err());
    // under-strands disagreeing on the orientation: 1->2 and 1<-2
    assert!(parse_pd("X[1,3,2,4] X[1,4,2,3]").is_err());
    // empty code: unknot
    let e = parse_pd("[]").unwrap();
    let be = pd_to_braid(&e);
    assert_eq!((be.strands, be.word.len(), be.components.clone()), (1, 0, vec![0]));
    // a single kink is the unknot
    let k = pd_to_braid(&parse_pd("X[1,1,2,2]").unwrap());
    assert_eq!((k.strands, k.word.len()), (1, 0));
    let k = pd_to_braid(&parse_pd("X[2,1,1,2]").unwrap());
    assert_eq!((k.strands, k.word.len()), (1, 0));
}

#[test]
fn split_diagrams() {
    // trefoil ⊔ Hopf link drawn side by side
    let s = "X[1,4,2,5] X[3,6,4,1] X[5,2,6,3] X[10,7,9,8] X[8,9,7,10]";
    let pd = pd_of(s);
    assert_eq!(pd.num_components(), 3);
    let b = pd_to_braid(&pd);
    check_braid_against_pd("3_1 + L2a1", &pd.crossings, &b);
}
