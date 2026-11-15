//! Foam_N pruning of the exact colored engine (`XOpts::prune_n`, rule `schur::Schur::killed_at`):
//! categorified validation by computation.
//!
//! * Pruned scan at N (killed summands dropped after every decomposition, Rickard complexes
//!   built by `crossing_cx(.., N)`) and the unpruned N-independent scan give identical gl_N
//!   Poincaré polynomials, unreduced and reduced: (2,2) words of length 1-5 (all sign patterns)
//!   at N = 2, 3; (2,1), (1,2), (3,1), (1,3) Hopf and T(2,4) at N = 2..4; (2,2) T(2,n), n <= 5, at
//!   N = 2..4; random mixed-sign 2-strand words with labels <= 2; (3,2) Hopf at N = 3, 4 (the
//!   unpruned scan takes about 40 s). A scan pruned at N closes correctly at every N' < N too.
//! * The truncated K_0 class (asserted inside the scan) equals `schur::decat_scan`'s class with
//!   the killed types removed; pruning drops summands in the tested cases.
//! * The closure check of d∘d components (`closure_kills_map`) is exercised on hand-made maps.
//! * New results (3-colored, (3,3)): Hopf link and trefoil at N = 3, 4 against the Hopf formula,
//!   Wang (tests/data/wang_ref.txt, k = 3), `schur::decat_chi`, both primes, and duality at
//!   N = 4 (uncolored gl_4); 3-colored T(2,5) and mixed-sign words at N = 4 against Wang,
//!   duality (knots, unreduced and reduced) and σσ^{-1} cancellation (links). Some of these have
//!   d∘d components in the ideal of killed types, whose closures are checked to vanish. N = 5
//!   is `#[ignore]` (Hom spaces of hours, see the test).

use krhom::colored::cclose::{ClosureDataC, EPoly, ERing, TypeBasisC};
use krhom::input::crossing::Gen;
use krhom::reference::cube::local_cx;
use krhom::uncolored::engine2::{khr_engine2_closure, Closure};
use krhom::algebra::field::{Field, F31, F32b};
use krhom::algebra::laurent::{LPoly, Poincare};
use krhom::combinat::schur::{decat_chi, decat_scan, ScanOpts, Schur, Table};
use krhom::colored::xengine::*;
use rand::{rngs::StdRng, Rng, SeedableRng};
use std::any::Any;
use std::collections::{BTreeMap, HashMap};
use std::sync::{Arc, Mutex, OnceLock};

// ------------------------------------------------------------------ helpers

fn gens(w: &[i32]) -> Vec<Gen> {
    w.iter().map(|&g| Gen { i: (g.unsigned_abs() - 1) as usize, positive: g > 0 }).collect()
}

type ScanKey = (String, Vec<usize>, Vec<i32>, Option<usize>);
type ScanCell = Arc<OnceLock<Arc<dyn Any + Send + Sync>>>;

/// One scan per (field, labels, word, prune_n), shared by all tests of this binary.
fn scan<F: Field>(labels: &[usize], word: &[i32], prune: Option<usize>) -> Arc<XScan<F>> {
    static CACHE: OnceLock<Mutex<HashMap<ScanKey, ScanCell>>> = OnceLock::new();
    let cell = CACHE.get_or_init(Default::default).lock().unwrap().entry((F::name(), labels.to_vec(), word.to_vec(), prune)).or_default().clone();
    let any = cell
        .get_or_init(|| {
            let t0 = std::time::Instant::now();
            let opts = XOpts { prune_n: prune, ..Default::default() };
            let s = xscan::<F>(labels, &gens(word), &opts).unwrap_or_else(|e| panic!("xscan {:?} {:?} prune {:?} over {}: {}", labels, word, prune, F::name(), e));
            let dt = t0.elapsed().as_secs_f64();
            if dt > 2.0 {
                eprintln!("[scan] {:?} {:?} prune {:?} over {}: {} final objects, {:.2}s (hom_space calls {}, max {:.2}s)", labels, word, prune, F::name(), s.data.objs.len(), dt, s.hom_calls, s.hom_max_secs);
            }
            Arc::new(s) as Arc<dyn Any + Send + Sync>
        })
        .clone();
    any.downcast::<XScan<F>>().unwrap()
}

fn close<F: Field>(labels: &[usize], word: &[i32], prune: Option<usize>, nn: usize, reduced: bool) -> Poincare {
    scan::<F>(labels, word, prune).close(&[nn], reduced).remove(&nn).unwrap()
}

/// Pruned (at N, and at every N0 in `also` with N0 >= N) and unpruned scans agree at N,
/// unreduced and reduced; colors > N: the pruned scan is refused and the unpruned gives 0.
fn compare(labels: &[usize], word: &[i32], nn: usize, also: &[usize]) -> Poincare {
    let un = close::<F31>(labels, word, None, nn, false);
    if labels.iter().any(|&a| a > nn) {
        assert_eq!(un, Poincare::default());
        let e = xscan::<F31>(labels, &gens(word), &XOpts { prune_n: Some(nn), ..Default::default() }).err().unwrap();
        assert!(e.contains("below a color"), "{}", e);
        return un;
    }
    let pr = close::<F31>(labels, word, Some(nn), nn, false);
    assert_eq!(pr, un, "{:?} {:?} N={}: pruned != unpruned", labels, word, nn);
    let unr = close::<F31>(labels, word, None, nn, true);
    let prr = close::<F31>(labels, word, Some(nn), nn, true);
    assert_eq!(prr, unr, "{:?} {:?} N={}: pruned != unpruned (reduced)", labels, word, nn);
    for &n0 in also.iter().filter(|&&n0| n0 > nn) {
        assert_eq!(close::<F31>(labels, word, Some(n0), nn, false), un, "{:?} {:?}: pruned at N0={} closed at N={}", labels, word, n0, nn);
    }
    un
}

fn one() -> Poincare {
    let mut p = Poincare::default();
    p.add_term(0, 0, 1);
    p
}

/// "t^0q^-6 + 2t^2q^-8 + ..." (the format of `Poincare`'s Display and of wang.py)
fn parse_poincare(s: &str) -> Poincare {
    let mut p = Poincare::default();
    let s = s.trim();
    if s == "0" || s.is_empty() {
        return p;
    }
    for term in s.split(" + ") {
        let term = term.trim();
        let ti = term.find("t^").unwrap_or_else(|| panic!("bad term {:?}", term));
        let k: i128 = if ti == 0 { 1 } else { term[..ti].parse().unwrap() };
        let qi = term.find("q^").unwrap();
        let t: i32 = term[ti + 2..qi].parse().unwrap();
        let q: i32 = term[qi + 2..].parse().unwrap();
        p.add_term(t, q, k);
    }
    p
}

/// tests/data/wang_ref.txt (r3scripts/wang.py k N m: T(2,2m+1) colored ∧^k at gl_N).
fn wang_ref() -> &'static BTreeMap<(usize, usize, usize), Poincare> {
    static R: OnceLock<BTreeMap<(usize, usize, usize), Poincare>> = OnceLock::new();
    R.get_or_init(|| {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/wang_ref.txt");
        let text = std::fs::read_to_string(&path).unwrap();
        let mut out = BTreeMap::new();
        let mut lines = text.lines();
        while let Some(l) = lines.next() {
            if let Some(rest) = l.strip_prefix("k N m = ") {
                let v: Vec<usize> = rest.split_whitespace().map(|x| x.parse().unwrap()).collect();
                let poly = parse_poincare(lines.next().unwrap());
                let rank: i128 = lines.next().unwrap().strip_prefix("rank ").unwrap().trim().parse().unwrap();
                assert_eq!(poly.total_rank(), rank, "wang_ref {:?}", v);
                out.insert((v[0], v[1], v[2]), poly);
            }
        }
        out
    })
}

/// R3 V4 closed formula for the closure of σ_1^2 with colors (a, b) (as tests/xengine.rs).
fn hopf_formula(nn: usize, a: usize, b: usize) -> Poincare {
    let (n, a, b) = (nn as i64, a as i64, b as i64);
    let m = a.min(b);
    let mut out = Poincare::default();
    if a > n || b > n {
        return out;
    }
    let jmax = m - 0.max(a + b - n);
    for j in 0..=jmax {
        let k = m - j;
        let parts = [k, a - k, b - k, n - a - b + k];
        assert!(parts.iter().all(|&x| x >= 0));
        let mut den = LPoly::one();
        for &x in &parts {
            den = den.mul(&LPoly::qfact(x));
        }
        let mult = LPoly::qfact(n).div_exact(&den);
        let s = -a * b - m * n + 2 * m * m - j * (n + 2);
        for (&e, &c) in &mult.c {
            out.add_term(2 * j as i32, (s + e as i64) as i32, c);
        }
    }
    out
}

fn chi_oracle(labels: &[usize], word: &[i32], nn: usize) -> LPoly {
    let lam: Vec<u8> = labels.iter().map(|&a| a as u8).collect();
    let mut s = Schur::new(&lam);
    let res = decat_scan(&mut s, word, &ScanOpts { peak: false, every_stage: false, marked_block: None });
    decat_chi(&mut s, &res, word, nn)
}

// ------------------------------------------------------------------ pruned = unpruned

#[test]
fn pruned_equals_unpruned_22_all_words() {
    // all 62 words of length 1..5 on (2,2), N = 2, 3 (at N = 2 only the identity type survives)
    let t = std::time::Instant::now();
    let mut words: Vec<Vec<i32>> = vec![];
    for len in 1..=5usize {
        for mask in 0..(1u32 << len) {
            words.push((0..len).map(|i| if mask >> i & 1 == 1 { -1 } else { 1 }).collect());
        }
    }
    assert_eq!(words.len(), 62);
    let mut dropped = 0;
    for w in &words {
        for nn in [2usize, 3] {
            compare(&[2, 2], w, nn, &[3]);
            dropped += scan::<F31>(&[2, 2], w, Some(nn)).dropped_summands;
        }
    }
    assert!(dropped > 0);
    // at N = 2 the 2-colored theory of any (2,2) closure is q^0 (∧^2 C^2 is a line)
    for w in &words {
        assert_eq!(close::<F31>(&[2, 2], w, Some(2), 2, false), one(), "{:?}", w);
    }
    eprintln!("(2,2) words of length <= 5 at N = 2, 3: {} summands dropped, {:.1}s", dropped, t.elapsed().as_secs_f64());
}

#[test]
fn pruned_equals_unpruned_hopf_and_t24() {
    let t = std::time::Instant::now();
    for labels in [[2usize, 1], [1, 2], [3, 1], [1, 3]] {
        for word in [vec![1, 1], vec![-1, -1], vec![1, 1, 1, 1]] {
            for nn in 2..=4usize {
                let p = compare(&labels, &word, nn, &[4]);
                if word == [1, 1] {
                    assert_eq!(p, hopf_formula(nn, labels[0], labels[1]), "Hopf {:?} N={}", labels, nn);
                }
            }
        }
    }
    eprintln!("Hopf and T(2,4), colors (2,1), (1,2), (3,1), (1,3), N = 2..4: {:.1}s", t.elapsed().as_secs_f64());
}

#[test]
fn pruned_equals_unpruned_22_torus() {
    let t = std::time::Instant::now();
    let r = wang_ref();
    for n in 1..=5usize {
        for nn in 2..=4usize {
            let p = compare(&[2, 2], &vec![1; n], nn, &[4]);
            if n % 2 == 1 {
                if let Some(w) = r.get(&(2, nn, n / 2)) {
                    assert_eq!(p, *w, "2-colored T(2,{}) N={} vs wang.py", n, nn);
                }
            }
        }
    }
    eprintln!("(2,2) T(2,n), n <= 5, N = 2..4: {:.1}s", t.elapsed().as_secs_f64());
}

#[test]
fn pruned_equals_unpruned_random_words() {
    // seeded random mixed-sign words, labels <= 2 (even length when the colors differ)
    let t = std::time::Instant::now();
    let mut rng = StdRng::seed_from_u64(20261110);
    let mut count = 0;
    for labels in [[1usize, 1], [2, 1], [1, 2], [2, 2]] {
        let mut done = 0;
        while done < 6 {
            let len = rng.gen_range(2..=7usize);
            if labels[0] != labels[1] && len % 2 == 1 {
                continue;
            }
            let w: Vec<i32> = (0..len).map(|_| if rng.gen_bool(0.5) { 1 } else { -1 }).collect();
            if w.iter().all(|&g| g == w[0]) {
                continue;
            }
            for nn in 1..=4usize {
                compare(&labels, &w, nn, &[4]);
            }
            done += 1;
            count += 1;
        }
    }
    eprintln!("{} random mixed-sign words, N = 1..4: {:.1}s", count, t.elapsed().as_secs_f64());
}

#[test]
fn pruned_equals_unpruned_32_hopf() {
    // the unpruned (3,2) Hopf scan takes about 40 s (86 Hom spaces); N = 3, 4, and the formula
    let t = std::time::Instant::now();
    for nn in 3..=4usize {
        let p = compare(&[3, 2], &[1, 1], nn, &[]);
        assert_eq!(p, hopf_formula(nn, 3, 2), "(3,2) Hopf N={}", nn);
        let sp = scan::<F31>(&[3, 2], &[1, 1], Some(nn));
        eprintln!("(3,2) Hopf N={}: {} (rank {}); pruned scan {:.1}s, {} hom_space calls (max {:.1}s), {} summands dropped", nn, p, p.total_rank(), sp.timings.0, sp.hom_calls, sp.hom_max_secs, sp.dropped_summands);
    }
    let su = scan::<F31>(&[3, 2], &[1, 1], None);
    eprintln!("(3,2) Hopf unpruned scan {:.1}s, {} hom_space calls (max {:.1}s); total {:.1}s", su.timings.0, su.hom_calls, su.hom_max_secs, t.elapsed().as_secs_f64());
}

// ------------------------------------------------------------------ K_0, pruned Rickard complexes, d∘d

#[test]
fn truncated_k0_class() {
    for (labels, word) in [(vec![2usize, 2], vec![1, 1, 1]), (vec![2, 2], vec![1, -1, 1, 1]), (vec![2, 1], vec![1, 1]), (vec![3, 1], vec![1, 1, 1, 1])] {
        let lam: Vec<u8> = labels.iter().map(|&a| a as u8).collect();
        let mut s = Schur::new(&lam);
        let res = decat_scan(&mut s, &word, &ScanOpts { peak: false, every_stage: false, marked_block: None });
        for nn in *labels.iter().max().unwrap()..=4usize {
            let sc = scan::<F31>(&labels, &word, Some(nn));
            let mut want: Vec<(String, i32, i128)> = res.k0.iter().filter(|&&(p, _, _)| !s.killed_at(p, nn)).map(|&(p, k, c)| (s.table(p).fmt(), k, c)).collect();
            want.sort();
            let mut got = sc.k0.clone();
            got.sort();
            assert_eq!(got, want, "{:?} {:?} N={}: truncated K_0 class", labels, word, nn);
            for (tab, _, _) in &sc.final_objects {
                let p = s.intern(Table::parse(tab));
                assert!(!s.killed_at(p, nn), "killed type {} in the pruned complex", tab);
            }
            for st in &sc.stages {
                assert!(st.objects as u128 >= st.lower_bound);
            }
        }
    }
}

#[test]
fn pruned_rickard_complexes_are_truncations() {
    // local_cx(.., n = N) keeps W_0..W_kmax with the same differentials as the untruncated one
    for (c, d) in [(1usize, 1usize), (2, 1), (1, 2), (2, 2), (3, 1), (1, 3), (3, 2)] {
        for positive in [true, false] {
            let full = local_cx::<F31>(Gen { i: 0, positive }, c, d, c + d);
            for nn in c.max(d)..c + d {
                let pr = local_cx::<F31>(Gen { i: 0, positive }, c, d, nn);
                let kept = pr.cx.terms.len();
                assert_eq!(kept, nn + 1 - c.max(d), "({},{}) N={}: W_k is kept iff max(c,d) + k <= N", c, d, nn);
                for k in 0..kept {
                    assert_eq!(pr.cx.terms[k].rungs, full.cx.terms[k].rungs);
                    assert_eq!(pr.fsizes[k], full.fsizes[k]);
                }
                for j in 0..kept.saturating_sub(1) {
                    assert!(pr.maps[j].sub(&full.maps[j]).is_zero(), "({},{}) N={}: differential {} differs", c, d, nn, j);
                }
            }
        }
    }
}

#[test]
fn closure_kills_map_detects_nonzero_maps() {
    // R^(2,2) -> R^(2,2), multiplication by e_{0,1} (degree 2): at N = 2 the closure is
    // H*(Gr(2,2))^⊗2 = k and the map closes to 0; at N = 3, 4 it does not. Multiplication by
    // e_{0,1} − e_{1,1} does not either; e_{0,2} at N = 2 (degree 4) closes to 0.
    let lam = vec![2u8, 2];
    let ring = ERing::new(&lam);
    let var = |i: usize, j: usize| {
        let mut ex = vec![0u32; ring.nvars];
        ex[ring.var(i, j)] = 1;
        EPoly::<F31>::mono(ring.pack(&ex), F31::one())
    };
    let tb = TypeBasisC { ins: vec![0], outs: vec![0], gram: vec![vec![EPoly::one()]] };
    let data = |t0: i32| ClosureDataC::<F31> { lam: lam.clone(), objs: vec![(0, 0, t0), (0, 1, 0)], types: BTreeMap::from([(0u32, tb.clone())]), edges: vec![], crossings: vec![] };
    let e01 = vec![vec![var(0, 1)]];
    assert!(closure_kills_map(&data(2), 0, 1, &e01, 2, None));
    assert!(!closure_kills_map(&data(2), 0, 1, &e01, 3, None));
    assert!(!closure_kills_map(&data(2), 0, 1, &e01, 4, None));
    let diff = vec![vec![var(0, 1).sub(&var(1, 1))]];
    assert!(!closure_kills_map(&data(2), 0, 1, &diff, 3, None));
    assert!(closure_kills_map(&data(4), 0, 1, &vec![vec![var(0, 2)]], 2, None));
    // reduced (block 0 marked): e_{1,1} on the unmarked block at N = 3 closes nontrivially
    assert!(!closure_kills_map(&data(2), 0, 1, &vec![vec![var(1, 1)]], 3, Some(0)));
}

#[test]
fn pruning_statistics() {
    // pruning is effective: (2,2) T(2,5) at N = 3 drops the E(2)F(2) summands; the stage
    // counters add up; d∘d components in the ideal (if any) are closed to 0 in `close`
    let un = scan::<F31>(&[2, 2], &[1, 1, 1, 1, 1], None);
    let pr = scan::<F31>(&[2, 2], &[1, 1, 1, 1, 1], Some(3));
    assert!(pr.dropped_summands > 0);
    assert!(pr.stages.iter().map(|s| s.dropped).sum::<usize>() > 0);
    assert!(pr.hom_calls <= un.hom_calls);
    assert!(un.stages.iter().all(|s| s.dropped == 0 && s.d2_ideal == 0) && un.final_d2.is_empty());
    let mut d2 = 0;
    for (labels, word) in [(vec![2usize, 2], vec![1, -1, 1, -1, 1]), (vec![2, 2], vec![1, 1, -1, -1, 1]), (vec![2, 1], vec![1, -1, 1, 1]), (vec![1, 1], vec![1, -1, 1, -1, -1, 1])] {
        for nn in *labels.iter().max().unwrap()..=3usize {
            let sc = scan::<F31>(&labels, &word, Some(nn));
            d2 += sc.stages.iter().map(|s| s.d2_ideal).sum::<usize>();
        }
    }
    eprintln!("(2,2) T(2,5): unpruned {} hom_space calls / {} objects, pruned at N=3 {} hom_space calls / {} objects ({} summands dropped); d∘d components in the ideal over mixed words: {}", un.hom_calls, un.data.objs.len(), pr.hom_calls, pr.data.objs.len(), pr.dropped_summands, d2);
    // a pruned scan refuses N above its N
    assert!(std::panic::catch_unwind(|| pr.close(&[4], false)).is_err());
}

// ------------------------------------------------------------------ 3 colors: new results

fn three_colored_check<F: Field>(word: &[i32], nn: usize) -> Poincare {
    let t = std::time::Instant::now();
    let sc = scan::<F>(&[3, 3], word, Some(nn));
    let p = sc.close(&[nn], false).remove(&nn).unwrap();
    let red = sc.close(&[nn], true).remove(&nn).unwrap();
    eprintln!(
        "(3,3) {:?} N={} over {}: {} (rank {}), reduced {} (rank {}); scan {:.1}s + closure data {:.1}s, {} hom_space calls (max {:.1}s), {} summands dropped, final objects {}, total {:.1}s",
        word,
        nn,
        F::name(),
        p,
        p.total_rank(),
        red,
        red.total_rank(),
        sc.timings.0,
        sc.timings.1,
        sc.hom_calls,
        sc.hom_max_secs,
        sc.dropped_summands,
        sc.data.objs.len(),
        t.elapsed().as_secs_f64()
    );
    assert_eq!(p.euler(), chi_oracle(&[3, 3], word, nn), "(3,3) {:?} N={}: χ != decat_chi", word, nn);
    if word.len() % 2 == 1 {
        assert_eq!(p.euler(), red.euler().mul(&LPoly::qbinom(nn as i64, 3)), "(3,3) {:?} N={}: χ(unred) = [N 3] χ(red)", word, nn);
    }
    p
}

#[test]
fn three_colored_hopf_and_trefoil_n3_n4() {
    let t = std::time::Instant::now();
    let r = wang_ref();
    for nn in [3usize, 4] {
        let h = three_colored_check::<F31>(&[1, 1], nn);
        assert_eq!(h, hopf_formula(nn, 3, 3), "3-colored Hopf N={} vs formula", nn);
        assert_eq!(three_colored_check::<F32b>(&[1, 1], nn), h, "3-colored Hopf N={}: two primes", nn);
        let tr = three_colored_check::<F31>(&[1, 1, 1], nn);
        assert_eq!(tr, r[&(3, nn, 1)], "3-colored trefoil N={} vs wang.py", nn);
        assert_eq!(three_colored_check::<F32b>(&[1, 1, 1], nn), tr, "3-colored trefoil N={}: two primes", nn);
    }
    assert_eq!(wang_ref()[&(3, 3, 1)].total_rank(), 1);
    assert_eq!(wang_ref()[&(3, 4, 1)].total_rank(), 10);
    // via duality (expected): ∧^3 gl_4 = uncolored gl_4 of the reversed word (knot: no shift)
    let un = khr_engine2_closure::<F31>(2, &gens(&[1, 1, 1]), &[4], &[], 7, false, Closure::Exact).poincare.remove(&4).unwrap();
    assert_eq!(close::<F31>(&[3, 3], &[1, 1, 1], Some(4), 4, false), un, "via duality (expected): 3-colored trefoil at N=4 vs uncolored gl_4");
    // the mirror (negative trefoil) at N = 4
    let neg = close::<F31>(&[3, 3], &[-1, -1, -1], Some(4), 4, false);
    let mut mir = Poincare::default();
    for (&(tt, q), &k) in &r[&(3, 4, 1)].c {
        mir.add_term(-tt, -q, k);
    }
    assert_eq!(neg, mir, "3-colored negative trefoil at N=4 = mirror");
    eprintln!("3-colored Hopf and trefoil, N = 3, 4, both primes: {:.1}s", t.elapsed().as_secs_f64());
}

#[test]
fn three_colored_mixed_words_n4() {
    // ∧^3 gl_4 (∧^3 C^4 = dual of C^4): knots against uncolored gl_4 of the reversed word (via
    // duality (expected), unreduced and reduced), links against the positive word with the same
    // exponent sum (σσ^{-1} cancellation), T(2,5) against Wang. Several of these scans leave
    // d∘d components in the ideal of killed types; `close` checks that each closes to 0.
    let t = std::time::Instant::now();
    let r = wang_ref();
    let (mut d2, mut final_d2) = (0usize, 0usize);
    let knots: Vec<Vec<i32>> = vec![
        vec![1, 1, 1, 1, 1],
        vec![1, 1, 1, -1, 1],
        vec![1, -1, 1, 1, 1],
        vec![1, -1, -1, 1, -1],
        vec![1, 1, -1, -1, -1],
        vec![-1, -1, 1, -1, -1],
        vec![1, 1, 1, 1, -1, 1, 1],
        vec![-1, 1, 1, 1, -1, 1, 1, -1, 1],
    ];
    for w in &knots {
        let rev: Vec<i32> = w.iter().rev().copied().collect();
        for reduced in [false, true] {
            let un = khr_engine2_closure::<F31>(2, &gens(&rev), &[4], &[], 7, reduced, Closure::Exact).poincare.remove(&4).unwrap();
            assert_eq!(close::<F31>(&[3, 3], w, Some(4), 4, reduced), un, "via duality (expected): 3-colored {:?} at N=4 (reduced {}) vs uncolored gl_4 of the reversed word", w, reduced);
        }
        let sc = scan::<F31>(&[3, 3], w, Some(4));
        d2 += sc.stages.iter().map(|s| s.d2_ideal).sum::<usize>();
        final_d2 += sc.final_d2.len();
    }
    assert_eq!(close::<F31>(&[3, 3], &[1; 5], Some(4), 4, false), r[&(3, 4, 2)], "3-colored T(2,5) at N=4 vs wang.py");
    for w in [vec![1, -1, 1, 1, 1, -1], vec![1, 1, 1, 1, -1, 1, 1, -1], vec![1, 1, -1, 1], vec![-1, -1, 1, -1, -1, -1]] {
        let net: i32 = w.iter().sum();
        let pos = vec![net.signum(); net.unsigned_abs() as usize];
        for reduced in [false, true] {
            assert_eq!(close::<F31>(&[3, 3], &w, Some(4), 4, reduced), close::<F31>(&[3, 3], &pos, Some(4), 4, reduced), "3-colored {:?} = {:?} at N=4 (reduced {})", w, pos, reduced);
        }
        let sc = scan::<F31>(&[3, 3], &w, Some(4));
        d2 += sc.stages.iter().map(|s| s.d2_ideal).sum::<usize>();
        final_d2 += sc.final_d2.len();
    }
    assert_eq!(close::<F31>(&[3, 3], &[1, -1, 1, 1, 1, -1], Some(4), 4, false), hopf_formula(4, 3, 3));
    // the second prime on a word with final d∘d components
    let w = [1, 1, 1, 1, -1, 1, 1, -1];
    assert!(!scan::<F32b>(&[3, 3], &w, Some(4)).final_d2.is_empty());
    assert_eq!(close::<F32b>(&[3, 3], &w, Some(4), 4, false), close::<F31>(&[3, 3], &w, Some(4), 4, false));
    assert!(d2 > 0 && final_d2 > 0, "no d∘d component in the ideal was exercised");
    eprintln!("3-colored mixed words at N = 4: {} d∘d components in the ideal over all stages, {} on final complexes (closed to 0), {:.1}s", d2, final_d2, t.elapsed().as_secs_f64());
}

#[test]
#[ignore]
fn three_colored_hopf_and_trefoil_n5() {
    // Expected Wang rank 46 (trefoil) at N = 5. Not feasible on a laptop: the pruned scan still
    // needs Hom spaces between 3-rung ladders of fibre 60 and the rank-30 type (5,1) (one such
    // `hom_space` call took 5 401 s in the unpruned cloud run). Cloud only.
    let r = wang_ref();
    let h = three_colored_check::<F31>(&[1, 1], 5);
    assert_eq!(h, hopf_formula(5, 3, 3), "3-colored Hopf N=5 vs formula");
    let tr = three_colored_check::<F31>(&[1, 1, 1], 5);
    assert_eq!(tr, r[&(3, 5, 1)], "3-colored trefoil N=5 vs wang.py");
    assert_eq!(three_colored_check::<F32b>(&[1, 1], 5), h);
    assert_eq!(three_colored_check::<F32b>(&[1, 1, 1], 5), tr);
}

#[test]
fn pruned_scan_refuses_unchecked_d2() {
    // a pruned complex is a complex only modulo the killed ideal; without the d∘d components
    // the closure cannot be certified, so the scan must refuse
    let opts = XOpts { check_d2: false, verbose: false, prune_n: Some(3), ..Default::default() };
    let r = xscan::<F31>(&[2, 2], &gens(&[1, 1, 1]), &opts);
    assert!(r.is_err(), "pruned scan with check_d2 = false must be refused");
}
