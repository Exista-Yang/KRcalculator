//! M4c A1, Stage 0 (src/library/leaves.rs): the 5 product closure types of (3,3,3) at N = 6 through
//! `ProductClosure` (product rule first, inner = the production chain), certified by
//! `certify_closure` inside the provider and again here; non-product types are refused by the
//! product rule.
use krhom::algebra::field::F31;
use krhom::library::homprov::*;
use krhom::library::leaves::*;
use krhom::combinat::schur::{Schur, Table};
use krhom::colored::sexpr::Expr;
use krhom::library::sindec::{canonical_ladder, provider_for};
use krhom::colored::smap::Arena;
use krhom::colored::tmpl::TmplLib;

const PRODUCTS: [&str; 5] = ["0,0,3/3,0,0/0,3,0", "0,3,0/0,0,3/3,0,0", "1,0,2/2,0,1/0,3,0", "1,2,0/0,0,3/2,1,0", "0,3,0/1,0,2/2,0,1"];

#[test]
fn stage0_factorizations() {
    let lam = [3u8, 3, 3];
    let mut s = Schur::new(&lam);
    for ts in PRODUCTS {
        let p = s.intern(Table::parse(ts));
        let w = canonical_ladder(&mut s, p);
        let f = product_factorization(&mut s, &w, p, 256).unwrap_or_else(|| panic!("{} is not found as a product", ts));
        assert!(f.canonical, "{}: factors are not canonical ladders", ts);
        assert!(f.swaps.len() <= 1, "{}: {} swaps", ts, f.swaps.len());
        assert_eq!(s.c_diag(f.a).mul(&s.c_diag(f.b)), s.c_diag(p));
    }
    // the 7 non-product failing types (L3 §3.2) are not products
    for ts in ["1,1,1/1,1,1/1,1,1", "1,0,2/1,2,0/1,1,1", "1,1,1/0,2,1/2,0,1", "0,2,1/2,0,1/1,1,1", "1,1,1/1,0,2/1,2,0", "0,1,2/2,0,1/1,2,0", "0,2,1/1,0,2/2,1,0"] {
        let p = s.intern(Table::parse(ts));
        let w = canonical_ladder(&mut s, p);
        assert!(product_factorization(&mut s, &w, p, 256).is_none(), "{} unexpectedly a product", ts);
    }
}

#[test]
fn stage0_product_closure_certified() {
    let lam = [3u8, 3, 3];
    let mut s = Schur::new(&lam);
    let mut ar = Arena::<F31>::new(&lam, TmplLib::new(None));
    let inner = provider_for::<F31>(&lam, Some(std::path::Path::new("tests/data/homprov/costs_333_m5c.txt")), 1.0);
    let mut pc = ProductClosure::new(inner);
    pc.product_first = true;
    let b = Budget::default();
    let r0 = ar.expr(Expr::new(&lam, &[]).unwrap());
    for ts in PRODUCTS {
        let p = s.intern(Table::parse(ts));
        let w = canonical_ladder(&mut s, p);
        let lp = ar.expr(Expr::new(&lam, &w).unwrap());
        let want = cdiag_degrees(&mut s, p);
        for outs in [false, true] {
            let (ms, _) = HomProvider::closure(&mut pc, &mut ar, &mut s, lp, p, None, outs, &b).unwrap_or_else(|e| panic!("stage0 {} outs={}: {:?}", ts, outs, e));
            certify_closure(&ar, &lam, r0, lp, None, &ms, &want, outs, 0xa1).unwrap();
        }
    }
    // the product rule refuses a non-product type
    let p = s.intern(Table::parse("1,1,1/1,1,1/1,1,1"));
    let w = canonical_ladder(&mut s, p);
    let lp = ar.expr(Expr::new(&lam, &w).unwrap());
    let e = pc.product_closure(&mut ar, &mut s, lp, p, true, &b).expect_err("non-product accepted");
    assert_eq!(e[0].kind, MissKind::Refused);
}

// =============================================================================================
// M4c A2: Lemma A leaves (Leaves backend alone)
// =============================================================================================
use krhom::combinat::dcoset::{top_labels, valid_word};
use krhom::algebra::field::{Field, F32b};
use krhom::combinat::schur::{flows, Tid};
use krhom::colored::smap::{points, MapNode};
use krhom::combinat::web::{Dir, Rung};
use krhom::colored::tmpl::Window;

/// every valid word of length ≤ len on lam (rungs of any size)
fn words(lam: &[u8], len: usize) -> Vec<Vec<Rung>> {
    let mut out = vec![vec![]];
    let mut frontier = vec![vec![]];
    for _ in 0..len {
        let mut nf = vec![];
        for w in &frontier {
            let top = top_labels(lam, w);
            for i in 0..lam.len() - 1 {
                for dir in [Dir::E, Dir::F] {
                    let giver = if dir == Dir::E { top[i + 1] } else { top[i] } as usize;
                    for r in 1..=giver {
                        let mut w2: Vec<Rung> = w.clone();
                        w2.push(Rung { i, dir, r });
                        if valid_word(lam, &w2) {
                            nf.push(w2);
                        }
                    }
                }
            }
        }
        out.extend(nf.iter().cloned());
        frontier = nf;
    }
    out
}

fn lp_degrees(p: &krhom::combinat::schur::Lp) -> Vec<i32> {
    let mut v = vec![];
    for (i, &c) in p.c.iter().enumerate() {
        for _ in 0..c {
            v.push(p.lo + i as i32);
        }
    }
    v
}

/// Schur-polynomial basis of R^λ in degree 2k: one partition per upright (≤ λ_u parts)
fn rlam_basis(lam: &[u8], k: usize) -> Vec<Vec<Vec<u8>>> {
    fn parts(k: usize, maxp: usize, rows: usize) -> Vec<Vec<u8>> {
        if k == 0 {
            return vec![vec![]];
        }
        if rows == 0 {
            return vec![];
        }
        let mut out = vec![];
        for f in (1..=k.min(maxp)).rev() {
            for mut rest in parts(k - f, f, rows - 1) {
                let mut v = vec![f as u8];
                v.append(&mut rest);
                out.push(v);
            }
        }
        out
    }
    fn rec(lam: &[u8], u: usize, k: usize, cur: &mut Vec<Vec<u8>>, out: &mut Vec<Vec<Vec<u8>>>) {
        if u == lam.len() {
            if k == 0 {
                out.push(cur.clone());
            }
            return;
        }
        for kk in 0..=k {
            for p in parts(kk, kk, lam[u] as usize) {
                cur.push(p);
                rec(lam, u + 1, k - kk, cur, out);
                cur.pop();
            }
        }
    }
    let mut out = vec![];
    rec(lam, 0, k, &mut vec![], &mut out);
    out
}

/// V0: for every word X with top λ (length ≤ len), the Lemma A leaves X -> R (Out) and R -> X
/// (In) to the identity type: n = c_id(char X)(1) maps with degree multiset c_id(char X); their
/// R^λ-span in each degree d equals the exact `sbim::hom_space` basis of Hom^d (stacked rank of
/// the leaves·monomials = dim Hom^d = rank of their union with the exact basis).
fn v0<F: Field>(lam: &[u8], len: usize, dspan: i32, max_unk: u128) -> usize {
    let mut s = Schur::new(lam);
    let mut ar = Arena::<F>::new(lam, TmplLib::new(None));
    let mut lv = Leaves::new();
    lv.use_elem = true;
    let id = s.intern(Table::identity(lam));
    let r0 = ar.expr(Expr::new(lam, &[]).unwrap());
    let c = identity_coset(lam);
    let mut checked = 0;
    for w in words(lam, len) {
        if top_labels(lam, &w) != lam {
            continue;
        }
        let ch = s.ladder_char(&w);
        let cid = ch.coef(id);
        if cid.is_zero() {
            continue;
        }
        let want = lp_degrees(&cid);
        let dmin = *want.iter().min().unwrap();
        let x = ar.expr(Expr::new(lam, &w).unwrap());
        for dir in [LeafDir::Out, LeafDir::In] {
            if (dmin..=dmin + dspan).any(|d| hom_unknowns(lam, &w, &[], d) > max_unk) {
                continue;
            }
            let mut st = LeafStats::default();
            let fams = lv.build_leaves(&mut ar, &mut s, x, &[id], dir, &mut st).unwrap_or_else(|e| panic!("V0 {:?} [{:?}] {:?}: {:?}", lam, w, dir, e));
            let f = fams.iter().find(|f| f.t == id).unwrap();
            let mut degs: Vec<i32> = f.maps.iter().map(|&m| ar.info(m).deg).collect();
            degs.sort_unstable();
            assert_eq!(degs, want, "V0 {:?} {:?} {:?}: leaf degrees", lam, w, dir);
            let (a, b) = if dir == LeafDir::Out { (x, r0) } else { (r0, x) };
            for &m in &f.maps {
                assert_eq!((ar.info(m).src, ar.info(m).tgt), (a, b));
                assert_genuine(&ar, m).unwrap();
            }
            for d in dmin..=dmin + dspan {
                // leaves times R^λ monomials of degree d - deg
                let mut span = vec![];
                for &m in &f.maps {
                    let dd = d - ar.info(m).deg;
                    if dd < 0 || dd % 2 != 0 {
                        continue;
                    }
                    for nus in rlam_basis(lam, (dd / 2) as usize) {
                        let mut g = m;
                        for (u, nu) in nus.iter().enumerate() {
                            if nu.is_empty() {
                                continue;
                            }
                            let dot = ar.add(MapNode::Dot { e: r0, level: 0, upright: u, nu: nu.clone() }, "v0 dot").unwrap();
                            g = ar.add(if dir == LeafDir::Out { MapNode::Comp(dot, g) } else { MapNode::Comp(g, dot) }, "v0").unwrap();
                        }
                        span.push(g);
                    }
                }
                let (la, lb) = if dir == LeafDir::Out { (w.clone(), vec![]) } else { (vec![], w.clone()) };
                let ids = ar.lib().try_hom_basis_cap(&Window::new(lam, &la), &Window::new(lam, &lb), d, u128::MAX).unwrap();
                let mut exact = vec![];
                for t in ids {
                    exact.push(ar.add(MapNode::Local { src: a, tgt: b, at: 0, first: 0, tmpl: t }, "v0 exact").unwrap());
                }
                let dim = exact.len();
                let pts = points::<F>(ar.k(), 0x0e0 + d as u64, dim.max(2) + 2, 0);
                let r_span = stacked_rank(&ar, &span, &pts, &[c]);
                let mut all = span.clone();
                all.extend(exact.iter().copied());
                let r_all = stacked_rank(&ar, &all, &pts, &[c]);
                let r_ex = stacked_rank(&ar, &exact, &pts, &[c]);
                assert_eq!(r_ex, dim, "V0 {:?} {:?}: exact basis not detected at {} points", lam, w, pts.len());
                assert_eq!((r_span, r_all), (dim, dim), "V0 {:?} [{:?}] {:?} degree {}: leaves span {} / union {} / Hom^d {}", lam, w, dir, d, r_span, r_all, dim);
            }
            checked += 1;
        }
    }
    checked
}

#[test]
fn v0_leaves_vs_hom_space() {
    // λ = (2) as (2,0): one root block and an empty upright to split into
    let mut tot = 0;
    for (lam, len, dspan) in [(vec![2u8, 0], 4, 4), (vec![1, 1], 4, 4), (vec![2, 1], 3, 4), (vec![1, 2, 1], 3, 2)] {
        let n = v0::<F31>(&lam, len, dspan, 30_000);
        eprintln!("[V0] {:?}: {} (word, direction) families, leaves·R^λ = sbim::hom_space in {} degrees each", lam, n, dspan + 1);
        assert!(n >= 4, "V0 {:?}: only {} cases", lam, n);
        tot += n;
    }
    eprintln!("[V0] total {} cases", tot);
}

/// SLL Ex 7.39: λ = (2) (as (2,0)), [s,∅,s] = F0.1 E0.1 -> [s] = ∅: two leaves of degrees ±1
#[test]
fn sll_ex_7_39() {
    let lam = [2u8, 0];
    let mut s = Schur::new(&lam);
    let mut ar = Arena::<F31>::new(&lam, TmplLib::new(None));
    let w = krhom::library::smapio::parse_rungs("F0.1,E0.1").unwrap();
    let x = ar.expr(Expr::new(&lam, &w).unwrap());
    let r0 = ar.expr(Expr::new(&lam, &[]).unwrap());
    let id = s.intern(Table::identity(&lam));
    let lv = Leaves::new();
    for dir in [LeafDir::Out, LeafDir::In] {
        let mut st = LeafStats::default();
        let fams = lv.build_leaves(&mut ar, &mut s, x, &[id], dir, &mut st).unwrap();
        let ms = &fams[0].maps;
        let mut degs: Vec<i32> = ms.iter().map(|&m| ar.info(m).deg).collect();
        degs.sort_unstable();
        assert_eq!(degs, vec![-1, 1], "Ex 7.39 {:?}", dir);
        certify_closure(&ar, &lam, r0, x, None, ms, &[-1, 1], dir == LeafDir::Out, 0x739).unwrap();
    }
}

/// SLL Ex 7.40/7.41 (S_3, I = (2,1), J = (1,2)): light leaves of every word (2,1) -> (1,2) to
/// every (I, J) coset; each family has c_t(char X)(1) maps, degree multiset c_t(char X), an
/// invertible block at rep(t) at two points, and is genuine. The paper's exact expressions are
/// not in the repository: the test covers all words of these I, J up to 5 rungs (and with an
/// empty third upright, where rungs into it are single EKLP steps), among them words with 4 and 5
/// light leaves in all (the sizes of Ex 7.40 and 7.41) and single-coset families of 4.
fn ex740(lam: &[u8], top: &[u8], len: usize) -> (std::collections::BTreeMap<usize, usize>, std::collections::BTreeSet<usize>) {
    let mut s = Schur::new(lam);
    let mut ar = Arena::<F31>::new(lam, TmplLib::new(None));
    let lv = Leaves::new();
    let mut sizes = std::collections::BTreeMap::new();
    let mut totals = std::collections::BTreeSet::new();
    for w in words(lam, len) {
        if top_labels(lam, &w) != top {
            continue;
        }
        let ch = s.ladder_char(&w);
        totals.insert(ch.t.iter().map(|(_, c)| c.at1() as usize).sum::<usize>());
        let x = ar.expr(Expr::new(lam, &w).unwrap());
        let ts: Vec<(Tid, krhom::combinat::schur::Lp)> = ch.t.clone();
        for (t, c) in ts {
            for dir in [LeafDir::Out, LeafDir::In] {
                let mut st = LeafStats::default();
                let fams = lv.build_leaves(&mut ar, &mut s, x, &[t], dir, &mut st).unwrap_or_else(|e| panic!("7.40 [{:?}] -> {}: {:?}", w, s.table(t).fmt(), e));
                let f = fams.iter().find(|f| f.t == t).unwrap();
                let mut degs: Vec<i32> = f.maps.iter().map(|&m| ar.info(m).deg).collect();
                degs.sort_unstable();
                assert_eq!(degs, lp_degrees(&c), "7.40 [{:?}] -> {} degrees", w, s.table(t).fmt());
                let rc = rep_rank_of(&s, lam, t);
                for pt in points::<F31>(ar.k(), 0x740, 2, 0) {
                    let mut ev = krhom::colored::smap::Eval::new(&ar, pt);
                    let rows: Vec<Vec<F31>> = f.maps.iter().map(|&m| ev.block(m, rc).a.clone()).collect();
                    assert_eq!(krhom::algebra::dense::rank(&rows, f.maps.len()), f.maps.len());
                }
                for &m in &f.maps {
                    assert_genuine(&ar, m).unwrap();
                }
                *sizes.entry(f.maps.len()).or_insert(0usize) += 1;
            }
        }
    }
    (sizes, totals)
}

#[test]
fn sll_ex_7_40_7_41() {
    let (a, ta) = ex740(&[2, 1], &[1, 2], 5);
    eprintln!("[7.40/7.41] (2,1) -> (1,2), words ≤ 5 rungs: family sizes (size: count) {:?}; leaves per word (all cosets) {:?}", a, ta);
    // with an empty third upright, rungs into it are single EKLP steps (splits)
    let (b, tb) = ex740(&[2, 1, 0], &[1, 2, 0], 4);
    eprintln!("[7.40/7.41] (2,1,0) -> (1,2,0), words ≤ 4 rungs: family sizes (size: count) {:?}; leaves per word {:?}", b, tb);
    assert!(a.contains_key(&4), "no family of 4 leaves");
    assert!(ta.contains(&4) && ta.contains(&5), "no word with 4 and 5 leaves in all");
}

/// G1(a) closure part with the Leaves backend alone: every type of (λ, λ), ins and outs
fn g1a_closure<F: Field>(lam: &[u8]) -> (usize, usize, std::collections::BTreeMap<String, usize>) {
    let mut s = Schur::new(lam);
    let mut ar = Arena::<F>::new(lam, TmplLib::new(None));
    let (_, order) = g1_targets(&mut s);
    let mut lv = Leaves::new();
    lv.use_elem = true;
    let r0 = ar.expr(Expr::new(lam, &[]).unwrap());
    let b = Budget::default();
    let (mut ok, mut tot) = (0, 0);
    let mut kinds = std::collections::BTreeMap::new();
    for p in order {
        let w = canonical_ladder(&mut s, p);
        let lp = ar.expr(Expr::new(lam, &w).unwrap());
        let want = cdiag_degrees(&mut s, p);
        for outs in [false, true] {
            tot += 1;
            match HomProvider::closure(&mut lv, &mut ar, &mut s, lp, p, None, outs, &b) {
                Ok((ms, _)) => {
                    certify_closure(&ar, lam, r0, lp, None, &ms, &want, outs, 0x9a1).unwrap();
                    ok += 1;
                    for r in &lv.last.steps {
                        *kinds.entry(format!("{:?}", r.kind)).or_insert(0usize) += 1;
                    }
                }
                Err(e) => eprintln!("[G1a closure] {:?} {} outs={}: {:?}", lam, s.table(p).fmt(), outs, e),
            }
        }
    }
    (ok, tot, kinds)
}

#[test]
fn g1a_closure_leaves_alone() {
    for lam in [vec![3u8, 3], vec![2, 2, 2]] {
        let t0 = std::time::Instant::now();
        let (a, n, ka) = g1a_closure::<F31>(&lam);
        let (b, m, _) = g1a_closure::<F32b>(&lam);
        eprintln!("[G1a closure, Leaves alone] {:?}: F31 {}/{}, F32b {}/{} ({:.1}s); step families (F31) {:?}", lam, a, n, b, m, t0.elapsed().as_secs_f64(), ka);
        let want = if lam.len() == 2 { 8 } else { 42 };
        assert_eq!((a, n, b, m), (want, want, want, want));
    }
}

#[allow(dead_code)]
fn _unused(_: u128) -> u128 {
    flows(&[1], &[])
}

/// validate_elem: every E2/E3 closed form with labels ≤ 3 is a nonzero multiple of the unique
/// exact map (both primes); a wrong window is rejected
#[test]
fn validate_elem_small() {
    use krhom::colored::smap::ElemKind::*;
    fn run<F: Field>() -> usize {
        let mut n = 0;
        for la in 0..=3u8 {
            for lb in 0..=3u8 {
                for k in 1..=3usize {
                    for (a, dirs) in [(vec![Rung::e(0, k), Rung::f(0, k)], 0), (vec![Rung::f(0, k), Rung::e(0, k)], 1)] {
                        let _ = dirs;
                        if !valid_word(&[la, lb], &a) {
                            continue;
                        }
                        for dir in [LeafDir::Out, LeafDir::In] {
                            if let Some((kind, _)) = elem_pattern(&[la, lb], &a, &[], dir) {
                                validate_elem::<F>(kind, &[la, lb], &a, 120_000).unwrap_or_else(|e| panic!("{}", e));
                                n += 1;
                            }
                        }
                    }
                }
            }
        }
        // a closed form on a window that is not its pattern is refused by the arena or rejected
        assert!(validate_elem::<F>(Eps, &[2, 2], &[Rung::e(0, 1), Rung::f(0, 1)], 120_000).is_err());
        n
    }
    let (a, b) = (run::<F31>(), run::<F32b>());
    eprintln!("[validate_elem] {} / {} closed forms validated (F31 / F32b)", a, b);
    assert!(a >= 30 && a == b);
}

/// A3 wiring: the production chain with Stage 0 and Leaves for closure (both orders) certifies
/// every closure type of (2,2,2); non-closure calls go to the inner chain
#[test]
fn m4c_chain_closure_222() {
    let lam = [2u8, 2, 2];
    for lf in [true, false] {
        let mut s = Schur::new(&lam);
        let mut ar = Arena::<F31>::new(&lam, TmplLib::new(None));
        let (mut prov, stamp) = provider_for_stamped_m4c::<F31>(&lam, None, 1.0, lf);
        assert!(stamp.contains(LEAVES_STAMP) && stamp.contains(STAGE0_STAMP_PF));
        let r0 = ar.expr(Expr::new(&lam, &[]).unwrap());
        let order = s.canonical_basis(&lam).order.clone();
        for p in order {
            let w = canonical_ladder(&mut s, p);
            let lp = ar.expr(Expr::new(&lam, &w).unwrap());
            let want = cdiag_degrees(&mut s, p);
            for outs in [false, true] {
                let (ms, _) = prov.closure(&mut ar, &mut s, lp, p, None, outs, &Budget::default()).unwrap_or_else(|e| panic!("m4c chain {} lf={}: {:?}", s.table(p).fmt(), lf, e));
                certify_closure(&ar, &lam, r0, lp, None, &ms, &want, outs, 0x3a).unwrap();
            }
        }
        eprintln!("[m4c chain] (2,2,2) closure, leaves_first = {}: all types certified ({})", lf, stamp);
    }
}
