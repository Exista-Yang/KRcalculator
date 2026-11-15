//! M4c B2: singlestep forms and the rex-chain planner (src/rexplan.rs), realized chains
//! (`leafgen::realize_chain`) checked against exact windows at small λ.
use krhom::combinat::dcoset::{reduced_for, top_labels};
use krhom::algebra::field::{Field, F31};
use krhom::library::homprov::hom_unknowns;
use krhom::library::leafgen::*;
use krhom::library::rexplan::*;
use krhom::combinat::schur::{Lp, Schur, Tid};
use krhom::colored::sexpr::{Expr, FlowIndex};
use krhom::colored::smap::*;
use krhom::colored::tmpl::*;
use krhom::combinat::web::Rung;
use std::collections::{BTreeMap, HashMap};

type F = F31;

#[test]
fn singlestep_forms() {
    // E0.1 on (2,1): split upright 1 (nothing new: r = whole), merge wall 2 -> [{2} +s2]
    assert_eq!(singlestep(&[2, 1], &[Rung::e(0, 1)]).unwrap(), vec![SStep::Up(2)]);
    // F0.1 on (2,1): split upright 0 as (1 | 1 moved), new wall 1, merge wall 2
    assert_eq!(singlestep(&[2, 1], &[Rung::f(0, 1)]).unwrap(), vec![SStep::Down(1), SStep::Up(2)]);
    // EKo (108) switchback on (1,2): E0.1 F0.1 = [{1} −s2 +s1 −s1? ...] has 4 steps or fewer
    let ss = singlestep(&[1, 2], &[Rung::e(0, 1), Rung::f(0, 1)]).unwrap();
    assert_eq!(ss.len(), 4, "{:?}", ss);
    // through a zero upright: merge then split = [+s −s]
    let ss = singlestep(&[1, 2], &[Rung::e(0, 2), Rung::f(0, 2)]).unwrap();
    assert_eq!(ss, vec![SStep::Up(1), SStep::Down(1)]);
    eprintln!("{}", singlestep_str(&[3, 3, 3], &[Rung::e(0, 3), Rung::f(1, 2)]));
}

#[test]
fn max_table_and_canonical_max() {
    let mut pl = Planner::new();
    for labels in [vec![2u8, 1], vec![1, 2], vec![2, 2], vec![1, 1, 1], vec![2, 0, 1], vec![2, 1, 2], vec![3, 0, 3], vec![3, 3]] {
        let k: u8 = labels.iter().sum();
        // every top composition of the same length
        let n = labels.len();
        let mut tops = vec![];
        fn comps(n: usize, k: u8, cur: &mut Vec<u8>, out: &mut Vec<Vec<u8>>) {
            if cur.len() == n - 1 {
                cur.push(k);
                out.push(cur.clone());
                cur.pop();
                return;
            }
            for a in 0..=k {
                cur.push(a);
                comps(n, k - a, cur, out);
                cur.pop();
            }
        }
        comps(n, k, &mut vec![], &mut tops);
        for top in tops {
            let cm = canonical_max(&labels, &top);
            assert!(pl.reduced_for_max(&labels, &cm), "canonical_max {:?} -> {:?} not reduced for the maximum", labels, top);
            let t = max_table(&labels, &top);
            assert_eq!(t.rows(), labels);
            assert_eq!(t.cols(), top);
        }
    }
}

/// All types q of (λ, μ) for μ a permutation of λ, with every reduced word of ≤ maxlen rungs.
fn reduced_words_by_type(lam: &[u8], maxlen: usize) -> (Schur, BTreeMap<Tid, Vec<Vec<Rung>>>) {
    let mut s = Schur::new(lam);
    let mut tops: Vec<Vec<u8>> = vec![];
    let mut p = lam.to_vec();
    p.sort();
    // all distinct permutations of lam
    fn perms(v: &mut Vec<u8>, k: usize, out: &mut Vec<Vec<u8>>) {
        if k == v.len() {
            if !out.contains(v) {
                out.push(v.clone());
            }
            return;
        }
        for i in k..v.len() {
            v.swap(k, i);
            perms(v, k + 1, out);
            v.swap(k, i);
        }
    }
    perms(&mut p, 0, &mut tops);
    let mut by: BTreeMap<Tid, Vec<Vec<Rung>>> = BTreeMap::new();
    for top in tops {
        for w in all_words(lam, &top, maxlen) {
            let ch = s.ladder_char(&w);
            let ex = s.expand(&ch);
            if ex.len() == 1 && ex[0].1 == Lp::one() {
                let q = ex[0].0;
                if reduced_for(&mut s, &w, q).ok() {
                    by.entry(q).or_default().push(w);
                }
            }
        }
    }
    (s, by)
}

/// Planner connectivity at small λ: from the canonical ladder of q to every reduced word of q
/// (≤ 4 rungs), every visited word reduced, every move well-formed.
#[test]
fn planner_connects_small() {
    for lam in [vec![1u8, 1, 1], vec![1, 2, 1], vec![2, 1, 2], vec![2, 2, 1], vec![2, 2, 2]] {
        let (mut s, by) = reduced_words_by_type(&lam, 4);
        let mut pl = Planner::new();
        let opts = PlanOpts { max_len: 8, exact_budget: 100_000, ..PlanOpts::default() };
        let (mut ok, mut fail, mut maxlen) = (0, 0, 0);
        let mut kinds: BTreeMap<String, usize> = BTreeMap::new();
        for (&q, ws) in &by {
            let lq = krhom::library::sindec::canonical_ladder(&mut s, q);
            for w in ws {
                match pl.plan(&mut s, q, w, &lq, &opts) {
                    Ok(ch) => {
                        ok += 1;
                        maxlen = maxlen.max(ch.len());
                        for (k, mv) in ch.moves.iter().enumerate() {
                            *kinds.entry(format!("{:?}", mv.kind)).or_insert(0) += 1;
                            let win = window_of(&lam, &ch.words[k], mv.at, mv.at + mv.src.len());
                            assert_eq!((win.labels.clone(), win.sub.clone(), win.first), (mv.labels.clone(), mv.src.clone(), mv.first), "move window mismatch");
                            assert!(reduced_for(&mut s, &ch.words[k + 1], q).ok());
                        }
                        assert_eq!(ch.words.last().unwrap(), &lq);
                    }
                    Err(e) => {
                        fail += 1;
                        eprintln!("{:?} q {}: [{}] -> [{}]: {}", lam, s.table(q).fmt(), krhom::library::smapio::rungs_str(w), krhom::library::smapio::rungs_str(&lq), e);
                    }
                }
            }
        }
        eprintln!("planner {:?}: {} types, {} chains found, {} failed, max length {}, moves {:?}, stats {:?}", lam, by.len(), ok, fail, maxlen, kinds, pl.stats);
        assert_eq!(fail, 0, "planner failed at {:?}", lam);
    }
}

/// Realized chains at small λ: genuine, degree 0, 1⊗ ↦ 1⊗ at every coset, scalar 1 at the top
/// cosets of q, and in the exact Hom⁰ window where affordable.
#[test]
fn realized_chains_small() {
    let mut checked = 0;
    for lam in [vec![1u8, 1, 1], vec![1, 2, 1], vec![2, 1, 2], vec![2, 2, 1]] {
        let (mut s, by) = reduced_words_by_type(&lam, 4);
        let mut pl = Planner::new();
        let opts = PlanOpts { max_len: 8, exact_budget: 100_000, ..PlanOpts::default() };
        let mut ar: Arena<F> = Arena::new(&lam, TmplLib::new(None));
        let mut cache = HashMap::new();
        let k: usize = lam.iter().map(|&a| a as usize).sum();
        for (&q, ws) in &by {
            let lq = krhom::library::sindec::canonical_ladder(&mut s, q);
            for w in ws.iter().take(6) {
                if *w == lq {
                    continue;
                }
                let ch = pl.plan(&mut s, q, w, &lq, &opts).unwrap();
                let (m, rep) = realize_chain(&mut ar, &ch, &mut cache, 70).unwrap_or_else(|e| panic!("{:?} [{}]: {}", lam, krhom::library::smapio::rungs_str(w), e));
                ar.assert_genuine(m).unwrap();
                assert_eq!(ar.info(m).deg, 0);
                let p = points::<F>(k, 4, 1, 0).pop().unwrap();
                let (ea, eb) = (ar.info(m).src, ar.info(m).tgt);
                let mut ev = Eval::new(&ar, p.clone());
                let lay = ev.layout(ea, eb);
                for &c in &lay.cosets {
                    let bl = ev.block(m, c);
                    let ones = Block { rows: bl.cols, cols: 1, a: vec![F::one(); bl.cols as usize] };
                    assert!(bl.mul(&ones).a.iter().all(|v| v.is_one()), "1⊗ not mapped to 1⊗");
                }
                let cq = krhom::library::leaves::rep_rank_of(&s, &lam, q);
                let bq = ev.block(m, cq);
                assert_eq!((bq.rows, bq.cols), (1, 1));
                assert!(bq.a[0].is_one());
                // exact window
                let u = hom_unknowns(&lam, w, &lq, 0);
                if u <= 20_000 {
                    let (wa, wb) = (Window::new(&lam, w), Window::new(&lam, &lq));
                    let ids = ar.lib().try_hom_basis_cap(&wa, &wb, 0, u128::MAX).unwrap();
                    let ys: Vec<Vec<F>> = (0..2).map(|i| points::<F>(k, 40 + i, 1, 0)[0].x.clone()).collect();
                    let flat = |wl: &WinLoc<F>| -> Vec<F> { wl.dense().into_iter().flatten().collect() };
                    let mut rows: Vec<Vec<F>> = ids.iter().map(|&t| ys.iter().flat_map(|y| flat(&ar.lib_ref().localize(t, y))).collect()).collect();
                    let r0 = krhom::algebra::dense::rank(&rows, rows.first().map_or(0, |r| r.len()));
                    // the chain localized on the whole word, in sorted flow order
                    let (fa, fb) = (FlowIndex::build(&Expr::new(&lam, w).unwrap()), FlowIndex::build(&Expr::new(&lam, &lq).unwrap()));
                    let mut v = vec![];
                    for y in &ys {
                        let mut d = vec![vec![F::zero(); fa.nflows()]; fb.nflows()];
                        let mut ev2 = Eval::new(&ar, Point { x: y.clone(), id: 0 });
                        for &c in &lay.cosets {
                            let bl = ev2.block(m, c);
                            for (ii, i) in fb.range(c).enumerate() {
                                for (jj, j) in fa.range(c).enumerate() {
                                    d[i as usize][j as usize] = bl.get(ii, jj);
                                }
                            }
                        }
                        v.extend(d.into_iter().flatten());
                    }
                    let n = v.len();
                    rows.push(v);
                    assert_eq!(krhom::algebra::dense::rank(&rows, n), r0, "{:?} chain [{}] -> [{}] not in the exact Hom^0 window", lam, krhom::library::smapio::rungs_str(w), krhom::library::smapio::rungs_str(&lq));
                }
                if checked < 3 {
                    eprintln!("{:?} [{}] -> [{}]: {:?}", lam, krhom::library::smapio::rungs_str(w), krhom::library::smapio::rungs_str(&lq), rep.lines);
                }
                checked += 1;
            }
        }
    }
    eprintln!("realized chains checked: {}", checked);
    assert!(checked > 10);
}

/// top_labels sanity for the window helpers
#[test]
fn window_splice_roundtrip() {
    let lam = [2u8, 1, 2];
    let w = vec![Rung::e(0, 1), Rung::e(1, 2), Rung::f(1, 1), Rung::f(0, 2)];
    for at in 0..w.len() {
        for end in at + 1..=w.len() {
            let win = window_of(&lam, &w, at, end);
            assert_eq!(splice(&w, &win, &win.sub), w);
            assert_eq!(top_labels(&win.labels, &win.sub).len(), win.labels.len());
        }
    }
}

/// B4 machinery at small λ: every rex-class pair target (X reduced for q, X ≠ L_q) of the G1
/// plan of (2,1,2) and (2,2,2): π, ι from rex chains, `certify_pairs`, π∘ι = id.
#[test]
fn rex_class_pairs_small() {
    let mut n = 0;
    for lam in [vec![2u8, 1, 2], vec![2, 2, 2]] {
        let mut s = Schur::new(&lam);
        let plan = krhom::library::sindec::pruned_plan(&mut s, None);
        for pp in &plan.pairs {
            let mut x = krhom::library::sindec::canonical_ladder(&mut s, pp.p);
            x.push(pp.rung);
            for tg in &pp.kept {
                let lq = krhom::library::sindec::canonical_ladder(&mut s, tg.q);
                if lq == x || !reduced_for(&mut s, &x, tg.q).ok() {
                    continue;
                }
                let mut ar: Arena<F> = Arena::new(&lam, TmplLib::new(None));
                let (xe, le) = (ar.expr(Expr::new(&lam, &x).unwrap()), ar.expr(Expr::new(&lam, &lq).unwrap()));
                let opts = PlanOpts { exact_budget: 100_000, max_len: 10, ..PlanOpts::default() };
                let (p, r) = krhom::library::leaves::rex_class_pairs(&mut ar, &mut s, xe, tg.q, le, &opts, 99).unwrap_or_else(|e| panic!("{:?} [{}] -> q {}: {}", lam, krhom::library::smapio::rungs_str(&x), s.table(tg.q).fmt(), e));
                assert_eq!(p.form, vec![vec![F::one()]]);
                assert!(r.theorem_assoc.is_empty());
                n += 1;
            }
        }
    }
    eprintln!("rex-class pair targets certified: {}", n);
    assert!(n > 0);
}
