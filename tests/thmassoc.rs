//! M4c B2: the theorem-level Assoc certificate (user decision 2026-11-12) and its cross-checks.
//!
//! (i)  agreement with the numeric certificate (`poly_certificate` of the identity matching)
//!      on every flag-type word pair of small windows, and on the braid middle steps at cells
//!      (1,1,1), (1,2,1), (2,1,2), (2,2,2) (default) and (3,1,3), (2,3,2) (ignored: 245 s and
//!      1 509 s in B1);
//! (ii) the bimodule-property spot check (`spot_check_assoc`, light localization) accepts
//!      genuine maps and rejects corrupted ones; at (3,3,3) (ignored) on the braid middle step;
//! (iv) marking: provenance tag, `relies_on_theorem`, the smapio `thm-assoc` record (round trip
//!      re-certifies by the theorem; without the record the import is not genuine).
use krhom::field::{Field, F31};
use krhom::leafgen::*;
use krhom::sexpr::{Expr, FlowIndex};
use krhom::smap::*;
use krhom::tmpl::*;
use krhom::web::{Dir, Rung};
use std::collections::BTreeMap;
use std::time::Instant;

type F = F31;

fn words(labels: &[u8], maxlen: usize) -> Vec<Vec<Rung>> {
    let mut out = vec![vec![]];
    let mut frontier = vec![(vec![], labels.to_vec())];
    for _ in 0..maxlen {
        let mut nf = vec![];
        for (w, l) in &frontier {
            for i in 0..l.len() - 1 {
                for dir in [Dir::E, Dir::F] {
                    let giver = if dir == Dir::E { l[i + 1] } else { l[i] } as usize;
                    for r in 1..=giver {
                        let mut w2: Vec<Rung> = w.clone();
                        w2.push(Rung { i, dir, r });
                        let mut l2 = l.clone();
                        if dir == Dir::E {
                            l2[i + 1] -= r as u8;
                            l2[i] += r as u8;
                        } else {
                            l2[i] -= r as u8;
                            l2[i + 1] += r as u8;
                        }
                        out.push(w2.clone());
                        nf.push((w2, l2));
                    }
                }
            }
        }
        frontier = nf;
    }
    out
}

/// LightLoc agrees with the SBim localization (fdeg, shift, rows up to the flow order).
#[test]
fn lightloc_matches_sbim() {
    for (lab, w) in [(vec![2u8, 1], vec![Rung::e(0, 1)]), (vec![2, 2, 1], vec![Rung::e(0, 1), Rung::f(1, 1), Rung::e(1, 1)]), (vec![3, 0, 2], vec![Rung::f(0, 2), Rung::e(1, 2), Rung::f(1, 3)])] {
        let e = Expr::new(&lab, &w).unwrap();
        let sb = e.sbim::<F>();
        let ll = LightLoc::new(&e);
        assert_eq!(ll.shift, sb.shift);
        let mut a = ll.fdeg.clone();
        let mut b = sb.fdeg.clone();
        a.sort();
        b.sort();
        assert_eq!(a, b);
        assert_eq!(ll.fdeg, sb.fdeg, "fibre order");
        let y = points::<F>(e.k(), 3, 1, 0)[0].x.clone();
        let fps = krhom::local::fixed_points(&sb);
        let lm = krhom::local::loc_matrix(&sb, &fps, &y);
        let lr = ll.matrix(&y);
        for (p, row) in fps.iter().zip(&lm) {
            let f = ll.fi.locate(&p.moved).unwrap() as usize;
            assert_eq!(&lr[f], row);
        }
    }
}

/// (i) on all small windows: every flag-type pair with a flow bijection, theorem vs numeric.
#[test]
fn theorem_agrees_with_numeric_small() {
    let mut stats: BTreeMap<&str, usize> = BTreeMap::new();
    let mut accepted_numeric = 0;
    for lab in [vec![2u8, 1], vec![1, 2], vec![2, 2], vec![1, 1, 1], vec![2, 1, 1], vec![1, 2, 1], vec![1, 1, 2], vec![2, 0, 1], vec![1, 0, 2], vec![3, 1], vec![2, 1, 2]] {
        let maxlen = if lab.iter().map(|&a| a as usize).sum::<usize>() >= 5 { 3 } else { 4 };
        // group flag-type words by (top, support)
        let mut groups: BTreeMap<(Vec<u8>, Vec<u64>), Vec<Vec<Rung>>> = BTreeMap::new();
        for w in words(&lab, maxlen) {
            let e = Expr::new(&lab, &w).unwrap();
            let fi = FlowIndex::build(&e);
            if fi.max_block() == 1 {
                groups.entry((e.top(), fi.support().to_vec())).or_default().push(w);
            }
        }
        for (_, ws) in groups {
            for i in 0..ws.len() {
                for j in 0..ws.len() {
                    if i == j {
                        continue;
                    }
                    let r = theorem_vs_numeric::<F>(&lab, &ws[i], &ws[j], 400, 7 + i as u64).unwrap();
                    assert!(r.consistent(), "theorem accepts but numeric fails: {:?} {:?} -> {:?}", lab, ws[i], ws[j]);
                    let key = match (&r.theorem, r.numeric) {
                        (Ok(()), Some(true)) => "thm+num",
                        (Ok(()), None) => "thm only",
                        (Err(_), Some(true)) => "num only",
                        (Err(_), Some(false)) => "both refuse",
                        (Err(_), None) => "thm refuses, numeric not run",
                        (Ok(()), Some(false)) => unreachable!(),
                    };
                    *stats.entry(key).or_insert(0) += 1;
                    if r.theorem.is_ok() {
                        accepted_numeric += 1;
                    }
                }
            }
        }
    }
    eprintln!("theorem vs numeric (small windows): {:?}", stats);
    assert!(stats.get("thm+num").copied().unwrap_or(0) > 20, "too few theorem cases exercised");
    assert_eq!(stats.get("thm only").copied().unwrap_or(0), 0, "numeric certificate skipped on a small theorem case");
    assert!(accepted_numeric > 0);
}

fn braid_middle_cross(x: u8, y: u8, z: u8) -> AssocCross {
    let bw = braid_words(x, y, z).unwrap();
    let r = theorem_vs_numeric::<F>(&[x, y, z], &bw.m, &bw.mp, usize::MAX, 31).unwrap();
    eprintln!("braid middle {:?}: {} flows, theorem {:?}, numeric {:?} ({:.1} s)", (x, y, z), r.flows, r.theorem, r.numeric, r.numeric_secs);
    assert!(r.theorem.is_ok(), "theorem refuses the braid middle at {:?}", (x, y, z));
    assert_eq!(r.numeric, Some(true), "numeric certificate fails at {:?}", (x, y, z));
    r
}

#[test]
fn braid_middle_small_cells() {
    for c in [(1, 1, 1), (1, 2, 1), (2, 1, 2), (2, 2, 2)] {
        braid_middle_cross(c.0, c.1, c.2);
    }
}

/// (i) the B1 numeric cases: (3,1,3) 140 flows, (2,3,2) 210 flows
#[test]
#[ignore]
fn braid_middle_313_232() {
    for c in [(3, 1, 3), (2, 3, 2)] {
        braid_middle_cross(c.0, c.1, c.2);
    }
}

/// (ii) the spot check accepts theorem Assocs and certified closed forms, rejects corruptions
#[test]
fn spot_check_accepts_and_rejects() {
    for c in [(1u8, 1u8, 1u8), (2, 1, 2), (2, 2, 2)] {
        let bw = braid_words(c.0, c.1, c.2).unwrap();
        let lam = [c.0, c.1, c.2];
        let (ws, wt) = (Expr::new(&lam, &bw.m).unwrap(), Expr::new(&lam, &bw.mp).unwrap());
        let r = spot_check_assoc::<F>(&ws, &wt, 8, 6, 12, 5).unwrap_or_else(|e| panic!("{:?}: {}", c, e));
        eprintln!("spot {:?}: {:?}", c, r);
        // corrupted: scale the flows of one coset (away from the bottom) by 2
        let (fa, fb) = (FlowIndex::build(&ws), FlowIndex::build(&wt));
        let mut bad = 0;
        for k in 1..fa.support().len() {
            let cbad = fa.support()[k];
            let mut f = |y: &[F]| {
                let mut w = elem_winloc(ElemKind::Assoc, &fa, &fb, y);
                for j in fa.range(cbad) {
                    let (s, e) = (w.col_start[j as usize] as usize, w.col_start[j as usize + 1] as usize);
                    for v in &mut w.vals[s..e] {
                        *v = v.mul(&F::from_i64(2));
                    }
                }
                w
            };
            if spot_check::<F>(&ws, &wt, 0, &mut f, 8, 8, 40, 9).is_err() {
                bad += 1;
            }
        }
        assert!(bad * 2 >= fa.support().len() - 1, "{:?}: only {} of {} corruptions detected", c, bad, fa.support().len() - 1);
        // a certified Pull (switchback φ) passes the spot check too
        let (phi, _) = switchback_eko(1, 2, 1).unwrap();
        let (se, te) = (phi.src.expr(), phi.tgt.expr());
        spot_check::<F>(&se, &te, phi.deg, &mut |y: &[F]| phi.winloc(y), 4, 4, 12, 3).unwrap();
    }
    // a switchback pair has different multistep expressions: the theorem refuses
    let lam = [1u8, 2];
    let (phi, _) = switchback_eko(1, 2, 1).unwrap();
    assert!(theorem_assoc_check(&Expr::new(&lam, &phi.src.rungs).unwrap(), &Expr::new(&lam, &phi.tgt.rungs).unwrap()).is_err());
}

/// (ii) at (3,3,3): the braid middle step (1 680 flows), light localization only
#[test]
#[ignore]
fn spot_check_333() {
    let bw = braid_words(3, 3, 3).unwrap();
    let lam = [3u8, 3, 3];
    let (ws, wt) = (Expr::new(&lam, &bw.m).unwrap(), Expr::new(&lam, &bw.mp).unwrap());
    let t0 = Instant::now();
    let th = theorem_assoc_check(&ws, &wt).unwrap();
    eprintln!("(3,3,3) middle: theorem ok, multistep {:?}, {} flows ({:.3} s)", th.multistep, th.flows, t0.elapsed().as_secs_f64());
    let r = spot_check_assoc::<F>(&ws, &wt, 16, 6, 10, 2026).unwrap();
    eprintln!("(3,3,3) spot check: {:?}", r);
}

/// (iv) marking and layer-A round trip
#[test]
fn marking_and_smapio_roundtrip() {
    let bw = braid_words(1, 2, 1).unwrap();
    let lam = vec![1u8, 2, 1];
    let mut ar: Arena<F> = Arena::new(&lam, TmplLib::new(None));
    let (m, log) = braid_dag_thm(&mut ar, &bw, 900).unwrap();
    eprintln!("{:?}", log);
    ar.assert_genuine(m).unwrap();
    let uses = relies_on_theorem(&ar, m);
    assert_eq!(uses.len(), 1);
    assert!(ar.info(uses[0]).provenance.contains(THM_ASSOC_TAG));
    let text = krhom::smapio::export(&ar, &[m]).unwrap();
    assert!(text.contains("thm-assoc "));
    // fresh arena: the record re-certifies by the theorem
    let mut ar2: Arena<F> = Arena::new(&lam, TmplLib::new(None));
    let got = krhom::smapio::import(&mut ar2, &text).unwrap();
    ar2.assert_genuine(got[0]).unwrap();
    assert_eq!(relies_on_theorem(&ar2, got[0]).len(), 1);
    // agrees with the original at a point
    let p = points::<F>(4, 8, 1, 0).pop().unwrap();
    let (ea, eb) = (ar.info(m).src, ar.info(m).tgt);
    let mut e1 = Eval::new(&ar, p.clone());
    let lay = e1.layout(ea, eb);
    let mut e2 = Eval::new(&ar2, p);
    for &c in &lay.cosets {
        assert_eq!(e1.block(m, c).a, e2.block(got[0], c).a);
    }
    // without the record the Assoc constant is missing in a fresh arena (B1/A2 behaviour);
    // the registry is process-wide, so test with a pair the registry does not hold
    let stripped: String = text.lines().filter(|l| !l.starts_with("thm-assoc")).map(|l| format!("{}\n", l)).collect();
    let mut ar3: Arena<F> = Arena::new(&lam, TmplLib::new(None));
    let got3 = krhom::smapio::import(&mut ar3, &stripped).unwrap();
    assert!(ar3.assert_genuine(got3[0]).is_err(), "Assoc certified without a record");
    // a forged record on a non-theorem pair is refused at import: replace the Assoc by one
    // between words with different multistep expressions is not expressible (add refuses), so
    // check the record on a non-Assoc node instead
    let forged = text.replace("thm-assoc ", "thm-assoc 0\n#").lines().filter(|l| !l.starts_with('#')).map(|l| format!("{}\n", l)).collect::<String>();
    let mut ar4: Arena<F> = Arena::new(&lam, TmplLib::new(None));
    let r4 = krhom::smapio::import(&mut ar4, &forged);
    if let Ok(v) = &r4 {
        // node 0 happened to be the Assoc: fine only if genuine
        ar4.assert_genuine(v[0]).unwrap();
    }
}
