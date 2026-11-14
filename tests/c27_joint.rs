//! C2-7: joint junction words (Leaves::joint_descent), Prop W checks W1/W2 and per-predecessor
//! final words; pure combinatorics (no maps are built here).
use krhom::dcoset::{junction_types, StepCache};
use krhom::leaves::{plan_steps, Leaves};
use krhom::schur::{Schur, Table};
use krhom::sindec::canonical_ladder;

fn plan_max(lam: &[u8], p: &str, rung: &str, q: &str, joint: bool) -> (u128, usize) {
    let mut s = Schur::new(lam);
    let id = s.intern(Table::identity(lam));
    let p = s.intern(Table::parse(p));
    let q = s.intern(Table::parse(q));
    let mut x = canonical_ladder(&mut s, p);
    x.push(krhom::smapio::parse_rung(rung).unwrap());
    let mut lv = Leaves::default();
    lv.joint = joint;
    let mut sc = StepCache::new();
    let js = junction_types(&mut s, &mut sc, &x, &[q]);
    let plan = lv.plan_words(&mut s, &mut sc, &x, &js, id);
    let finals = lv.final_pred_words(&mut s, &mut sc, &x, &plan, id);
    let steps = plan_steps(&mut s, &mut sc, &x, &plan);
    // a replaced final step is built into the per-predecessor word W' instead: it counts with the
    // replacement step's own W-aware score (review leaves #1), not as solved
    let mut wlc = std::collections::HashMap::new();
    let mut repl = std::collections::HashMap::new();
    for (&(qq, i), wa) in &finals {
        let (tp, wp) = plan[x.len() - 1][i].clone();
        let g = sc.step(&mut s, tp, &x[x.len() - 1]).into_iter().find(|e| e.0 == qq).map(|e| e.1).unwrap();
        let sw = krhom::leaves::step_score_w(&mut wlc, lam, &wp, &x[x.len() - 1], wa, &g);
        let e = repl.entry((qq, tp)).or_insert(0u128);
        *e = (*e).max(sw);
    }
    let mx = steps.iter().map(|st| if st.level == x.len() { repl.get(&(st.to, st.from)).copied().unwrap_or(st.score_w) } else { st.score_w }).max().unwrap();
    (mx, finals.len())
}

/// pair 1235 of the (2,2,2,2) prune-4 library (C2-5/C2-4 blocker): Lemma A step 83.7 M unknowns
/// with the default plan, ≤ 1e5 with joint words (C2-7 census, runs_c2/c27/steps_all*.txt)
#[test]
fn joint_words_pair_1235() {
    let lam = [2u8, 2, 2, 2];
    let (p, r, q) = ("0,1,0,1/1,0,1,0/1,1,0,0/0,1,0,1", "F1.1", "0,0,2,0/1,0,0,1/1,1,0,0/0,1,0,1");
    let (off, nf) = plan_max(&lam, p, r, q, false);
    assert_eq!(nf, 0, "per-predecessor finals are off without joint");
    assert!(off > 10_000_000, "default plan: {}", off);
    let (on, _) = plan_max(&lam, p, r, q, true);
    assert!(on <= 100_000, "joint plan: {}", on);
}

/// pair 1444 (rerouted call of C2-5's failure list): needs a per-predecessor final word
#[test]
fn per_predecessor_final_pair_1444() {
    let lam = [2u8, 2, 2, 2];
    let (on, nf) = plan_max(&lam, "0,1,0,1/1,0,1,0/1,1,0,0/0,0,1,1", "E1.2", "0,2,0,0/1,0,0,1/1,1,0,0/0,1,0,1", true);
    assert!(nf >= 1);
    assert!(on <= 100_000, "joint plan with per-predecessor finals: {}", on);
}

/// W1/W2: a suffix window (the step rung in the common suffix) is not W-certified, so the
/// W-aware score falls back to the whole words (the C2-5 failure of 1444's rerouted call)
#[test]
fn w12_rejects_suffix_window() {
    let lam = [2u8, 2, 2, 2];
    let mut s = Schur::new(&lam);
    let tp = s.intern(Table::parse("0,1,1,0/1,0,1,0/1,1,0,0/0,0,1,1"));
    let t = s.intern(Table::parse("0,2,0,0/1,0,1,0/1,1,0,0/0,1,0,1"));
    let pw = |x: &str| krhom::smapio::parse_rungs(x).unwrap();
    let wfrom = pw("E0.1,F0.1,E1.2,E0.1,F1.2,F0.1,E2.1");
    let rho = pw("E1.2")[0];
    let wto = pw("E1.1,E0.2,F1.1,F0.2,E2.1,E1.2");
    let mut sc = StepCache::new();
    let g = sc.step(&mut s, tp, &rho).into_iter().find(|x| x.0 == t).unwrap().1;
    let mut wl = std::collections::HashMap::new();
    assert!(!krhom::leaves::window_w12(&mut wl, &lam, &wfrom, &rho, &wto, &g));
    let sw = krhom::leaves::step_score_w(&mut wl, &lam, &wfrom, &rho, &wto, &g);
    let s0 = krhom::leaves::step_score(&lam, &wfrom, &rho, &wto, &g);
    assert!(s0 < 100_000 && sw > 1_000_000, "score {} score_w {}", s0, sw);
}
