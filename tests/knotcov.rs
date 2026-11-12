//! M5c: library jobs of a pruned scan of a braid word (`sindec::knot_jobs`).
use krhom::schur::{decat_scan, ScanOpts, Schur, Table};
use krhom::sindec::{knot_jobs, pruned_plan};
use std::collections::BTreeSet;

/// Unpruned: the lower-bound final types (+ identity) are exactly the types of the signed K_0
/// class of `decat_scan`; the lower bound is inside the upper bound; every job is a plan pair.
#[test]
fn knot_jobs_lb_matches_decat_scan() {
    for (lam, word) in [(vec![2u8, 2, 2], vec![1, -2, 1, -2]), (vec![1u8, 1, 1], vec![1, 2, 1, 2, 2]), (vec![2u8, 2, 2], vec![1, 1, 1, 2, -1, 2])] {
        let mut s = Schur::new(&lam);
        let plan = pruned_plan(&mut s, None);
        let kj = knot_jobs(&mut s, &plan, &word).unwrap();
        let r = decat_scan(&mut s, &word, &ScanOpts { peak: false, every_stage: false, marked_block: None });
        let id = s.intern(Table::identity(&lam));
        let k0: BTreeSet<_> = r.k0.iter().filter(|e| e.2 != 0).map(|e| e.0).filter(|&p| p != id).collect();
        assert_eq!(kj.closure_lb.iter().copied().collect::<BTreeSet<_>>(), k0, "{:?} {:?}", lam, word);
        assert!(kj.closure_lb.is_subset(&kj.closure));
        assert!(kj.pairs_lb.is_subset(&kj.pairs));
        assert!(kj.pairs.iter().all(|&k| k < plan.pairs.len() && !plan.pairs[k].p_killed));
    }
}

/// Pruned (3,3,3): 4_1 at N = 5 visits 89 pairs and 19 closure types (M5c_coverage.md), and the
/// pruned job sets are inside the level sets (every visited source is kept).
#[test]
fn knot_jobs_333_fig8() {
    let lam = [3u8, 3, 3];
    let mut s = Schur::new(&lam);
    let plan = pruned_plan(&mut s, Some(5));
    let kj = knot_jobs(&mut s, &plan, &[1, -2, 1, -2]).unwrap();
    assert_eq!((kj.pairs.len(), kj.closure.len()), (89, 19));
    assert_eq!((kj.pairs_lb.len(), kj.closure_lb.len()), (89, 15));
    for &q in &kj.closure {
        assert!(!s.killed_at(q, 5));
    }
}
