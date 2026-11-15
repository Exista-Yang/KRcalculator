//! M4c B4w: the Leaves pair route wired into the production provider (`leaves::leaf_pairs`,
//! `ClosureLeaves::pairs`, `provider_for_stamped_m4c`), through SIndec's full checks
//! (certify_pairs, dechom_postcheck, D14 scaling, Π∘Ι = id, kept multiset), layer-A save and a
//! replay with a provider that refuses every call (every map re-certified on load).
use krhom::algebra::field::{Field, F31};
use krhom::library::homprov::*;
use krhom::library::leaves::*;
use krhom::combinat::schur::{Schur, Tid};
use krhom::library::sindec::*;
use krhom::colored::smap::{Arena, ExprId, MapId};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

/// The production chain with every pair call refused (so that each pair target goes to Leaves).
struct RefusePairs<F: Field> {
    inner: Box<dyn HomProvider<F>>,
}
impl<F: Field> HomProvider<F> for RefusePairs<F> {
    fn name(&self) -> &'static str {
        "RefusePairs"
    }
    fn full(&mut self, ar: &mut Arena<F>, s: &mut Schur, a: ExprId, b: ExprId, d: i32, bu: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        self.inner.full(ar, s, a, b, d, bu)
    }
    fn pairs(&mut self, _: &mut Arena<F>, _: &mut Schur, _: ExprId, _: Tid, _: ExprId, _: Option<MapId>, sdeg: i32, mult: usize, _: &Budget) -> Result<(Pairs<F>, SearchStats), Vec<Missing>> {
        Err(vec![Missing { what: "test: pairs refused".into(), degree: sdeg, stratum: None, got: 0, want: mult, kind: MissKind::Refused }])
    }
    fn closure(&mut self, ar: &mut Arena<F>, s: &mut Schur, lp: ExprId, p: Tid, ep: Option<MapId>, outs: bool, b: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        self.inner.closure(ar, s, lp, p, ep, outs, b)
    }
    fn end0(&mut self, ar: &mut Arena<F>, s: &mut Schur, l: ExprId, b: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        self.inner.end0(ar, s, l, b)
    }
    fn failed_stats(&self) -> Option<SearchStats> {
        None
    }
}

/// Counts successful pair calls of the wrapped provider; refuses everything when `refuse_all`.
struct Count<F: Field> {
    inner: Option<Box<dyn HomProvider<F>>>,
    ok: Arc<AtomicUsize>,
    calls: Arc<AtomicUsize>,
}
fn refused<T>() -> Result<T, Vec<Missing>> {
    Err(vec![Missing { what: "replay: refused".into(), degree: 0, stratum: None, got: 0, want: 1, kind: MissKind::Refused }])
}
impl<F: Field> HomProvider<F> for Count<F> {
    fn name(&self) -> &'static str {
        "Count"
    }
    fn full(&mut self, ar: &mut Arena<F>, s: &mut Schur, a: ExprId, b: ExprId, d: i32, bu: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        self.calls.fetch_add(1, Ordering::SeqCst);
        match &mut self.inner {
            Some(p) => p.full(ar, s, a, b, d, bu),
            None => refused(),
        }
    }
    fn pairs(&mut self, ar: &mut Arena<F>, s: &mut Schur, x: ExprId, q: Tid, lq: ExprId, eq: Option<MapId>, sdeg: i32, mult: usize, b: &Budget) -> Result<(Pairs<F>, SearchStats), Vec<Missing>> {
        self.calls.fetch_add(1, Ordering::SeqCst);
        let r = match &mut self.inner {
            Some(p) => p.pairs(ar, s, x, q, lq, eq, sdeg, mult, b),
            None => refused(),
        };
        if r.is_ok() {
            self.ok.fetch_add(1, Ordering::SeqCst);
        }
        r
    }
    fn closure(&mut self, ar: &mut Arena<F>, s: &mut Schur, lp: ExprId, p: Tid, ep: Option<MapId>, outs: bool, b: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        self.calls.fetch_add(1, Ordering::SeqCst);
        match &mut self.inner {
            Some(pr) => pr.closure(ar, s, lp, p, ep, outs, b),
            None => refused(),
        }
    }
    fn end0(&mut self, ar: &mut Arena<F>, s: &mut Schur, l: ExprId, b: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        self.calls.fetch_add(1, Ordering::SeqCst);
        match &mut self.inner {
            Some(p) => p.end0(ar, s, l, b),
            None => refused(),
        }
    }
    fn failed_stats(&self) -> Option<SearchStats> {
        None
    }
}

fn run(lam: &[u8], max_pairs: usize) {
    let dir = std::env::temp_dir().join(format!("krhom_b4w_{}_{}", krhom::library::smapio::labels_str(lam).replace(',', ""), std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let (chain, cstamp) = provider_for_stamped::<F31>(lam, None, 1.0);
    let cl = ClosureLeaves::new(Box::new(RefusePairs { inner: chain }), Leaves::new(), false, true);
    let stamp = format!("{} test-refuse-pairs {}", cstamp, LEAF_PAIRS_STAMP);
    let (ok, calls) = (Arc::new(AtomicUsize::new(0)), Arc::new(AtomicUsize::new(0)));
    let prov = Count::<F31> { inner: Some(Box::new(cl)), ok: ok.clone(), calls: calls.clone() };
    let opts = SIndecOpts { cache_dir: Some(dir.clone()), ..SIndecOpts::default() };
    let mut sd: SIndec<F31> = SIndec::with_stamp(lam, Box::new(prov), stamp.clone(), opts.clone());
    let plan = pruned_plan(&mut sd.schur, None);
    let jobs: Vec<PlanPair> = plan.pairs.iter().filter(|pp| !pp.identity && !pp.all_killed && !pp.kept.is_empty()).take(max_pairs).cloned().collect();
    let mut failed = vec![];
    for pp in &jobs {
        if let Err(e) = sd.dec_rung(pp.p, pp.rung) {
            failed.push(format!("{} {:?}: {}", sd.schur.table(pp.p).fmt(), pp.rung, e));
        }
    }
    let st = sd.stats();
    eprintln!("[b4w] {:?}: {} jobs, {} failed, Leaves pair calls ok {}; stats saved {} loaded {}", lam, jobs.len(), failed.len(), ok.load(Ordering::SeqCst), st.layer_a_saved, st.layer_a_loaded);
    for f in &failed {
        eprintln!("[b4w]   FAILED {}", f);
    }
    assert!(failed.is_empty(), "{} jobs failed", failed.len());
    assert!(ok.load(Ordering::SeqCst) > 0);
    // replay: a fresh SIndec on the same cache with a provider refusing every call
    let (ok2, calls2) = (Arc::new(AtomicUsize::new(0)), Arc::new(AtomicUsize::new(0)));
    let none = Count::<F31> { inner: None, ok: ok2, calls: calls2.clone() };
    let mut sd2: SIndec<F31> = SIndec::with_stamp(lam, Box::new(none), stamp, opts);
    // type ids are per Schur instance: recompute the plan in the replay's Schur
    let plan2 = pruned_plan(&mut sd2.schur, None);
    let jobs2: Vec<PlanPair> = plan2.pairs.iter().filter(|pp| !pp.identity && !pp.all_killed && !pp.kept.is_empty()).take(max_pairs).cloned().collect();
    assert_eq!(jobs2.len(), jobs.len());
    for pp in &jobs2 {
        sd2.dec_rung(pp.p, pp.rung).unwrap_or_else(|e| panic!("replay {} {:?}: {}", sd2.schur.table(pp.p).fmt(), pp.rung, e));
    }
    let st2 = sd2.stats();
    eprintln!("[b4w] {:?} replay: loaded {}, rejected {}, provider calls {}", lam, st2.layer_a_loaded, st2.layer_a_rejected, calls2.load(Ordering::SeqCst));
    assert_eq!(st2.layer_a_rejected, 0);
    assert!(st2.layer_a_loaded > 0);
    let _ = std::fs::remove_dir_all(&dir);
}

/// (2,1,2): every non-identity pair job of the G1 plan with the production chain refusing every
/// pair call: all targets come from `leaf_pairs` (rex-class and Lemma A), pass SIndec's checks,
/// are saved to layer A and replay with no provider.
#[test]
fn leaf_pairs_through_sindec_212() {
    run(&[2, 1, 2], usize::MAX);
}

/// The production constructor carries the new stamp and asks Leaves for pairs.
#[test]
fn m4c_stamp() {
    let (_, st) = provider_for_stamped_m4c::<F31>(&[2, 2, 2], None, 1.0, false);
    assert!(st.contains(LEAF_PAIRS_STAMP) && st.contains(LEAVES_STAMP) && st.contains("leaves-last"));
    let (_, st3) = provider_for_stamped_m4c_a3::<F31>(&[2, 2, 2], None, 1.0, false);
    assert!(!st3.contains(LEAF_PAIRS_STAMP));
}
