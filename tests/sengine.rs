//! M6a-1a acceptance tests (SPEC_3C_GL6 rev1 §2 "M6a-1a"; SPEC_M4_M6a §6.3-6.4 T5; D14):
//! the localized scanning core on synthetic complexes with known minimal models.
//!
//! * known homology: P (D_min ⊕ cones) P⁻¹ eliminates to D_min (object multiset, K_0, filtered
//!   rank profile per coset, d² = 0), on graded multi-coset types;
//! * the Markowitz heap gives exactly the rescanning order;
//! * streaming (term by term, eligibility rule) = full-then-eliminate = the known minimal model
//!   of C ⊗ W, on 200 seeds; the variant that ignores eligibility is caught (Err(Stream)) or,
//!   with the run-time check off, gives a wrong complex on some seed (sensitivity);
//! * frozen structure: identical traces at 2 points (and in parallel), replay from point 0
//!   reproduces them, and a degenerate point is refused with PointDependent;
//! * memory accounting and the memory refusal; the pivot assertion; API misuse;
//! * the §5.1 library trait on a real (1,1) DecData and `gram_at_point`.
use krhom::field::{Field, F31};
use krhom::sengine::synth::{self, Params};
use krhom::sengine::*;
use krhom::smap::LocMap;
use std::collections::BTreeMap;
use std::sync::Arc;
use std::time::Instant;

type F = F31;

fn xval(i: u64) -> F {
    F::from_i64((i as i64).wrapping_mul(7_919).wrapping_add(123_457) % 2_000_000_000 + 2)
}

fn opts(streaming: bool, order: PivotOrder) -> EngOpts {
    EngOpts { streaming, order, check_d2: true, mem_point_bytes: 1 << 30, ..Default::default() }
}

fn c_params(seed: u64) -> Params {
    let s = seed as usize;
    Params {
        ntypes: 2 + s % 3,
        max_cosets: 3 + (seed % 4),
        max_dim: 1 + (seed % 3) as u32,
        n_min: s % 5,
        n_cones: 1 + (s / 5) % 4,
        hmin: 0,
        hmax: 1 + (s % 2) as i32,
        tmax: 1 + (s % 3) as i32,
        density: 0.5,
        conj: true,
        tid_base: 1,
    }
}
fn w_params(seed: u64) -> Params {
    let s = seed as usize;
    Params {
        ntypes: 2,
        max_cosets: 2 + (seed % 2),
        max_dim: 1 + ((seed / 3) % 2) as u32,
        n_min: 1 + s % 3,
        n_cones: (s / 3) % 3,
        hmin: 0,
        hmax: 1 + (s % 3) as i32, // m = 1..3
        tmax: 2,
        density: 0.5,
        conj: true,
        tid_base: 101,
    }
}

/// C ⊗ W at the point x into a fresh engine (C loaded raw, one crossing of |W levels| steps)
fn run_tensor(seed: u64, x: F, o: EngOpts, mode: Mode) -> Result<(PointScan<F>, SCx<F>, PEngineView), SErr> {
    let mut e = PEngine::<F>::new(mode, o);
    let c = synth::generate(seed, &c_params(seed), x.clone(), &mut e.reg);
    let x2 = x.clone();
    let w = synth::generate(seed ^ 0xabcdef, &w_params(seed), x, &mut e.reg);
    e.load(c.cx)?;
    // 3 of 4 seeds: every term in its own random basis (see synth::tensor_crossing)
    let twist = if seed % 4 != 0 { Some((seed, x2)) } else { None };
    synth::tensor_crossing(&mut e, &w.cx, twist)?;
    let reg = e.reg.clone();
    let (ps, cx) = e.finish(0)?;
    Ok((ps, cx, PEngineView { reg }))
}
struct PEngineView {
    reg: ShapeReg,
}

/// the known minimal model C_min ⊗ W_min at x (no pivots exist in it)
fn minimal_model(seed: u64, x: F) -> (SCx<F>, ShapeReg) {
    let mut e = PEngine::<F>::new(Mode::Record, opts(false, PivotOrder::Heap));
    let c = synth::generate(seed, &c_params(seed), x.clone(), &mut e.reg);
    let w = synth::generate(seed ^ 0xabcdef, &w_params(seed), x, &mut e.reg);
    if c.min.n_objs() == 0 || w.min.n_objs() == 0 {
        return (SCx::new(), e.reg.clone());
    }
    e.load(c.min).unwrap();
    synth::tensor_crossing(&mut e, &w.min, None).unwrap();
    assert_eq!(e.frozen().pivots(), 0, "the minimal model has no pivot");
    let reg = e.reg.clone();
    (e.finish(0).unwrap().1, reg)
}

type Inv = (Vec<SObj>, BTreeMap<(u32, i32), i64>, BTreeMap<(u64, i32, i32), usize>);
fn invariants(cx: &SCx<F>, reg: &ShapeReg) -> Inv {
    (cx.object_multiset(), cx.k0(), synth::rank_profile(cx, reg))
}

// ------------------------------------------------------------------ known homology

#[test]
fn known_homology_single_complexes() {
    let t = Instant::now();
    let (mut pivots, mut contractible, mut objs_in) = (0, 0, 0);
    for seed in 0..150u64 {
        let mut p = c_params(seed);
        if seed % 7 == 0 {
            p.n_min = 0; // contractible
        }
        let x = xval(seed);
        let mut e = PEngine::<F>::new(Mode::Record, opts(seed % 2 == 0, PivotOrder::Heap));
        let s = synth::generate(seed, &p, x, &mut e.reg);
        objs_in += s.cx.n_objs();
        synth::single_crossing(&mut e, &s.cx).unwrap();
        pivots += e.frozen().pivots();
        let reg = e.reg.clone();
        let (ps, cx) = e.finish(0).unwrap();
        assert_eq!(cx.bytes(), cx.recount_bytes());
        assert_eq!(invariants(&cx, &reg), invariants(&s.min, &reg), "seed {}", seed);
        assert!(cx.d2_nonzero(&mut reg.clone()).unwrap().is_empty());
        assert_eq!(ps.k0, s.min.k0());
        if p.n_min == 0 {
            assert_eq!(cx.n_objs(), 0, "contractible complex must eliminate to 0 (seed {})", seed);
            contractible += 1;
        }
        // a minimal complex has no pivot candidate left
        for (&a, outs) in cx.out() {
            for &b in outs.keys() {
                let (oa, ob) = (cx.objs()[&a], cx.objs()[&b]);
                assert!(!(oa.p == ob.p && oa.t == ob.t), "seed {}: candidate left", seed);
            }
        }
    }
    eprintln!("[known_homology] 150 complexes, {} objects in, {} pivots, {} contractible, {:.2} s", objs_in, pivots, contractible, t.elapsed().as_secs_f64());
}

// ------------------------------------------------------------------ heap = rescan

#[test]
fn markowitz_heap_equals_rescan_order() {
    let t = Instant::now();
    let mut piv = 0;
    for seed in 0..200u64 {
        let x = xval(seed);
        for streaming in [false, true] {
            let (a, _, _) = run_tensor(seed, x.clone(), opts(streaming, PivotOrder::Heap), Mode::Record).unwrap();
            let (b, _, _) = run_tensor(seed, x.clone(), opts(streaming, PivotOrder::Rescan), Mode::Record).unwrap();
            assert_eq!(a.frozen.rounds, b.frozen.rounds, "seed {} streaming {}: pivot sequences differ", seed, streaming);
            assert_eq!(a.trace, b.trace);
            piv += a.frozen.pivots();
        }
    }
    eprintln!("[heap=rescan] 400 scans, {} pivots, identical sequences, {:.2} s", piv, t.elapsed().as_secs_f64());
}

// ------------------------------------------------------------------ streaming

#[test]
fn streaming_equals_full_and_known_model() {
    let t = Instant::now();
    let (mut peak_s, mut peak_f, mut nonempty) = (0u64, 0u64, 0);
    for seed in 0..200u64 {
        let x = xval(seed);
        let (ps, cs, vs) = run_tensor(seed, x.clone(), opts(true, PivotOrder::Heap), Mode::Record).unwrap();
        let (pf, cf, vf) = run_tensor(seed, x.clone(), opts(false, PivotOrder::Heap), Mode::Record).unwrap();
        let (cm, rm) = minimal_model(seed, x);
        let (is, if_, im) = (invariants(&cs, &vs.reg), invariants(&cf, &vf.reg), invariants(&cm, &rm));
        assert_eq!(is, if_, "seed {}: streaming ≠ full", seed);
        assert_eq!(is, im, "seed {}: ≠ the known minimal model", seed);
        assert!(cs.d2_nonzero(&mut vs.reg.clone()).unwrap().is_empty());
        assert_eq!(cs.bytes(), cs.recount_bytes());
        assert!(ps.stages[0].rounds > 1 && pf.stages[0].rounds == 1);
        peak_s += ps.mem.peak_bytes;
        peak_f += pf.mem.peak_bytes;
        if cs.n_objs() > 0 {
            nonempty += 1;
        }
    }
    assert!(peak_s <= peak_f, "streaming must not raise the peak");
    eprintln!(
        "[streaming] 200 seeds ({} with nonempty homology): streaming = full = C_min ⊗ W_min; Σ peak bytes streaming {} vs full {} ({:.2}); {:.2} s",
        nonempty,
        peak_s,
        peak_f,
        peak_s as f64 / peak_f as f64,
        t.elapsed().as_secs_f64()
    );
}

#[test]
fn streaming_sensitivity_ignored_eligibility_is_caught() {
    let (mut caught, mut differs, mut early) = (0, 0, 0);
    for seed in 0..200u64 {
        let x = xval(seed);
        let (_, cf, vf) = run_tensor(seed, x.clone(), opts(true, PivotOrder::Heap), Mode::Record).unwrap();
        let reference = invariants(&cf, &vf.reg);
        // run-time check on: the engine must refuse
        let mut o = opts(true, PivotOrder::Heap);
        o.ignore_eligibility = true;
        match run_tensor(seed, x.clone(), o.clone(), Mode::Record) {
            Err(SErr::Stream(_)) => caught += 1,
            Err(e) => panic!("seed {}: unexpected {}", seed, e),
            Ok((_, c, v)) => assert_eq!(invariants(&c, &v.reg), reference, "seed {}: ignoring eligibility went undetected AND changed the result", seed),
        }
        // run-time check off: the result is wrong on some seeds (the generator is strong enough)
        o.check_dead = false;
        o.check_d2 = false;
        match run_tensor(seed, x, o, Mode::Record) {
            Err(_) => early += 1,
            Ok((_, c, mut v)) => {
                if invariants(&c, &v.reg) != reference || !c.d2_nonzero(&mut v.reg).unwrap().is_empty() {
                    differs += 1;
                }
            }
        }
    }
    eprintln!("[sensitivity] ignoring eligibility: refused on {} / 200 seeds; with the check off: wrong result on {}, other error on {}", caught, differs, early);
    assert!(caught >= 1, "the run-time streaming check never fired");
    assert!(differs >= 1, "generator too weak: ignoring eligibility never changed the result");
}

// ------------------------------------------------------------------ frozen structure

#[test]
fn identical_traces_at_two_points_and_replay() {
    let t = Instant::now();
    for seed in 0..60u64 {
        let (x0, x1) = (xval(1000 + seed), xval(2000 + seed));
        let o = opts(seed % 2 == 0, PivotOrder::Heap);
        let (p0, c0, _) = run_tensor(seed, x0, o.clone(), Mode::Record).unwrap();
        let (p1, _, _) = run_tensor(seed, x1.clone(), o.clone(), Mode::Record).unwrap();
        // structure decided independently at two generic points: identical
        assert_eq!(p0.frozen.pivot_sequence(), p1.frozen.pivot_sequence(), "seed {}", seed);
        p0.same_structure(&p1).unwrap();
        // replay of point 0's structure at point 1
        let (r1, c1, _) = run_tensor(seed, x1, o.clone(), Mode::Replay(p0.frozen.clone())).unwrap();
        p0.same_structure(&r1).unwrap();
        assert_eq!(*r1.frozen, *p0.frozen);
        assert_eq!(c0.object_multiset(), c1.object_multiset());
    }
    // parallel replays over 6 points against one frozen trace
    let seed = 17u64;
    let o = opts(true, PivotOrder::Heap);
    let (p0, _, _) = run_tensor(seed, xval(5), o.clone(), Mode::Record).unwrap();
    let fz = p0.frozen.clone();
    let pts: Vec<u64> = (10..16).collect();
    let res = par_map(&pts, |&i| run_tensor(seed, xval(i), o.clone(), Mode::Replay(fz.clone())).map(|r| r.0));
    for r in res {
        p0.same_structure(&r.unwrap()).unwrap();
    }
    eprintln!("[traces] 60 seeds × 2 points + 6 parallel replays identical, {:.2} s", t.elapsed().as_secs_f64());
}

#[test]
fn replay_refuses_a_degenerate_point() {
    let (mut refused, mut tested) = (0, 0);
    for seed in 0..60u64 {
        // a single complex with cones; the degenerate point kills one cone scalar
        let mut p = c_params(seed);
        p.n_cones = p.n_cones.max(1);
        let mk = |x: F, mode: Mode| -> Result<PointScan<F>, SErr> {
            let mut e = PEngine::<F>::new(mode, opts(false, PivotOrder::Heap));
            let s = synth::generate(seed, &p, x, &mut e.reg);
            synth::single_crossing(&mut e, &s.cx)?;
            Ok(e.finish(0)?.0)
        };
        let p0 = mk(xval(seed), Mode::Record).unwrap();
        let roots = {
            let mut reg = ShapeReg::new();
            synth::generate::<F>(seed, &p, xval(seed), &mut reg).cone_roots
        };
        for xs in roots {
            tested += 1;
            match mk(xs.clone(), Mode::Replay(p0.frozen.clone())) {
                Err(SErr::PointDependent { .. }) => refused += 1,
                Err(e) => panic!("seed {}: degenerate point gave {} instead of PointDependent", seed, e),
                Ok(_) => panic!("seed {}: degenerate point {} silently replayed", seed, xs),
            }
            // recording at the degenerate point gives a different structure (more objects)
            let pd = mk(xs, Mode::Record).unwrap();
            assert!(pd.objs.len() > p0.objs.len(), "seed {}: degenerate point did not change the minimal model", seed);
        }
        // the tensor case: C degenerate inside C ⊗ W, streaming
        let mut e = PEngine::<F>::new(Mode::Record, opts(true, PivotOrder::Heap));
        let c = synth::generate::<F>(seed, &c_params(seed), xval(seed), &mut e.reg);
        if let Some(xs) = c.cone_roots.first() {
            let (q0, _, _) = run_tensor(seed, xval(seed), opts(true, PivotOrder::Heap), Mode::Record).unwrap();
            let r = run_tensor(seed, xs.clone(), opts(true, PivotOrder::Heap), Mode::Replay(q0.frozen.clone()));
            assert!(matches!(r, Err(SErr::PointDependent { .. })), "seed {}: tensor replay at a degenerate point not refused", seed);
            tested += 1;
            refused += 1;
        }
    }
    eprintln!("[degenerate] {} degenerate replays, all refused with PointDependent ({})", tested, refused);
    assert!(tested > 50);
}

#[test]
fn replay_refuses_a_foreign_trace() {
    // a trace from another complex (or a truncated one) must be refused, not followed
    let o = opts(true, PivotOrder::Heap);
    let (pa, _, _) = run_tensor(3, xval(1), o.clone(), Mode::Record).unwrap();
    let r = run_tensor(4, xval(1), o.clone(), Mode::Replay(pa.frozen.clone()));
    assert!(matches!(r, Err(SErr::PointDependent { .. })));
    let mut cut = (*pa.frozen).clone();
    cut.rounds.pop();
    let r = run_tensor(3, xval(2), o, Mode::Replay(Arc::new(cut)));
    assert!(matches!(r, Err(SErr::PointDependent { .. })));
}

// ------------------------------------------------------------------ memory

#[test]
fn memory_accounting_and_refusal() {
    let seed = 11u64;
    let big = |budget: u64, streaming: bool| -> Result<PointScan<F>, SErr> {
        let mut o = opts(streaming, PivotOrder::Heap);
        o.mem_point_bytes = budget;
        o.check_d2 = false;
        let mut e = PEngine::<F>::new(Mode::Record, o);
        let cp = Params { ntypes: 4, max_cosets: 6, max_dim: 3, n_min: 8, n_cones: 10, hmax: 3, tmax: 4, ..Default::default() };
        let wp = Params { ntypes: 3, max_cosets: 4, max_dim: 2, n_min: 4, n_cones: 6, hmax: 3, tmax: 3, tid_base: 101, ..Default::default() };
        let c = synth::generate(seed, &cp, xval(1), &mut e.reg);
        let w = synth::generate(seed + 1, &wp, xval(1), &mut e.reg);
        e.set_extern_bytes(4096);
        e.load(c.cx)?;
        synth::tensor_crossing(&mut e, &w.cx, Some((seed, xval(1))))?;
        assert_eq!(e.cx().bytes(), e.cx().recount_bytes());
        Ok(e.finish(0)?.0)
    };
    let full = big(1 << 30, false).unwrap();
    let strm = big(1 << 30, true).unwrap();
    assert_eq!(full.objs.len(), strm.objs.len());
    let (pf, ps) = (full.mem.peak_bytes, strm.mem.peak_bytes);
    eprintln!("[memory] peak full {} B, streaming {} B; stage: {:?}", pf, ps, strm.stages[0]);
    assert!(ps > 4096 && ps <= pf);
    // a budget between the two peaks: full refused, streaming passes
    if ps < pf {
        let mid = (ps + pf) / 2;
        assert!(big(mid, true).is_ok());
        match big(mid, false) {
            Err(e @ SErr::Memory { .. }) => {
                let s = e.to_string();
                assert!(s.contains("per-point memory") && s.contains("run with streaming or on the cloud"), "{}", s);
                eprintln!("[memory] refusal: {}", s);
            }
            r => panic!("expected the memory refusal, got {:?}", r.map(|p| p.mem.peak_bytes)),
        }
    }
    // a tiny budget is refused for both
    assert!(matches!(big(20_000, true), Err(SErr::Memory { .. })));
    assert!(matches!(big(20_000, false), Err(SErr::Memory { .. })));
}

// ------------------------------------------------------------------ assertions and API

#[test]
fn pivot_assertion_and_api_misuse() {
    let mut e = PEngine::<F>::new(Mode::Record, opts(true, PivotOrder::Heap));
    let s = TypeShape::new(7, vec![0, 2], vec![1, 2], 0).unwrap();
    e.register(s.clone()).unwrap();
    assert!(matches!(e.add_obj(SObj { p: 7, h: 0, t: 0 }), Err(SErr::Usage(_))));
    e.begin_crossing().unwrap();
    e.begin_step().unwrap();
    let a = e.add_obj(SObj { p: 7, h: 0, t: 0 }).unwrap();
    let b = e.add_obj(SObj { p: 7, h: 1, t: 0 }).unwrap();
    let lay = e.layout(7, 7).unwrap();
    let mut m = LocMap::<F>::zeros(lay);
    m.block_mut(0)[0] = F::from_i64(2);
    m.block_mut(1).clone_from_slice(&[F::from_i64(3), F::zero(), F::zero(), F::from_i64(3)]);
    assert!(e.add_map(b, a, m.clone()).is_err(), "h must rise by 1");
    e.add_map(a, b, m).unwrap();
    e.end_step().unwrap();
    assert!(matches!(e.end_crossing(), Err(SErr::Pivot(_))), "non-scalar degree-0 component must be refused");
    // unregistered type, wrong layout
    let mut e = PEngine::<F>::new(Mode::Record, opts(true, PivotOrder::Heap));
    e.register(s).unwrap();
    e.register(TypeShape::new(8, vec![2], vec![1], 2).unwrap()).unwrap();
    assert!(e.register(TypeShape::new(8, vec![3], vec![1], 3).unwrap()).is_err());
    e.begin_crossing().unwrap();
    e.begin_step().unwrap();
    assert!(e.add_obj(SObj { p: 9, h: 0, t: 0 }).is_err());
    let a = e.add_obj(SObj { p: 7, h: 0, t: 1 }).unwrap();
    let b = e.add_obj(SObj { p: 8, h: 1, t: 0 }).unwrap();
    let wrong = LocMap::<F>::zeros(e.layout(7, 7).unwrap());
    assert!(matches!(e.add_map(a, b, wrong), Err(SErr::Shape(_))));
    assert!(e.end_crossing().is_err(), "end_crossing with an open step");
    // shapes
    assert!(TypeShape::new(1, vec![0, 1], vec![2, 2], 0).is_err(), "rep must have m = 1");
    assert!(TypeShape::new(1, vec![1, 0], vec![1, 1], 0).is_err());
}

// ------------------------------------------------------------------ the library trait

#[test]
fn library_trait_and_gram_at_point_on_11() {
    use krhom::sindec::{provider_for_stamped, SIndec, SIndecOpts};
    use krhom::smap::{points, Eval};
    let lam = [1u8, 1];
    let costs = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/homprov/costs_cloud.txt");
    let (prov, stamp) = provider_for_stamped::<F>(&lam, Some(&costs), 1.0);
    let mut sd = SIndec::with_stamp(&lam, prov, stamp, SIndecOpts::default());
    let order = sd.schur.canonical_basis(&lam).order.clone();
    let mut reps = vec![];
    for &p in &order {
        sd.ty(p).unwrap();
        if sd.schur.table(p).cols() == lam {
            sd.closure_gens(p).unwrap();
            reps.push((p, rep_coset(&sd.schur, &lam, p)));
        }
    }
    let dd = sd.freeze();
    let lib: &dyn SLib<F> = &*dd;
    assert_eq!(lib.lam(), &lam);
    let pts = points::<F>(lib.ar().k(), 0x5e9e, 2, 0);
    for &(p, rep) in &reps {
        let sh = lib_shape(lib, p, rep).unwrap();
        assert_eq!(sh.dim_at(rep), 1);
        let mut grams = vec![];
        for pt in &pts {
            let mut ev = Eval::new(lib.ar(), pt.clone());
            let g = gram_at_point(lib, &mut ev, p).unwrap();
            let n = g.ins.len();
            assert_eq!(g.outs.len(), n);
            assert_eq!(krhom::dense::rank(&g.g, n), n, "Gram of type {} singular at point {}", p, pt.id);
            grams.push(g);
        }
        eprintln!("[gram] (1,1) type {}: {} generators, degrees in {:?} out {:?}, rank S = {}, id coset {}", p, grams[0].ins.len(), grams[0].ins, grams[0].outs, sh.s(), grams[0].id_coset);
        if grams[0].ins.iter().chain(&grams[0].outs).any(|&d| d != 0) {
            assert_ne!(grams[0].g, grams[1].g, "Gram must depend on the point (degrees > 0 present)");
        }
    }
    // the trait object is usable from several threads over one Arc
    fn ss<T: Send + Sync + ?Sized>() {}
    ss::<dyn SLib<F>>();
}

// ------------------------------------------------------------------ measured costs

fn bench(seed: u64, cp: &Params, wp: &Params, streaming: bool, twist: bool) -> (PointScan<F>, f64) {
    let mut o = opts(streaming, PivotOrder::Heap);
    o.check_d2 = false;
    let mut e = PEngine::<F>::new(Mode::Record, o);
    let c = synth::generate(seed, cp, xval(1), &mut e.reg);
    let w = synth::generate(seed + 1, wp, xval(1), &mut e.reg);
    e.load(c.cx).unwrap();
    let t = Instant::now();
    synth::tensor_crossing(&mut e, &w.cx, if twist { Some((seed, xval(1))) } else { None }).unwrap();
    let _ = t;
    let s = e.stages[0].secs;
    (e.finish(0).unwrap().0, s)
}

fn cost_table(sizes: &[(usize, usize, usize, usize)], twist: bool) {
    eprintln!("| C objs (min+2·cones) | W objs | inserted | final | pivots | peak entries | Σ S final | c_e final | zero-block frac | peak bytes full | peak bytes stream | secs full | secs stream | replay secs |");
    for &(cm, cc, wm, wc) in sizes {
        let cp = Params { ntypes: 4, max_cosets: 6, max_dim: 3, n_min: cm, n_cones: cc, hmax: 3, tmax: 5, density: 0.3, ..Default::default() };
        let wp = Params { ntypes: 3, max_cosets: 4, max_dim: 2, n_min: wm, n_cones: wc, hmax: 3, tmax: 3, density: 0.3, tid_base: 101, ..Default::default() };
        let (f, sf) = bench(5, &cp, &wp, false, twist);
        let (s, ss) = bench(5, &cp, &wp, true, twist);
        assert_eq!(f.objs.len(), s.objs.len());
        // replay cost at a second point
        let mut o = opts(true, PivotOrder::Heap);
        o.check_d2 = false;
        let mut e = PEngine::<F>::new(Mode::Replay(s.frozen.clone()), o);
        let c = synth::generate(5, &cp, xval(1), &mut e.reg);
        let w = synth::generate(6, &wp, xval(1), &mut e.reg);
        e.load(c.cx).unwrap();
        let t = Instant::now();
        synth::tensor_crossing(&mut e, &w.cx, if twist { Some((5, xval(1))) } else { None }).unwrap();
        let _ = t;
        let sr = e.stages[0].secs;
        let st = &s.stages[0];
        eprintln!(
            "| {} | {} | {} | {} | {} | {} | {} | {:.2} | {:.2} | {} | {} | {:.3} | {:.3} | {:.3} |",
            cm + 2 * cc,
            wm + 2 * wc,
            st.inserted,
            st.objects,
            st.pivots,
            f.mem.peak_entries.max(s.mem.peak_entries),
            st.sum_s,
            st.c_e,
            st.zero_block_frac,
            f.mem.peak_bytes,
            s.mem.peak_bytes,
            sf,
            ss,
            sr
        );
    }
}

#[test]
fn measured_costs_small() {
    cost_table(&[(6, 6, 3, 3), (10, 12, 4, 5)], true);
}

#[test]
#[ignore]
fn measured_costs_large() {
    cost_table(&[(10, 20, 6, 8)], true);
    cost_table(&[(20, 40, 8, 12), (30, 60, 10, 16), (40, 80, 12, 20)], false);
}
