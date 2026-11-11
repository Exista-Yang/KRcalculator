//! MBOX test (e): the staged rank stage (`cclose::try_close_staged`) on synthetic data of 4_1
//! size stays within its memory estimate (`cclose::staged_mem_estimate`). The peak is measured
//! with a counting global allocator (bytes requested, all threads); this binary holds only these
//! tests, which take a lock so that they never overlap.
//!
//! Data: λ = (3,3,3); the final types of `schur::decat_scan` on 4_1 = σ_1σ_2^{-1}σ_1σ_2^{-1} kept
//! at N = 6 (not `killed_at` 6); per type the ins = outs degrees are the exponents of c_diag(b_p)
//! (the measured generator degrees, as xengine asserts); Gram entries are random sparse
//! homogeneous polynomials of weight (d_a + d_b)/2; objects follow the K_0 class (|c| copies of
//! each (type, q), h = 0 or 1 by sign, spread over h = −1..2) and P entries are random sparse
//! homogeneous of the bookkeeping weight between objects of adjacent h with |q_s − q_t| ≤ 4.

use krhom::cclose::*;
use krhom::field::{Field, F31};
use krhom::schur::{decat_scan, ScanOpts, Schur};
use rand::{rngs::StdRng, Rng, SeedableRng};
use std::alloc::{GlobalAlloc, Layout, System};
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

struct Counting;
static CUR: AtomicUsize = AtomicUsize::new(0);
static PEAK: AtomicUsize = AtomicUsize::new(0);

unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, l: Layout) -> *mut u8 {
        let p = System.alloc(l);
        if !p.is_null() {
            let c = CUR.fetch_add(l.size(), Ordering::Relaxed) + l.size();
            PEAK.fetch_max(c, Ordering::Relaxed);
        }
        p
    }
    unsafe fn dealloc(&self, p: *mut u8, l: Layout) {
        System.dealloc(p, l);
        CUR.fetch_sub(l.size(), Ordering::Relaxed);
    }
    unsafe fn realloc(&self, p: *mut u8, l: Layout, n: usize) -> *mut u8 {
        let q = System.realloc(p, l, n);
        if !q.is_null() {
            if n > l.size() {
                let c = CUR.fetch_add(n - l.size(), Ordering::Relaxed) + (n - l.size());
                PEAK.fetch_max(c, Ordering::Relaxed);
            } else {
                CUR.fetch_sub(l.size() - n, Ordering::Relaxed);
            }
        }
        q
    }
}

#[global_allocator]
static A: Counting = Counting;

static LOCK: Mutex<()> = Mutex::new(());

/// Peak bytes above the current level while running f.
fn measure<R>(f: impl FnOnce() -> R) -> (R, usize) {
    let base = CUR.load(Ordering::SeqCst);
    PEAK.store(base, Ordering::SeqCst);
    let r = f();
    (r, PEAK.load(Ordering::SeqCst) - base)
}

fn rand_poly(ring: &ERing, w: i32, nterms: usize, rng: &mut StdRng, cache: &mut BTreeMap<i32, Vec<u128>>) -> EPoly<F31> {
    if w < 0 {
        return EPoly::zero();
    }
    let monos = cache.entry(w).or_insert_with(|| ring.monos_of_weight(w as u32));
    let t = (0..nterms).map(|_| (monos[rng.gen_range(0..monos.len())], F31::from_i64(rng.gen_range(1..1_000_000_007)))).collect();
    EPoly::from_terms(t)
}

/// Synthetic closure data of 4_1 size at N = 6 (see the module doc); (data, pairing, #types).
fn synthetic_41(seed: u64, gram_terms: usize, max_w: i32) -> (ClosureDataC<F31>, Pairing<F31>, usize) {
    let lam = [3u8, 3, 3];
    let ring = ERing::new(&lam);
    let mut s = Schur::new(&lam);
    let word = [1, -2, 1, -2];
    let res = decat_scan(&mut s, &word, &ScanOpts { peak: false, every_stage: false, marked_block: None });
    let mut rng = StdRng::seed_from_u64(seed);
    let mut mc: BTreeMap<i32, Vec<u128>> = BTreeMap::new();
    let mut types = BTreeMap::new();
    let mut objs = vec![];
    for &(p, k, c) in &res.k0 {
        if s.killed_at(p, 6) {
            continue;
        }
        if !types.contains_key(&p) {
            let cd = s.c_diag(p);
            let mut degs = vec![];
            for (i, &m) in cd.c.iter().enumerate() {
                for _ in 0..m {
                    degs.push(cd.lo + i as i32);
                }
            }
            let gram: Vec<Vec<EPoly<F31>>> = degs
                .iter()
                .map(|&db| degs.iter().map(|&da| { let w = (da + db) / 2; if w > max_w { EPoly::zero() } else { rand_poly(&ring, w, gram_terms, &mut rng, &mut mc) } }).collect())
                .collect();
            types.insert(p, TypeBasisC { ins: degs.clone(), outs: degs, gram });
        }
        for j in 0..c.unsigned_abs() {
            let h = if c > 0 { 2 * (j as i32 % 2) } else { 1 - 2 * (j as i32 % 2) };
            objs.push((p, h, k));
        }
    }
    let ntypes = types.len();
    let data = ClosureDataC { lam: lam.to_vec(), objs: objs.clone(), types, edges: vec![], crossings: vec![] };
    let mut pairing = vec![];
    for (si, s_) in objs.iter().enumerate() {
        for (ti, t_) in objs.iter().enumerate() {
            if t_.1 != s_.1 + 1 || (s_.2 - t_.2).abs() > 4 || rng.gen_range(0..4) != 0 {
                continue;
            }
            let (ts, tt) = (&data.types[&s_.0], &data.types[&t_.0]);
            let p: Vec<Vec<EPoly<F31>>> = tt
                .outs
                .iter()
                .map(|&ob| ts.ins.iter().map(|&ia| { let d = s_.2 - t_.2 + ia + ob; if d % 2 != 0 || d / 2 > max_w { EPoly::zero() } else { rand_poly(&ring, d / 2, 2, &mut rng, &mut mc) } }).collect())
                .collect();
            pairing.push((si, ti, p));
        }
    }
    (data, pairing, ntypes)
}

fn run_case(marked: Option<usize>, max_w: i32, cap_mb: u64) {
    let _g = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let t = std::time::Instant::now();
    let (data, pairing, ntypes) = synthetic_41(41, 3, max_w);
    let ins: usize = data.types.values().map(|t| t.ins.len()).sum();
    let maxn = data.types.values().map(|t| t.ins.len()).max().unwrap();
    assert!((34..=40).contains(&ntypes), "4_1 kept types at N = 6: {}", ntypes);
    let est = staged_mem_estimate(&data, &pairing, 6, marked, None);
    let opts = StagedOpts { cap: Some((cap_mb as u128) << 20), rank_cap: None };
    let (pc, peak) = measure(|| try_close_staged(&data, &pairing, 6, marked, &opts).unwrap());
    let dims = try_staged_type_ranks(&data, 6, marked, &opts).unwrap();
    let chain: usize = data.objs.iter().map(|o| dims[&o.0].values().sum::<usize>()).sum();
    // with the rank cap (the uncapped ranks, as T0 would give them)
    let capo = StagedOpts { cap: opts.cap, rank_cap: Some(dims.clone()) };
    let est_c = staged_mem_estimate(&data, &pairing, 6, marked, Some(&dims));
    let (pc2, peak_c) = measure(|| try_close_staged(&data, &pairing, 6, marked, &capo).unwrap());
    assert_eq!(pc, pc2, "capped = uncapped");
    let old = exact_mem_estimate(&{ let mut d = data.clone(); d.edges.clear(); d }, 6, marked);
    eprintln!(
        "(e) 4_1-size synthetic (3,3,3) N=6 marked {:?}: {} types, {} ins (max {}), {} objects, {} pairing edges, chain dim {}; staged peak {:.1} MB <= estimate {:.1} MB (ratio {:.2}); capped peak {:.1} MB <= {:.1} MB; oracle exact_mem_estimate {:.1} MB; total rank {}; {:.1}s",
        marked,
        ntypes,
        ins,
        maxn,
        data.objs.len(),
        pairing.len(),
        chain,
        peak as f64 / 1048576.0,
        est as f64 / 1048576.0,
        peak as f64 / est as f64,
        peak_c as f64 / 1048576.0,
        est_c as f64 / 1048576.0,
        old as f64 / 1048576.0,
        pc.total_rank(),
        t.elapsed().as_secs_f64()
    );
    assert!((peak as u128) <= est, "staged peak {} above its estimate {}", peak, est);
    assert!((peak_c as u128) <= est_c, "capped staged peak {} above its estimate {}", peak_c, est_c);
}

#[test]
fn e_staged_rank_stage_41_size_reduced_within_estimate() {
    run_case(Some(0), 18, 4096);
}

#[test]
#[ignore = "MBOX (e): unreduced N = 6 at 4_1 size (dim A = 8000)"]
fn e_staged_rank_stage_41_size_unreduced_refused_before_allocating() {
    // Uncapped unreduced N = 6 at 4_1 size: the dense-charge model asks for terabytes (state
    // spaces up to 54 x 777 rows per degree), so the stage must refuse with Err before it
    // allocates (open problem: MOY caps from T0 plus a sparse fill model, R3 O4).
    let _g = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let (data, pairing, _) = synthetic_41(41, 3, 27);
    let est = staged_mem_estimate(&data, &pairing, 6, None, None);
    let opts = StagedOpts { cap: Some(14000u128 << 20), rank_cap: None };
    let (r, peak) = measure(|| try_close_staged(&data, &pairing, 6, None, &opts));
    eprintln!("(e) unreduced 4_1 size: estimate {:.0} MB -> {:?}; peak before refusal {:.1} MB", est as f64 / 1048576.0, r.as_ref().err(), peak as f64 / 1048576.0);
    assert!(r.is_err() && peak < (64 << 20));
}
