//! Milestone M1: scalable closure for the localized Soergel engine.
//!
//! The scanning (tangle level, N-independent) runs at a modest number of points.
//! For every indecomposable B_x in the final complex we choose R-bases of
//! Hom(R,B_x) and Hom(B_x,R) (graded Nakayama on multi-point data), interpolate the
//! Gram matrix G_x = h_b∘g_a and, for every final map f, the polynomial matrix F
//! with f∘g_a = Σ_c F_ca g'_c. The gl_N closure is then exact polynomial algebra:
//! τ_N(p) = coefficient of Π x_i^{N-1}. One scan serves every N.

use crate::input::crossing::Gen;
use crate::uncolored::engine::{monomials, mono_val, Cx, Obj, PointEngine};
use crate::algebra::field::Field;
use crate::combinat::hecke::Hecke;
use crate::algebra::laurent::Poincare;
use crate::algebra::linsys::Echelon;
use crate::algebra::sparse::SpVec;
use crate::uncolored::tower::{EndPolicy, TowerOpts};
use std::collections::{BTreeMap, HashMap};
use std::sync::Arc;

/// Exponent vector polynomial (dense map from exponent vectors to coefficients).
pub type EPoly<F> = HashMap<Vec<u32>, F>;

pub struct ScanData<F: Field> {
    pub objs: Vec<Obj>,
    /// (src index, tgt index, e-block matrix [tgt ebit][src ebit])
    pub edges: Vec<(usize, usize, Vec<Vec<F>>)>,
    /// per object type: (in-candidates, out-candidates) on the e-block
    pub gens: BTreeMap<usize, (Vec<(i32, Vec<F>)>, Vec<(i32, Vec<F>)>)>,
    pub max_objects: usize,
}

fn scan_point<F: Field>(n: usize, word: &[Gen], x: Vec<F>, hecke: Arc<Hecke>, start: usize) -> ScanData<F> {
    let mut pe = PointEngine::<F>::new(n, x, hecke);
    pe.ensure(start);
    let mut c = Cx::<F>::new();
    c.add_obj(Obj { x: start, q: 0, h: 0 });
    let mut maxo = 1;
    for g in word {
        c = pe.tensor_crossing(&c, *g);
        maxo = maxo.max(c.len());
        pe.simplify(&mut c);
    }
    let ids: Vec<usize> = c.objs.keys().copied().collect();
    let pos: HashMap<usize, usize> = ids.iter().enumerate().map(|(i, &id)| (id, i)).collect();
    let objs: Vec<Obj> = ids.iter().map(|id| c.objs[id]).collect();
    let mut gens = BTreeMap::new();
    let mut eb: HashMap<usize, Vec<u64>> = HashMap::new();
    for o in &objs {
        if !gens.contains_key(&o.x) {
            gens.insert(o.x, pe.closure_cands_e(o.x));
            pe.ensure(o.x);
            eb.insert(o.x, pe.ebits(o.x));
        }
    }
    let mut edges = vec![];
    for (&a, outs) in &c.out {
        for (&b, f) in outs {
            let xa = c.objs[&a].x;
            let xb = c.objs[&b].x;
            let m: Vec<Vec<F>> = eb[&xb].iter().map(|r| eb[&xa].iter().map(|cbit| f.cols[*cbit as usize].get(r).cloned().unwrap_or_else(F::zero)).collect()).collect();
            edges.push((pos[&a], pos[&b], m));
        }
    }
    ScanData { objs, edges, gens, max_objects: maxo }
}

fn gens_point<F: Field>(n: usize, x: Vec<F>, hecke: Arc<Hecke>, types: &[usize]) -> BTreeMap<usize, (Vec<(i32, Vec<F>)>, Vec<(i32, Vec<F>)>)> {
    let mut pe = PointEngine::<F>::new(n, x, hecke);
    types.iter().map(|&t| (t, pe.closure_cands_e(t))).collect()
}

fn par_map<T: Send + Sync, R: Send, Fun: Fn(&T) -> R + Sync>(items: &[T], f: Fun) -> Vec<R> {
    let nth = crate::num_threads();
    let chunks: Vec<Vec<usize>> = (0..nth).map(|t| (t..items.len()).step_by(nth).collect()).collect();
    let mut out: Vec<Option<R>> = (0..items.len()).map(|_| None).collect();
    let res: Vec<Vec<(usize, R)>> = std::thread::scope(|s| {
        let hs: Vec<_> = chunks.iter().map(|idx| { let f = &f; s.spawn(move || idx.iter().map(|&i| (i, f(&items[i]))).collect::<Vec<_>>()) }).collect();
        hs.into_iter().map(|h| h.join().unwrap()).collect()
    });
    for v in res {
        for (i, r) in v {
            out[i] = Some(r);
        }
    }
    out.into_iter().map(|x| x.unwrap()).collect()
}

/// Interpolation of homogeneous polynomials of weight k (degree 2k) from point values.
struct Interp<F: Field> {
    n: usize,
    cache: HashMap<usize, (Vec<Vec<u32>>, Vec<Vec<F>>)>,
}

impl<F: Field> Interp<F> {
    fn new(n: usize) -> Self {
        Interp { n, cache: HashMap::new() }
    }
    fn needed(&self, k: usize) -> usize {
        monomials(self.n, k).len()
    }
    /// Returns None if not enough points or values inconsistent with degree k.
    fn interpolate(&mut self, k: usize, pts: &[Vec<F>], vals: &[F]) -> Option<EPoly<F>> {
        let n = self.n;
        if !self.cache.contains_key(&k) {
            let ms = monomials(n, k);
            let m = ms.len();
            if pts.len() < m {
                return None;
            }
            let v: Vec<Vec<F>> = (0..m).map(|j| ms.iter().map(|mo| mono_val(mo, &pts[j])).collect()).collect();
            let inv = crate::soergel::local::inverse(&v)?;
            self.cache.insert(k, (ms, inv));
        }
        let (ms, inv) = &self.cache[&k];
        let m = ms.len();
        if pts.len() < m {
            return None;
        }
        let coef: Vec<F> = (0..m).map(|i| { let mut s = F::zero(); for j in 0..m { s.fma(&inv[i][j], &vals[j]); } s }).collect();
        // verify on remaining points
        for j in m..pts.len() {
            let mut s = F::zero();
            for (c, mo) in coef.iter().zip(ms.iter()) {
                if !c.is_zero() {
                    s.fma(c, &mono_val(mo, &pts[j]));
                }
            }
            if s != vals[j] {
                panic!("interpolation residual nonzero (degree bound violated) at weight {}", k);
            }
        }
        Some(ms.iter().cloned().zip(coef).filter(|(_, c)| !c.is_zero()).collect())
    }
}

fn padd<F: Field>(a: &mut EPoly<F>, b: &EPoly<F>, c: &F) {
    for (m, x) in b {
        let e = a.entry(m.clone()).or_insert_with(F::zero);
        e.fma(x, c);
    }
    a.retain(|_, v| !v.is_zero());
}
fn pmul<F: Field>(a: &EPoly<F>, b: &EPoly<F>) -> EPoly<F> {
    let mut out: EPoly<F> = HashMap::new();
    for (m1, x) in a {
        for (m2, y) in b {
            let m: Vec<u32> = m1.iter().zip(m2.iter()).map(|(p, q)| p + q).collect();
            out.entry(m).or_insert_with(F::zero).fma(x, y);
        }
    }
    out.retain(|_, v| !v.is_zero());
    out
}

/// Chosen R-bases for one object type.
pub struct TypeBasis<F: Field> {
    /// indices into in-candidates and their degrees (R-basis g_a of Hom(R,B_x))
    pub ins: Vec<(usize, i32)>,
    /// R-basis h_b of Hom(B_x,R)
    pub outs: Vec<(usize, i32)>,
    /// G[b][a] = h_b∘g_a polynomials
    pub gram: Vec<Vec<EPoly<F>>>,
}

/// Graded Nakayama selection: returns chosen candidate indices with degrees.
fn select_basis<F: Field>(cands: &[Vec<(i32, Vec<F>)>], pts: &[Vec<F>], target: &BTreeMap<i32, usize>, n: usize) -> Result<Vec<(usize, i32)>, String> {
    // cands[j][i] = (deg, vector) at point j
    let ncand = cands[0].len();
    let mut chosen: Vec<(usize, i32)> = vec![];
    let mut degs: Vec<i32> = cands[0].iter().map(|(d, _)| *d).collect();
    degs.sort();
    degs.dedup();
    let dimv = cands[0].first().map(|c| c.1.len()).unwrap_or(0);
    let len = pts.len() * dimv;
    for &d in &degs {
        let need = target.get(&d).copied().unwrap_or(0);
        let mut ech = Echelon::<F>::new(len.max(1));
        let mut nspan = 0;
        for &(ci, cd) in &chosen {
            if cd > d || (d - cd) % 2 != 0 {
                continue;
            }
            for mo in monomials(n, ((d - cd) / 2) as usize) {
                let mut v: SpVec<F> = vec![];
                for (j, p) in pts.iter().enumerate() {
                    let r = mono_val(&mo, p);
                    for (k, x) in cands[j][ci].1.iter().enumerate() {
                        let y = r.mul(x);
                        if !y.is_zero() {
                            v.push(((j * dimv + k) as u32, y));
                        }
                    }
                }
                if !ech.insert(&v) {
                    return Err(format!("not enough points for basis selection (degree {})", d));
                }
                nspan += 1;
            }
        }
        let _ = nspan;
        let mut got = 0;
        for i in 0..ncand {
            if cands[0][i].0 != d || got == need {
                continue;
            }
            let mut v: SpVec<F> = vec![];
            for j in 0..pts.len() {
                for (k, x) in cands[j][i].1.iter().enumerate() {
                    if !x.is_zero() {
                        v.push(((j * dimv + k) as u32, x.clone()));
                    }
                }
            }
            if ech.insert(&v) {
                chosen.push((i, d));
                got += 1;
            }
        }
        if got != need {
            return Err(format!("basis selection: degree {} found {} expected {}", d, got, need));
        }
    }
    Ok(chosen)
}

pub struct Engine2Result {
    pub poincare: BTreeMap<usize, Poincare>,
    pub points: usize,
    pub max_objects: usize,
    pub final_objects: usize,
    /// (scan, data preparation after the scan: bases + Gram + F, closure of all N)
    pub timings: (f64, f64, f64),
}

/// Which gl_N closure to run on the prepared data.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Closure {
    /// `close_exact`: one pairing matrix per (h, q) block of the whole closed complex (the oracle).
    Exact,
    /// Restriction tower (src/tower.rs), closing end variables in the given order.
    Tower(EndPolicy),
}

impl Closure {
    /// "exact", "tower" (right end), or "tower-left" / "tower-right" / "tower-alt".
    pub fn parse(s: &str) -> Option<Closure> {
        match s {
            "exact" => Some(Closure::Exact),
            "tower" => Some(Closure::Tower(EndPolicy::Right)),
            _ => s.strip_prefix("tower-").and_then(EndPolicy::parse).map(Closure::Tower),
        }
    }
}

/// The N-independent output of the scan: final complex, R-bases, Gram and edge polynomials.
pub struct ClosureData<F: Field> {
    pub n: usize,
    pub word: Vec<Gen>,
    pub objs: Vec<Obj>,
    /// (source, target) object indices of the final differential components
    pub edges: Vec<(usize, usize)>,
    /// fpolys[e][c][a]: f∘g_a = Σ_c F_ca g'_c
    pub fpolys: Vec<Vec<Vec<EPoly<F>>>>,
    pub bases: BTreeMap<usize, TypeBasis<F>>,
    pub points: usize,
    pub max_objects: usize,
    /// (scan, bases + Gram, F interpolation) seconds
    pub timings: (f64, f64, f64),
}

/// Truncation caps: N-1 per strand; the marked strand (bottom position 0) gets 0 in reduced mode.
pub fn closure_caps(n: usize, nn: usize, reduced: bool) -> Vec<u32> {
    assert!(nn >= 1);
    let mut top = vec![(nn - 1) as u32; n];
    if reduced {
        top[0] = 0;
    }
    top
}

/// Global q-shift of the closure (crossing shifts and the per-strand normalization).
pub fn closure_q_shift(n: usize, word: &[Gen], nn: usize, reduced: bool) -> i32 {
    let marked = if reduced { 1 } else { 0 };
    word.iter().map(|g| if g.positive { -((nn - 1) as i32) } else { (nn - 1) as i32 }).sum::<i32>() - ((n - marked) * (nn - 1)) as i32
}

/// gl_N homology (Poincaré polynomial) of the prepared data for one N.
pub fn close_with<F: Field>(data: &ClosureData<F>, nn: usize, reduced: bool, closure: Closure) -> Poincare {
    match closure {
        Closure::Exact => {
            let n = data.n;
            let top = closure_caps(n, nn, reduced);
            let global_q = closure_q_shift(n, &data.word, nn, reduced);
            close_exact(n, nn, &top, global_q, &data.objs, &data.edges, &data.fpolys, &data.bases)
        }
        Closure::Tower(policy) => close_tower_with(data, nn, reduced, TowerOpts::new(policy)).poincare,
    }
}

/// The tower closure with explicit options; returns per-stage statistics as well.
pub fn close_tower_with<F: Field>(data: &ClosureData<F>, nn: usize, reduced: bool, opts: TowerOpts) -> crate::uncolored::tower::TowerResult {
    let n = data.n;
    let top = closure_caps(n, nn, reduced);
    let global_q = closure_q_shift(n, &data.word, nn, reduced);
    let t = std::time::Instant::now();
    let r = crate::uncolored::tower::close_tower(n, &top, global_q, &data.objs, &data.edges, &data.fpolys, &data.bases, opts);
    let peak = r.stages.iter().map(|s| s.objects).max().unwrap_or(0);
    let alive: Vec<String> = r.stages.iter().map(|s| s.after_objects.to_string()).collect();
    let nonsplit: usize = r.stages.iter().map(|s| s.nonsplit).sum();
    eprintln!(
        "[tower N={}] {:?} end{}, objects after each stage {} (peak {} before elimination, {} non-split coords), {:.2}s",
        nn, opts.policy, if opts.unit_pivots_only { " (unit pivots only)" } else { "" }, alive.join(","), peak, nonsplit, t.elapsed().as_secs_f64()
    );
    r
}

/// Chain-level dimension dim Φ_N(B_x) of every type, computed by the tower decomposition.
pub fn tower_type_dims<F: Field>(data: &ClosureData<F>, nn: usize, reduced: bool, opts: TowerOpts) -> BTreeMap<usize, usize> {
    crate::uncolored::tower::type_dims(&closure_caps(data.n, nn, reduced), &data.bases, opts)
}

/// Chain-level dimension dim Φ_N(B_x) of every type, from the state spaces of `close_exact`.
pub fn exact_type_dims<F: Field>(data: &ClosureData<F>, nn: usize, reduced: bool) -> BTreeMap<usize, usize> {
    let top = closure_caps(data.n, nn, reduced);
    data.bases.iter().map(|(&x, tb)| (x, exact_state_space(data.n, &top, tb).0.values().map(|v| v.len()).sum())).collect()
}

/// Uncolored gl_N homology for every N in `ns` from one tangle-level scan (exact closure).
pub fn khr_engine2<F: Field>(n: usize, word: &[Gen], ns: &[usize], start_word: &[usize], seed: u64) -> Engine2Result {
    khr_engine2_mode::<F>(n, word, ns, start_word, seed, false)
}

/// As `khr_engine2`; with `reduced = true` the strand at bottom position 0 is
/// marked (reduced homology x_1^{N-1}·Φ_N, normalized so the reduced unknot is q^0).
pub fn khr_engine2_mode<F: Field>(n: usize, word: &[Gen], ns: &[usize], start_word: &[usize], seed: u64, reduced: bool) -> Engine2Result {
    khr_engine2_closure::<F>(n, word, ns, start_word, seed, reduced, Closure::Exact)
}

/// One scan and data preparation, then the chosen closure for every N in `ns` (sequentially).
pub fn khr_engine2_closure<F: Field>(n: usize, word: &[Gen], ns: &[usize], start_word: &[usize], seed: u64, reduced: bool, closure: Closure) -> Engine2Result {
    let data = prepare_closure_data::<F>(n, word, start_word, seed);
    let mut results = BTreeMap::new();
    let tc = std::time::Instant::now();
    for &nn in ns {
        let t1 = std::time::Instant::now();
        let p = close_with(&data, nn, reduced, closure);
        eprintln!("[engine2] N={} closed in {:.2}s, rank {}", nn, t1.elapsed().as_secs_f64(), p.total_rank());
        results.insert(nn, p);
    }
    Engine2Result {
        poincare: results,
        points: data.points,
        max_objects: data.max_objects,
        final_objects: data.objs.len(),
        timings: (data.timings.0, data.timings.1 + data.timings.2, tc.elapsed().as_secs_f64()),
    }
}

/// Scan at enough points, choose R-bases (graded Nakayama), interpolate Gram and edge
/// polynomials. N-independent; the scan data is released before returning.
pub fn prepare_closure_data<F: Field>(n: usize, word: &[Gen], start_word: &[usize], seed: u64) -> ClosureData<F> {
    use rand::{rngs::StdRng, Rng, SeedableRng};
    let t0 = std::time::Instant::now();
    let hecke = Arc::new(Hecke::new(n));
    let start = hecke.g.from_word(start_word);
    let mut rng = StdRng::seed_from_u64(seed);
    let mut newpt = || -> Vec<F> { (0..n).map(|_| F::from_i64(rng.gen_range(1..2_000_000_000))).collect() };
    let mut m = 24usize.max(monomials(n, 2).len() + 4);
    let mut pts: Vec<Vec<F>> = (0..m).map(|_| newpt()).collect();
    let mut scans: Vec<ScanData<F>> = par_map(&pts, |p| scan_point::<F>(n, word, p.clone(), hecke.clone(), start));
    let objs = scans[0].objs.clone();
    let max_objects = scans[0].max_objects;
    for sd in &scans {
        assert_eq!(sd.objs, objs, "final complexes differ between points");
    }
    let t_scan = t0.elapsed().as_secs_f64();
    let types: Vec<usize> = scans[0].gens.keys().copied().collect();
    // graded ranks of Hom(R,B_x) = h_{e,x}(v)
    let target = |x: usize| -> BTreeMap<i32, usize> { hecke.kl[x][0].c.iter().map(|(&e, &c)| (e, c as usize)).collect() };
    // ---- choose bases (may need more gens-only points) ----
    // Every type x uses the first np[x] points of gpts. A type gets more points only when
    // its own basis selection or Gram interpolation needs them (dense interpolation of
    // weight k needs C(k+n-1, n-1) points, plus 8 for the residual check), so the many
    // low-degree types are not evaluated at the hundreds of points the top types need.
    let mut gpts: Vec<Vec<F>> = pts.clone();
    let mut gdata: Vec<BTreeMap<usize, (Vec<(i32, Vec<F>)>, Vec<(i32, Vec<F>)>)>> = scans.iter().map(|s| s.gens.clone()).collect();
    let mut np: BTreeMap<usize, usize> = types.iter().map(|&x| (x, gpts.len())).collect();
    let t_bases_start = t0.elapsed().as_secs_f64();
    let mut bases: BTreeMap<usize, TypeBasis<F>> = BTreeMap::new();
    let mut interp = Interp::<F>::new(n);
    loop {
        let mut want: BTreeMap<usize, usize> = BTreeMap::new();
        for &x in &types {
            if bases.contains_key(&x) {
                continue;
            }
            let npx = np[&x];
            let px = &gpts[..npx];
            let incs: Vec<Vec<(i32, Vec<F>)>> = gdata[..npx].iter().map(|g| g[&x].0.clone()).collect();
            let outcs: Vec<Vec<(i32, Vec<F>)>> = gdata[..npx].iter().map(|g| g[&x].1.clone()).collect();
            let (ins, outs) = match (select_basis(&incs, px, &target(x), n), select_basis(&outcs, px, &target(x), n)) {
                (Ok(i), Ok(o)) => (i, o),
                _ => {
                    want.insert(x, 2 * npx);
                    continue;
                }
            };
            let mut need = 0;
            for &(_, ad) in &ins {
                for &(_, bd) in &outs {
                    need = need.max(interp.needed(((ad + bd) / 2) as usize));
                }
            }
            if need > npx {
                want.insert(x, need + 8);
                continue;
            }
            // Gram polynomials
            let mut gram = vec![];
            for &(bi, bd) in &outs {
                let mut row = vec![];
                for &(ai, ad) in &ins {
                    let k = ((ad + bd) / 2) as usize;
                    let vals: Vec<F> = (0..npx)
                        .map(|j| {
                            let a = &incs[j][ai].1;
                            let b = &outcs[j][bi].1;
                            let mut s = F::zero();
                            for (u, v) in a.iter().zip(b.iter()) { s.fma(u, v); }
                            s
                        })
                        .collect();
                    row.push(interp.interpolate(k, px, &vals).unwrap());
                }
                gram.push(row);
            }
            bases.insert(x, TypeBasis { ins, outs, gram });
        }
        if want.is_empty() {
            break;
        }
        // add gens-only points, evaluating at each only the types that want it
        let total = *want.values().max().unwrap();
        while gpts.len() < total {
            gpts.push(newpt());
            gdata.push(BTreeMap::new());
        }
        let lo = want.keys().map(|x| np[x]).min().unwrap();
        let jobs: Vec<(usize, Vec<usize>)> = (lo..total)
            .map(|j| (j, want.iter().filter(|&(x, &w)| np[x] <= j && j < w).map(|(&x, _)| x).collect::<Vec<usize>>()))
            .filter(|(_, ts)| !ts.is_empty())
            .collect();
        let more = par_map(&jobs, |(j, ts)| gens_point::<F>(n, gpts[*j].clone(), hecke.clone(), ts));
        for ((j, _), m) in jobs.iter().zip(more) {
            gdata[*j].extend(m);
        }
        for (x, w) in want {
            np.insert(x, w);
        }
    }
    let t_bases = t0.elapsed().as_secs_f64() - t_bases_start;
    eprintln!("[engine2] scan {:.2}s ({} points), bases+Gram {:.2}s ({} gens points)", t_scan, pts.len(), t_bases, gpts.len());
    drop(gdata);
    // ---- F matrices for each edge (may need more scan points) ----
    let deg_f = |e: &(usize, usize, Vec<Vec<F>>)| -> i32 { objs[e.0].q - objs[e.1].q };
    let mut maxk = 0usize;
    for e in &scans[0].edges {
        let bx = &bases[&objs[e.0].x];
        let by = &bases[&objs[e.1].x];
        for &(_, ad) in &bx.ins {
            for &(_, cd) in &by.ins {
                let dd = deg_f(e) + ad - cd;
                if dd >= 0 {
                    maxk = maxk.max((dd / 2) as usize);
                }
            }
        }
    }
    let need = monomials(n, maxk).len() + 4;
    let t_extra_start = t0.elapsed().as_secs_f64();
    let n_before_extra = pts.len();
    if need > pts.len() {
        let extra: Vec<Vec<F>> = (pts.len()..need).map(|_| newpt()).collect();
        let more = par_map(&extra, |p| scan_point::<F>(n, word, p.clone(), hecke.clone(), start));
        for sd in &more {
            assert_eq!(sd.objs, objs, "final complexes differ between points");
        }
        pts.extend(extra);
        scans.extend(more);
        m = pts.len();
    }
    let _ = m;
    // Extra scan points needed only for the F interpolation count as scan time.
    let t_extra = t0.elapsed().as_secs_f64() - t_extra_start;
    if pts.len() > n_before_extra {
        eprintln!("[engine2] extra scan for F: {} more points in {:.2}s", pts.len() - n_before_extra, t_extra);
    }
    let t_scan = t_scan + t_extra;
    let t_interp_start = t0.elapsed().as_secs_f64();
    // F polynomials: fmat[edge][c][a]
    let mut fpolys: Vec<Vec<Vec<EPoly<F>>>> = vec![];
    let mut interp2 = Interp::<F>::new(n);
    for (ei, e) in scans[0].edges.iter().enumerate() {
        let (xs, xt) = (objs[e.0].x, objs[e.1].x);
        let bx = &bases[&xs];
        let by = &bases[&xt];
        let mc = by.ins.len();
        let ma = bx.ins.len();
        // values per point
        let mut vals: Vec<Vec<Vec<F>>> = vec![vec![vec![]; ma]; mc];
        for (j, sd) in scans.iter().enumerate() {
            let fm = &sd.edges[ei].2;
            let gx = &sd.gens[&xs].0;
            let gy = &sd.gens[&xt].0;
            // basis vectors of y at this point
            let ybas: Vec<Vec<F>> = by.ins.iter().map(|(ci, _)| gy[*ci].1.clone()).collect();
            let rows: Vec<Vec<F>> = (0..ybas[0].len()).map(|r| (0..mc).map(|c| ybas[c][r].clone()).collect()).collect();
            for (a, (ai, _)) in bx.ins.iter().enumerate() {
                let g = &gx[*ai].1;
                let img: Vec<F> = fm.iter().map(|row| { let mut s = F::zero(); for (u, v) in row.iter().zip(g.iter()) { s.fma(u, v); } s }).collect();
                let sol = crate::algebra::dense::solve(&rows, mc, &img).expect("f∘g not in span of basis (point)");
                for c in 0..mc {
                    vals[c][a].push(sol[c].clone());
                }
            }
            let _ = j;
        }
        let mut fm = vec![];
        for c in 0..mc {
            let mut row = vec![];
            for a in 0..ma {
                let dd = deg_f(e) + bx.ins[a].1 - by.ins[c].1;
                if dd < 0 {
                    assert!(vals[c][a].iter().all(|v| v.is_zero()), "negative-degree component nonzero");
                    row.push(HashMap::new());
                } else {
                    row.push(interp2.interpolate((dd / 2) as usize, &pts, &vals[c][a]).expect("not enough points for F"));
                }
            }
            fm.push(row);
        }
        fpolys.push(fm);
    }
    let t_interp = t0.elapsed().as_secs_f64() - t_interp_start;
    eprintln!("[engine2] F interpolation {:.2}s", t_interp);
    // the per-point scan data is no longer needed: keep only the edge list
    let edges: Vec<(usize, usize)> = scans[0].edges.iter().map(|e| (e.0, e.1)).collect();
    drop(scans);
    ClosureData { n, word: word.to_vec(), objs, edges, fpolys, bases, points: pts.len(), max_objects, timings: (t_scan, t_bases, t_interp) }
}

/// Write prepared closure data to a text file (prime fields only), so that expensive
/// preparations can be closed again later without rescanning (`load_closure_data`).
pub fn save_closure_data<F: Field>(data: &ClosureData<F>, path: &str) -> std::io::Result<()> {
    use std::fmt::Write as _;
    assert!(F::characteristic() > 0, "save_closure_data supports prime fields only");
    let n = data.n;
    let mut s = String::new();
    let poly = |s: &mut String, p: &EPoly<F>| {
        let mut terms: Vec<(&Vec<u32>, &F)> = p.iter().collect();
        terms.sort_by(|a, b| a.0.cmp(b.0));
        write!(s, "{}", terms.len()).unwrap();
        for (e, c) in terms {
            assert_eq!(e.len(), n);
            for x in e {
                write!(s, " {}", x).unwrap();
            }
            write!(s, " {}", c).unwrap();
        }
        s.push('\n');
    };
    writeln!(s, "krhom-closure-data 1 {} {}", F::characteristic(), n).unwrap();
    writeln!(s, "{}", data.word.iter().map(|g| format!("{}", if g.positive { g.i as i64 + 1 } else { -(g.i as i64 + 1) })).collect::<Vec<_>>().join(" ")).unwrap();
    writeln!(s, "{} {} {} {} {}", data.points, data.max_objects, data.timings.0, data.timings.1, data.timings.2).unwrap();
    writeln!(s, "{}", data.objs.len()).unwrap();
    for o in &data.objs {
        writeln!(s, "{} {} {}", o.x, o.q, o.h).unwrap();
    }
    writeln!(s, "{}", data.edges.len()).unwrap();
    for (ei, &(a, b)) in data.edges.iter().enumerate() {
        let fm = &data.fpolys[ei];
        writeln!(s, "{} {} {} {}", a, b, fm.len(), fm.first().map(|r| r.len()).unwrap_or(0)).unwrap();
        for row in fm {
            for p in row {
                poly(&mut s, p);
            }
        }
    }
    writeln!(s, "{}", data.bases.len()).unwrap();
    for (x, tb) in &data.bases {
        let pairs = |v: &[(usize, i32)]| v.iter().map(|(i, d)| format!("{} {}", i, d)).collect::<Vec<_>>().join(" ");
        writeln!(s, "{} {} {}", x, tb.ins.len(), tb.outs.len()).unwrap();
        writeln!(s, "{}", pairs(&tb.ins)).unwrap();
        writeln!(s, "{}", pairs(&tb.outs)).unwrap();
        for row in &tb.gram {
            for p in row {
                poly(&mut s, p);
            }
        }
    }
    std::fs::write(path, s)
}

/// Read closure data written by `save_closure_data` (the same prime field).
pub fn load_closure_data<F: Field>(path: &str) -> Result<ClosureData<F>, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {}", path, e))?;
    let mut lines = text.lines();
    let mut next = || lines.next().ok_or_else(|| "unexpected end of file".to_string());
    let nums = |l: &str| -> Result<Vec<i64>, String> { l.split_whitespace().map(|t| t.parse::<i64>().map_err(|e| format!("{:?}: {}", t, e))).collect() };
    let head: Vec<String> = next()?.split_whitespace().map(|t| t.to_string()).collect();
    if head.len() != 4 || head[0] != "krhom-closure-data" || head[1] != "1" {
        return Err("not a krhom closure data file (version 1)".into());
    }
    if head[2].parse::<u64>().ok() != Some(F::characteristic()) {
        return Err(format!("field mismatch: file has characteristic {}, expected {}", head[2], F::characteristic()));
    }
    let n: usize = head[3].parse().map_err(|_| "bad strand count")?;
    let word: Vec<Gen> = nums(next()?)?.into_iter().map(|k| Gen { i: (k.unsigned_abs() - 1) as usize, positive: k > 0 }).collect();
    let meta: Vec<String> = next()?.split_whitespace().map(|t| t.to_string()).collect();
    if meta.len() != 5 {
        return Err("bad metadata line".into());
    }
    let points: usize = meta[0].parse().map_err(|_| "bad points")?;
    let max_objects: usize = meta[1].parse().map_err(|_| "bad max_objects")?;
    let tm = |k: usize| meta[k].parse::<f64>().map_err(|_| "bad timing".to_string());
    let timings = (tm(2)?, tm(3)?, tm(4)?);
    let parse_poly = |l: &str| -> Result<EPoly<F>, String> {
        let v = nums(l)?;
        let k = *v.first().ok_or("empty polynomial line")? as usize;
        if v.len() != 1 + k * (n + 1) {
            return Err(format!("bad polynomial line ({} numbers for {} terms)", v.len(), k));
        }
        let mut p: EPoly<F> = HashMap::new();
        for t in 0..k {
            let base = 1 + t * (n + 1);
            let e: Vec<u32> = v[base..base + n].iter().map(|&x| x as u32).collect();
            p.insert(e, F::from_i64(v[base + n]));
        }
        Ok(p)
    };
    let no = nums(next()?)?[0] as usize;
    let mut objs = Vec::with_capacity(no);
    for _ in 0..no {
        let v = nums(next()?)?;
        objs.push(Obj { x: v[0] as usize, q: v[1] as i32, h: v[2] as i32 });
    }
    let ne = nums(next()?)?[0] as usize;
    let (mut edges, mut fpolys) = (Vec::with_capacity(ne), Vec::with_capacity(ne));
    for _ in 0..ne {
        let v = nums(next()?)?;
        edges.push((v[0] as usize, v[1] as usize));
        let (r, c) = (v[2] as usize, v[3] as usize);
        let mut fm = Vec::with_capacity(r);
        for _ in 0..r {
            let mut row = Vec::with_capacity(c);
            for _ in 0..c {
                row.push(parse_poly(next()?)?);
            }
            fm.push(row);
        }
        fpolys.push(fm);
    }
    let nb = nums(next()?)?[0] as usize;
    let mut bases = BTreeMap::new();
    for _ in 0..nb {
        let v = nums(next()?)?;
        let (x, ni, nout) = (v[0] as usize, v[1] as usize, v[2] as usize);
        let pairs = |l: &str| -> Result<Vec<(usize, i32)>, String> { Ok(nums(l)?.chunks(2).map(|c| (c[0] as usize, c[1] as i32)).collect()) };
        let ins = pairs(next()?)?;
        let outs = pairs(next()?)?;
        if ins.len() != ni || outs.len() != nout {
            return Err(format!("type {}: bad basis lines", x));
        }
        let mut gram = Vec::with_capacity(nout);
        for _ in 0..nout {
            let mut row = Vec::with_capacity(ni);
            for _ in 0..ni {
                row.push(parse_poly(next()?)?);
            }
            gram.push(row);
        }
        bases.insert(x, TypeBasis { ins, outs, gram });
    }
    Ok(ClosureData { n, word, objs, edges, fpolys, bases, points, max_objects, timings })
}

/// State space of one type for close_exact: per degree d, representatives (a, m1) and
/// dual representatives (b, m2) of Φ_N(B_x) (rank of the pairing matrix).
type StateSpace = (BTreeMap<i32, Vec<(usize, Vec<u32>)>>, BTreeMap<i32, Vec<(usize, Vec<u32>)>>);

fn exact_state_space<F: Field>(n: usize, top: &[u32], tb: &TypeBasis<F>) -> StateSpace {
    let topdeg = 2 * top.iter().sum::<u32>() as i32;
    let coeff = |p: &EPoly<F>, m1: &[u32], m2: &[u32]| -> F {
        let mut e = Vec::with_capacity(n);
        for i in 0..n {
            let s = m1[i] + m2[i];
            if s > top[i] {
                return F::zero();
            }
            e.push(top[i] - s);
        }
        p.get(&e).cloned().unwrap_or_else(F::zero)
    };
    let mut reps = BTreeMap::new();
    let mut duals = BTreeMap::new();
    // A representative of degree d needs d >= min(ins); its dual partner has degree
    // topdeg - d, which needs topdeg - d >= min(outs). So d runs over [min(ins), topdeg - min(outs)].
    // (The old upper bound used max(outs) and silently dropped the top degrees when
    // min(outs) < min(ins); braid data never has that, synthetic data can.)
    let gmin = tb.ins.iter().map(|g| g.1).min().unwrap_or(0);
    let omin = tb.outs.iter().map(|g| g.1).min().unwrap_or(0);
    for d in gmin..=(topdeg - omin) {
        let mut span = vec![];
        for (a, &(_, ad)) in tb.ins.iter().enumerate() {
            let rd = d - ad;
            if rd < 0 || rd % 2 != 0 { continue; }
            for mo in monomials(n, (rd / 2) as usize) {
                if mo.iter().zip(top.iter()).all(|(&e, &t)| e <= t) { span.push((a, mo)); }
            }
        }
        let dd = topdeg - d;
        let mut dspan = vec![];
        for (b, &(_, bd)) in tb.outs.iter().enumerate() {
            let rd = dd - bd;
            if rd < 0 || rd % 2 != 0 { continue; }
            for mo in monomials(n, (rd / 2) as usize) {
                if mo.iter().zip(top.iter()).all(|(&e, &t)| e <= t) { dspan.push((b, mo)); }
            }
        }
        if span.is_empty() || dspan.is_empty() { continue; }
        let gram: Vec<Vec<F>> = span.iter().map(|(a, m1)| dspan.iter().map(|(b, m2)| coeff(&tb.gram[*b][*a], m1, m2)).collect()).collect();
        let rows = indep_rows(&gram);
        if rows.is_empty() { continue; }
        let sub: Vec<Vec<F>> = rows.iter().map(|&k| gram[k].clone()).collect();
        let cols = indep_rows(&transpose(&sub));
        reps.insert(d, rows.iter().map(|&k| span[k].clone()).collect::<Vec<_>>());
        duals.insert(d, cols.iter().map(|&c| dspan[c].clone()).collect::<Vec<_>>());
    }
    (reps, duals)
}

/// Exact closure. `top[i]` is the truncation cap of variable x_i (N-1 for an
/// unreduced strand, 0 for the marked strand in reduced mode); τ = coefficient of x^top.
fn close_exact<F: Field>(
    n: usize,
    nn: usize,
    top: &[u32],
    global_q: i32,
    objs: &[Obj],
    edges: &[(usize, usize)],
    fpolys: &[Vec<Vec<EPoly<F>>>],
    bases: &BTreeMap<usize, TypeBasis<F>>,
) -> Poincare {
    let top: Vec<u32> = top.to_vec();
    // state spaces per object type and degree: reps (a, mono), duals (b, mono).
    // (rank of d = rank of the pairing matrix, so no inverse Gram is needed)
    struct SS {
        reps: BTreeMap<i32, Vec<(usize, Vec<u32>)>>,
        duals: BTreeMap<i32, Vec<(usize, Vec<u32>)>>,
    }
    let mut ss: BTreeMap<usize, SS> = BTreeMap::new();
    let tss = std::time::Instant::now();
    let mut ssdim = 0usize;
    for (&x, tb) in bases {
        let (reps, duals) = exact_state_space(n, &top, tb);
        ssdim += reps.values().map(|v: &Vec<(usize, Vec<u32>)>| v.len()).sum::<usize>();
        ss.insert(x, SS { reps, duals });
    }
    eprintln!("[close N={}] state spaces: {} types, total dim {}, {:.2}s", nn, bases.len(), ssdim, tss.elapsed().as_secs_f64());
    let trk = std::time::Instant::now();
    let mut groups: BTreeMap<(i32, i32), Vec<(usize, i32, usize)>> = BTreeMap::new();
    for (i, o) in objs.iter().enumerate() {
        for (d, reps) in &ss[&o.x].reps {
            for k in 0..reps.len() {
                groups.entry((o.h, d + o.q + global_q)).or_default().push((i, *d, k));
            }
        }
    }
    // precompute P_ab = Σ_c F_ca G_y[b][c] per edge
    let mut out_edges: HashMap<usize, Vec<usize>> = HashMap::new();
    for (ei, &(s, _)) in edges.iter().enumerate() {
        out_edges.entry(s).or_default().push(ei);
    }
    let mut pab: Vec<Vec<Vec<EPoly<F>>>> = vec![];
    for (ei, &(_, t)) in edges.iter().enumerate() {
        let by = &bases[&objs[t].x];
        let ma = fpolys[ei].first().map(|r| r.len()).unwrap_or(0);
        let mut m = vec![];
        for a in 0..ma {
            let mut row = vec![];
            for b in 0..by.outs.len() {
                let mut acc: EPoly<F> = HashMap::new();
                for c in 0..by.ins.len() {
                    let prod = pmul(&fpolys[ei][c][a], &by.gram[b][c]);
                    padd(&mut acc, &prod, &F::one());
                }
                row.push(acc);
            }
            m.push(row);
        }
        pab.push(m);
    }
    // rank(d) = rank of the pairing matrix ⟨d(rep), dual⟩ (the target Gram is invertible),
    // so columns are indexed by target duals and no inverse Gram is needed.
    // Packed monomials in a u128 with `w` bits per variable. A row term is the sum of a
    // representative monomial m1 (exponents ≤ top) and a pruned term of P_ab (exponents
    // ≤ top), so every field of the sum is ≤ 2·max(top) < 2^w and never carries into the
    // next field; `fits` then tests each field separately and `toppack - sum` cannot borrow.
    let w = exact_pack_width(&top);
    assert!(n * w <= 128, "close_exact: {} strands with caps up to {} do not fit packed monomials ({} bits per variable)", n, top.iter().max().unwrap_or(&0), w);
    let fmask: u128 = (1u128 << w) - 1;
    let pack = |m: &[u32]| -> u128 { m.iter().enumerate().fold(0u128, |acc, (i, &e)| acc | ((e as u128) << (w * i))) };
    let toppack = pack(&top);
    let fits = |s: u128| -> bool { (0..n).all(|i| ((s >> (w * i)) & fmask) <= top[i] as u128) };
    // columns (duals) sorted by a lex term order on m2 (ascending) so that the
    // leftmost entry of a row m1·F is its leading term m1·lt(F): little fill-in.
    let lexkey = |m: &[u32]| -> u128 { m.iter().enumerate().fold(0u128, |acc, (i, &e)| acc | ((e as u128) << (w * (n - 1 - i)))) };
    let mut dual_lists: BTreeMap<(i32, i32), Vec<(u128, usize, usize, u128)>> = BTreeMap::new();
    for (i, o) in objs.iter().enumerate() {
        for (d, duals) in &ss[&o.x].duals {
            let key = (o.h, d + o.q + global_q);
            for (b, m2) in duals {
                dual_lists.entry(key).or_default().push((lexkey(m2), i, *b, pack(m2)));
            }
        }
    }
    let mut dual_groups: BTreeMap<(i32, i32), usize> = BTreeMap::new();
    let mut dual_index: HashMap<(u32, u32, u128), usize> = HashMap::new();
    for (key, mut list) in dual_lists {
        list.sort();
        for (idx, (_, i, b, m2p)) in list.iter().enumerate() {
            dual_index.insert((*i as u32, *b as u32, *m2p), idx);
        }
        dual_groups.insert(key, list.len());
    }
    // pruned packed P_ab terms
    let pab_packed: Vec<Vec<Vec<Vec<(u128, F)>>>> = pab
        .iter()
        .map(|m| m.iter().map(|row| row.iter().map(|p| p.iter().filter(|(t, _)| t.iter().zip(top.iter()).all(|(&e, &c)| e <= c)).map(|(t, c)| (pack(t), c.clone())).collect()).collect()).collect())
        .collect();
    let blocks: Vec<((i32, i32), Vec<(usize, i32, usize)>)> = groups.iter().map(|(k, v)| (*k, v.clone())).collect();
    let ranks: Vec<((i32, i32), usize)> = par_map(&blocks, |((h, q), elems)| {
        let tl = dual_groups.get(&(h + 1, *q)).copied().unwrap_or(0);
        if tl == 0 { return ((*h, *q), 0); }
        let mut ech = crate::algebra::linsys::RankEchelon::<F>::new(tl);
        for &(i, d, k) in elems {
            let (a, m1) = &ss[&objs[i].x].reps[&d][k];
            let m1p = pack(m1);
            let mut v: HashMap<usize, F> = HashMap::new();
            for &ei in out_edges.get(&i).map(|v| v.as_slice()).unwrap_or(&[]) {
                let t = edges[ei].1;
                let td = d + objs[i].q - objs[t].q;
                if !ss[&objs[t].x].duals.contains_key(&td) { continue; }
                for (b, terms) in pab_packed[ei][*a].iter().enumerate() {
                    for (tp, c) in terms {
                        let sum = tp + m1p;
                        if !fits(sum) { continue; }
                        let m2 = toppack - sum;
                        if let Some(&col) = dual_index.get(&(t as u32, b as u32, m2)) {
                            v.entry(col).or_insert_with(F::zero).add_assign(c);
                        }
                    }
                }
            }
            let mut sv: SpVec<F> = v.into_iter().filter(|(_, x)| !x.is_zero()).map(|(i, x)| (i as u32, x)).collect();
            sv.sort_by_key(|e| e.0);
            ech.insert(&sv);
        }
        ((*h, *q), ech.rank())
    });
    let rank: HashMap<(i32, i32), usize> = ranks.into_iter().collect();
    eprintln!("[close N={}] chain dim {}, ranks {:.2}s", nn, groups.values().map(|v| v.len()).sum::<usize>(), trk.elapsed().as_secs_f64());
    let mut pc = Poincare::default();
    for ((h, q), elems) in &groups {
        let r_out = rank.get(&(*h, *q)).copied().unwrap_or(0);
        let r_in = rank.get(&(h - 1, *q)).copied().unwrap_or(0);
        pc.add_term(*h, *q, elems.len() as i128 - r_out as i128 - r_in as i128);
    }
    pc
}

/// Bits per variable of the packed monomials of `close_exact`: enough for 2·max(top), the
/// largest exponent of a product of two monomials within the caps (at least 5 bits).
pub fn exact_pack_width(top: &[u32]) -> usize {
    let m = 2 * top.iter().copied().max().unwrap_or(0) as u64;
    ((64 - m.leading_zeros()) as usize).max(5)
}

fn transpose<F: Field>(m: &[Vec<F>]) -> Vec<Vec<F>> {
    if m.is_empty() { return vec![]; }
    (0..m[0].len()).map(|j| m.iter().map(|r| r[j].clone()).collect()).collect()
}

fn indep_rows<F: Field>(m: &[Vec<F>]) -> Vec<usize> {
    if m.is_empty() { return vec![]; }
    let mut ech = Echelon::<F>::new(m[0].len().max(1));
    let mut out = vec![];
    for (i, r) in m.iter().enumerate() {
        let v: SpVec<F> = r.iter().enumerate().filter(|(_, x)| !x.is_zero()).map(|(j, x)| (j as u32, x.clone())).collect();
        if ech.insert(&v) { out.push(i); }
    }
    out
}
