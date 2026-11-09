//! Milestone M1: scalable closure for the localized Soergel engine.
//!
//! The scanning (tangle level, N-independent) runs at a modest number of points.
//! For every indecomposable B_x in the final complex we choose R-bases of
//! Hom(R,B_x) and Hom(B_x,R) (graded Nakayama on multi-point data), interpolate the
//! Gram matrix G_x = h_b∘g_a and, for every final map f, the polynomial matrix F
//! with f∘g_a = Σ_c F_ca g'_c. The gl_N closure is then exact polynomial algebra:
//! τ_N(p) = coefficient of Π x_i^{N-1}. One scan serves every N.

use crate::crossing::Gen;
use crate::engine::{monomials, mono_val, Cx, Obj, PointEngine};
use crate::field::Field;
use crate::hecke::Hecke;
use crate::laurent::Poincare;
use crate::linsys::Echelon;
use crate::sparse::SpVec;
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
    let nth = std::thread::available_parallelism().map(|x| x.get()).unwrap_or(4);
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
            let inv = crate::local::inverse(&v)?;
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
struct TypeBasis<F: Field> {
    /// indices into in-candidates and their degrees
    ins: Vec<(usize, i32)>,
    outs: Vec<(usize, i32)>,
    /// G[b][a] polynomials
    gram: Vec<Vec<EPoly<F>>>,
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
    pub timings: (f64, f64, f64),
}

/// Uncolored gl_N homology for every N in `ns` from one tangle-level scan.
pub fn khr_engine2<F: Field>(n: usize, word: &[Gen], ns: &[usize], start_word: &[usize], seed: u64) -> Engine2Result {
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
    let mut gpts: Vec<Vec<F>> = pts.clone();
    let mut gdata: Vec<BTreeMap<usize, (Vec<(i32, Vec<F>)>, Vec<(i32, Vec<F>)>)>> = scans.iter().map(|s| s.gens.clone()).collect();
    let mut bases: BTreeMap<usize, TypeBasis<F>> = BTreeMap::new();
    let mut interp = Interp::<F>::new(n);
    loop {
        let mut ok = true;
        bases.clear();
        let mut need_pts = gpts.len();
        for &x in &types {
            let incs: Vec<Vec<(i32, Vec<F>)>> = gdata.iter().map(|g| g[&x].0.clone()).collect();
            let outcs: Vec<Vec<(i32, Vec<F>)>> = gdata.iter().map(|g| g[&x].1.clone()).collect();
            let ins = match select_basis(&incs, &gpts, &target(x), n) { Ok(b) => b, Err(_) => { ok = false; need_pts = need_pts * 2; break; } };
            let outs = match select_basis(&outcs, &gpts, &target(x), n) { Ok(b) => b, Err(_) => { ok = false; need_pts = need_pts * 2; break; } };
            // Gram polynomials
            let mut gram = vec![];
            for &(bi, bd) in &outs {
                let mut row = vec![];
                for &(ai, ad) in &ins {
                    let k = ((ad + bd) / 2) as usize;
                    if interp.needed(k) > gpts.len() {
                        ok = false;
                        need_pts = need_pts.max(interp.needed(k) + 8);
                        row.push(HashMap::new());
                        continue;
                    }
                    let vals: Vec<F> = (0..gpts.len())
                        .map(|j| {
                            let a = &incs[j][ai].1;
                            let b = &outcs[j][bi].1;
                            let mut s = F::zero();
                            for (u, v) in a.iter().zip(b.iter()) { s.fma(u, v); }
                            s
                        })
                        .collect();
                    row.push(interp.interpolate(k, &gpts, &vals).unwrap());
                }
                gram.push(row);
            }
            bases.insert(x, TypeBasis { ins, outs, gram });
        }
        if ok {
            break;
        }
        // add gens-only points
        let extra: Vec<Vec<F>> = (gpts.len()..need_pts).map(|_| newpt()).collect();
        let more = par_map(&extra, |p| gens_point::<F>(n, p.clone(), hecke.clone(), &types));
        gpts.extend(extra);
        gdata.extend(more);
        interp = Interp::<F>::new(n);
    }
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
                let sol = crate::dense::solve(&rows, mc, &img).expect("f∘g not in span of basis (point)");
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
    // ---- exact closure for each N ----
    let mut results = BTreeMap::new();
    let tc = std::time::Instant::now();
    for &nn in ns {
        let t1 = std::time::Instant::now();
        let p = close_exact(n, nn, word, &objs, &scans[0].edges.iter().map(|e| (e.0, e.1)).collect::<Vec<_>>(), &fpolys, &bases);
        eprintln!("[engine2] N={} closed in {:.2}s, rank {}", nn, t1.elapsed().as_secs_f64(), p.total_rank());
        results.insert(nn, p);
    }
    Engine2Result {
        poincare: results,
        points: pts.len(),
        max_objects,
        final_objects: objs.len(),
        timings: (t_scan, t_interp, tc.elapsed().as_secs_f64()),
    }
}

fn close_exact<F: Field>(
    n: usize,
    nn: usize,
    word: &[Gen],
    objs: &[Obj],
    edges: &[(usize, usize)],
    fpolys: &[Vec<Vec<EPoly<F>>>],
    bases: &BTreeMap<usize, TypeBasis<F>>,
) -> Poincare {
    let top = vec![(nn - 1) as u32; n];
    let topdeg = 2 * (n * (nn - 1)) as i32;
    let global_q: i32 = word.iter().map(|g| if g.positive { -((nn - 1) as i32) } else { (nn - 1) as i32 }).sum::<i32>() - (n * (nn - 1)) as i32;
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
    // state spaces per object type and degree: reps (a, mono), duals (b, mono), ginv
    struct SS<F: Field> {
        reps: BTreeMap<i32, Vec<(usize, Vec<u32>)>>,
        duals: BTreeMap<i32, Vec<(usize, Vec<u32>)>>,
        ginv: BTreeMap<i32, Vec<Vec<F>>>,
    }
    let mut ss: BTreeMap<usize, SS<F>> = BTreeMap::new();
    let tss = std::time::Instant::now();
    let mut ssdim = 0usize;
    for (&x, tb) in bases {
        let mut reps = BTreeMap::new();
        let mut duals = BTreeMap::new();
        let mut ginvs = BTreeMap::new();
        let gmin = tb.ins.iter().map(|g| g.1).min().unwrap_or(0);
        let gmax = tb.outs.iter().map(|g| g.1).max().unwrap_or(0);
        for d in gmin..=(topdeg - gmin.min(gmax)).max(topdeg) {
            let mut span = vec![];
            for (a, &(_, ad)) in tb.ins.iter().enumerate() {
                let rd = d - ad;
                if rd < 0 || rd % 2 != 0 { continue; }
                for mo in monomials(n, (rd / 2) as usize) {
                    if mo.iter().all(|&e| e <= (nn - 1) as u32) { span.push((a, mo)); }
                }
            }
            let dd = topdeg - d;
            let mut dspan = vec![];
            for (b, &(_, bd)) in tb.outs.iter().enumerate() {
                let rd = dd - bd;
                if rd < 0 || rd % 2 != 0 { continue; }
                for mo in monomials(n, (rd / 2) as usize) {
                    if mo.iter().all(|&e| e <= (nn - 1) as u32) { dspan.push((b, mo)); }
                }
            }
            if span.is_empty() || dspan.is_empty() { continue; }
            let gram: Vec<Vec<F>> = span.iter().map(|(a, m1)| dspan.iter().map(|(b, m2)| coeff(&tb.gram[*b][*a], m1, m2)).collect()).collect();
            let rows = indep_rows(&gram);
            if rows.is_empty() { continue; }
            let sub: Vec<Vec<F>> = rows.iter().map(|&k| gram[k].clone()).collect();
            let cols = indep_rows(&transpose(&sub));
            let g: Vec<Vec<F>> = sub.iter().map(|r| cols.iter().map(|&c| r[c].clone()).collect()).collect();
            reps.insert(d, rows.iter().map(|&k| span[k].clone()).collect::<Vec<_>>());
            duals.insert(d, cols.iter().map(|&c| dspan[c].clone()).collect::<Vec<_>>());
            ginvs.insert(d, crate::local::inverse(&g).unwrap());
        }
        ssdim += reps.values().map(|v: &Vec<(usize, Vec<u32>)>| v.len()).sum::<usize>();
        ss.insert(x, SS { reps, duals, ginv: ginvs });
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
    let gindex: HashMap<(usize, i32, usize), usize> = groups.values().flat_map(|v| v.iter().enumerate().map(|(k, x)| (*x, k))).collect();
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
    // packed monomials: 5 bits per variable (exponents ≤ 31), n ≤ 12
    assert!(n <= 12, "packed closure supports at most 12 strands");
    let pack = |m: &[u32]| -> u64 { m.iter().enumerate().fold(0u64, |acc, (i, &e)| acc | ((e as u64) << (5 * i))) };
    let toppack = pack(&top);
    let fits = |s: u64| -> bool { (0..n).all(|i| ((s >> (5 * i)) & 31) <= (nn - 1) as u64) };
    // columns (duals) sorted by a lex term order on m2 (ascending) so that the
    // leftmost entry of a row m1·F is its leading term m1·lt(F): little fill-in.
    let lexkey = |m: &[u32]| -> u64 { m.iter().enumerate().fold(0u64, |acc, (i, &e)| acc | ((e as u64) << (5 * (n - 1 - i)))) };
    let mut dual_lists: BTreeMap<(i32, i32), Vec<(u64, usize, usize, u64)>> = BTreeMap::new();
    for (i, o) in objs.iter().enumerate() {
        for (d, duals) in &ss[&o.x].duals {
            let key = (o.h, d + o.q + global_q);
            for (b, m2) in duals {
                dual_lists.entry(key).or_default().push((lexkey(m2), i, *b, pack(m2)));
            }
        }
    }
    let mut dual_groups: BTreeMap<(i32, i32), usize> = BTreeMap::new();
    let mut dual_index: HashMap<u128, usize> = HashMap::new();
    for (key, mut list) in dual_lists {
        list.sort();
        for (idx, (_, i, b, m2p)) in list.iter().enumerate() {
            dual_index.insert(((*i as u128) << 72) | ((*b as u128) << 64) | *m2p as u128, idx);
        }
        dual_groups.insert(key, list.len());
    }
    // pruned packed P_ab terms
    let pab_packed: Vec<Vec<Vec<Vec<(u64, F)>>>> = pab
        .iter()
        .map(|m| m.iter().map(|row| row.iter().map(|p| p.iter().filter(|(t, _)| t.iter().all(|&e| e <= (nn - 1) as u32)).map(|(t, c)| (pack(t), c.clone())).collect()).collect()).collect())
        .collect();
    let blocks: Vec<((i32, i32), Vec<(usize, i32, usize)>)> = groups.iter().map(|(k, v)| (*k, v.clone())).collect();
    let ranks: Vec<((i32, i32), usize)> = par_map(&blocks, |((h, q), elems)| {
        let tl = dual_groups.get(&(h + 1, *q)).copied().unwrap_or(0);
        if tl == 0 { return ((*h, *q), 0); }
        let mut ech = crate::linsys::RankEchelon::<F>::new(tl);
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
                        if let Some(&col) = dual_index.get(&(((t as u128) << 72) | ((b as u128) << 64) | m2 as u128)) {
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
    let _ = &gindex;
    eprintln!("[close N={}] chain dim {}, ranks {:.2}s", nn, groups.values().map(|v| v.len()).sum::<usize>(), trk.elapsed().as_secs_f64());
    let mut pc = Poincare::default();
    for ((h, q), elems) in &groups {
        let r_out = rank.get(&(*h, *q)).copied().unwrap_or(0);
        let r_in = rank.get(&(h - 1, *q)).copied().unwrap_or(0);
        pc.add_term(*h, *q, elems.len() as i128 - r_out as i128 - r_in as i128);
    }
    pc
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
