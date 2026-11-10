//! Bar-Natan style scanning engine for uncolored braid closures in the localized
//! Soergel calculus (milestone 6F). Each evaluation point runs the same
//! deterministic algorithm; the gl_N closure combines all points through
//! interpolation weights for τ_N.

use crate::crossing::Gen;
use crate::field::Field;
use crate::hecke::Hecke;
use crate::indec::{decompose, indecomposable, kl_mult, ops_to_lmat, Indec, LMat, Summand};
use crate::laurent::Poincare;
use crate::linsys::Echelon;
use crate::soergel::{double_leaves_all, Ctx, Op};
use crate::sparse::SpVec;
use std::collections::{BTreeMap, BTreeSet, HashMap};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Obj {
    pub x: usize,
    pub q: i32,
    pub h: i32,
}

/// Complex of indecomposables: objects and differential components (LMat at the point).
pub struct Cx<F: Field> {
    pub objs: BTreeMap<usize, Obj>,
    pub out: BTreeMap<usize, BTreeMap<usize, LMat<F>>>,
    pub inn: BTreeMap<usize, BTreeSet<usize>>,
    next_id: usize,
}

impl<F: Field> Cx<F> {
    pub fn new() -> Self {
        Cx { objs: BTreeMap::new(), out: BTreeMap::new(), inn: BTreeMap::new(), next_id: 0 }
    }
    pub fn add_obj(&mut self, o: Obj) -> usize {
        let id = self.next_id;
        self.next_id += 1;
        self.objs.insert(id, o);
        self.out.insert(id, BTreeMap::new());
        self.inn.insert(id, BTreeSet::new());
        id
    }
    pub fn add_map(&mut self, a: usize, b: usize, f: LMat<F>) {
        if f.is_zero() {
            return;
        }
        let e = self.out.get_mut(&a).unwrap();
        match e.get_mut(&b) {
            Some(g) => {
                g.add_scaled(&f, &F::one());
                if g.is_zero() {
                    e.remove(&b);
                    self.inn.get_mut(&b).unwrap().remove(&a);
                }
            }
            None => {
                e.insert(b, f);
                self.inn.get_mut(&b).unwrap().insert(a);
            }
        }
    }
    fn remove_obj(&mut self, a: usize) {
        if let Some(outs) = self.out.remove(&a) {
            for b in outs.keys() {
                self.inn.get_mut(b).unwrap().remove(&a);
            }
        }
        if let Some(ins) = self.inn.remove(&a) {
            for b in ins {
                self.out.get_mut(&b).unwrap().remove(&a);
            }
        }
        self.objs.remove(&a);
    }
    pub fn len(&self) -> usize {
        self.objs.len()
    }
}

/// Per-point data and caches.
pub struct PointEngine<F: Field> {
    pub ctx: Ctx<F>,
    pub hecke: std::sync::Arc<Hecke>,
    pub indecs: HashMap<usize, Indec<F>>,
    decomp: HashMap<(usize, usize), Vec<Summand<F>>>,
}

impl<F: Field> PointEngine<F> {
    pub fn new(n: usize, x: Vec<F>, hecke: std::sync::Arc<Hecke>) -> Self {
        let ctx = Ctx::new(n, x);
        PointEngine { ctx, hecke, indecs: HashMap::new(), decomp: HashMap::new() }
    }

    /// Compute the indecomposable B_x on demand.
    pub fn ensure(&mut self, x: usize) {
        if !self.indecs.contains_key(&x) {
            let ind = indecomposable(&mut self.ctx, x);
            self.indecs.insert(x, ind);
        }
    }

    fn decomposition(&mut self, x: usize, s: usize) -> Vec<Summand<F>> {
        if let Some(d) = self.decomp.get(&(x, s)) {
            return d.clone();
        }
        self.ensure(x);
        let mult0 = kl_mult(&self.hecke, x, s);
        for &(y, _, _) in &mult0 {
            self.ensure(y);
        }
        let ix = &self.indecs[&x];
        let mut word = ix.word.clone();
        word.push(s);
        let e = ix.e.tensor_letter(s);
        let mult: Vec<(usize, i32, usize)> = kl_mult(&self.hecke, x, s).into_iter().map(|(y, k, m)| (y, -k, m)).collect();
        let d = decompose(&mut self.ctx, &word, &e, &mult, &self.indecs);
        self.decomp.insert((x, s), d.clone());
        d
    }

    /// Tensor the complex with the crossing complex of g on the right (on top).
    pub fn tensor_crossing(&mut self, c: &Cx<F>, g: Gen) -> Cx<F> {
        let s = g.i;
        let mut nc = Cx::new();
        // F-terms: positive: R (q 0, h 0) -> B_s (q -1, h 1); negative: B_s (q 1, h -1) -> R (q 0, h 0)
        let (rq, rh, bq, bh) = if g.positive { (0, 0, -1, 1) } else { (0, 0, 1, -1) };
        // map old id -> new R-term id; old id -> list of (summand index, new id)
        for o in c.objs.values() {
            self.ensure(o.x);
        }
        let mut r_id: BTreeMap<usize, usize> = BTreeMap::new();
        let mut b_ids: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        let mut dec_of: BTreeMap<usize, Vec<Summand<F>>> = BTreeMap::new();
        for (&id, o) in &c.objs {
            r_id.insert(id, nc.add_obj(Obj { x: o.x, q: o.q + rq, h: o.h + rh }));
            let d = self.decomposition(o.x, s);
            let mut ids = vec![];
            for sm in &d {
                ids.push(nc.add_obj(Obj { x: sm.y, q: o.q + bq + sm.deg, h: o.h + bh }));
            }
            b_ids.insert(id, ids);
            dec_of.insert(id, d);
        }
        // old differentials ⊗ id with Koszul sign (-1)^{h_F}
        let sign_r = if rh % 2 == 0 { F::one() } else { F::one().neg() };
        let sign_b = if bh.rem_euclid(2) == 0 { F::one() } else { F::one().neg() };
        for (&a, outs) in &c.out {
            for (&b, f) in outs {
                nc.add_map(r_id[&a], r_id[&b], f.scale(&sign_r));
                let ft = f.tensor_letter(s);
                let da = &dec_of[&a];
                let db = &dec_of[&b];
                for (k, ska) in da.iter().enumerate() {
                    let fi = ft.compose(&ska.iota);
                    for (l, slb) in db.iter().enumerate() {
                        let comp = slb.pi.compose(&fi);
                        nc.add_map(b_ids[&a][k], b_ids[&b][l], comp.scale(&sign_b));
                    }
                }
            }
        }
        // crossing differential id_X ⊗ d_F
        for (&a, o) in &c.objs {
            let word = self.indecs[&o.x].word.clone();
            let e = self.indecs[&o.x].e.clone();
            let da = dec_of[&a].clone();
            if g.positive {
                // X -> X ⊗ B_s via DotStart at the end
                let dot = ops_to_lmat(&mut self.ctx, &word, &[(Op::DotStart(s), word.len())]).compose(&e);
                let dot = e.tensor_letter(s).compose(&dot);
                for (k, sm) in da.iter().enumerate() {
                    nc.add_map(r_id[&a], b_ids[&a][k], sm.pi.compose(&dot));
                }
            } else {
                let mut ws = word.clone();
                ws.push(s);
                let dot = ops_to_lmat(&mut self.ctx, &ws, &[(Op::DotEnd(s), word.len())]);
                let dot = e.compose(&dot);
                for (k, sm) in da.iter().enumerate() {
                    nc.add_map(b_ids[&a][k], r_id[&a], dot.compose(&sm.iota));
                }
            }
        }
        nc
    }

    /// Gaussian elimination of all degree-0 isomorphism components.
    pub fn simplify(&self, c: &mut Cx<F>) {
        loop {
            let mut found = None;
            'search: for (&a, outs) in &c.out {
                let oa = c.objs[&a];
                for (&b, f) in outs {
                    let ob = c.objs[&b];
                    if oa.x != ob.x || oa.q != ob.q {
                        continue;
                    }
                    let top = (1u64 << self.indecs[&oa.x].word.len()) - 1;
                    if let Some(v) = f.cols[top as usize].get(&top) {
                        if !v.is_zero() {
                            found = Some((a, b, v.clone()));
                            break 'search;
                        }
                    }
                }
            }
            let Some((a, b, cval)) = found else { break };
            let cinv = cval.inv();
            let ins: Vec<usize> = c.inn[&b].iter().copied().filter(|&k| k != a).collect();
            let outs: Vec<usize> = c.out[&a].keys().copied().filter(|&l| l != b).collect();
            let mut updates = vec![];
            for &k in &ins {
                let dkb = c.out[&k][&b].clone();
                for &l in &outs {
                    let dal = &c.out[&a][&l];
                    let m = dal.compose(&dkb).scale(&cinv.neg());
                    updates.push((k, l, m));
                }
            }
            c.remove_obj(a);
            c.remove_obj(b);
            for (k, l, m) in updates {
                c.add_map(k, l, m);
            }
        }
    }

    /// Fixed points of BS(rex x) with product e (sorted bits).
    pub fn ebits(&self, x: usize) -> Vec<u64> {
        let word = &self.indecs[&x].word;
        (0..(1u64 << word.len()))
            .filter(|&b| {
                let sub: Vec<usize> = word.iter().enumerate().filter(|(j, _)| b >> j & 1 == 1).map(|(_, &s)| s).collect();
                self.ctx.g.from_word(&sub) == 0
            })
            .collect()
    }

    /// Candidate generators restricted to the e-block: in-candidates as vectors over
    /// ebits(x), out-candidates as functionals over ebits(x).
    pub fn closure_cands_e(&mut self, x: usize) -> (Vec<(i32, Vec<F>)>, Vec<(i32, Vec<F>)>) {
        self.ensure(x);
        let eb = self.ebits(x);
        let (ins, outs) = self.closure_gens(x);
        let iv = ins
            .into_iter()
            .map(|(d, m)| (d, eb.iter().map(|b| m.cols[0].get(b).cloned().unwrap_or_else(F::zero)).collect()))
            .collect();
        let ov = outs
            .into_iter()
            .map(|(d, m)| (d, eb.iter().map(|b| m.cols[*b as usize].get(&0).cloned().unwrap_or_else(F::zero)).collect()))
            .collect();
        (iv, ov)
    }

    /// Generators of Hom(R, B_x) (as maps BS(∅) -> BS(word)) and of Hom(B_x, R), with degrees.
    pub fn closure_gens(&mut self, x: usize) -> (Vec<(i32, LMat<F>)>, Vec<(i32, LMat<F>)>) {
        self.ensure(x);
        let word = self.indecs[&x].word.clone();
        let e = self.indecs[&x].e.clone();
        let mut ins = vec![];
        let mut outs = vec![];
        let maxd = 2 * word.len() as i32 + 2;
        // one pass over all degrees (sorted by degree, then leaf), as degree-by-degree calls
        for (d, _, _, ops) in double_leaves_all(&self.ctx, &[], &word) {
            if d.abs() > maxd {
                continue;
            }
            let m = e.compose(&ops_to_lmat(&mut self.ctx, &[], &ops));
            if !m.is_zero() {
                ins.push((d, m));
            }
        }
        for (d, _, _, ops) in double_leaves_all(&self.ctx, &word, &[]) {
            if d.abs() > maxd {
                continue;
            }
            let m = ops_to_lmat(&mut self.ctx, &word, &ops).compose(&e);
            if !m.is_zero() {
                outs.push((d, m));
            }
        }
        (ins, outs)
    }
}

pub fn scalar_of<F: Field>(m: &LMat<F>) -> F {
    m.cols[0].get(&0).cloned().unwrap_or_else(F::zero)
}

/// Monomials of degree `deg` (each variable weight 1) in n variables.
pub fn monomials(n: usize, deg: usize) -> Vec<Vec<u32>> {
    let mut out = vec![];
    fn rec(i: usize, n: usize, left: usize, cur: &mut Vec<u32>, out: &mut Vec<Vec<u32>>) {
        if i == n - 1 {
            cur.push(left as u32);
            out.push(cur.clone());
            cur.pop();
            return;
        }
        for k in 0..=left {
            cur.push(k as u32);
            rec(i + 1, n, left - k, cur, out);
            cur.pop();
        }
    }
    rec(0, n, deg, &mut vec![], &mut out);
    out
}

pub fn mono_val<F: Field>(m: &[u32], x: &[F]) -> F {
    let mut v = F::one();
    for (i, &k) in m.iter().enumerate() {
        for _ in 0..k {
            v = v.mul(&x[i]);
        }
    }
    v
}

pub struct EngineResult {
    pub poincare: Poincare,
    pub final_objects: usize,
    pub max_objects: usize,
}

/// gl_N homology of the closure of an uncolored braid on n strands.
pub fn khr_engine<F: Field>(n: usize, nn: usize, word: &[Gen], seed: u64) -> EngineResult {
    khr_engine_start::<F>(n, nn, word, seed, &[])
}

/// As `khr_engine`, starting from the indecomposable B_x (x given as a word in the
/// simple reflections; empty = R) instead of R.
pub fn khr_engine_start<F: Field>(n: usize, nn: usize, word: &[Gen], seed: u64, start: &[usize]) -> EngineResult {
    use rand::{rngs::StdRng, Rng, SeedableRng};
    let t_start = std::time::Instant::now();
    let hecke = std::sync::Arc::new(Hecke::new(n));
    // τ_N weights: polynomials of degree D = n(N-1) in n variables (q-degree 2D)
    let dtop = n * (nn - 1);
    let monos = monomials(n, dtop);
    let mpts = monos.len();
    let mut rng = StdRng::seed_from_u64(seed);
    let pts: Vec<Vec<F>> = (0..mpts).map(|_| (0..n).map(|_| F::from_i64(rng.gen_range(1..2_000_000_000))).collect()).collect();
    // weights w: Σ_j w_j m(pt_j) = δ_{m = Π x^{N-1}}
    let target = vec![(nn - 1) as u32; n];
    let vmat: Vec<Vec<F>> = monos.iter().map(|m| pts.iter().map(|p| mono_val(m, p)).collect()).collect();
    let rhs: Vec<F> = monos.iter().map(|m| if *m == target { F::one() } else { F::zero() }).collect();
    let w = crate::dense::solve(&vmat, mpts, &rhs).expect("weights");
    // run the scanning at every point (in parallel)
    let nthreads = crate::num_threads();
    let chunks: Vec<Vec<usize>> = (0..nthreads).map(|t| (t..mpts).step_by(nthreads).collect()).collect();
    type PerPoint<F> = (usize, Vec<(Obj, usize)>, Vec<(usize, usize, Vec<Vec<F>>)>, usize, Vec<(usize, Vec<(i32, F)>, Vec<(i32, F)>, Vec<Vec<F>>)>);
    let results: Vec<Vec<PerPoint<F>>> = std::thread::scope(|scope| {
        let handles: Vec<_> = chunks
            .iter()
            .map(|idxs| {
                let pts = &pts;
                let hecke = hecke.clone();
                scope.spawn(move || {
                    let mut res = vec![];
                    for &j in idxs {
                        res.push(run_point::<F>(n, word, pts[j].clone(), hecke.clone(), j, start));
                    }
                    res
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    let mut per_point: Vec<Option<PerPoint<F>>> = (0..mpts).map(|_| None).collect();
    for r in results.into_iter().flatten() {
        let j = r.0;
        per_point[j] = Some(r);
    }
    let per_point: Vec<PerPoint<F>> = per_point.into_iter().map(|x| x.unwrap()).collect();
    let t_scan = t_start.elapsed().as_secs_f64();
    // consistency: same final objects at every point
    let objs0 = per_point[0].1.clone();
    for p in &per_point {
        assert_eq!(p.1, objs0, "different final complexes at different points");
    }
    let max_objects = per_point[0].3;
    // assemble closure: for each final object (x, q, h) with generators ins (deg), outs (deg)
    // values P[a][b] at each point = (out_b ∘ in_a) for the state space,
    // and for each edge: (out_b ∘ f ∘ in_a).
    let global_q: i32 = word.iter().map(|g| if g.positive { -((nn - 1) as i32) } else { (nn - 1) as i32 }).sum::<i32>() - (n * (nn - 1)) as i32;
    let ndeg = 2 * dtop as i32;
    let rmonos_cache: HashMap<i32, Vec<Vec<u32>>> = (0..=dtop as i32).map(|d| (d, monomials(n, d as usize))).collect();
    // For object index i: list of state-space basis elements per degree: (in-gen a, monomial r) and dual (out-gen b, monomial r')
    struct ObjSS<F: Field> {
        reps: BTreeMap<i32, Vec<(usize, Vec<u32>)>>,
        duals: BTreeMap<i32, Vec<(usize, Vec<u32>)>>,
        ginv: BTreeMap<i32, Vec<Vec<F>>>,
    }
    // monomial values per point: mval[deg][(monomial index)][j]
    let mval: HashMap<i32, Vec<Vec<F>>> = rmonos_cache
        .iter()
        .map(|(&d, ms)| (d, ms.iter().map(|m| pts.iter().map(|p| mono_val(m, p)).collect()).collect()))
        .collect();
    let mindex: HashMap<Vec<u32>, usize> = rmonos_cache.values().flat_map(|ms| ms.iter().enumerate().map(|(i, m)| (m.clone(), i))).collect();
    let mdeg = |m: &[u32]| -> i32 { m.iter().sum::<u32>() as i32 };
    let tau = |vals: &dyn Fn(usize) -> F, r: &[u32], r2: &[u32]| -> F {
        let a = &mval[&mdeg(r)][mindex[r]];
        let b = &mval[&mdeg(r2)][mindex[r2]];
        let mut s = F::zero();
        for j in 0..mpts {
            let v = vals(j);
            if v.is_zero() {
                continue;
            }
            s.fma(&w[j].mul(&v), &a[j].mul(&b[j]));
        }
        s
    };
    // block Gram: for fixed (a,b) values P(j), G[r][r2] = Σ_j w_j P(j) m_r(j) m_r2(j)
    let block = |pv: &[F], d1: i32, d2: i32| -> Vec<Vec<F>> {
        let m1 = &mval[&d1];
        let m2 = &mval[&d2];
        let wp: Vec<F> = (0..mpts).map(|j| w[j].mul(&pv[j])).collect();
        m1.iter()
            .map(|ra| {
                let t: Vec<F> = (0..mpts).map(|j| ra[j].mul(&wp[j])).collect();
                m2.iter().map(|rb| {
                    let mut s = F::zero();
                    for j in 0..mpts {
                        s.fma(&t[j], &rb[j]);
                    }
                    s
                }).collect()
            })
            .collect()
    };
    let mut ss: Vec<ObjSS<F>> = vec![];
    for (i, (o, _)) in objs0.iter().enumerate() {
        let (_, ins, outs, gram_pts) = {
            let p = &per_point[0].4[i];
            (p.0, p.1.clone(), p.2.clone(), ())
        };
        let _ = gram_pts;
        let mut reps = BTreeMap::new();
        let mut duals = BTreeMap::new();
        let mut ginvs = BTreeMap::new();
        for d in 0..=ndeg {
            // spanning set in degree d: (a, r) with deg(in_a) + 2|r| = d
            let mut span = vec![];
            for (a, (da, _)) in ins.iter().enumerate() {
                let rd = d - *da;
                if rd < 0 || rd % 2 != 0 || rd / 2 > dtop as i32 {
                    continue;
                }
                for r in &rmonos_cache[&(rd / 2)] {
                    span.push((a, r.clone()));
                }
            }
            let dd = ndeg - d;
            let mut dspan = vec![];
            for (b, (db, _)) in outs.iter().enumerate() {
                let rd = dd - *db;
                if rd < 0 || rd % 2 != 0 || rd / 2 > dtop as i32 {
                    continue;
                }
                for r in &rmonos_cache[&(rd / 2)] {
                    dspan.push((b, r.clone()));
                }
            }
            if span.is_empty() || dspan.is_empty() {
                continue;
            }
            let mut gram: Vec<Vec<F>> = vec![vec![F::zero(); dspan.len()]; span.len()];
            {
                // group span by generator a, dspan by generator b
                let mut sa: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
                for (k, (a, _)) in span.iter().enumerate() { sa.entry(*a).or_default().push(k); }
                let mut sb: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
                for (k, (b, _)) in dspan.iter().enumerate() { sb.entry(*b).or_default().push(k); }
                for (a, ka) in &sa {
                    for (b, kb) in &sb {
                        let pv: Vec<F> = (0..mpts).map(|j| per_point[j].4[i].3[*a][*b].clone()).collect();
                        if pv.iter().all(|x| x.is_zero()) { continue; }
                        let d1 = mdeg(&span[ka[0]].1);
                        let d2 = mdeg(&dspan[kb[0]].1);
                        let blk = block(&pv, d1, d2);
                        for &k1 in ka {
                            let r1 = mindex[&span[k1].1];
                            for &k2 in kb {
                                let r2 = mindex[&dspan[k2].1];
                                gram[k1][k2] = blk[r1][r2].clone();
                            }
                        }
                    }
                }
            }
            let rows = indep_rows(&gram);
            if rows.is_empty() {
                continue;
            }
            let sub: Vec<Vec<F>> = rows.iter().map(|&k| gram[k].clone()).collect();
            let cols = indep_rows(&transpose(&sub));
            let g: Vec<Vec<F>> = sub.iter().map(|r| cols.iter().map(|&c| r[c].clone()).collect()).collect();
            let gi = crate::local::inverse(&g).unwrap();
            reps.insert(d, rows.iter().map(|&k| span[k].clone()).collect::<Vec<_>>());
            duals.insert(d, cols.iter().map(|&c| dspan[c].clone()).collect::<Vec<_>>());
            ginvs.insert(d, gi);
        }
        let _ = o;
        ss.push(ObjSS { reps, duals, ginv: ginvs });
    }
    // chain groups
    let mut groups: BTreeMap<(i32, i32), Vec<(usize, i32, usize)>> = BTreeMap::new();
    for (i, (o, _)) in objs0.iter().enumerate() {
        for (d, reps) in &ss[i].reps {
            for k in 0..reps.len() {
                groups.entry((o.h, d + o.q + global_q)).or_default().push((i, *d, k));
            }
        }
    }
    let gindex: HashMap<(usize, i32, usize), usize> = groups.values().flat_map(|v| v.iter().enumerate().map(|(k, x)| (*x, k))).collect();
    // edges per point: per_point[j].2 = (src obj idx, tgt obj idx, values[a][b] of out_b∘f∘in_a)
    let edges0: Vec<(usize, usize)> = per_point[0].2.iter().map(|(a, b, _)| (*a, *b)).collect();
    let mut rank: HashMap<(i32, i32), usize> = HashMap::new();
    for ((h, q), elems) in &groups {
        let tgt_len = groups.get(&(h + 1, *q)).map(|v| v.len()).unwrap_or(0);
        if tgt_len == 0 {
            continue;
        }
        let mut ech = Echelon::<F>::new(tgt_len);
        for &(i, d, k) in elems {
            let (a, r) = &ss[i].reps[&d][k];
            let mut v: HashMap<usize, F> = HashMap::new();
            for (ei, &(src, tgt)) in edges0.iter().enumerate() {
                if src != i {
                    continue;
                }
                let deg_f = objs0[i].0.q - objs0[tgt].0.q;
                let td = d + deg_f;
                let Some(duals) = ss[tgt].duals.get(&td) else { continue };
                let ginv = &ss[tgt].ginv[&td];
                // pairing vector with duals of target
                let pv: Vec<F> = duals
                    .iter()
                    .map(|(b, r2)| tau(&|j: usize| per_point[j].2[ei].2[*a][*b].clone(), r, r2))
                    .collect();
                let m = ss[tgt].reps[&td].len();
                for s in 0..m {
                    let mut c = F::zero();
                    for (t, x) in pv.iter().enumerate() {
                        if !x.is_zero() {
                            c.fma(x, &ginv[t][s]);
                        }
                    }
                    if !c.is_zero() {
                        let gi = gindex[&(tgt, td, s)];
                        v.entry(gi).or_insert_with(F::zero).add_assign(&c);
                    }
                }
            }
            let mut sv: SpVec<F> = v.into_iter().filter(|(_, x)| !x.is_zero()).map(|(i, x)| (i as u32, x)).collect();
            sv.sort_by_key(|e| e.0);
            ech.insert(&sv);
        }
        rank.insert((*h, *q), ech.rank());
    }
    let mut pc = Poincare::default();
    for ((h, q), elems) in &groups {
        let r_out = rank.get(&(*h, *q)).copied().unwrap_or(0);
        let r_in = rank.get(&(h - 1, *q)).copied().unwrap_or(0);
        pc.add_term(*h, *q, elems.len() as i128 - r_out as i128 - r_in as i128);
    }
    eprintln!("[engine] points {} scan {:.2}s assemble {:.2}s", mpts, t_scan, t_start.elapsed().as_secs_f64() - t_scan);
    EngineResult { poincare: pc, final_objects: objs0.len(), max_objects }
}

fn transpose<F: Field>(m: &[Vec<F>]) -> Vec<Vec<F>> {
    if m.is_empty() {
        return vec![];
    }
    (0..m[0].len()).map(|j| m.iter().map(|r| r[j].clone()).collect()).collect()
}

fn indep_rows<F: Field>(m: &[Vec<F>]) -> Vec<usize> {
    if m.is_empty() {
        return vec![];
    }
    let mut ech = Echelon::<F>::new(m[0].len().max(1));
    let mut out = vec![];
    for (i, r) in m.iter().enumerate() {
        let v: SpVec<F> = r.iter().enumerate().filter(|(_, x)| !x.is_zero()).map(|(j, x)| (j as u32, x.clone())).collect();
        if ech.insert(&v) {
            out.push(i);
        }
    }
    out
}

/// Run the scanning at one point. Returns (point index, final objects with ids,
/// edges with pairing values, max complex size, closure generator data).
#[allow(clippy::type_complexity)]
fn run_point<F: Field>(
    n: usize,
    word: &[Gen],
    x: Vec<F>,
    hecke: std::sync::Arc<Hecke>,
    j: usize,
    start: &[usize],
) -> (usize, Vec<(Obj, usize)>, Vec<(usize, usize, Vec<Vec<F>>)>, usize, Vec<(usize, Vec<(i32, F)>, Vec<(i32, F)>, Vec<Vec<F>>)>) {
    let mut pe = PointEngine::<F>::new(n, x, hecke);
    let mut c = Cx::<F>::new();
    let x0 = pe.ctx.g.from_word(start);
    pe.ensure(x0);
    c.add_obj(Obj { x: x0, q: 0, h: 0 });
    let mut maxo = 1;
    for g in word {
        c = pe.tensor_crossing(&c, *g);
        maxo = maxo.max(c.len());
        pe.simplify(&mut c);
    }
    // final objects in deterministic order
    let ids: Vec<usize> = c.objs.keys().copied().collect();
    let pos: HashMap<usize, usize> = ids.iter().enumerate().map(|(i, &id)| (id, i)).collect();
    let objs: Vec<(Obj, usize)> = ids.iter().map(|id| (c.objs[id], 0)).collect();
    // closure generators per object type (computed once per type)
    let mut gens: HashMap<usize, (Vec<(i32, LMat<F>)>, Vec<(i32, LMat<F>)>)> = HashMap::new();
    for (o, _) in &objs {
        if !gens.contains_key(&o.x) {
            let g = pe.closure_gens(o.x);
            gens.insert(o.x, g);
        }
    }
    let mut cl = vec![];
    for (i, (o, _)) in objs.iter().enumerate() {
        let (ins, outs) = &gens[&o.x];
        let vals: Vec<Vec<F>> = ins.iter().map(|(_, a)| outs.iter().map(|(_, b)| scalar_of(&b.compose(a))).collect()).collect();
        cl.push((i, ins.iter().map(|(d, _)| (*d, F::zero())).collect(), outs.iter().map(|(d, _)| (*d, F::zero())).collect(), vals));
    }
    // fix degrees format: store degrees only
    let cl = cl
        .into_iter()
        .map(|(i, a, b, v)| (i, a, b, v))
        .collect::<Vec<_>>();
    let mut edges = vec![];
    for (&a, outs_m) in &c.out {
        for (&b, f) in outs_m {
            let (ins, _) = &gens[&c.objs[&a].x];
            let (_, outs) = &gens[&c.objs[&b].x];
            let vals: Vec<Vec<F>> = ins.iter().map(|(_, ia)| outs.iter().map(|(_, ob)| scalar_of(&ob.compose(&f.compose(ia)))).collect()).collect();
            edges.push((pos[&a], pos[&b], vals));
        }
    }
    (j, objs, edges, maxo, cl)
}
