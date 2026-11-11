//! Singular Bott–Samelson bimodules of ladder webs (untruncated), represented
//! as free modules over the bottom ring R_λ = ⊗_i Sym(e_{i,1..λ_i}) with an
//! explicit fibre basis; top edge variables act by polynomial matrices.
//! Bimodule maps are R_λ-linear polynomial matrices commuting with the top
//! actions. The gl_N state space of a closure is obtained by the universal
//! construction: Hom(R,B) modulo the radical of the pairing
//! (f, g) ↦ τ_N(g∘f), τ_N = top coefficient in ⊗ H*(Gr(λ_i, N)).

use crate::field::Field;
use crate::linsys::Echelon;
use crate::mpoly::{mono_exp, var_mono, Mono, MonoCache, PRing, Poly};
use crate::pmat::{PAcc, PMat, PVec};
use crate::sparse::SpVec;
use crate::symfun::GrTable;
use crate::web::Rung;
use std::collections::{BTreeMap, HashMap};
use std::sync::Arc;

#[derive(Clone, Debug)]
pub struct SBim<F: Field> {
    pub n: usize,
    pub ring: Arc<PRing>,
    /// variable index of e_{i,1} for upright i
    pub voff: Vec<usize>,
    pub bot_labels: Vec<usize>,
    pub top_labels: Vec<usize>,
    /// q-degrees of the fibre basis (unnormalized)
    pub fdeg: Vec<i32>,
    pub top: Vec<Vec<PMat<F>>>,
    /// normalization shift: -sum over rungs r(b-r)
    pub shift: i32,
    /// rungs added so far (bottom to top); fibre index is mixed radix over them
    pub rungs: Vec<Rung>,
}

pub fn bottom_ring(labels: &[usize]) -> (Arc<PRing>, Vec<usize>) {
    let mut w = vec![];
    let mut names = vec![];
    let mut voff = vec![];
    for (i, &l) in labels.iter().enumerate() {
        voff.push(w.len());
        for j in 1..=l {
            w.push(2 * j as i32);
            names.push(format!("e{}_{}", i, j));
        }
    }
    (PRing::new(w, names), voff)
}

impl<F: Field> SBim<F> {
    pub fn identity(n: usize, labels: &[usize]) -> Self {
        let (ring, voff) = bottom_ring(labels);
        let top = labels
            .iter()
            .enumerate()
            .map(|(i, &l)| (0..l).map(|j| PMat::scalar(1, &Poly::var(voff[i] + j))).collect())
            .collect();
        SBim { n, ring, voff, bot_labels: labels.to_vec(), top_labels: labels.to_vec(), fdeg: vec![0], top, shift: 0, rungs: vec![] }
    }

    pub fn fsize(&self) -> usize {
        self.fdeg.len()
    }

    /// Multiplication by the bottom variable e_{i,j} (1-based j) as a scalar.
    pub fn bot_var(&self, i: usize, j: usize) -> Poly<F> {
        Poly::var(self.voff[i] + j - 1)
    }

    /// Free rung extension B ⊗_{R_level} (rung bimodule).
    pub fn add_rung(&self, rung: &Rung) -> Self {
        let (rec, giv) = rung.ends();
        let a = self.top_labels[rec];
        let b = self.top_labels[giv];
        let r = rung.r;
        assert!(r <= b, "invalid rung {:?} on {:?}", rung, self.top_labels);
        if r == 0 {
            return self.clone();
        }
        let tab = GrTable::get(r, b);
        let fs = tab.fsize();
        let f0 = self.fsize();
        let c_ops = &self.top[giv];
        // c-monomial matrices
        let mut needed: Vec<Vec<u16>> = vec![];
        for t in tab.e_mul.iter().chain(tab.f_mul.iter()) {
            for row in t {
                for (_, cp) in row {
                    for (g, _) in cp {
                        needed.push(g.clone());
                    }
                }
            }
        }
        needed.sort();
        needed.dedup();
        let mut cm: HashMap<Vec<u16>, PMat<F>> = HashMap::new();
        cm.insert(vec![0; b], PMat::identity(f0));
        fn get<F: Field>(g: &Vec<u16>, c: &[PMat<F>], memo: &mut HashMap<Vec<u16>, PMat<F>>) -> PMat<F> {
            if let Some(m) = memo.get(g) {
                return m.clone();
            }
            let j = g.iter().position(|&k| k > 0).unwrap();
            let mut h = g.clone();
            h[j] -= 1;
            let prev = get(&h, c, memo);
            let m = c[j].matmul(&prev);
            memo.insert(g.clone(), m.clone());
            m
        }
        for g in &needed {
            get(g, c_ops, &mut cm);
        }
        let dim = f0 * fs;
        let build = |t: &Vec<Vec<(usize, Vec<(Vec<u16>, i64)>)>>| -> PMat<F> {
            let mut cols: Vec<PVec<F>> = Vec::with_capacity(dim);
            for beta in 0..f0 {
                for al in 0..fs {
                    let mut acc = PAcc::new();
                    for (bi, cp) in &t[al] {
                        for (g, k) in cp {
                            let col = &cm[g].cols[beta];
                            let shifted: PVec<F> = col.iter().map(|(row, p)| (*row * fs as u32 + *bi as u32, p.clone())).collect();
                            acc.add_scaled_const(&shifted, &F::from_i64(*k));
                        }
                    }
                    cols.push(acc.take());
                }
            }
            PMat { nrows: dim, ncols: dim, cols }
        };
        let e_ops: Vec<PMat<F>> = tab.e_mul.iter().map(|t| build(t)).collect();
        let f_ops: Vec<PMat<F>> = tab.f_mul.iter().map(|t| build(t)).collect();
        let mut top: Vec<Vec<PMat<F>>> = self.top.iter().map(|v| v.iter().map(|m| m.kron_id(fs)).collect()).collect();
        let k = a + r;
        let id = PMat::<F>::identity(dim);
        let mut merged = Vec::with_capacity(k);
        for j in 1..=k {
            let mut acc: Option<PMat<F>> = None;
            for l in 0..=j {
                let m = j - l;
                if l > a || m > r {
                    continue;
                }
                let term = match (l, m) {
                    (0, 0) => id.clone(),
                    (0, m) => e_ops[m - 1].clone(),
                    (l, 0) => top[rec][l - 1].clone(),
                    (l, m) => top[rec][l - 1].matmul(&e_ops[m - 1]),
                };
                acc = Some(match acc {
                    None => term,
                    Some(x) => x.add(&term),
                });
            }
            merged.push(acc.unwrap());
        }
        top[rec] = merged;
        top[giv] = f_ops;
        let mut tl = self.top_labels.clone();
        tl[rec] = k;
        tl[giv] = b - r;
        let mut fdeg = Vec::with_capacity(dim);
        for beta in 0..f0 {
            for al in 0..fs {
                fdeg.push(self.fdeg[beta] + 2 * tab.wt[al] as i32);
            }
        }
        SBim {
            n: self.n,
            ring: self.ring.clone(),
            voff: self.voff.clone(),
            bot_labels: self.bot_labels.clone(),
            top_labels: tl,
            fdeg,
            top,
            shift: self.shift - (r * (b - r)) as i32,
            rungs: {
                let mut v = self.rungs.clone();
                v.push(*rung);
                v
            },
        }
    }

    pub fn add_rungs(&self, rungs: &[Rung]) -> Self {
        let mut b = self.clone();
        for r in rungs {
            b = b.add_rung(r);
        }
        b
    }

    /// Unnormalized degree d of the element with fibre index beta and monomial m.
    pub fn elem_deg(&self, beta: usize, m: Mono) -> i32 {
        self.fdeg[beta] + self.ring.deg(m)
    }
}

/// Basis of bimodule maps B -> B' of unnormalized degree d
/// (entries Φ[β'][β] ∈ R of degree d + fdeg[β] - fdeg'[β']).
pub fn hom_space<F: Field>(b: &SBim<F>, bp: &SBim<F>, d: i32, mc: &mut MonoCache) -> Vec<PMat<F>> {
    assert_eq!(b.top_labels, bp.top_labels, "hom_space: top labels differ");
    assert_eq!(b.bot_labels, bp.bot_labels, "hom_space: bottom labels differ");
    let f = b.fsize();
    let fp = bp.fsize();
    // unknowns
    let mut unk: Vec<(u32, u32, Mono)> = vec![]; // (row beta', col beta, mono)
    let mut uidx: HashMap<(u32, u32, Mono), u32> = HashMap::new();
    for be in 0..f {
        for bep in 0..fp {
            let dd = d + b.fdeg[be] - bp.fdeg[bep];
            if dd < 0 {
                continue;
            }
            let (ms, _) = mc.get(dd).clone();
            for m in ms {
                uidx.insert((bep as u32, be as u32, m), unk.len() as u32);
                unk.push((bep as u32, be as u32, m));
            }
        }
    }
    if unk.is_empty() {
        return vec![];
    }
    // equations: for each top op pair (T, T'): T' Φ - Φ T = 0.
    // Key: (op index, row beta', col beta, mono) -> sparse row over unknowns
    let mut all_rows: Vec<SpVec<F>> = vec![];
    // group unknowns by (row, col)
    let mut by_entry: HashMap<(u32, u32), Vec<(Mono, u32)>> = HashMap::new();
    for (k, &(r, c, m)) in unk.iter().enumerate() {
        by_entry.entry((r, c)).or_default().push((m, k as u32));
    }
    // by column (for T'Φ): unknowns in column c
    let mut by_col: HashMap<u32, Vec<(u32, Mono, u32)>> = HashMap::new();
    let mut by_row: HashMap<u32, Vec<(u32, Mono, u32)>> = HashMap::new();
    for (k, &(r, c, m)) in unk.iter().enumerate() {
        by_col.entry(c).or_default().push((r, m, k as u32));
        by_row.entry(r).or_default().push((c, m, k as u32));
    }
    let tp_t: Vec<(PMat<F>, PMat<F>)> = b
        .top
        .iter()
        .zip(bp.top.iter())
        .flat_map(|(v, vp)| v.iter().zip(vp.iter()).map(|(t, tp)| (t.transpose(), tp.clone())))
        .collect();
    for (tt, tp) in &tp_t {
        // equations indexed by (beta', beta, mono)
        let mut eqs: HashMap<(u32, u32, Mono), Vec<(u32, F)>> = HashMap::new();
        // T'Φ: (T'Φ)[β'][β] = Σ_γ' T'[β'][γ'] Φ[γ'][β]
        for (gp, col) in tp.cols.iter().enumerate() {
            // T'[β'][γ'] entries are in column γ'
            for (bep, poly) in col {
                // all unknowns Φ[γ'][β] for any β
                if let Some(list) = by_row.get(&(gp as u32)) {
                    for &(be, m, k) in list {
                        for (mm, c) in &poly.t {
                            eqs.entry((*bep, be, m + mm)).or_default().push((k, c.clone()));
                        }
                    }
                }
            }
        }
        // ΦT: (ΦT)[β'][β] = Σ_γ Φ[β'][γ] T[γ][β]; tt = T^T so tt.cols[γ] lists (β, T[γ][β])
        for (g, col) in tt.cols.iter().enumerate() {
            for (be, poly) in col {
                if let Some(list) = by_col.get(&(g as u32)) {
                    for &(bep, m, k) in list {
                        for (mm, c) in &poly.t {
                            eqs.entry((bep, *be, m + mm)).or_default().push((k, c.neg()));
                        }
                    }
                }
            }
        }
        for (_, mut row) in eqs {
            row.sort_unstable_by_key(|e| e.0);
            let mut merged: SpVec<F> = Vec::with_capacity(row.len());
            for (k, c) in row {
                if let Some(last) = merged.last_mut() {
                    if last.0 == k {
                        last.1.add_assign(&c);
                        continue;
                    }
                }
                merged.push((k, c));
            }
            merged.retain(|e| !e.1.is_zero());
            if !merged.is_empty() {
                all_rows.push(merged);
            }
        }
    }
    let _ = by_entry;
    // equations in a deterministic order (the kernel basis does not depend on it)
    all_rows.sort_unstable_by(|x, y| x.len().cmp(&y.len()).then_with(|| x.iter().map(|e| e.0).cmp(y.iter().map(|e| e.0))));
    let ker = crate::linsys::kernel_lazy(unk.len(), &all_rows);
    ker.into_iter()
        .map(|v| {
            let mut cols: Vec<HashMap<u32, Vec<(Mono, F)>>> = vec![HashMap::new(); f];
            for (k, x) in v {
                let (r, c, m) = unk[k as usize];
                cols[c as usize].entry(r).or_default().push((m, x));
            }
            let cols = cols
                .into_iter()
                .map(|h| {
                    let mut v: PVec<F> = h
                        .into_iter()
                        .map(|(r, mut t)| {
                            t.sort_unstable_by_key(|e| e.0);
                            (r, Poly { t })
                        })
                        .collect();
                    v.sort_unstable_by_key(|e| e.0);
                    v
                })
                .collect();
            PMat { nrows: fp, ncols: f, cols }
        })
        .collect()
}

/// Reference implementation of `hom_space` (incremental `Echelon` over every equation); kept
/// for the equality test of the lazy solver.
pub fn hom_space_reference<F: Field>(b: &SBim<F>, bp: &SBim<F>, d: i32, mc: &mut MonoCache) -> Vec<PMat<F>> {
    assert_eq!(b.top_labels, bp.top_labels, "hom_space: top labels differ");
    assert_eq!(b.bot_labels, bp.bot_labels, "hom_space: bottom labels differ");
    let f = b.fsize();
    let fp = bp.fsize();
    // unknowns
    let mut unk: Vec<(u32, u32, Mono)> = vec![]; // (row beta', col beta, mono)
    let mut uidx: HashMap<(u32, u32, Mono), u32> = HashMap::new();
    for be in 0..f {
        for bep in 0..fp {
            let dd = d + b.fdeg[be] - bp.fdeg[bep];
            if dd < 0 {
                continue;
            }
            let (ms, _) = mc.get(dd).clone();
            for m in ms {
                uidx.insert((bep as u32, be as u32, m), unk.len() as u32);
                unk.push((bep as u32, be as u32, m));
            }
        }
    }
    if unk.is_empty() {
        return vec![];
    }
    // equations: for each top op pair (T, T'): T' Φ - Φ T = 0.
    // Key: (op index, row beta', col beta, mono) -> sparse row over unknowns
    let mut ech = Echelon::<F>::new(unk.len());
    // group unknowns by (row, col)
    let mut by_entry: HashMap<(u32, u32), Vec<(Mono, u32)>> = HashMap::new();
    for (k, &(r, c, m)) in unk.iter().enumerate() {
        by_entry.entry((r, c)).or_default().push((m, k as u32));
    }
    // by column (for T'Φ): unknowns in column c
    let mut by_col: HashMap<u32, Vec<(u32, Mono, u32)>> = HashMap::new();
    let mut by_row: HashMap<u32, Vec<(u32, Mono, u32)>> = HashMap::new();
    for (k, &(r, c, m)) in unk.iter().enumerate() {
        by_col.entry(c).or_default().push((r, m, k as u32));
        by_row.entry(r).or_default().push((c, m, k as u32));
    }
    let tp_t: Vec<(PMat<F>, PMat<F>)> = b
        .top
        .iter()
        .zip(bp.top.iter())
        .flat_map(|(v, vp)| v.iter().zip(vp.iter()).map(|(t, tp)| (t.transpose(), tp.clone())))
        .collect();
    for (tt, tp) in &tp_t {
        // equations indexed by (beta', beta, mono)
        let mut eqs: HashMap<(u32, u32, Mono), Vec<(u32, F)>> = HashMap::new();
        // T'Φ: (T'Φ)[β'][β] = Σ_γ' T'[β'][γ'] Φ[γ'][β]
        for (gp, col) in tp.cols.iter().enumerate() {
            // T'[β'][γ'] entries are in column γ'
            for (bep, poly) in col {
                // all unknowns Φ[γ'][β] for any β
                if let Some(list) = by_row.get(&(gp as u32)) {
                    for &(be, m, k) in list {
                        for (mm, c) in &poly.t {
                            eqs.entry((*bep, be, m + mm)).or_default().push((k, c.clone()));
                        }
                    }
                }
            }
        }
        // ΦT: (ΦT)[β'][β] = Σ_γ Φ[β'][γ] T[γ][β]; tt = T^T so tt.cols[γ] lists (β, T[γ][β])
        for (g, col) in tt.cols.iter().enumerate() {
            for (be, poly) in col {
                if let Some(list) = by_col.get(&(g as u32)) {
                    for &(bep, m, k) in list {
                        for (mm, c) in &poly.t {
                            eqs.entry((bep, *be, m + mm)).or_default().push((k, c.neg()));
                        }
                    }
                }
            }
        }
        for (_, mut row) in eqs {
            row.sort_unstable_by_key(|e| e.0);
            let mut merged: SpVec<F> = Vec::with_capacity(row.len());
            for (k, c) in row {
                if let Some(last) = merged.last_mut() {
                    if last.0 == k {
                        last.1.add_assign(&c);
                        continue;
                    }
                }
                merged.push((k, c));
            }
            merged.retain(|e| !e.1.is_zero());
            if !merged.is_empty() {
                ech.insert(&merged);
            }
        }
    }
    let _ = by_entry;
    let ker = ech.kernel();
    ker.into_iter()
        .map(|v| {
            let mut cols: Vec<HashMap<u32, Vec<(Mono, F)>>> = vec![HashMap::new(); f];
            for (k, x) in v {
                let (r, c, m) = unk[k as usize];
                cols[c as usize].entry(r).or_default().push((m, x));
            }
            let cols = cols
                .into_iter()
                .map(|h| {
                    let mut v: PVec<F> = h
                        .into_iter()
                        .map(|(r, mut t)| {
                            t.sort_unstable_by_key(|e| e.0);
                            (r, Poly { t })
                        })
                        .collect();
                    v.sort_unstable_by_key(|e| e.0);
                    v
                })
                .collect();
            PMat { nrows: fp, ncols: f, cols }
        })
        .collect()
}

/// τ_N on monomials of R_λ: top coefficient in ⊗ H*(Gr(λ_i,N)) (non-equivariant).
pub struct Tau<F: Field> {
    pub n: usize,
    pub labels: Vec<usize>,
    pub voff: Vec<usize>,
    cache: Vec<HashMap<Vec<u32>, F>>,
}

impl<F: Field> Tau<F> {
    pub fn new(n: usize, labels: &[usize], voff: &[usize]) -> Self {
        Tau { n, labels: labels.to_vec(), voff: voff.to_vec(), cache: vec![HashMap::new(); labels.len()] }
    }
    pub fn top_degree(&self) -> i32 {
        2 * self.labels.iter().map(|&l| (l * (self.n - l)) as i32).sum::<i32>()
    }
    fn tau_factor(&mut self, i: usize, exps: Vec<u32>) -> F {
        if let Some(x) = self.cache[i].get(&exps) {
            return x.clone();
        }
        let a = self.labels[i];
        let tab = GrTable::get(a, self.n);
        // vector over fibre basis, start at 1 (index 0)
        let mut v: Vec<F> = vec![F::zero(); tab.fsize()];
        v[0] = F::one();
        for (j, &e) in exps.iter().enumerate() {
            for _ in 0..e {
                let mut nv = vec![F::zero(); tab.fsize()];
                for (al, x) in v.iter().enumerate() {
                    if x.is_zero() {
                        continue;
                    }
                    for (be, cp) in &tab.e_mul[j][al] {
                        for (g, k) in cp {
                            if g.iter().all(|&z| z == 0) {
                                nv[*be].fma(x, &F::from_i64(*k));
                            }
                        }
                    }
                }
                v = nv;
            }
        }
        // top element: weight a(N-a)
        let topw = (a * (self.n - a)) as u32;
        let ti = tab.wt.iter().position(|&w| w == topw).unwrap();
        let x = v[ti].clone();
        self.cache[i].insert(exps, x.clone());
        x
    }
    pub fn tau_mono(&mut self, m: Mono) -> F {
        let mut r = F::one();
        for i in 0..self.labels.len() {
            let a = self.labels[i];
            let exps: Vec<u32> = (0..a).map(|j| mono_exp(m, self.voff[i] + j)).collect();
            let x = self.tau_factor(i, exps);
            if x.is_zero() {
                return F::zero();
            }
            r = r.mul(&x);
        }
        r
    }
    pub fn tau(&mut self, p: &Poly<F>) -> F {
        let mut s = F::zero();
        for (m, c) in &p.t {
            let t = self.tau_mono(*m);
            if !t.is_zero() {
                s.fma(c, &t);
            }
        }
        s
    }
}

/// Graded piece of the gl_N state space of the closure of B.
#[derive(Clone, Debug)]
pub struct StatePiece<F: Field> {
    /// representatives in Hom(R,B)_d (as F×1 matrices' single column)
    pub reps: Vec<PVec<F>>,
    /// dual elements in Hom(B,R)_{top-d} (rows, stored as 1×F matrices' columns transposed)
    pub duals: Vec<PMat<F>>,
    /// inverse of the Gram matrix G[s][t] = τ(dual_t · rep_s)
    pub ginv: Vec<Vec<F>>,
}

#[derive(Clone, Debug)]
pub struct StateSpace<F: Field> {
    pub pieces: BTreeMap<i32, StatePiece<F>>,
    pub top: i32,
}

pub fn pair<F: Field>(dual: &PMat<F>, x: &PVec<F>, tau: &mut Tau<F>) -> F {
    // dual is 1×F; x is a column vector
    let mut p = Poly::<F>::zero();
    for (be, xv) in x {
        let col = &dual.cols[*be as usize];
        if let Some((_, q)) = col.first() {
            p.add_assign(&q.mul(xv));
        }
    }
    tau.tau(&p)
}

fn invert<F: Field>(m: &[Vec<F>]) -> Vec<Vec<F>> {
    let n = m.len();
    let mut aug: Vec<Vec<F>> = m
        .iter()
        .enumerate()
        .map(|(i, r)| {
            let mut v = r.clone();
            v.extend((0..n).map(|j| if i == j { F::one() } else { F::zero() }));
            v
        })
        .collect();
    let piv = crate::dense::rref(&mut aug, 2 * n);
    assert!(piv.len() == n && piv.iter().enumerate().all(|(i, &p)| i == p), "singular Gram matrix");
    aug.into_iter().map(|r| r[n..].to_vec()).collect()
}

impl<F: Field> StateSpace<F> {
    pub fn compute(b: &SBim<F>, mc: &mut MonoCache) -> Self {
        assert_eq!(b.top_labels, b.bot_labels);
        let id = SBim::<F>::identity(b.n, &b.bot_labels);
        let mut tau = Tau::<F>::new(b.n, &b.bot_labels, &b.voff);
        let top = tau.top_degree();
        let maxf = *b.fdeg.iter().max().unwrap();
        let mut pieces = BTreeMap::new();
        for d in 0..=(top + maxf) {
            if d % 2 != 0 {
                continue;
            }
            let reps_m = hom_space(&id, b, d, mc);
            if reps_m.is_empty() {
                continue;
            }
            let duals = hom_space(b, &id, top - d, mc);
            if duals.is_empty() {
                continue;
            }
            let reps: Vec<PVec<F>> = reps_m.into_iter().map(|m| m.cols[0].clone()).collect();
            let gram: Vec<Vec<F>> = reps.iter().map(|x| duals.iter().map(|g| pair(g, x, &mut tau)).collect()).collect();
            // choose maximal independent rows and columns
            let rows_idx = independent_rows(&gram);
            if rows_idx.is_empty() {
                continue;
            }
            let sub: Vec<Vec<F>> = rows_idx.iter().map(|&i| gram[i].clone()).collect();
            let cols_idx = independent_rows(&transpose(&sub));
            let g2: Vec<Vec<F>> = sub.iter().map(|r| cols_idx.iter().map(|&j| r[j].clone()).collect()).collect();
            let ginv = invert(&g2);
            pieces.insert(
                d,
                StatePiece { reps: rows_idx.iter().map(|&i| reps[i].clone()).collect(), duals: cols_idx.iter().map(|&j| duals[j].clone()).collect(), ginv },
            );
        }
        StateSpace { pieces, top }
    }

    pub fn dim(&self) -> usize {
        self.pieces.values().map(|p| p.reps.len()).sum()
    }

    /// Coordinates of x ∈ Hom(R,B)_d in the chosen basis of the quotient.
    pub fn coords(&self, d: i32, x: &PVec<F>, tau: &mut Tau<F>) -> Vec<F> {
        let Some(p) = self.pieces.get(&d) else { return vec![] };
        let pv: Vec<F> = p.duals.iter().map(|g| pair(g, x, tau)).collect();
        // c = pv * ginv
        let n = p.reps.len();
        let mut c = vec![F::zero(); n];
        for (t, v) in pv.iter().enumerate() {
            if v.is_zero() {
                continue;
            }
            for s in 0..n {
                c[s].fma(v, &p.ginv[t][s]);
            }
        }
        c
    }
}

fn transpose<F: Field>(m: &[Vec<F>]) -> Vec<Vec<F>> {
    if m.is_empty() {
        return vec![];
    }
    (0..m[0].len()).map(|j| m.iter().map(|r| r[j].clone()).collect()).collect()
}

/// Indices of a maximal set of linearly independent rows (greedy).
fn independent_rows<F: Field>(m: &[Vec<F>]) -> Vec<usize> {
    if m.is_empty() {
        return vec![];
    }
    let ncols = m[0].len();
    let mut ech = Echelon::<F>::new(ncols);
    let mut out = vec![];
    for (i, r) in m.iter().enumerate() {
        let v: SpVec<F> = r.iter().enumerate().filter(|(_, x)| !x.is_zero()).map(|(j, x)| (j as u32, x.clone())).collect();
        if ech.insert(&v) {
            out.push(i);
        }
    }
    out
}

pub fn unit_mono(i: usize) -> Mono {
    var_mono(i)
}
