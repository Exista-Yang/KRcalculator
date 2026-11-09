//! OBSOLETE MODEL — kept only for the experiments documented in PLAN.md.
//! Truncated edge algebras do NOT model gl_N foams (see PLAN.md v2); use
//! `sbim` (untruncated singular Bott–Samelson bimodules + universal construction).
//!
//! Edge algebras of ladder webs (truncated singular Bott–Samelson bimodules).
//!
//! A ladder web W from bottom labels λ to top labels μ gives the commutative
//! algebra A_W = k[edge variables] / (vertex relations, facet truncations),
//! a bimodule over R_μ (top) and R_λ (bottom), where R_λ = ⊗_i H*_{GL_N}(Gr(λ_i,N)).
//! It is built incrementally: adding a rung that moves r units from an upright of
//! label b to an upright of label a is the free extension P ⊗_{k[c]} Gr(r,b)
//! (c = Chern classes of the b-edge being split), followed by the quotient by
//! the truncation relations of the merged (a+r)-edge.
//!
//! Basis elements are graded by q-degree (edge variable e_j has degree 2j).

use crate::field::Field;
use crate::poly::Mono;
use crate::sparse::{Acc, SpMat, SpVec};
use crate::symfun::GrTable;
use crate::web::{Ladder, Rung};
use std::collections::{BTreeMap, HashMap};

#[derive(Clone, Debug)]
pub struct Alg<F: Field> {
    pub n: usize,
    pub dim: usize,
    pub deg: Vec<i32>,
    pub top_labels: Vec<usize>,
    pub bot_labels: Vec<usize>,
    /// top[i][j-1] = multiplication by e_j of the top edge of upright i
    pub top: Vec<Vec<SpMat<F>>>,
    pub bot: Vec<Vec<SpMat<F>>>,
    /// multiplication by the equivariant parameters E_1..E_N
    pub eq: Vec<SpMat<F>>,
    /// q-shift of the bimodule (normalization so that closures are MOY-symmetric)
    pub shift: i32,
}

/// Result of extending a parent algebra P by a fibre with basis of size fsize:
/// alg = (P ⊗ fibre)/I, proj.cols[p*fsize + β] = image of p ⊗ β.
#[derive(Clone, Debug)]
pub struct Ext<F: Field> {
    pub alg: Alg<F>,
    pub parent_dim: usize,
    pub fsize: usize,
    pub fdeg: Vec<i32>,
    pub proj: SpMat<F>,
}

/// M ⊗ id_fsize, with index convention p*fsize + β.
pub fn tensor_id<F: Field>(m: &SpMat<F>, fsize: usize) -> SpMat<F> {
    let mut cols = Vec::with_capacity(m.ncols * fsize);
    for c in &m.cols {
        for b in 0..fsize {
            cols.push(c.iter().map(|(i, x)| (*i * fsize as u32 + b as u32, x.clone())).collect());
        }
    }
    SpMat { nrows: m.nrows * fsize, ncols: m.ncols * fsize, cols }
}

/// Operators for all c-monomials needed, computed by memoized products.
fn cmonomials<F: Field>(c_ops: &[SpMat<F>], dim: usize, needed: &[Mono]) -> HashMap<Mono, SpMat<F>> {
    let mut memo: HashMap<Mono, SpMat<F>> = HashMap::new();
    let b = c_ops.len();
    memo.insert(vec![0; b], SpMat::identity(dim));
    fn get<F: Field>(g: &Mono, c_ops: &[SpMat<F>], memo: &mut HashMap<Mono, SpMat<F>>) -> SpMat<F> {
        if let Some(m) = memo.get(g) {
            return m.clone();
        }
        let j = g.iter().position(|&k| k > 0).unwrap();
        let mut h = g.clone();
        h[j] -= 1;
        let prev = get(&h, c_ops, memo);
        let m = c_ops[j].matmul(&prev);
        memo.insert(g.clone(), m.clone());
        m
    }
    for g in needed {
        get(g, c_ops, &mut memo);
    }
    memo
}

/// Free extension data: operators of the fibre variables on P ⊗ Gr(r,b).
struct FreeExt<F: Field> {
    dim: usize,
    deg: Vec<i32>,
    fsize: usize,
    fdeg: Vec<i32>,
    e_ops: Vec<SpMat<F>>,
    f_ops: Vec<SpMat<F>>,
}

fn adjoin<F: Field>(p: &Alg<F>, c_ops: &[SpMat<F>], r: usize, b: usize) -> FreeExt<F> {
    assert_eq!(c_ops.len(), b);
    let tab = GrTable::get(r, b);
    let fsize = tab.fsize();
    let fdeg: Vec<i32> = tab.wt.iter().map(|&w| 2 * w as i32).collect();
    let dim = p.dim * fsize;
    let mut deg = Vec::with_capacity(dim);
    for pi in 0..p.dim {
        for bi in 0..fsize {
            deg.push(p.deg[pi] + fdeg[bi]);
        }
    }
    let mut needed: Vec<Mono> = Vec::new();
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
    let cm = cmonomials(c_ops, p.dim, &needed);
    let build = |t: &Vec<Vec<(usize, Vec<(Mono, i64)>)>>| -> SpMat<F> {
        let mut acc = Acc::new(dim);
        let mut cols: Vec<SpVec<F>> = Vec::with_capacity(dim);
        for pi in 0..p.dim {
            for ai in 0..fsize {
                for (bi, cp) in &t[ai] {
                    for (g, k) in cp {
                        let col = &cm[g].cols[pi];
                        let kf = F::from_i64(*k);
                        for (row, x) in col {
                            acc.add(*row * fsize as u32 + *bi as u32, &x.mul(&kf));
                        }
                    }
                }
                cols.push(acc.take());
            }
        }
        SpMat { nrows: dim, ncols: dim, cols }
    };
    let e_ops = tab.e_mul.iter().map(|t| build(t)).collect();
    let f_ops = tab.f_mul.iter().map(|t| build(t)).collect();
    FreeExt { dim, deg, fsize, fdeg, e_ops, f_ops }
}

/// Quotient of a graded space (basis degrees `deg`) by the span of the given
/// homogeneous vectors. Returns (kept old indices, projection matrix new x old).
pub fn quotient_by_span<F: Field>(deg: &[i32], gens: Vec<SpVec<F>>) -> (Vec<usize>, SpMat<F>) {
    let dim = deg.len();
    // group generator vectors by degree
    let mut by_deg: BTreeMap<i32, Vec<SpVec<F>>> = BTreeMap::new();
    for v in gens {
        if v.is_empty() {
            continue;
        }
        let d = deg[v[0].0 as usize];
        debug_assert!(v.iter().all(|(i, _)| deg[*i as usize] == d), "inhomogeneous relation");
        by_deg.entry(d).or_default().push(v);
    }
    let mut block: BTreeMap<i32, Vec<usize>> = BTreeMap::new();
    for (i, &d) in deg.iter().enumerate() {
        block.entry(d).or_default().push(i);
    }
    // pivot_rows[old index] = Some(reduced row as sparse over non-pivot old indices)
    let mut is_pivot = vec![false; dim];
    let mut pivot_row: HashMap<usize, Vec<(usize, F)>> = HashMap::new();
    for (d, vs) in by_deg {
        let coords = &block[&d];
        let nb = coords.len();
        let loc: HashMap<usize, usize> = coords.iter().enumerate().map(|(k, &i)| (i, k)).collect();
        // prefer pivots at large old indices
        let order: Vec<usize> = (0..nb).rev().collect();
        let mut rows: Vec<Vec<F>> = Vec::new();
        let mut rank_cap = false;
        for v in vs {
            if rank_cap {
                break;
            }
            let mut row = vec![F::zero(); nb];
            for (i, x) in v {
                row[loc[&(i as usize)]] = x;
            }
            rows.push(row);
            if rows.len() >= 2 * nb.max(8) {
                crate::dense::rref_with_order(&mut rows, &order);
                if rows.len() == nb {
                    rank_cap = true;
                }
            }
        }
        let piv = crate::dense::rref_with_order(&mut rows, &order);
        for (ri, &pc) in piv.iter().enumerate() {
            let old = coords[pc];
            is_pivot[old] = true;
            let mut entries = Vec::new();
            for (k, x) in rows[ri].iter().enumerate() {
                if k != pc && !x.is_zero() {
                    entries.push((coords[k], x.clone()));
                }
            }
            pivot_row.insert(old, entries);
        }
    }
    let keep: Vec<usize> = (0..dim).filter(|&i| !is_pivot[i]).collect();
    let mut newidx = vec![u32::MAX; dim];
    for (k, &i) in keep.iter().enumerate() {
        newidx[i] = k as u32;
    }
    let mut cols: Vec<SpVec<F>> = Vec::with_capacity(dim);
    for i in 0..dim {
        if !is_pivot[i] {
            cols.push(vec![(newidx[i], F::one())]);
        } else {
            let mut v: SpVec<F> = pivot_row[&i]
                .iter()
                .map(|(j, x)| {
                    debug_assert!(!is_pivot[*j]);
                    (newidx[*j], x.neg())
                })
                .collect();
            v.sort_by_key(|e| e.0);
            cols.push(v);
        }
    }
    let nk = keep.len();
    (keep, SpMat { nrows: nk, ncols: dim, cols })
}

/// Transform an operator M on the old space to the quotient: proj ∘ M ∘ incl.
fn descend<F: Field>(m: &SpMat<F>, keep: &[usize], proj: &SpMat<F>, acc: &mut Acc<F>) -> SpMat<F> {
    let cols = keep.iter().map(|&i| proj.mul_vec(&m.cols[i], acc)).collect();
    SpMat { nrows: keep.len(), ncols: keep.len(), cols }
}

impl<F: Field> Alg<F> {
    /// The ground ring (non-equivariant: E_j = 0).
    pub fn base(n: usize) -> Self {
        Alg {
            n,
            dim: 1,
            deg: vec![0],
            top_labels: vec![],
            bot_labels: vec![],
            top: vec![],
            bot: vec![],
            eq: (0..n).map(|_| SpMat::zero(1, 1)).collect(),
            shift: 0,
        }
    }

    /// Ground ring with constant equivariant parameters E_j = values[j-1]
    /// (an ungraded, filtered specialization).
    pub fn base_with_params(n: usize, values: &[F]) -> Self {
        let mut a = Self::base(n);
        a.eq = values.iter().map(|v| SpMat::scalar(1, v)).collect();
        a
    }

    /// Append a new upright of label `a` (identity strand) on the right.
    pub fn add_upright(&self, a: usize) -> Self {
        assert!(a <= self.n);
        let fe = adjoin(self, &self.eq, a, self.n);
        let tr = |m: &SpMat<F>| tensor_id(m, fe.fsize);
        let mut top: Vec<Vec<SpMat<F>>> = self.top.iter().map(|v| v.iter().map(&tr).collect()).collect();
        let mut bot: Vec<Vec<SpMat<F>>> = self.bot.iter().map(|v| v.iter().map(&tr).collect()).collect();
        top.push(fe.e_ops.clone());
        bot.push(fe.e_ops);
        let mut tl = self.top_labels.clone();
        tl.push(a);
        let mut bl = self.bot_labels.clone();
        bl.push(a);
        Alg { n: self.n, dim: fe.dim, deg: fe.deg, top_labels: tl, bot_labels: bl, top, bot, eq: self.eq.iter().map(&tr).collect(), shift: self.shift }
    }

    /// R_λ: identity web on uprights with the given labels.
    pub fn identity(n: usize, labels: &[usize]) -> Self {
        let mut a = Self::base(n);
        for &l in labels {
            a = a.add_upright(l);
        }
        a
    }

    pub fn identity_over(base: &Alg<F>, labels: &[usize]) -> Self {
        let mut a = base.clone();
        for &l in labels {
            a = a.add_upright(l);
        }
        a
    }

    /// Add a rung on top.
    pub fn add_rung(&self, rung: &Rung) -> Ext<F> {
        let n = self.n;
        let (rec, giv) = rung.ends();
        let a = self.top_labels[rec];
        let b = self.top_labels[giv];
        let r = rung.r;
        assert!(r <= b && a + r <= n, "invalid rung {:?} on {:?}", rung, self.top_labels);
        let fe = adjoin(self, &self.top[giv], r, b);
        let fsize = fe.fsize;
        let tr = |m: &SpMat<F>| tensor_id(m, fsize);
        let dim1 = fe.dim;
        let id1 = SpMat::<F>::identity(dim1);
        let mut top: Vec<Vec<SpMat<F>>> = self.top.iter().map(|v| v.iter().map(&tr).collect()).collect();
        let mut bot: Vec<Vec<SpMat<F>>> = self.bot.iter().map(|v| v.iter().map(&tr).collect()).collect();
        let mut eq: Vec<SpMat<F>> = self.eq.iter().map(&tr).collect();
        // merged edge e(a+r) = e(a) e(r)
        let k = a + r;
        let mut merged = Vec::with_capacity(k);
        for j in 1..=k {
            let mut acc: Option<SpMat<F>> = None;
            for l in 0..=j {
                let m = j - l;
                if l > a || m > r {
                    continue;
                }
                let term = match (l, m) {
                    (0, 0) => id1.clone(),
                    (0, m) => fe.e_ops[m - 1].clone(),
                    (l, 0) => top[rec][l - 1].clone(),
                    (l, m) => top[rec][l - 1].matmul(&fe.e_ops[m - 1]),
                };
                acc = Some(match acc {
                    None => term,
                    Some(x) => x.add(&term),
                });
            }
            merged.push(acc.unwrap());
        }
        top[rec] = merged;
        top[giv] = fe.f_ops.clone();
        let mut tl = self.top_labels.clone();
        tl[rec] = k;
        tl[giv] = b - r;

        // truncation relations of the merged edge (only needed if a > 0 and r > 0)
        let mut keep: Vec<usize> = (0..dim1).collect();
        let mut proj = SpMat::<F>::identity(dim1);
        let mut deg = fe.deg.clone();
        if a > 0 && r > 0 {
            let rels = truncation_relations(&top[rec], &eq, n, dim1);
            let mut gens = Vec::new();
            for m in rels {
                for c in m.cols {
                    if !c.is_empty() {
                        gens.push(c);
                    }
                }
            }
            if !gens.is_empty() {
                let (kp, pr) = quotient_by_span(&fe.deg, gens);
                keep = kp;
                proj = pr;
                deg = keep.iter().map(|&i| fe.deg[i]).collect();
                let mut acc = Acc::new(keep.len().max(1));
                let ds = |m: &SpMat<F>, acc: &mut Acc<F>| descend(m, &keep, &proj, acc);
                top = top.iter().map(|v| v.iter().map(|m| ds(m, &mut acc)).collect()).collect();
                bot = bot.iter().map(|v| v.iter().map(|m| ds(m, &mut acc)).collect()).collect();
                eq = eq.iter().map(|m| ds(m, &mut acc)).collect();
            }
        }
        let dim = keep.len();
        let alg = Alg {
            n,
            dim,
            deg,
            top_labels: tl,
            bot_labels: self.bot_labels.clone(),
            top,
            bot,
            eq,
            shift: self.shift - (r * (b - r)) as i32,
        };
        Ext { alg, parent_dim: self.dim, fsize, fdeg: fe.fdeg, proj }
    }

    /// Free (untruncated) rung extension: P ⊗_{k[c]} Gr(r,b), no truncation.
    /// Returns (algebra, fibre size). The shift is left unchanged.
    pub fn add_rung_free(&self, rung: &Rung) -> (Alg<F>, usize) {
        let (rec, giv) = rung.ends();
        let a = self.top_labels[rec];
        let b = self.top_labels[giv];
        let r = rung.r;
        assert!(r <= b, "invalid rung");
        let fe = adjoin(self, &self.top[giv], r, b);
        let fsize = fe.fsize;
        let tr = |m: &SpMat<F>| tensor_id(m, fsize);
        let id1 = SpMat::<F>::identity(fe.dim);
        let mut top: Vec<Vec<SpMat<F>>> = self.top.iter().map(|v| v.iter().map(&tr).collect()).collect();
        let bot: Vec<Vec<SpMat<F>>> = self.bot.iter().map(|v| v.iter().map(&tr).collect()).collect();
        let eq: Vec<SpMat<F>> = self.eq.iter().map(&tr).collect();
        let k = a + r;
        let mut merged = Vec::with_capacity(k);
        for j in 1..=k {
            let mut acc: Option<SpMat<F>> = None;
            for l in 0..=j {
                let m = j - l;
                if l > a || m > r {
                    continue;
                }
                let term = match (l, m) {
                    (0, 0) => id1.clone(),
                    (0, m) => fe.e_ops[m - 1].clone(),
                    (l, 0) => top[rec][l - 1].clone(),
                    (l, m) => top[rec][l - 1].matmul(&fe.e_ops[m - 1]),
                };
                acc = Some(match acc {
                    None => term,
                    Some(x) => x.add(&term),
                });
            }
            merged.push(acc.unwrap());
        }
        top[rec] = merged;
        top[giv] = fe.f_ops.clone();
        let mut tl = self.top_labels.clone();
        tl[rec] = k;
        tl[giv] = b - r;
        (Alg { n: self.n, dim: fe.dim, deg: fe.deg, top_labels: tl, bot_labels: self.bot_labels.clone(), top, bot, eq, shift: self.shift }, fsize)
    }

    /// Add several rungs; returns the composite extension (fibre index is
    /// mixed-radix: β = ((β1)*f2 + β2)*f3 + ...).
    pub fn add_rungs(&self, rungs: &[Rung]) -> Ext<F> {
        let mut cur = Ext {
            alg: self.clone(),
            parent_dim: self.dim,
            fsize: 1,
            fdeg: vec![0],
            proj: SpMat::identity(self.dim),
        };
        for rg in rungs {
            let next = cur.alg.add_rung(rg);
            cur = compose_ext(&cur, &next);
        }
        cur
    }

    pub fn ladder(lad: &Ladder) -> Self {
        let mut a = Self::identity(lad.n, &lad.bottom);
        for rg in &lad.rungs {
            a = a.add_rung(rg).alg;
        }
        a
    }

    pub fn ladder_over(base: &Alg<F>, lad: &Ladder) -> Self {
        let mut a = Self::identity_over(base, &lad.bottom);
        for rg in &lad.rungs {
            a = a.add_rung(rg).alg;
        }
        a
    }

    /// Basis-degree histogram.
    pub fn graded_dim(&self) -> BTreeMap<i32, usize> {
        let mut h = BTreeMap::new();
        for &d in &self.deg {
            *h.entry(d).or_insert(0) += 1;
        }
        h
    }
}

/// Relations h_m(e(k) - E) = 0 for m = N-k+1..N, as operators.
/// h(t) := E(t)/e(t) as a power series (signs absorbed: we use the
/// alternating convention h~_m = -sum_{i>=1} e_i h~_{m-i}).
pub fn truncation_relations<F: Field>(e: &[SpMat<F>], eq: &[SpMat<F>], n: usize, dim: usize) -> Vec<SpMat<F>> {
    let k = e.len();
    let id = SpMat::<F>::identity(dim);
    let mut h: Vec<SpMat<F>> = vec![id.clone()];
    for m in 1..=n {
        let mut acc = SpMat::<F>::zero(dim, dim);
        for i in 1..=m.min(k) {
            acc = acc.sub(&e[i - 1].matmul(&h[m - i]));
        }
        h.push(acc);
    }
    let eq_zero = eq.iter().all(|m| m.is_zero());
    let mut out = Vec::new();
    for m in (n + 1 - k)..=n {
        if m == 0 {
            continue;
        }
        let mut rho = h[m].clone();
        if !eq_zero {
            for i in 1..=m {
                rho = rho.add(&eq[i - 1].matmul(&h[m - i]));
            }
        }
        out.push(rho);
    }
    out
}

/// Compose extensions: first: P -> A1 (fibre f1), second: A1 -> A2 (fibre f2).
pub fn compose_ext<F: Field>(first: &Ext<F>, second: &Ext<F>) -> Ext<F> {
    let f1 = first.fsize;
    let f2 = second.fsize;
    let fsize = f1 * f2;
    let dim2 = second.alg.dim;
    let mut acc = Acc::new(dim2.max(1));
    let mut cols = Vec::with_capacity(first.parent_dim * fsize);
    for p in 0..first.parent_dim {
        for b1 in 0..f1 {
            let v = &first.proj.cols[p * f1 + b1];
            for b2 in 0..f2 {
                for (s, x) in v {
                    acc.axpy(x, &second.proj.cols[*s as usize * f2 + b2]);
                }
                cols.push(acc.take());
            }
        }
    }
    let mut fdeg = Vec::with_capacity(fsize);
    for b1 in 0..f1 {
        for b2 in 0..f2 {
            fdeg.push(first.fdeg[b1] + second.fdeg[b2]);
        }
    }
    Ext { alg: second.alg.clone(), parent_dim: first.parent_dim, fsize, fdeg, proj: SpMat { nrows: dim2, ncols: first.parent_dim * fsize, cols } }
}
