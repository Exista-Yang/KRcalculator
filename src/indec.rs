//! Indecomposable Soergel bimodules B_x (S_n) in the localized numeric calculus:
//! idempotents in End^0(BS(rex x)) and decompositions of B_x ⊗ B_s.

use crate::field::Field;
use crate::hecke::Hecke;
use crate::soergel::{double_leaves, Ctx, LVec, Op, Word};
use std::collections::HashMap;

/// Localized map BS(src) -> BS(tgt) at the context point: columns indexed by source bits.
#[derive(Clone, Debug)]
pub struct LMat<F: Field> {
    pub src: Word,
    pub tgt: Word,
    pub cols: Vec<LVec<F>>,
}

impl<F: Field> LMat<F> {
    pub fn identity(w: &[usize]) -> Self {
        let cols = (0..(1u64 << w.len())).map(|b| {
            let mut v = HashMap::new();
            v.insert(b, F::one());
            v
        }).collect();
        LMat { src: w.to_vec(), tgt: w.to_vec(), cols }
    }
    pub fn zero(src: &[usize], tgt: &[usize]) -> Self {
        LMat { src: src.to_vec(), tgt: tgt.to_vec(), cols: vec![HashMap::new(); 1 << src.len()] }
    }
    /// self ∘ a  (a: A -> B, self: B -> C)
    pub fn compose(&self, a: &LMat<F>) -> LMat<F> {
        assert_eq!(self.src, a.tgt, "compose: word mismatch");
        let cols = a
            .cols
            .iter()
            .map(|col| {
                let mut out: LVec<F> = HashMap::new();
                for (b, x) in col {
                    for (c, y) in &self.cols[*b as usize] {
                        out.entry(*c).or_insert_with(F::zero).fma(x, y);
                    }
                }
                out.retain(|_, v| !v.is_zero());
                out
            })
            .collect();
        LMat { src: a.src.clone(), tgt: self.tgt.clone(), cols }
    }
    pub fn add_scaled(&mut self, o: &LMat<F>, c: &F) {
        assert_eq!(self.src, o.src);
        assert_eq!(self.tgt, o.tgt);
        for (col, ocol) in self.cols.iter_mut().zip(o.cols.iter()) {
            for (b, x) in ocol {
                col.entry(*b).or_insert_with(F::zero).fma(x, c);
            }
            col.retain(|_, v| !v.is_zero());
        }
    }
    pub fn scale(&self, c: &F) -> LMat<F> {
        let mut z = LMat::zero(&self.src, &self.tgt);
        z.add_scaled(self, c);
        z
    }
    /// f ⊗ id_{B_s}: append letter s to both words.
    pub fn tensor_letter(&self, s: usize) -> LMat<F> {
        let ks = self.src.len();
        let kt = self.tgt.len();
        let mut cols = vec![HashMap::new(); 1 << (ks + 1)];
        for (b, col) in self.cols.iter().enumerate() {
            for last in 0..2u64 {
                let v: LVec<F> = col.iter().map(|(c, x)| (*c | (last << kt), x.clone())).collect();
                cols[b | ((last as usize) << ks)] = v;
            }
        }
        let mut src = self.src.clone();
        src.push(s);
        let mut tgt = self.tgt.clone();
        tgt.push(s);
        LMat { src, tgt, cols }
    }
    pub fn is_zero(&self) -> bool {
        self.cols.iter().all(|c| c.is_empty())
    }
    /// Flatten to a dense vector (for linear algebra on spans).
    pub fn flat(&self) -> Vec<F> {
        let nt = 1usize << self.tgt.len();
        let mut v = vec![F::zero(); self.cols.len() * nt];
        for (c, col) in self.cols.iter().enumerate() {
            for (b, x) in col {
                v[c * nt + *b as usize] = x.clone();
            }
        }
        v
    }
}

pub fn ops_to_lmat<F: Field>(ctx: &mut Ctx<F>, word: &[usize], ops: &[(Op, usize)]) -> LMat<F> {
    let (tw, cols) = ctx.map_columns(word, ops);
    LMat { src: word.to_vec(), tgt: tw, cols }
}

/// Coordinates of vectors in the span of a linearly independent family (via pivot columns).
pub struct SpanSolver<F: Field> {
    pivots: Vec<usize>,
    inv: Vec<Vec<F>>,
}

impl<F: Field> SpanSolver<F> {
    pub fn new(basis: &[Vec<F>]) -> Self {
        let d = basis.len();
        let dim = basis[0].len();
        // choose pivot coordinates: rref of basis rows
        let mut rows = basis.to_vec();
        let piv = crate::dense::rref(&mut rows, dim);
        assert_eq!(piv.len(), d, "basis not independent");
        let sub: Vec<Vec<F>> = (0..d).map(|i| piv.iter().map(|&p| basis[i][p].clone()).collect()).collect();
        // sub[i][j] = basis_i[pivot_j]; we need coords c with sum c_i basis_i[p_j] = v[p_j]
        let subt: Vec<Vec<F>> = (0..d).map(|j| (0..d).map(|i| sub[i][j].clone()).collect()).collect();
        let inv = crate::local::inverse(&subt).expect("pivot submatrix singular");
        SpanSolver { pivots: piv, inv }
    }
    pub fn coords(&self, v: &[F]) -> Vec<F> {
        let d = self.pivots.len();
        let rhs: Vec<F> = self.pivots.iter().map(|&p| v[p].clone()).collect();
        (0..d)
            .map(|i| {
                let mut s = F::zero();
                for j in 0..d {
                    s.fma(&self.inv[i][j], &rhs[j]);
                }
                s
            })
            .collect()
    }
}

// ---------- univariate polynomials over F ----------
fn ptrim<F: Field>(mut p: Vec<F>) -> Vec<F> {
    while p.last().map_or(false, |c| c.is_zero()) {
        p.pop();
    }
    p
}
fn pmul<F: Field>(a: &[F], b: &[F]) -> Vec<F> {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let mut r = vec![F::zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            r[i + j].fma(x, y);
        }
    }
    ptrim(r)
}
fn psub<F: Field>(a: &[F], b: &[F]) -> Vec<F> {
    let n = a.len().max(b.len());
    let mut r = vec![F::zero(); n];
    for i in 0..n {
        if i < a.len() {
            r[i].add_assign(&a[i]);
        }
        if i < b.len() {
            r[i].sub_assign(&b[i]);
        }
    }
    ptrim(r)
}
fn pdivrem<F: Field>(a: &[F], b: &[F]) -> (Vec<F>, Vec<F>) {
    let b = ptrim(b.to_vec());
    let mut r = ptrim(a.to_vec());
    if r.len() < b.len() {
        return (vec![], r);
    }
    let mut q = vec![F::zero(); r.len() - b.len() + 1];
    let inv = b.last().unwrap().inv();
    while r.len() >= b.len() && !r.is_empty() {
        let k = r.len() - b.len();
        let c = r.last().unwrap().mul(&inv);
        q[k] = c.clone();
        for (i, y) in b.iter().enumerate() {
            let t = c.mul(y);
            r[k + i].sub_assign(&t);
        }
        r = ptrim(r);
    }
    (ptrim(q), r)
}
/// extended gcd: returns (g, s, t) with s a + t b = g
fn pxgcd<F: Field>(a: &[F], b: &[F]) -> (Vec<F>, Vec<F>, Vec<F>) {
    let (mut r0, mut r1) = (ptrim(a.to_vec()), ptrim(b.to_vec()));
    let (mut s0, mut s1) = (vec![F::one()], vec![]);
    let (mut t0, mut t1) = (vec![], vec![F::one()]);
    while !r1.is_empty() {
        let (q, r) = pdivrem(&r0, &r1);
        let s2 = psub(&s0, &pmul(&q, &s1));
        let t2 = psub(&t0, &pmul(&q, &t1));
        r0 = r1;
        r1 = r;
        s0 = s1;
        s1 = s2;
        t0 = t1;
        t1 = t2;
    }
    (r0, s0, t0)
}

/// Finite-dimensional algebra given by a basis of localized maps (degree-0 endomorphisms).
pub struct EndAlg<F: Field> {
    pub basis: Vec<LMat<F>>,
    solver: SpanSolver<F>,
}

impl<F: Field> EndAlg<F> {
    pub fn new(basis: Vec<LMat<F>>) -> Self {
        let flat: Vec<Vec<F>> = basis.iter().map(|m| m.flat()).collect();
        let solver = SpanSolver::new(&flat);
        EndAlg { basis, solver }
    }
    pub fn coords(&self, m: &LMat<F>) -> Vec<F> {
        self.solver.coords(&m.flat())
    }
    pub fn elem(&self, c: &[F]) -> LMat<F> {
        let mut z = LMat::zero(&self.basis[0].src, &self.basis[0].tgt);
        for (b, x) in self.basis.iter().zip(c.iter()) {
            if !x.is_zero() {
                z.add_scaled(b, x);
            }
        }
        z
    }
}

/// Primitive idempotent e with χ(e) = 1 where χ(a) = a[top][top] (localized block at `top`).
pub fn primitive_idempotent<F: Field>(alg: &EndAlg<F>, top: u64, seed: u64) -> LMat<F> {
    let d = alg.basis.len();
    let chi = |m: &LMat<F>| -> F { m.cols[top as usize].get(&top).cloned().unwrap_or_else(F::zero) };
    // random element normalized to χ = 1
    let mut state = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
    let mut c: Vec<F> = (0..d)
        .map(|_| {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            F::from_i64(((state >> 33) % 1000003) as i64 + 1)
        })
        .collect();
    let a0 = alg.elem(&c);
    let x = chi(&a0);
    assert!(!x.is_zero(), "random element has χ = 0");
    let xi = x.inv();
    for v in c.iter_mut() {
        *v = v.mul(&xi);
    }
    let a = alg.elem(&c);
    // Krylov: powers of a in coordinates
    let mut pows: Vec<Vec<F>> = vec![];
    let id = LMat::identity(&a.src);
    let mut cur = id.clone();
    let mut minpoly: Vec<F>;
    loop {
        let cc = alg.coords(&cur);
        // check dependence of cc on previous
        let mut rows = pows.clone();
        rows.push(cc.clone());
        if crate::dense::rank(&rows, d) < rows.len() {
            // solve cc = sum λ_i pows_i
            let m = pows.len();
            let mat: Vec<Vec<F>> = (0..d).map(|r| (0..m).map(|i| pows[i][r].clone()).collect()).collect();
            let lam = crate::dense::solve(&mat, m, &cc).expect("krylov");
            minpoly = lam.iter().map(|l| l.neg()).collect();
            minpoly.push(F::one());
            break;
        }
        pows.push(cc);
        cur = a.compose(&cur);
    }
    // split off the (t-1)-primary part
    let tm1 = vec![F::one().neg(), F::one()];
    let mut q = minpoly.clone();
    let mut k = 0;
    loop {
        let (qq, r) = pdivrem(&q, &tm1);
        if r.is_empty() {
            q = qq;
            k += 1;
        } else {
            break;
        }
    }
    assert!(k > 0, "1 is not an eigenvalue");
    let mut pk = vec![F::one()];
    for _ in 0..k {
        pk = pmul(&pk, &tm1);
    }
    let (g, _s, t) = pxgcd(&pk, &q);
    assert_eq!(g.len(), 1, "gcd not constant");
    let ginv = g[0].inv();
    // e = t*q/g evaluated at a  (≡ 1 mod (t-1)^k, ≡ 0 mod q)
    let p: Vec<F> = pmul(&t, &q).into_iter().map(|c| c.mul(&ginv)).collect();
    let mut e = LMat::zero(&a.src, &a.tgt);
    let mut powm = id;
    for c in &p {
        if !c.is_zero() {
            e.add_scaled(&powm, c);
        }
        powm = a.compose(&powm);
    }
    e
}

/// Data of one indecomposable B_x at a point.
pub struct Indec<F: Field> {
    pub x: usize,
    pub word: Word,
    pub e: LMat<F>,
}

pub fn indecomposable<F: Field>(ctx: &mut Ctx<F>, x: usize) -> Indec<F> {
    let word = ctx.g.rex(x);
    let dls = double_leaves(ctx, &word, &word, 0);
    let basis: Vec<LMat<F>> = dls.iter().map(|(_, _, ops)| ops_to_lmat(ctx, &word, ops)).collect();
    let alg = EndAlg::new(basis);
    let top = (1u64 << word.len()) - 1;
    let e = primitive_idempotent(&alg, top, 12345 + x as u64);
    Indec { x, word, e }
}

/// Rank of the localized w-block of an endomorphism (fixed points with product w).
pub fn block_rank<F: Field>(ctx: &Ctx<F>, m: &LMat<F>, w: usize) -> usize {
    let word = &m.src;
    let idx: Vec<u64> = (0..(1u64 << word.len()))
        .filter(|&b| {
            let sub: Vec<usize> = word.iter().enumerate().filter(|(j, _)| b >> j & 1 == 1).map(|(_, &s)| s).collect();
            ctx.g.from_word(&sub) == w
        })
        .collect();
    let rows: Vec<Vec<F>> = idx.iter().map(|&c| idx.iter().map(|&r| m.cols[c as usize].get(&r).cloned().unwrap_or_else(F::zero)).collect()).collect();
    if rows.is_empty() {
        return 0;
    }
    crate::dense::rank(&rows, idx.len())
}

/// Summand data: B_y{shift} ⊂ X, with inclusion ι: BS(y) -> BS(X) and projection π, π∘ι = e_y.
#[derive(Clone)]
pub struct Summand<F: Field> {
    pub y: usize,
    /// degree of ι (as a map BS(rex y) -> BS(word of X)); π has degree -deg
    pub deg: i32,
    pub iota: LMat<F>,
    pub pi: LMat<F>,
}

/// Decompose (BS(word), e) into indecomposables with known multiplicities
/// mult[(y, deg)] = count. Uses local intersection forms of double leaves.
pub fn decompose<F: Field>(
    ctx: &mut Ctx<F>,
    word: &[usize],
    e: &LMat<F>,
    mult: &[(usize, i32, usize)],
    indecs: &HashMap<usize, Indec<F>>,
) -> Vec<Summand<F>> {
    let mut out = vec![];
    for &(y, deg, m) in mult {
        if m == 0 {
            continue;
        }
        let iy = &indecs[&y];
        let ytop = (1u64 << iy.word.len()) - 1;
        let dl_in = double_leaves(ctx, &iy.word, word, deg);
        let dl_out = double_leaves(ctx, word, &iy.word, -deg);
        let ins: Vec<LMat<F>> = dl_in.iter().map(|(_, _, ops)| e.compose(&ops_to_lmat(ctx, &iy.word, ops)).compose(&iy.e)).collect();
        let outs: Vec<LMat<F>> = dl_out.iter().map(|(_, _, ops)| iy.e.compose(&ops_to_lmat(ctx, word, ops)).compose(e)).collect();
        // Gram: G[p][i] = χ_y(π_p ∘ ι_i)
        let gram: Vec<Vec<F>> = outs
            .iter()
            .map(|p| ins.iter().map(|i| p.compose(i).cols[ytop as usize].get(&ytop).cloned().unwrap_or_else(F::zero)).collect())
            .collect();
        // pick m independent rows/cols
        let mut rows_sel = vec![];
        {
            let mut ech = crate::linsys::Echelon::<F>::new(ins.len().max(1));
            for (k, r) in gram.iter().enumerate() {
                let v: Vec<(u32, F)> = r.iter().enumerate().filter(|(_, x)| !x.is_zero()).map(|(j, x)| (j as u32, x.clone())).collect();
                if ech.insert(&v) {
                    rows_sel.push(k);
                }
            }
        }
        assert_eq!(rows_sel.len(), m, "multiplicity mismatch for y={} deg={}: rank {} expected {}", y, deg, rows_sel.len(), m);
        let sub: Vec<Vec<F>> = rows_sel.iter().map(|&k| gram[k].clone()).collect();
        let mut cols_sel = vec![];
        {
            let mut ech = crate::linsys::Echelon::<F>::new(m.max(1));
            for j in 0..ins.len() {
                let v: Vec<(u32, F)> = (0..m).filter(|&k| !sub[k][j].is_zero()).map(|k| (k as u32, sub[k][j].clone())).collect();
                if ech.insert(&v) {
                    cols_sel.push(j);
                }
            }
        }
        let g: Vec<Vec<F>> = (0..m).map(|a| cols_sel.iter().map(|&j| sub[a][j].clone()).collect()).collect();
        let ginv = crate::local::inverse(&g).expect("singular local intersection form");
        // ι'_b = Σ_j ι_{cols_sel[j]} ginv[j][b]  so that π_a ∘ ι'_b = δ_ab e_y
        for b in 0..m {
            let mut io = LMat::zero(&iy.word, word);
            for (j, &cj) in cols_sel.iter().enumerate() {
                if !ginv[j][b].is_zero() {
                    io.add_scaled(&ins[cj], &ginv[j][b]);
                }
            }
            out.push(Summand { y, deg, iota: io, pi: outs[rows_sel[b]].clone() });
        }
    }
    orthogonalize(&mut out, indecs);
    out
}

/// Make π_a ∘ ι_b = δ_ab e_{y_a}: with P = Π∘I = 1 + N (N of positive degree, nilpotent),
/// replace I by I ∘ P^{-1} = I ∘ Σ_k (-N)^k.
fn orthogonalize<F: Field>(sums: &mut Vec<Summand<F>>, indecs: &HashMap<usize, Indec<F>>) {
    let m = sums.len();
    // N[a][b] = π_a ∘ ι_b - δ_ab e_{y_a}
    let mut nmat: Vec<Vec<LMat<F>>> = vec![];
    for a in 0..m {
        let mut row = vec![];
        for b in 0..m {
            let mut c = sums[a].pi.compose(&sums[b].iota);
            if a == b {
                c.add_scaled(&indecs[&sums[a].y].e, &F::one().neg());
            }
            row.push(c);
        }
        nmat.push(row);
    }
    // Q = Σ_k (-N)^k as a matrix of maps BS(y_b) -> BS(y_a); start with identity (e_y on diagonal)
    let ident = |a: usize, b: usize| -> LMat<F> {
        if a == b { indecs[&sums[a].y].e.clone() } else { LMat::zero(&indecs[&sums[b].y].word, &indecs[&sums[a].y].word) }
    };
    let mut q: Vec<Vec<LMat<F>>> = (0..m).map(|a| (0..m).map(|b| ident(a, b)).collect()).collect();
    let mut term: Vec<Vec<LMat<F>>> = (0..m).map(|a| (0..m).map(|b| ident(a, b)).collect()).collect();
    for _ in 0..(2 * m + 4) {
        // term = (-N) ∘ term
        let mut next: Vec<Vec<LMat<F>>> = vec![];
        let mut all_zero = true;
        for a in 0..m {
            let mut row = vec![];
            for b in 0..m {
                let mut acc = LMat::zero(&indecs[&sums[b].y].word, &indecs[&sums[a].y].word);
                for c in 0..m {
                    if nmat[a][c].is_zero() || term[c][b].is_zero() {
                        continue;
                    }
                    acc.add_scaled(&nmat[a][c].compose(&term[c][b]), &F::one().neg());
                }
                if !acc.is_zero() {
                    all_zero = false;
                }
                row.push(acc);
            }
            next.push(row);
        }
        if all_zero {
            break;
        }
        for a in 0..m {
            for b in 0..m {
                q[a][b].add_scaled(&next[a][b], &F::one());
            }
        }
        term = next;
    }
    // new ι_b = Σ_a ι_a ∘ q[a][b]
    let old: Vec<LMat<F>> = sums.iter().map(|s| s.iota.clone()).collect();
    for b in 0..m {
        let mut acc = LMat::zero(&old[b].src, &old[b].tgt);
        for a in 0..m {
            if q[a][b].is_zero() {
                continue;
            }
            acc.add_scaled(&old[a].compose(&q[a][b]), &F::one());
        }
        sums[b].iota = acc;
    }
}

/// Multiplicities of b_x b_s in the KL basis as (y, degree-of-ι, count).
/// A summand v^k b_y corresponds to B_y{k} ⊂ B_x B_s; ι has degree -k... we return
/// the degree d of ι: BS(y) -> BS(x s) which is -k under ch(B{k}) = v^k ch(B).
pub fn kl_mult(h: &Hecke, x: usize, s: usize) -> Vec<(usize, i32, usize)> {
    let prod = h.mul_bs(&h.kl[x], s);
    let k = h.to_kl(&prod);
    let mut out = vec![];
    for (y, c) in k.iter().enumerate() {
        for (&e, &m) in &c.c {
            assert!(m > 0);
            out.push((y, e, m as usize));
        }
    }
    out
}
