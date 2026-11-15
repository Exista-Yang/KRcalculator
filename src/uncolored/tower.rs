//! Milestone M1b.1: the restriction-tower closure (design "algebraic-reduction").
//!
//! Notation. For the set of open variables S, A_S = F[x_i : i ∈ S]/(x_i^{cap_i+1}),
//! with cap_i = N-1 for an unreduced strand and cap = 0 for the marked strand of the
//! reduced theory. τ_i is the coefficient of x_i^{cap_i}.
//!
//! Objects. An object of a stage is X = A^c / ker(D ⊗ A) ≅ Im(D ⊗ A), where D is an
//! untruncated generator matrix over F[x_i : i ∈ S] (rows: relations, columns: the c
//! coordinates) and every coordinate t carries a q-degree deg[t]. Maps are matrices over
//! A on coordinates (v ↦ F v), homogeneous: deg F[t'][t] = deg_src[t] - deg_tgt[t'].
//!
//! * (P2) Stage n: Φ_N(B_x{q}[h]) ≅ Im(G_x ⊗ A_n), coordinate a of degree q + deg g_a,
//!   and Φ_N(f) = (v ↦ F v) with the polynomial matrix F of f in the R-bases.
//! * (P3) Closing the variable x_v: A_S = ⊕_{j ≤ cap} A_{S'} x_v^j (S' = S \ {v}) and
//!   τ_v is a perfect A_{S'}-valued pairing, so Im(D ⊗ A_S) restricted to A_{S'} is
//!   Im(M ⊗ A_{S'}) with the Hankel matrix M[(b,l)][(a,j)] = [x_v^{cap-l-j}] D_ba; the
//!   coordinate (a,j) (the x_v^j-component of coordinate a) has degree deg[a] + 2j, and a
//!   map F becomes the block-Toeplitz matrix T(F)[(c,j')][(a,j)] = [x_v^{j'-j}] F_ca.
//! * (P4) Unipotent column operations V and row operations U (exact division in the
//!   polynomial ring, so U and V are invertible over A) bring M to a block-diagonal form:
//!   1x1 pivots p (cyclic children Im(p ⊗ A)), junk coordinates (columns ≡ 0 mod I, in the
//!   kernel) and the connected components of the remainder. In the new coordinates
//!   v' = V^{-1} v a map is V_y^{-1} T(F) V_x restricted to the kept coordinates.
//! * (P5) Gaussian elimination in the category of graded A-modules: a component
//!   φ = [F]: X_i → X_j with matching coordinate degrees and invertible degree-0 part F_0
//!   is invertible as a matrix over A, maps ker_i into ker_j and induces a surjection; it
//!   is an isomorphism iff dim X_i = dim X_j, and then φ^{-1} = [F^{-1}].
//! * (P6) At stage 0 every object is a vector space; complete elimination leaves no map
//!   and the survivors are H(Φ_N(C)).
//!
//! Correctness never depends on how well the blocks split: a non-split block simply
//! stays one (larger) object.

use crate::uncolored::engine::Obj;
use crate::uncolored::engine2::{EPoly, TypeBasis};
use crate::algebra::field::Field;
use crate::algebra::laurent::Poincare;
use std::cmp::Reverse;
use std::collections::{BTreeMap, BTreeSet, BinaryHeap, HashMap, HashSet};
use std::rc::Rc;

/// Bits per exponent field of a packed monomial.
const BITS: usize = 6;
const FMASK: u64 = 63;
/// At most 10 variables fit in a u64 with 6-bit fields.
pub const MAX_VARS: usize = 10;
/// Products with at least this many term pairs try dense accumulation (`pmul_dense`).
const DENSE_MIN: usize = 256;

/// Packed polynomial: terms (monomial, coefficient) sorted by the packed key (a lex
/// order with the highest variable most significant), no zero coefficients.
pub type P<F> = Vec<(u64, F)>;

#[inline]
fn fld(m: u64, i: usize) -> u64 {
    (m >> (BITS * i)) & FMASK
}

#[inline]
fn mdeg(mut m: u64) -> u32 {
    let mut s = 0u32;
    while m != 0 {
        s += (m & FMASK) as u32;
        m >>= BITS;
    }
    s
}

/// Pack an exponent vector (every exponent ≤ 63, at most MAX_VARS variables).
pub fn pack(e: &[u32]) -> u64 {
    assert!(e.len() <= MAX_VARS, "tower closure supports at most {} strands", MAX_VARS);
    e.iter().enumerate().fold(0u64, |acc, (i, &x)| {
        assert!(x <= FMASK as u32, "exponent {} too large for packed monomials", x);
        acc | ((x as u64) << (BITS * i))
    })
}

/// Truncation caps: x_i^{cap_i+1} = 0.
#[derive(Clone, Debug)]
pub struct Caps {
    n: usize,
    cap: Vec<u32>,
    addc: u64,
    mask: u64,
}

impl Caps {
    pub fn new(cap: &[u32]) -> Caps {
        let n = cap.len();
        assert!(n <= MAX_VARS, "tower closure supports at most {} strands", MAX_VARS);
        let (mut addc, mut mask) = (0u64, 0u64);
        for (i, &c) in cap.iter().enumerate() {
            assert!(c <= 31, "truncation cap {} too large (N ≤ 32)", c);
            addc |= ((31 - c) as u64) << (BITS * i);
            mask |= 32u64 << (BITS * i);
        }
        Caps { n, cap: cap.to_vec(), addc, mask }
    }
    /// Within the caps, by one add and one mask: field e > cap iff e + (31 - cap) sets
    /// bit 5. Valid when every field e_i ≤ 32 + cap_i (no carry into the next field), in
    /// particular for a product of two truncated monomials (e_i ≤ 2·cap_i, cap_i ≤ 31).
    #[inline]
    fn fits_small(&self, m: u64) -> bool {
        (m + self.addc) & self.mask == 0
    }
    /// Within the caps, for any packed monomial.
    #[inline]
    fn fits(&self, m: u64) -> bool {
        (0..self.n).all(|i| fld(m, i) <= self.cap[i] as u64)
    }
}

// ---------------------------------------------------------------- polynomials

/// Sort, combine equal monomials, drop zeros. Compacts in place and releases the excess
/// capacity: callers pass uncombined product lists (up to |a|·|b| terms), and the result is
/// stored long term in matrices, so keeping that capacity multiplies the memory footprint.
fn norm<F: Field>(mut v: Vec<(u64, F)>) -> P<F> {
    // already strictly increasing (e.g. a single normalized contribution): keep as is
    if v.windows(2).all(|w| w[0].0 < w[1].0) {
        v.retain(|t| !t.1.is_zero());
        if v.capacity() > v.len() + v.len() / 4 + 4 {
            v.shrink_to_fit();
        }
        return v;
    }
    if v.len() >= 64 {
        if let Some(out) = norm_dense(&v) {
            return out;
        }
    }
    norm_sort(v)
}

/// `norm` by sorting.
fn norm_sort<F: Field>(mut v: Vec<(u64, F)>) -> P<F> {
    if v.len() > 1 {
        v.sort_unstable_by_key(|t| t.0);
        let mut k = 0usize;
        for i in 1..v.len() {
            if v[i].0 == v[k].0 {
                let c = v[i].1.clone();
                v[k].1.add_assign(&c);
            } else {
                if !v[k].1.is_zero() {
                    k += 1;
                }
                v.swap(k, i);
            }
        }
        if !v[k].1.is_zero() {
            k += 1;
        }
        v.truncate(k);
    } else {
        v.retain(|t| !t.1.is_zero());
    }
    if v.capacity() > v.len() + v.len() / 4 + 4 {
        v.shrink_to_fit();
    }
    v
}

fn nonzero_mod<F: Field>(p: &P<F>, caps: &Caps) -> bool {
    p.iter().any(|t| caps.fits(t.0))
}

fn trunc<F: Field>(p: &P<F>, caps: &Caps) -> P<F> {
    p.iter().filter(|t| caps.fits(t.0)).cloned().collect()
}

fn pscale<F: Field>(p: &P<F>, c: &F) -> P<F> {
    p.iter().map(|(m, x)| (*m, x.mul(c))).filter(|t| !t.1.is_zero()).collect()
}

/// a + cb · x^shift · b (merge of sorted term lists).
fn axpy<F: Field>(a: &P<F>, b: &P<F>, cb: &F, shift: u64) -> P<F> {
    let mut out = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0usize, 0usize);
    while i < a.len() && j < b.len() {
        let bm = b[j].0 + shift;
        if a[i].0 == bm {
            let c = a[i].1.add(&b[j].1.mul(cb));
            if !c.is_zero() {
                out.push((bm, c));
            }
            i += 1;
            j += 1;
        } else if a[i].0 < bm {
            out.push(a[i].clone());
            i += 1;
        } else {
            let c = b[j].1.mul(cb);
            if !c.is_zero() {
                out.push((bm, c));
            }
            j += 1;
        }
    }
    out.extend(a[i..].iter().cloned());
    for t in &b[j..] {
        let c = t.1.mul(cb);
        if !c.is_zero() {
            out.push((t.0 + shift, c));
        }
    }
    out
}

fn psub<F: Field>(a: &P<F>, b: &P<F>) -> P<F> {
    axpy(a, b, &F::one().neg(), 0)
}

/// Untruncated product of homogeneous polynomials (degree ≤ 63, so no field overflows).
fn pmul<F: Field>(a: &P<F>, b: &P<F>) -> P<F> {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    assert!(mdeg(a[0].0) + mdeg(b[0].0) <= FMASK as u32, "untruncated product of degree > 63");
    if a.len() * b.len() >= DENSE_MIN {
        if let Some(p) = pmul_dense(a, b, None) {
            return p;
        }
    }
    let mut v = Vec::with_capacity(a.len() * b.len());
    for (m1, c1) in a {
        for (m2, c2) in b {
            v.push((m1 + m2, c1.mul(c2)));
        }
    }
    norm(v)
}

/// Number of exponent fields in use by a set of packed monomials (from their bitwise or).
#[inline]
fn fields_used(or: u64) -> usize {
    (64 - or.leading_zeros() as usize).div_ceil(BITS)
}

/// Per-field minimum and maximum exponent of the terms, and whether all terms have the
/// same total degree.
fn exp_bounds<F: Field>(p: &[(u64, F)]) -> ([u32; MAX_VARS], [u32; MAX_VARS], bool) {
    let nf = fields_used(p.iter().fold(0u64, |a, t| a | t.0));
    let mut lo = [0u32; MAX_VARS];
    let mut hi = [0u32; MAX_VARS];
    lo[..nf].fill(u32::MAX);
    let d0 = mdeg(p[0].0);
    let mut hom = true;
    for (m, _) in p {
        let mut d = 0;
        for i in 0..nf {
            let e = fld(*m, i) as u32;
            lo[i] = lo[i].min(e);
            hi[i] = hi[i].max(e);
            d += e;
        }
        hom &= d == d0;
    }
    (lo, hi, hom)
}

/// A box lo ≤ e ≤ hi of exponent vectors addressed densely: fields constant on the box are
/// left out of the index, and when all monomials have the same total degree so is the
/// lowest varying field (it is implied). The index is mixed radix with the highest field
/// most significant, so increasing index is increasing packed order.
struct DenseBox {
    lo: [u32; MAX_VARS],
    hi: [u32; MAX_VARS],
    kept: Vec<usize>,
    stride: [usize; MAX_VARS],
    /// (lowest varying field, total degree) for homogeneous monomials
    dropped: Option<(usize, u32)>,
    /// the constant fields, packed
    base: u64,
    /// total degree of the constant fields
    base_deg: u32,
    size: usize,
}

impl DenseBox {
    /// None if the box has more than `limit` cells.
    fn new(lo: [u32; MAX_VARS], hi: [u32; MAX_VARS], total: Option<u32>, limit: usize) -> Option<DenseBox> {
        let varying: Vec<usize> = (0..MAX_VARS).filter(|&i| hi[i] > lo[i]).collect();
        let dropped = match (total, varying.first()) {
            (Some(t), Some(&i)) => Some((i, t)),
            _ => None,
        };
        let kept: Vec<usize> = varying.iter().copied().filter(|&i| dropped.map(|d| d.0) != Some(i)).collect();
        let mut stride = [0usize; MAX_VARS];
        let mut size = 1usize;
        for &i in &kept {
            stride[i] = size;
            size = size.checked_mul((hi[i] - lo[i] + 1) as usize)?;
            if size > limit {
                return None;
            }
        }
        let fixed = (0..MAX_VARS).filter(|i| !varying.contains(i));
        let base = fixed.clone().map(|i| (lo[i] as u64) << (BITS * i)).sum();
        let base_deg = fixed.map(|i| lo[i]).sum();
        Some(DenseBox { lo, hi, kept, stride, dropped, base, base_deg, size })
    }
    /// Index offset of m relative to the corner `l` (Σ over kept fields).
    #[inline]
    fn index(&self, m: u64, l: &[u32; MAX_VARS]) -> usize {
        self.kept.iter().map(|&i| (fld(m, i) as u32 - l[i]) as usize * self.stride[i]).sum()
    }
    fn decode(&self, idx: usize) -> u64 {
        let mut m = self.base;
        let mut d = self.base_deg;
        for &i in &self.kept {
            let e = self.lo[i] + ((idx / self.stride[i]) % (self.hi[i] - self.lo[i] + 1) as usize) as u32;
            m |= (e as u64) << (BITS * i);
            d += e;
        }
        if let Some((i, t)) = self.dropped {
            m |= ((t - d) as u64) << (BITS * i);
        }
        m
    }
    /// The nonzero cells in increasing packed order.
    fn collect<F: Field>(&self, table: Vec<F>) -> P<F> {
        let nz = table.iter().filter(|c| !c.is_zero()).count();
        let mut out = Vec::with_capacity(nz);
        for (idx, c) in table.into_iter().enumerate() {
            if !c.is_zero() {
                out.push((self.decode(idx), c));
            }
        }
        out
    }
}

/// Combine a long list of terms by direct addressing in the box they span (None if that
/// box is large compared with the list; `norm` then sorts).
fn norm_dense<F: Field>(v: &[(u64, F)]) -> Option<P<F>> {
    let (lo, hi, hom) = exp_bounds(v);
    let bx = DenseBox::new(lo, hi, if hom { Some(mdeg(v[0].0)) } else { None }, 2 * v.len() + 64)?;
    let mut table: Vec<F> = vec![F::zero(); bx.size];
    for (m, c) in v {
        table[bx.index(*m, &lo)].add_assign(c);
    }
    Some(bx.collect(table))
}

/// Product a·b (truncated to the caps if `caps` is given; a, b then truncated) by direct
/// addressing in the box of exponents the product can have, without building and sorting
/// all |a|·|b| products. None if that box is large compared with |a|·|b| (the caller then
/// sorts).
fn pmul_dense<F: Field>(a: &P<F>, b: &P<F>, caps: Option<&Caps>) -> Option<P<F>> {
    if a.is_empty() || b.is_empty() {
        return Some(vec![]);
    }
    let (la, ha, hom_a) = exp_bounds(a);
    let (lb, hb, hom_b) = exp_bounds(b);
    let mut lo = [0u32; MAX_VARS];
    let mut hi = [0u32; MAX_VARS];
    for i in 0..MAX_VARS {
        lo[i] = la[i] + lb[i];
        hi[i] = ha[i] + hb[i];
        if let Some(c) = caps {
            if i < c.n {
                if lo[i] > c.cap[i] {
                    return Some(vec![]);
                }
                hi[i] = hi[i].min(c.cap[i]);
            }
        }
    }
    let total = if hom_a && hom_b { Some(mdeg(a[0].0) + mdeg(b[0].0)) } else { None };
    let bx = DenseBox::new(lo, hi, total, 2 * a.len() * b.len() + 256)?;
    let ia: Vec<usize> = a.iter().map(|t| bx.index(t.0, &la)).collect();
    let ib: Vec<usize> = b.iter().map(|t| bx.index(t.0, &lb)).collect();
    let mut table: Vec<F> = vec![F::zero(); bx.size];
    match caps {
        None => {
            for (x, (_, c1)) in ia.iter().zip(a) {
                for (y, (_, c2)) in ib.iter().zip(b) {
                    table[x + y].fma(c1, c2);
                }
            }
        }
        Some(c) => {
            for (x, (m1, c1)) in ia.iter().zip(a) {
                for (y, (m2, c2)) in ib.iter().zip(b) {
                    if c.fits_small(m1 + m2) {
                        table[x + y].fma(c1, c2);
                    }
                }
            }
        }
    }
    Some(bx.collect(table))
}

/// a - q·b, untruncated (homogeneous q, b of total degree ≤ 63). Adding a fixed monomial
/// to every key preserves the order of packed monomials, so q·b is a merge of |q| shifted
/// copies of b: for a short factor this is a few linear merges instead of building and
/// sorting all |q|·|b| products.
fn sub_mul<F: Field>(a: &P<F>, q: &P<F>, b: &P<F>) -> P<F> {
    if q.is_empty() || b.is_empty() {
        return a.clone();
    }
    assert!(mdeg(q[0].0) + mdeg(b[0].0) <= FMASK as u32, "untruncated product of degree > 63");
    let (short, long) = if q.len() <= b.len() { (q, b) } else { (b, q) };
    if short.len() <= 4 {
        let mut acc = axpy(a, long, &short[0].1.neg(), short[0].0);
        for (m, c) in &short[1..] {
            acc = axpy(&acc, long, &c.neg(), *m);
        }
        acc
    } else {
        psub(a, &pmul(q, b))
    }
}

/// Push the terms of the truncated product a·b (a, b truncated) onto `out` (unnormalized).
#[inline]
fn pmul_trunc_into<F: Field>(a: &P<F>, b: &P<F>, caps: &Caps, out: &mut Vec<(u64, F)>) {
    for (m1, c1) in a {
        for (m2, c2) in b {
            let m = m1 + m2;
            if caps.fits_small(m) {
                out.push((m, c1.mul(c2)));
            }
        }
    }
}

fn pmul_trunc<F: Field>(a: &P<F>, b: &P<F>, caps: &Caps) -> P<F> {
    if a.len() * b.len() >= DENSE_MIN {
        if let Some(p) = pmul_dense(a, b, Some(caps)) {
            return p;
        }
    }
    let mut v = Vec::with_capacity(a.len() * b.len());
    pmul_trunc_into(a, b, caps, &mut v);
    norm(v)
}

#[inline]
fn mono_divides(d: u64, m: u64, n: usize) -> bool {
    (0..n).all(|i| fld(m, i) >= fld(d, i))
}

/// Exact division in the polynomial ring: Some(q) with a = q·p, or None.
fn pdiv_exact<F: Field>(a: &P<F>, p: &P<F>, n: usize) -> Option<P<F>> {
    if a.is_empty() {
        return Some(vec![]);
    }
    let (lm, lc) = p.last().expect("division by zero polynomial");
    let lci = lc.inv();
    if p.len() == 1 {
        let mut q = Vec::with_capacity(a.len());
        for (m, c) in a {
            if !mono_divides(*lm, *m, n) {
                return None;
            }
            q.push((m - lm, c.mul(&lci)));
        }
        return Some(q);
    }
    if mdeg(a[0].0) < mdeg(*lm) {
        return None;
    }
    let mut r: P<F> = a.clone();
    let mut q: Vec<(u64, F)> = vec![];
    while let Some((rm, rc)) = r.last().cloned() {
        if !mono_divides(*lm, rm, n) {
            return None;
        }
        let qm = rm - lm;
        let qc = rc.mul(&lci);
        r = axpy(&r, p, &qc.neg(), qm);
        q.push((qm, qc));
    }
    q.reverse();
    Some(q)
}

/// Coefficients of x_v^e for e = 0..=cap (field v cleared).
fn split_var<F: Field>(p: &P<F>, v: usize, cap: usize) -> Vec<P<F>> {
    let mut out: Vec<P<F>> = vec![vec![]; cap + 1];
    for (m, c) in p {
        let e = fld(*m, v) as usize;
        if e <= cap {
            out[e].push((m - ((e as u64) << (BITS * v)), c.clone()));
        }
    }
    out
}

fn const_term<F: Field>(p: &P<F>) -> F {
    match p.first() {
        Some((0, c)) => c.clone(),
        _ => F::zero(),
    }
}

fn monic<F: Field>(p: &P<F>) -> P<F> {
    let lci = p.last().unwrap().1.inv();
    pscale(p, &lci)
}

// ---------------------------------------------------------------- matrices

/// Sparse polynomial matrix: rows[r] = [(column, entry)] sorted by column, no zero entries.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct PMat<F: Field> {
    pub nr: usize,
    pub nc: usize,
    pub rows: Vec<Vec<(u32, P<F>)>>,
}

impl<F: Field> PMat<F> {
    fn is_zero(&self) -> bool {
        self.rows.iter().all(|r| r.is_empty())
    }
    fn nterms(&self) -> usize {
        self.rows.iter().map(|r| r.iter().map(|e| e.1.len()).sum::<usize>()).sum()
    }
    fn from_acc(nr: usize, nc: usize, acc: BTreeMap<(u32, u32), Vec<(u64, F)>>) -> PMat<F> {
        let mut rows: Vec<Vec<(u32, P<F>)>> = vec![vec![]; nr];
        for ((r, c), raw) in acc {
            let p = norm(raw);
            if !p.is_empty() {
                rows[r as usize].push((c, p));
            }
        }
        PMat { nr, nc, rows }
    }
}

/// Truncated product a·b over A (entries of a and b truncated).
fn mat_mul<F: Field>(a: &PMat<F>, b: &PMat<F>, caps: &Caps) -> PMat<F> {
    assert_eq!(a.nc, b.nr);
    let mut rows = Vec::with_capacity(a.nr);
    for row in &a.rows {
        let mut acc: BTreeMap<u32, Vec<(u64, F)>> = BTreeMap::new();
        for (k, pa) in row {
            for (j, pb) in &b.rows[*k as usize] {
                pmul_trunc_into(pa, pb, caps, acc.entry(*j).or_default());
            }
        }
        rows.push(acc.into_iter().map(|(j, raw)| (j, norm(raw))).filter(|e| !e.1.is_empty()).collect());
    }
    PMat { nr: a.nr, nc: b.nc, rows }
}

/// a + c·b.
fn mat_axpy<F: Field>(a: &PMat<F>, b: &PMat<F>, c: &F) -> PMat<F> {
    assert!(a.nr == b.nr && a.nc == b.nc);
    let mut rows = Vec::with_capacity(a.nr);
    for (ra, rb) in a.rows.iter().zip(b.rows.iter()) {
        let mut out: Vec<(u32, P<F>)> = Vec::with_capacity(ra.len() + rb.len());
        let (mut i, mut j) = (0usize, 0usize);
        while i < ra.len() || j < rb.len() {
            if j == rb.len() || (i < ra.len() && ra[i].0 < rb[j].0) {
                out.push(ra[i].clone());
                i += 1;
            } else if i == ra.len() || rb[j].0 < ra[i].0 {
                let p = pscale(&rb[j].1, c);
                if !p.is_empty() {
                    out.push((rb[j].0, p));
                }
                j += 1;
            } else {
                let p = axpy(&ra[i].1, &rb[j].1, c, 0);
                if !p.is_empty() {
                    out.push((ra[i].0, p));
                }
                i += 1;
                j += 1;
            }
        }
        rows.push(out);
    }
    PMat { nr: a.nr, nc: a.nc, rows }
}

/// Inverse of a dense constant matrix, None if singular.
fn dense_inverse<F: Field>(m: &[Vec<F>]) -> Option<Vec<Vec<F>>> {
    let n = m.len();
    let mut a: Vec<Vec<F>> = m.iter().enumerate().map(|(i, r)| { let mut r = r.clone(); r.extend((0..n).map(|j| if i == j { F::one() } else { F::zero() })); r }).collect();
    for col in 0..n {
        let piv = (col..n).find(|&r| !a[r][col].is_zero())?;
        a.swap(col, piv);
        let inv = a[col][col].inv();
        for x in a[col].iter_mut() {
            *x = x.mul(&inv);
        }
        for r in 0..n {
            if r != col && !a[r][col].is_zero() {
                let f = a[r][col].clone();
                let pr = a[col].clone();
                for (x, y) in a[r].iter_mut().zip(pr.iter()) {
                    x.sub_assign(&y.mul(&f));
                }
            }
        }
    }
    Some(a.into_iter().map(|r| r[n..].to_vec()).collect())
}

/// Rank of a constant matrix (all entries constants).
fn const_rank<F: Field>(d: &PMat<F>) -> usize {
    let mut rows: Vec<Vec<F>> = d.rows.iter().map(|r| { let mut v = vec![F::zero(); d.nc]; for (c, p) in r { v[*c as usize] = const_term(p); } v }).collect();
    let mut rank = 0;
    for col in 0..d.nc {
        let Some(piv) = (rank..rows.len()).find(|&r| !rows[r][col].is_zero()) else { continue };
        rows.swap(rank, piv);
        let inv = rows[rank][col].inv();
        let pr: Vec<F> = rows[rank].iter().map(|x| x.mul(&inv)).collect();
        for r in rank + 1..rows.len() {
            if !rows[r][col].is_zero() {
                let f = rows[r][col].clone();
                for (x, y) in rows[r].iter_mut().zip(pr.iter()) {
                    x.sub_assign(&y.mul(&f));
                }
            }
        }
        rank += 1;
    }
    rank
}

// ---------------------------------------------------------------- tower

/// Which end variable is closed at each stage.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EndPolicy {
    /// x_n, x_{n-1}, ..., x_1
    Right,
    /// x_1, x_2, ..., x_n
    Left,
    /// x_n, x_1, x_{n-1}, x_2, ...
    Alternate,
    /// Diagnostic only: middle strand first, then outward (x_{m}, x_{m-1}, x_{m+1}, ...).
    /// Valid by (P3) but non-planar: blocks split much worse than at an end.
    Middle,
}

impl EndPolicy {
    pub fn parse(s: &str) -> Option<EndPolicy> {
        match s {
            "right" | "r" => Some(EndPolicy::Right),
            "left" | "l" => Some(EndPolicy::Left),
            "alt" | "alternate" | "a" => Some(EndPolicy::Alternate),
            "middle" | "m" => Some(EndPolicy::Middle),
            _ => None,
        }
    }
    /// Variables in closing order; at every stage the closed one is an end of the
    /// interval of open strands.
    pub fn order(&self, n: usize) -> Vec<usize> {
        match self {
            EndPolicy::Right => (0..n).rev().collect(),
            EndPolicy::Left => (0..n).collect(),
            EndPolicy::Alternate => {
                let (mut lo, mut hi) = (0usize, n);
                let mut v = vec![];
                let mut right = true;
                while lo < hi {
                    if right {
                        hi -= 1;
                        v.push(hi);
                    } else {
                        v.push(lo);
                        lo += 1;
                    }
                    right = !right;
                }
                v
            }
            EndPolicy::Middle => {
                let m = n / 2;
                let mut v = vec![m];
                for d in 1..=n {
                    if m >= d {
                        v.push(m - d);
                    }
                    if m + d < n {
                        v.push(m + d);
                    }
                }
                v
            }
        }
    }
}

/// Decomposition of Im(M ⊗ A_{S'}) for the Hankel matrix M of one generator matrix D.
/// Restricted coordinates are indexed t = a·w + j (w = cap_v + 1).
struct Dec<F: Field> {
    /// children: interned generator matrix and its coordinates (restricted indices, in order)
    children: Vec<(Rc<PMat<F>>, Vec<u32>)>,
    /// restricted coordinate (new basis) -> (child, position); None = junk (in the kernel)
    coord: Vec<Option<(u32, u32)>>,
    /// kept t: column t of V (old restricted coordinate s -> entry), truncated
    vcol: Vec<Vec<(u32, P<F>)>>,
    /// old restricted coordinate s -> [(kept t, V^{-1}[t][s])], truncated
    vinv_by_s: Vec<Vec<(u32, P<F>)>>,
    nonsplit: usize,
}

#[derive(Clone)]
struct TObj<F: Field> {
    h: i32,
    deg: Vec<i32>,
    sdeg: Vec<i32>,
    d: Rc<PMat<F>>,
}

struct Edge<F: Field> {
    m: PMat<F>,
    ver: u64,
}

struct Cx<F: Field> {
    obj: Vec<TObj<F>>,
    alive: Vec<bool>,
    out: Vec<BTreeMap<u32, Edge<F>>>,
    inn: Vec<BTreeSet<u32>>,
}

impl<F: Field> Cx<F> {
    fn n_alive(&self) -> usize {
        self.alive.iter().filter(|&&a| a).count()
    }
    fn n_edges(&self) -> usize {
        self.out.iter().map(|o| o.len()).sum()
    }
    fn n_terms(&self) -> usize {
        self.out.iter().map(|o| o.values().map(|e| e.m.nterms()).sum::<usize>()).sum()
    }
    fn n_coords(&self) -> usize {
        (0..self.obj.len()).filter(|&i| self.alive[i]).map(|i| self.obj[i].deg.len()).sum()
    }
}

/// Per-stage statistics.
#[derive(Clone, Debug, Default)]
pub struct StageStats {
    pub var: usize,
    pub objects: usize,
    pub edges: usize,
    pub terms: usize,
    pub nonsplit: usize,
    pub after_objects: usize,
    pub after_coords: usize,
    pub after_edges: usize,
    pub after_terms: usize,
    /// eliminations: total, of non-cyclic objects (more than one coordinate), with a
    /// non-constant inverse (Neumann series), with dimensions certified by `gdim`
    /// (different generator matrices)
    pub elims: usize,
    pub elims_noncyclic: usize,
    pub elims_neumann: usize,
    pub elims_gdim: usize,
    pub secs: f64,
}

/// Options of the tower closure.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TowerOpts {
    pub policy: EndPolicy,
    /// Stress/diagnostic mode: only constant (unit) pivots, so that types which would
    /// split by exact division stay non-split blocks. Results must not change.
    pub unit_pivots_only: bool,
    /// Also run Gaussian elimination at stage n (before closing any variable). A minimal
    /// complex has no isomorphism components over A_n, so this only helps non-minimal
    /// input (synthetic complexes, tests).
    pub eliminate_top: bool,
}

impl TowerOpts {
    pub fn new(policy: EndPolicy) -> TowerOpts {
        TowerOpts { policy, unit_pivots_only: false, eliminate_top: false }
    }
}

pub struct Tower<F: Field> {
    n: usize,
    caps: Caps,
    order: Vec<usize>,
    unit_only: bool,
    counts: [usize; 4],
    intern: HashSet<Rc<PMat<F>>>,
    dec_cache: HashMap<(usize, usize), Rc<Dec<F>>>,
    gdim_cache: HashMap<(usize, u32), usize>,
    ver: u64,
}

impl<F: Field> Tower<F> {
    pub fn new(caps: &[u32], opts: TowerOpts) -> Tower<F> {
        let n = caps.len();
        Tower { n, caps: Caps::new(caps), order: opts.policy.order(n), unit_only: opts.unit_pivots_only, counts: [0; 4], intern: HashSet::new(), dec_cache: HashMap::new(), gdim_cache: HashMap::new(), ver: 0 }
    }

    fn intern(&mut self, m: PMat<F>) -> Rc<PMat<F>> {
        if let Some(r) = self.intern.get(&m) {
            return r.clone();
        }
        let r = Rc::new(m);
        self.intern.insert(r.clone());
        r
    }

    /// Interned generator matrix from Gram polynomials g[b][a] (exponent vectors).
    fn gram_matrix(&mut self, g: &[Vec<EPoly<F>>], ncols: usize) -> Rc<PMat<F>> {
        let mut rows = vec![];
        for row in g {
            let mut r = vec![];
            for (a, p) in row.iter().enumerate() {
                let mut terms: Vec<(u64, F)> = vec![];
                for (e, c) in p {
                    assert!(e.len() == self.n);
                    let total: u32 = e.iter().sum();
                    assert!(total <= FMASK as u32, "Gram entry of degree {} > 63", total);
                    terms.push((pack(e), c.clone()));
                }
                let p = norm(terms);
                if nonzero_mod(&p, &self.caps) {
                    r.push((a as u32, p));
                }
            }
            rows.push(r);
        }
        let mut m = PMat { nr: g.len(), nc: ncols, rows };
        normalize_rows(&mut m);
        self.intern(m)
    }

    fn decompose(&mut self, d: &Rc<PMat<F>>, v: usize) -> Rc<Dec<F>> {
        let key = (Rc::as_ptr(d) as usize, v);
        if let Some(x) = self.dec_cache.get(&key) {
            return x.clone();
        }
        let dec = Rc::new(self.decompose_raw(d, v));
        self.dec_cache.insert(key, dec.clone());
        dec
    }

    /// (P3)+(P4): Hankel matrix in x_v, 1x1 exact-division pivots, junk, components.
    fn decompose_raw(&mut self, d: &PMat<F>, v: usize) -> Dec<F> {
        let n = self.n;
        let caps = self.caps.clone();
        let cap = caps.cap[v] as usize;
        let w = cap + 1;
        let nr = d.nr * w;
        let nc = d.nc * w;
        let mut rows: Vec<BTreeMap<u32, P<F>>> = vec![BTreeMap::new(); nr];
        let mut cols: Vec<BTreeSet<u32>> = vec![BTreeSet::new(); nc];
        for (b, row) in d.rows.iter().enumerate() {
            for (a, poly) in row {
                let split = split_var(poly, v, cap);
                for (e, coef) in split.iter().enumerate() {
                    if coef.is_empty() || !nonzero_mod(coef, &caps) {
                        continue;
                    }
                    for l in 0..=(cap - e) {
                        let j = cap - e - l;
                        let u = (b * w + l) as u32;
                        let t = (*a as usize * w + j) as u32;
                        rows[u as usize].insert(t, coef.clone());
                        cols[t as usize].insert(u);
                    }
                }
            }
        }
        let mut row_act = vec![true; nr];
        let mut col_act = vec![true; nc];
        // V (columns) and V^{-1} (rows), both start as the identity
        let one: P<F> = vec![(0u64, F::one())];
        let mut vcols: Vec<BTreeMap<u32, P<F>>> = (0..nc).map(|t| BTreeMap::from([(t as u32, one.clone())])).collect();
        let mut vinv: Vec<BTreeMap<u32, P<F>>> = vcols.clone();
        let mut cyclic: Vec<(u32, P<F>)> = vec![];
        loop {
            let mut cands: Vec<(u32, usize, usize, u32, u32)> = vec![];
            for u in 0..nr {
                if !row_act[u] {
                    continue;
                }
                for (&t, p) in &rows[u] {
                    let mk = (rows[u].len() - 1) * (cols[t as usize].len() - 1);
                    cands.push((mdeg(p[0].0), p.len(), mk, u as u32, t));
                }
            }
            if cands.is_empty() {
                break;
            }
            cands.sort_unstable();
            let mut progress = false;
            for &(_, _, _, u, t) in &cands {
                let (uu, tt) = (u as usize, t as usize);
                if !row_act[uu] || !col_act[tt] {
                    continue;
                }
                let Some(p) = rows[uu].get(&t).cloned() else { continue };
                let unit = p.len() == 1 && p[0].0 == 0;
                if !unit && self.unit_only {
                    continue;
                }
                let mut qs: Vec<(u32, P<F>)> = vec![];
                let mut ok = true;
                if unit {
                    let ci = p[0].1.inv();
                    for (&t2, e) in &rows[uu] {
                        if t2 != t {
                            qs.push((t2, pscale(e, &ci)));
                        }
                    }
                } else {
                    for &u2 in &cols[tt] {
                        if u2 != u && pdiv_exact(&rows[u2 as usize][&t], &p, n).is_none() {
                            ok = false;
                            break;
                        }
                    }
                    if ok {
                        for (&t2, e) in &rows[uu] {
                            if t2 != t {
                                match pdiv_exact(e, &p, n) {
                                    Some(q) => qs.push((t2, q)),
                                    None => {
                                        ok = false;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
                if !ok {
                    continue;
                }
                // column ops col_t2 -= q·col_t (Schur update on the other rows of column t),
                // then row ops clear column t.
                let colents: Vec<(u32, P<F>)> = cols[tt].iter().filter(|&&u2| u2 != u).map(|&u2| (u2, rows[u2 as usize][&t].clone())).collect();
                for (u2, c) in &colents {
                    let r2 = &mut rows[*u2 as usize];
                    for (t2, q) in &qs {
                        let cur = r2.remove(t2).unwrap_or_default();
                        let new = sub_mul(&cur, q, c);
                        if !new.is_empty() && nonzero_mod(&new, &caps) {
                            r2.insert(*t2, new);
                            cols[*t2 as usize].insert(*u2);
                        } else {
                            cols[*t2 as usize].remove(u2);
                        }
                    }
                    r2.remove(&t);
                }
                for (t2, _) in &qs {
                    cols[*t2 as usize].remove(&u);
                }
                rows[uu].clear();
                cols[tt].clear();
                for (t2, q) in &qs {
                    let qt = trunc(q, &caps);
                    if qt.is_empty() {
                        continue;
                    }
                    let vt = vcols[tt].clone();
                    map_axpy(&mut vcols[*t2 as usize], &vt, &qt, &F::one().neg(), &caps);
                    let vt2 = vinv[*t2 as usize].clone();
                    map_axpy(&mut vinv[tt], &vt2, &qt, &F::one(), &caps);
                }
                row_act[uu] = false;
                col_act[tt] = false;
                cyclic.push((t, monic(&p)));
                progress = true;
            }
            if !progress {
                break;
            }
        }
        // remainder: connected components of the active entries; empty columns are junk
        let mut parent: Vec<usize> = (0..nc).collect();
        fn find(p: &mut Vec<usize>, x: usize) -> usize {
            let mut r = x;
            while p[r] != r {
                r = p[r];
            }
            let mut y = x;
            while p[y] != r {
                let z = p[y];
                p[y] = r;
                y = z;
            }
            r
        }
        for u in 0..nr {
            if !row_act[u] || rows[u].is_empty() {
                continue;
            }
            let first = *rows[u].keys().next().unwrap() as usize;
            let keys: Vec<u32> = rows[u].keys().copied().collect();
            for t in keys {
                let (a, b) = (find(&mut parent, first), find(&mut parent, t as usize));
                if a != b {
                    parent[a] = b;
                }
            }
        }
        let mut comp_cols: BTreeMap<usize, Vec<u32>> = BTreeMap::new();
        for t in 0..nc {
            if col_act[t] && !cols[t].is_empty() {
                let r = find(&mut parent, t);
                comp_cols.entry(r).or_default().push(t as u32);
            }
        }
        let mut comp_rows: BTreeMap<usize, Vec<u32>> = BTreeMap::new();
        for u in 0..nr {
            if row_act[u] && !rows[u].is_empty() {
                let first = *rows[u].keys().next().unwrap() as usize;
                let r = find(&mut parent, first);
                comp_rows.entry(r).or_default().push(u as u32);
            }
        }
        // children
        let mut children: Vec<(Rc<PMat<F>>, Vec<u32>)> = vec![];
        let mut coord: Vec<Option<(u32, u32)>> = vec![None; nc];
        for (t, p) in cyclic {
            let m = PMat { nr: 1, nc: 1, rows: vec![vec![(0u32, p)]] };
            coord[t as usize] = Some((children.len() as u32, 0));
            let rc = self.intern(m);
            children.push((rc, vec![t]));
        }
        let mut nonsplit = 0;
        for (r, cc) in comp_cols {
            let pos: HashMap<u32, u32> = cc.iter().enumerate().map(|(i, &t)| (t, i as u32)).collect();
            let rr = comp_rows.remove(&r).unwrap_or_default();
            let mrows: Vec<Vec<(u32, P<F>)>> = rr.iter().map(|&u| rows[u as usize].iter().map(|(t, p)| (pos[t], p.clone())).collect()).collect();
            let mut m = PMat { nr: rr.len(), nc: cc.len(), rows: mrows };
            normalize_rows(&mut m);
            for (i, &t) in cc.iter().enumerate() {
                coord[t as usize] = Some((children.len() as u32, i as u32));
            }
            nonsplit += cc.len();
            let rc = self.intern(m);
            children.push((rc, cc));
        }
        let mut vcol: Vec<Vec<(u32, P<F>)>> = vec![vec![]; nc];
        let mut vinv_by_s: Vec<Vec<(u32, P<F>)>> = vec![vec![]; nc];
        for t in 0..nc {
            if coord[t].is_some() {
                vcol[t] = std::mem::take(&mut vcols[t]).into_iter().collect();
                for (s, p) in std::mem::take(&mut vinv[t]) {
                    vinv_by_s[s as usize].push((t as u32, p));
                }
            }
        }
        Dec { children, coord, vcol, vinv_by_s, nonsplit }
    }

    /// Dimension of Im(D ⊗ A) over F at stage s (the variables order[..s] are closed).
    fn gdim(&mut self, d: &Rc<PMat<F>>, s: usize) -> usize {
        let mask = self.order[..s].iter().fold(0u32, |m, &v| m | (1 << v));
        self.gdim_mask(d, mask)
    }

    /// Dimension of Im(D ⊗ A_S) over F, S = the variables not in `closed` (bit mask).
    ///
    /// The dimension does not depend on the order in which the open variables are
    /// restricted, so this does not follow the tower's order through non-split blocks
    /// (whose recursive decomposition can cost far more than the tower itself, and is
    /// wasted when the certificate fails): it takes the tower's next variable if that
    /// splits D completely, else the first of the two ends of the open interval that does,
    /// else the candidate with the fewest non-split coordinates. Constant matrices and
    /// 1x1 monomials are counted directly.
    fn gdim_mask(&mut self, d: &Rc<PMat<F>>, closed: u32) -> usize {
        let key = (Rc::as_ptr(d) as usize, closed);
        if let Some(&x) = self.gdim_cache.get(&key) {
            return x;
        }
        let open: Vec<usize> = (0..self.n).filter(|&v| closed & (1 << v) == 0).collect();
        let r = if let Some(r) = self.gdim_direct(d, &open) {
            r
        } else {
            let mut cands: Vec<usize> = vec![];
            for v in [self.order.iter().copied().find(|&v| closed & (1 << v) == 0).unwrap(), open[open.len() - 1], open[0]] {
                if !cands.contains(&v) {
                    cands.push(v);
                }
            }
            // a cached complete splitting first, then compute in candidate order
            let mut pick: Option<(usize, Rc<Dec<F>>)> = cands.iter().find_map(|&v| {
                self.dec_cache.get(&(Rc::as_ptr(d) as usize, v)).filter(|dec| dec.nonsplit == 0).map(|dec| (v, dec.clone()))
            });
            if pick.is_none() {
                for &v in &cands {
                    let dec = self.decompose(d, v);
                    let better = match &pick {
                        None => true,
                        Some((_, best)) => dec.nonsplit < best.nonsplit,
                    };
                    let done = dec.nonsplit == 0;
                    if better {
                        pick = Some((v, dec));
                    }
                    if done {
                        break;
                    }
                }
            }
            let (v, dec) = pick.unwrap();
            let mut sum = 0;
            for (cd, _) in dec.children.iter() {
                sum += self.gdim_mask(cd, closed | (1 << v));
            }
            sum
        };
        self.gdim_cache.insert(key, r);
        r
    }

    /// dim Im(D ⊗ A_S) without decomposing, when D is constant (rank(D)·dim A_S) or a 1x1
    /// monomial x^e (Π_{i ∈ S} max(cap_i + 1 - e_i, 0)).
    fn gdim_direct(&self, d: &PMat<F>, open: &[usize]) -> Option<usize> {
        let dim_a: usize = open.iter().map(|&v| self.caps.cap[v] as usize + 1).product();
        if d.rows.iter().all(|r| r.iter().all(|(_, p)| p.len() == 1 && p[0].0 == 0)) {
            return Some(const_rank(d) * dim_a);
        }
        if d.nr == 1 && d.nc == 1 && d.rows[0].len() == 1 && d.rows[0][0].1.len() == 1 {
            let m = d.rows[0][0].1[0].0;
            return Some(open.iter().map(|&v| (self.caps.cap[v] as usize + 1).saturating_sub(fld(m, v) as usize)).product());
        }
        None
    }

    /// Restrict the complex from stage s (order[..s] closed) by the variable order[s].
    fn restrict(&mut self, cx: Cx<F>, s: usize, st: &mut StageStats) -> Cx<F> {
        let v = self.order[s];
        let caps = self.caps.clone();
        let cap = caps.cap[v] as usize;
        let w = cap + 1;
        let mut obj: Vec<TObj<F>> = vec![];
        let mut child_ids: Vec<Vec<u32>> = vec![vec![]; cx.obj.len()];
        let mut decs: Vec<Option<Rc<Dec<F>>>> = vec![None; cx.obj.len()];
        for i in 0..cx.obj.len() {
            if !cx.alive[i] {
                continue;
            }
            let o = &cx.obj[i];
            let dec = self.decompose(&o.d, v);
            st.nonsplit += dec.nonsplit;
            for (cd, coords) in &dec.children {
                let deg: Vec<i32> = coords.iter().map(|&t| o.deg[t as usize / w] + 2 * (t as usize % w) as i32).collect();
                let mut sdeg = deg.clone();
                sdeg.sort();
                child_ids[i].push(obj.len() as u32);
                obj.push(TObj { h: o.h, deg, sdeg, d: cd.clone() });
            }
            decs[i] = Some(dec);
        }
        let no = obj.len();
        let mut out: Vec<BTreeMap<u32, Edge<F>>> = (0..no).map(|_| BTreeMap::new()).collect();
        let mut inn: Vec<BTreeSet<u32>> = vec![BTreeSet::new(); no];
        for i in 0..cx.obj.len() {
            if !cx.alive[i] {
                continue;
            }
            let dx = decs[i].as_ref().unwrap();
            for (&j, e) in &cx.out[i] {
                let dy = decs[j as usize].as_ref().unwrap();
                for ((alpha, beta), m) in restrict_map(&e.m, dx, dy, v, w, &caps) {
                    let (a, b) = (child_ids[i][alpha as usize], child_ids[j as usize][beta as usize]);
                    out[a as usize].insert(b, Edge { m, ver: 0 });
                    inn[b as usize].insert(a);
                }
            }
        }
        Cx { obj, alive: vec![true; no], out, inn }
    }

    /// (P5) Inverse of the component i -> j if it is an isomorphism.
    fn iso_inverse(&mut self, cx: &Cx<F>, i: usize, j: usize, s: usize) -> Option<PMat<F>> {
        let (oi, oj) = (&cx.obj[i], &cx.obj[j]);
        if oi.sdeg != oj.sdeg {
            return None;
        }
        let f = &cx.out[i][&(j as u32)].m;
        let k = oi.deg.len();
        let (inv, pure) = inverse_over_a(f, &oi.deg, &oj.deg, &self.caps)?;
        let by_gdim = !Rc::ptr_eq(&oi.d, &oj.d);
        if by_gdim {
            let (di, dj) = (oi.d.clone(), oj.d.clone());
            if self.gdim(&di, s) != self.gdim(&dj, s) {
                return None;
            }
        }
        self.counts[0] += 1;
        if k > 1 {
            self.counts[1] += 1;
        }
        if !pure {
            self.counts[2] += 1;
        }
        if by_gdim {
            self.counts[3] += 1;
        }
        Some(inv)
    }

    /// Gaussian elimination (Markowitz order with lazy cost updates) at stage s.
    fn gauss(&mut self, cx: &mut Cx<F>, s: usize) -> usize {
        let cost = |cx: &Cx<F>, i: usize, j: usize| -> u64 { ((cx.inn[j].len() - 1) as u64) * ((cx.out[i].len() - 1) as u64) };
        let plausible = |cx: &Cx<F>, i: usize, j: usize| -> bool { cx.obj[i].sdeg == cx.obj[j].sdeg };
        let mut heap: BinaryHeap<Reverse<(u64, usize, u32, u32, u64)>> = BinaryHeap::new();
        for i in 0..cx.obj.len() {
            if !cx.alive[i] {
                continue;
            }
            for (&j, e) in &cx.out[i] {
                if plausible(cx, i, j as usize) {
                    heap.push(Reverse((cost(cx, i, j as usize), cx.obj[i].deg.len(), i as u32, j, e.ver)));
                }
            }
        }
        let mut nelim = 0;
        while let Some(Reverse((c0, _, i, j, ver))) = heap.pop() {
            let (iu, ju) = (i as usize, j as usize);
            if !cx.alive[iu] || !cx.alive[ju] {
                continue;
            }
            match cx.out[iu].get(&j) {
                Some(e) if e.ver == ver => {}
                _ => continue,
            }
            let c = cost(cx, iu, ju);
            if c > c0 {
                heap.push(Reverse((c, cx.obj[iu].deg.len(), i, j, ver)));
                continue;
            }
            let Some(inv) = self.iso_inverse(cx, iu, ju, s) else { continue };
            for (a, b) in self.eliminate(cx, iu, ju, &inv) {
                if plausible(cx, a, b) {
                    let e = &cx.out[a][&(b as u32)];
                    heap.push(Reverse((cost(cx, a, b), cx.obj[a].deg.len(), a as u32, b as u32, e.ver)));
                }
            }
            nelim += 1;
        }
        nelim
    }

    /// Cancel the isomorphism i -> j: d_ba -= d_bi φ^{-1} d_ja. Returns the modified edges.
    fn eliminate(&mut self, cx: &mut Cx<F>, i: usize, j: usize, inv: &PMat<F>) -> Vec<(usize, usize)> {
        let caps = self.caps.clone();
        let ins_j: Vec<u32> = cx.inn[j].iter().copied().filter(|&a| a as usize != i).collect();
        let outs_i: Vec<(u32, PMat<F>)> = cx.out[i].iter().filter(|(b, _)| **b as usize != j).map(|(b, e)| (*b, e.m.clone())).collect();
        let mut modified = vec![];
        for &a in &ins_j {
            let y = mat_mul(inv, &cx.out[a as usize][&(j as u32)].m, &caps);
            if y.is_zero() {
                continue;
            }
            for (b, dbi) in &outs_i {
                let delta = mat_mul(dbi, &y, &caps);
                if delta.is_zero() {
                    continue;
                }
                self.ver += 1;
                let ver = self.ver;
                let outa = &mut cx.out[a as usize];
                let new = match outa.get(b) {
                    Some(e) => mat_axpy(&e.m, &delta, &F::one().neg()),
                    None => mat_axpy(&PMat { nr: delta.nr, nc: delta.nc, rows: vec![vec![]; delta.nr] }, &delta, &F::one().neg()),
                };
                if new.is_zero() {
                    outa.remove(b);
                    cx.inn[*b as usize].remove(&a);
                } else {
                    outa.insert(*b, Edge { m: new, ver });
                    cx.inn[*b as usize].insert(a);
                    modified.push((a as usize, *b as usize));
                }
            }
        }
        for x in [i, j] {
            let outs: Vec<u32> = cx.out[x].keys().copied().collect();
            for b in outs {
                cx.inn[b as usize].remove(&(x as u32));
            }
            cx.out[x].clear();
            let ins: Vec<u32> = cx.inn[x].iter().copied().collect();
            for a in ins {
                cx.out[a as usize].remove(&(x as u32));
            }
            cx.inn[x].clear();
            cx.alive[x] = false;
        }
        modified.retain(|&(a, b)| cx.alive[a] && cx.alive[b] && cx.out[a].contains_key(&(b as u32)));
        modified.sort();
        modified.dedup();
        modified
    }
}

/// Inverse over A of a square homogeneous matrix F (columns: source coordinates of
/// degrees `ds`, rows: target coordinates of degrees `dt`) whose degree-0 part F_0 is
/// invertible: F^{-1} = Σ_m (-F_0^{-1} F_+)^m F_0^{-1}, a finite sum since F_+ raises the
/// degree. Returns (F^{-1}, F == F_0), or None if the degrees do not match or F_0 is singular.
fn inverse_over_a<F: Field>(f: &PMat<F>, ds: &[i32], dt: &[i32], caps: &Caps) -> Option<(PMat<F>, bool)> {
    let k = ds.len();
    if f.nc != k || f.nr != k || dt.len() != k {
        return None;
    }
    let mut f0 = vec![vec![F::zero(); k]; k];
    let mut nconst = 0;
    let mut pure = true;
    for (r, row) in f.rows.iter().enumerate() {
        for (c, p) in row {
            if dt[r] == ds[*c as usize] {
                debug_assert!(p.len() == 1 && p[0].0 == 0, "non-homogeneous map component");
                f0[r][*c as usize] = const_term(p);
                nconst += 1;
            } else {
                pure = false;
            }
        }
    }
    if nconst < k {
        return None;
    }
    let g = dense_inverse(&f0)?;
    // G = F_0^{-1} (rows: source coordinates, cols: target coordinates)
    let gm = PMat { nr: k, nc: k, rows: g.iter().map(|r| r.iter().enumerate().filter(|(_, x)| !x.is_zero()).map(|(c, x)| (c as u32, vec![(0u64, x.clone())])).collect()).collect() };
    if pure {
        return Some((gm, true));
    }
    let fplus = PMat { nr: k, nc: k, rows: f.rows.iter().enumerate().map(|(r, row)| row.iter().filter(|(c, _)| dt[r] != ds[*c as usize]).cloned().collect()).collect() };
    let nmat = mat_axpy(&PMat { nr: k, nc: k, rows: vec![vec![]; k] }, &mat_mul(&gm, &fplus, caps), &F::one().neg());
    let mut sum = gm.clone();
    let mut term = gm;
    // every factor raises the degree by ≥ 2, so at most (degree spread)/2 + 1 terms
    let spread = (ds.iter().max().unwrap() - ds.iter().min().unwrap()).max(0) as usize;
    for _ in 0..=spread / 2 + 1 {
        term = mat_mul(&nmat, &term, caps);
        if term.is_zero() {
            return Some((sum, false));
        }
        sum = mat_axpy(&sum, &term, &F::one());
    }
    panic!("Neumann series did not terminate (non-homogeneous map?)");
}

fn normalize_rows<F: Field>(m: &mut PMat<F>) {
    for row in m.rows.iter_mut() {
        if let Some((_, p)) = row.first() {
            let c = p.last().unwrap().1.inv();
            if !c.is_one() {
                for (_, p) in row.iter_mut() {
                    *p = pscale(p, &c);
                }
            }
        }
    }
}

/// map[s] += c·q·src[s] for sparse vectors of truncated polynomials.
fn map_axpy<F: Field>(map: &mut BTreeMap<u32, P<F>>, src: &BTreeMap<u32, P<F>>, q: &P<F>, c: &F, caps: &Caps) {
    for (s, p) in src {
        let prod = pmul_trunc(q, p, caps);
        if prod.is_empty() {
            continue;
        }
        let cur = map.remove(s).unwrap_or_default();
        let new = axpy(&cur, &prod, c, 0);
        if !new.is_empty() {
            map.insert(*s, new);
        }
    }
}

/// Components V_y^{-1} T(F) V_x of a restricted map, per (child of x, child of y).
fn restrict_map<F: Field>(f: &PMat<F>, dx: &Dec<F>, dy: &Dec<F>, v: usize, w: usize, caps: &Caps) -> Vec<((u32, u32), PMat<F>)> {
    let cap = w - 1;
    let mut fcol: Vec<Vec<(u32, Vec<P<F>>)>> = vec![vec![]; f.nc];
    for (c, row) in f.rows.iter().enumerate() {
        for (a, p) in row {
            fcol[*a as usize].push((c as u32, split_var(p, v, cap)));
        }
    }
    let mut acc: BTreeMap<(u32, u32), BTreeMap<(u32, u32), Vec<(u64, F)>>> = BTreeMap::new();
    for (t, coord) in dx.coord.iter().enumerate() {
        let Some((alpha, pos)) = *coord else { continue };
        // column t of T(F)·V_x, indexed by restricted coordinates of y
        let mut wcol: BTreeMap<u32, Vec<(u64, F)>> = BTreeMap::new();
        for (sx, vp) in &dx.vcol[t] {
            let (a2, j2) = (*sx as usize / w, *sx as usize % w);
            for (c, split) in &fcol[a2] {
                for e in 0..(w - j2) {
                    let coef = &split[e];
                    if coef.is_empty() {
                        continue;
                    }
                    let sy = *c * w as u32 + (j2 + e) as u32;
                    pmul_trunc_into(vp, coef, caps, wcol.entry(sy).or_default());
                }
            }
        }
        for (sy, raw) in wcol {
            let wp = norm(raw);
            if wp.is_empty() {
                continue;
            }
            for (t2, vi) in &dy.vinv_by_s[sy as usize] {
                let (beta, pos2) = dy.coord[*t2 as usize].unwrap();
                pmul_trunc_into(vi, &wp, caps, acc.entry((alpha, beta)).or_default().entry((pos2, pos)).or_default());
            }
        }
    }
    let mut res = vec![];
    for ((alpha, beta), ents) in acc {
        let m = PMat::from_acc(dy.children[beta as usize].1.len(), dx.children[alpha as usize].1.len(), ents);
        if !m.is_zero() {
            res.push(((alpha, beta), m));
        }
    }
    res
}

fn verbose() -> bool {
    std::env::var("KRHOM_VERBOSE").map(|v| v != "0" && !v.is_empty()).unwrap_or(false)
}

/// Result of a tower closure.
pub struct TowerResult {
    pub poincare: Poincare,
    pub stages: Vec<StageStats>,
}

/// dim Φ_N(B_x) = dim Im(G_x ⊗ A_n) for every type, by the tower decomposition.
pub fn type_dims<F: Field>(caps: &[u32], bases: &BTreeMap<usize, TypeBasis<F>>, opts: TowerOpts) -> BTreeMap<usize, usize> {
    let mut tw = Tower::<F>::new(caps, opts);
    let mut out = BTreeMap::new();
    for (&x, tb) in bases {
        let d = tw.gram_matrix(&tb.gram, tb.ins.len());
        out.insert(x, tw.gdim(&d, 0));
    }
    out
}

/// gl_N closure of the final complex by the restriction tower. `caps[i]` is the
/// truncation cap of x_i (N-1, or 0 for the marked strand of the reduced theory), and
/// `global_q` the global q-shift, both exactly as for `close_exact`.
#[allow(clippy::too_many_arguments)]
pub fn close_tower<F: Field>(
    n: usize,
    caps: &[u32],
    global_q: i32,
    objs: &[Obj],
    edges: &[(usize, usize)],
    fpolys: &[Vec<Vec<EPoly<F>>>],
    bases: &BTreeMap<usize, TypeBasis<F>>,
    opts: TowerOpts,
) -> TowerResult {
    assert_eq!(caps.len(), n);
    let mut tw = Tower::<F>::new(caps, opts);
    let mut tmat: BTreeMap<usize, Rc<PMat<F>>> = BTreeMap::new();
    for (&x, tb) in bases {
        let m = tw.gram_matrix(&tb.gram, tb.ins.len());
        tmat.insert(x, m);
    }
    // stage n: Φ_N(B_x) = Im(G_x ⊗ A_n); objects with G_x ≡ 0 mod I vanish
    let mut obj = vec![];
    let mut idx: Vec<Option<u32>> = vec![None; objs.len()];
    for (i, o) in objs.iter().enumerate() {
        let d = tmat[&o.x].clone();
        if d.is_zero() {
            continue;
        }
        let deg: Vec<i32> = bases[&o.x].ins.iter().map(|&(_, ad)| o.q + ad).collect();
        let mut sdeg = deg.clone();
        sdeg.sort();
        idx[i] = Some(obj.len() as u32);
        obj.push(TObj { h: o.h, deg, sdeg, d });
    }
    let no = obj.len();
    let mut out: Vec<BTreeMap<u32, Edge<F>>> = (0..no).map(|_| BTreeMap::new()).collect();
    let mut inn: Vec<BTreeSet<u32>> = vec![BTreeSet::new(); no];
    // Parallel edges (same source and target) add up, as in close_exact: collect the
    // truncated terms of all edges a -> b first, then combine them.
    let mut acc: BTreeMap<(u32, u32), Vec<Vec<Vec<(u64, F)>>>> = BTreeMap::new();
    for (ei, &(s, t)) in edges.iter().enumerate() {
        let (Some(a), Some(b)) = (idx[s], idx[t]) else { continue };
        let fm = &fpolys[ei];
        assert_eq!(fm.len(), obj[b as usize].deg.len());
        let ent = acc.entry((a, b)).or_insert_with(|| vec![vec![]; fm.len()]);
        for (r, row) in fm.iter().enumerate() {
            if ent[r].len() < row.len() {
                ent[r].resize(row.len(), vec![]);
            }
            for (c, p) in row.iter().enumerate() {
                ent[r][c].extend(p.iter().filter(|(e, _)| e.iter().zip(caps.iter()).all(|(x, c)| x <= c)).map(|(e, c)| (pack(e), c.clone())));
            }
        }
    }
    for ((a, b), ent) in acc {
        let mut rows = vec![];
        for row in ent {
            let mut r = vec![];
            for (c, terms) in row.into_iter().enumerate() {
                let p = norm(terms);
                if !p.is_empty() {
                    r.push((c as u32, p));
                }
            }
            rows.push(r);
        }
        let m = PMat { nr: obj[b as usize].deg.len(), nc: obj[a as usize].deg.len(), rows };
        if !m.is_zero() {
            out[a as usize].insert(b, Edge { m, ver: 0 });
            inn[b as usize].insert(a);
        }
    }
    let mut cx = Cx { obj, alive: vec![true; no], out, inn };
    let mut stages = vec![];
    if opts.eliminate_top {
        let t0 = std::time::Instant::now();
        let mut st = StageStats { var: usize::MAX, objects: cx.n_alive(), edges: cx.n_edges(), terms: cx.n_terms(), ..Default::default() };
        tw.counts = [0; 4];
        tw.gauss(&mut cx, 0);
        [st.elims, st.elims_noncyclic, st.elims_neumann, st.elims_gdim] = tw.counts;
        st.after_objects = cx.n_alive();
        st.after_coords = cx.n_coords();
        st.after_edges = cx.n_edges();
        st.after_terms = cx.n_terms();
        st.secs = t0.elapsed().as_secs_f64();
        if verbose() {
            eprintln!("[tower] stage n: {} objects -> {} ({} eliminations, {} Neumann, {} by gdim)", st.objects, st.after_objects, st.elims, st.elims_neumann, st.elims_gdim);
        }
        stages.push(st);
    }
    for s in 0..n {
        let t0 = std::time::Instant::now();
        let mut st = StageStats { var: tw.order[s], ..Default::default() };
        cx = tw.restrict(cx, s, &mut st);
        st.objects = cx.n_alive();
        st.edges = cx.n_edges();
        st.terms = cx.n_terms();
        tw.counts = [0; 4];
        tw.gauss(&mut cx, s + 1);
        [st.elims, st.elims_noncyclic, st.elims_neumann, st.elims_gdim] = tw.counts;
        st.after_objects = cx.n_alive();
        st.after_coords = cx.n_coords();
        st.after_edges = cx.n_edges();
        st.after_terms = cx.n_terms();
        st.secs = t0.elapsed().as_secs_f64();
        if verbose() {
            eprintln!(
                "[tower] close x{}: {} objects ({} non-split coords), {} maps, {} terms -> {} objects ({} coords), {} maps, {} terms; {} eliminations ({} non-cyclic, {} Neumann, {} by gdim), {:.3}s",
                st.var + 1, st.objects, st.nonsplit, st.edges, st.terms, st.after_objects, st.after_coords, st.after_edges, st.after_terms,
                st.elims, st.elims_noncyclic, st.elims_neumann, st.elims_gdim, st.secs
            );
        }
        stages.push(st);
    }
    // (P6) stage 0: complete elimination leaves no map
    assert_eq!(cx.n_edges(), 0, "tower: maps left at stage 0");
    let mut pc = Poincare::default();
    for i in 0..cx.obj.len() {
        if cx.alive[i] {
            let o = &cx.obj[i];
            assert!(o.deg.len() == 1 && const_rank(&o.d) == 1, "tower: stage-0 object is not one-dimensional");
            pc.add_term(o.h, o.deg[0] + global_q, 1);
        }
    }
    TowerResult { poincare: pc, stages }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra::field::F31;

    fn p(terms: &[(&[u32], i64)]) -> P<F31> {
        norm(terms.iter().map(|(e, c)| (pack(e), F31::from_i64(*c))).collect())
    }

    #[test]
    fn division_and_products() {
        // (x0 - x1)(x0 + 2 x1 + x2) / (x0 - x1)
        let a = p(&[(&[1, 0, 0], 1), (&[0, 1, 0], -1)]);
        let b = p(&[(&[1, 0, 0], 1), (&[0, 1, 0], 2), (&[0, 0, 1], 1)]);
        let ab = pmul(&a, &b);
        assert_eq!(pdiv_exact(&ab, &a, 3), Some(b.clone()));
        assert_eq!(pdiv_exact(&ab, &b, 3), Some(a.clone()));
        let c = p(&[(&[1, 0, 0], 1), (&[0, 0, 1], 1)]);
        assert_eq!(pdiv_exact(&ab, &c, 3), None);
        let caps = Caps::new(&[1, 1, 1]);
        let t = pmul_trunc(&a, &a, &caps);
        assert_eq!(t, p(&[(&[1, 1, 0], -2)]));
        let s = split_var(&ab, 1, 1);
        assert_eq!(s[1], p(&[(&[1, 0, 0], 1), (&[0, 0, 1], -1)]));
    }

    #[test]
    fn norm_and_fused_products() {
        use rand::{rngs::StdRng, Rng, SeedableRng};
        // norm: combine, cancel, sort, and no excess capacity
        let mut v: Vec<(u64, F31)> = vec![];
        for (e, c) in [(3u32, 1i64), (1, 2), (3, -1), (2, 5), (1, 1), (0, 0), (2, -5), (4, 7)] {
            v.push((pack(&[e, 4 - e]), F31::from_i64(c)));
        }
        v.reserve(1000);
        let r = norm(v);
        assert_eq!(r, vec![(pack(&[4, 0]), F31::from_i64(7)), (pack(&[1, 3]), F31::from_i64(3))]);
        assert!(r.capacity() < 16);
        assert!(norm::<F31>(vec![(5, F31::zero())]).is_empty());
        // sub_mul == psub(a, pmul(q, b)) for homogeneous polynomials of all lengths
        let mut rng = StdRng::seed_from_u64(9);
        let rand_poly = |rng: &mut StdRng, deg: usize, k: usize| -> P<F31> {
            let mut terms = vec![];
            for _ in 0..k {
                let mut e = [0u32; 4];
                for _ in 0..deg {
                    e[rng.gen_range(0..4)] += 1;
                }
                terms.push((pack(&e), F31::from_i64(rng.gen_range(-3..4))));
            }
            norm(terms)
        };
        for _ in 0..300 {
            let (dq, db) = (rng.gen_range(0..5), rng.gen_range(0..5));
            let (kq, kb, ka) = (rng.gen_range(0..9), rng.gen_range(0..9), rng.gen_range(0..9));
            let q = rand_poly(&mut rng, dq, kq);
            let b = rand_poly(&mut rng, db, kb);
            let a = rand_poly(&mut rng, dq + db, ka);
            assert_eq!(sub_mul(&a, &q, &b), psub(&a, &pmul(&q, &b)));
        }
    }

    #[test]
    fn dense_products_match_sorted() {
        use rand::{rngs::StdRng, Rng, SeedableRng};
        let naive = |a: &P<F31>, b: &P<F31>, caps: Option<&Caps>| -> P<F31> {
            let mut v = vec![];
            for (m1, c1) in a {
                for (m2, c2) in b {
                    if caps.is_none_or(|c| c.fits(m1 + m2)) {
                        v.push((m1 + m2, c1.mul(c2)));
                    }
                }
            }
            norm(v)
        };
        let mut rng = StdRng::seed_from_u64(17);
        let mut dense_used = 0;
        for trial in 0..400 {
            let nv = rng.gen_range(1..6usize);
            let caps = Caps::new(&(0..nv).map(|_| rng.gen_range(0..8u32)).collect::<Vec<_>>());
            let mut rand_poly = |deg: usize, k: usize, homogeneous: bool, truncated: bool| -> P<F31> {
                let mut terms = vec![];
                for _ in 0..k {
                    let d = if homogeneous { deg } else { rng.gen_range(0..=deg) };
                    let mut e = vec![0u32; nv];
                    for _ in 0..d {
                        e[rng.gen_range(0..nv)] += 1;
                    }
                    terms.push((pack(&e), F31::from_i64(rng.gen_range(-3..4))));
                }
                let p = norm(terms);
                if truncated { trunc(&p, &caps) } else { p }
            };
            let hom = trial % 4 != 0;
            let truncated = trial % 2 == 0;
            let (da, db) = (1 + trial % 7, 1 + (trial / 7) % 6);
            let a = rand_poly(da, 10 + trial % 40, hom, truncated);
            let b = rand_poly(db, 10 + (trial * 7) % 40, hom, truncated);
            if let Some(d) = pmul_dense(&a, &b, None) {
                assert_eq!(d, naive(&a, &b, None), "untruncated trial {}", trial);
                dense_used += 1;
            }
            if truncated {
                if let Some(d) = pmul_dense(&a, &b, Some(&caps)) {
                    assert_eq!(d, naive(&a, &b, Some(&caps)), "truncated trial {}", trial);
                }
                assert_eq!(pmul_trunc(&a, &b, &caps), naive(&a, &b, Some(&caps)));
            }
            assert_eq!(pmul(&a, &b), naive(&a, &b, None));
        }
        assert!(dense_used > 100, "dense path used {} times", dense_used);
        // norm (dense for long lists) == norm_sort, homogeneous or not
        let mut dense_norms = 0;
        for trial in 0..300 {
            let nv = rng.gen_range(1..7usize);
            let len = rng.gen_range(1..400usize);
            let deg = rng.gen_range(0..12u32);
            let hom = trial % 3 != 0;
            let v: Vec<(u64, F31)> = (0..len)
                .map(|_| {
                    let d = if hom { deg } else { rng.gen_range(0..=deg) };
                    let mut e = vec![0u32; nv];
                    for _ in 0..d {
                        e[rng.gen_range(0..nv)] += 1;
                    }
                    (pack(&e), F31::from_i64(rng.gen_range(-2..3)))
                })
                .collect();
            if v.len() >= 64 && norm_dense(&v).is_some() {
                dense_norms += 1;
            }
            assert_eq!(norm(v.clone()), norm_sort(v), "trial {}", trial);
        }
        assert!(dense_norms > 50, "dense norm used {} times", dense_norms);
    }

    #[test]
    fn end_orders() {
        assert_eq!(EndPolicy::Right.order(4), vec![3, 2, 1, 0]);
        assert_eq!(EndPolicy::Left.order(4), vec![0, 1, 2, 3]);
        assert_eq!(EndPolicy::Alternate.order(5), vec![4, 0, 3, 1, 2]);
        assert_eq!(EndPolicy::Middle.order(4), vec![2, 1, 3, 0]);
        assert_eq!(EndPolicy::Middle.order(5), vec![2, 1, 3, 0, 4]);
        assert_eq!(EndPolicy::Middle.order(2), vec![1, 0]);
    }

    #[test]
    fn caps_fast_test() {
        let caps = Caps::new(&[0, 3, 9]);
        for e0 in 0..3u32 {
            for e1 in 0..8u32 {
                for e2 in 0..19u32 {
                    let m = pack(&[e0, e1, e2]);
                    assert_eq!(caps.fits_small(m), caps.fits(m));
                    assert_eq!(caps.fits(m), e0 == 0 && e1 <= 3 && e2 <= 9);
                }
            }
        }
    }

    #[test]
    fn inverse_over_a_random() {
        use rand::{rngs::StdRng, Rng, SeedableRng};
        let mut rng = StdRng::seed_from_u64(3);
        let caps = Caps::new(&[3, 3, 3]);
        let mut checked = 0;
        for _ in 0..200 {
            // coordinates of degrees 0, 2, 2, 4, 6 (shuffled), F = F_0 + F_+ homogeneous
            let mut ds = vec![0, 2, 2, 4, 6];
            for i in (1..ds.len()).rev() {
                ds.swap(i, rng.gen_range(0..=i));
            }
            let k = ds.len();
            let mut dt = ds.clone();
            for i in (1..k).rev() {
                dt.swap(i, rng.gen_range(0..=i));
            }
            let rand_poly = |rng: &mut StdRng, w: usize| -> P<F31> {
                let mut terms = vec![];
                for _ in 0..3 {
                    let mut e = [0u32; 3];
                    for _ in 0..w {
                        e[rng.gen_range(0..3)] += 1;
                    }
                    terms.push((pack(&e), F31::from_i64(rng.gen_range(-5..6))));
                }
                trunc(&norm(terms), &caps)
            };
            let mut rows = vec![];
            for r in 0..k {
                let mut row = vec![];
                for c in 0..k {
                    let d = ds[c] - dt[r];
                    if d < 0 {
                        continue;
                    }
                    let p = rand_poly(&mut rng, (d / 2) as usize);
                    if !p.is_empty() {
                        row.push((c as u32, p));
                    }
                }
                rows.push(row);
            }
            let f = PMat { nr: k, nc: k, rows };
            let Some((inv, pure)) = inverse_over_a(&f, &ds, &dt, &caps) else { continue };
            assert!(!pure);
            let id = PMat { nr: k, nc: k, rows: (0..k).map(|i| vec![(i as u32, vec![(0u64, F31::one())])]).collect() };
            assert_eq!(mat_mul(&f, &inv, &caps), id);
            assert_eq!(mat_mul(&inv, &f, &caps), id);
            checked += 1;
        }
        assert!(checked > 50, "only {} invertible samples", checked);
    }

    #[test]
    fn bs_dimension() {
        // Φ_N(B_s) on 2 strands: Im((x0 - x1) ⊗ A_2) has dimension N(N-1)
        for nn in 2..=6u32 {
            for policy in [EndPolicy::Right, EndPolicy::Left] {
                let mut tw = Tower::<F31>::new(&[nn - 1, nn - 1], TowerOpts::new(policy));
                let d = tw.intern(PMat { nr: 1, nc: 1, rows: vec![vec![(0, p(&[(&[1, 0], 1), (&[0, 1], -1)]))]] });
                assert_eq!(tw.gdim(&d, 0), (nn * (nn - 1)) as usize);
            }
        }
    }
}
