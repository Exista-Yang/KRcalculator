//! CL1: the colored closure layer (SYNTHESIS §9.1; R3 A1-A7; C2 M1, M2).
//!
//! Ring. R^λ = ⊗_i Sym(λ_i) = F[e_{i,j} : 1 ≤ j ≤ λ_i] with e_{i,j} of weight j (q-degree 2j).
//! Monomials are packed in a u128, 6 bits per variable, at most 20 variables; every
//! product is checked for carries out of a field ([`emono_mul`]).
//!
//! Contract v2 ([`ClosureDataC`]). Per final object (type, h, q) with q N-independent; per
//! type the degrees of R^λ-bases g_a of Hom(R^λ, B) ("ins") and h_b of Hom(B, R^λ) ("outs")
//! and the Gram matrix G[b][a] = h_b ∘ g_a; per final differential component the matrix
//! F[c][a] with f ∘ g_a = Σ_c F_ca g'_c; per crossing (positive, m = min label).
//!
//! Closure (R3 A3, PROVEN there): Φ_N(B) = Im(G ⊗ A_λ) ⊂ A_λ^{outs}, A_λ = ⊗_i H*(Gr(λ_i, N)),
//! maps act by v ↦ F v on the ins-coordinates.
//! * [`close_exact_colored`]: the exact oracle. A_λ in the `symfun::GrTable(λ_i, N)`
//!   e-monomial box basis at c = 0 (structure constants by the Pieri rule in Schur
//!   coordinates, equal to GrTable at c = 0); τ_N = coefficient of the top class (the functional of
//!   `sbim::Tau`); ranks per (h, q) by `RankEchelon`. As in `engine2::close_exact`, a state
//!   space is a set of representatives (a, m) (m a box monomial) whose images m·G_{·a} are a
//!   basis of Im(G ⊗ A_λ) in that degree, together with a set of image coordinates (b, box
//!   element) on which Im(G ⊗ A_λ) projects isomorphically ("duals"). The rank of a
//!   differential block is the rank of the projections of m·(G' F)_{·a} onto the target duals.
//! * Reduced mode (R3 A7): the marked block m gets the one-dimensional algebra (substitute
//!   e_{m,*} = 0 in G and F) and its term −λ_m(N − λ_m) is dropped from the global shift.
//! * [`global_q`] (SYNTHESIS §2.3): −Σ_c ε_c m_c (N − m_c) − Σ_{unmarked blocks} λ_i (N − λ_i).
//! * [`close_delta_tower`]: the (Δ_λ G, F) adapter (C2 M1(c)). G and F are expanded in root
//!   variables (e_{i,j} ↦ e_j(roots of block i)), G is multiplied by Δ_λ (the product of the
//!   per-block Vandermondes), and the validated `tower::close_tower` runs with caps N−1 per
//!   root. Since Δ_a R_N ≅ H*(Gr(a, N)) ⊗ V_a with F acting as F ⊗ id (C2 M1(b), given
//!   freeness), the result is the colored Poincaré polynomial times Π_i gr(V_{λ_i}); it is
//!   divided by [`gr_v`] (normalization fixed on single types against the exact oracle,
//!   tests/cclose.rs: gr(V_a) = Π_{k ≤ a} [k]_{q^2} unbalanced, no extra shift). Freeness is
//!   certified per (a, N) by gate G4 (`src/bin/g4_free.rs`, whose certificate is in fact a
//!   proof for every (a, N): Δ_a R_N ≅ F[x] ⊗_{Sym_a} H*(Gr(a, N)), see there); the adapter
//!   refuses pairs not in [`G4_CERTIFIED`].
//!
//! Exact closure data of small closed ladders: [`closure_data_from_ladder`] chooses R^λ-bases
//! of Hom(R^λ, B) and Hom(B, R^λ) from `sbim::hom_space` by exact graded Nakayama, with a
//! freeness certificate (Hilbert function) and the one-point determinant criterion of C2 M2
//! (identity-coset values at generic points) as an independent check. [`solve_f`] computes F
//! for a bimodule map exactly. [`from_engine2`] converts engine2's uncolored `ClosureData`.

use crate::engine::Obj;
use crate::engine2::{ClosureData, TypeBasis};
use crate::field::Field;
use crate::laurent::{LPoly, Poincare};
use crate::linsys::{Echelon, RankEchelon};
use crate::mpoly::{MonoCache, Poly};
use crate::pmat::PMat;
use crate::sbim::{hom_space, SBim};
use crate::sparse::{Acc, SpVec};
use crate::tower::TowerOpts;
use crate::web::Ladder;
use std::collections::{BTreeMap, BTreeSet, HashMap};

/// Small multiplicative hasher (FxHash) for integer keys in the hot loops.
#[derive(Default, Clone, Copy)]
struct FxHasher(u64);

impl std::hash::Hasher for FxHasher {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.add(b as u64);
        }
    }
    #[inline]
    fn write_u8(&mut self, i: u8) {
        self.add(i as u64)
    }
    #[inline]
    fn write_u32(&mut self, i: u32) {
        self.add(i as u64)
    }
    #[inline]
    fn write_u64(&mut self, i: u64) {
        self.add(i)
    }
    #[inline]
    fn write_u128(&mut self, i: u128) {
        self.add(i as u64);
        self.add((i >> 64) as u64)
    }
    #[inline]
    fn write_usize(&mut self, i: usize) {
        self.add(i as u64)
    }
    #[inline]
    fn finish(&self) -> u64 {
        self.0
    }
}

impl FxHasher {
    #[inline]
    fn add(&mut self, x: u64) {
        self.0 = (self.0.rotate_left(5) ^ x).wrapping_mul(0x51_7c_c1_b7_27_22_0a_95);
    }
}

type FxMap<K, V> = HashMap<K, V, std::hash::BuildHasherDefault<FxHasher>>;

/// Composition (labels per bottom block), zeros allowed (re-exported from `schur`).
pub use crate::schur::Comp;

// ------------------------------------------------------------------ the ring R^λ

/// Bits per variable of a packed R^λ monomial.
pub const EBITS: usize = 6;
/// At most 20 e-variables (6·20 = 120 bits).
pub const EMAX_VARS: usize = 20;
const EFIELD: u128 = 63;
/// Bit 6k for k = 1..=20: a carry out of field k−1 shows up there.
const ECARRY: u128 = {
    let mut m = 0u128;
    let mut k = 1;
    while k <= EMAX_VARS {
        m |= 1u128 << (EBITS * k);
        k += 1;
    }
    m
};

/// Exponent of variable v in a packed monomial.
#[inline]
pub fn emono_exp(m: u128, v: usize) -> u32 {
    ((m >> (EBITS * v)) & EFIELD) as u32
}

/// Product of packed monomials; panics if some exponent exceeds 63.
#[inline]
pub fn emono_mul(a: u128, b: u128) -> u128 {
    let s = a + b;
    assert!((a ^ b ^ s) & ECARRY == 0, "R^lambda monomial exponent overflow (> 63)");
    s
}

/// The ring R^λ: variable layout and weights.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ERing {
    pub lam: Comp,
    /// variable index of e_{i,1}
    pub voff: Vec<usize>,
    pub nvars: usize,
    /// weight j of each variable (q-degree 2j)
    pub wt: Vec<u32>,
    /// block of each variable
    pub blk: Vec<usize>,
}

impl ERing {
    pub fn new(lam: &[u8]) -> ERing {
        let mut voff = vec![];
        let mut wt = vec![];
        let mut blk = vec![];
        for (i, &l) in lam.iter().enumerate() {
            voff.push(wt.len());
            for j in 1..=l as u32 {
                wt.push(j);
                blk.push(i);
            }
        }
        assert!(wt.len() <= EMAX_VARS, "R^lambda: {} e-variables > {}", wt.len(), EMAX_VARS);
        ERing { lam: lam.to_vec(), voff, nvars: wt.len(), wt, blk }
    }
    /// Variable index of e_{i,j} (j is 1-based).
    pub fn var(&self, i: usize, j: usize) -> usize {
        assert!(j >= 1 && j <= self.lam[i] as usize);
        self.voff[i] + j - 1
    }
    pub fn pack(&self, exps: &[u32]) -> u128 {
        assert_eq!(exps.len(), self.nvars);
        exps.iter().enumerate().fold(0u128, |acc, (v, &e)| {
            assert!(e as u128 <= EFIELD, "R^lambda exponent {} > 63", e);
            acc | ((e as u128) << (EBITS * v))
        })
    }
    pub fn unpack(&self, m: u128) -> Vec<u32> {
        (0..self.nvars).map(|v| emono_exp(m, v)).collect()
    }
    /// Weight Σ j·exp (q-degree = 2·weight).
    pub fn weight(&self, m: u128) -> u32 {
        (0..self.nvars).map(|v| self.wt[v] * emono_exp(m, v)).sum()
    }
    /// Exponents of block i's variables.
    pub fn block_exps(&self, m: u128, i: usize) -> Vec<u32> {
        (0..self.lam[i] as usize).map(|j| emono_exp(m, self.voff[i] + j)).collect()
    }
    /// All monomials of weight exactly w.
    pub fn monos_of_weight(&self, w: u32) -> Vec<u128> {
        fn rec(wt: &[u32], v: usize, left: u32, cur: u128, out: &mut Vec<u128>) {
            if v == wt.len() {
                if left == 0 {
                    out.push(cur);
                }
                return;
            }
            let mut k = 0u32;
            while k * wt[v] <= left {
                assert!(k as u128 <= EFIELD);
                rec(wt, v + 1, left - k * wt[v], cur | ((k as u128) << (EBITS * v)), out);
                k += 1;
            }
        }
        let mut out = vec![];
        rec(&self.wt, 0, w, 0, &mut out);
        out.sort_unstable();
        out
    }
}

/// Element of R^λ: terms (packed monomial, coefficient) sorted by monomial, no zeros.
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct EPoly<F> {
    pub t: Vec<(u128, F)>,
}

impl<F: Field> EPoly<F> {
    pub fn zero() -> Self {
        EPoly { t: vec![] }
    }
    pub fn constant(c: F) -> Self {
        Self::mono(0, c)
    }
    pub fn one() -> Self {
        Self::constant(F::one())
    }
    pub fn mono(m: u128, c: F) -> Self {
        if c.is_zero() { Self::zero() } else { EPoly { t: vec![(m, c)] } }
    }
    pub fn is_zero(&self) -> bool {
        self.t.is_empty()
    }
    /// Normalize an arbitrary term list (sort, combine, drop zeros).
    pub fn from_terms(mut v: Vec<(u128, F)>) -> Self {
        v.sort_unstable_by_key(|x| x.0);
        let mut out: Vec<(u128, F)> = Vec::with_capacity(v.len());
        for (m, c) in v {
            if let Some(last) = out.last_mut() {
                if last.0 == m {
                    last.1.add_assign(&c);
                    continue;
                }
            }
            out.push((m, c));
        }
        out.retain(|x| !x.1.is_zero());
        EPoly { t: out }
    }
    pub fn add(&self, o: &Self) -> Self {
        let mut out = Vec::with_capacity(self.t.len() + o.t.len());
        let (mut i, mut j) = (0, 0);
        while i < self.t.len() || j < o.t.len() {
            if j == o.t.len() || (i < self.t.len() && self.t[i].0 < o.t[j].0) {
                out.push(self.t[i].clone());
                i += 1;
            } else if i == self.t.len() || o.t[j].0 < self.t[i].0 {
                out.push(o.t[j].clone());
                j += 1;
            } else {
                let c = self.t[i].1.add(&o.t[j].1);
                if !c.is_zero() {
                    out.push((self.t[i].0, c));
                }
                i += 1;
                j += 1;
            }
        }
        EPoly { t: out }
    }
    pub fn add_assign(&mut self, o: &Self) {
        if !o.is_zero() {
            *self = self.add(o);
        }
    }
    pub fn scale(&self, c: &F) -> Self {
        if c.is_zero() {
            return Self::zero();
        }
        EPoly { t: self.t.iter().map(|(m, x)| (*m, x.mul(c))).collect() }
    }
    pub fn neg(&self) -> Self {
        self.scale(&F::one().neg())
    }
    pub fn sub(&self, o: &Self) -> Self {
        self.add(&o.neg())
    }
    pub fn mul_mono(&self, m: u128, c: &F) -> Self {
        if c.is_zero() {
            return Self::zero();
        }
        EPoly { t: self.t.iter().map(|(a, x)| (emono_mul(*a, m), x.mul(c))).collect() }
    }
    pub fn mul(&self, o: &Self) -> Self {
        if self.is_zero() || o.is_zero() {
            return Self::zero();
        }
        if o.t.len() == 1 {
            return self.mul_mono(o.t[0].0, &o.t[0].1);
        }
        if self.t.len() == 1 {
            return o.mul_mono(self.t[0].0, &self.t[0].1);
        }
        let mut h: HashMap<u128, F> = HashMap::with_capacity(self.t.len() * o.t.len());
        for (a, x) in &self.t {
            for (b, y) in &o.t {
                h.entry(emono_mul(*a, *b)).or_insert_with(F::zero).fma(x, y);
            }
        }
        Self::from_terms(h.into_iter().collect())
    }
    pub fn coeff(&self, m: u128) -> F {
        match self.t.binary_search_by_key(&m, |x| x.0) {
            Ok(i) => self.t[i].1.clone(),
            Err(_) => F::zero(),
        }
    }
    /// Some(w) if every term has weight w (None for 0 or inhomogeneous).
    pub fn weight_homog(&self, ring: &ERing) -> Option<u32> {
        let w = ring.weight(self.t.first()?.0);
        self.t.iter().all(|(m, _)| ring.weight(*m) == w).then_some(w)
    }
    /// Substitute e_{i,*} = 0 (block i).
    pub fn drop_block(&self, ring: &ERing, i: usize) -> Self {
        let mask: u128 = (0..ring.lam[i] as usize).fold(0u128, |acc, j| acc | (EFIELD << (EBITS * (ring.voff[i] + j))));
        EPoly { t: self.t.iter().filter(|(m, _)| m & mask == 0).cloned().collect() }
    }
    /// From an `mpoly::Poly` in the same variable order (8 bits per variable, e.g. the
    /// bottom ring of `sbim::SBim`).
    pub fn from_mpoly(p: &Poly<F>, nvars: usize) -> Self {
        let t = p
            .t
            .iter()
            .map(|(m, c)| {
                let mut e = 0u128;
                for v in 0..nvars {
                    let x = crate::mpoly::mono_exp(*m, v) as u128;
                    assert!(x <= EFIELD, "exponent {} > 63", x);
                    e |= x << (EBITS * v);
                }
                assert!(nvars >= 16 || *m >> (8 * nvars) == 0, "mpoly monomial uses more than {} variables", nvars);
                (e, c.clone())
            })
            .collect();
        Self::from_terms(t)
    }
    /// To an `mpoly::Poly` (at most 16 variables).
    pub fn to_mpoly(&self, nvars: usize) -> Poly<F> {
        assert!(nvars <= crate::mpoly::MAXVARS);
        let t: Vec<(u128, F)> = self
            .t
            .iter()
            .map(|(m, c)| {
                assert!(*m >> (EBITS * nvars) == 0);
                let mut e = 0u128;
                for v in 0..nvars {
                    e |= (emono_exp(*m, v) as u128) << (8 * v);
                }
                (e, c.clone())
            })
            .collect();
        let mut t = t;
        t.sort_unstable_by_key(|x| x.0);
        Poly { t }
    }
    /// From an exponent-vector polynomial (engine2's `EPoly`), variable v = exps[v].
    pub fn from_exps(p: &HashMap<Vec<u32>, F>) -> Self {
        Self::from_terms(
            p.iter()
                .filter(|(_, c)| !c.is_zero())
                .map(|(e, c)| {
                    assert!(e.len() <= EMAX_VARS);
                    let m = e.iter().enumerate().fold(0u128, |acc, (v, &x)| {
                        assert!(x as u128 <= EFIELD, "exponent {} > 63", x);
                        acc | ((x as u128) << (EBITS * v))
                    });
                    (m, c.clone())
                })
                .collect(),
        )
    }
}

// ------------------------------------------------------------------ contract v2

/// R^λ-bases of one type: degrees of g_a (ins) and h_b (outs) and G[b][a] = h_b ∘ g_a.
#[derive(Clone, Debug)]
pub struct TypeBasisC<F> {
    pub ins: Vec<i32>,
    pub outs: Vec<i32>,
    pub gram: Vec<Vec<EPoly<F>>>,
} // degrees, G[b][a]

/// N-independent closure data of a colored minimal complex (contract v2).
#[derive(Clone, Debug)]
pub struct ClosureDataC<F> {
    pub lam: Comp, // closable: final labels = lam
    pub objs: Vec<(u32 /*type*/, i32 /*h*/, i32 /*q, N-independent*/)>,
    pub types: BTreeMap<u32, TypeBasisC<F>>,
    pub edges: Vec<(usize, usize, Vec<Vec<EPoly<F>>>)>, // F[c][a]: rows = target ins, cols = source ins
    pub crossings: Vec<(bool /*positive*/, u8 /*m = min label*/)>,
}

/// Global q-shift at gl_N (SYNTHESIS §2.3, R3 A4):
/// −Σ_c ε_c m_c (N − m_c) − Σ_{unmarked blocks} λ_i (N − λ_i).
/// For labels 1 this is `engine2::closure_q_shift`.
pub fn global_q(lam: &[u8], crossings: &[(bool, u8)], nn: usize, marked: Option<usize>) -> i32 {
    let n = nn as i32;
    let cr: i32 = crossings
        .iter()
        .map(|&(pos, m)| {
            let m = m as i32;
            let x = m * (n - m);
            if pos { -x } else { x }
        })
        .sum();
    let bl: i32 = lam.iter().enumerate().filter(|&(i, _)| Some(i) != marked).map(|(_, &l)| (l as i32) * (n - l as i32)).sum();
    cr - bl
}

impl<F: Field> ClosureDataC<F> {
    pub fn ring(&self) -> ERing {
        ERing::new(&self.lam)
    }
    pub fn global_q(&self, nn: usize, marked: Option<usize>) -> i32 {
        global_q(&self.lam, &self.crossings, nn, marked)
    }
    /// Shapes, homological degrees and homogeneity of every entry.
    pub fn validate(&self) -> Result<(), String> {
        let ring = self.ring();
        let nv = ring.nvars;
        let fits = |p: &EPoly<F>| p.t.iter().all(|(m, _)| nv == EMAX_VARS || m >> (EBITS * nv) == 0);
        for (x, tb) in &self.types {
            if tb.gram.len() != tb.outs.len() {
                return Err(format!("type {}: gram has {} rows, {} outs", x, tb.gram.len(), tb.outs.len()));
            }
            for (b, row) in tb.gram.iter().enumerate() {
                if row.len() != tb.ins.len() {
                    return Err(format!("type {}: gram row {} has {} entries, {} ins", x, b, row.len(), tb.ins.len()));
                }
                for (a, p) in row.iter().enumerate() {
                    if !fits(p) {
                        return Err(format!("type {}: G[{}][{}] uses variables outside R^lambda", x, b, a));
                    }
                    let d = tb.ins[a] + tb.outs[b];
                    if p.t.iter().any(|(m, _)| 2 * ring.weight(*m) as i32 != d) {
                        return Err(format!("type {}: G[{}][{}] not homogeneous of degree {}", x, b, a, d));
                    }
                }
            }
        }
        for (i, o) in self.objs.iter().enumerate() {
            if !self.types.contains_key(&o.0) {
                return Err(format!("object {}: unknown type {}", i, o.0));
            }
        }
        for (e, (s, t, f)) in self.edges.iter().enumerate() {
            let (os, ot) = (self.objs[*s], self.objs[*t]);
            if ot.1 != os.1 + 1 {
                return Err(format!("edge {}: h {} -> {}", e, os.1, ot.1));
            }
            let (ts, tt) = (&self.types[&os.0], &self.types[&ot.0]);
            if f.len() != tt.ins.len() || f.iter().any(|r| r.len() != ts.ins.len()) {
                return Err(format!("edge {}: F is not {} x {}", e, tt.ins.len(), ts.ins.len()));
            }
            for (c, row) in f.iter().enumerate() {
                for (a, p) in row.iter().enumerate() {
                    if !fits(p) {
                        return Err(format!("edge {}: F[{}][{}] uses variables outside R^lambda", e, c, a));
                    }
                    let d = os.2 - ot.2 + ts.ins[a] - tt.ins[c];
                    if p.t.iter().any(|(m, _)| 2 * ring.weight(*m) as i32 != d) {
                        return Err(format!("edge {}: F[{}][{}] not homogeneous of degree {}", e, c, a, d));
                    }
                }
            }
        }
        Ok(())
    }
}

/// Engine2's uncolored data (λ = 1^n) in contract v2: e_{i,1} = x_i, crossings of label 1.
pub fn from_engine2<F: Field>(d: &ClosureData<F>) -> ClosureDataC<F> {
    let conv = |m: &Vec<Vec<HashMap<Vec<u32>, F>>>| -> Vec<Vec<EPoly<F>>> { m.iter().map(|r| r.iter().map(EPoly::from_exps).collect()).collect() };
    let types = d
        .bases
        .iter()
        .map(|(&x, tb)| (x as u32, TypeBasisC { ins: tb.ins.iter().map(|g| g.1).collect(), outs: tb.outs.iter().map(|g| g.1).collect(), gram: conv(&tb.gram) }))
        .collect();
    ClosureDataC {
        lam: vec![1; d.n],
        objs: d.objs.iter().map(|o| (o.x as u32, o.h, o.q)).collect(),
        types,
        edges: d.edges.iter().zip(d.fpolys.iter()).map(|(&(s, t), f)| (s, t, conv(f))).collect(),
        crossings: d.word.iter().map(|g| (g.positive, 1u8)).collect(),
    }
}

// ------------------------------------------------------------------ A_λ = ⊗ H*(Gr(λ_i, N))

/// One block of A_λ: H*(Gr(a, N)) in the GrTable e-monomial box basis at c = 0, or the
/// one-dimensional algebra (marked block in reduced mode, label 0).
struct BlockAlg<F: Field> {
    a: usize,
    /// box basis exponents α (|α| ≤ N − a)
    basis: Vec<Vec<u32>>,
    wt: Vec<u32>,
    top_wt: u32,
    top_idx: u32,
    /// normal form of e^β for every β of weight ≤ top_wt (key: 8 bits per exponent)
    nf: HashMap<u128, Vec<(u32, F)>>,
    /// the same as a dense table (index Σ_j β_j·dstride[j], β_j ≤ top_wt/(j+1)) when small:
    /// terms dterms[doff[i]..doff[i+1]]
    dense: Option<(Vec<u32>, Vec<u32>, Vec<(u32, F)>)>,
}

impl<F: Field> BlockAlg<F> {
    /// Build the dense lookup table from `nf` (if it has at most 4e6 entries).
    fn with_dense(mut self) -> Self {
        let caps: Vec<u64> = (0..self.a).map(|j| (self.top_wt / (j as u32 + 1)) as u64 + 1).collect();
        let size = caps.iter().try_fold(1u64, |acc, &c| acc.checked_mul(c));
        let Some(size) = size.filter(|&s| s <= 4_000_000) else { return self };
        let mut dstride = vec![0u32; self.a];
        let mut s = 1u64;
        for j in 0..self.a {
            dstride[j] = s as u32;
            s *= caps[j];
        }
        let mut entries: Vec<(u64, &Vec<(u32, F)>)> = self
            .nf
            .iter()
            .map(|(k, v)| {
                let idx: u64 = (0..self.a).map(|j| ((k >> (8 * j)) & 0xff) as u64 * dstride[j] as u64).sum();
                (idx, v)
            })
            .collect();
        entries.sort_unstable_by_key(|e| e.0);
        let mut doff = vec![0u32; size as usize + 1];
        let mut dterms = vec![];
        let mut it = entries.into_iter().peekable();
        for i in 0..size as usize {
            doff[i] = dterms.len() as u32;
            if let Some((idx, v)) = it.peek() {
                if *idx == i as u64 {
                    dterms.extend(v.iter().cloned());
                    it.next();
                }
            }
        }
        doff[size as usize] = dterms.len() as u32;
        self.dense = Some((dstride, doff, dterms));
        self
    }

    /// Normal form of e^β (β of weight ≤ top_wt).
    #[inline]
    fn lookup(&self, ex: &[u32]) -> &[(u32, F)] {
        match &self.dense {
            Some((dstride, doff, dterms)) => {
                let i: usize = ex.iter().zip(dstride.iter()).map(|(&e, &s)| e as usize * s as usize).sum();
                &dterms[doff[i] as usize..doff[i + 1] as usize]
            }
            None => &self.nf[&bkey(ex)],
        }
    }
}

#[inline]
fn bkey(exps: &[u32]) -> u128 {
    exps.iter().enumerate().fold(0u128, |acc, (j, &e)| acc | ((e as u128) << (8 * j)))
}

impl<F: Field> BlockAlg<F> {
    fn trivial(a: usize) -> Self {
        let mut nf = HashMap::new();
        nf.insert(0u128, vec![(0u32, F::one())]);
        BlockAlg { a, basis: vec![vec![0; a]], wt: vec![0], top_wt: 0, top_idx: 0, nf, dense: None }.with_dense()
    }
    /// H*(Gr(a, N)) = Sym_a / span{s_ν : ν_1 > N − a} in the e-monomial box basis
    /// {e^α : |α| ≤ N − a} (the basis of `symfun::GrTable(a, N)`, whose c = 0 structure
    /// constants these are; tests/cclose.rs compares them). Products are computed in Schur
    /// coordinates by the Pieri rule for e_j (vertical strips) with ν_1 > N − a dropped (the
    /// dropped span is the kernel, an ideal), then converted to the e-basis by the inverse of
    /// the unitriangular matrix S[α] = Schur coordinates of e^α. (GrTable::compute builds
    /// the full equivariant table and is far too slow beyond N ≈ 7.)
    fn grassmannian(a: usize, nn: usize) -> Self {
        assert!(a >= 1 && a <= nn && a <= 16);
        let s = nn - a;
        let top_wt = (a * s) as u32;
        assert!(top_wt <= 255, "H*(Gr({},{})): top weight {} > 255", a, nn, top_wt);
        // e-basis: α ∈ N^a with |α| ≤ N − a (same enumeration order as GrTable)
        let mut basis: Vec<Vec<u32>> = vec![];
        for tot in 0..=s {
            fn rec(i: usize, left: usize, cur: &mut Vec<u32>, out: &mut Vec<Vec<u32>>) {
                if i == cur.len() {
                    if left == 0 {
                        out.push(cur.clone());
                    }
                    return;
                }
                for k in 0..=left {
                    cur[i] = k as u32;
                    rec(i + 1, left - k, cur, out);
                }
                cur[i] = 0;
            }
            rec(0, tot, &mut vec![0; a], &mut basis);
        }
        let weight = |al: &Vec<u32>| -> u32 { al.iter().enumerate().map(|(j, &k)| (j as u32 + 1) * k).sum() };
        basis.sort_by_key(|al| (weight(al), al.clone())); // GrTable's order
        let wt: Vec<u32> = basis.iter().map(weight).collect();
        let tops: Vec<usize> = (0..wt.len()).filter(|&i| wt[i] == top_wt).collect();
        assert_eq!(tops.len(), 1, "H*(Gr({},{})): top degree not one-dimensional", a, nn);
        let zero = basis.iter().position(|al| al.iter().all(|&x| x == 0)).unwrap();
        // Schur box basis: partitions ν ⊂ a × (N − a), as non-increasing Vec<u8> of length a
        let mut parts: Vec<Vec<u8>> = vec![];
        fn prec(i: usize, maxv: usize, cur: &mut Vec<u8>, out: &mut Vec<Vec<u8>>) {
            if i == cur.len() {
                out.push(cur.clone());
                return;
            }
            for v in 0..=maxv {
                cur[i] = v as u8;
                prec(i + 1, v, cur, out);
            }
            cur[i] = 0;
        }
        prec(0, s, &mut vec![0; a], &mut parts);
        assert_eq!(parts.len(), basis.len());
        let pidx: HashMap<Vec<u8>, usize> = parts.iter().enumerate().map(|(i, p)| (p.clone(), i)).collect();
        // Pieri: s_ν · e_j = Σ s_μ over vertical strips μ/ν of size j (ℓ(μ) ≤ a), μ_1 ≤ N − a
        let pieri = |v: &BTreeMap<usize, F>, j: usize| -> BTreeMap<usize, F> {
            let mut out: BTreeMap<usize, F> = BTreeMap::new();
            for (&pi, c) in v {
                let nu = &parts[pi];
                for mask in 0u32..(1u32 << a) {
                    if mask.count_ones() as usize != j {
                        continue;
                    }
                    let mu: Vec<u8> = (0..a).map(|r| nu[r] + ((mask >> r) & 1) as u8).collect();
                    if mu.windows(2).any(|w| w[0] < w[1]) || mu[0] as usize > s {
                        continue;
                    }
                    out.entry(pidx[&mu]).or_insert_with(F::zero).add_assign(c);
                }
            }
            out.retain(|_, c| !c.is_zero());
            out
        };
        // Schur coordinates of e^β for every β of weight ≤ top_wt (increasing weight)
        let mut by_w: Vec<Vec<Vec<u32>>> = vec![vec![]; top_wt as usize + 1];
        fn rec(a: usize, j: usize, left: u32, cur: &mut Vec<u32>, w: u32, out: &mut Vec<Vec<Vec<u32>>>) {
            if j == a {
                out[w as usize].push(cur.clone());
                return;
            }
            let wj = j as u32 + 1;
            let mut k = 0u32;
            while k * wj <= left {
                cur[j] = k;
                rec(a, j + 1, left - k * wj, cur, w + k * wj, out);
                k += 1;
            }
            cur[j] = 0;
        }
        rec(a, 0, top_wt, &mut vec![0; a], 0, &mut by_w);
        let mut schur: HashMap<u128, BTreeMap<usize, F>> = HashMap::new();
        schur.insert(0, BTreeMap::from([(pidx[&vec![0u8; a]], F::one())]));
        for w in 1..=top_wt as usize {
            for be in &by_w[w] {
                let j = be.iter().position(|&x| x > 0).unwrap();
                let mut prev = be.clone();
                prev[j] -= 1;
                let v = pieri(&schur[&bkey(&prev)], j + 1);
                schur.insert(bkey(be), v);
            }
        }
        // S[α][ν] (rows: e-basis) and its inverse: Schur coordinates -> e-coordinates
        let dim = basis.len();
        let smat: Vec<Vec<F>> = basis
            .iter()
            .map(|al| {
                let mut row = vec![F::zero(); dim];
                for (pi, c) in &schur[&bkey(al)] {
                    row[*pi] = c.clone();
                }
                row
            })
            .collect();
        let sinv = crate::local::inverse(&smat).expect("e-monomial box basis is not a basis of H*(Gr(a,N))");
        let mut nf: HashMap<u128, Vec<(u32, F)>> = HashMap::new();
        for (key, v) in &schur {
            // c = w · S^{-1}
            let mut c = vec![F::zero(); dim];
            for (pi, x) in v {
                for (k, y) in sinv[*pi].iter().enumerate() {
                    if !y.is_zero() {
                        c[k].fma(x, y);
                    }
                }
            }
            nf.insert(*key, c.into_iter().enumerate().filter(|(_, x)| !x.is_zero()).map(|(k, x)| (k as u32, x)).collect());
        }
        debug_assert_eq!(nf[&0], vec![(zero as u32, F::one())]);
        BlockAlg { a, basis, wt, top_wt, top_idx: tops[0] as u32, nf, dense: None }.with_dense()
    }

    /// Multiplication by e_j (1-based) on the e-basis: e^{α+ε_j} in normal form.
    fn e_mul(&self, j: usize, idx: usize) -> Vec<(u32, F)> {
        let mut be = self.basis[idx].clone();
        be[j - 1] += 1;
        let w: u32 = be.iter().enumerate().map(|(i, &k)| (i as u32 + 1) * k).sum();
        if w > self.top_wt {
            return vec![];
        }
        self.nf[&bkey(&be)].clone()
    }
}

/// The c = 0 multiplication table of H*(Gr(a, N)) used by the exact oracle:
/// result[j−1][α] = e_j · e^α in the e-monomial box basis (`basis()` order), for comparison
/// with `symfun::GrTable(a, N)` at c = 0.
pub fn grassmannian_table<F: Field>(a: usize, nn: usize) -> (Vec<Vec<u32>>, Vec<Vec<Vec<(u32, F)>>>) {
    let b = BlockAlg::<F>::grassmannian(a, nn);
    let tab = (1..=a).map(|j| (0..b.basis.len()).map(|i| b.e_mul(j, i)).collect()).collect();
    (b.basis.clone(), tab)
}


/// A_λ = ⊗_i H*(Gr(λ_i, N)) (the marked block replaced by F in reduced mode), with normal
/// forms of R^λ monomials in the box basis and the Frobenius form τ_N.
pub struct GrAlg<F: Field> {
    pub lam: Comp,
    pub nn: usize,
    pub marked: Option<usize>,
    pub ring: ERing,
    blocks: Vec<BlockAlg<F>>,
    /// mixed-radix strides of the box multi-index key
    strides: Vec<u64>,
    pub dim: u64,
    pub top_weight: u32,
    top_key: u64,
    /// box basis by weight: (key, packed R^λ monomial)
    by_weight: Vec<Vec<(u64, u128)>>,
}

impl<F: Field> GrAlg<F> {
    pub fn new(lam: &[u8], nn: usize, marked: Option<usize>) -> Self {
        assert!(lam.iter().all(|&l| l as usize <= nn), "GrAlg: a label exceeds N = {}", nn);
        assert!(lam.iter().all(|&l| l <= 16), "GrAlg: labels above 16 are not supported");
        if let Some(m) = marked {
            assert!(m < lam.len() && lam[m] > 0, "marked block {} invalid for {:?}", m, lam);
        }
        let ring = ERing::new(lam);
        let blocks: Vec<BlockAlg<F>> = lam
            .iter()
            .enumerate()
            .map(|(i, &l)| if l == 0 || Some(i) == marked { BlockAlg::trivial(l as usize) } else { BlockAlg::grassmannian(l as usize, nn) })
            .collect();
        let mut strides = vec![];
        let mut dim: u64 = 1;
        for b in &blocks {
            strides.push(dim);
            dim = dim.checked_mul(b.basis.len() as u64).expect("A_lambda dimension overflows u64");
        }
        let top_weight: u32 = blocks.iter().map(|b| b.top_wt).sum();
        let top_key: u64 = blocks.iter().zip(strides.iter()).map(|(b, s)| b.top_idx as u64 * s).sum();
        assert!(dim <= 50_000_000, "A_lambda of dimension {} is too large for the exact oracle", dim);
        let mut by_weight: Vec<Vec<(u64, u128)>> = vec![vec![]; top_weight as usize + 1];
        let mut idx = vec![0usize; blocks.len()];
        loop {
            let mut key = 0u64;
            let mut mono = 0u128;
            let mut w = 0u32;
            for (i, b) in blocks.iter().enumerate() {
                key += idx[i] as u64 * strides[i];
                w += b.wt[idx[i]];
                if Some(i) != marked {
                    for (j, &e) in b.basis[idx[i]].iter().enumerate() {
                        mono |= (e as u128) << (EBITS * (ring.voff[i] + j));
                    }
                }
            }
            by_weight[w as usize].push((key, mono));
            let mut i = 0;
            loop {
                if i == blocks.len() {
                    for v in by_weight.iter_mut() {
                        v.sort_unstable();
                    }
                    return GrAlg { lam: lam.to_vec(), nn, marked, ring, blocks, strides, dim, top_weight, top_key, by_weight };
                }
                idx[i] += 1;
                if idx[i] < blocks[i].basis.len() {
                    break;
                }
                idx[i] = 0;
                i += 1;
            }
        }
    }

    /// Box basis elements of weight w: (key, packed R^λ monomial).
    pub fn box_of_weight(&self, w: u32) -> &[(u64, u128)] {
        self.by_weight.get(w as usize).map(|v| v.as_slice()).unwrap_or(&[])
    }

    /// Sort key of a box element: lex-descending in the per-block indices, block 0 most
    /// significant (for λ = 1^n: lex-descending exponent vectors).
    pub fn order_key(&self, key: u64) -> u64 {
        let mut r = 0u64;
        for (i, b) in self.blocks.iter().enumerate() {
            let size = b.basis.len() as u64;
            let idx = (key / self.strides[i]) % size;
            r = r * size + (size - 1 - idx);
        }
        r
    }

    /// Normal form of the product of two packed R^λ monomials (exponents added without
    /// packing, so no 6-bit overflow), as (box key, coefficient).
    pub fn nf_prod(&self, m1: u128, m2: u128, out: &mut Vec<(u64, F)>) {
        out.clear();
        out.push((0, F::one()));
        let mut ex = [0u32; 16];
        for (i, b) in self.blocks.iter().enumerate() {
            let a = b.a;
            if a == 0 {
                continue;
            }
            let off = self.ring.voff[i];
            let mut w = 0u32;
            for j in 0..a {
                let e = emono_exp(m1, off + j) + emono_exp(m2, off + j);
                ex[j] = e;
                w += e * (j as u32 + 1);
            }
            if w > b.top_wt {
                out.clear();
                return;
            }
            let v = b.lookup(&ex[..a]);
            if v.is_empty() {
                out.clear();
                return;
            }
            let s = self.strides[i];
            if v.len() == 1 {
                let (k, c) = (&v[0].0, &v[0].1);
                let add = *k as u64 * s;
                let unit = c.is_one();
                for t in out.iter_mut() {
                    t.0 += add;
                    if !unit {
                        t.1 = t.1.mul(c);
                    }
                }
            } else {
                let prev = std::mem::take(out);
                for (k0, c0) in &prev {
                    for (k, c) in v {
                        out.push((k0 + *k as u64 * s, c0.mul(c)));
                    }
                }
            }
        }
    }

    /// Normal form of p ∈ R^λ in A_λ (sorted by key, combined).
    pub fn nf(&self, p: &EPoly<F>) -> Vec<(u64, F)> {
        let mut acc: BTreeMap<u64, F> = BTreeMap::new();
        let mut buf = vec![];
        for (m, c) in &p.t {
            self.nf_prod(*m, 0, &mut buf);
            for (k, x) in &buf {
                acc.entry(*k).or_insert_with(F::zero).fma(c, x);
            }
        }
        acc.into_iter().filter(|(_, c)| !c.is_zero()).collect()
    }

    /// τ_N(p): coefficient of the top class (product of the per-block top classes e_a^{N−a}).
    pub fn tau(&self, p: &EPoly<F>) -> F {
        let mut s = F::zero();
        let mut buf = vec![];
        for (m, c) in &p.t {
            if self.ring.weight(*m) != self.top_weight {
                continue;
            }
            self.nf_prod(*m, 0, &mut buf);
            for (k, x) in &buf {
                if *k == self.top_key {
                    s.fma(c, x);
                }
            }
        }
        s
    }
}

// ------------------------------------------------------------------ exact oracle

/// State space of one type: representatives (a, box monomial) per degree d (degree of the
/// representative before the object shift) and dual image coordinates (b, box key).
struct TypeSS {
    reps: BTreeMap<i32, Vec<(u32, u128)>>,
    duals: BTreeMap<i32, Vec<(u32, u64)>>,
}

fn sorted_spvec<F: Field>(acc: HashMap<u32, F>) -> SpVec<F> {
    let mut v: SpVec<F> = acc.into_iter().filter(|(_, x)| !x.is_zero()).collect();
    v.sort_unstable_by_key(|e| e.0);
    v
}

fn type_state_space<F: Field>(alg: &GrAlg<F>, tb: &TypeBasisC<F>) -> TypeSS {
    let mut degs: BTreeSet<i32> = BTreeSet::new();
    for &da in &tb.ins {
        for w in 0..=alg.top_weight as i32 {
            degs.insert(da + 2 * w);
        }
    }
    let mut reps = BTreeMap::new();
    let mut duals = BTreeMap::new();
    let mut buf = vec![];
    for d in degs {
        let mut rows: Vec<(u32, u128)> = vec![];
        let mut vecs: Vec<SpVec<F>> = vec![];
        let mut col: HashMap<(u32, u64), u32> = HashMap::new();
        let mut cols: Vec<(u32, u64)> = vec![];
        for (a, &da) in tb.ins.iter().enumerate() {
            let rd = d - da;
            if rd < 0 || rd % 2 != 0 {
                continue;
            }
            for &(_, m1) in alg.box_of_weight((rd / 2) as u32) {
                let mut acc: HashMap<u32, F> = HashMap::new();
                for (b, row) in tb.gram.iter().enumerate() {
                    for (mono, c) in &row[a].t {
                        alg.nf_prod(m1, *mono, &mut buf);
                        for (k, x) in &buf {
                            let ci = *col.entry((b as u32, *k)).or_insert_with(|| {
                                cols.push((b as u32, *k));
                                (cols.len() - 1) as u32
                            });
                            acc.entry(ci).or_insert_with(F::zero).fma(c, x);
                        }
                    }
                }
                let v = sorted_spvec(acc);
                if !v.is_empty() {
                    rows.push((a as u32, m1));
                    vecs.push(v);
                }
            }
        }
        if rows.is_empty() {
            continue;
        }
        let mut ech = RankEchelon::<F>::new(cols.len());
        let mut sel = vec![];
        for (k, v) in vecs.iter().enumerate() {
            if ech.insert(v) {
                sel.push(k);
            }
        }
        if sel.is_empty() {
            continue;
        }
        // independent columns of the selected rows
        let r = sel.len();
        let mut tcols: Vec<SpVec<F>> = vec![vec![]; cols.len()];
        for (ri, &k) in sel.iter().enumerate() {
            for (c, x) in &vecs[k] {
                tcols[*c as usize].push((ri as u32, x.clone()));
            }
        }
        let mut ech2 = RankEchelon::<F>::new(r);
        let mut dsel = vec![];
        for (c, v) in tcols.iter().enumerate() {
            if ech2.insert(v) {
                dsel.push(cols[c]);
                if dsel.len() == r {
                    break;
                }
            }
        }
        assert_eq!(dsel.len(), r, "state space: row and column ranks differ");
        reps.insert(d, sel.iter().map(|&k| rows[k]).collect::<Vec<_>>());
        duals.insert(d, dsel);
    }
    TypeSS { reps, duals }
}

fn par_map<T: Sync, R: Send, Fun: Fn(&T) -> R + Sync>(items: &[T], f: Fun) -> Vec<R> {
    let nt = crate::num_threads().min(items.len()).max(1);
    if nt == 1 {
        return items.iter().map(&f).collect();
    }
    let next = std::sync::atomic::AtomicUsize::new(0);
    let mut out: Vec<Option<R>> = (0..items.len()).map(|_| None).collect();
    let results: Vec<Vec<(usize, R)>> = std::thread::scope(|s| {
        let hs: Vec<_> = (0..nt)
            .map(|_| {
                s.spawn(|| {
                    let mut v = vec![];
                    loop {
                        let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        if i >= items.len() {
                            break;
                        }
                        v.push((i, f(&items[i])));
                    }
                    v
                })
            })
            .collect();
        hs.into_iter().map(|h| h.join().unwrap()).collect()
    });
    for v in results {
        for (i, r) in v {
            out[i] = Some(r);
        }
    }
    out.into_iter().map(|x| x.unwrap()).collect()
}

/// Chain-level dimension dim Φ_N(B_p) of every type (exact oracle state spaces).
pub fn exact_type_dims_colored<F: Field>(data: &ClosureDataC<F>, nn: usize, marked: Option<usize>) -> BTreeMap<u32, usize> {
    if data.lam.iter().any(|&l| l as usize > nn) {
        return data.types.keys().map(|&x| (x, 0)).collect();
    }
    let alg = GrAlg::<F>::new(&data.lam, nn, marked);
    data.types.iter().map(|(&x, tb)| (x, type_state_space(&alg, tb).reps.values().map(|v| v.len()).sum())).collect()
}

/// Colored exact closure (the oracle): gl_N Poincaré polynomial of the closed complex.
/// `marked = Some(m)`: reduced theory with block m marked. A label > N gives 0.
pub fn close_exact_colored<F: Field>(data: &ClosureDataC<F>, nn: usize, marked: Option<usize>) -> Poincare {
    data.validate().unwrap_or_else(|e| panic!("close_exact_colored: invalid closure data: {}", e));
    if data.lam.iter().any(|&l| l as usize > nn) {
        return Poincare::default();
    }
    let alg = GrAlg::<F>::new(&data.lam, nn, marked);
    let gq = data.global_q(nn, marked);
    let tss = std::time::Instant::now();
    let ss: BTreeMap<u32, TypeSS> = data.types.iter().map(|(&x, tb)| (x, type_state_space(&alg, tb))).collect();
    let verbose = std::env::var("KRHOM_VERBOSE").map(|v| v != "0" && !v.is_empty()).unwrap_or(false);
    if verbose {
        let dim: usize = ss.values().map(|s| s.reps.values().map(|v| v.len()).sum::<usize>()).sum();
        eprintln!("[cclose N={}] state spaces: {} types, total dim {}, {:.2}s", nn, ss.len(), dim, tss.elapsed().as_secs_f64());
    }
    let trk = std::time::Instant::now();
    // groups (h, q) -> reps (object, degree, index)
    let mut groups: BTreeMap<(i32, i32), Vec<(usize, i32, usize)>> = BTreeMap::new();
    for (i, o) in data.objs.iter().enumerate() {
        for (d, reps) in &ss[&o.0].reps {
            for k in 0..reps.len() {
                groups.entry((o.1, d + o.2 + gq)).or_default().push((i, *d, k));
            }
        }
    }
    // dual columns per group
    // Columns sorted by the box element in lex-descending order of the per-block indices
    // (block 0 most significant), as `close_exact` does for λ = 1^n: the leftmost entry of a
    // row m1·P is then (mostly) its leading term, which keeps RankEchelon's fill-in small.
    let mut dual_lists: BTreeMap<(i32, i32), Vec<(u64, u32, u32, u64)>> = BTreeMap::new();
    for (j, o) in data.objs.iter().enumerate() {
        for (d, duals) in &ss[&o.0].duals {
            for &(b, key) in duals {
                dual_lists.entry((o.1, d + o.2 + gq)).or_default().push((alg.order_key(key), j as u32, b, key));
            }
        }
    }
    let mut dual_size: HashMap<(i32, i32), usize> = HashMap::new();
    let mut dual_index: FxMap<(u32, u32, u64), u32> = FxMap::default();
    for (g, mut list) in dual_lists {
        list.sort_unstable();
        for (c, e) in list.iter().enumerate() {
            dual_index.insert((e.1, e.2, e.3), c as u32);
        }
        dual_size.insert(g, list.len());
    }
    // P_ab = Σ_c G_t[b][c] F[c][a] per edge
    let mut out_edges: HashMap<usize, Vec<usize>> = HashMap::new();
    for (ei, (s, _, _)) in data.edges.iter().enumerate() {
        out_edges.entry(*s).or_default().push(ei);
    }
    let pab: Vec<Vec<Vec<EPoly<F>>>> = data
        .edges
        .iter()
        .map(|(s, t, f)| {
            let gt = &data.types[&data.objs[*t].0];
            let na = data.types[&data.objs[*s].0].ins.len();
            (0..na)
                .map(|a| {
                    (0..gt.outs.len())
                        .map(|b| {
                            let mut acc = EPoly::zero();
                            for c in 0..gt.ins.len() {
                                acc.add_assign(&gt.gram[b][c].mul(&f[c][a]));
                            }
                            acc
                        })
                        .collect()
                })
                .collect()
        })
        .collect();
    let blocks: Vec<((i32, i32), &Vec<(usize, i32, usize)>)> = groups.iter().map(|(k, v)| (*k, v)).collect();
    let ranks: Vec<usize> = par_map(&blocks, |((h, q), elems)| {
        let tl = dual_size.get(&(h + 1, *q)).copied().unwrap_or(0);
        if tl == 0 {
            return 0;
        }
        let mut ech = RankEchelon::<F>::new(tl);
        let mut acc = Acc::<F>::new(tl);
        let mut buf = vec![];
        for &(i, d, k) in elems.iter() {
            let (a, m1) = ss[&data.objs[i].0].reps[&d][k];
            for &ei in out_edges.get(&i).map(|v| v.as_slice()).unwrap_or(&[]) {
                let t = data.edges[ei].1;
                for (b, p) in pab[ei][a as usize].iter().enumerate() {
                    for (mono, c) in &p.t {
                        alg.nf_prod(m1, *mono, &mut buf);
                        for (key, x) in &buf {
                            if let Some(&col) = dual_index.get(&(t as u32, b as u32, *key)) {
                                acc.add(col, &c.mul(x));
                            }
                        }
                    }
                }
            }
            ech.insert(&acc.take());
        }
        ech.rank()
    });
    let rank: HashMap<(i32, i32), usize> = blocks.iter().map(|(k, _)| *k).zip(ranks).collect();
    if verbose {
        eprintln!("[cclose N={}] chain dim {}, ranks {:.2}s", nn, groups.values().map(|v| v.len()).sum::<usize>(), trk.elapsed().as_secs_f64());
    }
    let mut pc = Poincare::default();
    for ((h, q), elems) in &groups {
        let r_out = rank.get(&(*h, *q)).copied().unwrap_or(0);
        let r_in = rank.get(&(h - 1, *q)).copied().unwrap_or(0);
        pc.add_term(*h, *q, elems.len() as i128 - r_out as i128 - r_in as i128);
    }
    pc
}

// ------------------------------------------------------------------ (Δ_λ G, F) adapter

/// (a, N) pairs for which gate G4 certified that Δ_a R_N is free over H*(Gr(a, N)) of graded
/// rank t^{l(w_a)} π_a(t): output of `g4_free` (all 2 ≤ a ≤ 5, a ≤ N ≤ 10; certificate
/// I_N·Δ = 0 in R_N and staircase leading terms, i.e. dim Δ_a R_N = N!/(N−a)!; C2's direct
/// generator count agrees on every pair, see the module doc of src/bin/g4_free.rs). a ≤ 1 is
/// free trivially (Δ = 1). Pairs not listed here are refused by [`close_delta_tower`].
pub const G4_CERTIFIED: &[(u8, u8)] = &[
    (2, 2), (2, 3), (2, 4), (2, 5), (2, 6), (2, 7), (2, 8), (2, 9), (2, 10),
    (3, 3), (3, 4), (3, 5), (3, 6), (3, 7), (3, 8), (3, 9), (3, 10),
    (4, 4), (4, 5), (4, 6), (4, 7), (4, 8), (4, 9), (4, 10),
    (5, 5), (5, 6), (5, 7), (5, 8), (5, 9), (5, 10),
];

/// True if the adapter may use blocks of label a at gl_N.
pub fn g4_certified(a: usize, nn: usize) -> bool {
    a <= 1 || G4_CERTIFIED.iter().any(|&(x, y)| x as usize == a && y as usize == nn)
}

/// gr(V_a): the graded rank of Δ_a R_N over H*(Gr(a, N)) in the tower's coordinate degrees,
/// Π_{k=1}^{a} (1 + q^2 + … + q^{2(k−1)}). The tower's coordinate degrees do not see deg Δ
/// (the coordinate of v is v, the relation matrix is Δ G), so the t^{l(w_a)} of C2 M1(b)
/// (generator degrees in Δ_a R_N) cancels; the remaining shift is 0, as checked on single
/// types against the exact oracle (tests/cclose.rs).
pub fn gr_v(a: usize) -> LPoly {
    let mut r = LPoly::one();
    for k in 1..=a {
        let mut f = LPoly::zero();
        for j in 0..k {
            f.add_term(2 * j as i32, 1);
        }
        r = r.mul(&f);
    }
    r
}

/// Exact division of Laurent polynomials (d with leading coefficient ±1).
pub fn lpoly_div(p: &LPoly, d: &LPoly) -> Result<LPoly, String> {
    let (dmax, dlead) = d.c.iter().next_back().map(|(&e, &k)| (e, k)).ok_or("division by 0")?;
    if dlead != 1 && dlead != -1 {
        return Err("divisor needs a unit leading coefficient".into());
    }
    let dmin = d.min_deg().unwrap();
    let mut rem = p.clone();
    let mut quo = LPoly::zero();
    while let Some((&rmax, &rlead)) = rem.c.iter().next_back() {
        let rmin = rem.min_deg().unwrap();
        if rmax - dmax < rmin - dmin {
            return Err(format!("not divisible: remainder {}", rem));
        }
        let e = rmax - dmax;
        let k = rlead * dlead;
        quo.add_term(e, k);
        for (&de, &dk) in &d.c {
            rem.add_term(de + e, -dk * k);
        }
    }
    Ok(quo)
}

type RPoly<F> = HashMap<Vec<u32>, F>;

/// Expansion of R^λ in root variables for the adapter.
struct RootMap<F: Field> {
    ring: ERing,
    marked: Option<usize>,
    /// first root of each block (unmarked: λ_i roots; marked: one dummy root with cap 0)
    off: Vec<usize>,
    nroots: usize,
    caps: Vec<u32>,
    /// per block: β (block exponents) -> e^β(roots of the block), local exponents, truncated
    memo: Vec<HashMap<Vec<u32>, Vec<(Vec<u32>, F)>>>,
}

impl<F: Field> RootMap<F> {
    fn new(lam: &[u8], nn: usize, marked: Option<usize>) -> Self {
        let ring = ERing::new(lam);
        let mut off = vec![];
        let mut caps = vec![];
        for (i, &l) in lam.iter().enumerate() {
            off.push(caps.len());
            if Some(i) == marked {
                caps.push(0);
            } else {
                for _ in 0..l {
                    caps.push(nn as u32 - 1);
                }
            }
        }
        let nroots = caps.len();
        RootMap { memo: vec![HashMap::new(); lam.len()], ring, marked, off, nroots, caps }
    }
    fn local_mul(a: &[(Vec<u32>, F)], b: &[(Vec<u32>, F)], cap: u32) -> Vec<(Vec<u32>, F)> {
        let mut h: HashMap<Vec<u32>, F> = HashMap::new();
        for (ea, ca) in a {
            for (eb, cb) in b {
                let e: Vec<u32> = ea.iter().zip(eb.iter()).map(|(x, y)| x + y).collect();
                if e.iter().any(|&x| x > cap) {
                    continue;
                }
                h.entry(e).or_insert_with(F::zero).fma(ca, cb);
            }
        }
        let mut v: Vec<(Vec<u32>, F)> = h.into_iter().filter(|(_, c)| !c.is_zero()).collect();
        v.sort_unstable_by(|x, y| x.0.cmp(&y.0));
        v
    }
    /// e^β of block i in its λ_i local roots, exponents > cap dropped.
    fn block_poly(&mut self, i: usize, be: &[u32]) -> Vec<(Vec<u32>, F)> {
        if let Some(v) = self.memo[i].get(be) {
            return v.clone();
        }
        let a = self.ring.lam[i] as usize;
        let cap = self.caps[self.off[i]];
        let v = if be.iter().all(|&x| x == 0) {
            vec![(vec![0u32; a], F::one())]
        } else {
            let j = be.iter().position(|&x| x > 0).unwrap();
            let mut prev = be.to_vec();
            prev[j] -= 1;
            let pv = self.block_poly(i, &prev);
            // e_{j+1}(roots): all (j+1)-subsets
            let mut ej = vec![];
            for mask in 0u32..(1 << a) {
                if mask.count_ones() as usize == j + 1 {
                    ej.push(((0..a).map(|r| (mask >> r) & 1).collect::<Vec<u32>>(), F::one()));
                }
            }
            Self::local_mul(&pv, &ej, cap)
        };
        self.memo[i].insert(be.to_vec(), v.clone());
        v
    }
    fn expand(&mut self, p: &EPoly<F>) -> RPoly<F> {
        let mut out: RPoly<F> = HashMap::new();
        let nb = self.ring.lam.len();
        'terms: for (m, c) in &p.t {
            let mut parts: Vec<(usize, Vec<(Vec<u32>, F)>)> = vec![];
            for i in 0..nb {
                let be = self.ring.block_exps(*m, i);
                if Some(i) == self.marked {
                    if be.iter().any(|&x| x > 0) {
                        continue 'terms; // e_{m,*} = 0
                    }
                    continue;
                }
                if be.is_empty() {
                    continue;
                }
                let bp = self.block_poly(i, &be);
                if bp.is_empty() {
                    continue 'terms;
                }
                parts.push((i, bp));
            }
            let mut acc: Vec<(Vec<u32>, F)> = vec![(vec![0u32; self.nroots], c.clone())];
            for (i, bp) in &parts {
                let o = self.off[*i];
                let mut next = Vec::with_capacity(acc.len() * bp.len());
                for (e0, c0) in &acc {
                    for (e1, c1) in bp {
                        let mut e = e0.clone();
                        e[o..o + e1.len()].copy_from_slice(e1);
                        next.push((e, c0.mul(c1)));
                    }
                }
                acc = next;
            }
            for (e, x) in acc {
                out.entry(e).or_insert_with(F::zero).add_assign(&x);
            }
        }
        out.retain(|_, c| !c.is_zero());
        out
    }
    /// Δ_λ = Π over unmarked blocks Π_{p<q} (x_p − x_q), truncated.
    fn delta(&self) -> RPoly<F> {
        let mut d: RPoly<F> = HashMap::from([(vec![0u32; self.nroots], F::one())]);
        for (i, &l) in self.ring.lam.iter().enumerate() {
            if Some(i) == self.marked {
                continue;
            }
            let o = self.off[i];
            for p in 0..l as usize {
                for q in (p + 1)..l as usize {
                    let mut ep = vec![0u32; self.nroots];
                    ep[o + p] = 1;
                    let mut eq = vec![0u32; self.nroots];
                    eq[o + q] = 1;
                    let f: RPoly<F> = HashMap::from([(ep, F::one()), (eq, F::one().neg())]);
                    d = self.rmul(&d, &f);
                }
            }
        }
        d
    }
    fn rmul(&self, a: &RPoly<F>, b: &RPoly<F>) -> RPoly<F> {
        let mut h: RPoly<F> = HashMap::new();
        for (ea, ca) in a {
            for (eb, cb) in b {
                let e: Vec<u32> = ea.iter().zip(eb.iter()).map(|(x, y)| x + y).collect();
                if e.iter().zip(self.caps.iter()).any(|(x, c)| x > c) {
                    continue;
                }
                h.entry(e).or_insert_with(F::zero).fma(ca, cb);
            }
        }
        h.retain(|_, c| !c.is_zero());
        h
    }
}

/// The (Δ_λ G, F) adapter: `tower::close_tower` on root-variable data, divided by
/// Π_{unmarked i} gr(V_{λ_i}). Refuses (Err) if some unmarked block label a ≥ 2 has (a, N)
/// not certified by gate G4, if more than `tower::MAX_VARS` roots are needed, or if the
/// division is not exact.
pub fn close_delta_tower<F: Field>(data: &ClosureDataC<F>, nn: usize, marked: Option<usize>, opts: TowerOpts) -> Result<Poincare, String> {
    for (i, &l) in data.lam.iter().enumerate() {
        if Some(i) != marked && (l as usize) <= nn && !g4_certified(l as usize, nn) {
            return Err(format!(
                "(Delta G) adapter refused: freeness of Delta_a R_N over H*(Gr(a,N)) is not certified for (a, N) = ({}, {}) (gate G4); run src/bin/g4_free.rs and add the pair to cclose::G4_CERTIFIED",
                l, nn
            ));
        }
    }
    close_delta_tower_uncertified(data, nn, marked, opts)
}

/// As [`close_delta_tower`] without the G4 check; prints a loud warning for every
/// uncertified (a, N). For experiments only: the result is correct only if Δ_a R_N is free.
pub fn close_delta_tower_uncertified<F: Field>(data: &ClosureDataC<F>, nn: usize, marked: Option<usize>, opts: TowerOpts) -> Result<Poincare, String> {
    data.validate()?;
    if data.lam.iter().any(|&l| l as usize > nn) {
        return Ok(Poincare::default());
    }
    if let Some(m) = marked {
        if m >= data.lam.len() || data.lam[m] == 0 {
            return Err(format!("marked block {} invalid for {:?}", m, data.lam));
        }
    }
    for (i, &l) in data.lam.iter().enumerate() {
        if Some(i) != marked && !g4_certified(l as usize, nn) {
            eprintln!(
                "WARNING [cclose] (Delta G) adapter used with UNCERTIFIED (a, N) = ({}, {}): the result is correct only if Delta_a R_N is free over H*(Gr(a,N)) (gate G4)",
                l, nn
            );
        }
    }
    if nn > 32 {
        return Err(format!("(Delta G) adapter: N = {} > 32 (tower caps)", nn));
    }
    let mut rm = RootMap::<F>::new(&data.lam, nn, marked);
    if rm.nroots == 0 {
        return Err("(Delta G) adapter: no root variables (all labels 0)".into());
    }
    if rm.nroots > crate::tower::MAX_VARS {
        return Err(format!("(Delta G) adapter: {} root variables > {} (tower::MAX_VARS)", rm.nroots, crate::tower::MAX_VARS));
    }
    let delta = rm.delta();
    let mut bases: BTreeMap<usize, TypeBasis<F>> = BTreeMap::new();
    for (&x, tb) in &data.types {
        let gram: Vec<Vec<RPoly<F>>> = tb.gram.iter().map(|row| row.iter().map(|p| { let r = rm.expand(p); rm.rmul(&delta, &r) }).collect()).collect();
        bases.insert(
            x as usize,
            TypeBasis { ins: tb.ins.iter().enumerate().map(|(a, &d)| (a, d)).collect(), outs: tb.outs.iter().enumerate().map(|(b, &d)| (b, d)).collect(), gram },
        );
    }
    let objs: Vec<Obj> = data.objs.iter().map(|&(x, h, q)| Obj { x: x as usize, q, h }).collect();
    let edges: Vec<(usize, usize)> = data.edges.iter().map(|(s, t, _)| (*s, *t)).collect();
    let fpolys: Vec<Vec<Vec<RPoly<F>>>> = data.edges.iter().map(|(_, _, f)| f.iter().map(|row| row.iter().map(|p| rm.expand(p)).collect()).collect()).collect();
    let gq = data.global_q(nn, marked);
    let r = crate::tower::close_tower(rm.nroots, &rm.caps, gq, &objs, &edges, &fpolys, &bases, opts);
    let mut div = LPoly::one();
    for (i, &l) in data.lam.iter().enumerate() {
        if Some(i) != marked {
            div = div.mul(&gr_v(l as usize));
        }
    }
    let mut rows: BTreeMap<i32, LPoly> = BTreeMap::new();
    for (&(t, q), &k) in &r.poincare.c {
        rows.entry(t).or_default().add_term(q, k);
    }
    let mut out = Poincare::default();
    for (t, p) in rows {
        let quo = lpoly_div(&p, &div).map_err(|e| format!("(Delta G) adapter: t^{} row {} not divisible by gr(V) = {}: {}", t, p, div, e))?;
        for (&q, &k) in &quo.c {
            if k < 0 {
                return Err(format!("(Delta G) adapter: negative coefficient after division at t^{} q^{}", t, q));
            }
            out.add_term(t, q, k);
        }
    }
    Ok(out)
}

// ------------------------------------------------------------------ exact closure data of closed ladders

/// R^λ-bases of Hom(R^λ, B) ("ins": fsize×1 maps) and Hom(B, R^λ) ("outs": 1×fsize maps)
/// with their unnormalized degrees.
pub struct LadderBasis<F: Field> {
    pub ins: Vec<(i32, PMat<F>)>,
    pub outs: Vec<(i32, PMat<F>)>,
    /// number of identity-coset flows of the ladder (= the expected R^λ-rank of both Homs)
    pub rank: usize,
}

/// Number of flows of the ladder whose top configuration is the bottom one (identity coset):
/// the rank of Hom(R^λ, B) over R^λ (its localization is the identity-coset block).
pub fn identity_coset_flows<F: Field>(b: &SBim<F>) -> usize {
    let pts = crate::local::fixed_points(b);
    let off = crate::local::root_offsets(&b.bot_labels);
    let bottom: Vec<u64> = b.bot_labels.iter().enumerate().map(|(i, &l)| ((1u64 << l) - 1) << off[i]).collect();
    pts.iter().filter(|p| p.top == bottom).count()
}

/// Flatten a map R^λ → B (fsize×1) or B → R^λ (1×fsize) to (fibre index, polynomial).
fn flat<F: Field>(m: &PMat<F>, out: bool) -> Vec<(u32, &Poly<F>)> {
    if out {
        m.cols.iter().enumerate().filter_map(|(be, c)| c.first().map(|(_, p)| (be as u32, p))).collect()
    } else {
        m.cols[0].iter().map(|(be, p)| (*be, p)).collect()
    }
}

/// Graded Nakayama over R^λ on the exact Hom spaces of `hom_space`, degree by degree, until
/// `target` generators are found. At every processed degree (and two more after the last
/// generator) the R^λ-multiples of the chosen generators must be linearly independent and
/// span Hom_d (freeness certificate in those degrees).
fn nakayama<F: Field>(id: &SBim<F>, b: &SBim<F>, out: bool, target: usize, mc: &mut MonoCache) -> Result<Vec<(i32, PMat<F>)>, String> {
    let maxf = *b.fdeg.iter().max().unwrap();
    let k: i32 = b.bot_labels.iter().sum::<usize>() as i32;
    let (dmin, dmax) = if out { (-maxf, maxf + 4 * k * k + 8) } else { (0, maxf + 4 * k * k + 8) };
    let mut gens: Vec<(i32, PMat<F>)> = vec![];
    let mut extra = 0;
    let mut d = dmin;
    while d <= dmax {
        // coordinate space of degree d: (fibre index, monomial)
        let mut idx: HashMap<(u32, u128), u32> = HashMap::new();
        for be in 0..b.fsize() {
            let md = if out { d + b.fdeg[be] } else { d - b.fdeg[be] };
            if md < 0 {
                continue;
            }
            let ms = mc.get(md).0.clone();
            for m in ms {
                let n = idx.len() as u32;
                idx.insert((be as u32, m), n);
            }
        }
        let homs = if out { hom_space(b, id, d, mc) } else { hom_space(id, b, d, mc) };
        let tovec = |entries: &[(u32, Poly<F>)]| -> SpVec<F> {
            let mut v: SpVec<F> = vec![];
            for (be, p) in entries {
                for (m, c) in &p.t {
                    let i = idx.get(&(*be, *m)).unwrap_or_else(|| panic!("Nakayama: coordinate outside degree {}", d));
                    v.push((*i, c.clone()));
                }
            }
            v.sort_unstable_by_key(|e| e.0);
            v
        };
        let mut ech = Echelon::<F>::new(idx.len().max(1));
        let mut span = 0usize;
        for (dg, g) in &gens {
            let rd = d - dg;
            if rd < 0 {
                continue;
            }
            let ms = mc.get(rd).0.clone();
            for m in ms {
                let entries: Vec<(u32, Poly<F>)> = flat(g, out).into_iter().map(|(be, p)| (be, p.mul_mono(m, &F::one()))).collect();
                if !ech.insert(&tovec(&entries)) {
                    return Err(format!("Nakayama: R^lambda-multiples of the generators are dependent in degree {} (Hom not free)", d));
                }
                span += 1;
            }
        }
        let mut new = 0;
        for h in homs.iter() {
            let entries: Vec<(u32, Poly<F>)> = flat(h, out).into_iter().map(|(be, p)| (be, p.clone())).collect();
            if ech.insert(&tovec(&entries)) {
                gens.push((d, h.clone()));
                new += 1;
            }
        }
        if span + new != homs.len() {
            return Err(format!("Nakayama: degree {}: span {} + new {} != dim {}", d, span, new, homs.len()));
        }
        if gens.len() > target {
            return Err(format!("Nakayama: more than {} generators (degree {})", target, d));
        }
        if gens.len() == target {
            if new > 0 {
                extra = 0;
            }
            extra += 1;
            if extra > 4 {
                return Ok(gens);
            }
        }
        d += 1;
    }
    Err(format!("Nakayama: only {} of {} generators up to degree {}", gens.len(), target, dmax))
}

/// Exact R^λ-bases of Hom(R^λ, B) and Hom(B, R^λ) for a closed ladder bimodule B.
pub fn rlambda_bases<F: Field>(b: &SBim<F>) -> Result<LadderBasis<F>, String> {
    if b.top_labels != b.bot_labels {
        return Err(format!("ladder not closed: bottom {:?} top {:?}", b.bot_labels, b.top_labels));
    }
    let id = SBim::<F>::identity(b.n, &b.bot_labels);
    let mut mc = MonoCache::new(b.ring.clone());
    let rank = identity_coset_flows(b);
    let ins = nakayama(&id, b, false, rank, &mut mc)?;
    let outs = nakayama(&id, b, true, rank, &mut mc)?;
    Ok(LadderBasis { ins, outs, rank })
}

/// One-point determinant criterion (C2 M2): the identity-coset values of the in- (and out-)
/// generators at the point x (one root per cable) form nonsingular n×n matrices.
pub fn det_criterion<F: Field>(b: &SBim<F>, basis: &LadderBasis<F>, x: &[F]) -> (bool, bool) {
    use crate::local::{eval_bottom, fixed_points, inverse, loc_matrix, root_offsets};
    let pts = fixed_points(b);
    let off = root_offsets(&b.bot_labels);
    let bottom: Vec<u64> = b.bot_labels.iter().enumerate().map(|(i, &l)| ((1u64 << l) - 1) << off[i]).collect();
    let idp: Vec<usize> = (0..pts.len()).filter(|&p| pts[p].top == bottom).collect();
    let n = idp.len();
    let l = loc_matrix(b, &pts, x);
    let linv = match inverse(&l) {
        Some(m) => m,
        None => return (false, false),
    };
    let ev = |m: &PMat<F>, out: bool| -> Vec<F> {
        let mut v = vec![F::zero(); b.fsize()];
        for (be, p) in flat(m, out) {
            v[be as usize] = eval_bottom(b, p, x);
        }
        v
    };
    let ins: Vec<Vec<F>> = basis
        .ins
        .iter()
        .map(|(_, g)| {
            let gv = ev(g, false);
            idp.iter().map(|&p| (0..gv.len()).fold(F::zero(), |mut s, be| { s.fma(&l[p][be], &gv[be]); s })).collect()
        })
        .collect();
    let outs: Vec<Vec<F>> = basis
        .outs
        .iter()
        .map(|(_, h)| {
            let hv = ev(h, true);
            idp.iter().map(|&p| (0..hv.len()).fold(F::zero(), |mut s, be| { s.fma(&hv[be], &linv[be][p]); s })).collect()
        })
        .collect();
    let ok = |m: &Vec<Vec<F>>| m.len() == n && crate::dense::rank(m, n) == n;
    (ok(&ins), ok(&outs))
}

/// G[b][a] = h_b ∘ g_a in R^λ.
pub fn gram_of<F: Field>(basis: &LadderBasis<F>, nvars: usize) -> Vec<Vec<EPoly<F>>> {
    basis
        .outs
        .iter()
        .map(|(_, h)| {
            basis
                .ins
                .iter()
                .map(|(_, g)| {
                    let c = h.matmul(g);
                    EPoly::from_mpoly(&c.get(0, 0), nvars)
                })
                .collect()
        })
        .collect()
}

/// Exact F for a bimodule map f: B → B' of unnormalized degree `deg_f`:
/// f ∘ g_a = Σ_c F_ca g'_c with F_ca ∈ R^λ (unique, since the g'_c are an R^λ-basis).
pub fn solve_f<F: Field>(f: &PMat<F>, deg_f: i32, src: &LadderBasis<F>, tgt: &LadderBasis<F>, bp: &SBim<F>) -> Result<Vec<Vec<EPoly<F>>>, String> {
    let nvars = bp.ring.nvars;
    let mut mc = MonoCache::new(bp.ring.clone());
    let mut out = vec![vec![EPoly::zero(); src.ins.len()]; tgt.ins.len()];
    for (a, (da, ga)) in src.ins.iter().enumerate() {
        let lhs = f.matmul(ga); // fsize' × 1
        let d = deg_f + da;
        // unknowns (c, monomial of degree d − d'_c)
        let mut unk: Vec<(usize, u128)> = vec![];
        for (c, (dc, _)) in tgt.ins.iter().enumerate() {
            if d - dc < 0 {
                continue;
            }
            for &m in &mc.get(d - dc).0.clone() {
                unk.push((c, m));
            }
        }
        // coordinates (β', monomial)
        let mut idx: HashMap<(u32, u128), usize> = HashMap::new();
        let mut cols: Vec<Vec<((u32, u128), F)>> = vec![];
        for &(c, m) in &unk {
            let g = &tgt.ins[c].1;
            let v: Vec<((u32, u128), F)> = g.cols[0].iter().flat_map(|(be, p)| p.t.iter().map(move |(mm, x)| ((*be, mm + m), x.clone()))).collect();
            for (k, _) in &v {
                let n = idx.len();
                idx.entry(*k).or_insert(n);
            }
            cols.push(v);
        }
        let mut rhs_terms: Vec<((u32, u128), F)> = vec![];
        for (be, p) in &lhs.cols[0] {
            for (mm, x) in &p.t {
                rhs_terms.push(((*be, *mm), x.clone()));
            }
        }
        for (k, _) in &rhs_terms {
            let n = idx.len();
            idx.entry(*k).or_insert(n);
        }
        let nrow = idx.len();
        let mut rows = vec![vec![F::zero(); unk.len()]; nrow];
        for (u, v) in cols.iter().enumerate() {
            for (k, x) in v {
                rows[idx[k]][u].add_assign(x);
            }
        }
        let mut rhs = vec![F::zero(); nrow];
        for (k, x) in &rhs_terms {
            rhs[idx[k]].add_assign(x);
        }
        if unk.is_empty() {
            if rhs.iter().any(|x| !x.is_zero()) {
                return Err(format!("solve_f: f∘g_{} is not in the span (no unknowns)", a));
            }
            continue;
        }
        let sol = crate::dense::solve(&rows, unk.len(), &rhs).ok_or_else(|| format!("solve_f: f∘g_{} not in the R^lambda-span of the target basis", a))?;
        let mut per_c: Vec<Vec<(u128, F)>> = vec![vec![]; tgt.ins.len()];
        for (u, x) in sol.into_iter().enumerate() {
            if !x.is_zero() {
                per_c[unk[u].0].push((unk[u].1, x));
            }
        }
        for (c, terms) in per_c.into_iter().enumerate() {
            let p = Poly { t: { let mut t = terms; t.sort_unstable_by_key(|e| e.0); t } };
            out[c][a] = EPoly::from_mpoly(&p, nvars);
        }
    }
    Ok(out)
}

/// Closure data of a single closed ladder (one object at h = 0, q = MOY normalization shift
/// −Σ_rungs r(b−r)), from exact R^λ-bases (graded Nakayama) and Gram by composition.
/// The ladder's `n` only bounds the labels; the data are N-independent.
pub fn closure_data_from_ladder<F: Field>(lad: &Ladder) -> Result<(ClosureDataC<F>, LadderBasis<F>), String> {
    let b = SBim::<F>::identity(lad.n, &lad.bottom).add_rungs(&lad.rungs);
    let basis = rlambda_bases(&b)?;
    let gram = gram_of(&basis, b.ring.nvars);
    let tb = TypeBasisC { ins: basis.ins.iter().map(|g| g.0).collect(), outs: basis.outs.iter().map(|g| g.0).collect(), gram };
    let data = ClosureDataC { lam: lad.bottom.iter().map(|&l| l as u8).collect(), objs: vec![(0, 0, b.shift)], types: BTreeMap::from([(0u32, tb)]), edges: vec![], crossings: vec![] };
    data.validate()?;
    Ok((data, basis))
}
