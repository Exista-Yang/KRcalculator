//! MBOX: the closure layer at a fixed N ("box closure") for the localized multi-strand engine
//! (SPEC_3C_GL6 rev1 §2 MBOX, D5, D7, D14; R3_closure §2-6; R4_m6a §2; C_sound M1, M2, M5).
//!
//! **What is extracted.** The closure reads Gram entries G[b][a] and closed pairings
//! P[b][a] = (h^t_b ∘ d ∘ g^s_a) only through their images in A_λ = ⊗_i H*(Gr(λ_i, N)) (R3 §2).
//! An entry is a homogeneous f ∈ R^λ = ⊗_i Sym_{λ_i} of declared weight w (q-degree 2w). Write
//! f = Σ_ρ c'_ρ s_ρ(x | a) in factorial Schur functions (per block i with its own nodes a^(i));
//! then c'_ρ = 0 for |ρ| > w and c'_ρ = c_ρ, the ordinary Schur coefficient, for |ρ| = w, so
//! f mod I_N = Σ_{ρ ⊆ box, |ρ| = w} c_ρ s_ρ.
//!
//! **Points.** For μ = (μ^i) with μ^i ⊆ λ_i × (N − λ_i), the point z_μ has, in block i, the roots
//! a^(i)_{μ^i_j + λ_i + 1 − j} (j = 1..λ_i, 1-based nodes); [`point_coords`] gives the values of
//! the e_{i,j} there. By the vanishing property s_ρ(z_μ | a) = 0 unless ρ ⊆ μ, and
//! s_μ(z_μ | a) ≠ 0, so the matrix S_{μρ} = s_ρ(z_μ | a) = ⊗_i S^(i) is triangular for
//! containment, and so is S^{-1}: c_ρ = Σ_{μ ⊆ ρ} (S^{-1})_{ρμ} f(z_μ) needs only the down-set of ρ.
//! [`BoxGrid`] holds S^(i) (asserted triangular at construction) and its inverse.
//!
//! **Streaming.** [`BoxAcc`] accumulates, per entry of weight w, the coefficients of |ρ| = w
//! (outputs) and |ρ| = w + 1 (residuals, which must vanish: a polynomiality check that catches
//! an underestimated weight) as values arrive point by point, in any order, in shards. D14 /
//! C_sound M2 guards: a point added twice is refused; checkpoints are written atomically (temp
//! file + rename) with a header (prime, node-seed hash, library hash, word hash, prune_n, mode,
//! N, code stamp, entry list hash) and a payload hash, plus the done-point bitmap; merging
//! refuses mismatched headers and overlapping bitmaps; [`BoxAcc::finish`] refuses an incomplete
//! down-set and nonzero residuals.
//!
//! **Output.** Schur coordinates per entry. [`schur_to_epoly`] turns them into an element of R^λ
//! (dual Jacobi-Trudi per block) restricted to N' ≤ N (ρ^i ⊆ λ_i × (N' − λ_i)) and to the reduced
//! theory (ρ^m = ∅); its image in A_λ at N' is f mod I_{N'}.
//!
//! **Degrees (C_sound M1).** [`check_scaling`]: z_μ(λa) = λ·z_μ(a), so a value of weight w
//! satisfies v(λz) = λ^w v(z). Asserting this against the declared weight catches an
//! overestimated weight, which box extraction would turn into silent zeros.
//!
//! **Closure (D5, C_sound M5).** [`BoxClosureData`] wraps a `ClosureDataC` (types and objects;
//! Gram from the box; no F) with the closed pairing P, prune_n, the box N, valid_n =
//! min(prune_n, N_box), the extraction mode and the D² certificates. [`close_box`] refuses
//! N' > valid_n, a mode the data cannot serve, and pruned data without a [`D2Cert`] for that
//! (N', mode); it runs the staged rank stage on P (`cclose::try_close_staged`), optionally
//! capped by MOY ranks from an uncapped T0 ([`box_moy_check`]) on the same library hash.
//! [`closed_d2_freivalds`] runs the Freivalds test D_{h+1} Π^{-1} D_h u = 0
//! (`cclose::closed_d2_freivalds_raw`) and records a [`D2Cert`] bound to the pairing hash.

use crate::cclose::{closed_d2_freivalds_raw, global_q, try_close_staged, try_staged_type_ranks, validate_pairing, ClosureDataC, EPoly, ERing, Pairing, StagedOpts, TypeBasisC, TypeRanks};
use crate::field::{Field, Fp};
use crate::laurent::{LPoly, Poincare};
use rand::{rngs::StdRng, Rng, SeedableRng};
use std::collections::BTreeMap;
use std::io::{Read, Write};
use std::path::Path;
use std::sync::Arc;

/// Version stamp of this module's algorithms and file format (part of every checkpoint header).
pub const BOX_CODE_STAMP: &str = "boxclose-v2";

// ------------------------------------------------------------------ fields and hashes

/// A prime field whose elements serialize as u64 (checkpoints, hashes).
pub trait PrimeField: Field {
    fn to_u64(&self) -> u64;
    fn from_u64(x: u64) -> Self;
}

impl<const P: u64> PrimeField for Fp<P> {
    fn to_u64(&self) -> u64 {
        self.0
    }
    fn from_u64(x: u64) -> Self {
        Fp::new(x)
    }
}

/// FNV-1a 64-bit hash (stable across runs and builds; used for every stamp in this module).
#[derive(Clone, Copy, Debug)]
pub struct Fnv(pub u64);

impl Default for Fnv {
    fn default() -> Self {
        Fnv(0xcbf2_9ce4_8422_2325)
    }
}

impl Fnv {
    pub fn bytes(&mut self, b: &[u8]) -> &mut Self {
        for &x in b {
            self.0 ^= x as u64;
            self.0 = self.0.wrapping_mul(0x0000_0100_0000_01b3);
        }
        self
    }
    pub fn u64(&mut self, x: u64) -> &mut Self {
        self.bytes(&x.to_le_bytes())
    }
    pub fn i64(&mut self, x: i64) -> &mut Self {
        self.u64(x as u64)
    }
    pub fn u128(&mut self, x: u128) -> &mut Self {
        self.bytes(&x.to_le_bytes())
    }
    pub fn str(&mut self, s: &str) -> &mut Self {
        self.u64(s.len() as u64).bytes(s.as_bytes())
    }
    pub fn get(&self) -> u64 {
        self.0
    }
}

fn hash_epoly<F: PrimeField>(h: &mut Fnv, p: &EPoly<F>) {
    h.u64(p.t.len() as u64);
    for (m, c) in &p.t {
        h.u128(*m).u64(c.to_u64());
    }
}

/// Hash of a type library (ins, outs and Gram entries of every type; field included). T0
/// records and capped closures are bound to it (SPEC D14 "T0").
pub fn library_hash<F: PrimeField>(types: &BTreeMap<u32, TypeBasisC<F>>) -> u64 {
    let mut h = Fnv::default();
    h.str("library").u64(F::characteristic());
    for (x, tb) in types {
        h.u64(*x as u64).u64(tb.ins.len() as u64).u64(tb.outs.len() as u64);
        for d in tb.ins.iter().chain(tb.outs.iter()) {
            h.i64(*d as i64);
        }
        for row in &tb.gram {
            for p in row {
                hash_epoly(&mut h, p);
            }
        }
    }
    h.get()
}

/// Hash of a closed pairing (D² certificates are bound to it).
pub fn pairing_hash<F: PrimeField>(pairing: &[(usize, usize, Vec<Vec<EPoly<F>>>)]) -> u64 {
    let mut h = Fnv::default();
    h.str("pairing").u64(F::characteristic());
    for (s, t, p) in pairing {
        h.u64(*s as u64).u64(*t as u64).u64(p.len() as u64);
        for row in p {
            h.u64(row.len() as u64);
            for x in row {
                hash_epoly(&mut h, x);
            }
        }
    }
    h.get()
}

/// Hash of a list of braid words (checkpoint header field).
pub fn word_hash(words: &[Vec<i32>]) -> u64 {
    let mut h = Fnv::default();
    h.str("words");
    for w in words {
        h.u64(w.len() as u64);
        for &g in w {
            h.i64(g as i64);
        }
    }
    h.get()
}

// ------------------------------------------------------------------ nodes and points

/// Interpolation nodes: block i gets n nodes a^(i)_1..a^(i)_n, all distinct and nonzero across
/// the blocks, drawn from a seeded generator (a fresh seed per prime and per node sequence).
#[derive(Clone, Debug)]
pub struct BoxNodes<F> {
    pub lam: Vec<u8>,
    pub n: usize,
    pub seed: u64,
    /// scale factor applied to the drawn nodes (1 for the base sequence; λ for a scaling scan)
    pub scale: F,
    pub a: Vec<Vec<F>>,
}

impl<F: PrimeField> BoxNodes<F> {
    pub fn new(lam: &[u8], n: usize, seed: u64) -> Self {
        assert!(lam.iter().all(|&l| (l as usize) <= n), "BoxNodes: a label of {:?} exceeds N = {}", lam, n);
        let p = F::characteristic();
        assert!(p > (8 * n * lam.len().max(1)) as u64, "BoxNodes: field too small");
        let mut rng = StdRng::seed_from_u64(seed ^ 0x6d62_6f78_6e6f_6465);
        let mut seen = std::collections::HashSet::new();
        let a = lam
            .iter()
            .map(|_| {
                (0..n)
                    .map(|_| loop {
                        let x = rng.gen_range(1..p);
                        if seen.insert(x) {
                            break F::from_u64(x);
                        }
                    })
                    .collect()
            })
            .collect();
        BoxNodes { lam: lam.to_vec(), n, seed, scale: F::one(), a }
    }
    /// The node sequence λ·a (its points are λ·z_μ). λ must be nonzero.
    pub fn scaled(&self, l: &F) -> Self {
        assert!(!l.is_zero());
        BoxNodes { lam: self.lam.clone(), n: self.n, seed: self.seed, scale: self.scale.mul(l), a: self.a.iter().map(|v| v.iter().map(|x| x.mul(l)).collect()).collect() }
    }
    /// Hash of the node values (and prime, λ, N): the checkpoint's "node-seed hash".
    pub fn hash(&self) -> u64 {
        let mut h = Fnv::default();
        h.str("nodes").u64(F::characteristic()).u64(self.n as u64).u64(self.seed);
        for (l, v) in self.lam.iter().zip(self.a.iter()) {
            h.u64(*l as u64);
            for x in v {
                h.u64(x.to_u64());
            }
        }
        h.get()
    }
}

/// Partitions in the a × (n − a) box (non-increasing, length a), ordered by size, then
/// lexicographically. Index 0 is ∅.
pub fn box_parts(a: usize, n: usize) -> Vec<Vec<u8>> {
    let s = n.saturating_sub(a);
    let mut out = vec![];
    fn rec(i: usize, maxv: usize, cur: &mut Vec<u8>, out: &mut Vec<Vec<u8>>) {
        if i == cur.len() {
            out.push(cur.clone());
            return;
        }
        for v in 0..=maxv {
            cur[i] = v as u8;
            rec(i + 1, v, cur, out);
        }
        cur[i] = 0;
    }
    rec(0, s, &mut vec![0; a], &mut out);
    out.sort_by_key(|p| (p.iter().map(|&x| x as u32).sum::<u32>(), p.clone()));
    out
}

fn contained(r: &[u8], m: &[u8]) -> bool {
    r.iter().zip(m.iter()).all(|(x, y)| x <= y)
}

/// Values of e_1..e_a of block i at the point with partition mu (roots a_{mu_j + a + 1 − j}).
fn block_e_values<F: Field>(nodes: &[F], mu: &[u8]) -> Vec<F> {
    let a = mu.len();
    let roots: Vec<&F> = (0..a).map(|j| &nodes[mu[j] as usize + a - 1 - j]).collect();
    // elementary symmetric functions
    let mut e = vec![F::zero(); a + 1];
    e[0] = F::one();
    for r in roots {
        for k in (1..=a).rev() {
            let t = e[k - 1].mul(r);
            e[k].add_assign(&t);
        }
    }
    e[1..].to_vec()
}

/// Values of the R^λ variables e_{i,j} (ERing variable order) at z_μ, μ given per block.
pub fn point_coords<F: PrimeField>(nodes: &BoxNodes<F>, mu: &[Vec<u8>]) -> Vec<F> {
    assert_eq!(mu.len(), nodes.lam.len());
    let mut out = vec![];
    for (i, m) in mu.iter().enumerate() {
        assert_eq!(m.len(), nodes.lam[i] as usize, "point_coords: block {} needs a partition of length {}", i, nodes.lam[i]);
        assert!(m.first().map_or(true, |&x| x as usize + m.len() <= nodes.n), "point_coords: block {} partition {:?} outside the box", i, m);
        out.extend(block_e_values(&nodes.a[i], m));
    }
    out
}

/// Value of p ∈ R^λ at the coordinates `x` (ERing variable order).
pub fn eval_epoly<F: Field>(ring: &ERing, p: &EPoly<F>, x: &[F]) -> F {
    let mut s = F::zero();
    for (m, c) in &p.t {
        let mut t = c.clone();
        for v in 0..ring.nvars {
            let e = crate::cclose::emono_exp(*m, v);
            for _ in 0..e {
                t.mul_assign(&x[v]);
            }
        }
        s.add_assign(&t);
    }
    s
}

/// Factorial Schur function s_ρ(x | a) of one block at distinct roots x: det[(x_j | a)^{ρ_l + a − l}] /
/// det[(x_j | a)^{a − l}] with (x | a)^m = Π_{t ≤ m} (x − a_t).
fn factorial_schur<F: Field>(nodes: &[F], rho: &[u8], x: &[F]) -> F {
    let a = rho.len();
    let fpow = |xv: &F, m: usize| -> F {
        let mut r = F::one();
        for t in 0..m {
            r.mul_assign(&xv.sub(&nodes[t]));
        }
        r
    };
    let num: Vec<Vec<F>> = (0..a).map(|l| (0..a).map(|j| fpow(&x[j], rho[l] as usize + a - 1 - l)).collect()).collect();
    let den: Vec<Vec<F>> = (0..a).map(|l| (0..a).map(|j| fpow(&x[j], a - 1 - l)).collect()).collect();
    det(num).mul(&det(den).inv())
}

fn det<F: Field>(mut m: Vec<Vec<F>>) -> F {
    let n = m.len();
    let mut d = F::one();
    for k in 0..n {
        let Some(p) = (k..n).find(|&i| !m[i][k].is_zero()) else { return F::zero() };
        if p != k {
            m.swap(p, k);
            d = d.neg();
        }
        d.mul_assign(&m[k][k]);
        let inv = m[k][k].inv();
        for i in (k + 1)..n {
            if m[i][k].is_zero() {
                continue;
            }
            let f = m[i][k].mul(&inv).neg();
            for j in k..n {
                let t = m[k][j].mul(&f);
                m[i][j].add_assign(&t);
            }
        }
    }
    d
}

/// Dual Jacobi-Trudi: s_ρ = det[e_{ρ'_r − r + c}] in the variables of block i (offset `voff`).
fn jacobi_trudi<F: Field>(ring: &ERing, i: usize, rho: &[u8]) -> EPoly<F> {
    let a = rho.len();
    let l = rho.first().copied().unwrap_or(0) as usize;
    let conj: Vec<usize> = (0..l).map(|r| rho.iter().filter(|&&x| x as usize > r).count()).collect();
    let e = |k: i64| -> EPoly<F> {
        if k == 0 {
            EPoly::one()
        } else if k < 0 || k as usize > a {
            EPoly::zero()
        } else {
            let mut ex = vec![0u32; ring.nvars];
            ex[ring.var(i, k as usize)] = 1;
            EPoly::mono(ring.pack(&ex), F::one())
        }
    };
    let m: Vec<Vec<EPoly<F>>> = (0..l).map(|r| (0..l).map(|c| e(conj[r] as i64 - r as i64 + c as i64)).collect()).collect();
    fn pdet<F: Field>(m: &[Vec<EPoly<F>>], cols: &mut Vec<usize>, row: usize) -> EPoly<F> {
        if row == m.len() {
            return EPoly::one();
        }
        let mut s = EPoly::zero();
        for k in 0..cols.len() {
            let c = cols[k];
            if m[row][c].is_zero() {
                continue;
            }
            cols.remove(k);
            let sub = pdet(m, cols, row + 1);
            cols.insert(k, c);
            let t = m[row][c].mul(&sub);
            if k % 2 == 0 {
                s.add_assign(&t)
            } else {
                s = s.sub(&t)
            }
        }
        s
    }
    pdet(&m, &mut (0..l).collect(), 0)
}

/// The box at fixed N for λ and a mode: per block the partitions (only ∅ for the marked block),
/// the factorial Schur matrices S^(i) and their inverses, the point numbering and the
/// Jacobi-Trudi polynomials of the output basis.
pub struct BoxGrid<F> {
    pub nodes: BoxNodes<F>,
    pub lam: Vec<u8>,
    pub n: usize,
    pub marked: Option<usize>,
    pub ring: ERing,
    /// per block: partitions (index = per-block point index)
    pub parts: Vec<Vec<Vec<u8>>>,
    psize: Vec<Vec<u32>>,
    strides: Vec<u32>,
    pub npoints: u32,
    /// S^(i)[μ][ρ] = s_ρ(z_μ | a^(i)) (dense, per block)
    pub smat: Vec<Vec<Vec<F>>>,
    /// (S^(i))^{-1}[ρ][μ] (dense, per block)
    pub sinv: Vec<Vec<Vec<F>>>,
    /// per block, per μ: (ρ ⊇ μ, (S^(i))^{-1}[ρ][μ]) with nonzero coefficient
    sup: Vec<Vec<Vec<(u32, F)>>>,
    /// size |μ| of every point
    size_of: Vec<u32>,
    /// points of each size, and the position of a point within its size class
    by_size: Vec<Vec<u32>>,
    pos_in_size: Vec<u32>,
    /// largest size (top weight of the mode)
    pub top: u32,
    jt: Vec<Vec<EPoly<F>>>,
}

impl<F: PrimeField> BoxGrid<F> {
    /// The box of λ at N with the given nodes; `marked = Some(m)`: reduced mode (μ^m = ∅ only).
    pub fn new(nodes: BoxNodes<F>, marked: Option<usize>) -> Self {
        let lam = nodes.lam.clone();
        let n = nodes.n;
        if let Some(m) = marked {
            assert!(m < lam.len() && lam[m] > 0, "BoxGrid: invalid marked block {}", m);
        }
        let ring = ERing::new(&lam);
        let mut parts = vec![];
        let mut smat = vec![];
        let mut sinv = vec![];
        let mut sup = vec![];
        let mut jt = vec![];
        for (i, &l) in lam.iter().enumerate() {
            let a = l as usize;
            let ps = if Some(i) == marked { vec![vec![0u8; a]] } else { box_parts(a, n) };
            let k = ps.len();
            // S[μ][ρ] = s_ρ(z_μ | a)
            let s: Vec<Vec<F>> = ps
                .iter()
                .map(|mu| {
                    let roots: Vec<F> = (0..a).map(|j| nodes.a[i][mu[j] as usize + a - 1 - j].clone()).collect();
                    ps.iter().map(|rho| factorial_schur(&nodes.a[i], rho, &roots)).collect()
                })
                .collect();
            for (x, mu) in ps.iter().enumerate() {
                for (y, rho) in ps.iter().enumerate() {
                    let c = contained(rho, mu);
                    assert!(c || s[x][y].is_zero(), "factorial Schur vanishing fails: s_{:?}(z_{:?}) != 0", rho, mu);
                    if x == y {
                        assert!(!s[x][y].is_zero(), "factorial Schur: s_{:?}(z_{:?}) = 0 (degenerate nodes)", rho, mu);
                    }
                }
            }
            // inverse of the lower-triangular S (order by size is a linear extension of ⊆)
            let mut t = vec![vec![F::zero(); k]; k];
            for c in 0..k {
                // solve S t_c = e_c (column c of the inverse), forward substitution
                for r in 0..k {
                    let mut v = if r == c { F::one() } else { F::zero() };
                    for j in 0..r {
                        if !s[r][j].is_zero() && !t[j][c].is_zero() {
                            let q = s[r][j].mul(&t[j][c]);
                            v.sub_assign(&q);
                        }
                    }
                    t[r][c] = v.mul(&s[r][r].inv());
                }
            }
            let mut sp: Vec<Vec<(u32, F)>> = vec![vec![]; k];
            for (r, rho) in ps.iter().enumerate() {
                for (m, mu) in ps.iter().enumerate() {
                    if !t[r][m].is_zero() {
                        assert!(contained(mu, rho), "inverse factorial Schur matrix not triangular for containment");
                        sp[m].push((r as u32, t[r][m].clone()));
                    }
                }
            }
            jt.push(ps.iter().map(|rho| jacobi_trudi::<F>(&ring, i, rho)).collect());
            parts.push(ps);
            smat.push(s);
            sinv.push(t);
            sup.push(sp);
        }
        let psize: Vec<Vec<u32>> = parts.iter().map(|ps| ps.iter().map(|p| p.iter().map(|&x| x as u32).sum()).collect()).collect();
        let mut strides = vec![];
        let mut np: u64 = 1;
        for ps in &parts {
            strides.push(np as u32);
            np *= ps.len() as u64;
        }
        assert!(np < u32::MAX as u64, "BoxGrid: too many points");
        let npoints = np as u32;
        let top: u32 = psize.iter().map(|v| *v.iter().max().unwrap_or(&0)).sum();
        let mut size_of = vec![0u32; npoints as usize];
        let mut by_size = vec![vec![]; top as usize + 1];
        let mut pos_in_size = vec![0u32; npoints as usize];
        for id in 0..npoints {
            let s: u32 = (0..parts.len()).map(|i| psize[i][((id / strides[i]) % parts[i].len() as u32) as usize]).sum();
            size_of[id as usize] = s;
            pos_in_size[id as usize] = by_size[s as usize].len() as u32;
            by_size[s as usize].push(id);
        }
        BoxGrid { nodes, lam, n, marked, ring, parts, psize, strides, npoints, smat, sinv, sup, size_of, by_size, pos_in_size, top, jt }
    }

    /// Per-block partition indices of a point.
    pub fn index(&self, id: u32) -> Vec<u32> {
        (0..self.parts.len()).map(|i| (id / self.strides[i]) % self.parts[i].len() as u32).collect()
    }
    /// The point with the given per-block partition indices.
    pub fn id_of(&self, idx: &[u32]) -> u32 {
        idx.iter().zip(self.strides.iter()).map(|(a, b)| a * b).sum()
    }
    /// μ of a point, per block.
    pub fn mu(&self, id: u32) -> Vec<Vec<u8>> {
        self.index(id).iter().enumerate().map(|(i, &k)| self.parts[i][k as usize].clone()).collect()
    }
    /// The point of a multi-partition (None if outside the box or the mode).
    pub fn id_of_mu(&self, mu: &[Vec<u8>]) -> Option<u32> {
        let mut idx = vec![];
        for (i, m) in mu.iter().enumerate() {
            idx.push(self.parts[i].iter().position(|p| p == m)? as u32);
        }
        Some(self.id_of(&idx))
    }
    pub fn size(&self, id: u32) -> u32 {
        self.size_of[id as usize]
    }
    /// Number of box elements (ρ of the mode) of size s.
    pub fn count_of_size(&self, s: u32) -> usize {
        self.by_size.get(s as usize).map_or(0, |v| v.len())
    }
    /// The ρ of size s (point ids), in output order.
    pub fn of_size(&self, s: u32) -> &[u32] {
        self.by_size.get(s as usize).map(|v| v.as_slice()).unwrap_or(&[])
    }
    /// e-coordinates of the point (ERing order).
    pub fn coords(&self, id: u32) -> Vec<F> {
        point_coords(&self.nodes, &self.mu(id))
    }
    /// ρ ⊇ μ of size s with coefficient (S^{-1})_{ρμ}: (position of ρ in its size class, coefficient).
    pub fn supersets_of_size(&self, mu: u32, s: u32, out: &mut Vec<(u32, F)>) {
        out.clear();
        if s > self.top {
            return;
        }
        let idx = self.index(mu);
        let k = self.parts.len();
        // max remaining size after block i
        let mut rem = vec![0u32; k + 1];
        for i in (0..k).rev() {
            rem[i] = rem[i + 1] + *self.psize[i].iter().max().unwrap_or(&0);
        }
        fn rec<F: PrimeField>(g: &BoxGrid<F>, idx: &[u32], rem: &[u32], i: usize, left: u32, id: u32, c: F, out: &mut Vec<(u32, F)>) {
            if i == idx.len() {
                if left == 0 {
                    out.push((g.pos_in_size[id as usize], c));
                }
                return;
            }
            if left > rem[i] {
                return;
            }
            for (r, x) in &g.sup[i][idx[i] as usize] {
                let sz = g.psize[i][*r as usize];
                if sz > left {
                    continue;
                }
                rec(g, idx, rem, i + 1, left - sz, id + r * g.strides[i], c.mul(x), out);
            }
        }
        rec(self, &idx, &rem, 0, s, 0, F::one(), out);
    }
    /// All points of the mode with |μ| ≤ max_size (the down-set L(max_size)), by size.
    pub fn points_up_to(&self, max_size: u32) -> Vec<u32> {
        (0..=max_size.min(self.top)).flat_map(|s| self.of_size(s).iter().copied()).collect()
    }
}

/// All points of the box of λ at N in the mode, as per-block partitions (|μ| ≤ max_size).
pub fn box_points<F: PrimeField>(grid: &BoxGrid<F>, max_size: u32) -> Vec<Vec<Vec<u8>>> {
    grid.points_up_to(max_size).into_iter().map(|id| grid.mu(id)).collect()
}

/// R^λ element Σ_ρ c_ρ s_ρ from Schur coordinates (ρ given as point ids of the grid),
/// restricted to N' (ρ^i ⊆ λ_i × (N' − λ_i)) and, with `marked = Some(m)`, to ρ^m = ∅. Its image in
/// A_λ at gl_{N'} (in that mode) is the entry mod I_{N'}. A label > N' gives 0.
pub fn schur_to_epoly<F: PrimeField>(grid: &BoxGrid<F>, coords: &[(u32, F)], n_prime: usize, marked: Option<usize>) -> EPoly<F> {
    assert!(n_prime <= grid.n, "schur_to_epoly: N' = {} above the box N = {}", n_prime, grid.n);
    if grid.lam.iter().any(|&l| l as usize > n_prime) {
        return EPoly::zero();
    }
    let mut out = EPoly::zero();
    for (id, c) in coords {
        let idx = grid.index(*id);
        let mut ok = true;
        for (i, &k) in idx.iter().enumerate() {
            let p = &grid.parts[i][k as usize];
            let a = grid.lam[i] as usize;
            if p.first().map_or(false, |&x| x as usize > n_prime - a) || (Some(i) == marked && p.iter().any(|&x| x > 0)) {
                ok = false;
                break;
            }
        }
        if !ok {
            continue;
        }
        let mut t = EPoly::constant(c.clone());
        for (i, &k) in idx.iter().enumerate() {
            t = t.mul(&grid.jt[i][k as usize]);
        }
        out.add_assign(&t);
    }
    out
}

// ------------------------------------------------------------------ streaming accumulators

/// Checkpoint header (SPEC D14, C_sound M2). Merging requires equality of every field.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BoxHeader {
    pub prime: u64,
    pub node_seed_hash: u64,
    pub library_hash: u64,
    pub word_hash: u64,
    pub prune_n: Option<usize>,
    pub marked: Option<usize>,
    pub n: usize,
    pub lam: Vec<u8>,
    pub code_stamp: String,
    pub entries_hash: u64,
    /// what the entries mean (M6a-2a, review wave 3 box F2 / sengine F4)
    pub bind: BoxBinding,
}

/// The binding of a checkpoint to the meaning of its entries (M6a-2a): two shards whose weight
/// lists agree but whose entries are coordinates in different bases must never be summed.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct BoxBinding {
    /// closable stage (crossing count) the entries belong to (0: unbound, tests)
    pub stage: usize,
    /// hash of the closure skeleton: entry kinds, objects, edges, type degrees, crossings
    pub skeleton_hash: u64,
    /// engine version and the options that decide the frozen structure (pivot order, streaming)
    pub engine_stamp: String,
    /// hash of the frozen trace recorded at the record point (pivot sequence, patterns, values)
    pub trace_hash: u64,
    /// hash of the entry values at the record point (the library maps' content, in this basis)
    pub record_hash: u64,
}

/// The fields in which two headers differ (for refusal messages).
pub fn header_diff(a: &BoxHeader, b: &BoxHeader) -> String {
    let mut d = vec![];
    macro_rules! cmp {
        ($($f:ident).+) => {
            if a.$($f).+ != b.$($f).+ {
                d.push(format!("{} {:?} vs {:?}", stringify!($($f).+), a.$($f).+, b.$($f).+));
            }
        };
    }
    cmp!(prime);
    cmp!(node_seed_hash);
    cmp!(library_hash);
    cmp!(word_hash);
    cmp!(prune_n);
    cmp!(marked);
    cmp!(n);
    cmp!(lam);
    cmp!(code_stamp);
    cmp!(entries_hash);
    cmp!(bind.stage);
    cmp!(bind.skeleton_hash);
    cmp!(bind.engine_stamp);
    cmp!(bind.trace_hash);
    cmp!(bind.record_hash);
    d.join("; ")
}

/// Hash of an entry weight list (None: an entry that must be 0 at every point).
pub fn entries_hash(weights: &[Option<u32>]) -> u64 {
    let mut h = Fnv::default();
    h.str("entries").u64(weights.len() as u64);
    for w in weights {
        h.i64(w.map_or(-1, |x| x as i64));
    }
    h.get()
}

/// Streaming box accumulators of a list of entries (declared weights), with the done-point
/// bitmap and the checkpoint header.
pub struct BoxAcc<F> {
    pub header: BoxHeader,
    pub grid: Arc<BoxGrid<F>>,
    weights: Vec<Option<u32>>,
    off: Vec<usize>,
    nout: Vec<usize>,
    by_weight: BTreeMap<u32, Vec<usize>>,
    acc: Vec<F>,
    done: Vec<u64>,
    ndone: usize,
}

/// Extracted Schur coordinates per entry ((ρ point id, c_ρ), |ρ| = declared weight).
pub struct BoxResult<F> {
    pub grid: Arc<BoxGrid<F>>,
    pub coords: Vec<Vec<(u32, F)>>,
}

impl<F: PrimeField> BoxAcc<F> {
    /// Empty accumulators. The header's prime, node hash, mode, N, λ, code stamp and entries
    /// hash are filled from the grid and weights; library hash, word hash and prune_n come from
    /// the caller.
    pub fn new(grid: Arc<BoxGrid<F>>, weights: Vec<Option<u32>>, library_hash: u64, word_hash: u64, prune_n: Option<usize>) -> Self {
        Self::with_binding(grid, weights, library_hash, word_hash, prune_n, BoxBinding::default())
    }
    /// As [`BoxAcc::new`], with the binding of the entries (stage, skeleton, engine, trace).
    pub fn with_binding(grid: Arc<BoxGrid<F>>, weights: Vec<Option<u32>>, library_hash: u64, word_hash: u64, prune_n: Option<usize>, bind: BoxBinding) -> Self {
        let header = BoxHeader {
            prime: F::characteristic(),
            node_seed_hash: grid.nodes.hash(),
            library_hash,
            word_hash,
            prune_n,
            marked: grid.marked,
            n: grid.n,
            lam: grid.lam.clone(),
            code_stamp: BOX_CODE_STAMP.to_string(),
            entries_hash: entries_hash(&weights),
            bind,
        };
        let mut off = vec![];
        let mut nout = vec![];
        let mut tot = 0usize;
        let mut by_weight: BTreeMap<u32, Vec<usize>> = BTreeMap::new();
        for (e, w) in weights.iter().enumerate() {
            off.push(tot);
            match w {
                Some(w) => {
                    let o = grid.count_of_size(*w);
                    let r = grid.count_of_size(w + 1);
                    nout.push(o);
                    tot += o + r;
                    by_weight.entry(*w).or_default().push(e);
                }
                None => nout.push(0),
            }
        }
        let words = (grid.npoints as usize + 63) / 64;
        BoxAcc { header, grid, weights, off, nout, by_weight, acc: vec![F::zero(); tot], done: vec![0; words], ndone: 0 }
    }
    pub fn weights(&self) -> &[Option<u32>] {
        &self.weights
    }
    pub fn is_done(&self, id: u32) -> bool {
        (self.done[(id / 64) as usize] >> (id % 64)) & 1 == 1
    }
    pub fn points_done(&self) -> usize {
        self.ndone
    }
    /// Number of accumulator coefficients held (outputs + residuals).
    pub fn acc_len(&self) -> usize {
        self.acc.len()
    }

    /// Add the values of every entry at point `id`. Refuses a point already received (a
    /// double add would corrupt every c_ρ with ρ ⊇ μ), a wrong value count, and a nonzero value
    /// of an entry declared zero.
    pub fn add_point(&mut self, id: u32, vals: &[F]) -> Result<(), String> {
        if id >= self.grid.npoints {
            return Err(format!("point {} outside the box ({} points)", id, self.grid.npoints));
        }
        if self.is_done(id) {
            return Err(format!("point {} ({:?}) was already added", id, self.grid.mu(id)));
        }
        if vals.len() != self.weights.len() {
            return Err(format!("point {}: {} values for {} entries", id, vals.len(), self.weights.len()));
        }
        for (e, w) in self.weights.iter().enumerate() {
            if w.is_none() && !vals[e].is_zero() {
                return Err(format!("point {}: entry {} is declared 0 (negative or odd degree) but has a nonzero value", id, e));
            }
        }
        let mut sup_w = vec![];
        let mut sup_r = vec![];
        let ms = self.grid.size(id);
        for (&w, ents) in &self.by_weight {
            if w + 1 < ms {
                continue; // μ is not below any output or residual ρ of this weight
            }
            if ents.iter().all(|&e| vals[e].is_zero()) {
                continue;
            }
            self.grid.supersets_of_size(id, w, &mut sup_w);
            self.grid.supersets_of_size(id, w + 1, &mut sup_r);
            for &e in ents {
                let v = &vals[e];
                if v.is_zero() {
                    continue;
                }
                let (o, no) = (self.off[e], self.nout[e]);
                for (k, c) in &sup_w {
                    self.acc[o + *k as usize].fma(c, v);
                }
                for (k, c) in &sup_r {
                    self.acc[o + no + *k as usize].fma(c, v);
                }
            }
        }
        self.done[(id / 64) as usize] |= 1u64 << (id % 64);
        self.ndone += 1;
        Ok(())
    }

    /// Largest size of a point any output or residual needs (the down-set L(m)).
    pub fn needed_size(&self) -> Option<u32> {
        self.by_weight.keys().filter(|&&w| w <= self.grid.top).map(|&w| (w + 1).min(self.grid.top)).max()
    }

    /// Schur coordinates of every entry. Err if a point of the down-set {μ ⊆ ρ} of some output
    /// or residual ρ is missing, or if a residual (|ρ| = w + 1) is nonzero.
    pub fn finish(&self) -> Result<BoxResult<F>, String> {
        if let Some(m) = self.needed_size() {
            let missing: Vec<u32> = self.grid.points_up_to(m).into_iter().filter(|&id| !self.is_done(id)).collect();
            if !missing.is_empty() {
                return Err(format!("incomplete down-set: {} of the points with |mu| <= {} are missing (first {:?})", missing.len(), m, self.grid.mu(missing[0])));
            }
        }
        let mut bad = vec![];
        let mut coords = vec![];
        for (e, w) in self.weights.iter().enumerate() {
            let Some(w) = w else {
                coords.push(vec![]);
                continue;
            };
            let (o, no) = (self.off[e], self.nout[e]);
            let nr = self.grid.count_of_size(w + 1);
            if self.acc[o + no..o + no + nr].iter().any(|x| !x.is_zero()) {
                bad.push(e);
            }
            let ids = self.grid.of_size(*w);
            coords.push((0..no).filter(|&k| !self.acc[o + k].is_zero()).map(|k| (ids[k], self.acc[o + k].clone())).collect());
        }
        if !bad.is_empty() {
            return Err(format!("nonzero residual at weight w + 1 for {} entries (first: entry {} declared at weight {:?}): not homogeneous of the declared weight", bad.len(), bad[0], self.weights[bad[0]]));
        }
        Ok(BoxResult { grid: self.grid.clone(), coords })
    }

    /// Sum another shard into this one. Refuses a mismatched header (every field) and
    /// overlapping done-point bitmaps.
    pub fn merge(&mut self, o: &BoxAcc<F>) -> Result<(), String> {
        if o.header != self.header {
            return Err(format!("mismatched checkpoint header: {}", header_diff(&self.header, &o.header)));
        }
        if o.weights != self.weights || o.acc.len() != self.acc.len() || o.done.len() != self.done.len() {
            return Err("mismatched entry lists".to_string());
        }
        let overlap: u32 = self.done.iter().zip(o.done.iter()).map(|(a, b)| (a & b).count_ones()).sum();
        if overlap > 0 {
            return Err(format!("overlapping shards: {} points in both", overlap));
        }
        for (a, b) in self.acc.iter_mut().zip(o.acc.iter()) {
            a.add_assign(b);
        }
        for (a, b) in self.done.iter_mut().zip(o.done.iter()) {
            *a |= *b;
        }
        self.ndone += o.ndone;
        Ok(())
    }

    /// Content digest: header (every field) and payload (bitmap, weights, accumulators). Two
    /// accumulators with equal digests hold byte-identical checkpoints.
    pub fn digest(&self) -> u64 {
        let mut h = Fnv::default();
        h.str(&format!("{:?}", self.header)).bytes(&self.payload());
        h.get()
    }

    fn payload(&self) -> Vec<u8> {
        let mut v = Vec::with_capacity(8 * (self.done.len() + self.weights.len() + self.acc.len()));
        for w in &self.done {
            v.extend_from_slice(&w.to_le_bytes());
        }
        for w in &self.weights {
            v.extend_from_slice(&w.map_or(-1i64, |x| x as i64).to_le_bytes());
        }
        for x in &self.acc {
            v.extend_from_slice(&x.to_u64().to_le_bytes());
        }
        v
    }

    /// Write a checkpoint atomically: a temporary file in the same directory, fsync, rename.
    pub fn save(&self, path: &Path) -> Result<(), String> {
        let payload = self.payload();
        let h = &self.header;
        let opt = |x: Option<usize>| x.map_or("none".to_string(), |v| v.to_string());
        let head = format!(
            "krhom-boxacc v2\nprime={}\nnode_seed_hash={}\nlibrary_hash={}\nword_hash={}\nprune_n={}\nmarked={}\nn={}\nlam={}\ncode_stamp={}\nentries_hash={}\nstage={}\nskeleton_hash={}\nengine_stamp={}\ntrace_hash={}\nrecord_hash={}\nentries={}\nnpoints={}\nndone={}\npayload_len={}\npayload_hash={}\nend\n",
            h.prime,
            h.node_seed_hash,
            h.library_hash,
            h.word_hash,
            opt(h.prune_n),
            opt(h.marked),
            h.n,
            h.lam.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","),
            h.code_stamp,
            h.entries_hash,
            h.bind.stage,
            h.bind.skeleton_hash,
            h.bind.engine_stamp,
            h.bind.trace_hash,
            h.bind.record_hash,
            self.weights.len(),
            self.grid.npoints,
            self.ndone,
            payload.len(),
            Fnv::default().bytes(&payload).get()
        );
        let tmp = path.with_extension(format!("tmp{}", std::process::id()));
        let res = (|| -> std::io::Result<()> {
            let mut f = std::fs::File::create(&tmp)?;
            f.write_all(head.as_bytes())?;
            f.write_all(&payload)?;
            f.sync_all()?;
            drop(f);
            std::fs::rename(&tmp, path)?;
            if let Some(dir) = path.parent() {
                if let Ok(d) = std::fs::File::open(if dir.as_os_str().is_empty() { Path::new(".") } else { dir }) {
                    let _ = d.sync_all();
                }
            }
            Ok(())
        })();
        res.map_err(|e| {
            let _ = std::fs::remove_file(&tmp);
            format!("checkpoint {}: {}", path.display(), e)
        })
    }

    /// [`BoxAcc::load`], then require the whole header to equal `want` (library, word, prune_n,
    /// stage, skeleton, engine, trace and record bindings included): a resume or a merge from a
    /// stale file of another run with the same weight list is refused (review wave 3 box F2).
    pub fn load_expect(path: &Path, grid: Arc<BoxGrid<F>>, want: &BoxHeader) -> Result<Self, String> {
        let acc = Self::load(path, grid)?;
        if acc.header != *want {
            return Err(format!("checkpoint {}: header does not match this run: {}", path.display(), header_diff(&acc.header, want)));
        }
        Ok(acc)
    }

    /// Read a checkpoint written by [`BoxAcc::save`] against `grid`: Err if the file is
    /// truncated or corrupt (payload hash), or its N, λ, mode, prime or node hash differ from the
    /// grid's.
    pub fn load(path: &Path, grid: Arc<BoxGrid<F>>) -> Result<Self, String> {
        let mut raw = vec![];
        std::fs::File::open(path).and_then(|mut f| f.read_to_end(&mut raw)).map_err(|e| format!("checkpoint {}: {}", path.display(), e))?;
        let end = raw.windows(5).position(|w| w == b"\nend\n").ok_or("checkpoint: no header end")? + 5;
        let head = std::str::from_utf8(&raw[..end]).map_err(|_| "checkpoint: header not UTF-8")?;
        let mut kv: BTreeMap<&str, &str> = BTreeMap::new();
        let mut lines = head.lines();
        if lines.next() != Some("krhom-boxacc v2") {
            return Err("checkpoint: not a krhom-boxacc v2 file".into());
        }
        for l in lines {
            if l == "end" {
                break;
            }
            let (k, v) = l.split_once('=').ok_or_else(|| format!("checkpoint: bad header line {:?}", l))?;
            kv.insert(k, v);
        }
        let get = |k: &str| kv.get(k).copied().ok_or_else(|| format!("checkpoint: missing {}", k));
        let num = |k: &str| -> Result<u64, String> { get(k)?.parse::<u64>().map_err(|_| format!("checkpoint: bad {}", k)) };
        let opt = |k: &str| -> Result<Option<usize>, String> {
            let v = get(k)?;
            if v == "none" { Ok(None) } else { v.parse::<usize>().map(Some).map_err(|_| format!("checkpoint: bad {}", k)) }
        };
        let payload = &raw[end..];
        if payload.len() as u64 != num("payload_len")? || Fnv::default().bytes(payload).get() != num("payload_hash")? {
            return Err("checkpoint: payload truncated or corrupt (length or hash mismatch)".into());
        }
        let lam: Vec<u8> = get("lam")?.split(',').filter(|s| !s.is_empty()).map(|s| s.parse::<u8>().map_err(|_| "checkpoint: bad lam".to_string())).collect::<Result<_, _>>()?;
        let header = BoxHeader {
            prime: num("prime")?,
            node_seed_hash: num("node_seed_hash")?,
            library_hash: num("library_hash")?,
            word_hash: num("word_hash")?,
            prune_n: opt("prune_n")?,
            marked: opt("marked")?,
            n: num("n")? as usize,
            lam,
            code_stamp: get("code_stamp")?.to_string(),
            entries_hash: num("entries_hash")?,
            bind: BoxBinding { stage: num("stage")? as usize, skeleton_hash: num("skeleton_hash")?, engine_stamp: get("engine_stamp")?.to_string(), trace_hash: num("trace_hash")?, record_hash: num("record_hash")? },
        };
        if header.prime != F::characteristic() || header.node_seed_hash != grid.nodes.hash() || header.n != grid.n || header.lam != grid.lam || header.marked != grid.marked || num("npoints")? != grid.npoints as u64 {
            return Err("checkpoint: prime, nodes, N, lambda or mode differ from the grid".into());
        }
        if header.code_stamp != BOX_CODE_STAMP {
            return Err(format!("checkpoint: code stamp {} != {}", header.code_stamp, BOX_CODE_STAMP));
        }
        let ne = num("entries")? as usize;
        let words = (grid.npoints as usize + 63) / 64;
        let rd = |i: usize| u64::from_le_bytes(payload[8 * i..8 * i + 8].try_into().unwrap());
        if payload.len() < 8 * (words + ne) {
            return Err("checkpoint: payload too short".into());
        }
        let weights: Vec<Option<u32>> = (0..ne).map(|e| { let x = rd(words + e) as i64; if x < 0 { None } else { Some(x as u32) } }).collect();
        if entries_hash(&weights) != header.entries_hash {
            return Err("checkpoint: entries hash mismatch".into());
        }
        let mut acc = BoxAcc::with_binding(grid, weights, header.library_hash, header.word_hash, header.prune_n, header.bind.clone());
        if acc.header != header {
            return Err("checkpoint: header inconsistent with the grid".into());
        }
        if payload.len() != 8 * (words + ne + acc.acc.len()) {
            return Err("checkpoint: payload size does not match the entry list".into());
        }
        for i in 0..words {
            acc.done[i] = rd(i);
        }
        for k in 0..acc.acc.len() {
            acc.acc[k] = F::from_u64(rd(words + ne + k));
        }
        acc.ndone = acc.done.iter().map(|w| w.count_ones() as usize).sum();
        if acc.ndone as u64 != num("ndone")? {
            return Err("checkpoint: done count mismatch".into());
        }
        Ok(acc)
    }
}

// ------------------------------------------------------------------ degree scaling (C_sound M1)

/// Scaling law of a homogeneous value of weight w: v(λz) = λ^w v(z). Ok(true) if checked,
/// Ok(false) if the value is 0 at both points (no information), Err on a violation.
pub fn check_scaling<F: Field>(v_z: &F, v_lz: &F, lambda: &F, w: u32) -> Result<bool, String> {
    let mut lw = F::one();
    for _ in 0..w {
        lw.mul_assign(lambda);
    }
    if *v_lz != lw.mul(v_z) {
        return Err(format!("scaling law violated: v(lz) = {} != l^{} v(z) = {}", v_lz, w, lw.mul(v_z)));
    }
    Ok(!(v_z.is_zero() && v_lz.is_zero()))
}

/// Result of [`check_scaling_all`].
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ScalingReport {
    /// entries checked with a nonzero value
    pub checked: usize,
    /// entries with a declared weight that are 0 at both points (reported, not certified)
    pub zero_both: usize,
    /// entries declared 0 (negative or odd weight), verified 0 at both points (sengine review F1:
    /// reported apart from `zero_both`, which carries no information)
    pub declared_zero: usize,
}

/// [`check_scaling`] for every entry (declared weights; None: both values must be 0).
pub fn check_scaling_all<F: Field>(vz: &[F], vlz: &[F], weights: &[Option<u32>], lambda: &F) -> Result<ScalingReport, String> {
    let mut r = ScalingReport::default();
    for (e, w) in weights.iter().enumerate() {
        let ok = match w {
            Some(w) => check_scaling(&vz[e], &vlz[e], lambda, *w).map_err(|m| format!("entry {} (declared weight {}): {}", e, w, m))?,
            None => {
                if !vz[e].is_zero() || !vlz[e].is_zero() {
                    return Err(format!("entry {}: declared 0 but nonzero", e));
                }
                r.declared_zero += 1;
                continue;
            }
        };
        if ok {
            r.checked += 1;
        } else {
            r.zero_both += 1;
        }
    }
    Ok(r)
}

// ------------------------------------------------------------------ closure entries

/// One extracted closure entry: a Gram entry G[b][a] of a type, or a pairing entry P[b][a] of
/// an edge.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EntryKind {
    Gram { ty: u32, b: usize, a: usize },
    Pair { edge: usize, b: usize, a: usize },
}

/// The entry list of a closure skeleton: every Gram entry of every type, then every pairing
/// entry of every edge (s, t), with the declared weights from the bookkeeping (2w = ins + outs,
/// resp. q_s − q_t + ins_s + outs_t; None if negative or odd).
pub fn closure_entries<F: Field>(data: &ClosureDataC<F>, edges: &[(usize, usize)]) -> (Vec<EntryKind>, Vec<Option<u32>>) {
    let mut kinds = vec![];
    let mut ws = vec![];
    let wt = |d: i32| if d >= 0 && d % 2 == 0 { Some((d / 2) as u32) } else { None };
    for (&x, tb) in &data.types {
        for (b, &ob) in tb.outs.iter().enumerate() {
            for (a, &ia) in tb.ins.iter().enumerate() {
                kinds.push(EntryKind::Gram { ty: x, b, a });
                ws.push(wt(ia + ob));
            }
        }
    }
    for (ei, &(s, t)) in edges.iter().enumerate() {
        let (os, ot) = (data.objs[s], data.objs[t]);
        let (ts, tt) = (&data.types[&os.0], &data.types[&ot.0]);
        for (b, &ob) in tt.outs.iter().enumerate() {
            for (a, &ia) in ts.ins.iter().enumerate() {
                kinds.push(EntryKind::Pair { edge: ei, b, a });
                ws.push(wt(os.2 - ot.2 + ia + ob));
            }
        }
    }
    (kinds, ws)
}

/// Values of the entries at the coordinates `x` from exact R^λ data (the stand-in for a scan at a
/// box point: tests and the 2-strand regression).
pub fn eval_entries<F: Field>(data: &ClosureDataC<F>, pairing: &[(usize, usize, Vec<Vec<EPoly<F>>>)], kinds: &[EntryKind], x: &[F]) -> Vec<F> {
    let ring = data.ring();
    kinds
        .iter()
        .map(|k| match *k {
            EntryKind::Gram { ty, b, a } => eval_epoly(&ring, &data.types[&ty].gram[b][a], x),
            EntryKind::Pair { edge, b, a } => eval_epoly(&ring, &pairing[edge].2[b][a], x),
        })
        .collect()
}

// ------------------------------------------------------------------ box closure data (D5)

/// Stamps of box closure data.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BoxStamp {
    pub prime: u64,
    pub node_seed_hash: u64,
    pub library_hash: u64,
    pub word_hash: u64,
    pub code_stamp: String,
}

/// Freivalds certificate of D² = 0 on the closed complex at (N', mode), bound to the pairing.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct D2Cert {
    pub n: usize,
    pub marked: Option<usize>,
    /// random vectors per (h, q) block (false pass ≤ p^-vectors per block)
    pub vectors: usize,
    pub blocks: usize,
    pub seed: u64,
    pub pairing_hash: u64,
    pub library_hash: u64,
}

/// T0 record: uncapped graded dims of every type at (N', mode) equal to MOY, on this library.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct T0Record {
    pub n: usize,
    pub marked: Option<usize>,
    pub library_hash: u64,
    pub ranks: TypeRanks,
}

/// Box closure data (SPEC D5): the closure skeleton with Gram data from the box (`base.edges`
/// empty), the closed pairing P, and the validity guards.
#[derive(Clone, Debug)]
pub struct BoxClosureData<F> {
    pub base: ClosureDataC<F>,
    /// N of a pruned scan (None: unpruned)
    pub prune_n: Option<usize>,
    pub n_box: usize,
    /// min(prune_n, n_box)
    pub valid_n: usize,
    /// extraction mode: None (unreduced, serves both modes) or Some(m) (only μ^m = ∅ points:
    /// serves only the reduced theory with block m marked)
    pub mode: Option<usize>,
    pub pairing: Pairing<F>,
    pub d2_cert: BTreeMap<(usize, Option<usize>), D2Cert>,
    pub stamp: BoxStamp,
}

impl<F: PrimeField> BoxClosureData<F> {
    /// Wrap a skeleton and a pairing; valid_n = min(prune_n, n_box).
    pub fn new(mut base: ClosureDataC<F>, pairing: Pairing<F>, prune_n: Option<usize>, n_box: usize, mode: Option<usize>, stamp: BoxStamp) -> Result<Self, String> {
        base.edges.clear();
        base.validate()?;
        validate_pairing(&base, &pairing)?;
        let valid_n = prune_n.map_or(n_box, |p| p.min(n_box));
        Ok(BoxClosureData { base, prune_n, n_box, valid_n, mode, pairing, d2_cert: BTreeMap::new(), stamp })
    }

    /// Assemble from a box extraction of the entries of [`closure_entries`] (`edges` as there):
    /// Gram and pairing entries become [`schur_to_epoly`] at the box N in the grid's mode.
    pub fn from_box(skeleton: &ClosureDataC<F>, edges: &[(usize, usize)], kinds: &[EntryKind], res: &BoxResult<F>, prune_n: Option<usize>, stamp: BoxStamp) -> Result<Self, String> {
        let g = &res.grid;
        let mut types: BTreeMap<u32, TypeBasisC<F>> = skeleton.types.iter().map(|(&x, tb)| (x, TypeBasisC { ins: tb.ins.clone(), outs: tb.outs.clone(), gram: vec![vec![EPoly::zero(); tb.ins.len()]; tb.outs.len()] })).collect();
        let mut pairing: Pairing<F> = edges
            .iter()
            .map(|&(s, t)| {
                let (ts, tt) = (&skeleton.types[&skeleton.objs[s].0], &skeleton.types[&skeleton.objs[t].0]);
                (s, t, vec![vec![EPoly::zero(); ts.ins.len()]; tt.outs.len()])
            })
            .collect();
        if kinds.len() != res.coords.len() {
            return Err("from_box: entry list and result differ in length".into());
        }
        for (k, c) in kinds.iter().zip(res.coords.iter()) {
            let p = schur_to_epoly(g, c, g.n, g.marked);
            match *k {
                EntryKind::Gram { ty, b, a } => types.get_mut(&ty).ok_or("from_box: unknown type")?.gram[b][a] = p,
                EntryKind::Pair { edge, b, a } => pairing[edge].2[b][a] = p,
            }
        }
        let base = ClosureDataC { lam: skeleton.lam.clone(), objs: skeleton.objs.clone(), types, edges: vec![], crossings: skeleton.crossings.clone() };
        let mut stamp = stamp;
        stamp.library_hash = library_hash(&base.types);
        BoxClosureData::new(base, pairing, prune_n, g.n, g.marked, stamp)
    }

    fn check_request(&self, n_prime: usize, marked: Option<usize>) -> Result<(), String> {
        if n_prime > self.valid_n {
            return Err(format!("N' = {} above valid_n = {} (prune_n {:?}, box N {})", n_prime, self.valid_n, self.prune_n, self.n_box));
        }
        if let Some(m) = self.mode {
            if marked != Some(m) {
                return Err(format!("box data extracted in reduced mode (block {} marked) cannot serve marked = {:?}", m, marked));
            }
        }
        Ok(())
    }
}

/// Freivalds test of D² = 0 at (N', mode) on box data (k vectors per block); on success the
/// certificate is recorded in `data.d2_cert` and returned. Refuses N' > valid_n.
pub fn closed_d2_freivalds<F: PrimeField>(data: &mut BoxClosureData<F>, n_prime: usize, marked: Option<usize>, k: usize, seed: u64) -> Result<D2Cert, String> {
    data.check_request(n_prime, marked)?;
    let opts = StagedOpts::default();
    let blocks = closed_d2_freivalds_raw(&data.base, &data.pairing, n_prime, marked, &opts, k, seed)?;
    let cert = D2Cert { n: n_prime, marked, vectors: k, blocks, seed, pairing_hash: pairing_hash(&data.pairing), library_hash: library_hash(&data.base.types) };
    data.d2_cert.insert((n_prime, marked), cert.clone());
    Ok(cert)
}

/// The colored closure of box data at gl_{N'} (SPEC D5, C_sound M5). Refuses N' > valid_n, a
/// mode the extraction cannot serve, and pruned data without a D² certificate for (N', mode)
/// that matches the current pairing and library. With `t0`, the MOY rank cap is applied, after
/// checking that the T0 record is for this (N', mode) and this library hash.
pub fn close_box<F: PrimeField>(data: &BoxClosureData<F>, n_prime: usize, marked: Option<usize>, t0: Option<&T0Record>, cap: Option<u128>) -> Result<Poincare, String> {
    data.check_request(n_prime, marked)?;
    let lh = library_hash(&data.base.types);
    if data.prune_n.is_some() {
        let c = data.d2_cert.get(&(n_prime, marked)).ok_or_else(|| format!("pruned box data (prune_n {:?}) without a D² certificate at N' = {}, marked {:?}: run closed_d2_freivalds first", data.prune_n, n_prime, marked))?;
        if c.pairing_hash != pairing_hash(&data.pairing) || c.library_hash != lh {
            return Err("the D² certificate does not match the current pairing or library".into());
        }
    }
    let mut opts = StagedOpts { cap, rank_cap: None };
    if let Some(t) = t0 {
        if t.n != n_prime || t.marked != marked || t.library_hash != lh {
            return Err(format!("T0 record (N {}, marked {:?}, library {:x}) does not match (N {}, marked {:?}, library {:x})", t.n, t.marked, t.library_hash, n_prime, marked, lh));
        }
        opts.rank_cap = Some(t.ranks.clone());
    }
    try_close_staged(&data.base, &data.pairing, n_prime, marked, &opts)
}

/// MOY value tr_N(b_p) of the type with table `table` (in the block (λ, λ)), in the convention of
/// `schur::decat_chi` (the q-graded dimension of Φ_N(B_p) with the global shift of a 0-crossing
/// closure).
pub fn moy_type_trace(lam: &[u8], table: &str, n: usize) -> LPoly {
    let mut s = crate::schur::Schur::new(lam);
    let p = s.intern(crate::schur::Table::parse(table));
    let res = crate::schur::ScanResult { stages: vec![], k0: vec![(p, 0, 1)], gram: vec![] };
    crate::schur::decat_chi(&mut s, &res, &[], n)
}

/// T0 (SPEC D14): the uncapped graded dimension of Φ_{N'}(B_p) of every type of the box library
/// equals its MOY value `moy[p]` (both as q-polynomials with the 0-crossing global shift).
/// Returns the record (ranks per type and degree, library hash) that enables the rank cap.
pub fn box_moy_check<F: PrimeField>(data: &BoxClosureData<F>, n_prime: usize, marked: Option<usize>, moy: &BTreeMap<u32, LPoly>, cap: Option<u128>) -> Result<T0Record, String> {
    data.check_request(n_prime, marked)?;
    let ranks = try_staged_type_ranks(&data.base, n_prime, marked, &StagedOpts { cap, rank_cap: None })?;
    let gq = global_q(&data.base.lam, &[], n_prime, marked);
    for (x, r) in &ranks {
        let mut got = LPoly::zero();
        for (d, k) in r {
            got.add_term(d + gq, *k as i128);
        }
        let want = moy.get(x).ok_or_else(|| format!("T0: no MOY value for type {}", x))?;
        if &got != want {
            return Err(format!("T0 fails for type {} at N' = {} (marked {:?}): closure {:?} != MOY {:?}", x, n_prime, marked, got, want));
        }
    }
    Ok(T0Record { n: n_prime, marked, library_hash: library_hash(&data.base.types), ranks })
}
