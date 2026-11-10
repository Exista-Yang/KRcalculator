//! src/schur.rs — Hecke algebroid (q-Schur algebroid) of type A in the contingency-table model.
//! Conventions: EKLP arXiv:2401.03053 §3 (eqs. (44), (45)), Williamson arXiv:1010.1283v3 §2.3,
//! R1_singular_theory.md §1 and Appendix B. Rows = bottom blocks (lambda), columns = top blocks (mu).
//!
//! This is the decategorified oracle (milestone M3a). It is a port of the validated Python
//! reference `scripts/schur.py` (table model, split/merge/rung rules, canonical basis, Hom
//! formula), `scripts/t8.py` (decategorified scan, Grothendieck lower bound) and
//! `proto_D1/peak.py` (peak model), extended by minimal-flow ladders, Gram rows and the
//! decategorified closure (Euler characteristic) through `moy::closure_value`.
//!
//! * An element of the block (lambda, mu) is written in the standard basis h_p, p a table with
//!   row sums lambda and column sums mu; coefficients are Laurent polynomials in v (v = q).
//! * split (finer columns):  h_q -> sum_{p refines q} v^{lmax(q) - lmax(p)} h_p          (EKLP (44))
//! * merge (coarser columns): h_p -> v^{lmin(q) - lmin(p)} pi(LR q)/pi(LR p) h_q           (EKLP (45))
//! * rung E_i^(r): r from upright i+1 to i; F_i^(r): r from upright i to i+1 (as `web.rs`).
//! * b_p is the canonical (singular KL) basis: bar-invariant, h_p + vZ[v]-span of lower terms.
//!
//! All coefficient arithmetic is i128 with checked operations (overflow panics).

use crate::laurent::LPoly;
use crate::web::{Dir, Ladder, Rung};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BinaryHeap, HashMap, HashSet, VecDeque};
use std::hash::{BuildHasherDefault, Hasher};
use std::rc::Rc;

pub type Comp = Vec<u8>; // composition, zeros allowed, K = sum <= 64

const OVF: &str = "schur: i128 coefficient overflow";

// ---------------------------------------------------------------------------------------------
// small fast hasher for integer keys (no external dependency)

#[derive(Default, Clone, Copy)]
struct FxHasher {
    h: u64,
}
const FX_SEED: u64 = 0x51_7c_c1_b7_27_22_0a_95;
impl FxHasher {
    #[inline]
    fn add(&mut self, w: u64) {
        self.h = (self.h.rotate_left(5) ^ w).wrapping_mul(FX_SEED);
    }
}
impl Hasher for FxHasher {
    /// murmur3 fmix64 finalizer: the multiply-only state has weak low bits, which hashbrown uses
    /// for the bucket index (keys like (tid << 32) | rung_code would collide)
    fn finish(&self) -> u64 {
        let mut h = self.h;
        h ^= h >> 33;
        h = h.wrapping_mul(0xff51_afd7_ed55_8ccd);
        h ^= h >> 33;
        h = h.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
        h ^ (h >> 33)
    }
    fn write(&mut self, bytes: &[u8]) {
        for ch in bytes.chunks(8) {
            let mut b = [0u8; 8];
            b[..ch.len()].copy_from_slice(ch);
            self.add(u64::from_le_bytes(b));
        }
    }
    fn write_u8(&mut self, i: u8) {
        self.add(i as u64)
    }
    fn write_u16(&mut self, i: u16) {
        self.add(i as u64)
    }
    fn write_u32(&mut self, i: u32) {
        self.add(i as u64)
    }
    fn write_u64(&mut self, i: u64) {
        self.add(i)
    }
    fn write_usize(&mut self, i: usize) {
        self.add(i as u64)
    }
}
type FxBuild = BuildHasherDefault<FxHasher>;
type FxMap<K, V> = HashMap<K, V, FxBuild>;
type FxSet<K> = HashSet<K, FxBuild>;

fn c2(n: i64) -> i64 {
    n * (n - 1) / 2
}
fn lw(comp: &[u8]) -> i64 {
    comp.iter().map(|&a| c2(a as i64)).sum()
}
fn binom_u128(n: u64, k: u64) -> u128 {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut r: u128 = 1;
    for i in 0..k {
        r = r.checked_mul((n - i) as u128).expect("binomial overflow") / (i + 1) as u128;
    }
    r
}

// ---------------------------------------------------------------------------------------------
// tables

#[derive(Clone, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub struct Table {
    pub nr: u8,
    pub nc: u8,
    pub e: Vec<u8>,
} // row-major

impl Table {
    /// diag(lambda)
    pub fn identity(lam: &[u8]) -> Table {
        let n = lam.len();
        let mut e = vec![0u8; n * n];
        for i in 0..n {
            e[i * n + i] = lam[i];
        }
        Table { nr: n as u8, nc: n as u8, e }
    }
    pub fn get(&self, i: usize, j: usize) -> u8 {
        self.e[i * self.nc as usize + j]
    }
    /// row sums = lambda
    pub fn rows(&self) -> Comp {
        let m = self.nc as usize;
        (0..self.nr as usize).map(|i| self.e[i * m..(i + 1) * m].iter().sum()).collect()
    }
    /// column sums = mu
    pub fn cols(&self) -> Comp {
        let m = self.nc as usize;
        (0..m).map(|j| (0..self.nr as usize).map(|i| self.e[i * m + j]).sum()).collect()
    }
    /// l(p_min) = sum_{i>i', j<j'} M_ij M_i'j'
    pub fn lmin(&self) -> i64 {
        let n = self.nr as usize;
        let m = self.nc as usize;
        let mut s = 0i64;
        for i in 0..n {
            for j in 0..m {
                let a = self.e[i * m + j] as i64;
                if a == 0 {
                    continue;
                }
                for i2 in 0..i {
                    for j2 in j + 1..m {
                        s += a * self.e[i2 * m + j2] as i64;
                    }
                }
            }
        }
        s
    }
    /// l(p_max) = lmin + l(w_lambda) + l(w_mu) - sum C(M_ij, 2)
    pub fn lmax(&self) -> i64 {
        self.lmin() + lw(&self.rows()) + lw(&self.cols()) - self.e.iter().map(|&a| c2(a as i64)).sum::<i64>()
    }
    /// N(p) = prod_i lambda_i! / prod_j M_ij!  (number of left W_mu-cosets in p)
    pub fn ncosets(&self) -> u128 {
        let m = self.nc as usize;
        let mut r: u128 = 1;
        for i in 0..self.nr as usize {
            let mut part: u64 = 0;
            for j in 0..m {
                let a = self.e[i * m + j] as u64;
                part += a;
                r = r.checked_mul(binom_u128(part, a)).expect("ncosets overflow");
            }
        }
        r
    }
    /// prod_ij [M_ij]! (symmetric quantum factorials)
    pub fn pi_lr(&self) -> LPoly {
        let mut r = LPoly::one();
        for &a in &self.e {
            r = r.mul(&LPoly::qfact(a as i64));
        }
        r
    }
    /// "2,0/0,2"
    pub fn parse(s: &str) -> Table {
        let rows: Vec<Vec<u8>> = s
            .trim()
            .split('/')
            .map(|r| r.split(',').map(|x| x.trim().parse::<u8>().unwrap_or_else(|_| panic!("bad table {:?}", s))).collect())
            .collect();
        let nc = rows[0].len();
        assert!(rows.iter().all(|r| r.len() == nc), "ragged table {:?}", s);
        Table { nr: rows.len() as u8, nc: nc as u8, e: rows.concat() }
    }
    pub fn fmt(&self) -> String {
        let m = self.nc as usize;
        (0..self.nr as usize)
            .map(|i| self.e[i * m..(i + 1) * m].iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","))
            .collect::<Vec<_>>()
            .join("/")
    }
}

/// All tables with row sums lam and column sums mu, sorted by (lmax, e).
pub fn tables(lam: &[u8], mu: &[u8]) -> Vec<Table> {
    let n = lam.len();
    let m = mu.len();
    let mut out = vec![];
    if lam.iter().map(|&a| a as u32).sum::<u32>() != mu.iter().map(|&a| a as u32).sum::<u32>() || n == 0 || m == 0 {
        return out;
    }
    let mut e = vec![0u8; n * m];
    let mut colleft = mu.to_vec();
    fn rec(i: usize, j: usize, left: u8, lam: &[u8], m: usize, e: &mut Vec<u8>, colleft: &mut Vec<u8>, out: &mut Vec<Table>) {
        let n = lam.len();
        if i == n {
            if colleft.iter().all(|&c| c == 0) {
                out.push(Table { nr: n as u8, nc: m as u8, e: e.clone() });
            }
            return;
        }
        if j == m - 1 {
            if left <= colleft[j] {
                e[i * m + j] = left;
                colleft[j] -= left;
                let nl = if i + 1 < n { lam[i + 1] } else { 0 };
                rec(i + 1, 0, nl, lam, m, e, colleft, out);
                colleft[j] += left;
                e[i * m + j] = 0;
            }
            return;
        }
        let hi = left.min(colleft[j]);
        for x in 0..=hi {
            e[i * m + j] = x;
            colleft[j] -= x;
            rec(i, j + 1, left - x, lam, m, e, colleft, out);
            colleft[j] += x;
        }
        e[i * m + j] = 0;
    }
    rec(0, 0, lam[0], lam, m, &mut e, &mut colleft, &mut out);
    out.sort_by_cached_key(|t| (t.lmax(), t.e.clone()));
    out
}

// ---------------------------------------------------------------------------------------------
// dense Laurent polynomials

/// Dense Laurent polynomial, i128 coefficients, checked arithmetic (unsimplified counts reach 1e11+).
/// Invariant: `c` is empty (zero, then lo = 0) or has nonzero first and last entries.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Default)]
pub struct Lp {
    pub lo: i32,
    pub c: Vec<i128>,
}

impl Lp {
    pub fn zero() -> Lp {
        Lp { lo: 0, c: Vec::new() }
    }
    pub fn one() -> Lp {
        Lp { lo: 0, c: vec![1] }
    }
    pub fn mono(e: i32, k: i128) -> Lp {
        if k == 0 { Lp::zero() } else { Lp { lo: e, c: vec![k] } }
    }
    pub fn is_zero(&self) -> bool {
        self.c.is_empty()
    }
    pub fn at1(&self) -> i128 {
        self.c.iter().fold(0i128, |a, &b| a.checked_add(b).expect(OVF))
    }
    /// v -> v^{-1}
    pub fn bar(&self) -> Lp {
        if self.is_zero() {
            return Lp::zero();
        }
        let mut c = self.c.clone();
        c.reverse();
        Lp { lo: -(self.lo + self.c.len() as i32 - 1), c }
    }
    pub fn shift(&self, s: i32) -> Lp {
        if self.is_zero() { Lp::zero() } else { Lp { lo: self.lo + s, c: self.c.clone() } }
    }
    pub fn mul(&self, o: &Lp) -> Lp {
        let mut r = Lp::zero();
        r.add_mul(self, o, 1);
        r.normalize();
        r
    }
    /// self += k * v^shift * o
    pub fn add_scaled(&mut self, o: &Lp, k: i128, shift: i32) {
        if k == 0 || o.is_zero() {
            return;
        }
        let lo = o.lo + shift;
        let hi = lo + o.c.len() as i32 - 1;
        self.reserve(lo, hi);
        let off = (lo - self.lo) as usize;
        for (i, &a) in o.c.iter().enumerate() {
            let t = a.checked_mul(k).expect(OVF);
            self.c[off + i] = self.c[off + i].checked_add(t).expect(OVF);
        }
        self.normalize();
    }
    /// The unique bar-invariant m with self - m in vZ[v] (schur.py sym_part).
    pub fn sym_part(&self) -> Lp {
        let mut m = Lp::zero();
        for (i, &a) in self.c.iter().enumerate() {
            if a == 0 {
                continue;
            }
            let k = self.lo + i as i32;
            if k < 0 {
                m.add_scaled(&Lp::one(), a, k);
                m.add_scaled(&Lp::one(), a, -k);
            } else if k == 0 {
                m.add_scaled(&Lp::one(), a, 0);
            }
        }
        m
    }
    pub fn max_exp(&self) -> Option<i32> {
        if self.is_zero() { None } else { Some(self.lo + self.c.len() as i32 - 1) }
    }
    pub fn min_exp(&self) -> Option<i32> {
        if self.is_zero() { None } else { Some(self.lo) }
    }
    pub fn to_lpoly(&self) -> LPoly {
        let mut p = LPoly::zero();
        for (i, &k) in self.c.iter().enumerate() {
            p.add_term(self.lo + i as i32, k);
        }
        p
    }
    pub fn from_lpoly(p: &LPoly) -> Lp {
        let mut r = Lp::zero();
        if let (Some(lo), Some(hi)) = (p.min_deg(), p.max_deg()) {
            r.reserve(lo, hi);
            for (&e, &k) in &p.c {
                r.c[(e - lo) as usize] = k;
            }
            r.normalize();
        }
        r
    }
    /// coefficient of v^e
    pub fn coef(&self, e: i32) -> i128 {
        if self.is_zero() || e < self.lo || e > self.lo + self.c.len() as i32 - 1 {
            0
        } else {
            self.c[(e - self.lo) as usize]
        }
    }
    pub fn neg(&self) -> Lp {
        Lp { lo: self.lo, c: self.c.iter().map(|&a| a.checked_neg().expect(OVF)).collect() }
    }
    pub fn add(&self, o: &Lp) -> Lp {
        let mut r = self.clone();
        r.add_scaled(o, 1, 0);
        r
    }
    pub fn sub(&self, o: &Lp) -> Lp {
        let mut r = self.clone();
        r.add_scaled(o, -1, 0);
        r
    }
    /// all coefficients >= 0
    pub fn is_nonneg(&self) -> bool {
        self.c.iter().all(|&a| a >= 0)
    }
    /// "e:c e:c ..." (nonzero terms, ascending exponents); "0" never appears, zero is the empty string.
    pub fn fmt_ec(&self) -> String {
        self.c
            .iter()
            .enumerate()
            .filter(|(_, k)| **k != 0)
            .map(|(i, k)| format!("{}:{}", self.lo + i as i32, k))
            .collect::<Vec<_>>()
            .join(" ")
    }
    /// inverse of `fmt_ec` (whitespace-separated "e:c" tokens; repeated exponents add up)
    pub fn parse_ec(s: &str) -> Lp {
        let mut r = Lp::zero();
        for tok in s.split_whitespace() {
            let (e, k) = tok.split_once(':').unwrap_or_else(|| panic!("bad term {:?}", tok));
            r.add_scaled(&Lp::one(), k.parse::<i128>().expect("bad coefficient"), e.parse::<i32>().expect("bad exponent"));
        }
        r
    }
    // --- internals: un-normalized accumulation
    fn hi(&self) -> i32 {
        self.lo + self.c.len() as i32 - 1
    }
    fn reserve(&mut self, lo: i32, hi: i32) {
        if self.c.is_empty() {
            self.lo = lo;
            self.c = vec![0; (hi - lo + 1) as usize];
            return;
        }
        let shi = self.hi();
        if lo < self.lo {
            let k = (self.lo - lo) as usize;
            let mut nc = vec![0i128; k];
            nc.extend_from_slice(&self.c);
            self.c = nc;
            self.lo = lo;
        }
        if hi > shi {
            self.c.resize(self.c.len() + (hi - shi) as usize, 0);
        }
    }
    fn normalize(&mut self) {
        let mut a = 0;
        while a < self.c.len() && self.c[a] == 0 {
            a += 1;
        }
        if a == self.c.len() {
            self.c.clear();
            self.lo = 0;
            return;
        }
        let mut b = self.c.len();
        while self.c[b - 1] == 0 {
            b -= 1;
        }
        if a > 0 || b < self.c.len() {
            self.c = self.c[a..b].to_vec();
            self.lo += a as i32;
        }
    }
    /// self += sign * a * b, without normalizing
    fn add_mul(&mut self, a: &Lp, b: &Lp, sign: i128) {
        if a.is_zero() || b.is_zero() {
            return;
        }
        let lo = a.lo + b.lo;
        let hi = a.hi() + b.hi();
        self.reserve(lo, hi);
        let off = (lo - self.lo) as usize;
        for (i, &x) in a.c.iter().enumerate() {
            if x == 0 {
                continue;
            }
            let xs = x.checked_mul(sign).expect(OVF);
            for (j, &y) in b.c.iter().enumerate() {
                if y == 0 {
                    continue;
                }
                let t = xs.checked_mul(y).expect(OVF);
                self.c[off + i + j] = self.c[off + i + j].checked_add(t).expect(OVF);
            }
        }
    }
}

// ---------------------------------------------------------------------------------------------
// elements

pub type Tid = u32; // interned table id (per Schur instance)

/// Element of the block (lambda, mu) in the standard basis h_p.
#[derive(Clone, Debug, PartialEq)]
pub struct Elt {
    pub mu: Comp,
    pub t: Vec<(Tid, Lp)>,
} // sorted by Tid, no zero coefficients

impl Elt {
    pub fn zero(mu: &[u8]) -> Elt {
        Elt { mu: mu.to_vec(), t: vec![] }
    }
    pub fn is_zero(&self) -> bool {
        self.t.is_empty()
    }
    /// coefficient of h_p (zero if absent)
    pub fn coef(&self, p: Tid) -> Lp {
        match self.t.binary_search_by_key(&p, |e| e.0) {
            Ok(k) => self.t[k].1.clone(),
            Err(_) => Lp::zero(),
        }
    }
    /// v^s * self
    pub fn shift(&self, s: i32) -> Elt {
        Elt { mu: self.mu.clone(), t: self.t.iter().map(|(p, c)| (*p, c.shift(s))).collect() }
    }
    /// self + k v^s o (same block)
    pub fn add_scaled(&self, o: &Elt, k: i128, s: i32) -> Elt {
        if self.is_zero() && k == 1 && s == 0 {
            return o.clone();
        }
        assert!(self.is_zero() || o.is_zero() || self.mu == o.mu, "Elt::add: different blocks");
        let mu = if self.is_zero() { o.mu.clone() } else { self.mu.clone() };
        let mut out = Vec::with_capacity(self.t.len() + o.t.len());
        let (mut a, mut b) = (0, 0);
        while a < self.t.len() || b < o.t.len() {
            let ta = self.t.get(a).map(|e| e.0).unwrap_or(Tid::MAX);
            let tb = o.t.get(b).map(|e| e.0).unwrap_or(Tid::MAX);
            if ta < tb {
                out.push(self.t[a].clone());
                a += 1;
            } else if tb < ta {
                let mut c = Lp::zero();
                c.add_scaled(&o.t[b].1, k, s);
                if !c.is_zero() {
                    out.push((tb, c));
                }
                b += 1;
            } else {
                let mut c = self.t[a].1.clone();
                c.add_scaled(&o.t[b].1, k, s);
                if !c.is_zero() {
                    out.push((ta, c));
                }
                a += 1;
                b += 1;
            }
        }
        Elt { mu, t: out }
    }
    pub fn add(&self, o: &Elt) -> Elt {
        self.add_scaled(o, 1, 0)
    }
}

fn finish_acc(mu: Comp, acc: FxMap<Tid, Lp>) -> Elt {
    let mut t: Vec<(Tid, Lp)> = acc
        .into_iter()
        .filter_map(|(p, mut c)| {
            c.normalize();
            if c.is_zero() { None } else { Some((p, c)) }
        })
        .collect();
    t.sort_unstable_by_key(|e| e.0);
    Elt { mu, t }
}

// ---------------------------------------------------------------------------------------------
// canonical basis container

pub struct CanonicalBasis {
    pub mu: Comp,
    pub order: Vec<Tid>,                                            // all tables of (lambda, mu) by ascending lmax
    pub b: std::collections::HashMap<Tid, Elt>,                     // b_p, with coefficient 1 at p
    pub ladder_bfs: std::collections::HashMap<Tid, Vec<Rung>>,      // first found (schur.py order)
    pub ladder_minflow: std::collections::HashMap<Tid, Vec<Rung>>, // min flows, then length, then lexicographic
}

#[derive(Clone, Debug)]
struct TInfo {
    lmin: i64,
    lmax: i64,
    ncos: u128,
    cols: Comp,
}

fn rung_code(r: &Rung) -> u32 {
    assert!(r.i < 256 && r.r < 32768, "rung out of range");
    ((r.i as u32) << 16) | (if r.dir == Dir::F { 1 << 15 } else { 0 }) | r.r as u32
}
/// lexicographic key of a rung: (i, E < F, r), as Python tuples (i, 'E'|'F', r)
fn rung_lex(r: &Rung) -> (usize, u8, usize) {
    (r.i, if r.dir == Dir::E { 0 } else { 1 }, r.r)
}
fn rung_str(r: &Rung) -> String {
    format!("{}{}{}", r.i, if r.dir == Dir::E { 'E' } else { 'F' }, r.r)
}
fn parse_rung(tok: &str) -> Option<Rung> {
    let p = tok.find(['E', 'F'])?;
    let i = tok[..p].parse().ok()?;
    let r = tok[p + 1..].parse().ok()?;
    let dir = if &tok[p..p + 1] == "E" { Dir::E } else { Dir::F };
    Some(Rung { i, dir, r })
}
fn comp_str(c: &[u8]) -> String {
    c.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")
}
fn parse_comp(s: &str) -> Option<Comp> {
    s.trim().split(',').map(|x| x.trim().parse::<u8>().ok()).collect()
}
/// labels after applying a rung (panics if the giver is too small)
fn rung_labels(mu: &[u8], r: &Rung) -> Comp {
    let (rec, giv) = match r.dir {
        Dir::E => (r.i, r.i + 1),
        Dir::F => (r.i + 1, r.i),
    };
    assert!(giv < mu.len() && rec < mu.len(), "rung {:?} outside labels {:?}", r, mu);
    assert!(mu[giv] as usize >= r.r, "rung {:?} needs giver >= r on labels {:?}", r, mu);
    let mut l = mu.to_vec();
    l[giv] -= r.r as u8;
    l[rec] += r.r as u8;
    l
}

/// Refinements of column j into consecutive sub-columns with sums `parts` (split tables).
fn split_table(t: &Table, j: usize, parts: &[u8]) -> Vec<Table> {
    let n = t.nr as usize;
    let m = t.nc as usize;
    let k = parts.len();
    let nm = m + k - 1;
    let mut out = vec![];
    let mut pieces = vec![0u8; n * k];
    let mut left = parts.to_vec();
    fn rec(row: usize, l: usize, rem: u8, t: &Table, j: usize, k: usize, pieces: &mut Vec<u8>, left: &mut Vec<u8>, out: &mut Vec<Vec<u8>>) {
        let n = t.nr as usize;
        if row == n {
            if left.iter().all(|&x| x == 0) {
                out.push(pieces.clone());
            }
            return;
        }
        if l == k - 1 {
            if rem <= left[l] {
                pieces[row * k + l] = rem;
                left[l] -= rem;
                let nrem = if row + 1 < n { t.get(row + 1, j) } else { 0 };
                rec(row + 1, 0, nrem, t, j, k, pieces, left, out);
                left[l] += rem;
                pieces[row * k + l] = 0;
            }
            return;
        }
        for y in 0..=rem.min(left[l]) {
            pieces[row * k + l] = y;
            left[l] -= y;
            rec(row, l + 1, rem - y, t, j, k, pieces, left, out);
            left[l] += y;
        }
        pieces[row * k + l] = 0;
    }
    let mut raw = vec![];
    rec(0, 0, t.get(0, j), t, j, k, &mut pieces, &mut left, &mut raw);
    for pc in raw {
        let mut e = Vec::with_capacity(n * nm);
        for i in 0..n {
            for jj in 0..m {
                if jj == j {
                    e.extend_from_slice(&pc[i * k..(i + 1) * k]);
                } else {
                    e.push(t.get(i, jj));
                }
            }
        }
        out.push(Table { nr: n as u8, nc: nm as u8, e });
    }
    out
}

/// Merge columns j, j+1: (new table, v^{lmin(new) - lmin(old)} prod_i qbinom(M_ij + M_ij+1, M_ij)).
fn merge_table(t: &Table, j: usize, qbin: &[Vec<Lp>]) -> (Table, Lp) {
    let n = t.nr as usize;
    let m = t.nc as usize;
    assert!(j + 1 < m, "merge_cols: no column {}", j + 1);
    let mut e = Vec::with_capacity(n * (m - 1));
    let mut coef = Lp::one();
    for i in 0..n {
        for jj in 0..m {
            if jj == j {
                let (a, b) = (t.get(i, j), t.get(i, j + 1));
                e.push(a + b);
                if a > 0 && b > 0 {
                    coef = coef.mul(&qbin[(a + b) as usize][a as usize]);
                }
            } else if jj != j + 1 {
                e.push(t.get(i, jj));
            }
        }
    }
    let nt = Table { nr: n as u8, nc: (m - 1) as u8, e };
    let s = nt.lmin() - t.lmin();
    (nt, coef.shift(s as i32))
}

// ---------------------------------------------------------------------------------------------
// the algebroid

pub struct Schur {
    pub lam: Comp,
    tabs: Vec<Table>,
    index: FxMap<Table, Tid>,
    info: Vec<TInfo>,
    rcache: FxMap<u64, Vec<(Tid, Lp)>>,
    bases: HashMap<Comp, CanonicalBasis>,
    rank_cache: FxMap<Tid, u128>,
    tcache: FxMap<(Tid, Vec<u32>), Vec<(Tid, Lp)>>,
    qbin: Vec<Vec<Lp>>,
}

/// Heap entry of the minimal-flow search: ordered by (flows, length, lexicographic rung list).
struct MfEntry {
    flows: u128,
    path: Vec<Rung>,
    parent: Option<Rc<Elt>>,
}
impl MfEntry {
    fn key(&self) -> (u128, usize) {
        (self.flows, self.path.len())
    }
}
impl PartialEq for MfEntry {
    fn eq(&self, o: &Self) -> bool {
        self.cmp(o) == Ordering::Equal
    }
}
impl Eq for MfEntry {}
impl PartialOrd for MfEntry {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for MfEntry {
    // reversed: BinaryHeap pops the smallest (flows, len, lex)
    fn cmp(&self, o: &Self) -> Ordering {
        o.key().cmp(&self.key()).then_with(|| {
            let a: Vec<_> = o.path.iter().map(rung_lex).collect();
            let b: Vec<_> = self.path.iter().map(rung_lex).collect();
            a.cmp(&b)
        })
    }
}

impl Schur {
    pub fn new(lam: &[u8]) -> Schur {
        let k: usize = lam.iter().map(|&a| a as usize).sum();
        assert!(k <= 64, "K = sum of labels must be <= 64");
        // symmetric quantum binomials by Pascal: [n,k] = v^{-(n-k)}[n-1,k-1] + v^k [n-1,k]
        let mut qbin: Vec<Vec<Lp>> = vec![vec![Lp::one()]];
        for n in 1..=k {
            let mut row = vec![Lp::one(); n + 1];
            for kk in 1..n {
                let mut c = qbin[n - 1][kk - 1].shift(-((n - kk) as i32));
                c.add_scaled(&qbin[n - 1][kk], 1, kk as i32);
                row[kk] = c;
            }
            qbin.push(row);
        }
        Schur {
            lam: lam.to_vec(),
            tabs: vec![],
            index: FxMap::default(),
            info: vec![],
            rcache: FxMap::default(),
            bases: HashMap::new(),
            rank_cache: FxMap::default(),
            tcache: FxMap::default(),
            qbin,
        }
    }
    pub fn intern(&mut self, t: Table) -> Tid {
        if let Some(&id) = self.index.get(&t) {
            return id;
        }
        let id = self.tabs.len() as Tid;
        let info = TInfo { lmin: t.lmin(), lmax: t.lmax(), ncos: t.ncosets(), cols: t.cols() };
        self.tabs.push(t.clone());
        self.info.push(info);
        self.index.insert(t, id);
        id
    }
    pub fn table(&self, id: Tid) -> &Table {
        &self.tabs[id as usize]
    }
    /// cached l_max of an interned table
    pub fn lmax_of(&self, id: Tid) -> i64 {
        self.info[id as usize].lmax
    }
    /// cached l_min of an interned table
    pub fn lmin_of(&self, id: Tid) -> i64 {
        self.info[id as usize].lmin
    }
    /// cached N(p) of an interned table
    pub fn ncosets_of(&self, id: Tid) -> u128 {
        self.info[id as usize].ncos
    }
    /// h_{diag(lambda)}
    pub fn unit(&mut self) -> Elt {
        let lam = self.lam.clone();
        let d = self.intern(Table::identity(&lam));
        Elt { mu: lam, t: vec![(d, Lp::one())] }
    }
    /// EKLP (44): refine column j into consecutive parts.
    pub fn split_col(&mut self, x: &Elt, j: usize, parts: &[u8]) -> Elt {
        assert_eq!(parts.iter().map(|&a| a as u32).sum::<u32>(), x.mu[j] as u32, "split_col: parts must sum to the column");
        let mut mu = x.mu[..j].to_vec();
        mu.extend_from_slice(parts);
        mu.extend_from_slice(&x.mu[j + 1..]);
        let mut acc: FxMap<Tid, Lp> = FxMap::default();
        for (p, c) in &x.t {
            let t = self.tabs[*p as usize].clone();
            let lq = self.info[*p as usize].lmax;
            for s in split_table(&t, j, parts) {
                let sid = self.intern(s);
                let e = (lq - self.info[sid as usize].lmax) as i32;
                acc.entry(sid).or_default().add_mul(c, &Lp::mono(e, 1), 1);
            }
        }
        finish_acc(mu, acc)
    }
    /// EKLP (45): merge columns j, j+1.
    pub fn merge_cols(&mut self, x: &Elt, j: usize) -> Elt {
        let mut mu = x.mu[..j].to_vec();
        mu.push(x.mu[j] + x.mu[j + 1]);
        mu.extend_from_slice(&x.mu[j + 2..]);
        let mut acc: FxMap<Tid, Lp> = FxMap::default();
        for (p, c) in &x.t {
            let t = self.tabs[*p as usize].clone();
            let (nt, coef) = merge_table(&t, j, &self.qbin);
            let nid = self.intern(nt);
            acc.entry(nid).or_default().add_mul(c, &coef, 1);
        }
        finish_acc(mu, acc)
    }
    /// image of h_p under a rung (split then merge), cached per (table, rung)
    fn ensure_rung_image(&mut self, p: Tid, r: &Rung) -> u64 {
        let key = ((p as u64) << 32) | rung_code(r) as u64;
        if self.rcache.contains_key(&key) {
            return key;
        }
        let t = self.tabs[p as usize].clone();
        let mu = &self.info[p as usize].cols;
        let lq = self.info[p as usize].lmax;
        let (sj, parts, mj) = match r.dir {
            // E: split column i+1 as (r, b-r), merge columns i, i+1
            Dir::E => (r.i + 1, vec![r.r as u8, mu[r.i + 1] - r.r as u8], r.i),
            // F: split column i as (a-r, r), merge columns i+1, i+2
            Dir::F => (r.i, vec![mu[r.i] - r.r as u8, r.r as u8], r.i + 1),
        };
        let mut img: Vec<(Tid, Lp)> = vec![];
        for s in split_table(&t, sj, &parts) {
            let e = (lq - s.lmax()) as i32;
            let (nt, coef) = merge_table(&s, mj, &self.qbin);
            let nid = self.intern(nt);
            img.push((nid, coef.shift(e)));
        }
        img.sort_unstable_by_key(|e| e.0);
        // distinct split tables give distinct merged tables (the r-part is recorded in the result)
        debug_assert!(img.windows(2).all(|w| w[0].0 != w[1].0));
        self.rcache.insert(key, img);
        key
    }
    /// E_i^(r): split col i+1 as (r, b-r), merge cols i, i+1;
    /// F_i^(r): split col i as (a-r, r), merge cols i+1, i+2
    pub fn rung(&mut self, x: &Elt, r: &Rung) -> Elt {
        if r.r == 0 {
            return x.clone();
        }
        let mu = rung_labels(&x.mu, r);
        let keys: Vec<u64> = x.t.iter().map(|(p, _)| self.ensure_rung_image(*p, r)).collect();
        let mut acc: FxMap<Tid, Lp> = FxMap::default();
        acc.reserve(x.t.len() * 2);
        for ((_, c), key) in x.t.iter().zip(keys) {
            for (q, k) in &self.rcache[&key] {
                acc.entry(*q).or_default().add_mul(c, k, 1);
            }
        }
        finish_acc(mu, acc)
    }
    pub fn mul_rungs(&mut self, x: &Elt, rs: &[Rung]) -> Elt {
        let mut y = x.clone();
        for r in rs {
            y = self.rung(&y, r);
        }
        y
    }
    /// unit * rungs
    pub fn ladder_char(&mut self, rs: &[Rung]) -> Elt {
        let u = self.unit();
        self.mul_rungs(&u, rs)
    }

    /// leading tables (maximal lmax) and whether the element is clean (unique leading table, coefficient 1)
    fn leading(&self, x: &Elt) -> (Vec<Tid>, bool) {
        let mut best = i64::MIN;
        let mut tops: Vec<Tid> = vec![];
        let mut clean = false;
        for (p, c) in &x.t {
            let l = self.info[*p as usize].lmax;
            if l > best {
                best = l;
                tops.clear();
                tops.push(*p);
                clean = *c == Lp::one();
            } else if l == best {
                tops.push(*p);
            }
        }
        let clean = clean && tops.len() == 1;
        tops.sort_unstable();
        (tops, clean)
    }

    /// Breadth-first search of schur.py: ladders from 1_lambda over rungs, extending only clean
    /// states, deduplicated on the set of leading tables; records the first ladder for each table
    /// of (lambda, mu) together with its character.
    fn bfs_ladders(&mut self, mu: &[u8], targets: &FxSet<Tid>) -> FxMap<Tid, (Vec<Rung>, Elt)> {
        let n = self.lam.len();
        let mut found: FxMap<Tid, (Vec<Rung>, Elt)> = FxMap::default();
        let mut seen: FxSet<Vec<Tid>> = FxSet::default();
        let mut queue: VecDeque<(Vec<Rung>, Elt)> = VecDeque::new();
        let start = self.unit();
        let (tops, clean) = self.leading(&start);
        seen.insert(tops.clone());
        if clean {
            if start.mu == mu && targets.contains(&tops[0]) {
                found.insert(tops[0], (vec![], start.clone()));
            }
            queue.push_back((vec![], start));
        }
        'outer: while let Some((w, x)) = queue.pop_front() {
            if found.len() == targets.len() {
                break;
            }
            let cur = x.mu.clone();
            for i in 0..n.saturating_sub(1) {
                for dir in [Dir::E, Dir::F] {
                    let giver = if dir == Dir::E { cur[i + 1] } else { cur[i] } as usize;
                    for r in 1..=giver {
                        let rg = Rung { i, dir, r };
                        let y = self.rung(&x, &rg);
                        if y.is_zero() {
                            continue;
                        }
                        let (tops, clean) = self.leading(&y);
                        if !seen.insert(tops.clone()) {
                            continue;
                        }
                        if !clean {
                            continue;
                        }
                        let mut w2 = w.clone();
                        w2.push(rg);
                        if y.mu == mu && targets.contains(&tops[0]) && !found.contains_key(&tops[0]) {
                            found.insert(tops[0], (w2.clone(), y.clone()));
                            if found.len() == targets.len() {
                                break 'outer;
                            }
                        }
                        queue.push_back((w2, y));
                    }
                }
            }
        }
        assert_eq!(found.len(), targets.len(), "canonical basis BFS: not every table of {:?} -> {:?} was reached", self.lam, mu);
        found
    }

    /// Uniform-cost search over the same clean states, ordered by (flows, length, lexicographic
    /// rung list); the first pop of each table of (lambda, mu) gives its minimal-flow ladder.
    /// Unclean states are dropped (not marked), so they never block a clean state.
    /// Entries with flows above `bound` are not pushed (memory). If every target is still found,
    /// the pops form a prefix of the unbounded search's pops (flows never decrease along a path
    /// and pops come in priority order), so the result is that of the unbounded search; if not,
    /// the caller retries without a bound (deduplication could hide a ladder below the bound).
    fn minflow_ladders(&mut self, mu: &[u8], targets: &FxSet<Tid>, bound: u128) -> Option<FxMap<Tid, Vec<Rung>>> {
        let n = self.lam.len();
        let mut found: FxMap<Tid, Vec<Rung>> = FxMap::default();
        let mut seen: FxSet<Tid> = FxSet::default();
        let mut heap: BinaryHeap<MfEntry> = BinaryHeap::new();
        heap.push(MfEntry { flows: 1, path: vec![], parent: None });
        while let Some(ent) = heap.pop() {
            let y = match &ent.parent {
                None => self.unit(),
                Some(par) => {
                    let r = *ent.path.last().unwrap();
                    self.rung(par, &r)
                }
            };
            drop(ent.parent);
            if y.is_zero() {
                continue;
            }
            let (tops, clean) = self.leading(&y);
            if !clean || !seen.insert(tops[0]) {
                continue;
            }
            if y.mu == mu && targets.contains(&tops[0]) {
                found.insert(tops[0], ent.path.clone());
                if found.len() == targets.len() {
                    break;
                }
            }
            let cur = y.mu.clone();
            let yr = Rc::new(y);
            for i in 0..n.saturating_sub(1) {
                for dir in [Dir::E, Dir::F] {
                    let giver = if dir == Dir::E { cur[i + 1] } else { cur[i] } as usize;
                    for r in 1..=giver {
                        let f = ent.flows.checked_mul(binom_u128(giver as u64, r as u64)).expect("flows overflow");
                        if f > bound {
                            continue;
                        }
                        let mut p = ent.path.clone();
                        p.push(Rung { i, dir, r });
                        heap.push(MfEntry { flows: f, path: p, parent: Some(yr.clone()) });
                    }
                }
            }
        }
        if found.len() == targets.len() { Some(found) } else { None }
    }

    /// Soergel's algorithm: peel b_q (q below p, descending lmax) off the character x with
    /// leading term exactly h_p, subtracting sym_part(x[q]) b_q; asserts the result.
    fn peel(&self, p: Tid, x: &Elt, b: &HashMap<Tid, Elt>) -> Elt {
        let mut acc: FxMap<Tid, Lp> = x.t.iter().cloned().collect();
        let mut heap: BinaryHeap<(i64, Tid)> = x.t.iter().filter(|e| e.0 != p).map(|e| (self.info[e.0 as usize].lmax, e.0)).collect();
        let lp = self.info[p as usize].lmax;
        while let Some((_, q)) = heap.pop() {
            let m = match acc.get_mut(&q) {
                Some(c) => {
                    c.normalize();
                    c.sym_part()
                }
                None => continue,
            };
            if m.is_zero() {
                continue;
            }
            assert!(self.info[q as usize].lmax < lp, "peel: non-leading table above p");
            let bq = b.get(&q).unwrap_or_else(|| panic!("peel: b_q missing for {}", self.tabs[q as usize].fmt()));
            for (r, c) in &bq.t {
                let present = acc.contains_key(r);
                acc.entry(*r).or_default().add_mul(&m, c, -1);
                if !present && *r != q {
                    heap.push((self.info[*r as usize].lmax, *r));
                }
            }
        }
        let out = finish_acc(x.mu.clone(), acc);
        for (q, c) in &out.t {
            if *q == p {
                assert_eq!(*c, Lp::one(), "peel: coefficient of h_p is not 1");
            } else {
                assert!(c.min_exp().unwrap() > 0, "peel: coefficient of {} in b_{} not in vZ[v]: {}", self.tabs[*q as usize].fmt(), self.tabs[p as usize].fmt(), c.fmt_ec());
            }
        }
        out
    }

    /// Canonical basis of (lambda, mu) peeled from the given ladders (one per table), in
    /// ascending lmax. Used for the check that the basis does not depend on the ladders (T-c2).
    pub fn basis_from_ladders(&mut self, mu: &[u8], ladders: &HashMap<Tid, Vec<Rung>>) -> HashMap<Tid, Elt> {
        let order: Vec<Tid> = tables(&self.lam, mu).into_iter().map(|t| self.intern(t)).collect();
        let mut b: HashMap<Tid, Elt> = HashMap::new();
        for &p in &order {
            let lad = ladders.get(&p).unwrap_or_else(|| panic!("no ladder for {}", self.tabs[p as usize].fmt()));
            let x = self.ladder_char(lad);
            let (tops, clean) = self.leading(&x);
            assert!(clean && tops == [p], "ladder for {} does not have leading term h_p", self.tabs[p as usize].fmt());
            let bp = self.peel(p, &x, &b);
            b.insert(p, bp);
        }
        b
    }

    fn ensure_basis(&mut self, mu: &[u8]) {
        if self.bases.contains_key(mu) {
            return;
        }
        let order: Vec<Tid> = tables(&self.lam, mu).into_iter().map(|t| self.intern(t)).collect();
        assert!(!order.is_empty(), "empty block {:?} -> {:?}", self.lam, mu);
        let targets: FxSet<Tid> = order.iter().copied().collect();
        let found = self.bfs_ladders(mu, &targets);
        let mut b: HashMap<Tid, Elt> = HashMap::new();
        for &p in &order {
            let (_, x) = &found[&p];
            let bp = self.peel(p, x, &b);
            b.insert(p, bp);
        }
        let ladder_bfs: HashMap<Tid, Vec<Rung>> = found.into_iter().map(|(p, (w, _))| (p, w)).collect();
        let bound = ladder_bfs.values().map(|w| flows(&self.lam, w)).max().unwrap_or(1);
        let mf = match self.minflow_ladders(mu, &targets, bound) {
            Some(m) => m,
            None => self.minflow_ladders(mu, &targets, u128::MAX).unwrap_or_else(|| {
                panic!("minimal-flow search: not every table of {:?} -> {:?} was reached", self.lam, mu)
            }),
        };
        let ladder_minflow: HashMap<Tid, Vec<Rung>> = mf.into_iter().collect();
        self.bases.insert(mu.to_vec(), CanonicalBasis { mu: mu.to_vec(), order, b, ladder_bfs, ladder_minflow });
    }
    /// lazy, cached
    pub fn canonical_basis(&mut self, mu: &[u8]) -> &CanonicalBasis {
        self.ensure_basis(mu);
        &self.bases[mu]
    }
    /// x = sum m_q b_q; panics unless all m_q in N[v,v^-1]
    pub fn expand(&mut self, x: &Elt) -> Vec<(Tid, Lp)> {
        if x.is_zero() {
            return vec![];
        }
        self.ensure_basis(&x.mu);
        let basis = &self.bases[&x.mu].b;
        let mut acc: FxMap<Tid, Lp> = x.t.iter().cloned().collect();
        let mut heap: BinaryHeap<(i64, Tid)> = x.t.iter().map(|e| (self.info[e.0 as usize].lmax, e.0)).collect();
        let mut out = vec![];
        while let Some((_, t)) = heap.pop() {
            let mut m = match acc.remove(&t) {
                Some(c) => c,
                None => continue,
            };
            m.normalize();
            if m.is_zero() {
                continue;
            }
            assert!(m.is_nonneg(), "expand: negative multiplicity {} of b_{}", m.fmt_ec(), self.tabs[t as usize].fmt());
            let bt = &basis[&t];
            for (q, c) in &bt.t {
                if *q == t {
                    continue;
                }
                let present = acc.contains_key(q);
                acc.entry(*q).or_default().add_mul(&m, c, -1);
                if !present {
                    heap.push((self.info[*q as usize].lmax, *q));
                }
            }
            out.push((t, m));
        }
        out.sort_unstable_by_key(|e| e.0);
        out
    }
    /// EKLP Thm 4.16 (y self-dual): gd Hom(x, y) = sum_M x[M] y[M] gd(R_M), degrees dlo..=dhi
    pub fn hom_gd(&self, x: &Elt, y: &Elt, dlo: i32, dhi: i32) -> Vec<i128> {
        assert!(dhi >= dlo);
        let mut out = vec![0i128; (dhi - dlo + 1) as usize];
        if x.mu != y.mu {
            return out;
        }
        let (mut a, mut b) = (0, 0);
        while a < x.t.len() && b < y.t.len() {
            let (ta, tb) = (x.t[a].0, y.t[b].0);
            if ta < tb {
                a += 1;
                continue;
            }
            if tb < ta {
                b += 1;
                continue;
            }
            let cc = x.t[a].1.mul(&y.t[b].1);
            a += 1;
            b += 1;
            let lo = match cc.min_exp() {
                Some(l) => l,
                None => continue,
            };
            if dhi < lo {
                continue;
            }
            let dmax = (dhi - lo) as usize;
            // gd(R_M) = prod_ij prod_{k=1..M_ij} 1/(1 - v^{2k}), up to degree dmax
            let mut ser = vec![0i128; dmax + 1];
            ser[0] = 1;
            for &x in &self.tabs[ta as usize].e {
                for k in 1..=x as usize {
                    let st = 2 * k;
                    for d in st..=dmax {
                        ser[d] = ser[d].checked_add(ser[d - st]).expect(OVF);
                    }
                }
            }
            for (i, &k) in cc.c.iter().enumerate() {
                if k == 0 {
                    continue;
                }
                let e = lo + i as i32;
                for d in dlo.max(e)..=dhi {
                    let s = ser[(d - e) as usize];
                    if s != 0 {
                        let idx = (d - dlo) as usize;
                        out[idx] = out[idx].checked_add(k.checked_mul(s).expect(OVF)).expect(OVF);
                    }
                }
            }
        }
        out
    }
    /// rank over R^lambda = sum_q c_q(b_p)(1) N(q)
    pub fn rank(&mut self, mu: &[u8], p: Tid) -> u128 {
        if let Some(&r) = self.rank_cache.get(&p) {
            return r;
        }
        self.ensure_basis(mu);
        let bp = &self.bases[mu].b[&p];
        let mut r: u128 = 0;
        for (q, c) in &bp.t {
            let a = c.at1();
            assert!(a >= 0);
            r = r.checked_add((a as u128).checked_mul(self.info[*q as usize].ncos).expect("rank overflow")).expect("rank overflow");
        }
        self.rank_cache.insert(p, r);
        r
    }
    /// coefficient of h_diag in b_p (mu = lambda)
    pub fn c_diag(&mut self, p: Tid) -> Lp {
        let lam = self.lam.clone();
        let d = self.intern(Table::identity(&lam));
        self.ensure_basis(&lam);
        self.bases[&lam].b[&p].coef(d)
    }
    /// max_{p,q} c_q(b_p)(1)
    pub fn kappa(&mut self, mu: &[u8]) -> i128 {
        self.ensure_basis(mu);
        self.bases[mu].b.values().flat_map(|e| e.t.iter().map(|(_, c)| c.at1())).max().unwrap_or(0)
    }
    /// S(p,q) = sum_r c_r(b_p)(1) c_r(b_q)(1) N(r)
    pub fn s_storage(&mut self, mu: &[u8], p: Tid, q: Tid) -> u128 {
        self.ensure_basis(mu);
        let cb = &self.bases[mu];
        let (bp, bq) = (&cb.b[&p], &cb.b[&q]);
        let mut s: u128 = 0;
        for (r, c) in &bp.t {
            let d = bq.coef(*r).at1();
            if d == 0 {
                continue;
            }
            let t = (c.at1() as u128).checked_mul(d as u128).and_then(|x| x.checked_mul(self.info[*r as usize].ncos)).expect("S overflow");
            s = s.checked_add(t).expect("S overflow");
        }
        s
    }
    /// expand(b_p * [rungs]), cached per (p, rungs)
    pub fn term_mult(&mut self, p: Tid, rs: &[Rung]) -> Vec<(Tid, Lp)> {
        let key = (p, rs.iter().map(rung_code).collect::<Vec<u32>>());
        if let Some(v) = self.tcache.get(&key) {
            return v.clone();
        }
        let mu = self.info[p as usize].cols.clone();
        self.ensure_basis(&mu);
        let bp = self.bases[&mu].b[&p].clone();
        let y = self.mul_rungs(&bp, rs);
        let e = self.expand(&y);
        self.tcache.insert(key, e.clone());
        e
    }
    /// header "schur-cb v1 lam=.. mu=.."; then CB, LAD (bfs) and LADMF (minimal-flow) lines
    pub fn save_basis(&self, mu: &[u8], path: &std::path::Path) -> std::io::Result<()> {
        use std::io::Write;
        let cb = self.bases.get(mu).ok_or_else(|| std::io::Error::new(std::io::ErrorKind::NotFound, "basis not computed"))?;
        let lam = comp_str(&self.lam);
        let mut f = std::io::BufWriter::new(std::fs::File::create(path)?);
        writeln!(f, "schur-cb v1 lam={} mu={}", lam, comp_str(mu))?;
        let ord = |q: &Tid| (self.info[*q as usize].lmax, self.tabs[*q as usize].e.clone());
        for p in &cb.order {
            let mut ent: Vec<&(Tid, Lp)> = cb.b[p].t.iter().collect();
            ent.sort_by_key(|e| ord(&e.0));
            let body: Vec<String> = ent.iter().map(|(q, c)| format!("{} {}", self.tabs[*q as usize].fmt(), c.fmt_ec())).collect();
            writeln!(f, "CB {}|{}|{}", lam, self.tabs[*p as usize].fmt(), body.join("|"))?;
        }
        for (tag, lads) in [("LAD", &cb.ladder_bfs), ("LADMF", &cb.ladder_minflow)] {
            for p in &cb.order {
                let w: Vec<String> = lads[p].iter().map(rung_str).collect();
                writeln!(f, "{} {}|{}|{}", tag, lam, self.tabs[*p as usize].fmt(), w.join(" "))?;
            }
        }
        f.flush()
    }
    pub fn load_basis(&mut self, path: &std::path::Path) -> std::io::Result<()> {
        let bad = |m: String| std::io::Error::new(std::io::ErrorKind::InvalidData, m);
        let text = std::fs::read_to_string(path)?;
        let mut lines = text.lines();
        let header = lines.next().ok_or_else(|| bad("empty basis file".into()))?;
        let hp: Vec<&str> = header.split_whitespace().collect();
        if hp.len() != 4 || hp[0] != "schur-cb" || hp[1] != "v1" {
            return Err(bad(format!("bad header {:?}", header)));
        }
        let lam = hp[2].strip_prefix("lam=").and_then(parse_comp).ok_or_else(|| bad("bad lam".into()))?;
        let mu = hp[3].strip_prefix("mu=").and_then(parse_comp).ok_or_else(|| bad("bad mu".into()))?;
        if lam != self.lam {
            return Err(bad(format!("basis file is for lam {:?}, this instance has {:?}", lam, self.lam)));
        }
        let mut b: HashMap<Tid, Elt> = HashMap::new();
        let mut lb: HashMap<Tid, Vec<Rung>> = HashMap::new();
        let mut lm: HashMap<Tid, Vec<Rung>> = HashMap::new();
        for line in lines {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let (tag, rest) = line.split_once(' ').ok_or_else(|| bad(format!("bad line {:?}", line)))?;
            let f: Vec<&str> = rest.split('|').collect();
            if f.len() < 3 || parse_comp(f[0]).as_deref() != Some(&lam[..]) {
                return Err(bad(format!("bad line {:?}", line)));
            }
            let p = self.intern(Table::parse(f[1]));
            match tag {
                "CB" => {
                    let mut acc: FxMap<Tid, Lp> = FxMap::default();
                    for ent in &f[2..] {
                        let ent = ent.trim();
                        if ent.is_empty() {
                            continue;
                        }
                        let (tq, poly) = ent.split_once(' ').ok_or_else(|| bad(format!("bad entry {:?}", ent)))?;
                        let q = self.intern(Table::parse(tq));
                        acc.entry(q).or_default().add_scaled(&Lp::parse_ec(poly), 1, 0);
                    }
                    b.insert(p, finish_acc(mu.clone(), acc));
                }
                "LAD" | "LADMF" => {
                    let w: Option<Vec<Rung>> = f[2].split_whitespace().map(parse_rung).collect();
                    let w = w.ok_or_else(|| bad(format!("bad ladder {:?}", line)))?;
                    if tag == "LAD" {
                        lb.insert(p, w);
                    } else {
                        lm.insert(p, w);
                    }
                }
                _ => return Err(bad(format!("unknown tag {:?}", tag))),
            }
        }
        let order: Vec<Tid> = tables(&self.lam, &mu).into_iter().map(|t| self.intern(t)).collect();
        for p in &order {
            if !b.contains_key(p) || !lb.contains_key(p) || !lm.contains_key(p) {
                return Err(bad(format!("basis file misses table {}", self.tabs[*p as usize].fmt())));
            }
        }
        if b.len() != order.len() {
            return Err(bad("basis file has tables outside the block".into()));
        }
        self.rank_cache.clear();
        self.tcache.clear();
        self.bases.insert(mu.clone(), CanonicalBasis { mu, order, b, ladder_bfs: lb, ladder_minflow: lm });
        Ok(())
    }
}

// ---------------------------------------------------------------------------------------------
// crossings, flows, closure weights

pub struct Term {
    pub rungs: Vec<Rung>,
    pub hdeg: i32,
    pub dq: i32,
}
/// = t8.py cross_terms; zero rungs dropped. Untruncated (all k = 0..min(c,d)).
/// Positive: hdeg = k, dq = -k; negative: hdeg = -k, dq = k (N-independent part of convention A).
pub fn crossing_terms(c: u8, d: u8, positive: bool, i: usize) -> Vec<Term> {
    let (c, d) = (c as usize, d as usize);
    let m = c.min(d);
    let mut out = vec![];
    for k in 0..=m {
        let rs = if positive {
            if c >= d { vec![Rung::e(i, k), Rung::f(i, c - d + k)] } else { vec![Rung::f(i, k), Rung::e(i, d - c + k)] }
        } else if d >= c {
            vec![Rung::e(i, d - c + k), Rung::f(i, k)]
        } else {
            vec![Rung::f(i, c - d + k), Rung::e(i, k)]
        };
        let k = k as i32;
        out.push(Term {
            rungs: rs.into_iter().filter(|r| r.r > 0).collect(),
            hdeg: if positive { k } else { -k },
            dq: if positive { -k } else { k },
        });
    }
    out
}
/// prod over rungs of C(giver, r)  (fibre dimension of the ladder)
pub fn flows(bottom: &[u8], rs: &[Rung]) -> u128 {
    let mut lab = bottom.to_vec();
    let mut f: u128 = 1;
    for r in rs {
        if r.r == 0 {
            continue;
        }
        let giver = match r.dir {
            Dir::E => lab[r.i + 1],
            Dir::F => lab[r.i],
        };
        f = f.checked_mul(binom_u128(giver as u64, r.r as u64)).expect("flows overflow");
        lab = rung_labels(&lab, r);
    }
    f
}
/// e-monomials of weight w in R^lambda (e_j of block i has weight j, j = 1..lambda_i)
pub fn emono_count(lam: &[u8], w: usize) -> u128 {
    let mut ser = vec![0u128; w + 1];
    ser[0] = 1;
    for &a in lam {
        for j in 1..=a as usize {
            for d in j..=w {
                ser[d] = ser[d].checked_add(ser[d - j]).expect("emono overflow");
            }
        }
    }
    ser[w]
}

// ---------------------------------------------------------------------------------------------
// decategorified scan

pub struct ScanOpts {
    pub peak: bool,
    pub every_stage: bool,
    pub marked_block: Option<usize>,
}
pub struct StageStats {
    pub step: usize,
    pub labels: Comp,
    pub unsimplified_objects: u128,
    pub unsimplified_ldim: u128,
    pub lb_objects: u128,
    pub lb_types: usize,
    pub lb_ldim: u128,
    pub lb_sum_s: u128,
    pub max_rank_present: u128,
    pub kappa_present: i128,
    pub peak_before_next: Option<(u128, u128, u128)>, // (objects, ldim, sum S) of LB (x) next crossing, expanded
}
pub struct GramRow {
    pub p: Tid,
    pub weight: i32,
    pub unknowns: u128,
    pub unknowns_reduced: u128,
}
pub struct ScanResult {
    pub stages: Vec<StageStats>,
    pub k0: Vec<(Tid, i32, i128)>, // signed multiset (type, q-shift, sum_h (-1)^h mult): the K_0 class
    pub gram: Vec<GramRow>,        // for the final LB types (block (lambda, lambda))
}

type PeakCache = FxMap<(Tid, Vec<u32>), (u128, u128, u128)>;

fn stage_stats(
    s: &mut Schur,
    x: &BTreeMap<i32, Elt>,
    labels: &[u8],
    step: usize,
    next: Option<(u8, u8, bool, usize)>,
    pcache: &mut PeakCache,
) -> (StageStats, BTreeMap<(Tid, i32), i128>) {
    let mut uobj: u128 = 0;
    let mut uldim: u128 = 0;
    let mut signed: BTreeMap<(Tid, i32), i128> = BTreeMap::new();
    for (h, e) in x {
        let sg: i128 = if h.rem_euclid(2) == 0 { 1 } else { -1 };
        for (p, m) in s.expand(e) {
            let n = m.at1() as u128;
            uobj = uobj.checked_add(n).expect("count overflow");
            let rk = s.rank(labels, p);
            uldim = uldim.checked_add(n.checked_mul(rk).expect("count overflow")).expect("count overflow");
            for (i, &c) in m.c.iter().enumerate() {
                if c != 0 {
                    let ent = signed.entry((p, m.lo + i as i32)).or_insert(0);
                    *ent = ent.checked_add(sg * c).expect(OVF);
                }
            }
        }
    }
    signed.retain(|_, c| *c != 0);
    let mut lb_objects: u128 = 0;
    let mut lb_ldim: u128 = 0;
    let mut lb_sum_s: u128 = 0;
    let mut types: Vec<Tid> = signed.keys().map(|k| k.0).collect();
    types.dedup();
    let mut sself: FxMap<Tid, u128> = FxMap::default();
    for &p in &types {
        let v = s.s_storage(labels, p, p);
        sself.insert(p, v);
    }
    for (&(p, _), &c) in &signed {
        let a = c.unsigned_abs();
        lb_objects += a;
        lb_ldim += a * s.rank(labels, p);
        lb_sum_s += a * sself[&p];
    }
    let mut max_rank_present = 0;
    let mut kappa_present = 0;
    for &p in &types {
        max_rank_present = max_rank_present.max(s.rank(labels, p));
        let k = s.canonical_basis(labels).b[&p].t.iter().map(|(_, c)| c.at1()).max().unwrap_or(0);
        kappa_present = kappa_present.max(k);
    }
    let peak_before_next = next.map(|(c, d, pos, i)| {
        let terms = crossing_terms(c, d, pos, i);
        let mut nl = labels.to_vec();
        nl.swap(i, i + 1);
        let (mut a, mut b, mut cs) = (0u128, 0u128, 0u128);
        for (&(p, _), &c) in &signed {
            let mult = c.unsigned_abs();
            for t in &terms {
                let key = (p, t.rungs.iter().map(rung_code).collect::<Vec<u32>>());
                let tr = match pcache.get(&key) {
                    Some(v) => *v,
                    None => {
                        let e = s.term_mult(p, &t.rungs);
                        let (mut x1, mut x2, mut x3) = (0u128, 0u128, 0u128);
                        for (q, m) in &e {
                            let k = m.at1() as u128;
                            x1 += k;
                            x2 += k * s.rank(&nl, *q);
                            x3 += k * s.s_storage(&nl, *q, *q);
                        }
                        pcache.insert(key, (x1, x2, x3));
                        (x1, x2, x3)
                    }
                };
                a += mult * tr.0;
                b += mult * tr.1;
                cs += mult * tr.2;
            }
        }
        (a, b, cs)
    });
    let st = StageStats {
        step,
        labels: labels.to_vec(),
        unsimplified_objects: uobj,
        unsimplified_ldim: uldim,
        lb_objects,
        lb_types: types.len(),
        lb_ldim,
        lb_sum_s,
        max_rank_present,
        kappa_present,
        peak_before_next,
    };
    (st, signed)
}

fn crossing_at(labels: &[u8], g: i32) -> (u8, u8, bool, usize) {
    assert!(g != 0, "braid generator 0");
    let i = (g.unsigned_abs() - 1) as usize;
    assert!(i + 1 < labels.len(), "generator {} needs {} strands", g, i + 2);
    (labels[i], labels[i + 1], g > 0, i)
}

/// Decategorified scan (t8.py `scan` + stats, decat_scan.py lower bound, peak.py peak model).
///
/// * Stage statistics are computed in the block (lambda, current labels).
/// * The final stage is always recorded. With `every_stage` or `peak`, every stage from step 0
///   (the identity complex) to the end is recorded; `peak_before_next` of the stage after j
///   crossings is the expanded LB complex tensored with crossing j+1 (peak.py), None at the end.
/// * `max_rank_present` and `kappa_present` are taken over the LB types (decat_scan.py).
/// * `gram` is filled when the final labels equal lambda; `marked_block` (if any) is removed
///   from lambda for `unknowns_reduced` (otherwise unknowns_reduced = unknowns).
pub fn decat_scan(s: &mut Schur, word: &[i32], opts: &ScanOpts) -> ScanResult {
    let lam = s.lam.clone();
    let mut labels = lam.clone();
    let mut x: BTreeMap<i32, Elt> = BTreeMap::new();
    x.insert(0, s.unit());
    let record_all = opts.every_stage || opts.peak;
    let mut stages = vec![];
    let mut pcache: PeakCache = FxMap::default();
    let mut final_signed = BTreeMap::new();
    for step in 0..=word.len() {
        if step > 0 {
            let (c, d, pos, i) = crossing_at(&labels, word[step - 1]);
            let terms = crossing_terms(c, d, pos, i);
            let mut y: BTreeMap<i32, Elt> = BTreeMap::new();
            for (h, e) in &x {
                for t in &terms {
                    let z = s.mul_rungs(e, &t.rungs);
                    let slot = y.entry(h + t.hdeg).or_insert_with(|| Elt::zero(&z.mu));
                    *slot = slot.add_scaled(&z, 1, t.dq);
                }
            }
            y.retain(|_, e| !e.is_zero());
            x = y;
            labels.swap(i, i + 1);
        }
        if record_all || step == word.len() {
            let next = if opts.peak && step < word.len() { Some(crossing_at(&labels, word[step])) } else { None };
            let (st, signed) = stage_stats(s, &x, &labels, step, next, &mut pcache);
            if step == word.len() {
                final_signed = signed;
            }
            stages.push(st);
        }
    }
    let k0: Vec<(Tid, i32, i128)> = final_signed.iter().map(|(&(p, k), &c)| (p, k, c)).collect();
    let mut gram = vec![];
    if labels == lam {
        let lw_lam = lw(&lam);
        let red: Comp = match opts.marked_block {
            Some(mb) => lam.iter().enumerate().filter(|(i, _)| *i != mb).map(|(_, &a)| a).collect(),
            None => lam.clone(),
        };
        let mut types: Vec<Tid> = k0.iter().map(|e| e.0).collect();
        types.sort_unstable();
        types.dedup();
        types.sort_by_key(|&p| (s.lmax_of(p), s.table(p).e.clone()));
        for p in types {
            let w = s.c_diag(p).max_exp().expect("c_diag(b_p) = 0");
            assert_eq!(w as i64, s.lmax_of(p) - lw_lam, "Gram weight != l(p_max) - l(w_lambda) for {}", s.table(p).fmt());
            gram.push(GramRow { p, weight: w, unknowns: emono_count(&lam, w as usize), unknowns_reduced: emono_count(&red, w as usize) });
        }
    }
    ScanResult { stages, k0, gram }
}

// ---------------------------------------------------------------------------------------------
// decategorified closure

fn closure_checked(lam: &[u8], rs: &[Rung], n: usize) -> LPoly {
    let mut lab: Vec<usize> = lam.iter().map(|&a| a as usize).collect();
    if lab.iter().any(|&a| a > n) {
        return LPoly::zero();
    }
    for r in rs {
        match r.apply(&lab, n) {
            Some(l) => lab = l,
            None => return LPoly::zero(), // a level with a label > N: that weight space is zero
        }
    }
    let lad = Ladder { n, bottom: lam.iter().map(|&a| a as usize).collect(), rungs: rs.to_vec() };
    crate::moy::closure_value(&lad)
}

fn trace_b(s: &mut Schur, p: Tid, n: usize, inv: bool, memo: &mut FxMap<Tid, LPoly>) -> LPoly {
    if let Some(v) = memo.get(&p) {
        return v.clone();
    }
    let lam = s.lam.clone();
    let lad = s.canonical_basis(&lam).ladder_minflow[&p].clone();
    let ch = s.ladder_char(&lad);
    let mult = s.expand(&ch);
    let mut v = closure_checked(&lam, &lad, n);
    for (q, m) in mult {
        if q == p {
            assert_eq!(m, Lp::one(), "ladder character of b_p has multiplicity != 1 at p");
            continue;
        }
        let tq = trace_b(s, q, n, inv, memo);
        let mq = if inv { m.bar() } else { m }.to_lpoly();
        v = v.sub(&mq.mul(&tq));
    }
    memo.insert(p, v.clone());
    v
}

/// Decategorified closure: chi(KhR_N) = q^{-sum_c eps_c m_c (N - m_c)} * sum_{(p,k)} signed * q^k * tr_N(b_p),
/// tr_N(b_p) = moy::closure_value(L_p) - sum_{q != p} m_q(q) tr_N(b_q), char(L_p) = sum m_q b_q,
/// with closure_value := 0 when L_p has an intermediate label > N.
pub fn decat_chi(s: &mut Schur, res: &ScanResult, word: &[i32], n_gl: usize) -> LPoly {
    decat_chi_conv(s, res, word, n_gl, false)
}

/// `decat_chi` with a choice of the identification of v: `v_to_qinv = false` is v -> q (the
/// convention of §2.3), `true` is the alternative v -> q^{-1} (kept only to show that test T-n
/// discriminates between the two).
pub fn decat_chi_conv(s: &mut Schur, res: &ScanResult, word: &[i32], n_gl: usize, v_to_qinv: bool) -> LPoly {
    let lam = s.lam.clone();
    let mut labels = lam.clone();
    let mut glob: i64 = 0;
    for &g in word {
        let (c, d, pos, i) = crossing_at(&labels, g);
        let m = c.min(d) as i64;
        let eps = if pos { 1 } else { -1 };
        glob -= eps * m * (n_gl as i64 - m);
        labels.swap(i, i + 1);
    }
    assert_eq!(labels, lam, "decat_chi: the braid closure needs matching colors");
    if lam.iter().any(|&a| a as usize > n_gl) {
        return LPoly::zero();
    }
    let mut memo: FxMap<Tid, LPoly> = FxMap::default();
    let mut total = LPoly::zero();
    for &(p, k, c) in &res.k0 {
        assert_eq!(s.table(p).cols(), lam, "decat_chi: K_0 class not in the block (lambda, lambda)");
        let t = trace_b(s, p, n_gl, v_to_qinv, &mut memo);
        total = total.add(&t.shift(if v_to_qinv { -k } else { k }).scale(c));
    }
    total.shift(glob as i32)
}
