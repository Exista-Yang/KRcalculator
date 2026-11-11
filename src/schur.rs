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
    /// "2,0/0,2" (panics on malformed input; see [`Table::try_parse`])
    pub fn parse(s: &str) -> Table {
        Table::try_parse(s).unwrap_or_else(|e| panic!("{}", e))
    }
    /// "2,0/0,2": rows separated by '/', entries by ','; Err on bad entries, ragged rows or
    /// more than 255 rows or columns.
    pub fn try_parse(s: &str) -> Result<Table, String> {
        let mut rows: Vec<Vec<u8>> = vec![];
        for r in s.trim().split('/') {
            let row: Option<Vec<u8>> = r.split(',').map(|x| x.trim().parse::<u8>().ok()).collect();
            rows.push(row.ok_or_else(|| format!("bad table {:?}", s))?);
        }
        let nc = rows[0].len();
        if rows.iter().any(|r| r.len() != nc) {
            return Err(format!("ragged table {:?}", s));
        }
        if rows.len() > 255 || nc > 255 {
            return Err(format!("table {:?} too large", s));
        }
        Ok(Table { nr: rows.len() as u8, nc: nc as u8, e: rows.concat() })
    }
    pub fn fmt(&self) -> String {
        let m = self.nc as usize;
        (0..self.nr as usize)
            .map(|i| self.e[i * m..(i + 1) * m].iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","))
            .collect::<Vec<_>>()
            .join("/")
    }
    /// Two-sided cell weight ν(p) of the canonical basis element b_p (a partition of K, parts in
    /// decreasing order): the RSK shape of the table with its columns reversed. See
    /// [`Schur::cell_weight`] for the rule and its status.
    ///
    /// RSK of a nonnegative integer matrix M (Knuth): the two-line array of the pairs (i, j),
    /// each repeated M_ij times, in lexicographic order; the j's are row-inserted (Schensted:
    /// an entry x bumps the leftmost entry > x of a row into the next row). The shape does not
    /// depend on whether rows or columns are reversed (Schützenberger: RSK(R M R) has the shape
    /// of RSK(M)), nor on transposing.
    pub fn cell_weight(&self) -> Comp {
        let (n, m) = (self.nr as usize, self.nc as usize);
        let mut rows: Vec<Vec<u8>> = vec![];
        for i in 0..n {
            for j in 0..m {
                // column-reversed table: entry (i, j) of M' is M[i][m-1-j]
                for _ in 0..self.e[i * m + (m - 1 - j)] {
                    let mut x = j as u8;
                    let mut r = 0;
                    loop {
                        if r == rows.len() {
                            rows.push(vec![x]);
                            break;
                        }
                        let pos = rows[r].partition_point(|&y| y <= x);
                        if pos == rows[r].len() {
                            rows[r].push(x);
                            break;
                        }
                        std::mem::swap(&mut rows[r][pos], &mut x);
                        r += 1;
                    }
                }
            }
        }
        rows.iter().map(|r| r.len() as u8).collect()
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
    /// self + o, or None if a coefficient overflows i128
    pub fn checked_add(&self, o: &Lp) -> Option<Lp> {
        if o.is_zero() {
            return Some(self.clone());
        }
        if self.is_zero() {
            return Some(o.clone());
        }
        let lo = self.lo.min(o.lo);
        let hi = self.hi().max(o.hi());
        let mut r = Lp { lo, c: vec![0; (hi - lo + 1) as usize] };
        for x in [self, o] {
            let off = (x.lo - lo) as usize;
            for (i, &a) in x.c.iter().enumerate() {
                r.c[off + i] = r.c[off + i].checked_add(a)?;
            }
        }
        r.normalize();
        Some(r)
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
    /// inverse of `fmt_ec` (whitespace-separated "e:c" tokens; repeated exponents add up);
    /// panics on malformed input, see [`Lp::try_parse_ec`]
    pub fn parse_ec(s: &str) -> Lp {
        Lp::try_parse_ec(s).unwrap_or_else(|e| panic!("{}", e))
    }
    /// inverse of `fmt_ec`; Err on a malformed token, a coefficient overflow or an exponent
    /// outside ±2^20 (dense storage)
    pub fn try_parse_ec(s: &str) -> Result<Lp, String> {
        Lp::try_parse_ec_in(s, -(1 << 20), 1 << 20)
    }
    /// As [`Lp::try_parse_ec`], with every exponent required to lie in [lo, hi] (the dense
    /// storage then has at most min(hi − lo, 2^20) + 1 coefficients).
    pub fn try_parse_ec_in(s: &str, lo: i32, hi: i32) -> Result<Lp, String> {
        let mut r = Lp::zero();
        for tok in s.split_whitespace() {
            let (e, k) = tok.split_once(':').ok_or_else(|| format!("bad term {:?}", tok))?;
            let k = k.parse::<i128>().map_err(|_| format!("bad coefficient in {:?}", tok))?;
            let e = e.parse::<i32>().map_err(|_| format!("bad exponent in {:?}", tok))?;
            if e < lo || e > hi {
                return Err(format!("exponent out of range [{}, {}] in {:?}", lo, hi, tok));
            }
            if let (Some(a), Some(b)) = (r.min_exp(), r.max_exp()) {
                if b.max(e) as i64 - a.min(e) as i64 > 1 << 20 {
                    return Err(format!("exponent range too wide at {:?}", tok));
                }
            }
            r.coef(e).checked_add(k).ok_or_else(|| format!("coefficient overflow at {:?}", tok))?;
            r.add_scaled(&Lp::one(), k, e);
        }
        Ok(r)
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

/// sum of the absolute values of the coefficients (saturating)
fn l1_norm(c: &Lp) -> u128 {
    c.c.iter().fold(0u128, |s, &a| s.saturating_add(a.unsigned_abs()))
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

/// (i, dir, r) packed so that the numeric order is the order of `rung_lex`
fn rung_code(r: &Rung) -> u32 {
    assert!(r.i < 256 && r.r < 32768, "rung out of range");
    ((r.i as u32) << 16) | (if r.dir == Dir::F { 1 << 15 } else { 0 }) | r.r as u32
}
/// inverse of `rung_code`
fn rung_decode(c: u32) -> Rung {
    Rung { i: (c >> 16) as usize, dir: if c & (1 << 15) != 0 { Dir::F } else { Dir::E }, r: (c & 0x7fff) as usize }
}
/// lexicographic key of a rung: (i, E < F, r), as Python tuples (i, 'E'|'F', r); the order
/// of `rung_code` (tested)
#[cfg(test)]
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
/// labels after applying a rung, or None if the rung leaves the labels or its giver is too small
fn try_rung_labels(mu: &[u8], r: &Rung) -> Option<Comp> {
    let next = r.i.checked_add(1)?;
    let (rec, giv) = match r.dir {
        Dir::E => (r.i, next),
        Dir::F => (next, r.i),
    };
    if giv >= mu.len() || rec >= mu.len() || (mu[giv] as usize) < r.r || mu[rec] as usize + r.r > u8::MAX as usize {
        return None;
    }
    let mut l = mu.to_vec();
    l[giv] -= r.r as u8;
    l[rec] += r.r as u8;
    Some(l)
}

/// Ok if the block (lambda, mu) can be reached by ladders from 1_lambda: rungs preserve the
/// number of uprights and the total label K, so mu must have the length and the sum of lambda
/// (every such mu is reachable; zero labels are allowed).
pub fn check_block(lam: &[u8], mu: &[u8]) -> Result<(), String> {
    if mu.len() != lam.len() {
        return Err(format!(
            "block {:?} -> {:?}: mu must have as many entries as lambda ({}), since ladders preserve the number of uprights",
            lam,
            mu,
            lam.len()
        ));
    }
    let (a, b) = (lam.iter().map(|&x| x as u32).sum::<u32>(), mu.iter().map(|&x| x as u32).sum::<u32>());
    if a != b {
        return Err(format!("block {:?} -> {:?}: the labels of lambda and mu must have the same sum ({} != {})", lam, mu, a, b));
    }
    Ok(())
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
/// The path is stored as rung codes (`rung_code`), whose numeric order is the lexicographic
/// order (i, E < F, r) of `rung_lex`: 4 bytes per rung instead of 24, and comparisons need no
/// allocation.
struct MfEntry {
    flows: u128,
    path: Vec<u32>,
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
        o.key().cmp(&self.key()).then_with(|| o.path.cmp(&self.path))
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
    /// [`Schur::rung`], or None if a coefficient could overflow i128. Every coefficient of every
    /// partial sum is bounded by |rung(x)|_1 ≤ Σ_p |x_p|_1 · |rung(h_p)|_1 (|·|_1 = sum of the
    /// absolute values of all coefficients), which is checked first.
    fn try_rung(&mut self, x: &Elt, r: &Rung) -> Option<Elt> {
        if r.r == 0 {
            return Some(x.clone());
        }
        let mut bound: u128 = 0;
        for (p, c) in &x.t {
            let key = self.ensure_rung_image(*p, r);
            let img = self.rcache[&key].iter().fold(0u128, |s, (_, k)| s.saturating_add(l1_norm(k)));
            bound = bound.saturating_add(l1_norm(c).saturating_mul(img));
        }
        if bound > i128::MAX as u128 {
            return None;
        }
        Some(self.rung(x, r))
    }
    /// [`Schur::ladder_char`] of an untrusted ladder of valid rungs: None if a coefficient could
    /// overflow i128 (e.g. a long ladder read from a file).
    fn try_ladder_char(&mut self, rs: &[Rung]) -> Option<Elt> {
        let mut y = self.unit();
        for r in rs {
            y = self.try_rung(&y, r)?;
        }
        Some(y)
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
                    let r = rung_decode(*ent.path.last().unwrap());
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
                found.insert(tops[0], ent.path.iter().map(|&c| rung_decode(c)).collect());
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
                        let mut p = Vec::with_capacity(ent.path.len() + 1);
                        p.extend_from_slice(&ent.path);
                        p.push(rung_code(&Rung { i, dir, r }));
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
        if let Err(e) = check_block(&self.lam, mu) {
            panic!("basis_from_ladders: {}", e);
        }
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
        // fail before the BFS, which would otherwise exhaust every reachable state first
        if let Err(e) = check_block(&self.lam, mu) {
            panic!("canonical_basis: {}", e);
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
    /// lazy, cached. Panics (before any search) unless `check_block(lambda, mu)` is Ok.
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
    /// Two-sided cell weight ν(p) of b_p: the RSK shape of the table p with its columns reversed
    /// ([`Table::cell_weight`]).
    ///
    /// **Rule (conjectured, validated by computation).** Under skew Howe duality the q-Schur
    /// algebra S(n, K) of the n uprights acts on Λ^K(C^n ⊗ C^N) = ⊕_ν V_ν(gl_n) ⊗ V_ν'(gl_N),
    /// ν ⊢ K with at most n rows and ν_1 ≤ N; its canonical basis is cellular with two-sided
    /// cells indexed by partitions ν, and b_p annihilates V_μ unless ν(p) ⊴ μ. So the kernel of
    /// the gl_N quotient (the image of S(n, K) in Foam_N) is the span of the b_p with
    /// ν(p)_1 > N, and since the classes of indecomposables are a positive basis, B_p ↦ 0 in
    /// Foam_N (Karoubi completed) exactly for those p ([`Schur::killed_at`]). The convention
    /// "RSK of the column-reversed table" puts the identity diag(λ) into the smallest cell
    /// sort(λ) and the full merge-split (anti-diagonal) into the largest one (K); for the
    /// Rickard terms p_k of a (c, d) crossing it gives (max(c, d) + k, min(c, d) − k), whose
    /// first part is the intermediate label of W_k, so killed_at(p_k, N) is exactly the label
    /// pruning of `crossing::crossing_cx(.., N)`. Hand-checked: diag(3,3) → (3,3), the
    /// anti-diagonal (3,3) table → (6), diag(3,3,3) → (3,3,3), and for (2,2) at N = 2 only the
    /// identity survives. Decategorified check (unit test `kernel_of_gl_n_quotient`): the MOY
    /// trace pairing tr_N(b_p · y) vanishes for every y exactly on the killed p, and the kept b_p
    /// stay linearly independent modulo it, on all blocks of (2,2), (3,3), (2,1), (3,2), (1,1,1),
    /// (2,1,2) and the end blocks of (1,1,1,1), (2,2,2), (3,3,3) (`#[ignore]`: (2,2,2,2)), at
    /// every N from the largest label to K. Categorified check: `xengine::XOpts::prune_n`
    /// (pruned and unpruned scans give the same gl_N homology, tests/ntrunc.rs).
    ///
    /// Use this (Greene: ν(p)_1 = the largest entry sum along a down-left chain of p), never
    /// "the minimal-flow ladder of p has a label > N": the two agree on every 2- and 3-strand
    /// block tested, but not in general. Counterexample (review, 2026-11-10): λ = (3,3,3,1), end
    /// block, p = 1,0,2,0/0,1,1,1/2,1,0,0/0,1,0,0 has cell (7,3) and <b_p, b_{p^T}>_6 = 0
    /// (killed at N = 6), while its minimal-flow ladder never exceeds label 6.
    pub fn cell_weight(&self, p: Tid) -> Comp {
        self.tabs[p as usize].cell_weight()
    }
    /// True if B_p maps to 0 in the gl_N theory: a boundary label (row or column sum) > N, or
    /// ν(p)_1 > N for the cell weight ν(p) = [`Schur::cell_weight`] (the first condition is
    /// implied by the second, since ν(p) dominates both margins; it is kept as a cheap guard).
    pub fn killed_at(&self, p: Tid, n: usize) -> bool {
        let t = &self.tabs[p as usize];
        if t.rows().iter().chain(t.cols().iter()).any(|&a| a as usize > n) {
            return true;
        }
        t.cell_weight().first().map_or(false, |&a| a as usize > n)
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
    /// Load a basis written by `save_basis` (replacing a cached basis of the same block).
    ///
    /// Returns `Err(InvalidData)` (never panics) on malformed input, and verifies the basis:
    /// * the header's lambda equals this instance's, and `check_block(lambda, mu)` holds;
    /// * every table (of a CB, LAD or LADMF line and of every CB entry) lies in the block
    ///   (lambda, mu) (row and column sums compared in u32, so u8 entries cannot wrap around);
    ///   every table of the block has exactly one CB, LAD and LADMF line;
    /// * every exponent of a CB entry lies in [0, L], L the largest lmax of the block (a loose
    ///   bound: the coefficient of h_q in b_p has degree ≤ lmax(p) − lmax(q), as for
    ///   Kazhdan–Lusztig polynomials, checked on computed bases of up to 5 strands); this
    ///   bounds the dense storage of a coefficient before any other check; repeated entries of
    ///   one table are added with overflow checks;
    /// * b_p has coefficient exactly 1 at h_p and every other coefficient is nonzero, in vZ[v]
    ///   and sits at a table q with lmax(q) < lmax(p);
    /// * every ladder is a valid rung sequence on lambda ending at mu, and its character has the
    ///   unique leading term h_p (coefficient 1); a character whose coefficients could
    ///   overflow i128 (e.g. a long ladder) is an error;
    /// * b_p equals the canonical element peeled from the character of its LAD ladder with the
    ///   loaded lower b_q (`peel`, in ascending lmax). By induction on lmax this holds only for the
    ///   true canonical basis, so it subsumes a bar-invariance check: a corrupted b_p is
    ///   rejected even if every expansion stays non-negative. The cost is that of the peeling
    ///   step of `canonical_basis` (no ladder search).
    pub fn load_basis(&mut self, path: &std::path::Path) -> std::io::Result<()> {
        let bad = |m: String| std::io::Error::new(std::io::ErrorKind::InvalidData, format!("{}: {}", path.display(), m));
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
        check_block(&lam, &mu).map_err(bad)?;
        // the tables of the block; every exponent of a b_p lies in [0, largest lmax of the block]
        let block = tables(&lam, &mu);
        let emax = i32::try_from(block.iter().map(|t| t.lmax()).max().unwrap_or(0)).map_err(|_| bad("block too large".into()))?;
        // parse and check that every table lies in the block before interning it (row and
        // column sums in u32: u8 sums could wrap around, e.g. 128 + 130 = 2 mod 256)
        let in_block = |s: &str| -> Result<Table, std::io::Error> {
            let t = Table::try_parse(s).map_err(bad)?;
            let (n, m) = (t.nr as usize, t.nc as usize);
            let sums_ok = n == lam.len()
                && m == mu.len()
                && (0..n).all(|i| t.e[i * m..(i + 1) * m].iter().map(|&x| x as u32).sum::<u32>() == lam[i] as u32)
                && (0..m).all(|j| (0..n).map(|i| t.e[i * m + j] as u32).sum::<u32>() == mu[j] as u32);
            if !sums_ok {
                return Err(bad(format!("table {} is outside the block {:?} -> {:?}", s.trim(), lam, mu)));
            }
            Ok(t)
        };
        let mut raw_b: Vec<(Table, Vec<(Table, Lp)>)> = vec![];
        let mut raw_l: Vec<(bool, Table, Vec<Rung>)> = vec![];
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
            let pt = in_block(f[1])?;
            match tag {
                "CB" => {
                    let mut ents = vec![];
                    for ent in &f[2..] {
                        let ent = ent.trim();
                        if ent.is_empty() {
                            continue;
                        }
                        let (tq, poly) = ent.split_once(' ').ok_or_else(|| bad(format!("bad entry {:?}", ent)))?;
                        ents.push((in_block(tq)?, Lp::try_parse_ec_in(poly, 0, emax).map_err(bad)?));
                    }
                    raw_b.push((pt, ents));
                }
                "LAD" | "LADMF" => {
                    if f.len() != 3 {
                        return Err(bad(format!("bad ladder line {:?}", line)));
                    }
                    let w: Option<Vec<Rung>> = f[2].split_whitespace().map(parse_rung).collect();
                    let w = w.ok_or_else(|| bad(format!("bad ladder {:?}", line)))?;
                    let mut lab = lam.clone();
                    for r in &w {
                        lab = try_rung_labels(&lab, r).ok_or_else(|| bad(format!("invalid rung {} in {:?}", rung_str(r), line)))?;
                    }
                    if lab != mu {
                        return Err(bad(format!("ladder ends at {:?}, not at mu: {:?}", lab, line)));
                    }
                    raw_l.push((tag == "LAD", pt, w));
                }
                _ => return Err(bad(format!("unknown tag {:?}", tag))),
            }
        }
        let mut b: HashMap<Tid, Elt> = HashMap::new();
        let mut lb: HashMap<Tid, Vec<Rung>> = HashMap::new();
        let mut lm: HashMap<Tid, Vec<Rung>> = HashMap::new();
        for (pt, ents) in raw_b {
            let p = self.intern(pt);
            let mut acc: FxMap<Tid, Lp> = FxMap::default();
            for (qt, c) in ents {
                let q = self.intern(qt);
                let slot = acc.entry(q).or_default();
                *slot = slot
                    .checked_add(&c)
                    .ok_or_else(|| bad(format!("b_{}: coefficient overflow in the repeated entries of h_{}", self.tabs[p as usize].fmt(), self.tabs[q as usize].fmt())))?;
            }
            if b.insert(p, finish_acc(mu.clone(), acc)).is_some() {
                return Err(bad(format!("two CB lines for table {}", self.tabs[p as usize].fmt())));
            }
        }
        for (is_bfs, pt, w) in raw_l {
            let p = self.intern(pt);
            let dup = if is_bfs { lb.insert(p, w) } else { lm.insert(p, w) };
            if dup.is_some() {
                return Err(bad(format!("two {} lines for table {}", if is_bfs { "LAD" } else { "LADMF" }, self.tabs[p as usize].fmt())));
            }
        }
        let order: Vec<Tid> = block.into_iter().map(|t| self.intern(t)).collect();
        for p in &order {
            if !b.contains_key(p) || !lb.contains_key(p) || !lm.contains_key(p) {
                return Err(bad(format!("basis file misses table {}", self.tabs[*p as usize].fmt())));
            }
        }
        if b.len() != order.len() || lb.len() != order.len() || lm.len() != order.len() {
            return Err(bad("basis file has tables outside the block".into()));
        }
        // b_p = h_p + sum_{lmax(q) < lmax(p)} vZ[v] h_q
        for &p in &order {
            let lp = self.info[p as usize].lmax;
            for (q, c) in &b[&p].t {
                if *q == p {
                    if *c != Lp::one() {
                        return Err(bad(format!("b_{}: coefficient of h_p is {}, not 1", self.tabs[p as usize].fmt(), c.fmt_ec())));
                    }
                } else if c.min_exp().map_or(true, |e| e <= 0) {
                    return Err(bad(format!("b_{}: coefficient {} of h_{} is not in vZ[v]", self.tabs[p as usize].fmt(), c.fmt_ec(), self.tabs[*q as usize].fmt())));
                } else if self.info[*q as usize].lmax >= lp {
                    return Err(bad(format!("b_{}: term h_{} is not below p (lmax)", self.tabs[p as usize].fmt(), self.tabs[*q as usize].fmt())));
                }
            }
            if b[&p].coef(p) != Lp::one() {
                return Err(bad(format!("b_{}: coefficient of h_p is not 1", self.tabs[p as usize].fmt())));
            }
        }
        // ladders lead to h_p; b_p is the canonical element peeled from its LAD ladder
        for &p in &order {
            let pf = self.tabs[p as usize].fmt();
            for (tag, lads) in [("LADMF", &lm), ("LAD", &lb)] {
                let x = self.try_ladder_char(&lads[&p]).ok_or_else(|| bad(format!("{} ladder of {}: its character overflows i128 coefficients", tag, pf)))?;
                let (tops, clean) = self.leading(&x);
                if !clean || tops != [p] {
                    return Err(bad(format!("{} ladder of {} does not have the leading term h_p", tag, pf)));
                }
                if tag == "LAD" && self.peel(p, &x, &b) != b[&p] {
                    return Err(bad(format!("b_{} is not the canonical basis element (it differs from the one peeled from its ladder)", pf)));
                }
            }
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

const CNT: &str = "schur: u128 count overflow";

/// acc + a * b with checked u128 arithmetic (release builds do not check overflow)
#[inline]
fn add_mul_u128(acc: u128, a: u128, b: u128) -> u128 {
    acc.checked_add(a.checked_mul(b).expect(CNT)).expect(CNT)
}

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
                    *ent = ent.checked_add(sg.checked_mul(c).expect(OVF)).expect(OVF);
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
        lb_objects = add_mul_u128(lb_objects, a, 1);
        lb_ldim = add_mul_u128(lb_ldim, a, s.rank(labels, p));
        lb_sum_s = add_mul_u128(lb_sum_s, a, sself[&p]);
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
                            x1 = add_mul_u128(x1, k, 1);
                            x2 = add_mul_u128(x2, k, s.rank(&nl, *q));
                            x3 = add_mul_u128(x3, k, s.s_storage(&nl, *q, *q));
                        }
                        pcache.insert(key, (x1, x2, x3));
                        (x1, x2, x3)
                    }
                };
                a = add_mul_u128(a, mult, tr.0);
                b = add_mul_u128(b, mult, tr.1);
                cs = add_mul_u128(cs, mult, tr.2);
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

/// Checks the input of [`decat_scan`]: every generator of the word acts on two existing
/// strands, and the marked block (if any) is a block of lambda (0-based) with a nonzero label.
/// `decat_scan` panics with this message; front ends call it first to report a usage error.
pub fn check_scan(lam: &[u8], word: &[i32], opts: &ScanOpts) -> Result<(), String> {
    for (k, &g) in word.iter().enumerate() {
        if g == 0 {
            return Err(format!("word position {}: braid generator 0", k + 1));
        }
        let i = g.unsigned_abs() as usize - 1;
        if i + 1 >= lam.len() {
            return Err(format!("word position {}: generator {} needs {} strands, lambda {:?} has {}", k + 1, g, i + 2, lam, lam.len()));
        }
    }
    if let Some(mb) = opts.marked_block {
        if mb >= lam.len() {
            return Err(format!("marked block {} out of range: lambda {:?} has blocks 0..{}", mb, lam, lam.len()));
        }
        if lam[mb] == 0 {
            return Err(format!("marked block {} of lambda {:?} has label 0", mb, lam));
        }
    }
    Ok(())
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
/// * Panics unless [`check_scan`] accepts the input (e.g. a marked block out of range).
pub fn decat_scan(s: &mut Schur, word: &[i32], opts: &ScanOpts) -> ScanResult {
    if let Err(e) = check_scan(&s.lam, word, opts) {
        panic!("decat_scan: {}", e);
    }
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stage_counts_are_checked() {
        assert_eq!(add_mul_u128(5, 3, 4), 17);
        assert_eq!(add_mul_u128(u128::MAX - 1, 1, 1), u128::MAX);
        // release builds do not check overflow: the stage counters must
        assert!(std::panic::catch_unwind(|| add_mul_u128(u128::MAX, 1, 1)).is_err());
        assert!(std::panic::catch_unwind(|| add_mul_u128(0, u128::MAX, 2)).is_err());
    }

    // ----------------------------------------------------------------- cell weights (Foam_N pruning)

    #[test]
    fn cell_weight_hand_checked() {
        let w = |s: &str| Table::parse(s).cell_weight();
        assert_eq!(w("3,0/0,3"), vec![3, 3]); // diag(3,3): the identity, smallest cell
        assert_eq!(w("0,3/3,0"), vec![6]); // anti-diagonal: the full merge-split
        assert_eq!(w("3,0,0/0,3,0/0,0,3"), vec![3, 3, 3]);
        assert_eq!(w("2,1/1,2"), vec![4, 2]);
        assert_eq!(w("1,2/2,1"), vec![5, 1]);
        assert_eq!(w("1,1/1,1"), vec![3, 1]);
        // permutation matrices: the transpose of the RSK shape of w (identity -> (1^n), w0 -> (n))
        assert_eq!(w("1,0,0/0,1,0/0,0,1"), vec![1, 1, 1]);
        assert_eq!(w("0,0,1/0,1,0/1,0,0"), vec![3]);
        assert_eq!(w("0,1,0/1,0,0/0,0,1"), vec![2, 1]);
        // reversing rows instead of columns, or transposing, gives the same shape
        for s in ["2,1,0/0,1,2/1,0,1", "1,0,2/2,1,0", "0,2/1,1/2,0", "3,1,0,1/0,2,2,0"] {
            let t = Table::parse(s);
            let (n, m) = (t.nr as usize, t.nc as usize);
            let rev_rows = Table { nr: t.nr, nc: t.nc, e: (0..n).flat_map(|i| t.e[(n - 1 - i) * m..(n - i) * m].to_vec()).collect() };
            let rev_both = Table { nr: t.nr, nc: t.nc, e: rev_rows.e.chunks(m).flat_map(|r| r.iter().rev().copied().collect::<Vec<_>>()).collect() };
            let tr = Table { nr: t.nc, nc: t.nr, e: (0..m).flat_map(|j| (0..n).map(move |i| (i, j))).map(|(i, j)| t.e[i * m + j]).collect() };
            // column reversal of rev_both = row reversal of t
            assert_eq!(rev_both.cell_weight(), t.cell_weight(), "{}", s);
            assert_eq!(tr.cell_weight(), t.cell_weight(), "{} transposed", s);
            let k: u32 = t.e.iter().map(|&a| a as u32).sum();
            assert_eq!(t.cell_weight().iter().map(|&a| a as u32).sum::<u32>(), k);
        }
    }

    #[test]
    fn killed_at_examples() {
        // (2,2) at N = 2: only the identity type survives (M6.0's 2-colored gl_2 results are q^0)
        let mut s = Schur::new(&[2, 2]);
        let ids: Vec<Tid> = tables(&[2, 2], &[2, 2]).into_iter().map(|t| s.intern(t)).collect();
        let alive: Vec<String> = ids.iter().filter(|&&p| !s.killed_at(p, 2)).map(|&p| s.table(p).fmt()).collect();
        assert_eq!(alive, vec!["2,0/0,2".to_string()]);
        assert_eq!(ids.iter().filter(|&&p| !s.killed_at(p, 3)).count(), 2); // (2,2), (3,1)
        assert_eq!(ids.iter().filter(|&&p| !s.killed_at(p, 4)).count(), 3);
        // a boundary label > N kills
        let t = s.intern(Table::parse("2,0/2,0"));
        assert!(s.killed_at(t, 3) && !s.killed_at(t, 4));
        // (3,3) at N = 3, 4, 5, 6: 1, 2, 3, 4 surviving types
        let mut s = Schur::new(&[3, 3]);
        let ids: Vec<Tid> = tables(&[3, 3], &[3, 3]).into_iter().map(|t| s.intern(t)).collect();
        for (n, alive) in [(3usize, 1usize), (4, 2), (5, 3), (6, 4)] {
            assert_eq!(ids.iter().filter(|&&p| !s.killed_at(p, n)).count(), alive, "N={}", n);
        }
        // Rickard terms: killed_at(p_k, N) iff the web W_k has a label > N (crossing_cx's pruning)
        for c in 1..=4u8 {
            for d in 1..=4u8 {
                for positive in [true, false] {
                    let mut s = Schur::new(&[c, d]);
                    for t in crossing_terms(c, d, positive, 0) {
                        let ch = s.ladder_char(&t.rungs);
                        let e = s.expand(&ch);
                        assert!(e.len() == 1 && e[0].1 == Lp::one(), "Rickard term not indecomposable");
                        let p = e[0].0;
                        let mut lab = vec![c as usize, d as usize];
                        let mut top = c.max(d) as usize;
                        for r in &t.rungs {
                            lab = r.apply(&lab, 64).unwrap();
                            top = top.max(lab[0]).max(lab[1]);
                        }
                        assert_eq!(s.cell_weight(p)[0] as usize, top, "({},{}) {} {:?}", c, d, positive, t.rungs);
                        for n in c.max(d) as usize..=(c + d) as usize {
                            assert_eq!(s.killed_at(p, n), top > n);
                        }
                    }
                }
            }
        }
    }

    const PR: u64 = 2_147_483_647;

    fn pow_mod(mut b: u64, mut e: u64) -> u64 {
        let mut r = 1u64;
        b %= PR;
        while e > 0 {
            if e & 1 == 1 {
                r = r * b % PR;
            }
            b = b * b % PR;
            e >>= 1;
        }
        r
    }

    /// value of a Laurent polynomial at v = x, mod PR
    fn eval_mod(p: &LPoly, x: u64) -> u64 {
        let xi = pow_mod(x, PR - 2);
        let mut acc = 0u64;
        for (&e, &k) in &p.c {
            let pw = if e >= 0 { pow_mod(x, e as u64) } else { pow_mod(xi, (-e) as u64) };
            let km = k.rem_euclid(PR as i128) as u64;
            acc = (acc + km * pw % PR) % PR;
        }
        acc
    }

    fn rank_mod(mut rows: Vec<Vec<u64>>) -> usize {
        let mut rank = 0;
        let ncols = rows.first().map_or(0, |r| r.len());
        for c in 0..ncols {
            let Some(piv) = (rank..rows.len()).find(|&i| rows[i][c] != 0) else { continue };
            rows.swap(rank, piv);
            let inv = pow_mod(rows[rank][c], PR - 2);
            for i in 0..rows.len() {
                if i != rank && rows[i][c] != 0 {
                    let f = rows[i][c] * inv % PR;
                    for j in c..ncols {
                        let sub = f * rows[rank][j] % PR;
                        rows[i][j] = (rows[i][j] + PR - sub) % PR;
                    }
                }
            }
            rank += 1;
        }
        rank
    }

    fn compositions(k: u8, parts: usize) -> Vec<Comp> {
        if parts == 1 {
            return vec![vec![k]];
        }
        let mut out = vec![];
        for a in 0..=k {
            for mut rest in compositions(k - a, parts - 1) {
                rest.insert(0, a);
                out.push(rest);
            }
        }
        out
    }

    /// Decategorified check of the cell rule on the block (lam, mu) at gl_N: the trace pairing
    /// T[p][q] = tr_N(b_p · ch(L_q^rev)) (L_q^rev = the minimal-flow ladder of q read top to
    /// bottom, a ladder mu -> lam; these characters span 1_mu S 1_lam) vanishes identically on
    /// the killed p, and the kept rows have full rank (evaluated at v = 3 and v = 5 mod p).
    /// Since the MOY closure is the gl_N quantum trace and the trace form is nondegenerate on
    /// the semisimple gl_N quotient, this says: the kernel of the gl_N quotient on the block is
    /// exactly the span of the killed b_p. Returns (killed, kept, killed types whose canonical
    /// ladder never leaves the labels <= N).
    fn kernel_check(s: &mut Schur, mu: &[u8], n: usize, memo: &mut FxMap<Tid, LPoly>) -> (usize, usize, usize) {
        let lam = s.lam.clone();
        let ids: Vec<Tid> = tables(&lam, mu).into_iter().map(|t| s.intern(t)).collect();
        let cb = s.canonical_basis(mu);
        let bs: Vec<Elt> = ids.iter().map(|p| cb.b[p].clone()).collect();
        let lads: Vec<Vec<Rung>> = ids.iter().map(|p| cb.ladder_minflow[p].clone()).collect();
        let revs: Vec<Vec<Rung>> = lads.iter().map(|l| l.iter().rev().map(|r| Rung { i: r.i, dir: if r.dir == Dir::E { Dir::F } else { Dir::E }, r: r.r }).collect()).collect();
        let (mut killed, mut kept, mut invisible) = (0, 0, 0);
        let mut kept_rows: Vec<Vec<LPoly>> = vec![];
        for (ip, &p) in ids.iter().enumerate() {
            let mut row = vec![];
            for y in &revs {
                let x = s.mul_rungs(&bs[ip], y);
                let mut v = LPoly::zero();
                for (r, m) in s.expand(&x) {
                    let t = trace_b(s, r, n, false, memo);
                    v = v.add(&m.to_lpoly().mul(&t));
                }
                row.push(v);
            }
            // does the canonical (minimal-flow) ladder of p stay within the labels <= N?
            let mut lab: Vec<usize> = lam.iter().map(|&a| a as usize).collect();
            let mut within = lab.iter().all(|&a| a <= n);
            for r in &lads[ip] {
                lab = r.apply(&lab, 64).unwrap();
                within &= lab.iter().all(|&a| a <= n);
            }
            if s.killed_at(p, n) {
                killed += 1;
                assert!(row.iter().all(|v| v.is_zero()), "{:?} -> {:?} N={}: killed type {} (cell {:?}) pairs nontrivially", lam, mu, n, s.table(p).fmt(), s.cell_weight(p));
                if within {
                    invisible += 1;
                }
            } else {
                // B_p is a summand of its ladder, which is 0 in Foam_N if it has a label > N
                assert!(within, "{:?} -> {:?} N={}: kept type {} (cell {:?}) has a ladder through a label > N", lam, mu, n, s.table(p).fmt(), s.cell_weight(p));
                kept += 1;
                kept_rows.push(row);
            }
        }
        let rk = [3u64, 5].iter().map(|&x| rank_mod(kept_rows.iter().map(|r| r.iter().map(|v| eval_mod(v, x)).collect()).collect())).max().unwrap_or(0);
        assert_eq!(rk, kept, "{:?} -> {:?} N={}: the kept b_p are dependent modulo the gl_N kernel (rank {} < {})", lam, mu, n, rk, kept);
        (killed, kept, invisible)
    }

    #[test]
    fn kernel_of_gl_n_quotient() {
        let t0 = std::time::Instant::now();
        // (lambda, N range, all blocks mu or only mu = lambda)
        let cases: Vec<(Vec<u8>, std::ops::RangeInclusive<usize>, bool)> = vec![
            (vec![2, 2], 2..=4, true),
            (vec![3, 3], 3..=6, true),
            (vec![2, 1], 2..=3, true),
            (vec![3, 2], 3..=5, true),
            (vec![1, 1, 1], 1..=3, true),
            (vec![2, 1, 2], 2..=5, true),
            (vec![1, 1, 1, 1], 1..=4, false),
            (vec![2, 2, 2], 2..=6, false),
            (vec![3, 3, 3], 3..=6, false),
        ];
        let (mut nk, mut nalive, mut ninv, mut nblocks) = (0, 0, 0, 0);
        let mut report = vec![];
        for (lam, ns, all_mu) in cases {
            let k: u8 = lam.iter().sum();
            let mus = if all_mu { compositions(k, lam.len()) } else { vec![lam.clone()] };
            for n in ns {
                let mut s = Schur::new(&lam);
                let mut memo: FxMap<Tid, LPoly> = FxMap::default();
                let (mut a, mut b, mut c) = (0, 0, 0);
                for mu in &mus {
                    let (x, y, z) = kernel_check(&mut s, mu, n, &mut memo);
                    a += x;
                    b += y;
                    c += z;
                    nblocks += 1;
                }
                report.push(format!("{:?} N={}: {} killed, {} kept, {} killed with all ladder labels <= N", lam, n, a, b, c));
                nk += a;
                nalive += b;
                ninv += c;
            }
        }
        for l in &report {
            eprintln!("{}", l);
        }
        eprintln!("kernel_of_gl_n_quotient: {} blocks, {} killed, {} kept, {} invisible killed, {:.1}s", nblocks, nk, nalive, ninv, t0.elapsed().as_secs_f64());
        assert!(nk > 0 && nalive > 0);
    }

    #[test]
    #[ignore]
    fn kernel_of_gl_n_quotient_2222() {
        // (2,2,2,2) -> (2,2,2,2) (282 types) at N = 2..5; slow (282^2 products per N)
        let t0 = std::time::Instant::now();
        for n in 2..=5usize {
            let mut s = Schur::new(&[2, 2, 2, 2]);
            let mut memo: FxMap<Tid, LPoly> = FxMap::default();
            let (a, b, c) = kernel_check(&mut s, &[2, 2, 2, 2], n, &mut memo);
            eprintln!("[2, 2, 2, 2] N={}: {} killed, {} kept, {} killed with all ladder labels <= N ({:.1}s)", n, a, b, c, t0.elapsed().as_secs_f64());
        }
    }

    #[test]
    fn rung_codes_order_like_python_tuples() {
        // the minimal-flow heap compares paths of rung codes: same order as (i, 'E'|'F', r)
        let mut rs = vec![];
        for i in [0usize, 1, 2, 7, 255] {
            for dir in [Dir::E, Dir::F] {
                for r in [1usize, 2, 3, 9, 300] {
                    rs.push(Rung { i, dir, r });
                }
            }
        }
        for a in &rs {
            assert_eq!(rung_decode(rung_code(a)), *a);
            for b in &rs {
                assert_eq!(rung_code(a).cmp(&rung_code(b)), rung_lex(a).cmp(&rung_lex(b)), "{:?} {:?}", a, b);
            }
        }
    }
}
