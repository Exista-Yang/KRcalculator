//! Multivariate polynomials over a field with packed monomials (8 bits per
//! variable, at most 16 variables), graded by per-variable weights.

use crate::field::Field;
use std::collections::HashMap;
use std::sync::Arc;

pub type Mono = u128;

pub const MAXVARS: usize = 16;

#[inline]
pub fn mono_exp(m: Mono, i: usize) -> u32 {
    ((m >> (8 * i)) & 0xff) as u32
}
#[inline]
pub fn var_mono(i: usize) -> Mono {
    1u128 << (8 * i)
}
#[inline]
pub fn mono_mul(a: Mono, b: Mono) -> Mono {
    // caller guarantees no byte overflow (exponents < 256)
    a + b
}

/// Polynomial ring data: variable weights (q-degree of each variable).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PRing {
    pub nvars: usize,
    pub w: Vec<i32>,
    pub names: Vec<String>,
}

impl PRing {
    pub fn new(w: Vec<i32>, names: Vec<String>) -> Arc<PRing> {
        assert!(w.len() <= MAXVARS, "too many variables");
        Arc::new(PRing { nvars: w.len(), w, names })
    }
    #[inline]
    pub fn deg(&self, m: Mono) -> i32 {
        let mut d = 0;
        let mut x = m;
        let mut i = 0;
        while x != 0 {
            d += (x & 0xff) as i32 * self.w[i];
            x >>= 8;
            i += 1;
        }
        d
    }
    /// All monomials of exactly degree d (cached by caller).
    pub fn monos_of_deg(&self, d: i32) -> Vec<Mono> {
        let mut out = vec![];
        if d < 0 {
            return out;
        }
        fn rec(w: &[i32], i: usize, left: i32, cur: Mono, out: &mut Vec<Mono>) {
            if i == w.len() {
                if left == 0 {
                    out.push(cur);
                }
                return;
            }
            let mut k = 0i32;
            while k * w[i] <= left {
                rec(w, i + 1, left - k * w[i], cur + (k as u128) * var_mono(i), out);
                k += 1;
                if w[i] == 0 {
                    break;
                }
            }
        }
        rec(&self.w, 0, d, 0, &mut out);
        out.sort_unstable();
        out
    }
}

/// Cache of monomial bases per degree with index maps.
pub struct MonoCache {
    pub ring: Arc<PRing>,
    pub by_deg: HashMap<i32, (Vec<Mono>, HashMap<Mono, u32>)>,
}

impl MonoCache {
    pub fn new(ring: Arc<PRing>) -> Self {
        MonoCache { ring, by_deg: HashMap::new() }
    }
    pub fn get(&mut self, d: i32) -> &(Vec<Mono>, HashMap<Mono, u32>) {
        let ring = self.ring.clone();
        self.by_deg.entry(d).or_insert_with(|| {
            let v = ring.monos_of_deg(d);
            let m = v.iter().enumerate().map(|(i, &x)| (x, i as u32)).collect();
            (v, m)
        })
    }
}

/// Sparse polynomial: terms sorted by monomial, nonzero coefficients.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct Poly<F: Field> {
    pub t: Vec<(Mono, F)>,
}

impl<F: Field> Poly<F> {
    pub fn zero() -> Self {
        Poly { t: vec![] }
    }
    pub fn constant(c: F) -> Self {
        if c.is_zero() { Self::zero() } else { Poly { t: vec![(0, c)] } }
    }
    pub fn one() -> Self {
        Self::constant(F::one())
    }
    pub fn mono(m: Mono, c: F) -> Self {
        if c.is_zero() { Self::zero() } else { Poly { t: vec![(m, c)] } }
    }
    pub fn var(i: usize) -> Self {
        Self::mono(var_mono(i), F::one())
    }
    pub fn is_zero(&self) -> bool {
        self.t.is_empty()
    }
    fn normalize(mut v: Vec<(Mono, F)>) -> Self {
        v.sort_unstable_by_key(|x| x.0);
        let mut out: Vec<(Mono, F)> = Vec::with_capacity(v.len());
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
        Poly { t: out }
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
        Poly { t: out }
    }
    pub fn add_assign(&mut self, o: &Self) {
        if o.is_zero() {
            return;
        }
        *self = self.add(o);
    }
    pub fn neg(&self) -> Self {
        Poly { t: self.t.iter().map(|(m, c)| (*m, c.neg())).collect() }
    }
    pub fn sub(&self, o: &Self) -> Self {
        self.add(&o.neg())
    }
    pub fn scale(&self, c: &F) -> Self {
        if c.is_zero() {
            return Self::zero();
        }
        Poly { t: self.t.iter().map(|(m, x)| (*m, x.mul(c))).collect() }
    }
    pub fn mul_mono(&self, m: Mono, c: &F) -> Self {
        if c.is_zero() {
            return Self::zero();
        }
        Poly { t: self.t.iter().map(|(a, x)| (a + m, x.mul(c))).collect() }
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
        let mut h: HashMap<Mono, F> = HashMap::with_capacity(self.t.len() * o.t.len());
        for (a, x) in &self.t {
            for (b, y) in &o.t {
                h.entry(a + b).or_insert_with(F::zero).fma(x, y);
            }
        }
        Self::normalize(h.into_iter().collect())
    }
    /// Coefficient of a monomial.
    pub fn coeff(&self, m: Mono) -> F {
        match self.t.binary_search_by_key(&m, |x| x.0) {
            Ok(i) => self.t[i].1.clone(),
            Err(_) => F::zero(),
        }
    }
}

pub fn fmt_poly<F: Field>(p: &Poly<F>, ring: &PRing) -> String {
    if p.is_zero() {
        return "0".into();
    }
    let mut parts = vec![];
    for (m, c) in &p.t {
        let mut s = format!("{}", c);
        for i in 0..ring.nvars {
            let e = mono_exp(*m, i);
            if e > 0 {
                s += &format!("*{}", ring.names[i]);
                if e > 1 {
                    s += &format!("^{}", e);
                }
            }
        }
        parts.push(s);
    }
    parts.join(" + ")
}
