//! Small multivariate polynomials with i64 coefficients (used for precomputing
//! universal structure constants; not performance critical).

use std::collections::HashMap;

pub type Mono = Vec<u16>;

#[derive(Clone, Debug, PartialEq)]
pub struct MPoly {
    pub nvars: usize,
    pub t: HashMap<Mono, i64>,
}

impl MPoly {
    pub fn zero(nvars: usize) -> Self {
        MPoly { nvars, t: HashMap::new() }
    }
    pub fn constant(nvars: usize, c: i64) -> Self {
        let mut p = Self::zero(nvars);
        if c != 0 {
            p.t.insert(vec![0; nvars], c);
        }
        p
    }
    pub fn one(nvars: usize) -> Self {
        Self::constant(nvars, 1)
    }
    pub fn var(nvars: usize, i: usize) -> Self {
        let mut m = vec![0; nvars];
        m[i] = 1;
        Self::monomial(m)
    }
    pub fn monomial(m: Mono) -> Self {
        let nvars = m.len();
        let mut p = Self::zero(nvars);
        p.t.insert(m, 1);
        p
    }
    pub fn add_term(&mut self, m: Mono, c: i64) {
        if c == 0 {
            return;
        }
        let e = self.t.entry(m.clone()).or_insert(0);
        *e = e.checked_add(c).expect("MPoly overflow");
        if *e == 0 {
            self.t.remove(&m);
        }
    }
    pub fn add(&self, o: &MPoly) -> MPoly {
        let mut r = self.clone();
        for (m, &c) in &o.t {
            r.add_term(m.clone(), c);
        }
        r
    }
    pub fn mul(&self, o: &MPoly) -> MPoly {
        let mut r = MPoly::zero(self.nvars);
        for (m1, &c1) in &self.t {
            for (m2, &c2) in &o.t {
                let m: Mono = m1.iter().zip(m2.iter()).map(|(a, b)| a + b).collect();
                r.add_term(m, c1.checked_mul(c2).expect("MPoly overflow"));
            }
        }
        r
    }
    pub fn pow(&self, k: u32) -> MPoly {
        let mut r = MPoly::one(self.nvars);
        for _ in 0..k {
            r = r.mul(self);
        }
        r
    }
    pub fn is_zero(&self) -> bool {
        self.t.is_empty()
    }
}

/// All exponent vectors m in N^w.len() with sum_i w[i]*m[i] == d.
pub fn monomials_of_weight(w: &[u32], d: u32) -> Vec<Mono> {
    let mut out = Vec::new();
    let mut cur = vec![0u16; w.len()];
    fn rec(w: &[u32], i: usize, left: u32, cur: &mut Vec<u16>, out: &mut Vec<Mono>) {
        if i == w.len() {
            if left == 0 {
                out.push(cur.clone());
            }
            return;
        }
        let mut k = 0u32;
        while k * w[i] <= left {
            cur[i] = k as u16;
            rec(w, i + 1, left - k * w[i], cur, out);
            k += 1;
        }
        cur[i] = 0;
    }
    rec(w, 0, d, &mut cur, &mut out);
    out
}
