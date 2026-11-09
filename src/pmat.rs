//! Sparse matrices with polynomial entries (column-major).

use crate::field::Field;
use crate::mpoly::{Mono, Poly};
use std::collections::HashMap;

pub type PVec<F> = Vec<(u32, Poly<F>)>;

#[derive(Clone, Debug)]
pub struct PMat<F: Field> {
    pub nrows: usize,
    pub ncols: usize,
    pub cols: Vec<PVec<F>>,
}

pub fn pvec_add<F: Field>(a: &PVec<F>, b: &PVec<F>) -> PVec<F> {
    let mut out = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() || j < b.len() {
        if j == b.len() || (i < a.len() && a[i].0 < b[j].0) {
            out.push(a[i].clone());
            i += 1;
        } else if i == a.len() || b[j].0 < a[i].0 {
            out.push(b[j].clone());
            j += 1;
        } else {
            let p = a[i].1.add(&b[j].1);
            if !p.is_zero() {
                out.push((a[i].0, p));
            }
            i += 1;
            j += 1;
        }
    }
    out
}

pub fn pvec_scale<F: Field>(a: &PVec<F>, p: &Poly<F>) -> PVec<F> {
    a.iter().map(|(i, x)| (*i, x.mul(p))).filter(|(_, x)| !x.is_zero()).collect()
}

/// Accumulate sum of polynomial-vectors efficiently.
pub struct PAcc<F: Field> {
    h: HashMap<u32, HashMap<Mono, F>>,
}
impl<F: Field> PAcc<F> {
    pub fn new() -> Self {
        PAcc { h: HashMap::new() }
    }
    pub fn add_scaled(&mut self, v: &PVec<F>, p: &Poly<F>) {
        for (i, x) in v {
            let e = self.h.entry(*i).or_default();
            for (a, c) in &x.t {
                for (b, d) in &p.t {
                    e.entry(a + b).or_insert_with(F::zero).fma(c, d);
                }
            }
        }
    }
    pub fn add_scaled_const(&mut self, v: &PVec<F>, k: &F) {
        for (i, x) in v {
            let e = self.h.entry(*i).or_default();
            for (a, c) in &x.t {
                e.entry(*a).or_insert_with(F::zero).fma(c, k);
            }
        }
    }
    pub fn take(&mut self) -> PVec<F> {
        let mut out: PVec<F> = Vec::with_capacity(self.h.len());
        for (i, m) in self.h.drain() {
            let mut t: Vec<(Mono, F)> = m.into_iter().filter(|(_, c)| !c.is_zero()).collect();
            if t.is_empty() {
                continue;
            }
            t.sort_unstable_by_key(|x| x.0);
            out.push((i, Poly { t }));
        }
        out.sort_unstable_by_key(|x| x.0);
        out
    }
}

impl<F: Field> PMat<F> {
    pub fn zero(nrows: usize, ncols: usize) -> Self {
        PMat { nrows, ncols, cols: vec![vec![]; ncols] }
    }
    pub fn identity(n: usize) -> Self {
        PMat { nrows: n, ncols: n, cols: (0..n).map(|i| vec![(i as u32, Poly::one())]).collect() }
    }
    pub fn scalar(n: usize, p: &Poly<F>) -> Self {
        if p.is_zero() {
            return Self::zero(n, n);
        }
        PMat { nrows: n, ncols: n, cols: (0..n).map(|i| vec![(i as u32, p.clone())]).collect() }
    }
    pub fn is_zero(&self) -> bool {
        self.cols.iter().all(|c| c.is_empty())
    }
    pub fn mul_vec(&self, v: &PVec<F>) -> PVec<F> {
        let mut acc = PAcc::new();
        for (j, p) in v {
            acc.add_scaled(&self.cols[*j as usize], p);
        }
        acc.take()
    }
    pub fn matmul(&self, o: &PMat<F>) -> PMat<F> {
        assert_eq!(self.ncols, o.nrows);
        PMat { nrows: self.nrows, ncols: o.ncols, cols: o.cols.iter().map(|c| self.mul_vec(c)).collect() }
    }
    pub fn add(&self, o: &PMat<F>) -> PMat<F> {
        assert_eq!((self.nrows, self.ncols), (o.nrows, o.ncols));
        PMat { nrows: self.nrows, ncols: self.ncols, cols: self.cols.iter().zip(o.cols.iter()).map(|(a, b)| pvec_add(a, b)).collect() }
    }
    pub fn scale(&self, p: &Poly<F>) -> PMat<F> {
        PMat { nrows: self.nrows, ncols: self.ncols, cols: self.cols.iter().map(|c| pvec_scale(c, p)).collect() }
    }
    pub fn neg(&self) -> PMat<F> {
        self.scale(&Poly::constant(F::one().neg()))
    }
    pub fn sub(&self, o: &PMat<F>) -> PMat<F> {
        self.add(&o.neg())
    }
    /// self ⊗ I_k with index convention row*k + b.
    pub fn kron_id(&self, k: usize) -> PMat<F> {
        let mut cols = Vec::with_capacity(self.ncols * k);
        for c in &self.cols {
            for b in 0..k {
                cols.push(c.iter().map(|(i, p)| (*i * k as u32 + b as u32, p.clone())).collect());
            }
        }
        PMat { nrows: self.nrows * k, ncols: self.ncols * k, cols }
    }
    pub fn transpose(&self) -> PMat<F> {
        let mut cols: Vec<PVec<F>> = vec![vec![]; self.nrows];
        for (j, c) in self.cols.iter().enumerate() {
            for (i, p) in c {
                cols[*i as usize].push((j as u32, p.clone()));
            }
        }
        PMat { nrows: self.ncols, ncols: self.nrows, cols }
    }
    pub fn get(&self, r: usize, c: usize) -> Poly<F> {
        match self.cols[c].binary_search_by_key(&(r as u32), |x| x.0) {
            Ok(i) => self.cols[c][i].1.clone(),
            Err(_) => Poly::zero(),
        }
    }
}
