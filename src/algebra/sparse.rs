//! Sparse vectors and column-major sparse matrices over a field.

use crate::algebra::field::Field;

/// Sparse vector: sorted (index, nonzero value) pairs.
pub type SpVec<F> = Vec<(u32, F)>;

/// Sparse accumulator for building linear combinations of sparse vectors.
pub struct Acc<F: Field> {
    vals: Vec<F>,
    used: Vec<bool>,
    idx: Vec<u32>,
}

impl<F: Field> Acc<F> {
    pub fn new(n: usize) -> Self {
        Acc { vals: vec![F::zero(); n], used: vec![false; n], idx: Vec::new() }
    }
    pub fn len(&self) -> usize {
        self.vals.len()
    }
    #[inline]
    pub fn add(&mut self, i: u32, v: &F) {
        let iu = i as usize;
        if !self.used[iu] {
            self.used[iu] = true;
            self.idx.push(i);
            self.vals[iu] = v.clone();
        } else {
            self.vals[iu].add_assign(v);
        }
    }
    /// self += c * v
    #[inline]
    pub fn axpy(&mut self, c: &F, v: &SpVec<F>) {
        for (i, x) in v {
            let iu = *i as usize;
            if !self.used[iu] {
                self.used[iu] = true;
                self.idx.push(*i);
                self.vals[iu] = c.mul(x);
            } else {
                self.vals[iu].fma(c, x);
            }
        }
    }
    /// Extract the accumulated vector (sorted, zeros removed) and reset.
    pub fn take(&mut self) -> SpVec<F> {
        self.idx.sort_unstable();
        let mut out = Vec::with_capacity(self.idx.len());
        for &i in &self.idx {
            let iu = i as usize;
            self.used[iu] = false;
            let v = std::mem::replace(&mut self.vals[iu], F::zero());
            if !v.is_zero() {
                out.push((i, v));
            }
        }
        self.idx.clear();
        out
    }
}

pub fn spvec_scale<F: Field>(v: &SpVec<F>, c: &F) -> SpVec<F> {
    if c.is_zero() {
        return vec![];
    }
    v.iter().map(|(i, x)| (*i, x.mul(c))).collect()
}

/// Column-major sparse matrix.
#[derive(Clone, Debug)]
pub struct SpMat<F: Field> {
    pub nrows: usize,
    pub ncols: usize,
    pub cols: Vec<SpVec<F>>,
}

impl<F: Field> SpMat<F> {
    pub fn zero(nrows: usize, ncols: usize) -> Self {
        SpMat { nrows, ncols, cols: vec![Vec::new(); ncols] }
    }
    pub fn identity(n: usize) -> Self {
        SpMat { nrows: n, ncols: n, cols: (0..n).map(|i| vec![(i as u32, F::one())]).collect() }
    }
    pub fn scalar(n: usize, c: &F) -> Self {
        if c.is_zero() {
            return Self::zero(n, n);
        }
        SpMat { nrows: n, ncols: n, cols: (0..n).map(|i| vec![(i as u32, c.clone())]).collect() }
    }
    pub fn nnz(&self) -> usize {
        self.cols.iter().map(|c| c.len()).sum()
    }
    pub fn is_zero(&self) -> bool {
        self.cols.iter().all(|c| c.is_empty())
    }
    pub fn mul_vec(&self, v: &SpVec<F>, acc: &mut Acc<F>) -> SpVec<F> {
        debug_assert!(acc.len() >= self.nrows);
        for (j, x) in v {
            acc.axpy(x, &self.cols[*j as usize]);
        }
        acc.take()
    }
    /// self * other
    pub fn matmul(&self, other: &SpMat<F>) -> SpMat<F> {
        assert_eq!(self.ncols, other.nrows, "matmul dimension mismatch");
        let mut acc = Acc::new(self.nrows);
        let cols = other.cols.iter().map(|c| self.mul_vec(c, &mut acc)).collect();
        SpMat { nrows: self.nrows, ncols: other.ncols, cols }
    }
    /// a*self + b*other
    pub fn lincomb(&self, a: &F, other: &SpMat<F>, b: &F) -> SpMat<F> {
        assert_eq!(self.nrows, other.nrows);
        assert_eq!(self.ncols, other.ncols);
        let mut acc = Acc::new(self.nrows);
        let mut cols = Vec::with_capacity(self.ncols);
        for j in 0..self.ncols {
            acc.axpy(a, &self.cols[j]);
            acc.axpy(b, &other.cols[j]);
            cols.push(acc.take());
        }
        SpMat { nrows: self.nrows, ncols: self.ncols, cols }
    }
    pub fn add(&self, other: &SpMat<F>) -> SpMat<F> {
        self.lincomb(&F::one(), other, &F::one())
    }
    pub fn sub(&self, other: &SpMat<F>) -> SpMat<F> {
        self.lincomb(&F::one(), other, &F::one().neg())
    }
    pub fn scale(&self, c: &F) -> SpMat<F> {
        SpMat { nrows: self.nrows, ncols: self.ncols, cols: self.cols.iter().map(|v| spvec_scale(v, c)).collect() }
    }
    pub fn transpose(&self) -> SpMat<F> {
        let mut cols: Vec<SpVec<F>> = vec![Vec::new(); self.nrows];
        for (j, c) in self.cols.iter().enumerate() {
            for (i, x) in c {
                cols[*i as usize].push((j as u32, x.clone()));
            }
        }
        SpMat { nrows: self.ncols, ncols: self.nrows, cols }
    }
    pub fn to_dense(&self) -> Vec<Vec<F>> {
        let mut d = vec![vec![F::zero(); self.ncols]; self.nrows];
        for (j, c) in self.cols.iter().enumerate() {
            for (i, x) in c {
                d[*i as usize][j] = x.clone();
            }
        }
        d
    }
    pub fn eq_mat(&self, other: &SpMat<F>) -> bool {
        self.nrows == other.nrows && self.ncols == other.ncols && self.cols == other.cols
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra::field::F31;
    #[test]
    fn matmul_identity() {
        let mut m = SpMat::<F31>::zero(3, 3);
        m.cols[0] = vec![(1, F31::from_i64(2))];
        m.cols[2] = vec![(0, F31::from_i64(5)), (2, F31::from_i64(-1))];
        let id = SpMat::identity(3);
        assert!(m.matmul(&id).eq_mat(&m));
        assert!(id.matmul(&m).eq_mat(&m));
        let m2 = m.matmul(&m);
        // m^2 e2 = m(5 e0 - e2) = 10 e1 - (5 e0 - e2) = -5e0 + 10 e1 + e2
        assert_eq!(m2.cols[2], vec![(0, F31::from_i64(-5)), (1, F31::from_i64(10)), (2, F31::from_i64(1))]);
        assert!(m.transpose().transpose().eq_mat(&m));
    }
}
