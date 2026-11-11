//! Sparse linear systems over a field: kernel and rank via incremental
//! row reduction.

use crate::field::Field;
use crate::sparse::{Acc, SpVec};
use std::collections::HashMap;

/// Incremental echelon basis of a row space (rows sparse over columns 0..ncols).
pub struct Echelon<F: Field> {
    pub ncols: usize,
    /// pivot column -> row (normalized: leading coefficient 1 at pivot,
    /// pivot is the smallest column index in the row)
    pub rows: HashMap<u32, SpVec<F>>,
    acc: Acc<F>,
}

impl<F: Field> Echelon<F> {
    pub fn new(ncols: usize) -> Self {
        Echelon { ncols, rows: HashMap::new(), acc: Acc::new(ncols.max(1)) }
    }
    pub fn rank(&self) -> usize {
        self.rows.len()
    }
    /// Reduce v against the current basis (eliminating pivot columns).
    pub fn reduce(&mut self, v: &SpVec<F>) -> SpVec<F> {
        if v.is_empty() {
            return vec![];
        }
        // dense-ish reduction using a heap of columns would be ideal; we do
        // repeated sweeps in increasing column order.
        let mut cur = v.clone();
        loop {
            let mut changed = false;
            let mut k = 0;
            while k < cur.len() {
                let (c, x) = (cur[k].0, cur[k].1.clone());
                if let Some(r) = self.rows.get(&c) {
                    // cur -= x * r
                    let f = x.neg();
                    self.acc.axpy(&F::one(), &cur);
                    self.acc.axpy(&f, r);
                    cur = self.acc.take();
                    changed = true;
                    // restart from the position of c (entries before c unchanged)
                    k = cur.partition_point(|e| e.0 < c);
                    continue;
                }
                k += 1;
            }
            if !changed {
                break;
            }
        }
        cur
    }
    /// Insert a row; returns true if it increased the rank.
    pub fn insert(&mut self, v: &SpVec<F>) -> bool {
        let r = self.reduce(v);
        if r.is_empty() {
            return false;
        }
        let inv = r[0].1.inv();
        let r: SpVec<F> = r.into_iter().map(|(i, x)| (i, x.mul(&inv))).collect();
        self.rows.insert(r[0].0, r);
        true
    }
    /// Fully reduce all rows (RREF).
    pub fn rref(&mut self) {
        let mut piv: Vec<u32> = self.rows.keys().copied().collect();
        piv.sort_unstable();
        // back substitution from the largest pivot downwards
        for &p in piv.iter().rev() {
            let row = self.rows.remove(&p).unwrap();
            // reduce entries other than the pivot
            let head = vec![row[0].clone()];
            let tail: SpVec<F> = row[1..].to_vec();
            let red = self.reduce(&tail);
            let mut full = head;
            full.extend(red);
            self.rows.insert(p, full);
        }
    }
    /// Kernel basis of the linear map whose rows are the inserted equations.
    pub fn kernel(&mut self) -> Vec<SpVec<F>> {
        self.rref();
        let mut out = vec![];
        for free in 0..self.ncols as u32 {
            if self.rows.contains_key(&free) {
                continue;
            }
            let mut v: SpVec<F> = vec![(free, F::one())];
            for (&p, row) in &self.rows {
                if let Ok(i) = row.binary_search_by_key(&free, |e| e.0) {
                    v.push((p, row[i].1.neg()));
                }
            }
            v.sort_unstable_by_key(|e| e.0);
            out.push(v);
        }
        out
    }
}

/// Rank of a set of sparse vectors.
pub fn rank_of<F: Field>(vs: &[SpVec<F>], ncols: usize) -> usize {
    let mut e = Echelon::new(ncols);
    for v in vs {
        e.insert(v);
    }
    e.rank()
}

/// Fast incremental rank computation: pivot rows stored by leading column, rows
/// reduced in one left-to-right pass with a dense accumulator.
pub struct RankEchelon<F: Field> {
    ncols: usize,
    /// pivot row for each column (leading coefficient normalized to 1), entries strictly right of the pivot
    piv: Vec<Option<Vec<(u32, F)>>>,
    acc: Vec<F>,
    rank: usize,
}

impl<F: Field> RankEchelon<F> {
    pub fn new(ncols: usize) -> Self {
        RankEchelon { ncols, piv: vec![None; ncols], acc: vec![F::zero(); ncols], rank: 0 }
    }
    pub fn rank(&self) -> usize {
        self.rank
    }
    /// Insert a sparse row; returns true if the rank increased.
    pub fn insert(&mut self, v: &SpVec<F>) -> bool {
        if v.is_empty() {
            return false;
        }
        let mut lo = usize::MAX;
        let mut hi = 0usize;
        for (i, x) in v {
            let i = *i as usize;
            self.acc[i] = x.clone();
            lo = lo.min(i);
            hi = hi.max(i);
        }
        let mut col = lo;
        let mut lead: Option<usize> = None;
        while col <= hi && col < self.ncols {
            if !self.acc[col].is_zero() {
                if let Some(p) = &self.piv[col] {
                    let c = self.acc[col].neg();
                    self.acc[col] = F::zero();
                    for (j, y) in p {
                        let j = *j as usize;
                        self.acc[j].fma(&c, y);
                        if j > hi {
                            hi = j;
                        }
                    }
                } else if lead.is_none() {
                    lead = Some(col);
                }
            }
            col += 1;
        }
        match lead {
            None => {
                for c in lo..=hi.min(self.ncols - 1) {
                    self.acc[c] = F::zero();
                }
                false
            }
            Some(l) => {
                let inv = self.acc[l].inv();
                let mut row = vec![];
                for c in (l + 1)..=hi.min(self.ncols - 1) {
                    if !self.acc[c].is_zero() {
                        row.push((c as u32, self.acc[c].mul(&inv)));
                    }
                }
                for c in lo..=hi.min(self.ncols - 1) {
                    self.acc[c] = F::zero();
                }
                self.piv[l] = Some(row);
                self.rank += 1;
                true
            }
        }
    }
}

/// Kernel of the linear map whose rows are `eqs` (columns 0..ncols): the reduced basis, one
/// vector per free column (1 there, 0 at the other free columns, sorted), i.e. exactly
/// `Echelon::kernel` after inserting every row (the reduced kernel basis is determined by the
/// row space and the column order).
///
/// Rows are reduced with a dense accumulator. After `run` consecutive redundant rows, when the
/// current kernel has dimension <= `KER_CHECK_MAX`, the kernel of the rows inserted so far is
/// evaluated on every remaining row; violated rows are inserted and the check repeats; when no
/// remaining row is violated the scan stops. This is exact: kernel(all) ⊆ kernel(inserted)
/// always, and a basis of kernel(inserted) annihilating every row gives the reverse inclusion.
pub fn kernel_lazy<F: Field>(ncols: usize, eqs: &[SpVec<F>]) -> Vec<SpVec<F>> {
    const KER_CHECK_MAX: usize = 128;
    let mut piv: Vec<Option<Vec<(u32, F)>>> = vec![None; ncols];
    let mut acc: Vec<F> = vec![F::zero(); ncols.max(1)];
    let mut rank = 0usize;
    let mut insert = |v: &SpVec<F>, piv: &mut Vec<Option<Vec<(u32, F)>>>, acc: &mut Vec<F>| -> bool {
        if v.is_empty() {
            return false;
        }
        let mut lo = usize::MAX;
        let mut hi = 0usize;
        for (i, x) in v {
            let i = *i as usize;
            acc[i].add_assign(x);
            lo = lo.min(i);
            hi = hi.max(i);
        }
        let mut lead: Option<usize> = None;
        let mut col = lo;
        while col <= hi {
            if !acc[col].is_zero() {
                if let Some(p) = &piv[col] {
                    let c = acc[col].neg();
                    acc[col] = F::zero();
                    for (j, y) in p {
                        let j = *j as usize;
                        acc[j].fma(&c, y);
                        if j > hi {
                            hi = j;
                        }
                    }
                } else if lead.is_none() {
                    lead = Some(col);
                }
            }
            col += 1;
        }
        let r = match lead {
            None => false,
            Some(l) => {
                let inv = acc[l].inv();
                let mut row = vec![];
                for c in (l + 1)..=hi {
                    if !acc[c].is_zero() {
                        row.push((c as u32, acc[c].mul(&inv)));
                    }
                }
                piv[l] = Some(row);
                true
            }
        };
        for c in lo..=hi {
            acc[c] = F::zero();
        }
        r
    };
    // dense kernel basis of the current pivot rows (back substitution)
    let kernel_dense = |piv: &Vec<Option<Vec<(u32, F)>>>| -> Vec<(usize, Vec<F>)> {
        let mut out = vec![];
        for f in 0..ncols {
            if piv[f].is_some() {
                continue;
            }
            let mut x = vec![F::zero(); ncols];
            x[f] = F::one();
            for p in (0..f).rev() {
                if let Some(row) = &piv[p] {
                    let mut s = F::zero();
                    for (j, y) in row {
                        let xj = &x[*j as usize];
                        if !xj.is_zero() {
                            s.fma(y, xj);
                        }
                    }
                    x[p] = s.neg();
                }
            }
            out.push((f, x));
        }
        out
    };
    let violated = |v: &SpVec<F>, ker: &[(usize, Vec<F>)]| -> bool {
        ker.iter().any(|(_, x)| {
            let mut s = F::zero();
            for (j, y) in v {
                s.fma(y, &x[*j as usize]);
            }
            !s.is_zero()
        })
    };
    let mut done = vec![false; eqs.len()];
    let mut run_max = 64usize;
    let mut run = 0usize;
    let mut i = 0usize;
    let mut finished_early: Option<Vec<(usize, Vec<F>)>> = None;
    while i < eqs.len() {
        if !done[i] {
            if insert(&eqs[i], &mut piv, &mut acc) {
                rank += 1;
                run = 0;
            } else {
                run += 1;
            }
            done[i] = true;
        }
        i += 1;
        if run >= run_max && ncols - rank <= KER_CHECK_MAX && i < eqs.len() {
            run = 0;
            loop {
                let ker = kernel_dense(&piv);
                let bad: Vec<usize> = (i..eqs.len()).filter(|&j| !done[j] && violated(&eqs[j], &ker)).collect();
                if bad.is_empty() {
                    finished_early = Some(ker);
                    break;
                }
                for j in bad {
                    if insert(&eqs[j], &mut piv, &mut acc) {
                        rank += 1;
                    }
                    done[j] = true;
                }
            }
            if finished_early.is_some() {
                break;
            }
            run_max *= 2;
        }
    }
    let _ = rank;
    let ker = finished_early.unwrap_or_else(|| kernel_dense(&piv));
    ker.into_iter()
        .map(|(_, x)| x.into_iter().enumerate().filter(|(_, v)| !v.is_zero()).map(|(j, v)| (j as u32, v)).collect())
        .collect()
}
