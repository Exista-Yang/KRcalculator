//! Sparse linear systems over a field: kernel and rank via incremental
//! row reduction.

use crate::algebra::field::Field;
use crate::algebra::sparse::{Acc, SpVec};
use std::collections::HashMap;
use std::sync::atomic::{AtomicUsize, Ordering};

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
    let (ker, tr) = kernel_lazy_traced(ncols, eqs);
    LAZY_CALLS.fetch_add(1, Ordering::Relaxed);
    if tr.early_stop {
        LAZY_EARLY_STOPS.fetch_add(1, Ordering::Relaxed);
    }
    if tr.reinserted > 0 {
        LAZY_REINSERTS.fetch_add(1, Ordering::Relaxed);
    }
    ker
}

/// What one `kernel_lazy` call did (branch coverage of the early stop; tests and diagnostics)
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct KernelTrace {
    /// number of kernel checks started (after a run of redundant rows)
    pub checks: usize,
    /// the scan stopped before the last row (every remaining row annihilated the kernel)
    pub early_stop: bool,
    /// rows left unscanned by the early stop
    pub skipped: usize,
    /// check rounds that found violated rows, and the rows inserted by them
    pub reinsert_rounds: usize,
    pub reinserted: usize,
}

static LAZY_CALLS: AtomicUsize = AtomicUsize::new(0);
static LAZY_EARLY_STOPS: AtomicUsize = AtomicUsize::new(0);
static LAZY_REINSERTS: AtomicUsize = AtomicUsize::new(0);

/// process-wide counters of `kernel_lazy`: (calls, calls that stopped early, calls in which a
/// check inserted violated rows)
pub fn kernel_lazy_counters() -> (usize, usize, usize) {
    (LAZY_CALLS.load(Ordering::Relaxed), LAZY_EARLY_STOPS.load(Ordering::Relaxed), LAZY_REINSERTS.load(Ordering::Relaxed))
}

/// `kernel_lazy` with its branch trace
pub fn kernel_lazy_traced<F: Field>(ncols: usize, eqs: &[SpVec<F>]) -> (Vec<SpVec<F>>, KernelTrace) {
    const KER_CHECK_MAX: usize = 128;
    let mut piv: Vec<Option<Vec<(u32, F)>>> = vec![None; ncols];
    let mut acc: Vec<F> = vec![F::zero(); ncols.max(1)];
    let mut rank = 0usize;
    // M5c: sparse sweep. The accumulator's nonzero columns are visited in increasing order through
    // a min-heap (a column enters it once, when first touched) instead of a dense scan of every
    // column between the row's first and last entry: the same pivots are applied in the same
    // order with the same arithmetic, so the stored pivot rows (and the kernel) are identical;
    // the cost is O(touched · log) instead of O(span).
    let mut mark: Vec<bool> = vec![false; ncols.max(1)];
    let mut heap: std::collections::BinaryHeap<std::cmp::Reverse<u32>> = std::collections::BinaryHeap::new();
    let mut touched: Vec<u32> = vec![];
    let mut insert = |v: &SpVec<F>, piv: &mut Vec<Option<Vec<(u32, F)>>>, acc: &mut Vec<F>| -> bool {
        if v.is_empty() {
            return false;
        }
        for (i, x) in v {
            let iu = *i as usize;
            acc[iu].add_assign(x);
            if !mark[iu] {
                mark[iu] = true;
                heap.push(std::cmp::Reverse(*i));
                touched.push(*i);
            }
        }
        let mut lead: Option<usize> = None;
        let mut rest: Vec<u32> = vec![];
        while let Some(std::cmp::Reverse(c32)) = heap.pop() {
            let col = c32 as usize;
            if acc[col].is_zero() {
                continue;
            }
            if let Some(p) = &piv[col] {
                let c = acc[col].neg();
                acc[col] = F::zero();
                for (j, y) in p {
                    let ju = *j as usize;
                    acc[ju].fma(&c, y);
                    if !mark[ju] {
                        mark[ju] = true;
                        heap.push(std::cmp::Reverse(*j));
                        touched.push(*j);
                    }
                }
            } else if lead.is_none() {
                lead = Some(col);
            } else {
                // final value: every pivot left of col has been applied (columns pop in order)
                rest.push(c32);
            }
        }
        let r = match lead {
            None => false,
            Some(l) => {
                let inv = acc[l].inv();
                let mut row = Vec::with_capacity(rest.len());
                for c in rest {
                    let x = &acc[c as usize];
                    if !x.is_zero() {
                        row.push((c, x.mul(&inv)));
                    }
                }
                piv[l] = Some(row);
                true
            }
        };
        for c in touched.drain(..) {
            acc[c as usize] = F::zero();
            mark[c as usize] = false;
        }
        r
    };
    // back substitution for the free column f into the zeroed buffer x (pivots p > f stay 0)
    let back = |piv: &Vec<Option<Vec<(u32, F)>>>, f: usize, x: &mut Vec<F>| {
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
    };
    // dense kernel basis of the current pivot rows, used only for the check (dimension <= KER_CHECK_MAX)
    let kernel_dense = |piv: &Vec<Option<Vec<(u32, F)>>>| -> Vec<(usize, Vec<F>)> {
        let mut out = vec![];
        for f in 0..ncols {
            if piv[f].is_some() {
                continue;
            }
            let mut x = vec![F::zero(); ncols];
            back(piv, f, &mut x);
            out.push((f, x));
        }
        out
    };
    // the final kernel, built sparsely: one reused dense buffer, each vector pushed sparse at once
    // (peak memory O(ncols + nnz(kernel)) instead of dim · ncols)
    let kernel_sparse = |piv: &Vec<Option<Vec<(u32, F)>>>| -> Vec<SpVec<F>> {
        let mut out = vec![];
        let mut x = vec![F::zero(); ncols.max(1)];
        for f in 0..ncols {
            if piv[f].is_some() {
                continue;
            }
            back(piv, f, &mut x);
            let mut v: SpVec<F> = vec![];
            for (j, xj) in x.iter_mut().enumerate().take(f + 1) {
                if !xj.is_zero() {
                    v.push((j as u32, std::mem::replace(xj, F::zero())));
                }
            }
            out.push(v);
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
    let mut tr = KernelTrace::default();
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
            tr.checks += 1;
            loop {
                let ker = kernel_dense(&piv);
                let bad: Vec<usize> = (i..eqs.len()).filter(|&j| !done[j] && violated(&eqs[j], &ker)).collect();
                if bad.is_empty() {
                    tr.early_stop = true;
                    tr.skipped = (i..eqs.len()).filter(|&j| !done[j]).count();
                    break;
                }
                tr.reinsert_rounds += 1;
                for j in bad {
                    if insert(&eqs[j], &mut piv, &mut acc) {
                        rank += 1;
                    }
                    tr.reinserted += 1;
                    done[j] = true;
                }
            }
            if tr.early_stop {
                break;
            }
            run_max *= 2;
        }
    }
    let _ = rank;
    (kernel_sparse(&piv), tr)
}
