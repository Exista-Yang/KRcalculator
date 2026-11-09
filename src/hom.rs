//! Hom spaces of edge-algebra bimodules.

use crate::algebra::Alg;
use crate::dense;
use crate::field::Field;
use crate::laurent::LPoly;
use crate::sparse::{SpMat, SpVec};
use std::collections::{BTreeMap, HashMap};

/// Kernel of a family of homogeneous operators restricted to each degree block.
/// `ops` are (operator, degree of operator). Returns basis per degree.
pub fn joint_kernel<F: Field>(deg: &[i32], ops: &[(SpMat<F>, i32)]) -> BTreeMap<i32, Vec<SpVec<F>>> {
    let mut block: BTreeMap<i32, Vec<usize>> = BTreeMap::new();
    for (i, &d) in deg.iter().enumerate() {
        block.entry(d).or_default().push(i);
    }
    let mut out = BTreeMap::new();
    for (&d, coords) in &block {
        let nb = coords.len();
        let mut rows: Vec<Vec<F>> = Vec::new();
        for (op, _od) in ops {
            let mut rowmap: HashMap<u32, Vec<F>> = HashMap::new();
            for (k, &c) in coords.iter().enumerate() {
                for (t, x) in &op.cols[c] {
                    rowmap.entry(*t).or_insert_with(|| vec![F::zero(); nb])[k] = x.clone();
                }
            }
            rows.extend(rowmap.into_values());
            if rows.len() > 2 * nb.max(8) {
                dense::rref(&mut rows, nb);
            }
        }
        let ker = dense::kernel(&rows, nb);
        let vecs: Vec<SpVec<F>> = ker
            .into_iter()
            .map(|v| v.into_iter().enumerate().filter(|(_, x)| !x.is_zero()).map(|(k, x)| (coords[k] as u32, x)).collect())
            .collect();
        if !vecs.is_empty() {
            out.insert(d, vecs);
        }
    }
    out
}

/// The operators top_i,j - bot_i,j whose joint kernel is Hom(R, A) (closure).
pub fn closure_ops<F: Field>(a: &Alg<F>) -> Vec<(SpMat<F>, i32)> {
    assert_eq!(a.top_labels, a.bot_labels, "closure needs matching labels");
    let mut ops = Vec::new();
    for i in 0..a.top_labels.len() {
        for j in 0..a.top_labels[i] {
            let d = a.top[i][j].sub(&a.bot[i][j]);
            if !d.is_zero() {
                ops.push((d, 2 * (j as i32 + 1)));
            }
        }
    }
    ops
}

/// Basis (per unshifted degree) of the closure state space Hom(R, A).
pub fn closure_space<F: Field>(a: &Alg<F>) -> BTreeMap<i32, Vec<SpVec<F>>> {
    joint_kernel(&a.deg, &closure_ops(a))
}

/// q-shift turning unshifted algebra degrees into closure q-degrees.
pub fn closure_shift<F: Field>(a: &Alg<F>) -> i32 {
    let n = a.n as i32;
    a.shift - a.bot_labels.iter().map(|&l| l as i32 * (n - l as i32)).sum::<i32>()
}

pub fn closure_graded_dim<F: Field>(a: &Alg<F>) -> LPoly {
    let sh = closure_shift(a);
    let mut p = LPoly::zero();
    for (d, v) in closure_space(a) {
        p.add_term(d + sh, v.len() as i128);
    }
    p
}
