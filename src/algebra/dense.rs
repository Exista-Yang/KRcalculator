//! Dense linear algebra over a field: row reduction, rank, kernel.

use crate::algebra::field::Field;

/// In-place reduced row echelon form. Columns are scanned in the order given by
/// `col_order` (a permutation of 0..ncols), so pivots are preferably chosen
/// early in that order. Returns the pivot column of each nonzero row; rows are
/// reordered so that the first `rank` rows are the pivot rows, the rest removed.
pub fn rref_with_order<F: Field>(rows: &mut Vec<Vec<F>>, col_order: &[usize]) -> Vec<usize> {
    let nrows = rows.len();
    let mut pivots = Vec::new();
    let mut r = 0;
    for &c in col_order {
        if r >= nrows {
            break;
        }
        // find pivot
        let mut p = None;
        for i in r..nrows {
            if !rows[i][c].is_zero() {
                p = Some(i);
                break;
            }
        }
        let Some(p) = p else { continue };
        rows.swap(r, p);
        let inv = rows[r][c].inv();
        if !inv.is_one() {
            for x in rows[r].iter_mut() {
                if !x.is_zero() {
                    *x = x.mul(&inv);
                }
            }
        }
        let prow = rows[r].clone();
        let nz: Vec<usize> = (0..prow.len()).filter(|&k| !prow[k].is_zero()).collect();
        for i in 0..nrows {
            if i == r || rows[i][c].is_zero() {
                continue;
            }
            let f = rows[i][c].neg();
            let row = &mut rows[i];
            for &k in &nz {
                row[k].fma(&f, &prow[k]);
            }
        }
        pivots.push(c);
        r += 1;
    }
    rows.truncate(r);
    pivots
}

pub fn rref<F: Field>(rows: &mut Vec<Vec<F>>, ncols: usize) -> Vec<usize> {
    let order: Vec<usize> = (0..ncols).collect();
    rref_with_order(rows, &order)
}

pub fn rank<F: Field>(rows: &[Vec<F>], ncols: usize) -> usize {
    let mut r = rows.to_vec();
    rref(&mut r, ncols).len()
}

/// Basis of the kernel {x : A x = 0} of the matrix given by `rows` (each of length ncols).
pub fn kernel<F: Field>(rows: &[Vec<F>], ncols: usize) -> Vec<Vec<F>> {
    let mut r = rows.to_vec();
    let piv = rref(&mut r, ncols);
    let mut is_piv = vec![false; ncols];
    for &p in &piv {
        is_piv[p] = true;
    }
    let mut basis = Vec::new();
    for free in 0..ncols {
        if is_piv[free] {
            continue;
        }
        let mut v = vec![F::zero(); ncols];
        v[free] = F::one();
        for (i, &p) in piv.iter().enumerate() {
            v[p] = r[i][free].neg();
        }
        basis.push(v);
    }
    basis
}

/// Solve A x = b for x (one solution) if consistent. `rows` is A (m x n).
pub fn solve<F: Field>(rows: &[Vec<F>], ncols: usize, b: &[F]) -> Option<Vec<F>> {
    let mut aug: Vec<Vec<F>> = rows
        .iter()
        .zip(b.iter())
        .map(|(r, x)| {
            let mut v = r.clone();
            v.push(x.clone());
            v
        })
        .collect();
    let piv = rref(&mut aug, ncols + 1);
    if piv.last() == Some(&ncols) {
        return None;
    }
    let mut x = vec![F::zero(); ncols];
    for (i, &p) in piv.iter().enumerate() {
        x[p] = aug[i][ncols].clone();
    }
    Some(x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra::field::{F31, Q};
    #[test]
    fn kernel_small() {
        let f = |x: i64| Q::from_i64(x);
        let a = vec![vec![f(1), f(2), f(3)], vec![f(2), f(4), f(6)]];
        let k = kernel(&a, 3);
        assert_eq!(k.len(), 2);
        for v in &k {
            for row in &a {
                let mut s = Q::zero();
                for j in 0..3 {
                    s.fma(&row[j], &v[j]);
                }
                assert!(s.is_zero());
            }
        }
        let g = |x: i64| F31::from_i64(x);
        assert_eq!(rank(&[vec![g(1), g(1)], vec![g(1), g(-1)]], 2), 2);
    }
}
