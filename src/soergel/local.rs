//! Equivariant localization of ladder (singular Bott–Samelson) bimodules at
//! torus fixed points, evaluated numerically at a point of T = k[x_1..x_K]
//! (K = total label). A fixed point records, for each rung, which subset of the
//! giver's roots moves; loc_p(r ⊗ β) = r(roots) · Π_rungs e^α(moved roots).
//! Bimodule maps become block-diagonal: blocks indexed by the top configuration.

use crate::algebra::field::Field;
use crate::algebra::mpoly::mono_exp;
use crate::algebra::pmat::PMat;
use crate::soergel::sbim::SBim;
use crate::algebra::symfun::GrTable;

/// Roots are numbered 0..K; bottom upright i owns roots off[i]..off[i]+λ_i.
#[derive(Clone, Debug)]
pub struct FixedPoint {
    /// root bitmask on each top upright
    pub top: Vec<u64>,
    /// chosen (moved) root bitmask at each rung
    pub moved: Vec<u64>,
}

pub fn root_offsets(labels: &[usize]) -> Vec<usize> {
    let mut off = vec![];
    let mut k = 0;
    for &l in labels {
        off.push(k);
        k += l;
    }
    off
}

fn subsets_of(mask: u64, r: usize) -> Vec<u64> {
    let bits: Vec<u32> = (0..64).filter(|&i| mask >> i & 1 == 1).collect();
    let mut out = vec![];
    let nb = bits.len();
    if r > nb {
        return out;
    }
    let mut idx: Vec<usize> = (0..r).collect();
    loop {
        out.push(idx.iter().fold(0u64, |m, &i| m | 1u64 << bits[i]));
        // next combination
        let mut i = r;
        loop {
            if i == 0 {
                return out;
            }
            i -= 1;
            if idx[i] != i + nb - r {
                break;
            }
            if i == 0 && idx[0] == nb - r {
                return out;
            }
        }
        idx[i] += 1;
        for j in i + 1..r {
            idx[j] = idx[j - 1] + 1;
        }
    }
}

/// Fixed points in the same order as the fibre basis enumeration is NOT implied;
/// fixed points are enumerated by successive subset choices.
pub fn fixed_points<F: Field>(b: &SBim<F>) -> Vec<FixedPoint> {
    let off = root_offsets(&b.bot_labels);
    let start: Vec<u64> = b.bot_labels.iter().enumerate().map(|(i, &l)| ((1u64 << l) - 1) << off[i]).collect();
    let mut pts = vec![FixedPoint { top: start, moved: vec![] }];
    for rg in &b.rungs {
        let (rec, giv) = rg.ends();
        let mut next = vec![];
        for p in &pts {
            for c in subsets_of(p.top[giv], rg.r) {
                let mut t = p.top.clone();
                t[giv] &= !c;
                t[rec] |= c;
                let mut mv = p.moved.clone();
                mv.push(c);
                next.push(FixedPoint { top: t, moved: mv });
            }
        }
        pts = next;
    }
    pts
}

/// Elementary symmetric polynomials e_0..e_k of the roots in `mask`.
pub fn esym<F: Field>(mask: u64, x: &[F]) -> Vec<F> {
    let mut e = vec![F::one()];
    for i in 0..64 {
        if mask >> i & 1 == 1 {
            let mut ne = vec![F::zero(); e.len() + 1];
            for (j, c) in e.iter().enumerate() {
                ne[j].add_assign(c);
                ne[j + 1].fma(c, &x[i]);
            }
            e = ne;
        }
    }
    e
}

/// Value of a bottom-ring polynomial at the point (e_{i,j} = e_j(roots of upright i)).
pub fn eval_bottom<F: Field>(b: &SBim<F>, p: &crate::algebra::mpoly::Poly<F>, x: &[F]) -> F {
    let off = root_offsets(&b.bot_labels);
    let es: Vec<Vec<F>> = b.bot_labels.iter().enumerate().map(|(i, &l)| esym(((1u64 << l) - 1) << off[i], x)).collect();
    let mut s = F::zero();
    for (m, c) in &p.t {
        let mut v = c.clone();
        for (i, &l) in b.bot_labels.iter().enumerate() {
            for j in 0..l {
                for _ in 0..mono_exp(*m, b.voff[i] + j) {
                    v = v.mul(&es[i][j + 1]);
                }
            }
        }
        s.add_assign(&v);
    }
    s
}

/// Localization matrix L[p][β] at the numeric point x (|P| × |F|).
pub fn loc_matrix<F: Field>(b: &SBim<F>, pts: &[FixedPoint], x: &[F]) -> Vec<Vec<F>> {
    // fibre index is mixed radix over rungs: β = ((β_1)*f_2 + β_2)*f_3 + ...
    let mut tabs = vec![];
    let mut labels = b.bot_labels.clone();
    for rg in &b.rungs {
        let (rec, giv) = rg.ends();
        tabs.push(GrTable::get(rg.r, labels[giv]));
        labels[giv] -= rg.r;
        labels[rec] += rg.r;
    }
    let fs: Vec<usize> = tabs.iter().map(|t| t.fsize()).collect();
    let total: usize = fs.iter().product();
    assert_eq!(total, b.fsize());
    let mut out = Vec::with_capacity(pts.len());
    for p in pts {
        // per rung: values of basis monomials e^α(moved)
        let vals: Vec<Vec<F>> = tabs
            .iter()
            .zip(p.moved.iter())
            .map(|(t, &c)| {
                let e = esym(c, x);
                t.basis
                    .iter()
                    .map(|al| {
                        let mut v = F::one();
                        for (j, &k) in al.iter().enumerate() {
                            for _ in 0..k {
                                v = v.mul(&e[j + 1]);
                            }
                        }
                        v
                    })
                    .collect()
            })
            .collect();
        let mut row = vec![F::zero(); total];
        for (beta, slot) in row.iter_mut().enumerate() {
            let mut rem = beta;
            let mut v = F::one();
            for k in (0..fs.len()).rev() {
                let a = rem % fs[k];
                rem /= fs[k];
                v = v.mul(&vals[k][a]);
            }
            *slot = v;
        }
        out.push(row);
    }
    out
}

/// Numeric evaluation of a polynomial matrix (bottom ring entries) at x.
pub fn eval_pmat<F: Field>(b: &SBim<F>, m: &PMat<F>, x: &[F]) -> Vec<Vec<F>> {
    let mut d = vec![vec![F::zero(); m.ncols]; m.nrows];
    for (j, col) in m.cols.iter().enumerate() {
        for (i, p) in col {
            d[*i as usize][j] = eval_bottom(b, p, x);
        }
    }
    d
}

pub fn matmul<F: Field>(a: &[Vec<F>], b: &[Vec<F>]) -> Vec<Vec<F>> {
    let n = a.len();
    let m = if b.is_empty() { 0 } else { b[0].len() };
    let k = b.len();
    let mut c = vec![vec![F::zero(); m]; n];
    for i in 0..n {
        for t in 0..k {
            if a[i][t].is_zero() {
                continue;
            }
            for j in 0..m {
                if !b[t][j].is_zero() {
                    c[i][j].fma(&a[i][t], &b[t][j]);
                }
            }
        }
    }
    c
}

pub fn inverse<F: Field>(a: &[Vec<F>]) -> Option<Vec<Vec<F>>> {
    let n = a.len();
    let mut aug: Vec<Vec<F>> = a
        .iter()
        .enumerate()
        .map(|(i, r)| {
            let mut v = r.clone();
            v.extend((0..n).map(|j| if i == j { F::one() } else { F::zero() }));
            v
        })
        .collect();
    let piv = crate::algebra::dense::rref(&mut aug, 2 * n);
    if piv.len() < n || piv.iter().take(n).enumerate().any(|(i, &p)| p != i) {
        return None;
    }
    Some(aug.into_iter().map(|r| r[n..].to_vec()).collect())
}

/// C2-7: the columns `cols` of a⁻¹ (n × k, column order of `cols`) by Gaussian elimination of
/// [a | E_cols] (forward elimination with row pivoting, then back substitution): n³/3 + n²k
/// instead of the n³-type rref of [a | I] in `inverse`. None if a is singular.
pub fn inverse_cols<F: Field>(a: &[Vec<F>], cols: &[usize]) -> Option<Vec<Vec<F>>> {
    let n = a.len();
    let k = cols.len();
    let mut m: Vec<Vec<F>> = a.to_vec();
    let mut r: Vec<Vec<F>> = (0..n).map(|i| cols.iter().map(|&c| if c == i { F::one() } else { F::zero() }).collect()).collect();
    for c in 0..n {
        let piv = (c..n).find(|&i| !m[i][c].is_zero())?;
        m.swap(c, piv);
        r.swap(c, piv);
        let inv = m[c][c].inv();
        let (top, bot) = m.split_at_mut(c + 1);
        let (rtop, rbot) = r.split_at_mut(c + 1);
        let prow = &top[c];
        let prr = &rtop[c];
        for (row, rr) in bot.iter_mut().zip(rbot.iter_mut()) {
            if row[c].is_zero() {
                continue;
            }
            let f = row[c].mul(&inv).neg();
            for j in c..n {
                if !prow[j].is_zero() {
                    row[j].fma(&f, &prow[j]);
                }
            }
            for j in 0..k {
                if !prr[j].is_zero() {
                    rr[j].fma(&f, &prr[j]);
                }
            }
        }
    }
    // back substitution: x_i = (r_i - Σ_{j>i} m_ij x_j) / m_ii
    let mut x: Vec<Vec<F>> = vec![vec![F::zero(); k]; n];
    for i in (0..n).rev() {
        let mut v = r[i].clone();
        for j in i + 1..n {
            if m[i][j].is_zero() {
                continue;
            }
            let f = m[i][j].neg();
            for t in 0..k {
                if !x[j][t].is_zero() {
                    v[t].fma(&f, &x[j][t]);
                }
            }
        }
        let inv = m[i][i].inv();
        for t in 0..k {
            v[t] = v[t].mul(&inv);
        }
        x[i] = v;
    }
    Some(x)
}

/// Localize a bimodule map Φ: B -> B' at x: f = L' Φ L^{-1} (|P'| × |P|).
pub fn localize_map<F: Field>(b: &SBim<F>, bp: &SBim<F>, phi: &PMat<F>, x: &[F]) -> (Vec<FixedPoint>, Vec<FixedPoint>, Vec<Vec<F>>) {
    let p = fixed_points(b);
    let pp = fixed_points(bp);
    let l = loc_matrix(b, &p, x);
    let lp = loc_matrix(bp, &pp, x);
    let linv = inverse(&l).expect("localization matrix singular at this point");
    let ph = eval_pmat(b, phi, x);
    let f = matmul(&matmul(&lp, &ph), &linv);
    (p, pp, f)
}
