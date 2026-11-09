//! Universal structure constants for relative Grassmannians.
//!
//! Gr(r,b) over the base k[c_1..c_b] is the ring k[e_1..e_r, f_1..f_{b-r}] in which
//! c(t) = e(t) f(t), i.e. c_j = sum_i e_i f_{j-i}. It is free over k[c] with basis
//! { e^alpha : alpha in N^r, |alpha| <= b-r } (e-monomials, unitriangular to Schur
//! polynomials in a box). We tabulate multiplication by e_j and by f_j in this basis,
//! with coefficients polynomials in c (integers).

use crate::dense::rref;
use crate::field::{Field, Fp};
use crate::poly::{monomials_of_weight, MPoly, Mono};
use std::collections::HashMap;
use std::sync::{Mutex, OnceLock};

/// Polynomial in c_1..c_b with integer coefficients: list of (exponent vector, coeff).
pub type CPoly = Vec<(Mono, i64)>;

#[derive(Clone, Debug)]
pub struct GrTable {
    pub r: usize,
    pub b: usize,
    /// basis exponent vectors alpha (length r)
    pub basis: Vec<Mono>,
    /// weight (= half the q-degree) of each basis element
    pub wt: Vec<u32>,
    /// e_mul[j-1][a] = list of (beta index, coefficient polynomial in c)
    pub e_mul: Vec<Vec<Vec<(usize, CPoly)>>>,
    /// f_mul[j-1][a] for j = 1..b-r
    pub f_mul: Vec<Vec<Vec<(usize, CPoly)>>>,
}

fn weight(m: &[u16]) -> u32 {
    m.iter().enumerate().map(|(i, &k)| (i as u32 + 1) * k as u32).sum()
}

impl GrTable {
    pub fn fsize(&self) -> usize {
        self.basis.len()
    }

    pub fn compute(r: usize, b: usize) -> GrTable {
        assert!(r <= b);
        let s = b - r;
        // fibre basis
        let mut basis: Vec<Mono> = Vec::new();
        for tot in 0..=s {
            // all alpha with |alpha| = tot
            let mut cur = vec![0u16; r];
            fn rec(i: usize, left: usize, cur: &mut Vec<u16>, out: &mut Vec<Mono>) {
                if i == cur.len() {
                    if left == 0 {
                        out.push(cur.clone());
                    }
                    return;
                }
                for k in 0..=left {
                    cur[i] = k as u16;
                    rec(i + 1, left - k, cur, out);
                }
                cur[i] = 0;
            }
            if r == 0 {
                if tot == 0 {
                    basis.push(vec![]);
                }
            } else {
                rec(0, tot, &mut cur, &mut basis);
            }
        }
        basis.sort_by_key(|m| (weight(m), m.clone()));
        let wt: Vec<u32> = basis.iter().map(|m| weight(m)).collect();
        let bindex: HashMap<Mono, usize> = basis.iter().cloned().enumerate().map(|(i, m)| (m, i)).collect();
        let _ = bindex;

        // polynomial ring in e_1..e_r, f_1..f_s
        let nv = r + s;
        let vw: Vec<u32> = (1..=r as u32).chain(1..=s as u32).collect();
        let evar = |i: usize| -> MPoly {
            // e_i, i in 0..=r
            if i == 0 { MPoly::one(nv) } else if i <= r { MPoly::var(nv, i - 1) } else { MPoly::zero(nv) }
        };
        let fvar = |i: usize| -> MPoly {
            if i == 0 { MPoly::one(nv) } else if i <= s { MPoly::var(nv, r + i - 1) } else { MPoly::zero(nv) }
        };
        let cpoly: Vec<MPoly> = (0..=b)
            .map(|j| {
                let mut p = MPoly::zero(nv);
                for i in 0..=j {
                    p = p.add(&evar(i).mul(&fvar(j - i)));
                }
                p
            })
            .collect();
        let cw: Vec<u32> = (1..=b as u32).collect();

        // embed alpha (e-exponents) into a monomial of the big ring
        let emono = |a: &Mono| -> Mono {
            let mut m = vec![0u16; nv];
            m[..r].copy_from_slice(a);
            m
        };

        // targets: (kind, j, a) -> polynomial
        let mut targets: Vec<((u8, usize, usize), MPoly)> = Vec::new();
        for j in 1..=r {
            for (ai, a) in basis.iter().enumerate() {
                let mut m = emono(a);
                m[j - 1] += 1;
                targets.push(((0, j, ai), MPoly::monomial(m)));
            }
        }
        for j in 1..=s {
            for (ai, a) in basis.iter().enumerate() {
                let mut m = emono(a);
                m[r + j - 1] += 1;
                targets.push(((1, j, ai), MPoly::monomial(m)));
            }
        }
        // group by weight
        let mut by_w: HashMap<u32, Vec<usize>> = HashMap::new();
        for (ti, (_, p)) in targets.iter().enumerate() {
            let w = weight_big(p.t.keys().next().unwrap(), &vw);
            by_w.entry(w).or_default().push(ti);
        }
        let mut cmemo: HashMap<Mono, MPoly> = HashMap::new();
        let mut results: Vec<Vec<(usize, CPoly)>> = vec![Vec::new(); targets.len()];
        for (&d, tis) in &by_w {
            // columns: pairs (gamma, beta)
            let mut pairs: Vec<(Mono, usize)> = Vec::new();
            for (bi, &bw) in wt.iter().enumerate() {
                if bw > d {
                    continue;
                }
                for g in monomials_of_weight(&cw, d - bw) {
                    pairs.push((g, bi));
                }
            }
            let rows_m = monomials_of_weight(&vw, d);
            assert_eq!(rows_m.len(), pairs.len(), "Gr({},{}) weight {}: not square", r, b, d);
            let ridx: HashMap<Mono, usize> = rows_m.iter().cloned().enumerate().map(|(i, m)| (m, i)).collect();
            let np = pairs.len();
            // Expand the columns once (integer coefficients).
            let mut colpolys: Vec<Vec<(usize, i64)>> = Vec::with_capacity(np);
            for (g, bi) in pairs.iter() {
                let cg = cmemo
                    .entry(g.clone())
                    .or_insert_with(|| {
                        let mut p = MPoly::one(nv);
                        for (j, &k) in g.iter().enumerate() {
                            p = p.mul(&cpoly[j + 1].pow(k as u32));
                        }
                        p
                    })
                    .clone();
                let full = cg.mul(&MPoly::monomial(emono(&basis[*bi])));
                colpolys.push(full.t.iter().map(|(m, &c)| (ridx[m], c)).collect());
            }
            let tpolys: Vec<Vec<(usize, i64)>> =
                tis.iter().map(|&ti| targets[ti].1.t.iter().map(|(m, &c)| (ridx[m], c)).collect()).collect();
            let sol1 = solve_mod::<2147483647>(&colpolys, &tpolys, rows_m.len(), r, b);
            let sol2 = solve_mod::<1000000007>(&colpolys, &tpolys, rows_m.len(), r, b);
            assert_eq!(sol1, sol2, "Gr({},{}): structure constants not small integers", r, b);
            for (k, &ti) in tis.iter().enumerate() {
                let mut res: HashMap<usize, CPoly> = HashMap::new();
                for pc in 0..np {
                    let xi = sol1[k][pc];
                    if xi == 0 {
                        continue;
                    }
                    res.entry(pairs[pc].1).or_default().push((pairs[pc].0.clone(), xi));
                }
                let mut v: Vec<(usize, CPoly)> = res.into_iter().collect();
                v.sort_by_key(|x| x.0);
                for (_, cp) in v.iter_mut() {
                    cp.sort();
                }
                results[ti] = v;
            }
        }
        let mut e_mul = vec![vec![Vec::new(); basis.len()]; r];
        let mut f_mul = vec![vec![Vec::new(); basis.len()]; s];
        for (ti, ((kind, j, ai), _)) in targets.iter().enumerate() {
            if *kind == 0 {
                e_mul[j - 1][*ai] = results[ti].clone();
            } else {
                f_mul[j - 1][*ai] = results[ti].clone();
            }
        }
        GrTable { r, b, basis, wt, e_mul, f_mul }
    }

    /// Cached access.
    pub fn get(r: usize, b: usize) -> &'static GrTable {
        static CACHE: OnceLock<Mutex<HashMap<(usize, usize), &'static GrTable>>> = OnceLock::new();
        let m = CACHE.get_or_init(|| Mutex::new(HashMap::new()));
        let mut g = m.lock().unwrap();
        if let Some(t) = g.get(&(r, b)) {
            return t;
        }
        let t: &'static GrTable = Box::leak(Box::new(GrTable::compute(r, b)));
        g.insert((r, b), t);
        t
    }
}


/// Solve sum_pc x[pc] * col[pc] = target (for each target) modulo P, returning
/// symmetric representatives.
fn solve_mod<const P: u64>(
    cols: &[Vec<(usize, i64)>],
    targets: &[Vec<(usize, i64)>],
    nrows: usize,
    r: usize,
    b: usize,
) -> Vec<Vec<i64>> {
    let np = cols.len();
    let nt = targets.len();
    let mut mat = vec![vec![Fp::<P>(0); np + nt]; nrows];
    for (pc, c) in cols.iter().enumerate() {
        for &(i, x) in c {
            mat[i][pc] = Fp::<P>::from_i64(x);
        }
    }
    for (k, c) in targets.iter().enumerate() {
        for &(i, x) in c {
            mat[i][np + k] = Fp::<P>::from_i64(x);
        }
    }
    let piv = rref(&mut mat, np + nt);
    assert!(piv.len() == np && piv.iter().enumerate().all(|(i, &p)| i == p), "Gr({},{}): singular", r, b);
    (0..nt).map(|k| (0..np).map(|pc| mat[pc][np + k].signed()).collect()).collect()
}

fn weight_big(m: &[u16], w: &[u32]) -> u32 {
    m.iter().zip(w.iter()).map(|(&k, &x)| k as u32 * x).sum()
}

/// Binomial coefficient.
pub fn binom(n: usize, k: usize) -> usize {
    if k > n {
        return 0;
    }
    let mut r = 1usize;
    for i in 0..k {
        r = r * (n - i) / (i + 1);
    }
    r
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn table_sizes() {
        for b in 0..=6 {
            for r in 0..=b {
                let t = GrTable::get(r, b);
                assert_eq!(t.fsize(), binom(b, r));
            }
        }
    }
    #[test]
    fn projective_line() {
        // Gr(1,2): basis 1, e1 ; e1*e1 = c1 e1 - c2
        let t = GrTable::get(1, 2);
        assert_eq!(t.basis, vec![vec![0], vec![1]]);
        let r = &t.e_mul[0][1];
        // expect (0, -c2), (1, c1)
        assert_eq!(r.len(), 2);
        assert_eq!(r[0], (0, vec![(vec![0, 1], -1)]));
        assert_eq!(r[1], (1, vec![(vec![1, 0], 1)]));
        // f1 = c1 - e1
        let f = &t.f_mul[0][0];
        assert_eq!(f, &vec![(0, vec![(vec![1, 0], 1)]), (1, vec![(vec![0, 0], -1)])]);
    }
}
