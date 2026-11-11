//! Helpers shared by the integration tests: Poincaré polynomial parsing, the Wang references
//! (tests/data/wang_ref.txt) and the closed Hopf-link formula.
#![allow(dead_code)]

use krhom::laurent::{LPoly, Poincare};
use std::collections::BTreeMap;
use std::sync::OnceLock;

/// "t^0q^-6 + 2t^2q^-8 + ..." (the format of `Poincare`'s Display and of wang.py)
pub fn parse_poincare(s: &str) -> Poincare {
    let mut p = Poincare::default();
    let s = s.trim();
    if s == "0" || s.is_empty() {
        return p;
    }
    for term in s.split(" + ") {
        let term = term.trim();
        let ti = term.find("t^").unwrap_or_else(|| panic!("bad term {:?}", term));
        let k: i128 = if ti == 0 { 1 } else { term[..ti].parse().unwrap() };
        let qi = term.find("q^").unwrap();
        let t: i32 = term[ti + 2..qi].parse().unwrap();
        let q: i32 = term[qi + 2..].parse().unwrap();
        p.add_term(t, q, k);
    }
    p
}

/// tests/data/wang_ref.txt: output of r3scripts/wang.py k N m (T(2,2m+1) colored ∧^k at gl_N,
/// krhom conventions), blocks "k N m = .." followed by the Poincaré polynomial.
pub fn wang_ref() -> &'static BTreeMap<(usize, usize, usize), Poincare> {
    static R: OnceLock<BTreeMap<(usize, usize, usize), Poincare>> = OnceLock::new();
    R.get_or_init(|| {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/wang_ref.txt");
        let text = std::fs::read_to_string(&path).unwrap();
        let mut out = BTreeMap::new();
        let mut lines = text.lines();
        while let Some(l) = lines.next() {
            if let Some(rest) = l.strip_prefix("k N m = ") {
                let v: Vec<usize> = rest.split_whitespace().map(|x| x.parse().unwrap()).collect();
                let poly = parse_poincare(lines.next().unwrap());
                let rank: i128 = lines.next().unwrap().strip_prefix("rank ").unwrap().trim().parse().unwrap();
                assert_eq!(poly.total_rank(), rank, "wang_ref {:?}", v);
                out.insert((v[0], v[1], v[2]), poly);
            }
        }
        out
    })
}

/// R3 V4 (Wang arXiv:2211.08409 Thm 3.5 / Lemma 3.4, converted; `hopf_ours.py`): closure of σ_1^2
/// with colors (a, b), m = min(a, b):
/// Σ_{j=0}^{m − max(0, a+b−N)} t^{2j} q^{−ab − mN + 2m² − j(N+2)} [N; m−j, a−m+j, b−m+j, N−a−b+m−j].
pub fn hopf_formula(nn: usize, a: usize, b: usize) -> Poincare {
    let (n, a, b) = (nn as i64, a as i64, b as i64);
    let m = a.min(b);
    let mut out = Poincare::default();
    if a > n || b > n {
        return out;
    }
    let jmax = m - 0.max(a + b - n);
    for j in 0..=jmax {
        let k = m - j;
        let parts = [k, a - k, b - k, n - a - b + k];
        assert!(parts.iter().all(|&x| x >= 0));
        let mut den = LPoly::one();
        for &x in &parts {
            den = den.mul(&LPoly::qfact(x));
        }
        let mult = LPoly::qfact(n).div_exact(&den);
        let s = -a * b - m * n + 2 * m * m - j * (n + 2);
        for (&e, &c) in &mult.c {
            out.add_term(2 * j as i32, (s + e as i64) as i32, c);
        }
    }
    out
}

