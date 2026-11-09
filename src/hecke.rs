//! Hecke algebra of S_n (Soergel's normalization) and Kazhdan–Lusztig basis.
//!
//! (H_s + v)(H_s - v^{-1}) = 0, b_s = H_s + v, b_x = H_x + Σ_{y<x} h_{y,x} H_y with
//! h_{y,x} ∈ v Z[v]. For Soergel bimodules ch(B_x) = b_x and
//! grk_R Hom(B, B') = Σ_w c_w(ch B) c_w(ch B') (bilinear), with v ↔ q-degree 1.

use crate::laurent::LPoly;
use std::collections::HashMap;

#[derive(Clone, Debug)]
pub struct SymGroup {
    pub n: usize,
    /// elements in one-line notation, sorted by length then lexicographically
    pub elems: Vec<Vec<u8>>,
    pub index: HashMap<Vec<u8>, usize>,
    pub len: Vec<usize>,
    /// right multiplication by s_i: rmul[x][i] = index of x·s_i
    pub rmul: Vec<Vec<usize>>,
    /// left multiplication by s_i
    pub lmul: Vec<Vec<usize>>,
}

fn inversions(p: &[u8]) -> usize {
    let mut c = 0;
    for i in 0..p.len() {
        for j in i + 1..p.len() {
            if p[i] > p[j] {
                c += 1;
            }
        }
    }
    c
}

impl SymGroup {
    pub fn new(n: usize) -> Self {
        let mut elems: Vec<Vec<u8>> = vec![];
        let mut cur: Vec<u8> = (0..n as u8).collect();
        fn perms(k: usize, cur: &mut Vec<u8>, out: &mut Vec<Vec<u8>>) {
            if k == cur.len() {
                out.push(cur.clone());
                return;
            }
            for i in k..cur.len() {
                cur.swap(k, i);
                perms(k + 1, cur, out);
                cur.swap(k, i);
            }
        }
        perms(0, &mut cur, &mut elems);
        elems.sort_by_key(|p| (inversions(p), p.clone()));
        let index: HashMap<Vec<u8>, usize> = elems.iter().cloned().enumerate().map(|(i, p)| (p, i)).collect();
        let len: Vec<usize> = elems.iter().map(|p| inversions(p)).collect();
        let nn = n.saturating_sub(1);
        let rmul = elems
            .iter()
            .map(|p| {
                (0..nn)
                    .map(|i| {
                        let mut q = p.clone();
                        q.swap(i, i + 1);
                        index[&q]
                    })
                    .collect()
            })
            .collect();
        let lmul = elems
            .iter()
            .map(|p| {
                (0..nn)
                    .map(|i| {
                        // s_i · p : swap values i and i+1
                        let q: Vec<u8> = p.iter().map(|&v| if v == i as u8 { i as u8 + 1 } else if v == i as u8 + 1 { i as u8 } else { v }).collect();
                        index[&q]
                    })
                    .collect()
            })
            .collect();
        SymGroup { n, elems, index, len, rmul, lmul }
    }
    pub fn size(&self) -> usize {
        self.elems.len()
    }
    pub fn identity(&self) -> usize {
        0
    }
    /// A reduced expression (list of simple reflection indices) for x, lexicographically
    /// smallest in the sense of always peeling the smallest right descent.
    pub fn rex(&self, x: usize) -> Vec<usize> {
        let mut w = x;
        let mut out = vec![];
        while self.len[w] > 0 {
            let i = (0..self.n - 1).find(|&i| self.len[self.rmul[w][i]] < self.len[w]).unwrap();
            out.push(i);
            w = self.rmul[w][i];
        }
        out.reverse();
        out
    }
    pub fn from_word(&self, word: &[usize]) -> usize {
        let mut w = 0;
        for &i in word {
            w = self.rmul[w][i];
        }
        w
    }
    pub fn right_descent(&self, x: usize, i: usize) -> bool {
        self.len[self.rmul[x][i]] < self.len[x]
    }
    /// Bruhat order y ≤ x (subword criterion via recursion).
    pub fn bruhat_le(&self, y: usize, x: usize) -> bool {
        if self.len[y] > self.len[x] {
            return false;
        }
        if self.len[x] == 0 {
            return y == 0;
        }
        let i = (0..self.n - 1).find(|&i| self.right_descent(x, i)).unwrap();
        let xs = self.rmul[x][i];
        if self.right_descent(y, i) {
            self.bruhat_le(self.rmul[y][i], xs)
        } else {
            self.bruhat_le(y, xs)
        }
    }
}

/// Element of the Hecke algebra in the standard basis.
pub type HElem = Vec<LPoly>;

pub struct Hecke {
    pub g: SymGroup,
    /// KL basis elements in the standard basis
    pub kl: Vec<HElem>,
    /// mu(y, x)
    pub mu: HashMap<(usize, usize), i128>,
}

impl Hecke {
    pub fn zero(&self) -> HElem {
        vec![LPoly::zero(); self.g.size()]
    }
    /// h · b_s (right multiplication)
    pub fn mul_bs(&self, h: &HElem, s: usize) -> HElem {
        let mut out = self.zero();
        let v = LPoly::mono(1, 1);
        let vinv = LPoly::mono(-1, 1);
        for (y, c) in h.iter().enumerate() {
            if c.is_zero() {
                continue;
            }
            let ys = self.g.rmul[y][s];
            out[ys] = out[ys].add(c);
            if self.g.len[ys] > self.g.len[y] {
                out[y] = out[y].add(&c.mul(&v));
            } else {
                out[y] = out[y].add(&c.mul(&vinv));
            }
        }
        out
    }
    pub fn new(n: usize) -> Self {
        let g = SymGroup::new(n);
        let size = g.size();
        let mut kl: Vec<HElem> = vec![vec![]; size];
        let mut e = vec![LPoly::zero(); size];
        e[0] = LPoly::one();
        kl[0] = e;
        let mut mu: HashMap<(usize, usize), i128> = HashMap::new();
        let mut h = Hecke { g, kl: vec![], mu: HashMap::new() };
        // elements are sorted by length; build b_x from b_{x'} b_s with x = x's > x'
        for x in 1..size {
            let i = (0..n - 1).find(|&i| h.g.right_descent(x, i)).unwrap();
            let xp = h.g.rmul[x][i];
            let mut b = h.mul_bs(&kl[xp], i);
            // subtract mu(y, x') b_y for y < x', y s < y
            for y in 0..size {
                if h.g.len[y] >= h.g.len[xp] || kl[y].is_empty() {
                    continue;
                }
                if !h.g.right_descent(y, i) {
                    continue;
                }
                let m = *mu.get(&(y, xp)).unwrap_or(&0);
                if m != 0 {
                    for z in 0..size {
                        if !kl[y][z].is_zero() {
                            b[z] = b[z].sub(&kl[y][z].scale(m));
                        }
                    }
                }
            }
            // record mu(z, x) = coefficient of v^1 in h_{z,x}
            for z in 0..size {
                if z != x {
                    let c = b[z].c.get(&1).copied().unwrap_or(0);
                    if c != 0 {
                        mu.insert((z, x), c);
                    }
                }
            }
            kl[x] = b;
        }
        h.kl = kl;
        h.mu = mu;
        h
    }
    /// KL polynomial P_{y,x}(q) in Soergel's normalization h_{y,x} = v^{l(x)-l(y)} P_{y,x}(v^{-2}).
    pub fn h(&self, y: usize, x: usize) -> LPoly {
        self.kl[x][y].clone()
    }
    /// Character of a Bott–Samelson object b_{s1}...b_{sk}.
    pub fn bs_char(&self, word: &[usize]) -> HElem {
        let mut e = self.zero();
        e[0] = LPoly::one();
        for &s in word {
            e = self.mul_bs(&e, s);
        }
        e
    }
    /// Graded rank (in v) of Hom(B, B') given characters.
    pub fn hom_grk(&self, a: &HElem, b: &HElem) -> LPoly {
        let mut r = LPoly::zero();
        for (x, c) in a.iter().enumerate() {
            if !c.is_zero() && !b[x].is_zero() {
                r = r.add(&c.mul(&b[x]));
            }
        }
        r
    }
    /// Expand an element in the KL basis (assumes it is a Z[v,v^-1]-combination of b_x).
    pub fn to_kl(&self, h: &HElem) -> Vec<LPoly> {
        let mut rem = h.clone();
        let size = self.g.size();
        let mut out = vec![LPoly::zero(); size];
        for x in (0..size).rev() {
            if rem[x].is_zero() {
                continue;
            }
            let c = rem[x].clone();
            out[x] = c.clone();
            for z in 0..size {
                if !self.kl[x][z].is_zero() {
                    rem[z] = rem[z].sub(&self.kl[x][z].mul(&c));
                }
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn s3_s4_kl() {
        let h3 = Hecke::new(3);
        assert_eq!(h3.g.size(), 6);
        // in S_3 all KL polynomials are trivial: h_{y,x} = v^{l(x)-l(y)} for y ≤ x
        for x in 0..6 {
            for y in 0..6 {
                let le = h3.g.bruhat_le(y, x);
                let expect = if le { LPoly::mono((h3.g.len[x] - h3.g.len[y]) as i32, 1) } else { LPoly::zero() };
                assert_eq!(h3.h(y, x), expect, "S3 h({},{})", y, x);
            }
        }
        let h4 = Hecke::new(4);
        // 3412 and 4231 have P_{e,x} = 1 + q, i.e. h_{e,x} = v^4 + v^2
        for p in [vec![2u8, 3, 0, 1], vec![3u8, 1, 2, 0]] {
            let x = h4.g.index[&p];
            let l = h4.g.len[x] as i32;
            assert_eq!(h4.h(0, x), LPoly::mono(l, 1).add(&LPoly::mono(l - 2, 1)), "{:?}", p);
        }
        // b_x b_s expands positively in KL basis
        for x in 0..24 {
            for s in 0..3 {
                let p = h4.mul_bs(&h4.kl[x], s);
                let k = h4.to_kl(&p);
                for c in &k {
                    assert!(c.c.values().all(|&v| v >= 0));
                }
            }
        }
    }
}
