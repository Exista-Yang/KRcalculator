//! Independent oracle: MOY evaluation of closed ladder webs via quantum skew
//! Howe duality. U_q(gl_m) acts on (∧_q C^m)^{⊗N} (rows j = 1..N), with
//! Δ(E) = E⊗K + 1⊗E and Δ(F) = F⊗1 + K^{-1}⊗F; a basis vector is an N×m 0/1
//! matrix, stored as N row bitmasks. Column i corresponds to upright i, whose
//! label is the column sum. The closure of a ladder is the gl_N quantum trace
//! sum_b <b|W|b> q^{<2ρ, wt(b)>}.

use crate::laurent::LPoly;
use crate::web::{Dir, Ladder};
use std::collections::HashMap;

type State = Vec<u32>;
type Vect = HashMap<State, LPoly>;

fn bit(x: u32, i: usize) -> i32 {
    ((x >> i) & 1) as i32
}

fn apply_e_once(v: &Vect, i: usize, dir: Dir) -> Vect {
    let mut out: Vect = HashMap::new();
    for (s, c) in v {
        let nrows = s.len();
        for j in 0..nrows {
            let row = s[j];
            let (from, to) = match dir {
                Dir::E => (i + 1, i),
                Dir::F => (i, i + 1),
            };
            if bit(row, from) == 1 && bit(row, to) == 0 {
                let mut ns = s.clone();
                ns[j] = (row & !(1 << from)) | (1 << to);
                let h = |x: u32| bit(x, i) - bit(x, i + 1);
                let e: i32 = match dir {
                    Dir::E => (j + 1..nrows).map(|k| h(s[k])).sum(),
                    Dir::F => -(0..j).map(|k| h(s[k])).sum::<i32>(),
                };
                let term = c.shift(e);
                let ent = out.entry(ns).or_insert_with(LPoly::zero);
                *ent = ent.add(&term);
            }
        }
    }
    out.retain(|_, c| !c.is_zero());
    out
}

fn apply_rung(v: &Vect, i: usize, dir: Dir, r: usize) -> Vect {
    let mut w = v.clone();
    for _ in 0..r {
        w = apply_e_once(&w, i, dir);
    }
    let fact = LPoly::qfact(r as i64);
    w.into_iter().map(|(s, c)| (s, c.div_exact(&fact))).filter(|(_, c)| !c.is_zero()).collect()
}

/// All states with given column sums.
pub fn states(n: usize, labels: &[usize]) -> Vec<State> {
    let m = labels.len();
    // choose a subset of rows for each column
    let mut out = vec![vec![0u32; n]];
    for (col, &k) in labels.iter().enumerate() {
        let mut next = Vec::new();
        for s in &out {
            for subset in subsets(n, k) {
                let mut ns = s.clone();
                for j in subset {
                    ns[j] |= 1 << col;
                }
                next.push(ns);
            }
        }
        out = next;
    }
    let _ = m;
    out
}

fn subsets(n: usize, k: usize) -> Vec<Vec<usize>> {
    let mut out = Vec::new();
    let mut cur = Vec::new();
    fn rec(start: usize, n: usize, k: usize, cur: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
        if cur.len() == k {
            out.push(cur.clone());
            return;
        }
        for j in start..n {
            cur.push(j);
            rec(j + 1, n, k, cur, out);
            cur.pop();
        }
    }
    rec(0, n, k, &mut cur, &mut out);
    out
}

/// MOY value of the closure (each upright closed to the right) of a ladder
/// whose top labels equal its bottom labels.
pub fn closure_value(lad: &Ladder) -> LPoly {
    assert_eq!(lad.top(), lad.bottom, "closure needs top == bottom");
    let n = lad.n;
    let mut total = LPoly::zero();
    for s in states(n, &lad.bottom) {
        let mut v: Vect = HashMap::new();
        v.insert(s.clone(), LPoly::one());
        for rg in &lad.rungs {
            v = apply_rung(&v, rg.i, rg.dir, rg.r);
            if v.is_empty() {
                break;
            }
        }
        if let Some(c) = v.get(&s) {
            let w: i32 = (0..n).map(|j| (n as i32 + 1 - 2 * (j as i32 + 1)) * s[j].count_ones() as i32).sum();
            total = total.add(&c.shift(w));
        }
    }
    total
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::web::Rung;
    #[test]
    fn circles_and_thetas() {
        for n in 1..=5usize {
            for a in 0..=n {
                let lad = Ladder::new(n, vec![a]);
                assert_eq!(closure_value(&lad), LPoly::qbinom(n as i64, a as i64));
            }
            // F then E on (1,1): theta = [2][n choose 2]
            if n >= 2 {
                let mut lad = Ladder::new(n, vec![1, 1]);
                lad.push(Rung::f(0, 1)).push(Rung::e(0, 1));
                assert_eq!(closure_value(&lad), LPoly::qint(2).mul(&LPoly::qbinom(n as i64, 2)));
            }
            // (a,0) -> all moved -> back: circle
            for a in 0..=n {
                let mut lad = Ladder::new(n, vec![a, 0]);
                if a > 0 {
                    lad.push(Rung::f(0, a)).push(Rung::e(0, a));
                }
                assert_eq!(closure_value(&lad), LPoly::qbinom(n as i64, a as i64));
            }
        }
    }
    #[test]
    fn ef_relation() {
        // E F 1_(a,b) = F E 1_(a,b) + [a-b] 1_(a,b)  (closures)
        let n = 4;
        for a in 0..=n {
            for b in 0..=n {
                let mut ef = Ladder::new(n, vec![a, b]);
                ef.push(Rung::f(0, 1)).push(Rung::e(0, 1));
                let mut fe = Ladder::new(n, vec![a, b]);
                fe.push(Rung::e(0, 1)).push(Rung::f(0, 1));
                if !ef.is_valid() || !fe.is_valid() {
                    continue;
                }
                let id = Ladder::new(n, vec![a, b]);
                let lhs = closure_value(&ef);
                let rhs = closure_value(&fe).add(&LPoly::qint(a as i64 - b as i64).mul(&closure_value(&id)));
                assert_eq!(lhs, rhs, "a={} b={}", a, b);
            }
        }
    }
}
