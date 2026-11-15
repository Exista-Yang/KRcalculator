//! Crossing complexes (colored Rickard complexes, ETW Definition 3.3) in ladder form.
//!
//! A braid generator acts on uprights (i, i+1) with current labels (c, d).
//! Positive σ_i: over-strand is the left one (ETW Figure 5). Terms k = 0..min(c,d):
//!   c >= d: W_k = F_i^{(c-d+k)} E_i^{(k)} 1_{(c,d)}   (rungs bottom to top: E(k) then F(c-d+k))
//!   c <  d: W_k = E_i^{(d-c+k)} F_i^{(k)} 1_{(c,d)}   (mirror: F(k) then E(d-c+k))
//! Negative σ_i^{-1}: vertical flip of the positive crossing on (d,c).
//! Homological degree: +k (positive), -k (negative).

use crate::combinat::web::Rung;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Gen {
    pub i: usize,
    pub positive: bool,
}

#[derive(Clone, Debug)]
pub struct CrossingTerm {
    pub k: usize,
    pub rungs: Vec<Rung>,
    pub hdeg: i32,
}

#[derive(Clone, Debug)]
pub struct CrossingCx {
    pub g: Gen,
    pub c: usize,
    pub d: usize,
    pub m: usize,
    /// x = m(N-m), m = min(c,d)
    pub x: i32,
    pub terms: Vec<CrossingTerm>,
}

fn nz(v: Vec<Rung>) -> Vec<Rung> {
    v.into_iter().filter(|r| r.r > 0).collect()
}

pub fn crossing_cx(g: Gen, c: usize, d: usize, n: usize) -> CrossingCx {
    let i = g.i;
    let m = c.min(d);
    let x = (m * (n - m)) as i32;
    let mut terms = vec![];
    for k in 0..=m {
        let rungs = if g.positive {
            if c >= d { nz(vec![Rung::e(i, k), Rung::f(i, c - d + k)]) } else { nz(vec![Rung::f(i, k), Rung::e(i, d - c + k)]) }
        } else {
            // vertical flip of the positive crossing on (d, c)
            if d >= c { nz(vec![Rung::e(i, d - c + k), Rung::f(i, k)]) } else { nz(vec![Rung::f(i, c - d + k), Rung::e(i, k)]) }
        };
        // webs with a label outside 0..=N are zero objects in Foam_N
        let mut lab = vec![c, d];
        let mut valid = true;
        for r in &rungs {
            let rr = Rung { i: 0, dir: r.dir, r: r.r };
            match rr.apply(&lab, n) {
                Some(l) => lab = l,
                None => {
                    valid = false;
                    break;
                }
            }
        }
        if valid {
            terms.push(CrossingTerm { k, rungs, hdeg: if g.positive { k as i32 } else { -(k as i32) } });
        }
    }
    CrossingCx { g, c, d, m, x, terms }
}

/// Grading convention for the q-shift of term k (in the convention where a
/// shift s multiplies graded dimensions by q^s, webs MOY-normalized).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum QConv {
    /// positive: q^{-(x+k)}, negative: q^{x+k}  (KR / Lewark–Lobb direction)
    A,
    /// positive: q^{x+k}, negative: q^{-(x+k)}
    B,
    /// literal ETW (3-2): positive q^{k-x}, negative q^{x-k}
    Etw,
    /// ETW (3-2) read with ETW's shift convention (q^m in ETW = q^{-m} here):
    /// positive q^{x-k}, negative q^{k-x}
    EtwRev,
}

pub fn qshift(cx: &CrossingCx, k: usize, conv: QConv) -> i32 {
    let s = cx.x + k as i32;
    let pos = cx.g.positive;
    match (conv, pos) {
        (QConv::A, true) | (QConv::B, false) => -s,
        (QConv::A, false) | (QConv::B, true) => s,
        (QConv::Etw, true) => k as i32 - cx.x,
        (QConv::Etw, false) => cx.x - k as i32,
        (QConv::EtwRev, true) => cx.x - k as i32,
        (QConv::EtwRev, false) => k as i32 - cx.x,
    }
}
