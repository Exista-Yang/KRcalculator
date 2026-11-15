//! Colored (∧^a) braid closures via cabling: each a-colored strand becomes a
//! parallel strands, a colored crossing becomes the a·b cabled crossing, and a
//! thick-edge insertion B_{w_0(S_a)} is placed at the bottom of each colored block.
//! The resulting homology is the colored homology times Π [a]! over the colored strand
//! positions (one factor per thick-edge insertion, as q-graded spaces in homological degree
//! 0), shifted by `cable_crossing_shift` per crossing (no per-strand shift); see
//! `khr_colored_cable`.

use crate::input::crossing::Gen;

#[derive(Clone, Debug)]
pub struct Cabled {
    pub strands: usize,
    pub word: Vec<Gen>,
    /// start object as a word in simple reflections: product of w_0 of each block
    pub start: Vec<usize>,
    /// a-values of each colored strand (bottom positions)
    pub labels: Vec<usize>,
}

/// Longest element of S_a on positions off..off+a as a reduced word.
fn w0_word(off: usize, a: usize) -> Vec<usize> {
    let mut w = vec![];
    for i in 0..a {
        for j in (0..a - 1 - i).rev() {
            let _ = j;
        }
        let _ = i;
    }
    // standard reduced word: (s_1)(s_2 s_1)(s_3 s_2 s_1)...
    for k in 1..a {
        for j in (0..k).rev() {
            w.push(off + j);
        }
    }
    w
}

pub fn cable(labels: &[usize], word: &[Gen]) -> Cabled {
    let mut lab = labels.to_vec();
    let mut out = vec![];
    for g in word {
        let i = g.i;
        let off: usize = lab[..i].iter().sum();
        let (a, b) = (lab[i], lab[i + 1]);
        // the a-block at off..off+a crosses the b-block at off+a..off+a+b.
        // Positive σ: the left (a) block goes over. As a braid: move each strand of the
        // left block rightwards past all strands of the right block.
        // cabled word: for k = a-1 down to 0: strand at position off+k moves right b steps
        for k in (0..a).rev() {
            for j in 0..b {
                out.push(Gen { i: off + k + j, positive: g.positive });
            }
        }
        lab.swap(i, i + 1);
    }
    let mut start = vec![];
    let mut off = 0;
    for &a in labels {
        start.extend(w0_word(off, a));
        off += a;
    }
    Cabled { strands: labels.iter().sum(), word: out, start, labels: labels.to_vec() }
}

// ------------------------------------------------------------------ cabling as a colored engine

use crate::input::crossing::{crossing_cx, qshift, QConv};
use crate::uncolored::engine2::{khr_engine2_closure, Closure};
use crate::algebra::field::Field;
use crate::algebra::laurent::{LPoly, Poincare};
use std::collections::BTreeMap;

/// Labels (left, right) of the two uprights at every crossing of `word`, starting from the
/// bottom `labels`; errors unless the top labels equal the bottom labels (braid closure).
pub fn crossing_labels(labels: &[usize], word: &[Gen]) -> Result<Vec<(usize, usize)>, String> {
    let mut lab = labels.to_vec();
    let mut out = Vec::with_capacity(word.len());
    for g in word {
        if g.i + 1 >= lab.len() {
            return Err(format!("generator {} needs at least {} strands", g.i + 1, g.i + 2));
        }
        out.push((lab[g.i], lab[g.i + 1]));
        lab.swap(g.i, g.i + 1);
    }
    if lab != labels {
        return Err(format!("colors do not close up: bottom {:?}, top {:?}", labels, lab));
    }
    Ok(out)
}

/// (t, q) shift of one cabled crossing against the colored crossing, for labels (a, b), in the
/// conventions of `cube::khr_full_cube` with `QConv::A` and of the uncolored Soergel engine:
/// with m = min(a,b) and T = ab - m, the cabled crossing (thick-edge insertions B_{w0} on
/// both blocks) has the homology of the colored crossing times t^T q^{-(N T + m(m-1))}
/// (positive) or t^{-T} q^{N T + m(m-1)} (negative). Heuristic origin: colored term k sits in
/// cable degree T + k, and the uncolored shifts q^{-(N-1)} per crossing and q^{-1} per
/// resolved crossing differ from the colored q^{-(m(N-m)+k)} by exactly this amount; the same
/// q-power is the ratio of the two braidings on the Cartan component of ∧^a ⊗ ∧^b (content
/// formula), a decategorified sketch only. Σ T over a closed braid is always even (crossings
/// between two components come in pairs, T = a(a-1) on one component). Verified against the
/// full cube for (a,b) = (1,1), (2,1), (1,2), (2,2), (3,1) and both signs (tests/acceptance.rs);
/// (3,2) and m = 3 are unverified (cables too large). `khr_colored_cable` re-checks the range.
pub fn cable_crossing_shift(a: usize, b: usize, positive: bool, nn: usize) -> (i32, i32) {
    let m = a.min(b) as i32;
    let t = (a * b) as i32 - m;
    let q = nn as i32 * t + m * (m - 1);
    if positive { (t, -q) } else { (-t, q) }
}

/// Sum of `cable_crossing_shift` over the crossings: cable homology / Π[a]! = colored · t^t q^q.
pub fn cable_total_shift(labels: &[usize], word: &[Gen], nn: usize) -> Result<(i32, i32), String> {
    let mut s = (0, 0);
    for (g, (a, b)) in word.iter().zip(crossing_labels(labels, word)?) {
        let (t, q) = cable_crossing_shift(a, b, g.positive, nn);
        s.0 += t;
        s.1 += q;
    }
    Ok(s)
}

/// Π over the colored strands (bottom positions, i.e. thick-edge insertions) of [a]!.
pub fn cable_divisor(labels: &[usize]) -> LPoly {
    labels.iter().fold(LPoly::one(), |acc, &a| acc.mul(&LPoly::qfact(a as i64)))
}

/// Exact division of every homological degree of `p` by `d`; errors if some degree is not
/// divisible or the quotient has a negative coefficient.
pub fn divide_poincare(p: &Poincare, d: &LPoly) -> Result<Poincare, String> {
    let mut by_t: BTreeMap<i32, LPoly> = BTreeMap::new();
    for (&(t, q), &k) in &p.c {
        by_t.entry(t).or_insert_with(LPoly::zero).add_term(q, k);
    }
    let mut out = Poincare::default();
    for (t, lp) in by_t {
        let quo = div_checked(&lp, d).ok_or_else(|| format!("homological degree {} is not divisible by {}", t, d))?;
        for (&q, &k) in &quo.c {
            if k < 0 {
                return Err(format!("negative coefficient at t^{} q^{} after division", t, q));
            }
            out.add_term(t, q, k);
        }
    }
    Ok(out)
}

/// Exact division p / d (d with leading coefficient ±1), or None if d does not divide p.
fn div_checked(p: &LPoly, d: &LPoly) -> Option<LPoly> {
    let (&dmax, &dlead) = d.c.iter().next_back()?;
    if dlead != 1 && dlead != -1 {
        return None;
    }
    let dmin = d.min_deg()?;
    let mut rem = p.clone();
    let mut quo = LPoly::zero();
    loop {
        let Some((&rmax, &rlead)) = rem.c.iter().next_back() else { break };
        let rmin = rem.min_deg()?;
        if rmax - dmax < rmin - dmin {
            return None;
        }
        let (e, k) = (rmax - dmax, rlead * dlead);
        quo.add_term(e, k);
        for (&de, &dk) in &d.c {
            rem.add_term(de + e, -dk * k);
        }
    }
    Some(quo)
}

/// p · t^dt q^dq.
pub fn shift_poincare(p: &Poincare, dt: i32, dq: i32) -> Poincare {
    let mut out = Poincare::default();
    for (&(t, q), &k) in &p.c {
        out.add_term(t + dt, q + dq, k);
    }
    out
}

/// Homological range [lo, hi] of the colored chain complex of the closure in gl_N (only the
/// Rickard terms with all labels in 0..=N count, as in the full cube).
pub fn colored_hdeg_range(labels: &[usize], word: &[Gen], nn: usize) -> Result<(i32, i32), String> {
    let (mut lo, mut hi) = (0, 0);
    for (g, (a, b)) in word.iter().zip(crossing_labels(labels, word)?) {
        let cx = crossing_cx(*g, a, b, nn);
        lo += cx.terms.iter().map(|t| t.hdeg).min().unwrap_or(0);
        hi += cx.terms.iter().map(|t| t.hdeg).max().unwrap_or(0);
    }
    Ok((lo, hi))
}

/// Graded Euler characteristic of the colored closure from the MOY evaluation of every
/// resolution, Σ_states (-1)^h q^{Σ shifts} <closure of the web> (`moy::closure_value`; the
/// same q-shifts as `khr_full_cube` with `QConv::A`). Independent of both homology engines.
/// `None` if the number of resolutions exceeds `max_states`.
pub fn colored_euler_moy(labels: &[usize], word: &[Gen], nn: usize, max_states: usize) -> Result<Option<LPoly>, String> {
    use crate::combinat::web::Ladder;
    let labs = crossing_labels(labels, word)?;
    if labels.iter().any(|&a| a > nn) {
        return Ok(Some(LPoly::zero()));
    }
    let cxs: Vec<_> = word.iter().zip(&labs).map(|(g, &(a, b))| crossing_cx(*g, a, b, nn)).collect();
    let mut count: usize = 1;
    for cx in &cxs {
        count = count.saturating_mul(cx.terms.len());
        if count > max_states {
            return Ok(None);
        }
    }
    // depth-first over the resolutions, sharing rung prefixes
    let mut total = LPoly::zero();
    let mut stack: Vec<(usize, Vec<crate::combinat::web::Rung>, i32, i32)> = vec![(0, vec![], 0, 0)];
    while let Some((l, rungs, h, q)) = stack.pop() {
        if l == cxs.len() {
            let lad = Ladder { n: nn, bottom: labels.to_vec(), rungs };
            let v = crate::combinat::moy::closure_value(&lad).shift(q);
            total = total.add(&if h.rem_euclid(2) == 0 { v } else { v.scale(-1) });
            continue;
        }
        for t in &cxs[l].terms {
            let mut r = rungs.clone();
            r.extend(t.rungs.iter().copied());
            stack.push((l + 1, r, h + t.hdeg, q + qshift(&cxs[l], t.k, QConv::A)));
        }
    }
    Ok(Some(total))
}

/// Whether cabling is the default (affordable) method for these strand labels. Measured
/// (Apple M5, KRHOM_THREADS=2, 2026-11-09): cables with at most 4 strands, and 5-strand cables
/// whose colors are at most 2, finish in seconds; the (3,2) Hopf link (5 strands, start object
/// B_{w0(S_3)} s_4) reached 6 GB without finishing, and 6-strand cables (3-colored single
/// crossing, 2-colored figure-eight) reached 0.4-5 GB within 12-60 s without finishing the scan.
pub fn cable_default(labels: &[usize]) -> bool {
    let strands: usize = labels.iter().sum();
    strands <= 4 || (strands == 5 && labels.iter().all(|&a| a <= 2))
}

/// Colored homology of a braid closure computed by cabling.
pub struct ColoredCable {
    pub cabled: Cabled,
    /// homology of the cable (uncolored engine normalization), for the N that were computed
    pub raw: BTreeMap<usize, Poincare>,
    /// colored homology in krhom's conventions (equal to `khr_full_cube` with `QConv::A`)
    pub colored: BTreeMap<usize, Poincare>,
    /// engine2 timings (scan, data preparation, closure of all N), seconds
    pub timings: (f64, f64, f64),
    pub final_objects: usize,
}

/// Colored gl_N homology of the closure of `word` with bottom `labels` (one ∧^a label per
/// strand position, constant along components) for every N in `ns`, by cabling: the cable of
/// `cable` with one B_{w0(S_a)} per colored strand is run through the uncolored Soergel engine
/// (one tangle-level scan for all N, then `closure`), every homological degree is divided by
/// Π[a]! and the per-crossing shifts `cable_crossing_shift` are undone. N < max label gives 0
/// (∧^a C^N = 0) without running the engine for that N. Errors (instead of returning a wrong
/// answer) if a degree is not divisible by Π[a]! or the result leaves the homological range of
/// the colored complex; both would indicate a wrong normalization.
pub fn khr_colored_cable<F: Field>(labels: &[usize], word: &[Gen], ns: &[usize], seed: u64, closure: Closure) -> Result<ColoredCable, String> {
    if labels.iter().any(|&a| a == 0) {
        return Err("colors must be at least 1".into());
    }
    crossing_labels(labels, word)?;
    let cb = cable(labels, word);
    let amax = labels.iter().copied().max().unwrap_or(1);
    let run: Vec<usize> = ns.iter().copied().filter(|&nn| nn >= amax).collect();
    let mut raw = BTreeMap::new();
    let mut timings = (0.0, 0.0, 0.0);
    let mut final_objects = 0;
    if !run.is_empty() {
        let r = khr_engine2_closure::<F>(cb.strands, &cb.word, &run, &cb.start, seed, false, closure);
        raw = r.poincare;
        timings = r.timings;
        final_objects = r.final_objects;
    }
    let div = cable_divisor(labels);
    let mut colored = BTreeMap::new();
    for &nn in ns {
        let Some(p) = raw.get(&nn) else {
            colored.insert(nn, Poincare::default());
            continue;
        };
        let (dt, dq) = cable_total_shift(labels, word, nn)?;
        let c = shift_poincare(&divide_poincare(p, &div).map_err(|e| format!("N={}: {}", nn, e))?, -dt, -dq);
        let (lo, hi) = colored_hdeg_range(labels, word, nn)?;
        if let Some(&(t, q)) = c.c.keys().find(|&&(t, _)| t < lo || t > hi) {
            return Err(format!("N={}: t^{} q^{} lies outside the colored complex range [{}, {}]", nn, t, q, lo, hi));
        }
        colored.insert(nn, c);
    }
    Ok(ColoredCable { cabled: cb, raw, colored, timings, final_objects })
}
