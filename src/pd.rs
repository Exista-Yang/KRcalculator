//! Planar diagram (PD) codes and their conversion to braid words (Vogel's algorithm).
//!
//! # PD convention (KnotTheory / Knot Atlas)
//! A crossing `X[i,j,k,l]` lists the four edge labels **counterclockwise**, starting
//! from the incoming under-strand: the under-strand runs `i -> k`, the over-strand
//! joins `j` and `l`. Edges are oriented along increasing labels, with wraparound per
//! component. Here the orientation of each component is obtained by walking along it;
//! its direction is fixed by the under-passages (`i -> k` is unambiguous). Only for a
//! component that never passes under (necessarily a split unknot) the label order is
//! used.
//!
//! Sign: with the under-strand pointing north (`i` = S, `j` = E, `k` = N, `l` = W),
//! the crossing is positive (right-handed) iff the over-strand runs `l -> j`
//! (west to east). This is KnotTheory's `PositiveQ`: `j - l == 1 || l - j > 1`
//! (and the kinks `i == j`, `k == l` are positive).
//!
//! # Braid convention ([`crate::crossing::Gen`])
//! Strands oriented upward, word read bottom to top; `Gen { i, positive: true }` is a
//! positive crossing in which the strand at bottom position `i` (the left one) passes
//! over the strand at position `i + 1`.
//!
//! # Algorithm
//! P. Vogel, *Representation of links by braids: a new algorithm*, Comment. Math.
//! Helv. 65 (1990). While some face of the diagram contains two edges that belong to
//! different Seifert circles and both have the face on the same side, a Reidemeister II
//! move pushes one edge over the other across that face. This preserves the number of
//! Seifert circles and the writhe; when no such face is left, the Seifert circles are
//! coherently nested and the braid is read off around the common axis. The braid then
//! has as many strands as the input diagram has Seifert circles. Afterwards
//! [`BraidWord::simplify`] applies cyclic free reductions, braid relations and Markov
//! destabilizations.

use crate::crossing::Gen;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

/// A PD code in the KnotTheory convention (see the module documentation).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Pd {
    /// `X[i,j,k,l]`: edge labels counterclockwise from the incoming under-strand.
    pub crossings: Vec<[usize; 4]>,
}

/// A braid word together with the link component of every strand.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BraidWord {
    pub strands: usize,
    pub word: Vec<Gen>,
    /// `components[p]` = link component of the strand at bottom position `p`.
    /// For braids obtained from a PD code, components are numbered as in the PD code
    /// (the component containing the smallest edge label is component 0, ...).
    pub components: Vec<usize>,
}

impl fmt::Display for Pd {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "PD[")?;
        for (n, x) in self.crossings.iter().enumerate() {
            if n > 0 {
                write!(f, ", ")?;
            }
            write!(f, "X[{},{},{},{}]", x[0], x[1], x[2], x[3])?;
        }
        write!(f, "]")
    }
}

/// Parse a PD code. Accepted forms include `[[1,4,2,5],[3,6,4,1],[5,2,6,3]]`,
/// `X[1,4,2,5] X[3,6,4,1] X[5,2,6,3]` and `PD[X[1,4,2,5], ...]`: all integers are read
/// in order and grouped into crossings of four. The code is validated (every label
/// occurs exactly twice, consistent orientations, planarity).
pub fn parse_pd(s: &str) -> Result<Pd, String> {
    let mut nums: Vec<usize> = vec![];
    let mut cur: Option<usize> = None;
    for ch in s.chars() {
        if let Some(d) = ch.to_digit(10) {
            let v = cur.unwrap_or(0).checked_mul(10).and_then(|v| v.checked_add(d as usize));
            cur = Some(v.ok_or_else(|| "edge label too large".to_string())?);
        } else {
            if let Some(v) = cur.take() {
                nums.push(v);
            }
            if !(ch.is_whitespace() || "[](){},;_XxPDpd".contains(ch)) {
                return Err(format!("unexpected character {ch:?} in PD code"));
            }
        }
    }
    if let Some(v) = cur {
        nums.push(v);
    }
    if nums.len() % 4 != 0 {
        return Err(format!("PD code has {} labels, not a multiple of 4", nums.len()));
    }
    let crossings = nums.chunks(4).map(|c| [c[0], c[1], c[2], c[3]]).collect();
    let pd = Pd { crossings };
    pd.validate()?;
    Ok(pd)
}

impl Pd {
    /// Check that the code describes an oriented planar link diagram.
    pub fn validate(&self) -> Result<(), String> {
        Diagram::from_pd(self).map(|_| ())
    }

    /// Number of link components (0 for the empty code).
    pub fn num_components(&self) -> usize {
        self.diagram().ncomp
    }

    /// Signs of the crossings (`true` = positive / right-handed).
    pub fn crossing_signs(&self) -> Vec<bool> {
        self.diagram().pos
    }

    pub fn writhe(&self) -> i64 {
        self.crossing_signs().iter().map(|&p| if p { 1 } else { -1 }).sum()
    }

    /// Link component of every edge label (components numbered by smallest label).
    pub fn label_components(&self) -> BTreeMap<usize, usize> {
        let d = self.diagram();
        let mut m = BTreeMap::new();
        for (c, x) in self.crossings.iter().enumerate() {
            for s in 0..4 {
                m.insert(x[s], d.comp[d.x[c][s]]);
            }
        }
        m
    }

    /// Number of Seifert circles of the diagram (= strands produced by Vogel's algorithm).
    pub fn seifert_circles(&self) -> usize {
        if self.crossings.is_empty() {
            return 1;
        }
        self.diagram().seifert_circles().1.len()
    }

    fn diagram(&self) -> Diagram {
        Diagram::from_pd(self).unwrap_or_else(|e| panic!("invalid PD code: {e}"))
    }
}

/// Convert a PD code to a braid whose closure is the same oriented link:
/// Vogel's algorithm followed by [`BraidWord::simplify`].
pub fn pd_to_braid(pd: &Pd) -> BraidWord {
    pd_to_braid_vogel(pd).simplify()
}

/// Vogel's algorithm only (no simplification). The number of strands equals the number
/// of Seifert circles of the diagram and the writhe equals that of the diagram.
pub fn pd_to_braid_vogel(pd: &Pd) -> BraidWord {
    let d = pd.diagram();
    if d.x.is_empty() {
        return BraidWord { strands: 1, word: vec![], components: vec![0] };
    }
    // split diagrams: braid each connected piece and juxtapose the braids
    let mut out = BraidWord { strands: 0, word: vec![], components: vec![] };
    for mut piece in d.connected_pieces() {
        piece.vogel();
        let b = piece.read_braid();
        let off = out.strands;
        out.word.extend(b.word.iter().map(|g| Gen { i: g.i + off, positive: g.positive }));
        out.components.extend(b.components);
        out.strands += b.strands;
    }
    out
}

/// PD code of the closure of a braid (closing arcs run around the right-hand side).
/// Edges are labelled 1, 2, ... consecutively along each component in the direction of
/// the orientation; components are numbered in order of their first bottom position.
/// Panics if some strand position is never touched by a crossing (such a split unknot
/// cannot be expressed as a PD code).
pub fn braid_to_pd(b: &BraidWord) -> Pd {
    let n = b.strands;
    let mut next_id = n;
    let mut cur: Vec<usize> = (0..n).collect();
    let mut xs: Vec<[usize; 4]> = vec![];
    let mut pos: Vec<bool> = vec![];
    for g in &b.word {
        assert!(g.i + 1 < n, "generator {} out of range for {} strands", g.i, n);
        let (l, r) = (cur[g.i], cur[g.i + 1]);
        let (ne, nw) = (next_id, next_id + 1);
        next_id += 2;
        // SW = l, SE = r, NE = ne (continues l), NW = nw (continues r)
        if g.positive {
            // left over right: under strand SE -> NW
            xs.push([r, ne, nw, l]);
        } else {
            // right over left: under strand SW -> NE
            xs.push([l, r, ne, nw]);
        }
        pos.push(g.positive);
        cur[g.i] = nw;
        cur[g.i + 1] = ne;
    }
    let mut rep: Vec<usize> = (0..next_id).collect();
    for p in 0..n {
        assert!(cur[p] != p, "braid_to_pd: strand position {p} has no crossings");
        rep[cur[p]] = p;
    }
    let xs: Vec<[usize; 4]> = xs.iter().map(|x| x.map(|e| rep[e])).collect();
    let mut next = vec![usize::MAX; next_id];
    for (c, x) in xs.iter().enumerate() {
        next[x[0]] = x[2];
        if pos[c] {
            next[x[3]] = x[1];
        } else {
            next[x[1]] = x[3];
        }
    }
    let mut label = vec![0usize; next_id];
    let mut lab = 0;
    for p in 0..n {
        if label[p] != 0 {
            continue;
        }
        let mut e = p;
        loop {
            lab += 1;
            label[e] = lab;
            e = next[e];
            if e == p {
                break;
            }
        }
    }
    Pd { crossings: xs.iter().map(|x| x.map(|e| label[e])).collect() }
}

impl BraidWord {
    /// Braid from signed 1-based generators (`k > 0`: σ_k, `k < 0`: σ_|k|^{-1}).
    /// Components are the cycles of the permutation, numbered by smallest bottom position.
    pub fn from_signed(strands: usize, w: &[i64]) -> BraidWord {
        let word: Vec<Gen> = w
            .iter()
            .map(|&k| {
                assert!(k != 0 && (k.unsigned_abs() as usize) < strands, "bad generator {k}");
                Gen { i: k.unsigned_abs() as usize - 1, positive: k > 0 }
            })
            .collect();
        let mut b = BraidWord { strands, word, components: vec![0; strands] };
        let perm = b.permutation();
        let mut comp = vec![usize::MAX; strands];
        let mut nc = 0;
        for p in 0..strands {
            if comp[p] == usize::MAX {
                let mut q = p;
                while comp[q] == usize::MAX {
                    comp[q] = nc;
                    q = perm[q];
                }
                nc += 1;
            }
        }
        b.components = comp;
        b
    }

    /// Signed 1-based generators.
    pub fn signed(&self) -> Vec<i64> {
        self.word.iter().map(|g| if g.positive { g.i as i64 + 1 } else { -(g.i as i64 + 1) }).collect()
    }

    /// `perm[p]` = top position of the strand that starts at bottom position `p`.
    pub fn permutation(&self) -> Vec<usize> {
        let mut at: Vec<usize> = (0..self.strands).collect();
        for g in &self.word {
            at.swap(g.i, g.i + 1);
        }
        let mut perm = vec![0; self.strands];
        for (pos, &s) in at.iter().enumerate() {
            perm[s] = pos;
        }
        perm
    }

    /// Components of the strands at the top (must equal `components` for a closure).
    pub fn top_components(&self) -> Vec<usize> {
        let mut c = self.components.clone();
        for g in &self.word {
            c.swap(g.i, g.i + 1);
        }
        c
    }

    pub fn writhe(&self) -> i64 {
        self.word.iter().map(|g| if g.positive { 1 } else { -1 }).sum()
    }

    /// Number of components of the closure (cycles of the permutation).
    pub fn num_components(&self) -> usize {
        let perm = self.permutation();
        let mut seen = vec![false; self.strands];
        let mut nc = 0;
        for p in 0..self.strands {
            if !seen[p] {
                nc += 1;
                let mut q = p;
                while !seen[q] {
                    seen[q] = true;
                    q = perm[q];
                }
            }
        }
        nc
    }

    /// Reduce strands and length while keeping the closure (as an oriented link, with
    /// the component labels carried along): cyclic free reduction (also across
    /// commuting letters), Markov destabilization at either end, and braid relations
    /// that lower the number of occurrences of an end generator to one; then a seeded
    /// (deterministic) random search, see [`BraidWord::simplify_search`].
    pub fn simplify(&self) -> BraidWord {
        self.simplify_greedy().simplify_search(20_000)
    }

    /// The deterministic part of [`BraidWord::simplify`].
    pub fn simplify_greedy(&self) -> BraidWord {
        let mut w = Work { n: self.strands, w: self.word.clone(), comps: self.components.clone() };
        w.greedy();
        BraidWord { strands: w.n, word: w.w, components: w.comps }
    }

    /// Random walk on braid words representing the same closure: conjugation (rotation
    /// and by single generators), commutations, braid relations of length three and
    /// insertion of cancelling pairs, accepted by a Metropolis rule that favours few
    /// occurrences of an end generator and short words; Markov destabilizations are
    /// applied whenever possible. Stops after `patience` steps without improvement.
    /// Deterministic (fixed seed). Returns the best braid found (fewest strands, then
    /// shortest).
    pub fn simplify_search(&self, patience: usize) -> BraidWord {
        let mut cur = Work { n: self.strands, w: self.word.clone(), comps: self.components.clone() };
        cur.greedy();
        let ncomp = BraidWord { strands: cur.n, word: cur.w.clone(), components: cur.comps.clone() }.num_components();
        let mut best = cur.clone();
        let mut rng = XorShift(0x2545_F491_4F6C_DD1D);
        let max_len = 2 * cur.w.len() + 12;
        let temp = 40.0;
        let mut since = 0;
        while since < patience && best.n > ncomp.max(1) && cur.n >= 2 {
            since += 1;
            let len = cur.w.len();
            let mut cand = cur.clone();
            let mv = rng.below(100);
            let ok = match mv {
                0..=34 if len >= 3 => {
                    cand.rotate(rng.below(len as u64) as usize);
                    cand.relation_at_front()
                }
                35..=64 if len >= 2 => {
                    cand.rotate(rng.below(len as u64) as usize);
                    if cand.w[0].i.abs_diff(cand.w[1].i) >= 2 {
                        cand.w.swap(0, 1);
                        true
                    } else {
                        false
                    }
                }
                65..=79 if len >= 1 => {
                    cand.rotate(1 + rng.below(len as u64) as usize);
                    true
                }
                80..=94 => {
                    let g = Gen { i: rng.below(cand.n as u64 - 1) as usize, positive: rng.below(2) == 0 };
                    cand.comps.swap(g.i, g.i + 1);
                    cand.w.insert(0, g);
                    cand.w.push(inv(g));
                    true
                }
                _ => {
                    let g = Gen { i: rng.below(cand.n as u64 - 1) as usize, positive: rng.below(2) == 0 };
                    let k = rng.below(len as u64 + 1) as usize;
                    cand.w.insert(k, g);
                    cand.w.insert(k + 1, inv(g));
                    true
                }
            };
            if !ok {
                continue;
            }
            cand.free_reduce_linear();
            cand.destab_both();
            if cand.w.len() > max_len {
                continue;
            }
            let (fc, fo) = (cand.score(), cur.score());
            if fc <= fo || (rng.below(1 << 20) as f64) < ((fo - fc) / temp).exp() * (1u64 << 20) as f64 {
                cur = cand;
                if since % 64 == 0 {
                    cur.greedy();
                }
                if (cur.n, cur.w.len()) < (best.n, best.w.len()) {
                    cur.greedy();
                    best = cur.clone();
                    since = 0;
                }
            }
        }
        best.greedy();
        BraidWord { strands: best.n, word: best.w, components: best.comps }
    }
}

/// xorshift64* (deterministic pseudo-random numbers for the braid search)
struct XorShift(u64);

impl XorShift {
    fn below(&mut self, n: u64) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        (self.0.wrapping_mul(0x2545_F491_4F6C_DD1D) >> 11) % n
    }
}

// ---------------------------------------------------------------------------
// braid simplification
// ---------------------------------------------------------------------------

#[derive(Clone)]
struct Work {
    n: usize,
    w: Vec<Gen>,
    comps: Vec<usize>,
}

fn inv(g: Gen) -> Gen {
    Gen { i: g.i, positive: !g.positive }
}

impl Work {
    fn greedy(&mut self) {
        loop {
            self.free_reduce_cyclic();
            if self.destab_both() {
                continue;
            }
            if self.try_destab_by_relations() {
                continue;
            }
            self.flip();
            let d = self.try_destab_by_relations();
            self.flip();
            if d {
                continue;
            }
            break;
        }
    }

    /// Destabilize at either end as long as possible.
    fn destab_both(&mut self) -> bool {
        let mut any = false;
        loop {
            if self.destab_last() {
                any = true;
                continue;
            }
            if self.n < 2 {
                return any;
            }
            self.flip();
            let d = self.destab_last();
            self.flip();
            if !d {
                return any;
            }
            any = true;
        }
    }

    /// Search objective: strands, then occurrences of the rarer end generator, then length.
    fn score(&self) -> f64 {
        if self.n < 2 {
            return self.w.len() as f64;
        }
        let first = self.w.iter().filter(|g| g.i == 0).count();
        let last = self.w.iter().filter(|g| g.i == self.n - 2).count();
        (self.n * 10_000 + first.min(last) * 100 + self.w.len()) as f64
    }

    /// Apply a length-three braid relation to the first three letters, if one applies:
    /// σ_i^a σ_j^a σ_i^a = σ_j^a σ_i^a σ_j^a and σ_i^a σ_j^b σ_i^{-a} = σ_j^{-a} σ_i^b σ_j^a
    /// (|i - j| = 1).
    fn relation_at_front(&mut self) -> bool {
        if self.w.len() < 3 {
            return false;
        }
        let (x, y, z) = (self.w[0], self.w[1], self.w[2]);
        if x.i != z.i || x.i.abs_diff(y.i) != 1 {
            return false;
        }
        let rep = if x.positive == y.positive && y.positive == z.positive {
            [y, x, y]
        } else if z.positive != x.positive {
            [Gen { i: y.i, positive: !x.positive }, Gen { i: x.i, positive: y.positive }, Gen { i: y.i, positive: x.positive }]
        } else {
            return false;
        };
        self.w[..3].copy_from_slice(&rep);
        true
    }

    /// Conjugate: move the first `r` letters to the end (bottom components follow).
    fn rotate(&mut self, r: usize) {
        if self.w.is_empty() {
            return;
        }
        let r = r % self.w.len();
        for g in &self.w[..r] {
            self.comps.swap(g.i, g.i + 1);
        }
        self.w.rotate_left(r);
    }

    /// Conjugate by the half twist: σ_i -> σ_{n-1-i} (1-based), components reversed.
    fn flip(&mut self) {
        for g in &mut self.w {
            g.i = self.n - 2 - g.i;
        }
        self.comps.reverse();
    }

    /// Cancel σ^ε ... σ^{-ε} when every letter in between commutes with σ.
    fn free_reduce_linear(&mut self) -> bool {
        let mut any = false;
        'outer: loop {
            for a in 0..self.w.len() {
                for b in a + 1..self.w.len() {
                    if self.w[b] == inv(self.w[a]) {
                        self.w.remove(b);
                        self.w.remove(a);
                        any = true;
                        continue 'outer;
                    }
                    if self.w[b].i.abs_diff(self.w[a].i) < 2 {
                        break;
                    }
                }
            }
            return any;
        }
    }

    fn free_reduce_cyclic(&mut self) -> bool {
        let mut any = false;
        loop {
            let mut changed = self.free_reduce_linear();
            if !changed {
                for _ in 0..self.w.len() {
                    self.rotate(1);
                    if self.free_reduce_linear() {
                        changed = true;
                        break;
                    }
                }
            }
            if !changed {
                return any;
            }
            any = true;
        }
    }

    /// Markov destabilization if σ_{n-1} (1-based) occurs exactly once.
    fn destab_last(&mut self) -> bool {
        if self.n < 2 {
            return false;
        }
        let l = self.n - 2;
        let occ: Vec<usize> = (0..self.w.len()).filter(|&k| self.w[k].i == l).collect();
        if occ.len() != 1 {
            return false;
        }
        self.rotate(occ[0] + 1);
        let g = self.w.pop().unwrap();
        debug_assert_eq!(g.i, l);
        let c = self.comps.pop().unwrap();
        debug_assert_eq!(c, self.comps[l]);
        self.n -= 1;
        true
    }

    /// One braid relation lowering the number of occurrences of the last generator σ_L:
    /// σ_L^a Y σ_{L-1}^b Z σ_L^c with Y, Z free of σ_{L-1}, σ_L and either c = -a
    /// (→ Y σ_{L-1}^{-a} σ_L^b σ_{L-1}^a Z) or a = b = c (→ Y σ_{L-1}^a σ_L^a σ_{L-1}^a Z).
    fn reduce_last_once(&mut self) -> bool {
        if self.n < 3 {
            return false;
        }
        let l = self.n - 2;
        let count = self.w.iter().filter(|g| g.i == l).count();
        if count < 2 {
            return false;
        }
        for _ in 0..count {
            let i0 = self.w.iter().position(|g| g.i == l).unwrap();
            self.rotate(i0);
            let j = 1 + self.w[1..].iter().position(|g| g.i == l).unwrap();
            let mids: Vec<usize> = (1..j).filter(|&m| self.w[m].i == l - 1).collect();
            if mids.len() == 1 {
                let m = mids[0];
                let (a, b, c) = (self.w[0].positive, self.w[m].positive, self.w[j].positive);
                let rep = if c != a {
                    Some([
                        Gen { i: l - 1, positive: !a },
                        Gen { i: l, positive: b },
                        Gen { i: l - 1, positive: a },
                    ])
                } else if a == b && b == c {
                    Some([Gen { i: l - 1, positive: a }, Gen { i: l, positive: a }, Gen { i: l - 1, positive: a }])
                } else {
                    None
                };
                if let Some(rep) = rep {
                    let mut nw: Vec<Gen> = self.w[1..m].to_vec();
                    nw.extend_from_slice(&rep);
                    nw.extend_from_slice(&self.w[m + 1..j]);
                    nw.extend_from_slice(&self.w[j + 1..]);
                    self.w = nw;
                    return true;
                }
            }
            self.rotate(j);
        }
        false
    }

    /// Apply relations to bring the last generator down to one occurrence, then
    /// destabilize; on failure the braid is left unchanged.
    fn try_destab_by_relations(&mut self) -> bool {
        let saved = self.clone();
        let limit = 4 * self.w.len() + 16;
        for _ in 0..limit {
            self.free_reduce_cyclic();
            if self.destab_last() {
                return true;
            }
            if !self.reduce_last_once() {
                break;
            }
        }
        *self = saved;
        false
    }
}

// ---------------------------------------------------------------------------
// oriented planar diagrams
// ---------------------------------------------------------------------------

/// Oriented planar diagram. Slots of a crossing are counterclockwise, slot 0 is the
/// incoming under-strand, slot 2 the outgoing under-strand; the over-strand enters at
/// slot 3 (positive crossing) or slot 1 (negative crossing).
#[derive(Clone, Debug)]
struct Diagram {
    x: Vec<[usize; 4]>,
    pos: Vec<bool>,
    /// edge -> (crossing, slot) where it ends / starts
    head: Vec<(usize, usize)>,
    tail: Vec<(usize, usize)>,
    /// edge -> link component
    comp: Vec<usize>,
    ncomp: usize,
}

impl Diagram {
    fn from_pd(pd: &Pd) -> Result<Diagram, String> {
        let mut labels: Vec<usize> = pd.crossings.iter().flatten().copied().collect();
        labels.sort_unstable();
        labels.dedup();
        let ne = labels.len();
        let x: Vec<[usize; 4]> =
            pd.crossings.iter().map(|c| c.map(|l| labels.binary_search(&l).unwrap())).collect();
        let mut occ: Vec<Vec<(usize, usize)>> = vec![vec![]; ne];
        for (c, xc) in x.iter().enumerate() {
            for s in 0..4 {
                occ[xc[s]].push((c, s));
            }
        }
        for (e, o) in occ.iter().enumerate() {
            if o.len() != 2 {
                return Err(format!("edge label {} occurs {} times (expected 2)", labels[e], o.len()));
            }
        }
        let other = |e: usize, o: (usize, usize)| if occ[e][0] == o { occ[e][1] } else { occ[e][0] };
        let mut head = vec![(usize::MAX, 0); ne];
        let mut tail = vec![(usize::MAX, 0); ne];
        let mut comp = vec![usize::MAX; ne];
        let mut ncomp = 0;
        for e0 in 0..ne {
            if comp[e0] != usize::MAX {
                continue;
            }
            // walk the component, tentatively entering e0's crossing at occ[e0][0]
            let mut steps: Vec<(usize, (usize, usize), (usize, usize))> = vec![];
            let (mut e, mut h) = (e0, occ[e0][0]);
            loop {
                steps.push((e, other(e, h), h));
                let p = (h.0, h.1 ^ 2);
                let e2 = x[p.0][p.1];
                h = other(e2, p);
                e = e2;
                if e == e0 && h == occ[e0][0] {
                    break;
                }
                if steps.len() > ne {
                    return Err("malformed PD code (component walk does not close)".into());
                }
            }
            let fwd = steps.iter().filter(|s| s.2 .1 == 0).count();
            let bwd = steps.iter().filter(|s| s.2 .1 == 2).count();
            if fwd > 0 && bwd > 0 {
                return Err(format!(
                    "inconsistent orientation of the component through label {}: \
                     under-strands point both ways",
                    labels[e0]
                ));
            }
            let reverse = if fwd + bwd > 0 {
                bwd > 0
            } else {
                // never passes under: use the label order
                let m = steps.len();
                let (mut inc, mut dec) = (0, 0);
                for k in 0..m {
                    let (a, b) = (labels[steps[k].0], labels[steps[(k + 1) % m].0]);
                    if b == a + 1 {
                        inc += 1;
                    }
                    if a == b + 1 {
                        dec += 1;
                    }
                }
                dec > inc
            };
            for &(e, t, h) in &steps {
                if comp[e] != usize::MAX {
                    return Err("malformed PD code (edge visited twice)".into());
                }
                comp[e] = ncomp;
                if reverse {
                    head[e] = t;
                    tail[e] = h;
                } else {
                    head[e] = h;
                    tail[e] = t;
                }
            }
            ncomp += 1;
        }
        let mut pos = vec![];
        for (c, xc) in x.iter().enumerate() {
            let p = head[xc[3]] == (c, 3);
            let ok = head[xc[0]] == (c, 0)
                && tail[xc[2]] == (c, 2)
                && if p { tail[xc[1]] == (c, 1) } else { head[xc[1]] == (c, 1) && tail[xc[3]] == (c, 3) };
            if !ok {
                return Err(format!("crossing {} has an inconsistent orientation", c + 1));
            }
            pos.push(p);
        }
        let d = Diagram { x, pos, head, tail, comp, ncomp };
        // planarity: every connected piece is a sphere (F = V + 2)
        let pieces = d.connected_pieces().len();
        let nf = d.faces().0.len();
        if nf != d.x.len() + 2 * pieces {
            return Err(format!(
                "PD code is not planar ({} faces, expected {})",
                nf,
                d.x.len() + 2 * pieces
            ));
        }
        Ok(d)
    }

    fn over_in(&self, c: usize) -> usize {
        if self.pos[c] { 3 } else { 1 }
    }

    /// Next edge along the Seifert circle (oriented smoothing).
    fn seifert_next(&self, e: usize) -> usize {
        let (c, s) = self.head[e];
        let out = if s == 0 {
            if self.pos[c] { 1 } else { 3 }
        } else {
            2
        };
        self.x[c][out]
    }

    /// Next edge along the link.
    fn link_next(&self, e: usize) -> usize {
        let (c, s) = self.head[e];
        self.x[c][s ^ 2]
    }

    /// (circle of each edge, edges of each circle in order)
    fn seifert_circles(&self) -> (Vec<usize>, Vec<Vec<usize>>) {
        let ne = self.head.len();
        let mut id = vec![usize::MAX; ne];
        let mut circles = vec![];
        for e0 in 0..ne {
            if id[e0] != usize::MAX {
                continue;
            }
            let k = circles.len();
            let mut cyc = vec![];
            let mut e = e0;
            loop {
                id[e] = k;
                cyc.push(e);
                e = self.seifert_next(e);
                if e == e0 {
                    break;
                }
            }
            circles.push(cyc);
        }
        (id, circles)
    }

    /// Faces as cycles of darts `(edge, forward)`, traversed with the face on the left;
    /// so `forward == true` means the face lies to the left of the oriented edge.
    /// Also returns `face_of[e][forward as usize]`.
    fn faces(&self) -> (Vec<Vec<(usize, bool)>>, Vec<[usize; 2]>) {
        let ne = self.head.len();
        let mut face_of = vec![[usize::MAX; 2]; ne];
        let mut faces = vec![];
        for e in 0..ne {
            for d in [false, true] {
                if face_of[e][d as usize] != usize::MAX {
                    continue;
                }
                let fid = faces.len();
                let mut f = vec![];
                let (mut ce, mut cd) = (e, d);
                loop {
                    face_of[ce][cd as usize] = fid;
                    f.push((ce, cd));
                    let (c, s) = if cd { self.head[ce] } else { self.tail[ce] };
                    // arriving through slot s, the face on the left continues out of
                    // the clockwise-next slot s - 1
                    let s2 = (s + 3) % 4;
                    let e2 = self.x[c][s2];
                    ce = e2;
                    cd = self.tail[e2] == (c, s2);
                    if face_of[ce][cd as usize] != usize::MAX {
                        debug_assert_eq!((ce, cd), (e, d));
                        break;
                    }
                }
                faces.push(f);
            }
        }
        (faces, face_of)
    }

    /// Connected pieces of the diagram (as separate diagrams, component ids kept).
    fn connected_pieces(&self) -> Vec<Diagram> {
        let nx = self.x.len();
        let mut parent: Vec<usize> = (0..nx).collect();
        fn find(p: &mut [usize], mut a: usize) -> usize {
            while p[a] != a {
                p[a] = p[p[a]];
                a = p[a];
            }
            a
        }
        for e in 0..self.head.len() {
            let (a, b) = (find(&mut parent, self.head[e].0), find(&mut parent, self.tail[e].0));
            parent[a] = b;
        }
        let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for c in 0..nx {
            let r = find(&mut parent, c);
            groups.entry(r).or_default().push(c);
        }
        if groups.len() <= 1 {
            return if nx == 0 { vec![] } else { vec![self.clone()] };
        }
        let mut out: Vec<(usize, Diagram)> = vec![];
        for (_, cs) in groups {
            let mut cmap = vec![usize::MAX; nx];
            for (k, &c) in cs.iter().enumerate() {
                cmap[c] = k;
            }
            let mut emap = vec![usize::MAX; self.head.len()];
            let mut edges = vec![];
            for e in 0..self.head.len() {
                if cmap[self.head[e].0] != usize::MAX {
                    emap[e] = edges.len();
                    edges.push(e);
                }
            }
            let d = Diagram {
                x: cs.iter().map(|&c| self.x[c].map(|e| emap[e])).collect(),
                pos: cs.iter().map(|&c| self.pos[c]).collect(),
                head: edges.iter().map(|&e| (cmap[self.head[e].0], self.head[e].1)).collect(),
                tail: edges.iter().map(|&e| (cmap[self.tail[e].0], self.tail[e].1)).collect(),
                comp: edges.iter().map(|&e| self.comp[e]).collect(),
                ncomp: self.ncomp,
            };
            out.push((cs[0], d));
        }
        out.sort_by_key(|p| p.0);
        out.into_iter().map(|p| p.1).collect()
    }

    /// A face with two edges on different Seifert circles that both have the face on
    /// the same side: (e1, e2, face is on the left of both).
    fn find_defect(seif: &[usize], faces: &[Vec<(usize, bool)>]) -> Option<(usize, usize, bool)> {
        let mut best: Option<(usize, (usize, usize, bool))> = None;
        for f in faces {
            for a in 0..f.len() {
                for b in a + 1..f.len() {
                    let ((e1, d1), (e2, d2)) = (f[a], f[b]);
                    if d1 == d2 && seif[e1] != seif[e2] {
                        // prefer small faces (fewer new edges near long boundaries)
                        if best.as_ref().is_none_or(|bb| f.len() < bb.0) {
                            best = Some((f.len(), (e1, e2, d1)));
                        }
                    }
                }
            }
        }
        best.map(|b| b.1)
    }

    /// Vogel move: Reidemeister II pushing a finger of e1 over e2 across their common face.
    fn vogel_move(&mut self, e1: usize, e2: usize, left: bool) {
        let (t1, h1) = (self.tail[e1], self.head[e1]);
        let (t2, h2) = (self.tail[e2], self.head[e2]);
        let ne = self.head.len();
        let (e1a, e1b, e1c) = (e1, ne, ne + 1);
        let (e2a, e2b, e2c) = (e2, ne + 2, ne + 3);
        self.head.extend([(0, 0); 4]);
        self.tail.extend([(0, 0); 4]);
        let (c1, c2) = (self.comp[e1], self.comp[e2]);
        self.comp.extend([c1, c1, c2, c2]);
        let p = self.x.len();
        let q = p + 1;
        if left {
            // face to the left of both: e1 runs east along the bottom, e2 west along the
            // top; e1 rises through e2 at P (x = 0) and comes back down at Q (x = 1)
            self.x.push([e2b, e1b, e2c, e1a]); // P: positive
            self.pos.push(true);
            self.x.push([e2a, e1b, e2b, e1c]); // Q: negative
            self.pos.push(false);
            self.tail[e1a] = t1;
            self.head[e1a] = (p, 3);
            self.tail[e1b] = (p, 1);
            self.head[e1b] = (q, 1);
            self.tail[e1c] = (q, 3);
            self.head[e1c] = h1;
            self.tail[e2a] = t2;
            self.head[e2a] = (q, 0);
            self.tail[e2b] = (q, 2);
            self.head[e2b] = (p, 0);
            self.tail[e2c] = (p, 2);
            self.head[e2c] = h2;
        } else {
            // face to the right of both: e1 runs west along the bottom, e2 east along
            // the top; e1 rises through e2 at p (x = 1), comes down at q (x = 0)
            self.x.push([e2b, e1a, e2c, e1b]); // negative
            self.pos.push(false);
            self.x.push([e2a, e1c, e2b, e1b]); // positive
            self.pos.push(true);
            self.tail[e1a] = t1;
            self.head[e1a] = (p, 1);
            self.tail[e1b] = (p, 3);
            self.head[e1b] = (q, 3);
            self.tail[e1c] = (q, 1);
            self.head[e1c] = h1;
            self.tail[e2a] = t2;
            self.head[e2a] = (q, 0);
            self.tail[e2b] = (q, 2);
            self.head[e2b] = (p, 0);
            self.tail[e2c] = (p, 2);
            self.head[e2c] = h2;
        }
        self.x[h1.0][h1.1] = e1c;
        self.x[h2.0][h2.1] = e2c;
    }

    /// Apply Vogel moves until the Seifert circles are coherently nested.
    fn vogel(&mut self) {
        let s0 = self.seifert_circles().1.len();
        let limit = 64 + 8 * s0 * s0 + 4 * self.x.len();
        for _ in 0..limit {
            let (seif, circles) = self.seifert_circles();
            let (faces, _) = self.faces();
            assert_eq!(circles.len(), s0, "Vogel move changed the number of Seifert circles");
            assert_eq!(faces.len(), self.x.len() + 2, "Vogel move broke planarity");
            match Self::find_defect(&seif, &faces) {
                None => return,
                Some((e1, e2, left)) => self.vogel_move(e1, e2, left),
            }
        }
        panic!("Vogel's algorithm did not terminate after {limit} moves");
    }

    /// Read the braid off a connected diagram whose Seifert circles are coherently
    /// nested.
    fn read_braid(&self) -> BraidWord {
        let (seif, circles) = self.seifert_circles();
        let s = circles.len();
        let nx = self.x.len();
        let mut adj: Vec<Vec<usize>> = vec![vec![]; s];
        let mut ends = vec![(0usize, 0usize); nx];
        for c in 0..nx {
            let a = seif[self.x[c][0]];
            let b = seif[self.x[c][self.over_in(c)]];
            assert_ne!(a, b, "crossing joins a Seifert circle to itself");
            ends[c] = (a, b);
            if !adj[a].contains(&b) {
                adj[a].push(b);
                adj[b].push(a);
            }
        }
        // the Seifert graph must be a path C_0 - C_1 - ... - C_{s-1}
        assert!(adj.iter().all(|a| a.len() <= 2), "Seifert circles are not nested");
        let start = (0..s).find(|&k| adj[k].len() <= 1).expect("Seifert graph is a cycle");
        let mut order = vec![start];
        let mut posn = vec![usize::MAX; s];
        posn[start] = 0;
        while order.len() < s {
            let cur = *order.last().unwrap();
            let nxt = adj[cur].iter().copied().find(|&k| posn[k] == usize::MAX).expect("Seifert graph disconnected");
            posn[nxt] = order.len();
            order.push(nxt);
        }
        // cut points: an edge on each circle, consecutive ones sharing a face
        let (faces, face_of) = self.faces();
        let mut cut = vec![circles[order[0]][0]];
        for p in 0..s - 1 {
            let e = cut[p];
            let target = order[p + 1];
            let found = (0..2).find_map(|side| {
                faces[face_of[e][side]].iter().find(|&&(e2, _)| seif[e2] == target).map(|&(e2, _)| e2)
            });
            cut.push(found.expect("no face joins consecutive Seifert circles"));
        }
        // order of crossings along each circle, starting after its cut point
        let mut succ: Vec<Vec<usize>> = vec![vec![]; nx];
        let mut indeg = vec![0usize; nx];
        for &e0 in &cut {
            let mut seq = vec![];
            let mut e = e0;
            loop {
                seq.push(self.head[e].0);
                e = self.seifert_next(e);
                if e == e0 {
                    break;
                }
            }
            for w in seq.windows(2) {
                succ[w[0]].push(w[1]);
                indeg[w[1]] += 1;
            }
        }
        let gen_of = |c: usize| {
            let (pa, pb) = (posn[ends[c].0], posn[ends[c].1]);
            assert_eq!(pa.abs_diff(pb), 1, "crossing between non-adjacent Seifert circles");
            Gen { i: pa.min(pb), positive: self.pos[c] }
        };
        let mut ready: BTreeSet<(usize, usize)> =
            (0..nx).filter(|&c| indeg[c] == 0).map(|c| (gen_of(c).i, c)).collect();
        let mut word = vec![];
        while let Some(&(i, c)) = ready.iter().next() {
            ready.remove(&(i, c));
            word.push(gen_of(c));
            for &d in &succ[c] {
                indeg[d] -= 1;
                if indeg[d] == 0 {
                    ready.insert((gen_of(d).i, d));
                }
            }
        }
        assert_eq!(word.len(), nx, "crossing order along Seifert circles is cyclic");
        let b = BraidWord { strands: s, word, components: cut.iter().map(|&e| self.comp[e]).collect() };
        // consistency: following the link from cut point p, the next cut point is perm[p]
        let perm = b.permutation();
        let mut cut_pos = vec![usize::MAX; self.head.len()];
        for (p, &e) in cut.iter().enumerate() {
            cut_pos[e] = p;
        }
        for p in 0..s {
            let mut e = self.link_next(cut[p]);
            while cut_pos[e] == usize::MAX {
                e = self.link_next(e);
            }
            assert_eq!(cut_pos[e], perm[p], "braid readout inconsistent with the diagram");
        }
        b
    }
}
