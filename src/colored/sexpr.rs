//! Ladder expressions with zero labels (no truncation) and their flows, grouped by top coset.
//!
//! A flow of an expression is the vector of moved root masks, one per rung, in global root
//! numbering (`local::fixed_points`). `FlowIndex` stores the flows packed (binding layout, spec
//! F10): sorted by (top coset rank, lexicographic moved masks), all masks in one `Vec<u64>`.
use crate::combinat::coset::{subsets, Conf, CosetSpace, Mask, Roots};
use crate::combinat::web::{Dir, Rung};

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Expr {
    pub bottom: Vec<u8>,
    pub rungs: Vec<Rung>,
}

/// (receiver, giver) uprights of a rung
#[inline]
pub fn rung_ends(r: &Rung) -> (usize, usize) {
    match r.dir {
        Dir::E => (r.i, r.i + 1),
        Dir::F => (r.i + 1, r.i),
    }
}

impl Expr {
    /// giver >= r at every rung; r > 0; uprights in range; labels stay <= 255
    pub fn new(bottom: &[u8], rungs: &[Rung]) -> Result<Expr, String> {
        let k: usize = bottom.iter().map(|&a| a as usize).sum();
        if k > 64 {
            return Err(format!("Expr: K = {} > 64", k));
        }
        let mut lab = bottom.to_vec();
        for (t, r) in rungs.iter().enumerate() {
            if r.r == 0 {
                return Err(format!("Expr: rung {} ({:?}) has r = 0", t, r));
            }
            if r.i + 1 >= lab.len() {
                return Err(format!("Expr: rung {} ({:?}) out of range for {} uprights", t, r, lab.len()));
            }
            let (rec, giv) = rung_ends(r);
            if (lab[giv] as usize) < r.r {
                return Err(format!("Expr: rung {} ({:?}) moves {} but the giver has label {} (labels {:?})", t, r, r.r, lab[giv], lab));
            }
            if lab[rec] as usize + r.r > 255 {
                return Err(format!("Expr: rung {} ({:?}) overflows label {}", t, r, lab[rec]));
            }
            lab[giv] -= r.r as u8;
            lab[rec] += r.r as u8;
        }
        Ok(Expr { bottom: bottom.to_vec(), rungs: rungs.to_vec() })
    }
    /// labels after 0..=len rungs
    pub fn levels(&self) -> Vec<Vec<u8>> {
        let mut out = Vec::with_capacity(self.rungs.len() + 1);
        let mut lab = self.bottom.clone();
        out.push(lab.clone());
        for r in &self.rungs {
            let (rec, giv) = rung_ends(r);
            lab[giv] -= r.r as u8;
            lab[rec] += r.r as u8;
            out.push(lab.clone());
        }
        out
    }
    /// labels after `level` rungs
    pub fn level(&self, level: usize) -> Vec<u8> {
        let mut lab = self.bottom.clone();
        for r in &self.rungs[..level] {
            let (rec, giv) = rung_ends(r);
            lab[giv] -= r.r as u8;
            lab[rec] += r.r as u8;
        }
        lab
    }
    pub fn top(&self) -> Vec<u8> {
        self.level(self.rungs.len())
    }
    /// -sum_rungs r(b - r), b = giver label (= SBim::shift)
    pub fn shift(&self) -> i32 {
        let mut lab = self.bottom.clone();
        let mut s = 0i32;
        for r in &self.rungs {
            let (rec, giv) = rung_ends(r);
            let b = lab[giv] as i32;
            s -= r.r as i32 * (b - r.r as i32);
            lab[giv] -= r.r as u8;
            lab[rec] += r.r as u8;
        }
        s
    }
    /// = schur::flows(bottom, rungs)
    pub fn nflows(&self) -> u128 {
        crate::combinat::schur::flows(&self.bottom, &self.rungs)
    }
    /// exact bridge, K <= 16
    pub fn sbim<F: crate::algebra::field::Field>(&self) -> crate::soergel::sbim::SBim<F> {
        let lab: Vec<usize> = self.bottom.iter().map(|&a| a as usize).collect();
        let k: usize = lab.iter().sum();
        assert!(k <= crate::algebra::mpoly::MAXVARS, "Expr::sbim: K = {} > {} (mpoly::MAXVARS)", k, crate::algebra::mpoly::MAXVARS);
        crate::soergel::sbim::SBim::<F>::identity(2 * k.max(1), &lab).add_rungs(&self.rungs)
    }
    pub fn concat(&self, more: &[Rung]) -> Result<Expr, String> {
        let mut r = self.rungs.clone();
        r.extend_from_slice(more);
        Expr::new(&self.bottom, &r)
    }
    /// total label K
    pub fn k(&self) -> usize {
        self.bottom.iter().map(|&a| a as usize).sum()
    }
}

#[derive(Clone, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub struct Flow(pub Vec<Mask>); // owned copy of one flow (templates, tests, oracles)

/// Packed flow index (binding layout, F10): flows sorted by (top coset rank, lexicographic moved
/// masks); flow f (global id) occupies masks[f*nr .. (f+1)*nr]. No hash maps: heap memory is
/// 8·nr bytes per flow + 12 bytes per coset of the support.
#[derive(Clone, Debug)]
pub struct FlowIndex {
    pub expr: Expr,
    pub top: CosetSpace, // block (bottom, top labels)
    pub nr: usize,       // rungs = masks per flow
    pub cosets: Vec<u64>, // support, ascending rank
    pub start: Vec<u32>, // len = cosets.len() + 1; flows at cosets[i] are start[i]..start[i+1]
    pub masks: Vec<Mask>, // len = nr * nflows
}

impl FlowIndex {
    /// forward enumeration as local::fixed_points, then one sort
    pub fn build(expr: &Expr) -> FlowIndex {
        let roots = Roots::new(&expr.bottom);
        let top = CosetSpace::new(&expr.bottom, &expr.top());
        let nr = expr.rungs.len();
        let nf = expr.nflows();
        assert!(nf < u32::MAX as u128, "FlowIndex: {} flows exceed u32", nf);
        let nf = nf as usize;
        let mut flat: Vec<Mask> = Vec::with_capacity(nf * nr);
        let mut keys: Vec<(u64, u32)> = Vec::with_capacity(nf);
        let mut conf = roots.bottom().0;
        let mut cur: Vec<Mask> = vec![0; nr];
        // precomputed subset lists would cost memory; enumerate per node
        fn rec(expr: &Expr, top: &CosetSpace, level: usize, conf: &mut Vec<Mask>, cur: &mut Vec<Mask>, flat: &mut Vec<Mask>, keys: &mut Vec<(u64, u32)>) {
            if level == expr.rungs.len() {
                let id = keys.len() as u32;
                keys.push((top.rank(&Conf(conf.clone())), id));
                flat.extend_from_slice(cur);
                return;
            }
            let r = &expr.rungs[level];
            let (rec_u, giv) = rung_ends(r);
            for s in subsets(conf[giv], r.r) {
                conf[giv] &= !s;
                conf[rec_u] |= s;
                cur[level] = s;
                rec(expr, top, level + 1, conf, cur, flat, keys);
                conf[rec_u] &= !s;
                conf[giv] |= s;
            }
        }
        rec(expr, &top, 0, &mut conf, &mut cur, &mut flat, &mut keys);
        assert_eq!(keys.len(), nf, "FlowIndex: enumeration count != schur::flows");
        keys.sort_unstable_by(|a, b| a.0.cmp(&b.0).then_with(|| flat[a.1 as usize * nr..(a.1 as usize + 1) * nr].cmp(&flat[b.1 as usize * nr..(b.1 as usize + 1) * nr])));
        let mut masks: Vec<Mask> = Vec::with_capacity(nf * nr);
        let mut cosets: Vec<u64> = vec![];
        let mut start: Vec<u32> = vec![];
        for (g, &(c, id)) in keys.iter().enumerate() {
            if cosets.last() != Some(&c) {
                cosets.push(c);
                start.push(g as u32);
            }
            masks.extend_from_slice(&flat[id as usize * nr..(id as usize + 1) * nr]);
        }
        start.push(nf as u32);
        drop(flat);
        drop(keys);
        cosets.shrink_to_fit();
        start.shrink_to_fit();
        masks.shrink_to_fit();
        FlowIndex { expr: expr.clone(), top, nr, cosets, start, masks }
    }
    pub fn nflows(&self) -> usize {
        *self.start.last().unwrap() as usize
    }
    fn coset_pos(&self, c: u64) -> Option<usize> {
        self.cosets.binary_search(&c).ok()
    }
    /// global ids of the flows at c; empty outside the support
    pub fn range(&self, c: u64) -> std::ops::Range<u32> {
        match self.coset_pos(c) {
            Some(i) => self.start[i]..self.start[i + 1],
            None => 0..0,
        }
    }
    pub fn flow(&self, f: u32) -> &[Mask] {
        let f = f as usize;
        &self.masks[f * self.nr..(f + 1) * self.nr]
    }
    /// row index inside coset c (binary search)
    pub fn find(&self, c: u64, moved: &[Mask]) -> Option<u32> {
        let i = self.coset_pos(c)?;
        let (lo, hi) = (self.start[i] as usize, self.start[i + 1] as usize);
        let nr = self.nr;
        let (mut a, mut b) = (lo, hi);
        while a < b {
            let mid = (a + b) / 2;
            match self.masks[mid * nr..(mid + 1) * nr].cmp(moved) {
                std::cmp::Ordering::Less => a = mid + 1,
                std::cmp::Ordering::Greater => b = mid,
                std::cmp::Ordering::Equal => return Some((mid - lo) as u32),
            }
        }
        None
    }
    /// ascending
    pub fn support(&self) -> &[u64] {
        &self.cosets
    }
    /// configuration reached by replaying the moved masks `moved[..level]` from the bottom
    pub fn replay(&self, moved: &[Mask], level: usize) -> Conf {
        let mut c = Roots::new(&self.expr.bottom).bottom().0;
        for (r, &m) in self.expr.rungs[..level].iter().zip(moved) {
            let (rec, giv) = rung_ends(r);
            c[giv] &= !m;
            c[rec] |= m;
        }
        Conf(c)
    }
    /// configuration after `level` rungs
    pub fn conf_at(&self, f: u32, level: usize) -> Conf {
        self.replay(self.flow(f), level)
    }
    /// top configuration of a moved-mask vector (no validity check)
    pub fn top_conf(&self, moved: &[Mask]) -> Conf {
        self.replay(moved, self.nr)
    }
    /// global id of a flow given by its moved masks (None if it is not a flow)
    pub fn locate(&self, moved: &[Mask]) -> Option<u32> {
        if moved.len() != self.nr {
            return None;
        }
        let tc = self.top_conf(moved);
        if tc.0.iter().zip(&self.top.mu).any(|(m, &mu)| m.count_ones() as u8 != mu) {
            return None;
        }
        let c = self.top.rank(&tc);
        let i = self.coset_pos(c)?;
        self.find(c, moved).map(|r| self.start[i] + r)
    }
    /// coset (rank) of global flow f
    pub fn coset_of_flow(&self, f: u32) -> u64 {
        let i = self.start.partition_point(|&s| s <= f) - 1;
        self.cosets[i]
    }
    pub fn max_block(&self) -> usize {
        self.start.windows(2).map(|w| (w[1] - w[0]) as usize).max().unwrap_or(0)
    }
    /// heap bytes (A16)
    pub fn bytes(&self) -> usize {
        let expr = self.expr.bottom.capacity() + self.expr.rungs.capacity() * std::mem::size_of::<Rung>();
        let cs = self.top.roots.lam.capacity() + self.top.roots.off.capacity() * std::mem::size_of::<usize>() + self.top.mu.capacity();
        expr + cs + self.cosets.capacity() * 8 + self.start.capacity() * 4 + self.masks.capacity() * 8
    }
}

/// Flows per coset against the singular Deodhar coefficient c_{table(c)}(ladder_char)(1),
/// and every coset of every table with nonzero coefficient present.
pub fn deodhar_check(fi: &FlowIndex, s: &mut crate::combinat::schur::Schur) -> Result<(), String> {
    if s.lam != fi.expr.bottom {
        return Err(format!("deodhar_check: Schur lambda {:?} != bottom {:?}", s.lam, fi.expr.bottom));
    }
    let ch = s.ladder_char(&fi.expr.rungs);
    let mut expected: std::collections::BTreeMap<u64, i128> = std::collections::BTreeMap::new();
    let mut total: i128 = 0;
    for (tid, c) in &ch.t {
        let t = s.table(*tid).clone();
        let m = c.at1();
        if m < 0 {
            return Err(format!("deodhar_check: negative coefficient {} at table {}", m, t.fmt()));
        }
        if m == 0 {
            continue;
        }
        let cos = fi.top.cosets_of(&t);
        if cos.len() as u128 != t.ncosets() {
            return Err(format!("deodhar_check: cosets_of({}) has {} cosets, N(t) = {}", t.fmt(), cos.len(), t.ncosets()));
        }
        for c in cos {
            expected.insert(c, m);
            total += m;
        }
    }
    if total as usize != fi.nflows() {
        return Err(format!("deodhar_check {:?} {:?}: sum of coefficients {} != flows {}", fi.expr.bottom, fi.expr.rungs, total, fi.nflows()));
    }
    for (i, &c) in fi.cosets.iter().enumerate() {
        let n = (fi.start[i + 1] - fi.start[i]) as i128;
        match expected.get(&c) {
            Some(&m) if m == n => {}
            Some(&m) => {
                let t = fi.top.table(&fi.top.unrank(c));
                return Err(format!("deodhar_check {:?} {:?}: coset {} (table {}) has {} flows, coefficient {}", fi.expr.bottom, fi.expr.rungs, c, t.fmt(), n, m));
            }
            None => {
                let t = fi.top.table(&fi.top.unrank(c));
                return Err(format!("deodhar_check {:?} {:?}: coset {} (table {}) has {} flows but coefficient 0", fi.expr.bottom, fi.expr.rungs, c, t.fmt(), n));
            }
        }
    }
    if expected.len() != fi.cosets.len() {
        return Err(format!("deodhar_check {:?} {:?}: {} cosets expected, {} present", fi.expr.bottom, fi.expr.rungs, expected.len(), fi.cosets.len()));
    }
    Ok(())
}
