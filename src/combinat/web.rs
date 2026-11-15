//! Ladder webs: m upward uprights with labels in 0..=N, connected by rungs.
//!
//! Convention (skew Howe, Queffelec–Rose): E_i^{(r)} moves r units from upright
//! i+1 to upright i, F_i^{(r)} moves r units from upright i to upright i+1
//! (uprights 0-based, rung i joins uprights i and i+1). Rungs are listed bottom
//! to top.

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Dir {
    E,
    F,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Rung {
    pub i: usize,
    pub dir: Dir,
    pub r: usize,
}

impl Rung {
    pub fn e(i: usize, r: usize) -> Self {
        Rung { i, dir: Dir::E, r }
    }
    pub fn f(i: usize, r: usize) -> Self {
        Rung { i, dir: Dir::F, r }
    }
    /// (receiver, giver) uprights
    pub fn ends(&self) -> (usize, usize) {
        match self.dir {
            Dir::E => (self.i, self.i + 1),
            Dir::F => (self.i + 1, self.i),
        }
    }
    /// Apply to labels; None if some label leaves 0..=n.
    pub fn apply(&self, labels: &[usize], n: usize) -> Option<Vec<usize>> {
        let (rec, giv) = self.ends();
        if giv >= labels.len() || rec >= labels.len() {
            return None;
        }
        if labels[giv] < self.r || labels[rec] + self.r > n {
            return None;
        }
        let mut l = labels.to_vec();
        l[giv] -= self.r;
        l[rec] += self.r;
        Some(l)
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Ladder {
    pub n: usize,
    pub bottom: Vec<usize>,
    pub rungs: Vec<Rung>,
}

impl Ladder {
    pub fn new(n: usize, bottom: Vec<usize>) -> Self {
        Ladder { n, bottom, rungs: vec![] }
    }
    pub fn push(&mut self, r: Rung) -> &mut Self {
        self.rungs.push(r);
        self
    }
    /// Labels after each rung (len = rungs+1). Panics if invalid.
    pub fn levels(&self) -> Vec<Vec<usize>> {
        let mut out = vec![self.bottom.clone()];
        for r in &self.rungs {
            let next = r.apply(out.last().unwrap(), self.n).unwrap_or_else(|| panic!("invalid rung {:?} on {:?}", r, out.last()));
            out.push(next);
        }
        out
    }
    pub fn top(&self) -> Vec<usize> {
        self.levels().pop().unwrap()
    }
    pub fn is_valid(&self) -> bool {
        let mut l = self.bottom.clone();
        if l.iter().any(|&x| x > self.n) {
            return false;
        }
        for r in &self.rungs {
            match r.apply(&l, self.n) {
                Some(x) => l = x,
                None => return false,
            }
        }
        true
    }
    /// Vertical reflection: read from top to bottom (E <-> F).
    pub fn flip(&self) -> Ladder {
        let rungs = self
            .rungs
            .iter()
            .rev()
            .map(|r| Rung { i: r.i, dir: if r.dir == Dir::E { Dir::F } else { Dir::E }, r: r.r })
            .collect();
        Ladder { n: self.n, bottom: self.top(), rungs }
    }
    /// Stack `other` on top of self.
    pub fn compose(&self, other: &Ladder) -> Ladder {
        assert_eq!(self.top(), other.bottom, "compose: label mismatch");
        let mut l = self.clone();
        l.rungs.extend(other.rungs.iter().copied());
        l
    }
}

/// Random closed ladder (top labels == bottom labels) for testing.
pub fn random_closed_ladder<R: rand::Rng>(rng: &mut R, n: usize, bottom: Vec<usize>, steps: usize) -> Ladder {
    let m = bottom.len();
    let mut lad = Ladder::new(n, bottom.clone());
    let mut cur = bottom.clone();
    if m < 2 {
        return lad;
    }
    for _ in 0..steps {
        for _try in 0..20 {
            let i = rng.gen_range(0..m - 1);
            let dir = if rng.gen_bool(0.5) { Dir::E } else { Dir::F };
            let rmax = n;
            let r = rng.gen_range(1..=rmax);
            let rung = Rung { i, dir, r };
            if let Some(l) = rung.apply(&cur, n) {
                cur = l;
                lad.rungs.push(rung);
                break;
            }
        }
    }
    // restore labels greedily from left to right
    for i in 0..m - 1 {
        if cur[i] < bottom[i] {
            // need to move bottom[i]-cur[i] from the right: do it in steps
            let need = bottom[i] - cur[i];
            let mut moved = 0;
            while moved < need {
                // pull from i+1, pushing into i+1 from further right if empty
                if cur[i + 1] == 0 {
                    // find j>i+1 with cur[j]>0 and shift left one step at a time
                    let j = (i + 2..m).find(|&j| cur[j] > 0).expect("cannot restore");
                    for k in (i + 1..j).rev() {
                        let amt = cur[k + 1].min(n - cur[k]);
                        let rg = Rung::e(k, amt.max(1).min(cur[k + 1]));
                        cur = rg.apply(&cur, n).expect("restore");
                        lad.rungs.push(rg);
                    }
                    continue;
                }
                let amt = (need - moved).min(cur[i + 1]);
                let rg = Rung::e(i, amt);
                cur = rg.apply(&cur, n).expect("restore e");
                lad.rungs.push(rg);
                moved += amt;
            }
        } else if cur[i] > bottom[i] {
            let amt = cur[i] - bottom[i];
            // push to i+1; may exceed n, so move in chunks
            let mut left = amt;
            while left > 0 {
                let room = n - cur[i + 1];
                if room == 0 {
                    // push i+1 rightwards first
                    let j = (i + 2..m).find(|&j| cur[j] < n).expect("cannot restore (full)");
                    for k in (i + 1..j).rev() {
                        let amt2 = cur[k].min(n - cur[k + 1]);
                        let rg = Rung::f(k, amt2);
                        cur = rg.apply(&cur, n).expect("restore f2");
                        lad.rungs.push(rg);
                    }
                    continue;
                }
                let a = left.min(room);
                let rg = Rung::f(i, a);
                cur = rg.apply(&cur, n).expect("restore f");
                lad.rungs.push(rg);
                left -= a;
            }
        }
    }
    lad.rungs.retain(|r| r.r > 0);
    debug_assert!(lad.is_valid());
    lad
}
