//! Localized Soergel calculus for uncolored braids (Soergel bimodules of S_n),
//! evaluated numerically at a generic point x ∈ k^n.
//!
//! A Bott–Samelson word (s_{i_1},...,s_{i_k}) localizes to k^{2^k}: fixed points are
//! bit strings (bit j = 1 iff letter j swaps). Elementary maps are exact bimodule
//! maps on small local webs (from `sbim::hom_space`), localized; tensoring with
//! identities is local application. Light leaves / double leaves are sequences of
//! local operations.

use crate::field::Field;
use crate::hecke::SymGroup;
use crate::local::{fixed_points, localize_map, FixedPoint};
use crate::mpoly::MonoCache;
use crate::pmat::PMat;
use crate::sbim::{hom_space, SBim};
use crate::web::Rung;
use std::collections::{HashMap, HashSet, VecDeque};

pub type Word = Vec<usize>;
/// Sparse localized vector: fixed point bits -> value.
pub type LVec<F> = HashMap<u64, F>;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Op {
    /// remove letter s at position
    DotEnd(usize),
    /// insert letter s at position
    DotStart(usize),
    /// letters (s,s) at pos,pos+1 -> s at pos
    Merge(usize),
    /// s at pos -> (s,s)
    Split(usize),
    /// (a,b) -> (b,a), |a-b| >= 2
    Comm(usize, usize),
    /// (a,b,a) -> (b,a,b), |a-b| = 1
    Braid(usize, usize),
}

impl Op {
    pub fn flip(&self) -> Op {
        match *self {
            Op::DotEnd(s) => Op::DotStart(s),
            Op::DotStart(s) => Op::DotEnd(s),
            Op::Merge(s) => Op::Split(s),
            Op::Split(s) => Op::Merge(s),
            Op::Comm(a, b) => Op::Comm(b, a),
            Op::Braid(a, b) => Op::Braid(b, a),
        }
    }
    /// (source letters, target letters)
    pub fn letters(&self) -> (Word, Word) {
        match *self {
            Op::DotEnd(s) => (vec![s], vec![]),
            Op::DotStart(s) => (vec![], vec![s]),
            Op::Merge(s) => (vec![s, s], vec![s]),
            Op::Split(s) => (vec![s], vec![s, s]),
            Op::Comm(a, b) => (vec![a, b], vec![b, a]),
            Op::Braid(a, b) => (vec![a, b, a], vec![b, a, b]),
        }
    }
    /// degree of the map
    pub fn degree(&self) -> i32 {
        match self {
            Op::DotEnd(_) | Op::DotStart(_) => 1,
            Op::Merge(_) | Op::Split(_) => -1,
            _ => 0,
        }
    }
}

/// Local template: elementary map on a small number of strands with local words.
struct Template<F: Field> {
    strands: usize,
    src: SBim<F>,
    tgt: SBim<F>,
    phi: PMat<F>,
    /// local letters (in local strand numbering)
    src_word: Word,
    tgt_word: Word,
}

fn bs_sbim<F: Field>(strands: usize, word: &[usize]) -> SBim<F> {
    let mut b = SBim::<F>::identity(1, &vec![1; strands]);
    // N is irrelevant for untruncated bimodules; use N = strands to keep labels valid
    b.n = strands.max(2);
    for &i in word {
        b = b.add_rung(&Rung::e(i, 1)).add_rung(&Rung::f(i, 1));
    }
    b
}

fn make_template<F: Field>(strands: usize, src_word: Word, tgt_word: Word, deg: i32) -> Template<F> {
    let src = bs_sbim::<F>(strands, &src_word);
    let tgt = bs_sbim::<F>(strands, &tgt_word);
    let mut mc = MonoCache::new(src.ring.clone());
    let u = deg + src.shift - tgt.shift;
    let h = hom_space(&src, &tgt, u, &mut mc);
    assert_eq!(h.len(), 1, "elementary map {:?}->{:?} not unique (dim {})", src_word, tgt_word, h.len());
    Template { strands, src, tgt, phi: h.into_iter().next().unwrap(), src_word, tgt_word }
}

/// Convert a local FixedPoint of a Bott–Samelson SBim to bits.
fn fp_bits(fp: &FixedPoint, word: &[usize], strands: usize) -> u64 {
    // replay: track root on each upright
    let mut pi: Vec<u64> = (0..strands).map(|u| 1u64 << u).collect();
    let mut bits = 0u64;
    for (j, &i) in word.iter().enumerate() {
        let a = pi[i];
        let moved = fp.moved[2 * j + 1];
        if moved == a {
            bits |= 1 << j;
            pi.swap(i, i + 1);
        }
    }
    bits
}

pub struct Ctx<F: Field> {
    pub n: usize,
    pub g: SymGroup,
    pub x: Vec<F>,
    templates: HashMap<u8, Template<F>>,
    cache: HashMap<(u8, Vec<u8>), Vec<(u64, u64, F)>>,
    /// memoized `rex_path` (point-independent breadth-first search over reduced words)
    rex_paths: std::cell::RefCell<HashMap<(Word, Word), Vec<(Op, usize)>>>,
}

fn kind_of(op: &Op) -> (u8, Vec<usize>) {
    // returns template id and local strand list (global uprights)
    match *op {
        Op::DotEnd(s) => (0, vec![s, s + 1]),
        Op::DotStart(s) => (1, vec![s, s + 1]),
        Op::Merge(s) => (2, vec![s, s + 1]),
        Op::Split(s) => (3, vec![s, s + 1]),
        Op::Comm(a, b) => {
            if a < b { (4, vec![a, a + 1, b, b + 1]) } else { (5, vec![b, b + 1, a, a + 1]) }
        }
        Op::Braid(a, b) => {
            if b == a + 1 { (6, vec![a, a + 1, a + 2]) } else { (7, vec![b, b + 1, b + 2]) }
        }
    }
}

impl<F: Field> Ctx<F> {
    pub fn new(n: usize, x: Vec<F>) -> Self {
        assert_eq!(x.len(), n);
        let mut t = HashMap::new();
        t.insert(0, make_template::<F>(2, vec![0], vec![], 1));
        t.insert(1, make_template::<F>(2, vec![], vec![0], 1));
        t.insert(2, make_template::<F>(2, vec![0, 0], vec![0], -1));
        t.insert(3, make_template::<F>(2, vec![0], vec![0, 0], -1));
        if n >= 4 {
            t.insert(4, make_template::<F>(4, vec![0, 2], vec![2, 0], 0));
            t.insert(5, make_template::<F>(4, vec![2, 0], vec![0, 2], 0));
        }
        if n >= 3 {
            t.insert(6, make_template::<F>(3, vec![0, 1, 0], vec![1, 0, 1], 0));
            t.insert(7, make_template::<F>(3, vec![1, 0, 1], vec![0, 1, 0], 0));
        }
        Ctx { n, g: SymGroup::new(n), x, templates: t, cache: HashMap::new(), rex_paths: Default::default() }
    }

    /// Localized local matrix entries (tgt bits, src bits, value) for template k at local roots.
    fn local(&mut self, k: u8, roots: &[u8]) -> Vec<(u64, u64, F)> {
        let key = (k, roots.to_vec());
        if let Some(v) = self.cache.get(&key) {
            return v.clone();
        }
        let t = &self.templates[&k];
        let xl: Vec<F> = roots.iter().map(|&r| self.x[r as usize].clone()).collect();
        let (p, pp, f) = localize_map(&t.src, &t.tgt, &t.phi, &xl);
        let sb: Vec<u64> = p.iter().map(|fp| fp_bits(fp, &t.src_word, t.strands)).collect();
        let tb: Vec<u64> = pp.iter().map(|fp| fp_bits(fp, &t.tgt_word, t.strands)).collect();
        let mut out = vec![];
        for (a, row) in f.iter().enumerate() {
            for (b, v) in row.iter().enumerate() {
                if !v.is_zero() {
                    out.push((tb[a], sb[b], v.clone()));
                }
            }
        }
        self.cache.insert(key, out.clone());
        out
    }

    /// Root sitting on each upright after the letters word[..len] with the given bits.
    fn perm_after(&self, word: &[usize], bits: u64, len: usize) -> Vec<u8> {
        let mut pi: Vec<u8> = (0..self.n as u8).collect();
        for j in 0..len {
            if bits >> j & 1 == 1 {
                pi.swap(word[j], word[j] + 1);
            }
        }
        pi
    }

    /// Apply an op at position pos to a localized vector on BS(word). Returns new word and vector.
    pub fn apply(&mut self, word: &[usize], v: &LVec<F>, op: Op, pos: usize) -> (Word, LVec<F>) {
        let (sl, tl) = op.letters();
        assert_eq!(&word[pos..pos + sl.len()], &sl[..], "op {:?} does not match word {:?} at {}", op, word, pos);
        let mut nw: Word = word[..pos].to_vec();
        nw.extend(tl.iter().copied());
        nw.extend(word[pos + sl.len()..].iter().copied());
        let (k, strands) = kind_of(&op);
        let ls = sl.len();
        let lt = tl.len();
        let mut out: LVec<F> = HashMap::new();
        for (&bits, val) in v {
            let pi = self.perm_after(word, bits, pos);
            let roots: Vec<u8> = strands.iter().map(|&u| pi[u]).collect();
            let loc_src = (bits >> pos) & ((1u64 << ls) - 1);
            let prefix = bits & ((1u64 << pos) - 1);
            let suffix = bits >> (pos + ls);
            let entries = self.local(k, &roots);
            for (tb, sb, c) in entries {
                if sb != loc_src {
                    continue;
                }
                let nb = prefix | (tb << pos) | (suffix << (pos + lt));
                out.entry(nb).or_insert_with(F::zero).fma(val, &c);
            }
        }
        out.retain(|_, c| !c.is_zero());
        (nw, out)
    }

    pub fn apply_seq(&mut self, word: &[usize], v: &LVec<F>, ops: &[(Op, usize)]) -> (Word, LVec<F>) {
        let mut w = word.to_vec();
        let mut cur = v.clone();
        for &(op, pos) in ops {
            let (nw, nv) = self.apply(&w, &cur, op, pos);
            w = nw;
            cur = nv;
        }
        (w, cur)
    }

    /// Element of S_n given by a word.
    pub fn elem(&self, word: &[usize]) -> usize {
        self.g.from_word(word)
    }

    /// Shortest sequence of braid moves from reduced word u to reduced word v (memoized).
    pub fn rex_path(&self, u: &[usize], v: &[usize]) -> Vec<(Op, usize)> {
        if u == v {
            return vec![];
        }
        let key = (u.to_vec(), v.to_vec());
        if let Some(p) = self.rex_paths.borrow().get(&key) {
            return p.clone();
        }
        let path = self.rex_path_bfs(u, v);
        self.rex_paths.borrow_mut().insert(key, path.clone());
        path
    }

    fn rex_path_bfs(&self, u: &[usize], v: &[usize]) -> Vec<(Op, usize)> {
        let mut prev: HashMap<Word, (Word, Op, usize)> = HashMap::new();
        let mut seen: HashSet<Word> = HashSet::new();
        let mut q = VecDeque::new();
        seen.insert(u.to_vec());
        q.push_back(u.to_vec());
        while let Some(w) = q.pop_front() {
            if w == v {
                break;
            }
            for j in 0..w.len() {
                if j + 1 < w.len() && (w[j] as i64 - w[j + 1] as i64).abs() >= 2 {
                    let mut w2 = w.clone();
                    w2.swap(j, j + 1);
                    if seen.insert(w2.clone()) {
                        prev.insert(w2.clone(), (w.clone(), Op::Comm(w[j], w[j + 1]), j));
                        q.push_back(w2);
                    }
                }
                if j + 2 < w.len() && w[j] == w[j + 2] && (w[j] as i64 - w[j + 1] as i64).abs() == 1 {
                    let w2: Word = w[..j].iter().copied().chain([w[j + 1], w[j], w[j + 1]]).chain(w[j + 3..].iter().copied()).collect();
                    if seen.insert(w2.clone()) {
                        prev.insert(w2.clone(), (w.clone(), Op::Braid(w[j], w[j + 1]), j));
                        q.push_back(w2);
                    }
                }
            }
        }
        let mut path = vec![];
        let mut cur = v.to_vec();
        while cur != u {
            let (p, op, j) = prev.get(&cur).expect("no rex path").clone();
            path.push((op, j));
            cur = p;
        }
        path.reverse();
        path
    }

    /// A reduced word for w ending with s (requires ws < w).
    pub fn rex_ending_with(&self, w: usize, s: usize) -> Word {
        let ws = self.g.rmul[w][s];
        let mut r = self.g.rex(ws);
        r.push(s);
        r
    }

    /// Operations of the light leaf LL_{word, e}: BS(word) -> BS(rex(w_e)), and its degree.
    pub fn light_leaf(&self, word: &[usize], e: u64) -> (Vec<(Op, usize)>, usize, i32) {
        let mut ops = vec![];
        let mut w = 0usize; // current element
        let mut deg = 0;
        for (i, &s) in word.iter().enumerate() {
            let p = self.g.len[w];
            let r = self.g.rex(w);
            let up = self.g.len[self.g.rmul[w][s]] > self.g.len[w];
            let bit = e >> i & 1 == 1;
            if up {
                if !bit {
                    ops.push((Op::DotEnd(s), p));
                    deg += 1;
                } else {
                    let mut r2 = r.clone();
                    r2.push(s);
                    let ws = self.g.rmul[w][s];
                    for (op, j) in self.rex_path(&r2, &self.g.rex(ws)) {
                        ops.push((op, j));
                    }
                    w = ws;
                }
            } else {
                let r2 = self.rex_ending_with(w, s);
                for (op, j) in self.rex_path(&r, &r2) {
                    ops.push((op, j));
                }
                ops.push((Op::Merge(s), p - 1));
                if !bit {
                    deg -= 1;
                    for (op, j) in self.rex_path(&r2, &self.g.rex(w)) {
                        ops.push((op, j));
                    }
                } else {
                    ops.push((Op::DotEnd(s), p - 1));
                    let ws = self.g.rmul[w][s];
                    let r3: Word = r2[..r2.len() - 1].to_vec();
                    for (op, j) in self.rex_path(&r3, &self.g.rex(ws)) {
                        ops.push((op, j));
                    }
                    w = ws;
                }
            }
        }
        (ops, w, deg)
    }

    /// End element w_e and degree of the light leaf LL_{word, e}, without building its
    /// operations (the same recursion as `light_leaf`).
    pub fn light_leaf_wd(&self, word: &[usize], e: u64) -> (usize, i32) {
        let mut w = 0usize;
        let mut deg = 0;
        for (i, &s) in word.iter().enumerate() {
            let up = self.g.len[self.g.rmul[w][s]] > self.g.len[w];
            let bit = e >> i & 1 == 1;
            if bit {
                w = self.g.rmul[w][s];
            } else if up {
                deg += 1;
            } else {
                deg -= 1;
            }
        }
        (w, deg)
    }

    /// Double leaf LL_{y,f}^flip ∘ LL_{x,e} as a list of operations on BS(x).
    pub fn double_leaf_ops(&self, x: &[usize], e: u64, y: &[usize], f: u64) -> Option<(Vec<(Op, usize)>, i32)> {
        let (o1, w1, d1) = self.light_leaf(x, e);
        let (o2, w2, d2) = self.light_leaf(y, f);
        if w1 != w2 {
            return None;
        }
        let mut ops = o1;
        for &(op, j) in o2.iter().rev() {
            ops.push((op.flip(), j));
        }
        Some((ops, d1 + d2))
    }

    /// Localized matrix (as columns per source bits) of an op sequence on BS(word).
    pub fn map_columns(&mut self, word: &[usize], ops: &[(Op, usize)]) -> (Word, Vec<LVec<F>>) {
        let k = word.len();
        let mut cols = vec![];
        let mut tw = vec![];
        for b in 0..(1u64 << k) {
            let mut v = HashMap::new();
            v.insert(b, F::one());
            let (w, out) = self.apply_seq(word, &v, ops);
            tw = w;
            cols.push(out);
        }
        (tw, cols)
    }
}

/// All double leaves Hom(BS(x), BS(y)) of a given degree (op sequences).
pub fn double_leaves<F: Field>(ctx: &Ctx<F>, x: &[usize], y: &[usize], deg: i32) -> Vec<(u64, u64, Vec<(Op, usize)>)> {
    double_leaves_filtered(ctx, x, y, Some(deg)).into_iter().map(|(_, e, f, ops)| (e, f, ops)).collect()
}

/// All double leaves Hom(BS(x), BS(y)) of every degree, as (degree, e, f, ops), sorted by
/// (degree, e, f): the concatenation of `double_leaves` over increasing degrees.
pub fn double_leaves_all<F: Field>(ctx: &Ctx<F>, x: &[usize], y: &[usize]) -> Vec<(i32, u64, u64, Vec<(Op, usize)>)> {
    double_leaves_filtered(ctx, x, y, None)
}

/// Matching uses only the end element and degree of each light leaf (cheap); the
/// operations are built once per leaf and only for leaves that occur in a double leaf.
fn double_leaves_filtered<F: Field>(ctx: &Ctx<F>, x: &[usize], y: &[usize], deg: Option<i32>) -> Vec<(i32, u64, u64, Vec<(Op, usize)>)> {
    let lx: Vec<(usize, i32)> = (0..(1u64 << x.len())).map(|e| ctx.light_leaf_wd(x, e)).collect();
    let mut by_w: HashMap<usize, Vec<(u64, i32)>> = HashMap::new();
    for f in 0..(1u64 << y.len()) {
        let (w, d) = ctx.light_leaf_wd(y, f);
        by_w.entry(w).or_default().push((f, d));
    }
    let mut ops_y: HashMap<u64, Vec<(Op, usize)>> = HashMap::new();
    let mut out = vec![];
    for (e, &(w1, d1)) in lx.iter().enumerate() {
        let Some(fs) = by_w.get(&w1) else { continue };
        let mut ops_x: Option<Vec<(Op, usize)>> = None;
        for &(f, d2) in fs {
            if deg.is_some_and(|d| d != d1 + d2) {
                continue;
            }
            let ox = ops_x.get_or_insert_with(|| ctx.light_leaf(x, e as u64).0);
            let oy = ops_y.entry(f).or_insert_with(|| ctx.light_leaf(y, f).0);
            let mut ops = ox.clone();
            for &(op, j) in oy.iter().rev() {
                ops.push((op.flip(), j));
            }
            out.push((d1 + d2, e as u64, f, ops));
        }
    }
    out.sort_by_key(|t| (t.0, t.1, t.2));
    out
}
