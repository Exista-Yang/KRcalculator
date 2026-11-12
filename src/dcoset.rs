//! src/dcoset.rs — double-coset combinatorics for the M4c light-leaves backend (milestone A1).
//!
//! Pure combinatorics (no bimodule maps): it runs in seconds and is the input of the G-L0 census
//! (src/bin/m4c.rs) and of the later `leaves` provider. Conventions as in schur.rs: a (λ, μ)
//! double coset is a contingency `Table` (rows = bottom blocks λ, columns = top blocks μ, zero
//! labels allowed); a rung E_i^(r) / F_i^(r) is split-then-merge of columns (schur.rs).
//!
//! Paper keys (docs/colored_design/m4c/L1_eklp.md): SLL = Elias-Ko-Libedinsky-Patimo, singular
//! light leaves, arXiv 2401.03053; CORE = arXiv 2402.08673 (reduced expressions of core double
//! cosets); SUB = Bruhat order on double cosets (cited in SLL, Def 2.4).
//!
//! Contents
//! * [`conf_type`], [`rep_rank`]: configuration -> table, representative configuration of a table.
//! * [`step_gf_table`] / [`step_gf`]: the generating function G_{t'->t}(v) of one rung, written
//!   directly from SLL Def 3.5 (descending defect ℓ(q̄) − ℓ(p̄); ascending defect
//!   ℓ(q_) − ℓ(p_) − ℓ(LR q) + ℓ(LR p) and polynomial factor π⁺(LR q)/π⁺(LR p) as a product of
//!   Gaussian binomials). It is an implementation independent of `Schur::rung`; gate G-L1 checks
//!   that the Deodhar sums Σ_{t'} c_{t'}(X')·G_{t'->t} equal `Schur::ladder_char`.
//! * [`forward`], [`reach_back`], [`junction_types`]: types per level of a word, the types from
//!   which a target set is reachable through the remaining rungs, their intersection, and the
//!   degree windows of paths through each (level, type).
//! * [`min_perm`], [`bruhat_le`], [`reduced_for`]: minimal coset representatives in S_K, the
//!   Bruhat order of double cosets (p ≤ q iff p_ ≤ q_, SUB Thm 2.16), and the reducedness test of
//!   a junction word (c_t(char W) = 1 exactly and every table of char W is ≤ t).
//! * [`junction_word`]: prefix-extension choice of reduced junction words (C_m4c m3), with the
//!   minimal-flow canonical ladder and a bounded BFS as fallbacks.
//! * [`psi`], [`grass`]: CORE's squashed permutation ψ(p) and SLL Thm 2.6 Grassmannian reduction
//!   at rung level (z, m', m'' with t' = z.m', t = z.m'', same z).
use crate::coset::{Conf, CosetSpace};
use crate::schur::{Lp, Schur, Table, Tid};
use crate::web::{Dir, Rung};
use std::collections::{BTreeMap, HashMap, VecDeque};

// ---------------------------------------------------------------------------------------------
// configurations

/// The table (double coset) of a top configuration: M_ij = |B_i ∩ T_j|.
pub fn conf_type(cs: &CosetSpace, conf: &Conf) -> Table {
    cs.table(conf)
}

/// Rank of the representative configuration of a table (`CosetSpace::rep`), as `rep_rank` of
/// homprov (private there).
pub fn rep_rank(cs: &CosetSpace, t: &Table) -> u64 {
    cs.rank(&cs.rep(t))
}

// ---------------------------------------------------------------------------------------------
// one rung on tables (SLL Def 3.5, written independently of schur.rs)

fn c2(n: i64) -> i64 {
    n * (n - 1) / 2
}
fn l_lr(t: &Table) -> i64 {
    t.e.iter().map(|&a| c2(a as i64)).sum()
}

/// Gaussian binomial [n choose k] as a polynomial in v (exponents 0, 2, 4, ...): π⁺(n)/(π⁺(k)π⁺(n−k))
fn gauss_v(n: usize, k: usize) -> Lp {
    let g = crate::homprov::gauss_binom(n, k);
    let mut out = Lp::zero();
    for (j, &c) in g.iter().enumerate() {
        if c != 0 {
            out = out.add(&Lp::mono(2 * j as i32, c as i128));
        }
    }
    out
}

/// (receiver, giver) uprights of a rung
pub fn rung_rg(r: &Rung) -> (usize, usize) {
    match r.dir {
        Dir::E => (r.i, r.i + 1),
        Dir::F => (r.i + 1, r.i),
    }
}

/// labels after a rung, None if invalid
pub fn apply_labels(l: &[u8], r: &Rung) -> Option<Vec<u8>> {
    let (rec, giv) = rung_rg(r);
    if giv >= l.len() || rec >= l.len() || (l[giv] as usize) < r.r {
        return None;
    }
    let mut v = l.to_vec();
    v[giv] -= r.r as u8;
    v[rec] += r.r as u8;
    Some(v)
}

/// One term of a rung step: the moved vector x (x_i roots of row i go from the giver column to
/// the receiver column), the target table, the refined (split) table, and v^def·poly.
#[derive(Clone, Debug)]
pub struct StepTerm {
    pub x: Vec<u8>,
    pub to: Table,
    pub split: Table,
    pub def_desc: i64,
    pub def_asc: i64,
    pub gf: Lp,
}

/// All steps of the rung r from the table t (columns = uprights). The split refines the giver
/// column g into (moving part, staying part) in the order of the rung (E: moving part first, as
/// schur.rs "split col i+1 as (r, b−r)"; F: staying part first), then the moving part merges
/// into the receiver column. Defects per SLL Def 3.5: descending ℓ(t̄) − ℓ(ū); ascending
/// ℓ(t_) − ℓ(u_) − ℓ(LR t) + ℓ(LR u), polynomial factor π⁺(LR t)/π⁺(LR u) = Π_i [t_{i,rec}, x_i]_{v²}.
pub fn step_terms(t: &Table, r: &Rung) -> Vec<StepTerm> {
    let (n, m) = (t.nr as usize, t.nc as usize);
    let (rec, giv) = rung_rg(r);
    assert!(giv < m && rec < m, "step_terms: rung {:?} outside {} columns", r, m);
    let col: Vec<u8> = (0..n).map(|i| t.get(i, giv)).collect();
    let mut out = vec![];
    // enumerate x with x_i <= col_i, Σ x = r
    let mut x = vec![0u8; n];
    fn rec_x(i: usize, left: usize, col: &[u8], x: &mut Vec<u8>, acc: &mut Vec<Vec<u8>>) {
        if i == col.len() {
            if left == 0 {
                acc.push(x.clone());
            }
            return;
        }
        let hi = (col[i] as usize).min(left);
        for v in 0..=hi {
            x[i] = v as u8;
            rec_x(i + 1, left - v, col, x, acc);
        }
        x[i] = 0;
    }
    let mut xs = vec![];
    rec_x(0, r.r, &col, &mut x, &mut xs);
    for x in xs {
        // split table u: columns 0..m with the giver column replaced by two columns
        let (first, second): (Vec<u8>, Vec<u8>) = match r.dir {
            Dir::E => (x.clone(), (0..n).map(|i| col[i] - x[i]).collect()),
            Dir::F => ((0..n).map(|i| col[i] - x[i]).collect(), x.clone()),
        };
        let mut ue = Vec::with_capacity(n * (m + 1));
        for i in 0..n {
            for j in 0..m {
                if j == giv {
                    ue.push(first[i]);
                    ue.push(second[i]);
                } else {
                    ue.push(t.get(i, j));
                }
            }
        }
        let u = Table { nr: n as u8, nc: (m + 1) as u8, e: ue };
        // target table
        let mut te = t.e.clone();
        for i in 0..n {
            te[i * m + giv] -= x[i];
            te[i * m + rec] += x[i];
        }
        let to = Table { nr: n as u8, nc: m as u8, e: te };
        // merged (coarse) table q after the merge = `to`; the merged pair of u is (x, t[.,rec])
        let def_desc = t.lmax() - u.lmax();
        let def_asc = to.lmin() - u.lmin() - l_lr(&to) + l_lr(&u);
        let mut poly = Lp::one();
        for i in 0..n {
            let a = to.get(i, rec) as usize;
            if x[i] > 0 && a > x[i] as usize {
                poly = poly.mul(&gauss_v(a, x[i] as usize));
            }
        }
        let gf = poly.shift((def_desc + def_asc) as i32);
        out.push(StepTerm { x, to, split: u, def_desc, def_asc, gf });
    }
    out
}

/// G_{t->·} of one rung on tables: (target table, G(v)), one entry per target (distinct moved
/// vectors give distinct targets).
pub fn step_gf_table(t: &Table, r: &Rung) -> Vec<(Table, Lp)> {
    step_terms(t, r).into_iter().map(|s| (s.to, s.gf)).collect()
}

/// Cache of rung steps on interned tables.
#[derive(Default)]
pub struct StepCache {
    map: HashMap<(Tid, (usize, bool, usize)), Vec<(Tid, Lp)>>,
}
impl StepCache {
    pub fn new() -> StepCache {
        StepCache::default()
    }
    /// G_{t->·}(v) for the rung r, as interned tables
    pub fn step(&mut self, s: &mut Schur, t: Tid, r: &Rung) -> Vec<(Tid, Lp)> {
        let key = (t, (r.i, r.dir == Dir::E, r.r));
        if let Some(v) = self.map.get(&key) {
            return v.clone();
        }
        let v = step_gf(s, t, r);
        self.map.insert(key, v.clone());
        v
    }
}

/// G_{t->·}(v) for the rung r on an interned table (sorted by target Tid).
pub fn step_gf(s: &mut Schur, t: Tid, r: &Rung) -> Vec<(Tid, Lp)> {
    if r.r == 0 {
        return vec![(t, Lp::one())];
    }
    let tab = s.table(t).clone();
    let mut v: Vec<(Tid, Lp)> = step_gf_table(&tab, r).into_iter().map(|(to, g)| (s.intern(to), g)).collect();
    v.sort_by_key(|e| e.0);
    v
}

/// Forward Deodhar sums of a word: entry k = Σ_paths v^def h_{t_k} after k rungs (entry 0 = the
/// identity table with coefficient 1). By SLL Cor 3.7 the last entry is the character
/// `Schur::ladder_char` (gate G-L1).
pub fn forward(s: &mut Schur, sc: &mut StepCache, w: &[Rung]) -> Vec<BTreeMap<Tid, Lp>> {
    let lam = s.lam.clone();
    let id = s.intern(Table::identity(&lam));
    let mut cur: BTreeMap<Tid, Lp> = BTreeMap::new();
    cur.insert(id, Lp::one());
    let mut out = vec![cur.clone()];
    for r in w {
        let mut next: BTreeMap<Tid, Lp> = BTreeMap::new();
        for (t, c) in &cur {
            for (u, g) in sc.step(s, *t, r) {
                let e = next.entry(u).or_insert_with(Lp::zero);
                *e = e.add(&c.mul(&g));
            }
        }
        next.retain(|_, c| !c.is_zero());
        out.push(next.clone());
        cur = next;
    }
    out
}

/// Backward reachability: entry k = the types at level k (every table of the level's block that
/// the forward pass meets) from which some target is reached through rungs k.. of w, with the
/// generating function of the path suffixes (Σ over suffix paths of v^def, summed over targets).
pub fn reach_back(s: &mut Schur, sc: &mut StepCache, w: &[Rung], targets: &[Tid]) -> Vec<BTreeMap<Tid, Lp>> {
    let fw = forward(s, sc, w);
    let m = w.len();
    let mut out: Vec<BTreeMap<Tid, Lp>> = vec![BTreeMap::new(); m + 1];
    for &t in targets {
        if fw[m].contains_key(&t) {
            out[m].insert(t, Lp::one());
        }
    }
    for k in (0..m).rev() {
        let mut here: BTreeMap<Tid, Lp> = BTreeMap::new();
        let keys: Vec<Tid> = fw[k].keys().copied().collect();
        for t in keys {
            let mut acc = Lp::zero();
            for (u, g) in sc.step(s, t, &w[k]) {
                if let Some(b) = out[k + 1].get(&u) {
                    acc = acc.add(&g.mul(b));
                }
            }
            if !acc.is_zero() {
                here.insert(t, acc);
            }
        }
        out[k] = here;
    }
    out
}

/// One (level, type) on a path to the targets: forward and backward generating functions and the
/// degree window [dmin, dmax] of whole paths through it.
#[derive(Clone, Debug)]
pub struct Junction {
    pub level: usize,
    pub t: Tid,
    pub fwd: Lp,
    pub back: Lp,
    pub dmin: i32,
    pub dmax: i32,
}

/// The junction types of a word for a target set: per level, the types with nonzero forward and
/// backward generating functions. Σ_t fwd·back at any level = Σ_targets c_target(char w).
pub fn junction_types(s: &mut Schur, sc: &mut StepCache, w: &[Rung], targets: &[Tid]) -> Vec<Vec<Junction>> {
    let fw = forward(s, sc, w);
    let bk = reach_back(s, sc, w, targets);
    let mut out = vec![];
    for k in 0..=w.len() {
        let mut lv = vec![];
        for (t, b) in &bk[k] {
            if let Some(f) = fw[k].get(t) {
                let p = f.mul(b);
                if p.is_zero() {
                    continue;
                }
                lv.push(Junction { level: k, t: *t, fwd: f.clone(), back: b.clone(), dmin: p.min_exp().unwrap(), dmax: p.max_exp().unwrap() });
            }
        }
        out.push(lv);
    }
    out
}

// ---------------------------------------------------------------------------------------------
// permutations, minimal representatives, Bruhat order

/// Minimal element p_ of the double coset of a table, as a permutation of 0..K in one-line
/// notation: position a (top positions: top block 0 first, ascending) ↦ root (bottom block i owns
/// roots off_i..off_i+λ_i). Top block j takes, for i = 0.., the next M_ij unused roots of B_i
/// (the witness of `CosetSpace::rep`). Its length is `Table::lmin` (tested).
pub fn min_perm(t: &Table) -> Vec<u8> {
    let (n, m) = (t.nr as usize, t.nc as usize);
    let rows = t.rows();
    let mut next = vec![0usize; n];
    let mut off = 0usize;
    for i in 0..n {
        next[i] = off;
        off += rows[i] as usize;
    }
    let mut out = Vec::with_capacity(off);
    for j in 0..m {
        for i in 0..n {
            for _ in 0..t.get(i, j) {
                out.push(next[i] as u8);
                next[i] += 1;
            }
        }
    }
    out
}

/// number of inversions
pub fn perm_len(w: &[u8]) -> i64 {
    let mut s = 0i64;
    for a in 0..w.len() {
        for b in a + 1..w.len() {
            if w[a] > w[b] {
                s += 1;
            }
        }
    }
    s
}

/// (u ∘ v)(a) = u(v(a))
pub fn perm_comp(u: &[u8], v: &[u8]) -> Vec<u8> {
    v.iter().map(|&a| u[a as usize]).collect()
}

pub fn perm_inv(u: &[u8]) -> Vec<u8> {
    let mut o = vec![0u8; u.len()];
    for (a, &b) in u.iter().enumerate() {
        o[b as usize] = a as u8;
    }
    o
}

/// Bruhat order on S_K (tableau criterion): u ≤ w iff for all i, j:
/// #{a ≤ i : u(a) ≥ j} ≤ #{a ≤ i : w(a) ≥ j}.
pub fn perm_bruhat_le(u: &[u8], w: &[u8]) -> bool {
    let k = u.len();
    assert_eq!(k, w.len());
    for j in 0..k as u8 {
        let (mut cu, mut cw) = (0i32, 0i32);
        for i in 0..k {
            cu += (u[i] >= j) as i32;
            cw += (w[i] >= j) as i32;
            if cu > cw {
                return false;
            }
        }
    }
    true
}

/// Bruhat order of double cosets of the same (λ, μ): p ≤ q iff p_ ≤ q_ (SUB Thm 2.16).
pub fn bruhat_le(p: &Table, q: &Table) -> bool {
    assert_eq!((p.rows(), p.cols()), (q.rows(), q.cols()), "bruhat_le: different blocks");
    perm_bruhat_le(&min_perm(p), &min_perm(q))
}

/// Table of the double coset (I rows, J columns) containing the permutation w (positions of J ↦
/// roots of I): M_ij = #{a ∈ J_j : w(a) ∈ I_i}.
pub fn table_of_perm(w: &[u8], rows: &[u8], cols: &[u8]) -> Table {
    let (n, m) = (rows.len(), cols.len());
    let mut blk = vec![0usize; w.len()];
    let mut off = 0usize;
    for (i, &a) in rows.iter().enumerate() {
        for b in off..off + a as usize {
            blk[b] = i;
        }
        off += a as usize;
    }
    let mut e = vec![0u8; n * m];
    let mut pos = 0usize;
    for (j, &b) in cols.iter().enumerate() {
        for _ in 0..b {
            e[blk[w[pos] as usize] * m + j] += 1;
            pos += 1;
        }
    }
    Table { nr: n as u8, nc: m as u8, e }
}

// ---------------------------------------------------------------------------------------------
// reducedness of junction words

/// Reducedness of a word W for the type t (D1 §1.1): `coef_one` = c_t(char W) is exactly 1
/// (so W has exactly one flow at every configuration of type t), `support_le` = every table of
/// char W is ≤ t in the Bruhat order (W has no flows at configurations of type ≰ t).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Reduced {
    pub coef_one: bool,
    pub support_le: bool,
}
impl Reduced {
    pub fn ok(&self) -> bool {
        self.coef_one && self.support_le
    }
}

/// Reducedness test of a word for t (uses `Schur::ladder_char`).
pub fn reduced_for(s: &mut Schur, w: &[Rung], t: Tid) -> Reduced {
    let lam = s.lam.clone();
    if !valid_word(&lam, w) {
        return Reduced { coef_one: false, support_le: false };
    }
    let ch = s.ladder_char(w);
    let tt = s.table(t).clone();
    if ch.mu != tt.cols() {
        return Reduced { coef_one: false, support_le: false };
    }
    let coef_one = ch.coef(t) == Lp::one();
    let tp = min_perm(&tt);
    let support_le = ch.t.iter().all(|(u, _)| perm_bruhat_le(&min_perm(s.table(*u)), &tp));
    Reduced { coef_one, support_le }
}

pub fn valid_word(lam: &[u8], w: &[Rung]) -> bool {
    let mut l = lam.to_vec();
    for r in w {
        if r.r == 0 {
            return false;
        }
        match apply_labels(&l, r) {
            Some(x) => l = x,
            None => return false,
        }
    }
    true
}

pub fn top_labels(lam: &[u8], w: &[Rung]) -> Vec<u8> {
    let mut l = lam.to_vec();
    for r in w {
        l = apply_labels(&l, r).expect("top_labels: invalid word");
    }
    l
}

/// How a junction word was chosen.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum WordKind {
    /// W_t = W_{t'}·ρ for the junction t' at the previous level (C_m4c m3): the top step is Id
    Prefix(Tid),
    /// the minimal-flow canonical ladder of t (reduced)
    Canonical,
    /// a word found by the bounded BFS (the canonical ladder is not reduced in the strict sense)
    Bfs,
    /// no reduced word found (canonical ladder returned, flagged for the census: C_m4c m2)
    NotReduced,
}

/// Reduced junction word for t at a level: prefix extensions W_{t'}·ρ (t' over `prev`, the
/// junction words of the previous level) that are reduced for t, minimal flows first; otherwise
/// the minimal-flow canonical ladder if it is reduced; otherwise a bounded BFS (words of at most
/// `bfs_len` rungs, at most `bfs_cap` words) for a reduced word of minimal flows.
pub fn junction_word(s: &mut Schur, t: Tid, prev: &[(Tid, Vec<Rung>)], rho: Option<&Rung>, bfs_len: usize, bfs_cap: usize) -> (Vec<Rung>, WordKind) {
    let lam = s.lam.clone();
    if let Some(r) = rho {
        let mut best: Option<(u128, Vec<Rung>, Tid)> = None;
        for (tp, wp) in prev {
            let mut w2 = wp.clone();
            w2.push(*r);
            if !valid_word(&lam, &w2) {
                continue;
            }
            if reduced_for(s, &w2, t).ok() {
                let f = crate::schur::flows(&lam, &w2);
                if best.as_ref().map_or(true, |b| f < b.0) {
                    best = Some((f, w2, *tp));
                }
            }
        }
        if let Some((_, w, tp)) = best {
            return (w, WordKind::Prefix(tp));
        }
    }
    let can = crate::sindec::canonical_ladder(s, t);
    if reduced_for(s, &can, t).ok() {
        return (can, WordKind::Canonical);
    }
    match bfs_reduced(s, t, bfs_len, bfs_cap) {
        Some(w) => (w, WordKind::Bfs),
        None => (can, WordKind::NotReduced),
    }
}

/// Breadth-first search over valid words (length ≤ max_len, ≤ cap words visited) for a word
/// reduced for t; minimal flows among the reduced words of the shortest length found.
pub fn bfs_reduced(s: &mut Schur, t: Tid, max_len: usize, cap: usize) -> Option<Vec<Rung>> {
    let lam = s.lam.clone();
    let mu = s.table(t).cols();
    let n = lam.len();
    let mut q: VecDeque<Vec<Rung>> = VecDeque::new();
    q.push_back(vec![]);
    let mut seen = 0usize;
    let mut best: Option<(u128, Vec<Rung>)> = None;
    let mut best_len = usize::MAX;
    while let Some(w) = q.pop_front() {
        if w.len() > best_len {
            break;
        }
        seen += 1;
        if seen > cap {
            break;
        }
        let top = top_labels(&lam, &w);
        if top == mu && reduced_for(s, &w, t).ok() {
            let f = crate::schur::flows(&lam, &w);
            if best.as_ref().map_or(true, |b| f < b.0) {
                best = Some((f, w.clone()));
                best_len = w.len();
            }
            continue;
        }
        if w.len() >= max_len {
            continue;
        }
        for i in 0..n.saturating_sub(1) {
            for dir in [Dir::E, Dir::F] {
                let giver = if dir == Dir::E { top[i + 1] } else { top[i] } as usize;
                for r in 1..=giver {
                    let mut w2 = w.clone();
                    w2.push(Rung { i, dir, r });
                    q.push_back(w2);
                }
            }
        }
    }
    best.map(|b| b.1)
}

// ---------------------------------------------------------------------------------------------
// cells, ψ, Grassmannian reduction

/// Nonzero cells of a table in row-major order: (row, column, size).
pub fn cells(t: &Table) -> Vec<(usize, usize, u8)> {
    let (n, m) = (t.nr as usize, t.nc as usize);
    let mut out = vec![];
    for i in 0..n {
        for j in 0..m {
            let a = t.get(i, j);
            if a > 0 {
                out.push((i, j, a));
            }
        }
    }
    out
}

/// LR(p): λ refined by the cells, read by rows (CORE Prop 3.1); RR(p): μ refined by the cells,
/// read by columns (Prop 3.3). Compositions of K.
pub fn lr(t: &Table) -> Vec<u8> {
    cells(t).iter().map(|c| c.2).collect()
}
pub fn rr(t: &Table) -> Vec<u8> {
    let (n, m) = (t.nr as usize, t.nc as usize);
    let mut out = vec![];
    for j in 0..m {
        for i in 0..n {
            let a = t.get(i, j);
            if a > 0 {
                out.push(a);
            }
        }
    }
    out
}

/// ψ(p) ∈ S_{#cells} (CORE Def 3.7 / Thm 3.8): psi[k] = position, in column-major order, of the
/// k-th nonzero cell in row-major order. Σ over inverted pairs (k<l, psi[k]>psi[l]) of
/// size_k·size_l = ℓ(p_) = `Table::lmin` (tested).
pub fn psi(t: &Table) -> Vec<u8> {
    let cs = cells(t);
    let mut colmajor: Vec<(usize, usize)> = cs.iter().map(|c| (c.1, c.0)).collect();
    colmajor.sort();
    cs.iter().map(|c| colmajor.iter().position(|&x| x == (c.1, c.0)).unwrap() as u8).collect()
}

/// Weighted length of ψ: Σ over inverted pairs of the product of the cell sizes.
pub fn psi_weighted_len(t: &Table) -> i64 {
    let cs = cells(t);
    let p = psi(t);
    let mut s = 0i64;
    for a in 0..p.len() {
        for b in a + 1..p.len() {
            if p[a] > p[b] {
                s += cs[a].2 as i64 * cs[b].2 as i64;
            }
        }
    }
    s
}

/// SLL Thm 2.6 at rung level (D1 §1.7). For a step t' -> t of the rung ρ on uprights K = {i, i+1}:
/// Q = the (λ, μ'_K) coset containing t' and t (columns i, i+1 merged), L = RR(Q),
/// z = the (λ, L) coset with z_ = Q_, m' = the (L, μ') coset of Q_^{-1} t'_, m'' = the (L, μ)
/// coset of Q_^{-1} t_. Then t'_ = Q_ m'_ and t_ = Q_ m''_ with lengths adding (tested), and the
/// local step m' -> m'' of the same rung has the same generating function (Lemma 3.9; tested).
#[derive(Clone, Debug)]
pub struct Grass {
    pub q: Table,
    pub l: Vec<u8>,
    pub z: Table,
    pub m_from: Table,
    pub m_to: Table,
    /// the cells of Q in the merged column (row totals inside K): the local bottom of the window
    pub local_lam: Vec<u8>,
    pub k_uprights: (usize, usize),
}

/// merge columns j, j+1 of a table
pub fn merge_columns(t: &Table, j: usize) -> Table {
    let (n, m) = (t.nr as usize, t.nc as usize);
    let mut e = Vec::with_capacity(n * (m - 1));
    for i in 0..n {
        for jj in 0..m {
            if jj == j {
                e.push(t.get(i, j) + t.get(i, j + 1));
            } else if jj != j + 1 {
                e.push(t.get(i, jj));
            }
        }
    }
    Table { nr: n as u8, nc: (m - 1) as u8, e }
}

pub fn grass(t_from: &Table, t_to: &Table, r: &Rung) -> Grass {
    let k = r.i;
    let q = merge_columns(t_from, k);
    assert_eq!(q, merge_columns(t_to, k), "grass: t and t' are not in the same Q");
    let l = rr(&q);
    let qmin = min_perm(&q);
    let qinv = perm_inv(&qmin);
    let lam = t_from.rows();
    let z = table_of_perm(&qmin, &lam, &l);
    let m_from = table_of_perm(&perm_comp(&qinv, &min_perm(t_from)), &l, &t_from.cols());
    let m_to = table_of_perm(&perm_comp(&qinv, &min_perm(t_to)), &l, &t_to.cols());
    let local_lam: Vec<u8> = (0..q.nr as usize).map(|i| q.get(i, k)).filter(|&a| a > 0).collect();
    Grass { q, l, z, m_from, m_to, local_lam, k_uprights: (k, k + 1) }
}
