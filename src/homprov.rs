//! M4b: Hom provider, stacked-point certification and the G1 target enumerator
//! (SPEC_M4_M6a §4; design note `colored_design/M4b_design.md`).
use crate::coset::CosetSpace;
use crate::field::Field;
use crate::pmat::PMat;
use crate::schur::{Elt, Schur, Tid};
use crate::sexpr::Expr;
use crate::smap::{points, Arena, Eval, ExprId, MapId, MapNode, Point};
use crate::tmpl::{Template, TmplId, TmplSource, Window};
use crate::web::{Dir, Rung};
use std::collections::{BTreeMap, HashMap};
use std::time::Instant;

// ------------------------------------------------------------------ G1 targets

pub struct G1Pair {
    pub mu: Vec<u8>,
    pub p: Tid,
    pub rung: Rung,
    pub targets: Vec<(Tid, i32, usize)>,
    pub identity: bool,
}

fn apply_rung(l: &[u8], r: &Rung) -> Option<Vec<u8>> {
    let (rec, giv) = if r.dir == Dir::E { (r.i, r.i + 1) } else { (r.i + 1, r.i) };
    if giv >= l.len() || (l[giv] as usize) < r.r {
        return None;
    }
    let mut v = l.to_vec();
    v[giv] -= r.r as u8;
    v[rec] += r.r as u8;
    Some(v)
}

/// (block, rung) pairs of every crossing term on lambda (both signs, every position): the first
/// rung from lambda, the second rung from the intermediate block (probe P1's work list)
pub fn g1_work(lam: &[u8]) -> Vec<(Vec<u8>, Rung)> {
    let n = lam.len();
    let mut work: Vec<(Vec<u8>, Rung)> = vec![];
    for i in 0..n.saturating_sub(1) {
        for pos in [true, false] {
            for t in crate::schur::crossing_terms(lam[i], lam[i + 1], pos, i) {
                let mut l = lam.to_vec();
                for r in &t.rungs {
                    if !work.contains(&(l.clone(), *r)) {
                        work.push((l.clone(), *r));
                    }
                    l = apply_rung(&l, r).expect("crossing term rung invalid");
                }
            }
        }
    }
    work
}

/// G1 targets of lambda = s.lam: every (type, rung) pair of the blocks reached by one crossing
/// term (all types of each block, in canonical-basis order), with its summand targets
/// (q, s, mult) from `term_mult`; `identity` marks the xengine identity shortcut (a single
/// summand of multiplicity 1 whose canonical ladder is the word L_p R). Second component: the
/// closure types (all types of the block (lambda, lambda)).
pub fn g1_targets(s: &mut Schur) -> (Vec<G1Pair>, Vec<Tid>) {
    let lam = s.lam.clone();
    let mut out = vec![];
    for (mu, r) in g1_work(&lam) {
        let order = s.canonical_basis(&mu).order.clone();
        for p in order {
            let lp = s.canonical_basis(&mu).ladder_minflow[&p].clone();
            let mut x = lp.clone();
            x.push(r);
            let mult = s.term_mult(p, &[r]);
            let mut targets = vec![];
            for (q, m) in &mult {
                for (i, &c) in m.c.iter().enumerate() {
                    if c != 0 {
                        assert!(c > 0, "negative multiplicity");
                        targets.push((*q, m.lo + i as i32, c as usize));
                    }
                }
            }
            let identity = if mult.len() == 1 && mult[0].1 == crate::schur::Lp::one() {
                let q = mult[0].0;
                let mq = s.table(q).cols();
                s.canonical_basis(&mq).ladder_minflow[&q] == x
            } else {
                false
            };
            out.push(G1Pair { mu: mu.clone(), p, rung: r, targets, identity });
        }
    }
    let order = s.canonical_basis(&lam).order.clone();
    (out, order)
}

// ------------------------------------------------------------------ generator degrees and M

/// R^lambda-generator degrees of Hom(src, tgt): hom_gd(char src, char tgt) on [-4K², 4K²]
/// multiplied by prod_{i, j <= lambda_i} (1 - v^{2j}); asserted finite with nonnegative
/// coefficients (§4.3 item 1)
pub fn hom_gens(s: &mut Schur, src_char: &Elt, tgt_char: &Elt) -> BTreeMap<i32, i128> {
    let lam = s.lam.clone();
    let k: i32 = lam.iter().map(|&a| a as i32).sum::<i32>().max(1);
    let (lo, hi) = (-4 * k * k, 4 * k * k);
    let mut g = s.hom_gd(src_char, tgt_char, lo, hi);
    for &a in &lam {
        for j in 1..=a as usize {
            let st = 2 * j;
            for d in (st..g.len()).rev() {
                g[d] -= g[d - st];
            }
        }
    }
    let mut out = BTreeMap::new();
    for (i, &c) in g.iter().enumerate() {
        assert!(c >= 0, "hom_gens: negative generator count {} in degree {} (Hom not free?)", c, lo + i as i32);
        if c != 0 {
            out.insert(lo + i as i32, c);
        }
    }
    // finiteness: no generator in the top quarter of the window (the series must have ended)
    if let Some((&dmax, _)) = out.iter().next_back() {
        assert!(dmax < hi - 2 * k * k, "hom_gens: generator degrees reach {} (window up to {}): not a finite polynomial", dmax, hi);
    }
    out
}

/// M(src, tgt, d) = max over generator degrees d_j <= d with d - d_j even of
/// emono_count(lambda, (d - d_j)/2) (§4.3 item 2; 0 if no generator qualifies)
pub fn m_points(lam: &[u8], gens: &BTreeMap<i32, i128>, d: i32) -> usize {
    let mut m = 0u128;
    for &dj in gens.keys() {
        if dj <= d && (d - dj) % 2 == 0 {
            m = m.max(crate::schur::emono_count(lam, ((d - dj) / 2) as usize));
        }
    }
    m as usize
}

// ------------------------------------------------------------------ words: commutation and windows

/// Do adjacent rungs r1 (below) and r2 (above) commute as bimodules (degree-0 isomorphism
/// [r1, r2] ≅ [r2, r1])? Disjoint upright pairs always commute; rungs sharing one upright commute
/// when both merge into it or both split from it (associativity of a three-part merge or split).
pub fn rungs_commute(r1: &Rung, r2: &Rung) -> bool {
    let (p1, p2) = (r1.i as i64, r2.i as i64);
    if (p1 - p2).abs() >= 2 {
        return true;
    }
    if p1 == p2 {
        return false;
    }
    // shared upright s; on its left pair a rung F merges into s and E splits from s; on its right
    // pair E merges into s and F splits from s
    let (left, right) = if p2 == p1 + 1 { (r1, r2) } else { (r2, r1) };
    left.dir != right.dir
}

/// All words reachable from `w` by swapping adjacent commuting rungs (breadth first, capped)
pub fn commutation_class(w: &[Rung], cap: usize) -> Vec<Vec<Rung>> {
    let mut seen: std::collections::HashSet<Vec<Rung>> = std::collections::HashSet::new();
    let mut out = vec![w.to_vec()];
    seen.insert(w.to_vec());
    let mut k = 0;
    while k < out.len() && out.len() < cap {
        let cur = out[k].clone();
        k += 1;
        for i in 0..cur.len().saturating_sub(1) {
            if rungs_commute(&cur[i], &cur[i + 1]) {
                let mut nw = cur.clone();
                nw.swap(i, i + 1);
                if seen.insert(nw.clone()) {
                    out.push(nw);
                    if out.len() >= cap {
                        break;
                    }
                }
            }
        }
    }
    out
}

/// The non-common window of two words with equal bottom and top: common prefix length `at`,
/// common suffix length, window uprights first..first+labels.len(), window words (relative
/// uprights) and fibres
#[derive(Clone, Debug)]
pub struct WinSplit {
    pub at: usize,
    pub suffix: usize,
    pub first: usize,
    pub labels: Vec<u8>,
    pub a: Vec<Rung>,
    pub b: Vec<Rung>,
    pub fa: u128,
    pub fb: u128,
}

pub fn win_split(lam: &[u8], a: &[Rung], b: &[Rung]) -> WinSplit {
    let mut pre = 0;
    while pre < a.len() && pre < b.len() && a[pre] == b[pre] {
        pre += 1;
    }
    let mut suf = 0;
    while suf < a.len() - pre && suf < b.len() - pre && a[a.len() - 1 - suf] == b[b.len() - 1 - suf] {
        suf += 1;
    }
    let (wa, wb) = (&a[pre..a.len() - suf], &b[pre..b.len() - suf]);
    let mut lv = lam.to_vec();
    for r in &a[..pre] {
        lv = apply_rung(&lv, r).expect("win_split: invalid word");
    }
    if wa.is_empty() && wb.is_empty() {
        return WinSplit { at: pre, suffix: suf, first: 0, labels: vec![], a: vec![], b: vec![], fa: 1, fb: 1 };
    }
    let lo = wa.iter().chain(wb).map(|r| r.i).min().unwrap();
    let hi = wa.iter().chain(wb).map(|r| r.i + 1).max().unwrap();
    let labels = lv[lo..=hi].to_vec();
    let sh = |rs: &[Rung]| -> Vec<Rung> { rs.iter().map(|r| Rung { i: r.i - lo, dir: r.dir, r: r.r }).collect() };
    let (ra, rb) = (sh(wa), sh(wb));
    let fa = crate::schur::flows(&labels, &ra);
    let fb = crate::schur::flows(&labels, &rb);
    WinSplit { at: pre, suffix: suf, first: lo, labels, a: ra, b: rb, fa, fb }
}

// ------------------------------------------------------------------ decategorified reduction planner

/// A window of a word: rungs at..at+len, all on uprights first..first+w (relative rungs)
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct WinPos {
    pub at: usize,
    pub len: usize,
    pub first: usize,
    pub w: usize,
}

/// labels after the first `k` rungs of `w` from `lam`
pub fn level_of(lam: &[u8], w: &[Rung], k: usize) -> Vec<u8> {
    let mut l = lam.to_vec();
    for r in &w[..k] {
        l = apply_rung(&l, r).expect("level_of: invalid word");
    }
    l
}

/// relative window word of a position
pub fn window_word(w: &[Rung], p: &WinPos) -> Vec<Rung> {
    w[p.at..p.at + p.len].iter().map(|r| Rung { i: r.i - p.first, dir: r.dir, r: r.r }).collect()
}

/// Schur instances per window labels (decategorified planner)
#[derive(Default)]
pub struct SchurCache {
    pub by_lab: std::collections::HashMap<Vec<u8>, Schur>,
}
impl SchurCache {
    pub fn get(&mut self, lab: &[u8]) -> &mut Schur {
        self.by_lab.entry(lab.to_vec()).or_insert_with(|| Schur::new(lab))
    }
    /// decomposition of a window word into canonical ladders of its block: (canonical word, shift, count)
    pub fn expand_word(&mut self, lab: &[u8], rs: &[Rung]) -> Vec<(Vec<Rung>, i32, usize, Tid)> {
        let s = self.get(lab);
        let ch = s.ladder_char(rs);
        let mut out = vec![];
        for (t, m) in s.expand(&ch) {
            let mu = s.table(t).cols();
            let lt = s.canonical_basis(&mu).ladder_minflow[&t].clone();
            for (i, &c) in m.c.iter().enumerate() {
                if c != 0 {
                    out.push((lt.clone(), m.lo + i as i32, c as usize, t));
                }
            }
        }
        out
    }
}

// ------------------------------------------------------------------ exact-call size estimates

/// histogram of the fibre degrees of a word: entry k counts fibre basis elements of
/// unnormalized degree 2k (product over rungs of the Gaussian binomial [giver choose r]_q)
pub fn fdeg_hist(labels: &[u8], rungs: &[Rung]) -> Vec<u128> {
    let mut h: Vec<u128> = vec![1];
    let mut lab = labels.to_vec();
    for r in rungs {
        let giv = if r.dir == Dir::E { r.i + 1 } else { r.i };
        let b = lab[giv] as usize;
        let g = gauss_binom(b, r.r);
        let mut out = vec![0u128; h.len() + g.len() - 1];
        for (i, &x) in h.iter().enumerate() {
            if x == 0 {
                continue;
            }
            for (j, &y) in g.iter().enumerate() {
                out[i + j] += x * y;
            }
        }
        h = out;
        lab = apply_rung(&lab, r).expect("fdeg_hist: invalid word");
    }
    h
}

/// coefficients of the Gaussian binomial [n choose k]_q (partitions in a k × (n-k) box by size)
pub fn gauss_binom(n: usize, k: usize) -> Vec<u128> {
    if k > n {
        return vec![0];
    }
    // standard recursion: [n,k] = [n-1,k-1] + q^k [n-1,k]
    fn rec(n: usize, k: usize, memo: &mut std::collections::HashMap<(usize, usize), Vec<u128>>) -> Vec<u128> {
        if k == 0 || k == n {
            return vec![1];
        }
        if let Some(v) = memo.get(&(n, k)) {
            return v.clone();
        }
        let a = rec(n - 1, k - 1, memo);
        let b = rec(n - 1, k, memo);
        let mut out = vec![0u128; (k * (n - k)) + 1];
        for (i, &x) in a.iter().enumerate() {
            out[i] += x;
        }
        for (i, &x) in b.iter().enumerate() {
            out[i + k] += x;
        }
        memo.insert((n, k), out.clone());
        out
    }
    rec(n, k, &mut std::collections::HashMap::new())
}

/// normalization shift of a word: -sum_rungs r (b - r), b = giver label
pub fn word_shift(labels: &[u8], rungs: &[Rung]) -> i32 {
    let mut lab = labels.to_vec();
    let mut s = 0i32;
    for r in rungs {
        let giv = if r.dir == Dir::E { r.i + 1 } else { r.i };
        s -= r.r as i32 * (lab[giv] as i32 - r.r as i32);
        lab = apply_rung(&lab, r).expect("word_shift: invalid word");
    }
    s
}

/// Number of unknowns of `sbim::hom_space(src, tgt, u)` for the normalized degree `d`:
/// Σ_{β, β'} #monos(u + fdeg_src[β] - fdeg_tgt[β']) with u = d - shift(tgt) + shift(src)
pub fn hom_unknowns(labels: &[u8], src: &[Rung], tgt: &[Rung], d: i32) -> u128 {
    let u = d - word_shift(labels, tgt) + word_shift(labels, src);
    let hs = fdeg_hist(labels, src);
    let ht = fdeg_hist(labels, tgt);
    let mut memo: std::collections::HashMap<i32, u128> = std::collections::HashMap::new();
    let mut tot: u128 = 0;
    for (i, &a) in hs.iter().enumerate() {
        if a == 0 {
            continue;
        }
        for (j, &b) in ht.iter().enumerate() {
            if b == 0 {
                continue;
            }
            let dd = u + 2 * i as i32 - 2 * j as i32;
            if dd < 0 || dd % 2 != 0 {
                continue;
            }
            let m = *memo.entry(dd).or_insert_with(|| crate::schur::emono_count(labels, (dd / 2) as usize));
            tot = tot.saturating_add(a.saturating_mul(b).saturating_mul(m));
        }
    }
    tot
}

// ------------------------------------------------------------------ interface types (§4.2)

#[derive(Clone, Debug)]
pub struct Budget {
    pub l_max: usize,
    pub c_max: usize,
    pub t_max_s: f64,
    pub m_max_bytes: u64,
    /// 0 = 4·max(flows(src), flows(tgt)) at call time
    pub flows_max: u64,
    pub index_max_bytes: u64,
}
impl Default for Budget {
    fn default() -> Self {
        Budget { l_max: 6, c_max: 4096, t_max_s: 600.0, m_max_bytes: 1 << 30, flows_max: 0, index_max_bytes: 256 << 20 }
    }
}

#[derive(Clone, Debug, Default)]
pub struct SearchStats {
    pub candidates: usize,
    pub evaluated: usize,
    pub secs: f64,
    pub tier_hits: [usize; 7],
    pub max_flows: u64,
    pub points: usize,
    pub visited_words: usize,
    pub indexed_words: usize,
    pub peak_index_bytes: u64,
    pub peak_arena_bytes: u64,
    /// computed (not cached) whole-ladder hom_space calls of the Exact backend
    pub exact_calls: usize,
    pub exact_secs: f64,
    /// computed (not cached) exact window calls made by Search (splits, commutations, rex moves,
    /// leaf isomorphisms, S6 windows); counted apart from the Exact backend (review tests F5)
    pub window_calls: usize,
    pub window_secs: f64,
    pub detect: &'static str, /* "witness" | "orbits" | "more-points" */
}
impl SearchStats {
    /// add the counters of another call (Chain)
    pub fn absorb(&mut self, o: &SearchStats) {
        self.candidates += o.candidates;
        self.evaluated += o.evaluated;
        self.secs += o.secs;
        for k in 0..7 {
            self.tier_hits[k] += o.tier_hits[k];
        }
        self.max_flows = self.max_flows.max(o.max_flows);
        self.points = self.points.max(o.points);
        self.visited_words += o.visited_words;
        self.indexed_words += o.indexed_words;
        self.peak_index_bytes = self.peak_index_bytes.max(o.peak_index_bytes);
        self.peak_arena_bytes = self.peak_arena_bytes.max(o.peak_arena_bytes);
        self.exact_calls += o.exact_calls;
        self.exact_secs += o.exact_secs;
        self.window_calls += o.window_calls;
        self.window_secs += o.window_secs;
        if !o.detect.is_empty() {
            self.detect = o.detect;
        }
    }
}

/// Why a provider call stopped without a certified answer. Exhausted = candidate stream ended with
/// budget left (completeness signal); Time/Memory = budget overrun (engineering); Refused = an
/// admission rule declined the call (Exact).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MissKind {
    Exhausted,
    Time,
    Memory,
    Refused,
}
#[derive(Clone, Debug)]
pub struct Missing {
    pub what: String,
    pub degree: i32,
    pub stratum: Option<String>,
    pub got: usize,
    pub want: usize,
    pub kind: MissKind,
}

/// form[j][i] = top scalar of pis[j]∘iotas[i], rank = mult
#[derive(Clone, Debug)]
pub struct Pairs<F> {
    pub iotas: Vec<MapId>,
    pub pis: Vec<MapId>,
    pub form: Vec<Vec<F>>,
}

pub trait HomProvider<F: Field> {
    fn name(&self) -> &'static str;
    /// certified k-basis of Hom^d(src, tgt) (|S| = Hom formula, stacked rank, fresh residual)
    fn full(&mut self, ar: &mut Arena<F>, s: &mut Schur, src: ExprId, tgt: ExprId, d: i32, b: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>>;
    /// `mult` pairs (iota: L_q -> X degree s, pi: X -> L_q degree -s) whose intersection form at a
    /// top-stratum coset of q has rank `mult`; e_q applied on both sides when present
    #[allow(clippy::too_many_arguments)]
    fn pairs(&mut self, ar: &mut Arena<F>, s: &mut Schur, x: ExprId, q: Tid, lq: ExprId, eq: Option<MapId>, sdeg: i32, mult: usize, b: &Budget) -> Result<(Pairs<F>, SearchStats), Vec<Missing>>;
    /// n = c_diag(b_p)(1) maps R^lambda -> L_p (ins) or L_p -> R^lambda (outs) with degree multiset
    /// c_diag(b_p) and nonsingular identity-coset matrix at two points (C2 M2)
    #[allow(clippy::too_many_arguments)]
    fn closure(&mut self, ar: &mut Arena<F>, s: &mut Schur, lp: ExprId, p: Tid, ep: Option<MapId>, outs: bool, b: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>>;
    /// a spanning set of End^0(L) (certified by stacked rank = Hom formula)
    fn end0(&mut self, ar: &mut Arena<F>, s: &mut Schur, l: ExprId, b: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>>;
    /// the counters of the last call if it failed (None after a success; Chain absorbs them so
    /// that the time and candidates of failed providers are reported). Additive to the SPEC
    /// interface.
    fn failed_stats(&self) -> Option<SearchStats> {
        None
    }
}

fn miss(what: &str, degree: i32, stratum: Option<String>, got: usize, want: usize, kind: MissKind) -> Vec<Missing> {
    vec![Missing { what: what.to_string(), degree, stratum, got, want, kind }]
}

// ------------------------------------------------------------------ small helpers

fn eid<F: Field>(ar: &mut Arena<F>, w: &[Rung]) -> Result<ExprId, String> {
    let lam = ar.lam().to_vec();
    let e = Expr::new(&lam, w)?;
    Ok(ar.expr(e))
}
fn word<F: Field>(ar: &Arena<F>, e: ExprId) -> Vec<Rung> {
    ar.expr_of(e).rungs.clone()
}
fn top_of(lam: &[u8], w: &[Rung]) -> Vec<u8> {
    level_of(lam, w, w.len())
}
/// Some(c) if m = c·id (zero gives Some(0)), None otherwise
pub fn pmat_scalar<F: Field>(m: &PMat<F>) -> Option<F> {
    if m.is_zero() {
        return Some(F::zero());
    }
    if m.nrows != m.ncols {
        return None;
    }
    let mut c: Option<F> = None;
    for (j, col) in m.cols.iter().enumerate() {
        if col.len() != 1 || col[0].0 as usize != j {
            return None;
        }
        let p = &col[0].1;
        if p.t.len() != 1 || p.t[0].0 != 0 {
            return None;
        }
        match &c {
            None => c = Some(p.t[0].1.clone()),
            Some(x) if *x == p.t[0].1 => {}
            Some(_) => return None,
        }
    }
    c
}
fn invert<F: Field>(a: &[Vec<F>]) -> Option<Vec<Vec<F>>> {
    crate::local::inverse(a)
}
/// coset rank of the witness rep(t) of an interned table in its block
fn rep_rank(s: &Schur, lam: &[u8], t: Tid) -> u64 {
    let tab = s.table(t).clone();
    let cs = CosetSpace::new(lam, &tab.cols());
    cs.rank(&cs.rep(&tab))
}
/// coset rank of rep(p) for a type of (lambda, lambda) (the top-stratum coset of L_p)
fn identity_coset_top(lam: &[u8], s: &Schur, p: Tid) -> u64 {
    rep_rank(s, lam, p)
}
fn identity_coset(lam: &[u8]) -> u64 {
    let cs = CosetSpace::new(lam, lam);
    cs.rank(&cs.roots.bottom())
}
fn elapsed(t0: &Instant) -> f64 {
    t0.elapsed().as_secs_f64()
}

// ------------------------------------------------------------------ certification (§4.3)

/// genuineness (§4.3 item 9): every leaf of the DAG is an exact template (hom_space, Rickard or a
/// certified added template), an Elem with a recorded constant, or a Dot
pub fn assert_genuine<F: Field>(ar: &Arena<F>, m: MapId) -> Result<(), String> {
    ar.assert_genuine(m)
}

/// tables t with c_t(char src)(1) · c_t(char tgt)(1) > 0, with both coefficients
fn common_tables(chs: &Elt, cht: &Elt) -> Vec<(Tid, i128, i128)> {
    let mut out = vec![];
    for (t, c) in &chs.t {
        let a = c.at1();
        let b = cht.coef(*t).at1();
        if a > 0 && b > 0 {
            out.push((*t, a, b));
        }
    }
    out
}

/// The detecting set of §4.3 item 3 with its point count.
#[derive(Clone, Debug)]
pub struct DetectPlan {
    pub cosets: Vec<u64>,
    pub npoints: usize,
    pub kind: &'static str,
    pub dim: usize,
    pub m: usize,
    pub rows_per_point: usize,
}

/// witness cosets rep(t) of the tables t with c_t(char src)(1) c_t(char tgt)(1) > 0
pub fn detecting_cosets<F: Field>(ar: &Arena<F>, s: &mut Schur, src: ExprId, tgt: ExprId) -> Vec<u64> {
    let lam = ar.lam().to_vec();
    assert_eq!(s.lam, lam, "detecting_cosets: Schur lambda != arena lambda");
    let (ws, wt) = (word(ar, src), word(ar, tgt));
    let chs = s.ladder_char(&ws);
    let cht = s.ladder_char(&wt);
    let mut out: Vec<u64> = common_tables(&chs, &cht).iter().map(|&(t, _, _)| rep_rank(s, &lam, t)).collect();
    out.sort_unstable();
    out
}

/// dimension, M, detecting set and number of points for Hom^d(src, tgt) (§4.3 items 1-3)
pub fn detect_plan<F: Field>(ar: &Arena<F>, s: &mut Schur, src: ExprId, tgt: ExprId, d: i32) -> DetectPlan {
    let lam = ar.lam().to_vec();
    assert_eq!(s.lam, lam, "detect_plan: Schur lambda != arena lambda");
    let (ws, wt) = (word(ar, src), word(ar, tgt));
    let chs = s.ladder_char(&ws);
    let cht = s.ladder_char(&wt);
    let dim = s.hom_gd(&chs, &cht, d, d)[0] as usize;
    let gens = hom_gens(s, &chs, &cht);
    let m = m_points(&lam, &gens, d);
    let ts = common_tables(&chs, &cht);
    let ew: usize = ts.iter().map(|&(_, a, b)| (a * b) as usize).sum();
    let mu = top_of(&lam, &ws);
    let cs = CosetSpace::new(&lam, &mu);
    let mut wit: Vec<u64> = ts.iter().map(|&(t, _, _)| cs.rank(&cs.rep(s.table(t)))).collect();
    wit.sort_unstable();
    if (m + 2) * ew >= dim {
        return DetectPlan { cosets: wit, npoints: m + 2, kind: "witness", dim, m, rows_per_point: ew };
    }
    let orb_rows: usize = ts.iter().map(|&(t, a, b)| (a * b) as usize * s.ncosets_of(t) as usize).sum();
    if (m + 2) * orb_rows >= dim {
        let mut cos: Vec<u64> = vec![];
        for &(t, _, _) in &ts {
            cos.extend(cs.cosets_of(s.table(t)));
        }
        cos.sort_unstable();
        return DetectPlan { cosets: cos, npoints: m + 2, kind: "orbits", dim, m, rows_per_point: orb_rows };
    }
    let np = dim.div_ceil(ew.max(1)) + 2;
    DetectPlan { cosets: wit, npoints: np, kind: "more-points", dim, m, rows_per_point: ew }
}

/// rank of the maps' blocks at the cosets, stacked over the points (one row per map)
pub fn stacked_rank<F: Field>(ar: &Arena<F>, ms: &[MapId], pts: &[Point<F>], cosets: &[u64]) -> usize {
    if ms.is_empty() {
        return 0;
    }
    let rows = Eval::stacked(ar, ms, pts, cosets);
    let n = rows[0].len();
    crate::dense::rank(&rows, n)
}

/// §4.3 item 4: |S| = [v^d] hom_gd, stacked rank |S| at the detecting set (counting bound asserted
/// first), rank still |S| with 2 fresh points appended, every map genuine with the right ends and
/// degree. Returns the detecting plan used.
pub fn certify_full_plan<F: Field>(ar: &Arena<F>, s: &mut Schur, src: ExprId, tgt: ExprId, d: i32, ms: &[MapId], seed: u64) -> Result<DetectPlan, String> {
    let plan = detect_plan(ar, s, src, tgt, d);
    certify_full_at(ar, s, src, tgt, d, ms, seed, plan)
}

/// `certify_full_plan` at a given detecting plan (tests force a plan, B4b); the counting bound
/// rows_per_point · npoints >= dim is asserted before any ranking
#[allow(clippy::too_many_arguments)]
pub fn certify_full_at<F: Field>(ar: &Arena<F>, s: &mut Schur, src: ExprId, tgt: ExprId, d: i32, ms: &[MapId], seed: u64, plan: DetectPlan) -> Result<DetectPlan, String> {
    for &m in ms {
        let i = ar.info(m);
        if (i.src, i.tgt, i.deg) != (src, tgt, d) {
            return Err(format!("certify_full: map {} is {}->{} of degree {}, expected {}->{} degree {}", m, i.src, i.tgt, i.deg, src, tgt, d));
        }
        assert_genuine(ar, m)?;
    }
    if ms.len() != plan.dim {
        return Err(format!("certify_full: {} maps, dim Hom^{} = {}", ms.len(), d, plan.dim));
    }
    if plan.dim == 0 {
        return Ok(plan);
    }
    if plan.rows_per_point * plan.npoints < plan.dim {
        return Err(format!("certify_full: counting bound violated: {} rows per point x {} points < dim {} ({})", plan.rows_per_point, plan.npoints, plan.dim, plan.kind));
    }
    let pts = points::<F>(ar.k(), seed, plan.npoints + 2, 0);
    let r = stacked_rank(ar, ms, &pts[..plan.npoints], &plan.cosets);
    if r != plan.dim {
        // Deterministic fallback (a too-small detecting set only causes false negatives): the
        // full W_lambda-orbits of the common tables, with M + 2 points and then twice as many.
        // Stacked rank = |S| at any coset set and points proves independence, so the fallback
        // cannot certify a dependent family. Measured: hom60 case (1,3) F0(1)E0(3)F0(2) ->
        // E0(3)F0(2)F0(1) degree 6, witness rank 141 < 145 at 18 points.
        let first = format!("stacked rank {} < dim {} at {} points ({} detecting cosets, {})", r, plan.dim, plan.npoints, plan.cosets.len(), plan.kind);
        let (cos, orb_rows, np0) = orbit_fallback(ar, s, src, tgt, &plan);
        for (k, np) in [np0, 2 * np0].into_iter().enumerate() {
            let pts = points::<F>(ar.k(), seed ^ 0x0b17_0000 ^ k as u64, np + 2, 0);
            if stacked_rank(ar, ms, &pts[..np], &cos) == plan.dim {
                if stacked_rank(ar, ms, &pts, &cos) != plan.dim {
                    return Err(format!("certify_full: residual rank != {} (orbit fallback)", plan.dim));
                }
                return Ok(DetectPlan { cosets: cos, npoints: np, kind: if k == 0 { "orbits-fallback" } else { "orbits-fallback-2x" }, dim: plan.dim, m: plan.m, rows_per_point: orb_rows });
            }
        }
        return Err(format!("certify_full: {}; orbit fallback ({} cosets, {} and {} points) also below dim", first, cos.len(), np0, 2 * np0));
    }
    // residual: the 2 fresh points appended (columns only; the rank cannot drop, kept as specified)
    let r2 = stacked_rank(ar, ms, &pts, &plan.cosets);
    if r2 != plan.dim {
        return Err(format!("certify_full: residual rank {} != {}", r2, plan.dim));
    }
    Ok(plan)
}

/// the orbit fallback of `certify_full_plan`: the full W_λ-orbits of the common tables
/// (deduplicated), their rows per point, and the first point count max(M + 2, ⌈dim/rows⌉ + 2)
/// (the second attempt doubles it)
pub fn orbit_fallback<F: Field>(ar: &Arena<F>, s: &mut Schur, src: ExprId, tgt: ExprId, plan: &DetectPlan) -> (Vec<u64>, usize, usize) {
    let (ws, wt) = (word(ar, src), word(ar, tgt));
    let lam = ar.lam().to_vec();
    let (chs, cht) = (s.ladder_char(&ws), s.ladder_char(&wt));
    let ts = common_tables(&chs, &cht);
    let mu = top_of(&lam, &ws);
    let cs = CosetSpace::new(&lam, &mu);
    let mut cos: Vec<u64> = vec![];
    for &(t, _, _) in &ts {
        cos.extend(cs.cosets_of(s.table(t)));
    }
    cos.sort_unstable();
    cos.dedup();
    let orb_rows: usize = ts.iter().map(|&(t, a, b)| (a * b) as usize * s.ncosets_of(t) as usize).sum();
    let np0 = (plan.m + 2).max(plan.dim.div_ceil(orb_rows.max(1)) + 2);
    (cos, orb_rows, np0)
}

pub fn certify_full<F: Field>(ar: &Arena<F>, s: &mut Schur, src: ExprId, tgt: ExprId, d: i32, ms: &[MapId], seed: u64) -> Result<(), String> {
    certify_full_plan(ar, s, src, tgt, d, ms, seed).map(|_| ())
}

/// §4.3 item 5: a single map certifies a 1-dimensional space when it is nonzero at two points
pub fn certify_one<F: Field>(ar: &Arena<F>, m: MapId, seed: u64) -> Result<(), String> {
    assert_genuine(ar, m)?;
    for p in points::<F>(ar.k(), seed, 2, 0) {
        if Eval::new(ar, p).is_zero(m) {
            return Err(format!("certify_one: map {} vanishes at a point", m));
        }
    }
    Ok(())
}

/// form[j][i] = (block(π_j, c) · block(ι_i, c))[0][0] at the coset c (L_q has one flow there)
pub fn form_at<F: Field>(ar: &Arena<F>, iotas: &[MapId], pis: &[MapId], c: u64, pt: &Point<F>) -> Result<Vec<Vec<F>>, String> {
    let mut ev = Eval::new(ar, pt.clone());
    let ib: Vec<std::rc::Rc<crate::smap::Block<F>>> = iotas.iter().map(|&m| ev.block(m, c)).collect();
    let pb: Vec<std::rc::Rc<crate::smap::Block<F>>> = pis.iter().map(|&m| ev.block(m, c)).collect();
    let mut out = vec![vec![F::zero(); iotas.len()]; pis.len()];
    for (j, p) in pb.iter().enumerate() {
        if p.rows != 1 {
            return Err(format!("form: L_q has {} flows at the coset {} (expected 1)", p.rows, c));
        }
        for (i, b) in ib.iter().enumerate() {
            if b.cols != 1 || b.rows != p.cols {
                return Err(format!("form: block shapes {}x{} and {}x{}", p.rows, p.cols, b.rows, b.cols));
            }
            let mut v = F::zero();
            for k in 0..b.rows as usize {
                v.fma(&p.a[k], &b.a[k]);
            }
            out[j][i] = v;
        }
    }
    Ok(out)
}

/// independent columns I and rows J of a form of rank `want` (None if the rank differs)
fn select_form<F: Field>(g: &[Vec<F>], ni: usize, want: usize) -> Option<(Vec<usize>, Vec<usize>)> {
    if g.is_empty() || ni == 0 {
        return if want == 0 { Some((vec![], vec![])) } else { None };
    }
    let mut gr = g.to_vec();
    let icols = crate::dense::rref(&mut gr, ni);
    if icols.len() != want {
        return None;
    }
    let mut gt: Vec<Vec<F>> = (0..ni).map(|i| (0..g.len()).map(|j| g[j][i].clone()).collect()).collect();
    let jrows = crate::dense::rref(&mut gt, g.len());
    Some((icols, jrows))
}

/// apply the idempotent of q on both sides (ι ↦ ι∘e_q, π ↦ e_q∘π)
fn with_idem<F: Field>(ar: &mut Arena<F>, eq: Option<MapId>, iotas: &mut [MapId], pis: &mut [MapId]) -> Result<(), String> {
    if let Some(e) = eq {
        for m in iotas.iter_mut() {
            *m = ar.add(MapNode::Comp(*m, e), "iota e_q")?;
        }
        for m in pis.iter_mut() {
            *m = ar.add(MapNode::Comp(e, *m), "e_q pi")?;
        }
    }
    Ok(())
}

/// ends and degree of a map
fn check_info<F: Field>(ar: &Arena<F>, m: MapId, src: ExprId, tgt: ExprId, deg: i32, what: &str) -> Result<(), String> {
    let i = ar.info(m);
    if (i.src, i.tgt, i.deg) != (src, tgt, deg) {
        return Err(format!("{}: map {} is {}->{} of degree {}, expected {}->{} degree {}", what, m, i.src, i.tgt, i.deg, src, tgt, deg));
    }
    Ok(())
}

/// pairs certification (§4.3 item 6): every ι is L_q -> X of degree sdeg and every π is X -> L_q
/// of degree -sdeg (review F3), rank `mult`, the same matrix at a second point, genuine maps;
/// returns exactly `mult` ι and `mult` π with an invertible form
#[allow(clippy::too_many_arguments)]
pub fn certify_pairs<F: Field>(ar: &Arena<F>, x: ExprId, lq: ExprId, sdeg: i32, iotas: &[MapId], pis: &[MapId], c: u64, mult: usize, seed: u64) -> Result<Pairs<F>, String> {
    for &m in iotas {
        check_info(ar, m, lq, x, sdeg, "pairs (iota)")?;
    }
    for &m in pis {
        check_info(ar, m, x, lq, -sdeg, "pairs (pi)")?;
    }
    let pts = points::<F>(ar.k(), seed, 2, 0);
    let g1 = form_at(ar, iotas, pis, c, &pts[0])?;
    let (ii, jj) = select_form(&g1, iotas.len(), mult).ok_or_else(|| format!("form rank != {}", mult))?;
    let io: Vec<MapId> = ii.iter().map(|&i| iotas[i]).collect();
    let pi: Vec<MapId> = jj.iter().map(|&j| pis[j]).collect();
    let f1 = form_at(ar, &io, &pi, c, &pts[0])?;
    let f2 = form_at(ar, &io, &pi, c, &pts[1])?;
    if f1 != f2 {
        return Err("point-dependent map: the intersection form differs at a second point".into());
    }
    if crate::dense::rank(&f1, mult) != mult {
        return Err("selected form is singular".into());
    }
    for &m in io.iter().chain(&pi) {
        assert_genuine(ar, m)?;
    }
    Ok(Pairs { iotas: io, pis: pi, form: f1 })
}

/// row-major square matrix product
fn sq_mul<F: Field>(a: &[F], b: &[F], n: usize) -> Vec<F> {
    let mut out = vec![F::zero(); n * n];
    for i in 0..n {
        for k in 0..n {
            let x = &a[i * n + k];
            if x.is_zero() {
                continue;
            }
            for j in 0..n {
                out[i * n + j].fma(x, &b[k * n + j]);
            }
        }
    }
    out
}

/// The e_q / e_p contract (SPEC §5.2 item 6, checked cheaply here): e is a genuine degree-0
/// endomorphism of L, b_p occurs in char L exactly once in degree 0, and at two points, at the
/// witness coset rep(t) of every table t of char L, block(e) is idempotent of rank c_t(b_p)(1);
/// at rep(p) (one flow) its scalar is 1. With these ranks Krull-Schmidt gives Im e ≅ B_p.
pub fn check_idempotent<F: Field>(ar: &Arena<F>, s: &mut Schur, l: ExprId, p: Tid, e: MapId, seed: u64) -> Result<(), String> {
    let lam = ar.lam().to_vec();
    check_info(ar, e, l, l, 0, "idempotent")?;
    assert_genuine(ar, e)?;
    let mu = s.table(p).cols();
    let ch = s.ladder_char(&word(ar, l));
    let mult = s.expand(&ch).into_iter().find(|t| t.0 == p).map(|t| t.1).unwrap_or_else(crate::schur::Lp::zero);
    if mult != crate::schur::Lp::one() {
        return Err(format!("idempotent: b_p occurs in char L with multiplicity {:?}, expected 1", mult));
    }
    let bp = s.canonical_basis(&mu).b[&p].clone();
    let tables: Vec<Tid> = ch.t.iter().map(|t| t.0).collect();
    for pt in points::<F>(ar.k(), seed, 2, 0) {
        let mut ev = Eval::new(ar, pt);
        for &t in &tables {
            let c = rep_rank(s, &lam, t);
            let b = ev.block(e, c);
            if b.rows != b.cols {
                return Err(format!("idempotent: block at rep({}) is {}x{}", s.table(t).fmt(), b.rows, b.cols));
            }
            let n = b.rows as usize;
            if sq_mul(&b.a, &b.a, n) != b.a {
                return Err(format!("idempotent: e^2 != e at rep({})", s.table(t).fmt()));
            }
            let rows: Vec<Vec<F>> = (0..n).map(|i| b.a[i * n..(i + 1) * n].to_vec()).collect();
            let r = if n == 0 { 0 } else { crate::dense::rank(&rows, n) } as i128;
            let want = bp.coef(t).at1();
            if r != want {
                return Err(format!("idempotent: rank {} at rep({}), c_t(b_p)(1) = {}", r, s.table(t).fmt(), want));
            }
            if t == p && (n != 1 || b.a[0] != F::one()) {
                return Err(format!("idempotent: top-coset scalar at rep({}) is not 1 ({} flows)", s.table(p).fmt(), n));
            }
        }
    }
    Ok(())
}

/// Preconditions of `pairs`, checked by every backend before any work (review F1, F2): the ends
/// have the top of q, L_q is indecomposable of type q (char L_q = b_q) unless an idempotent e_q
/// is given (then `check_idempotent`), and `mult` = [v^sdeg] of the multiplicity of b_q in char X.
/// Err = Refused.
#[allow(clippy::too_many_arguments)]
pub fn pairs_precheck<F: Field>(ar: &Arena<F>, s: &mut Schur, x: ExprId, q: Tid, lq: ExprId, eq: Option<MapId>, sdeg: i32, mult: usize) -> Result<(), String> {
    let lam = ar.lam().to_vec();
    if s.lam != lam {
        return Err(format!("Schur lambda {:?} != arena lambda {:?}", s.lam, lam));
    }
    let (xw, qw) = (word(ar, x), word(ar, lq));
    let mu = s.table(q).cols();
    if top_of(&lam, &qw) != mu || top_of(&lam, &xw) != mu {
        return Err(format!("tops: L_q {:?}, X {:?}, q is a table of {:?}", top_of(&lam, &qw), top_of(&lam, &xw), mu));
    }
    match eq {
        None => {
            let chq = s.ladder_char(&qw);
            let ex = s.expand(&chq);
            if ex.len() != 1 || ex[0].0 != q || ex[0].1 != crate::schur::Lp::one() {
                let parts: Vec<String> = ex.iter().map(|(t, c)| format!("{}{{{}}}", s.table(*t).fmt(), c.fmt_ec())).collect();
                return Err(format!("L_q decomposable: e_q required (char L_q = {} for q = {})", parts.join(" + "), s.table(q).fmt()));
            }
        }
        Some(e) => check_idempotent(ar, s, lq, q, e, 0x5eed_e000).map_err(|m| format!("e_q rejected: {}", m))?,
    }
    let chx = s.ladder_char(&xw);
    let m = s.expand(&chx).into_iter().find(|t| t.0 == q).map(|t| t.1.coef(sdeg)).unwrap_or(0);
    if m != mult as i128 {
        return Err(format!("mult {} != [v^{}] of the multiplicity of b_q in char X = {}", mult, sdeg, m));
    }
    Ok(())
}

/// Preconditions of `closure` (review F1/F4): L_p has top λ and p is a type of (λ, λ); without
/// an idempotent L_p must be indecomposable of type p (char L_p = b_p), with one it must pass
/// `check_idempotent`. Err = Refused.
pub fn closure_precheck<F: Field>(ar: &Arena<F>, s: &mut Schur, lp: ExprId, p: Tid, ep: Option<MapId>) -> Result<(), String> {
    let lam = ar.lam().to_vec();
    if s.lam != lam {
        return Err(format!("Schur lambda {:?} != arena lambda {:?}", s.lam, lam));
    }
    let w = word(ar, lp);
    if top_of(&lam, &w) != lam || s.table(p).cols() != lam {
        return Err(format!("closure: L_p has top {:?} and p is a table of {:?}, expected {:?}", top_of(&lam, &w), s.table(p).cols(), lam));
    }
    match ep {
        None => {
            let ch = s.ladder_char(&w);
            let ex = s.expand(&ch);
            if ex.len() != 1 || ex[0].0 != p || ex[0].1 != crate::schur::Lp::one() {
                let parts: Vec<String> = ex.iter().map(|(t, c)| format!("{}{{{}}}", s.table(*t).fmt(), c.fmt_ec())).collect();
                return Err(format!("L_p decomposable: e_p required (char L_p = {} for p = {})", parts.join(" + "), s.table(p).fmt()));
            }
            Ok(())
        }
        Some(e) => check_idempotent(ar, s, lp, p, e, 0x5eed_e001).map_err(|m| format!("e_p rejected: {}", m)),
    }
}

/// closure check (§4.3 item 7): every map is R^λ -> L_p (ins) or L_p -> R^λ (outs) (review F3),
/// degree multiset = c_diag(b_p), identity-coset matrix of rank n at two points.
///
/// Without e_p the identity-coset blocks have exactly n = c_diag(b_p)(1) entries (n×n matrix).
/// With e_p (SPEC §4.3 item 7, review F4) every map must be e_p∘g (ins) or g∘e_p (outs) as a DAG
/// node, so it lies in Hom(R^λ, B_p) resp. Hom(B_p, R^λ) with B_p = Im e_p; its block has N >= n
/// entries in the n-dimensional image (resp. coimage) of the block of e_p, and rank n of the n×N
/// matrix is the same as nonsingularity in reduced coordinates of Im e_p. The determinant argument
/// then applies to the free rank-n module Hom(R^λ, B_p) (degree multiset c_diag(b_p)).
#[allow(clippy::too_many_arguments)]
pub fn certify_closure<F: Field>(ar: &Arena<F>, lam: &[u8], r0: ExprId, lp: ExprId, ep: Option<MapId>, ms: &[MapId], want_degs: &[i32], outs: bool, seed: u64) -> Result<(), String> {
    let (a, b) = if outs { (lp, r0) } else { (r0, lp) };
    for &m in ms {
        let i = ar.info(m);
        if (i.src, i.tgt) != (a, b) {
            return Err(format!("closure: map {} is {}->{}, expected {}->{} ({})", m, i.src, i.tgt, a, b, if outs { "outs" } else { "ins" }));
        }
        if let Some(e) = ep {
            let wrapped = match ar.node(m) {
                MapNode::Comp(g, h) => if outs { *h == e } else { *g == e },
                _ => false,
            };
            if !wrapped {
                return Err(format!("closure: map {} is not {} (not in the image of e_p)", m, if outs { "g∘e_p" } else { "e_p∘g" }));
            }
        }
    }
    let mut degs: Vec<i32> = ms.iter().map(|&m| ar.info(m).deg).collect();
    degs.sort_unstable();
    let mut w = want_degs.to_vec();
    w.sort_unstable();
    if degs != w {
        return Err(format!("closure: degrees {:?} != c_diag {:?}", degs, w));
    }
    let c = identity_coset(lam);
    for p in points::<F>(ar.k(), seed, 2, 0) {
        let mut ev = Eval::new(ar, p);
        let rows: Vec<Vec<F>> = ms.iter().map(|&m| ev.block(m, c).a.clone()).collect();
        let n = ms.len();
        let width = if ep.is_some() { rows.first().map_or(n, |r| r.len()) } else { n };
        if rows.iter().any(|r| r.len() != width) || width < n {
            return Err(format!("closure: identity-coset blocks are not of size {} ({})", width, if outs { "outs" } else { "ins" }));
        }
        if crate::dense::rank(&rows, width) != n {
            return Err("closure: identity-coset matrix singular".into());
        }
    }
    for &m in ms {
        assert_genuine(ar, m)?;
    }
    Ok(())
}

/// c_diag(b_p) as a sorted degree list
fn cdiag_degrees(s: &mut Schur, p: Tid) -> Vec<i32> {
    let cd = s.c_diag(p);
    let mut want = vec![];
    for (i, &c) in cd.c.iter().enumerate() {
        for _ in 0..c {
            want.push(cd.lo + i as i32);
        }
    }
    want
}

// ------------------------------------------------------------------ cost table (§4.5)

/// one measured exact call
#[derive(Clone, Debug)]
pub struct CostRec {
    pub labels: Vec<u8>,
    pub fs: u128,
    pub ft: u128,
    pub unknowns: u128,
    pub deg: i32,
    pub secs: f64,
    pub rss_mb: f64,
    /// killed by the watchdog: `secs` is a lower bound
    pub censored: bool,
}

/// Generation of the exact solver whose timings a cost table holds (`sbim::hom_space` with
/// `linsys::kernel_lazy`); bump when the solver's cost changes. Tables with another stamp are
/// rejected (review tests finding 3).
pub const COST_ALG: &str = "hs-kl1";

/// The machine class a cost table was measured on: CPU count and model name (not the host name,
/// which is random in cloud containers).
pub fn host_stamp() -> String {
    let n = std::thread::available_parallelism().map(|x| x.get()).unwrap_or(1);
    let model = std::fs::read_to_string("/proc/cpuinfo")
        .ok()
        .and_then(|t| t.lines().find(|l| l.starts_with("model name")).and_then(|l| l.split_once(':')).map(|(_, m)| m.trim().to_string()))
        .or_else(|| {
            std::process::Command::new("sysctl").args(["-n", "machdep.cpu.brand_string"]).output().ok().and_then(|o| String::from_utf8(o.stdout).ok()).map(|x| x.trim().to_string()).filter(|x| !x.is_empty())
        })
        .unwrap_or_else(|| "unknown".into());
    let model: String = model.chars().map(|c| if c.is_ascii_alphanumeric() || c == '.' || c == '@' { c } else { '_' }).collect();
    format!("{}cpu-{}", n, model)
}

/// Kind of an exact call for the cost buckets (M5c): 0 = one-sided (one end is the empty word,
/// fibre 1: the whole-ladder closure calls Hom(R, L_p), Hom(L_p, R)), 1 = two-sided (pair and
/// window calls). Their costs at equal fibre product differ by orders of magnitude, so one
/// censored closure call must not close the bucket of pair calls (M5b open problem 3).
pub fn call_kind(fs: u128, ft: u128) -> u8 {
    if fs <= 1 || ft <= 1 { 0 } else { 1 }
}

/// measured hom_space / rlambda_bases costs keyed by (sorted labels, call kind, fibre product
/// rounded up to a power of 2, degree) -> max seconds observed; built and saved (text) in M4b
#[derive(Clone, Debug)]
pub struct CostTable {
    pub recs: Vec<CostRec>,
    map: BTreeMap<(Vec<u8>, u8, u32, i32), f64>,
    /// solver generation and machine class of the measurements (header stamp)
    pub alg: String,
    pub host: String,
}
impl Default for CostTable {
    fn default() -> Self {
        CostTable { recs: vec![], map: BTreeMap::new(), alg: COST_ALG.to_string(), host: host_stamp() }
    }
}

fn bucket(fp: u128) -> u32 {
    let mut b = 0u32;
    while (1u128 << b) < fp.max(1) {
        b += 1;
    }
    b
}

impl CostTable {
    pub fn new() -> CostTable {
        CostTable::default()
    }
    fn key(labels: &[u8], kind: u8, fp: u128, deg: i32) -> (Vec<u8>, u8, u32, i32) {
        let mut l = labels.to_vec();
        l.sort_unstable();
        (l, kind, bucket(fp), deg)
    }
    pub fn record(&mut self, r: CostRec) {
        let k = CostTable::key(&r.labels, call_kind(r.fs, r.ft), r.fs * r.ft, r.deg);
        let e = self.map.entry(k).or_insert(0.0);
        // a censored record is only a lower bound (killed at its limit): the bucket prediction
        // becomes +infinity, so it can never admit a call (review tests finding 3: a bucket
        // whose maximum was a censored 60 s record admitted every window in it)
        *e = if r.censored { f64::INFINITY } else { e.max(r.secs) };
        self.recs.push(r);
    }
    pub fn len(&self) -> usize {
        self.recs.len()
    }
    pub fn is_empty(&self) -> bool {
        self.recs.is_empty()
    }
    /// records per sorted label class
    pub fn classes(&self) -> BTreeMap<Vec<u8>, usize> {
        let mut m = BTreeMap::new();
        for r in &self.recs {
            let mut l = r.labels.clone();
            l.sort_unstable();
            *m.entry(l).or_insert(0) += 1;
        }
        m
    }
    /// maximum over records of exactly this call shape (labels in order, fibres, unknowns,
    /// degree); used by Search's window admission before the bucket prediction
    pub fn predict_exact(&self, labels: &[u8], fs: u128, ft: u128, unknowns: u128, deg: i32) -> Option<f64> {
        // completed measurements give the time; a watchdog kill is a lower bound, which only
        // matters (infinite prediction) if no completed measurement reaches it
        let (mut done, mut bound): (Option<f64>, Option<f64>) = (None, None);
        for r in &self.recs {
            if r.labels == labels && r.fs == fs && r.ft == ft && r.unknowns == unknowns && r.deg == deg {
                let slot = if r.censored { &mut bound } else { &mut done };
                *slot = Some(slot.map_or(r.secs, |x: f64| x.max(r.secs)));
            }
        }
        match (done, bound) {
            (Some(t), Some(b)) if b > t => Some(f64::INFINITY),
            (Some(t), _) => Some(t),
            (None, Some(_)) => Some(f64::INFINITY),
            (None, None) => None,
        }
    }
    /// bucket maximum over both call kinds (the pre-M5c key; kept for reports and tests)
    pub fn predict(&self, labels: &[u8], fibre_product: u64, deg: i32) -> Option<f64> {
        let a = self.map.get(&CostTable::key(labels, 0, fibre_product as u128, deg)).cloned();
        let b = self.map.get(&CostTable::key(labels, 1, fibre_product as u128, deg)).cloned();
        match (a, b) {
            (Some(x), Some(y)) => Some(x.max(y)),
            (x, y) => x.or(y),
        }
    }
    /// bucket maximum of the calls of the same kind (`call_kind(fs, ft)`) as a call fs x ft
    pub fn predict_kind(&self, labels: &[u8], fs: u128, ft: u128, deg: i32) -> Option<f64> {
        self.map.get(&CostTable::key(labels, call_kind(fs, ft), fs.saturating_mul(ft), deg)).cloned()
    }
    /// prediction for one call: the exact-shape record if there is one, else the bucket of its kind
    pub fn predict_call(&self, labels: &[u8], fs: u128, ft: u128, unknowns: u128, deg: i32) -> Option<f64> {
        self.predict_exact(labels, fs, ft, unknowns, deg).or_else(|| self.predict_kind(labels, fs, ft, deg))
    }
    pub fn save(&self, path: &std::path::Path) -> std::io::Result<()> {
        use std::io::Write;
        let tmp = path.with_extension(format!("{}.tmp", std::process::id()));
        {
            let mut w = std::io::BufWriter::new(std::fs::File::create(&tmp)?);
            writeln!(w, "krhom-costs v2 alg={} host={} fields=labels,fs,ft,unknowns,deg,secs,rss_mb,censored", self.alg, self.host)?;
            for r in &self.recs {
                let l: Vec<String> = r.labels.iter().map(|x| x.to_string()).collect();
                writeln!(w, "{} {} {} {} {} {:.6} {:.1} {}", l.join(","), r.fs, r.ft, r.unknowns, r.deg, r.secs, r.rss_mb, if r.censored { 1 } else { 0 })?;
            }
            w.flush()?;
        }
        std::fs::rename(&tmp, path)
    }
    pub fn load(path: &std::path::Path) -> std::io::Result<CostTable> {
        let text = std::fs::read_to_string(path)?;
        let bad = |l: &str| std::io::Error::new(std::io::ErrorKind::InvalidData, format!("cost table: bad line {:?}", l));
        let mut lines = text.lines();
        let h = lines.next().unwrap_or("");
        let mismatch = |m: String| std::io::Error::new(std::io::ErrorKind::InvalidData, m);
        if h.starts_with("krhom-costs v1 ") {
            return Err(mismatch(format!("{}: unstamped v1 cost table (no solver / host stamp; recalibrate with `g1 calib`)", path.display())));
        }
        if !h.starts_with("krhom-costs v2 ") {
            return Err(bad(h));
        }
        let field = |k: &str| h.split_whitespace().find_map(|f| f.strip_prefix(k)).unwrap_or("").to_string();
        let (alg, host) = (field("alg="), field("host="));
        if alg != COST_ALG {
            return Err(mismatch(format!("{}: cost table solver stamp {:?} != {:?} (recalibrate with `g1 calib`)", path.display(), alg, COST_ALG)));
        }
        let here = host_stamp();
        if host != here {
            if std::env::var("KRHOM_COSTS_ANY_HOST").is_err() {
                return Err(mismatch(format!("{}: cost table host stamp {:?} != this machine {:?} (recalibrate, or set KRHOM_COSTS_ANY_HOST=1)", path.display(), host, here)));
            }
            eprintln!("[costs] warning: {} was measured on {:?}, this machine is {:?} (KRHOM_COSTS_ANY_HOST)", path.display(), host, here);
        }
        let mut t = CostTable { alg, host, ..CostTable::new() };
        for l in lines {
            if l.trim().is_empty() {
                continue;
            }
            let f: Vec<&str> = l.split_whitespace().collect();
            if f.len() != 8 {
                return Err(bad(l));
            }
            let labels: Vec<u8> = f[0].split(',').map(|x| x.parse::<u8>().map_err(|_| bad(l))).collect::<Result<_, _>>()?;
            let p = |x: &str| x.parse::<u128>().map_err(|_| bad(l));
            t.record(CostRec {
                labels,
                fs: p(f[1])?,
                ft: p(f[2])?,
                unknowns: p(f[3])?,
                deg: f[4].parse().map_err(|_| bad(l))?,
                secs: f[5].parse().map_err(|_| bad(l))?,
                rss_mb: f[6].parse().map_err(|_| bad(l))?,
                censored: f[7] == "1",
            });
        }
        Ok(t)
    }
}

// ------------------------------------------------------------------ exact window calls (shared)

/// one exact window call through the template library, with its cost in the stats
/// counts computed (not cached) calls into (calls, secs)
fn exact_call<F: Field>(ar: &mut Arena<F>, labels: &[u8], a: &[Rung], b: &[Rung], d: i32, cnt: (&mut usize, &mut f64)) -> Result<Vec<TmplId>, String> {
    let before = ar.lib_ref().stats().computed;
    let t0 = Instant::now();
    let trace = std::env::var("KRHOM_TRACE").is_ok() && !ar.lib_ref().has_hom(&Window::new(labels, a), &Window::new(labels, b), d);
    if trace {
        eprintln!("[exact] start {:?} [{}] -> [{}] degree {} (fibres {} x {}, unknowns {})", labels, wstr(a), wstr(b), d, crate::schur::flows(labels, a), crate::schur::flows(labels, b), hom_unknowns(labels, a, b, d));
    }
    let ids = ar.lib().try_hom_basis_cap(&Window::new(labels, a), &Window::new(labels, b), d, u128::MAX)?;
    if ar.lib_ref().stats().computed > before {
        *cnt.0 += 1;
        *cnt.1 += elapsed(&t0);
    }
    if trace {
        eprintln!("[exact] done {:.2}s dim {}", elapsed(&t0), ids.len());
    }
    // KRHOM_WINLOG=PATH: append every exact call (computed or cached) in the `g1 calib --windows`
    // format (labels, source word, target word, degree), for recalibrating the cost table
    if let Ok(path) = std::env::var("KRHOM_WINLOG") {
        use std::io::Write as _;
        let code = |rs: &[Rung]| -> String {
            if rs.is_empty() {
                return "-".into();
            }
            rs.iter().map(|r| format!("{}{}.{}", if r.dir == crate::web::Dir::E { "E" } else { "F" }, r.i, r.r)).collect::<Vec<_>>().join(",")
        };
        let lab: Vec<String> = labels.iter().map(|x| x.to_string()).collect();
        if let Ok(mut f) = std::fs::OpenOptions::new().create(true).append(true).open(&path) {
            let _ = writeln!(f, "{} {} {} {}", lab.join(","), code(a), code(b), d);
        }
    }
    Ok(ids)
}

// ------------------------------------------------------------------ the Exact backend (§4.5)

/// defaults: 1 600, 40 000, max_label 2, costs None, max_call_s 10.0 (§4.5; the unknown cap was
/// 5 000 for the old hom_space solver; with `linsys::kernel_lazy` every Hom^d (d = 0..4) between
/// canonical ladders of (2,1,2) and (2,2,2) with fibre product <= 1 600 takes <= 7.4 s on the
/// cloud box, at most 39 165 unknowns)
#[derive(Clone, Debug)]
pub struct Exact {
    pub max_fibre_product: u64,
    pub max_unknowns: u64,
    pub max_label: u8,
    pub costs: Option<CostTable>,
    pub max_call_s: f64,
    /// counters of the last call if it failed (`failed_stats`)
    pub failed: Option<SearchStats>,
}
impl Default for Exact {
    fn default() -> Self {
        Exact { max_fibre_product: 1600, max_unknowns: 40_000, max_label: 2, costs: None, max_call_s: 10.0, failed: None }
    }
}

impl Exact {
    /// Admission (§4.5): (max label <= max_label, fibre product <= cap, unknowns <= cap) or a
    /// cost-table prediction <= max_call_s; a prediction above max_call_s vetoes the label rule
    /// (per-call ceiling). The prediction is the record of exactly this call shape if measured,
    /// else the bucket of the same call kind (M5c: a censored closure call no longer vetoes pair
    /// calls of the same fibre product). Ok(rule) or Err(reason).
    pub fn admit(&self, labels: &[u8], fs: u128, ft: u128, unknowns: u128, deg: i32) -> Result<&'static str, String> {
        let fp = fs.saturating_mul(ft);
        let pred = self.costs.as_ref().and_then(|c| c.predict_call(labels, fs, ft, unknowns, deg));
        let maxl = labels.iter().copied().max().unwrap_or(0);
        let label_rule = maxl <= self.max_label && fp <= self.max_fibre_product as u128 && unknowns <= self.max_unknowns as u128;
        match pred {
            Some(t) if t <= self.max_call_s => Ok("cost-table"),
            Some(t) if label_rule => Err(format!("label rule vetoed: cost table predicts {:.1}s > {:.1}s", t, self.max_call_s)),
            _ if label_rule => Ok("label-rule"),
            Some(t) => Err(format!("cost table predicts {:.1}s > {:.1}s (labels {:?}, fibres {}x{}, degree {})", t, self.max_call_s, labels, fs, ft, deg)),
            None => Err(format!("not admitted: labels {:?} (max {} > {} or caps: fibre product {} / {}, unknowns {} / {}) and no cost-table entry", labels, maxl, self.max_label, fp, self.max_fibre_product, unknowns, self.max_unknowns)),
        }
    }
    /// exact basis of Hom^d(src, tgt) on the whole ladders, as whole-ladder Local nodes
    fn basis<F: Field>(&self, ar: &mut Arena<F>, src: ExprId, tgt: ExprId, d: i32, st: &mut SearchStats) -> Result<Vec<MapId>, (MissKind, String)> {
        let lam = ar.lam().to_vec();
        let (ws, wt) = (word(ar, src), word(ar, tgt));
        let (fs, ft) = (crate::schur::flows(&lam, &ws), crate::schur::flows(&lam, &wt));
        let unk = hom_unknowns(&lam, &ws, &wt, d);
        self.admit(&lam, fs, ft, unk, d).map_err(|e| (MissKind::Refused, e))?;
        let ids = exact_call(ar, &lam, &ws, &wt, d, (&mut st.exact_calls, &mut st.exact_secs)).map_err(|e| (MissKind::Refused, e))?;
        let mut out = vec![];
        for t in ids {
            out.push(ar.add(MapNode::Local { src, tgt, at: 0, first: 0, tmpl: t }, "exact").map_err(|e| (MissKind::Exhausted, e))?);
        }
        Ok(out)
    }
}

impl<F: Field> HomProvider<F> for Exact {
    fn name(&self) -> &'static str {
        "Exact"
    }
    fn full(&mut self, ar: &mut Arena<F>, s: &mut Schur, src: ExprId, tgt: ExprId, d: i32, _b: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        let t0 = Instant::now();
        let mark = ar.mark();
        let mut st = SearchStats::default();
        self.failed = None;
        let ms = match self.basis(ar, src, tgt, d, &mut st) {
            Ok(v) => v,
            Err((k, e)) => {
                ar.commit(mark, &[]);
                st.secs = elapsed(&t0);
                self.failed = Some(st);
                return Err(miss(&format!("Exact full: {}", e), d, None, 0, 0, k));
            }
        };
        match certify_full_plan(ar, s, src, tgt, d, &ms, 0x5eed_0001) {
            Ok(plan) => {
                st.detect = plan.kind;
                st.points = plan.npoints;
                st.candidates = ms.len();
                st.evaluated = ms.len();
                let keep = ar.commit(mark, &ms);
                st.secs = elapsed(&t0);
                Ok((keep, st))
            }
            Err(e) => {
                ar.commit(mark, &[]);
                st.secs = elapsed(&t0);
                self.failed = Some(st);
                Err(miss(&format!("Exact full: {}", e), d, None, ms.len(), 0, MissKind::Exhausted))
            }
        }
    }
    fn pairs(&mut self, ar: &mut Arena<F>, s: &mut Schur, x: ExprId, q: Tid, lq: ExprId, eq: Option<MapId>, sdeg: i32, mult: usize, _b: &Budget) -> Result<(Pairs<F>, SearchStats), Vec<Missing>> {
        let t0 = Instant::now();
        let mark = ar.mark();
        let mut st = SearchStats::default();
        self.failed = None;
        let strat = Some(s.table(q).fmt());
        let lam = ar.lam().to_vec();
        let r = (|| -> Result<Pairs<F>, (MissKind, String)> {
            pairs_precheck(ar, s, x, q, lq, eq, sdeg, mult).map_err(|e| (MissKind::Refused, e))?;
            let mut io = self.basis(ar, lq, x, sdeg, &mut st)?;
            let mut pi = self.basis(ar, x, lq, -sdeg, &mut st)?;
            st.candidates = io.len() + pi.len();
            st.evaluated = st.candidates;
            with_idem(ar, eq, &mut io, &mut pi).map_err(|e| (MissKind::Exhausted, e))?;
            certify_pairs(ar, x, lq, sdeg, &io, &pi, rep_rank(s, &lam, q), mult, 0x5eed_0002).map_err(|e| (MissKind::Exhausted, e))
        })();
        match r {
            Ok(p) => {
                let mut keep = p.iotas.clone();
                keep.extend(&p.pis);
                let k = ar.commit(mark, &keep);
                st.secs = elapsed(&t0);
                st.points = 2;
                Ok((Pairs { iotas: k[..mult].to_vec(), pis: k[mult..].to_vec(), form: p.form }, st))
            }
            Err((kind, e)) => {
                ar.commit(mark, &[]);
                st.secs = elapsed(&t0);
                self.failed = Some(st);
                Err(miss(&format!("Exact pairs: {}", e), sdeg, strat, 0, mult, kind))
            }
        }
    }
    fn closure(&mut self, ar: &mut Arena<F>, s: &mut Schur, lp: ExprId, p: Tid, ep: Option<MapId>, outs: bool, _b: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        let t0 = Instant::now();
        let mark = ar.mark();
        let mut st = SearchStats::default();
        self.failed = None;
        let want = cdiag_degrees(s, p);
        let lam = ar.lam().to_vec();
        let r = (|| -> Result<Vec<MapId>, (MissKind, String)> {
            closure_precheck(ar, s, lp, p, ep).map_err(|e| (MissKind::Refused, e))?;
            let r0 = eid(ar, &[]).map_err(|e| (MissKind::Exhausted, e))?;
            let mut sel = Selector::new(identity_coset(&lam), ar.k(), want.len(), 0x5eed_0003, ep.is_some());
            let mut degs: Vec<i32> = want.clone();
            degs.dedup();
            for d in degs {
                let need = want.iter().filter(|&&x| x == d).count();
                let (a, b) = if outs { (lp, r0) } else { (r0, lp) };
                let mut cands = self.basis(ar, a, b, d, &mut st)?;
                if let Some(e) = ep {
                    for m in cands.iter_mut() {
                        *m = ar.add(if outs { MapNode::Comp(*m, e) } else { MapNode::Comp(e, *m) }, "e_p").map_err(|e| (MissKind::Exhausted, e))?;
                    }
                }
                let got = sel.offer_many(ar, &cands, need);
                st.candidates += cands.len();
                st.evaluated += cands.len();
                if got < need {
                    return Err((MissKind::Exhausted, format!("closure: degree {}: {} of {} generators", d, got, need)));
                }
            }
            let ms = sel.chosen.clone();
            certify_closure(ar, &lam, r0, lp, ep, &ms, &want, outs, 0x5eed_0004).map_err(|e| (MissKind::Exhausted, e))?;
            Ok(ms)
        })();
        match r {
            Ok(ms) => {
                let k = ar.commit(mark, &ms);
                st.secs = elapsed(&t0);
                Ok((k, st))
            }
            Err((kind, e)) => {
                ar.commit(mark, &[]);
                st.secs = elapsed(&t0);
                self.failed = Some(st);
                Err(miss(&format!("Exact closure: {}", e), 0, Some(s.table(p).fmt()), 0, want.len(), kind))
            }
        }
    }
    fn end0(&mut self, ar: &mut Arena<F>, s: &mut Schur, l: ExprId, b: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        self.full(ar, s, l, l, 0, b)
    }
    fn failed_stats(&self) -> Option<SearchStats> {
        self.failed.clone()
    }
}

/// greedy degree-by-degree determinant selection on the identity coset (closure)
struct Selector<F: Field> {
    coset: u64,
    pts: Vec<Point<F>>,
    rows: Vec<Vec<F>>,
    chosen: Vec<MapId>,
    n: usize,
    /// block length: n without e_p; with e_p the flows of L_p at the identity coset (set by the
    /// first candidate)
    width: Option<usize>,
}
impl<F: Field> Selector<F> {
    fn new(coset: u64, k: usize, n: usize, seed: u64, with_ep: bool) -> Self {
        Selector { coset, pts: points::<F>(k, seed, 1, 0), rows: vec![], chosen: vec![], n, width: if with_ep { None } else { Some(n) } }
    }
    /// offer candidates; keep those that raise the rank, up to `need` of them; returns kept count
    fn offer_many(&mut self, ar: &Arena<F>, cands: &[MapId], need: usize) -> usize {
        let mut got = 0;
        let mut ev = Eval::new(ar, self.pts[0].clone());
        for &m in cands {
            if got >= need {
                break;
            }
            if self.chosen.len() >= self.n {
                break;
            }
            let v = ev.block(m, self.coset).a.clone();
            let w = *self.width.get_or_insert(v.len());
            if v.len() != w || w < self.n {
                continue;
            }
            let mut rows = self.rows.clone();
            rows.push(v.clone());
            if crate::dense::rank(&rows, w) == rows.len() {
                self.rows.push(v);
                self.chosen.push(m);
                got += 1;
            }
        }
        got
    }
}

// ------------------------------------------------------------------ the Search provider (§4.4)

/// windows of a word: contiguous rungs at..at+len (2 <= len <= lmax) spanning exactly w adjacent
/// uprights, 2 <= w <= wmax
pub fn windows_of(w: &[Rung], wmax: usize, lmax: usize) -> Vec<WinPos> {
    let mut out = vec![];
    for ww in 2..=wmax {
        for at in 0..w.len() {
            for len in 2..=lmax {
                if at + len > w.len() {
                    break;
                }
                let lo = w[at..at + len].iter().map(|r| r.i).min().unwrap();
                let hi = w[at..at + len].iter().map(|r| r.i + 1).max().unwrap();
                if hi - lo + 1 != ww {
                    continue;
                }
                out.push(WinPos { at, len, first: lo, w: ww });
            }
        }
    }
    out
}

fn shift_word(rs: &[Rung], by: isize) -> Vec<Rung> {
    rs.iter().map(|r| Rung { i: (r.i as isize + by) as usize, dir: r.dir, r: r.r }).collect()
}

/// distinct segments of a word: (level, upright) where an upright with a nonzero label starts a
/// run untouched by rungs (level 0, or the rung below touches the upright); `interior` keeps
/// segments that start above the bottom and end below the top
pub fn segments(lam: &[u8], w: &[Rung], interior: bool) -> Vec<(usize, usize, u8)> {
    let n = lam.len();
    let lv = Expr::new(lam, w).expect("segments: invalid word").levels();
    let touches = |r: &Rung, u: usize| r.i == u || r.i + 1 == u;
    let mut out = vec![];
    for (level, labels) in lv.iter().enumerate() {
        for u in 0..n {
            let starts = level == 0 || touches(&w[level - 1], u);
            if !starts || labels[u] == 0 {
                continue;
            }
            if interior {
                let ends_below_top = w[level..].iter().any(|r| touches(r, u));
                if level == 0 || !ends_below_top {
                    continue;
                }
            }
            out.push((level, u, labels[u]));
        }
    }
    out
}

/// partitions of k with at most `rows` parts (descending order of parts, lexicographic)
fn partitions(k: usize, rows: usize) -> Vec<Vec<u8>> {
    fn rec(k: usize, maxp: usize, rows: usize, cur: &mut Vec<u8>, out: &mut Vec<Vec<u8>>) {
        if k == 0 {
            out.push(cur.clone());
            return;
        }
        if cur.len() == rows {
            return;
        }
        for p in (1..=maxp.min(k)).rev() {
            cur.push(p as u8);
            rec(k - p, p, rows, cur, out);
            cur.pop();
        }
    }
    let mut out = vec![];
    rec(k, k, rows, &mut vec![], &mut out);
    out
}

/// dot monomials of total degree 2k on the given segments: products of at most `maxf` Schur
/// dots (segment index, partition), ordered by number of factors, then segment, then partition
fn dot_monomials(segs: &[(usize, usize, u8)], k: usize, maxf: usize, cap: usize) -> Vec<Vec<(usize, Vec<u8>)>> {
    let mut out: Vec<Vec<(usize, Vec<u8>)>> = vec![];
    if k == 0 {
        return vec![vec![]];
    }
    fn rec(segs: &[(usize, usize, u8)], k: usize, from: usize, left: usize, cur: &mut Vec<(usize, Vec<u8>)>, out: &mut Vec<Vec<(usize, Vec<u8>)>>, cap: usize) {
        if out.len() >= cap {
            return;
        }
        if k == 0 {
            out.push(cur.clone());
            return;
        }
        if left == 0 {
            return;
        }
        for si in from..segs.len() {
            let rows = segs[si].2 as usize;
            for part in 1..=k {
                for nu in partitions(part, rows) {
                    cur.push((si, nu));
                    rec(segs, k - part, si + 1, left - 1, cur, out, cap);
                    cur.pop();
                    if out.len() >= cap {
                        return;
                    }
                }
            }
        }
    }
    for f in 1..=maxf {
        let mut layer = vec![];
        rec(segs, k, 0, f, &mut vec![], &mut layer, cap);
        for m in layer {
            if m.len() == f && out.len() < cap {
                out.push(m);
            }
        }
    }
    out
}

/// an isomorphism move of a word (iso_reach)
#[derive(Clone, Debug)]
enum IsoMove {
    Swap(usize),
    Rex(WinPos, Vec<Rung>),
}

/// one local summand of a window: canonical word of the window block, shift, ι (L_t -> W) and
/// π (W -> L_t) templates with π_k ∘ ι_l = δ_kl
#[derive(Clone, Debug)]
struct Piece {
    word: Vec<Rung>,
    s: i32,
    iota: TmplId,
    pi: TmplId,
}

#[derive(Clone, Debug)]
struct Leaf {
    w: Vec<Rung>,
    tword: Vec<Rung>,
    io: MapId,
    pi: MapId,
    shift: i32,
}

/// stop reason inside a call
type Stop = (MissKind, String);

struct Ctx {
    t0: Instant,
    b: Budget,
    st: SearchStats,
    stuck: usize,
    /// stuck branches whose only cause was window admission (subset of `stuck`)
    stuck_adm: usize,
    /// tiers or degrees skipped because admission refused their exact window
    tier_refused: usize,
}
impl Ctx {
    /// the kind of a call whose candidate stream ended: Refused when every stuck branch and every
    /// skipped tier was caused by admission (and there was at least one), Exhausted otherwise
    /// (review F5: Exhausted is the completeness signal, admission is not)
    fn end_kind(&self) -> MissKind {
        if self.stuck_adm + self.tier_refused > 0 && self.stuck == self.stuck_adm {
            MissKind::Refused
        } else {
            MissKind::Exhausted
        }
    }
    fn end_note(&self) -> String {
        format!("{} stuck branches, {} by admission, {} tiers refused", self.stuck, self.stuck_adm, self.tier_refused)
    }
    fn stuck_by(&mut self, k: MissKind) {
        self.stuck += 1;
        if k == MissKind::Refused {
            self.stuck_adm += 1;
        }
    }
    fn check<F: Field>(&mut self, ar: &Arena<F>) -> Result<(), Stop> {
        let by = ar.bytes();
        self.st.peak_arena_bytes = self.st.peak_arena_bytes.max(by.total() as u64);
        self.st.peak_index_bytes = self.st.peak_index_bytes.max(by.indices as u64);
        if elapsed(&self.t0) > self.b.t_max_s {
            return Err((MissKind::Time, format!("time budget {:.0}s exceeded", self.b.t_max_s)));
        }
        if by.total() as u64 > self.b.m_max_bytes {
            return Err((MissKind::Memory, format!("arena {} bytes > {}", by.total(), self.b.m_max_bytes)));
        }
        if self.st.candidates > self.b.c_max {
            return Err((MissKind::Exhausted, format!("candidate budget {} exhausted", self.b.c_max)));
        }
        Ok(())
    }
}

/// Search provider: tiers S0..S6 (a sub-range can be selected), window templates in the TmplLib,
/// per-window decomposition and commutation caches (template ids: one Search per Arena)
pub struct Search {
    lo: usize,
    hi: usize,
    /// measured costs for window admission (otherwise the unknown-count cap)
    pub costs: Option<CostTable>,
    pub max_window_unknowns: u128,
    /// fibre-product cap of Search windows without a measured record (6 400; was
    /// tmpl::HOM_WINDOW_CAP = 1 600 for the old hom_space solver)
    pub max_window_fibre: u128,
    pub max_window_s: f64,
    /// unknown cap of commutation windows (two rungs, degree 0)
    pub max_comm_unknowns: u128,
    /// unknown cap of degree-0 windows between isomorphic words (rex moves, leaf isomorphisms)
    pub max_iso_unknowns: u128,
    /// maximal window width (uprights) and length (rungs) of a local decomposition
    pub win_w: usize,
    pub win_len: usize,
    sc: SchurCache,
    splits: HashMap<(Vec<u8>, Vec<Rung>), Result<Vec<Piece>, String>>,
    comms: HashMap<(Vec<u8>, Vec<Rung>), Option<(TmplId, TmplId)>>,
    /// `TmplLib::instance` of the library the TmplIds in `splits` / `comms` belong to; the caches
    /// are cleared when a call comes with another library (review F3)
    lib_instance: u64,
    /// the last `choose_window` refused at least one decomposable window by admission
    last_win_refused: bool,
    /// counters of the last call if it failed (`failed_stats`)
    failed: Option<SearchStats>,
    /// one line per decision of the last call (tiers tried, outcomes)
    pub trace: Vec<String>,
    /// exact window calls computed by Search (not cached) and their seconds (cumulative)
    pub window_calls: usize,
    pub window_secs: f64,
}

impl Search {
    fn note(&mut self, m: String) {
        if std::env::var("KRHOM_TRACE").is_ok() {
            eprintln!("[search S{}-S{}] {}", self.lo, self.hi, m);
        }
        self.trace.push(m);
    }
    pub fn tiers(lo: usize, hi: usize) -> Search {
        assert!(lo <= hi && hi <= 6, "Search::tiers: 0 <= lo <= hi <= 6");
        Search {
            lo,
            hi,
            costs: None,
            // measured with linsys::kernel_lazy (cloud box): (2,4,0) [F1(1) F0(1) E0(1)] <->
            // [F1(2) E0(2) F0(2) E1(1)] degree 0, 32 x 72, 8 762 unknowns: 0.77 / 0.97 s; whole
            // (2,2,2) ladders with 39 165 unknowns: 7.4 s; (1,5,0) [F1(2) E0(3) F0(2) E1(1)] ->
            // [E0(4) F1(1) F0(3)] degree 1, 120 x 50, 33 134 unknowns: 5.6 s; (5,1,0) 60 x 90,
            // 21 300 unknowns: 1.8 s (the old caps were 2 500 / 1 600)
            max_window_unknowns: 40_000,
            max_window_fibre: 6400,
            max_window_s: 120.0,
            max_comm_unknowns: 20000,
            max_iso_unknowns: 6000,
            win_w: 3,
            win_len: 4,
            sc: SchurCache::default(),
            splits: HashMap::new(),
            comms: HashMap::new(),
            lib_instance: 0,
            last_win_refused: false,
            failed: None,
            trace: vec![],
            window_calls: 0,
            window_secs: 0.0,
        }
    }
    fn on(&self, k: usize) -> bool {
        self.lo <= k && k <= self.hi
    }
    /// bind the TmplId caches to the Arena's template library (cleared on a change)
    fn bind<F: Field>(&mut self, ar: &Arena<F>) {
        let id = ar.lib_ref().instance();
        if id != self.lib_instance {
            self.splits.clear();
            self.comms.clear();
            self.lib_instance = id;
        }
    }
    /// window admission as Search decides it (cache hit, measured shape, caps, bucket) (tests)
    pub fn window_admitted<F: Field>(&self, ar: &Arena<F>, labels: &[u8], a: &[Rung], b: &[Rung], d: i32) -> bool {
        self.win_ok(ar, labels, a, b, d)
    }
    /// number of cached window decompositions and commutations (tests)
    pub fn cached_windows(&self) -> (usize, usize) {
        (self.splits.len(), self.comms.len())
    }
    /// predicted admission of an exact window call (cache hits are always admitted)
    fn win_ok<F: Field>(&self, ar: &Arena<F>, labels: &[u8], a: &[Rung], b: &[Rung], d: i32) -> bool {
        if ar.lib_ref().has_hom(&Window::new(labels, a), &Window::new(labels, b), d) {
            return true;
        }
        let (fa, fb) = (crate::schur::flows(labels, a), crate::schur::flows(labels, b));
        let fp = fa * fb;
        let unk = hom_unknowns(labels, a, b, d);
        // a measured call of exactly this shape decides (also above the 1 600 fibre cap of
        // `hom_basis`: "M4b decides")
        if let Some(t) = self.costs.as_ref().and_then(|c| c.predict_exact(labels, fa, fb, unk, d)) {
            return t <= self.max_window_s;
        }
        if fp > self.max_window_fibre {
            return false;
        }
        if let Some(c) = self.costs.as_ref() {
            // the bucket maximum (a different call shape may dominate it) can admit, not veto
            if let Some(t) = c.predict_kind(labels, fa, fb, d) {
                if t <= self.max_window_s {
                    return true;
                }
            }
        }
        unk <= self.max_window_unknowns
    }
    /// admission of a degree-0 window between isomorphic words (rex moves, leaf isomorphisms):
    /// the window rule, or the larger isomorphism cap when no measured record refuses it
    fn iso_ok<F: Field>(&self, ar: &Arena<F>, labels: &[u8], a: &[Rung], b: &[Rung]) -> bool {
        if self.win_ok(ar, labels, a, b, 0) {
            return true;
        }
        let (fa, fb) = (crate::schur::flows(labels, a), crate::schur::flows(labels, b));
        let unk = hom_unknowns(labels, a, b, 0);
        if self.costs.as_ref().and_then(|c| c.predict_exact(labels, fa, fb, unk, 0)).is_some() {
            return false; // measured and too slow
        }
        unk <= self.max_iso_unknowns
    }
    fn win_cost<F: Field>(&self, ar: &Arena<F>, labels: &[u8], a: &[Rung], b: &[Rung], d: i32) -> u128 {
        if ar.lib_ref().has_hom(&Window::new(labels, a), &Window::new(labels, b), d) {
            0
        } else {
            hom_unknowns(labels, a, b, d)
        }
    }
    /// one exact window call, counted in the call's SearchStats (window_calls / window_secs) and
    /// in the cumulative window counters
    fn xcall<F: Field>(&mut self, ar: &mut Arena<F>, labels: &[u8], a: &[Rung], b: &[Rung], d: i32, st: &mut SearchStats) -> Result<Vec<TmplId>, String> {
        let (mut n, mut t) = (0usize, 0f64);
        let r = exact_call(ar, labels, a, b, d, (&mut n, &mut t));
        self.window_calls += n;
        self.window_secs += t;
        st.window_calls += n;
        st.window_secs += t;
        r
    }
    fn exact<F: Field>(&mut self, ar: &mut Arena<F>, labels: &[u8], a: &[Rung], b: &[Rung], d: i32, st: &mut SearchStats) -> Result<Vec<TmplId>, Stop> {
        if !self.win_ok(ar, labels, a, b, d) {
            return Err((MissKind::Refused, format!("window {:?} {:?} -> {:?} degree {} not admitted", labels, a, b, d)));
        }
        self.xcall(ar, labels, a, b, d, st).map_err(|e| (MissKind::Refused, e))
    }

    // ---------------------------------------------------------- commutation isomorphisms

    fn comm_tmpls<F: Field>(&mut self, ar: &mut Arena<F>, labels: &[u8], rel: &[Rung], st: &mut SearchStats) -> Result<(TmplId, TmplId, F), Stop> {
        let key = (labels.to_vec(), rel.to_vec());
        let sw = vec![rel[1], rel[0]];
        let ids = match self.comms.get(&key) {
            Some(Some(x)) => *x,
            Some(None) => return Err((MissKind::Exhausted, format!("commutation {:?} {:?}: no unique degree-0 map", labels, rel))),
            None => {
                // commutation windows (two rungs, degree 0, one-dimensional) are cheap whatever
                // their unknown count (measured: 2 936 unknowns in 0.09 s); own cap
                let ok = |me: &Search, a: &[Rung], b: &[Rung]| me.win_ok(ar, labels, a, b, 0) || hom_unknowns(labels, a, b, 0) <= me.max_comm_unknowns;
                if !ok(self, rel, &sw) || !ok(self, &sw, rel) {
                    return Err((MissKind::Refused, format!("commutation window {:?} [{}] not admitted", labels, wstr(rel))));
                }
                let mut f = vec![];
                let mut b = vec![];
                for (x, y, out) in [(rel, &sw[..], &mut f), (&sw[..], rel, &mut b)] {
                    *out = self.xcall(ar, labels, x, y, 0, st).map_err(|e| (MissKind::Refused, e))?;
                }
                let v = if f.len() == 1 && b.len() == 1 { Some((f[0], b[0])) } else { None };
                self.comms.insert(key.clone(), v);
                match v {
                    Some(x) => x,
                    None => return Err((MissKind::Exhausted, format!("commutation {:?} {:?}: Hom^0 dims {} / {}", labels, rel, f.len(), b.len()))),
                }
            }
        };
        let lib = ar.lib_ref();
        let prod = lib.get(ids.1).phi.matmul(&lib.get(ids.0).phi);
        match pmat_scalar(&prod) {
            Some(c) if !c.is_zero() => Ok((ids.0, ids.1, c)),
            _ => {
                self.comms.insert(key, None);
                Err((MissKind::Exhausted, format!("commutation {:?} {:?}: not an isomorphism", labels, rel)))
            }
        }
    }
    /// (w -> swapped, swapped -> w) exact inverses, for swapping rungs i, i+1 of w
    fn swap_maps<F: Field>(&mut self, ar: &mut Arena<F>, w: &[Rung], i: usize, st: &mut SearchStats) -> Result<(MapId, MapId, Vec<Rung>), Stop> {
        let lam = ar.lam().to_vec();
        let (r1, r2) = (w[i], w[i + 1]);
        let lo = r1.i.min(r2.i);
        let hi = r1.i.max(r2.i) + 1;
        let labels = level_of(&lam, w, i)[lo..=hi].to_vec();
        let rel = shift_word(&[r1, r2], -(lo as isize));
        let (tf, tb, c) = self.comm_tmpls(ar, &labels, &rel, st)?;
        let mut sw = w.to_vec();
        sw.swap(i, i + 1);
        let e = |ar: &mut Arena<F>, x: &[Rung]| eid(ar, x).map_err(|e| (MissKind::Exhausted, e));
        let (a, b) = (e(ar, w)?, e(ar, &sw)?);
        let f = ar.add(MapNode::Local { src: a, tgt: b, at: i, first: lo, tmpl: tf }, "commute").map_err(|e| (MissKind::Exhausted, e))?;
        let g = ar.add(MapNode::Local { src: b, tgt: a, at: i, first: lo, tmpl: tb }, "commute back").map_err(|e| (MissKind::Exhausted, e))?;
        let g = ar.add(MapNode::Lin(vec![(c.inv(), g)]), "commute inverse").map_err(|e| (MissKind::Exhausted, e))?;
        Ok((f, g, sw))
    }
    /// (from -> to, to -> from), mutually inverse, for commutation-equivalent words
    fn iso_chain<F: Field>(&mut self, ar: &mut Arena<F>, from: &[Rung], to: &[Rung], st: &mut SearchStats) -> Result<(MapId, MapId), Stop> {
        let ef = eid(ar, from).map_err(|e| (MissKind::Exhausted, e))?;
        if from == to {
            let id = ar.add(MapNode::Id(ef), "id").map_err(|e| (MissKind::Exhausted, e))?;
            return Ok((id, id));
        }
        // breadth-first search over swaps
        let mut prev: HashMap<Vec<Rung>, (Vec<Rung>, usize)> = HashMap::new();
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(from.to_vec());
        prev.insert(from.to_vec(), (vec![], usize::MAX));
        let mut found = false;
        while let Some(cur) = queue.pop_front() {
            if cur == to {
                found = true;
                break;
            }
            if prev.len() > 5000 {
                break;
            }
            for i in 0..cur.len().saturating_sub(1) {
                if rungs_commute(&cur[i], &cur[i + 1]) {
                    let mut nw = cur.clone();
                    nw.swap(i, i + 1);
                    if !prev.contains_key(&nw) {
                        prev.insert(nw.clone(), (cur.clone(), i));
                        queue.push_back(nw);
                    }
                }
            }
        }
        if !found {
            return Err((MissKind::Exhausted, "iso_chain: words are not commutation-equivalent".into()));
        }
        let mut steps = vec![];
        let mut cur = to.to_vec();
        while cur != from {
            let (p, i) = prev[&cur].clone();
            steps.push((p.clone(), i));
            cur = p;
        }
        steps.reverse();
        let mut fwd: Option<MapId> = None;
        let mut bwd: Option<MapId> = None;
        for (w, i) in steps {
            let (f, g, _) = self.swap_maps(ar, &w, i, st)?;
            fwd = Some(match fwd {
                None => f,
                Some(x) => ar.add(MapNode::Comp(f, x), "chain").map_err(|e| (MissKind::Exhausted, e))?,
            });
            bwd = Some(match bwd {
                None => g,
                Some(x) => ar.add(MapNode::Comp(x, g), "chain back").map_err(|e| (MissKind::Exhausted, e))?,
            });
        }
        Ok((fwd.unwrap(), bwd.unwrap()))
    }

    // ---------------------------------------------------------- local window decompositions

    /// exact decomposition of a window word into canonical ladders of its block (xengine's
    /// decompose_direct on the window), orthogonal; π' are certified added templates
    fn split<F: Field>(&mut self, ar: &mut Arena<F>, labels: &[u8], rel: &[Rung], st: &mut SearchStats) -> Result<Vec<Piece>, Stop> {
        let key = (labels.to_vec(), rel.to_vec());
        if let Some(r) = self.splits.get(&key) {
            return r.clone().map_err(|e| (MissKind::Exhausted, e));
        }
        let r = self.split_compute(ar, labels, rel, st);
        match &r {
            Ok(v) => {
                self.splits.insert(key, Ok(v.clone()));
            }
            Err((MissKind::Exhausted, e)) => {
                self.splits.insert(key, Err(e.clone()));
            }
            Err(_) => {}
        }
        r
    }
    fn split_compute<F: Field>(&mut self, ar: &mut Arena<F>, labels: &[u8], rel: &[Rung], st: &mut SearchStats) -> Result<Vec<Piece>, Stop> {
        let parts = self.sc.expand_word(labels, rel);
        let fw = crate::schur::flows(labels, rel) as usize;
        let mut io_t: Vec<Option<TmplId>> = vec![];
        let mut io_m: Vec<PMat<F>> = vec![];
        let mut pi_m: Vec<PMat<F>> = vec![];
        let mut words: Vec<(Vec<Rung>, i32)> = vec![];
        for (lw, sdeg, cnt, _) in &parts {
            // candidate ι (L_t -> W, degree s) and π (W -> L_t, degree -s) as exact PMats: the exact
            // space when admitted, otherwise dot multiples of the family extremes
            let fam: Vec<i32> = parts.iter().filter(|p| p.0 == *lw).map(|p| p.1).collect();
            let (s_lo, s_hi) = (*fam.iter().min().unwrap(), *fam.iter().max().unwrap());
            let (ic_m, ic_t) = self.local_cands(ar, labels, rel, lw, *sdeg, s_lo, true, st)?;
            let (pc_m, _) = self.local_cands(ar, labels, rel, lw, *sdeg, s_hi, false, st)?;
            let g: Vec<Vec<F>> = pc_m
                .iter()
                .map(|pj| ic_m.iter().map(|ii| pmat_scalar(&pj.matmul(ii)).ok_or_else(|| "π∘ι is not a scalar".to_string())).collect::<Result<Vec<F>, String>>())
                .collect::<Result<_, _>>()
                .map_err(|e| (MissKind::Exhausted, format!("window {:?} {:?}: {}", labels, rel, e)))?;
            let (ic, jr) = select_form(&g, ic_m.len(), *cnt).ok_or_else(|| (MissKind::Exhausted, format!("window {:?} {:?}: local form rank != {} for {:?}{{{}}} ({} / {} candidates)", labels, rel, cnt, lw, sdeg, ic_m.len(), pc_m.len())))?;
            let sub: Vec<Vec<F>> = jr.iter().map(|&j| ic.iter().map(|&i| g[j][i].clone()).collect()).collect();
            let inv = invert(&sub).ok_or_else(|| (MissKind::Exhausted, "local form minor singular".to_string()))?;
            for a in 0..*cnt {
                let mut pm = PMat::zero(crate::schur::flows(labels, lw) as usize, fw);
                for bb in 0..*cnt {
                    if !inv[a][bb].is_zero() {
                        pm = pm.add(&pc_m[jr[bb]].scale(&crate::mpoly::Poly::constant(inv[a][bb].clone())));
                    }
                }
                io_t.push(ic_t[ic[a]]);
                io_m.push(ic_m[ic[a]].clone());
                pi_m.push(pm);
                words.push((lw.clone(), *sdeg));
            }
        }
        let total: usize = words.iter().map(|(w, _)| crate::schur::flows(labels, w) as usize).sum();
        if total != fw {
            return Err((MissKind::Exhausted, format!("window {:?} {:?}: summand ranks {} != {}", labels, rel, total, fw)));
        }
        // orthogonalize: Π∘Ι = I + N, N nilpotent (nonzero only from lower to higher shift)
        let n = words.len();
        let mut nb: Vec<Vec<Option<PMat<F>>>> = vec![vec![None; n]; n];
        for t in 0..n {
            for s2 in 0..n {
                let m = pi_m[t].matmul(&io_m[s2]);
                if t == s2 {
                    if pmat_scalar(&m).map(|c| c.is_one()) != Some(true) {
                        return Err((MissKind::Exhausted, "window: π∘ι != id on a summand".into()));
                    }
                } else if !m.is_zero() {
                    if words[s2].1 <= words[t].1 {
                        return Err((MissKind::Exhausted, "window: nonzero π_t∘ι_s with s_s <= s_t".into()));
                    }
                    nb[t][s2] = Some(m.neg());
                }
            }
        }
        let fsz: Vec<usize> = words.iter().map(|(w, _)| crate::schur::flows(labels, w) as usize).collect();
        let mut minv: Vec<Vec<Option<PMat<F>>>> = vec![vec![None; n]; n];
        for t in 0..n {
            minv[t][t] = Some(PMat::identity(fsz[t]));
        }
        let mut pw = minv.clone();
        let mut iter = 0;
        loop {
            let mut next: Vec<Vec<Option<PMat<F>>>> = vec![vec![None; n]; n];
            let mut any = false;
            for t in 0..n {
                for s2 in 0..n {
                    let mut acc: Option<PMat<F>> = None;
                    for k in 0..n {
                        if let (Some(a), Some(b)) = (&nb[t][k], &pw[k][s2]) {
                            let m = a.matmul(b);
                            acc = Some(match acc {
                                None => m,
                                Some(x) => x.add(&m),
                            });
                        }
                    }
                    if let Some(m) = acc {
                        if !m.is_zero() {
                            next[t][s2] = Some(m);
                            any = true;
                        }
                    }
                }
            }
            if !any {
                break;
            }
            iter += 1;
            if iter > n {
                return Err((MissKind::Exhausted, "window: Π∘Ι not unipotent".into()));
            }
            for t in 0..n {
                for s2 in 0..n {
                    if let Some(z) = &next[t][s2] {
                        minv[t][s2] = Some(match minv[t][s2].take() {
                            None => z.clone(),
                            Some(w) => w.add(z),
                        });
                    }
                }
            }
            pw = next;
        }
        let old = pi_m.clone();
        for t in 0..n {
            let mut acc = PMat::zero(fsz[t], fw);
            for s2 in 0..n {
                if let Some(m) = &minv[t][s2] {
                    acc = acc.add(&m.matmul(&old[s2]));
                }
            }
            pi_m[t] = acc;
        }
        for t in 0..n {
            for s2 in 0..n {
                let c = pmat_scalar(&pi_m[t].matmul(&io_m[s2])).ok_or_else(|| (MissKind::Exhausted, "window: Π∘Ι not scalar".to_string()))?;
                if (t == s2 && !c.is_one()) || (t != s2 && !c.is_zero()) {
                    return Err((MissKind::Exhausted, "window: Π∘Ι != id after orthogonalization".into()));
                }
            }
        }
        let mut out = vec![];
        for t in 0..n {
            let (lw, sdeg) = words[t].clone();
            let tm = Template::new(Window::new(labels, rel), Window::new(labels, &lw), -sdeg, pi_m[t].clone(), TmplSource::Added { certified: false });
            let id = ar.lib().add_template(tm);
            ar.certify_template(id, 0x5eed_0100 + t as u64).map_err(|e| (MissKind::Exhausted, format!("window projection certificate: {}", e)))?;
            let iota = match io_t[t] {
                Some(x) => x,
                None => {
                    let tm = Template::new(Window::new(labels, &lw), Window::new(labels, rel), sdeg, io_m[t].clone(), TmplSource::Added { certified: false });
                    let i2 = ar.lib().add_template(tm);
                    ar.certify_template(i2, 0x5eed_0180 + t as u64).map_err(|e| (MissKind::Exhausted, format!("window inclusion certificate: {}", e)))?;
                    i2
                }
            };
            out.push(Piece { word: lw, s: sdeg, iota, pi: id });
        }
        Ok(out)
    }

    /// Candidate local maps for one summand (L_t{s}) of the window word W: ι (L_t -> W, degree s)
    /// when `iota`, else π (W -> L_t, degree -s). The exact basis when admitted; otherwise dot
    /// multiples (interior segments of W and of L_t, products of up to two Schur dots) of the
    /// family's extreme map (ι at s_lo or π at s_hi, `ext`). Second component: the hom_space
    /// template ids when the maps are unmodified exact templates.
    #[allow(clippy::too_many_arguments)]
    fn local_cands<F: Field>(&mut self, ar: &mut Arena<F>, labels: &[u8], w: &[Rung], lt: &[Rung], s: i32, ext: i32, iota: bool, st: &mut SearchStats) -> Result<(Vec<PMat<F>>, Vec<Option<TmplId>>), Stop> {
        let (a, b, d) = if iota { (lt, w, s) } else { (w, lt, -s) };
        if self.win_ok(ar, labels, a, b, d) {
            let ts = self.exact(ar, labels, a, b, d, st)?;
            let lib = ar.lib_ref();
            return Ok((ts.iter().map(|&t| lib.get(t).phi.clone()).collect(), ts.iter().map(|&t| Some(t)).collect()));
        }
        let gap = if iota { s - ext } else { ext - s };
        if gap <= 0 || gap % 2 != 0 {
            return Err((MissKind::Refused, format!("window {:?} {:?} -> {:?} degree {}: not admitted and no extreme to multiply", labels, a, b, d)));
        }
        let (ea, eb, ed) = if iota { (lt, w, ext) } else { (w, lt, -ext) };
        let ex_ts = self.exact(ar, labels, ea, eb, ed, st)?;
        let exm: Vec<PMat<F>> = ex_ts.iter().map(|&t| ar.lib_ref().get(t).phi.clone()).collect();
        let k = (gap / 2) as usize;
        let mut out = vec![];
        for (word_, on_w) in [(w, true), (lt, false)] {
            let segs = segments(labels, word_, true);
            for mono in dot_monomials(&segs, k, 2, 80) {
                let mut dm: Option<PMat<F>> = None;
                for (si, nu) in &mono {
                    let (level, u, _) = segs[*si];
                    let dp = dot_pmat::<F>(labels, word_, level, u, nu);
                    dm = Some(match dm {
                        None => dp,
                        Some(x) => dp.matmul(&x),
                    });
                }
                let Some(dm) = dm else { continue };
                for e in &exm {
                    // ι: Dot_W ∘ ι_lo or ι_lo ∘ Dot_Lt; π: π_hi ∘ Dot_W or Dot_Lt ∘ π_hi
                    let m = match (iota, on_w) {
                        (true, true) => dm.matmul(e),
                        (true, false) => e.matmul(&dm),
                        (false, true) => e.matmul(&dm),
                        (false, false) => dm.matmul(e),
                    };
                    if !m.is_zero() {
                        out.push(m);
                    }
                }
            }
        }
        let n = out.len();
        Ok((out, vec![None; n]))
    }

    /// the cheapest admissible decomposable window over the commutation class of z
    fn choose_window<F: Field>(&mut self, ar: &Arena<F>, z: &[Rung]) -> Option<(Vec<Rung>, WinPos, Vec<u8>, Vec<Rung>)> {
        let lam = ar.lam().to_vec();
        let n = lam.len();
        let mut best: Option<(u128, Vec<Rung>, WinPos, Vec<u8>, Vec<Rung>)> = None;
        self.last_win_refused = false;
        for zz in commutation_class(z, 64) {
            for p in windows_of(&zz, self.win_w.min(n), self.win_len) {
                let lv = level_of(&lam, &zz, p.at);
                let labels = lv[p.first..p.first + p.w].to_vec();
                let rel = window_word(&zz, &p);
                let key = (labels.clone(), rel.clone());
                let cached = self.splits.get(&key).map(|r| r.is_ok());
                if cached == Some(false) {
                    continue;
                }
                let parts = self.sc.expand_word(&labels, &rel);
                if parts.len() == 1 && parts[0].1 == 0 && parts[0].2 == 1 {
                    continue;
                }
                let mut cost = 0u128;
                let mut ok = true;
                if cached != Some(true) {
                    for (lw, sdeg, _, _) in &parts {
                        let fam: Vec<i32> = parts.iter().filter(|p| p.0 == *lw).map(|p| p.1).collect();
                        let (s_lo, s_hi) = (*fam.iter().min().unwrap(), *fam.iter().max().unwrap());
                        // ι at s (or at the family's lowest degree, then dots); π at -s (or at -s_hi)
                        let io_d = if self.win_ok(ar, &labels, lw, &rel, *sdeg) { *sdeg } else { s_lo };
                        let pi_d = if self.win_ok(ar, &labels, &rel, lw, -*sdeg) { -*sdeg } else { -s_hi };
                        if !self.win_ok(ar, &labels, lw, &rel, io_d) || !self.win_ok(ar, &labels, &rel, lw, pi_d) {
                            if std::env::var("KRHOM_TRACE").is_ok() {
                                eprintln!("[choose_window] refused {:?} [{}] summand [{}]{{{}}}: iota d {} ({} x {}, unknowns {}), pi d {} (unknowns {})", labels, wstr(&rel), wstr(lw), sdeg, io_d, crate::schur::flows(&labels, lw), crate::schur::flows(&labels, &rel), hom_unknowns(&labels, lw, &rel, io_d), pi_d, hom_unknowns(&labels, &rel, lw, pi_d));
                            }
                            ok = false;
                            self.last_win_refused = true;
                            break;
                        }
                        cost = cost.max(self.win_cost(ar, &labels, lw, &rel, io_d)).max(self.win_cost(ar, &labels, &rel, lw, pi_d));
                    }
                }
                if ok && best.as_ref().map_or(true, |b| cost < b.0) {
                    best = Some((cost, zz.clone(), p, labels, rel));
                }
            }
        }
        best.map(|(_, a, b, c, d)| (a, b, c, d))
    }

    /// local decomposition tree of x (§2 of the design note): leaves (word, type word, ι, π, shift)
    /// with π_a∘ι_b = δ_ab; with `want` = (type word, shift) only branches that can reach it
    fn tree<F: Field>(&mut self, ar: &mut Arena<F>, x: ExprId, want: Option<(&[Rung], i32)>, ctx: &mut Ctx) -> Result<Vec<Leaf>, Stop> {
        let lam = ar.lam().to_vec();
        let xw = word(ar, x);
        let id = ar.add(MapNode::Id(x), "tree root").map_err(|e| (MissKind::Exhausted, e))?;
        let flows_max = if ctx.b.flows_max == 0 { 4 * crate::schur::flows(&lam, &xw) as u64 } else { ctx.b.flows_max };
        let mut stack: Vec<(Vec<Rung>, MapId, MapId, i32, usize)> = vec![(xw, id, id, 0, 0)];
        let mut leaves = vec![];
        while let Some((z, io, pi, sh, depth)) = stack.pop() {
            ctx.check(ar)?;
            ctx.st.visited_words += 1;
            let fz = crate::schur::flows(&lam, &z) as u64;
            ctx.st.max_flows = ctx.st.max_flows.max(fz);
            let ex = self.sc.expand_word(&lam, &z);
            if let Some((qw, s)) = want {
                if !ex.iter().any(|e| e.0 == qw && e.1 == s - sh) {
                    continue;
                }
            }
            if ex.len() == 1 && ex[0].1 == 0 && ex[0].2 == 1 {
                leaves.push(Leaf { w: z, tword: ex[0].0.clone(), io, pi, shift: sh });
                continue;
            }
            if depth >= ctx.b.l_max {
                ctx.stuck += 1;
                self.note(format!("tree: depth {} reached at [{}]", depth, wstr(&z)));
                continue;
            }
            let Some((zz, p, labels, rel)) = self.choose_window(ar, &z) else {
                ctx.stuck_by(if self.last_win_refused { MissKind::Refused } else { MissKind::Exhausted });
                self.note(format!("tree: no admissible decomposable window in [{}]", wstr(&z)));
                continue;
            };
            let pieces = match self.split(ar, &labels, &rel, &mut ctx.st) {
                Ok(v) => v,
                Err((k, e)) => {
                    if k == MissKind::Time || k == MissKind::Memory {
                        return Err((k, e));
                    }
                    ctx.stuck_by(k);
                    self.note(format!("tree: split of {:?} [{}] failed: {}", labels, wstr(&rel), e));
                    continue;
                }
            };
            let (c, cinv) = match self.iso_chain(ar, &z, &zz, &mut ctx.st) {
                Ok(x) => x,
                Err((k, e)) => {
                    if k == MissKind::Time || k == MissKind::Memory {
                        return Err((k, e));
                    }
                    ctx.stuck_by(k);
                    self.note(format!("tree: commutation [{}] -> [{}] failed: {}", wstr(&z), wstr(&zz), e));
                    continue;
                }
            };
            let ioz = ar.add(MapNode::Comp(io, cinv), "tree").map_err(|e| (MissKind::Exhausted, e))?;
            let piz = ar.add(MapNode::Comp(c, pi), "tree").map_err(|e| (MissKind::Exhausted, e))?;
            let ez = eid(ar, &zz).map_err(|e| (MissKind::Exhausted, e))?;
            for pc in pieces.iter().rev() {
                let mut child: Vec<Rung> = zz[..p.at].to_vec();
                child.extend(shift_word(&pc.word, p.first as isize));
                child.extend_from_slice(&zz[p.at + p.len..]);
                if crate::schur::flows(&lam, &child) as u64 > flows_max {
                    continue;
                }
                let ec = eid(ar, &child).map_err(|e| (MissKind::Exhausted, e))?;
                let li = ar.add(MapNode::Local { src: ec, tgt: ez, at: p.at, first: p.first, tmpl: pc.iota }, "split iota").map_err(|e| (MissKind::Exhausted, e))?;
                let lp = ar.add(MapNode::Local { src: ez, tgt: ec, at: p.at, first: p.first, tmpl: pc.pi }, "split pi").map_err(|e| (MissKind::Exhausted, e))?;
                let ioc = ar.add(MapNode::Comp(ioz, li), "tree").map_err(|e| (MissKind::Exhausted, e))?;
                let pic = ar.add(MapNode::Comp(lp, piz), "tree").map_err(|e| (MissKind::Exhausted, e))?;
                stack.push((child, ioc, pic, sh + pc.s, depth + 1));
            }
        }
        Ok(leaves)
    }

    /// one degree-0 isomorphism move of a word, as mutually inverse maps (w -> w', w' -> w)
    fn iso_move<F: Field>(&mut self, ar: &mut Arena<F>, w: &[Rung], mv: &IsoMove, st: &mut SearchStats) -> Result<(MapId, MapId, Vec<Rung>), Stop> {
        match mv {
            IsoMove::Swap(i) => self.swap_maps(ar, w, *i, st),
            IsoMove::Rex(p, lw) => {
                let lam = ar.lam().to_vec();
                let labels = level_of(&lam, w, p.at)[p.first..p.first + p.w].to_vec();
                let rel = window_word(w, p);
                let f = self.xcall(ar, &labels, &rel, lw, 0, st).map_err(|e| (MissKind::Refused, e))?;
                let g = self.xcall(ar, &labels, lw, &rel, 0, st).map_err(|e| (MissKind::Refused, e))?;
                if f.len() != 1 || g.len() != 1 {
                    return Err((MissKind::Exhausted, format!("rex move {:?} [{}] -> [{}]: Hom^0 dims {} / {}", labels, wstr(&rel), wstr(lw), f.len(), g.len())));
                }
                let c = pmat_scalar(&ar.lib_ref().get(g[0]).phi.matmul(&ar.lib_ref().get(f[0]).phi)).filter(|c| !c.is_zero()).ok_or_else(|| (MissKind::Exhausted, "rex move: not an isomorphism".to_string()))?;
                let mut nw: Vec<Rung> = w[..p.at].to_vec();
                nw.extend(shift_word(lw, p.first as isize));
                nw.extend_from_slice(&w[p.at + p.len..]);
                let (a, b) = (eid(ar, w).map_err(|e| (MissKind::Exhausted, e))?, eid(ar, &nw).map_err(|e| (MissKind::Exhausted, e))?);
                let fm = ar.add(MapNode::Local { src: a, tgt: b, at: p.at, first: p.first, tmpl: f[0] }, "rex").map_err(|e| (MissKind::Exhausted, e))?;
                let gm = ar.add(MapNode::Local { src: b, tgt: a, at: p.at, first: p.first, tmpl: g[0] }, "rex back").map_err(|e| (MissKind::Exhausted, e))?;
                let gm = ar.add(MapNode::Lin(vec![(c.inv(), gm)]), "rex inverse").map_err(|e| (MissKind::Exhausted, e))?;
                Ok((fm, gm, nw))
            }
        }
    }
    /// words reachable from w by isomorphism moves (commutations; rewriting an indecomposable
    /// window to its canonical local word when both degree-0 window maps are admitted), breadth
    /// first up to `cap` words; parent pointers for the path
    fn iso_reach<F: Field>(&mut self, ar: &Arena<F>, w: &[Rung], cap: usize) -> HashMap<Vec<Rung>, Option<(Vec<Rung>, IsoMove)>> {
        let lam = ar.lam().to_vec();
        let n = lam.len();
        let mut prev: HashMap<Vec<Rung>, Option<(Vec<Rung>, IsoMove)>> = HashMap::new();
        prev.insert(w.to_vec(), None);
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(w.to_vec());
        while let Some(cur) = queue.pop_front() {
            if prev.len() >= cap {
                break;
            }
            let mut moves: Vec<(IsoMove, Vec<Rung>)> = vec![];
            for i in 0..cur.len().saturating_sub(1) {
                if rungs_commute(&cur[i], &cur[i + 1]) {
                    let mut nw = cur.clone();
                    nw.swap(i, i + 1);
                    moves.push((IsoMove::Swap(i), nw));
                }
            }
            for p in windows_of(&cur, self.win_w.min(n), self.win_len) {
                let labels = level_of(&lam, &cur, p.at)[p.first..p.first + p.w].to_vec();
                let rel = window_word(&cur, &p);
                let ex = self.sc.expand_word(&labels, &rel);
                if ex.len() != 1 || ex[0].1 != 0 || ex[0].2 != 1 || ex[0].0 == rel {
                    continue;
                }
                let lw = ex[0].0.clone();
                if !self.iso_ok(ar, &labels, &rel, &lw) || !self.iso_ok(ar, &labels, &lw, &rel) {
                    continue;
                }
                let mut nw: Vec<Rung> = cur[..p.at].to_vec();
                nw.extend(shift_word(&lw, p.first as isize));
                nw.extend_from_slice(&cur[p.at + p.len..]);
                moves.push((IsoMove::Rex(p, lw), nw));
            }
            for (mv, nw) in moves {
                if !prev.contains_key(&nw) {
                    prev.insert(nw.clone(), Some((cur.clone(), mv)));
                    queue.push_back(nw);
                }
            }
        }
        prev
    }
    /// mutually inverse maps (w -> y, y -> w) along the recorded moves from w to y
    fn iso_along<F: Field>(&mut self, ar: &mut Arena<F>, prev: &HashMap<Vec<Rung>, Option<(Vec<Rung>, IsoMove)>>, y: &[Rung], st: &mut SearchStats) -> Result<(MapId, MapId), Stop> {
        let mut steps = vec![];
        let mut cur = y.to_vec();
        while let Some(Some((p, mv))) = prev.get(&cur) {
            steps.push((p.clone(), mv.clone()));
            cur = p.clone();
        }
        steps.reverse();
        let e0 = eid(ar, &cur).map_err(|e| (MissKind::Exhausted, e))?;
        let id = ar.add(MapNode::Id(e0), "id").map_err(|e| (MissKind::Exhausted, e))?;
        let (mut fwd, mut bwd) = (id, id);
        for (w, mv) in steps {
            let (f, g, _) = self.iso_move(ar, &w, &mv, st)?;
            fwd = ar.add(MapNode::Comp(f, fwd), "iso path").map_err(|e| (MissKind::Exhausted, e))?;
            bwd = ar.add(MapNode::Comp(bwd, g), "iso path back").map_err(|e| (MissKind::Exhausted, e))?;
        }
        Ok((fwd, bwd))
    }

    /// degree-0 maps (leaf word -> type word, type word -> leaf word) with nonzero top scalar at
    /// the coset c (both are indecomposable of the same type: any such map is an isomorphism)
    fn leaf_iso<F: Field>(&mut self, ar: &mut Arena<F>, lw: &[Rung], qw: &[Rung], c: u64, st: &mut SearchStats) -> Result<(MapId, MapId), Stop> {
        let lam = ar.lam().to_vec();
        if commutation_class(qw, 500).iter().any(|w| w == lw) {
            return self.iso_chain(ar, lw, qw, st);
        }
        // meet in the middle over isomorphism moves (commutations and window rewrites)
        let ra = self.iso_reach(ar, lw, 400);
        let rb = self.iso_reach(ar, qw, 400);
        let mut common: Vec<&Vec<Rung>> = ra.keys().filter(|w| rb.contains_key(*w)).collect();
        common.sort_by_key(|w| (w.len(), wstr(w)));
        if let Some(m) = common.first().map(|w| (*w).clone()) {
            let (l2m, m2l) = self.iso_along(ar, &ra, &m, st)?;
            let (q2m, m2q) = self.iso_along(ar, &rb, &m, st)?;
            let f = ar.add(MapNode::Comp(m2q, l2m), "leaf iso").map_err(|e| (MissKind::Exhausted, e))?;
            let g = ar.add(MapNode::Comp(m2l, q2m), "leaf iso").map_err(|e| (MissKind::Exhausted, e))?;
            return Ok((f, g));
        }
        // best commutation-equivalent pair with an admissible window
        let mut best: Option<(u128, Vec<Rung>, Vec<Rung>, WinSplit)> = None;
        // cheapest window refused by admission (M5c: reported as Refused, and traced for
        // `libbuild calib --windows`, instead of Exhausted)
        let mut refused: Option<(u128, WinSplit)> = None;
        for a in commutation_class(lw, 64) {
            for b in commutation_class(qw, 64) {
                let ws = win_split(&lam, &a, &b);
                let cost = ws.fa * ws.fb;
                if best.as_ref().map_or(true, |x| cost < x.0) {
                    if self.iso_ok(ar, &ws.labels, &ws.a, &ws.b) && self.iso_ok(ar, &ws.labels, &ws.b, &ws.a) {
                        best = Some((cost, a.clone(), b.clone(), ws));
                    } else {
                        let u = hom_unknowns(&ws.labels, &ws.a, &ws.b, 0).max(hom_unknowns(&ws.labels, &ws.b, &ws.a, 0));
                        if refused.as_ref().map_or(true, |r| u < r.0) {
                            refused = Some((u, ws));
                        }
                    }
                }
            }
        }
        let Some((_, a, b, ws)) = best else {
            if let Some((u, ws)) = refused {
                if std::env::var("KRHOM_TRACE").is_ok() {
                    eprintln!("[leaf_iso] refused {:?} [{}] <-> [{}] degree 0 ({} x {}, unknowns {})", ws.labels, wstr(&ws.a), wstr(&ws.b), ws.fa, ws.fb, u);
                }
                return Err((MissKind::Refused, format!("leaf [{}] -> type word [{}]: no admissible window (cheapest refused: {:?} [{}] <-> [{}], {} unknowns)", wstr(lw), wstr(qw), ws.labels, wstr(&ws.a), wstr(&ws.b), u)));
            }
            return Err((MissKind::Exhausted, format!("leaf [{}] -> type word [{}]: no admissible window", wstr(lw), wstr(qw))));
        };
        let fwd = self.xcall(ar, &ws.labels, &ws.a, &ws.b, 0, st).map_err(|e| (MissKind::Refused, e))?;
        let bwd = self.xcall(ar, &ws.labels, &ws.b, &ws.a, 0, st).map_err(|e| (MissKind::Refused, e))?;
        let (ea, eb) = (eid(ar, &a).map_err(|e| (MissKind::Exhausted, e))?, eid(ar, &b).map_err(|e| (MissKind::Exhausted, e))?);
        let pt = points::<F>(ar.k(), 0x5eed_0200, 1, 0).remove(0);
        let pick = |ar: &mut Arena<F>, ts: &[TmplId], s: ExprId, t: ExprId| -> Result<MapId, Stop> {
            for &tm in ts {
                let m = ar.add(MapNode::Local { src: s, tgt: t, at: ws.at, first: ws.first, tmpl: tm }, "leaf iso").map_err(|e| (MissKind::Exhausted, e))?;
                let blk = Eval::new(ar, pt.clone()).block(m, c);
                if blk.rows == 1 && blk.cols == 1 && !blk.a[0].is_zero() {
                    return Ok(m);
                }
            }
            Err((MissKind::Exhausted, "leaf iso: no degree-0 window map with nonzero top scalar".into()))
        };
        let phi = pick(ar, &fwd, ea, eb)?;
        let psi = pick(ar, &bwd, eb, ea)?;
        let (l2a, a2l) = self.iso_chain(ar, lw, &a, st)?;
        let (b2q, q2b) = self.iso_chain(ar, &b, qw, st)?;
        let f1 = ar.add(MapNode::Comp(phi, l2a), "leaf iso").map_err(|e| (MissKind::Exhausted, e))?;
        let f = ar.add(MapNode::Comp(b2q, f1), "leaf iso").map_err(|e| (MissKind::Exhausted, e))?;
        let g1 = ar.add(MapNode::Comp(psi, q2b), "leaf iso").map_err(|e| (MissKind::Exhausted, e))?;
        let g = ar.add(MapNode::Comp(a2l, g1), "leaf iso").map_err(|e| (MissKind::Exhausted, e))?;
        Ok((f, g))
    }

    /// (ι: L_q -> X, π: X -> L_q) for every leaf of type q at shift s of the tree of x
    fn tree_pairs<F: Field>(&mut self, ar: &mut Arena<F>, x: ExprId, qw: &[Rung], s: i32, c: u64, ctx: &mut Ctx) -> Result<(Vec<MapId>, Vec<MapId>), Stop> {
        let leaves = self.tree(ar, x, Some((qw, s)), ctx)?;
        let mut io = vec![];
        let mut pi = vec![];
        for lf in leaves {
            if lf.tword != qw || lf.shift != s {
                continue;
            }
            if lf.w == qw {
                io.push(lf.io);
                pi.push(lf.pi);
                continue;
            }
            match self.leaf_iso(ar, &lf.w, qw, c, &mut ctx.st) {
                Ok((f, g)) => {
                    io.push(ar.add(MapNode::Comp(lf.io, g), "leaf").map_err(|e| (MissKind::Exhausted, e))?);
                    pi.push(ar.add(MapNode::Comp(f, lf.pi), "leaf").map_err(|e| (MissKind::Exhausted, e))?);
                }
                Err(e) => {
                    ctx.stuck_by(e.0);
                    self.note(format!("leaf iso failed: {}", e.1));
                }
            }
        }
        ctx.st.candidates += io.len() + pi.len();
        Ok((io, pi))
    }

    /// dot monomial maps of degree 2k on the segments of expression e (composition of Dot nodes)
    fn dot_maps<F: Field>(&self, ar: &mut Arena<F>, e: ExprId, k: usize, interior: bool, cap: usize) -> Result<Vec<MapId>, Stop> {
        let lam = ar.lam().to_vec();
        let w = word(ar, e);
        let segs = segments(&lam, &w, interior);
        let mut out = vec![];
        for mono in dot_monomials(&segs, k, 3, cap) {
            let mut m: Option<MapId> = None;
            for (si, nu) in mono {
                let (level, u, _) = segs[si];
                let d = ar.add(MapNode::Dot { e, level, upright: u, nu }, "dot").map_err(|e| (MissKind::Exhausted, e))?;
                m = Some(match m {
                    None => d,
                    Some(x) => ar.add(MapNode::Comp(d, x), "dots").map_err(|e| (MissKind::Exhausted, e))?,
                });
            }
            if let Some(m) = m {
                out.push(m);
            } else {
                out.push(ar.add(MapNode::Id(e), "dot 1").map_err(|e| (MissKind::Exhausted, e))?);
            }
        }
        Ok(out)
    }

    /// exact maps on the best non-common window of commutation-equivalent words of a and b
    /// (S6): Local(φ) composed with the commutation isomorphisms, degree d
    fn window_maps<F: Field>(&mut self, ar: &mut Arena<F>, a: ExprId, b: ExprId, d: i32, st: &mut SearchStats) -> Result<Vec<MapId>, Stop> {
        let lam = ar.lam().to_vec();
        let (wa, wb) = (word(ar, a), word(ar, b));
        let mut best: Option<(u128, Vec<Rung>, Vec<Rung>, WinSplit)> = None;
        let mut refused = 0usize;
        for x in commutation_class(&wa, 64) {
            for y in commutation_class(&wb, 64) {
                let mut ws = win_split(&lam, &x, &y);
                if ws.labels.is_empty() {
                    // equal words: the only window of Hom^d(X, X) is the whole ladder
                    ws = WinSplit { at: 0, suffix: 0, first: 0, labels: lam.clone(), a: x.clone(), b: y.clone(), fa: crate::schur::flows(&lam, &x), fb: crate::schur::flows(&lam, &y) };
                }
                let cost = hom_unknowns(&ws.labels, &ws.a, &ws.b, d);
                if best.as_ref().map_or(true, |z| cost < z.0) {
                    if self.win_ok(ar, &ws.labels, &ws.a, &ws.b, d) {
                        best = Some((cost, x.clone(), y.clone(), ws));
                    } else {
                        refused += 1;
                    }
                }
            }
        }
        let Some((_, x, y, ws)) = best else {
            return Err((if refused > 0 { MissKind::Refused } else { MissKind::Exhausted }, format!("S6: no admissible window ({} refused by admission)", refused)));
        };
        let ts = self.exact(ar, &ws.labels, &ws.a, &ws.b, d, st)?;
        if ts.is_empty() {
            return Ok(vec![]);
        }
        let (ex, ey) = (eid(ar, &x).map_err(|e| (MissKind::Exhausted, e))?, eid(ar, &y).map_err(|e| (MissKind::Exhausted, e))?);
        let (a2x, _) = self.iso_chain(ar, &wa, &x, st)?;
        let (y2b, _) = self.iso_chain(ar, &y, &wb, st)?;
        let mut out = vec![];
        for t in ts {
            let l = ar.add(MapNode::Local { src: ex, tgt: ey, at: ws.at, first: ws.first, tmpl: t }, "S6 window").map_err(|e| (MissKind::Exhausted, e))?;
            let m1 = ar.add(MapNode::Comp(l, a2x), "S6").map_err(|e| (MissKind::Exhausted, e))?;
            out.push(ar.add(MapNode::Comp(y2b, m1), "S6").map_err(|e| (MissKind::Exhausted, e))?);
        }
        Ok(out)
    }
}

/// s_ν as a polynomial in e_1..e_l (variables 0..l-1) by the dual Jacobi-Trudi determinant
pub fn schur_e_poly<F: Field>(nu: &[u8], l: usize) -> crate::mpoly::Poly<F> {
    use crate::mpoly::Poly;
    let nu: Vec<usize> = nu.iter().map(|&x| x as usize).filter(|&x| x > 0).collect();
    if nu.is_empty() {
        return Poly::one();
    }
    let conj: Vec<usize> = (1..=nu[0]).map(|i| nu.iter().filter(|&&x| x >= i).count()).collect();
    let n = conj.len();
    let ent = |i: usize, j: usize| -> Poly<F> {
        let k = conj[i] as isize - i as isize + j as isize;
        if k == 0 {
            Poly::one()
        } else if k < 0 || k as usize > l {
            Poly::zero()
        } else {
            Poly::var(k as usize - 1)
        }
    };
    // permutation expansion (n <= 8 here)
    fn rec<F: Field>(i: usize, n: usize, used: &mut Vec<bool>, sign: bool, acc: Poly<F>, ent: &dyn Fn(usize, usize) -> Poly<F>, out: &mut Poly<F>) {
        if acc.is_zero() {
            return;
        }
        if i == n {
            *out = if sign { out.sub(&acc) } else { out.add(&acc) };
            return;
        }
        let mut inv = 0;
        for j in (0..n).rev() {
            if used[j] {
                inv += 1;
                continue;
            }
            let e = ent(i, j);
            if !e.is_zero() {
                used[j] = true;
                rec(i + 1, n, used, sign ^ (inv % 2 == 1), acc.mul(&e), ent, out);
                used[j] = false;
            }
        }
    }
    let mut out = Poly::zero();
    rec(0, n, &mut vec![false; n], false, Poly::one(), &ent, &mut out);
    out
}

/// exact PMat of the dot s_ν on upright u at `level` of the word w over `labels`
/// (tensor_exact of a one-upright scalar template on the prefix, identity on the suffix)
pub fn dot_pmat<F: Field>(labels: &[u8], w: &[Rung], level: usize, u: usize, nu: &[u8]) -> PMat<F> {
    let lv = level_of(labels, w, level);
    let pre: crate::sbim::SBim<F> = Expr::new(labels, &w[..level]).expect("dot_pmat: prefix").sbim();
    let win = Window::new(&[lv[u]], &[]);
    let deg: i32 = 2 * nu.iter().map(|&x| x as i32).sum::<i32>();
    let t = Template::new(win.clone(), win, deg, PMat::scalar(1, &schur_e_poly::<F>(nu, lv[u] as usize)), TmplSource::Added { certified: false });
    let fq = crate::schur::flows(&lv, &w[level..]) as usize;
    crate::tmpl::tensor_exact(&pre, u, &t, fq)
}

pub fn wstr(rs: &[Rung]) -> String {
    rs.iter().map(|r| format!("{}{}({})", if r.dir == Dir::E { "E" } else { "F" }, r.i, r.r)).collect::<Vec<_>>().join(" ")
}

/// incremental dense echelon over stacked rows (full / end0)
struct DenseEch<F: Field> {
    rows: Vec<(usize, Vec<F>)>,
}
impl<F: Field> DenseEch<F> {
    fn new() -> Self {
        DenseEch { rows: vec![] }
    }
    fn insert(&mut self, mut v: Vec<F>) -> bool {
        for (p, r) in &self.rows {
            if !v[*p].is_zero() {
                let f = v[*p].neg();
                for (x, y) in v.iter_mut().zip(r) {
                    if !y.is_zero() {
                        x.fma(&f, y);
                    }
                }
            }
        }
        match v.iter().position(|x| !x.is_zero()) {
            None => false,
            Some(p) => {
                let inv = v[p].inv();
                for x in v.iter_mut() {
                    *x = x.mul(&inv);
                }
                // keep earlier rows reduced at the new pivot
                for (_, r) in self.rows.iter_mut() {
                    if !r[p].is_zero() {
                        let f = r[p].neg();
                        for (x, y) in r.iter_mut().zip(&v) {
                            if !y.is_zero() {
                                x.fma(&f, y);
                            }
                        }
                    }
                }
                self.rows.push((p, v));
                true
            }
        }
    }
}

impl Search {
    fn pairs_inner<F: Field>(&mut self, ar: &mut Arena<F>, s: &mut Schur, x: ExprId, q: Tid, lq: ExprId, eq: Option<MapId>, sdeg: i32, mult: usize, ctx: &mut Ctx) -> Result<Pairs<F>, Stop> {
        pairs_precheck(ar, s, x, q, lq, eq, sdeg, mult).map_err(|e| (MissKind::Refused, e))?;
        let lam = ar.lam().to_vec();
        let (xw, qw) = (word(ar, x), word(ar, lq));
        let c = rep_rank(s, &lam, q);
        let chx = s.ladder_char(&xw);
        let mq = s.expand(&chx).into_iter().find(|e| e.0 == q).map(|e| e.1).ok_or_else(|| (MissKind::Refused, "q is not a summand of X".to_string()))?;
        let (s_lo, s_hi) = (mq.lo, mq.lo + mq.c.len() as i32 - 1);
        let seed = 0x5eed_1000 ^ (sdeg as i64 as u64);
        let finish = |ar: &mut Arena<F>, mut io: Vec<MapId>, mut pi: Vec<MapId>| -> Result<Pairs<F>, String> {
            with_idem(ar, eq, &mut io, &mut pi)?;
            certify_pairs(ar, x, lq, sdeg, &io, &pi, c, mult, seed)
        };
        // S0: identity
        if self.on(0) && xw == qw && sdeg == 0 && mult == 1 {
            let id = ar.add(MapNode::Id(x), "S0 identity").map_err(|e| (MissKind::Exhausted, e))?;
            ctx.st.candidates += 1;
            ctx.st.points = 2;
            if let Ok(p) = finish(ar, vec![id], vec![id]) {
                ctx.st.tier_hits[0] += 1;
                return Ok(p);
            }
        }
        // S1: extremes (one leaf of the local decomposition tree at s_lo or s_hi)
        if self.on(1) && (sdeg == s_lo || sdeg == s_hi) {
            let (io, pi) = self.tree_pairs(ar, x, &qw, sdeg, c, ctx)?;
            self.note(format!("S1: {} leaves at the extreme {}", io.len(), sdeg));
            match finish(ar, io, pi) {
                Ok(p) => {
                    ctx.st.tier_hits[1] += 1;
                    return Ok(p);
                }
                Err(e) => self.note(format!("S1 failed: {}", e)),
            }
        }
        // S2: dot multiples of the extremes
        if self.on(2) && s_lo < sdeg && sdeg < s_hi && (sdeg - s_lo) % 2 == 0 && (s_hi - sdeg) % 2 == 0 {
            let (io_lo, _) = self.tree_pairs(ar, x, &qw, s_lo, c, ctx)?;
            let (_, pi_hi) = self.tree_pairs(ar, x, &qw, s_hi, c, ctx)?;
            if let (Some(&il), Some(&ph)) = (io_lo.first(), pi_hi.first()) {
                let (ki, kp) = (((sdeg - s_lo) / 2) as usize, ((s_hi - sdeg) / 2) as usize);
                let cap = ctx.b.c_max.min(400);
                let mut ic = vec![];
                for d in self.dot_maps(ar, x, ki, false, cap)? {
                    ic.push(ar.add(MapNode::Comp(d, il), "S2 dot iota").map_err(|e| (MissKind::Exhausted, e))?);
                }
                for d in self.dot_maps(ar, lq, ki, false, cap)? {
                    ic.push(ar.add(MapNode::Comp(il, d), "S2 iota dot").map_err(|e| (MissKind::Exhausted, e))?);
                }
                let mut pc = vec![];
                for d in self.dot_maps(ar, x, kp, false, cap)? {
                    pc.push(ar.add(MapNode::Comp(ph, d), "S2 pi dot").map_err(|e| (MissKind::Exhausted, e))?);
                }
                for d in self.dot_maps(ar, lq, kp, false, cap)? {
                    pc.push(ar.add(MapNode::Comp(d, ph), "S2 dot pi").map_err(|e| (MissKind::Exhausted, e))?);
                }
                ctx.check(ar)?;
                // greedy: keep candidates raising the rank of their top-coset vectors
                let pt = points::<F>(ar.k(), seed ^ 0x77, 1, 0).remove(0);
                let mut ev = Eval::new(ar, pt);
                let mut keep_i = vec![];
                let mut ech = DenseEch::new();
                for &m in &ic {
                    ctx.st.evaluated += 1;
                    if ech.insert(ev.block(m, c).a.clone()) {
                        keep_i.push(m);
                    }
                }
                let mut keep_p = vec![];
                let mut ech = DenseEch::new();
                for &m in &pc {
                    ctx.st.evaluated += 1;
                    if ech.insert(ev.block(m, c).a.clone()) {
                        keep_p.push(m);
                    }
                }
                drop(ev);
                ctx.st.candidates += ic.len() + pc.len();
                self.note(format!("S2: {} / {} independent dot multiples (iota / pi)", keep_i.len(), keep_p.len()));
                match finish(ar, keep_i, keep_p) {
                    Ok(p) => {
                        ctx.st.tier_hits[2] += 1;
                        return Ok(p);
                    }
                    Err(e) => self.note(format!("S2 failed: {}", e)),
                }
            }
        }
        // S3: the decomposition tree at the target degree
        if self.on(3) {
            let (io, pi) = self.tree_pairs(ar, x, &qw, sdeg, c, ctx)?;
            self.note(format!("S3: {} leaves at {}", io.len(), sdeg));
            match finish(ar, io, pi) {
                Ok(p) => {
                    ctx.st.tier_hits[3] += 1;
                    return Ok(p);
                }
                Err(e) => self.note(format!("S3 failed: {}", e)),
            }
        }
        // S6: exact windows after commutation canonicalization
        if self.on(6) {
            let r = (|| -> Result<Pairs<F>, Stop> {
                let io = self.window_maps(ar, lq, x, sdeg, &mut ctx.st)?;
                let pi = self.window_maps(ar, x, lq, -sdeg, &mut ctx.st)?;
                ctx.st.candidates += io.len() + pi.len();
                finish(ar, io, pi).map_err(|e| (MissKind::Exhausted, e))
            })();
            match r {
                Ok(p) => {
                    ctx.st.tier_hits[6] += 1;
                    return Ok(p);
                }
                Err(e) => {
                    if e.0 == MissKind::Refused {
                        ctx.tier_refused += 1;
                    }
                    self.note(format!("S6 failed: {}", e.1))
                }
            }
        }
        let k = ctx.end_kind();
        Err((k, format!("tiers S{}-S{} {} ({})", self.lo, self.hi, if k == MissKind::Refused { "refused by admission" } else { "exhausted" }, ctx.end_note())))
    }

    fn closure_inner<F: Field>(&mut self, ar: &mut Arena<F>, s: &mut Schur, lp: ExprId, p: Tid, ep: Option<MapId>, outs: bool, ctx: &mut Ctx, depth: usize) -> Result<Vec<MapId>, Stop> {
        let mut seeds: Option<Vec<(MapId, MapId, ExprId)>> = None;
        let lam = ar.lam().to_vec();
        let want = cdiag_degrees(s, p);
        if depth == 0 {
            closure_precheck(ar, s, lp, p, ep).map_err(|e| (MissKind::Refused, e))?;
        }
        let r0 = eid(ar, &[]).map_err(|e| (MissKind::Exhausted, e))?;
        let mut sel: Selector<F> = Selector::new(identity_coset(&lam), ar.k(), want.len(), 0x5eed_2000, ep.is_some());
        let mut degs = want.clone();
        degs.dedup();
        let wrap = |ar: &mut Arena<F>, ms: Vec<MapId>| -> Result<Vec<MapId>, Stop> {
            match ep {
                None => Ok(ms),
                Some(e) => ms.into_iter().map(|m| ar.add(if outs { MapNode::Comp(m, e) } else { MapNode::Comp(e, m) }, "e_p").map_err(|e| (MissKind::Exhausted, e))).collect(),
            }
        };
        let d0 = degs[0];
        if self.on(0) && word(ar, lp).is_empty() && want == vec![0] {
            let id = ar.add(MapNode::Id(lp), "S0 closure identity").map_err(|e| (MissKind::Exhausted, e))?;
            let id = wrap(ar, vec![id])?;
            ctx.st.candidates += 1;
            if sel.offer_many(ar, &id, 1) == 1 {
                ctx.st.tier_hits[0] += 1;
                let ms = sel.chosen.clone();
                certify_closure(ar, &lam, r0, lp, ep, &ms, &want, outs, 0x5eed_2001).map_err(|e| (MissKind::Exhausted, e))?;
                return Ok(ms);
            }
        }
        for d in degs {
            ctx.check(ar)?;
            let need = want.iter().filter(|&&x| x == d).count();
            let mut got = 0;
            // S2: dot multiples of chosen generators of lower degree (interior segments of L_p)
            if self.on(2) && d > d0 {
                let chosen: Vec<MapId> = sel.chosen.clone();
                for g in chosen {
                    let dg = ar.info(g).deg;
                    if (d - dg) % 2 != 0 || d < dg {
                        continue;
                    }
                    let k = ((d - dg) / 2) as usize;
                    let dots = self.dot_maps(ar, lp, k, true, 200)?;
                    let mut cs = vec![];
                    for dt in dots {
                        cs.push(ar.add(if outs { MapNode::Comp(g, dt) } else { MapNode::Comp(dt, g) }, "S2 closure").map_err(|e| (MissKind::Exhausted, e))?);
                    }
                    // a dot need not commute with e_p: re-wrap so the candidate lies in the
                    // image of e_p again (review F4)
                    let cs = wrap(ar, cs)?;
                    ctx.st.candidates += cs.len();
                    ctx.st.evaluated += cs.len();
                    got += sel.offer_many(ar, &cs, need - got);
                    if got >= need {
                        break;
                    }
                }
                if got >= need {
                    ctx.st.tier_hits[2] += 1;
                    continue;
                }
            }
            // S1 (lowest degree) / S6: exact window [] <-> L_p at degree d
            let tier = if d == d0 { 1 } else { 6 };
            if self.on(tier) {
                let (a, b) = if outs { (lp, r0) } else { (r0, lp) };
                match self.window_maps(ar, a, b, d, &mut ctx.st) {
                    Ok(ms) => {
                        let ms = wrap(ar, ms)?;
                        ctx.st.candidates += ms.len();
                        ctx.st.evaluated += ms.len();
                        got += sel.offer_many(ar, &ms, need - got);
                    }
                    Err(e) => {
                        if e.0 == MissKind::Refused {
                            ctx.tier_refused += 1;
                        }
                        self.note(format!("closure degree {}: {}", d, e.1))
                    }
                }
                if got >= need {
                    ctx.st.tier_hits[tier] += 1;
                    continue;
                }
            }
            // S4: double-leaf shape through a smaller type, with dots on interior segments of L_p
            if self.on(4) {
                if seeds.is_none() {
                    seeds = Some(self.closure_seeds(ar, s, p, lp, outs, ctx, depth)?);
                }
                let sd = seeds.clone().unwrap();
                'seeds: for (left, right, ey) in sd {
                    let dg = ar.info(left).deg + ar.info(right).deg;
                    if d < dg || (d - dg) % 2 != 0 {
                        continue;
                    }
                    let k = ((d - dg) / 2) as usize;
                    // dots on the intermediate word Y (interior segments), between the two halves
                    let dots = self.dot_maps(ar, ey, k, true, 300)?;
                    let mut cs = vec![];
                    for dt in dots {
                        let a = ar.add(MapNode::Comp(dt, right), "S4 closure").map_err(|e| (MissKind::Exhausted, e))?;
                        cs.push(ar.add(MapNode::Comp(left, a), "S4 closure").map_err(|e| (MissKind::Exhausted, e))?);
                    }
                    let cs = wrap(ar, cs)?;
                    ctx.st.candidates += cs.len();
                    ctx.st.evaluated += cs.len();
                    got += sel.offer_many(ar, &cs, need - got);
                    if got >= need {
                        break 'seeds;
                    }
                    ctx.check(ar)?;
                }
                if got >= need {
                    ctx.st.tier_hits[4] += 1;
                    continue;
                }
            }
            let k = ctx.end_kind();
            return Err((k, format!("closure degree {}: {} of {} generators ({})", d, got, need, ctx.end_note())));
        }
        let ms = sel.chosen.clone();
        certify_closure(ar, &lam, r0, lp, ep, &ms, &want, outs, 0x5eed_2001).map_err(|e| (MissKind::Exhausted, e))?;
        Ok(ms)
    }

    /// S4 closure seeds of type p: π_leaf ∘ (id ⊗ unit of [R, R^∨]) ∘ g' (ins) or
    /// g' ∘ (id ⊗ counit) ∘ ι_leaf (outs), for a smaller type p' with p a summand of
    /// b_{p'}·[R, R^∨], g' running over the closure generators of p' (recursive)
    #[allow(clippy::too_many_arguments)]
    fn closure_seeds<F: Field>(&mut self, ar: &mut Arena<F>, s: &mut Schur, p: Tid, lp: ExprId, outs: bool, ctx: &mut Ctx, depth: usize) -> Result<Vec<(MapId, MapId, ExprId)>, Stop> {
        if depth >= 4 {
            return Ok(vec![]);
        }
        let lam = ar.lam().to_vec();
        let n = lam.len();
        let pw = word(ar, lp);
        let rank_p = s.rank(&lam, p);
        let order = s.canonical_basis(&lam).order.clone();
        let mut opts: Vec<(u128, i32, Tid, Rung, Rung, i32)> = vec![];
        for &p2 in &order {
            let r2 = s.rank(&lam, p2);
            if r2 >= rank_p {
                continue;
            }
            for i in 0..n.saturating_sub(1) {
                for dir in [Dir::E, Dir::F] {
                    let giv = if dir == Dir::E { lam[i + 1] } else { lam[i] };
                    for r in 1..=giv as usize {
                        let rr = Rung { i, dir, r };
                        let rv = Rung { i, dir: if dir == Dir::E { Dir::F } else { Dir::E }, r };
                        for (q, m) in s.term_mult(p2, &[rr, rv]) {
                            if q != p {
                                continue;
                            }
                            for (k, &c) in m.c.iter().enumerate() {
                                if c > 0 {
                                    let t = m.lo + k as i32;
                                    opts.push((u128::MAX - r2, t.abs(), p2, rr, rv, t));
                                }
                            }
                        }
                    }
                }
            }
        }
        opts.sort_by(|a, b| (a.0, a.1, wstr(&[a.3])).cmp(&(b.0, b.1, wstr(&[b.3]))));
        let mut seeds = vec![];
        for (_, _, p2, rr, rv, t) in opts.into_iter().take(3) {
            ctx.check(ar)?;
            let w2 = s.canonical_basis(&lam).ladder_minflow[&p2].clone();
            let e2 = eid(ar, &w2).map_err(|e| (MissKind::Exhausted, e))?;
            let g2s = match self.closure_inner(ar, s, e2, p2, None, outs, ctx, depth + 1) {
                Ok(v) => v,
                Err((k, e)) => {
                    if k == MissKind::Time || k == MissKind::Memory {
                        return Err((k, e));
                    }
                    self.note(format!("S4 closure: generators of {} failed: {}", s.table(p2).fmt(), e));
                    continue;
                }
            };
            // the closed window [R, R^∨] on uprights i, i+1 at the top (labels lambda)
            let i = rr.i;
            let wl = vec![lam[i], lam[i + 1]];
            let rel = vec![Rung { i: 0, dir: rr.dir, r: rr.r }, Rung { i: 0, dir: rv.dir, r: rv.r }];
            let dmin = {
                let sch = self.sc.get(&wl);
                let u = sch.unit();
                let cw = sch.ladder_char(&rel);
                let gd = sch.hom_gd(&u, &cw, -60, 60);
                match gd.iter().position(|&x| x > 0) {
                    Some(k) => -60 + k as i32,
                    None => continue,
                }
            };
            let units = if outs { self.exact(ar, &wl, &rel, &[], dmin, &mut ctx.st) } else { self.exact(ar, &wl, &[], &rel, dmin, &mut ctx.st) };
            let Ok(units) = units else { continue };
            let mut yw = w2.clone();
            yw.extend_from_slice(&[rr, rv]);
            let ey = eid(ar, &yw).map_err(|e| (MissKind::Exhausted, e))?;
            let leaves = self.tree(ar, ey, Some((&pw, t)), ctx)?;
            let qc = identity_coset_top(&lam, s, p);
            for lf in leaves {
                if lf.tword != pw || lf.shift != t {
                    continue;
                }
                let (pi_f, io_f) = if lf.w == pw {
                    (lf.pi, lf.io)
                } else {
                    match self.leaf_iso(ar, &lf.w, &pw, qc, &mut ctx.st) {
                        Ok((f, g)) => (
                            ar.add(MapNode::Comp(f, lf.pi), "S4 leaf").map_err(|e| (MissKind::Exhausted, e))?,
                            ar.add(MapNode::Comp(lf.io, g), "S4 leaf").map_err(|e| (MissKind::Exhausted, e))?,
                        ),
                        Err(e) => {
                            self.note(format!("S4 closure leaf iso: {}", e.1));
                            continue;
                        }
                    }
                };
                for &u in &units {
                    let un = if outs {
                        ar.add(MapNode::Local { src: ey, tgt: e2, at: w2.len(), first: i, tmpl: u }, "S4 counit")
                    } else {
                        ar.add(MapNode::Local { src: e2, tgt: ey, at: w2.len(), first: i, tmpl: u }, "S4 unit")
                    }
                    .map_err(|e| (MissKind::Exhausted, e))?;
                    for &g2 in &g2s {
                        // candidate = left ∘ Dot_Y ∘ right
                        let (left, right) = if outs {
                            (ar.add(MapNode::Comp(g2, un), "S4").map_err(|e| (MissKind::Exhausted, e))?, io_f)
                        } else {
                            (pi_f, ar.add(MapNode::Comp(un, g2), "S4").map_err(|e| (MissKind::Exhausted, e))?)
                        };
                        seeds.push((left, right, ey));
                    }
                }
            }
        }
        seeds.sort_by_key(|&(l, r, _)| ar.info(l).deg + ar.info(r).deg);
        self.note(format!("S4 closure seeds for {}: {} (degrees {:?})", s.table(p).fmt(), seeds.len(), seeds.iter().map(|&(l, r, _)| ar.info(l).deg + ar.info(r).deg).collect::<Vec<_>>()));
        Ok(seeds)
    }

    fn full_inner<F: Field>(&mut self, ar: &mut Arena<F>, s: &mut Schur, src: ExprId, tgt: ExprId, d: i32, ctx: &mut Ctx) -> Result<Vec<MapId>, Stop> {
        let plan = detect_plan(ar, s, src, tgt, d);
        ctx.st.detect = plan.kind;
        ctx.st.points = plan.npoints;
        let dim = plan.dim;
        if dim == 0 {
            return Ok(vec![]);
        }
        if plan.rows_per_point * plan.npoints < dim {
            return Err((MissKind::Exhausted, "counting bound violated".into()));
        }
        let pts = points::<F>(ar.k(), 0x5eed_3000, plan.npoints, 0);
        let mut ech: DenseEch<F> = DenseEch::new();
        let mut kept: Vec<MapId> = vec![];
        // every candidate ranked, for the orbit fallback below
        let offered: std::cell::RefCell<Vec<MapId>> = std::cell::RefCell::new(vec![]);
        let offer = |ar: &Arena<F>, ms: &[MapId], kept: &mut Vec<MapId>, ech: &mut DenseEch<F>, ctx: &mut Ctx| {
            for &m in ms {
                if kept.len() >= dim {
                    break;
                }
                offered.borrow_mut().push(m);
                ctx.st.evaluated += 1;
                let row = Eval::stacked(ar, &[m], &pts, &plan.cosets).remove(0);
                if ech.insert(row) {
                    kept.push(m);
                }
            }
        };
        // S0: identity
        if self.on(0) && src == tgt && d == 0 {
            let id = ar.add(MapNode::Id(src), "S0 identity").map_err(|e| (MissKind::Exhausted, e))?;
            offer(ar, &[id], &mut kept, &mut ech, ctx);
            if kept.len() == dim {
                ctx.st.tier_hits[0] += 1;
                return Ok(kept);
            }
        }
        // S2: dots on src (endomorphisms) when src == tgt
        if self.on(2) && src == tgt && d > 0 && d % 2 == 0 {
            let ms = self.dot_maps(ar, src, (d / 2) as usize, false, ctx.b.c_max.min(2000))?;
            ctx.st.candidates += ms.len();
            offer(ar, &ms, &mut kept, &mut ech, ctx);
            if kept.len() == dim {
                ctx.st.tier_hits[2] += 1;
                return Ok(kept);
            }
        }
        // S3/S4: leaves of the trees of src and tgt with equal types, dots in between
        if self.on(3) || self.on(4) {
            let ls = self.tree(ar, src, None, ctx)?;
            let lt = self.tree(ar, tgt, None, ctx)?;
            for a in &ls {
                for bl in &lt {
                    ctx.check(ar)?;
                    if a.tword != bl.tword {
                        continue;
                    }
                    let k2 = d - (bl.shift - a.shift);
                    if k2 < 0 || k2 % 2 != 0 {
                        continue;
                    }
                    let iso = if a.w == bl.w {
                        let e = eid(ar, &a.w).map_err(|e| (MissKind::Exhausted, e))?;
                        ar.add(MapNode::Id(e), "id").map_err(|e| (MissKind::Exhausted, e))?
                    } else if commutation_class(&bl.w, 500).iter().any(|w| *w == a.w) {
                        self.iso_chain(ar, &a.w, &bl.w, &mut ctx.st)?.0
                    } else {
                        continue;
                    };
                    let ea = eid(ar, &a.w).map_err(|e| (MissKind::Exhausted, e))?;
                    let dots = self.dot_maps(ar, ea, (k2 / 2) as usize, false, 200)?;
                    let mut cs = vec![];
                    for dt in dots {
                        let m1 = ar.add(MapNode::Comp(dt, a.pi), "S4").map_err(|e| (MissKind::Exhausted, e))?;
                        let m2 = ar.add(MapNode::Comp(iso, m1), "S4").map_err(|e| (MissKind::Exhausted, e))?;
                        cs.push(ar.add(MapNode::Comp(bl.io, m2), "S4").map_err(|e| (MissKind::Exhausted, e))?);
                    }
                    ctx.st.candidates += cs.len();
                    offer(ar, &cs, &mut kept, &mut ech, ctx);
                    if kept.len() == dim {
                        ctx.st.tier_hits[4] += 1;
                        return Ok(kept);
                    }
                }
            }
        }
        // S6: exact windows (degree d, and degree d - 2k followed by dots on tgt)
        if self.on(6) {
            for k in 0..=(d.max(0) / 2) {
                let dd = d - 2 * k;
                match self.window_maps(ar, src, tgt, dd, &mut ctx.st) {
                    Ok(ws) => {
                        let mut cs = vec![];
                        if k == 0 {
                            cs = ws;
                        } else {
                            let dots = self.dot_maps(ar, tgt, k as usize, false, 100)?;
                            for w in &ws {
                                for &dt in &dots {
                                    cs.push(ar.add(MapNode::Comp(dt, *w), "S6 dots").map_err(|e| (MissKind::Exhausted, e))?);
                                }
                            }
                        }
                        ctx.st.candidates += cs.len();
                        offer(ar, &cs, &mut kept, &mut ech, ctx);
                    }
                    Err(e) => {
                        if e.0 == MissKind::Refused {
                            ctx.tier_refused += 1;
                        }
                        self.note(format!("S6 degree {}: {}", dd, e.1))
                    }
                }
                if kept.len() == dim {
                    ctx.st.tier_hits[6] += 1;
                    return Ok(kept);
                }
                ctx.check(ar)?;
            }
        }
        // the planned detecting set is not faithful in general (hom60 (1,3) degree 6: witness rank
        // 141 < 145): before reporting Exhausted, re-rank every candidate at the orbit cosets with
        // the fallback point counts, as certify_full_plan does (review F6)
        let all = offered.into_inner();
        let (cos, _, np0) = orbit_fallback(ar, s, src, tgt, &plan);
        for (k, np) in [np0, 2 * np0].into_iter().enumerate() {
            ctx.check(ar)?;
            let pts = points::<F>(ar.k(), 0x5eed_3100 ^ k as u64, np, 0);
            let mut ech: DenseEch<F> = DenseEch::new();
            let mut kept2: Vec<MapId> = vec![];
            for &m in &all {
                if kept2.len() >= dim {
                    break;
                }
                let row = Eval::stacked(ar, &[m], &pts, &cos).remove(0);
                if ech.insert(row) {
                    kept2.push(m);
                }
            }
            self.note(format!("full Hom^{}: orbit re-rank ({} cosets, {} points): {} of {} ({} candidates)", d, cos.len(), np, kept2.len(), dim, all.len()));
            if kept2.len() == dim {
                ctx.st.detect = if k == 0 { "orbits-fallback" } else { "orbits-fallback-2x" };
                ctx.st.points = np;
                return Ok(kept2);
            }
        }
        let k = ctx.end_kind();
        Err((k, format!("full Hom^{}: rank {} of {} ({}; orbit re-rank also below)", d, kept.len(), dim, ctx.end_note())))
    }
}

fn run_call<F: Field, T>(ar: &mut Arena<F>, b: &Budget, f: impl FnOnce(&mut Arena<F>, &mut Ctx) -> Result<T, Stop>, keep_of: impl Fn(&T) -> Vec<MapId>, rebuild: impl Fn(T, Vec<MapId>) -> T) -> Result<(T, SearchStats), (Stop, SearchStats)> {
    let mark = ar.mark();
    let n_expr0 = ar.n_exprs();
    let mut ctx = Ctx { t0: Instant::now(), b: b.clone(), st: SearchStats::default(), stuck: 0, stuck_adm: 0, tier_refused: 0 };
    let r = f(ar, &mut ctx);
    ctx.st.indexed_words = ar.n_exprs().saturating_sub(n_expr0);
    ctx.st.secs = elapsed(&ctx.t0);
    match r {
        Ok(v) => {
            let keep = keep_of(&v);
            let k = ar.commit(mark, &keep);
            Ok((rebuild(v, k), ctx.st))
        }
        Err(e) => {
            ar.commit(mark, &[]);
            Err((e, ctx.st))
        }
    }
}

impl<F: Field> HomProvider<F> for Search {
    fn name(&self) -> &'static str {
        "Search"
    }
    fn full(&mut self, ar: &mut Arena<F>, s: &mut Schur, src: ExprId, tgt: ExprId, d: i32, b: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        self.trace.clear();
        self.bind(ar);
        self.failed = None;
        let r = run_call(ar, b, |ar, ctx| {
            let ms = self.full_inner(ar, s, src, tgt, d, ctx)?;
            certify_full(ar, s, src, tgt, d, &ms, 0x5eed_3001).map_err(|e| (MissKind::Exhausted, e))?;
            Ok(ms)
        }, |v| v.clone(), |_, k| k);
        r.map_err(|((k, e), st)| {
            let m = miss(&format!("Search full: {} ({} candidates)", e, st.candidates), d, None, 0, 0, k);
            self.failed = Some(st);
            m
        })
    }
    fn pairs(&mut self, ar: &mut Arena<F>, s: &mut Schur, x: ExprId, q: Tid, lq: ExprId, eq: Option<MapId>, sdeg: i32, mult: usize, b: &Budget) -> Result<(Pairs<F>, SearchStats), Vec<Missing>> {
        self.trace.clear();
        self.bind(ar);
        self.failed = None;
        let strat = Some(s.table(q).fmt());
        let r = run_call(
            ar,
            b,
            |ar, ctx| self.pairs_inner(ar, s, x, q, lq, eq, sdeg, mult, ctx),
            |p: &Pairs<F>| {
                let mut v = p.iotas.clone();
                v.extend(&p.pis);
                v
            },
            |p, k| Pairs { iotas: k[..p.iotas.len()].to_vec(), pis: k[p.iotas.len()..].to_vec(), form: p.form },
        );
        r.map_err(|((k, e), st)| {
            let m = miss(&format!("Search pairs: {} ({} candidates)", e, st.candidates), sdeg, strat, 0, mult, k);
            self.failed = Some(st);
            m
        })
    }
    fn closure(&mut self, ar: &mut Arena<F>, s: &mut Schur, lp: ExprId, p: Tid, ep: Option<MapId>, outs: bool, b: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        self.trace.clear();
        self.bind(ar);
        self.failed = None;
        let strat = Some(s.table(p).fmt());
        let n = cdiag_degrees(s, p).len();
        let r = run_call(ar, b, |ar, ctx| self.closure_inner(ar, s, lp, p, ep, outs, ctx, 0), |v| v.clone(), |_, k| k);
        r.map_err(|((k, e), st)| {
            let m = miss(&format!("Search closure ({}): {} ({} candidates)", if outs { "outs" } else { "ins" }, e, st.candidates), 0, strat, 0, n, k);
            self.failed = Some(st);
            m
        })
    }
    fn end0(&mut self, ar: &mut Arena<F>, s: &mut Schur, l: ExprId, b: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        self.full(ar, s, l, l, 0, b)
    }
    fn failed_stats(&self) -> Option<SearchStats> {
        self.failed.clone()
    }
}

// ------------------------------------------------------------------ Chain

/// first Ok wins; Missing lists merged; counters of every attempted provider accumulated
pub struct Chain<F>(pub Vec<Box<dyn HomProvider<F>>>);

impl<F: Field> Chain<F> {
    /// §4.5 order: max label <= 2: [Exact, Search S0-S6]; max label >= 3: [Search S0-S2, Exact
    /// (cost table only), Search S3-S6]
    pub fn for_labels(max_label: u8, exact: Exact) -> Chain<F> {
        let mut s_all = Search::tiers(0, 6);
        s_all.costs = exact.costs.clone();
        if max_label <= 2 {
            return Chain(vec![Box::new(exact), Box::new(s_all)]);
        }
        let mut s_lo = Search::tiers(0, 2);
        s_lo.costs = exact.costs.clone();
        let mut s_hi = Search::tiers(3, 6);
        s_hi.costs = exact.costs.clone();
        let mut ex = exact;
        ex.max_label = 0; // cost table only
        Chain(vec![Box::new(s_lo), Box::new(ex), Box::new(s_hi)])
    }
}

macro_rules! chain_call {
    ($self:ident, $call:expr) => {{
        let mut all: Vec<Missing> = vec![];
        let mut acc = SearchStats::default();
        for p in $self.0.iter_mut() {
            let name = p.name();
            match $call(p) {
                Ok((v, st)) => {
                    let mut st2 = acc.clone();
                    st2.absorb(&st);
                    st2.detect = if st.detect.is_empty() { acc.detect } else { st.detect };
                    return Ok((v, st2));
                }
                Err(ms) => {
                    if let Some(st) = p.failed_stats() {
                        acc.absorb(&st);
                    }
                    for mut m in ms {
                        m.what = format!("[{}] {}", name, m.what);
                        all.push(m);
                    }
                }
            }
        }
        let _ = acc;
        Err(all)
    }};
}

impl<F: Field> HomProvider<F> for Chain<F> {
    fn name(&self) -> &'static str {
        "Chain"
    }
    fn full(&mut self, ar: &mut Arena<F>, s: &mut Schur, src: ExprId, tgt: ExprId, d: i32, b: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        chain_call!(self, |p: &mut Box<dyn HomProvider<F>>| p.full(ar, s, src, tgt, d, b))
    }
    fn pairs(&mut self, ar: &mut Arena<F>, s: &mut Schur, x: ExprId, q: Tid, lq: ExprId, eq: Option<MapId>, sdeg: i32, mult: usize, b: &Budget) -> Result<(Pairs<F>, SearchStats), Vec<Missing>> {
        chain_call!(self, |p: &mut Box<dyn HomProvider<F>>| p.pairs(ar, s, x, q, lq, eq, sdeg, mult, b))
    }
    fn closure(&mut self, ar: &mut Arena<F>, s: &mut Schur, lp: ExprId, p: Tid, ep: Option<MapId>, outs: bool, b: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        chain_call!(self, |pr: &mut Box<dyn HomProvider<F>>| pr.closure(ar, s, lp, p, ep, outs, b))
    }
    fn end0(&mut self, ar: &mut Arena<F>, s: &mut Schur, l: ExprId, b: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        chain_call!(self, |p: &mut Box<dyn HomProvider<F>>| p.end0(ar, s, l, b))
    }
    fn failed_stats(&self) -> Option<SearchStats> {
        let mut acc = SearchStats::default();
        let mut any = false;
        for p in &self.0 {
            if let Some(st) = p.failed_stats() {
                acc.absorb(&st);
                any = true;
            }
        }
        if any { Some(acc) } else { None }
    }
}
