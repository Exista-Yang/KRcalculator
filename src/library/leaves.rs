//! src/library/leaves.rs — M4c light-leaves provider. Milestone A1: Stage 0 only (M4C_PLAN §3.2, D2 §2.1,
//! C_m4c "do it first").
//!
//! **Stage 0 (product closure).** Let p be a (λ, λ) type whose canonical ladder L_p is,
//! after adjacent commutations of rungs (`homprov::rungs_commute`), a product A·B of two words
//! with top λ that are indecomposable of (λ, λ) types a, b with c_diag(b_p) = c_diag(b_a)·c_diag(b_b)
//! (L3 §3.2). The closure generators of p are then
//! - ins (R -> L_p):  iso ∘ Comp(Kron(α_i, B), β_j)   with α_i: R -> A, β_j: R -> B;
//! - outs (L_p -> R): Comp(β'_j, Kron(α'_i, B)) ∘ iso  with α'_i: A -> R, β'_j: B -> R;
//!
//! where α, β are the closure generators of a and b from the inner provider (the production
//! chain, which certifies those types at level 4/5), and iso is the chain of degree-0
//! commutation maps L_p ⇄ A·B, each an exact `hom_space` template on the 2-rung window (one-
//! dimensional; genuine by construction). Degrees add, so the degree multiset is
//! c_diag(a)·c_diag(b) = c_diag(p).
//!
//! Nothing is accepted on this argument: the result is accepted only through
//! `homprov::certify_closure` (ends, degree multiset = c_diag(b_p), identity-coset matrix of rank
//! n at two points, genuineness), and sindec re-certifies it after every layer-A load.
use crate::algebra::field::Field;
use crate::library::homprov::*;
use crate::combinat::schur::{Lp, Schur, Table, Tid};
use crate::colored::sexpr::Expr;
use crate::colored::smap::{Arena, ExprId, MapId, MapNode};
use crate::colored::tmpl::Window;
use crate::combinat::web::Rung;
use std::collections::{HashMap, VecDeque};
use std::time::Instant;

/// provider stamp suffix of Stage 0 (layer-A headers)
pub const STAGE0_STAMP: &str = "stage0=ProductClosure(v1,inner-first)";
/// stamp suffix with `product_first`
pub const STAGE0_STAMP_PF: &str = "stage0=ProductClosure(v1,product-first)";

/// A product factorization of a closure word (Stage 0).
#[derive(Clone, Debug)]
pub struct Factorization {
    /// the commuted word A·B
    pub word: Vec<Rung>,
    /// |A|
    pub split: usize,
    pub a: Tid,
    pub b: Tid,
    /// positions i of the adjacent swaps (i, i+1) taking the original word to `word`, in order
    pub swaps: Vec<usize>,
    /// A and B are the canonical ladders of a and b
    pub canonical: bool,
}

fn indecomposable_type(s: &mut Schur, w: &[Rung]) -> Option<Tid> {
    let ch = s.ladder_char(w);
    let ex = s.expand(&ch);
    if ex.len() == 1 && ex[0].1 == Lp::one() { Some(ex[0].0) } else { None }
}

/// Commutation class of w with the swap path to every word (breadth first, at most `cap` words).
fn class_with_paths(w: &[Rung], cap: usize) -> Vec<(Vec<Rung>, Vec<usize>)> {
    let mut seen: HashMap<Vec<Rung>, usize> = HashMap::new();
    let mut out: Vec<(Vec<Rung>, Vec<usize>)> = vec![(w.to_vec(), vec![])];
    seen.insert(w.to_vec(), 0);
    let mut q = VecDeque::from([0usize]);
    while let Some(k) = q.pop_front() {
        if out.len() >= cap {
            break;
        }
        let (cur, path) = out[k].clone();
        for i in 0..cur.len().saturating_sub(1) {
            if rungs_commute(&cur[i], &cur[i + 1]) {
                let mut nw = cur.clone();
                nw.swap(i, i + 1);
                if !seen.contains_key(&nw) {
                    let mut p2 = path.clone();
                    p2.push(i);
                    seen.insert(nw.clone(), out.len());
                    out.push((nw, p2));
                    q.push_back(out.len() - 1);
                }
            }
        }
    }
    out
}

/// Stage 0 test: a factorization of the word w of the (λ, λ) type p as A·B (after commutations)
/// with indecomposable non-identity factors of types a, b and c_diag(p) = c_diag(a)·c_diag(b).
/// Words with canonical factors are preferred, then fewer swaps.
pub fn product_factorization(s: &mut Schur, w: &[Rung], p: Tid, cap: usize) -> Option<Factorization> {
    let lam = s.lam.clone();
    let id = s.intern(Table::identity(&lam));
    let cdp = s.c_diag(p);
    let mut best: Option<Factorization> = None;
    for (word, swaps) in class_with_paths(w, cap) {
        let mut lv = lam.clone();
        for k in 1..word.len() {
            lv = crate::combinat::dcoset::apply_labels(&lv, &word[k - 1]).expect("product_factorization: invalid word");
            if lv != lam {
                continue;
            }
            let (aw, bw) = (&word[..k], &word[k..]);
            let (Some(a), Some(b)) = (indecomposable_type(s, aw), indecomposable_type(s, bw)) else { continue };
            if a == id || b == id {
                continue;
            }
            if s.c_diag(a).mul(&s.c_diag(b)) != cdp {
                continue;
            }
            let canonical = crate::library::sindec::canonical_ladder(s, a) == aw && crate::library::sindec::canonical_ladder(s, b) == bw;
            let f = Factorization { word: word.clone(), split: k, a, b, swaps: swaps.clone(), canonical };
            let better = match &best {
                None => true,
                Some(o) => (f.canonical && !o.canonical) || (f.canonical == o.canonical && f.swaps.len() < o.swaps.len()),
            };
            if better {
                best = Some(f);
            }
        }
    }
    best
}

/// Stage 0 provider: the inner provider first (Leaves last, C_m4c M4); for a closure call that
/// the inner provider does not certify and whose word is a product, the product closure.
/// `full`, `pairs` and `end0` go to the inner provider unchanged.
pub struct ProductClosure<F: Field> {
    pub inner: Box<dyn HomProvider<F>>,
    /// try the product rule before the inner provider
    pub product_first: bool,
    /// unknown cap of the 2-rung commutation windows
    pub comm_cap: u128,
    pub seed: u64,
    last_fail: Option<SearchStats>,
}

impl<F: Field> ProductClosure<F> {
    pub fn new(inner: Box<dyn HomProvider<F>>) -> Self {
        ProductClosure { inner, product_first: false, comm_cap: 200_000, seed: 0x57a6e0, last_fail: None }
    }

    fn miss(what: String, kind: MissKind) -> Vec<Missing> {
        vec![Missing { what, degree: 0, stratum: None, got: 0, want: 1, kind }]
    }

    /// one commutation map: the Local node of the degree-0 template on the 2-rung window of
    /// `from` at position i, from -> swapped (from and to interned)
    fn swap_map(&self, ar: &mut Arena<F>, from: &[Rung], i: usize, forward: bool) -> Result<(MapId, Vec<Rung>), Vec<Missing>> {
        let lam = ar.lam().to_vec();
        let (r1, r2) = (from[i], from[i + 1]);
        let lo = r1.i.min(r2.i);
        let hi = r1.i.max(r2.i) + 1;
        let labels = level_of(&lam, from, i)[lo..=hi].to_vec();
        let sh = |r: Rung| Rung { i: r.i - lo, dir: r.dir, r: r.r };
        let rel = vec![sh(r1), sh(r2)];
        let sw = vec![sh(r2), sh(r1)];
        let mut to = from.to_vec();
        to.swap(i, i + 1);
        let (a, b, wa, wb) = if forward { (from.to_vec(), to.clone(), rel, sw) } else { (to.clone(), from.to_vec(), sw, rel) };
        let u = hom_unknowns(&labels, &wa, &wb, 0);
        if u > self.comm_cap {
            return Err(Self::miss(format!("stage0: commutation window {:?} [{}] has {} unknowns > cap {}", labels, wstr(&wa), u, self.comm_cap), MissKind::Refused));
        }
        let ts = ar.lib().try_hom_basis_cap(&Window::new(&labels, &wa), &Window::new(&labels, &wb), 0, self.comm_cap).map_err(|e| Self::miss(format!("stage0: commutation window: {}", e), MissKind::Refused))?;
        if ts.len() != 1 {
            return Err(Self::miss(format!("stage0: commutation window {:?} [{}] -> [{}]: Hom^0 has dimension {}", labels, wstr(&wa), wstr(&wb), ts.len()), MissKind::Exhausted));
        }
        let ea = ar.expr(Expr::new(&lam, &a).map_err(|e| Self::miss(e, MissKind::Exhausted))?);
        let eb = ar.expr(Expr::new(&lam, &b).map_err(|e| Self::miss(e, MissKind::Exhausted))?);
        let m = ar.add(MapNode::Local { src: ea, tgt: eb, at: i, first: lo, tmpl: ts[0] }, "stage0 commute").map_err(|e| Self::miss(e, MissKind::Exhausted))?;
        Ok((m, to))
    }

    /// The product closure (Stage 0) of the word of `lp`, certified.
    pub fn product_closure(&mut self, ar: &mut Arena<F>, s: &mut Schur, lp: ExprId, p: Tid, outs: bool, b: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        let t0 = Instant::now();
        let w = ar.expr_of(lp).rungs.clone();
        let Some(fac) = product_factorization(s, &w, p, 256) else {
            return Err(Self::miss(format!("stage0: {} is not a product type", s.table(p).fmt()), MissKind::Refused));
        };
        let mark = ar.mark();
        match self.product_inner(ar, s, lp, p, outs, b, &fac, t0) {
            Ok(r) => Ok(r),
            Err(e) => {
                ar.commit(mark, &[]);
                Err(e)
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn product_inner(&mut self, ar: &mut Arena<F>, s: &mut Schur, lp: ExprId, p: Tid, outs: bool, b: &Budget, fac: &Factorization, t0: Instant) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        let lam = ar.lam().to_vec();
        let w = ar.expr_of(lp).rungs.clone();
        let mut st = SearchStats::default();
        let mark = ar.mark();
        let (aw, bw) = (fac.word[..fac.split].to_vec(), fac.word[fac.split..].to_vec());
        let ea = ar.expr(Expr::new(&lam, &aw).map_err(|e| Self::miss(e, MissKind::Exhausted))?);
        let eb = ar.expr(Expr::new(&lam, &bw).map_err(|e| Self::miss(e, MissKind::Exhausted))?);
        let (fa, sa) = self.inner.closure(ar, s, ea, fac.a, None, outs, b)?;
        st.absorb(&sa);
        let (fb, sb) = self.inner.closure(ar, s, eb, fac.b, None, outs, b)?;
        st.absorb(&sb);
        // commutation chain between w and A·B: outs needs w -> A·B, ins needs A·B -> w
        let mut iso: Option<MapId> = None;
        let mut cur = w.clone();
        for &i in &fac.swaps {
            let (m, next) = self.swap_map(ar, &cur, i, outs)?;
            st.window_calls += 1;
            iso = Some(match iso {
                None => m,
                // outs: maps compose forward (w -> ... -> A·B); ins: backward (A·B -> ... -> w)
                Some(prev) => {
                    let n = if outs { MapNode::Comp(m, prev) } else { MapNode::Comp(prev, m) };
                    ar.add(n, "stage0 iso").map_err(|e| Self::miss(e, MissKind::Exhausted))?
                }
            });
            cur = next;
        }
        assert_eq!(cur, fac.word, "stage0: swap path does not reach A·B");
        let mut maps = vec![];
        let add = |ar: &mut Arena<F>, n: MapNode<F>| ar.add(n, "stage0 product").map_err(|e| Self::miss(e, MissKind::Exhausted));
        for &ai in &fa {
            let k = add(ar, MapNode::Kron(ai, bw.clone()))?;
            for &bj in &fb {
                let m = if outs { add(ar, MapNode::Comp(bj, k))? } else { add(ar, MapNode::Comp(k, bj))? };
                let m = match iso {
                    None => m,
                    Some(c) => {
                        if outs {
                            add(ar, MapNode::Comp(m, c))?
                        } else {
                            add(ar, MapNode::Comp(c, m))?
                        }
                    }
                };
                maps.push(m);
            }
        }
        let kept = ar.commit(mark, &maps);
        let r0 = ar.expr(Expr::new(&lam, &[]).map_err(|e| Self::miss(e, MissKind::Exhausted))?);
        let want = cdiag_degrees(s, p);
        certify_closure(ar, &lam, r0, lp, None, &kept, &want, outs, self.seed).map_err(|e| Self::miss(format!("stage0 certify_closure {}: {}", s.table(p).fmt(), e), MissKind::Exhausted))?;
        st.candidates += kept.len();
        st.evaluated += kept.len();
        st.secs = t0.elapsed().as_secs_f64();
        st.detect = "stage0-product";
        Ok((kept, st))
    }
}

/// c_diag(b_p) as a sorted degree list (as homprov's private helper)
pub fn cdiag_degrees(s: &mut Schur, p: Tid) -> Vec<i32> {
    let cd = s.c_diag(p);
    let mut want = vec![];
    for (i, &c) in cd.c.iter().enumerate() {
        for _ in 0..c {
            want.push(cd.lo + i as i32);
        }
    }
    want
}

fn wstr(w: &[Rung]) -> String {
    crate::library::smapio::rungs_str(w)
}

impl<F: Field> HomProvider<F> for ProductClosure<F> {
    fn name(&self) -> &'static str {
        "Stage0"
    }
    fn full(&mut self, ar: &mut Arena<F>, s: &mut Schur, src: ExprId, tgt: ExprId, d: i32, b: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        self.last_fail = None;
        let r = self.inner.full(ar, s, src, tgt, d, b);
        if r.is_err() {
            self.last_fail = self.inner.failed_stats();
        }
        r
    }
    fn pairs(&mut self, ar: &mut Arena<F>, s: &mut Schur, x: ExprId, q: Tid, lq: ExprId, eq: Option<MapId>, sdeg: i32, mult: usize, b: &Budget) -> Result<(Pairs<F>, SearchStats), Vec<Missing>> {
        self.last_fail = None;
        let r = self.inner.pairs(ar, s, x, q, lq, eq, sdeg, mult, b);
        if r.is_err() {
            self.last_fail = self.inner.failed_stats();
        }
        r
    }
    fn closure(&mut self, ar: &mut Arena<F>, s: &mut Schur, lp: ExprId, p: Tid, ep: Option<MapId>, outs: bool, b: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        self.last_fail = None;
        if ep.is_some() {
            let r = self.inner.closure(ar, s, lp, p, ep, outs, b);
            if r.is_err() {
                self.last_fail = self.inner.failed_stats();
            }
            return r;
        }
        let w = ar.expr_of(lp).rungs.clone();
        let is_prod = product_factorization(s, &w, p, 256).is_some();
        if self.product_first && is_prod {
            if let Ok(r) = self.product_closure(ar, s, lp, p, outs, b) {
                return Ok(r);
            }
        }
        match self.inner.closure(ar, s, lp, p, ep, outs, b) {
            Ok(r) => Ok(r),
            Err(mut e) => {
                let inner_fail = self.inner.failed_stats();
                if !is_prod || self.product_first {
                    self.last_fail = inner_fail;
                    return Err(e);
                }
                match self.product_closure(ar, s, lp, p, outs, b) {
                    Ok(r) => Ok(r),
                    Err(e2) => {
                        e.extend(e2);
                        self.last_fail = inner_fail;
                        Err(e)
                    }
                }
            }
        }
    }
    fn end0(&mut self, ar: &mut Arena<F>, s: &mut Schur, l: ExprId, b: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        self.last_fail = None;
        let r = self.inner.end0(ar, s, l, b);
        if r.is_err() {
            self.last_fail = self.inner.failed_stats();
        }
        r
    }
    fn failed_stats(&self) -> Option<SearchStats> {
        self.last_fail.clone()
    }
}

// =============================================================================================
// M4c A2: localized light leaves (Lemma A, D1 §1.3; M4C_PLAN §3.1-§3.2)
// =============================================================================================
//
// **Lemma A.** X = X'·ρ. For every junction type t' at the top level of X' a family
// LL_{t'} ⊂ Hom(X', W_{t'}) (Out) with an invertible block at rep(t'), and for every step t' -> t
// a family S_{t'->t} ⊂ Hom(W_{t'}·ρ, W_t) with degree multiset G_{t'->t} and an invertible
// top-to-top matrix at rep(t) (rows = maps, columns = the flows of W_{t'}·ρ at rep(t) whose
// configuration after |W_{t'}| rungs has type t'). Then LL_t = { S ∘ Kron(ℓ, [ρ]) } has an
// invertible block at rep(t) and the degrees add. In: the dual families (columns), maps
// W_t -> X built as Kron(ℓ, [ρ]) ∘ S with S ⊂ Hom(W_t, W_{t'}·ρ).
//
// **Step families** (D1 §2.3, M4C_PLAN §3.2), in order: `Id` when W_t = W_{t'}·ρ; then exact
// `hom_space` bases on the differing window (common prefix and suffix stripped: id_P ⊗ s ⊗ id_Q,
// not a priori complete, A1 report §5 (a)); then exact bases of the whole words. Selection is
// greedy by ascending degree over every degree of supp G (C_m4c m4) on the top-to-top columns;
// a family is accepted only with degree multiset exactly G and full top-to-top rank.
// The budget rule (M4C_PLAN §0.3): no exact call above `budget` (10⁵) unknowns; up to `hard_cap`
// (4·10⁵) only with an uncensored exact-shape cost record. Otherwise the step is Refused.
//
// Internal rank tests are diagnostics; acceptance is `certify_closure` only.

use crate::combinat::coset::CosetSpace;
use crate::combinat::dcoset::{junction_types, junction_word, StepCache};
use crate::colored::smap::{points, Eval};

/// direction of a leaf family: Out = maps X -> W_t, In = maps W_t -> X
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum LeafDir {
    Out,
    In,
}

/// how a step family was realized
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum StepKind {
    Id,
    Elem,
    Window,
    Whole,
}

/// one leaf family at the top level of X: the junction type, its junction word, the maps
#[derive(Clone, Debug)]
pub struct LeafFamily {
    pub t: Tid,
    pub word: Vec<Rung>,
    pub wid: ExprId,
    pub maps: Vec<MapId>,
}

/// record of one step family (for reports and the census comparison)
#[derive(Clone, Debug)]
pub struct StepRec {
    pub level: usize,
    pub from: Tid,
    pub to: Tid,
    pub kind: StepKind,
    pub degs: Vec<i32>,
    /// max unknowns over the exact calls made for this step (0 for Id)
    pub unknowns: u128,
    pub secs: f64,
}

/// one leaves-stage failure
#[derive(Clone, Debug)]
pub struct LeafFail {
    pub kind: MissKind,
    pub what: String,
}
fn lfail(kind: MissKind, what: String) -> LeafFail {
    LeafFail { kind, what }
}

/// cost counters of a Leaves call
#[derive(Clone, Debug, Default)]
pub struct LeafStats {
    pub steps: Vec<StepRec>,
    pub exact_calls: usize,
    pub exact_secs: f64,
    pub max_unknowns: u128,
    pub max_word_flows: u128,
    pub maps_built: usize,
}

/// the M4c Leaves provider (A2: closure only; `full`, `pairs`, `end0` refuse)
pub struct Leaves {
    /// no exact call above this many unknowns (M4C_PLAN §0.3)
    pub budget: u128,
    /// calls in (budget, hard_cap] only with an uncensored exact-shape cost record
    pub hard_cap: u128,
    /// one-sided calls (an empty end) on windows with a label > 3 above this many unknowns only
    /// with a cost record
    pub one_sided_cap: u128,
    pub costs: Option<CostTable>,
    /// try the differing window before the whole words
    pub use_window: bool,
    /// E2/E3 closed forms for one-dimensional cup/cap steps (validated by `validate_elem`). Off by
    /// default: Elem constants are recorded per arena and are not part of layer A, so a replay
    /// (or another session) would reject the loaded maps as not genuine (A2 finding)
    pub use_elem: bool,
    /// BFS limits of `junction_word`
    pub bfs_len: usize,
    pub bfs_cap: usize,
    /// alternative reduced junction words (≤ alt_len rungs, ≤ alt_cap words visited) when a step
    /// into the default word exceeds the budget; 0 = off
    pub alt_len: usize,
    pub alt_cap: usize,
    pub seed: u64,
    /// counters of the last call
    pub last: LeafStats,
    /// a log line per step when set (also with KRHOM_TRACE)
    pub verbose: bool,
    /// C2-7 joint junction words (`joint_descent`; off by default; stamp suffix `+joint` only when on)
    pub joint: bool,
    /// suffix length, suffixes visited per prefix, candidate pairs fully evaluated, sweeps
    pub joint_suffix: usize,
    pub joint_cap: usize,
    pub joint_keep: usize,
    pub joint_sweeps: usize,
    /// C2-7: Leaves pairs also for targets with an idempotent e_q (`pairs_leaves_eq`)
    pub eq_pairs: bool,
    /// C2-7: per-step target words with rex-chain detours at every level (`step_words`); with
    /// `joint` only; off = final level only
    pub detours: bool,
    /// chains planned by `step_words` ((type, from word, to word) -> chain), reused by `build_leaves`
    chain_memo: std::cell::RefCell<HashMap<(Tid, Vec<Rung>, Vec<Rung>), crate::library::rexplan::Chain>>,
    last_fail: Option<SearchStats>,
}

/// provider stamp suffix of the Leaves backend (layer-A headers)
pub const LEAVES_STAMP: &str = "leaves=Leaves(v1,closure,threeup,prefix-ext+costalt5,window-whole,budget=1e5/4e5/os2e4l4)";

impl Default for Leaves {
    fn default() -> Self {
        Leaves { budget: 100_000, hard_cap: 400_000, one_sided_cap: 20_000, costs: None, use_window: true, use_elem: false, bfs_len: 6, bfs_cap: 20_000, alt_len: 5, alt_cap: 400_000, seed: 0x1ea7e5, last: LeafStats::default(), verbose: false, joint: false, joint_suffix: 4, joint_cap: 60_000, joint_keep: 40, joint_sweeps: 3, eq_pairs: false, detours: false, chain_memo: Default::default(), last_fail: None }
    }
}

/// C2-7: the words P·A reduced for t with A of at most `max_len` rungs (breadth first from the
/// top labels of P, at most `cap` suffixes visited), by ascending flows.
pub fn suffix_words(s: &mut Schur, p: &[Rung], t: Tid, max_len: usize, cap: usize) -> Vec<Vec<Rung>> {
    use crate::combinat::web::Dir;
    let lam = s.lam.clone();
    if !crate::combinat::dcoset::valid_word(&lam, p) {
        return vec![];
    }
    let mu = s.table(t).cols();
    let n = lam.len();
    let mut q: VecDeque<Vec<Rung>> = VecDeque::new();
    q.push_back(p.to_vec());
    let mut seen = 0usize;
    let mut out = vec![];
    while let Some(w) = q.pop_front() {
        seen += 1;
        if seen > cap {
            break;
        }
        let top = crate::combinat::dcoset::top_labels(&lam, &w);
        if top == mu && crate::combinat::dcoset::reduced_for(s, &w, t).ok() {
            out.push(w.clone());
        }
        if w.len() >= p.len() + max_len {
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
    out.sort_by_key(|w| crate::combinat::schur::flows(&lam, w));
    out
}

/// C2-7: the local type of a word (W1): the table t of char w with w reduced for t
/// (`reduced_for`: coefficient 1, support ≤ t); ∅ gives the identity type.
pub fn local_type(s: &mut Schur, w: &[Rung]) -> Option<Tid> {
    if w.is_empty() {
        let lam = s.lam.clone();
        return Some(s.intern(Table::identity(&lam)));
    }
    let ch = s.ladder_char(w);
    let cands: Vec<Tid> = ch.t.iter().filter(|(_, c)| *c == Lp::one()).map(|x| x.0).collect();
    cands.into_iter().find(|&t| crate::combinat::dcoset::reduced_for(s, w, t).ok())
}

/// C2-7 Prop W checks W1 + W2 (D_grass §1.3) of the window of the step W_{t'}·ρ -> W_t: the
/// window has no common suffix (it ends in ρ: A'·ρ_loc -> B'), A' and B' are reduced for local
/// types a', b on the window labels ν (W1), and b occurs in step_gf(a', ρ_loc) with
/// G_{a'->b} = G (W2). W3/W4 (flow-level) are not checked a priori: the runtime rank test of
/// `step_family` decides (a failure is Exhausted, never a wrong basis). Planning only.
pub fn window_w12(local: &mut HashMap<Vec<u8>, Schur>, lam: &[u8], wfrom: &[Rung], rho: &Rung, wto: &[Rung], g: &Lp) -> bool {
    let mut a = wfrom.to_vec();
    a.push(*rho);
    let ws = win_split(lam, &a, wto);
    if ws.suffix > 0 || ws.a.is_empty() || ws.labels.is_empty() {
        return false;
    }
    let s = local.entry(ws.labels.clone()).or_insert_with(|| Schur::new(&ws.labels));
    let ap = &ws.a[..ws.a.len() - 1];
    let rl = ws.a[ws.a.len() - 1];
    let (Some(ta), Some(tb)) = (local_type(s, ap), local_type(s, &ws.b)) else { return false };
    let mut sc = StepCache::new();
    sc.step(s, ta, &rl).into_iter().any(|(t, gl)| t == tb && gl == *g)
}

/// C2-7 W-aware step score (joint planning): the window score only when the window passes
/// W1 + W2 (`window_w12`), else the whole-word score (what `step_family` falls back to when the
/// window family has too small a rank).
pub fn step_score_w(local: &mut HashMap<Vec<u8>, Schur>, lam: &[u8], wfrom: &[Rung], rho: &Rung, wto: &[Rung], g: &Lp) -> u128 {
    let mut a = wfrom.to_vec();
    a.push(*rho);
    if a == wto {
        return 0;
    }
    if crate::combinat::schur::flows(lam, wto) == 0 || crate::combinat::dcoset::top_labels(lam, &a) != crate::combinat::dcoset::top_labels(lam, wto) {
        return u128::MAX;
    }
    let ws = win_split(lam, &a, wto);
    let proper = ws.at > 0 || ws.suffix > 0;
    let wok = proper && window_w12(local, lam, wfrom, rho, wto, g);
    let mut mx = 0u128;
    for (d, _) in degs_of(g) {
        let u = if wok { hom_unknowns(&ws.labels, &ws.a, &ws.b, d).max(hom_unknowns(&ws.labels, &ws.b, &ws.a, d)) } else { hom_unknowns(lam, &a, wto, d).max(hom_unknowns(lam, wto, &a, d)) };
        mx = mx.max(u);
    }
    mx
}

/// One step of a Lemma A plan (C2-7 census: every step of the executor's own plan).
#[derive(Clone, Debug)]
pub struct PlanStep {
    pub level: usize,
    pub from: Tid,
    pub to: Tid,
    pub wfrom: Vec<Rung>,
    pub rho: Rung,
    pub wto: Vec<Rung>,
    pub g: Lp,
    /// `step_score` (what the executor's budget sees, both directions, max over supp G)
    pub score: u128,
    /// window of W_{t'}·ρ vs W_t (common prefix length, common suffix length, labels)
    pub win_at: usize,
    pub win_suffix: usize,
    pub win_labels: Vec<u8>,
    /// per degree of supp G: (d, whole max of both directions, window max of both directions or
    /// u128::MAX when there is no proper window)
    pub degs: Vec<(i32, u128, u128)>,
    /// W1 + W2 of the window (`window_w12`) and the W-aware score (`step_score_w`)
    pub w12: bool,
    pub score_w: u128,
}

/// Every step (level k, every junction type t' at k−1 with a step to every t at k) of a plan.
pub fn plan_steps(s: &mut Schur, sc: &mut StepCache, w: &[Rung], plan: &[Vec<(Tid, Vec<Rung>)>]) -> Vec<PlanStep> {
    let lam = s.lam.clone();
    let mut out = vec![];
    let mut wl: HashMap<Vec<u8>, Schur> = HashMap::new();
    for k in 1..=w.len() {
        for (t, wt) in &plan[k] {
            for (tp, wp) in &plan[k - 1] {
                let Some(g) = sc.step(s, *tp, &w[k - 1]).into_iter().find(|e| e.0 == *t).map(|e| e.1) else { continue };
                let score = step_score(&lam, wp, &w[k - 1], wt, &g);
                let mut a = wp.clone();
                a.push(w[k - 1]);
                let (mut at, mut suffix, mut labels, mut degs) = (0, 0, vec![], vec![]);
                if a != *wt && crate::combinat::schur::flows(&lam, wt) > 0 && crate::combinat::dcoset::top_labels(&lam, &a) == crate::combinat::dcoset::top_labels(&lam, wt) {
                    let ws = win_split(&lam, &a, wt);
                    at = ws.at;
                    suffix = ws.suffix;
                    labels = ws.labels.clone();
                    for (d, _) in degs_of(&g) {
                        let whole = hom_unknowns(&lam, &a, wt, d).max(hom_unknowns(&lam, wt, &a, d));
                        let win = if ws.at > 0 || ws.suffix > 0 { hom_unknowns(&ws.labels, &ws.a, &ws.b, d).max(hom_unknowns(&ws.labels, &ws.b, &ws.a, d)) } else { u128::MAX };
                        degs.push((d, whole, win));
                    }
                }
                let w12 = window_w12(&mut wl, &lam, wp, &w[k - 1], wt, &g);
                let score_w = step_score_w(&mut wl, &lam, wp, &w[k - 1], wt, &g);
                out.push(PlanStep { level: k, from: *tp, to: *t, wfrom: wp.clone(), rho: w[k - 1], wto: wt.clone(), g, score, win_at: at, win_suffix: suffix, win_labels: labels, degs, w12, score_w });
            }
        }
    }
    out
}

/// identity coset of the block (λ, λ) (as homprov's private helper `identity_coset`)
pub fn identity_coset(lam: &[u8]) -> u64 {
    let cs = CosetSpace::new(lam, lam);
    cs.rank(&cs.roots.bottom())
}

/// coset rank of rep(t) of an interned table in its block (as homprov's private `rep_rank`)
pub fn rep_rank_of(s: &Schur, lam: &[u8], t: Tid) -> u64 {
    let tab = s.table(t).clone();
    let cs = CosetSpace::new(lam, &tab.cols());
    cs.rank(&cs.rep(&tab))
}

fn degs_of(g: &Lp) -> Vec<(i32, usize)> {
    (0..g.c.len()).filter(|&i| g.c[i] != 0).map(|i| (g.lo + i as i32, g.c[i] as usize)).collect()
}

/// Estimated unknowns of the step family W_{t'}·ρ -> W_t (max over supp G and both directions;
/// the differing window when the words share a prefix or suffix, else the whole words; 0 for Id)
pub fn step_score(lam: &[u8], wfrom: &[Rung], rho: &Rung, wto: &[Rung], g: &Lp) -> u128 {
    let mut a = wfrom.to_vec();
    a.push(*rho);
    if a == wto {
        return 0;
    }
    if crate::combinat::schur::flows(lam, wto) == 0 || crate::combinat::dcoset::top_labels(lam, &a) != crate::combinat::dcoset::top_labels(lam, wto) {
        return u128::MAX;
    }
    let ws = win_split(lam, &a, wto);
    let mut mx = 0u128;
    for (d, _) in degs_of(g) {
        let whole = hom_unknowns(lam, &a, wto, d).max(hom_unknowns(lam, wto, &a, d));
        let u = if ws.at > 0 || ws.suffix > 0 { whole.min(hom_unknowns(&ws.labels, &ws.a, &ws.b, d).max(hom_unknowns(&ws.labels, &ws.b, &ws.a, d))) } else { whole };
        mx = mx.max(u);
    }
    mx
}

/// The E2/E3 closed form (smap `ElemKind`) of a 2-upright step window, if any: Out: window
/// a -> ∅ (Eps, EpsL, Trace, TraceL); In: ∅ -> a (Eta, EtaL, Iota, IotaL). Returns the kind and
/// its constant key (E2: (a, b) labels; E3: (c, rung size)).
pub fn elem_pattern(labels: &[u8], a: &[Rung], b: &[Rung], dir: LeafDir) -> Option<(crate::colored::smap::ElemKind, (u8, u8))> {
    use crate::colored::smap::ElemKind::*;
    use crate::combinat::web::Dir;
    if labels.len() != 2 || !b.is_empty() || a.len() != 2 || a[0].i != 0 || a[1].i != 0 || a[0].r != a[1].r || a[0].dir == a[1].dir {
        return None;
    }
    let k = a[0].r as u8;
    let (l0, l1) = (labels[0], labels[1]);
    let out = dir == LeafDir::Out;
    if a[0].dir == Dir::E {
        // [E0.k, F0.k]
        if l0 == 0 {
            return Some((if out { TraceL } else { IotaL }, (l1, k)));
        }
        if k == l1 {
            return Some((if out { Eps } else { Eta }, (l0, l1)));
        }
    } else {
        // [F0.k, E0.k]
        if l1 == 0 {
            return Some((if out { Trace } else { Iota }, (l0, k)));
        }
        if k == l0 {
            return Some((if out { EpsL } else { EtaL }, (l0, l1)));
        }
    }
    None
}

/// `validate_elem` (M4C_PLAN §3.4, C_m4c M5): the closed form `kind` on the 2-upright window
/// (labels, a) with scalar 1 is a nonzero multiple of the exact `hom_space` map, which must be the
/// unique one up to scalars (Hom^deg one-dimensional). Compared blockwise at every coset of the
/// window at two points (exact comparison; refused above `cap` unknowns). Ok(ratio) on success.
pub fn validate_elem<F: Field>(kind: crate::colored::smap::ElemKind, labels: &[u8], a: &[Rung], cap: u128) -> Result<F, String> {
    use crate::colored::smap::ElemKind::*;
    let out = matches!(kind, Eps | EpsL | Trace | TraceL);
    let mut ar = Arena::<F>::new(labels, crate::colored::tmpl::TmplLib::new(None));
    let ew = ar.expr(Expr::new(labels, a)?);
    let e0 = ar.expr(Expr::new(labels, &[])?);
    let (src, tgt) = if out { (ew, e0) } else { (e0, ew) };
    let el = ar.add(MapNode::Elem { kind, src, tgt, at: 0, upright: 0, scalar: F::one() }, "validate_elem")?;
    let d = ar.info(el).deg;
    let (wa, wb): (Vec<Rung>, Vec<Rung>) = if out { (a.to_vec(), vec![]) } else { (vec![], a.to_vec()) };
    let u = hom_unknowns(labels, &wa, &wb, d);
    if u > cap {
        return Err(format!("validate_elem {:?} {:?}: {} unknowns > cap {}", kind, labels, u, cap));
    }
    let ids = ar.lib().try_hom_basis_cap(&Window::new(labels, &wa), &Window::new(labels, &wb), d, u128::MAX)?;
    if ids.len() != 1 {
        return Err(format!("validate_elem {:?} {:?}: Hom^{} has dimension {}", kind, labels, d, ids.len()));
    }
    let tm = ar.add(MapNode::Local { src, tgt, at: 0, first: 0, tmpl: ids[0] }, "validate_elem exact")?;
    let cos: Vec<u64> = ar.flows(ew).support().to_vec();
    let mut ratio: Option<F> = None;
    for pt in points::<F>(ar.k(), 0xe1e7, 2, 0) {
        let mut ev = Eval::new(&ar, pt);
        for &c in &cos {
            let (be, bt) = (ev.block(el, c), ev.block(tm, c));
            for (x, y) in be.a.iter().zip(&bt.a) {
                if y.is_zero() {
                    if !x.is_zero() {
                        return Err(format!("validate_elem {:?} {:?}: not proportional at coset {}", kind, labels, c));
                    }
                    continue;
                }
                let r = x.mul(&y.inv());
                match &ratio {
                    None => ratio = Some(r),
                    Some(q) if *q == r => {}
                    Some(_) => return Err(format!("validate_elem {:?} {:?}: not proportional at coset {}", kind, labels, c)),
                }
            }
        }
    }
    match ratio {
        Some(r) if !r.is_zero() => Ok(r),
        _ => Err(format!("validate_elem {:?} {:?}: closed form vanishes", kind, labels)),
    }
}

impl Leaves {
    pub fn new() -> Leaves {
        Leaves::default()
    }

    fn trace(&self) -> bool {
        self.verbose || std::env::var("KRHOM_TRACE").is_ok()
    }

    /// budget rule: Ok(()) if the exact call may be made
    fn admit(&self, labels: &[u8], a: &[Rung], b: &[Rung], d: i32, unknowns: u128) -> Result<(), String> {
        // one-sided calls (one end the empty word, `call_kind` 0) on windows with a label > 3
        // (outside the measured (3,3,3) class) can cost far more per unknown (A2:
        // Hom^18(∅, F0.3 E0.3) on (3,6), 42 872 unknowns, > 25 min and 3.3 GB, while one-sided
        // calls on labels ≤ 3 up to 83 641 unknowns took ≤ 66 s): above `one_sided_cap` they need
        // a finite exact-shape cost record
        let (fs, ft) = (crate::combinat::schur::flows(labels, a), crate::combinat::schur::flows(labels, b));
        if call_kind(fs, ft) == 0 && labels.iter().any(|&l| l > 3) && unknowns > self.one_sided_cap {
            let pred = self.costs.as_ref().and_then(|c| c.predict_exact(labels, fs, ft, unknowns, d));
            return match pred {
                Some(t) if t.is_finite() && unknowns <= self.hard_cap => Ok(()),
                _ => Err(format!("one-sided step call {:?} [{}] -> [{}] d={}: {} unknowns > {} and no uncensored exact-shape cost record", labels, wstr(a), wstr(b), d, unknowns, self.one_sided_cap)),
            };
        }
        if unknowns <= self.budget {
            return Ok(());
        }
        if unknowns <= self.hard_cap {
            if let Some(c) = &self.costs {
                let (fs, ft) = (crate::combinat::schur::flows(labels, a), crate::combinat::schur::flows(labels, b));
                if let Some(t) = c.predict_exact(labels, fs, ft, unknowns, d) {
                    if t.is_finite() {
                        return Ok(());
                    }
                }
            }
            return Err(format!("step call {:?} [{}] -> [{}] d={}: {} unknowns > budget {} and no uncensored exact-shape cost record", labels, wstr(a), wstr(b), d, unknowns, self.budget));
        }
        Err(format!("step call {:?} [{}] -> [{}] d={}: {} unknowns > hard cap {}", labels, wstr(a), wstr(b), d, unknowns, self.hard_cap))
    }

    /// exact basis of Hom^d on a window (labels, a -> b) as Local nodes src -> tgt at (at, first)
    #[allow(clippy::too_many_arguments)]
    fn exact_local<F: Field>(&self, ar: &mut Arena<F>, labels: &[u8], a: &[Rung], b: &[Rung], d: i32, src: ExprId, tgt: ExprId, at: usize, first: usize, st: &mut LeafStats) -> Result<(Vec<MapId>, u128), LeafFail> {
        let u = hom_unknowns(labels, a, b, d);
        self.admit(labels, a, b, d, u).map_err(|e| lfail(MissKind::Refused, e))?;
        let (wa, wb) = (Window::new(labels, a), Window::new(labels, b));
        let cached = ar.lib_ref().has_hom(&wa, &wb, d);
        if self.trace() && !cached {
            eprintln!("[leaves] exact {:?} [{}] -> [{}] d={} ({} unknowns)", labels, wstr(a), wstr(b), d, u);
        }
        let t0 = Instant::now();
        let ids = ar.lib().try_hom_basis_cap(&wa, &wb, d, u128::MAX).map_err(|e| lfail(MissKind::Refused, format!("step call: {}", e)))?;
        if !cached {
            st.exact_calls += 1;
            st.exact_secs += t0.elapsed().as_secs_f64();
        }
        st.max_unknowns = st.max_unknowns.max(u);
        let mut out = vec![];
        for t in ids {
            out.push(ar.add(MapNode::Local { src, tgt, at, first, tmpl: t }, "leaves step").map_err(|e| lfail(MissKind::Exhausted, e))?);
        }
        Ok((out, u))
    }

    /// The step family S_{t'->t} (Out: W_{t'}·ρ -> W_t; In: W_t -> W_{t'}·ρ) with degree multiset
    /// G and an invertible top-to-top matrix at rep(t). `wfrom` = W_{t'}, `wto` = W_t.
    #[allow(clippy::too_many_arguments)]
    pub fn step_family<F: Field>(&self, ar: &mut Arena<F>, s: &mut Schur, wfrom: &[Rung], rho: &Rung, wto: &[Rung], tp: Tid, t: Tid, g: &Lp, dir: LeafDir, level: usize, st: &mut LeafStats) -> Result<Vec<MapId>, LeafFail> {
        let lam = ar.lam().to_vec();
        let t0 = Instant::now();
        let mut a = wfrom.to_vec();
        a.push(*rho);
        let ea = ar.expr(Expr::new(&lam, &a).map_err(|e| lfail(MissKind::Exhausted, e))?);
        let eb = ar.expr(Expr::new(&lam, wto).map_err(|e| lfail(MissKind::Exhausted, e))?);
        if a == wto {
            if *g != Lp::one() {
                return Err(lfail(MissKind::Exhausted, format!("step {} -> {}: W_t = W_t'·ρ but G = {:?} != 1", s.table(tp).fmt(), s.table(t).fmt(), g)));
            }
            let m = ar.add(MapNode::Id(ea), "leaves step id").map_err(|e| lfail(MissKind::Exhausted, e))?;
            st.steps.push(StepRec { level, from: tp, to: t, kind: StepKind::Id, degs: vec![0], unknowns: 0, secs: 0.0 });
            return Ok(vec![m]);
        }
        let need = degs_of(g);
        let n: usize = need.iter().map(|x| x.1).sum();
        // top-to-top columns: flows of W_{t'}·ρ at rep(t) whose configuration after |W_{t'}|
        // rungs has type t'
        let c = rep_rank_of(s, &lam, t);
        let cols: Vec<usize> = {
            let fi = ar.index(ea);
            let lv = level_of(&lam, &a, wfrom.len());
            let cs = CosetSpace::new(&lam, &lv);
            let tpt = s.table(tp).clone();
            fi.range(c).enumerate().filter(|(_, f)| cs.table(&fi.conf_at(*f, wfrom.len())) == tpt).map(|(j, _)| j).collect()
        };
        if cols.len() != n {
            return Err(lfail(MissKind::Exhausted, format!("step {} -> {}: {} top-to-top columns at rep(t), G(1) = {}", s.table(tp).fmt(), s.table(t).fmt(), cols.len(), n)));
        }
        let (src, tgt) = match dir {
            LeafDir::Out => (ea, eb),
            LeafDir::In => (eb, ea),
        };
        let ws = win_split(&lam, &a, wto);
        let whole_is_window = ws.at == 0 && ws.suffix == 0;
        let mut kinds: Vec<StepKind> = vec![];
        if self.use_window && !whole_is_window {
            kinds.push(StepKind::Window);
        }
        kinds.push(StepKind::Whole);
        let pt = points::<F>(ar.k(), self.seed ^ 0x57e9, 1, 0).remove(0);
        // E2/E3 closed form for a one-dimensional cup/cap step (validated once per arena)
        if n == 1 && self.use_elem {
            if let Some((ek, key)) = elem_pattern(&ws.labels, &ws.a, &ws.b, dir) {
                let valid = ar.elem_certified(ek, key).is_some() || match validate_elem::<F>(ek, &ws.labels, &ws.a, 120_000) {
                    Ok(_) => {
                        ar.certify_elem(ek, key, F::one());
                        true
                    }
                    Err(e) => {
                        if self.trace() {
                            eprintln!("[leaves] {}", e);
                        }
                        false
                    }
                };
                if valid {
                    if let Ok(m) = ar.add(MapNode::Elem { kind: ek, src, tgt, at: ws.at, upright: ws.first, scalar: F::one() }, "leaves step elem") {
                        let d = ar.info(m).deg;
                        if d == need[0].0 {
                            let mut ev = Eval::new(ar, pt.clone());
                            let blk = ev.block(m, c);
                            if !blk.a[cols[0]].is_zero() {
                                st.steps.push(StepRec { level, from: tp, to: t, kind: StepKind::Elem, degs: vec![d], unknowns: 0, secs: t0.elapsed().as_secs_f64() });
                                return Ok(vec![m]);
                            }
                        }
                    }
                }
            }
        }
        let mut cands: Vec<(i32, MapId)> = vec![];
        let mut maxu = 0u128;
        let mut last_err = String::new();
        for kind in kinds {
            for &(d, _) in &need {
                let r = match kind {
                    StepKind::Window => {
                        let (la, lb) = match dir {
                            LeafDir::Out => (&ws.a, &ws.b),
                            LeafDir::In => (&ws.b, &ws.a),
                        };
                        self.exact_local(ar, &ws.labels, la, lb, d, src, tgt, ws.at, ws.first, st)
                    }
                    _ => {
                        let (la, lb) = match dir {
                            LeafDir::Out => (&a, &wto.to_vec()),
                            LeafDir::In => (&wto.to_vec(), &a),
                        };
                        self.exact_local(ar, &lam, la, lb, d, src, tgt, 0, 0, st)
                    }
                };
                match r {
                    Ok((ms, u)) => {
                        maxu = maxu.max(u);
                        cands.extend(ms.into_iter().map(|m| (d, m)));
                    }
                    Err(e) => {
                        if kind == StepKind::Whole {
                            return Err(e);
                        }
                        last_err = e.what;
                    }
                }
            }
            // greedy by ascending degree on the top-to-top columns
            let mut ev = Eval::new(ar, pt.clone());
            let mut rows: Vec<Vec<F>> = vec![];
            let mut chosen: Vec<MapId> = vec![];
            let mut ok = true;
            for &(d, k) in &need {
                let mut got = 0;
                for &(dd, m) in &cands {
                    if dd != d || got >= k {
                        continue;
                    }
                    let blk = ev.block(m, c);
                    let v: Vec<F> = cols.iter().map(|&j| blk.a[j].clone()).collect();
                    rows.push(v);
                    if crate::algebra::dense::rank(&rows, n) == rows.len() {
                        chosen.push(m);
                        got += 1;
                    } else {
                        rows.pop();
                    }
                }
                if got < k {
                    ok = false;
                    last_err = format!("{:?}: degree {}: {} of {} top-to-top generators", kind, d, got, k);
                    break;
                }
            }
            if ok {
                let degs: Vec<i32> = need.iter().flat_map(|&(d, k)| std::iter::repeat(d).take(k)).collect();
                let secs = t0.elapsed().as_secs_f64();
                if self.trace() {
                    eprintln!("[leaves] level {} step {} -> {} ({:?}, {:?}): {:?} degs {:?} max unknowns {} ({:.2}s)", level, s.table(tp).fmt(), s.table(t).fmt(), dir, kind, wstr(&a), degs, maxu, secs);
                }
                st.steps.push(StepRec { level, from: tp, to: t, kind, degs, unknowns: maxu, secs });
                return Ok(chosen);
            }
        }
        Err(lfail(MissKind::Exhausted, format!("step {} -> {} ([{}] -> [{}]): {}", s.table(tp).fmt(), s.table(t).fmt(), wstr(&a), wstr(wto), last_err)))
    }

    /// Junction words of every level (pure combinatorics, before any map is built): the default
    /// rule (`junction_word`: prefix extension, canonical ladder, BFS; ∅ for the final identity
    /// type), then, if some step exceeds the budget, coordinate descent over the reduced words of
    /// each junction type (≤ alt_len rungs) on the estimated step unknowns (`step_score`) of the
    /// steps into and out of it (objective: max, then flows), at most 4 sweeps.
    pub fn plan_words(&self, s: &mut Schur, sc: &mut StepCache, w: &[Rung], js: &[Vec<crate::combinat::dcoset::Junction>], id: Tid) -> Vec<Vec<(Tid, Vec<Rung>)>> {
        let lam = s.lam.clone();
        let m = w.len();
        let mut plan: Vec<Vec<(Tid, Vec<Rung>)>> = vec![vec![(id, vec![])]];
        for k in 1..=m {
            let prev = plan[k - 1].clone();
            let mut lv = vec![];
            for j in &js[k] {
                let wt = if k == m && j.t == id { vec![] } else { junction_word(s, j.t, &prev, Some(&w[k - 1]), self.bfs_len, self.bfs_cap).0 };
                lv.push((j.t, wt));
            }
            plan.push(lv);
        }
        if self.alt_len == 0 {
            return plan;
        }
        // step score of (level k type index a) -> (level k+1 type index b)
        let step = |s: &mut Schur, sc: &mut StepCache, plan: &Vec<Vec<(Tid, Vec<Rung>)>>, k: usize, a: usize, b: usize| -> Option<u128> {
            let (tp, wp) = &plan[k][a];
            let (t, wt) = &plan[k + 1][b];
            let g = sc.step(s, *tp, &w[k]).into_iter().find(|x| x.0 == *t).map(|x| x.1)?;
            Some(step_score(&lam, wp, &w[k], wt, &g))
        };
        let local = |s: &mut Schur, sc: &mut StepCache, plan: &Vec<Vec<(Tid, Vec<Rung>)>>, k: usize, j: usize| -> u128 {
            let mut mx = 0u128;
            for a in 0..plan[k - 1].len() {
                if let Some(x) = step(s, sc, plan, k - 1, a, j) {
                    mx = mx.max(x);
                }
            }
            if k < m {
                for b in 0..plan[k + 1].len() {
                    if let Some(x) = step(s, sc, plan, k, j, b) {
                        mx = mx.max(x);
                    }
                }
            }
            mx
        };
        let mut alts: HashMap<Tid, Vec<Vec<Rung>>> = HashMap::new();
        for _sweep in 0..4 {
            let mut changed = false;
            for k in 1..m {
                for j in 0..plan[k].len() {
                    let s0 = local(s, sc, &plan, k, j);
                    if s0 <= self.budget {
                        continue;
                    }
                    let t = plan[k][j].0;
                    let cands = alts.entry(t).or_insert_with(|| crate::combinat::dcoset::reduced_words(s, t, self.alt_len, self.alt_cap)).clone();
                    let orig = plan[k][j].1.clone();
                    let mut best = (s0, crate::combinat::schur::flows(&lam, &orig), orig.clone());
                    for c in cands {
                        plan[k][j].1 = c.clone();
                        let sa = local(s, sc, &plan, k, j);
                        let fa = crate::combinat::schur::flows(&lam, &c);
                        if (sa, fa) < (best.0, best.1) {
                            best = (sa, fa, c);
                        }
                    }
                    if best.2 != orig {
                        changed = true;
                        if self.trace() {
                            eprintln!("[leaves] plan level {} {}: [{}] ({} unknowns) -> [{}] ({} unknowns)", k, s.table(t).fmt(), wstr(&orig), s0, wstr(&best.2), best.0);
                        }
                    }
                    plan[k][j].1 = best.2;
                }
            }
            if !changed {
                break;
            }
        }
        if self.joint {
            self.joint_descent(s, sc, w, &mut plan, id);
        }
        plan
    }

    /// C2-7 joint junction words (B5 §2.1 "adopt"; off by default, `Leaves::joint`). For every
    /// step t' -ρ-> t of the plan whose `step_score` exceeds the budget, choose W_{t'} and W_t
    /// **together** as P·A and P·B with a common prefix P (P over the prefixes of the current
    /// words of t' and t, A and B suffix words of ≤ `joint_suffix` rungs with P·A reduced for t'
    /// and P·B reduced for t), so that the step window shrinks to A·ρ -> B on the labels after P.
    /// A candidate pair is adopted only when it lowers the plan objective (max step score over
    /// every step into or out of t' and t, then flows): the plan structure (one word per junction
    /// type per level) is unchanged, so no detour maps are needed. The identity type keeps ∅;
    /// a final non-identity type q may change its word (pairs_leaves composes the exit chain
    /// W_q -> L_q, which must then be found by the planner). Pure combinatorics; nothing is
    /// certified here (acceptance is unchanged: step families must reach degree multiset G and
    /// full top-to-top rank, then certify_pairs / certify_closure).
    pub fn joint_descent(&self, s: &mut Schur, sc: &mut StepCache, w: &[Rung], plan: &mut Vec<Vec<(Tid, Vec<Rung>)>>, id: Tid) {
        let lam = s.lam.clone();
        let m = w.len();
        let wl: std::cell::RefCell<HashMap<Vec<u8>, Schur>> = std::cell::RefCell::new(HashMap::new());
        let step = |s: &mut Schur, sc: &mut StepCache, plan: &Vec<Vec<(Tid, Vec<Rung>)>>, k: usize, a: usize, b: usize| -> Option<u128> {
            let (tp, wp) = &plan[k][a];
            let (t, wt) = &plan[k + 1][b];
            let g = sc.step(s, *tp, &w[k]).into_iter().find(|x| x.0 == *t).map(|x| x.1)?;
            Some(step_score_w(&mut wl.borrow_mut(), &lam, wp, &w[k], wt, &g))
        };
        // max step score over every step touching (level k, index j)
        let local = |s: &mut Schur, sc: &mut StepCache, plan: &Vec<Vec<(Tid, Vec<Rung>)>>, k: usize, j: usize| -> u128 {
            let mut mx = 0u128;
            if k > 0 {
                for a in 0..plan[k - 1].len() {
                    if let Some(x) = step(s, sc, plan, k - 1, a, j) {
                        mx = mx.max(x);
                    }
                }
            }
            if k < m {
                for b in 0..plan[k + 1].len() {
                    if let Some(x) = step(s, sc, plan, k, j, b) {
                        mx = mx.max(x);
                    }
                }
            }
            mx
        };
        let mut suff_cache: HashMap<(Vec<Rung>, Tid), Vec<Vec<Rung>>> = HashMap::new();
        for _sweep in 0..self.joint_sweeps {
            let mut changed = false;
            for k in 0..m {
                for a in 0..plan[k].len() {
                    for b in 0..plan[k + 1].len() {
                        let Some(s0) = step(s, sc, plan, k, a, b) else { continue };
                        if s0 <= self.budget {
                            continue;
                        }
                        let (tp, t) = (plan[k][a].0, plan[k + 1][b].0);
                        let fix_a = k == 0 || tp == id;
                        let fix_b = t == id;
                        let (wa0, wb0) = (plan[k][a].1.clone(), plan[k + 1][b].1.clone());
                        let obj0 = local(s, sc, plan, k, a).max(local(s, sc, plan, k + 1, b));
                        let f0 = crate::combinat::schur::flows(&lam, &wa0) + crate::combinat::schur::flows(&lam, &wb0);
                        // prefixes
                        let mut prefs: Vec<Vec<Rung>> = vec![];
                        for src in [&wa0, &wb0] {
                            for l in 0..=src.len() {
                                let p = src[..l].to_vec();
                                if !prefs.contains(&p) {
                                    prefs.push(p);
                                }
                            }
                        }
                        // candidate pairs ranked by the score of this step alone
                        let mut cands: Vec<(u128, u128, Vec<Rung>, Vec<Rung>)> = vec![];
                        let g = sc.step(s, tp, &w[k]).into_iter().find(|x| x.0 == t).map(|x| x.1).unwrap();
                        for p in &prefs {
                            let la: Vec<Vec<Rung>> = if fix_a {
                                if wa0.starts_with(p) { vec![wa0.clone()] } else { vec![] }
                            } else {
                                suff_cache.entry((p.clone(), tp)).or_insert_with(|| suffix_words(s, p, tp, self.joint_suffix, self.joint_cap)).clone()
                            };
                            if la.is_empty() {
                                continue;
                            }
                            let lb: Vec<Vec<Rung>> = if fix_b {
                                if wb0.starts_with(p) { vec![wb0.clone()] } else { vec![] }
                            } else {
                                suff_cache.entry((p.clone(), t)).or_insert_with(|| suffix_words(s, p, t, self.joint_suffix, self.joint_cap)).clone()
                            };
                            for wa in &la {
                                for wb in &lb {
                                    let sc1 = step_score_w(&mut wl.borrow_mut(), &lam, wa, &w[k], wb, &g);
                                    if sc1 < s0 {
                                        cands.push((sc1, crate::combinat::schur::flows(&lam, wa) + crate::combinat::schur::flows(&lam, wb), wa.clone(), wb.clone()));
                                    }
                                }
                            }
                        }
                        cands.sort_by(|x, y| (x.0, x.1).cmp(&(y.0, y.1)));
                        cands.dedup_by(|x, y| x.2 == y.2 && x.3 == y.3);
                        let mut best = (obj0, f0, wa0.clone(), wb0.clone());
                        for (_, _, wa, wb) in cands.into_iter().take(self.joint_keep) {
                            plan[k][a].1 = wa.clone();
                            plan[k + 1][b].1 = wb.clone();
                            let o = local(s, sc, plan, k, a).max(local(s, sc, plan, k + 1, b));
                            let f = crate::combinat::schur::flows(&lam, &wa) + crate::combinat::schur::flows(&lam, &wb);
                            if (o, f) < (best.0, best.1) {
                                best = (o, f, wa, wb);
                            }
                        }
                        plan[k][a].1 = best.2.clone();
                        plan[k + 1][b].1 = best.3.clone();
                        if best.2 != wa0 || best.3 != wb0 {
                            changed = true;
                            if self.trace() {
                                eprintln!("[leaves] joint level {} {} -> {}: [{}]/[{}] ({}) -> [{}]/[{}] ({})", k + 1, s.table(tp).fmt(), s.table(t).fmt(), wstr(&wa0), wstr(&wb0), obj0, wstr(&best.2), wstr(&best.3), best.0);
                            }
                        }
                    }
                }
            }
            if !changed {
                break;
            }
        }
    }

    /// C2-7 per-predecessor final words (B5 §2.1 item 5, C_B5 F2; with `joint`): for every
    /// predecessor t'_i (index i at level m−1) of a final non-identity type q, the word
    /// W'_{q,i} reduced for q that minimizes the W-aware score of its own final step (candidates:
    /// the plan word W_q and P·B, P over the prefixes of W_{t'_i}·ρ, |B| ≤ `joint_suffix`), kept
    /// only when it is strictly better than W_q and a degree-0 rex chain W'_{q,i} ⇄ W_q is found
    /// (B2 planner, `leaf_pair_plan_opts`). `build_leaves` composes predecessor i's final leaves
    /// with that chain, so every leaf still lands in W_q (the chains are degree-0 maps between
    /// reduced words of q with top scalar 1). Returns (q, i) -> W'_{q,i}.
    pub fn final_pred_words(&self, s: &mut Schur, sc: &mut StepCache, w: &[Rung], plan: &[Vec<(Tid, Vec<Rung>)>], id: Tid) -> HashMap<(Tid, usize), Vec<Rung>> {
        let m = w.len();
        self.step_words(s, sc, w, plan, id, false).into_iter().filter(|((k, _, _), _)| *k == m).map(|((_, t, i), v)| ((t, i), v)).collect()
    }

    /// C2-7 per-step target words ("detours", B5 §2.1 item 4; `final_pred_words` is the final
    /// level): for every step t'_i -ρ-> t at level k (all levels when `all_levels`, else only the
    /// final level; never into the identity type) whose W-aware score exceeds the budget, the word
    /// W'_{t,i} reduced for t minimizing that step's score (candidates P·B, P over the prefixes of
    /// W_{t'_i}·ρ, |B| ≤ `joint_suffix`), kept only when strictly better and degree-0 rex chains
    /// W'_{t,i} ⇄ W_t are found (B2 planner, `leaf_pair_plan_opts`). `build_leaves` composes the
    /// step family into W'_{t,i} with the chain, so every leaf at t still lands in W_t. Returns
    /// (k, t, i) -> W'_{t,i}.
    pub fn step_words(&self, s: &mut Schur, sc: &mut StepCache, w: &[Rung], plan: &[Vec<(Tid, Vec<Rung>)>], id: Tid, all_levels: bool) -> HashMap<(usize, Tid, usize), Vec<Rung>> {
        let mut out = HashMap::new();
        let m = w.len();
        if !self.joint || m == 0 {
            return out;
        }
        let lam = s.lam.clone();
        let mut wl: HashMap<Vec<u8>, Schur> = HashMap::new();
        let opts = leaf_pair_plan_opts();
        let k0 = if all_levels { 1 } else { m };
        for k in k0..=m {
            for (q, wq) in plan[k].clone() {
                if q == id {
                    continue;
                }
                for (i, (tp, wp)) in plan[k - 1].iter().enumerate() {
                    let Some(g) = sc.step(s, *tp, &w[k - 1]).into_iter().find(|x| x.0 == q).map(|x| x.1) else { continue };
                    let s0 = step_score_w(&mut wl, &lam, wp, &w[k - 1], &wq, &g);
                    if s0 <= self.budget {
                        continue;
                    }
                    let mut a = wp.clone();
                    a.push(w[k - 1]);
                    let mut cands: Vec<(u128, u128, Vec<Rung>)> = vec![];
                    for l in 0..=a.len() {
                        for c in suffix_words(s, &a[..l], q, self.joint_suffix, self.joint_cap) {
                            let sc1 = step_score_w(&mut wl, &lam, wp, &w[k - 1], &c, &g);
                            if sc1 < s0 {
                                cands.push((sc1, crate::combinat::schur::flows(&lam, &c), c));
                            }
                        }
                    }
                    cands.sort_by(|x, y| (x.0, x.1).cmp(&(y.0, y.1)));
                    cands.dedup_by(|x, y| x.2 == y.2);
                    if self.trace() {
                        eprintln!("[leaves] level {} step {} -> {} ({}): {} better candidate words{}", k, s.table(*tp).fmt(), s.table(q).fmt(), s0, cands.len(), cands.first().map(|c| format!(", best [{}] {}", wstr(&c.2), c.0)).unwrap_or_default());
                    }
                    let mut pl = crate::library::rexplan::Planner::new();
                    for (sc1, _, c) in cands.into_iter().take(8) {
                        let (Ok(c1), Ok(c2)) = (pl.plan(s, q, &c, &wq, &opts), pl.plan(s, q, &wq, &c, &opts)) else { continue };
                        {
                            let mut memo = self.chain_memo.borrow_mut();
                            memo.insert((q, c.clone(), wq.clone()), c1);
                            memo.insert((q, wq.clone(), c.clone()), c2);
                        }
                        {
                            if self.trace() {
                                eprintln!("[leaves] level {} word for predecessor {} of {}: [{}] ({}) instead of [{}] ({})", k, s.table(*tp).fmt(), s.table(q).fmt(), wstr(&c), sc1, wstr(&wq), s0);
                            }
                            out.insert((k, q, i), c);
                            break;
                        }
                    }
                }
            }
        }
        out
    }

    /// Lemma A, level by level (D1 §1.3/§2.3): the leaf families of X for the junction types at
    /// its top level that lie on paths to `targets`. Final words: ∅ for the identity type,
    /// otherwise `junction_word` (prefix extension, canonical ladder, BFS). Intermediate
    /// families are checked (diagnostic) to be square with an invertible block at rep(t).
    pub fn build_leaves<F: Field>(&self, ar: &mut Arena<F>, s: &mut Schur, x: ExprId, targets: &[Tid], dir: LeafDir, st: &mut LeafStats) -> Result<Vec<LeafFamily>, LeafFail> {
        let lam = ar.lam().to_vec();
        let w = ar.expr_of(x).rungs.clone();
        let id = s.intern(Table::identity(&lam));
        let mut sc = StepCache::new();
        let js = junction_types(s, &mut sc, &w, targets);
        if js[w.len()].is_empty() {
            return Err(lfail(MissKind::Exhausted, "build_leaves: no target type is reached".into()));
        }
        let plan = self.plan_words(s, &mut sc, &w, &js, id);
        let finals = self.step_words(s, &mut sc, &w, &plan, id, self.detours);
        let mut chain_cache = HashMap::new();
        let e0 = ar.expr(Expr::new(&lam, &[]).map_err(|e| lfail(MissKind::Exhausted, e))?);
        let id0 = ar.add(MapNode::Id(e0), "leaves base").map_err(|e| lfail(MissKind::Exhausted, e))?;
        let mut cur: Vec<LeafFamily> = vec![LeafFamily { t: id, word: vec![], wid: e0, maps: vec![id0] }];
        let pt = points::<F>(ar.k(), self.seed ^ 0xd1a9, 1, 0).remove(0);
        for k in 1..=w.len() {
            let rho = w[k - 1];
            let mut next: Vec<LeafFamily> = vec![];
            for j in &js[k] {
                let t = j.t;
                let wt = plan[k].iter().find(|x| x.0 == t).map(|x| x.1.clone()).expect("planned junction word");
                if !crate::combinat::dcoset::reduced_for(s, &wt, t).ok() {
                    return Err(lfail(MissKind::Exhausted, format!("build_leaves: no reduced junction word for {} at level {}", s.table(t).fmt(), k)));
                }
                st.max_word_flows = st.max_word_flows.max(crate::combinat::schur::flows(&lam, &wt));
                let wid = ar.expr(Expr::new(&lam, &wt).map_err(|e| lfail(MissKind::Exhausted, e))?);
                let mut maps = vec![];
                for (fi, f) in cur.iter().enumerate() {
                    let Some(g) = sc.step(s, f.t, &rho).into_iter().find(|x| x.0 == t).map(|x| x.1) else { continue };
                    // C2-7: per-predecessor final word W'_{q,i}, then the rex chain to / from W_q
                    let alt = finals.get(&(k, t, fi)).cloned();
                    let fam = match &alt {
                        None => self.step_family(ar, s, &f.word, &rho, &wt, f.t, t, &g, dir, k, st)?,
                        Some(wa) => {
                            let fam = self.step_family(ar, s, &f.word, &rho, wa, f.t, t, &g, dir, k, st)?;
                            let mut pl = crate::library::rexplan::Planner::new();
                            let opts = leaf_pair_plan_opts();
                            let (from, to) = match dir {
                                LeafDir::Out => (wa.clone(), wt.clone()),
                                LeafDir::In => (wt.clone(), wa.clone()),
                            };
                            let memo = self.chain_memo.borrow().get(&(t, from.clone(), to.clone())).cloned();
                            let ch = match memo {
                                Some(c) => c,
                                None => pl.plan(s, t, &from, &to, &opts).map_err(|e| lfail(MissKind::Refused, format!("final chain: {}", e)))?,
                            };
                            let (r, _) = crate::library::leafgen::realize_chain_opts(ar, &ch, &mut chain_cache, self.seed ^ 0xf1a1 ^ fi as u64, true).map_err(|e| lfail(MissKind::Refused, format!("final chain: {}", e)))?;
                            let mut out = vec![];
                            for sm in fam {
                                let n = match dir {
                                    LeafDir::Out => MapNode::Comp(r, sm),
                                    LeafDir::In => MapNode::Comp(sm, r),
                                };
                                out.push(ar.add(n, "leaves final chain").map_err(|e| lfail(MissKind::Exhausted, e))?);
                            }
                            out
                        }
                    };
                    for &l in &f.maps {
                        let kr = ar.add(MapNode::Kron(l, vec![rho]), "leaves kron").map_err(|e| lfail(MissKind::Exhausted, e))?;
                        for &sm in &fam {
                            let is_id = matches!(ar.node(sm), MapNode::Id(_));
                            let m = if is_id {
                                kr
                            } else {
                                let n = match dir {
                                    LeafDir::Out => MapNode::Comp(sm, kr),
                                    LeafDir::In => MapNode::Comp(kr, sm),
                                };
                                ar.add(n, "leaves leaf").map_err(|e| lfail(MissKind::Exhausted, e))?
                            };
                            maps.push(m);
                        }
                    }
                }
                // diagnostic: square and invertible at rep(t)
                let c = rep_rank_of(s, &lam, t);
                let want = j.fwd.at1() as usize;
                if maps.len() != want {
                    return Err(lfail(MissKind::Exhausted, format!("build_leaves: level {} type {}: {} leaves, c_t(X)(1) = {}", k, s.table(t).fmt(), maps.len(), want)));
                }
                {
                    let mut ev = Eval::new(ar, pt.clone());
                    let rows: Vec<Vec<F>> = maps.iter().map(|&m| ev.block(m, c).a.clone()).collect();
                    if rows.iter().any(|r| r.len() != want) || crate::algebra::dense::rank(&rows, want) != want {
                        return Err(lfail(MissKind::Exhausted, format!("build_leaves: level {} type {}: block at rep(t) not invertible ({} leaves)", k, s.table(t).fmt(), want)));
                    }
                }
                st.maps_built += maps.len();
                next.push(LeafFamily { t, word: wt, wid, maps });
            }
            cur = next;
        }
        Ok(cur)
    }

    /// closure generators (outs: L_p -> R, ins: R -> L_p) through Lemma A to the identity type,
    /// certified by `certify_closure`
    pub fn closure_leaves<F: Field>(&mut self, ar: &mut Arena<F>, s: &mut Schur, lp: ExprId, p: Tid, outs: bool) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        let t0 = Instant::now();
        let lam = ar.lam().to_vec();
        let want = cdiag_degrees(s, p);
        let mut st = LeafStats::default();
        let mark = ar.mark();
        let mk = |kind: MissKind, what: String| vec![Missing { what, degree: 0, stratum: None, got: 0, want: want.len(), kind }];
        let r = (|| -> Result<Vec<MapId>, LeafFail> {
            closure_precheck(ar, s, lp, p, None).map_err(|e| lfail(MissKind::Refused, e))?;
            let id = s.intern(Table::identity(&lam));
            let fams = self.build_leaves(ar, s, lp, &[id], if outs { LeafDir::Out } else { LeafDir::In }, &mut st)?;
            let f = fams.into_iter().find(|f| f.t == id).ok_or_else(|| lfail(MissKind::Exhausted, "no identity family".into()))?;
            Ok(f.maps)
        })();
        let ms = match r {
            Ok(ms) => ms,
            Err(e) => {
                ar.commit(mark, &[]);
                self.last = st;
                let ss = SearchStats { secs: t0.elapsed().as_secs_f64(), exact_calls: self.last.exact_calls, exact_secs: self.last.exact_secs, ..SearchStats::default() };
                self.last_fail = Some(ss);
                return Err(mk(e.kind, format!("Leaves closure {} ({}): {}", s.table(p).fmt(), if outs { "outs" } else { "ins" }, e.what)));
            }
        };
        let kept = ar.commit(mark, &ms);
        let r0 = ar.expr(Expr::new(&lam, &[]).expect("empty word"));
        if let Err(e) = certify_closure(ar, &lam, r0, lp, None, &kept, &want, outs, self.seed ^ 0xce57) {
            self.last = st;
            self.last_fail = Some(SearchStats { secs: t0.elapsed().as_secs_f64(), ..SearchStats::default() });
            return Err(mk(MissKind::Exhausted, format!("Leaves closure {}: certify_closure: {}", s.table(p).fmt(), e)));
        }
        let ss = SearchStats { candidates: kept.len(), evaluated: kept.len(), secs: t0.elapsed().as_secs_f64(), exact_calls: st.exact_calls, exact_secs: st.exact_secs, max_flows: st.max_word_flows as u64, detect: "leaves", ..SearchStats::default() };
        self.last = st;
        Ok((kept, ss))
    }
}

fn refuse_leaves(what: &str) -> Vec<Missing> {
    vec![Missing { what: format!("Leaves: {} not implemented (M4c A4/A5)", what), degree: 0, stratum: None, got: 0, want: 1, kind: MissKind::Refused }]
}

impl<F: Field> HomProvider<F> for Leaves {
    fn name(&self) -> &'static str {
        "Leaves"
    }
    fn full(&mut self, _: &mut Arena<F>, _: &mut Schur, _: ExprId, _: ExprId, _: i32, _: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        self.last_fail = None;
        Err(refuse_leaves("full"))
    }
    fn pairs(&mut self, _: &mut Arena<F>, _: &mut Schur, _: ExprId, _: Tid, _: ExprId, _: Option<MapId>, _: i32, _: usize, _: &Budget) -> Result<(Pairs<F>, SearchStats), Vec<Missing>> {
        self.last_fail = None;
        Err(refuse_leaves("pairs"))
    }
    fn closure(&mut self, ar: &mut Arena<F>, s: &mut Schur, lp: ExprId, p: Tid, ep: Option<MapId>, outs: bool, _b: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        self.last_fail = None;
        if ep.is_some() {
            return Err(refuse_leaves("closure with e_p"));
        }
        self.closure_leaves(ar, s, lp, p, outs)
    }
    fn end0(&mut self, _: &mut Arena<F>, _: &mut Schur, _: ExprId, _: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        self.last_fail = None;
        Err(refuse_leaves("end0"))
    }
    fn failed_stats(&self) -> Option<SearchStats> {
        self.last_fail.clone()
    }
}

// ------------------------------------------------------------------ production wiring (A3)

/// The M4c production chain (A3 closure wiring + B4w pair wiring), built around the frozen
/// production chain `sindec::provider_for_stamped` without editing sindec.rs:
/// - closure calls: Stage 0 product closure first for product types (A1 report §4: before
///   Exact), then the production chain, then Leaves (`leaves_first = false`, M4C_PLAN §0.2
///   "Leaves last"), or Leaves first (`leaves_first = true`);
/// - pair calls (B4w): the production chain first; on failure, `leaf_pairs` (rex-class targets:
///   `rex_class_pairs`; other targets: `pairs_leaves`, Lemma A within the Leaves budget plus
///   rex-chain exit/entry), Leaves last;
/// - `full`, `end0` and closure with e_p: the chain alone.
///
/// The stamp extends the chain stamp (`STAGE0_STAMP_PF`, `LEAVES_STAMP(order)`,
/// `LEAF_PAIRS_STAMP`), so layer A built with it is separate from the M5c layer A.
pub fn provider_for_stamped_m4c<F: Field>(lam: &[u8], costs: Option<&std::path::Path>, slowdown: f64, leaves_first: bool) -> (Box<dyn HomProvider<F>>, String) {
    let (chain, stamp) = crate::library::sindec::provider_for_stamped::<F>(lam, costs, slowdown);
    let mut pc = ProductClosure::new(chain);
    pc.product_first = true;
    let mut lv = Leaves::new();
    lv.costs = costs.and_then(|p| CostTable::load(p).ok());
    let order = if leaves_first { "leaves-first" } else { "leaves-last" };
    // C2-7 (off by default): KRHOM_LEAVES_C27=1 turns on joint junction words (`joint_descent`)
    // and Leaves pairs for targets with e_q; the stamp then carries the extra token C27_STAMP
    // (a route token: `sindec::route_header` drops it, so with --stamp-route layer A is shared)
    let c27 = std::env::var("KRHOM_LEAVES_C27").map(|v| v == "1").unwrap_or(false);
    lv.joint = c27;
    lv.eq_pairs = c27;
    lv.detours = c27 && std::env::var("KRHOM_LEAVES_C27_DETOURS").map(|v| v != "0").unwrap_or(true);
    (
        Box::new(ClosureLeaves { inner: Box::new(pc), leaves: lv, leaves_first, leaf_pairs: true, last_fail: None }),
        format!("{} {} {}({}) {}{}", stamp, STAGE0_STAMP_PF, LEAVES_STAMP, order, LEAF_PAIRS_STAMP, if c27 { format!(" {}", C27_STAMP) } else { String::new() }),
    )
}

/// The A3 constructor (closure wiring only, A3 stamp): Leaves is never asked for pairs.
pub fn provider_for_stamped_m4c_a3<F: Field>(lam: &[u8], costs: Option<&std::path::Path>, slowdown: f64, leaves_first: bool) -> (Box<dyn HomProvider<F>>, String) {
    let (chain, stamp) = crate::library::sindec::provider_for_stamped::<F>(lam, costs, slowdown);
    let mut pc = ProductClosure::new(chain);
    pc.product_first = true;
    let mut lv = Leaves::new();
    lv.costs = costs.and_then(|p| CostTable::load(p).ok());
    let order = if leaves_first { "leaves-first" } else { "leaves-last" };
    (Box::new(ClosureLeaves { inner: Box::new(pc), leaves: lv, leaves_first, leaf_pairs: false, last_fail: None }), format!("{} {} {}({})", stamp, STAGE0_STAMP_PF, LEAVES_STAMP, order))
}

/// C2-7 route token (joint junction words W1/W2-aware, Leaves pairs with e_q); see
/// `provider_for_stamped_m4c`
pub const C27_STAMP: &str = "c27=joint(sfx4,keep40,w12)+eqpairs+detours";

/// provider stamp suffix of the B4w pair route (layer-A headers)
pub const LEAF_PAIRS_STAMP: &str = "leafpairs=LeafPairs(v1,last,rex-chain|lemmaA+exit-entry,plan=exact1e5/len12/visit6e4,assoc-as-template,persistable-only)";

/// Planner options of the production pair route (those of the B2 `gp0 --b4` runs).
pub fn leaf_pair_plan_opts() -> crate::library::rexplan::PlanOpts {
    crate::library::rexplan::PlanOpts { exact_budget: 100_000, max_len: 12, max_visit: 60_000, ..crate::library::rexplan::PlanOpts::default() }
}

/// Elem nodes reachable from `roots` whose constant is not persisted in layer A (every Elem node
/// except the theorem-level Assoc nodes, which smapio marks and re-checks on import). A layer-A
/// file containing one would be rejected on load in a fresh session (A2 finding), so the
/// production pair route refuses such maps.
pub fn unpersistable_elems<F: Field>(ar: &Arena<F>, roots: &[MapId]) -> Vec<MapId> {
    let mut out = vec![];
    let mut stack = roots.to_vec();
    let mut seen = std::collections::HashSet::new();
    while let Some(m) = stack.pop() {
        if !seen.insert(m) {
            continue;
        }
        match ar.node(m) {
            MapNode::Elem { kind, .. } => {
                let thm = *kind == crate::colored::smap::ElemKind::Assoc && ar.elem_windows(m).map(|(a, b)| crate::library::leafgen::thm_assoc_registered(a, b)).unwrap_or(false);
                if !thm {
                    out.push(m);
                }
            }
            MapNode::Comp(g, f) => {
                stack.push(*g);
                stack.push(*f);
            }
            MapNode::Lin(t) => stack.extend(t.iter().map(|x| x.1)),
            MapNode::Kron(f, _) => stack.push(*f),
            _ => {}
        }
    }
    out.sort_unstable();
    out
}

/// Route of a `leaf_pairs` result.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum LeafPairRoute {
    /// X is a reduced expression of q (s = 0, mult 1): π, ι are realized rex chains
    RexClass,
    /// Lemma A leaves at q plus rex-chain exit/entry
    LemmaA,
}

/// B4w: the Leaves pair route of the production provider. Rex-class targets (X reduced for q,
/// s = 0, mult 1) go to `rex_class_pairs`; the others to `pairs_leaves` (whose Lemma A steps obey
/// the Leaves budget rule: ≤ 1e5 unknowns, ≤ 4e5 with an uncensored record, the whole-word
/// fallback included). Everything runs inside a mark/commit; the result is refused when a map
/// uses an Elem constant that layer A cannot persist. Acceptance is `certify_pairs` here and
/// then sindec's own checks (dechom_postcheck, D14 scaling, Π∘Ι = id, kept multiset, layer-A
/// save and re-certified load), exactly as for any provider.
#[allow(clippy::too_many_arguments)]
pub fn leaf_pairs<F: Field>(lv: &Leaves, ar: &mut Arena<F>, s: &mut Schur, x: ExprId, q: Tid, lq: ExprId, eq: Option<MapId>, sdeg: i32, mult: usize) -> Result<(Pairs<F>, LeafPairRoute, LeafStats), String> {
    if eq.is_some() && !lv.eq_pairs {
        return Err("Leaves pairs: target with an idempotent e_q (not implemented)".into());
    }
    let xw = ar.expr_of(x).rungs.clone();
    let rex = eq.is_none() && sdeg == 0 && mult == 1 && crate::combinat::dcoset::reduced_for(s, &xw, q).ok();
    let opts = leaf_pair_plan_opts();
    let mark = ar.mark();
    let mut st = LeafStats::default();
    let r = if rex {
        rex_class_pairs_opts(ar, s, x, q, lq, &opts, lv.seed ^ 0xb4b4, true).map(|(p, _)| (p, LeafPairRoute::RexClass))
    } else {
        pairs_leaves_eq(lv, ar, s, x, q, lq, eq, sdeg, mult, &opts, &mut st, true).map(|(p, _)| (p, LeafPairRoute::LemmaA))
    };
    let (p, route) = match r {
        Ok(v) => v,
        Err(e) => {
            ar.commit(mark, &[]);
            return Err(format!("Leaves pairs ({}): {}", if rex { "rex class" } else { "Lemma A" }, e));
        }
    };
    let roots: Vec<MapId> = p.iotas.iter().chain(&p.pis).copied().collect();
    let bad = unpersistable_elems(ar, &roots);
    if !bad.is_empty() {
        ar.commit(mark, &[]);
        return Err(format!("Leaves pairs: {} Elem node(s) with constants not persisted in layer A", bad.len()));
    }
    let kept = ar.commit(mark, &roots);
    let n = p.iotas.len();
    Ok((Pairs { iotas: kept[..n].to_vec(), pis: kept[n..].to_vec(), form: p.form }, route, st))
}

/// The inner provider plus Leaves for closure calls without e_p (before or after the inner
/// provider); every other call goes to the inner provider alone, so its failure reports are
/// unchanged.
pub struct ClosureLeaves<F: Field> {
    pub inner: Box<dyn HomProvider<F>>,
    pub leaves: Leaves,
    pub leaves_first: bool,
    /// B4w: pair calls the inner provider fails go to `leaf_pairs` (Leaves last)
    pub leaf_pairs: bool,
    last_fail: Option<SearchStats>,
}

impl<F: Field> ClosureLeaves<F> {
    pub fn new(inner: Box<dyn HomProvider<F>>, leaves: Leaves, leaves_first: bool, leaf_pairs: bool) -> Self {
        ClosureLeaves { inner, leaves, leaves_first, leaf_pairs, last_fail: None }
    }
}

impl<F: Field> HomProvider<F> for ClosureLeaves<F> {
    fn name(&self) -> &'static str {
        "ClosureLeaves"
    }
    fn full(&mut self, ar: &mut Arena<F>, s: &mut Schur, src: ExprId, tgt: ExprId, d: i32, b: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        let r = self.inner.full(ar, s, src, tgt, d, b);
        self.last_fail = if r.is_err() { self.inner.failed_stats() } else { None };
        r
    }
    fn pairs(&mut self, ar: &mut Arena<F>, s: &mut Schur, x: ExprId, q: Tid, lq: ExprId, eq: Option<MapId>, sdeg: i32, mult: usize, b: &Budget) -> Result<(Pairs<F>, SearchStats), Vec<Missing>> {
        let r = self.inner.pairs(ar, s, x, q, lq, eq, sdeg, mult, b);
        self.last_fail = if r.is_err() { self.inner.failed_stats() } else { None };
        let errs = match r {
            Ok(v) => return Ok(v),
            Err(e) if !self.leaf_pairs || (eq.is_some() && !self.leaves.eq_pairs) => return Err(e),
            Err(e) => e,
        };
        let mut acc = self.last_fail.clone().unwrap_or_default();
        let t0 = Instant::now();
        let r = leaf_pairs(&self.leaves, ar, s, x, q, lq, eq, sdeg, mult);
        let secs = t0.elapsed().as_secs_f64();
        match r {
            Ok((p, route, st)) => {
                if std::env::var_os("KRHOM_TRACE").is_some() {
                    eprintln!("[leafpairs] q {} s {} mult {}: certified by {:?} in {:.1} s (inner: {})", s.table(q).fmt(), sdeg, mult, route, secs, errs.first().map(|m| m.what.as_str()).unwrap_or(""));
                }
                let ss = SearchStats { candidates: p.pis.len(), evaluated: p.pis.len(), secs, exact_calls: st.exact_calls, exact_secs: st.exact_secs, detect: if route == LeafPairRoute::RexClass { "leaves-rex" } else { "leaves-lemmaA" }, ..SearchStats::default() };
                acc.absorb(&ss);
                self.last_fail = None;
                Ok((p, acc))
            }
            Err(e) => {
                if std::env::var_os("KRHOM_TRACE").is_some() {
                    eprintln!("[leafpairs] q {} s {} mult {}: refused in {:.1} s: {}", s.table(q).fmt(), sdeg, mult, secs, e);
                }
                acc.absorb(&SearchStats { secs, ..SearchStats::default() });
                self.last_fail = Some(acc);
                let mut errs = errs;
                errs.iter_mut().for_each(|m| m.what = format!("[inner] {}", m.what));
                errs.push(Missing { what: format!("[Leaves] {}", e), degree: sdeg, stratum: None, got: 0, want: mult, kind: MissKind::Refused });
                Err(errs)
            }
        }
    }
    fn end0(&mut self, ar: &mut Arena<F>, s: &mut Schur, l: ExprId, b: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        let r = self.inner.end0(ar, s, l, b);
        self.last_fail = if r.is_err() { self.inner.failed_stats() } else { None };
        r
    }
    fn closure(&mut self, ar: &mut Arena<F>, s: &mut Schur, lp: ExprId, p: Tid, ep: Option<MapId>, outs: bool, b: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        self.last_fail = None;
        if ep.is_some() {
            let r = self.inner.closure(ar, s, lp, p, ep, outs, b);
            self.last_fail = if r.is_err() { self.inner.failed_stats() } else { None };
            return r;
        }
        let mut errs: Vec<Missing> = vec![];
        let mut acc = SearchStats::default();
        for k in 0..2 {
            let use_leaves = (k == 0) == self.leaves_first;
            let r = if use_leaves { self.leaves.closure_leaves(ar, s, lp, p, outs) } else { self.inner.closure(ar, s, lp, p, ep, outs, b) };
            match r {
                Ok((v, st)) => {
                    let mut st2 = acc.clone();
                    st2.absorb(&st);
                    return Ok((v, st2));
                }
                Err(ms) => {
                    let fs = if use_leaves { HomProvider::<F>::failed_stats(&self.leaves) } else { self.inner.failed_stats() };
                    if let Some(st) = fs {
                        acc.absorb(&st);
                    }
                    let name = if use_leaves { "Leaves" } else { "inner" };
                    errs.extend(ms.into_iter().map(|mut m| {
                        m.what = format!("[{}] {}", name, m.what);
                        m
                    }));
                }
            }
        }
        self.last_fail = Some(acc);
        Err(errs)
    }
    fn failed_stats(&self) -> Option<SearchStats> {
        self.last_fail.clone()
    }
}

// =============================================================================================
// M4c B2/B4: rex-class pairs (X ≅ B_q: X is itself a reduced expression of q, s = 0, mult 1)
// =============================================================================================

/// Report of a rex-class pair built from rex chains.
#[derive(Clone, Debug)]
pub struct RexPairReport {
    pub pi_chain: Vec<String>,
    pub iota_chain: Vec<String>,
    pub templates: usize,
    pub exact_moves: usize,
    pub max_exact_unknowns: u128,
    pub theorem_assoc: Vec<MapId>,
    pub secs: f64,
}

/// π: X → L_q and ι: L_q → X as realized rex chains (src/library/rexplan.rs, `leafgen::realize_chain`),
/// certified by `certify_pairs` (form rank 1 at rep(q) at two points, genuine) and checked:
/// π∘ι = id_{L_q} at every coset of the support of L_q (End⁰(B_q) = k and the top scalar is 1).
/// Every map that uses the theorem-level Assoc certificate is listed in `theorem_assoc`.
pub fn rex_class_pairs<F: Field>(ar: &mut Arena<F>, s: &mut Schur, x: ExprId, q: Tid, lq: ExprId, opts: &crate::library::rexplan::PlanOpts, seed: u64) -> Result<(Pairs<F>, RexPairReport), String> {
    rex_class_pairs_opts(ar, s, x, q, lq, opts, seed, false)
}

/// `rex_class_pairs` with `persistable` (`leafgen::realize_chain_opts`).
#[allow(clippy::too_many_arguments)]
pub fn rex_class_pairs_opts<F: Field>(ar: &mut Arena<F>, s: &mut Schur, x: ExprId, q: Tid, lq: ExprId, opts: &crate::library::rexplan::PlanOpts, seed: u64, persistable: bool) -> Result<(Pairs<F>, RexPairReport), String> {
    let t0 = Instant::now();
    let lam = ar.lam().to_vec();
    let (xw, lw) = (ar.expr_of(x).rungs.clone(), ar.expr_of(lq).rungs.clone());
    if !crate::combinat::dcoset::reduced_for(s, &xw, q).ok() {
        return Err("rex_class_pairs: X is not a reduced expression of q".into());
    }
    let mut pl = crate::library::rexplan::Planner::new();
    let ch_pi = pl.plan(s, q, &xw, &lw, opts)?;
    let ch_io = pl.plan(s, q, &lw, &xw, opts)?;
    let mut cache = HashMap::new();
    let (pi, rp) = crate::library::leafgen::realize_chain_opts(ar, &ch_pi, &mut cache, seed, persistable)?;
    let (io, ri) = crate::library::leafgen::realize_chain_opts(ar, &ch_io, &mut cache, seed + 1000, persistable)?;
    let c = rep_rank_of(s, &lam, q);
    let pairs = certify_pairs(ar, x, lq, 0, &[io], &[pi], c, 1, seed ^ 0xb4)?;
    // π∘ι = id on L_q at every coset of its support (two points)
    let comp = ar.add(MapNode::Comp(pi, io), "rex pair pi.iota")?;
    for pt in points::<F>(ar.k(), seed ^ 0x1d, 2, 0) {
        let mut ev = Eval::new(ar, pt);
        let lay = ev.layout(lq, lq);
        for &cc in &lay.cosets {
            let b = ev.block(comp, cc);
            let n = b.rows as usize;
            for i in 0..n {
                for j in 0..n {
                    let want = if i == j { F::one() } else { F::zero() };
                    if b.a[i * n + j] != want {
                        return Err(format!("rex_class_pairs: pi.iota != id at coset {}", cc));
                    }
                }
            }
        }
    }
    let mut thm = crate::library::leafgen::relies_on_theorem(ar, pi);
    thm.extend(crate::library::leafgen::relies_on_theorem(ar, io));
    let exact: Vec<&crate::library::rexplan::Move> = ch_pi.moves.iter().chain(&ch_io.moves).filter(|m| m.kind == crate::library::rexplan::MoveKind::Exact).collect();
    Ok((
        pairs,
        RexPairReport {
            pi_chain: rp.lines,
            iota_chain: ri.lines,
            templates: rp.templates + ri.templates,
            exact_moves: exact.len(),
            max_exact_unknowns: exact.iter().map(|m| m.unknowns).max().unwrap_or(0),
            theorem_assoc: thm,
            secs: t0.elapsed().as_secs_f64(),
        },
    ))
}

/// M4c B2/B4: pairs through Lemma A (A4 form) with a rex-chain exit/entry: π = exit ∘ LL_q
/// (Out leaves X → W_q of degree −s, then the chain W_q → L_q), ι = LL^∨_q ∘ entry (the chain
/// L_q → W_q, then In leaves W_q → X of degree s); `certify_pairs` selects `mult` pairs with an
/// invertible form at rep(q) (two points, genuine). W_q is the final junction word of
/// `plan_words`; when it differs from L_q the chains come from src/library/rexplan.rs.
#[allow(clippy::too_many_arguments)]
pub fn pairs_leaves<F: Field>(lv: &Leaves, ar: &mut Arena<F>, s: &mut Schur, x: ExprId, q: Tid, lq: ExprId, sdeg: i32, mult: usize, opts: &crate::library::rexplan::PlanOpts, st: &mut LeafStats) -> Result<(Pairs<F>, Vec<MapId>), String> {
    pairs_leaves_opts(lv, ar, s, x, q, lq, sdeg, mult, opts, st, false)
}

/// `pairs_leaves` with `persistable` (`leafgen::realize_chain_opts`).
#[allow(clippy::too_many_arguments)]
pub fn pairs_leaves_opts<F: Field>(lv: &Leaves, ar: &mut Arena<F>, s: &mut Schur, x: ExprId, q: Tid, lq: ExprId, sdeg: i32, mult: usize, opts: &crate::library::rexplan::PlanOpts, st: &mut LeafStats, persistable: bool) -> Result<(Pairs<F>, Vec<MapId>), String> {
    pairs_leaves_eq(lv, ar, s, x, q, lq, None, sdeg, mult, opts, st, persistable)
}

/// `pairs_leaves_opts` for a target with an idempotent e_q (C2-7, `Leaves::eq_pairs`): the
/// Lemma A leaves and exit / entry chains as for an indecomposable target, then ι ↦ ι∘e_q,
/// π ↦ e_q∘π (as the production backends' `with_idem`) before `certify_pairs`. e_q is the
/// identity at rep(q) (the extra summands of L_q have types < q), so the form at rep(q) is that
/// of the unwrapped leaves; acceptance is unchanged (certify_pairs here, then sindec's checks with
/// e_q on the diagonal of Π∘Ι and Σιπ = e).
#[allow(clippy::too_many_arguments)]
pub fn pairs_leaves_eq<F: Field>(lv: &Leaves, ar: &mut Arena<F>, s: &mut Schur, x: ExprId, q: Tid, lq: ExprId, eq: Option<MapId>, sdeg: i32, mult: usize, opts: &crate::library::rexplan::PlanOpts, st: &mut LeafStats, persistable: bool) -> Result<(Pairs<F>, Vec<MapId>), String> {
    let lam = ar.lam().to_vec();
    let lw = ar.expr_of(lq).rungs.clone();
    let fo = lv.build_leaves(ar, s, x, &[q], LeafDir::Out, st).map_err(|e| e.what)?;
    let fi = lv.build_leaves(ar, s, x, &[q], LeafDir::In, st).map_err(|e| e.what)?;
    let fo = fo.into_iter().find(|f| f.t == q).ok_or("pairs_leaves: no Out family at q")?;
    let fi = fi.into_iter().find(|f| f.t == q).ok_or("pairs_leaves: no In family at q")?;
    let mut pis: Vec<MapId> = fo.maps.iter().copied().filter(|&m| ar.info(m).deg == -sdeg).collect();
    let mut ios: Vec<MapId> = fi.maps.iter().copied().filter(|&m| ar.info(m).deg == sdeg).collect();
    if fo.word != lw || fi.word != lw {
        let mut pl = crate::library::rexplan::Planner::new();
        let mut cache = HashMap::new();
        if fo.word != lw {
            let ch = pl.plan(s, q, &fo.word, &lw, opts)?;
            let (exit, _) = crate::library::leafgen::realize_chain_opts(ar, &ch, &mut cache, lv.seed ^ 0xe1, persistable)?;
            pis = pis.into_iter().map(|m| ar.add(MapNode::Comp(exit, m), "leaves pi exit")).collect::<Result<_, _>>()?;
        }
        if fi.word != lw {
            let ch = pl.plan(s, q, &lw, &fi.word, opts)?;
            let (entry, _) = crate::library::leafgen::realize_chain_opts(ar, &ch, &mut cache, lv.seed ^ 0xe2, persistable)?;
            ios = ios.into_iter().map(|m| ar.add(MapNode::Comp(m, entry), "leaves iota entry")).collect::<Result<_, _>>()?;
        }
    }
    if let Some(e) = eq {
        ios = ios.into_iter().map(|m| ar.add(MapNode::Comp(m, e), "leaves iota e_q")).collect::<Result<_, _>>()?;
        pis = pis.into_iter().map(|m| ar.add(MapNode::Comp(e, m), "leaves e_q pi")).collect::<Result<_, _>>()?;
    }
    let c = rep_rank_of(s, &lam, q);
    let p = certify_pairs(ar, x, lq, sdeg, &ios, &pis, c, mult, lv.seed ^ 0xa4)?;
    let mut thm = vec![];
    for &m in p.iotas.iter().chain(&p.pis) {
        thm.extend(crate::library::leafgen::relies_on_theorem(ar, m));
    }
    Ok((p, thm))
}
