//! src/leaves.rs — M4c light-leaves provider. Milestone A1: Stage 0 only (M4C_PLAN §3.2, D2 §2.1,
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
use crate::field::Field;
use crate::homprov::*;
use crate::schur::{Lp, Schur, Table, Tid};
use crate::sexpr::Expr;
use crate::smap::{Arena, ExprId, MapId, MapNode};
use crate::tmpl::Window;
use crate::web::Rung;
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
            lv = crate::dcoset::apply_labels(&lv, &word[k - 1]).expect("product_factorization: invalid word");
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
            let canonical = crate::sindec::canonical_ladder(s, a) == aw && crate::sindec::canonical_ladder(s, b) == bw;
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
    crate::smapio::rungs_str(w)
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
