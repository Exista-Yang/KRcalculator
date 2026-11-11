//! M6.0: the exact small colored engine (SYNTHESIS §9.2; design note
//! `colored_design/M6_0_design.md`).
//!
//! Colored (∧^a per strand) gl_N homology of 2-strand braid closures through the singular Soergel
//! path, exactly over a prime field:
//!
//! * **Objects** are (type p, homological degree h, q-shift t), p a table of the block
//!   (λ, current labels) of `schur::Schur(λ)`. On two uprights every canonical basis element is
//!   the character of a ladder with at most two rungs (M3a test T-h2), so B_p is represented by
//!   that ladder L_p itself (idempotent e_p = id). The object (p, h, t) is q^t L_p, K_0 class
//!   v^t b_p.
//! * **Maps** are exact `PMat`s over R^λ between the ladders. A component B_p{t} → B_q{t'} has
//!   normalized degree e = t − t' (SYNTHESIS §2.3); every component is checked to be homogeneous
//!   of that degree.
//! * **Crossings.** C ⊗ X with X the untruncated colored Rickard complex (`cube::local_cx` with
//!   n = c + d): d_C ⊗ id = `PMat::kron_id`, id ⊗ d = `cube::apply_local_poly` with the Koszul sign
//!   (−1)^h. Each L_p ⊗ W_k is decomposed into canonical ladders B_q{s} with multiplicities from
//!   `schur::term_mult`: inclusions ι and projections π from `sbim::hom_space`, scalar
//!   intersection forms π∘ι ∈ End^0(L_q) = k (asserted to be scalars, of rank [v^s] m_q), summand
//!   ranks adding up to the rank of the tensor product (asserted), and Π∘Ι = id across summands
//!   (unitriangular, inverted exactly by a Neumann series; asserted). By default the two rungs of
//!   W_k are decomposed one at a time (L_p ⊗ R_1 = ⊕ L_r{s_1}, then L_r ⊗ R_2 = ⊕ L_q{s_2}, with
//!   ι = (ι_1 ⊗ id) ∘ ι_2, π = π_2 ∘ (π_1 ⊗ id)): the Hom spaces are then between ladders with at
//!   most three rungs instead of four (fibres up to 90 instead of 900 for colors 3), and the
//!   resulting multiplicities are asserted to equal `term_mult(p, W_k)`. A tensor product that is
//!   literally the canonical ladder of its only summand (multiplicity 1) needs no Hom space
//!   (ι = π = id). `XOpts::rungwise = false` decomposes L_p ⊗ W_k in one step (cross-check).
//!   Then Gaussian elimination with scalar pivots in Markowitz order. The whole tensor complex is
//!   built before eliminating, so every pivot is eligible (the streaming rule of SYNTHESIS §1.4 is
//!   not needed).
//! * **Invariants asserted per crossing:** d∘d = 0 (option), the K_0 class equals the
//!   decategorified scan (`schur` `mul_rungs`/`expand`, as `schur::decat_scan`), and the number of
//!   objects is at least the Grothendieck lower bound.
//! * **Closure data, exact:** per final type, R^λ-bases of Hom(R^λ, L_p) and Hom(L_p, R^λ) by
//!   exact graded Nakayama (`cclose::rlambda_bases`), checked by the one-point determinant
//!   criterion at two seeded points and against the degrees of c_diag(b_p); Gram by composition;
//!   F per final component by an exact solve (`cclose::solve_f`). The closure is the CL1 exact
//!   oracle `cclose::close_exact_colored` per N (unreduced, or reduced with block 0 marked) with
//!   `cclose::global_q`. One scan serves every N; labels > N give 0.
//!
//! **Pruning (`XOpts::prune_n = Some(N)`, a scan for gl_N' with N' <= N only).** A type q is
//! killed at N when `schur::Schur::killed_at(q, N)` (cell weight ν(q)_1 > N, the rule documented
//! at `Schur::cell_weight`): B_q maps to 0 in Foam_N. The killed types span a two-sided ideal I_N
//! (the killed b_q span one, with positive structure constants), so the scan can run in the
//! quotient SSBim / I_N, which the gl_N closure factors through when the rule is right:
//! * after every decomposition L_p ⊗ (rungs) the killed summands are dropped (their ι, π, i.e.
//!   their Hom spaces, are never computed; rung by rung, a killed intermediate summand is not
//!   tensored further); the kept summands keep Π∘Ι = id among themselves, which is a
//!   decomposition modulo I_N; ranks are asserted as Σ kept + Σ killed (`schur::rank`);
//! * the Rickard complexes are built by `crossing_cx(.., N)`, which drops the terms W_k with a
//!   label > N. This is consistent: W_k = B_{p_k} with ν(p_k) = (max(c,d) + k, min(c,d) − k),
//!   so the dropped terms are exactly the killed ones (asserted per crossing), and the kept
//!   terms W_0..W_kmax with their differentials are a stupid truncation of the untruncated
//!   complex, isomorphic to it modulo I_N;
//! * the K_0 assertion compares with the decategorified class with the killed types removed
//!   (the projection along the killed span); the lower bound is that of the truncated class;
//! * d∘d is no longer 0 on the nose: dropping a killed y between kept x and z leaves
//!   −Σ_y d_{yz} d_{xy} ∈ I_N. With `check_d2` the nonzero components are counted per stage
//!   (`StageInfo::d2_ideal`), and those of the final complex are kept (`XScan::final_d2`) and
//!   checked to close to 0 at every requested N before the closure (`closure_kills_map`);
//! * the closure is valid for every N' <= N (I_N ⊆ I_N'), and `XScan::close` refuses N' > N.

use crate::cclose::{close_exact_colored, det_criterion, gram_of, rlambda_bases, solve_f, ClosureDataC, EPoly, LadderBasis, TypeBasisC};
use crate::crossing::Gen;
use crate::cube::{apply_local_poly, local_cx, LocalCx};
use crate::field::Field;
use crate::laurent::Poincare;
use crate::mpoly::{MonoCache, Poly};
use crate::pmat::{PAcc, PMat, PVec};
use crate::sbim::{hom_space, SBim};
use crate::schur::{crossing_terms, flows, Elt, Lp, Schur, Tid};
use crate::web::{Dir, Rung};
use rand::{rngs::StdRng, Rng, SeedableRng};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::rc::Rc;
use std::time::Instant;

/// Largest label accepted by the exact engine.
pub const MAX_LABEL: usize = 3;

/// Options of a scan.
#[derive(Clone, Debug)]
pub struct XOpts {
    /// assert d∘d = 0 on the minimal complex after every crossing
    pub check_d2: bool,
    /// per-stage statistics on stderr (default: the environment variable KRHOM_VERBOSE)
    pub verbose: bool,
    /// seed of the points of the determinant criterion
    pub seed: u64,
    /// decompose L_p ⊗ W_k rung by rung (L_p ⊗ R_1 first, then each summand ⊗ R_2); false:
    /// decompose L_p ⊗ W_k in one step (larger Hom spaces; kept as a cross-check)
    pub rungwise: bool,
    /// Foam_N pruning for a scan that is closed only at gl_N' with N' <= N (module doc,
    /// "Pruning"): every summand B_q with `Schur::killed_at(q, N)` is dropped right after each
    /// decomposition (its ι, π are never computed), and the Rickard complexes are built by
    /// `crossing_cx(.., N)` (terms with a label > N dropped). None: the N-independent scan.
    pub prune_n: Option<usize>,
}

impl Default for XOpts {
    fn default() -> Self {
        let verbose = std::env::var("KRHOM_VERBOSE").map(|v| v != "0" && !v.is_empty()).unwrap_or(false);
        XOpts { check_d2: true, verbose, seed: 7, rungwise: true, prune_n: None }
    }
}

/// Statistics of one crossing step.
#[derive(Clone, Debug)]
pub struct StageInfo {
    /// number of crossings tensored so far
    pub step: usize,
    /// labels at the top of the stage
    pub labels: Vec<usize>,
    /// objects of C ⊗ X after decomposition, before elimination
    pub unsimplified: usize,
    /// objects and components of the minimal complex
    pub objects: usize,
    pub components: usize,
    /// Grothendieck lower bound Σ |signed multiplicities| of the K_0 class (pruned scan: of
    /// the truncated class, killed types removed)
    pub lower_bound: u128,
    pub secs: f64,
    /// pruned scan: summands of C ⊗ X killed at N and dropped (counted with multiplicity)
    pub dropped: usize,
    /// pruned scan: nonzero components of d∘d on the minimal complex (maps through dropped
    /// summands, zero modulo the ideal of killed types; always 0 for an unpruned scan, where
    /// d∘d = 0 is asserted)
    pub d2_ideal: usize,
}

/// Result of a scan: N-independent closure data and statistics.
pub struct XScan<F: Field> {
    pub labels: Vec<usize>,
    pub word: Vec<Gen>,
    /// contract-v2 closure data (type ids = schur table ids, object q = t, normalized degrees)
    pub data: ClosureDataC<F>,
    pub stages: Vec<StageInfo>,
    /// final objects (table, h, t), sorted
    pub final_objects: Vec<(String, i32, i32)>,
    /// K_0 class of the final complex: (table, t, Σ_h (−1)^h multiplicity) (pruned scan: the
    /// decategorified class with the killed types removed, as asserted after every crossing)
    pub k0: Vec<(String, i32, i128)>,
    /// `hom_space` calls of the scan (decompositions): count, total and maximal seconds
    pub hom_calls: usize,
    pub hom_total_secs: f64,
    pub hom_max_secs: f64,
    /// largest fibre (rank over R^λ) of a tensor product that was decomposed by Hom spaces
    pub max_x_fsize: usize,
    /// number of tensor products decomposed by Hom spaces and intersection forms
    pub decompositions: usize,
    /// number of tensor products that are literally the canonical ladder of their only summand
    pub identity_decompositions: usize,
    /// every type used by the scan: (table, canonical ladder, rank over R^λ)
    pub ladders: Vec<(String, String, usize)>,
    /// canonical ladder (rungs, bottom to top) of every type of the closure data
    pub type_ladders: BTreeMap<u32, Vec<Rung>>,
    /// canonical ladder (rungs, bottom to top) of every type used by the scan (intermediate
    /// stages and decomposition summands included), keyed by table; bottom labels = `labels`
    pub scan_ladders: BTreeMap<String, Vec<Rung>>,
    /// seconds: (scan, closure data)
    pub timings: (f64, f64),
    /// the N of a pruned scan (`XOpts::prune_n`); `close` then accepts only N' <= N
    pub prune_n: Option<usize>,
    /// pruned scan: summands dropped in the decompositions computed (with multiplicity), and
    /// decompositions whose summands were all killed (no Hom space computed)
    pub dropped_summands: usize,
    pub all_killed_decompositions: usize,
    /// pruned scan: the nonzero components of d∘d on the final complex, as closure data F
    /// (source index, target index into `data.objs`, F[c][a]); `close` checks that the closure
    /// at every requested N kills each of them (so the closed complex is a complex)
    pub final_d2: Vec<(usize, usize, Vec<Vec<EPoly<F>>>)>,
}

// ------------------------------------------------------------------ internal structures

/// Canonical ladder L_p of a type (= B_p on two uprights).
struct Lad<F: Field> {
    rungs: Vec<Rung>,
    b: SBim<F>,
}

/// One summand B_q{s} of L_p ⊗ W_k: ι: L_q → X of degree s, π: X → L_q of degree −s, with
/// Π∘Ι = id over all summands of X.
struct Summand<F: Field> {
    q: Tid,
    s: i32,
    iota: PMat<F>,
    pi: PMat<F>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct XObj {
    p: Tid,
    h: i32,
    t: i32,
}

/// A complex of canonical ladders with exact components.
struct XCx<F: Field> {
    objs: BTreeMap<usize, XObj>,
    out: BTreeMap<usize, BTreeMap<usize, PMat<F>>>,
    inn: BTreeMap<usize, BTreeSet<usize>>,
    next: usize,
}

impl<F: Field> XCx<F> {
    fn new() -> Self {
        XCx { objs: BTreeMap::new(), out: BTreeMap::new(), inn: BTreeMap::new(), next: 0 }
    }
    fn add_obj(&mut self, o: XObj) -> usize {
        let id = self.next;
        self.next += 1;
        self.objs.insert(id, o);
        self.out.insert(id, BTreeMap::new());
        self.inn.insert(id, BTreeSet::new());
        id
    }
    /// x → y += m (zero results are removed)
    fn add_map(&mut self, x: usize, y: usize, m: PMat<F>) {
        if m.is_zero() {
            return;
        }
        let outs = self.out.get_mut(&x).unwrap();
        let sum = match outs.remove(&y) {
            Some(old) => old.add(&m),
            None => m,
        };
        if sum.is_zero() {
            self.inn.get_mut(&y).unwrap().remove(&x);
        } else {
            outs.insert(y, sum);
            self.inn.get_mut(&y).unwrap().insert(x);
        }
    }
    fn remove_obj(&mut self, a: usize) {
        if let Some(outs) = self.out.remove(&a) {
            for y in outs.keys() {
                self.inn.get_mut(y).unwrap().remove(&a);
            }
        }
        if let Some(ins) = self.inn.remove(&a) {
            for x in ins {
                self.out.get_mut(&x).unwrap().remove(&a);
            }
        }
        self.objs.remove(&a);
    }
    fn components(&self) -> usize {
        self.out.values().map(|m| m.len()).sum()
    }
    /// signed multiset (type, t) -> Σ (−1)^h
    fn k0(&self) -> BTreeMap<(Tid, i32), i128> {
        let mut m: BTreeMap<(Tid, i32), i128> = BTreeMap::new();
        for o in self.objs.values() {
            *m.entry((o.p, o.t)).or_insert(0) += if o.h.rem_euclid(2) == 0 { 1 } else { -1 };
        }
        m.retain(|_, c| *c != 0);
        m
    }
}

/// Some(c) if m = c·id (c may be zero; a zero matrix of any shape gives Some(0)), None otherwise.
fn scalar_of<F: Field>(m: &PMat<F>) -> Option<F> {
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
            Some(x) => {
                if *x != p.t[0].1 {
                    return None;
                }
            }
        }
    }
    c
}

fn scale_const<F: Field>(m: &PMat<F>, c: &F) -> PMat<F> {
    m.scale(&Poly::constant(c.clone()))
}

/// Is every entry of m: src → tgt homogeneous of normalized degree e?
fn homogeneous<F: Field>(m: &PMat<F>, src: &SBim<F>, tgt: &SBim<F>, e: i32) -> bool {
    let u = e - tgt.shift + src.shift;
    m.cols.iter().enumerate().all(|(be, col)| col.iter().all(|(bp, p)| p.t.iter().all(|(mono, _)| src.ring.deg(*mono) == u + src.fdeg[be] - tgt.fdeg[*bp as usize])))
}

/// Block matrix of maps between summands: entry [tau][sigma]: L_sigma → L_tau.
type Blocks<F> = Vec<Vec<Option<PMat<F>>>>;

fn blocks_mul<F: Field>(a: &Blocks<F>, b: &Blocks<F>) -> Blocks<F> {
    let n = a.len();
    let mut out: Blocks<F> = vec![vec![None; n]; n];
    for t in 0..n {
        for s in 0..n {
            let mut acc: Option<PMat<F>> = None;
            for r in 0..n {
                if let (Some(x), Some(y)) = (&a[t][r], &b[r][s]) {
                    let z = x.matmul(y);
                    acc = Some(match acc {
                        None => z,
                        Some(w) => w.add(&z),
                    });
                }
            }
            out[t][s] = acc.filter(|m| !m.is_zero());
        }
    }
    out
}

/// id_{L_p} ⊗ D for a local map D: W_s → W_t on the two uprights (0, 1) at the top of P
/// (fibre index of P ⊗ W = β_P · fsize(W) + α).
fn id_tensor_local<F: Field>(p: &SBim<F>, lmap: &PMat<F>, c: usize, d: usize, fs: usize, ft: usize) -> PMat<F> {
    let fp = p.fsize();
    let mut cols: Vec<PVec<F>> = Vec::with_capacity(fp * fs);
    for beta in 0..fp * fs {
        let (bp, bl) = (beta / fs, beta % fs);
        let unit: PVec<F> = vec![(bp as u32, Poly::one())];
        let mut acc = PAcc::new();
        for (blp, lpoly) in &lmap.cols[bl] {
            let w = apply_local_poly(p, 0, c, d, lpoly, &unit);
            let placed: PVec<F> = w.into_iter().map(|(x, poly)| (x * ft as u32 + *blp, poly)).collect();
            acc.add_scaled_const(&placed, &F::one());
        }
        cols.push(acc.take());
    }
    PMat { nrows: fp * ft, ncols: fp * fs, cols }
}

/// Are a and b the same SBim in the same fibre basis: labels, rungs, fibre degrees, shift and
/// the matrices of every top operator?
fn same_sbim<F: Field>(a: &SBim<F>, b: &SBim<F>) -> bool {
    a.bot_labels == b.bot_labels
        && a.top_labels == b.top_labels
        && a.rungs == b.rungs
        && a.fsize() == b.fsize()
        && a.fdeg == b.fdeg
        && a.shift == b.shift
        && a.top.len() == b.top.len()
        && a.top.iter().zip(&b.top).all(|(x, y)| x.len() == y.len() && x.iter().zip(y).all(|(m, n)| m.nrows == n.nrows && m.ncols == n.ncols && m.sub(n).is_zero()))
}

fn rung_key(r: &Rung) -> (u8, usize, usize) {
    (if r.dir == Dir::E { 0 } else { 1 }, r.i, r.r)
}

// ------------------------------------------------------------------ the engine

struct XEngine<F: Field> {
    lam: Vec<usize>,
    schur: Schur,
    base: SBim<F>,
    mc: MonoCache,
    lads: HashMap<Tid, Rc<Lad<F>>>,
    locals: HashMap<(usize, usize, bool), Rc<LocalCx<F>>>,
    decomp: HashMap<(Tid, bool, usize), Rc<Vec<Summand<F>>>>,
    wdecomp: HashMap<(Tid, Vec<(u8, usize, usize)>), Rc<Vec<Summand<F>>>>,
    /// decompose crossing terms rung by rung (default) or in one step
    rungwise: bool,
    /// re-check Π∘Ι = id on composite decompositions
    check_orth: bool,
    direct_decomps: usize,
    identity_decomps: usize,
    dloc: HashMap<(Tid, bool, usize), Rc<Vec<(usize, usize, PMat<F>)>>>,
    /// decategorified scan: homological degree -> element of the current block
    decat: BTreeMap<i32, Elt>,
    hom_calls: usize,
    hom_total: f64,
    hom_max: f64,
    max_x: usize,
    verbose: bool,
    /// Foam_N pruning (`XOpts::prune_n`)
    prune: Option<usize>,
    /// killed summands (with multiplicity) of L_p ⊗ W_k, per (p, positive, k)
    dropped: HashMap<(Tid, bool, usize), usize>,
    /// killed summands dropped in the current stage, and in the decompositions computed
    stage_dropped: usize,
    dropped_total: usize,
    all_killed_decomps: usize,
}

impl<F: Field> XEngine<F> {
    fn new(lam: &[usize], verbose: bool, prune: Option<usize>) -> Self {
        let lam8: Vec<u8> = lam.iter().map(|&a| a as u8).collect();
        let k: usize = lam.iter().sum();
        let mut schur = Schur::new(&lam8);
        let unit = schur.unit();
        let base = SBim::<F>::identity(2 * k, lam);
        let mc = MonoCache::new(base.ring.clone());
        XEngine {
            lam: lam.to_vec(),
            schur,
            base,
            mc,
            lads: HashMap::new(),
            locals: HashMap::new(),
            decomp: HashMap::new(),
            wdecomp: HashMap::new(),
            rungwise: true,
            check_orth: true,
            direct_decomps: 0,
            identity_decomps: 0,
            dloc: HashMap::new(),
            decat: BTreeMap::from([(0, unit)]),
            hom_calls: 0,
            hom_total: 0.0,
            hom_max: 0.0,
            max_x: 0,
            verbose,
            prune,
            dropped: HashMap::new(),
            stage_dropped: 0,
            dropped_total: 0,
            all_killed_decomps: 0,
        }
    }

    /// Is B_q zero in the gl_N theory of a pruned scan?
    fn killed(&self, q: Tid) -> bool {
        self.prune.is_some_and(|n| self.schur.killed_at(q, n))
    }

    /// Σ over the killed summands B_q{s} of L_p ⊗ (rungs) of rank(B_q), with multiplicity
    /// (0 for an unpruned scan): with the kept summands it adds up to the rank of L_p ⊗ (rungs).
    fn killed_rank(&mut self, p: Tid, rs: &[Rung]) -> usize {
        let Some(n) = self.prune else { return 0 };
        let mut r: u128 = 0;
        for (q, m) in self.schur.term_mult(p, rs) {
            if self.schur.killed_at(q, n) {
                let mu = self.schur.table(q).cols();
                r += m.at1() as u128 * self.schur.rank(&mu, q);
            }
        }
        r as usize
    }

    fn identity_type(&mut self) -> Tid {
        let u = self.schur.unit();
        assert_eq!(u.t.len(), 1);
        u.t[0].0
    }

    /// The canonical ladder of type p: a ladder with at most two rungs whose character is b_p
    /// (minimal flows, then fewest rungs, then lexicographic).
    fn ladder(&mut self, p: Tid) -> Rc<Lad<F>> {
        if let Some(l) = self.lads.get(&p) {
            return l.clone();
        }
        let mu8 = self.schur.table(p).cols();
        let mu: Vec<usize> = mu8.iter().map(|&a| a as usize).collect();
        let bp = self.schur.canonical_basis(&mu8).b[&p].clone();
        let k: usize = self.lam.iter().sum();
        let lab0 = self.lam.clone();
        let mut cands: Vec<Vec<Rung>> = vec![vec![]];
        if lab0.len() == 2 {
            for d1 in [Dir::E, Dir::F] {
                for r1 in 1..=k {
                    let a = Rung { i: 0, dir: d1, r: r1 };
                    let Some(l1) = a.apply(&lab0, k) else { continue };
                    cands.push(vec![a]);
                    for d2 in [Dir::E, Dir::F] {
                        for r2 in 1..=k {
                            let b = Rung { i: 0, dir: d2, r: r2 };
                            if b.apply(&l1, k).is_some() {
                                cands.push(vec![a, b]);
                            }
                        }
                    }
                }
            }
        }
        let lam8: Vec<u8> = self.lam.iter().map(|&a| a as u8).collect();
        let mut best: Option<((u128, usize, Vec<(u8, usize, usize)>), Vec<Rung>)> = None;
        for rs in cands {
            let mut lab = lab0.clone();
            for r in &rs {
                lab = r.apply(&lab, k).unwrap();
            }
            if lab != mu {
                continue;
            }
            let ch = self.schur.ladder_char(&rs);
            if ch != bp {
                continue;
            }
            let key = (flows(&lam8, &rs), rs.len(), rs.iter().map(rung_key).collect::<Vec<_>>());
            if best.as_ref().map_or(true, |(bk, _)| key < *bk) {
                best = Some((key, rs));
            }
        }
        let (_, rungs) = best.unwrap_or_else(|| panic!("type {}: b_p is not the character of a ladder with at most two rungs (T-h2)", self.schur.table(p).fmt()));
        let b = self.base.add_rungs(&rungs);
        let rk = self.schur.rank(&mu8, p);
        assert_eq!(b.fsize() as u128, rk, "type {}: fsize of the canonical ladder != rank of B_p", self.schur.table(p).fmt());
        let l = Rc::new(Lad { rungs, b });
        self.lads.insert(p, l.clone());
        l
    }

    /// The untruncated Rickard complex of a crossing on labels (c, d).
    fn local(&mut self, c: usize, d: usize, positive: bool) -> Rc<LocalCx<F>> {
        if let Some(l) = self.locals.get(&(c, d, positive)) {
            return l.clone();
        }
        // pruned scan: crossing_cx(.., N) drops the terms W_k whose web has a label > N; the
        // kept terms are W_0..W_kmax (the intermediate label grows with k) and their SBims, hence
        // the unique differentials between them, are those of the untruncated complex
        // (`SBim::identity(n, ..)` records n but neither add_rungs nor hom_space uses it)
        let n = self.prune.map_or(c + d, |nn| nn.min(c + d));
        let loc = local_cx::<F>(Gen { i: 0, positive }, c, d, n);
        let st = crossing_terms(c as u8, d as u8, positive, 0);
        assert!(loc.cx.terms.len() <= st.len() && !loc.cx.terms.is_empty());
        for (k, a) in loc.cx.terms.iter().enumerate() {
            let b = &st[k];
            assert_eq!(a.k, k, "pruned Rickard complex of ({},{}) is not W_0..W_kmax", c, d);
            assert_eq!(a.rungs, b.rungs, "crossing.rs and schur.rs disagree on term {} of ({},{})", k, c, d);
            assert_eq!(a.hdeg, b.hdeg);
        }
        match self.prune {
            None => assert_eq!(loc.cx.terms.len(), c.min(d) + 1, "Rickard complex pruned"),
            Some(nn) => {
                // the label pruning of crossing.rs is the cell rule: W_k = B_{p_k} is dropped iff
                // p_k is killed at N (schur::Schur::cell_weight(p_k) = (max(c,d) + k, min(c,d) − k))
                let mut sc = Schur::new(&[c as u8, d as u8]);
                for (k, t) in st.iter().enumerate() {
                    let ch = sc.ladder_char(&t.rungs);
                    let e = sc.expand(&ch);
                    assert!(e.len() == 1 && e[0].1 == Lp::one(), "Rickard term {} of ({},{}) is not indecomposable", k, c, d);
                    assert_eq!(sc.killed_at(e[0].0, nn), k >= loc.cx.terms.len(), "Rickard term {} of ({},{}): label pruning at N={} != cell rule (cell {:?})", k, c, d, nn, sc.cell_weight(e[0].0));
                }
            }
        }
        let l = Rc::new(loc);
        self.locals.insert((c, d, positive), l.clone());
        l
    }

    /// Hom space b → bp in normalized degree e (timed).
    fn hom_norm(&mut self, b: &SBim<F>, bp: &SBim<F>, e: i32) -> Vec<PMat<F>> {
        let u = e - bp.shift + b.shift;
        let t0 = Instant::now();
        let h = hom_space(b, bp, u, &mut self.mc);
        let dt = t0.elapsed().as_secs_f64();
        self.hom_calls += 1;
        self.hom_total += dt;
        self.hom_max = self.hom_max.max(dt);
        if self.verbose && dt > 5.0 {
            eprintln!("[xengine]   hom_space fibres {} -> {} degree {}: dim {} in {:.1}s", b.fsize(), bp.fsize(), e, h.len(), dt);
        }
        h
    }

    /// Decomposition of L_p ⊗ W_k (term k of the crossing at the top of p), cached.
    fn decompose(&mut self, p: Tid, positive: bool, k: usize) -> Rc<Vec<Summand<F>>> {
        if let Some(d) = self.decomp.get(&(p, positive, k)) {
            return d.clone();
        }
        let lp = self.ladder(p);
        let mu = lp.b.top_labels.clone();
        let loc = self.local(mu[0], mu[1], positive);
        let wr = loc.cx.terms[k].rungs.clone();
        let r = self.decompose_word(p, &wr);
        // the multiplicities are those of the Hecke algebroid: b_p * [W_k] = Σ m_q b_q
        let mut got: BTreeMap<(Tid, i32), i128> = BTreeMap::new();
        for sm in r.iter() {
            *got.entry((sm.q, sm.s)).or_insert(0) += 1;
        }
        // (pruned scan: the kept part of b_p * [W_k]; the killed part is dropped)
        let mut want: BTreeMap<(Tid, i32), i128> = BTreeMap::new();
        let mut dropped = 0usize;
        for (q, m) in self.schur.term_mult(p, &wr) {
            if self.killed(q) {
                dropped += m.at1() as usize;
                continue;
            }
            for (i, &c) in m.c.iter().enumerate() {
                if c != 0 {
                    want.insert((q, m.lo + i as i32), c);
                }
            }
        }
        assert_eq!(got, want, "decomposition of {} ⊗ W_{}: summands != schur::term_mult", self.schur.table(p).fmt(), k);
        self.dropped.insert((p, positive, k), dropped);
        self.dropped_total += dropped;
        self.decomp.insert((p, positive, k), r.clone());
        r
    }

    /// Decomposition of L_p ⊗ (rungs), cached. Rung by rung (default): L_p ⊗ R_1 = ⊕ L_r{s_1},
    /// then L_r ⊗ (R_2 ...) recursively, with ι = (ι_1 ⊗ id) ∘ ι_2 and π = π_2 ∘ (π_1 ⊗ id); in the
    /// direct mode the whole word is decomposed at once.
    fn decompose_word(&mut self, p: Tid, rs: &[Rung]) -> Rc<Vec<Summand<F>>> {
        let key = (p, rs.iter().map(rung_key).collect::<Vec<_>>());
        if let Some(d) = self.wdecomp.get(&key) {
            return d.clone();
        }
        let r = if rs.is_empty() {
            debug_assert!(!self.killed(p));
            let n = self.ladder(p).b.fsize();
            Rc::new(vec![Summand { q: p, s: 0, iota: PMat::identity(n), pi: PMat::identity(n) }])
        } else if rs.len() == 1 || !self.rungwise {
            self.decompose_direct(p, rs)
        } else {
            let first = self.decompose_word(p, &rs[..1]);
            let lp = self.ladder(p);
            let f1 = lp.b.add_rungs(&rs[..1]);
            let x = lp.b.add_rungs(rs);
            let fr = x.fsize() / f1.fsize();
            let mut out = vec![];
            for s1 in first.iter() {
                let rest = self.decompose_word(s1.q, &rs[1..]);
                let i1 = s1.iota.kron_id(fr);
                let p1 = s1.pi.kron_id(fr);
                for s2 in rest.iter() {
                    out.push(Summand { q: s2.q, s: s1.s + s2.s, iota: i1.matmul(&s2.iota), pi: s2.pi.matmul(&p1) });
                }
            }
            // pruned scan: a killed summand B_r of L_p ⊗ R_1 has only killed summands in
            // B_r ⊗ R_2 ... (the killed span is a two-sided ideal with positive structure
            // constants), so the kept summands of L_p ⊗ (rungs) all come from kept B_r
            let total: usize = out.iter().map(|sm| self.lads[&sm.q].b.fsize()).sum();
            let killed = self.killed_rank(p, rs);
            assert_eq!(total + killed, x.fsize(), "rung-by-rung decomposition: summand ranks {} + killed {} != rank {}", total, killed, x.fsize());
            if self.prune.is_some() {
                let mut got: BTreeMap<(Tid, i32), i128> = BTreeMap::new();
                for sm in &out {
                    *got.entry((sm.q, sm.s)).or_insert(0) += 1;
                }
                let mut want: BTreeMap<(Tid, i32), i128> = BTreeMap::new();
                for (q, m) in self.schur.term_mult(p, rs) {
                    if !self.killed(q) {
                        for (i, &c) in m.c.iter().enumerate() {
                            if c != 0 {
                                want.insert((q, m.lo + i as i32), c);
                            }
                        }
                    }
                }
                assert_eq!(got, want, "pruned rung-by-rung decomposition of {}: kept summands != schur::term_mult", self.schur.table(p).fmt());
            }
            if self.check_orth {
                for (t, a) in out.iter().enumerate() {
                    for (u, b) in out.iter().enumerate() {
                        let c = scalar_of(&a.pi.matmul(&b.iota)).expect("Π∘Ι != id (composite)");
                        assert!(if t == u { c.is_one() } else { c.is_zero() }, "Π∘Ι != id (composite)");
                    }
                }
            }
            Rc::new(out)
        };
        self.wdecomp.insert(key, r.clone());
        r
    }

    /// Decomposition of X = L_p ⊗ (rungs) by Hom spaces and intersection forms.
    fn decompose_direct(&mut self, p: Tid, rs: &[Rung]) -> Rc<Vec<Summand<F>>> {
        let t0 = Instant::now();
        let lp = self.ladder(p);
        let mult = self.schur.term_mult(p, rs);
        // pruned scan: the killed summands are dropped; only the kept ones get ι, π (Hom
        // spaces). Any family of kept summands with Π∘Ι = id is a decomposition in the quotient
        // by the ideal of killed types: its complement is a summand with the killed multiset
        // (Krull-Schmidt), hence killed.
        let killed_rank = self.killed_rank(p, rs);
        let kept: Vec<(Tid, Lp)> = mult.iter().filter(|(q, _)| !self.killed(*q)).cloned().collect();
        if kept.is_empty() {
            let mu8: Vec<u8> = lp.b.top_labels.iter().map(|&a| a as u8).collect();
            let xsize = lp.b.fsize() as u128 * flows(&mu8, rs);
            assert_eq!(killed_rank as u128, xsize, "all summands killed: ranks {} != rank {}", killed_rank, xsize);
            self.all_killed_decomps += 1;
            return Rc::new(vec![]);
        }
        // X is literally the canonical ladder of its only summand B_q (multiplicity 1, shift 0):
        // ι = π = id, no Hom space needed
        if mult.len() == 1 && mult[0].1 == crate::schur::Lp::one() {
            let q = mult[0].0;
            let lq = self.ladder(q);
            let mut w = lp.rungs.clone();
            w.extend_from_slice(rs);
            if w == lq.rungs {
                // ι = π = id is only valid if X and L_q are the same singular Bott-Samelson
                // bimodule in the same fibre basis; this holds because both are built by
                // `SBim::add_rungs` from the same rung word, and is asserted here so that a change
                // of the add_rungs layout cannot silently corrupt the transported components
                let x = lp.b.add_rungs(rs);
                assert!(same_sbim(&x, &lq.b), "identity decomposition of {} ⊗ {:?}: L_p ⊗ R and the canonical ladder of {} differ as SBims", self.schur.table(p).fmt(), rs, self.schur.table(q).fmt());
                let n = lq.b.fsize();
                self.identity_decomps += 1;
                return Rc::new(vec![Summand { q, s: 0, iota: PMat::identity(n), pi: PMat::identity(n) }]);
            }
        }
        let x = lp.b.add_rungs(rs);
        self.max_x = self.max_x.max(x.fsize());
        let word = || rs.iter().map(|r| format!("{}{}({})", if r.dir == Dir::E { "E" } else { "F" }, r.i, r.r)).collect::<Vec<_>>().join(" ");
        let mut sums: Vec<Summand<F>> = vec![];
        for (q, m) in kept {
            let lq = self.ladder(q);
            assert_eq!(lq.b.top_labels, x.top_labels);
            for (i, &cnt) in m.c.iter().enumerate() {
                if cnt == 0 {
                    continue;
                }
                assert!(cnt > 0);
                let cnt = cnt as usize;
                let s = m.lo + i as i32;
                let iotas = self.hom_norm(&lq.b, &x, s);
                let pis = self.hom_norm(&x, &lq.b, -s);
                // intersection form G[j][i] = π_j ∘ ι_i ∈ End^0(L_q) = k
                let g: Vec<Vec<F>> = pis
                    .iter()
                    .map(|pj| {
                        iotas
                            .iter()
                            .map(|ii| {
                                let c = pj.matmul(ii);
                                scalar_of(&c).unwrap_or_else(|| panic!("π∘ι is not a scalar: End^0 of type {} is not k", self.schur.table(q).fmt()))
                            })
                            .collect()
                    })
                    .collect();
                let mut gr = g.clone();
                let icols = crate::dense::rref(&mut gr, iotas.len());
                assert_eq!(icols.len(), cnt, "decomposition of {} ⊗ [{}]: rank of the intersection form for {} at v^{} is {}, multiplicity {}", self.schur.table(p).fmt(), word(), self.schur.table(q).fmt(), s, icols.len(), cnt);
                let mut gt: Vec<Vec<F>> = (0..iotas.len()).map(|i| (0..pis.len()).map(|j| g[j][i].clone()).collect()).collect();
                let jrows = crate::dense::rref(&mut gt, pis.len());
                assert_eq!(jrows.len(), cnt);
                // A = G[J][I]^{-1}; π'_a = Σ_b A[a][b] π_{J[b]}
                let sub: Vec<Vec<F>> = jrows.iter().map(|&j| icols.iter().map(|&i| g[j][i].clone()).collect()).collect();
                let inv = invert(&sub).expect("intersection form submatrix not invertible");
                for a in 0..cnt {
                    let mut pi = PMat::zero(lq.b.fsize(), x.fsize());
                    for b in 0..cnt {
                        if !inv[a][b].is_zero() {
                            pi = pi.add(&scale_const(&pis[jrows[b]], &inv[a][b]));
                        }
                    }
                    sums.push(Summand { q, s, iota: iotas[icols[a]].clone(), pi });
                }
            }
        }
        let total: usize = sums.iter().map(|sm| self.lads[&sm.q].b.fsize()).sum();
        assert_eq!(total + killed_rank, x.fsize(), "decomposition of {} ⊗ [{}]: summand ranks {} + killed {} != rank {}", self.schur.table(p).fmt(), word(), total, killed_rank, x.fsize());
        // Π∘Ι = I + N with N nilpotent; Π <- (I + N)^{-1} Π
        let n = sums.len();
        let mut nb: Blocks<F> = vec![vec![None; n]; n];
        for t in 0..n {
            for s in 0..n {
                let m = sums[t].pi.matmul(&sums[s].iota);
                if t == s {
                    let c = scalar_of(&m).expect("π∘ι not scalar");
                    assert!(c.is_one(), "π∘ι != id on a summand");
                } else {
                    if !m.is_zero() {
                        assert!(sums[s].s > sums[t].s, "nonzero π_τ∘ι_σ with s_σ <= s_τ");
                    }
                    nb[t][s] = if m.is_zero() { None } else { Some(m.neg()) };
                }
            }
        }
        let mut minv: Blocks<F> = vec![vec![None; n]; n];
        for (t, row) in minv.iter_mut().enumerate() {
            row[t] = Some(PMat::identity(self.lads[&sums[t].q].b.fsize()));
        }
        let mut pw = minv.clone();
        let mut iter = 0;
        loop {
            pw = blocks_mul(&nb, &pw);
            if pw.iter().all(|r| r.iter().all(|e| e.is_none())) {
                break;
            }
            iter += 1;
            assert!(iter <= n, "Π∘Ι is not unipotent");
            for t in 0..n {
                for s in 0..n {
                    if let Some(z) = &pw[t][s] {
                        minv[t][s] = Some(match minv[t][s].take() {
                            None => z.clone(),
                            Some(w) => w.add(z),
                        })
                        .filter(|m| !m.is_zero());
                    }
                }
            }
        }
        if iter > 0 {
            let old: Vec<PMat<F>> = sums.iter().map(|sm| sm.pi.clone()).collect();
            for t in 0..n {
                let mut acc = PMat::zero(old[t].nrows, x.fsize());
                for s in 0..n {
                    if let Some(m) = &minv[t][s] {
                        acc = acc.add(&m.matmul(&old[s]));
                    }
                }
                sums[t].pi = acc;
            }
            for t in 0..n {
                for s in 0..n {
                    let m = sums[t].pi.matmul(&sums[s].iota);
                    let c = scalar_of(&m).expect("Π∘Ι != id after orthogonalization");
                    assert!(if t == s { c.is_one() } else { c.is_zero() }, "Π∘Ι != id after orthogonalization");
                }
            }
        }
        self.direct_decomps += 1;
        if self.verbose {
            eprintln!(
                "[xengine]   decomposed {} ⊗ [{}] (fibre {}) into {} summands in {:.2}s",
                self.schur.table(p).fmt(),
                word(),
                x.fsize(),
                n,
                t0.elapsed().as_secs_f64()
            );
        }
        Rc::new(sums)
    }

    /// Transported components of id_{L_p} ⊗ d_j between the summands, cached (no sign).
    fn dloc(&mut self, p: Tid, positive: bool, j: usize) -> Rc<Vec<(usize, usize, PMat<F>)>> {
        if let Some(d) = self.dloc.get(&(p, positive, j)) {
            return d.clone();
        }
        let lp = self.ladder(p);
        let (c, d) = (lp.b.top_labels[0], lp.b.top_labels[1]);
        let loc = self.local(c, d, positive);
        let (ks, kt) = if positive { (j, j + 1) } else { (j + 1, j) };
        let src = self.decompose(p, positive, ks);
        let tgt = self.decompose(p, positive, kt);
        let phi = id_tensor_local(&lp.b, &loc.maps[j], c, d, loc.fsizes[ks], loc.fsizes[kt]);
        let mut out = vec![];
        for (si, sg) in src.iter().enumerate() {
            let a = phi.matmul(&sg.iota);
            if a.is_zero() {
                continue;
            }
            for (ti, tu) in tgt.iter().enumerate() {
                let comp = tu.pi.matmul(&a);
                if !comp.is_zero() {
                    out.push((si, ti, comp));
                }
            }
        }
        let r = Rc::new(out);
        self.dloc.insert((p, positive, j), r.clone());
        r
    }

    /// C ⊗ (Rickard complex of the crossing at the top of C), decomposed and transported.
    fn tensor_crossing(&mut self, cx: &XCx<F>, labels: &[usize], positive: bool) -> XCx<F> {
        let (c, d) = (labels[0], labels[1]);
        let loc = self.local(c, d, positive);
        let nt = loc.cx.terms.len();
        let mut out = XCx::new();
        let mut place: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
        for (&o, ob) in &cx.objs {
            for k in 0..nt {
                let sums = self.decompose(ob.p, positive, k);
                self.stage_dropped += self.dropped[&(ob.p, positive, k)];
                let hk = loc.cx.terms[k].hdeg;
                let dq = if positive { -(k as i32) } else { k as i32 };
                let ids: Vec<usize> = sums.iter().map(|sm| out.add_obj(XObj { p: sm.q, h: ob.h + hk, t: ob.t + dq + sm.s })).collect();
                place.insert((o, k), ids);
            }
        }
        // d_C ⊗ id
        for (&o, outs) in &cx.out {
            let po = cx.objs[&o].p;
            for (&o2, f) in outs {
                let po2 = cx.objs[&o2].p;
                for k in 0..nt {
                    let phi = f.kron_id(loc.fsizes[k]);
                    let ss = self.decompose(po, positive, k);
                    let ts = self.decompose(po2, positive, k);
                    for (si, sg) in ss.iter().enumerate() {
                        let a = phi.matmul(&sg.iota);
                        if a.is_zero() {
                            continue;
                        }
                        for (ti, tu) in ts.iter().enumerate() {
                            let comp = tu.pi.matmul(&a);
                            out.add_map(place[&(o, k)][si], place[&(o2, k)][ti], comp);
                        }
                    }
                }
            }
        }
        // (−1)^h id ⊗ d
        for (&o, ob) in &cx.objs {
            for j in 0..nt.saturating_sub(1) {
                let (ks, kt) = if positive { (j, j + 1) } else { (j + 1, j) };
                let comps = self.dloc(ob.p, positive, j);
                let neg = ob.h.rem_euclid(2) == 1;
                for (si, ti, m) in comps.iter() {
                    out.add_map(place[&(o, ks)][*si], place[&(o, kt)][*ti], if neg { m.neg() } else { m.clone() });
                }
            }
        }
        out
    }

    fn check_degrees(&mut self, cx: &XCx<F>) {
        for (&a, outs) in &cx.out {
            let oa = cx.objs[&a];
            let la = self.ladder(oa.p);
            for (&b, m) in outs {
                let ob = cx.objs[&b];
                assert_eq!(ob.h, oa.h + 1, "component does not raise h by 1");
                let lb = self.ladder(ob.p);
                assert!(homogeneous(m, &la.b, &lb.b, oa.t - ob.t), "component of total q-degree != 0");
            }
        }
    }

    /// Gaussian elimination with scalar pivots in Markowitz order.
    ///
    /// Cost note: every pivot choice rescans all components, so elimination is quadratic in the
    /// size of the complex. That is irrelevant at M6.0 sizes (at most about 130 objects) but must
    /// not be carried over to sengine, which needs an incrementally updated pivot queue.
    fn simplify(&self, cx: &mut XCx<F>) {
        loop {
            let mut best: Option<(usize, usize, usize)> = None;
            for (&a, outs) in &cx.out {
                let oa = cx.objs[&a];
                for &b in outs.keys() {
                    let ob = cx.objs[&b];
                    if oa.p != ob.p || oa.t != ob.t {
                        continue;
                    }
                    let cost = (cx.inn[&b].len() - 1) * (outs.len() - 1);
                    if best.map_or(true, |x| (cost, a, b) < x) {
                        best = Some((cost, a, b));
                    }
                }
            }
            let Some((_, a, b)) = best else { break };
            let c = scalar_of(&cx.out[&a][&b]).expect("degree-0 component between equal types is not a scalar");
            assert!(!c.is_zero());
            let f = c.inv().neg();
            let xs: Vec<(usize, PMat<F>)> = cx.inn[&b].iter().filter(|&&x| x != a).map(|&x| (x, cx.out[&x][&b].clone())).collect();
            let ys: Vec<(usize, PMat<F>)> = cx.out[&a].iter().filter(|(y, _)| **y != b).map(|(&y, m)| (y, m.clone())).collect();
            for (x, g) in &xs {
                for (y, h) in &ys {
                    let m = scale_const(&h.matmul(g), &f);
                    cx.add_map(*x, *y, m);
                }
            }
            cx.remove_obj(a);
            cx.remove_obj(b);
        }
    }

    /// The nonzero components x -> z of d∘d.
    fn d2_components(&self, cx: &XCx<F>) -> Vec<(usize, usize, PMat<F>)> {
        let mut out = vec![];
        for (&x, outs) in &cx.out {
            let mut acc: BTreeMap<usize, PMat<F>> = BTreeMap::new();
            for (&y, f) in outs {
                for (&z, g) in &cx.out[&y] {
                    let m = g.matmul(f);
                    let e = acc.remove(&z);
                    acc.insert(z, match e {
                        None => m,
                        Some(w) => w.add(&m),
                    });
                }
            }
            out.extend(acc.into_iter().filter(|(_, m)| !m.is_zero()).map(|(z, m)| (x, z, m)));
        }
        out
    }

    fn decat_step(&mut self, c: usize, d: usize, positive: bool) {
        let terms = crossing_terms(c as u8, d as u8, positive, 0);
        let mut y: BTreeMap<i32, Elt> = BTreeMap::new();
        let x = std::mem::take(&mut self.decat);
        for (h, e) in &x {
            for t in &terms {
                let z = self.schur.mul_rungs(e, &t.rungs);
                let slot = y.entry(h + t.hdeg).or_insert_with(|| Elt::zero(&z.mu));
                *slot = slot.add_scaled(&z, 1, t.dq);
            }
        }
        y.retain(|_, e| !e.is_zero());
        self.decat = y;
    }

    fn decat_k0(&mut self) -> BTreeMap<(Tid, i32), i128> {
        let mut signed: BTreeMap<(Tid, i32), i128> = BTreeMap::new();
        let x = self.decat.clone();
        for (h, e) in &x {
            let sg: i128 = if h.rem_euclid(2) == 0 { 1 } else { -1 };
            for (p, m) in self.schur.expand(e) {
                for (i, &c) in m.c.iter().enumerate() {
                    if c != 0 {
                        *signed.entry((p, m.lo + i as i32)).or_insert(0) += sg * c;
                    }
                }
            }
        }
        signed.retain(|_, c| *c != 0);
        signed
    }

    /// Exact contract-v2 closure data of the final complex.
    /// (`d2`: nonzero components of d∘d of a pruned scan, returned as closure data F matrices.)
    #[allow(clippy::type_complexity)]
    fn closure_data(&mut self, cx: &XCx<F>, crossings: Vec<(bool, u8)>, seed: u64, d2: &[(usize, usize, PMat<F>)]) -> Result<(ClosureDataC<F>, Vec<(usize, usize, Vec<Vec<EPoly<F>>>)>), String> {
        let lam8: Vec<u8> = self.lam.iter().map(|&a| a as u8).collect();
        let k: usize = self.lam.iter().sum();
        let mut rng = StdRng::seed_from_u64(seed);
        let types: BTreeSet<Tid> = cx.objs.values().map(|o| o.p).collect();
        let mut bases: HashMap<Tid, LadderBasis<F>> = HashMap::new();
        let mut tbs: BTreeMap<u32, TypeBasisC<F>> = BTreeMap::new();
        for &p in &types {
            let lad = self.ladder(p);
            if lad.b.top_labels != self.lam {
                return Err(format!("final type {} is not closable", self.schur.table(p).fmt()));
            }
            let basis = rlambda_bases(&lad.b).map_err(|e| format!("type {}: {}", self.schur.table(p).fmt(), e))?;
            for _ in 0..2 {
                let x: Vec<F> = (0..k).map(|_| F::from_i64(rng.gen_range(1..2_000_000_000))).collect();
                let (a, b) = det_criterion(&lad.b, &basis, &x);
                if !(a && b) {
                    return Err(format!("type {}: one-point determinant criterion fails (ins {}, outs {})", self.schur.table(p).fmt(), a, b));
                }
            }
            let sh = lad.b.shift;
            let mut ins: Vec<i32> = basis.ins.iter().map(|g| g.0 + sh).collect();
            let mut outs: Vec<i32> = basis.outs.iter().map(|g| g.0 - sh).collect();
            let cd = self.schur.c_diag(p);
            let mut want: Vec<i32> = vec![];
            for (i, &c) in cd.c.iter().enumerate() {
                for _ in 0..c {
                    want.push(cd.lo + i as i32);
                }
            }
            let gram = gram_of(&basis, lad.b.ring.nvars);
            let tb = TypeBasisC { ins: ins.clone(), outs: outs.clone(), gram };
            ins.sort();
            outs.sort();
            want.sort();
            if ins != want || outs != want {
                return Err(format!("type {}: generator degrees ins {:?} outs {:?} != c_diag {:?}", self.schur.table(p).fmt(), ins, outs, want));
            }
            tbs.insert(p, tb);
            bases.insert(p, basis);
        }
        let ids: Vec<usize> = cx.objs.keys().copied().collect();
        let index: HashMap<usize, usize> = ids.iter().enumerate().map(|(i, &o)| (o, i)).collect();
        let objs: Vec<(u32, i32, i32)> = ids.iter().map(|o| { let x = cx.objs[o]; (x.p, x.h, x.t) }).collect();
        let mut edges = vec![];
        for (&a, outs) in &cx.out {
            let oa = cx.objs[&a];
            let la = self.ladder(oa.p);
            for (&b, f) in outs {
                let ob = cx.objs[&b];
                let lb = self.ladder(ob.p);
                let u = (oa.t - ob.t) - lb.b.shift + la.b.shift;
                let fm = solve_f(f, u, &bases[&oa.p], &bases[&ob.p], &lb.b).map_err(|e| format!("component {} -> {}: {}", a, b, e))?;
                edges.push((index[&a], index[&b], fm));
            }
        }
        let mut d2f = vec![];
        for (a, b, f) in d2 {
            let (oa, ob) = (cx.objs[a], cx.objs[b]);
            let (la, lb) = (self.ladder(oa.p), self.ladder(ob.p));
            let u = (oa.t - ob.t) - lb.b.shift + la.b.shift;
            let fm = solve_f(f, u, &bases[&oa.p], &bases[&ob.p], &lb.b).map_err(|e| format!("d∘d component {} -> {}: {}", a, b, e))?;
            d2f.push((index[a], index[b], fm));
        }
        let data = ClosureDataC { lam: lam8, objs, types: tbs, edges, crossings };
        data.validate()?;
        Ok((data, d2f))
    }
}

/// Inverse of a square matrix over F (None if singular).
fn invert<F: Field>(a: &[Vec<F>]) -> Option<Vec<Vec<F>>> {
    let n = a.len();
    let mut aug: Vec<Vec<F>> = a
        .iter()
        .enumerate()
        .map(|(i, r)| {
            let mut v = r.clone();
            v.extend((0..n).map(|j| if i == j { F::one() } else { F::zero() }));
            v
        })
        .collect();
    let order: Vec<usize> = (0..n).collect();
    let piv = crate::dense::rref_with_order(&mut aug, &order);
    if piv.len() != n {
        return None;
    }
    Some(aug.into_iter().map(|r| r[n..].to_vec()).collect())
}

fn check_input(labels: &[usize], word: &[Gen]) -> Result<(), String> {
    if labels.is_empty() || labels.len() > 2 {
        return Err(format!("the exact singular engine handles braids on at most 2 strands (got {})", labels.len()));
    }
    if let Some(&a) = labels.iter().find(|&&a| a == 0 || a > MAX_LABEL) {
        return Err(format!("colors must lie in 1..={} for the exact singular engine (got {})", MAX_LABEL, a));
    }
    for g in word {
        if g.i + 1 >= labels.len() {
            return Err(format!("generator {} needs {} strands", g.i + 1, g.i + 2));
        }
    }
    let mut lab = labels.to_vec();
    for g in word {
        lab.swap(g.i, g.i + 1);
    }
    if lab != labels {
        return Err(format!("colors do not close up: bottom {:?}, top {:?}", labels, lab));
    }
    Ok(())
}

/// The N-independent scan of a colored braid on at most 2 strands, with exact closure data.
pub fn xscan<F: Field>(labels: &[usize], word: &[Gen], opts: &XOpts) -> Result<XScan<F>, String> {
    check_input(labels, word)?;
    if let Some(n) = opts.prune_n {
        if labels.iter().any(|&a| a > n) {
            return Err(format!("prune_n = {} is below a color {:?}: the gl_N homology is 0 for every N <= {}", n, labels, n));
        }
        if !opts.check_d2 {
            // a pruned complex has d∘d = 0 only modulo I_N; its final components are certified
            // (closure_kills_map) only when check_d2 collects them, so refuse the uncertified path
            return Err(format!("prune_n = {} requires check_d2: d∘d vanishes only modulo the killed ideal, and the closure certifies the final components only when they are collected", n));
        }
    }
    let t0 = Instant::now();
    let mut eng = XEngine::<F>::new(labels, opts.verbose, opts.prune_n);
    eng.rungwise = opts.rungwise;
    eng.check_orth = opts.check_d2;
    let p0 = eng.identity_type();
    eng.ladder(p0);
    let mut cx = XCx::<F>::new();
    cx.add_obj(XObj { p: p0, h: 0, t: 0 });
    let mut labs = labels.to_vec();
    let mut stages = vec![];
    let mut crossings = vec![];
    let mut last_d2: Vec<(usize, usize, PMat<F>)> = vec![];
    for (step, g) in word.iter().enumerate() {
        let ts = Instant::now();
        let (c, d) = (labs[0], labs[1]);
        crossings.push((g.positive, c.min(d) as u8));
        eng.decat_step(c, d, g.positive);
        eng.stage_dropped = 0;
        let mut ncx = eng.tensor_crossing(&cx, &labs, g.positive);
        let unsimplified = ncx.objs.len();
        eng.check_degrees(&ncx);
        eng.simplify(&mut ncx);
        let mut d2_ideal = 0;
        if opts.check_d2 {
            let d2 = eng.d2_components(&ncx);
            if !d2.is_empty() {
                if eng.prune.is_none() {
                    return Err(format!("d∘d != 0 from object {} to object {}", d2[0].0, d2[0].1));
                }
                // pruned: d∘d = −Σ_{killed y} d_{y·} d_{·y} lies in the ideal of killed types;
                // its components on the final complex are checked to close to 0 (`XScan::close`)
                d2_ideal = d2.len();
            }
            last_d2 = d2;
        }
        labs.swap(0, 1);
        // pruned scan: the class of the pruned complex is the decategorified class with the
        // killed types removed (the killed span is an ideal; dropping it is the projection)
        let mut want = eng.decat_k0();
        want.retain(|(p, _), _| !eng.killed(*p));
        let got = ncx.k0();
        if want != got {
            let fmt = |m: &BTreeMap<(Tid, i32), i128>, s: &Schur| m.iter().map(|((p, t), c)| format!("{}{{{}}}:{}", s.table(*p).fmt(), t, c)).collect::<Vec<_>>().join(" ");
            return Err(format!("crossing {}: K_0 class {} != decategorified scan {}", step + 1, fmt(&got, &eng.schur), fmt(&want, &eng.schur)));
        }
        let lb: u128 = want.values().map(|c| c.unsigned_abs()).sum();
        if (ncx.objs.len() as u128) < lb {
            return Err(format!("crossing {}: {} objects < lower bound {}", step + 1, ncx.objs.len(), lb));
        }
        let info = StageInfo {
            step: step + 1,
            labels: labs.clone(),
            unsimplified,
            objects: ncx.objs.len(),
            components: ncx.components(),
            lower_bound: lb,
            secs: ts.elapsed().as_secs_f64(),
            dropped: eng.stage_dropped,
            d2_ideal,
        };
        if opts.verbose {
            eprintln!(
                "[xengine] crossing {}/{} ({}{}): {} objects -> {} (LB {}), {} components, {:.2}s{}",
                step + 1,
                word.len(),
                if g.positive { "+" } else { "-" },
                format!("({},{})", c, d),
                unsimplified,
                info.objects,
                lb,
                info.components,
                info.secs,
                if eng.prune.is_some() { format!(" [pruned at N={}: {} summands dropped, d∘d components in the ideal {}]", eng.prune.unwrap(), info.dropped, d2_ideal) } else { String::new() }
            );
        }
        stages.push(info);
        cx = ncx;
    }
    let t_scan = t0.elapsed().as_secs_f64();
    let t1 = Instant::now();
    if word.is_empty() || !opts.check_d2 {
        last_d2.clear();
    }
    let (data, final_d2) = eng.closure_data(&cx, crossings, opts.seed, &last_d2)?;
    let t_data = t1.elapsed().as_secs_f64();
    if opts.verbose {
        eprintln!("[xengine] closure data: {} objects, {} types, {} edges in {:.2}s", data.objs.len(), data.types.len(), data.edges.len(), t_data);
    }
    let mut final_objects: Vec<(String, i32, i32)> = cx.objs.values().map(|o| (eng.schur.table(o.p).fmt(), o.h, o.t)).collect();
    final_objects.sort();
    let k0 = cx.k0().into_iter().map(|((p, t), c)| (eng.schur.table(p).fmt(), t, c)).collect();
    let mut ladders: Vec<(String, String, usize)> = eng
        .lads
        .iter()
        .map(|(p, l)| {
            let w: Vec<String> = l.rungs.iter().map(|r| format!("{}{}({})", if r.dir == Dir::E { "E" } else { "F" }, r.i, r.r)).collect();
            (eng.schur.table(*p).fmt(), w.join(" "), l.b.fsize())
        })
        .collect();
    ladders.sort();
    let type_ladders: BTreeMap<u32, Vec<Rung>> = data.types.keys().map(|&p| (p, eng.lads[&p].rungs.clone())).collect();
    let scan_ladders: BTreeMap<String, Vec<Rung>> = eng.lads.iter().map(|(&p, l)| (eng.schur.table(p).fmt(), l.rungs.clone())).collect();
    Ok(XScan {
        labels: labels.to_vec(),
        word: word.to_vec(),
        data,
        stages,
        final_objects,
        k0,
        hom_calls: eng.hom_calls,
        hom_total_secs: eng.hom_total,
        hom_max_secs: eng.hom_max,
        max_x_fsize: eng.max_x,
        decompositions: eng.direct_decomps,
        identity_decompositions: eng.identity_decomps,
        ladders,
        type_ladders,
        scan_ladders,
        timings: (t_scan, t_data),
        prune_n: opts.prune_n,
        dropped_summands: eng.dropped_total,
        all_killed_decompositions: eng.all_killed_decomps,
        final_d2,
    })
}

impl<F: Field> XScan<F> {
    /// gl_N Poincaré polynomials for every N in `ns` (reduced: block 0 marked). Colors > N give 0.
    /// A pruned scan (`prune_n = Some(N0)`) is valid for N <= N0 only (panics otherwise); there
    /// the closure of every final d∘d component is first checked to vanish (panics otherwise).
    pub fn close(&self, ns: &[usize], reduced: bool) -> BTreeMap<usize, Poincare> {
        let marked = if reduced { Some(0) } else { None };
        ns.iter()
            .map(|&nn| {
                if let Some(n0) = self.prune_n {
                    assert!(nn <= n0, "a scan pruned at N = {} gives the gl_N homology only for N <= {} (requested N = {})", n0, n0, nn);
                    for (a, b, f) in &self.final_d2 {
                        assert!(closure_kills_map(&self.data, *a, *b, f, nn, marked), "pruned scan at N = {}: the d∘d component {} -> {} does not close to 0 at gl_{}", n0, a, b, nn);
                    }
                }
                (nn, close_exact_colored(&self.data, nn, marked))
            })
            .collect()
    }
}

/// Does the closure at gl_N (`marked`: reduced) send the map f: B_x → B_w between two final
/// objects of `data` (closure-data F matrix, rows = ins of w, columns = ins of x; normalized
/// degree t_x − t_w) to 0? The two objects are placed in homological degrees 0 and 1 with the
/// single edge f; the closed map is 0 iff the Poincaré polynomial equals that without the edge.
pub fn closure_kills_map<F: Field>(data: &ClosureDataC<F>, x: usize, w: usize, f: &[Vec<EPoly<F>>], nn: usize, marked: Option<usize>) -> bool {
    let (ox, ow) = (data.objs[x], data.objs[w]);
    let mut types = BTreeMap::new();
    types.insert(ox.0, data.types[&ox.0].clone());
    types.insert(ow.0, data.types[&ow.0].clone());
    let two = |edges| ClosureDataC { lam: data.lam.clone(), objs: vec![(ox.0, 0, ox.2), (ow.0, 1, ow.2)], types: types.clone(), edges, crossings: data.crossings.clone() };
    let with = two(vec![(0, 1, f.to_vec())]);
    with.validate().unwrap_or_else(|e| panic!("closure_kills_map: {}", e));
    close_exact_colored(&with, nn, marked) == close_exact_colored(&two(vec![]), nn, marked)
}

/// Colored gl_N homology of the closure of a braid on at most 2 strands (labels per bottom
/// position in 1..=3), for every N in `ns`, by the exact singular engine (one scan for all N).
/// `reduced`: the block at bottom position 0 is marked (reduced unknot = q^0).
pub fn khr_singular_exact<F: Field>(labels: &[usize], word: &[Gen], ns: &[usize], reduced: bool) -> Result<BTreeMap<usize, Poincare>, String> {
    let scan = xscan::<F>(labels, word, &XOpts::default())?;
    Ok(scan.close(ns, reduced))
}
