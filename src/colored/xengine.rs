//! M6.0: the exact small colored engine (SYNTHESIS §9.2; design note
//! `docs/dev/colored_design/M6_0_design.md`).
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
//!   Then Gaussian elimination with scalar pivots in Markowitz order (a lazy heap reproducing
//!   the full-rescan order).
//! * **Memory (M2S′, SPEC_3C_GL6_rev1 §2; measurements in gl6/M2S_step0.md).** By default
//!   (`XOpts::stream`) C ⊗ X is built slice by slice in the homological degree of C and every
//!   eligible pivot is eliminated after each slice (`tensor_crossing_stream`: a → b only once
//!   out(a) and in(b) are complete); `stream = false` builds the whole tensor complex first (the
//!   cross-check). Components are stored compactly (`CMat`, `XOpts::compact`: interned u32
//!   monomials and u32 coefficients, F31/F32b) and decoded on use. d_C ⊗ id is transported rung
//!   by rung through the composite decompositions (`XEngine::transport`,
//!   `XOpts::rung_transport`). `XOpts::close_stages` keeps the closure data of every stage
//!   (`XScan::close_stage`).
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

use crate::colored::cclose::{close_exact_colored, det_criterion, try_close_exact_colored, gram_of, rlambda_bases, solve_f, ClosureDataC, EPoly, LadderBasis, TypeBasisC};
use crate::input::crossing::Gen;
use crate::reference::cube::{apply_local_poly, local_cx, LocalCx};
use crate::algebra::field::Field;
use crate::algebra::laurent::Poincare;
use crate::algebra::mpoly::{MonoCache, Poly};
use crate::algebra::pmat::{PAcc, PMat, PVec};
use crate::soergel::sbim::{hom_space, SBim};
use crate::combinat::schur::{crossing_terms, flows, Elt, Lp, Schur, Tid};
use crate::combinat::web::{Dir, Rung};
use rand::{rngs::StdRng, Rng, SeedableRng};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::rc::Rc;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};
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
    /// build C ⊗ X slice by slice in h and eliminate eligible pivots after every slice
    /// (`tensor_crossing_stream`); false: the full build, then one elimination (cross-check)
    pub stream: bool,
    /// keep the closure data of every stage m >= 2 (`XScan::close_stage`)
    pub close_stages: bool,
    /// store the components of the complexes compactly (`CMat`; F31 and F32b only)
    pub compact: bool,
    /// streaming build: transport d_C ⊗ id rung by rung (`transport_view`, M2S′ step 4)
    pub rung_transport: bool,
    /// with `close_stages`: close every kept stage at these N right away (unreduced and
    /// reduced) and print the Poincaré polynomials on stderr, so that a long scan reports its
    /// stage results before it ends
    pub stage_close_ns: Vec<usize>,
}

impl Default for XOpts {
    fn default() -> Self {
        let verbose = std::env::var("KRHOM_VERBOSE").map(|v| v != "0" && !v.is_empty()).unwrap_or(false);
        XOpts { check_d2: true, verbose, seed: 7, rungwise: true, prune_n: None, stream: true, close_stages: false, compact: true, rung_transport: true, stage_close_ns: vec![] }
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
    /// measurements (M2S′ step 0): components and Poly terms of C ⊗ X before `simplify`
    /// (full build; streaming: the largest unsimplified part held at once) and of the minimal
    /// complex after it
    pub mem: StageMem,
}

/// Memory and time measurements of one stage (M2S′ step 0).
#[derive(Clone, Debug, Default)]
pub struct StageMem {
    /// components / nonzero entries / Poly terms of the complex before simplify
    pub comps_before: usize,
    pub entries_before: usize,
    pub terms_before: usize,
    /// estimated bytes of the stored components before simplify
    pub bytes_before: usize,
    /// the same after simplify (the minimal complex)
    pub entries_after: usize,
    pub terms_after: usize,
    pub bytes_after: usize,
    /// estimated bytes of the caches (decompositions incl. composite ones, dloc) after the stage
    pub bytes_caches: usize,
    /// RSS (kB) after the tensor product (before simplify), after the stage, and the peak RSS
    /// of the stage (VmHWM, reset at the start of the stage where the kernel allows it)
    pub rss_tensor_kb: usize,
    pub rss_end_kb: usize,
    pub rss_peak_kb: usize,
    /// seconds: decompositions (Hom spaces, orthogonalization, composites), transport (the
    /// rest of the tensor product: kron_id, products with ι and π, dloc), elimination, d∘d
    pub secs_decomp: f64,
    pub secs_transport: f64,
    pub secs_elim: f64,
    pub secs_d2: f64,
    /// streaming: number of partial simplify calls and the largest object count held at once
    pub stream_passes: usize,
    pub peak_objects: usize,
}

/// (VmRSS, VmHWM) of this process in kB (0 where /proc is unavailable).
pub fn rss_kb() -> (usize, usize) {
    let s = std::fs::read_to_string("/proc/self/status").unwrap_or_default();
    let get = |k: &str| s.lines().find(|l| l.starts_with(k)).and_then(|l| l.split_whitespace().nth(1)).and_then(|x| x.parse().ok()).unwrap_or(0);
    (get("VmRSS:"), get("VmHWM:"))
}

/// Reset the peak RSS (VmHWM) of this process; false if the kernel refuses.
pub fn reset_hwm() -> bool {
    std::fs::write("/proc/self/clear_refs", "5").is_ok()
}

/// (nonzero entries, Poly terms, estimated heap bytes) of a PMat: a `PVec` entry (u32, Poly) is
/// 32 B, a term (u128 Mono, F) 32 B, a column header 24 B.
pub fn pmat_size<F: Field>(m: &PMat<F>) -> (usize, usize, usize) {
    let mut e = 0;
    let mut t = 0;
    let mut b = 24 * m.cols.len() + 48;
    for c in &m.cols {
        e += c.len();
        b += 32 * c.capacity();
        for (_, p) in c {
            t += p.t.len();
            b += 32 * p.t.capacity();
        }
    }
    (e, t, b)
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
    /// `XOpts::close_stages`: closure data (and final d∘d components, as `final_d2`) of the
    /// complex after crossing m, for every m >= 2 whose top labels close up, read by
    /// `close_stage`; the crossings of the data are `crossings[..m]`
    #[allow(clippy::type_complexity)]
    pub stage_data: BTreeMap<usize, (ClosureDataC<F>, Vec<(usize, usize, Vec<Vec<EPoly<F>>>)>)>,
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

// ------------------------------------------------------------------ compact component store

/// Interned monomials of the compact store (index -> Mono, Mono -> index), shared by all
/// threads; indices are never reused, so a packed component stays valid.
type MonoTab = (Vec<crate::algebra::mpoly::Mono>, HashMap<crate::algebra::mpoly::Mono, u32>);
fn monos() -> &'static std::sync::RwLock<MonoTab> {
    static M: std::sync::OnceLock<std::sync::RwLock<MonoTab>> = std::sync::OnceLock::new();
    M.get_or_init(|| std::sync::RwLock::new((vec![], HashMap::new())))
}

/// Representative in [0, p) of a coefficient of F31 or F32b (None for any other field).
fn coef_u32<F: Field>(c: &F) -> Option<u32> {
    let a = c as &dyn std::any::Any;
    if let Some(x) = a.downcast_ref::<crate::algebra::field::F31>() {
        return Some(x.0 as u32);
    }
    a.downcast_ref::<crate::algebra::field::F32b>().map(|x| x.0 as u32)
}

fn packable<F: Field>() -> bool {
    coef_u32(&F::one()).is_some()
}

/// A component in the compact store (M2S′ step 2): column-compressed, every entry a row index
/// and a term range, every term (interned monomial index << 32 | coefficient representative),
/// i.e. 8 B per entry and 8 B per term instead of 32 B per entry and 32 B per term (u128 Mono +
/// Fp(u64)) plus vector headers. Only for F31 and F32b (both < 2^32); other fields keep the PMat.
#[derive(Clone, Debug)]
pub enum CMat<F: Field> {
    Plain(PMat<F>),
    Packed { nrows: u32, ncols: u32, colptr: Vec<u32>, rows: Vec<u32>, tptr: Vec<u32>, terms: Vec<u64> },
}

impl<F: Field> CMat<F> {
    /// Encode (compact: pack if the field allows it).
    pub fn new(m: PMat<F>, compact: bool) -> Self {
        if !compact || !packable::<F>() {
            return CMat::Plain(m);
        }
        let ne: usize = m.cols.iter().map(|c| c.len()).sum();
        let mut colptr = Vec::with_capacity(m.ncols + 1);
        let mut rows = Vec::with_capacity(ne);
        let mut tptr = Vec::with_capacity(ne + 1);
        let mut terms = Vec::new();
        colptr.push(0u32);
        tptr.push(0u32);
        // monomial indices: read lock, then a write lock only for new monomials
        let mut idx: Vec<u32> = Vec::with_capacity(m.cols.iter().flat_map(|c| c.iter()).map(|(_, p)| p.t.len()).sum());
        let mut missing = false;
        {
            let mt = monos().read().unwrap();
            for c in &m.cols {
                for (_, p) in c {
                    for (mono, _) in &p.t {
                        match mt.1.get(mono) {
                            Some(&i) => idx.push(i),
                            None => {
                                missing = true;
                                idx.push(u32::MAX);
                            }
                        }
                    }
                }
            }
        }
        if missing {
            let mut mt = monos().write().unwrap();
            let mut it = 0;
            for c in &m.cols {
                for (_, p) in c {
                    for (mono, _) in &p.t {
                        if idx[it] == u32::MAX {
                            idx[it] = match mt.1.get(mono) {
                                Some(&i) => i,
                                None => {
                                    let i = mt.0.len() as u32;
                                    mt.0.push(*mono);
                                    mt.1.insert(*mono, i);
                                    i
                                }
                            };
                        }
                        it += 1;
                    }
                }
            }
        }
        let mut it = 0;
        for c in &m.cols {
            for (r, p) in c {
                rows.push(*r);
                for (_, cf) in &p.t {
                    terms.push(((idx[it] as u64) << 32) | coef_u32(cf).unwrap() as u64);
                    it += 1;
                }
                tptr.push(terms.len() as u32);
            }
            colptr.push(rows.len() as u32);
        }
        assert!(terms.len() < u32::MAX as usize && rows.len() < u32::MAX as usize, "compact store: component too large");
        terms.shrink_to_fit();
        CMat::Packed { nrows: m.nrows as u32, ncols: m.ncols as u32, colptr, rows, tptr, terms }
    }
    /// Decode to a PMat.
    pub fn get(&self) -> PMat<F> {
        match self {
            CMat::Plain(m) => m.clone(),
            CMat::Packed { nrows, ncols, colptr, rows, tptr, terms } => {
                let mt = monos().read().unwrap();
                let mut cols = Vec::with_capacity(*ncols as usize);
                for j in 0..*ncols as usize {
                    let (e0, e1) = (colptr[j] as usize, colptr[j + 1] as usize);
                    let mut col: PVec<F> = Vec::with_capacity(e1 - e0);
                    for e in e0..e1 {
                        let t: Vec<(crate::algebra::mpoly::Mono, F)> = terms[tptr[e] as usize..tptr[e + 1] as usize].iter().map(|&w| (mt.0[(w >> 32) as usize], F::from_i64((w & 0xffff_ffff) as i64))).collect();
                        col.push((rows[e], Poly { t }));
                    }
                    cols.push(col);
                }
                PMat { nrows: *nrows as usize, ncols: *ncols as usize, cols }
            }
        }
    }
    /// (entries, terms, estimated heap bytes)
    pub fn size(&self) -> (usize, usize, usize) {
        match self {
            CMat::Plain(m) => pmat_size(m),
            CMat::Packed { colptr, rows, tptr, terms, .. } => (rows.len(), terms.len(), 4 * (colptr.capacity() + rows.capacity() + tptr.capacity()) + 8 * terms.capacity() + 120),
        }
    }
}

/// A complex of canonical ladders with exact components (stored as `CMat`).
struct XCx<F: Field> {
    objs: BTreeMap<usize, XObj>,
    out: BTreeMap<usize, BTreeMap<usize, CMat<F>>>,
    inn: BTreeMap<usize, BTreeSet<usize>>,
    next: usize,
    /// pack the components (`XOpts::compact`)
    compact: bool,
}

impl<F: Field> XCx<F> {
    /// (components, entries, terms, bytes) of the stored components (+ 64 B per map node)
    fn size(&self) -> (usize, usize, usize, usize) {
        let (mut c, mut e, mut t, mut b) = (0, 0, 0, 0);
        for outs in self.out.values() {
            for m in outs.values() {
                let (e1, t1, b1) = m.size();
                c += 1;
                e += e1;
                t += t1;
                b += b1 + 64;
            }
        }
        (c, e, t, b)
    }
    fn new(compact: bool) -> Self {
        XCx { objs: BTreeMap::new(), out: BTreeMap::new(), inn: BTreeMap::new(), next: 0, compact }
    }
    /// the component x → y, decoded
    fn map(&self, x: usize, y: usize) -> PMat<F> {
        self.out[&x][&y].get()
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
        let compact = self.compact;
        let outs = self.out.get_mut(&x).unwrap();
        let sum = match outs.remove(&y) {
            Some(old) => old.get().add(&m),
            None => m,
        };
        if sum.is_zero() {
            self.inn.get_mut(&y).unwrap().remove(&x);
        } else {
            outs.insert(y, CMat::new(sum, compact));
            self.inn.get_mut(&y).unwrap().insert(x);
        }
    }
    /// replace the component x → y (None or a zero matrix removes it)
    fn set_map(&mut self, x: usize, y: usize, m: Option<CMat<F>>) {
        match m {
            Some(m) => {
                self.out.get_mut(&x).unwrap().insert(y, m);
                self.inn.get_mut(&y).unwrap().insert(x);
            }
            None => {
                self.out.get_mut(&x).unwrap().remove(&y);
                self.inn.get_mut(&y).unwrap().remove(&x);
            }
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

/// a · b for PMats, by sorting the term products of each column instead of hashing them
/// (same result as `PMat::matmul`, in the same canonical form: rows ascending, monomials
/// ascending, no zero entries or terms); `c`: also scale by the constant c.
fn mm<F: Field>(a: &PMat<F>, b: &PMat<F>, c: Option<&F>) -> PMat<F> {
    assert_eq!(a.ncols, b.nrows);
    let mut buf: Vec<(u32, crate::algebra::mpoly::Mono, F)> = vec![];
    let cols = b
        .cols
        .iter()
        .map(|col| {
            buf.clear();
            for (k, pk) in col {
                for (i, pa) in &a.cols[*k as usize] {
                    for (ma, ca) in &pa.t {
                        for (mb, cb) in &pk.t {
                            buf.push((*i, ma + mb, ca.mul(cb)));
                        }
                    }
                }
            }
            buf.sort_unstable_by_key(|x| (x.0, x.1));
            let mut out: PVec<F> = vec![];
            let mut cur: Option<(u32, Vec<(crate::algebra::mpoly::Mono, F)>)> = None;
            let mut j = 0;
            while j < buf.len() {
                let (r, m) = (buf[j].0, buf[j].1);
                let mut v = buf[j].2.clone();
                j += 1;
                while j < buf.len() && buf[j].0 == r && buf[j].1 == m {
                    v.add_assign(&buf[j].2);
                    j += 1;
                }
                if let Some(cc) = c {
                    v.mul_assign(cc);
                }
                if v.is_zero() {
                    continue;
                }
                match &mut cur {
                    Some((cr, t)) if *cr == r => t.push((m, v)),
                    _ => {
                        if let Some((cr, t)) = cur.take() {
                            out.push((cr, Poly { t }));
                        }
                        cur = Some((r, vec![(m, v)]));
                    }
                }
            }
            if let Some((cr, t)) = cur {
                out.push((cr, Poly { t }));
            }
            out
        })
        .collect();
    PMat { nrows: a.nrows, ncols: b.ncols, cols }
}

/// Key of the composite decomposition cache: (type, rung word).
type WKey = (Tid, Vec<(u8, usize, usize)>);

/// Map `f` over `jobs` on `crate::num_threads()` scoped threads (in order; sequential for
/// fewer than two jobs or one thread).
fn par_map<T: Sync, R: Send>(jobs: &[T], f: impl Fn(&T) -> R + Sync) -> Vec<R> {
    let nt = crate::num_threads().min(jobs.len());
    if nt <= 1 {
        return jobs.iter().map(&f).collect();
    }
    let next = AtomicUsize::new(0);
    let slots: Vec<Mutex<Option<R>>> = (0..jobs.len()).map(|_| Mutex::new(None)).collect();
    std::thread::scope(|sc| {
        for _ in 0..nt {
            sc.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::Relaxed);
                if i >= jobs.len() {
                    break;
                }
                let r = f(&jobs[i]);
                *slots[i].lock().unwrap() = Some(r);
            });
        }
    });
    slots.into_iter().map(|m| m.into_inner().unwrap().unwrap()).collect()
}

/// Transported components of f ⊗ id_W between the summands of L_p ⊗ W and L_p2 ⊗ W (W = `rs`)
/// from the composite decomposition cache `wd`, which must hold every decomposition used
/// (`XEngine::decompose_word(p, rs)` and `(p2, rs)` fill it): (source index, target index,
/// π_t ∘ (f ⊗ id) ∘ ι_s), nonzero only. Rung by rung (M2S′ step 4): with ι = (ι_1 ⊗ id) ∘ ι_2
/// and π = π_2 ∘ (π_1 ⊗ id), π_t (f ⊗ id_{R_1 R_2}) ι_s = π_2 ((π_1 (f ⊗ id_{R_1}) ι_1) ⊗
/// id_{R_2}) ι_2, so f is first transported through R_1 (between the summands B_r, B_r' of
/// L_p ⊗ R_1, L_p2 ⊗ R_1) and each nonzero result then through the rest; the large f ⊗ id_W is
/// never formed. `rungwise = false`: in one step through the summands of L_p ⊗ W.
fn transport_view<F: Field>(wd: &HashMap<WKey, Arc<Vec<Summand<F>>>>, rungwise: bool, p: Tid, p2: Tid, f: &PMat<F>, rs: &[Rung]) -> Vec<(usize, usize, PMat<F>)> {
    let get = |q: Tid, w: &[Rung]| wd.get(&(q, w.iter().map(rung_key).collect::<Vec<_>>())).unwrap_or_else(|| panic!("transport: decomposition of type {} ⊗ {:?} not cached", q, w)).clone();
    let ss = get(p, rs);
    let ts = get(p2, rs);
    let mut out = vec![];
    if ss.is_empty() || ts.is_empty() {
        return out;
    }
    if rs.len() <= 1 || !rungwise {
        let fw = ss[0].iota.nrows / f.ncols;
        assert_eq!(fw * f.ncols, ss[0].iota.nrows);
        let phi = if fw == 1 { f.clone() } else { f.kron_id(fw) };
        for (si, sg) in ss.iter().enumerate() {
            let a = mm(&phi, &sg.iota, None);
            if a.is_zero() {
                continue;
            }
            for (ti, tu) in ts.iter().enumerate() {
                let comp = mm(&tu.pi, &a, None);
                if !comp.is_zero() {
                    out.push((si, ti, comp));
                }
            }
        }
        return out;
    }
    let first_s = get(p, &rs[..1]);
    let first_t = get(p2, &rs[..1]);
    let offs = |first: &[Summand<F>]| -> Vec<usize> {
        let mut o = Vec::with_capacity(first.len());
        let mut acc = 0;
        for sm in first {
            o.push(acc);
            acc += get(sm.q, &rs[1..]).len();
        }
        o
    };
    let os = offs(&first_s);
    let ot = offs(&first_t);
    for (i, j, g) in transport_view(wd, rungwise, p, p2, f, &rs[..1]) {
        for (u, v, h) in transport_view(wd, rungwise, first_s[i].q, first_t[j].q, &g, &rs[1..]) {
            out.push((os[i] + u, ot[j] + v, h));
        }
    }
    out
}

// ------------------------------------------------------------------ the engine

struct XEngine<F: Field> {
    lam: Vec<usize>,
    schur: Schur,
    base: SBim<F>,
    mc: MonoCache,
    lads: HashMap<Tid, Rc<Lad<F>>>,
    locals: HashMap<(usize, usize, bool), Rc<LocalCx<F>>>,
    decomp: HashMap<(Tid, bool, usize), Arc<Vec<Summand<F>>>>,
    wdecomp: HashMap<WKey, Arc<Vec<Summand<F>>>>,
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
    /// seconds spent in decompositions (cache misses of `decompose`)
    t_decomp: f64,
    /// objects inserted by streaming tensor products (all stages)
    stream_inserted: usize,
    /// compact component store (`XOpts::compact`)
    compact: bool,
    /// transport d_C ⊗ id rung by rung in the streaming build (`XOpts::rung_transport`)
    transport_rungwise: bool,
}

impl<F: Field> XEngine<F> {
    /// estimated bytes of the decomposition and dloc caches
    fn cache_bytes(&self) -> usize {
        let mut b = 0;
        for v in self.wdecomp.values() {
            for sm in v.iter() {
                b += pmat_size(&sm.iota).2 + pmat_size(&sm.pi).2;
            }
        }
        for v in self.dloc.values() {
            for (_, _, m) in v.iter() {
                b += pmat_size(m).2;
            }
        }
        b
    }

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
            t_decomp: 0.0,
            stream_inserted: 0,
            compact: true,
            transport_rungwise: true,
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
    fn decompose(&mut self, p: Tid, positive: bool, k: usize) -> Arc<Vec<Summand<F>>> {
        if let Some(d) = self.decomp.get(&(p, positive, k)) {
            return d.clone();
        }
        let t0 = Instant::now();
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
        self.t_decomp += t0.elapsed().as_secs_f64();
        r
    }

    /// Decomposition of L_p ⊗ (rungs), cached. Rung by rung (default): L_p ⊗ R_1 = ⊕ L_r{s_1},
    /// then L_r ⊗ (R_2 ...) recursively, with ι = (ι_1 ⊗ id) ∘ ι_2 and π = π_2 ∘ (π_1 ⊗ id); in the
    /// direct mode the whole word is decomposed at once.
    fn decompose_word(&mut self, p: Tid, rs: &[Rung]) -> Arc<Vec<Summand<F>>> {
        let key = (p, rs.iter().map(rung_key).collect::<Vec<_>>());
        if let Some(d) = self.wdecomp.get(&key) {
            return d.clone();
        }
        let r = if rs.is_empty() {
            debug_assert!(!self.killed(p));
            let n = self.ladder(p).b.fsize();
            Arc::new(vec![Summand { q: p, s: 0, iota: PMat::identity(n), pi: PMat::identity(n) }])
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
            Arc::new(out)
        };
        self.wdecomp.insert(key, r.clone());
        r
    }

    /// Decomposition of X = L_p ⊗ (rungs) by Hom spaces and intersection forms.
    fn decompose_direct(&mut self, p: Tid, rs: &[Rung]) -> Arc<Vec<Summand<F>>> {
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
            return Arc::new(vec![]);
        }
        // X is literally the canonical ladder of its only summand B_q (multiplicity 1, shift 0):
        // ι = π = id, no Hom space needed
        if mult.len() == 1 && mult[0].1 == crate::combinat::schur::Lp::one() {
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
                return Arc::new(vec![Summand { q, s: 0, iota: PMat::identity(n), pi: PMat::identity(n) }]);
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
                let icols = crate::algebra::dense::rref(&mut gr, iotas.len());
                assert_eq!(icols.len(), cnt, "decomposition of {} ⊗ [{}]: rank of the intersection form for {} at v^{} is {}, multiplicity {}", self.schur.table(p).fmt(), word(), self.schur.table(q).fmt(), s, icols.len(), cnt);
                let mut gt: Vec<Vec<F>> = (0..iotas.len()).map(|i| (0..pis.len()).map(|j| g[j][i].clone()).collect()).collect();
                let jrows = crate::algebra::dense::rref(&mut gt, pis.len());
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
        Arc::new(sums)
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
        let mut out = XCx::new(self.compact);
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
                    let phi = f.get().kron_id(loc.fsizes[k]);
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

    /// Streaming C ⊗ (Rickard complex), simplified (`XOpts::stream`). The objects o of C are
    /// inserted slice by slice in h(o): slice s brings o ⊗ W_k for every o with h(o) = s, the
    /// components id ⊗ d inside each o ⊗ X, and d_C ⊗ id from slice s − 1. After slice s, an
    /// object of slice ≤ s has all its incoming components (they come from slices s − 1 and s)
    /// and an object of slice ≤ s − 1 all its outgoing ones; a pivot a → b is eliminated only if
    /// out(a) and in(b) are complete (SPEC_M4_M6a §6.4). Gaussian elimination of a → b uses only
    /// d_{x b} and d_{a y}; components into a or out of b inserted later are dropped, which is
    /// what the elimination of a → b does to them in the full build. C is consumed slice by
    /// slice. Returns the complex (simplified after the last slice), and the largest
    /// (objects, components, entries, terms, bytes) held at once before a partial pass.
    #[allow(clippy::type_complexity)]
    fn tensor_crossing_stream(&mut self, mut cx: XCx<F>, labels: &[usize], positive: bool, mem: &mut StageMem) -> XCx<F> {
        let (c, d) = (labels[0], labels[1]);
        let loc = self.local(c, d, positive);
        let nt = loc.cx.terms.len();
        let mut out = XCx::new(self.compact);
        let mut place: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
        let mut slice_of: HashMap<usize, i32> = HashMap::new();
        let mut by_h: BTreeMap<i32, Vec<usize>> = BTreeMap::new();
        for (&o, ob) in &cx.objs {
            by_h.entry(ob.h).or_default().push(o);
        }
        let hs: Vec<i32> = by_h.keys().copied().collect();
        for (si, &hs_cur) in hs.iter().enumerate() {
            let tt = Instant::now();
            let td0 = self.t_decomp;
            for &o in &by_h[&hs_cur] {
                let ob = cx.objs[&o];
                for k in 0..nt {
                    let sums = self.decompose(ob.p, positive, k);
                    self.stage_dropped += self.dropped[&(ob.p, positive, k)];
                    let hk = loc.cx.terms[k].hdeg;
                    let dq = if positive { -(k as i32) } else { k as i32 };
                    let ids: Vec<usize> = sums.iter().map(|sm| out.add_obj(XObj { p: sm.q, h: ob.h + hk, t: ob.t + dq + sm.s })).collect();
                    self.stream_inserted += ids.len();
                    for &id in &ids {
                        slice_of.insert(id, hs_cur);
                    }
                    place.insert((o, k), ids);
                }
                // (−1)^h id ⊗ d
                for j in 0..nt.saturating_sub(1) {
                    let (ks, kt) = if positive { (j, j + 1) } else { (j + 1, j) };
                    let comps = self.dloc(ob.p, positive, j);
                    let neg = ob.h.rem_euclid(2) == 1;
                    for (sx, tx, m) in comps.iter() {
                        out.add_map(place[&(o, ks)][*sx], place[&(o, kt)][*tx], if neg { m.neg() } else { m.clone() });
                    }
                }
            }
            // d_C ⊗ id from slice s − 1 into slice s (sources eliminated as b are skipped),
            // transported in parallel (`transport_view`) in chunks, added in job order
            if si > 0 && hs[si - 1] == hs_cur - 1 {
                let mut jobs: Vec<(usize, usize, usize, Arc<PMat<F>>)> = vec![];
                for &o in &by_h[&hs[si - 1]] {
                    let po = cx.objs[&o].p;
                    let outs = cx.out.remove(&o).unwrap_or_default();
                    for (o2, f) in outs {
                        let po2 = cx.objs[&o2].p;
                        let f = Arc::new(f.get());
                        for k in 0..nt {
                            self.decompose(po, positive, k);
                            self.decompose(po2, positive, k);
                            if place[&(o, k)].iter().any(|x| out.objs.contains_key(x)) {
                                jobs.push((o, o2, k, f.clone()));
                            }
                        }
                    }
                }
                let rw = self.transport_rungwise && self.rungwise;
                let words: Vec<Vec<Rung>> = loc.cx.terms.iter().map(|t| t.rungs.clone()).collect();
                let chunk = 8 * crate::num_threads().max(1);
                for part in jobs.chunks(chunk) {
                    let wd = &self.wdecomp;
                    let objs = &cx.objs;
                    let res = par_map(part, |(o, o2, k, f)| transport_view(wd, rw, objs[o].p, objs[o2].p, f, &words[*k]));
                    for ((o, o2, k, _), comps) in part.iter().zip(res) {
                        let (srcs, tgts) = (&place[&(*o, *k)], &place[&(*o2, *k)]);
                        for (sx, tx, comp) in comps {
                            if out.objs.contains_key(&srcs[sx]) {
                                out.add_map(srcs[sx], tgts[tx], comp);
                            }
                        }
                    }
                }
            }
            // C's slice s − 1 is no longer needed (its components into slice s are transported)
            if si > 0 {
                for &o in &by_h[&hs[si - 1]] {
                    cx.out.remove(&o);
                    for k in 0..nt {
                        place.remove(&(o, k));
                    }
                }
            }
            let dd = self.t_decomp - td0;
            mem.secs_decomp += dd;
            mem.secs_transport += tt.elapsed().as_secs_f64() - dd;
            let sz = out.size();
            if sz.3 > mem.bytes_before {
                (mem.comps_before, mem.entries_before, mem.terms_before, mem.bytes_before) = sz;
                mem.rss_tensor_kb = rss_kb().0;
            }
            mem.peak_objects = mem.peak_objects.max(out.objs.len());
            self.check_degrees(&out);
            let te = Instant::now();
            let last = si + 1 == hs.len();
            if last {
                self.simplify(&mut out, None);
            } else {
                // out complete: every slice before s (its successors in C have h <= s)
                let out_done = |x: usize| slice_of[&x] < hs_cur;
                let in_done = |x: usize| slice_of[&x] <= hs_cur;
                self.simplify(&mut out, Some((&out_done, &in_done)));
            }
            mem.secs_elim += te.elapsed().as_secs_f64();
            mem.stream_passes += 1;
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
                assert!(homogeneous(&m.get(), &la.b, &lb.b, oa.t - ob.t), "component of total q-degree != 0");
            }
        }
    }

    /// Gaussian elimination with scalar pivots in Markowitz order.
    ///
    /// Cost note: every pivot choice rescans all components, so elimination is quadratic in the
    /// size of the complex. That is irrelevant at M6.0 sizes (at most about 130 objects) but must
    /// not be carried over to sengine, which needs an incrementally updated pivot queue.
    /// `elig`: (out(a) complete, in(b) complete) for a streaming pass; None: every pivot.
    #[allow(clippy::type_complexity)]
    fn simplify(&self, cx: &mut XCx<F>, elig: Option<(&dyn Fn(usize) -> bool, &dyn Fn(usize) -> bool)>) {
        // Markowitz order (cost, a, b) with an incrementally updated lazy heap: every edge whose
        // cost (|in(b)| − 1)(|out(a)| − 1) may have changed is pushed again with its current
        // cost, and a popped entry is used only if its cost is still current, so the pivot
        // sequence is that of the full rescan (smallest (cost, a, b) among eligible pivots).
        use std::cmp::Reverse;
        use std::collections::BinaryHeap;
        let ok_a = |a: usize| elig.map_or(true, |(f, _)| f(a));
        let ok_b = |b: usize| elig.map_or(true, |(_, g)| g(b));
        let is_piv = |cx: &XCx<F>, a: usize, b: usize| {
            let (oa, ob) = (cx.objs[&a], cx.objs[&b]);
            oa.p == ob.p && oa.t == ob.t && ok_a(a) && ok_b(b)
        };
        let cost = |cx: &XCx<F>, a: usize, b: usize| (cx.inn[&b].len() - 1) * (cx.out[&a].len() - 1);
        let mut heap: BinaryHeap<Reverse<(usize, usize, usize)>> = BinaryHeap::new();
        for (&a, outs) in &cx.out {
            for &b in outs.keys() {
                if is_piv(cx, a, b) {
                    heap.push(Reverse((cost(cx, a, b), a, b)));
                }
            }
        }
        while let Some(Reverse((c0, a, b))) = heap.pop() {
            if !cx.objs.contains_key(&a) || !cx.objs.contains_key(&b) || !cx.out[&a].contains_key(&b) {
                continue;
            }
            if cost(cx, a, b) != c0 {
                continue; // a current entry was pushed when the cost changed
            }
            let c = scalar_of(&cx.map(a, b)).expect("degree-0 component between equal types is not a scalar");
            assert!(!c.is_zero());
            let f = c.inv().neg();
            let xs: Vec<(usize, PMat<F>)> = cx.inn[&b].iter().filter(|&&x| x != a).map(|&x| (x, cx.map(x, b))).collect();
            let ys: Vec<(usize, PMat<F>)> = cx.out[&a].iter().filter(|(y, _)| **y != b).map(|(&y, m)| (y, m.get())).collect();
            // objects whose out (resp. in) set may change
            let mut xo: BTreeSet<usize> = cx.inn[&a].iter().chain(cx.inn[&b].iter()).copied().collect();
            let mut yi: BTreeSet<usize> = cx.out[&a].keys().chain(cx.out[&b].keys()).copied().collect();
            // fill-in x → y += f · d_ay d_xb: product, sum with the old component and encoding
            // in parallel (each (x, y) once per pivot), stored in pair order
            let pairs: Vec<(usize, usize)> = (0..xs.len()).flat_map(|i| (0..ys.len()).map(move |j| (i, j))).collect();
            let compact = cx.compact;
            for part in pairs.chunks(256) {
                let jobs: Vec<(usize, usize, Option<CMat<F>>)> = part.iter().map(|&(i, j)| (i, j, cx.out[&xs[i].0].get(&ys[j].0).cloned())).collect();
                let work = |(i, j, old): &(usize, usize, Option<CMat<F>>)| -> Option<Option<CMat<F>>> {
                    let m = mm(&ys[*j].1, &xs[*i].1, Some(&f));
                    if m.is_zero() {
                        return None; // unchanged
                    }
                    let sum = match old {
                        Some(o) => o.get().add(&m),
                        None => m,
                    };
                    Some(if sum.is_zero() { None } else { Some(CMat::new(sum, compact)) })
                };
                let res: Vec<Option<Option<CMat<F>>>> = if part.len() >= 4 { par_map(&jobs, work) } else { jobs.iter().map(work).collect() };
                for ((i, j, _), r) in jobs.iter().zip(res) {
                    if let Some(m) = r {
                        cx.set_map(xs[*i].0, ys[*j].0, m);
                    }
                }
            }
            cx.remove_obj(a);
            cx.remove_obj(b);
            xo.remove(&a);
            xo.remove(&b);
            yi.remove(&a);
            yi.remove(&b);
            for &x in &xo {
                for &v in cx.out[&x].keys() {
                    if is_piv(cx, x, v) {
                        heap.push(Reverse((cost(cx, x, v), x, v)));
                    }
                }
            }
            for &y in &yi {
                for &u in &cx.inn[&y] {
                    if is_piv(cx, u, y) {
                        heap.push(Reverse((cost(cx, u, y), u, y)));
                    }
                }
            }
        }
    }

    /// The nonzero components x -> z of d∘d.
    fn d2_components(&self, cx: &XCx<F>) -> Vec<(usize, usize, PMat<F>)> {
        let mut out = vec![];
        for (&x, outs) in &cx.out {
            let mut acc: BTreeMap<usize, PMat<F>> = BTreeMap::new();
            for (&y, f) in outs {
                let f = f.get();
                for (&z, g) in &cx.out[&y] {
                    let m = g.get().matmul(&f);
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
                let fm = solve_f(&f.get(), u, &bases[&oa.p], &bases[&ob.p], &lb.b).map_err(|e| format!("component {} -> {}: {}", a, b, e))?;
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
    let piv = crate::algebra::dense::rref_with_order(&mut aug, &order);
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
    eng.compact = opts.compact;
    eng.transport_rungwise = opts.rung_transport;
    let p0 = eng.identity_type();
    eng.ladder(p0);
    let mut cx = XCx::<F>::new(opts.compact);
    cx.add_obj(XObj { p: p0, h: 0, t: 0 });
    let mut labs = labels.to_vec();
    let mut stages = vec![];
    let mut crossings = vec![];
    let mut last_d2: Vec<(usize, usize, PMat<F>)> = vec![];
    let mut stage_data = BTreeMap::new();
    for (step, g) in word.iter().enumerate() {
        let ts = Instant::now();
        let (c, d) = (labs[0], labs[1]);
        crossings.push((g.positive, c.min(d) as u8));
        eng.decat_step(c, d, g.positive);
        eng.stage_dropped = 0;
        reset_hwm();
        let mut mem = StageMem::default();
        let mut ncx;
        let unsimplified;
        if opts.stream {
            let n_before = eng.stream_inserted;
            ncx = eng.tensor_crossing_stream(std::mem::replace(&mut cx, XCx::new(opts.compact)), &labs, g.positive, &mut mem);
            unsimplified = eng.stream_inserted - n_before;
        } else {
            let td0 = eng.t_decomp;
            let tt = Instant::now();
            ncx = eng.tensor_crossing(&cx, &labs, g.positive);
            drop(std::mem::replace(&mut cx, XCx::new(opts.compact)));
            mem.secs_decomp = eng.t_decomp - td0;
            mem.secs_transport = tt.elapsed().as_secs_f64() - mem.secs_decomp;
            unsimplified = ncx.objs.len();
            (mem.comps_before, mem.entries_before, mem.terms_before, mem.bytes_before) = ncx.size();
            mem.rss_tensor_kb = rss_kb().0;
            mem.peak_objects = unsimplified;
            eng.check_degrees(&ncx);
            let te = Instant::now();
            eng.simplify(&mut ncx, None);
            mem.secs_elim = te.elapsed().as_secs_f64();
        }
        let (_, ea, ta, ba) = ncx.size();
        (mem.entries_after, mem.terms_after, mem.bytes_after) = (ea, ta, ba);
        let mut d2_ideal = 0;
        let t2 = Instant::now();
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
        mem.secs_d2 = t2.elapsed().as_secs_f64();
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
        let mut info = StageInfo {
            step: step + 1,
            labels: labs.clone(),
            unsimplified,
            objects: ncx.objs.len(),
            components: ncx.components(),
            lower_bound: lb,
            secs: ts.elapsed().as_secs_f64(),
            dropped: eng.stage_dropped,
            d2_ideal,
            mem,
        };
        let (r, h) = rss_kb();
        info.mem.rss_end_kb = r;
        info.mem.rss_peak_kb = h;
        info.mem.bytes_caches = eng.cache_bytes();
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
            let m = &info.mem;
            eprintln!(
                "[xengine]   mem: before simplify {} comps / {} entries / {} terms / {:.1} MB; after {} entries / {} terms / {:.1} MB; caches {:.1} MB; RSS tensor {:.0} MB end {:.0} MB peak {:.0} MB; secs decomp {:.2} transport {:.2} elim {:.2} d2 {:.2}",
                m.comps_before, m.entries_before, m.terms_before, m.bytes_before as f64 / 1e6, m.entries_after, m.terms_after, m.bytes_after as f64 / 1e6, m.bytes_caches as f64 / 1e6,
                m.rss_tensor_kb as f64 / 1e3, m.rss_end_kb as f64 / 1e3, m.rss_peak_kb as f64 / 1e3, m.secs_decomp, m.secs_transport, m.secs_elim, m.secs_d2
            );
            if opts.stream {
                eprintln!("[xengine]   stream: {} passes, at most {} objects held", m.stream_passes, m.peak_objects);
            }
        }
        if opts.close_stages && step + 1 >= 2 && step + 1 < word.len() && labs == labels {
            let tc = Instant::now();
            let d2 = if opts.check_d2 { last_d2.clone() } else { vec![] };
            let sd = eng.closure_data(&ncx, crossings.clone(), opts.seed, &d2).map_err(|e| format!("closure data of stage {}: {}", step + 1, e))?;
            if opts.verbose {
                eprintln!("[xengine]   stage {} closure data: {} objects, {} edges in {:.2}s", step + 1, sd.0.objs.len(), sd.0.edges.len(), tc.elapsed().as_secs_f64());
            }
            if !opts.stage_close_ns.is_empty() {
                for marked in [None, Some(0)] {
                    for (nn, p) in close_guarded(&sd.0, &sd.1, opts.prune_n, &opts.stage_close_ns, marked) {
                        eprintln!("[xengine] stage {} N = {}{}: rank {}: {}", step + 1, nn, if marked.is_some() { " reduced" } else { "" }, p.total_rank(), p);
                    }
                }
            }
            stage_data.insert(step + 1, sd);
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
        stage_data,
    })
}

impl<F: Field> XScan<F> {
    /// gl_N Poincaré polynomials for every N in `ns` (reduced: block 0 marked). Colors > N give 0.
    /// A pruned scan (`prune_n = Some(N0)`) is valid for N <= N0 only (panics otherwise); there
    /// the closure of every final d∘d component is first checked to vanish (panics otherwise).
    pub fn close(&self, ns: &[usize], reduced: bool) -> BTreeMap<usize, Poincare> {
        close_guarded(&self.data, &self.final_d2, self.prune_n, ns, if reduced { Some(0) } else { None })
    }

    /// As [`XScan::close`], with Err instead of a panic when the exact oracle refuses (memory
    /// cap, unsupported shape) or a pruned scan is asked for N above its pruning level.
    pub fn try_close(&self, ns: &[usize], reduced: bool) -> Result<BTreeMap<usize, Poincare>, String> {
        try_close_guarded(&self.data, &self.final_d2, self.prune_n, ns, if reduced { Some(0) } else { None })
    }

    /// Stages whose closure data were kept (`XOpts::close_stages`).
    pub fn closed_stages(&self) -> Vec<usize> {
        self.stage_data.keys().copied().collect()
    }

    /// gl_N Poincaré polynomials of the closure of the first m crossings (the braid word
    /// `word[..m]`), from the stage data kept with `XOpts::close_stages`; m = word.len() is the
    /// final complex. The guards of `close` apply: a pruned scan refuses N > prune_n, and the
    /// stage-m d∘d components must close to 0. Err if stage m was not kept.
    pub fn close_stage(&self, m: usize, ns: &[usize], reduced: bool) -> Result<BTreeMap<usize, Poincare>, String> {
        self.close_stage_marked(m, ns, if reduced { Some(0) } else { None })
    }

    /// `close_stage` with an arbitrary marked block (None: unreduced; Some(i): the block at
    /// bottom position i is marked). For a knot the reduced homology does not depend on the
    /// marked block (the base-point check: Some(1) = Some(0)).
    pub fn close_stage_marked(&self, m: usize, ns: &[usize], marked: Option<usize>) -> Result<BTreeMap<usize, Poincare>, String> {
        if m == self.word.len() {
            if let Some(n0) = self.prune_n {
                if let Some(&nn) = ns.iter().find(|&&nn| nn > n0) {
                    return Err(format!("a scan pruned at N = {} gives the gl_N homology only for N <= {} (requested N = {})", n0, n0, nn));
                }
            }
            return Ok(close_guarded(&self.data, &self.final_d2, self.prune_n, ns, marked));
        }
        if let Some(n0) = self.prune_n {
            if let Some(&nn) = ns.iter().find(|&&nn| nn > n0) {
                return Err(format!("a scan pruned at N = {} gives the gl_N homology only for N <= {} (requested N = {})", n0, n0, nn));
            }
        }
        let (data, d2) = self.stage_data.get(&m).ok_or_else(|| format!("stage {} has no closure data (close_stages off, m < 2, or labels do not close up)", m))?;
        Ok(close_guarded(data, d2, self.prune_n, ns, marked))
    }
}

/// `XScan::close` on given closure data and final d∘d components.
#[allow(clippy::type_complexity)]
fn close_guarded<F: Field>(data: &ClosureDataC<F>, final_d2: &[(usize, usize, Vec<Vec<EPoly<F>>>)], prune_n: Option<usize>, ns: &[usize], marked: Option<usize>) -> BTreeMap<usize, Poincare> {
    ns.iter()
        .map(|&nn| {
            if let Some(n0) = prune_n {
                assert!(nn <= n0, "a scan pruned at N = {} gives the gl_N homology only for N <= {} (requested N = {})", n0, n0, nn);
                for (a, b, f) in final_d2 {
                    assert!(closure_kills_map(data, *a, *b, f, nn, marked), "pruned scan at N = {}: the d∘d component {} -> {} does not close to 0 at gl_{}", n0, a, b, nn);
                }
            }
            (nn, close_exact_colored(data, nn, marked))
        })
        .collect()
}

/// As [`close_guarded`], with Err instead of a panic when the exact oracle refuses (memory
/// estimate above the cap, an unsupported A_λ shape such as top weight > 255) or when a pruned
/// scan is asked for N above its pruning level. The pruned-scan d∘d guard still panics (it is a
/// certification check, unchanged).
#[allow(clippy::type_complexity)]
fn try_close_guarded<F: Field>(data: &ClosureDataC<F>, final_d2: &[(usize, usize, Vec<Vec<EPoly<F>>>)], prune_n: Option<usize>, ns: &[usize], marked: Option<usize>) -> Result<BTreeMap<usize, Poincare>, String> {
    let mut out = BTreeMap::new();
    for &nn in ns {
        if let Some(n0) = prune_n {
            if nn > n0 {
                return Err(format!("a scan pruned at N = {} gives the gl_N homology only for N <= {} (requested N = {})", n0, n0, nn));
            }
            for (a, b, f) in final_d2 {
                assert!(closure_kills_map(data, *a, *b, f, nn, marked), "pruned scan at N = {}: the d∘d component {} -> {} does not close to 0 at gl_{}", n0, a, b, nn);
            }
        }
        out.insert(nn, try_close_exact_colored(data, nn, marked)?);
    }
    Ok(out)
}

impl<F: Field> XScan<F> {
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
/// `reduced`: the block at bottom position 0 is marked (reduced unknot = q^0). Err if the scan
/// fails or the exact closure refuses (memory estimate above `KRHOM_MEM_CAP_MB`, or an A_λ
/// shape the oracle does not support, e.g. top weight a(N − a) > 255); never panics on those.
pub fn khr_singular_exact<F: Field>(labels: &[usize], word: &[Gen], ns: &[usize], reduced: bool) -> Result<BTreeMap<usize, Poincare>, String> {
    let scan = xscan::<F>(labels, word, &XOpts::default())?;
    scan.try_close(ns, reduced)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algebra::field::F31;

    #[test]
    fn sort_product_equals_matmul() {
        let mut rng = StdRng::seed_from_u64(5);
        let rnd = |rng: &mut StdRng, nr: usize, nc: usize| {
            let mut m = PMat::<F31>::zero(nr, nc);
            for col in m.cols.iter_mut() {
                for i in 0..nr {
                    if rng.gen_bool(0.5) {
                        let mut t: Vec<(crate::algebra::mpoly::Mono, F31)> = (0..rng.gen_range(1..4)).map(|_| ((rng.gen_range(0..3u128)) << (8 * rng.gen_range(0..4)), F31::from_i64(rng.gen_range(-3..4)))).collect();
                        t.sort_by_key(|x| x.0);
                        t.dedup_by_key(|x| x.0);
                        t.retain(|x| !x.1.is_zero());
                        if !t.is_empty() {
                            col.push((i as u32, Poly { t }));
                        }
                    }
                }
            }
            m
        };
        for _ in 0..200 {
            let (n, k, m) = (rng.gen_range(1..6), rng.gen_range(1..6), rng.gen_range(1..6));
            let a = rnd(&mut rng, n, k);
            let b = rnd(&mut rng, k, m);
            let c = F31::from_i64(rng.gen_range(1..5));
            assert_eq!(mm(&a, &b, None).cols, a.matmul(&b).cols);
            assert_eq!(mm(&a, &b, Some(&c)).cols, scale_const(&a.matmul(&b), &c).cols);
        }
    }
}
