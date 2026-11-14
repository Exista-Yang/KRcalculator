//! M5a: singular indecomposables and merge decompositions with Foam_N pruning built in
//! (SPEC_M4_M6a §5, SPEC_3C_GL6 rev1 D3/D4/D8/D11/D14 and §2 "M5a", design note
//! gl6/R2_m5prune.md).
//!
//! * **Types.** L_p = the minimal-flow canonical ladder of p (`schur::CanonicalBasis::
//!   ladder_minflow`, [`canonical_ladder`]). A kept type whose ladder is decomposable
//!   (char L_p ≠ b_p) needs an idempotent e_p; M5a does not build idempotents (M5i does): `ty`
//!   returns [`DecErr::NeedsIdempotent`] for such a type, and [`idem_needed`] lists them.
//! * **Pruning** (D3). `SIndecOpts::prune_n` is fixed per builder; `killed(q)` is
//!   `Schur::killed_at(q, N)` and nothing else. Killed summands get no ι, π, no ladder and no
//!   further rungs. A frozen [`DecData`] built at prune_n serves every n' ≤ prune_n through
//!   [`DecData::restrict`].
//! * **Decompositions** rung by rung: one rung by Hom pairs from the provider (through the
//!   [`DecHom`] adapter, the only code that touches `HomProvider`), re-certified here (form at
//!   rep(q) at 2 points, genuineness, ends, degrees, D14 degree scaling), normalized with the
//!   inverse intersection form and orthogonalized among the kept summands by a Neumann series;
//!   composites by Kron and composition. Release asserts: the kept multiset equals the kept part
//!   of `term_mult` (D4), rank bookkeeping, Π∘Ι = id at `check_points` points on every coset of
//!   the support of X; Σ ιπ = e_X when nothing is killed.
//! * **Rickard components** π_τ ∘ (id_{L_p} ⊗ d_j) ∘ ι_σ, nonzero ones only, with the
//!   cell-rule assert on the Rickard terms.
//! * **Layer A** (R2 §4.1): certified provider results (pairs and closure generators) on disk,
//!   one file per target, keyed by λ, field, table strings and canonical-ladder rung words; the
//!   header carries the stamps (format, ladder function, `TMPL_ALG`, provider stamp) and a
//!   mismatch is rejected (recomputed). Loaded maps go through the same re-certification as
//!   fresh ones.
//!
//! Degree-scaling law (D14, measured in M5a): for a map m of normalized degree d between
//! expressions with normalization shifts sh(src), sh(tgt), every localized entry at every coset
//! satisfies v(λz) = λ^u v(z) with 2u = d − sh(tgt) + sh(src) (the unnormalized polynomial
//! degree, roots of degree 2). [`check_scaling`] asserts it at the top-stratum coset.

use crate::coset::CosetSpace;
use crate::field::Field;
use crate::homprov::{assert_genuine, certify_closure, check_idempotent, form_at, Budget, Chain, CostRec, CostTable, Exact, HomProvider, MissKind, Search, COST_ALG};
use crate::schur::{crossing_terms, flows, Lp, Schur, Tid};
use crate::sexpr::Expr;
use crate::smap::{points, Arena, Eval, ExprId, MapId, MapNode, Point, Reduced};
use crate::tmpl::{TmplLib, TMPL_ALG};
use crate::web::Rung;
use std::collections::{BTreeMap, HashMap};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

/// format stamp of layer-A files
pub const SINDEC_VERSION: &str = "krhom-sindec v2";
/// stamp of the ladder function (`canonical_ladder` = schur's minimal-flow ladder: min flows,
/// then length, then lexicographic)
pub const LADDER_STAMP: &str = "schur.ladder_minflow.v1";
/// M5i: the idempotent of a decomposable canonical ladder (`SIndec::build_idem`: e_p = id − Σ
/// ι_r π_r from certified pairs of the extra summands); enters the layer-A keys of maps that
/// contain an e_q / e_p (eq=, ep=)
pub const M5I_STAMP: &str = "m5i-peel1";

// ------------------------------------------------------------------ types

pub struct SType {
    pub p: Tid,
    pub lad: ExprId,
    pub e: Option<MapId>,
    pub rank: u128,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SSummand {
    pub q: Tid,
    pub s: i32,
    pub iota: MapId,
    pub pi: MapId,
}

#[derive(Clone, Debug)]
pub struct Decomp {
    pub p: Tid,
    pub rungs: Vec<Rung>,
    pub x: ExprId,
    pub summands: Vec<SSummand>,
    /// ι = π = id (X is the canonical ladder of its only summand)
    pub identity: bool,
    /// every summand is killed (no provider call); M6a drops the object
    pub all_killed: bool,
    /// killed summands (q, s, multiplicity), dropped
    pub killed: Vec<(Tid, i32, usize)>,
}

#[derive(Clone, Debug)]
pub struct ClosureGens {
    pub p: Tid,
    pub ins: Vec<MapId>,
    pub outs: Vec<MapId>,
}

#[derive(Clone, Debug)]
pub struct SIndecOpts {
    pub cache_dir: Option<PathBuf>,
    pub budget: Budget,
    /// points of the Π∘Ι and form checks (2)
    pub check_points: usize,
    pub seed: u64,
    /// Foam_N pruning (None = unpruned), fixed for the life of the builder
    pub prune_n: Option<u8>,
    /// M5c route stamp: a layer-A file is accepted if its header equals the current one up to
    /// the cost-table token (`costs=ALG:RECSrecs:HASH`, compared as `costs=ALG`); the provider
    /// route (chain, caps, slowdown, solver generation, ladder and template stamps, key) must still
    /// match. Sound because every loaded map is re-certified (`target_pairs`: ends, degrees, D14
    /// scaling, then Π∘Ι, form and genuineness in `dec_one`; `closure_gens`: `certify_closure`),
    /// whichever table found it. Files are always written with the full stamp. Default false (the
    /// M5b behaviour: exact header match).
    pub layer_a_route_stamp: bool,
}
impl Default for SIndecOpts {
    fn default() -> Self {
        SIndecOpts { cache_dir: None, budget: Budget::default(), check_points: 2, seed: 0x5d5e_c001, prune_n: None, layer_a_route_stamp: false }
    }
}

#[derive(Clone, Debug, Default)]
pub struct DecStats {
    pub decomps: usize,
    pub identity: usize,
    pub all_killed: usize,
    pub killed_summands: usize,
    pub targets: usize,
    pub sum_mult: usize,
    pub provider_calls: usize,
    pub provider_secs: f64,
    pub max_x_flows: u64,
    pub idempotents: usize,
    pub layer_a_loaded: usize,
    /// of the loaded: accepted by the route stamp only (another cost table, M5c)
    pub layer_a_route: usize,
    pub layer_a_saved: usize,
    pub layer_a_rejected: usize,
    pub rickard: usize,
    pub closure: usize,
    /// entries checked by the degree-scaling check (D14), and of those the ones 0 at both points
    pub scaling_entries: usize,
    pub scaling_zero: usize,
    /// seconds spent in M5's own checks (forms, Π∘Ι, scaling)
    pub check_secs: f64,
}

/// Errors of sindec (SPEC_3C_GL6 rev1 §2 M5a).
#[derive(Clone, Debug, PartialEq)]
pub enum DecErr {
    /// the type is killed at N (B_p = 0 in Foam_N)
    Killed { p: String, n: u8 },
    /// a kept type with a decomposable canonical ladder (M5i)
    NeedsIdempotent { p: String },
    /// a DecHom guard declined the call before the provider was asked
    Refused(String),
    /// the provider failed; engineering = Time | Memory | Refused (not Exhausted)
    Provider { engineering: bool, text: String },
    /// a certification check of M5 failed
    Check(String),
    /// frozen data: the item was not precomputed
    NotPrecomputed(String),
    /// invalid request (e.g. restrict above prune_n)
    Invalid(String),
    /// layer-A cache error
    Cache(String),
}
impl DecErr {
    /// the same error with a context prefix on its text (variant, hence the failure kind, kept)
    pub fn ctx(self, c: &str) -> DecErr {
        match self {
            DecErr::Refused(t) => DecErr::Refused(format!("{}: {}", c, t)),
            DecErr::Provider { engineering, text } => DecErr::Provider { engineering, text: format!("{}: {}", c, text) },
            DecErr::Check(t) => DecErr::Check(format!("{}: {}", c, t)),
            DecErr::Cache(t) => DecErr::Cache(format!("{}: {}", c, t)),
            e => e,
        }
    }
}
impl std::fmt::Display for DecErr {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DecErr::Killed { p, n } => write!(f, "type {} is killed at N = {}", p, n),
            DecErr::NeedsIdempotent { p } => write!(f, "type {}: canonical ladder decomposable, idempotent required (M5i)", p),
            DecErr::Refused(s) => write!(f, "DecHom refused: {}", s),
            DecErr::Provider { engineering, text } => write!(f, "provider failed ({}): {}", if *engineering { "engineering" } else { "exhausted" }, text),
            DecErr::Check(s) => write!(f, "check failed: {}", s),
            DecErr::NotPrecomputed(s) => write!(f, "not precomputed: {}", s),
            DecErr::Invalid(s) => write!(f, "invalid request: {}", s),
            DecErr::Cache(s) => write!(f, "layer-A cache: {}", s),
        }
    }
}
impl From<DecErr> for String {
    fn from(e: DecErr) -> String {
        e.to_string()
    }
}

// ------------------------------------------------------------------ ladders and helpers

/// The canonical ladder of q (the ladder function shared by `ty`, the DecHom guard and
/// `idem_needed`): schur's minimal-flow ladder of q in its block (λ, μ(q)).
pub fn canonical_ladder(s: &mut Schur, q: Tid) -> Vec<Rung> {
    if let Some(w) = forced_ladder(s, q) {
        return w;
    }
    let mu = s.table(q).cols();
    s.canonical_basis(&mu).ladder_minflow[&q].clone()
}

/// M5i validation only: `KRHOM_FORCE_LADDER="TABLE=WORD;..."` replaces the canonical ladder of
/// the type TABLE (smapio table format) of the λ given by the word by WORD (smapio rung format),
/// e.g. a decomposable word, so that the engine runs through e_p on a case with a known answer.
/// The word must end at the block of the type (asserted).
pub fn forced_ladder(s: &mut Schur, q: Tid) -> Option<Vec<Rung>> {
    let v = std::env::var("KRHOM_FORCE_LADDER").ok()?;
    let name = s.table(q).fmt();
    for item in v.split(';') {
        let Some((t, w)) = item.split_once('=') else { continue };
        if t.trim() == name {
            let w = crate::smapio::parse_rungs(w.trim()).unwrap_or_else(|e| panic!("KRHOM_FORCE_LADDER: {}", e));
            let top = Expr::new(&s.lam, &w).map(|e| e.top()).unwrap_or_else(|e| panic!("KRHOM_FORCE_LADDER {}: {}", name, e));
            assert_eq!(top, s.table(q).cols(), "KRHOM_FORCE_LADDER {}: word ends at {:?}", name, top);
            return Some(w);
        }
    }
    None
}

/// expand(ladder_char(word)) == [(q, 1)]: the word is indecomposable of type q
pub fn is_indecomposable_of(s: &mut Schur, w: &[Rung], q: Tid) -> bool {
    let ch = s.ladder_char(w);
    let ex = s.expand(&ch);
    ex.len() == 1 && ex[0].0 == q && ex[0].1 == Lp::one()
}

/// the smallest N at which q is kept: killed_at(q, N) ⇔ N < kill_level(q)
pub fn kill_level(s: &Schur, q: Tid) -> usize {
    let t = s.table(q);
    let m = t.rows().iter().chain(t.cols().iter()).map(|&a| a as usize).max().unwrap_or(0);
    let c = t.cell_weight().first().map_or(0, |&a| a as usize);
    m.max(c)
}

fn rkey(rs: &[Rung]) -> Vec<(usize, u8, usize)> {
    rs.iter().map(|r| (r.i, if r.dir == crate::web::Dir::E { 0 } else { 1 }, r.r)).collect()
}

fn coset_of_rep(lam: &[u8], s: &Schur, t: Tid) -> u64 {
    let tab = s.table(t).clone();
    let cs = CosetSpace::new(lam, &tab.cols());
    cs.rank(&cs.rep(&tab))
}

fn fnv64(s: &str) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for b in s.bytes() {
        h ^= b as u64;
        h = h.wrapping_mul(0x0100_0000_01b3);
    }
    h
}

fn word_of<F: Field>(ar: &Arena<F>, e: ExprId) -> Vec<Rung> {
    ar.expr_of(e).rungs.clone()
}

/// λ^e in F (e may be negative)
fn fpow<F: Field>(l: &F, e: i32) -> F {
    let b = if e < 0 { l.inv() } else { l.clone() };
    let mut r = F::one();
    for _ in 0..e.unsigned_abs() {
        r = r.mul(&b);
    }
    r
}

/// The scaling factor of the D14 check (large multiplicative order in both primes).
pub const SCALE_LAMBDA: i64 = 1_234_567;

/// D14 degree-scaling check of one map at the coset c: every entry satisfies
/// v(λz) = λ^u v(z) with 2u = deg − sh(tgt) + sh(src), where `deg` is the declared normalized
/// degree. Returns (entries checked, entries 0 at both points); Err on a violation (also when the
/// unnormalized degree is odd and an entry is nonzero).
pub fn check_scaling<F: Field>(ar: &Arena<F>, m: MapId, c: u64, deg: i32, z: &Point<F>) -> Result<(usize, usize), String> {
    let info = ar.info(m);
    let u2 = deg - ar.expr_of(info.tgt).shift() + ar.expr_of(info.src).shift();
    let lam = F::from_i64(SCALE_LAMBDA);
    let lz = Point { id: z.id ^ 0x8000_0000, x: z.x.iter().map(|v| v.mul(&lam)).collect() };
    let b1 = Eval::new(ar, z.clone()).block(m, c);
    let b2 = Eval::new(ar, lz).block(m, c);
    let mut zero = 0;
    let f = if u2 % 2 == 0 { Some(fpow(&lam, u2 / 2)) } else { None };
    for (a, b) in b1.a.iter().zip(&b2.a) {
        if a.is_zero() && b.is_zero() {
            zero += 1;
            continue;
        }
        match &f {
            Some(f) if a.mul(f) == *b => {}
            _ => return Err(format!("degree scaling: map {} (declared degree {}) at coset {}: an entry does not scale by λ^{}/2", m, deg, c, u2)),
        }
    }
    Ok((b1.a.len(), zero))
}

// ------------------------------------------------------------------ DecHom adapter (R2 §5.2)

#[derive(Clone, Debug)]
pub struct PairTarget {
    pub x: ExprId,
    pub q: Tid,
    pub lq: ExprId,
    pub eq: Option<MapId>,
    pub s: i32,
    pub mult: usize,
}
#[derive(Clone, Debug)]
pub struct ProvOut {
    pub maps_in: Vec<MapId>,
    pub maps_out: Vec<MapId>,
    pub secs: f64,
    pub summary: String,
}
#[derive(Clone, Debug)]
pub struct ProvFail {
    /// Time | Memory | Refused (or a guard), not Exhausted
    pub engineering: bool,
    /// a DecHom guard (precondition or output check) failed, not the provider
    pub guard: bool,
    pub text: String,
}

pub trait DecHom<F: Field> {
    fn pairs(&mut self, ar: &mut Arena<F>, s: &mut Schur, t: &PairTarget, b: &Budget) -> Result<ProvOut, ProvFail>;
    #[allow(clippy::too_many_arguments)]
    fn closure(&mut self, ar: &mut Arena<F>, s: &mut Schur, lp: ExprId, p: Tid, ep: Option<MapId>, outs: bool, b: &Budget) -> Result<ProvOut, ProvFail>;
    fn end0(&mut self, ar: &mut Arena<F>, s: &mut Schur, l: ExprId, b: &Budget) -> Result<ProvOut, ProvFail>;
}

fn guard(text: String) -> ProvFail {
    ProvFail { engineering: true, guard: true, text }
}

/// The D4 guards of `DecHom::pairs` before any provider work (review F1/F2): lq is the
/// canonical ladder of q; without e_q, L_q is indecomposable of type q; with e_q, the e_q
/// contract (`check_idempotent`); mult = [v^s] of the multiplicity of b_q in char X; the tops
/// agree.
pub fn dechom_precheck<F: Field>(ar: &Arena<F>, s: &mut Schur, t: &PairTarget) -> Result<(), String> {
    let lam = ar.lam().to_vec();
    if s.lam != lam {
        return Err(format!("Schur lambda {:?} != arena lambda {:?}", s.lam, lam));
    }
    let (xw, qw) = (word_of(ar, t.x), word_of(ar, t.lq));
    let canon = canonical_ladder(s, t.q);
    if qw != canon {
        return Err(format!("wrong lq: L_q word {} is not the canonical ladder {} of q = {}", crate::smapio::rungs_str(&qw), crate::smapio::rungs_str(&canon), s.table(t.q).fmt()));
    }
    let mu = s.table(t.q).cols();
    if Expr::new(&lam, &xw).map(|e| e.top()).ok() != Some(mu.clone()) {
        return Err(format!("X does not end at the block {:?} of q", mu));
    }
    match t.eq {
        None => {
            if !is_indecomposable_of(s, &qw, t.q) {
                return Err(format!("L_q decomposable: e_q required (q = {})", s.table(t.q).fmt()));
            }
        }
        Some(e) => check_idempotent(ar, s, t.lq, t.q, e, 0x5eed_d0c0).map_err(|m| format!("e_q rejected: {}", m))?,
    }
    let chx = s.ladder_char(&xw);
    let m = s.expand(&chx).into_iter().find(|u| u.0 == t.q).map(|u| u.1.coef(t.s)).unwrap_or(0);
    if m != t.mult as i128 {
        return Err(format!("wrong mult: {} != [v^{}] of the multiplicity of b_q in char X = {}", t.mult, t.s, m));
    }
    Ok(())
}

/// output checks of `DecHom::pairs`: mult maps each way, ι: lq → X of degree s and π: X → lq
/// of degree −s (read from `ar.info`), and the D14 scaling law at the top-stratum coset rep(q)
pub fn dechom_postcheck<F: Field>(ar: &Arena<F>, s: &Schur, t: &PairTarget, ins: &[MapId], outs: &[MapId], seed: u64) -> Result<(usize, usize), String> {
    if ins.len() != t.mult || outs.len() != t.mult {
        return Err(format!("provider returned {} ι and {} π, mult {}", ins.len(), outs.len(), t.mult));
    }
    for &m in ins {
        let i = ar.info(m);
        if (i.src, i.tgt, i.deg) != (t.lq, t.x, t.s) {
            return Err(format!("wrong degree or ends: ι {} is {}->{} of degree {}, expected {}->{} degree {}", m, i.src, i.tgt, i.deg, t.lq, t.x, t.s));
        }
    }
    for &m in outs {
        let i = ar.info(m);
        if (i.src, i.tgt, i.deg) != (t.x, t.lq, -t.s) {
            return Err(format!("wrong degree or ends: π {} is {}->{} of degree {}, expected {}->{} degree {}", m, i.src, i.tgt, i.deg, t.x, t.lq, -t.s));
        }
    }
    let c = coset_of_rep(ar.lam(), s, t.q);
    let z = points::<F>(ar.k(), seed ^ 0x5ca1e, 1, 0).remove(0);
    let (mut n, mut zr) = (0, 0);
    for &m in ins {
        let (a, b) = check_scaling(ar, m, c, t.s, &z)?;
        n += a;
        zr += b;
    }
    for &m in outs {
        let (a, b) = check_scaling(ar, m, c, -t.s, &z)?;
        n += a;
        zr += b;
    }
    Ok((n, zr))
}

fn prov_fail(ms: Vec<crate::homprov::Missing>) -> ProvFail {
    // R2 §5.2: Missing.kind matched with a wildcard; a failure is a completeness failure iff some
    // provider of the chain ended Exhausted, otherwise (Time | Memory | Refused only) engineering
    let engineering = !ms.iter().any(|m| m.kind == MissKind::Exhausted);
    let text = ms.iter().map(|m| format!("{} [{:?}, degree {}, {}/{}]", m.what, m.kind, m.degree, m.got, m.want)).collect::<Vec<_>>().join("; ");
    ProvFail { engineering, guard: false, text }
}

impl<F: Field, P: HomProvider<F> + ?Sized> DecHom<F> for P {
    fn pairs(&mut self, ar: &mut Arena<F>, s: &mut Schur, t: &PairTarget, b: &Budget) -> Result<ProvOut, ProvFail> {
        dechom_precheck(ar, s, t).map_err(guard)?;
        let t0 = Instant::now();
        let r = HomProvider::pairs(self, ar, s, t.x, t.q, t.lq, t.eq, t.s, t.mult, b);
        let secs = t0.elapsed().as_secs_f64();
        match r {
            Ok((p, st)) => {
                dechom_postcheck(ar, s, t, &p.iotas, &p.pis, 0x5eed_d0c1).map_err(guard)?;
                Ok(ProvOut { maps_in: p.iotas, maps_out: p.pis, secs, summary: format!("{} tiers {:?} exact {} detect {}", self.name(), st.tier_hits, st.exact_calls, st.detect) })
            }
            Err(ms) => Err(prov_fail(ms)),
        }
    }
    fn closure(&mut self, ar: &mut Arena<F>, s: &mut Schur, lp: ExprId, p: Tid, ep: Option<MapId>, outs: bool, b: &Budget) -> Result<ProvOut, ProvFail> {
        let w = word_of(ar, lp);
        if w != canonical_ladder(s, p) {
            return Err(guard(format!("closure: L_p is not the canonical ladder of {}", s.table(p).fmt())));
        }
        if ep.is_none() && !is_indecomposable_of(s, &w, p) {
            return Err(guard(format!("closure: L_p decomposable: e_p required ({})", s.table(p).fmt())));
        }
        let t0 = Instant::now();
        match HomProvider::closure(self, ar, s, lp, p, ep, outs, b) {
            Ok((ms, st)) => Ok(ProvOut { maps_in: if outs { vec![] } else { ms.clone() }, maps_out: if outs { ms } else { vec![] }, secs: t0.elapsed().as_secs_f64(), summary: format!("{} tiers {:?}", self.name(), st.tier_hits) }),
            Err(ms) => Err(prov_fail(ms)),
        }
    }
    fn end0(&mut self, ar: &mut Arena<F>, s: &mut Schur, l: ExprId, b: &Budget) -> Result<ProvOut, ProvFail> {
        let t0 = Instant::now();
        match HomProvider::end0(self, ar, s, l, b) {
            Ok((ms, st)) => Ok(ProvOut { maps_in: ms, maps_out: vec![], secs: t0.elapsed().as_secs_f64(), summary: format!("{} tiers {:?}", self.name(), st.tier_hits) }),
            Err(ms) => Err(prov_fail(ms)),
        }
    }
}

/// A layer-A header with the cost-table token reduced to its solver generation
/// (`costs=hs-kl1:273recs:b8c4…` -> `costs=hs-kl1`), without the C2-7 route token `c27=…`, for
/// `SIndecOpts::layer_a_route_stamp`.
pub fn route_header(h: &str) -> String {
    h.split(' ')
        // C2-7: the `c27=` token only steers which Leaves words / routes are attempted
        .filter(|t| !t.starts_with("c27="))
        .map(|t| match t.strip_prefix("costs=") {
            Some(v) => format!("costs={}", v.split(':').next().unwrap_or("")),
            None => t.to_string(),
        })
        .collect::<Vec<_>>()
        .join(" ")
}

// ------------------------------------------------------------------ provider construction (D8)

/// per-call time limit of Exact at label 3 (admission by cost records only)
pub const EXACT_L3_MAX_CALL_S: f64 = 600.0;

/// A copy of the cost table with every time multiplied by `slowdown` (machine factor).
pub fn scaled_costs(t: &CostTable, slowdown: f64) -> CostTable {
    let mut out = CostTable::new();
    for r in &t.recs {
        out.record(CostRec { secs: r.secs * slowdown, ..r.clone() });
    }
    out
}

/// The provider of λ and its stamp (D8). Max label ≤ 2: `Chain::for_labels` (Exact first).
/// Max label ≥ 3: an own chain [Search S0-S2, Exact { max_label: 3 } with the scaled cost
/// table, Search S3-S6] (`Chain::for_labels` closes Exact for labels ≥ 3).
pub fn provider_for_stamped<F: Field>(lam: &[u8], costs: Option<&Path>, slowdown: f64) -> (Box<dyn HomProvider<F>>, String) {
    let maxl = lam.iter().copied().max().unwrap_or(0);
    let (table, tstamp) = match costs {
        None => (None, "none".to_string()),
        Some(p) => match CostTable::load(p) {
            Ok(t) => {
                let st = format!("{}:{}recs:{:016x}", t.alg, t.len(), fnv64(&std::fs::read_to_string(p).unwrap_or_default()));
                (Some(scaled_costs(&t, slowdown)), st)
            }
            Err(e) => {
                eprintln!("[sindec] cost table {} not used: {}", p.display(), e);
                (None, "rejected".to_string())
            }
        },
    };
    let stamp = format!("homprov[{}] maxl={} costs={} slowdown={}", COST_ALG, maxl, tstamp, slowdown);
    if maxl <= 2 {
        let ex = Exact { costs: table, ..Exact::default() };
        return (Box::new(Chain::<F>::for_labels(maxl, ex)), format!("{} chain=Exact,S0-6", stamp));
    }
    let mut s_lo = Search::tiers(0, 2);
    s_lo.costs = table.clone();
    // Exact at label 3 only through measured cost records (the label rule's caps are kept); its
    // per-call ceiling is the libbuild watchdog wall (600 s): the (3,3,3) closure calls of the
    // Markov corner measure 16-340 s and are otherwise vetoed (M5b_measure.md)
    let ex = Exact { max_label: 3, costs: table.clone(), max_call_s: EXACT_L3_MAX_CALL_S, ..Exact::default() };
    let mut s_hi = Search::tiers(3, 6);
    s_hi.costs = table;
    (Box::new(Chain::<F>(vec![Box::new(s_lo), Box::new(ex), Box::new(s_hi)])), format!("{} chain=S0-2,Exact(l3,{}s),S3-6", stamp, EXACT_L3_MAX_CALL_S))
}

pub fn provider_for<F: Field>(lam: &[u8], costs: Option<&Path>, slowdown: f64) -> Box<dyn HomProvider<F>> {
    provider_for_stamped(lam, costs, slowdown).0
}

// ------------------------------------------------------------------ the builder

pub struct SIndec<F: Field> {
    pub schur: Schur,
    pub ar: Arena<F>,
    prov: Box<dyn HomProvider<F>>,
    prov_stamp: String,
    opts: SIndecOpts,
    lam: Vec<u8>,
    types: HashMap<Tid, SType>,
    /// every ladder type built (kept or not): (L_q, e_q, rank B_q); killed types appear here only
    /// as extra summands of a decomposable canonical ladder (M5i)
    tall: HashMap<Tid, (ExprId, Option<MapId>, u128)>,
    decs: HashMap<(Tid, Vec<(usize, u8, usize)>), Arc<Decomp>>,
    rick: HashMap<(Tid, usize, bool, usize), Arc<Vec<(usize, usize, MapId)>>>,
    clos: HashMap<Tid, Arc<ClosureGens>>,
    stats: DecStats,
    /// per idempotent type: (table, construction, seconds)
    pub idem_log: Vec<(String, String, f64)>,
    /// (p, rung) source-idempotent decompositions in progress (commutation reroute; cycle guard)
    busy: std::collections::HashSet<(Tid, (usize, u8, usize))>,
    /// per source-idempotent decomposition: (p, rung, route) for the report
    pub route_log: Vec<(String, String, String)>,
    /// per fresh provider call: (p table, q table or "closure-ins|outs", s, mult, seconds,
    /// provider summary with tier hits), for the report
    pub target_log: Vec<(String, String, i32, usize, f64, String)>,
    /// when set, every `target_log` entry is also appended at once to this file (one JSON line),
    /// so that a job killed by a watchdog keeps the timings of its finished targets
    pub target_log_file: Option<PathBuf>,
}

impl<F: Field> SIndec<F> {
    /// Builder with its own Arena (template library in `opts.cache_dir`/tmpl when given). The
    /// provider stamp enters every layer-A header.
    pub fn new(lam: &[u8], provider: Box<dyn HomProvider<F>>, opts: SIndecOpts) -> Self {
        SIndec::with_stamp(lam, provider, "unstamped".into(), opts)
    }
    pub fn with_stamp(lam: &[u8], provider: Box<dyn HomProvider<F>>, prov_stamp: String, opts: SIndecOpts) -> Self {
        let tdir = opts.cache_dir.as_ref().map(|d| d.join("tmpl"));
        let lib = TmplLib::new(tdir);
        if let Some(d) = &opts.cache_dir {
            let _ = std::fs::create_dir_all(d.join("layer_a"));
        }
        SIndec {
            schur: Schur::new(lam),
            ar: Arena::new(lam, lib),
            prov: provider,
            prov_stamp,
            opts,
            lam: lam.to_vec(),
            types: HashMap::new(),
            tall: HashMap::new(),
            decs: HashMap::new(),
            rick: HashMap::new(),
            clos: HashMap::new(),
            stats: DecStats::default(),
            idem_log: vec![],
            busy: std::collections::HashSet::new(),
            route_log: vec![],
            target_log: vec![],
            target_log_file: None,
        }
    }
    pub fn prune_n(&self) -> Option<u8> {
        self.opts.prune_n
    }
    pub fn stats(&self) -> DecStats {
        self.stats.clone()
    }
    pub fn killed(&self, q: Tid) -> bool {
        self.opts.prune_n.is_some_and(|n| self.schur.killed_at(q, n as usize))
    }
    fn tfmt(&self, q: Tid) -> String {
        self.schur.table(q).fmt()
    }
    fn pts(&self, n: usize, salt: u64) -> Vec<Point<F>> {
        points::<F>(self.ar.k(), self.opts.seed ^ salt, n, 0)
    }

    /// canonical ladder of p, interned, with its idempotent e_p when the ladder is decomposable
    /// (M5i); Err for a killed p
    pub fn ty(&mut self, p: Tid) -> Result<&SType, DecErr> {
        if let Some(n) = self.opts.prune_n {
            if self.schur.killed_at(p, n as usize) {
                return Err(DecErr::Killed { p: self.tfmt(p), n });
            }
        }
        if self.types.contains_key(&p) {
            return Ok(&self.types[&p]);
        }
        let (lad, e, rank) = self.ty_full(p)?;
        if let Some(n) = self.opts.prune_n {
            // R2 §2.1 (debug check): on <= 3 strands a kept type has a minimal-flow ladder with labels <= N
            if self.lam.len() <= 3 {
                let w = word_of(&self.ar, lad);
                let ex = Expr::new(&self.lam, &w).map_err(DecErr::Invalid)?;
                debug_assert!(ex.levels().iter().all(|l| l.iter().all(|&a| a <= n)), "kept type {} has a ladder label > N = {}", self.tfmt(p), n);
            }
        }
        self.types.insert(p, SType { p, lad, e, rank });
        Ok(&self.types[&p])
    }

    /// (L_q, e_q, rank B_q) of any type, killed or not (no pruning check): the canonical ladder,
    /// and for a decomposable one the idempotent of [`SIndec::build_idem`]
    fn ty_full(&mut self, p: Tid) -> Result<(ExprId, Option<MapId>, u128), DecErr> {
        if let Some(&t) = self.tall.get(&p) {
            return Ok(t);
        }
        let w = canonical_ladder(&mut self.schur, p);
        let mu = self.schur.table(p).cols();
        let rank = self.schur.rank(&mu, p);
        let lad = self.ar.expr(Expr::new(&self.lam, &w).map_err(DecErr::Invalid)?);
        let e = if is_indecomposable_of(&mut self.schur, &w, p) {
            assert_eq!(flows(&self.lam, &w), rank, "type {}: flows of the canonical ladder != rank of B_p", self.tfmt(p));
            None
        } else {
            if std::env::var("KRHOM_NO_M5I").is_ok() {
                return Err(DecErr::NeedsIdempotent { p: self.tfmt(p) });
            }
            Some(self.build_idem(p, lad, &w)?)
        };
        self.tall.insert(p, (lad, e, rank));
        Ok((lad, e, rank))
    }

    /// M5i: the idempotent e_p of B_p inside a decomposable canonical ladder L_p (char L_p =
    /// b_p + Σ v^s b_r), C2_2_report §2 step 1. For every extra summand (r, s, cnt): certified
    /// pairs (ι_r: L_r → L_p of degree s, π_r of degree −s) from the provider (layer A, kind=idem;
    /// e_r passed when L_r is itself decomposable), normalized by the inverse form at rep(r),
    /// orthogonalized by the Neumann series; e_p = id − Σ ι_r π_r. Checks: Π∘Ι = δ (e_r on the
    /// diagonal) at `check_points` points on every coset; `check_idempotent` (e² = e and rank =
    /// c_t(b_p)(1) at rep(t) of every table of char L_p, top scalar 1) at 2 points; e² = e on
    /// every coset at a third point with Σ_c rank = rank B_p.
    fn build_idem(&mut self, p: Tid, lad: ExprId, w: &[Rung]) -> Result<MapId, DecErr> {
        let t0 = Instant::now();
        // prefix propagation: L_p = L_{p'} ⊗ R with L_{p'} the (decomposable) canonical ladder of a
        // type p' and b_{p'}·R = b_p: then e_p = Kron(e_{p'}, R) (its image is B_{p'} ⊗ R ≅ B_p)
        if let Some((pp, e)) = self.idem_by_prefix(p, w)? {
            let what = format!("M5i e_p of {} (prefix {})", self.tfmt(p), self.tfmt(pp));
            self.idem_checks(p, lad, e, &what)?;
            self.stats.idempotents += 1;
            self.idem_log.push((self.tfmt(p), format!("Kron(e of {}, {})", self.tfmt(pp), crate::smapio::rung_str(&w[w.len() - 1])), t0.elapsed().as_secs_f64()));
            return Ok(e);
        }
        let ch = self.schur.ladder_char(w);
        let ex = self.schur.expand(&ch);
        let mut tg: Vec<(Tid, i32, usize)> = vec![];
        for (q, m) in &ex {
            if *q == p {
                if *m != Lp::one() {
                    return Err(DecErr::Invalid(format!("M5i: b_p occurs in char L_p of {} with multiplicity {}", self.tfmt(p), m.fmt_ec())));
                }
                continue;
            }
            for (i, &c) in m.c.iter().enumerate() {
                if c != 0 {
                    if c < 0 {
                        return Err(DecErr::Invalid("M5i: negative multiplicity in char L_p".into()));
                    }
                    tg.push((*q, m.lo + i as i32, c as usize));
                }
            }
        }
        if tg.is_empty() {
            return Err(DecErr::Invalid(format!("M5i: L_p of {} has no extra summand", self.tfmt(p))));
        }
        tg.sort_by(|a, b| self.schur.lmax_of(b.0).cmp(&self.schur.lmax_of(a.0)).then(a.0.cmp(&b.0)).then(a.1.cmp(&b.1)));
        let mut sums: Vec<SSummand> = vec![];
        let mut desc = vec![];
        for (q, s, cnt) in tg {
            let (lq, eq, _) = self.ty_full(q)?;
            desc.push(format!("{}{{v^{}:{}}}", self.tfmt(q), s, cnt));
            let t = PairTarget { x: lad, q, lq, eq, s, mult: cnt };
            let key = self.idem_key(p, &t);
            let (ins, outs) = self.target_pairs_keyed(key, p, &t)?;
            let g = self.form_twice(&ins, &outs, q, s, p)?;
            let inv = crate::local::inverse(&g).ok_or_else(|| DecErr::Check(format!("M5i: form of {} at v^{} in L_p of {} singular", self.tfmt(q), s, self.tfmt(p))))?;
            for a in 0..cnt {
                let terms: Vec<(F, MapId)> = (0..cnt).filter(|&b| !inv[a][b].is_zero()).map(|b| (inv[a][b].clone(), outs[b])).collect();
                let pi = if terms.len() == 1 && terms[0].0.is_one() { terms[0].1 } else { self.ar.add(MapNode::Lin(terms), "m5i pi'").map_err(DecErr::Check)? };
                sums.push(SSummand { q, s, iota: ins[a], pi });
            }
        }
        let sums = self.orthogonalize(lad, sums)?;
        let what = format!("M5i e_p of {}", self.tfmt(p));
        self.check_pi_iota(lad, &sums, None, false, &what)?;
        let id = self.id_map(lad)?;
        let mut terms: Vec<(F, MapId)> = vec![(F::one(), id)];
        let mone = F::zero().sub(&F::one());
        for sm in &sums {
            let ip = self.ar.add(MapNode::Comp(sm.iota, sm.pi), "m5i iota pi").map_err(DecErr::Check)?;
            terms.push((mone.clone(), ip));
        }
        let e = self.ar.add(MapNode::Lin(terms), "m5i e_p").map_err(DecErr::Check)?;
        self.idem_checks(p, lad, e, &what)?;
        self.stats.idempotents += 1;
        self.idem_log.push((self.tfmt(p), desc.join(" + "), t0.elapsed().as_secs_f64()));
        Ok(e)
    }

    /// prefix propagation of [`SIndec::build_idem`]: Some((p', Kron(e_{p'}, R))) when w = L_{p'}·R
    /// with L_{p'} the canonical ladder of a decomposable type p' and term_mult(p', R) = b_p
    fn idem_by_prefix(&mut self, p: Tid, w: &[Rung]) -> Result<Option<(Tid, MapId)>, DecErr> {
        if w.len() < 2 || std::env::var("KRHOM_M5I_NO_PREFIX").is_ok() {
            return Ok(None);
        }
        let (w1, r) = (&w[..w.len() - 1], w[w.len() - 1]);
        let ch = self.schur.ladder_char(w1);
        let cands: Vec<Tid> = self.schur.expand(&ch).into_iter().map(|t| t.0).collect();
        for pp in cands {
            if canonical_ladder(&mut self.schur, pp) != w1 || is_indecomposable_of(&mut self.schur, w1, pp) {
                continue;
            }
            let tm = self.schur.term_mult(pp, &[r]);
            if tm.len() != 1 || tm[0].0 != p || tm[0].1 != Lp::one() {
                continue;
            }
            let (_, epp, _) = self.ty_full(pp)?;
            let Some(epp) = epp else { continue };
            let e = self.ar.add(MapNode::Kron(epp, vec![r]), "m5i Kron(e_p', R)").map_err(DecErr::Check)?;
            return Ok(Some((pp, e)));
        }
        Ok(None)
    }

    /// the checks of an idempotent e_p of L_p (spec §5.2 step 6): `check_idempotent` at 2 points
    /// (e² = e and rank c_t(b_p)(1) at rep(t) of every table of char L_p, top scalar 1), then e² = e
    /// on every coset at a third point with Σ_c rank = rank B_p
    fn idem_checks(&mut self, p: Tid, lad: ExprId, e: MapId, what: &str) -> Result<(), DecErr> {
        let tc = Instant::now();
        check_idempotent(&self.ar, &mut self.schur, lad, p, e, self.opts.seed ^ 0x1de3).map_err(|m| DecErr::Check(format!("{}: {}", what, m)))?;
        let pt = self.pts(3, 0x1de3_e5a1).remove(2);
        let mut ev = Eval::new(&self.ar, pt);
        let cos = self.ar.flows(lad).cosets.clone();
        let mut tot: u128 = 0;
        for &c in &cos {
            let b = ev.block(e, c);
            if b.mul(&b) != *b {
                return Err(DecErr::Check(format!("{}: e² != e at coset {} (third point)", what, c)));
            }
            let n = b.rows as usize;
            if n > 0 {
                tot += crate::dense::rank(&b.to_rows(), n) as u128;
            }
        }
        let mu = self.schur.table(p).cols();
        let rk = self.schur.rank(&mu, p);
        if tot != rk {
            return Err(DecErr::Check(format!("{}: Σ_c rank(e_c) = {} != rank B_p = {}", what, tot, rk)));
        }
        self.stats.check_secs += tc.elapsed().as_secs_f64();
        Ok(())
    }

    /// the form π_j∘ι_i at rep(q) at check_points (≥ 2) points, equal at all of them
    fn form_twice(&mut self, ins: &[MapId], outs: &[MapId], q: Tid, s: i32, p: Tid) -> Result<Vec<Vec<F>>, DecErr> {
        let tc = Instant::now();
        let c = coset_of_rep(&self.lam, &self.schur, q);
        let pts = self.pts(self.opts.check_points.max(2), 0xf0f0);
        let g = form_at(&self.ar, ins, outs, c, &pts[0]).map_err(DecErr::Check)?;
        for pt in &pts[1..] {
            if form_at(&self.ar, ins, outs, c, pt).map_err(DecErr::Check)? != g {
                return Err(DecErr::Check(format!("form of {} at v^{} in {} differs at a second point", self.tfmt(q), s, self.tfmt(p))));
            }
        }
        for &m in ins.iter().chain(outs) {
            assert_genuine(&self.ar, m).map_err(DecErr::Check)?;
        }
        self.stats.check_secs += tc.elapsed().as_secs_f64();
        Ok(g)
    }

    fn killed_rank(&mut self, p: Tid, rs: &[Rung]) -> u128 {
        let Some(n) = self.opts.prune_n else { return 0 };
        let mut r: u128 = 0;
        for (q, m) in self.schur.term_mult(p, rs) {
            if self.schur.killed_at(q, n as usize) {
                let mu = self.schur.table(q).cols();
                r += m.at1() as u128 * self.schur.rank(&mu, q);
            }
        }
        r
    }

    /// kept and killed parts of term_mult(p, rs) as (q, s) -> multiplicity
    fn split_mult(&mut self, p: Tid, rs: &[Rung]) -> (BTreeMap<(Tid, i32), usize>, BTreeMap<(Tid, i32), usize>) {
        let (mut kept, mut kill) = (BTreeMap::new(), BTreeMap::new());
        for (q, m) in self.schur.term_mult(p, rs) {
            let k = self.killed(q);
            for (i, &c) in m.c.iter().enumerate() {
                if c != 0 {
                    assert!(c > 0, "negative multiplicity");
                    let e = if k { &mut kill } else { &mut kept };
                    e.insert((q, m.lo + i as i32), c as usize);
                }
            }
        }
        (kept, kill)
    }

    pub fn dec_rung(&mut self, p: Tid, r: Rung) -> Result<Arc<Decomp>, DecErr> {
        self.dec_word(p, &[r])
    }

    /// Decomposition of L_p ⊗ (rungs), rung by rung, cached.
    pub fn dec_word(&mut self, p: Tid, rs: &[Rung]) -> Result<Arc<Decomp>, DecErr> {
        let key = (p, rkey(rs));
        if let Some(d) = self.decs.get(&key) {
            return Ok(d.clone());
        }
        let (lad, ep) = {
            let t = self.ty(p)?;
            (t.lad, t.e)
        };
        let d = if rs.is_empty() {
            let m = match ep {
                Some(e) => e,
                None => self.id_map(lad)?,
            };
            Decomp { p, rungs: vec![], x: lad, summands: vec![SSummand { q: p, s: 0, iota: m, pi: m }], identity: ep.is_none(), all_killed: false, killed: vec![] }
        } else if rs.len() == 1 {
            self.dec_one(p, rs[0])?
        } else {
            self.dec_composite(p, rs)?
        };
        let d = Arc::new(d);
        self.decs.insert(key, d.clone());
        Ok(d)
    }

    fn id_map(&mut self, e: ExprId) -> Result<MapId, DecErr> {
        self.ar.add(MapNode::Id(e), "sindec id").map_err(DecErr::Check)
    }

    fn x_of(&mut self, p: Tid, rs: &[Rung]) -> Result<(ExprId, Vec<Rung>), DecErr> {
        let lad = self.ty(p)?.lad;
        let mut w = word_of(&self.ar, lad);
        w.extend_from_slice(rs);
        let e = Expr::new(&self.lam, &w).map_err(DecErr::Invalid)?;
        Ok((self.ar.expr(e), w))
    }

    fn killed_list(kill: &BTreeMap<(Tid, i32), usize>) -> Vec<(Tid, i32, usize)> {
        kill.iter().map(|(&(q, s), &c)| (q, s, c)).collect()
    }

    /// the source idempotent Kron(e_p, rs) of a decomposable L_p (None otherwise)
    fn src_idem(&mut self, p: Tid, rs: &[Rung]) -> Result<Option<MapId>, DecErr> {
        let Some(e) = self.ty(p)?.e else { return Ok(None) };
        if rs.is_empty() {
            return Ok(Some(e));
        }
        Ok(Some(self.ar.add(MapNode::Kron(e, rs.to_vec()), "m5i Kron(e_p, R)").map_err(DecErr::Check)?))
    }

    /// rank of B_p ⊗ rs: flows of X when L_p = B_p, else Σ term_mult ranks
    fn src_rank(&mut self, p: Tid, rs: &[Rung], xflows: u128) -> Result<u128, DecErr> {
        if self.ty(p)?.e.is_none() {
            return Ok(xflows);
        }
        let mut r: u128 = 0;
        for (q, m) in self.schur.term_mult(p, rs) {
            let mu = self.schur.table(q).cols();
            r += m.at1() as u128 * self.schur.rank(&mu, q);
        }
        Ok(r)
    }

    fn dec_one(&mut self, p: Tid, r: Rung) -> Result<Decomp, DecErr> {
        let rs = [r];
        let (x, xw) = self.x_of(p, &rs)?;
        let xflows = flows(&self.lam, &xw);
        let mult = self.schur.term_mult(p, &rs);
        let (kept, kill) = self.split_mult(p, &rs);
        let killed_rank = self.killed_rank(p, &rs);
        let srank = self.src_rank(p, &rs, xflows)?;
        let ep = self.ty(p)?.e;
        self.stats.killed_summands += kill.values().sum::<usize>();
        if kept.is_empty() {
            assert_eq!(killed_rank, srank, "all summands of {} ⊗ {} killed: killed ranks {} != rank {}", self.tfmt(p), crate::smapio::rung_str(&r), killed_rank, srank);
            self.stats.all_killed += 1;
            return Ok(Decomp { p, rungs: rs.to_vec(), x, summands: vec![], identity: false, all_killed: true, killed: Self::killed_list(&kill) });
        }
        // identity shortcut: X is literally the canonical ladder of its only summand (e_p = e_q = None)
        if ep.is_none() && mult.len() == 1 && mult[0].1 == Lp::one() {
            let q = mult[0].0;
            if canonical_ladder(&mut self.schur, q) == xw && self.ty(q)?.e.is_none() {
                let lq = self.ty(q)?.lad;
                assert_eq!(lq, x, "identity shortcut: interned ladder of q != X");
                let id = self.id_map(x)?;
                self.stats.identity += 1;
                self.stats.decomps += 1;
                return Ok(Decomp { p, rungs: rs.to_vec(), x, summands: vec![SSummand { q, s: 0, iota: id, pi: id }], identity: true, all_killed: false, killed: vec![] });
            }
        }
        self.stats.max_x_flows = self.stats.max_x_flows.max(xflows as u64);
        // source idempotent (M5i): the provider sees X = L_p ⊗ R with the multiplicities of char X
        let ex = self.src_idem(p, &rs)?;
        let xmult: BTreeMap<(Tid, i32), usize> = if ex.is_some() {
            let chx = self.schur.ladder_char(&xw);
            let mut m = BTreeMap::new();
            for (q, c) in self.schur.expand(&chx) {
                for (i, &v) in c.c.iter().enumerate() {
                    if v > 0 {
                        m.insert((q, c.lo + i as i32), v as usize);
                    }
                }
            }
            m
        } else {
            BTreeMap::new()
        };
        // targets in SPEC order: l_max(q) descending, then s ascending
        let mut tg: Vec<(Tid, i32, usize)> = kept.iter().map(|(&(q, s), &c)| (q, s, c)).collect();
        tg.sort_by(|a, b| self.schur.lmax_of(b.0).cmp(&self.schur.lmax_of(a.0)).then(a.0.cmp(&b.0)).then(a.1.cmp(&b.1)));
        // M5i: X equal to the canonical ladder of a decomposable type xq (its extra-summand pairs
        // are those of e_xq); a decomposable source: the commutation reroute family
        let xq: Option<Tid> = {
            let ch = self.schur.ladder_char(&xw);
            let mut f = None;
            for (q, _) in self.schur.expand(&ch) {
                if canonical_ladder(&mut self.schur, q) == xw && !is_indecomposable_of(&mut self.schur, &xw, q) {
                    f = Some(q);
                }
            }
            f
        };
        let mut fam_tried = ex.is_some() && xq.is_none();
        let mut fam = if fam_tried { self.reroute_family(p, r).map_err(|e| e.ctx(&format!("M5i reroute of {} ⊗ {}", self.tfmt(p), crate::smapio::rung_str(&r))))? } else { None };
        let mut sums: Vec<SSummand> = vec![];
        for (q, s, cnt) in tg {
            let (lq, eq) = {
                let t = self.ty(q)?;
                (t.lad, t.e)
            };
            let mx = match ex {
                None => cnt,
                Some(_) => *xmult.get(&(q, s)).ok_or_else(|| DecErr::Check(format!("M5i: b_q of {} at v^{} in term_mult but not in char X", self.tfmt(q), s)))?,
            };
            if mx < cnt {
                return Err(DecErr::Check(format!("M5i: multiplicity {} in char X < {} in term_mult", mx, cnt)));
            }
            let t = PairTarget { x, q, lq, eq, s, mult: mx };
            let comm_q = if let (Some(_), false, 0, 1) = (eq, lq == x, s, mx) { self.comm_chain(&xw, &word_of(&self.ar, lq))? } else { None };
            let (mut ins, mut outs) = if let (Some(e), true, 0, 1) = (eq, lq == x, s, mx) {
                // X is the (decomposable) canonical ladder of q itself: ι = π = e_q
                (vec![e], vec![e])
            } else if let (Some(e), Some((c, cb))) = (eq, comm_q) {
                // X is commutation-equivalent to L_q (decomposable): ι = c⁻¹∘e_q, π = e_q∘c
                let io = self.ar.add(MapNode::Comp(cb, e), "m5i comm iota").map_err(DecErr::Check)?;
                let pi = self.ar.add(MapNode::Comp(e, c), "m5i comm pi").map_err(DecErr::Check)?;
                (vec![io], vec![pi])
            } else if let Some(fam) = fam.as_ref() {
                // source idempotent: the family of X from the commutation reroute
                let (a, b) = fam.get(&(q, s)).cloned().unwrap_or_default();
                if a.len() != mx || b.len() != mx {
                    return Err(DecErr::Check(format!("M5i reroute of {} ⊗ {}: {} pairs for {} at v^{}, char X has {}", self.tfmt(p), crate::smapio::rung_str(&r), a.len(), self.tfmt(q), s, mx)));
                }
                (a, b)
            } else if let Some(xq) = xq {
                // X = L_{xq} (decomposable): the extra-summand pairs of e_{xq} (layer A kind=idem)
                let key = self.idem_key(xq, &t);
                self.target_pairs_keyed(key, p, &t)?
            } else {
                let what = format!("target {} at v^{} of {} ⊗ {}", self.tfmt(q), s, self.tfmt(p), crate::smapio::rung_str(&r));
                match self.target_pairs(p, &rs, &t) {
                    Ok(x) => x,
                    Err(err) if !fam_tried && !matches!(err, DecErr::Refused(_)) && std::env::var("KRHOM_NO_REROUTE_FALLBACK").is_err() => {
                        // C2-4: a target the provider cannot reach (an idempotent target has no leaf
                        // in Search's tree; a non-local splitting leaves stuck branches): the
                        // commutation reroute of X, decomposed from an earlier prefix type
                        fam_tried = true;
                        fam = self.reroute_family(p, r).map_err(|e| e.ctx(&format!("{} (reroute after: {})", what, err)))?;
                        let Some(f) = fam.as_ref() else { return Err(err.ctx(&format!("{} (no reroute)", what))) };
                        let (a, b) = f.get(&(q, s)).cloned().unwrap_or_default();
                        if a.len() != mx || b.len() != mx {
                            return Err(DecErr::Check(format!("M5i reroute: {}: {} pairs, char X has {}", what, a.len(), mx)));
                        }
                        (a, b)
                    }
                    Err(err) => return Err(err.ctx(&what)),
                }
            };
            if let Some(e) = ex {
                // apply the source idempotent: ι ↦ Kron(e_p, R)∘ι, π ↦ π∘Kron(e_p, R)
                for m in ins.iter_mut() {
                    *m = self.ar.add(MapNode::Comp(e, *m), "m5i E iota").map_err(DecErr::Check)?;
                }
                for m in outs.iter_mut() {
                    *m = self.ar.add(MapNode::Comp(*m, e), "m5i pi E").map_err(DecErr::Check)?;
                }
            }
            // M5's own re-certification (R2 §2.2 item 4): form at rep(q) at 2 points
            let g = self.form_twice(&ins, &outs, q, s, p)?;
            let (ins, outs, g) = if ex.is_some() {
                // cnt independent columns I and rows J of the form E-projected form (rank = cnt)
                let mut gr = g.clone();
                let icols = crate::dense::rref(&mut gr, mx);
                if icols.len() != cnt {
                    return Err(DecErr::Check(format!("M5i: projected form of {} at v^{} in {} has rank {} != {}", self.tfmt(q), s, self.tfmt(p), icols.len(), cnt)));
                }
                let mut gt: Vec<Vec<F>> = (0..mx).map(|i| (0..mx).map(|j| g[j][i].clone()).collect()).collect();
                let jrows = crate::dense::rref(&mut gt, mx);
                let gs: Vec<Vec<F>> = jrows.iter().map(|&j| icols.iter().map(|&i| g[j][i].clone()).collect()).collect();
                (icols.iter().map(|&i| ins[i]).collect::<Vec<_>>(), jrows.iter().map(|&j| outs[j]).collect::<Vec<_>>(), gs)
            } else {
                (ins, outs, g)
            };
            let inv = crate::local::inverse(&g).ok_or_else(|| DecErr::Check(format!("form of {} at v^{} singular (rank < {})", self.tfmt(q), s, cnt)))?;
            // π'_a = Σ_b inv[a][b] π_b
            for a in 0..cnt {
                let terms: Vec<(F, MapId)> = (0..cnt).filter(|&b| !inv[a][b].is_zero()).map(|b| (inv[a][b].clone(), outs[b])).collect();
                let pi = if terms.len() == 1 && terms[0].0.is_one() { terms[0].1 } else { self.ar.add(MapNode::Lin(terms), "sindec pi'").map_err(DecErr::Check)? };
                sums.push(SSummand { q, s, iota: ins[a], pi });
            }
            self.stats.targets += 1;
            self.stats.sum_mult += cnt;
        }
        // rank bookkeeping
        let mut total: u128 = 0;
        for sm in &sums {
            total += self.types[&sm.q].rank;
        }
        assert_eq!(total + killed_rank, srank, "decomposition of {} ⊗ {}: summand ranks {} + killed {} != rank {}", self.tfmt(p), crate::smapio::rung_str(&r), total, killed_rank, srank);
        let sums = self.orthogonalize(x, sums)?;
        let full = kill.is_empty();
        self.check_pi_iota(x, &sums, ex, full, &format!("{} ⊗ {}", self.tfmt(p), crate::smapio::rung_str(&r)))?;
        self.stats.decomps += 1;
        Ok(Decomp { p, rungs: rs.to_vec(), x, summands: sums, identity: false, all_killed: false, killed: Self::killed_list(&kill) })
    }

    /// the commutation isomorphism of rungs i, i+1 of w (commuting rungs): (w → w', w' → w),
    /// mutually inverse, from the degree-0 window templates (as homprov's Search)
    fn swap_iso(&mut self, w: &[Rung], i: usize) -> Result<(MapId, MapId, Vec<Rung>), DecErr> {
        let lam = self.lam.clone();
        let (r1, r2) = (w[i], w[i + 1]);
        let lo = r1.i.min(r2.i);
        let hi = r1.i.max(r2.i) + 1;
        let labels = crate::homprov::level_of(&lam, w, i)[lo..=hi].to_vec();
        let sh = |r: Rung| Rung { i: r.i - lo, ..r };
        let rel = vec![sh(r1), sh(r2)];
        let swr = vec![rel[1], rel[0]];
        let (wa, wb) = (crate::tmpl::Window::new(&labels, &rel), crate::tmpl::Window::new(&labels, &swr));
        let tf = self.ar.lib().try_hom_basis_cap(&wa, &wb, 0, u128::MAX).map_err(DecErr::Check)?;
        let tb = self.ar.lib().try_hom_basis_cap(&wb, &wa, 0, u128::MAX).map_err(DecErr::Check)?;
        if tf.len() != 1 || tb.len() != 1 {
            return Err(DecErr::Check(format!("commutation {:?} {}: Hom^0 dims {} / {}", labels, crate::smapio::rungs_str(&rel), tf.len(), tb.len())));
        }
        let mut sw = w.to_vec();
        sw.swap(i, i + 1);
        let a = self.ar.expr(Expr::new(&lam, w).map_err(DecErr::Invalid)?);
        let b = self.ar.expr(Expr::new(&lam, &sw).map_err(DecErr::Invalid)?);
        let f = self.ar.add(MapNode::Local { src: a, tgt: b, at: i, first: lo, tmpl: tf[0] }, "m5i commute").map_err(DecErr::Check)?;
        let g = self.ar.add(MapNode::Local { src: b, tgt: a, at: i, first: lo, tmpl: tb[0] }, "m5i commute back").map_err(DecErr::Check)?;
        // g∘f = c·id (one-dimensional End^0 of the window word): normalize
        let gf = self.ar.add(MapNode::Comp(g, f), "m5i commute gf").map_err(DecErr::Check)?;
        let pt = self.pts(1, 0xc0aa).remove(0);
        let mut ev = Eval::new(&self.ar, pt);
        let cos = self.ar.flows(a).cosets.clone();
        let mut c: Option<F> = None;
        for &cc in &cos {
            let blk = ev.block(gf, cc);
            let n = blk.rows as usize;
            if n == 0 {
                continue;
            }
            let c0 = blk.a[0].clone();
            let mut want = crate::smap::Block::<F>::zeros(n, n);
            for k in 0..n {
                want.a[k * n + k] = c0.clone();
            }
            if *blk != want {
                return Err(DecErr::Check(format!("commutation {:?} {}: g∘f is not scalar", labels, crate::smapio::rungs_str(&rel))));
            }
            if c.is_none() {
                c = Some(c0);
            } else if c.as_ref() != Some(&c0) {
                return Err(DecErr::Check("commutation: g∘f scalar differs between cosets".into()));
            }
        }
        let c = c.filter(|x| !x.is_zero()).ok_or_else(|| DecErr::Check("commutation: g∘f = 0".into()))?;
        let g = self.ar.add(MapNode::Lin(vec![(c.inv(), g)]), "m5i commute inverse").map_err(DecErr::Check)?;
        Ok((f, g, sw))
    }

    /// (from → to, to → from) mutually inverse commutation isomorphisms when `to` is in the
    /// commutation class of `from` (breadth first over adjacent commuting swaps, ≤ 4 096 words)
    fn comm_chain(&mut self, from: &[Rung], to: &[Rung]) -> Result<Option<(MapId, MapId)>, DecErr> {
        if from.len() != to.len() || from == to {
            return Ok(None);
        }
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
            if prev.len() > 4096 {
                break;
            }
            for i in 0..cur.len().saturating_sub(1) {
                if crate::homprov::rungs_commute(&cur[i], &cur[i + 1]) {
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
            return Ok(None);
        }
        let mut steps = vec![];
        let mut cur = to.to_vec();
        while cur != from {
            let (p, i) = prev[&cur].clone();
            steps.push((p.clone(), i));
            cur = p;
        }
        steps.reverse();
        let (mut fwd, mut bwd): (Option<MapId>, Option<MapId>) = (None, None);
        for (w, i) in steps {
            let (f, g, _) = self.swap_iso(&w, i)?;
            fwd = Some(match fwd {
                None => f,
                Some(x) => self.ar.add(MapNode::Comp(f, x), "m5i chain").map_err(DecErr::Check)?,
            });
            bwd = Some(match bwd {
                None => g,
                Some(x) => self.ar.add(MapNode::Comp(x, g), "m5i chain back").map_err(DecErr::Check)?,
            });
        }
        Ok(Some((fwd.unwrap(), bwd.unwrap())))
    }

    /// M5i, decomposable source (C2-4): the decomposition of X = L_p·R through the commutation
    /// class: R commutes with the tail w[j..] of w = L_p (j < |w| minimal), so X ≅ Y = w[..j]·R·w[j..];
    /// w[..j] is the canonical ladder of a type a, and Y is decomposed rung by rung from a
    /// (`dec_word(a, R·w[j..])`, the ordinary machinery, which never meets the non-local splitting
    /// of L_p). Returns per kept (q, s) the family (ι: L_q → X, π: X → L_q) transported through the
    /// commutation isomorphism, Π∘Ι = δ: a basis of the multiplicity space of b_q⟨s⟩ in char X.
    /// None if R does not commute past the last rung of w or w[..j] is not a canonical ladder.
    #[allow(clippy::type_complexity)]
    fn reroute_family(&mut self, p: Tid, r: Rung) -> Result<Option<BTreeMap<(Tid, i32), (Vec<MapId>, Vec<MapId>)>>, DecErr> {
        if std::env::var("KRHOM_M5I_NO_REROUTE").is_ok() {
            return Ok(None);
        }
        let lad = self.ty(p)?.lad;
        let w = word_of(&self.ar, lad);
        let mut j = w.len();
        while j > 0 && crate::homprov::rungs_commute(&w[j - 1], &r) {
            j -= 1;
        }
        if j == w.len() {
            return Ok(None);
        }
        let w1 = w[..j].to_vec();
        let ch = self.schur.ladder_char(&w1);
        let a = self.schur.expand(&ch).into_iter().map(|t| t.0).find(|&q| canonical_ladder(&mut self.schur, q) == w1);
        let Some(a) = a else { return Ok(None) };
        let key = (p, (r.i, r.dir as u8, r.r));
        if !self.busy.insert(key) {
            return Err(DecErr::Check(format!("M5i reroute cycle at {} ⊗ {}", self.tfmt(p), crate::smapio::rung_str(&r))));
        }
        let mut rs2 = vec![r];
        rs2.extend_from_slice(&w[j..]);
        let d = self.dec_word(a, &rs2);
        self.busy.remove(&key);
        let d = d.map_err(|e| e.ctx(&format!("rerouted dec_word({}, {})", self.tfmt(a), crate::smapio::rungs_str(&rs2))))?;
        // commutation chain X = w·R → Y (R moved from position |w| to j)
        let mut cur = w.clone();
        cur.push(r);
        let (mut fwd, mut bwd): (Option<MapId>, Option<MapId>) = (None, None);
        for i in (j..w.len()).rev() {
            let (f, g, sw) = self.swap_iso(&cur, i)?;
            fwd = Some(match fwd {
                None => f,
                Some(x) => self.ar.add(MapNode::Comp(f, x), "m5i chain").map_err(DecErr::Check)?,
            });
            bwd = Some(match bwd {
                None => g,
                Some(x) => self.ar.add(MapNode::Comp(x, g), "m5i chain back").map_err(DecErr::Check)?,
            });
            cur = sw;
        }
        let (c, cb) = (fwd.unwrap(), bwd.unwrap());
        if self.ar.info(c).tgt != d.x {
            return Err(DecErr::Check("M5i reroute: the commuted word is not the decomposed X".into()));
        }
        let mut fam: BTreeMap<(Tid, i32), (Vec<MapId>, Vec<MapId>)> = BTreeMap::new();
        for sm in &d.summands {
            let io = self.ar.add(MapNode::Comp(cb, sm.iota), "m5i reroute iota").map_err(DecErr::Check)?;
            let pi = self.ar.add(MapNode::Comp(sm.pi, c), "m5i reroute pi").map_err(DecErr::Check)?;
            let e = fam.entry((sm.q, sm.s)).or_default();
            e.0.push(io);
            e.1.push(pi);
        }
        self.route_log.push((self.tfmt(p), crate::smapio::rung_str(&r), format!("dec_word({}, {})", self.tfmt(a), crate::smapio::rungs_str(&rs2))));
        Ok(Some(fam))
    }

    /// layer-A lookup or a provider call for one target
    fn target_pairs(&mut self, p: Tid, rs: &[Rung], t: &PairTarget) -> Result<(Vec<MapId>, Vec<MapId>), DecErr> {
        let key = self.pairs_key(p, rs, t);
        self.target_pairs_keyed(key, p, t)
    }

    fn target_pairs_keyed(&mut self, key: String, p: Tid, t: &PairTarget) -> Result<(Vec<MapId>, Vec<MapId>), DecErr> {
        if let Some(v) = self.layer_a_load(&key, 2 * t.mult)? {
            let (ins, outs) = (v[..t.mult].to_vec(), v[t.mult..].to_vec());
            // loaded maps pass the same output checks as fresh ones (ends, degrees, D14 scaling)
            let tc = Instant::now();
            let (n, z) = dechom_postcheck(&self.ar, &self.schur, t, &ins, &outs, self.opts.seed).map_err(|e| DecErr::Cache(format!("loaded target rejected: {}", e)))?;
            self.stats.scaling_entries += n;
            self.stats.scaling_zero += z;
            self.stats.check_secs += tc.elapsed().as_secs_f64();
            return Ok((ins, outs));
        }
        let b = self.opts.budget.clone();
        let mut r = DecHom::pairs(&mut *self.prov, &mut self.ar, &mut self.schur, t, &b);
        self.stats.provider_calls += 1;
        if let Err(f) = &r {
            // M5i: an extra summand of a decomposable canonical ladder is a global (not local)
            // splitting, which Search's windows do not see; Exact with raised caps
            // (KRHOM_M5I_EXACT_UNK unknowns, default 1.3 M) is the fallback, its output certified
            // as any provider's (DecHom output checks, then M5's own)
            let cap: u64 = if self.prov.name() == "NoProvider" { 0 } else { std::env::var("KRHOM_M5I_EXACT_UNK").ok().and_then(|v| v.parse().ok()).unwrap_or(1_300_000) };
            let m5i = key.starts_with("kind=idem") || t.eq.is_some() || self.types.get(&p).is_some_and(|x| x.e.is_some());
            if m5i && !f.guard && cap > 0 {
                let (xw, qw) = (word_of(&self.ar, t.x), word_of(&self.ar, t.lq));
                let unk = crate::homprov::hom_unknowns(&self.lam, &qw, &xw, t.s).max(crate::homprov::hom_unknowns(&self.lam, &xw, &qw, -t.s));
                if unk <= cap as u128 {
                    eprintln!("[sindec] M5i: provider chain failed for {} -> {} ({}); Exact fallback ({} unknowns)", self.tfmt(t.q), self.tfmt(p), f.text.chars().take(120).collect::<String>(), unk);
                    let mut ex = Exact { max_label: 9, max_fibre_product: u64::MAX, max_unknowns: u64::MAX, max_call_s: 1e9, ..Exact::default() };
                    r = DecHom::pairs(&mut ex, &mut self.ar, &mut self.schur, t, &b);
                    self.stats.provider_calls += 1;
                }
            }
        }
        match r {
            Ok(o) => {
                self.stats.provider_secs += o.secs;
                // DecHom ran the output checks; repeat the scaling check at M5's own point for the counts
                let (n, z) = dechom_postcheck(&self.ar, &self.schur, t, &o.maps_in, &o.maps_out, self.opts.seed).map_err(DecErr::Check)?;
                self.stats.scaling_entries += n;
                self.stats.scaling_zero += z;
                let e = (self.tfmt(p), self.tfmt(t.q), t.s, t.mult, o.secs, o.summary.clone());
                self.log_target(e);
                let mut all = o.maps_in.clone();
                all.extend(&o.maps_out);
                self.layer_a_save(&key, &all)?;
                Ok((o.maps_in, o.maps_out))
            }
            Err(f) if f.guard => Err(DecErr::Refused(f.text)),
            Err(f) => Err(DecErr::Provider { engineering: f.engineering, text: f.text }),
        }
    }

    /// Neumann orthogonalization among the kept summands (R2 §2.2 item 5): N_ts = π_t∘ι_s
    /// (t ≠ s) is nonzero only for s_s > s_t (asserted); Π' = Σ_k (−N)^k Π.
    fn orthogonalize(&mut self, x: ExprId, sums: Vec<SSummand>) -> Result<Vec<SSummand>, DecErr> {
        let n = sums.len();
        if n <= 1 {
            return Ok(sums);
        }
        let tc = Instant::now();
        let pts = self.pts(self.opts.check_points.max(2), 0x0a0a);
        let cos = self.ar.flows(x).cosets.clone();
        // N_ts nonzero pattern by block products at the points (no nodes for zero entries)
        let mut nz = vec![vec![false; n]; n];
        for pt in &pts {
            let mut ev = Eval::new(&self.ar, pt.clone());
            for &c in &cos {
                let ib: Vec<_> = sums.iter().map(|sm| ev.block(sm.iota, c)).collect();
                let pb: Vec<_> = sums.iter().map(|sm| ev.block(sm.pi, c)).collect();
                for t in 0..n {
                    if pb[t].rows == 0 {
                        continue;
                    }
                    for s in 0..n {
                        if t != s && !nz[t][s] && ib[s].cols > 0 && !pb[t].mul(&ib[s]).is_zero() {
                            nz[t][s] = true;
                        }
                    }
                }
            }
        }
        for t in 0..n {
            for s in 0..n {
                if nz[t][s] && sums[s].s <= sums[t].s {
                    return Err(DecErr::Check(format!("nonzero π_τ∘ι_σ with s_σ = {} <= s_τ = {} (not unitriangular)", sums[s].s, sums[t].s)));
                }
            }
        }
        self.stats.check_secs += tc.elapsed().as_secs_f64();
        if !nz.iter().any(|r| r.iter().any(|&b| b)) {
            return Ok(sums);
        }
        let mark = self.ar.mark();
        let mut nm: Vec<Vec<Option<MapId>>> = vec![vec![None; n]; n];
        for t in 0..n {
            for s in 0..n {
                if nz[t][s] {
                    nm[t][s] = Some(self.ar.add(MapNode::Comp(sums[t].pi, sums[s].iota), "sindec N").map_err(DecErr::Check)?);
                }
            }
        }
        // M - I = Σ_{k>=1} (−N)^k, accumulated per entry as Lin terms
        let mut macc: Vec<Vec<Vec<(F, MapId)>>> = vec![vec![vec![]; n]; n];
        let mut pw: Vec<Vec<Option<MapId>>> = (0..n).map(|t| (0..n).map(|s| nm[t][s]).collect()).collect(); // N^1
        let mut sign = F::zero().sub(&F::one());
        let mut iter = 0;
        loop {
            let mut any = false;
            for t in 0..n {
                for s in 0..n {
                    if let Some(m) = pw[t][s] {
                        macc[t][s].push((sign.clone(), m));
                        any = true;
                    }
                }
            }
            if !any {
                break;
            }
            iter += 1;
            if iter > n {
                return Err(DecErr::Check("Π∘Ι is not unipotent".into()));
            }
            // pw <- N · pw
            let mut next: Vec<Vec<Option<MapId>>> = vec![vec![None; n]; n];
            for t in 0..n {
                for s in 0..n {
                    let mut terms = vec![];
                    for u in 0..n {
                        if let (Some(a), Some(b)) = (nm[t][u], pw[u][s]) {
                            terms.push((F::one(), self.ar.add(MapNode::Comp(a, b), "sindec N^k").map_err(DecErr::Check)?));
                        }
                    }
                    if terms.is_empty() {
                        continue;
                    }
                    let m = if terms.len() == 1 { terms[0].1 } else { self.ar.add(MapNode::Lin(terms), "sindec N^k").map_err(DecErr::Check)? };
                    let zero = pts.iter().all(|pt| Eval::new(&self.ar, pt.clone()).is_zero(m));
                    if !zero {
                        next[t][s] = Some(m);
                    }
                }
            }
            pw = next;
            sign = F::zero().sub(&sign);
        }
        let mut out = sums.clone();
        for t in 0..n {
            let mut terms: Vec<(F, MapId)> = vec![(F::one(), sums[t].pi)];
            for s in 0..n {
                if macc[t][s].is_empty() {
                    continue;
                }
                let m = if macc[t][s].len() == 1 && macc[t][s][0].0.is_one() {
                    macc[t][s][0].1
                } else {
                    self.ar.add(MapNode::Lin(macc[t][s].clone()), "sindec M").map_err(DecErr::Check)?
                };
                terms.push((F::one(), self.ar.add(MapNode::Comp(m, sums[s].pi), "sindec M pi").map_err(DecErr::Check)?));
            }
            if terms.len() > 1 {
                out[t].pi = self.ar.add(MapNode::Lin(terms), "sindec Pi'").map_err(DecErr::Check)?;
            }
        }
        // drop the scratch nodes not reachable from the new projections
        let keep: Vec<MapId> = out.iter().map(|sm| sm.pi).collect();
        let kept = self.ar.commit(mark, &keep);
        for (t, k) in kept.into_iter().enumerate() {
            out[t].pi = k;
        }
        Ok(out)
    }

    /// Π∘Ι = id (e_q on the diagonal for an idempotent type q) at check_points points on every
    /// coset of the support of X; with `full` (nothing killed), also Σ ι_t π_t = e_X (identity,
    /// or the source idempotent `ex` = Kron(e_p, R) of a decomposable L_p)
    fn check_pi_iota(&mut self, x: ExprId, sums: &[SSummand], ex: Option<MapId>, full: bool, what: &str) -> Result<(), DecErr> {
        let tc = Instant::now();
        let pts = self.pts(self.opts.check_points.max(2), 0x7e57);
        let cos = self.ar.flows(x).cosets.clone();
        let n = sums.len();
        let eqs: Vec<Option<MapId>> = sums.iter().map(|sm| self.tall.get(&sm.q).and_then(|t| t.1)).collect();
        for pt in &pts {
            let mut ev = Eval::new(&self.ar, pt.clone());
            for &c in &cos {
                let ib: Vec<_> = sums.iter().map(|sm| ev.block(sm.iota, c)).collect();
                let pb: Vec<_> = sums.iter().map(|sm| ev.block(sm.pi, c)).collect();
                for t in 0..n {
                    for s in 0..n {
                        if pb[t].rows == 0 || ib[s].cols == 0 {
                            continue;
                        }
                        let m = pb[t].mul(&ib[s]);
                        let ok = if t == s {
                            match eqs[t] {
                                None => m == crate::smap::Block::identity(m.rows as usize),
                                Some(e) => m == *ev.block(e, c),
                            }
                        } else {
                            m.is_zero()
                        };
                        if !ok {
                            return Err(DecErr::Check(format!("{}: Π∘Ι != id (summands {} {}) at coset {}", what, t, s, c)));
                        }
                    }
                }
                if full {
                    let nx = ev.block(sums[0].iota, c).rows as usize;
                    let mut acc = crate::smap::Block::<F>::zeros(nx, nx);
                    for t in 0..n {
                        if ib[t].cols == 0 {
                            continue;
                        }
                        let m = ib[t].mul(&pb[t]);
                        for (a, b) in acc.a.iter_mut().zip(&m.a) {
                            a.add_assign(b);
                        }
                    }
                    let want = match ex {
                        None => crate::smap::Block::identity(nx),
                        Some(e) => (*ev.block(e, c)).clone(),
                    };
                    if acc != want {
                        return Err(DecErr::Check(format!("{}: Σ ιπ != e_X at coset {}", what, c)));
                    }
                }
            }
        }
        self.stats.check_secs += tc.elapsed().as_secs_f64();
        Ok(())
    }

    fn dec_composite(&mut self, p: Tid, rs: &[Rung]) -> Result<Decomp, DecErr> {
        let first = self.dec_word(p, &rs[..1])?;
        let rest = rs[1..].to_vec();
        let (x, xw) = self.x_of(p, rs)?;
        let mut out: Vec<SSummand> = vec![];
        for s1 in first.summands.iter() {
            let d2 = self.dec_word(s1.q, &rest)?;
            let (i1, p1) = if first.identity {
                (None, None)
            } else {
                (Some(self.ar.add(MapNode::Kron(s1.iota, rest.clone()), "sindec kron iota").map_err(DecErr::Check)?), Some(self.ar.add(MapNode::Kron(s1.pi, rest.clone()), "sindec kron pi").map_err(DecErr::Check)?))
            };
            for s2 in d2.summands.iter() {
                let (iota, pi) = match (i1, p1, d2.identity) {
                    (None, None, _) => (s2.iota, s2.pi),
                    (Some(a), Some(b), true) => (a, b),
                    (Some(a), Some(b), false) => (self.ar.add(MapNode::Comp(a, s2.iota), "sindec iota").map_err(DecErr::Check)?, self.ar.add(MapNode::Comp(s2.pi, b), "sindec pi").map_err(DecErr::Check)?),
                    _ => unreachable!(),
                };
                if self.ar.info(iota).tgt != x || self.ar.info(pi).src != x {
                    return Err(DecErr::Check("composite: ends differ from X".into()));
                }
                out.push(SSummand { q: s2.q, s: s1.s + s2.s, iota, pi });
            }
        }
        // the kept multiset equals the kept part of term_mult (D4, release assert)
        let (kept, kill) = self.split_mult(p, rs);
        let mut got: BTreeMap<(Tid, i32), usize> = BTreeMap::new();
        for sm in &out {
            *got.entry((sm.q, sm.s)).or_insert(0) += 1;
        }
        assert_eq!(got, kept, "dec_word {} ⊗ {}: composite kept multiset != kept part of term_mult", self.tfmt(p), crate::smapio::rungs_str(rs));
        let xflows = flows(&self.lam, &xw);
        let xflows = self.src_rank(p, rs, xflows)?;
        let total: u128 = out.iter().map(|sm| self.types[&sm.q].rank).sum();
        let killed_rank = self.killed_rank(p, rs);
        assert_eq!(total + killed_rank, xflows, "dec_word {} ⊗ {}: ranks {} + killed {} != rank {}", self.tfmt(p), crate::smapio::rungs_str(rs), total, killed_rank, xflows);
        let all_killed = out.is_empty();
        if all_killed {
            self.stats.all_killed += 1;
        } else {
            let ex = self.src_idem(p, rs)?;
            self.check_pi_iota(x, &out, ex, kill.is_empty(), &format!("{} ⊗ {}", self.tfmt(p), crate::smapio::rungs_str(rs)))?;
        }
        let identity = out.len() == 1 && self.ar.info(out[0].iota).src == x && matches!(self.ar.node(out[0].iota), MapNode::Id(_));
        Ok(Decomp { p, rungs: rs.to_vec(), x, summands: out, identity, all_killed, killed: Self::killed_list(&kill) })
    }

    /// Rickard components π_τ ∘ (id_{L_p} ⊗ d_j) ∘ ι_σ of the crossing at strands (i, i+1) at the
    /// top of p, between the summands of L_p ⊗ W_ks and L_p ⊗ W_kt ((ks, kt) = (j, j+1) for a
    /// positive crossing, (j+1, j) for a negative one); nonzero ones only. Err if p or one of the
    /// two terms is killed (M6a truncates).
    pub fn rickard_comp_at(&mut self, p: Tid, i: usize, positive: bool, j: usize) -> Result<Arc<Vec<(usize, usize, MapId)>>, DecErr> {
        if let Some(r) = self.rick.get(&(p, i, positive, j)) {
            return Ok(r.clone());
        }
        let lad = self.ty(p)?.lad;
        let mu = self.schur.table(p).cols();
        if i + 1 >= mu.len() {
            return Err(DecErr::Invalid(format!("crossing position {} out of range", i)));
        }
        let (c, d) = (mu[i], mu[i + 1]);
        let terms = crossing_terms(c, d, positive, i);
        if j + 1 >= terms.len() {
            return Err(DecErr::Invalid(format!("Rickard ({},{}) has no d_{}", c, d, j)));
        }
        let (ks, kt) = if positive { (j, j + 1) } else { (j + 1, j) };
        // cell rule (R2 §2.5): the term type p_k is killed at N iff max(c, d) + k > N
        if let Some(n) = self.opts.prune_n {
            let mut sc = Schur::new(&[c, d]);
            let t0 = crossing_terms(c, d, positive, 0);
            for (k, t) in t0.iter().enumerate() {
                let ch = sc.ladder_char(&t.rungs);
                let e = sc.expand(&ch);
                assert!(e.len() == 1 && e[0].1 == Lp::one(), "Rickard term {} of ({},{}) is not indecomposable", k, c, d);
                assert_eq!(sc.killed_at(e[0].0, n as usize), c.max(d) as usize + k > n as usize, "Rickard term {} of ({},{}) at N = {}: cell rule != label rule", k, c, d, n);
            }
            for k in [ks, kt] {
                if c.max(d) as usize + k > n as usize {
                    return Err(DecErr::Killed { p: format!("Rickard term W_{} of ({},{})", k, c, d), n });
                }
            }
        }
        let src = self.dec_word(p, &terms[ks].rungs)?;
        let tgt = self.dec_word(p, &terms[kt].rungs)?;
        let (_, ids) = self.ar.lib().rickard(c, d, positive);
        let at = word_of(&self.ar, lad).len();
        let phi = self.ar.add(MapNode::Local { src: src.x, tgt: tgt.x, at, first: i, tmpl: ids[j] }, "sindec rickard d_j").map_err(DecErr::Check)?;
        let pts = self.pts(2, 0x41c4);
        let mark = self.ar.mark();
        let mut cand = vec![];
        for (si, sg) in src.summands.iter().enumerate() {
            let a = self.ar.add(MapNode::Comp(phi, sg.iota), "sindec d iota").map_err(DecErr::Check)?;
            for (ti, tu) in tgt.summands.iter().enumerate() {
                let m = self.ar.add(MapNode::Comp(tu.pi, a), "sindec rickard comp").map_err(DecErr::Check)?;
                let nonzero = pts.iter().any(|pt| !Eval::new(&self.ar, pt.clone()).is_zero(m));
                if nonzero {
                    cand.push((si, ti, m));
                }
            }
        }
        let keep: Vec<MapId> = cand.iter().map(|x| x.2).collect();
        let kept = self.ar.commit(mark, &keep);
        let out: Vec<(usize, usize, MapId)> = cand.iter().zip(kept).map(|(&(a, b, _), m)| (a, b, m)).collect();
        // degree of each component: 1 + s_σ − s_τ
        for &(si, ti, m) in &out {
            let want = 1 + src.summands[si].s - tgt.summands[ti].s;
            if self.ar.info(m).deg != want {
                return Err(DecErr::Check(format!("Rickard component degree {} != {}", self.ar.info(m).deg, want)));
            }
        }
        self.stats.rickard += 1;
        let r = Arc::new(out);
        self.rick.insert((p, i, positive, j), r.clone());
        Ok(r)
    }

    /// SPEC §5.1 form (2 strands: the crossing at position 0 of the top labels (c, d) of p)
    pub fn rickard_comp(&mut self, p: Tid, c: u8, d: u8, positive: bool, j: usize) -> Result<Arc<Vec<(usize, usize, MapId)>>, DecErr> {
        let mu = self.schur.table(p).cols();
        if mu.len() < 2 || (mu[0], mu[1]) != (c, d) {
            return Err(DecErr::Invalid(format!("rickard_comp: top labels {:?} of p do not start with ({}, {})", mu, c, d)));
        }
        self.rickard_comp_at(p, 0, positive, j)
    }

    /// closure generators (ins: R^λ → L_p, outs: L_p → R^λ) of a kept type of (λ, λ), layer A
    pub fn closure_gens(&mut self, p: Tid) -> Result<Arc<ClosureGens>, DecErr> {
        if let Some(c) = self.clos.get(&p) {
            return Ok(c.clone());
        }
        let (lad, ep) = {
            let t = self.ty(p)?;
            (t.lad, t.e)
        };
        if self.schur.table(p).cols() != self.lam {
            return Err(DecErr::Invalid(format!("closure_gens: {} is not a type of (λ, λ)", self.tfmt(p))));
        }
        let r0 = self.ar.expr(Expr::new(&self.lam, &[]).map_err(DecErr::Invalid)?);
        let cd = self.schur.c_diag(p);
        let mut want = vec![];
        for (i, &c) in cd.c.iter().enumerate() {
            for _ in 0..c {
                want.push(cd.lo + i as i32);
            }
        }
        let mut res = [vec![], vec![]];
        for (k, outs) in [false, true].into_iter().enumerate() {
            let key = self.closure_key(p, outs);
            let ms = match self.layer_a_load(&key, want.len())? {
                // M5i: a loaded generator carries its own (imported) copy of e_p; re-wrap with this
                // builder's e_p so that certify_closure sees Comp(e_p, g) / Comp(g, e_p) (it then
                // re-checks the family numerically: rank n at the identity coset at 2 points)
                Some(v) => match ep {
                    None => v,
                    Some(e) => v.into_iter().map(|m| self.ar.add(if outs { MapNode::Comp(m, e) } else { MapNode::Comp(e, m) }, "m5i rewrap e_p")).collect::<Result<Vec<_>, _>>().map_err(DecErr::Check)?,
                },
                None => {
                    let b = self.opts.budget.clone();
                    let mut r = DecHom::closure(&mut *self.prov, &mut self.ar, &mut self.schur, lad, p, ep, outs, &b);
                    self.stats.provider_calls += 1;
                    if let (Err(f), Some(_)) = (&r, ep) {
                        // M5i: closure through e_p, Exact fallback with raised caps (as for pairs)
                        let cap: u64 = if self.prov.name() == "NoProvider" { 0 } else { std::env::var("KRHOM_M5I_EXACT_UNK").ok().and_then(|v| v.parse().ok()).unwrap_or(1_300_000) };
                        let w = word_of(&self.ar, lad);
                        let dmax = want.iter().copied().max().unwrap_or(0);
                        let unk = (0..=dmax.max(0)).map(|d| crate::homprov::hom_unknowns(&self.lam, &[], &w, d).max(crate::homprov::hom_unknowns(&self.lam, &w, &[], d))).max().unwrap_or(0);
                        if !f.guard && cap > 0 && unk <= cap as u128 {
                            eprintln!("[sindec] M5i: closure of {} ({}) failed through the chain ({}); Exact fallback (≤ {} unknowns per degree)", self.tfmt(p), if outs { "outs" } else { "ins" }, f.text.chars().take(120).collect::<String>(), unk);
                            let mut ex = Exact { max_label: 9, max_fibre_product: u64::MAX, max_unknowns: u64::MAX, max_call_s: 1e9, ..Exact::default() };
                            r = DecHom::closure(&mut ex, &mut self.ar, &mut self.schur, lad, p, ep, outs, &b);
                            self.stats.provider_calls += 1;
                        }
                    }
                    match r {
                        Ok(o) => {
                            self.stats.provider_secs += o.secs;
                            let e = (self.tfmt(p), format!("closure-{}", if outs { "outs" } else { "ins" }), 0, want.len(), o.secs, o.summary.clone());
                            self.log_target(e);
                            let ms = if outs { o.maps_out } else { o.maps_in };
                            self.layer_a_save(&key, &ms)?;
                            ms
                        }
                        Err(f) if f.guard => return Err(DecErr::Refused(f.text)),
                        Err(f) => return Err(DecErr::Provider { engineering: f.engineering, text: f.text }),
                    }
                }
            };
            certify_closure(&self.ar, &self.lam, r0, lad, ep, &ms, &want, outs, self.opts.seed ^ 0xc105e).map_err(DecErr::Check)?;
            res[k] = ms;
        }
        self.stats.closure += 1;
        let [ins, outs] = res;
        let g = Arc::new(ClosureGens { p, ins, outs });
        self.clos.insert(p, g.clone());
        Ok(g)
    }

    fn log_target(&mut self, e: (String, String, i32, usize, f64, String)) {
        if let Some(f) = &self.target_log_file {
            use std::io::Write;
            let esc = |x: &str| x.replace('\\', "\\\\").replace('"', "\\\"");
            let line = format!("{{\"p\":\"{}\",\"q\":\"{}\",\"s\":{},\"mult\":{},\"secs\":{:.3},\"summary\":\"{}\"}}\n", esc(&e.0), esc(&e.1), e.2, e.3, e.4, esc(&e.5));
            if let Ok(mut fh) = std::fs::OpenOptions::new().create(true).append(true).open(f) {
                let _ = fh.write_all(line.as_bytes());
            }
        }
        self.target_log.push(e);
    }

    /// visit every kept-source pair of a plan (dec_rung), in order; returns the failures
    pub fn prebuild(&mut self, plan: &Plan) -> Vec<(usize, DecErr)> {
        let mut fails = vec![];
        for (k, pp) in plan.pairs.iter().enumerate() {
            if pp.p_killed {
                continue;
            }
            if let Err(e) = self.dec_rung(pp.p, pp.rung) {
                fails.push((k, e));
            }
        }
        fails
    }

    // ---------------------------------------------------------- layer A

    fn layer_a_dir(&self) -> Option<PathBuf> {
        self.opts.cache_dir.as_ref().map(|d| d.join("layer_a"))
    }
    /// the stamps of a layer-A header (format, ladder function, TMPL_ALG, provider)
    pub fn layer_a_stamps(&self) -> String {
        format!("{} ladfn={} tmpl_alg={} prov={}", SINDEC_VERSION, LADDER_STAMP, TMPL_ALG, self.prov_stamp)
    }
    fn pairs_key(&mut self, p: Tid, rs: &[Rung], t: &PairTarget) -> String {
        let lp = canonical_ladder(&mut self.schur, p);
        let lq = canonical_ladder(&mut self.schur, t.q);
        format!(
            "kind=pairs;lam={};field={};p={};Lp={};R={};q={};Lq={};s={};mult={};eq={}",
            crate::smapio::labels_str(&self.lam),
            F::name(),
            self.tfmt(p),
            crate::smapio::rungs_str(&lp),
            crate::smapio::rungs_str(rs),
            self.tfmt(t.q),
            crate::smapio::rungs_str(&lq),
            t.s,
            t.mult,
            if t.eq.is_some() { M5I_STAMP } else { "none" }
        )
    }
    /// layer-A key of an extra summand of a decomposable canonical ladder (M5i): pairs
    /// L_q ⇄ X = L_p
    fn idem_key(&mut self, p: Tid, t: &PairTarget) -> String {
        let lp = canonical_ladder(&mut self.schur, p);
        let lq = canonical_ladder(&mut self.schur, t.q);
        format!(
            "kind=idem;lam={};field={};p={};Lp={};q={};Lq={};s={};mult={};eq={}",
            crate::smapio::labels_str(&self.lam),
            F::name(),
            self.tfmt(p),
            crate::smapio::rungs_str(&lp),
            self.tfmt(t.q),
            crate::smapio::rungs_str(&lq),
            t.s,
            t.mult,
            if t.eq.is_some() { M5I_STAMP } else { "none" }
        )
    }
    fn closure_key(&mut self, p: Tid, outs: bool) -> String {
        let lp = canonical_ladder(&mut self.schur, p);
        let ep = if self.types.get(&p).is_some_and(|t| t.e.is_some()) { M5I_STAMP } else { "none" };
        format!("kind=closure;lam={};field={};p={};Lp={};outs={};ep={}", crate::smapio::labels_str(&self.lam), F::name(), self.tfmt(p), crate::smapio::rungs_str(&lp), outs, ep)
    }
    /// path of the layer-A file of a key (named by the hash of the key only; stamps are in the
    /// header so that a stale file is found and rejected)
    pub fn layer_a_path(&self, key: &str) -> Option<PathBuf> {
        self.layer_a_dir().map(|d| d.join(format!("a-{:016x}.txt", fnv64(key))))
    }
    fn layer_a_header(&self, key: &str) -> String {
        format!("{} key={}", self.layer_a_stamps(), key)
    }
    /// Ok(None) if absent or stale (rejected, counted); Ok(Some(maps)) if loaded
    fn layer_a_load(&mut self, key: &str, want: usize) -> Result<Option<Vec<MapId>>, DecErr> {
        let Some(path) = self.layer_a_path(key) else { return Ok(None) };
        let Ok(text) = std::fs::read_to_string(&path) else { return Ok(None) };
        let (h, body) = text.split_once('\n').unwrap_or((text.as_str(), ""));
        let want_h = self.layer_a_header(key);
        let by_route = h != want_h;
        if by_route && !(self.opts.layer_a_route_stamp && route_header(h) == route_header(&want_h)) {
            self.stats.layer_a_rejected += 1;
            return Ok(None);
        }
        match crate::smapio::import(&mut self.ar, body) {
            Ok(v) if v.len() == want => {
                self.stats.layer_a_loaded += 1;
                self.stats.layer_a_route += by_route as usize;
                Ok(Some(v))
            }
            Ok(v) => {
                self.stats.layer_a_rejected += 1;
                eprintln!("[sindec] layer A {}: {} maps, expected {} (rejected)", path.display(), v.len(), want);
                Ok(None)
            }
            Err(e) => {
                self.stats.layer_a_rejected += 1;
                eprintln!("[sindec] layer A {}: {} (rejected)", path.display(), e);
                Ok(None)
            }
        }
    }
    fn layer_a_save(&mut self, key: &str, maps: &[MapId]) -> Result<(), DecErr> {
        let Some(path) = self.layer_a_path(key) else { return Ok(()) };
        let body = crate::smapio::export(&self.ar, maps).map_err(DecErr::Cache)?;
        let tmp = path.with_extension(format!("{}.tmp", std::process::id()));
        std::fs::write(&tmp, format!("{}\n{}", self.layer_a_header(key), body)).map_err(|e| DecErr::Cache(e.to_string()))?;
        std::fs::rename(&tmp, &path).map_err(|e| DecErr::Cache(e.to_string()))?;
        // the templates the export references by call must be on disk for a later import
        let _ = self.ar.lib_ref().save();
        self.stats.layer_a_saved += 1;
        Ok(())
    }
    /// bytes of the layer-A files on disk
    pub fn layer_a_bytes(&self) -> u64 {
        let Some(d) = self.layer_a_dir() else { return 0 };
        std::fs::read_dir(d).map(|it| it.filter_map(|e| e.ok()).filter_map(|e| e.metadata().ok()).map(|m| m.len()).sum()).unwrap_or(0)
    }

    /// immutable snapshot for parallel scans
    pub fn freeze(self) -> Arc<DecData<F>> {
        let mut level = HashMap::new();
        let mut names = HashMap::new();
        let mut ids: Vec<Tid> = self.types.keys().copied().collect();
        for d in self.decs.values() {
            ids.push(d.p);
            ids.extend(d.summands.iter().map(|s| s.q));
            ids.extend(d.killed.iter().map(|k| k.0));
        }
        for q in ids {
            level.insert(q, kill_level(&self.schur, q));
            names.insert(q, self.schur.table(q).fmt());
        }
        Arc::new(DecData { ar: Arc::new(self.ar), lam: self.lam, prune_n: self.opts.prune_n, types: self.types.into_iter().map(|(k, v)| (k, Arc::new(v))).collect(), decs: self.decs, rick: self.rick, clos: self.clos, level, names })
    }
}

// ------------------------------------------------------------------ frozen data

/// Frozen data: read-only lookups that return Err(NotPrecomputed) instead of computing.
/// Send + Sync (no interior mutability).
pub struct DecData<F: Field> {
    ar: Arc<Arena<F>>,
    lam: Vec<u8>,
    prune_n: Option<u8>,
    types: HashMap<Tid, Arc<SType>>,
    decs: HashMap<(Tid, Vec<(usize, u8, usize)>), Arc<Decomp>>,
    rick: HashMap<(Tid, usize, bool, usize), Arc<Vec<(usize, usize, MapId)>>>,
    clos: HashMap<Tid, Arc<ClosureGens>>,
    level: HashMap<Tid, usize>,
    names: HashMap<Tid, String>,
}

impl<F: Field> DecData<F> {
    pub fn ar(&self) -> &Arena<F> {
        &self.ar
    }
    pub fn lam(&self) -> &[u8] {
        &self.lam
    }
    pub fn prune_n(&self) -> Option<u8> {
        self.prune_n
    }
    fn name(&self, p: Tid) -> String {
        self.names.get(&p).cloned().unwrap_or_else(|| format!("tid {}", p))
    }
    pub fn ty(&self, p: Tid) -> Result<&SType, DecErr> {
        self.types.get(&p).map(|t| t.as_ref()).ok_or_else(|| DecErr::NotPrecomputed(format!("type {}", self.name(p))))
    }
    pub fn dec_word(&self, p: Tid, rs: &[Rung]) -> Result<Arc<Decomp>, DecErr> {
        self.decs.get(&(p, rkey(rs))).cloned().ok_or_else(|| DecErr::NotPrecomputed(format!("decomposition of {} ⊗ {}", self.name(p), crate::smapio::rungs_str(rs))))
    }
    pub fn rickard_comp_at(&self, p: Tid, i: usize, positive: bool, j: usize) -> Result<Arc<Vec<(usize, usize, MapId)>>, DecErr> {
        self.rick.get(&(p, i, positive, j)).cloned().ok_or_else(|| DecErr::NotPrecomputed(format!("Rickard components of {} at {} ({}) d_{}", self.name(p), i, positive, j)))
    }
    pub fn rickard_comp(&self, p: Tid, _c: u8, _d: u8, positive: bool, j: usize) -> Result<Arc<Vec<(usize, usize, MapId)>>, DecErr> {
        self.rickard_comp_at(p, 0, positive, j)
    }
    pub fn closure_gens(&self, p: Tid) -> Result<Arc<ClosureGens>, DecErr> {
        self.clos.get(&p).cloned().ok_or_else(|| DecErr::NotPrecomputed(format!("closure generators of {}", self.name(p))))
    }
    /// reduced coordinates of B_p at the coset c: identity for e_p = None; for an idempotent
    /// type (M5i) P (flows × m), Q (m × flows) with QP = id_m and PQ = block(e_p) (smap::reduced_of,
    /// pivot columns of e_p: the basis is point dependent, every engine quantity is invariant)
    pub fn reduced(&self, ev: &mut Eval<F>, p: Tid, c: u64) -> Result<Reduced<F>, DecErr> {
        let t = self.ty(p)?;
        match t.e {
            None => {
                let n = self.ar.flows(t.lad).range(c).len();
                Ok(Reduced { p: crate::smap::Block::identity(n), q: crate::smap::Block::identity(n), m: n })
            }
            Some(e) => crate::smap::reduced_of(&ev.block(e, c)).map_err(|m| DecErr::Check(format!("reduced coordinates of {} at coset {}: {}", self.name(p), c, m))),
        }
    }
    /// number of precomputed decompositions and Rickard lists
    pub fn sizes(&self) -> (usize, usize, usize, usize) {
        (self.types.len(), self.decs.len(), self.rick.len(), self.clos.len())
    }
    /// all precomputed (p, rungs) keys
    pub fn dec_keys(&self) -> Vec<(Tid, Vec<Rung>)> {
        let mut v: Vec<(Tid, Vec<Rung>)> = self.decs.values().map(|d| (d.p, d.rungs.clone())).collect();
        v.sort_by(|a, b| a.0.cmp(&b.0).then(rkey(&a.1).cmp(&rkey(&b.1))));
        v
    }
    fn killed_at(&self, q: Tid, n: usize) -> bool {
        self.level.get(&q).is_some_and(|&l| n < l)
    }
    /// The data for gl_N' with N' = n ≤ prune_n (R2 §4.4): the same Arena, types killed at n
    /// removed, summands killed at n dropped from every decomposition (Π∘Ι = id holds for the
    /// subfamily), Rickard lists filtered and re-indexed. Err for n > prune_n.
    pub fn restrict(&self, n: u8) -> Result<DecData<F>, DecErr> {
        if let Some(n0) = self.prune_n {
            if n > n0 {
                return Err(DecErr::Invalid(format!("restrict({}) of data pruned at N = {}", n, n0)));
            }
        }
        let nn = n as usize;
        let types: HashMap<Tid, Arc<SType>> = self.types.iter().filter(|(q, _)| !self.killed_at(**q, nn)).map(|(k, v)| (*k, v.clone())).collect();
        let mut decs = HashMap::new();
        let mut remap: HashMap<(Tid, Vec<(usize, u8, usize)>), Vec<Option<usize>>> = HashMap::new();
        for (k, d) in &self.decs {
            if self.killed_at(d.p, nn) {
                continue;
            }
            let mut idx = vec![];
            let mut sums = vec![];
            let mut killed = d.killed.clone();
            for sm in &d.summands {
                if self.killed_at(sm.q, nn) {
                    idx.push(None);
                    killed.push((sm.q, sm.s, 1));
                } else {
                    idx.push(Some(sums.len()));
                    sums.push(sm.clone());
                }
            }
            let all_killed = sums.is_empty();
            decs.insert(k.clone(), Arc::new(Decomp { p: d.p, rungs: d.rungs.clone(), x: d.x, summands: sums, identity: d.identity && !all_killed, all_killed, killed }));
            remap.insert(k.clone(), idx);
        }
        let mut rick = HashMap::new();
        for (&(p, i, positive, j), comps) in &self.rick {
            if self.killed_at(p, nn) {
                continue;
            }
            let mu: Vec<u8> = {
                let lad = self.types[&p].lad;
                Expr::new(&self.lam, &self.ar.expr_of(lad).rungs).map(|e| e.top()).unwrap_or_default()
            };
            let (c, d) = (mu[i], mu[i + 1]);
            let (ks, kt) = if positive { (j, j + 1) } else { (j + 1, j) };
            if c.max(d) as usize + ks.max(kt) > nn {
                continue; // a term killed at n: M6a truncates
            }
            let terms = crossing_terms(c, d, positive, i);
            let ms = &remap[&(p, rkey(&terms[ks].rungs))];
            let mt = &remap[&(p, rkey(&terms[kt].rungs))];
            let v: Vec<(usize, usize, MapId)> = comps.iter().filter_map(|&(a, b, m)| Some((ms[a]?, mt[b]?, m))).collect();
            rick.insert((p, i, positive, j), Arc::new(v));
        }
        let clos = self.clos.iter().filter(|(q, _)| !self.killed_at(**q, nn)).map(|(k, v)| (*k, v.clone())).collect();
        Ok(DecData { ar: self.ar.clone(), lam: self.lam.clone(), prune_n: Some(n), types, decs, rick, clos, level: self.level.clone(), names: self.names.clone() })
    }
}

// ------------------------------------------------------------------ planning (schur only)

#[derive(Clone, Debug)]
pub struct PlanTarget {
    pub q: Tid,
    pub s: i32,
    pub mult: usize,
    /// max(kill_level(p), kill_level(q)): the smallest N at which the target is needed
    pub min_level: usize,
}

#[derive(Clone, Debug)]
pub struct PlanPair {
    pub mu: Vec<u8>,
    pub p: Tid,
    pub rung: Rung,
    pub p_killed: bool,
    pub kept: Vec<PlanTarget>,
    pub killed: Vec<PlanTarget>,
    pub all_killed: bool,
    pub identity: bool,
}

#[derive(Clone, Debug)]
pub struct Plan {
    pub lam: Vec<u8>,
    pub n: Option<u8>,
    pub pairs: Vec<PlanPair>,
}

impl Plan {
    /// (pairs, pairs with a kept source, non-identity kept-source pairs with a provider call,
    /// all-killed pairs, kept targets, Σ mult of kept targets, killed targets)
    pub fn counts(&self) -> (usize, usize, usize, usize, usize, usize, usize) {
        let kp: Vec<&PlanPair> = self.pairs.iter().filter(|p| !p.p_killed).collect();
        let calls = kp.iter().filter(|p| !p.identity && !p.all_killed).count();
        let ak = kp.iter().filter(|p| p.all_killed).count();
        let real: Vec<&&PlanPair> = kp.iter().filter(|p| !p.identity && !p.all_killed).collect();
        let kt: usize = real.iter().map(|p| p.kept.len()).sum();
        let sm: usize = real.iter().map(|p| p.kept.iter().map(|t| t.mult).sum::<usize>()).sum();
        let kd: usize = kp.iter().map(|p| p.killed.len()).sum();
        (self.pairs.len(), kp.len(), calls, ak, kt, sm, kd)
    }
    pub fn summary(&self) -> String {
        let (a, b, c, d, e, f, g) = self.counts();
        format!(
            "plan lam={} N={}: {} pairs, {} with kept source, {} non-identity with provider targets, {} all-killed; kept targets {} (Σ mult {}), killed targets {}",
            crate::smapio::labels_str(&self.lam),
            self.n.map_or("none".to_string(), |n| n.to_string()),
            a,
            b,
            c,
            d,
            e,
            f,
            g
        )
    }
}

/// The G1 work list of λ (every (block type, rung) pair of one crossing term, both signs, every
/// position) with kept and killed targets at N (R2 §6). Schur only.
pub fn pruned_plan(s: &mut Schur, n: Option<u8>) -> Plan {
    let lam = s.lam.clone();
    let mut pairs = vec![];
    for (mu, r) in crate::homprov::g1_work(&lam) {
        let order = s.canonical_basis(&mu).order.clone();
        for p in order {
            let lvp = kill_level(s, p);
            let p_killed = n.is_some_and(|nn| s.killed_at(p, nn as usize));
            let mut lp = canonical_ladder(s, p);
            let mult = s.term_mult(p, &[r]);
            let (mut kept, mut killed) = (vec![], vec![]);
            for (q, m) in &mult {
                let kq = n.is_some_and(|nn| s.killed_at(*q, nn as usize));
                let lvq = kill_level(s, *q);
                for (i, &c) in m.c.iter().enumerate() {
                    if c != 0 {
                        let t = PlanTarget { q: *q, s: m.lo + i as i32, mult: c as usize, min_level: lvp.max(lvq) };
                        if kq { killed.push(t) } else { kept.push(t) }
                    }
                }
            }
            lp.push(r);
            let identity = mult.len() == 1 && mult[0].1 == Lp::one() && !killed.iter().any(|_| true) && canonical_ladder(s, mult[0].0) == lp;
            let all_killed = kept.is_empty();
            pairs.push(PlanPair { mu: mu.clone(), p, rung: r, p_killed, kept, killed, all_killed, identity });
        }
    }
    Plan { lam, n, pairs }
}

/// One decomposable canonical ladder over the visited blocks of λ.
#[derive(Clone, Debug)]
pub struct IdemCase {
    pub mu: Vec<u8>,
    pub q: Tid,
    pub table: String,
    pub ladder: String,
    pub flows: u128,
    /// char L_q expanded: (table, multiplicity)
    pub char_expansion: Vec<(String, String)>,
    pub kill_level: usize,
}

/// The visited blocks of λ: λ and the blocks of the G1 work list (both crossing signs).
pub fn visited_blocks(lam: &[u8]) -> Vec<Vec<u8>> {
    let mut bl: Vec<Vec<u8>> = vec![lam.to_vec()];
    for (mu, _) in crate::homprov::g1_work(lam) {
        if !bl.contains(&mu) {
            bl.push(mu);
        }
    }
    bl
}

/// Every type of the visited blocks of λ whose canonical ladder is decomposable (char L_q ≠ b_q),
/// with its kill level; `idem_needed(s, n)` keeps those kept at n (the idempotents M5i must
/// build before level n).
pub fn decomposable_ladders(s: &mut Schur) -> Vec<IdemCase> {
    let lam = s.lam.clone();
    let mut out = vec![];
    for mu in visited_blocks(&lam) {
        let order = s.canonical_basis(&mu).order.clone();
        for q in order {
            let w = canonical_ladder(s, q);
            if is_indecomposable_of(s, &w, q) {
                continue;
            }
            let ch = s.ladder_char(&w);
            let ex: Vec<(String, String)> = s.expand(&ch).iter().map(|(t, c)| (s.table(*t).fmt(), c.fmt_ec())).collect();
            out.push(IdemCase { mu: mu.clone(), q, table: s.table(q).fmt(), ladder: crate::smapio::rungs_str(&w), flows: flows(&lam, &w), char_expansion: ex, kill_level: kill_level(s, q) });
        }
    }
    out
}

pub fn idem_needed(s: &mut Schur, n: Option<u8>) -> Vec<IdemCase> {
    decomposable_ladders(s).into_iter().filter(|c| n.is_none_or(|nn| c.kill_level <= nn as usize)).collect()
}

// ------------------------------------------------------------------ knot coverage (schur only)

/// Library jobs a pruned localized scan of a braid word needs (M5c): the kept-source pairs of
/// `pruned_plan` visited rung by rung and the closure types of the final chain objects, as an
/// upper bound (no cancellation) and a lower bound (types with a nonzero K_0 class).
#[derive(Clone, Debug, Default)]
pub struct KnotJobs {
    /// upper bound: indices into `plan.pairs`
    pub pairs: std::collections::BTreeSet<usize>,
    /// upper bound: non-identity kept types of (λ, λ) present after the last crossing
    pub closure: std::collections::BTreeSet<Tid>,
    /// lower bound: pairs visited from types whose signed K_0 class is nonzero
    pub pairs_lb: std::collections::BTreeSet<usize>,
    /// lower bound: non-identity final types with a nonzero K_0 class
    pub closure_lb: std::collections::BTreeSet<Tid>,
    /// (upper, lower) number of kept types after each crossing
    pub stage_types: Vec<(usize, usize)>,
}

/// Jobs a pruned scan (`plan.n`) of the braid word (generators ±i, 1-based, positive =
/// positive crossing) visits, starting from the identity type of (λ, λ): for every type p present
/// before a crossing at strands (i, i+1) and every Rickard term W_k that is not killed by the cell
/// rule (max(c, d) + k ≤ N), the rungs of W_k are applied one at a time; each (type, rung) step
/// is a plan pair, its kept summand types (`term_mult`, killed types dropped) feed the next rung.
///
/// Upper bound: types are united over all homological degrees (no cancellation). Lower bound:
/// only types whose signed class Σ_h (−1)^h [object] in K_0 (with q-shifts, as `decat_scan`) is
/// nonzero, which every complex homotopy equivalent to the stage contains (Gaussian elimination
/// cannot remove them); the steps of a crossing from those types are needed whatever the scan
/// cancels. Err if a visited step is not a plan pair (e.g. unequal labels: `g1_work` only has the
/// blocks reached by one crossing from λ).
pub fn knot_jobs(s: &mut Schur, plan: &Plan, word: &[i32]) -> Result<KnotJobs, String> {
    use std::collections::BTreeSet;
    let lam = plan.lam.clone();
    let idx: HashMap<(Tid, Rung), usize> = plan.pairs.iter().enumerate().map(|(k, pp)| ((pp.p, pp.rung), k)).collect();
    let killed = |s: &Schur, q: Tid| plan.n.is_some_and(|n| s.killed_at(q, n as usize));
    let id = s.intern(crate::schur::Table::identity(&lam));
    let mut cur: BTreeSet<Tid> = [id].into_iter().collect();
    let mut cls: BTreeMap<Tid, Lp> = [(id, Lp::one())].into_iter().collect();
    let mut out = KnotJobs::default();
    for &g in word {
        if g == 0 || g.unsigned_abs() as usize >= lam.len() {
            return Err(format!("generator {} out of range for {} strands", g, lam.len()));
        }
        let i = g.unsigned_abs() as usize - 1;
        let positive = g > 0;
        let mut next: BTreeSet<Tid> = BTreeSet::new();
        let mut ncls: BTreeMap<Tid, Lp> = BTreeMap::new();
        for &p in &cur {
            let cp = cls.get(&p).cloned().unwrap_or_else(Lp::zero);
            let mu = s.table(p).cols();
            let (c, d) = (mu[i], mu[i + 1]);
            for (k, t) in crossing_terms(c, d, positive, i).iter().enumerate() {
                if plan.n.is_some_and(|n| c.max(d) as usize + k > n as usize) {
                    continue;
                }
                let mut layer: BTreeSet<Tid> = [p].into_iter().collect();
                for r in &t.rungs {
                    let mut nl = BTreeSet::new();
                    for &q in &layer {
                        let k = *idx.get(&(q, *r)).ok_or_else(|| format!("step {} ⊗ {} is not a plan pair", s.table(q).fmt(), crate::smapio::rung_str(r)))?;
                        out.pairs.insert(k);
                        if !cp.is_zero() {
                            out.pairs_lb.insert(k);
                        }
                        for (q2, _) in s.term_mult(q, &[*r]) {
                            if !killed(s, q2) {
                                nl.insert(q2);
                            }
                        }
                    }
                    layer = nl;
                }
                next.extend(layer);
                if !cp.is_zero() {
                    let sign = if t.hdeg.rem_euclid(2) == 0 { 1 } else { -1 };
                    let tm = if t.rungs.is_empty() { vec![(p, Lp::one())] } else { s.term_mult(p, &t.rungs) };
                    for (q2, m) in tm {
                        if !killed(s, q2) {
                            let e = ncls.entry(q2).or_insert_with(Lp::zero);
                            e.add_scaled(&cp.mul(&m), sign, t.dq);
                        }
                    }
                }
            }
        }
        ncls.retain(|_, v| !v.is_zero());
        out.stage_types.push((next.len(), ncls.len()));
        cur = next;
        cls = ncls;
    }
    for q in cur {
        if !canonical_ladder(s, q).is_empty() {
            out.closure.insert(q);
            if cls.contains_key(&q) {
                out.closure_lb.insert(q);
            }
        }
    }
    Ok(out)
}
