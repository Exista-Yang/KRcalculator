//! M6a-1a: the core of the localized multi-strand scanning engine (SPEC_M4_M6a §6.3-6.4,
//! SPEC_3C_GL6 rev1 §2 "M6a-1a" and D14), without the decomposition library.
//!
//! * [`SCx`]: a complex of localized objects at one evaluation point. Objects are (type, h, t);
//!   every component is a [`LocMap`] over an `Arc<Layout>` shared per (src type, tgt type)
//!   (a type determines its top block). Blocks are in reduced coordinates (= flow coordinates
//!   while every e_p = id, which holds for every kept type at N ≤ 6, M5a report).
//! * [`TypeShape`] / [`ShapeReg`]: per type the support cosets, the block sizes m_c and the
//!   top-stratum coset rep(p) (where m_c = 1); layouts are built from shapes and cached.
//! * [`PEngine`]: one scan at one point. Crossings are inserted step by step (one step = one
//!   Rickard term, in the order in which the term differential runs: ascending terms for a
//!   positive crossing, descending for a negative one, so the id ⊗ d_W arrows always run from
//!   step s to step s+1). Gaussian elimination with scalar pivots in Markowitz order kept in a
//!   lazy, versioned min-heap ([`PivotOrder::Heap`]; [`PivotOrder::Rescan`] is the quadratic
//!   reference order of xengine and must give the identical pivot sequence). With streaming on,
//!   a pivot a → b is eligible only once out(a) is complete, i.e. step(a) < the current step
//!   (in(b) is always complete); after the last step every pivot is eligible (SPEC §6.4).
//!   Arrows from objects that died as pivot targets (and into objects that died as pivot
//!   sources) are discarded; any other arrow touching a dead object is an error
//!   ([`SErr::Stream`]) — this is the run-time form of the eligibility rule.
//! * Frozen structure (D14): [`Mode::Record`] (point 0) decides the pivots from the values and
//!   records a [`FrozenTrace`]: the sparsity-pattern hash after every insertion step, every
//!   elimination round's pivot list and resulting pattern hash, and a determinism hash per
//!   crossing. [`Mode::Replay`] (every other point) never decides structure: it applies the
//!   recorded pivots, asserts that each one is still a nonzero scalar, and compares every
//!   pattern hash; any difference is [`SErr::PointDependent`] (detected, never ignored).
//! * Memory accounting: an incremental byte counter of the complex (LocMap data plus fixed
//!   per-object and per-edge overheads), the previous complex while a crossing is built, the
//!   layouts, and external bytes declared by the caller (the Eval memo). It is checked after
//!   every arrow insertion against `mem_point_bytes`: [`SErr::Memory`] refuses before the
//!   budget is exceeded by more than one LocMap.
//! * [`PointScan`] (the per-point output for the box driver) and [`gram_at_point`] /
//!   [`GramAtPoint`] (G[b][a] = h_b ∘ g_a at the identity coset).
//! * [`SLib`]: the read-only library interface of SPEC_M4_M6a §5.1 (frozen `DecData`), which
//!   M6a-1b wires into the crossing step; implemented here for `sindec::DecData`.
//! * [`synth`]: synthetic localized complexes with known minimal models (tests, measurements).
use crate::field::Field;
use crate::schur::{Schur, Tid};
use crate::sindec::{ClosureGens, DecData, DecErr, Decomp, SType};
use crate::smap::{gemm_acc, Arena, Eval, Layout, LocMap, MapId, Reduced};
use crate::web::Rung;
use std::cmp::Reverse;
use std::collections::{BTreeMap, BTreeSet, BinaryHeap, HashMap};
use std::sync::Arc;
use std::time::Instant;

// ------------------------------------------------------------------ errors

/// Errors of the scanning engine. `Memory` and `PointDependent` are the two refusals of the
/// spec; the others are failed assertions.
#[derive(Clone, Debug, PartialEq)]
pub enum SErr {
    /// per-point memory above the budget (SPEC §6.3 step 4)
    Memory { bytes: u64, budget: u64, at: String },
    /// a replayed point needs a different structure than point 0 (D14)
    PointDependent { stage: usize, round: usize, what: String },
    /// a pivot candidate is not a nonzero scalar c·I at every coset
    Pivot(String),
    /// the streaming invariant was violated (arrow at an object eliminated too early)
    Stream(String),
    /// unknown type / mismatched shape or layout
    Shape(String),
    /// d∘d ≠ 0
    D2(String),
    /// library error (M6a-1b)
    Lib(String),
    /// API misuse (steps not opened, ...)
    Usage(String),
}
impl std::fmt::Display for SErr {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SErr::Memory { bytes, budget, at } => write!(
                f,
                "per-point memory {:.1} MB above the budget {:.1} MB at {}: run with streaming or on the cloud",
                *bytes as f64 / 1e6,
                *budget as f64 / 1e6,
                at
            ),
            SErr::PointDependent { stage, round, what } => write!(f, "point-dependent scan at stage {} round {}: {}", stage, round, what),
            SErr::Pivot(s) => write!(f, "pivot assertion: {}", s),
            SErr::Stream(s) => write!(f, "streaming invariant: {}", s),
            SErr::Shape(s) => write!(f, "shape: {}", s),
            SErr::D2(s) => write!(f, "d∘d ≠ 0: {}", s),
            SErr::Lib(s) => write!(f, "library: {}", s),
            SErr::Usage(s) => write!(f, "usage: {}", s),
        }
    }
}
impl From<SErr> for String {
    fn from(e: SErr) -> String {
        e.to_string()
    }
}
impl From<DecErr> for SErr {
    fn from(e: DecErr) -> SErr {
        SErr::Lib(e.to_string())
    }
}

// ------------------------------------------------------------------ the library interface

/// What the engine reads from the frozen library: exactly the read-only part of SPEC_M4_M6a
/// §5.1 (`DecData`), with M5a's 3-strand `rickard_comp_at`. Send + Sync so that points can be
/// scanned in parallel over one `Arc`.
pub trait SLib<F: Field>: Send + Sync {
    fn ar(&self) -> &Arena<F>;
    fn lam(&self) -> &[u8];
    fn prune_n(&self) -> Option<u8>;
    fn ty(&self, p: Tid) -> Result<&SType, DecErr>;
    fn dec_word(&self, p: Tid, rs: &[Rung]) -> Result<Arc<Decomp>, DecErr>;
    fn rickard_comp_at(&self, p: Tid, i: usize, positive: bool, j: usize) -> Result<Arc<Vec<(usize, usize, MapId)>>, DecErr>;
    fn closure_gens(&self, p: Tid) -> Result<Arc<ClosureGens>, DecErr>;
    /// reduced coordinates of B_p at the coset c (identity while e_p = None); callers memoize per point
    fn reduced(&self, ev: &mut Eval<F>, p: Tid, c: u64) -> Result<Reduced<F>, DecErr>;
}

impl<F: Field> SLib<F> for DecData<F> {
    fn ar(&self) -> &Arena<F> {
        DecData::ar(self)
    }
    fn lam(&self) -> &[u8] {
        DecData::lam(self)
    }
    fn prune_n(&self) -> Option<u8> {
        DecData::prune_n(self)
    }
    fn ty(&self, p: Tid) -> Result<&SType, DecErr> {
        DecData::ty(self, p)
    }
    fn dec_word(&self, p: Tid, rs: &[Rung]) -> Result<Arc<Decomp>, DecErr> {
        DecData::dec_word(self, p, rs)
    }
    fn rickard_comp_at(&self, p: Tid, i: usize, positive: bool, j: usize) -> Result<Arc<Vec<(usize, usize, MapId)>>, DecErr> {
        DecData::rickard_comp_at(self, p, i, positive, j)
    }
    fn closure_gens(&self, p: Tid) -> Result<Arc<ClosureGens>, DecErr> {
        DecData::closure_gens(self, p)
    }
    fn reduced(&self, ev: &mut Eval<F>, p: Tid, c: u64) -> Result<Reduced<F>, DecErr> {
        DecData::reduced(self, ev, p, c)
    }
}

// ------------------------------------------------------------------ shapes and layouts

/// The localized shape of a type: support cosets (ascending ranks), block sizes m_c (reduced
/// coordinates) and the top-stratum coset rep(p), where m_c = 1.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TypeShape {
    pub p: Tid,
    pub cosets: Vec<u64>,
    pub dims: Vec<u32>,
    pub rep: u64,
}

impl TypeShape {
    pub fn new(p: Tid, cosets: Vec<u64>, dims: Vec<u32>, rep: u64) -> Result<TypeShape, SErr> {
        if cosets.len() != dims.len() || cosets.is_empty() {
            return Err(SErr::Shape(format!("type {}: {} cosets, {} dims", p, cosets.len(), dims.len())));
        }
        if !cosets.windows(2).all(|w| w[0] < w[1]) {
            return Err(SErr::Shape(format!("type {}: cosets not strictly ascending", p)));
        }
        if dims.iter().any(|&d| d == 0) {
            return Err(SErr::Shape(format!("type {}: empty block in the support", p)));
        }
        match cosets.binary_search(&rep) {
            Ok(i) if dims[i] == 1 => {}
            Ok(i) => return Err(SErr::Shape(format!("type {}: rep coset {} has m = {} ≠ 1", p, rep, dims[i]))),
            Err(_) => return Err(SErr::Shape(format!("type {}: rep coset {} outside the support", p, rep))),
        }
        Ok(TypeShape { p, cosets, dims, rep })
    }
    /// index of rep(p) in `cosets`
    pub fn rep_index(&self) -> usize {
        self.cosets.binary_search(&self.rep).unwrap()
    }
    /// localized rank Σ m_c
    pub fn rank(&self) -> u64 {
        self.dims.iter().map(|&d| d as u64).sum()
    }
    /// End block size S(p) = Σ m_c²
    pub fn s(&self) -> u64 {
        self.dims.iter().map(|&d| d as u64 * d as u64).sum()
    }
    pub fn dim_at(&self, c: u64) -> u32 {
        self.cosets.binary_search(&c).map_or(0, |i| self.dims[i])
    }
}

/// Layout of the maps src → tgt: the joint support, rows = tgt blocks, cols = src blocks.
pub fn layout_of(src: &TypeShape, tgt: &TypeShape) -> Layout {
    let (a, b) = (&src.cosets, &tgt.cosets);
    let (mut i, mut j) = (0, 0);
    let (mut cos, mut rows, mut cols) = (vec![], vec![], vec![]);
    while i < a.len() && j < b.len() {
        if a[i] < b[j] {
            i += 1;
        } else if b[j] < a[i] {
            j += 1;
        } else {
            cos.push(a[i]);
            cols.push(src.dims[i]);
            rows.push(tgt.dims[j]);
            i += 1;
            j += 1;
        }
    }
    Layout::new(cos, rows, cols)
}

/// rep(p): the coset of the witness configuration of p's table (as sindec's form checks).
pub fn rep_coset(schur: &Schur, lam: &[u8], p: Tid) -> u64 {
    let tab = schur.table(p).clone();
    let cs = crate::coset::CosetSpace::new(lam, &tab.cols());
    cs.rank(&cs.rep(&tab))
}

/// The shape of a library type: flows of its canonical ladder per coset (e_p = None). `rep` is
/// `rep_coset(&builder.schur, lam, p)`, taken before `freeze` (Tids are per Schur instance).
pub fn lib_shape<F: Field, L: SLib<F> + ?Sized>(lib: &L, p: Tid, rep: u64) -> Result<TypeShape, SErr> {
    let t = lib.ty(p)?;
    if t.e.is_some() {
        return Err(SErr::Lib(format!("type {}: idempotent reduced coordinates are M5i", p)));
    }
    let fi = lib.ar().flows(t.lad);
    let dims: Vec<u32> = (0..fi.cosets.len()).map(|i| fi.start[i + 1] - fi.start[i]).collect();
    TypeShape::new(p, fi.cosets.clone(), dims, rep)
}

fn layout_bytes(l: &Layout) -> u64 {
    (l.cosets.len() * 8 + l.rows.len() * 12) as u64 + std::mem::size_of::<Layout>() as u64 + 32
}

/// Registered shapes and the per-scan layout cache keyed by (src type, tgt type).
#[derive(Clone, Debug, Default)]
pub struct ShapeReg {
    shapes: HashMap<Tid, Arc<TypeShape>>,
    lays: HashMap<(Tid, Tid), Arc<Layout>>,
    bytes: u64,
}

impl ShapeReg {
    pub fn new() -> Self {
        Self::default()
    }
    /// registering the same type twice requires the same shape
    pub fn register(&mut self, s: TypeShape) -> Result<(), SErr> {
        if let Some(o) = self.shapes.get(&s.p) {
            if **o != s {
                return Err(SErr::Shape(format!("type {} registered with two shapes", s.p)));
            }
            return Ok(());
        }
        self.bytes += (s.cosets.len() * 12) as u64 + 64;
        self.shapes.insert(s.p, Arc::new(s));
        Ok(())
    }
    pub fn contains(&self, p: Tid) -> bool {
        self.shapes.contains_key(&p)
    }
    pub fn shape(&self, p: Tid) -> Result<&Arc<TypeShape>, SErr> {
        self.shapes.get(&p).ok_or_else(|| SErr::Shape(format!("type {} not registered", p)))
    }
    pub fn layout(&mut self, src: Tid, tgt: Tid) -> Result<Arc<Layout>, SErr> {
        if let Some(l) = self.lays.get(&(src, tgt)) {
            return Ok(l.clone());
        }
        let l = Arc::new(layout_of(self.shape(src)?, self.shape(tgt)?));
        self.bytes += layout_bytes(&l);
        self.lays.insert((src, tgt), l.clone());
        Ok(l)
    }
    /// shapes and layouts
    pub fn bytes(&self) -> u64 {
        self.bytes
    }
    pub fn n_layouts(&self) -> usize {
        self.lays.len()
    }
}

// ------------------------------------------------------------------ the complex

/// An object B_p ⟨t⟩ [h] of the complex.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct SObj {
    pub p: Tid,
    pub h: i32,
    pub t: i32,
}

/// fixed overheads of the byte counter (BTreeMap/BTreeSet entries, measured order of magnitude)
pub const OBJ_BYTES: u64 = 160;
pub const EDGE_BYTES: u64 = 96;

/// A complex of localized objects at one point (SPEC §6.3 "State"). Mutate only through the
/// methods: they keep `inn` and the byte counter consistent.
#[derive(Clone, Debug)]
pub struct SCx<F> {
    objs: BTreeMap<usize, SObj>,
    out: BTreeMap<usize, BTreeMap<usize, LocMap<F>>>,
    inn: BTreeMap<usize, BTreeSet<usize>>,
    next: usize,
    bytes: u64,
    /// number of stored components (kept with the byte counter)
    nedges: u64,
}

impl<F: Field> Default for SCx<F> {
    fn default() -> Self {
        Self::new()
    }
}

impl<F: Field> SCx<F> {
    pub fn new() -> Self {
        SCx { objs: BTreeMap::new(), out: BTreeMap::new(), inn: BTreeMap::new(), next: 0, bytes: 0, nedges: 0 }
    }
    pub fn objs(&self) -> &BTreeMap<usize, SObj> {
        &self.objs
    }
    pub fn out(&self) -> &BTreeMap<usize, BTreeMap<usize, LocMap<F>>> {
        &self.out
    }
    pub fn inn(&self) -> &BTreeMap<usize, BTreeSet<usize>> {
        &self.inn
    }
    pub fn get(&self, x: usize, y: usize) -> Option<&LocMap<F>> {
        self.out.get(&x).and_then(|m| m.get(&y))
    }
    pub fn add_obj(&mut self, o: SObj) -> usize {
        let id = self.next;
        self.next += 1;
        self.objs.insert(id, o);
        self.out.insert(id, BTreeMap::new());
        self.inn.insert(id, BTreeSet::new());
        self.bytes += OBJ_BYTES;
        id
    }
    /// x → y += m; a zero sum removes the component. The layout must be the one of (type x, type y).
    pub fn add_map(&mut self, x: usize, y: usize, m: LocMap<F>) {
        assert!(self.objs.contains_key(&x) && self.objs.contains_key(&y), "SCx::add_map: unknown object");
        let outs = self.out.get_mut(&x).unwrap();
        match outs.get_mut(&y) {
            Some(old) => {
                assert!(Arc::ptr_eq(&old.lay, &m.lay) || *old.lay == *m.lay, "SCx::add_map: layouts differ");
                for (a, b) in old.data.iter_mut().zip(&m.data) {
                    a.add_assign(b);
                }
                if old.is_zero() {
                    let gone = outs.remove(&y).unwrap();
                    self.bytes -= gone.bytes() as u64 + EDGE_BYTES;
                    self.nedges -= 1;
                    self.inn.get_mut(&y).unwrap().remove(&x);
                }
            }
            None => {
                if m.is_zero() {
                    return;
                }
                self.bytes += m.bytes() as u64 + EDGE_BYTES;
                self.nedges += 1;
                outs.insert(y, m);
                self.inn.get_mut(&y).unwrap().insert(x);
            }
        }
    }
    fn take_map(&mut self, x: usize, y: usize) -> Option<LocMap<F>> {
        let m = self.out.get_mut(&x)?.remove(&y)?;
        self.inn.get_mut(&y).unwrap().remove(&x);
        self.bytes -= m.bytes() as u64 + EDGE_BYTES;
        self.nedges -= 1;
        Some(m)
    }
    pub fn remove_obj(&mut self, a: usize) {
        if let Some(outs) = self.out.remove(&a) {
            for (y, m) in outs {
                self.inn.get_mut(&y).unwrap().remove(&a);
                self.bytes -= m.bytes() as u64 + EDGE_BYTES;
                self.nedges -= 1;
            }
        }
        if let Some(ins) = self.inn.remove(&a) {
            for x in ins {
                let m = self.out.get_mut(&x).unwrap().remove(&a).unwrap();
                self.bytes -= m.bytes() as u64 + EDGE_BYTES;
                self.nedges -= 1;
            }
        }
        if self.objs.remove(&a).is_some() {
            self.bytes -= OBJ_BYTES;
        }
    }
    pub fn n_objs(&self) -> usize {
        self.objs.len()
    }
    pub fn components(&self) -> usize {
        debug_assert_eq!(self.nedges as usize, self.out.values().map(|m| m.len()).sum::<usize>());
        self.nedges as usize
    }
    /// stored block entries (all LocMap lengths)
    pub fn stored_entries(&self) -> u64 {
        self.out.values().flat_map(|m| m.values()).map(|l| l.data.len() as u64).sum()
    }
    /// fraction of zero blocks among the stored blocks (§6.8)
    pub fn zero_block_frac(&self) -> f64 {
        let (mut z, mut n) = (0usize, 0usize);
        for l in self.out.values().flat_map(|m| m.values()) {
            let k = l.lay.cosets.len();
            n += k;
            z += k - l.nonzero_blocks();
        }
        if n == 0 {
            0.0
        } else {
            z as f64 / n as f64
        }
    }
    /// the incremental byte counter (LocMap data + OBJ_BYTES per object + EDGE_BYTES per edge)
    pub fn bytes(&self) -> u64 {
        self.bytes
    }
    /// the same quantity recomputed from scratch (tests)
    pub fn recount_bytes(&self) -> u64 {
        self.objs.len() as u64 * OBJ_BYTES + self.out.values().flat_map(|m| m.values()).map(|l| l.bytes() as u64 + EDGE_BYTES).sum::<u64>()
    }
    /// signed K_0 class: (type, t) → Σ (−1)^h
    pub fn k0(&self) -> BTreeMap<(Tid, i32), i64> {
        let mut m: BTreeMap<(Tid, i32), i64> = BTreeMap::new();
        for o in self.objs.values() {
            *m.entry((o.p, o.t)).or_insert(0) += if o.h.rem_euclid(2) == 0 { 1 } else { -1 };
        }
        m.retain(|_, c| *c != 0);
        m
    }
    /// sorted object multiset
    pub fn object_multiset(&self) -> Vec<SObj> {
        let mut v: Vec<SObj> = self.objs.values().copied().collect();
        v.sort();
        v
    }
    /// Σ_objects S(type) (the decategorified End size of the complex)
    pub fn sum_s(&self, reg: &ShapeReg) -> u64 {
        self.objs.values().map(|o| reg.shape(o.p).map_or(0, |s| s.s())).sum()
    }
    /// hash of the objects (id, p, h, t) and the sparsity pattern (x, y)
    pub fn pattern_hash(&self) -> u64 {
        let mut h = Fnv::new();
        h.word(self.objs.len() as u64);
        for (&id, o) in &self.objs {
            h.word(id as u64);
            h.word(o.p as u64);
            h.word(o.h as i64 as u64);
            h.word(o.t as i64 as u64);
        }
        h.word(0xed6e);
        for (&x, outs) in &self.out {
            for &y in outs.keys() {
                h.word(x as u64);
                h.word(y as u64);
            }
        }
        h.finish()
    }
    /// all nonzero components x → z of d∘d (blockwise per coset)
    pub fn d2_nonzero(&self, reg: &mut ShapeReg) -> Result<Vec<(usize, usize)>, SErr> {
        let mut bad = vec![];
        for (&x, outs) in &self.out {
            let mut acc: BTreeMap<usize, LocMap<F>> = BTreeMap::new();
            for (&y, f) in outs {
                for (&z, g) in &self.out[&y] {
                    let lay = reg.layout(self.objs[&x].p, self.objs[&z].p)?;
                    let e = acc.entry(z).or_insert_with(|| LocMap::zeros(lay.clone()));
                    compose_acc(e, g, f, F::one());
                }
            }
            for (z, m) in acc {
                if !m.is_zero() {
                    bad.push((x, z));
                }
            }
        }
        Ok(bad)
    }
}

/// dst(x → z) += alpha · g(y → z) ∘ f(x → y), blockwise over dst's cosets (zero where g or f
/// has no block). No allocation.
pub fn compose_acc<F: Field>(dst: &mut LocMap<F>, g: &LocMap<F>, f: &LocMap<F>, alpha: F) {
    let lay = dst.lay.clone();
    for (i, &c) in lay.cosets.iter().enumerate() {
        let (Some(ig), Some(jf)) = (g.lay.index(c), f.lay.index(c)) else { continue };
        let (m, k, n) = (lay.rows[i] as usize, g.lay.cols[ig] as usize, lay.cols[i] as usize);
        debug_assert!(g.lay.rows[ig] as usize == m && f.lay.cols[jf] as usize == n && f.lay.rows[jf] as usize == k);
        let o = lay.off[i] as usize;
        let len = m * n;
        gemm_acc(&mut dst.data[o..o + len], g.block(ig), f.block(jf), m, k, n, alpha.clone());
    }
}

/// multiply-adds of compose_acc(dst, g, f)
pub fn compose_cost<F: Field>(dst: &LocMap<F>, g: &LocMap<F>, f: &LocMap<F>) -> u64 {
    let mut t = 0u64;
    for (i, &c) in dst.lay.cosets.iter().enumerate() {
        if let (Some(ig), Some(_)) = (g.lay.index(c), f.lay.index(c)) {
            t += dst.lay.rows[i] as u64 * g.lay.cols[ig] as u64 * dst.lay.cols[i] as u64;
        }
    }
    t
}

/// Some(c) if every block of m is c·I (m: p → p over layout(p, p)); the scalar is read at the
/// rep coset (one entry); cost Σ_c m_c².
pub fn scalar_of<F: Field>(m: &LocMap<F>, shape: &TypeShape) -> Option<F> {
    let i0 = m.lay.index(shape.rep)?;
    let c = m.block(i0)[0].clone();
    for i in 0..m.lay.cosets.len() {
        let (r, k) = (m.lay.rows[i] as usize, m.lay.cols[i] as usize);
        if r != k {
            return None;
        }
        let b = m.block(i);
        for a in 0..r {
            for j in 0..r {
                let v = &b[a * r + j];
                if a == j {
                    if *v != c {
                        return None;
                    }
                } else if !v.is_zero() {
                    return None;
                }
            }
        }
    }
    if m.lay.cosets.len() != shape.cosets.len() {
        return None;
    }
    Some(c)
}

// ------------------------------------------------------------------ hashing

/// FNV-1a over 64-bit words (stable across runs and platforms)
#[derive(Clone, Copy, Debug)]
pub struct Fnv(u64);
impl Fnv {
    pub fn new() -> Self {
        Fnv(0xcbf2_9ce4_8422_2325)
    }
    pub fn word(&mut self, w: u64) {
        for b in w.to_le_bytes() {
            self.0 ^= b as u64;
            self.0 = self.0.wrapping_mul(0x0100_0000_01b3);
        }
    }
    pub fn finish(&self) -> u64 {
        self.0
    }
}
impl Default for Fnv {
    fn default() -> Self {
        Self::new()
    }
}

// ------------------------------------------------------------------ frozen structure

/// One elimination round: the stage (crossing index), the pivots in order, and the pattern
/// hash after the round.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Round {
    pub stage: usize,
    pub pivots: Vec<(usize, usize)>,
    pub pattern: u64,
}

/// The structure decided at point 0 (D14): replayed, never re-decided, at every other point.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct FrozenTrace {
    /// pattern hash after every insertion step and every `load`, in call order
    pub checkpoints: Vec<u64>,
    /// every elimination round, in call order
    pub rounds: Vec<Round>,
    /// determinism record per crossing: hash of (sorted objects, pivot sequence, pattern)
    pub stages: Vec<u64>,
}

impl FrozenTrace {
    pub fn pivots(&self) -> usize {
        self.rounds.iter().map(|r| r.pivots.len()).sum()
    }
    /// all pivots in order
    pub fn pivot_sequence(&self) -> Vec<(usize, usize)> {
        self.rounds.iter().flat_map(|r| r.pivots.iter().copied()).collect()
    }
    /// approximate bytes (stored per point set, not per point)
    pub fn bytes(&self) -> u64 {
        (self.checkpoints.len() * 8 + self.stages.len() * 8 + self.rounds.iter().map(|r| 40 + r.pivots.len() * 16).sum::<usize>()) as u64
    }
}

#[derive(Clone, Debug)]
pub enum Mode {
    /// point 0: values decide the pivots; the trace is recorded
    Record,
    /// other points: the recorded structure is applied and checked
    Replay(Arc<FrozenTrace>),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PivotOrder {
    /// lazy versioned min-heap (production)
    Heap,
    /// full rescan per pivot (xengine's order; reference for tests)
    Rescan,
    /// lazy-update Markowitz (production for dense candidate graphs): a popped pair whose cost
    /// grew is re-pushed at its current cost; only fill-in candidates are pushed after a pivot.
    /// Deterministic (structure only), so the frozen trace is as well defined as with `Heap`.
    Lazy,
    /// the M6a-1b implementation of `Heap` (BTreeMap/HashMap lookups per push; reference)
    HeapV1,
}

#[derive(Clone, Debug)]
pub struct EngOpts {
    /// eliminate eligible pivots after every insertion step (§6.4)
    pub streaming: bool,
    /// per-point budget: complex + previous complex + layouts + declared external bytes (Eval memo)
    pub mem_point_bytes: u64,
    pub order: PivotOrder,
    /// d∘d = 0 after every crossing (unpruned scans only, R4 §8)
    pub check_d2: bool,
    /// Err on an arrow at an object that died the wrong way (always on in production)
    pub check_dead: bool,
    /// TEST ONLY: every pivot eligible at every step (the sensitivity test of §6.4)
    pub ignore_eligibility: bool,
    /// the coset whose blocks `PointScan` extracts (the identity coset for closure data)
    pub id_coset: Option<u64>,
    /// record the location string of the memory peak (a format! per new peak)
    pub peak_detail: bool,
}
impl Default for EngOpts {
    fn default() -> Self {
        EngOpts {
            streaming: true,
            mem_point_bytes: 512 << 20,
            order: PivotOrder::Heap,
            check_d2: false,
            check_dead: true,
            ignore_eligibility: false,
            id_coset: None,
            peak_detail: true,
        }
    }
}

/// Per-crossing statistics.
#[derive(Clone, Debug, Default)]
pub struct StageStat {
    pub stage: usize,
    /// objects inserted in this crossing (the unsimplified size)
    pub inserted: usize,
    pub objects: usize,
    pub components: usize,
    pub pivots: usize,
    pub rounds: usize,
    pub stored_entries: u64,
    pub sum_s: u64,
    /// stored entries / Σ S of the objects (§6.8)
    pub c_e: f64,
    pub zero_block_frac: f64,
    /// bytes after the crossing, and the peak during it
    pub bytes: u64,
    pub peak_bytes: u64,
    pub secs: f64,
    pub hash: u64,
}

#[derive(Clone, Debug, Default)]
pub struct MemStats {
    pub peak_bytes: u64,
    pub peak_at: String,
    /// peak of the complex alone (entries + overheads)
    pub peak_cx_bytes: u64,
    pub peak_entries: u64,
    /// objects and components of the complex, and the previous complex's bytes, at the peak
    pub peak_objs: u64,
    pub peak_edges: u64,
    pub peak_prev_bytes: u64,
}

/// The output of one scan at one point (what the box driver stores per point).
#[derive(Clone, Debug)]
pub struct PointScan<F> {
    pub point_id: u32,
    /// final objects by id
    pub objs: Vec<(usize, SObj)>,
    /// final components restricted to `id_coset`: (src id, tgt id, rows, cols, row-major block);
    /// components without a block there are listed with an empty block
    pub id_blocks: Vec<(usize, usize, u32, u32, Vec<F>)>,
    /// the determinism record (one hash per crossing)
    pub trace: Vec<u64>,
    /// the structure (recorded at point 0, the replayed copy elsewhere)
    pub frozen: Arc<FrozenTrace>,
    pub stages: Vec<StageStat>,
    pub mem: MemStats,
    pub k0: BTreeMap<(Tid, i32), i64>,
    pub secs: f64,
}

impl<F: Field> PointScan<F> {
    /// identical final objects and determinism traces (the §6.5 step 1 requirement)
    pub fn same_structure(&self, o: &PointScan<F>) -> Result<(), String> {
        if self.objs != o.objs {
            return Err(format!("points {} and {}: final objects differ", self.point_id, o.point_id));
        }
        if self.trace != o.trace {
            return Err(format!("points {} and {}: determinism traces differ", self.point_id, o.point_id));
        }
        Ok(())
    }
    /// bytes of the extracted values
    pub fn sample_bytes(&self) -> u64 {
        self.id_blocks.iter().map(|b| (b.4.len() * std::mem::size_of::<F>()) as u64 + 24).sum()
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Died {
    Source,
    Target,
}

/// One scan at one point (SPEC §6.3 `PEngine`): the complex, the shapes and layouts, the mode
/// (record or replay), streaming eligibility, memory accounting and statistics.
pub struct PEngine<F: Field> {
    pub reg: ShapeReg,
    cx: SCx<F>,
    prev: Option<Arc<SCx<F>>>,
    prev_bytes: u64,
    pub opts: EngOpts,
    mode: Mode,
    rec: FrozenTrace,
    chk_i: usize,
    round_i: usize,
    stage: usize,
    in_crossing: bool,
    step_open: bool,
    cur_step: Option<u32>,
    steps: HashMap<usize, u32>,
    dead: HashMap<usize, Died>,
    stage_pivots: Vec<(usize, usize)>,
    stage_rounds: usize,
    stage_peak: u64,
    inserted: usize,
    extern_bytes: u64,
    pub mem: MemStats,
    pub stages: Vec<StageStat>,
    t_stage: Instant,
    t0: Instant,
    /// profile of the current crossing: pivot updates, multiply-adds, seconds in apply_pivot
    pub prof: (u64, u64, f64),
    pub prof2: (f64, u64, f64),
    /// edges created by the last pivot (fill-in), for the candidate index of the heap
    new_edges: Vec<(usize, usize)>,
    track_new: bool,
    /// cumulative profile of the whole scan: (apply_pivot secs, pattern-hash secs, updates, multiply-adds)
    pub tprof: (f64, f64, u64, u64),
    /// composition plans per (x type, pivot type, y type): blocks of d_xy += d_ay ∘ d_xb
    plans: HashMap<(Tid, Tid, Tid), Arc<CompPlan>>,
    /// freed component buffers by length (reused by fill-in edges)
    pool: HashMap<u32, Vec<Vec<F>>>,
    pool_bytes: u64,
}

/// Blockwise composition dst(x → y) += g(a → y) ∘ f(x → a) over the registry layouts:
/// per coset of dst where g and f have blocks, (dst off, g off, f off, m, k, n).
#[derive(Debug)]
pub struct CompPlan {
    pub blocks: Vec<(u32, u32, u32, u32, u32, u32)>,
    pub madds: u64,
}

impl CompPlan {
    pub fn new(dst: &Layout, g: &Layout, f: &Layout) -> CompPlan {
        let mut blocks = vec![];
        let mut madds = 0;
        for (i, &c) in dst.cosets.iter().enumerate() {
            let (Some(ig), Some(jf)) = (g.index(c), f.index(c)) else { continue };
            let (m, k, n) = (dst.rows[i], g.cols[ig], dst.cols[i]);
            assert!(g.rows[ig] == m && f.cols[jf] == n && f.rows[jf] == k, "CompPlan: block shapes differ");
            if m == 0 || k == 0 || n == 0 {
                continue;
            }
            madds += m as u64 * k as u64 * n as u64;
            blocks.push((dst.off[i], g.off[ig], f.off[jf], m, k, n));
        }
        CompPlan { blocks, madds }
    }
    /// dst += alpha · g ∘ f (no allocation)
    #[inline]
    pub fn apply<F: Field>(&self, dst: &mut [F], g: &[F], f: &[F], alpha: &F) {
        for &(od, og, of, m, k, n) in &self.blocks {
            let (m, k, n) = (m as usize, k as usize, n as usize);
            let (od, og, of) = (od as usize, og as usize, of as usize);
            for i in 0..m {
                let crow = &mut dst[od + i * n..od + (i + 1) * n];
                for t in 0..k {
                    let x = &g[og + i * k + t];
                    if x.is_zero() {
                        continue;
                    }
                    let s = x.mul(alpha);
                    let brow = &f[of + t * n..of + (t + 1) * n];
                    for (cj, bj) in crow.iter_mut().zip(brow) {
                        if !bj.is_zero() {
                            cj.fma(&s, bj);
                        }
                    }
                }
            }
        }
    }
}

impl<F: Field> PEngine<F> {
    pub fn new(mode: Mode, opts: EngOpts) -> Self {
        PEngine {
            reg: ShapeReg::new(),
            cx: SCx::new(),
            prev: None,
            prev_bytes: 0,
            opts,
            mode,
            rec: FrozenTrace::default(),
            chk_i: 0,
            round_i: 0,
            stage: 0,
            in_crossing: false,
            step_open: false,
            cur_step: None,
            steps: HashMap::new(),
            dead: HashMap::new(),
            stage_pivots: vec![],
            stage_rounds: 0,
            stage_peak: 0,
            inserted: 0,
            extern_bytes: 0,
            mem: MemStats::default(),
            stages: vec![],
            t_stage: Instant::now(),
            t0: Instant::now(),
            prof: (0, 0, 0.0),
            prof2: (0.0, 0, 0.0),
            new_edges: vec![],
            track_new: false,
            tprof: (0.0, 0.0, 0, 0),
            plans: HashMap::new(),
            pool: HashMap::new(),
            pool_bytes: 0,
        }
    }
    pub fn is_replay(&self) -> bool {
        matches!(self.mode, Mode::Replay(_))
    }
    pub fn cx(&self) -> &SCx<F> {
        &self.cx
    }
    pub fn register(&mut self, s: TypeShape) -> Result<(), SErr> {
        self.reg.register(s)
    }
    pub fn layout(&mut self, src: Tid, tgt: Tid) -> Result<Arc<Layout>, SErr> {
        self.reg.layout(src, tgt)
    }
    /// the recorded (or replayed) structure so far
    pub fn frozen(&self) -> &FrozenTrace {
        &self.rec
    }
    /// bytes counted against the budget right now
    pub fn bytes(&self) -> u64 {
        self.cx.bytes() + self.prev_bytes + self.reg.bytes() + self.extern_bytes + self.pool_bytes
    }
    /// declare memory held outside the engine for this point (Eval memo, evaluated library blocks)
    pub fn set_extern_bytes(&mut self, b: u64) {
        self.extern_bytes = b;
    }

    fn check_mem(&mut self, at: &dyn Fn() -> String) -> Result<(), SErr> {
        let b = self.bytes();
        if b > self.stage_peak {
            self.stage_peak = b;
        }
        if b > self.mem.peak_bytes {
            self.mem.peak_bytes = b;
            if self.opts.peak_detail {
                self.mem.peak_at = at();
            }
            self.mem.peak_objs = self.cx.objs.len() as u64;
            self.mem.peak_edges = self.cx.nedges;
            self.mem.peak_prev_bytes = self.prev_bytes;
        }
        if self.cx.bytes() > self.mem.peak_cx_bytes {
            self.mem.peak_cx_bytes = self.cx.bytes();
        }
        if b > self.opts.mem_point_bytes {
            return Err(SErr::Memory { bytes: b, budget: self.opts.mem_point_bytes, at: at() });
        }
        Ok(())
    }

    fn checkpoint(&mut self, what: &str) -> Result<(), SErr> {
        let th = Instant::now();
        let h = self.cx.pattern_hash();
        self.tprof.1 += th.elapsed().as_secs_f64();
        match &self.mode {
            Mode::Record => self.rec.checkpoints.push(h),
            Mode::Replay(fz) => {
                let want = fz.checkpoints.get(self.chk_i).copied();
                if want != Some(h) {
                    return Err(SErr::PointDependent {
                        stage: self.stage,
                        round: self.round_i,
                        what: format!("sparsity pattern after {} (checkpoint {}) differs from point 0", what, self.chk_i),
                    });
                }
                self.rec.checkpoints.push(h);
            }
        }
        self.chk_i += 1;
        Ok(())
    }

    /// Replace the complex without elimination (the start object, or a given complex in tests).
    pub fn load(&mut self, cx: SCx<F>) -> Result<(), SErr> {
        if self.in_crossing {
            return Err(SErr::Usage("load inside a crossing".into()));
        }
        for o in cx.objs.values() {
            self.reg.shape(o.p)?;
        }
        self.cx = cx;
        self.checkpoint("load")?;
        self.check_mem(&|| "load".into())
    }

    /// Start a crossing: the current complex becomes the (read-only) previous complex, returned
    /// to the caller, who builds C ⊗ W term by term into the engine. It is counted against the
    /// budget until `end_crossing`.
    pub fn begin_crossing(&mut self) -> Result<Arc<SCx<F>>, SErr> {
        if self.in_crossing {
            return Err(SErr::Usage("begin_crossing twice".into()));
        }
        let old = std::mem::take(&mut self.cx);
        self.prev_bytes = old.bytes();
        let a = Arc::new(old);
        self.prev = Some(a.clone());
        self.in_crossing = true;
        self.cur_step = None;
        self.steps.clear();
        self.dead.clear();
        self.stage_pivots.clear();
        self.stage_rounds = 0;
        self.stage_peak = 0;
        self.inserted = 0;
        self.t_stage = Instant::now();
        Ok(a)
    }

    /// Open the next insertion step (one Rickard term). Arrows of the term differential must run
    /// from step s to step s+1.
    pub fn begin_step(&mut self) -> Result<u32, SErr> {
        if !self.in_crossing || self.step_open {
            return Err(SErr::Usage("begin_step outside a crossing or with a step open".into()));
        }
        let s = self.cur_step.map_or(0, |s| s + 1);
        if s == 0 {
            // stage seconds = insertion + elimination (not the caller's construction of C ⊗ W)
            self.t_stage = Instant::now();
        }
        self.cur_step = Some(s);
        self.step_open = true;
        Ok(s)
    }

    pub fn add_obj(&mut self, o: SObj) -> Result<usize, SErr> {
        if !self.step_open {
            return Err(SErr::Usage("add_obj without an open step".into()));
        }
        self.reg.shape(o.p)?;
        let id = self.cx.add_obj(o);
        self.steps.insert(id, self.cur_step.unwrap());
        self.inserted += 1;
        Ok(id)
    }

    /// Is the object alive in the complex being built?
    pub fn alive(&self, id: usize) -> bool {
        self.cx.objs.contains_key(&id)
    }

    /// x → y += m. Returns Ok(false) if the arrow was discarded because an endpoint was
    /// eliminated legitimately (x as a pivot target, y as a pivot source); Err(Stream) if it was
    /// eliminated the other way (its out- or in-set was incomplete: the eligibility rule failed).
    pub fn add_map(&mut self, x: usize, y: usize, m: LocMap<F>) -> Result<bool, SErr> {
        let (ax, ay) = (self.alive(x), self.alive(y));
        if !ax || !ay {
            for (id, alive, ok) in [(x, ax, Died::Target), (y, ay, Died::Source)] {
                if alive {
                    continue;
                }
                match self.dead.get(&id) {
                    None => return Err(SErr::Usage(format!("add_map: unknown object {}", id))),
                    Some(d) if *d != ok && self.opts.check_dead => {
                        return Err(SErr::Stream(format!(
                            "arrow {} → {} at object {} which was eliminated as pivot {} before its {} arrows were complete",
                            x,
                            y,
                            id,
                            if *d == Died::Source { "source" } else { "target" },
                            if *d == Died::Source { "outgoing" } else { "incoming" }
                        )))
                    }
                    _ => {}
                }
            }
            return Ok(false);
        }
        let (ox, oy) = (self.cx.objs[&x], self.cx.objs[&y]);
        if oy.h != ox.h + 1 {
            return Err(SErr::Usage(format!("add_map {} → {}: h {} → {}", x, y, ox.h, oy.h)));
        }
        let lay = self.reg.layout(ox.p, oy.p)?;
        if !Arc::ptr_eq(&lay, &m.lay) && *lay != *m.lay {
            return Err(SErr::Shape(format!("add_map {} → {}: layout is not the one of ({}, {})", x, y, ox.p, oy.p)));
        }
        self.cx.add_map(x, y, m);
        let stage = self.stage;
        self.check_mem(&|| format!("stage {} insertion", stage))?;
        Ok(true)
    }

    /// Close the step: pattern checkpoint, memory check, and with streaming the elimination of
    /// every eligible pivot (step(a) < this step).
    pub fn end_step(&mut self) -> Result<(), SErr> {
        if !self.step_open {
            return Err(SErr::Usage("end_step without an open step".into()));
        }
        self.step_open = false;
        let s = self.cur_step.unwrap();
        self.checkpoint(&format!("step {} of stage {}", s, self.stage))?;
        let stage = self.stage;
        self.check_mem(&|| format!("stage {} step {}", stage, s))?;
        if self.opts.streaming {
            self.eliminate_round(false)?;
        }
        Ok(())
    }

    /// Close the crossing: eliminate everything, optional d∘d check, determinism hash, stats.
    pub fn end_crossing(&mut self) -> Result<&StageStat, SErr> {
        if !self.in_crossing || self.step_open {
            return Err(SErr::Usage("end_crossing outside a crossing or with a step open".into()));
        }
        self.eliminate_round(true)?;
        if std::env::var("KRHOM_SPROF").is_ok() {
            eprintln!("[sprof] stage {}: {} pivots, {} updates, {:.3e} mult-adds, {:.2}s in apply_pivot, {:.2}s stage; rounds {:.2}s, push loop {:.2}s, touched {}, edges {}", self.stage, self.stage_pivots.len(), self.prof.0, self.prof.1 as f64, self.prof.2, self.t_stage.elapsed().as_secs_f64(), self.prof2.0, self.prof2.2, self.prof2.1, self.cx.components());
            self.prof2 = (0.0, 0, 0.0);
            self.prof = (0, 0, 0.0);
        }
        if self.opts.check_d2 {
            let bad = self.cx.d2_nonzero(&mut self.reg)?;
            if let Some(&(x, z)) = bad.first() {
                return Err(SErr::D2(format!("stage {}: {} nonzero components, first {} → {}", self.stage, bad.len(), x, z)));
            }
        }
        // determinism record
        let mut h = Fnv::new();
        for o in self.cx.object_multiset() {
            h.word(o.p as u64);
            h.word(o.h as i64 as u64);
            h.word(o.t as i64 as u64);
        }
        h.word(0x9175);
        for &(a, b) in &self.stage_pivots {
            h.word(a as u64);
            h.word(b as u64);
        }
        h.word(self.cx.pattern_hash());
        let hash = h.finish();
        match &self.mode {
            Mode::Record => self.rec.stages.push(hash),
            Mode::Replay(fz) => {
                if fz.stages.get(self.stage) != Some(&hash) {
                    return Err(SErr::PointDependent { stage: self.stage, round: self.round_i, what: "determinism hash of the crossing differs from point 0".into() });
                }
                self.rec.stages.push(hash);
            }
        }
        self.prev = None;
        self.prev_bytes = 0;
        let stage = self.stage;
        self.check_mem(&|| format!("end of stage {}", stage))?;
        let entries = self.cx.stored_entries();
        if entries > self.mem.peak_entries {
            self.mem.peak_entries = entries;
        }
        let sum_s = self.cx.sum_s(&self.reg);
        let st = StageStat {
            stage: self.stage,
            inserted: self.inserted,
            objects: self.cx.n_objs(),
            components: self.cx.components(),
            pivots: self.stage_pivots.len(),
            rounds: self.stage_rounds,
            stored_entries: entries,
            sum_s,
            c_e: if sum_s == 0 { 0.0 } else { entries as f64 / sum_s as f64 },
            zero_block_frac: self.cx.zero_block_frac(),
            bytes: self.bytes(),
            peak_bytes: self.stage_peak,
            secs: self.t_stage.elapsed().as_secs_f64(),
            hash,
        };
        self.stages.push(st);
        self.stage += 1;
        self.in_crossing = false;
        self.cur_step = None;
        Ok(self.stages.last().unwrap())
    }

    /// Finish the scan at this point: the per-point output and the final complex.
    pub fn finish(self, point_id: u32) -> Result<(PointScan<F>, SCx<F>), SErr> {
        if self.in_crossing {
            return Err(SErr::Usage("finish inside a crossing".into()));
        }
        if let Mode::Replay(fz) = &self.mode {
            if self.rec != **fz {
                return Err(SErr::PointDependent { stage: self.stage, round: self.round_i, what: "the replay did not consume the whole frozen trace".into() });
            }
        }
        let mut id_blocks = vec![];
        if let Some(c) = self.opts.id_coset {
            for (&x, outs) in &self.cx.out {
                for (&y, m) in outs {
                    match m.lay.index(c) {
                        Some(i) => id_blocks.push((x, y, m.lay.rows[i], m.lay.cols[i], m.block(i).to_vec())),
                        None => id_blocks.push((x, y, 0, 0, vec![])),
                    }
                }
            }
        }
        let ps = PointScan {
            point_id,
            objs: self.cx.objs.iter().map(|(&i, &o)| (i, o)).collect(),
            id_blocks,
            trace: self.rec.stages.clone(),
            frozen: Arc::new(self.rec),
            stages: self.stages,
            mem: self.mem,
            k0: self.cx.k0(),
            secs: self.t0.elapsed().as_secs_f64(),
        };
        Ok((ps, self.cx))
    }

    // -------------------------------------------------------------- elimination

    fn eligible(&self, a: usize, all: bool) -> bool {
        if all || self.opts.ignore_eligibility {
            return true;
        }
        match (self.steps.get(&a), self.cur_step) {
            (Some(&s), Some(cur)) => s < cur,
            _ => true,
        }
    }

    fn is_candidate(&self, a: usize, b: usize, all: bool) -> bool {
        let (oa, ob) = (self.cx.objs[&a], self.cx.objs[&b]);
        oa.p == ob.p && oa.t == ob.t && self.eligible(a, all)
    }

    fn cost(&self, a: usize, b: usize) -> u64 {
        (self.cx.inn[&b].len() as u64 - 1) * (self.cx.out[&a].len() as u64 - 1)
    }

    fn eliminate_round(&mut self, all: bool) -> Result<(), SErr> {
        let tr = Instant::now();
        let r = self.eliminate_round_inner(all);
        self.prof2.0 += tr.elapsed().as_secs_f64();
        r
    }
    fn eliminate_round_inner(&mut self, all: bool) -> Result<(), SErr> {
        let round_no = self.round_i;
        let pivots = match self.mode.clone() {
            Mode::Record => match self.opts.order {
                PivotOrder::Heap => self.markowitz_heap_fast(all)?,
                PivotOrder::HeapV1 => self.markowitz_heap(all)?,
                PivotOrder::Rescan => self.markowitz_rescan(all)?,
                PivotOrder::Lazy => self.markowitz_lazy(all)?,
            },
            Mode::Replay(fz) => {
                let Some(r) = fz.rounds.get(round_no) else {
                    return Err(SErr::PointDependent { stage: self.stage, round: round_no, what: "more elimination rounds than at point 0".into() });
                };
                if r.stage != self.stage {
                    return Err(SErr::PointDependent { stage: self.stage, round: round_no, what: format!("round belongs to stage {} at point 0", r.stage) });
                }
                for &(a, b) in &r.pivots {
                    if !(self.alive(a) && self.alive(b) && self.cx.get(a, b).is_some() && self.is_candidate(a, b, all)) {
                        return Err(SErr::PointDependent { stage: self.stage, round: round_no, what: format!("recorded pivot {} → {} is not a component here", a, b) });
                    }
                    self.apply_pivot(a, b).map_err(|e| match e {
                        SErr::Pivot(s) => SErr::PointDependent { stage: self.stage, round: round_no, what: format!("recorded pivot {} → {}: {}", a, b, s) },
                        e => e,
                    })?;
                }
                r.pivots.clone()
            }
        };
        let th = Instant::now();
        let pattern = self.cx.pattern_hash();
        self.tprof.1 += th.elapsed().as_secs_f64();
        if let Mode::Replay(fz) = &self.mode {
            if fz.rounds[round_no].pattern != pattern {
                return Err(SErr::PointDependent { stage: self.stage, round: round_no, what: "sparsity pattern after the round differs from point 0".into() });
            }
        }
        self.stage_pivots.extend_from_slice(&pivots);
        self.stage_rounds += 1;
        self.rec.rounds.push(Round { stage: self.stage, pivots, pattern });
        self.round_i += 1;
        Ok(())
    }

    fn markowitz_lazy(&mut self, all: bool) -> Result<Vec<(usize, usize)>, SErr> {
        let mut heap: BinaryHeap<Reverse<(u64, usize, usize)>> = BinaryHeap::new();
        for (&a, outs) in &self.cx.out {
            for &b in outs.keys() {
                if self.is_candidate(a, b, all) {
                    heap.push(Reverse((self.cost(a, b), a, b)));
                }
            }
        }
        self.track_new = true;
        let mut done = vec![];
        while let Some(Reverse((c, a, b))) = heap.pop() {
            if !self.alive(a) || !self.alive(b) || self.cx.get(a, b).is_none() {
                continue;
            }
            let cur = self.cost(a, b);
            if cur > c {
                heap.push(Reverse((cur, a, b)));
                continue;
            }
            self.new_edges.clear();
            if let Err(e) = self.apply_pivot(a, b) {
                self.track_new = false;
                return Err(e);
            }
            done.push((a, b));
            let ne = std::mem::take(&mut self.new_edges);
            for &(x, y) in &ne {
                if self.alive(x) && self.alive(y) && self.cx.get(x, y).is_some() && self.is_candidate(x, y, all) {
                    heap.push(Reverse((self.cost(x, y), x, y)));
                }
            }
            self.new_edges = ne;
        }
        self.track_new = false;
        Ok(done)
    }

    /// Reference order: rescan every candidate per pivot, minimum of (cost, a, b).
    fn markowitz_rescan(&mut self, all: bool) -> Result<Vec<(usize, usize)>, SErr> {
        let mut done = vec![];
        loop {
            let mut best: Option<(u64, usize, usize)> = None;
            for (&a, outs) in &self.cx.out {
                for &b in outs.keys() {
                    if !self.is_candidate(a, b, all) {
                        continue;
                    }
                    let k = (self.cost(a, b), a, b);
                    if best.map_or(true, |x| k < x) {
                        best = Some(k);
                    }
                }
            }
            let Some((_, a, b)) = best else { break };
            self.apply_pivot(a, b)?;
            done.push((a, b));
        }
        Ok(done)
    }

    /// `Heap` order (= the rescan order: minimum of (cost, a, b) over all candidates) with
    /// dense per-object arrays: degrees, versions, (type, t) keys, eligibility (fixed within a
    /// round) and candidate neighbour lists; after a pivot only the touched objects' candidate
    /// pairs are re-pushed; stale entries are skipped by version.
    fn markowitz_heap_fast(&mut self, all: bool) -> Result<Vec<(usize, usize)>, SErr> {
        let n = self.cx.next;
        let mut key: Vec<(Tid, i32)> = vec![(Tid::MAX, i32::MIN); n];
        let mut alive = vec![false; n];
        let mut elig = vec![false; n];
        let mut indeg = vec![0u32; n];
        let mut outdeg = vec![0u32; n];
        let mut ver = vec![0u32; n];
        for (&id, o) in &self.cx.objs {
            key[id] = (o.p, o.t);
            alive[id] = true;
            elig[id] = self.eligible(id, all);
            indeg[id] = self.cx.inn[&id].len() as u32;
            outdeg[id] = self.cx.out[&id].len() as u32;
        }
        let cost = |indeg: &[u32], outdeg: &[u32], a: usize, b: usize| (indeg[b].max(1) as u64 - 1) * (outdeg[a].max(1) as u64 - 1);
        let mut c_out: Vec<Vec<u32>> = vec![vec![]; n];
        let mut c_in: Vec<Vec<u32>> = vec![vec![]; n];
        let mut heap: BinaryHeap<Reverse<(u64, u32, u32, u32, u32)>> = BinaryHeap::new();
        for (&x, outs) in &self.cx.out {
            for &y in outs.keys() {
                if key[x] == key[y] {
                    c_out[x].push(y as u32);
                    c_in[y].push(x as u32);
                    if elig[x] {
                        heap.push(Reverse((cost(&indeg, &outdeg, x, y), x as u32, y as u32, 0, 0)));
                    }
                }
            }
        }
        self.track_new = true;
        let mut done = vec![];
        let tpre = Instant::now();
        let mut npush = 0u64;
        let (mut npush_in, mut npop, mut maxheap, mut tap) = (0u64, 0u64, 0usize, 0.0f64);
        let mut compact_at = (2 * heap.len()).max(1 << 16);
        while let Some(Reverse((_, a, b, va, vb))) = heap.pop() {
            npop += 1;
            maxheap = maxheap.max(heap.len());
            if heap.len() > compact_at {
                // drop stale entries (dead objects, old versions) and duplicates: O(n), and the
                // pop order of the valid entries is unchanged
                let mut v = std::mem::take(&mut heap).into_vec();
                v.retain(|Reverse((_, x, y, vx, vy))| alive[*x as usize] && alive[*y as usize] && ver[*x as usize] == *vx && ver[*y as usize] == *vy);
                v.sort_unstable();
                v.dedup();
                heap = BinaryHeap::from(v);
                compact_at = (2 * heap.len()).max(1 << 16);
            }
            let (a, b) = (a as usize, b as usize);
            if !alive[a] || !alive[b] || ver[a] != va || ver[b] != vb {
                continue;
            }
            if self.cx.get(a, b).is_none() {
                continue;
            }
            self.new_edges.clear();
            let tq = Instant::now();
            let touched = match self.apply_pivot(a, b) {
                Ok(t) => t,
                Err(e) => {
                    self.track_new = false;
                    return Err(e);
                }
            };
            tap += tq.elapsed().as_secs_f64();
            done.push((a, b));
            alive[a] = false;
            alive[b] = false;
            for &(x, y) in &self.new_edges {
                if key[x] == key[y] {
                    c_out[x].push(y as u32);
                    c_in[y].push(x as u32);
                }
            }
            let tpush = Instant::now();
            for &t in &touched {
                ver[t] += 1;
                indeg[t] = self.cx.inn[&t].len() as u32;
                outdeg[t] = self.cx.out[&t].len() as u32;
            }
            for &t in &touched {
                let vt = ver[t];
                c_out[t].retain(|&y| alive[y as usize]);
                c_in[t].retain(|&x| alive[x as usize]);
                if elig[t] {
                    for &y in &c_out[t] {
                        npush += 1;
                        heap.push(Reverse((cost(&indeg, &outdeg, t, y as usize), t as u32, y, vt, ver[y as usize])));
                    }
                }
                for &x in &c_in[t] {
                    if elig[x as usize] {
                        npush_in += 1;
                        heap.push(Reverse((cost(&indeg, &outdeg, x as usize, t), x, t as u32, ver[x as usize], vt)));
                    }
                }
            }
            self.prof2.2 += tpush.elapsed().as_secs_f64();
            self.prof2.1 += touched.len() as u64;
        }
        self.track_new = false;
        if std::env::var("KRHOM_SPROF2").is_ok() {
            eprintln!("[sprof2] round: {} pivots, {}+{} pushes, {} pops, max heap {}, heap loop {:.2}s (apply {:.2}s)", done.len(), npush, npush_in, npop, maxheap, tpre.elapsed().as_secs_f64(), tap);
        }
        Ok(done)
    }

    /// Lazy versioned min-heap: entries (cost, a, b, ver a, ver b); after a pivot every object
    /// whose in- or out-set changed gets a new version and its candidate pairs are re-pushed;
    /// stale entries are skipped. Equivalent to the rescan order.
    fn markowitz_heap(&mut self, all: bool) -> Result<Vec<(usize, usize)>, SErr> {
        let mut ver: HashMap<usize, u32> = HashMap::new();
        let mut heap: BinaryHeap<Reverse<(u64, usize, usize, u32, u32)>> = BinaryHeap::new();
        let tf = Instant::now();
        for (&a, outs) in &self.cx.out {
            for &b in outs.keys() {
                if self.is_candidate(a, b, all) {
                    heap.push(Reverse((self.cost(a, b), a, b, 0, 0)));
                }
            }
        }
        // candidates join objects of equal (type, t): per object the neighbours of equal key
        // (built here, extended by the fill-in edges of every pivot); the re-push after a pivot
        // scans only these (stale entries are skipped by the edge test): the same pairs as a
        // scan of all neighbours
        let maxid = self.cx.objs.keys().next_back().map_or(0, |&m| m + 1);
        let mut keyv: Vec<(Tid, i32)> = vec![(Tid::MAX, i32::MIN); maxid];
        for (&id, o) in &self.cx.objs {
            keyv[id] = (o.p, o.t);
        }
        let mut c_out: Vec<Vec<usize>> = vec![vec![]; maxid];
        let mut c_in: Vec<Vec<usize>> = vec![vec![]; maxid];
        for (&x, outs) in &self.cx.out {
            for &y in outs.keys() {
                if keyv[x] == keyv[y] {
                    c_out[x].push(y);
                    c_in[y].push(x);
                }
            }
        }
        self.track_new = true;
        let mut npush = 0u64;
        let tpre = Instant::now();
        let mut done = vec![];
        while let Some(Reverse((_, a, b, va, vb))) = heap.pop() {
            if !self.alive(a) || !self.alive(b) {
                continue;
            }
            if ver.get(&a).copied().unwrap_or(0) != va || ver.get(&b).copied().unwrap_or(0) != vb {
                continue;
            }
            if self.cx.get(a, b).is_none() {
                continue;
            }
            self.new_edges.clear();
            let touched = match self.apply_pivot(a, b) {
                Ok(t) => t,
                Err(e) => {
                    self.track_new = false;
                    return Err(e);
                }
            };
            done.push((a, b));
            for &(x, y) in &self.new_edges {
                if keyv[x] == keyv[y] {
                    c_out[x].push(y);
                    c_in[y].push(x);
                }
            }
            let tpush = Instant::now();
            for &t in &touched {
                *ver.entry(t).or_insert(0) += 1;
            }
            for &t in &touched {
                let vt = ver[&t];
                c_out[t].retain(|y| self.cx.out[&t].contains_key(y));
                c_in[t].retain(|x| self.cx.inn[&t].contains(x));
                for &y in &c_out[t] {
                    if self.is_candidate(t, y, all) {
                        npush += 1;
                        heap.push(Reverse((self.cost(t, y), t, y, vt, ver.get(&y).copied().unwrap_or(0))));
                    }
                }
                for &x in &c_in[t] {
                    if self.is_candidate(x, t, all) {
                        heap.push(Reverse((self.cost(x, t), x, t, ver.get(&x).copied().unwrap_or(0), vt)));
                    }
                }
            }
            self.prof2.2 += tpush.elapsed().as_secs_f64();
            self.prof2.1 += touched.len() as u64;
        }
        self.track_new = false;
        if std::env::var("KRHOM_SPROF2").is_ok() {
            eprintln!("[sprof2] round: {} pivots, {} pushes (out side), heap loop {:.2}s, heap len at end {}", done.len(), npush, tpre.elapsed().as_secs_f64(), heap.len());
        }
        Ok(done)
    }

    /// Gaussian elimination of the scalar pivot a → b: d_xy −= d_ay c⁻¹ d_xb for x ∈ in(b)∖a,
    /// y ∈ out(a)∖b, blockwise per coset; a and b are removed. Returns the objects whose in- or
    /// out-set changed.
    fn apply_pivot(&mut self, a: usize, b: usize) -> Result<Vec<usize>, SErr> {
        let tp = Instant::now();
        let (u0, m0) = (self.prof.0, self.prof.1);
        let r = self.apply_pivot_inner(a, b);
        let dt = tp.elapsed().as_secs_f64();
        self.prof.2 += dt;
        self.tprof.0 += dt;
        self.tprof.2 += self.prof.0 - u0;
        self.tprof.3 += self.prof.1 - m0;
        r
    }
    fn apply_pivot_inner(&mut self, a: usize, b: usize) -> Result<Vec<usize>, SErr> {
        let (oa, ob) = (self.cx.objs[&a], self.cx.objs[&b]);
        if ob.h != oa.h + 1 || oa.p != ob.p || oa.t != ob.t {
            return Err(SErr::Pivot(format!("{} → {} is not a degree-0 component between equal types", a, b)));
        }
        let shape = self.reg.shape(oa.p)?.clone();
        let dab = self.cx.get(a, b).ok_or_else(|| SErr::Pivot(format!("{} → {} is not a component", a, b)))?;
        let c = scalar_of(dab, &shape).ok_or_else(|| SErr::Pivot(format!("component {} → {} (type {}) is not c·I at every coset (End^0 = k violated)", a, b, oa.p)))?;
        if c.is_zero() {
            return Err(SErr::Pivot(format!("component {} → {} has scalar 0 at rep but is nonzero", a, b)));
        }
        let alpha = c.inv().neg();
        let mut touched: Vec<usize> = Vec::with_capacity(self.cx.inn[&a].len() + self.cx.inn[&b].len() + self.cx.out[&a].len() + self.cx.out[&b].len());
        touched.extend(self.cx.inn[&a].iter().copied());
        touched.extend(self.cx.inn[&b].iter().copied());
        touched.extend(self.cx.out[&a].keys().copied());
        touched.extend(self.cx.out[&b].keys().copied());
        touched.sort_unstable();
        touched.dedup();
        touched.retain(|&t| t != a && t != b);
        let xs: Vec<usize> = self.cx.inn[&b].iter().copied().filter(|&x| x != a).collect();
        let ys: Vec<usize> = self.cx.out[&a].keys().copied().filter(|&y| y != b).collect();
        if !xs.is_empty() && !ys.is_empty() {
            let mut day: Vec<(usize, LocMap<F>)> = Vec::with_capacity(ys.len());
            for &y in &ys {
                day.push((y, self.cx.take_map(a, y).unwrap()));
            }
            let ptypes: Vec<Tid> = day.iter().map(|(y, _)| self.cx.objs[y].p).collect();
            // (x type, y type) → (layout, plan), local to this pivot (few types)
            let mut local: Vec<((Tid, Tid), Arc<Layout>, Arc<CompPlan>)> = vec![];
            let mut sel: Vec<usize> = vec![0; day.len()];
            for &x in &xs {
                let dxb = self.cx.take_map(x, b).unwrap();
                let px = self.cx.objs[&x].p;
                for (j, ((_, g), &py)) in day.iter().zip(&ptypes).enumerate() {
                    sel[j] = match local.iter().position(|e| e.0 == (px, py)) {
                        Some(i) => i,
                        None => {
                            let lay = self.reg.layout(px, py)?;
                            let plan = match self.plans.get(&(px, oa.p, py)) {
                                Some(p) => p.clone(),
                                None => {
                                    let p = Arc::new(CompPlan::new(&lay, &g.lay, &dxb.lay));
                                    self.plans.insert((px, oa.p, py), p.clone());
                                    p
                                }
                            };
                            local.push(((px, py), lay, plan));
                            local.len() - 1
                        }
                    };
                }
                let xo = self.cx.out.get_mut(&x).unwrap();
                let mut zeroed: Vec<usize> = vec![];
                for (j, (y, g)) in day.iter().enumerate() {
                    let (_, lay, plan) = &local[sel[j]];
                    let dst = match xo.entry(*y) {
                        std::collections::btree_map::Entry::Occupied(e) => e.into_mut(),
                        std::collections::btree_map::Entry::Vacant(e) => {
                            if self.track_new {
                                self.new_edges.push((x, *y));
                            }
                            let n = lay.len as usize;
                            let data = match self.pool.get_mut(&lay.len).and_then(|v| v.pop()) {
                                Some(mut d) => {
                                    self.pool_bytes -= (d.capacity() * std::mem::size_of::<F>()) as u64;
                                    d.clear();
                                    d.resize(n, F::zero());
                                    d
                                }
                                None => vec![F::zero(); n],
                            };
                            let z = LocMap { lay: lay.clone(), data };
                            self.cx.bytes += z.bytes() as u64 + EDGE_BYTES;
                            self.cx.nedges += 1;
                            self.cx.inn.get_mut(y).unwrap().insert(x);
                            e.insert(z)
                        }
                    };
                    self.prof.0 += 1;
                    self.prof.1 += plan.madds;
                    plan.apply(&mut dst.data, &g.data, &dxb.data, &alpha);
                    if dst.is_zero() {
                        zeroed.push(*y);
                    }
                }
                for y in zeroed {
                    if let Some(m) = self.cx.take_map(x, y) {
                        self.recycle(m);
                    }
                }
                self.recycle(dxb);
            }
            for (_, g) in day {
                self.recycle(g);
            }
            let stage = self.stage;
            self.check_mem(&|| format!("stage {} elimination", stage))?;
        }
        self.cx.remove_obj(a);
        self.cx.remove_obj(b);
        if self.in_crossing {
            self.dead.insert(a, Died::Source);
            self.dead.insert(b, Died::Target);
        }
        Ok(touched)
    }

    /// keep a freed component buffer for reuse by fill-in edges (bounded; counted as external)
    fn recycle(&mut self, m: LocMap<F>) {
        let b = (m.data.capacity() * std::mem::size_of::<F>()) as u64;
        if self.pool_bytes + b > POOL_BYTES {
            return;
        }
        self.pool_bytes += b;
        self.pool.entry(m.lay.len).or_default().push(m.data);
    }
}

/// cap of the buffer pool of `PEngine` (freed component data reused by fill-in edges)
pub const POOL_BYTES: u64 = 16 << 20;

// ------------------------------------------------------------------ closure Gram at a point

/// The Gram matrix of one final type at one point: G[b][a] = (h_b ∘ g_a) at the identity coset,
/// with g = ins (r0 → L_p) and h = outs (L_p → r0); degrees from the Arena.
#[derive(Clone, Debug, PartialEq)]
pub struct GramAtPoint<F> {
    pub p: Tid,
    pub point_id: u32,
    pub id_coset: u64,
    pub ins: Vec<i32>,
    pub outs: Vec<i32>,
    pub g: Vec<Vec<F>>,
}

/// The identity coset of λ: the single coset of the empty ladder r0 (source of every `ins`).
pub fn identity_coset<F: Field>(ar: &Arena<F>, gens: &ClosureGens) -> Result<u64, SErr> {
    let r0 = if let Some(&m) = gens.ins.first() {
        ar.info(m).src
    } else if let Some(&m) = gens.outs.first() {
        ar.info(m).tgt
    } else {
        return Err(SErr::Lib(format!("type {}: no closure generators", gens.p)));
    };
    let fi = ar.flows(r0);
    if fi.cosets.len() != 1 || fi.nflows() != 1 {
        return Err(SErr::Lib(format!("closure generators of {} do not start at the empty ladder", gens.p)));
    }
    Ok(fi.cosets[0])
}

pub fn gram_at_point<F: Field, L: SLib<F> + ?Sized>(lib: &L, ev: &mut Eval<F>, p: Tid) -> Result<GramAtPoint<F>, SErr> {
    let gens = lib.closure_gens(p)?;
    let ar = lib.ar();
    let c = identity_coset(ar, &gens)?;
    let red = lib.reduced(ev, p, c)?;
    let gs: Vec<Vec<F>> = gens
        .ins
        .iter()
        .map(|&m| {
            let b = ev.block(m, c); // flows × 1
            let mut v = vec![F::zero(); red.m];
            gemm_acc(&mut v, &red.q.a, &b.a, red.m, b.rows as usize, 1, F::one());
            v
        })
        .collect();
    let hs: Vec<Vec<F>> = gens
        .outs
        .iter()
        .map(|&m| {
            let b = ev.block(m, c); // 1 × flows
            let mut v = vec![F::zero(); red.m];
            gemm_acc(&mut v, &b.a, &red.p.a, 1, b.cols as usize, red.m, F::one());
            v
        })
        .collect();
    let g = hs.iter().map(|h| gs.iter().map(|g| h.iter().zip(g).fold(F::zero(), |acc, (x, y)| acc.add(&x.mul(y)))).collect()).collect();
    Ok(GramAtPoint {
        p,
        point_id: ev.point().id,
        id_coset: c,
        ins: gens.ins.iter().map(|&m| ar.info(m).deg).collect(),
        outs: gens.outs.iter().map(|&m| ar.info(m).deg).collect(),
        g,
    })
}

// ------------------------------------------------------------------ parallel points

/// Map over items on KRHOM_THREADS scoped threads (order preserved).
pub fn par_map<T: Sync, R: Send, Fun: Fn(&T) -> R + Sync>(items: &[T], f: Fun) -> Vec<R> {
    let nt = crate::num_threads().min(items.len()).max(1);
    if nt == 1 {
        return items.iter().map(&f).collect();
    }
    let next = std::sync::atomic::AtomicUsize::new(0);
    let mut out: Vec<Option<R>> = (0..items.len()).map(|_| None).collect();
    let slots: Vec<std::sync::Mutex<&mut Option<R>>> = out.iter_mut().map(std::sync::Mutex::new).collect();
    std::thread::scope(|s| {
        for _ in 0..nt {
            s.spawn(|| loop {
                let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                if i >= items.len() {
                    break;
                }
                let r = f(&items[i]);
                **slots[i].lock().unwrap() = Some(r);
            });
        }
    });
    drop(slots);
    out.into_iter().map(|r| r.unwrap()).collect()
}

// ------------------------------------------------------------------ synthetic complexes

/// Synthetic localized complexes with known minimal models (tests and measurements).
///
/// Model: types are random shapes; a component x → y has degree e = t_x − t_y ≥ 0, and degree 0
/// only between equal types as c·I (Hom^0(B_p, B_q) = δ_pq k, Hom^{<0} = 0). A complex is
/// P (D_min ⊕ cones) P⁻¹ with D_min minimal (arrows of degree > 0 from even h to h+1 only, so
/// d² = 0), cones c(x)·I between equal objects in adjacent h, and P = I + N per h, N of degree
/// ≥ 0 and strictly increasing in (−t, index) (unitriangular). Every value is α + β·x with
/// (α, β) drawn from the structure seed in a fixed order, so the structure is the same at
/// every generic point x and degenerates exactly where a cone scalar vanishes.
pub mod synth {
    use super::*;
    use rand::{rngs::StdRng, Rng, SeedableRng};

    #[derive(Clone, Debug)]
    pub struct Params {
        pub ntypes: usize,
        /// cosets are drawn from 0..max_cosets (< 64)
        pub max_cosets: u64,
        pub max_dim: u32,
        pub n_min: usize,
        pub n_cones: usize,
        pub hmin: i32,
        pub hmax: i32,
        pub tmax: i32,
        /// probability of a minimal arrow / of an N entry
        pub density: f64,
        pub conj: bool,
        /// type ids are tid_base + 0..ntypes
        pub tid_base: Tid,
    }
    impl Default for Params {
        fn default() -> Self {
            Params { ntypes: 3, max_cosets: 5, max_dim: 3, n_min: 4, n_cones: 3, hmin: 0, hmax: 2, tmax: 3, density: 0.5, conj: true, tid_base: 1 }
        }
    }

    /// A generated complex at one point, its minimal model, and the cone roots.
    pub struct Synth<F: Field> {
        pub shapes: Vec<TypeShape>,
        pub cx: SCx<F>,
        pub min: SCx<F>,
        /// points x where some cone scalar vanishes (degenerate points)
        pub cone_roots: Vec<F>,
    }

    struct Vals<F> {
        rng: StdRng,
        x: F,
    }
    impl<F: Field> Vals<F> {
        fn v(&mut self) -> F {
            let a = F::from_i64(self.rng.gen_range(-1_000_000i64..1_000_000));
            let b = F::from_i64(self.rng.gen_range(-1_000_000i64..1_000_000));
            a.add(&b.mul(&self.x))
        }
        fn ab(&mut self) -> (F, F) {
            let a = F::from_i64(self.rng.gen_range(1i64..1_000_000));
            let b = F::from_i64(self.rng.gen_range(1i64..1_000_000));
            (a, b)
        }
    }

    pub fn shapes(seed: u64, p: &Params) -> Vec<TypeShape> {
        let mut rng = StdRng::seed_from_u64(seed ^ 0x5ba9e5);
        (0..p.ntypes)
            .map(|i| {
                let mut cos: Vec<u64> = (0..p.max_cosets).filter(|_| rng.gen_bool(0.6)).collect();
                if cos.is_empty() {
                    cos.push(rng.gen_range(0..p.max_cosets));
                }
                let rep = cos[rng.gen_range(0..cos.len())];
                let dims: Vec<u32> = cos.iter().map(|&c| if c == rep { 1 } else { rng.gen_range(1..=p.max_dim) }).collect();
                TypeShape::new(p.tid_base + i as Tid, cos, dims, rep).unwrap()
            })
            .collect()
    }

    /// dense inverse (Gauss-Jordan); panics if singular
    fn inverse<F: Field>(a: &[Vec<F>]) -> Vec<Vec<F>> {
        let n = a.len();
        let mut m: Vec<Vec<F>> = a.iter().enumerate().map(|(i, r)| {
            let mut row = r.clone();
            row.extend((0..n).map(|j| if i == j { F::one() } else { F::zero() }));
            row
        }).collect();
        for c in 0..n {
            let p = (c..n).find(|&i| !m[i][c].is_zero()).expect("synth: singular P");
            m.swap(c, p);
            let inv = m[c][c].inv();
            for v in m[c].iter_mut() {
                *v = v.mul(&inv);
            }
            let pr = m[c].clone();
            for i in 0..n {
                if i != c && !m[i][c].is_zero() {
                    let f = m[i][c].clone();
                    for (v, w) in m[i].iter_mut().zip(&pr) {
                        v.sub_assign(&f.mul(w));
                    }
                }
            }
        }
        m.into_iter().map(|r| r[n..].to_vec()).collect()
    }
    fn matmul<F: Field>(a: &[Vec<F>], b: &[Vec<F>]) -> Vec<Vec<F>> {
        let (n, k, m) = (a.len(), b.len(), b.first().map_or(0, |r| r.len()));
        let mut c = vec![vec![F::zero(); m]; n];
        for i in 0..n {
            for t in 0..k {
                if a[i][t].is_zero() {
                    continue;
                }
                for j in 0..m {
                    c[i][j].fma(&a[i][t], &b[t][j]);
                }
            }
        }
        c
    }

    /// Generate the complex of `seed` at the point x (types must be registered in `reg`).
    pub fn generate<F: Field>(seed: u64, p: &Params, x: F, reg: &mut ShapeReg) -> Synth<F> {
        let shapes = shapes(seed, p);
        for s in &shapes {
            reg.register(s.clone()).unwrap();
        }
        let mut rs = StdRng::seed_from_u64(seed);
        let mut vals = Vals { rng: StdRng::seed_from_u64(seed ^ 0x7a1e5), x: x.clone() };
        let ty = |rs: &mut StdRng| shapes[rs.gen_range(0..shapes.len())].p;
        // objects: minimal ones, then cone pairs
        let mut objs: Vec<SObj> = vec![];
        for _ in 0..p.n_min {
            objs.push(SObj { p: ty(&mut rs), h: rs.gen_range(p.hmin..=p.hmax), t: rs.gen_range(0..=p.tmax) });
        }
        let nmin = objs.len();
        let mut cones = vec![];
        let mut cone_roots = vec![];
        for _ in 0..p.n_cones {
            let o = SObj { p: ty(&mut rs), h: rs.gen_range(p.hmin..p.hmax.max(p.hmin + 1)), t: rs.gen_range(0..=p.tmax) };
            objs.push(o);
            objs.push(SObj { h: o.h + 1, ..o });
            let (a, b) = vals.ab();
            cone_roots.push(a.neg().mul(&b.inv()));
            cones.push((objs.len() - 2, objs.len() - 1, a.add(&b.mul(&x))));
        }
        let n = objs.len();
        let shp = |o: &SObj| shapes.iter().find(|s| s.p == o.p).unwrap();
        let allcos: BTreeSet<u64> = shapes.iter().flat_map(|s| s.cosets.iter().copied()).collect();
        // structure: arrows of D_min and entries of N, with values drawn now (fixed order)
        // arrow values per coset are drawn for every joint coset in ascending order
        let mut dmin: Vec<(usize, usize, BTreeMap<u64, Vec<F>>)> = vec![];
        for i in 0..nmin {
            for j in 0..nmin {
                let (a, b) = (objs[i], objs[j]);
                if a.h.rem_euclid(2) == 0 && b.h == a.h + 1 && a.t > b.t && rs.gen_bool(p.density) {
                    let mut blocks = BTreeMap::new();
                    for &c in &allcos {
                        let (dc, dr) = (shp(&a).dim_at(c), shp(&b).dim_at(c));
                        if dc > 0 && dr > 0 {
                            blocks.insert(c, (0..dc * dr).map(|_| vals.v()).collect());
                        }
                    }
                    dmin.push((i, j, blocks));
                }
            }
        }
        let key = |i: usize| (-objs[i].t, i);
        let mut nmap: Vec<(usize, usize, BTreeMap<u64, Vec<F>>)> = vec![];
        if p.conj {
            for i in 0..n {
                for j in 0..n {
                    let (a, b) = (objs[i], objs[j]);
                    if a.h != b.h || key(i) >= key(j) {
                        continue;
                    }
                    let scalar = a.t == b.t;
                    if scalar && a.p != b.p {
                        continue;
                    }
                    if !rs.gen_bool(p.density) {
                        continue;
                    }
                    let mut blocks = BTreeMap::new();
                    let s = vals.v();
                    for &c in &allcos {
                        let (dc, dr) = (shp(&a).dim_at(c), shp(&b).dim_at(c));
                        if dc > 0 && dr > 0 {
                            let blk: Vec<F> = if scalar {
                                (0..dr * dc).map(|k| if k / dc == k % dc { s.clone() } else { F::zero() }).collect()
                            } else {
                                (0..dc * dr).map(|_| vals.v()).collect()
                            };
                            blocks.insert(c, blk);
                        }
                    }
                    nmap.push((i, j, blocks));
                }
            }
        }
        // assemble per coset and conjugate
        let mut cx = SCx::new();
        for o in &objs {
            cx.add_obj(*o);
        }
        let mut min = SCx::new();
        for o in &objs[..nmin] {
            min.add_obj(*o);
        }
        let mut comp: BTreeMap<(usize, usize), BTreeMap<u64, Vec<F>>> = BTreeMap::new();
        for &c in &allcos {
            let mut off = vec![usize::MAX; n];
            let mut tot = 0;
            for i in 0..n {
                let d = shp(&objs[i]).dim_at(c) as usize;
                if d > 0 {
                    off[i] = tot;
                    tot += d;
                }
            }
            if tot == 0 {
                continue;
            }
            let dim = |i: usize| shp(&objs[i]).dim_at(c) as usize;
            let mut d0 = vec![vec![F::zero(); tot]; tot];
            for (i, j, bl) in &dmin {
                if let Some(b) = bl.get(&c) {
                    let (dc, dr) = (dim(*i), dim(*j));
                    for r in 0..dr {
                        for k in 0..dc {
                            d0[off[*j] + r][off[*i] + k] = b[r * dc + k].clone();
                        }
                    }
                }
            }
            for (i, j, s) in &cones {
                if off[*i] != usize::MAX {
                    for r in 0..dim(*i) {
                        d0[off[*j] + r][off[*i] + r] = s.clone();
                    }
                }
            }
            let mut pm = vec![vec![F::zero(); tot]; tot];
            for (k, row) in pm.iter_mut().enumerate() {
                row[k] = F::one();
            }
            for (i, j, bl) in &nmap {
                if let Some(b) = bl.get(&c) {
                    let (dc, dr) = (dim(*i), dim(*j));
                    for r in 0..dr {
                        for k in 0..dc {
                            pm[off[*j] + r][off[*i] + k] = b[r * dc + k].clone();
                        }
                    }
                }
            }
            let d = if p.conj { matmul(&matmul(&pm, &d0), &inverse(&pm)) } else { d0.clone() };
            for i in 0..n {
                for j in 0..n {
                    if off[i] == usize::MAX || off[j] == usize::MAX || objs[j].h != objs[i].h + 1 {
                        continue;
                    }
                    let (dc, dr) = (dim(i), dim(j));
                    let blk: Vec<F> = (0..dr * dc).map(|k| d[off[j] + k / dc][off[i] + k % dc].clone()).collect();
                    if blk.iter().any(|v| !v.is_zero()) {
                        comp.entry((i, j)).or_default().insert(c, blk);
                    }
                }
            }
        }
        for ((i, j), bl) in comp {
            let lay = reg.layout(objs[i].p, objs[j].p).unwrap();
            let mut m = LocMap::zeros(lay.clone());
            for (k, &c) in lay.cosets.iter().enumerate() {
                if let Some(b) = bl.get(&c) {
                    m.block_mut(k).clone_from_slice(b);
                }
            }
            cx.add_map(i, j, m);
        }
        for (i, j, bl) in &dmin {
            let lay = reg.layout(objs[*i].p, objs[*j].p).unwrap();
            let mut m = LocMap::zeros(lay.clone());
            for (k, &c) in lay.cosets.iter().enumerate() {
                if let Some(b) = bl.get(&c) {
                    m.block_mut(k).clone_from_slice(b);
                }
            }
            min.add_map(*i, *j, m);
        }
        Synth { shapes, cx, min, cone_roots }
    }

    /// type id of the product (p from C, q from W)
    pub fn ptid(p: Tid, q: Tid) -> Tid {
        assert!(q < 1000 && p < 4_000_000);
        1_000_000 + p * 1000 + q
    }
    /// product shape: cosets (c1, c2) ↦ 64·c1 + c2, dims multiply
    pub fn product_shape(a: &TypeShape, b: &TypeShape) -> TypeShape {
        let (mut cos, mut dims) = (vec![], vec![]);
        for (i, &c1) in a.cosets.iter().enumerate() {
            for (j, &c2) in b.cosets.iter().enumerate() {
                assert!(c2 < 64);
                cos.push(c1 * 64 + c2);
                dims.push(a.dims[i] * b.dims[j]);
            }
        }
        TypeShape::new(ptid(a.p, b.p), cos, dims, a.rep * 64 + b.rep).unwrap()
    }

    /// f ⊗ id (left = true) or id ⊗ g (left = false), scaled by sign, over layout(src, tgt)
    fn kron<F: Field>(lay: Arc<Layout>, m: &LocMap<F>, other_src: &TypeShape, other_tgt: &TypeShape, left: bool, sign: &F) -> LocMap<F> {
        let mut out = LocMap::zeros(lay.clone());
        for (k, &cc) in lay.cosets.iter().enumerate() {
            let (c1, c2) = (cc / 64, cc % 64);
            let (cm, co) = if left { (c1, c2) } else { (c2, c1) };
            let Some(i) = m.lay.index(cm) else { continue };
            let (mr, mc) = (m.lay.rows[i] as usize, m.lay.cols[i] as usize);
            let (ns, nt) = (other_src.dim_at(co) as usize, other_tgt.dim_at(co) as usize);
            if ns == 0 {
                continue;
            }
            assert_eq!(ns, nt);
            let blk = m.block(i);
            let cols = lay.cols[k] as usize;
            let o = out.block_mut(k);
            for r in 0..mr {
                for c in 0..mc {
                    let v = blk[r * mc + c].mul(sign);
                    if v.is_zero() {
                        continue;
                    }
                    for e in 0..ns {
                        let (row, col) = if left { (r * ns + e, c * ns + e) } else { (e * mr + r, e * mc + c) };
                        o[row * cols + col] = v.clone();
                    }
                }
            }
        }
        out
    }

    /// Build prev ⊗ W as one crossing: one step per h-level of W (ascending), d_prev ⊗ id
    /// within a step, (−1)^{h} id ⊗ d_W from step k−1 to step k. Product shapes are registered
    /// on the fly. With `twist = Some((seed, x))` every term is written in its own random basis
    /// (a unitriangular automorphism Q_k of C ⊗ W_k, applied as Q d Q⁻¹ blockwise over terms):
    /// an isomorphic complex whose within-term pivots differ from term to term, as the summand
    /// bases of real decompositions do. Without a twist every within-term pivot of C is
    /// mirrored in the next term, which makes early elimination harmless (see the report).
    pub fn tensor_crossing<F: Field>(eng: &mut PEngine<F>, w: &SCx<F>, twist: Option<(u64, F)>) -> Result<(), SErr> {
        let prev = eng.begin_crossing()?;
        let wh: BTreeSet<i32> = w.objs().values().map(|o| o.h).collect();
        let (h0, h1) = (*wh.iter().next().unwrap(), *wh.iter().last().unwrap());
        let mut full: SCx<F> = SCx::new();
        let mut term: Vec<u32> = vec![];
        let mut place: HashMap<(usize, usize), usize> = HashMap::new();
        for k in h0..=h1 {
            for (&wi, ow) in w.objs().iter().filter(|(_, o)| o.h == k) {
                for (&o, oc) in prev.objs() {
                    let (sc, sw) = (eng.reg.shape(oc.p)?.clone(), eng.reg.shape(ow.p)?.clone());
                    let pt = ptid(oc.p, ow.p);
                    if !eng.reg.contains(pt) {
                        eng.register(product_shape(&sc, &sw))?;
                    }
                    place.insert((o, wi), full.add_obj(SObj { p: pt, h: oc.h + ow.h, t: oc.t + ow.t }));
                    term.push((k - h0) as u32);
                }
            }
        }
        for (&o, outs) in prev.out() {
            for (&o2, f) in outs {
                for (&wi, ow) in w.objs() {
                    let (p1, p2) = (prev.objs()[&o].p, prev.objs()[&o2].p);
                    let lay = eng.layout(ptid(p1, ow.p), ptid(p2, ow.p))?;
                    let sw = eng.reg.shape(ow.p)?.clone();
                    full.add_map(place[&(o, wi)], place[&(o2, wi)], kron(lay, f, &sw, &sw, true, &F::one()));
                }
            }
        }
        for (&wj, outs) in w.out() {
            for (&wi, g) in outs {
                let (pj, pi) = (w.objs()[&wj].p, w.objs()[&wi].p);
                for (&o, oc) in prev.objs() {
                    let lay = eng.layout(ptid(oc.p, pj), ptid(oc.p, pi))?;
                    let so = eng.reg.shape(oc.p)?.clone();
                    let sign = if oc.h.rem_euclid(2) == 0 { F::one() } else { F::one().neg() };
                    full.add_map(place[&(o, wj)], place[&(o, wi)], kron(lay, g, &so, &so, false, &sign));
                }
            }
        }
        drop(prev);
        let full = match twist {
            Some((seed, x)) => twist_terms(&full, &mut eng.reg, &term, seed, x),
            None => full,
        };
        insert_terms(eng, &full, &term)?;
        eng.end_crossing()?;
        Ok(())
    }

    /// Insert a complex whose objects are tagged with terms 0..n into the open crossing, one step
    /// per term: the objects of term k, the arrows inside term k and from term k−1 to term k.
    /// Panics if an arrow skips or descends terms.
    pub fn insert_terms<F: Field>(eng: &mut PEngine<F>, full: &SCx<F>, term: &[u32]) -> Result<(), SErr> {
        let nt = term.iter().copied().max().map_or(0, |m| m + 1);
        let mut ids: HashMap<usize, usize> = HashMap::new();
        for k in 0..nt {
            eng.begin_step()?;
            for (&i, &o) in full.objs() {
                if term[i] == k {
                    ids.insert(i, eng.add_obj(o)?);
                }
            }
            for (&x, outs) in full.out() {
                for (&y, m) in outs {
                    if term[y] != k {
                        continue;
                    }
                    assert!(term[x] == k || term[x] + 1 == k, "insert_terms: arrow {} → {} from term {} to {}", x, y, term[x], k);
                    eng.add_map(ids[&x], ids[&y], m.clone())?;
                }
            }
            eng.end_step()?;
        }
        Ok(())
    }

    /// Q d Q⁻¹ with Q = I + N block diagonal over terms: N between objects of the same term and
    /// h, strictly increasing in (−t, id), degree ≥ 0 (degree 0: scalar between equal types).
    pub fn twist_terms<F: Field>(cx: &SCx<F>, reg: &mut ShapeReg, term: &[u32], seed: u64, x: F) -> SCx<F> {
        let objs: Vec<SObj> = cx.objs().values().copied().collect();
        assert!(cx.objs().keys().enumerate().all(|(i, &k)| i == k), "twist_terms: ids must be 0..n");
        let n = objs.len();
        let mut rs = StdRng::seed_from_u64(seed ^ 0x7715);
        let mut vals = Vals { rng: StdRng::seed_from_u64(seed ^ 0x7716), x };
        let shp: Vec<Arc<TypeShape>> = objs.iter().map(|o| reg.shape(o.p).unwrap().clone()).collect();
        let allcos: BTreeSet<u64> = shp.iter().flat_map(|s| s.cosets.iter().copied()).collect();
        let key = |i: usize| (-objs[i].t, i);
        let mut nmap: Vec<(usize, usize, BTreeMap<u64, Vec<F>>)> = vec![];
        for i in 0..n {
            for j in 0..n {
                let (a, b) = (objs[i], objs[j]);
                if term[i] != term[j] || a.h != b.h || key(i) >= key(j) {
                    continue;
                }
                let scalar = a.t == b.t;
                if scalar && a.p != b.p {
                    continue;
                }
                if !rs.gen_bool(0.5) {
                    continue;
                }
                let s = vals.v();
                let mut blocks = BTreeMap::new();
                for &c in &allcos {
                    let (dc, dr) = (shp[i].dim_at(c), shp[j].dim_at(c));
                    if dc > 0 && dr > 0 {
                        let blk: Vec<F> = if scalar {
                            (0..dr * dc).map(|k| if k / dc == k % dc { s.clone() } else { F::zero() }).collect()
                        } else {
                            (0..dc * dr).map(|_| vals.v()).collect()
                        };
                        blocks.insert(c, blk);
                    }
                }
                nmap.push((i, j, blocks));
            }
        }
        let mut out = SCx::new();
        for o in &objs {
            out.add_obj(*o);
        }
        let mut comp: BTreeMap<(usize, usize), BTreeMap<u64, Vec<F>>> = BTreeMap::new();
        for &c in &allcos {
            let mut off = vec![usize::MAX; n];
            let mut tot = 0;
            for i in 0..n {
                let d = shp[i].dim_at(c) as usize;
                if d > 0 {
                    off[i] = tot;
                    tot += d;
                }
            }
            if tot == 0 {
                continue;
            }
            let dim = |i: usize| shp[i].dim_at(c) as usize;
            let mut d0 = vec![vec![F::zero(); tot]; tot];
            for (&i, outs) in cx.out() {
                for (&j, m) in outs {
                    let Some(k) = m.lay.index(c) else { continue };
                    let b = m.block(k);
                    let (dc, dr) = (dim(i), dim(j));
                    for r in 0..dr {
                        for q in 0..dc {
                            d0[off[j] + r][off[i] + q] = b[r * dc + q].clone();
                        }
                    }
                }
            }
            let mut pm = vec![vec![F::zero(); tot]; tot];
            for (k, row) in pm.iter_mut().enumerate() {
                row[k] = F::one();
            }
            for (i, j, bl) in &nmap {
                if let Some(b) = bl.get(&c) {
                    let (dc, dr) = (dim(*i), dim(*j));
                    for r in 0..dr {
                        for q in 0..dc {
                            pm[off[*j] + r][off[*i] + q] = b[r * dc + q].clone();
                        }
                    }
                }
            }
            let d = matmul(&matmul(&pm, &d0), &inverse(&pm));
            for i in 0..n {
                for j in 0..n {
                    if off[i] == usize::MAX || off[j] == usize::MAX || objs[j].h != objs[i].h + 1 {
                        continue;
                    }
                    let (dc, dr) = (dim(i), dim(j));
                    let blk: Vec<F> = (0..dr * dc).map(|k| d[off[j] + k / dc][off[i] + k % dc].clone()).collect();
                    if blk.iter().any(|v| !v.is_zero()) {
                        comp.entry((i, j)).or_default().insert(c, blk);
                    }
                }
            }
        }
        for ((i, j), bl) in comp {
            let lay = reg.layout(objs[i].p, objs[j].p).unwrap();
            let mut m = LocMap::zeros(lay.clone());
            for (k, &c) in lay.cosets.iter().enumerate() {
                if let Some(b) = bl.get(&c) {
                    m.block_mut(k).clone_from_slice(b);
                }
            }
            out.add_map(i, j, m);
        }
        out
    }

    /// Insert a whole complex as one crossing of one step (then eliminate).
    pub fn single_crossing<F: Field>(eng: &mut PEngine<F>, c: &SCx<F>) -> Result<(), SErr> {
        eng.begin_crossing()?;
        eng.begin_step()?;
        let mut ids = HashMap::new();
        for (&i, &o) in c.objs() {
            ids.insert(i, eng.add_obj(o)?);
        }
        for (&x, outs) in c.out() {
            for (&y, m) in outs {
                eng.add_map(ids[&x], ids[&y], m.clone())?;
            }
        }
        eng.end_step()?;
        eng.end_crossing()?;
        Ok(())
    }

    /// Iso invariants of a complex of this model: for every coset c, h and threshold s, the rank
    /// of d_h restricted to the subcomplex of objects with t ≤ s (d and every isomorphism of the
    /// model preserve this filtration).
    pub fn rank_profile<F: Field>(cx: &SCx<F>, reg: &ShapeReg) -> BTreeMap<(u64, i32, i32), usize> {
        let mut out = BTreeMap::new();
        let mut cos: BTreeSet<u64> = BTreeSet::new();
        for o in cx.objs().values() {
            cos.extend(reg.shape(o.p).unwrap().cosets.iter().copied());
        }
        let hs: BTreeSet<i32> = cx.objs().values().map(|o| o.h).collect();
        for &c in &cos {
            for &h in &hs {
                let src: Vec<usize> = cx.objs().iter().filter(|(_, o)| o.h == h && reg.shape(o.p).unwrap().dim_at(c) > 0).map(|(&i, _)| i).collect();
                let tgt: Vec<usize> = cx.objs().iter().filter(|(_, o)| o.h == h + 1 && reg.shape(o.p).unwrap().dim_at(c) > 0).map(|(&i, _)| i).collect();
                if src.is_empty() || tgt.is_empty() {
                    continue;
                }
                let mut toff = HashMap::new();
                let mut nr = 0usize;
                for &y in &tgt {
                    toff.insert(y, nr);
                    nr += reg.shape(cx.objs()[&y].p).unwrap().dim_at(c) as usize;
                }
                let ts: BTreeSet<i32> = src.iter().map(|i| cx.objs()[i].t).collect();
                for &s in &ts {
                    // columns = sources with t ≤ s; rows of the transposed matrix
                    let mut rows: Vec<Vec<F>> = vec![];
                    for &x in src.iter().filter(|&&x| cx.objs()[&x].t <= s) {
                        let dx = reg.shape(cx.objs()[&x].p).unwrap().dim_at(c) as usize;
                        for k in 0..dx {
                            let mut col = vec![F::zero(); nr];
                            for (&y, m) in &cx.out()[&x] {
                                let Some(&o) = toff.get(&y) else { continue };
                                let Some(i) = m.lay.index(c) else { continue };
                                let b = m.block(i);
                                let cols = m.lay.cols[i] as usize;
                                for r in 0..m.lay.rows[i] as usize {
                                    col[o + r] = b[r * cols + k].clone();
                                }
                            }
                            rows.push(col);
                        }
                    }
                    out.insert((c, h, s), crate::dense::rank(&rows, nr));
                }
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::F31;

    #[test]
    fn scalar_of_detects_non_scalar() {
        let s = TypeShape::new(1, vec![0, 3], vec![1, 2], 0).unwrap();
        let lay = Arc::new(layout_of(&s, &s));
        let mut m = LocMap::<F31>::zeros(lay);
        m.block_mut(0)[0] = F31::from_i64(5);
        m.block_mut(1).clone_from_slice(&[F31::from_i64(5), F31::zero(), F31::zero(), F31::from_i64(5)]);
        assert_eq!(scalar_of(&m, &s), Some(F31::from_i64(5)));
        m.block_mut(1)[1] = F31::one();
        assert_eq!(scalar_of(&m, &s), None);
    }
}
