//! Global (point-independent) maps between ladder expressions as an immutable DAG; evaluation one
//! top coset at a time at a numeric point, memoized per point.
//!
//! * `Arena`: interned expressions (with lazily built packed `FlowIndex`es), map nodes with
//!   their `MapInfo` (ends, normalized degree, certified flag, provenance), the template library
//!   and the recorded closed-form constants. No interior mutability: `Send + Sync`.
//! * `Eval`: one point, one thread; memoizes blocks per (map, coset), localized window maps per
//!   (template or closed form, window root masks) and layouts per (src, tgt).
//! * Closed forms (§3.4.3) are normalized so that each equals the unique exact map (constant 1):
//!   Eta = e(B, A) at S = B, EtaL = e(A, B) at S = A, Trace = Σ_S m(S)/e(S, A∖S).
use crate::coset::{low_mask, Conf, Mask};
use crate::field::Field;
use crate::sexpr::{rung_ends, Expr, FlowIndex};
use crate::tmpl::{TmplId, TmplLib, WinLoc};
use crate::web::Rung;
use std::collections::HashMap;
use std::rc::Rc;

/// C2-7: templates whose source window has at least this many fixed points are localized per
/// coset on the needed columns only (`Eval::window_block`, `TmplLib::localize_cols`)
pub const LOC_COLS_MIN: usize = 48;
use std::sync::Arc;

pub type ExprId = u32;
pub type MapId = u32;

// ------------------------------------------------------------------ blocks and layouts

/// One coset block, dense row-major: rows = tgt flows at c, cols = src flows at c.
#[derive(Clone, Debug, PartialEq)]
pub struct Block<F> {
    pub rows: u32,
    pub cols: u32,
    pub a: Vec<F>,
} // a.len() == rows * cols

impl<F: Field> Block<F> {
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Block { rows: rows as u32, cols: cols as u32, a: vec![F::zero(); rows * cols] }
    }
    pub fn identity(n: usize) -> Self {
        let mut b = Self::zeros(n, n);
        for i in 0..n {
            b.a[i * n + i] = F::one();
        }
        b
    }
    pub fn get(&self, i: usize, j: usize) -> F {
        self.a[i * self.cols as usize + j].clone()
    }
    /// self · o, dimension-checked
    pub fn mul(&self, o: &Block<F>) -> Block<F> {
        assert_eq!(self.cols, o.rows, "Block::mul: {}x{} · {}x{}", self.rows, self.cols, o.rows, o.cols);
        let mut c = Block::zeros(self.rows as usize, o.cols as usize);
        gemm_acc(&mut c.a, &self.a, &o.a, self.rows as usize, self.cols as usize, o.cols as usize, F::one());
        c
    }
    /// tests and oracles only
    pub fn to_rows(&self) -> Vec<Vec<F>> {
        (0..self.rows as usize).map(|i| self.a[i * self.cols as usize..(i + 1) * self.cols as usize].to_vec()).collect()
    }
    pub fn is_zero(&self) -> bool {
        self.a.iter().all(|x| x.is_zero())
    }
}

/// c[m×n] += alpha · a[m×k] · b[k×n] on row-major slices, without allocation. Every block product
/// of M4a-M6a goes through it.
pub fn gemm_acc<F: Field>(c: &mut [F], a: &[F], b: &[F], m: usize, k: usize, n: usize, alpha: F) {
    assert!(c.len() >= m * n && a.len() >= m * k && b.len() >= k * n, "gemm_acc: slice too short");
    if alpha.is_zero() {
        return;
    }
    let one = alpha.is_one();
    for i in 0..m {
        let crow = &mut c[i * n..(i + 1) * n];
        for t in 0..k {
            let x = &a[i * k + t];
            if x.is_zero() {
                continue;
            }
            let s = if one { x.clone() } else { x.mul(&alpha) };
            let brow = &b[t * n..(t + 1) * n];
            for (cj, bj) in crow.iter_mut().zip(brow) {
                if !bj.is_zero() {
                    cj.fma(&s, bj);
                }
            }
        }
    }
}

/// Block shapes of all maps between two flow spaces over the cosets of their joint support. Computed
/// once per (src, tgt): expressions in M4a/M5, (src type, tgt type, top block) in M6a. Shared by every
/// map with these ends.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Layout {
    pub cosets: Vec<u64>,
    pub rows: Vec<u32>,
    pub cols: Vec<u32>,
    pub off: Vec<u32>,
    pub len: u32,
}

impl Layout {
    /// offsets = prefix sums; asserts len < 2^32
    pub fn new(cosets: Vec<u64>, rows: Vec<u32>, cols: Vec<u32>) -> Layout {
        assert!(cosets.len() == rows.len() && rows.len() == cols.len(), "Layout::new: lengths differ");
        assert!(cosets.windows(2).all(|w| w[0] < w[1]), "Layout::new: cosets not strictly ascending");
        let mut off = Vec::with_capacity(cosets.len());
        let mut s: u64 = 0;
        for (r, c) in rows.iter().zip(&cols) {
            off.push(s as u32);
            s += *r as u64 * *c as u64;
            assert!(s < 1u64 << 32, "Layout::new: total length >= 2^32");
        }
        Layout { cosets, rows, cols, off, len: s as u32 }
    }
    /// cosets in both supports; flow counts
    pub fn of(src: &FlowIndex, tgt: &FlowIndex) -> Layout {
        let (a, b) = (src.support(), tgt.support());
        let (mut i, mut j) = (0, 0);
        let (mut cos, mut rows, mut cols) = (vec![], vec![], vec![]);
        while i < a.len() && j < b.len() {
            if a[i] < b[j] {
                i += 1;
            } else if b[j] < a[i] {
                j += 1;
            } else {
                cos.push(a[i]);
                cols.push(src.start[i + 1] - src.start[i]);
                rows.push(tgt.start[j + 1] - tgt.start[j]);
                i += 1;
                j += 1;
            }
        }
        Layout::new(cos, rows, cols)
    }
    /// binary search in `cosets`
    pub fn index(&self, c: u64) -> Option<usize> {
        self.cosets.binary_search(&c).ok()
    }
    pub fn block_len(&self, i: usize) -> usize {
        self.rows[i] as usize * self.cols[i] as usize
    }
}

/// A map localized at one point (binding layout, F10): one contiguous vector over a shared layout.
/// Zero blocks are stored as zeros; their fraction is measured (§6.8).
#[derive(Clone, Debug)]
pub struct LocMap<F> {
    pub lay: Arc<Layout>,
    pub data: Vec<F>,
}

impl<F: Field> LocMap<F> {
    pub fn zeros(lay: Arc<Layout>) -> Self {
        let n = lay.len as usize;
        LocMap { lay, data: vec![F::zero(); n] }
    }
    /// i indexes lay.cosets
    pub fn block(&self, i: usize) -> &[F] {
        let o = self.lay.off[i] as usize;
        &self.data[o..o + self.lay.block_len(i)]
    }
    pub fn block_mut(&mut self, i: usize) -> &mut [F] {
        let o = self.lay.off[i] as usize;
        let l = self.lay.block_len(i);
        &mut self.data[o..o + l]
    }
    pub fn is_zero(&self) -> bool {
        self.data.iter().all(|x| x.is_zero())
    }
    pub fn nonzero_blocks(&self) -> usize {
        (0..self.lay.cosets.len()).filter(|&i| self.block(i).iter().any(|x| !x.is_zero())).count()
    }
    /// 8·len + O(1); the shared layout is not counted
    pub fn bytes(&self) -> usize {
        std::mem::size_of::<F>() * self.data.capacity() + std::mem::size_of::<Self>()
    }
}

// ------------------------------------------------------------------ points

#[derive(Clone, Debug)]
pub struct Point<F> {
    pub id: u32,
    pub x: Vec<F>,
}

/// The 32-byte RNG seed of point `id` of a point set: (seed, id) written into disjoint bytes,
/// so distinct (seed, id) pairs never share a generator (no XOR aliasing; M4a review finding g).
pub fn point_seed(seed: u64, id: u32) -> [u8; 32] {
    let mut s = [0u8; 32];
    s[..8].copy_from_slice(&seed.to_le_bytes());
    s[8..12].copy_from_slice(&id.to_le_bytes());
    s[12..20].copy_from_slice(b"krhompts");
    s
}

/// `count` points with ids first_id.., each with K pairwise distinct coordinates drawn uniformly
/// from 1..2·10^9 (redrawn on a collision); point `id` depends only on (seed, id), injectively
/// (`point_seed`).
pub fn points<F: Field>(k: usize, seed: u64, count: usize, first_id: u32) -> Vec<Point<F>> {
    use rand::{rngs::StdRng, Rng, SeedableRng};
    (0..count)
        .map(|i| {
            let id = first_id + i as u32;
            let mut rng = StdRng::from_seed(point_seed(seed, id));
            let mut x: Vec<F> = Vec::with_capacity(k);
            while x.len() < k {
                let v = F::from_i64(rng.gen_range(1..2_000_000_000i64));
                if !x.contains(&v) {
                    x.push(v);
                }
            }
            Point { id, x }
        })
        .collect()
}

// ------------------------------------------------------------------ nodes

/// Closed-form elementary maps (R1 §2.4), all on a 2-upright window (u, u+1) starting at rung `at`
/// (Assoc: on the window between the common prefix and the longest common rung suffix).
/// Normalized so that each equals the unique exact map (constant 1, F8, §3.4.3).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum ElemKind {
    /// E2 on window [E_u^(b), F_u^(b)] over labels (a, b) (all of upright u+1 merged and split
    /// back): Eps (window -> empty, degree +ab) and Eta (empty -> window, degree +ab, value e(B,A));
    /// mirrored window [F_u^(a), E_u^(a)] (all of upright u): EpsL, EtaL (value e(A,B))
    Eps,
    Eta,
    EpsL,
    EtaL,
    /// E3 digon [F_u^(b), E_u^(b)] over (c, 0), and mirrored [E_u^(b), F_u^(b)] over (0, c):
    /// Iota (empty -> digon, degree -ab), Trace (digon -> empty, degree -ab, Σ_S m(S)/e(S, A∖S))
    Iota,
    Trace,
    IotaL,
    TraceL,
    /// E4: identity by matching flows with equal intermediate configurations outside the
    /// window and equal top configuration; only for window pairs whose flows are in bijection
    /// with their top configurations (asserted). Labels key (K_w, w).
    Assoc,
}

#[derive(Clone, Debug)]
pub enum MapNode<F> {
    Id(ExprId),
    /// E1: multiplication by s_nu (Schur polynomial; nu = partition, empty = 1) of the roots on
    /// upright `upright` after `level` rungs; degree 2|nu|
    Dot { e: ExprId, level: usize, upright: usize, nu: Vec<u8> },
    Elem { kind: ElemKind, src: ExprId, tgt: ExprId, at: usize, upright: usize, scalar: F },
    /// id_prefix (x) template (x) id_suffix: src = P·A·Q, tgt = P·A'·Q, A at rungs at..at+|A|
    Local { src: ExprId, tgt: ExprId, at: usize, first: usize, tmpl: TmplId },
    Comp(MapId, MapId), // Comp(g, f) = g ∘ f
    Lin(Vec<(F, MapId)>), // constant coefficients only
    Kron(MapId, Vec<Rung>), // f (x) id_rungs
}

#[derive(Clone, Debug)]
pub struct MapInfo {
    pub src: ExprId,
    pub tgt: ExprId,
    pub deg: i32,
    pub certified: bool,
    pub provenance: String,
}

/// window record of a Local or Elem node
#[derive(Clone, Copy, Debug)]
struct WinRec {
    at: u32,
    first: u32,
    w: u32,
    la: u32,
    lb: u32,
    ws: u32,
    wt: u32,
    elem: Option<(ElemKind, (u8, u8))>,
}

/// window expressions, first upright, labels key and normalized degree of a closed form
struct ElemSpec {
    ws: Expr,
    wt: Expr,
    first: usize,
    key: (u8, u8),
    deg: i32,
}

fn shift_rungs(rs: &[Rung], by: isize) -> Vec<Rung> {
    rs.iter().map(|r| Rung { i: (r.i as isize + by) as usize, dir: r.dir, r: r.r }).collect()
}

/// window of a closed form at (at, u) of src -> tgt
fn elem_spec(kind: ElemKind, src: &Expr, tgt: &Expr, at: usize, u: usize) -> Result<ElemSpec, String> {
    use ElemKind::*;
    let n = src.bottom.len();
    let name = format!("{:?}", kind);
    if at > src.rungs.len() || at > tgt.rungs.len() {
        return Err(format!("{}: at = {} beyond the expressions", name, at));
    }
    if kind == Assoc {
        if src.rungs[..at] != tgt.rungs[..at] {
            return Err(format!("Assoc: prefixes of src and tgt differ below rung {}", at));
        }
        let (sa, ta) = (&src.rungs[at..], &tgt.rungs[at..]);
        let mut s = 0;
        while s < sa.len() && s < ta.len() && sa[sa.len() - 1 - s] == ta[ta.len() - 1 - s] {
            s += 1;
        }
        let (a, b) = (&sa[..sa.len() - s], &ta[..ta.len() - s]);
        if a.is_empty() && b.is_empty() {
            return Err("Assoc: src and tgt are equal (empty window)".into());
        }
        let lo = a.iter().chain(b).map(|r| r.i).min().unwrap();
        let hi = a.iter().chain(b).map(|r| r.i + 1).max().unwrap();
        if u > lo {
            return Err(format!("Assoc: upright {} is right of the window's first upright {}", u, lo));
        }
        let w = hi + 1 - u;
        if u + w > n {
            return Err("Assoc: window out of range".into());
        }
        let lv = src.level(at);
        let wl = lv[u..u + w].to_vec();
        let ws = Expr::new(&wl, &shift_rungs(a, -(u as isize))).map_err(|e| format!("Assoc window: {}", e))?;
        let wt = Expr::new(&wl, &shift_rungs(b, -(u as isize))).map_err(|e| format!("Assoc window: {}", e))?;
        let kw: usize = wl.iter().map(|&x| x as usize).sum();
        return Ok(ElemSpec { ws, wt, first: u, key: (kw as u8, w as u8), deg: 0 });
    }
    if u + 1 >= n {
        return Err(format!("{}: upright {} has no right neighbour", name, u));
    }
    let lv = src.level(at);
    let (la, lb) = (lv[u], lv[u + 1]);
    // the 2-rung window and the side (src or tgt) that carries it
    let (rungs, on_src, key, deg): (Vec<Rung>, bool, (u8, u8), i32) = match kind {
        Eps | Eta => {
            if lb == 0 {
                return Err(format!("{}: label of upright {} is 0", name, u + 1));
            }
            let b = lb as usize;
            (vec![Rung::e(0, b), Rung::f(0, b)], kind == Eps, (la, lb), la as i32 * lb as i32)
        }
        EpsL | EtaL => {
            if la == 0 {
                return Err(format!("{}: label of upright {} is 0", name, u));
            }
            let a = la as usize;
            (vec![Rung::f(0, a), Rung::e(0, a)], kind == EpsL, (la, lb), la as i32 * lb as i32)
        }
        Iota | Trace | IotaL | TraceL => {
            let mirror = matches!(kind, IotaL | TraceL);
            let on_src = matches!(kind, Trace | TraceL);
            let side = if on_src { src } else { tgt };
            let (c, zero) = if mirror { (lb, la) } else { (la, lb) };
            if zero != 0 {
                return Err(format!("{}: digon needs labels {} (got {:?})", name, if mirror { "(0, c)" } else { "(c, 0)" }, (la, lb)));
            }
            let r0 = side.rungs.get(at).ok_or_else(|| format!("{}: no rung at {}", name, at))?;
            let b = r0.r;
            if b == 0 || b > c as usize {
                return Err(format!("{}: digon rung {:?} invalid for c = {}", name, r0, c));
            }
            let rs = if mirror { vec![Rung::e(0, b), Rung::f(0, b)] } else { vec![Rung::f(0, b), Rung::e(0, b)] };
            let a = c as i32 - b as i32;
            (rs, on_src, (c, b as u8), -(a * b as i32))
        }
        Assoc => unreachable!(),
    };
    let wl = vec![la, lb];
    let wwin = Expr::new(&wl, &rungs).map_err(|e| format!("{} window: {}", name, e))?;
    let wempty = Expr::new(&wl, &[]).unwrap();
    let (ws, wt) = if on_src { (wwin, wempty) } else { (wempty, wwin) };
    Ok(ElemSpec { ws, wt, first: u, key, deg })
}

/// src = P·A·Q, tgt = P·A'·Q with A = ws (shifted by `first`) at rungs at..; window bottom labels
/// match; window tops agree
fn check_window(src: &Expr, tgt: &Expr, at: usize, first: usize, ws: &Expr, wt: &Expr) -> Result<(), String> {
    let w = ws.bottom.len();
    if src.bottom != tgt.bottom {
        return Err(format!("window map: bottoms differ ({:?} vs {:?})", src.bottom, tgt.bottom));
    }
    if ws.bottom != wt.bottom {
        return Err("window map: window bottoms differ".into());
    }
    if first + w > src.bottom.len() {
        return Err(format!("window map: window uprights {}..{} out of range", first, first + w));
    }
    let (la, lb) = (ws.rungs.len(), wt.rungs.len());
    if at + la > src.rungs.len() || at + lb > tgt.rungs.len() {
        return Err(format!("window map: window at rung {} longer than src/tgt", at));
    }
    if src.rungs[..at] != tgt.rungs[..at] {
        return Err(format!("window map: prefixes of src and tgt differ below rung {}", at));
    }
    let lv = src.level(at);
    if lv[first..first + w] != ws.bottom[..] {
        return Err(format!("window map: labels {:?} at rung {} on uprights {}..{} != window labels {:?}", &lv[first..first + w], at, first, first + w, ws.bottom));
    }
    if shift_rungs(&ws.rungs, first as isize) != src.rungs[at..at + la] {
        return Err(format!("window map: src rungs {:?} at {} do not match the window {:?} (first upright {})", &src.rungs[at..at + la], at, ws.rungs, first));
    }
    if shift_rungs(&wt.rungs, first as isize) != tgt.rungs[at..at + lb] {
        return Err(format!("window map: tgt rungs {:?} at {} do not match the window {:?} (first upright {})", &tgt.rungs[at..at + lb], at, wt.rungs, first));
    }
    if src.rungs[at + la..] != tgt.rungs[at + lb..] {
        return Err("window map: suffixes of src and tgt differ".into());
    }
    if ws.top() != wt.top() {
        return Err("window map: window tops differ".into());
    }
    Ok(())
}

/// e(S, T)(y) = Π_{s ∈ S, t ∈ T} (y_s − y_t)
pub fn e_pair<F: Field>(s: Mask, t: Mask, y: &[F]) -> F {
    let mut v = F::one();
    let mut a = s;
    while a != 0 {
        let i = a.trailing_zeros() as usize;
        a &= a - 1;
        let mut b = t;
        while b != 0 {
            let j = b.trailing_zeros() as usize;
            b &= b - 1;
            v = v.mul(&y[i].sub(&y[j]));
        }
    }
    v
}

/// Closed-form window map (constant 1) at window roots y, CSC over the sorted window flows.
pub fn elem_winloc<F: Field>(kind: ElemKind, ws: &FlowIndex, wt: &FlowIndex, y: &[F]) -> WinLoc<F> {
    use ElemKind::*;
    let (ns, nt) = (ws.nflows(), wt.nflows());
    let lab = &ws.expr.bottom;
    let mut cols: Vec<Vec<(u32, F)>> = vec![vec![]; ns];
    match kind {
        Eps | EpsL | Eta | EtaL => {
            let (a, b) = (lab[0] as usize, lab[1] as usize);
            let am = low_mask(a);
            let bm = low_mask(b) << a;
            let s = if matches!(kind, Eps | Eta) { bm } else { am };
            match kind {
                Eps | EpsL => {
                    for j in 0..ns as u32 {
                        if ws.flow(j)[1] == s {
                            cols[j as usize].push((0, F::one()));
                        }
                    }
                }
                _ => {
                    let v = if kind == Eta { e_pair(bm, am, y) } else { e_pair(am, bm, y) };
                    for i in 0..nt as u32 {
                        if wt.flow(i)[1] == s {
                            cols[0].push((i, v.clone()));
                        }
                    }
                }
            }
        }
        Iota | IotaL => {
            for i in 0..nt as u32 {
                cols[0].push((i, F::one()));
            }
        }
        Trace | TraceL => {
            let c = if kind == Trace { lab[0] } else { lab[1] } as usize;
            let all = low_mask(c); // the digon's roots are the window roots 0..c in both cases
            for j in 0..ns as u32 {
                let s = ws.flow(j)[0];
                cols[j as usize].push((0, e_pair(s, all & !s, y).inv()));
            }
        }
        Assoc => {
            for j in 0..ns as u32 {
                let tc = ws.top_conf(ws.flow(j));
                let r = wt.range(wt.top.rank(&tc));
                assert_eq!(r.len(), 1, "Assoc: top configuration without a unique target flow");
                cols[j as usize].push((r.start, F::one()));
            }
        }
    }
    let mut col_start = vec![0u32];
    let (mut rows, mut vals) = (vec![], vec![]);
    for c in cols {
        for (i, v) in c {
            rows.push(i);
            vals.push(v);
        }
        col_start.push(rows.len() as u32);
    }
    WinLoc { nrows: nt as u32, col_start, rows, vals }
}

// ------------------------------------------------------------------ the arena

/// Send + Sync (A16): no Rc, RefCell or Cell inside; immutable once frozen into DecData (M5).
pub struct Arena<F: Field> {
    lam: Vec<u8>,
    k: usize,
    exprs: Vec<Expr>,
    expr_ids: HashMap<Expr, ExprId>,
    idx: Vec<Option<FlowIndex>>,
    idx_gen: Vec<u64>,
    stamp: u64,
    nodes: Vec<MapNode<F>>,
    infos: Vec<MapInfo>,
    wrec: Vec<Option<WinRec>>,
    wins: Vec<FlowIndex>,
    win_ids: HashMap<Expr, u32>,
    lib: TmplLib<F>,
    consts: HashMap<(ElemKind, (u8, u8)), F>,
    /// Assoc constants per window pair (src window, tgt window) (M4a review finding e)
    assoc: HashMap<(Expr, Expr), F>,
}

#[derive(Clone, Copy, Debug)]
pub struct ArenaMark {
    pub nodes: u32,
    pub exprs: u32,
    pub stamp: u64,
    pub indices: u32,
}

#[derive(Clone, Debug, Default)]
pub struct ArenaBytes {
    pub exprs: usize,
    pub indices: usize,
    pub indexed_exprs: usize,
    pub nodes: usize,
    pub templates: usize,
}

impl ArenaBytes {
    pub fn total(&self) -> usize {
        self.exprs + self.indices + self.nodes + self.templates
    }
}

fn children<F>(n: &MapNode<F>) -> Vec<MapId> {
    match n {
        MapNode::Comp(g, f) => vec![*g, *f],
        MapNode::Lin(t) => t.iter().map(|x| x.1).collect(),
        MapNode::Kron(f, _) => vec![*f],
        _ => vec![],
    }
}

impl<F: Field> Arena<F> {
    pub fn new(lam: &[u8], lib: TmplLib<F>) -> Self {
        let k: usize = lam.iter().map(|&a| a as usize).sum();
        assert!(k <= 64, "Arena: K = {} > 64", k);
        Arena {
            lam: lam.to_vec(),
            k,
            exprs: vec![],
            expr_ids: HashMap::new(),
            idx: vec![],
            idx_gen: vec![],
            stamp: 0,
            nodes: vec![],
            infos: vec![],
            wrec: vec![],
            wins: vec![],
            win_ids: HashMap::new(),
            lib,
            consts: HashMap::new(),
            assoc: HashMap::new(),
        }
    }
    pub fn lam(&self) -> &[u8] {
        &self.lam
    }
    pub fn k(&self) -> usize {
        self.k
    }
    /// interns only; no FlowIndex is built (panics unless the bottom is lambda)
    pub fn expr(&mut self, e: Expr) -> ExprId {
        assert_eq!(e.bottom, self.lam, "Arena::expr: bottom {:?} != lambda {:?}", e.bottom, self.lam);
        if let Some(&id) = self.expr_ids.get(&e) {
            return id;
        }
        let id = self.exprs.len() as ExprId;
        self.expr_ids.insert(e.clone(), id);
        self.exprs.push(e);
        self.idx.push(None);
        self.idx_gen.push(0);
        id
    }
    /// builds the FlowIndex on first call
    pub fn index(&mut self, id: ExprId) -> &FlowIndex {
        let i = id as usize;
        if self.idx[i].is_none() {
            self.idx[i] = Some(FlowIndex::build(&self.exprs[i]));
            self.stamp += 1;
            self.idx_gen[i] = self.stamp;
        }
        self.idx[i].as_ref().unwrap()
    }
    /// indices of every expression reachable from m
    pub fn ensure_indices(&mut self, m: MapId) {
        let mut stack = vec![m];
        let mut seen = std::collections::HashSet::new();
        while let Some(x) = stack.pop() {
            if !seen.insert(x) {
                continue;
            }
            let (s, t) = (self.infos[x as usize].src, self.infos[x as usize].tgt);
            self.index(s);
            self.index(t);
            stack.extend(children(&self.nodes[x as usize]));
        }
    }
    /// frees it; rebuilt identically on demand
    pub fn drop_index(&mut self, id: ExprId) {
        self.idx[id as usize] = None;
    }
    pub fn has_index(&self, id: ExprId) -> bool {
        self.idx[id as usize].is_some()
    }
    pub fn expr_of(&self, id: ExprId) -> &Expr {
        &self.exprs[id as usize]
    }
    pub fn n_exprs(&self) -> usize {
        self.exprs.len()
    }
    pub fn n_nodes(&self) -> usize {
        self.nodes.len()
    }
    /// panics "index not built" (Eval relies on ensure_indices)
    pub fn flows(&self, id: ExprId) -> &FlowIndex {
        self.idx[id as usize].as_ref().unwrap_or_else(|| panic!("index not built for expression {}", id))
    }
    pub fn lib(&mut self) -> &mut TmplLib<F> {
        &mut self.lib
    }
    pub fn lib_ref(&self) -> &TmplLib<F> {
        &self.lib
    }
    /// gives the template library back (to reuse it in another arena)
    pub fn into_lib(self) -> TmplLib<F> {
        self.lib
    }
    pub fn node(&self, m: MapId) -> &MapNode<F> {
        &self.nodes[m as usize]
    }
    fn win_id(&mut self, e: Expr) -> u32 {
        if let Some(&i) = self.win_ids.get(&e) {
            return i;
        }
        let i = self.wins.len() as u32;
        self.wins.push(FlowIndex::build(&e));
        self.win_ids.insert(e, i);
        i
    }
    fn check_expr(&self, e: ExprId) -> Result<(), String> {
        if (e as usize) < self.exprs.len() { Ok(()) } else { Err(format!("unknown expression id {}", e)) }
    }
    fn check_map(&self, m: MapId) -> Result<(), String> {
        if (m as usize) < self.nodes.len() { Ok(()) } else { Err(format!("unknown map id {}", m)) }
    }
    /// checks compatibility (equal bottoms/tops, window match, Lin terms with equal src/tgt/deg),
    /// computes the degree, builds the indices of src and tgt; `certified` = all leaves certified
    /// (Templates are; Elem iff its (kind, labels) is recorded; Dot always; Comp/Lin/Kron: all children)
    pub fn add(&mut self, n: MapNode<F>, provenance: &str) -> Result<MapId, String> {
        let mut wrec: Option<WinRec> = None;
        let (src, tgt, deg, certified) = match &n {
            MapNode::Id(e) => {
                self.check_expr(*e)?;
                (*e, *e, 0, true)
            }
            MapNode::Dot { e, level, upright, nu } => {
                self.check_expr(*e)?;
                let ex = &self.exprs[*e as usize];
                if *level > ex.rungs.len() {
                    return Err(format!("Dot: level {} > {} rungs", level, ex.rungs.len()));
                }
                if *upright >= ex.bottom.len() {
                    return Err(format!("Dot: upright {} out of range", upright));
                }
                if nu.windows(2).any(|w| w[0] < w[1]) {
                    return Err(format!("Dot: {:?} is not a partition", nu));
                }
                let d: i32 = nu.iter().map(|&x| x as i32).sum();
                (*e, *e, 2 * d, true)
            }
            MapNode::Elem { kind, src, tgt, at, upright, scalar } => {
                self.check_expr(*src)?;
                self.check_expr(*tgt)?;
                if scalar.is_zero() {
                    return Err(format!("Elem {:?}: zero scalar", kind));
                }
                let (se, te) = (self.exprs[*src as usize].clone(), self.exprs[*tgt as usize].clone());
                let sp = elem_spec(*kind, &se, &te, *at, *upright)?;
                check_window(&se, &te, *at, sp.first, &sp.ws, &sp.wt)?;
                let (la, lb, w) = (sp.ws.rungs.len(), sp.wt.rungs.len(), sp.ws.bottom.len());
                let ws = self.win_id(sp.ws);
                let wt = self.win_id(sp.wt);
                if *kind == ElemKind::Assoc {
                    let (a, b) = (&self.wins[ws as usize], &self.wins[wt as usize]);
                    if a.support() != b.support() || a.max_block() != 1 || b.max_block() != 1 {
                        return Err(format!(
                            "Assoc: window flows are not in bijection with equal top configurations ({:?} -> {:?})",
                            a.expr.rungs, b.expr.rungs
                        ));
                    }
                    // degree 0 needs equal normalization shifts (M4a review finding e)
                    if a.expr.shift() != b.expr.shift() {
                        return Err(format!("Assoc: window shifts differ ({} vs {}) for {:?} -> {:?}", a.expr.shift(), b.expr.shift(), a.expr.rungs, b.expr.rungs));
                    }
                }
                let rec = WinRec { at: *at as u32, first: sp.first as u32, w: w as u32, la: la as u32, lb: lb as u32, ws, wt, elem: Some((*kind, sp.key)) };
                let cert = self.elem_rec_certified(&rec);
                wrec = Some(rec);
                (*src, *tgt, sp.deg, cert)
            }
            MapNode::Local { src, tgt, at, first, tmpl } => {
                self.check_expr(*src)?;
                self.check_expr(*tgt)?;
                if *tmpl as usize >= self.lib.len() {
                    return Err(format!("Local: unknown template {}", tmpl));
                }
                let t = self.lib.get(*tmpl);
                let (ws, wt) = (t.src.expr(), t.tgt.expr());
                let deg = t.deg;
                let (se, te) = (&self.exprs[*src as usize], &self.exprs[*tgt as usize]);
                check_window(se, te, *at, *first, &ws, &wt)?;
                let (la, lb, w) = (ws.rungs.len(), wt.rungs.len(), ws.bottom.len());
                let wsi = self.win_id(ws);
                let wti = self.win_id(wt);
                wrec = Some(WinRec { at: *at as u32, first: *first as u32, w: w as u32, la: la as u32, lb: lb as u32, ws: wsi, wt: wti, elem: None });
                (*src, *tgt, deg, self.lib.is_certified(*tmpl))
            }
            MapNode::Comp(g, f) => {
                self.check_map(*g)?;
                self.check_map(*f)?;
                let (ig, i_f) = (&self.infos[*g as usize], &self.infos[*f as usize]);
                if i_f.tgt != ig.src {
                    return Err(format!("Comp: tgt of f (expr {}) != src of g (expr {})", i_f.tgt, ig.src));
                }
                (i_f.src, ig.tgt, ig.deg + i_f.deg, ig.certified && i_f.certified)
            }
            MapNode::Lin(terms) => {
                if terms.is_empty() {
                    return Err("Lin: no terms".into());
                }
                for (_, m) in terms {
                    self.check_map(*m)?;
                }
                let i0 = &self.infos[terms[0].1 as usize];
                let mut cert = true;
                for (_, m) in terms {
                    let i = &self.infos[*m as usize];
                    if i.src != i0.src || i.tgt != i0.tgt {
                        return Err(format!("Lin: terms with different ends ({}->{} vs {}->{})", i0.src, i0.tgt, i.src, i.tgt));
                    }
                    if i.deg != i0.deg {
                        return Err(format!("Lin: terms with mixed degrees ({} vs {})", i0.deg, i.deg));
                    }
                    cert &= i.certified;
                }
                (i0.src, i0.tgt, i0.deg, cert)
            }
            MapNode::Kron(f, rs) => {
                self.check_map(*f)?;
                let i_f = self.infos[*f as usize].clone();
                let se = self.exprs[i_f.src as usize].concat(rs).map_err(|e| format!("Kron: {}", e))?;
                let te = self.exprs[i_f.tgt as usize].concat(rs).map_err(|e| format!("Kron: {}", e))?;
                let s = self.expr(se);
                let t = self.expr(te);
                (s, t, i_f.deg, i_f.certified)
            }
        };
        {
            let (se, te) = (&self.exprs[src as usize], &self.exprs[tgt as usize]);
            if se.bottom != te.bottom || se.top() != te.top() {
                return Err(format!("map ends with different bottoms or tops: {:?} -> {:?}", se, te));
            }
        }
        self.index(src);
        self.index(tgt);
        self.nodes.push(n);
        self.infos.push(MapInfo { src, tgt, deg, certified, provenance: provenance.to_string() });
        self.wrec.push(wrec);
        Ok((self.nodes.len() - 1) as MapId)
    }
    /// `MapInfo::certified` is kept consistent with `assert_genuine` (M4a review finding d): it is
    /// computed at `add` and recomputed for every node whenever a constant or a template is
    /// certified later (`certify_elem`, `certify_assoc`, `certify_template`). Certification only
    /// turns leaves from uncertified to certified, so a template certified directly through
    /// `lib()` (bypassing `certify_template`) can only leave the flag conservatively false.
    /// Decisions about genuineness use `assert_genuine`, which walks the DAG.
    pub fn info(&self, m: MapId) -> &MapInfo {
        &self.infos[m as usize]
    }
    /// is the closed form of a window record certified (E2/E3 by (kind, labels); Assoc per window)?
    fn elem_rec_certified(&self, w: &WinRec) -> bool {
        match w.elem {
            Some((ElemKind::Assoc, _)) => self.assoc.contains_key(&(self.wins[w.ws as usize].expr.clone(), self.wins[w.wt as usize].expr.clone())),
            Some((kind, lab)) => self.consts.contains_key(&(kind, lab)),
            None => true,
        }
    }
    /// recompute every node's certified flag (children precede parents)
    fn refresh_certified(&mut self) {
        for m in 0..self.nodes.len() {
            let c = match &self.nodes[m] {
                MapNode::Id(_) | MapNode::Dot { .. } => true,
                MapNode::Elem { .. } => self.wrec[m].map(|w| self.elem_rec_certified(&w)).unwrap_or(false),
                MapNode::Local { tmpl, .. } => self.lib.is_certified(*tmpl),
                MapNode::Comp(g, f) => self.infos[*g as usize].certified && self.infos[*f as usize].certified,
                MapNode::Lin(t) => t.iter().all(|(_, x)| self.infos[*x as usize].certified),
                MapNode::Kron(f, _) => self.infos[*f as usize].certified,
            };
            self.infos[m].certified = c;
        }
    }
    /// records a validated constant (A7). Assoc constants are per window: use `certify_assoc`.
    pub fn certify_elem(&mut self, kind: ElemKind, labels: (u8, u8), scalar: F) {
        assert!(kind != ElemKind::Assoc, "certify_elem: Assoc constants are recorded per window pair (certify_assoc)");
        self.consts.insert((kind, labels), scalar);
        self.refresh_certified();
    }
    pub fn elem_certified(&self, kind: ElemKind, labels: (u8, u8)) -> Option<F> {
        self.consts.get(&(kind, labels)).cloned()
    }
    /// records a validated Assoc constant for one window pair (window expressions as in
    /// `elem_windows`) (A8)
    pub fn certify_assoc(&mut self, ws: &Expr, wt: &Expr, scalar: F) {
        self.assoc.insert((ws.clone(), wt.clone()), scalar);
        self.refresh_certified();
    }
    pub fn assoc_certified(&self, ws: &Expr, wt: &Expr) -> Option<F> {
        self.assoc.get(&(ws.clone(), wt.clone())).cloned()
    }
    /// certify an added template of the library (`TmplLib::certify_added`) and refresh the flags
    pub fn certify_template(&mut self, t: TmplId, seed: u64) -> Result<crate::tmpl::CertReport, String> {
        let r = self.lib.certify_added(t, seed)?;
        self.refresh_certified();
        Ok(r)
    }
    /// (kind, labels key) of an Elem node (Assoc: (K_w, w), informational only)
    pub fn elem_key(&self, m: MapId) -> Option<(ElemKind, (u8, u8))> {
        self.wrec[m as usize].and_then(|w| w.elem)
    }
    /// window expressions (src window, tgt window) of a Local or Elem node
    pub fn elem_windows(&self, m: MapId) -> Option<(&Expr, &Expr)> {
        self.wrec[m as usize].map(|w| (&self.wins[w.ws as usize].expr, &self.wins[w.wt as usize].expr))
    }
    /// M4a form of M4b's assert_genuine: Err naming the first leaf that is not certified (checked
    /// now, not at add time): an Elem whose (kind, labels) (Assoc: window pair) has no recorded
    /// constant, or a Local whose template is an uncertified added template
    pub fn assert_genuine(&self, m: MapId) -> Result<(), String> {
        let mut stack = vec![m];
        let mut seen = std::collections::HashSet::new();
        while let Some(x) = stack.pop() {
            if !seen.insert(x) {
                continue;
            }
            match &self.nodes[x as usize] {
                MapNode::Elem { .. } => {
                    let w = self.wrec[x as usize].ok_or_else(|| format!("map {}: Elem node {} without window record", m, x))?;
                    if !self.elem_rec_certified(&w) {
                        let (kind, lab) = w.elem.unwrap();
                        if kind == ElemKind::Assoc {
                            return Err(format!(
                                "map {} is not genuine: closed form Assoc on the window {:?} -> {:?} has no recorded constant (node {})",
                                m, self.wins[w.ws as usize].expr.rungs, self.wins[w.wt as usize].expr.rungs, x
                            ));
                        }
                        return Err(format!("map {} is not genuine: closed form {:?} with labels {:?} has no recorded constant (node {})", m, kind, lab, x));
                    }
                }
                MapNode::Local { tmpl, .. } => {
                    if !self.lib.is_certified(*tmpl) {
                        return Err(format!("map {} is not genuine: template {} is an uncertified added template (node {})", m, tmpl, x));
                    }
                }
                _ => {}
            }
            stack.extend(children(&self.nodes[x as usize]));
        }
        Ok(())
    }
    /// Scratch regions (used by every provider call, §4.4): `commit` drops everything created after
    /// `mark` (nodes, expressions, indices) except the nodes reachable from `keep`, with their
    /// expressions and indices. It re-interns them and returns their new ids in the order of `keep`.
    /// Ids below the mark are unchanged; templates are always kept.
    pub fn mark(&self) -> ArenaMark {
        ArenaMark { nodes: self.nodes.len() as u32, exprs: self.exprs.len() as u32, stamp: self.stamp, indices: self.idx.iter().filter(|x| x.is_some()).count() as u32 }
    }
    pub fn commit(&mut self, mark: ArenaMark, keep: &[MapId]) -> Vec<MapId> {
        let (n0, e0) = (mark.nodes as usize, mark.exprs as usize);
        let (nn, ne) = (self.nodes.len(), self.exprs.len());
        assert!(n0 <= nn && e0 <= ne, "commit: stale mark");
        let mut keepn = vec![false; nn - n0];
        let mut stack: Vec<MapId> = keep.to_vec();
        while let Some(m) = stack.pop() {
            let mu = m as usize;
            assert!(mu < nn, "commit: unknown map {}", m);
            if mu < n0 || keepn[mu - n0] {
                continue;
            }
            keepn[mu - n0] = true;
            stack.extend(children(&self.nodes[mu]));
        }
        let mut nmap = vec![u32::MAX; nn - n0];
        let mut next = n0 as u32;
        for (i, &k) in keepn.iter().enumerate() {
            if k {
                nmap[i] = next;
                next += 1;
            }
        }
        let mut keepe = vec![false; ne - e0];
        for (i, &k) in keepn.iter().enumerate() {
            if k {
                let inf = &self.infos[n0 + i];
                for e in [inf.src, inf.tgt] {
                    if e as usize >= e0 {
                        keepe[e as usize - e0] = true;
                    }
                }
            }
        }
        let mut emap = vec![u32::MAX; ne - e0];
        let mut nexte = e0 as u32;
        for (i, &k) in keepe.iter().enumerate() {
            if k {
                emap[i] = nexte;
                nexte += 1;
            }
        }
        let me = |e: ExprId| -> ExprId { if (e as usize) < e0 { e } else { emap[e as usize - e0] } };
        let mm = |m: MapId| -> MapId { if (m as usize) < n0 { m } else { nmap[m as usize - n0] } };
        // expressions
        let old_exprs: Vec<Expr> = self.exprs.drain(e0..).collect();
        let old_idx: Vec<Option<FlowIndex>> = self.idx.drain(e0..).collect();
        let old_gen: Vec<u64> = self.idx_gen.drain(e0..).collect();
        for ((i, e), (ix, g)) in old_exprs.into_iter().enumerate().zip(old_idx.into_iter().zip(old_gen)) {
            if keepe[i] {
                self.expr_ids.insert(e.clone(), emap[i]);
                self.exprs.push(e);
                self.idx.push(ix);
                self.idx_gen.push(g);
            } else {
                self.expr_ids.remove(&e);
            }
        }
        // nodes
        let old_nodes: Vec<MapNode<F>> = self.nodes.drain(n0..).collect();
        let old_infos: Vec<MapInfo> = self.infos.drain(n0..).collect();
        let old_wrec: Vec<Option<WinRec>> = self.wrec.drain(n0..).collect();
        for (i, ((nd, mut inf), wr)) in old_nodes.into_iter().zip(old_infos).zip(old_wrec).enumerate() {
            if !keepn[i] {
                continue;
            }
            let nd = match nd {
                MapNode::Id(e) => MapNode::Id(me(e)),
                MapNode::Dot { e, level, upright, nu } => MapNode::Dot { e: me(e), level, upright, nu },
                MapNode::Elem { kind, src, tgt, at, upright, scalar } => MapNode::Elem { kind, src: me(src), tgt: me(tgt), at, upright, scalar },
                MapNode::Local { src, tgt, at, first, tmpl } => MapNode::Local { src: me(src), tgt: me(tgt), at, first, tmpl },
                MapNode::Comp(g, f) => MapNode::Comp(mm(g), mm(f)),
                MapNode::Lin(t) => MapNode::Lin(t.into_iter().map(|(c, m)| (c, mm(m))).collect()),
                MapNode::Kron(f, rs) => MapNode::Kron(mm(f), rs),
            };
            inf.src = me(inf.src);
            inf.tgt = me(inf.tgt);
            self.nodes.push(nd);
            self.infos.push(inf);
            self.wrec.push(wr);
        }
        // indices of old expressions built after the mark and used by no surviving node
        let mut used = vec![false; e0];
        for inf in &self.infos {
            for e in [inf.src, inf.tgt] {
                if (e as usize) < e0 {
                    used[e as usize] = true;
                }
            }
        }
        for e in 0..e0 {
            if self.idx[e].is_some() && self.idx_gen[e] > mark.stamp && !used[e] {
                self.idx[e] = None;
            }
        }
        self.exprs.shrink_to_fit();
        self.idx.shrink_to_fit();
        self.idx_gen.shrink_to_fit();
        self.nodes.shrink_to_fit();
        self.infos.shrink_to_fit();
        self.wrec.shrink_to_fit();
        keep.iter().map(|&m| mm(m)).collect()
    }
    pub fn bytes(&self) -> ArenaBytes {
        let rs = std::mem::size_of::<Rung>();
        let eb = |e: &Expr| e.bottom.capacity() + e.rungs.capacity() * rs + std::mem::size_of::<Expr>();
        let exprs: usize = self.exprs.iter().map(|e| 2 * eb(e) + 8).sum::<usize>() + self.idx.capacity() * std::mem::size_of::<Option<FlowIndex>>() + self.idx_gen.capacity() * 8;
        let indices: usize = self.idx.iter().flatten().map(|f| f.bytes()).sum();
        let indexed_exprs = self.idx.iter().filter(|x| x.is_some()).count();
        let mut nodes = self.nodes.capacity() * std::mem::size_of::<MapNode<F>>() + self.infos.capacity() * std::mem::size_of::<MapInfo>() + self.wrec.capacity() * std::mem::size_of::<Option<WinRec>>();
        for (n, i) in self.nodes.iter().zip(&self.infos) {
            nodes += i.provenance.capacity();
            nodes += match n {
                MapNode::Dot { nu, .. } => nu.capacity(),
                MapNode::Lin(t) => t.capacity() * std::mem::size_of::<(F, MapId)>(),
                MapNode::Kron(_, r) => r.capacity() * rs,
                _ => 0,
            };
        }
        let templates = self.lib.bytes() + self.wins.iter().map(|w| w.bytes()).sum::<usize>();
        ArenaBytes { exprs, indices, indexed_exprs, nodes, templates }
    }
}

// ------------------------------------------------------------------ evaluation

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
enum WKey {
    T(TmplId),
    E(ElemKind, u32, u32),
}

/// s_nu from e_0..e_l of a root set (dual Jacobi-Trudi: det(e_{nu'_i - i + j}))
pub fn schur_from_e<F: Field>(nu: &[u8], e: &[F]) -> F {
    let nu: Vec<usize> = nu.iter().map(|&x| x as usize).filter(|&x| x > 0).collect();
    if nu.is_empty() {
        return F::one();
    }
    let l = nu[0];
    let conj: Vec<usize> = (1..=l).map(|i| nu.iter().filter(|&&x| x >= i).count()).collect();
    let n = conj.len();
    let get = |k: isize| -> F { if k < 0 || k as usize >= e.len() { F::zero() } else { e[k as usize].clone() } };
    let mut m: Vec<Vec<F>> = (0..n).map(|i| (0..n).map(|j| get(conj[i] as isize - i as isize + j as isize)).collect()).collect();
    // determinant by elimination
    let mut det = F::one();
    for c in 0..n {
        let Some(p) = (c..n).find(|&r| !m[r][c].is_zero()) else { return F::zero() };
        if p != c {
            m.swap(p, c);
            det = det.neg();
        }
        det = det.mul(&m[c][c]);
        let inv = m[c][c].inv();
        for r in c + 1..n {
            if m[r][c].is_zero() {
                continue;
            }
            let f = m[r][c].mul(&inv).neg();
            for k in c..n {
                let v = m[c][k].clone();
                m[r][k].fma(&f, &v);
            }
        }
    }
    det
}

fn remap(m: Mask, tab: &[u8; 64]) -> Mask {
    let mut out = 0u64;
    let mut b = m;
    while b != 0 {
        let i = b.trailing_zeros() as usize;
        b &= b - 1;
        debug_assert!(tab[i] != u8::MAX, "remap: root outside the window");
        out |= 1u64 << tab[i];
    }
    out
}

/// Per point and per thread (not Send): memo (MapId, coset) -> Rc<Block>, template/closed-form memo
/// (TmplId or Elem, root ids) -> Rc<WinLoc>, layout memo (src, tgt) -> Arc<Layout>.
pub struct Eval<'a, F: Field> {
    ar: &'a Arena<F>,
    pt: Point<F>,
    blocks: HashMap<(MapId, u64), Rc<Block<F>>>,
    wins: HashMap<(WKey, Vec<Mask>), Rc<WinLoc<F>>>,
    lays: HashMap<(ExprId, ExprId), Arc<Layout>>,
}

impl<'a, F: Field> Eval<'a, F> {
    pub fn new(ar: &'a Arena<F>, pt: Point<F>) -> Self {
        assert_eq!(pt.x.len(), ar.k, "Eval: point has {} coordinates, K = {}", pt.x.len(), ar.k);
        Eval { ar, pt, blocks: HashMap::new(), wins: HashMap::new(), lays: HashMap::new() }
    }
    pub fn point(&self) -> &Point<F> {
        &self.pt
    }
    /// zero block of the right shape if outside support
    pub fn block(&mut self, m: MapId, c: u64) -> Rc<Block<F>> {
        if let Some(b) = self.blocks.get(&(m, c)) {
            return b.clone();
        }
        let ar = self.ar;
        let info = ar.info(m);
        let (sfi, tfi) = (ar.flows(info.src), ar.flows(info.tgt));
        let (rows, cols) = (tfi.range(c).len(), sfi.range(c).len());
        let blk = if rows == 0 || cols == 0 {
            Block::zeros(rows, cols)
        } else {
            match &ar.nodes[m as usize] {
                MapNode::Id(_) => Block::identity(rows),
                MapNode::Dot { level, upright, nu, .. } => {
                    let mut b: Block<F> = Block::zeros(rows, cols);
                    for (j, f) in sfi.range(c).enumerate() {
                        let conf = sfi.conf_at(f, *level);
                        let e = crate::local::esym(conf.0[*upright], &self.pt.x);
                        b.a[j * cols + j] = schur_from_e(nu, &e);
                    }
                    b
                }
                MapNode::Comp(g, f) => {
                    let bf = self.block(*f, c);
                    let bg = self.block(*g, c);
                    assert_eq!(bg.cols, bf.rows, "Comp: block dimensions differ at coset {}", c);
                    let mut b: Block<F> = Block::zeros(rows, cols);
                    gemm_acc(&mut b.a, &bg.a, &bf.a, rows, bf.rows as usize, cols, F::one());
                    b
                }
                MapNode::Lin(terms) => {
                    let mut b: Block<F> = Block::zeros(rows, cols);
                    for (cf, mi) in terms {
                        if cf.is_zero() {
                            continue;
                        }
                        let bi = self.block(*mi, c);
                        for (x, y) in b.a.iter_mut().zip(&bi.a) {
                            if !y.is_zero() {
                                x.fma(cf, y);
                            }
                        }
                    }
                    b
                }
                MapNode::Kron(f, _) => self.kron_block(m, *f, c, rows, cols),
                MapNode::Local { .. } | MapNode::Elem { .. } => self.window_block(m, c, rows, cols),
            }
        };
        let rc = Rc::new(blk);
        self.blocks.insert((m, c), rc.clone());
        rc
    }
    fn kron_block(&mut self, m: MapId, f: MapId, c: u64, rows: usize, cols: usize) -> Block<F> {
        let ar = self.ar;
        let info = ar.info(m);
        let (sfi, tfi) = (ar.flows(info.src), ar.flows(info.tgt));
        let fi = ar.info(f);
        let (fs, ft) = (ar.flows(fi.src), ar.flows(fi.tgt));
        let nf = fs.nr;
        let mut b: Block<F> = Block::zeros(rows, cols);
        let mut tau: Vec<Mask> = vec![0; tfi.nr];
        for (col, g) in sfi.range(c).enumerate() {
            let flow = sfi.flow(g);
            let (sigma, w) = flow.split_at(nf);
            let cp = fs.top.rank(&fs.top_conf(sigma));
            let colp = fs.find(cp, sigma).expect("Kron: source prefix flow missing") as usize;
            let bf = self.block(f, cp);
            let r0 = ft.range(cp).start;
            for rp in 0..bf.rows as usize {
                let v = &bf.a[rp * bf.cols as usize + colp];
                if v.is_zero() {
                    continue;
                }
                let t = ft.flow(r0 + rp as u32);
                tau[..t.len()].copy_from_slice(t);
                tau[t.len()..].copy_from_slice(w);
                let row = tfi.find(c, &tau).expect("Kron: target flow missing") as usize;
                b.a[row * cols + col].add_assign(v);
            }
        }
        b
    }
    fn win_map(&mut self, key: WKey, wmasks: &[Mask], y: &[F], ws: u32, wt: u32) -> Rc<WinLoc<F>> {
        let k = (key, wmasks.to_vec());
        if let Some(w) = self.wins.get(&k) {
            return w.clone();
        }
        let ar = self.ar;
        let wl = match key {
            WKey::T(t) => ar.lib.localize(t, y),
            WKey::E(kind, _, _) => elem_winloc(kind, &ar.wins[ws as usize], &ar.wins[wt as usize], y),
        };
        let rc = Rc::new(wl);
        self.wins.insert(k, rc.clone());
        rc
    }
    fn window_block(&mut self, m: MapId, c: u64, rows: usize, cols: usize) -> Block<F> {
        let ar = self.ar;
        let info = ar.info(m);
        let (sfi, tfi) = (ar.flows(info.src), ar.flows(info.tgt));
        let wr = ar.wrec[m as usize].expect("window record missing");
        let (key, scalar) = match &ar.nodes[m as usize] {
            MapNode::Local { tmpl, .. } => (WKey::T(*tmpl), F::one()),
            MapNode::Elem { kind, scalar, .. } => (WKey::E(*kind, wr.ws, wr.wt), scalar.clone()),
            _ => unreachable!(),
        };
        let (at, first, w, la, lb) = (wr.at as usize, wr.first as usize, wr.w as usize, wr.la as usize, wr.lb as usize);
        let (ws, wt) = (&ar.wins[wr.ws as usize], &ar.wins[wr.wt as usize]);
        let mut b: Block<F> = Block::zeros(rows, cols);
        let mut tau: Vec<Mask> = vec![0; tfi.nr];
        let mut local: Vec<Mask> = vec![0; la];
        let mut l2g: Vec<u8> = Vec::with_capacity(64);
        let mut y: Vec<F> = Vec::with_capacity(64);
        // C2-7: for templates on large windows, localize only the source columns this coset
        // needs (`TmplLib::localize_cols`, one per distinct window configuration) instead of the
        // full n × n inverse per configuration; KRHOM_LOC_FULL=1 restores the full localization
        let mut partial: HashMap<Vec<Mask>, Rc<WinLoc<F>>> = HashMap::new();
        if let WKey::T(t) = key {
            if ar.lib.src_fixed_points(t) >= LOC_COLS_MIN && std::env::var_os("KRHOM_LOC_FULL").is_none() {
                let mut need: HashMap<Vec<Mask>, (Vec<u32>, Vec<F>)> = HashMap::new();
                for g in sfi.range(c) {
                    let sigma = sfi.flow(g);
                    let conf = sfi.conf_at(g, at);
                    let wmasks = &conf.0[first..first + w];
                    if self.wins.contains_key(&(key, wmasks.to_vec())) {
                        continue;
                    }
                    let mut g2l = [u8::MAX; 64];
                    l2g.clear();
                    for &mk in wmasks {
                        let mut bits = mk;
                        while bits != 0 {
                            let r = bits.trailing_zeros() as usize;
                            bits &= bits - 1;
                            g2l[r] = l2g.len() as u8;
                            l2g.push(r as u8);
                        }
                    }
                    for (k, &mk) in sigma[at..at + la].iter().enumerate() {
                        local[k] = remap(mk, &g2l);
                    }
                    let j = ws.locate(&local).unwrap_or_else(|| panic!("window map {}: source window flow {:?} not found", m, local));
                    let e = need.entry(wmasks.to_vec()).or_insert_with(|| (vec![], l2g.iter().map(|&gr| self.pt.x[gr as usize].clone()).collect()));
                    e.0.push(j);
                }
                for (wm, (mut js, yy)) in need {
                    js.sort_unstable();
                    js.dedup();
                    partial.insert(wm, Rc::new(ar.lib.localize_cols(t, &yy, &js)));
                }
            }
        }
        for (col, g) in sfi.range(c).enumerate() {
            let sigma = sfi.flow(g);
            let conf = sfi.conf_at(g, at);
            let wmasks = &conf.0[first..first + w];
            let mut g2l = [u8::MAX; 64];
            l2g.clear();
            for &mk in wmasks {
                let mut bits = mk;
                while bits != 0 {
                    let r = bits.trailing_zeros() as usize;
                    bits &= bits - 1;
                    g2l[r] = l2g.len() as u8;
                    l2g.push(r as u8);
                }
            }
            let mut l2g_tab = [u8::MAX; 64];
            for (l, &gr) in l2g.iter().enumerate() {
                l2g_tab[l] = gr;
            }
            for (k, &mk) in sigma[at..at + la].iter().enumerate() {
                local[k] = remap(mk, &g2l);
            }
            let j = ws.locate(&local).unwrap_or_else(|| panic!("window map {}: source window flow {:?} not found", m, local)) as usize;
            y.clear();
            y.extend(l2g.iter().map(|&gr| self.pt.x[gr as usize].clone()));
            let wl = match partial.get(wmasks) {
                Some(p) => p.clone(),
                None => self.win_map(key, wmasks, &y, wr.ws, wr.wt),
            };
            tau[..at].copy_from_slice(&sigma[..at]);
            tau[at + lb..].copy_from_slice(&sigma[at + la..]);
            for kk in wl.col_start[j] as usize..wl.col_start[j + 1] as usize {
                let i = wl.rows[kk];
                let v = &wl.vals[kk];
                if v.is_zero() {
                    continue;
                }
                let ta = wt.flow(i);
                for (k, &mk) in ta.iter().enumerate() {
                    tau[at + k] = remap(mk, &l2g_tab);
                }
                let row = tfi.find(c, &tau).unwrap_or_else(|| panic!("window map {}: target flow {:?} missing at coset {}", m, tau, c)) as usize;
                let val = if scalar.is_one() { v.clone() } else { v.mul(&scalar) };
                b.a[row * cols + col].add_assign(&val);
            }
        }
        b
    }
    /// Layout::of, memoized
    pub fn layout(&mut self, src: ExprId, tgt: ExprId) -> Arc<Layout> {
        if let Some(l) = self.lays.get(&(src, tgt)) {
            return l.clone();
        }
        let l = Arc::new(Layout::of(self.ar.flows(src), self.ar.flows(tgt)));
        self.lays.insert((src, tgt), l.clone());
        l
    }
    /// over layout(src, tgt)
    pub fn all(&mut self, m: MapId) -> LocMap<F> {
        let info = self.ar.info(m);
        let lay = self.layout(info.src, info.tgt);
        let mut lm = LocMap::zeros(lay.clone());
        for (i, &c) in lay.cosets.iter().enumerate() {
            let b = self.block(m, c);
            lm.block_mut(i).clone_from_slice(&b.a);
        }
        lm
    }
    /// single-map zero test (all cosets)
    pub fn is_zero(&mut self, m: MapId) -> bool {
        let info = self.ar.info(m);
        let lay = self.layout(info.src, info.tgt);
        lay.cosets.iter().all(|&c| self.block(m, c).is_zero())
    }
    /// one flattened row per map
    pub fn stacked(ar: &Arena<F>, ms: &[MapId], pts: &[Point<F>], cosets: &[u64]) -> Vec<Vec<F>> {
        let mut out: Vec<Vec<F>> = vec![vec![]; ms.len()];
        for p in pts {
            let mut ev = Eval::new(ar, p.clone());
            for (r, &m) in ms.iter().enumerate() {
                for &c in cosets {
                    let b = ev.block(m, c);
                    out[r].extend_from_slice(&b.a);
                }
            }
        }
        out
    }
    /// the whole map as CSC over the sorted global flows (rows = tgt flow ids, cols = src flow ids)
    pub fn winloc(&mut self, m: MapId) -> WinLoc<F> {
        let info = self.ar.info(m);
        let (sfi, tfi) = (self.ar.flows(info.src), self.ar.flows(info.tgt));
        let lay = self.layout(info.src, info.tgt);
        let mut cols: Vec<Vec<(u32, F)>> = vec![vec![]; sfi.nflows()];
        for &c in &lay.cosets {
            let b = self.block(m, c);
            let (r0, c0) = (tfi.range(c).start, sfi.range(c).start);
            for i in 0..b.rows as usize {
                for j in 0..b.cols as usize {
                    let v = &b.a[i * b.cols as usize + j];
                    if !v.is_zero() {
                        cols[c0 as usize + j].push((r0 + i as u32, v.clone()));
                    }
                }
            }
        }
        let mut col_start = vec![0u32];
        let (mut rows, mut vals) = (vec![], vec![]);
        for c in cols {
            for (i, v) in c {
                rows.push(i);
                vals.push(v);
            }
            col_start.push(rows.len() as u32);
        }
        WinLoc { nrows: tfi.nflows() as u32, col_start, rows, vals }
    }
    pub fn memo_bytes(&self) -> usize {
        let fs = std::mem::size_of::<F>();
        let b: usize = self.blocks.values().map(|x| x.a.capacity() * fs + 48).sum::<usize>() + self.blocks.len() * 32;
        let w: usize = self.wins.iter().map(|(k, v)| k.1.capacity() * 8 + v.bytes() + 64).sum();
        let l: usize = self.lays.values().map(|l| l.cosets.capacity() * 8 + (l.rows.capacity() + l.cols.capacity() + l.off.capacity()) * 4 + 64).sum();
        b + w + l
    }
    pub fn clear(&mut self) {
        self.blocks.clear();
        self.wins.clear();
        self.lays.clear();
    }
}

// ------------------------------------------------------------------ reduced coordinates

/// Reduced coordinates of an idempotent's image at one coset: P (flows × m), Q (m × flows),
/// Q P = I_m, P Q = E_c (D2's U_c, V_c).
#[derive(Clone, Debug)]
pub struct Reduced<F> {
    pub p: Block<F>,
    pub q: Block<F>,
    pub m: usize,
}

/// asserts E^2 = E
pub fn reduced_of<F: Field>(e_c: &Block<F>) -> Result<Reduced<F>, String> {
    let n = e_c.rows as usize;
    if e_c.cols as usize != n {
        return Err(format!("reduced_of: block {}x{} is not square", e_c.rows, e_c.cols));
    }
    if e_c.mul(e_c) != *e_c {
        return Err("reduced_of: E^2 != E".into());
    }
    let mut r = e_c.to_rows();
    let piv = crate::dense::rref(&mut r, n);
    let m = piv.len();
    let mut p = Block::zeros(n, m);
    for (k, &c) in piv.iter().enumerate() {
        for i in 0..n {
            p.a[i * m + k] = e_c.get(i, c);
        }
    }
    let mut q = Block::zeros(m, n);
    for (k, row) in r.iter().enumerate() {
        for j in 0..n {
            q.a[k * n + j] = row[j].clone();
        }
    }
    if q.mul(&p) != Block::identity(m) || p.mul(&q) != *e_c {
        return Err("reduced_of: factorization check failed".into());
    }
    Ok(Reduced { p, q, m })
}

/// Conf helper: the configuration of a coset rank in a FlowIndex's block
pub fn coset_conf(fi: &FlowIndex, c: u64) -> Conf {
    fi.top.unrank(c)
}

/// (receiver, giver) of a rung (re-exported for callers of this module)
pub fn ends(r: &Rung) -> (usize, usize) {
    rung_ends(r)
}
