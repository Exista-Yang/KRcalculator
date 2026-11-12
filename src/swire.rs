//! M6a-1b: the localized multi-strand engine (sengine) wired to the decomposition library
//! (sindec) and to the box closure (boxclose). SPEC_3C_GL6 rev1 §2 "M6a-1b".
//!
//! * [`WLib`]: a frozen `DecData` plus everything the engine needs that `DecData` does not export
//!   (type shapes with rep(p), table names, the identity type and coset, closure-generator
//!   degrees). Built by [`build_wlib`] from an `SIndec` builder by a breadth-first walk over the
//!   crossing terms of a set of words (every type of every unsimplified stage). All library-side
//!   asserts (ι/π ends and degrees, Rickard ends and degrees, c_diag of the closure generators)
//!   run once here.
//! * [`scan_point`]: one scan of a braid word at one point (record or replay): crossings via
//!   `dec_word` and `rickard_comp_at`; d_C ⊗ id rung word by rung word, blockwise per coset (the
//!   Kron rule on localized blocks); (−1)^h id ⊗ d; K_0 = the decategorified scan and #objects ≥ LB
//!   (optional, point independent); d∘d (`EngOpts::check_d2`, unpruned). The closure entries
//!   (Gram and pairing P = H^y · d_id · G^x) of the requested stages are extracted at the
//!   identity coset.
//! * [`box_word`]: the one-word box loop (record at a generic point, replay at every box point,
//!   `BoxAcc` per stage and mode, `close_box`, Freivalds `D2Cert`), plus [`scaling_check`] (D14:
//!   scans at z and λz).
use crate::boxclose::*;
use crate::cclose::{ClosureDataC, TypeBasisC};
use crate::field::Field;
use crate::laurent::Poincare;
use crate::cclose::EPoly;
use crate::schur::{crossing_terms, Elt, Schur, Table, Tid};
use crate::sengine::*;
use crate::sindec::{DecData, SIndec};
use crate::smap::{gemm_acc, Eval, LocMap, MapId, Point};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::sync::Arc;
use std::time::Instant;

pub const SWIRE_STAMP: &str = "swire-v1";

// ------------------------------------------------------------------ helpers

/// Generator index i and sign of a braid letter ±(i+1).
pub fn letter(g: i32) -> (usize, bool) {
    ((g.unsigned_abs() - 1) as usize, g > 0)
}

/// Kept crossing terms at prune N: the term k of a (c, d) crossing has a web edge of label
/// max(c, d) + k, killed above N (the terms kept are a prefix 0..K).
pub fn kept_terms(c: u8, d: u8, positive: bool, i: usize, prune: Option<u8>) -> Vec<usize> {
    let nt = crossing_terms(c, d, positive, i).len();
    (0..nt).filter(|&k| prune.is_none_or(|n| c.max(d) as usize + k <= n as usize)).collect()
}

// ------------------------------------------------------------------ the wired library

/// Closure data of one closable type: generator maps and degrees.
#[derive(Clone, Debug)]
pub struct WClosure {
    pub ins: Vec<MapId>,
    pub outs: Vec<MapId>,
    pub ins_deg: Vec<i32>,
    pub outs_deg: Vec<i32>,
}

pub struct WLib<F: Field> {
    pub dd: Arc<DecData<F>>,
    pub lam: Vec<u8>,
    pub prune_n: Option<u8>,
    pub shapes: BTreeMap<Tid, TypeShape>,
    pub names: BTreeMap<Tid, String>,
    pub by_name: HashMap<String, Tid>,
    pub id_tid: Tid,
    /// the identity coset of λ (source coset of every closure generator)
    pub id_coset: Option<u64>,
    pub closure: BTreeMap<Tid, WClosure>,
    /// types of the walk without closure generators in the library (error only if final)
    pub no_closure: BTreeMap<Tid, String>,
    pub words: Vec<Vec<i32>>,
    pub hash: u64,
    pub build_secs: f64,
    pub stats: crate::sindec::DecStats,
    /// point-independent Kron flow maps, keyed by (L_a, L_b, X_a, X_b, coset)
    pub kron: std::sync::RwLock<HashMap<(u32, u32, u32, u32, u64), Arc<KronIdx>>>,
}

impl<F: Field> WLib<F> {
    pub fn kron_idx(&self, la: u32, lb: u32, xa: u32, xb: u32, c: u64) -> Arc<KronIdx> {
        let key = (la, lb, xa, xb, c);
        if let Some(k) = self.kron.read().unwrap().get(&key) {
            return k.clone();
        }
        let k = Arc::new(kron_idx(self.dd.ar(), la, lb, xa, xb, c));
        self.kron.write().unwrap().insert(key, k.clone());
        k
    }
}

/// Build the wired library for `words` (every type of every unsimplified stage of every word,
/// every kept term, every Rickard component of adjacent kept terms; closure generators of every
/// type reached at a stage whose labels are λ, when `closure`).
pub fn build_wlib<F: Field>(mut sd: SIndec<F>, words: &[Vec<i32>], closure: bool) -> Result<WLib<F>, String> {
    let t0 = Instant::now();
    let lam = sd.schur.lam.clone();
    let prune = sd.prune_n();
    let id = sd.schur.intern(Table::identity(&lam));
    sd.ty(id).map_err(String::from)?;
    let mut all: BTreeSet<Tid> = BTreeSet::new();
    all.insert(id);
    let mut closable: BTreeSet<Tid> = BTreeSet::new();
    let mut seen_step: BTreeSet<(Tid, usize, bool)> = BTreeSet::new();
    // per-term (dq, summand s) for degree asserts of the Rickard components
    for w in words {
        let mut cur: BTreeSet<Tid> = [id].into_iter().collect();
        let mut labels = lam.clone();
        for &g in w {
            let (i, pos) = letter(g);
            if i + 1 >= labels.len() {
                return Err(format!("letter {} needs {} strands", g, i + 2));
            }
            let (c, d) = (labels[i], labels[i + 1]);
            let terms = crossing_terms(c, d, pos, i);
            let kept = kept_terms(c, d, pos, i, prune);
            let mut next = BTreeSet::new();
            for &p in &cur {
                let cols = sd.schur.table(p).cols();
                if (cols[i], cols[i + 1]) != (c, d) {
                    return Err(format!("type {} has top labels {:?}, expected {:?}", sd.schur.table(p).fmt(), cols, labels));
                }
                let mut sums: BTreeMap<usize, Vec<(Tid, i32)>> = BTreeMap::new();
                for &k in &kept {
                    let dec = sd.dec_word(p, &terms[k].rungs).map_err(|e| format!("dec_word({} ⊗ {}): {}", sd.schur.table(p).fmt(), crate::smapio::rungs_str(&terms[k].rungs), String::from(e)))?;
                    for sm in &dec.summands {
                        next.insert(sm.q);
                    }
                    sums.insert(k, dec.summands.iter().map(|s| (s.q, s.s)).collect());
                }
                if seen_step.insert((p, i, pos)) {
                    for j in 0..terms.len().saturating_sub(1) {
                        let (ks, kt) = if pos { (j, j + 1) } else { (j + 1, j) };
                        if !kept.contains(&ks) || !kept.contains(&kt) {
                            continue;
                        }
                        let comps = sd.rickard_comp_at(p, i, pos, j).map_err(|e| format!("rickard_comp_at({}, {}, {}, {}): {}", sd.schur.table(p).fmt(), i, pos, j, String::from(e)))?;
                        for &(si, ti, m) in comps.iter() {
                            let (qs, ss) = sums[&ks][si];
                            let (qt, st) = sums[&kt][ti];
                            let info = sd.ar.info(m).clone();
                            let (ls, lt) = (sd.ty(qs).map_err(String::from)?.lad, sd.ty(qt).map_err(String::from)?.lad);
                            if info.src != ls || info.tgt != lt {
                                return Err(format!("Rickard component of {} ({}, {}, d_{}): ends are not the canonical ladders", sd.schur.table(p).fmt(), i, pos, j));
                            }
                            let want = (terms[ks].dq + ss) - (terms[kt].dq + st);
                            if info.deg != want {
                                return Err(format!("Rickard component of {} d_{} {}→{}: degree {} != t_src − t_tgt = {}", sd.schur.table(p).fmt(), j, si, ti, info.deg, want));
                            }
                        }
                    }
                }
            }
            labels.swap(i, i + 1);
            cur = next;
            all.extend(cur.iter().copied());
            if labels == lam {
                closable.extend(cur.iter().copied());
            }
        }
    }
    let mut cl: BTreeMap<Tid, WClosure> = BTreeMap::new();
    let mut no_cl = BTreeMap::new();
    if closure {
        for &p in &closable {
            match sd.closure_gens(p) {
                Ok(g) => {
                    let ins_deg: Vec<i32> = g.ins.iter().map(|&m| sd.ar.info(m).deg).collect();
                    let outs_deg: Vec<i32> = g.outs.iter().map(|&m| sd.ar.info(m).deg).collect();
                    // generator degrees = the exponents of c_diag(b_p) (as xengine)
                    let cd = sd.schur.c_diag(p);
                    let mut want: Vec<i32> = vec![];
                    for (k, &c) in cd.c.iter().enumerate() {
                        for _ in 0..c {
                            want.push(cd.lo + k as i32);
                        }
                    }
                    let (mut a, mut b) = (ins_deg.clone(), outs_deg.clone());
                    a.sort();
                    b.sort();
                    want.sort();
                    if a != want || b != want {
                        return Err(format!("type {}: closure generator degrees ins {:?} outs {:?} != c_diag {:?}", sd.schur.table(p).fmt(), a, b, want));
                    }
                    cl.insert(p, WClosure { ins: g.ins.clone(), outs: g.outs.clone(), ins_deg, outs_deg });
                }
                Err(e) => {
                    no_cl.insert(p, String::from(e));
                }
            }
        }
    }
    for &p in &all {
        let lad = sd.ty(p).map_err(String::from)?.lad;
        sd.ar.index(lad);
    }
    for &p in cl.keys() {
        let c = &cl[&p];
        for &m in c.ins.iter().chain(c.outs.iter()) {
            sd.ar.ensure_indices(m);
        }
    }
    let mut reps = BTreeMap::new();
    let mut names = BTreeMap::new();
    for &p in &all {
        reps.insert(p, rep_coset(&sd.schur, &lam, p));
        names.insert(p, sd.schur.table(p).fmt());
    }
    let stats = sd.stats();
    let dd = sd.freeze();
    let mut shapes = BTreeMap::new();
    for &p in &all {
        shapes.insert(p, lib_shape(&*dd, p, reps[&p]).map_err(String::from)?);
    }
    // degree and end asserts of every decomposition (ι: L_q → X of degree s, π: X → L_q of −s)
    for (p, rs) in dd.dec_keys() {
        let dec = dd.dec_word(p, &rs).map_err(String::from)?;
        let ar = dd.ar();
        if !ar.has_index(dec.x) {
            return Err(format!("decomposition of {}: no flow index of X", names[&p]));
        }
        for sm in &dec.summands {
            let lq = dd.ty(sm.q).map_err(String::from)?.lad;
            let (ii, pi) = (ar.info(sm.iota), ar.info(sm.pi));
            if ii.src != lq || ii.tgt != dec.x || pi.src != dec.x || pi.tgt != lq || ii.deg != sm.s || pi.deg != -sm.s {
                return Err(format!("decomposition of {}: summand {} has wrong ends or degrees (ι {} π {}, s {})", names[&p], names.get(&sm.q).cloned().unwrap_or_default(), ii.deg, pi.deg, sm.s));
            }
        }
    }
    let id_coset = match cl.values().next() {
        Some(c) => {
            let g = dd.closure_gens(*cl.keys().next().unwrap()).map_err(String::from)?;
            let _ = c;
            Some(identity_coset(dd.ar(), &g).map_err(String::from)?)
        }
        None => None,
    };
    let mut h = crate::boxclose::Fnv::default();
    h.str(SWIRE_STAMP).u64(F::characteristic()).str(&format!("{:?}{:?}", lam, prune));
    for (p, n) in &names {
        h.u64(*p as u64).str(n);
    }
    for (p, c) in &cl {
        h.u64(*p as u64).str(&format!("{:?}{:?}", c.ins_deg, c.outs_deg));
    }
    for w in words {
        h.str(&format!("{:?}", w));
    }
    let by_name = names.iter().map(|(p, n)| (n.clone(), *p)).collect();
    Ok(WLib { dd, lam, prune_n: prune, shapes, names, by_name, id_tid: id, id_coset, closure: cl, no_closure: no_cl, words: words.to_vec(), hash: h.get(), build_secs: t0.elapsed().as_secs_f64(), stats, kron: Default::default() })
}

// ------------------------------------------------------------------ the decategorified oracle

/// The decategorified scan (K_0 class per stage), in its own Schur instance; types are matched
/// to the library by table.
pub struct Decat {
    s: Schur,
    x: BTreeMap<i32, Elt>,
}

impl Decat {
    pub fn new(lam: &[u8]) -> Decat {
        let mut s = Schur::new(lam);
        let mut x = BTreeMap::new();
        x.insert(0, s.unit());
        Decat { s, x }
    }
    pub fn step(&mut self, c: u8, d: u8, positive: bool, i: usize) {
        let terms = crossing_terms(c, d, positive, i);
        let mut y: BTreeMap<i32, Elt> = BTreeMap::new();
        let x = std::mem::take(&mut self.x);
        for (h, e) in &x {
            for t in &terms {
                let z = self.s.mul_rungs(e, &t.rungs);
                let slot = y.entry(h + t.hdeg).or_insert_with(|| Elt::zero(&z.mu));
                *slot = slot.add_scaled(&z, 1, t.dq);
            }
        }
        y.retain(|_, e| !e.is_zero());
        self.x = y;
    }
    /// signed K_0 class (table, t) → Σ (−1)^h, with the types killed at `prune` removed
    pub fn k0(&mut self, prune: Option<u8>) -> BTreeMap<(String, i32), i64> {
        let mut signed: BTreeMap<(String, i32), i64> = BTreeMap::new();
        let x = self.x.clone();
        for (h, e) in &x {
            let sg: i64 = if h.rem_euclid(2) == 0 { 1 } else { -1 };
            for (p, m) in self.s.expand(e) {
                if prune.is_some_and(|n| self.s.killed_at(p, n as usize)) {
                    continue;
                }
                let name = self.s.table(p).fmt();
                for (k, &c) in m.c.iter().enumerate() {
                    if c != 0 {
                        *signed.entry((name.clone(), m.lo + k as i32)).or_insert(0) += sg * c as i64;
                    }
                }
            }
        }
        signed.retain(|_, c| *c != 0);
        signed
    }
}

// ------------------------------------------------------------------ one scan at one point

#[derive(Clone, Debug)]
pub struct ScanOpts {
    pub eng: EngOpts,
    /// stages (1-based crossing counts) whose closure entries are extracted (labels must be λ)
    pub close_stages: Vec<usize>,
    /// K_0 = decategorified scan and #objects ≥ LB after every crossing (point independent)
    pub k0_check: bool,
}

/// Closure entries of one stage at one point.
#[derive(Clone, Debug)]
pub struct StageOut<F> {
    pub stage: usize,
    /// skeleton: objects, types (degrees; Gram zero), crossings; edges empty
    pub skel: ClosureDataC<F>,
    pub edges: Vec<(usize, usize)>,
    pub kinds: Vec<EntryKind>,
    pub weights: Vec<Option<u32>>,
    pub vals: Vec<F>,
}

#[derive(Clone, Debug)]
pub struct ScanOut<F> {
    pub ps: PointScan<F>,
    pub stages: Vec<StageOut<F>>,
    /// per crossing: (objects, LB) from the K_0 check (empty without it)
    pub lb: Vec<(usize, u64)>,
    pub eval_bytes: usize,
    /// seconds: construction of C ⊗ W (Kron, π K ι, Rickard), engine (insertion + elimination)
    pub t_build: f64,
    pub t_engine: f64,
}

/// Point-independent flow map of f ⊗ id_W at a coset c of X_a = L_a ⊗ W: per source flow g of
/// X_a at c the coset cp and column of its L_a prefix, and for every row rp of f's block at cp the
/// row of X_b = L_b ⊗ W at c.
pub struct KronIdx {
    pub rows: usize,
    pub cols: usize,
    /// per column g: (cp, colp, offset into rowmap); rowmap[off + rp] = row
    pub colinfo: Vec<(u64, u32, u32)>,
    pub rowmap: Vec<u32>,
}

fn kron_idx<F: Field>(ar: &crate::smap::Arena<F>, la: u32, lb: u32, xa: u32, xb: u32, c: u64) -> KronIdx {
    let (fs, ft) = (ar.flows(la), ar.flows(lb));
    let (sfi, tfi) = (ar.flows(xa), ar.flows(xb));
    let nf = fs.nr;
    let (rows, cols) = (tfi.range(c).len(), sfi.range(c).len());
    let mut colinfo = Vec::with_capacity(cols);
    let mut rowmap = vec![];
    let mut tau = vec![0; tfi.nr];
    for g in sfi.range(c) {
        let flow = sfi.flow(g);
        let (sigma, w) = flow.split_at(nf);
        let cp = fs.top.rank(&fs.top_conf(sigma));
        let colp = fs.find(cp, sigma).expect("Kron: source prefix flow missing");
        let off = rowmap.len() as u32;
        for t in ft.range(cp) {
            let tf = ft.flow(t);
            tau[..tf.len()].copy_from_slice(tf);
            tau[tf.len()..].copy_from_slice(w);
            rowmap.push(tfi.find(c, &tau).expect("Kron: target flow missing"));
        }
        colinfo.push((cp, colp, off));
    }
    KronIdx { rows, cols, colinfo, rowmap }
}

/// A library map evaluated on the engine's layout of (src type, tgt type), times `scale`.
fn lib_map<F: Field>(ev: &mut Eval<F>, eng: &mut PEngine<F>, m: MapId, qs: Tid, qt: Tid, scale: &F) -> Result<LocMap<F>, String> {
    let lay = eng.layout(qs, qt).map_err(String::from)?;
    let mut lm = LocMap::zeros(lay.clone());
    for (i, &c) in lay.cosets.iter().enumerate() {
        let b = ev.block(m, c);
        if b.a.len() != lay.block_len(i) {
            return Err(format!("library map block at coset {}: {}x{} vs layout {}x{}", c, b.rows, b.cols, lay.rows[i], lay.cols[i]));
        }
        for (d, s) in lm.block_mut(i).iter_mut().zip(b.a.iter()) {
            *d = s.mul(scale);
        }
    }
    Ok(lm)
}

/// One scan of `word` at the point x (Arena coordinates: the roots of block 0, then block 1, ...).
pub fn scan_point<F: Field>(lib: &WLib<F>, word: &[i32], x: Vec<F>, point_id: u32, mode: Mode, opts: &ScanOpts) -> Result<ScanOut<F>, String> {
    let t0s = Instant::now();
    let dd = &*lib.dd;
    let ar = dd.ar();
    let mut ev = Eval::new(ar, Point { id: point_id, x });
    let mut eng: PEngine<F> = PEngine::new(mode, opts.eng.clone());
    for s in lib.shapes.values() {
        eng.register(s.clone()).map_err(String::from)?;
    }
    let mut cx0 = SCx::new();
    cx0.add_obj(SObj { p: lib.id_tid, h: 0, t: 0 });
    eng.load(cx0).map_err(String::from)?;
    let mut labels = lib.lam.clone();
    let mut decat = if opts.k0_check { Some(Decat::new(&lib.lam)) } else { None };
    let mut lbs = vec![];
    let mut outs = vec![];
    let mut crossings = vec![];
    let (mut t_build, mut t_engine) = (0.0f64, 0.0f64);
    let one = F::one();
    let mone = F::zero().sub(&one);
    for (st, &g) in word.iter().enumerate() {
        let (i, pos) = letter(g);
        let (c, d) = (labels[i], labels[i + 1]);
        let terms = crossing_terms(c, d, pos, i);
        let kept = kept_terms(c, d, pos, i, lib.prune_n);
        let order: Vec<usize> = if pos { kept.clone() } else { kept.iter().rev().copied().collect() };
        let prev = eng.begin_crossing().map_err(String::from)?;
        let mut place: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
        for (si, &k) in order.iter().enumerate() {
            let rs = &terms[k].rungs;
            eng.begin_step().map_err(String::from)?;
            let mut decs = HashMap::new();
            for (&o, ob) in prev.objs() {
                let dec = dd.dec_word(ob.p, rs).map_err(String::from)?;
                let ids = dec.summands.iter().map(|sm| eng.add_obj(SObj { p: sm.q, h: ob.h + terms[k].hdeg, t: ob.t + terms[k].dq + sm.s })).collect::<Result<Vec<_>, _>>().map_err(String::from)?;
                place.insert((o, k), ids);
                decs.insert(ob.p, dec);
            }
            // d_C ⊗ id inside the term
            for (&a, aouts) in prev.out() {
                let pa = prev.objs()[&a].p;
                let da = decs[&pa].clone();
                if da.summands.is_empty() {
                    continue;
                }
                for (&b, f) in aouts {
                    let pb = prev.objs()[&b].p;
                    let db = decs[&pb].clone();
                    if db.summands.is_empty() {
                        continue;
                    }
                    let (la, lb) = (dd.ty(pa).map_err(String::from)?.lad, dd.ty(pb).map_err(String::from)?.lad);
                    // per coset of X_a ∩ X_b: K_c, then π K ι per summand pair
                    let (sa, sb) = (ar.flows(da.x).support().to_vec(), ar.flows(db.x).support().to_vec());
                    let mut per: Vec<Vec<BTreeMap<u64, (usize, usize, Vec<F>)>>> = vec![vec![BTreeMap::new(); db.summands.len()]; da.summands.len()];
                    for &cc in sa.iter().filter(|cc| sb.binary_search(cc).is_ok()) {
                        let ki_ = lib.kron_idx(la, lb, da.x, db.x, cc);
                        let kr = ki_.rows;
                        // nonzero test of the f blocks this coset reads
                        if !ki_.colinfo.iter().any(|&(cp, _, _)| f.lay.index(cp).is_some_and(|i| f.block(i).iter().any(|v| !v.is_zero()))) {
                            continue;
                        }
                        for (ia, sg) in da.summands.iter().enumerate() {
                            let ib = ev.block(sg.iota, cc);
                            if ib.cols == 0 {
                                continue;
                            }
                            let m_s = ib.cols as usize;
                            // KI = (f ⊗ id) ι, sparse in f: row(rp) += f[rp, colp] · ι[g, :]
                            let mut ki = vec![F::zero(); kr * m_s];
                            let mut any = false;
                            for (g, &(cp, colp, off)) in ki_.colinfo.iter().enumerate() {
                                let Some(fi) = f.lay.index(cp) else { continue };
                                let irow = &ib.a[g * m_s..(g + 1) * m_s];
                                if irow.iter().all(|v| v.is_zero()) {
                                    continue;
                                }
                                let (fr, fc) = (f.lay.rows[fi] as usize, f.lay.cols[fi] as usize);
                                let blk = f.block(fi);
                                for rp in 0..fr {
                                    let v = &blk[rp * fc + colp as usize];
                                    if v.is_zero() {
                                        continue;
                                    }
                                    any = true;
                                    let row = ki_.rowmap[off as usize + rp] as usize;
                                    for (d, x) in ki[row * m_s..(row + 1) * m_s].iter_mut().zip(irow) {
                                        if !x.is_zero() {
                                            d.fma(v, x);
                                        }
                                    }
                                }
                            }
                            if !any {
                                continue;
                            }
                            let nzr: Vec<usize> = (0..kr).filter(|&r| ki[r * m_s..(r + 1) * m_s].iter().any(|v| !v.is_zero())).collect();
                            for (it, tu) in db.summands.iter().enumerate() {
                                let pb_ = ev.block(tu.pi, cc);
                                if pb_.rows == 0 {
                                    continue;
                                }
                                let mt = pb_.rows as usize;
                                let mut m = vec![F::zero(); mt * m_s];
                                // m = π · KI over the nonzero rows of KI only
                                for &r in &nzr {
                                    let kir = &ki[r * m_s..(r + 1) * m_s];
                                    for i in 0..mt {
                                        let pv = &pb_.a[i * kr + r];
                                        if pv.is_zero() {
                                            continue;
                                        }
                                        for (d, x) in m[i * m_s..(i + 1) * m_s].iter_mut().zip(kir) {
                                            if !x.is_zero() {
                                                d.fma(pv, x);
                                            }
                                        }
                                    }
                                }
                                if m.iter().any(|v| !v.is_zero()) {
                                    per[ia][it].insert(cc, (pb_.rows as usize, m_s, m));
                                }
                            }
                        }
                    }
                    for (ia, sg) in da.summands.iter().enumerate() {
                        for (it, tu) in db.summands.iter().enumerate() {
                            if per[ia][it].is_empty() {
                                continue;
                            }
                            let lay = eng.layout(sg.q, tu.q).map_err(String::from)?;
                            let mut lm = LocMap::zeros(lay.clone());
                            for (cc, (r, cl, m)) in std::mem::take(&mut per[ia][it]) {
                                let li = lay.index(cc).ok_or_else(|| format!("coset {} of π K ι outside the layout of the summand types", cc))?;
                                if lay.rows[li] as usize != r || lay.cols[li] as usize != cl {
                                    return Err(format!("π K ι block at coset {}: {}x{} vs layout {}x{}", cc, r, cl, lay.rows[li], lay.cols[li]));
                                }
                                lm.block_mut(li).clone_from_slice(&m);
                            }
                            let (x, y) = (place[&(a, k)][ia], place[&(b, k)][it]);
                            eng.add_map(x, y, lm).map_err(String::from)?;
                        }
                    }
                }
            }
            // (−1)^h id ⊗ d from the previously inserted term
            if si > 0 {
                let kp = order[si - 1];
                let j = kp.min(k);
                for (&o, ob) in prev.objs() {
                    let comps = dd.rickard_comp_at(ob.p, i, pos, j).map_err(String::from)?;
                    let dp = dd.dec_word(ob.p, &terms[kp].rungs).map_err(String::from)?;
                    let dk = &decs[&ob.p];
                    let sgn = if ob.h.rem_euclid(2) == 1 { &mone } else { &one };
                    for &(s_i, t_i, m) in comps.iter() {
                        let lm = lib_map(&mut ev, &mut eng, m, dp.summands[s_i].q, dk.summands[t_i].q, sgn)?;
                        if lm.is_zero() {
                            continue;
                        }
                        let (x, y) = (place[&(o, kp)][s_i], place[&(o, k)][t_i]);
                        eng.add_map(x, y, lm).map_err(String::from)?;
                    }
                }
            }
            eng.set_extern_bytes(ev.memo_bytes() as u64);
            let te = Instant::now();
            eng.end_step().map_err(String::from)?;
            t_engine += te.elapsed().as_secs_f64();
        }
        drop(prev);
        let te = Instant::now();
        eng.end_crossing().map_err(String::from)?;
        t_engine += te.elapsed().as_secs_f64();
        labels.swap(i, i + 1);
        crossings.push((pos, c.min(d)));
        if let Some(dc) = decat.as_mut() {
            dc.step(c, d, pos, i);
            let want = dc.k0(lib.prune_n);
            let got: BTreeMap<(String, i32), i64> = eng.cx().k0().into_iter().map(|((p, t), v)| ((lib.names.get(&p).cloned().unwrap_or_else(|| format!("tid {}", p)), t), v)).collect();
            if want != got {
                return Err(format!("crossing {}: K_0 class {:?} != decategorified scan {:?}", st + 1, got, want));
            }
            let lb: u64 = want.values().map(|c| c.unsigned_abs()).sum();
            let n = eng.cx().n_objs();
            if (n as u64) < lb {
                return Err(format!("crossing {}: {} objects < lower bound {}", st + 1, n, lb));
            }
            lbs.push((n, lb));
        }
        if opts.close_stages.contains(&(st + 1)) {
            if labels != lib.lam {
                return Err(format!("stage {} is not closable (labels {:?})", st + 1, labels));
            }
            outs.push(extract_stage(lib, &mut ev, eng.cx(), st + 1, &crossings)?);
        }
    }
    let eval_bytes = ev.memo_bytes();
    let (ps, _) = eng.finish(point_id).map_err(String::from)?;
    t_build = t0s.elapsed().as_secs_f64() - t_engine;
    Ok(ScanOut { ps, stages: outs, lb: lbs, eval_bytes, t_build, t_engine })
}

/// Closure skeleton and entry values of the current complex at the identity coset.
fn extract_stage<F: Field>(lib: &WLib<F>, ev: &mut Eval<F>, cx: &SCx<F>, stage: usize, crossings: &[(bool, u8)]) -> Result<StageOut<F>, String> {
    let cid = lib.id_coset.ok_or("no closure data in the library")?;
    let ids: Vec<usize> = cx.objs().keys().copied().collect();
    let index: HashMap<usize, usize> = ids.iter().enumerate().map(|(i, &o)| (o, i)).collect();
    let objs: Vec<(u32, i32, i32)> = ids.iter().map(|o| { let x = cx.objs()[o]; (x.p, x.h, x.t) }).collect();
    let mut types: BTreeMap<u32, TypeBasisC<F>> = BTreeMap::new();
    // generator vectors at the identity coset (reduced coordinates = flows)
    let mut gv: HashMap<Tid, (Vec<Vec<F>>, Vec<Vec<F>>)> = HashMap::new();
    for &(p, _, _) in &objs {
        if types.contains_key(&p) {
            continue;
        }
        let c = lib.closure.get(&p).ok_or_else(|| format!("final type {} has no closure generators: {}", lib.names[&p], lib.no_closure.get(&p).cloned().unwrap_or_default()))?;
        let g: Vec<Vec<F>> = c.ins.iter().map(|&m| ev.block(m, cid).a.clone()).collect();
        let h: Vec<Vec<F>> = c.outs.iter().map(|&m| ev.block(m, cid).a.clone()).collect();
        types.insert(p, TypeBasisC { ins: c.ins_deg.clone(), outs: c.outs_deg.clone(), gram: vec![vec![EPoly::zero(); c.ins.len()]; c.outs.len()] });
        gv.insert(p, (g, h));
    }
    let mut edges = vec![];
    let mut dblk: Vec<(usize, usize, Vec<F>)> = vec![];
    for (&a, outs) in cx.out() {
        for (&b, m) in outs {
            edges.push((index[&a], index[&b]));
            match m.lay.index(cid) {
                Some(i) => dblk.push((m.lay.rows[i] as usize, m.lay.cols[i] as usize, m.block(i).to_vec())),
                None => dblk.push((0, 0, vec![])),
            }
        }
    }
    let skel = ClosureDataC { lam: lib.lam.clone(), objs: objs.clone(), types, edges: vec![], crossings: crossings.to_vec() };
    let (kinds, weights) = closure_entries(&skel, &edges);
    let dot = |a: &[F], b: &[F]| a.iter().zip(b).fold(F::zero(), |s, (x, y)| s.add(&x.mul(y)));
    let mut vals = Vec::with_capacity(kinds.len());
    // D·g per (edge, a) memo
    let mut dg: HashMap<(usize, usize), Vec<F>> = HashMap::new();
    for k in &kinds {
        let v = match *k {
            EntryKind::Gram { ty, b, a } => {
                let (g, h) = &gv[&ty];
                dot(&h[b], &g[a])
            }
            EntryKind::Pair { edge, b, a } => {
                let (s, t) = edges[edge];
                let (r, c, blk) = &dblk[edge];
                if *r == 0 {
                    F::zero()
                } else {
                    let x = dg.entry((edge, a)).or_insert_with(|| {
                        let g = &gv[&objs[s].0].0[a];
                        let mut y = vec![F::zero(); *r];
                        gemm_acc(&mut y, blk, g, *r, *c, 1, F::one());
                        y
                    });
                    dot(&gv[&objs[t].0].1[b], x)
                }
            }
        };
        vals.push(v);
    }
    Ok(StageOut { stage, skel, edges, kinds, weights, vals })
}

// ------------------------------------------------------------------ box points

/// Arena coordinates of the box point z_μ: the roots a_{μ_j + λ_i + 1 − j} of every block.
pub fn box_x<F: PrimeField>(nodes: &BoxNodes<F>, mu: &[Vec<u8>]) -> Vec<F> {
    let mut x = vec![];
    for (i, m) in mu.iter().enumerate() {
        let a = m.len();
        for j in 0..a {
            x.push(nodes.a[i][m[j] as usize + a - 1 - j].clone());
        }
    }
    x
}

/// A generic record point (smap's point generator).
pub fn record_x<F: Field>(k: usize, seed: u64) -> Vec<F> {
    crate::smap::points::<F>(k, seed, 1, 0).remove(0).x
}

/// Per-stage closures of one word.
#[derive(Clone, Debug)]
pub struct StageResult {
    pub stage: usize,
    /// (N', marked) → Poincaré polynomial
    pub closed: BTreeMap<(usize, Option<usize>), Poincare>,
    pub entries: usize,
    pub objects: usize,
    pub d2_blocks: BTreeMap<(usize, Option<usize>), usize>,
}

#[derive(Clone, Debug, Default)]
pub struct PointCost {
    pub point_id: u32,
    pub secs: f64,
    pub peak_bytes: u64,
    pub eval_bytes: usize,
    /// cumulative engine seconds after each crossing, and the engine peak bytes per crossing
    pub stage_secs: Vec<f64>,
    pub stage_peak: Vec<u64>,
}

#[derive(Clone, Debug)]
pub struct BoxWordOut {
    pub stages: Vec<StageResult>,
    pub record: Vec<StageStat>,
    pub record_secs: f64,
    pub lb: Vec<(usize, u64)>,
    pub points: usize,
    pub costs: Vec<PointCost>,
    pub scan_secs: f64,
    pub close_secs: f64,
    pub final_types: BTreeSet<String>,
}

#[derive(Clone, Debug)]
pub struct BoxOpts {
    pub n_box: usize,
    pub seed: u64,
    /// closures computed: N' in this list, unreduced (None) and reduced (Some(0))
    pub n_close: Vec<usize>,
    pub modes: Vec<Option<usize>>,
    pub stages: Vec<usize>,
    pub eng: EngOpts,
    /// check_d2 at every box point (otherwise only at the record point, when eng.check_d2)
    pub d2_all_points: bool,
    pub freivalds_k: usize,
}

/// The one-word box loop: record at a generic point, replay at every needed box point of the
/// unreduced grid (the reduced grid's points are the μ^0 = ∅ subset), accumulate the entries of
/// every stage per mode, assemble `BoxClosureData`, record Freivalds `D2Cert`s and close.
pub fn box_word<F: PrimeField>(lib: &WLib<F>, word: &[i32], o: &BoxOpts) -> Result<BoxWordOut, String> {
    let k: usize = lib.lam.iter().map(|&a| a as usize).sum();
    let t0 = Instant::now();
    let mut ropts = ScanOpts { eng: o.eng.clone(), close_stages: o.stages.clone(), k0_check: true };
    ropts.eng.id_coset = lib.id_coset;
    let rec = scan_point(lib, word, record_x::<F>(k, o.seed ^ 0x7265_636f_7264), u32::MAX, Mode::Record, &ropts)?;
    let record_secs = t0.elapsed().as_secs_f64();
    let frozen = rec.ps.frozen.clone();
    let nodes = BoxNodes::<F>::new(&lib.lam, o.n_box, o.seed);
    let grids: Vec<(Option<usize>, Arc<BoxGrid<F>>)> = o.modes.iter().map(|&m| (m, Arc::new(BoxGrid::new(nodes.clone(), m)))).collect();
    let whash = word_hash(&[word.to_vec()]);
    // accumulators per (stage, mode)
    let mut accs: Vec<Vec<BoxAcc<F>>> = rec
        .stages
        .iter()
        .map(|so| grids.iter().map(|(_, g)| BoxAcc::new(g.clone(), so.weights.clone(), lib.hash, whash, lib.prune_n.map(|n| n as usize))).collect())
        .collect();
    let ug = BoxGrid::new(nodes.clone(), None);
    let mut maxsize = 0;
    for a in accs.iter().flatten() {
        maxsize = maxsize.max(a.needed_size().unwrap_or(0));
    }
    let pts = ug.points_up_to(maxsize);
    let mut popts = ScanOpts { eng: o.eng.clone(), close_stages: o.stages.clone(), k0_check: false };
    popts.eng.id_coset = lib.id_coset;
    popts.eng.check_d2 = o.eng.check_d2 && o.d2_all_points;
    let ts = Instant::now();
    let res: Vec<Result<(u32, Vec<Vec<F>>, PointCost), String>> = par_map(&pts, |&id| {
        let tp = Instant::now();
        let mu = ug.mu(id);
        let out = scan_point(lib, word, box_x(&nodes, &mu), id, Mode::Replay(frozen.clone()), &popts).map_err(|e| format!("box point {:?}: {}", mu, e))?;
        out.ps.same_structure(&rec.ps).map_err(|e| format!("box point {:?}: {}", mu, e))?;
        for (a, b) in out.stages.iter().zip(rec.stages.iter()) {
            if a.skel.objs != b.skel.objs || a.edges != b.edges {
                return Err(format!("box point {:?}: closure skeleton of stage {} differs from the record point", mu, a.stage));
            }
        }
        let mut acc = 0.0;
        let stage_secs = out.ps.stages.iter().map(|s| { acc += s.secs; acc }).collect();
        let stage_peak = out.ps.stages.iter().map(|s| s.peak_bytes).collect();
        let cost = PointCost { point_id: id, secs: tp.elapsed().as_secs_f64(), peak_bytes: out.ps.mem.peak_bytes, eval_bytes: out.eval_bytes, stage_secs, stage_peak };
        Ok((id, out.stages.into_iter().map(|s| s.vals).collect(), cost))
    });
    let scan_secs = ts.elapsed().as_secs_f64();
    let mut costs = vec![];
    for r in res {
        let (id, vals, cost) = r?;
        costs.push(cost);
        let mu = ug.mu(id);
        for (si, v) in vals.iter().enumerate() {
            for (gi, (m, g)) in grids.iter().enumerate() {
                let gid = match m {
                    None => Some(id),
                    Some(mk) => if mu[*mk].iter().all(|&x| x == 0) { g.id_of_mu(&mu) } else { None },
                };
                let Some(gid) = gid else { continue };
                let acc = &mut accs[si][gi];
                if g.size(gid) <= acc.needed_size().unwrap_or(0) {
                    acc.add_point(gid, v)?;
                }
            }
        }
    }
    let tc = Instant::now();
    let mut stages = vec![];
    for (si, so) in rec.stages.iter().enumerate() {
        let mut sr = StageResult { stage: so.stage, closed: BTreeMap::new(), entries: so.kinds.len(), objects: so.skel.objs.len(), d2_blocks: BTreeMap::new() };
        for (gi, (m, g)) in grids.iter().enumerate() {
            let res = accs[si][gi].finish().map_err(|e| format!("stage {} mode {:?}: {}", so.stage, m, e))?;
            let stamp = BoxStamp { prime: F::characteristic(), node_seed_hash: g.nodes.hash(), library_hash: lib.hash, word_hash: whash, code_stamp: BOX_CODE_STAMP.into() };
            let mut bd = BoxClosureData::from_box(&so.skel, &so.edges, &so.kinds, &res, lib.prune_n.map(|n| n as usize), stamp)?;
            for &nn in &o.n_close {
                for marked in [None, Some(0)] {
                    if m.is_some() && marked != *m {
                        continue;
                    }
                    if m.is_none() && marked.is_some() && o.modes.contains(&marked) {
                        continue; // served by the reduced extraction
                    }
                    let cert = closed_d2_freivalds(&mut bd, nn, marked, o.freivalds_k, o.seed ^ (nn as u64 * 7919))?;
                    sr.d2_blocks.insert((nn, marked), cert.blocks);
                    let p = close_box(&bd, nn, marked, None, None).map_err(|e| format!("stage {} N' {} marked {:?}: {}", so.stage, nn, marked, e))?;
                    sr.closed.insert((nn, marked), p);
                }
            }
        }
        stages.push(sr);
    }
    let final_types = rec.stages.last().map(|s| s.skel.objs.iter().map(|o| lib.names[&o.0].clone()).collect()).unwrap_or_default();
    Ok(BoxWordOut { stages, record: rec.ps.stages.clone(), record_secs, lb: rec.lb.clone(), points: pts.len(), costs, scan_secs, close_secs: tc.elapsed().as_secs_f64(), final_types })
}

/// D14 scaling check: scans at z_μ and λ·z_μ (scaled nodes) for the given points; every entry of
/// every stage must satisfy v(λz) = λ^w v(z). Returns (checked, zero at both points).
pub fn scaling_check<F: PrimeField>(lib: &WLib<F>, word: &[i32], o: &BoxOpts, lambda: &F, point_ids: &[u32]) -> Result<(usize, usize), String> {
    let k: usize = lib.lam.iter().map(|&a| a as usize).sum();
    let mut ropts = ScanOpts { eng: o.eng.clone(), close_stages: o.stages.clone(), k0_check: false };
    ropts.eng.id_coset = lib.id_coset;
    let rec = scan_point(lib, word, record_x::<F>(k, o.seed ^ 0x7265_636f_7264), u32::MAX, Mode::Record, &ropts)?;
    let nodes = BoxNodes::<F>::new(&lib.lam, o.n_box, o.seed);
    let sn = nodes.scaled(lambda);
    let grid = BoxGrid::new(nodes.clone(), None);
    let (mut checked, mut zero) = (0, 0);
    for &id in point_ids {
        let mu = grid.mu(id);
        let a = scan_point(lib, word, box_x(&nodes, &mu), id, Mode::Replay(rec.ps.frozen.clone()), &ropts)?;
        let b = scan_point(lib, word, box_x(&sn, &mu), id, Mode::Replay(rec.ps.frozen.clone()), &ropts)?;
        for (sa, sb) in a.stages.iter().zip(b.stages.iter()) {
            let r = check_scaling_all(&sa.vals, &sb.vals, &sa.weights, lambda).map_err(|e| format!("stage {} point {:?}: {}", sa.stage, mu, e))?;
            checked += r.checked;
            zero += r.zero_both;
        }
    }
    Ok((checked, zero))
}

// ------------------------------------------------------------------ replay-only provider

/// A provider that refuses every call: the library must come entirely from layer A.
pub struct NoProvider;
fn refuse(what: &str) -> Vec<crate::homprov::Missing> {
    vec![crate::homprov::Missing { what: format!("replay: {} not in layer A", what), degree: 0, stratum: None, got: 0, want: 1, kind: crate::homprov::MissKind::Refused }]
}
impl<F: Field> crate::homprov::HomProvider<F> for NoProvider {
    fn name(&self) -> &'static str {
        "NoProvider"
    }
    fn full(&mut self, _: &mut crate::smap::Arena<F>, _: &mut Schur, _: crate::smap::ExprId, _: crate::smap::ExprId, _: i32, _: &crate::homprov::Budget) -> Result<(Vec<MapId>, crate::homprov::SearchStats), Vec<crate::homprov::Missing>> {
        Err(refuse("full"))
    }
    fn pairs(&mut self, _: &mut crate::smap::Arena<F>, _: &mut Schur, _: crate::smap::ExprId, _: Tid, _: crate::smap::ExprId, _: Option<MapId>, _: i32, _: usize, _: &crate::homprov::Budget) -> Result<(crate::homprov::Pairs<F>, crate::homprov::SearchStats), Vec<crate::homprov::Missing>> {
        Err(refuse("pairs"))
    }
    fn closure(&mut self, _: &mut crate::smap::Arena<F>, _: &mut Schur, _: crate::smap::ExprId, _: Tid, _: Option<MapId>, _: bool, _: &crate::homprov::Budget) -> Result<(Vec<MapId>, crate::homprov::SearchStats), Vec<crate::homprov::Missing>> {
        Err(refuse("closure"))
    }
    fn end0(&mut self, _: &mut crate::smap::Arena<F>, _: &mut Schur, _: crate::smap::ExprId, _: &crate::homprov::Budget) -> Result<(Vec<MapId>, crate::homprov::SearchStats), Vec<crate::homprov::Missing>> {
        Err(refuse("end0"))
    }
}
