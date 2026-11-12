//! Exact window maps (sbim::hom_space, cube::local_cx), localized per window root tuple,
//! applied locally inside larger ladders; polynomiality certificates; disk cache.
//!
//! A window is a ladder on w adjacent uprights (labels at the window bottom, rungs relative to
//! the window's first upright). Window roots are numbered as `local::root_offsets(labels)`; at a
//! global flow they are the roots of the window uprights at the window's bottom level, ascending
//! inside each upright (`smap::Eval`). Window flows are sorted as in `FlowIndex`, which is the
//! row/column order of every `WinLoc`.
use crate::field::Field;
use crate::local::{esym, eval_pmat, fixed_points, inverse, loc_matrix, matmul, FixedPoint};
use crate::mpoly::{mono_exp, MonoCache, Poly};
use crate::pmat::{PAcc, PMat, PVec};
use crate::sbim::{hom_space, SBim};
use crate::sexpr::{Expr, Flow, FlowIndex};
use crate::web::{Dir, Rung};
use std::collections::HashMap;
use std::sync::Arc;
use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Window {
    pub labels: Vec<u8>,
    pub rungs: Vec<Rung>,
} // rungs relative to the window's first upright

impl Window {
    pub fn new(labels: &[u8], rungs: &[Rung]) -> Window {
        Window { labels: labels.to_vec(), rungs: rungs.to_vec() }
    }
    /// the window as an expression (panics if invalid)
    pub fn expr(&self) -> Expr {
        Expr::new(&self.labels, &self.rungs).unwrap_or_else(|e| panic!("window {:?}: {}", self, e))
    }
    /// canonical text key, e.g. "1,1|E0.1,F0.1"
    pub fn key(&self) -> String {
        let l: Vec<String> = self.labels.iter().map(|x| x.to_string()).collect();
        let r: Vec<String> = self.rungs.iter().map(|r| format!("{}{}.{}", if r.dir == Dir::E { 'E' } else { 'F' }, r.i, r.r)).collect();
        format!("{}|{}", l.join(","), r.join(","))
    }
}

pub type TmplId = u32;

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum TmplSource {
    HomSpace { index: usize },
    Rickard { c: u8, d: u8, positive: bool, j: usize },
    /// added through `TmplLib::add_template` (tests, idempotents): uncertified until
    /// `TmplLib::certify_added` has checked it (degree, exact bimodule check, polynomiality
    /// certificate with block diagonality); `Arena` and `assert_genuine` treat it as uncertified
    /// before that (M4a review finding b)
    Added { certified: bool },
}

#[derive(Clone, Debug)]
pub struct Template<F: Field> {
    pub src: Window,
    pub tgt: Window,
    pub deg: i32,     // normalized
    pub phi: PMat<F>, // exact over R^{src.labels}, unnormalized degree
    pub source: TmplSource,
    pub src_flows: Vec<Flow>, // window-local moved masks, sorted as in FlowIndex
    pub tgt_flows: Vec<Flow>, // (point-independent; computed once)
}

impl<F: Field> Template<F> {
    /// computes the sorted window flows; checks shapes (equal bottoms and tops, PMat size)
    pub fn new(src: Window, tgt: Window, deg: i32, phi: PMat<F>, source: TmplSource) -> Template<F> {
        let se = src.expr();
        let te = tgt.expr();
        assert_eq!(src.labels, tgt.labels, "template: window bottoms differ");
        assert_eq!(se.top(), te.top(), "template: window tops differ");
        let sf = FlowIndex::build(&se);
        let tf = FlowIndex::build(&te);
        assert_eq!((phi.nrows, phi.ncols), (tf.nflows(), sf.nflows()), "template: PMat shape != (tgt, src) fibre sizes");
        let src_flows = (0..sf.nflows() as u32).map(|f| Flow(sf.flow(f).to_vec())).collect();
        let tgt_flows = (0..tf.nflows() as u32).map(|f| Flow(tf.flow(f).to_vec())).collect();
        Template { src, tgt, deg, phi, source, src_flows, tgt_flows }
    }
}

/// Localized window map at one root tuple, sparse by column (CSC over the sorted window flows):
/// the nonzeros of column j (source flow j) are (rows[k], vals[k]) for k in col_start[j]..col_start[j+1].
/// Window maps are block diagonal over window top configurations, so columns hold few entries.
#[derive(Clone, Debug, PartialEq)]
pub struct WinLoc<F> {
    pub nrows: u32,
    pub col_start: Vec<u32>,
    pub rows: Vec<u32>,
    pub vals: Vec<F>,
}

impl<F: Field> WinLoc<F> {
    /// rows = tgt flows, cols = src flows; tests and oracles only
    pub fn dense(&self) -> Vec<Vec<F>> {
        let nc = self.ncols();
        let mut d = vec![vec![F::zero(); nc]; self.nrows as usize];
        for j in 0..nc {
            for k in self.col_start[j] as usize..self.col_start[j + 1] as usize {
                d[self.rows[k] as usize][j] = self.vals[k].clone();
            }
        }
        d
    }
    pub fn ncols(&self) -> usize {
        self.col_start.len() - 1
    }
    /// CSC of a dense matrix (rows x ncols)
    pub fn from_dense(m: &[Vec<F>], ncols: usize) -> WinLoc<F> {
        let mut col_start = vec![0u32];
        let mut rows = vec![];
        let mut vals = vec![];
        for j in 0..ncols {
            for (i, r) in m.iter().enumerate() {
                if !r[j].is_zero() {
                    rows.push(i as u32);
                    vals.push(r[j].clone());
                }
            }
            col_start.push(rows.len() as u32);
        }
        WinLoc { nrows: m.len() as u32, col_start, rows, vals }
    }
    pub fn bytes(&self) -> usize {
        self.col_start.capacity() * 4 + self.rows.capacity() * 4 + self.vals.capacity() * std::mem::size_of::<F>()
    }
}

/// localization data of one window (its SBim, fixed points, fixed point -> sorted flow id);
/// shared by every template with this window as source or target (M5c: a hom space of
/// dimension k used to store k clones of both whole-ladder SBims, several GB for the (3,3,3)
/// whole-ladder Exact pair calls of the Memory pairs)
#[derive(Debug)]
struct Side<F: Field> {
    b: SBim<F>,
    fp: Vec<FixedPoint>,
    pos: Vec<u32>,
}

/// per-template data for localization (source and target window sides, shared)
#[derive(Clone, Debug)]
struct Aux<F: Field> {
    src: Arc<Side<F>>,
    tgt: Arc<Side<F>>,
}

fn positions(fps: &[FixedPoint], fi: &FlowIndex) -> Vec<u32> {
    fps.iter().map(|p| fi.locate(&p.moved).expect("fixed point is not a flow of the window")).collect()
}

#[derive(Clone, Debug)]
enum Call {
    Hom(Window, Window, i32, Vec<TmplId>),
    Rick(u8, u8, bool, Vec<TmplId>),
}

/// counters of the template library
#[derive(Clone, Debug, Default)]
pub struct LibStats {
    pub computed: usize,
    pub loaded: usize,
    pub rejected: usize,
    pub compute_secs: f64,
}

/// Send + Sync (A16): no Rc, RefCell or Cell inside.
pub struct TmplLib<F: Field> {
    cache_dir: Option<PathBuf>,
    tmpls: Vec<Template<F>>,
    aux: Vec<Aux<F>>,
    /// one localization side per window (shared by the templates' `Aux`)
    sides: HashMap<Window, Arc<Side<F>>>,
    hom: HashMap<(Window, Window, i32), Vec<TmplId>>,
    rick: HashMap<(u8, u8, bool), (Vec<Window>, Vec<TmplId>)>,
    calls: Vec<Call>,
    stats: LibStats,
    /// identity of this library (process-unique nonce): template ids are only meaningful within
    /// the library that issued them; caches of TmplIds outside the library (homprov::Search) are
    /// bound to it (M4b review F3)
    instance: u64,
}

/// process-unique library identities
static LIB_INSTANCE: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(1);

/// does the first line of a cache file equal `header`? (a stale or foreign file is not "cached")
fn header_ok(path: &Path, header: &str) -> bool {
    let Ok(f) = std::fs::File::open(path) else { return false };
    let mut l = String::new();
    if std::io::BufReader::new(f).read_line(&mut l).is_err() {
        return false;
    }
    l.trim_end_matches(['\n', '\r']) == header
}

/// admission bound of `hom_basis`: fsize(src) * fsize(tgt)
pub const HOM_WINDOW_CAP: u128 = 1600;

fn fnv64(s: &str) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for b in s.bytes() {
        h ^= b as u64;
        h = h.wrapping_mul(0x0100_0000_01b3);
    }
    h
}

/// Version stamp of the algorithms whose output the cache stores (`sbim::hom_space` basis order
/// and normalization, `cube::local_cx`). Bump it whenever either changes: files with another
/// stamp are rejected and recomputed (M4a review finding f).
pub const TMPL_ALG: &str = "hs1.lcx1";

/// cache header line for a call
pub fn cache_header<F: Field>(kind: &str, key: &str) -> String {
    format!("krhom-tmpl v2 alg={} field={} kind={} key={}", TMPL_ALG, F::name(), kind, key)
}

/// per-process counter of temporary cache file names
static TMP_COUNTER: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

fn hom_key(src: &Window, tgt: &Window, deg: i32) -> String {
    format!("src={};tgt={};deg={}", src.key(), tgt.key(), deg)
}
fn rick_key(c: u8, d: u8, positive: bool) -> String {
    format!("c={};d={};positive={}", c, d, positive)
}

fn write_cache_file<F: Field>(path: &Path, header: &str, maps: &[&PMat<F>]) -> std::io::Result<()> {
    // per-process, per-call temporary name: concurrent writers never share a tmp file
    let n = TMP_COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let tmp = path.with_extension(format!("{}.{}.tmp", std::process::id(), n));
    {
        let mut w = std::io::BufWriter::new(std::fs::File::create(&tmp)?);
        writeln!(w, "{}", header)?;
        writeln!(w, "maps {}", maps.len())?;
        for m in maps {
            let nt: usize = m.cols.iter().map(|c| c.iter().map(|(_, p)| p.t.len()).sum::<usize>()).sum();
            writeln!(w, "map {} {} {}", m.nrows, m.ncols, nt)?;
            for (j, col) in m.cols.iter().enumerate() {
                for (i, p) in col {
                    for (mono, c) in &p.t {
                        writeln!(w, "{} {} {:x} {}", j, i, mono, c)?;
                    }
                }
            }
        }
        writeln!(w, "end")?;
        w.flush()?;
    }
    let r = std::fs::rename(&tmp, path);
    if r.is_err() {
        let _ = std::fs::remove_file(&tmp);
    }
    r
}

/// Read a cache file; Err if the header differs from `header` or the body is malformed.
pub fn read_cache_file<F: Field>(path: &Path, header: &str) -> Result<Vec<PMat<F>>, String> {
    let f = std::fs::File::open(path).map_err(|e| format!("{}: {}", path.display(), e))?;
    let mut lines = std::io::BufReader::new(f).lines();
    let mut next = || -> Result<String, String> { lines.next().ok_or_else(|| "unexpected end of file".to_string())?.map_err(|e| e.to_string()) };
    let h = next()?;
    if h != header {
        return Err(format!("{}: header mismatch: {:?} != {:?}", path.display(), h, header));
    }
    let l = next()?;
    let nm: usize = l.strip_prefix("maps ").and_then(|x| x.trim().parse().ok()).ok_or_else(|| format!("bad line {:?}", l))?;
    let mut out = vec![];
    for _ in 0..nm {
        let l = next()?;
        let p: Vec<usize> = l.strip_prefix("map ").ok_or_else(|| format!("bad line {:?}", l))?.split_whitespace().map(|x| x.parse::<usize>().map_err(|e| e.to_string())).collect::<Result<_, _>>()?;
        if p.len() != 3 {
            return Err(format!("bad line {:?}", l));
        }
        let (nrows, ncols, nt) = (p[0], p[1], p[2]);
        let mut cols: Vec<Vec<(u32, Vec<(u128, F)>)>> = vec![vec![]; ncols];
        for _ in 0..nt {
            let l = next()?;
            let t: Vec<&str> = l.split_whitespace().collect();
            if t.len() != 4 {
                return Err(format!("bad term line {:?}", l));
            }
            let j: usize = t[0].parse().map_err(|_| format!("bad term line {:?}", l))?;
            let i: u32 = t[1].parse().map_err(|_| format!("bad term line {:?}", l))?;
            let mono = u128::from_str_radix(t[2], 16).map_err(|_| format!("bad term line {:?}", l))?;
            let c: i64 = t[3].parse().map_err(|_| format!("bad coefficient in {:?}", l))?;
            if j >= ncols || i as usize >= nrows {
                return Err(format!("term out of range {:?}", l));
            }
            let col = &mut cols[j];
            if col.last().map(|e| e.0) != Some(i) {
                col.push((i, vec![]));
            }
            col.last_mut().unwrap().1.push((mono, F::from_i64(c)));
        }
        let cols: Vec<PVec<F>> = cols.into_iter().map(|c| c.into_iter().map(|(i, t)| (i, Poly { t })).collect()).collect();
        out.push(PMat { nrows, ncols, cols });
    }
    if next()? != "end" {
        return Err("missing end line".into());
    }
    Ok(out)
}

impl<F: Field> TmplLib<F> {
    /// cache file names carry F::name() (through the hashed key and the header)
    pub fn new(cache_dir: Option<PathBuf>) -> Self {
        if let Some(d) = &cache_dir {
            let _ = std::fs::create_dir_all(d);
        }
        TmplLib { cache_dir, tmpls: vec![], aux: vec![], sides: HashMap::new(), hom: HashMap::new(), rick: HashMap::new(), calls: vec![], stats: LibStats::default(), instance: LIB_INSTANCE.fetch_add(1, std::sync::atomic::Ordering::Relaxed) }
    }
    /// process-unique identity of this library (TmplIds are valid only here)
    pub fn instance(&self) -> u64 {
        self.instance
    }
    pub fn get(&self, t: TmplId) -> &Template<F> {
        &self.tmpls[t as usize]
    }
    pub fn len(&self) -> usize {
        self.tmpls.len()
    }
    pub fn is_empty(&self) -> bool {
        self.tmpls.is_empty()
    }
    pub fn stats(&self) -> LibStats {
        self.stats.clone()
    }
    pub fn cache_dir(&self) -> Option<&Path> {
        self.cache_dir.as_deref()
    }
    /// path of the cache file of a call (kind "hom" or "rickard")
    pub fn cache_path(&self, kind: &str, key: &str) -> Option<PathBuf> {
        let full = format!("{}|{}|{}", F::name(), kind, key);
        self.cache_dir.as_ref().map(|d| d.join(format!("tmpl-{:016x}.txt", fnv64(&full))))
    }
    /// cache path of a hom_basis call
    pub fn hom_cache_path(&self, src: &Window, tgt: &Window, deg: i32) -> Option<PathBuf> {
        self.cache_path("hom", &hom_key(src, tgt, deg))
    }
    /// cache path of a rickard call
    pub fn rickard_cache_path(&self, c: u8, d: u8, positive: bool) -> Option<PathBuf> {
        self.cache_path("rickard", &rick_key(c, d, positive))
    }
    fn load(&mut self, kind: &str, key: &str) -> Option<Vec<PMat<F>>> {
        let path = self.cache_path(kind, key)?;
        if !path.exists() {
            return None;
        }
        match read_cache_file::<F>(&path, &cache_header::<F>(kind, key)) {
            Ok(v) => {
                self.stats.loaded += 1;
                Some(v)
            }
            Err(_) => {
                self.stats.rejected += 1;
                None
            }
        }
    }
    fn store(&self, kind: &str, key: &str, maps: &[&PMat<F>]) -> std::io::Result<()> {
        if let Some(path) = self.cache_path(kind, key) {
            write_cache_file(&path, &cache_header::<F>(kind, key), maps)?;
        }
        Ok(())
    }
    /// the shared localization side of a window (built once; `b` = its SBim if the caller has it)
    fn side(&mut self, w: &Window, b: Option<&SBim<F>>) -> Arc<Side<F>> {
        if let Some(s) = self.sides.get(w) {
            return s.clone();
        }
        let b = b.cloned().unwrap_or_else(|| w.expr().sbim());
        let fp = fixed_points(&b);
        let pos = positions(&fp, &FlowIndex::build(&w.expr()));
        let s = Arc::new(Side { b, fp, pos });
        self.sides.insert(w.clone(), s.clone());
        s
    }
    fn push(&mut self, t: Template<F>, sb: Option<&SBim<F>>, tb: Option<&SBim<F>>) -> TmplId {
        let src = self.side(&t.src, sb);
        let tgt = self.side(&t.tgt, tb);
        assert_eq!(src.b.fsize(), t.phi.ncols);
        assert_eq!(tgt.b.fsize(), t.phi.nrows);
        self.aux.push(Aux { src, tgt });
        self.tmpls.push(t);
        (self.tmpls.len() - 1) as TmplId
    }
    /// add a template built elsewhere (tests: dots, idempotents); not written to the cache.
    /// Whatever source the caller gives, the template is recorded as `TmplSource::Added
    /// { certified: false }`: it is not genuine until `certify_added` passes.
    pub fn add_template(&mut self, mut t: Template<F>) -> TmplId {
        t.source = TmplSource::Added { certified: false };
        self.push(t, None, None)
    }
    /// is the template a certified leaf (hom_space, Rickard, or an added template that passed
    /// `certify_added`)?
    pub fn is_certified(&self, t: TmplId) -> bool {
        match &self.tmpls[t as usize].source {
            TmplSource::HomSpace { .. } | TmplSource::Rickard { .. } => true,
            TmplSource::Added { certified } => *certified,
        }
    }
    /// Certify an added template: (1) every entry of φ is homogeneous of the unnormalized degree
    /// of `deg`; (2) φ commutes exactly with every top operator of the window (bimodule map);
    /// (3) the polynomiality certificate of its localization passes (including block
    /// diagonality). On success the template becomes a certified leaf. Templates from hom_space
    /// or Rickard need no certification (Ok with an empty report).
    pub fn certify_added(&mut self, t: TmplId, seed: u64) -> Result<CertReport, String> {
        if !matches!(self.tmpls[t as usize].source, TmplSource::Added { .. }) {
            return Ok(CertReport::default());
        }
        let tm = &self.tmpls[t as usize];
        let a = &self.aux[t as usize];
        let u = tm.deg - a.tgt.b.shift + a.src.b.shift;
        for (b, col) in tm.phi.cols.iter().enumerate() {
            for (bp, p) in col {
                for (mono, _) in &p.t {
                    let want = u + a.src.b.fdeg[b] - a.tgt.b.fdeg[*bp as usize];
                    if a.src.b.ring.deg(*mono) != want {
                        return Err(format!("added template {}: entry ({}, {}) has a term of degree {} != {} (normalized degree {})", t, bp, b, a.src.b.ring.deg(*mono), want, tm.deg));
                    }
                }
            }
        }
        for (up, (os, ot)) in a.src.b.top.iter().zip(&a.tgt.b.top).enumerate() {
            for (j, (ms, mt)) in os.iter().zip(ot).enumerate() {
                if !tm.phi.matmul(ms).sub(&mt.matmul(&tm.phi)).is_zero() {
                    return Err(format!("added template {}: does not commute with the top operator e_{} of upright {} (not a bimodule map)", t, j + 1, up));
                }
            }
        }
        let (src, tgt, deg) = (tm.src.clone(), tm.tgt.clone(), tm.deg);
        let rep = {
            let me: &TmplLib<F> = self;
            poly_certificate::<F>(&src, &tgt, deg, &mut |y: &[F]| me.localize(t, y), seed)?
        };
        self.tmpls[t as usize].source = TmplSource::Added { certified: true };
        Ok(rep)
    }
    /// is the basis of Hom^deg(src, tgt) already available (in memory or in the disk cache)?
    pub fn has_hom(&self, src: &Window, tgt: &Window, deg: i32) -> bool {
        // a file counts only if its header (format version, TMPL_ALG, field, full key) matches:
        // a stale file would be rejected by `load` and recomputed (M4b review, tests finding 4)
        self.hom.contains_key(&(src.clone(), tgt.clone(), deg))
            || self.hom_cache_path(src, tgt, deg).map(|p| header_ok(&p, &cache_header::<F>("hom", &hom_key(src, tgt, deg)))).unwrap_or(false)
    }
    /// exact basis of Hom^deg(src, tgt) (normalized degree), computed once, cached on disk
    pub fn hom_basis(&mut self, src: &Window, tgt: &Window, deg: i32) -> Vec<TmplId> {
        self.try_hom_basis(src, tgt, deg).unwrap_or_else(|e| panic!("hom_basis: {}", e))
    }
    /// as `hom_basis`; Err for invalid or non-admitted windows (fsize(src)·fsize(tgt) > 1 600)
    pub fn try_hom_basis(&mut self, src: &Window, tgt: &Window, deg: i32) -> Result<Vec<TmplId>, String> {
        self.try_hom_basis_cap(src, tgt, deg, HOM_WINDOW_CAP)
    }
    /// as `try_hom_basis` with another fibre-product cap (the M4b Exact backend admits calls
    /// through its own rules: label caps or a measured cost table)
    pub fn try_hom_basis_cap(&mut self, src: &Window, tgt: &Window, deg: i32, cap: u128) -> Result<Vec<TmplId>, String> {
        let k3 = (src.clone(), tgt.clone(), deg);
        if let Some(v) = self.hom.get(&k3) {
            return Ok(v.clone());
        }
        let se = Expr::new(&src.labels, &src.rungs)?;
        let te = Expr::new(&tgt.labels, &tgt.rungs)?;
        if src.labels != tgt.labels || se.top() != te.top() {
            return Err(format!("windows {} and {} have different bottoms or tops", src.key(), tgt.key()));
        }
        if se.k() > crate::mpoly::MAXVARS {
            return Err(format!("window K = {} > {}", se.k(), crate::mpoly::MAXVARS));
        }
        let prod = se.nflows() * te.nflows();
        if prod > cap {
            return Err(format!("window {} -> {}: fsize product {} > {} (not admitted)", src.key(), tgt.key(), prod, cap));
        }
        let key = hom_key(src, tgt, deg);
        let sb: SBim<F> = se.sbim();
        let tb: SBim<F> = te.sbim();
        let maps = match self.load("hom", &key) {
            Some(v) => v,
            None => {
                let t0 = Instant::now();
                let u = deg - tb.shift + sb.shift;
                let mut mc = MonoCache::new(sb.ring.clone());
                let v = hom_space(&sb, &tb, u, &mut mc);
                self.stats.computed += 1;
                self.stats.compute_secs += t0.elapsed().as_secs_f64();
                let refs: Vec<&PMat<F>> = v.iter().collect();
                let _ = self.store("hom", &key, &refs);
                v
            }
        };
        let mut ids = vec![];
        for (index, phi) in maps.into_iter().enumerate() {
            let t = Template::new(src.clone(), tgt.clone(), deg, phi, TmplSource::HomSpace { index });
            ids.push(self.push(t, Some(&sb), Some(&tb)));
        }
        self.hom.insert(k3, ids.clone());
        self.calls.push(Call::Hom(src.clone(), tgt.clone(), deg, ids.clone()));
        Ok(ids)
    }
    /// untruncated Rickard complex of a crossing on labels (c, d): term windows and d_j templates
    /// (positive: d_j from term j to j+1; negative: from term j+1 to term j; normalized degree 1)
    pub fn rickard(&mut self, c: u8, d: u8, positive: bool) -> (Vec<Window>, Vec<TmplId>) {
        if let Some(v) = self.rick.get(&(c, d, positive)) {
            return v.clone();
        }
        let terms = crate::schur::crossing_terms(c, d, positive, 0);
        let wins: Vec<Window> = terms.iter().map(|t| Window::new(&[c, d], &t.rungs)).collect();
        let key = rick_key(c, d, positive);
        let maps = match self.load("rickard", &key) {
            Some(v) => v,
            None => {
                let t0 = Instant::now();
                let loc = crate::cube::local_cx::<F>(crate::crossing::Gen { i: 0, positive }, c as usize, d as usize, (c + d) as usize);
                assert_eq!(loc.cx.terms.len(), terms.len(), "Rickard complex of ({},{}) pruned", c, d);
                for (k, (a, b)) in loc.cx.terms.iter().zip(&terms).enumerate() {
                    assert_eq!(a.rungs, b.rungs, "cube::local_cx and schur::crossing_terms disagree on term {} of ({},{})", k, c, d);
                    assert_eq!(a.hdeg, b.hdeg);
                }
                for (j, &u) in loc.udeg.iter().enumerate() {
                    let (s, t) = if positive { (j, j + 1) } else { (j + 1, j) };
                    assert_eq!(u + loc.shifts[t] - loc.shifts[s], 1, "Rickard map {} of ({},{}) not of normalized degree 1", j, c, d);
                }
                self.stats.computed += 1;
                self.stats.compute_secs += t0.elapsed().as_secs_f64();
                let refs: Vec<&PMat<F>> = loc.maps.iter().collect();
                let _ = self.store("rickard", &key, &refs);
                loc.maps
            }
        };
        assert_eq!(maps.len(), wins.len().saturating_sub(1), "Rickard ({},{}): wrong number of maps", c, d);
        // each d_j is the unique map of its degree (M4a review finding c): dim Hom^1 = 1 and
        // Hom^{<1} = 0 between the term characters (schur::hom_gd), and d_j != 0
        {
            let mut sch = crate::schur::Schur::new(&[c, d]);
            let chars: Vec<crate::schur::Elt> = wins.iter().map(|w| sch.ladder_char(&w.rungs)).collect();
            let kk = (c as i32 + d as i32).max(1);
            let lo = -4 * kk * kk - 4;
            for (j, phi) in maps.iter().enumerate() {
                let (s, t) = if positive { (j, j + 1) } else { (j + 1, j) };
                let gd = sch.hom_gd(&chars[s], &chars[t], lo, 1);
                let n = gd.len();
                assert!(gd[..n - 1].iter().all(|&x| x == 0) && gd[n - 1] == 1, "Rickard ({},{}) {}: d_{} is not the unique map of degree 1 (dim Hom^<=1 = {:?})", c, d, if positive { "+" } else { "-" }, j, &gd[n.saturating_sub(4)..]);
                assert!(!phi.is_zero(), "Rickard ({},{}) {}: d_{} is zero", c, d, if positive { "+" } else { "-" }, j);
            }
        }
        let mut ids = vec![];
        for (j, phi) in maps.into_iter().enumerate() {
            let (s, t) = if positive { (j, j + 1) } else { (j + 1, j) };
            let tm = Template::new(wins[s].clone(), wins[t].clone(), 1, phi, TmplSource::Rickard { c, d, positive, j });
            ids.push(self.push(tm, None, None));
        }
        self.rick.insert((c, d, positive), (wins.clone(), ids.clone()));
        self.calls.push(Call::Rick(c, d, positive, ids.clone()));
        (wins, ids)
    }
    /// localization at window roots y (upright by upright, ascending global root id inside an upright);
    /// not memoized here (the caller, `Eval`, memoizes per point)
    pub fn localize(&self, t: TmplId, y: &[F]) -> WinLoc<F> {
        let tm = &self.tmpls[t as usize];
        let a = &self.aux[t as usize];
        let ls = loc_matrix(&a.src.b, &a.src.fp, y);
        let lt = loc_matrix(&a.tgt.b, &a.tgt.fp, y);
        let linv = inverse(&ls).expect("localization matrix singular at this point (roots not distinct?)");
        let ph = eval_pmat(&a.src.b, &tm.phi, y);
        let f = matmul(&matmul(&lt, &ph), &linv); // rows tgt fixed points, cols src fixed points
        let ns = a.src.fp.len();
        let mut src_of_sorted = vec![0usize; ns];
        for (p, &j) in a.src.pos.iter().enumerate() {
            src_of_sorted[j as usize] = p;
        }
        let mut col_start = Vec::with_capacity(ns + 1);
        col_start.push(0u32);
        let mut rows = vec![];
        let mut vals = vec![];
        let mut col: Vec<(u32, F)> = vec![];
        for j in 0..ns {
            let p = src_of_sorted[j];
            col.clear();
            for (q, r) in f.iter().enumerate() {
                if !r[p].is_zero() {
                    col.push((a.tgt.pos[q], r[p].clone()));
                }
            }
            col.sort_unstable_by_key(|e| e.0);
            for (i, v) in col.drain(..) {
                rows.push(i);
                vals.push(v);
            }
            col_start.push(rows.len() as u32);
        }
        WinLoc { nrows: a.tgt.fp.len() as u32, col_start, rows, vals }
    }
    /// header "krhom-tmpl v1 field=.. "
    pub fn save(&self) -> std::io::Result<()> {
        if self.cache_dir.is_none() {
            return Ok(());
        }
        for call in &self.calls {
            match call {
                Call::Hom(s, t, d, ids) => {
                    let refs: Vec<&PMat<F>> = ids.iter().map(|&i| &self.tmpls[i as usize].phi).collect();
                    self.store("hom", &hom_key(s, t, *d), &refs)?;
                }
                Call::Rick(c, d, p, ids) => {
                    let refs: Vec<&PMat<F>> = ids.iter().map(|&i| &self.tmpls[i as usize].phi).collect();
                    self.store("rickard", &rick_key(*c, *d, *p), &refs)?;
                }
            }
        }
        Ok(())
    }
    /// approximate heap bytes of the templates (PMats, flows, window SBims)
    pub fn bytes(&self) -> usize {
        let mut b = 0usize;
        for t in &self.tmpls {
            b += pmat_bytes(&t.phi);
            b += t.src_flows.iter().chain(&t.tgt_flows).map(|f| f.0.capacity() * 8 + 24).sum::<usize>();
        }
        for sd in self.sides.values() {
            b += sd.b.top.iter().flat_map(|v| v.iter()).map(pmat_bytes).sum::<usize>();
            b += sd.fp.len() * 64 + sd.pos.len() * 4;
        }
        b
    }
}

fn pmat_bytes<F: Field>(m: &PMat<F>) -> usize {
    m.cols.iter().map(|c| 24 + c.iter().map(|(_, p)| 32 + p.t.len() * (16 + std::mem::size_of::<F>())).sum::<usize>()).sum()
}

/// Exact id_P (x) phi (x) id_Q as a PMat (for tests and certificates): phi acts on the uprights
/// first..first+w of the top of `prefix` (substitution of the prefix's top operators, generalizing
/// cube::apply_local_poly to w >= 2 uprights), then kron_id(suffix_fsize).
pub fn tensor_exact<F: Field>(prefix: &SBim<F>, first: usize, t: &Template<F>, suffix_fsize: usize) -> PMat<F> {
    let w = t.src.labels.len();
    let wl: Vec<usize> = t.src.labels.iter().map(|&a| a as usize).collect();
    assert_eq!(&prefix.top_labels[first..first + w], &wl[..], "tensor_exact: window labels != prefix top labels");
    let mut voff = vec![];
    let mut k = 0;
    for &l in &wl {
        voff.push(k);
        k += l;
    }
    let fp = prefix.fsize();
    let (fa, fa2) = (t.phi.ncols, t.phi.nrows);
    // apply a window polynomial to a prefix vector through the prefix's top operators
    let apply = |poly: &Poly<F>, v: &PVec<F>| -> PVec<F> {
        let mut acc = PAcc::new();
        for (m, coef) in &poly.t {
            let mut x = v.clone();
            for j in 0..w {
                for kk in 0..wl[j] {
                    for _ in 0..mono_exp(*m, voff[j] + kk) {
                        x = prefix.top[first + j][kk].mul_vec(&x);
                    }
                }
            }
            acc.add_scaled_const(&x, coef);
        }
        acc.take()
    };
    let mut cols: Vec<PVec<F>> = vec![vec![]; fp * fa];
    for bp in 0..fp {
        let unit: PVec<F> = vec![(bp as u32, Poly::one())];
        for ba in 0..fa {
            let mut acc = PAcc::new();
            for (ba2, poly) in &t.phi.cols[ba] {
                let img = apply(poly, &unit);
                let placed: PVec<F> = img.into_iter().map(|(r, p)| ((r as usize * fa2 + *ba2 as usize) as u32, p)).collect();
                acc.add_scaled_const(&placed, &F::one());
            }
            cols[bp * fa + ba] = acc.take();
        }
    }
    let m = PMat { nrows: fp * fa2, ncols: fp * fa, cols };
    if suffix_fsize == 1 { m } else { m.kron_id(suffix_fsize) }
}

/// Polynomiality certificate (C3 H3, §3.4.5).
#[derive(Clone, Debug, Default)]
pub struct CertReport {
    pub points: usize,
    pub max_monos: usize,
    pub distinct_degrees: usize,
    pub entries: usize,
    pub redraws: usize,
    pub secs: f64,
}

/// The function y ↦ WinLoc (in sorted window-flow order) is the localization of a polynomial
/// bimodule map of normalized degree `deg`: Φ(y) = L_tgt(y)^{-1} f(y) L_src(y) has, at each entry
/// (β', β), a polynomial in the window bottom ring of degree D = u + fdeg_src[β] - fdeg_tgt[β'].
pub fn poly_certificate<F: Field>(src: &Window, tgt: &Window, deg: i32, f: &mut dyn FnMut(&[F]) -> WinLoc<F>, seed: u64) -> Result<CertReport, String> {
    let t0 = Instant::now();
    let sb: SBim<F> = src.expr().sbim();
    let tb: SBim<F> = tgt.expr().sbim();
    let u = deg - tb.shift + sb.shift;
    let sfp = fixed_points(&sb);
    let tfp = fixed_points(&tb);
    let sfi = FlowIndex::build(&src.expr());
    let tfi = FlowIndex::build(&tgt.expr());
    if sfi.top.mu != tfi.top.mu {
        return Err(format!("certificate: window tops differ ({:?} vs {:?})", sfi.top.mu, tfi.top.mu));
    }
    let spos = positions(&sfp, &sfi);
    let tpos = positions(&tfp, &tfi);
    // window top coset of every sorted window flow (right-linearity = block diagonality)
    let scos: Vec<u64> = (0..sfi.nflows() as u32).map(|f| sfi.coset_of_flow(f)).collect();
    let tcos: Vec<u64> = (0..tfi.nflows() as u32).map(|f| tfi.coset_of_flow(f)).collect();
    let (ns, nt) = (sfp.len(), tfp.len());
    let ring = sb.ring.clone();
    let kw: usize = src.labels.iter().map(|&a| a as usize).sum();
    // entries grouped by degree
    let mut by_d: std::collections::BTreeMap<i32, Vec<(usize, usize)>> = Default::default();
    let mut neg: Vec<(usize, usize, i32)> = vec![];
    for bp in 0..nt {
        for b in 0..ns {
            let d = u + sb.fdeg[b] - tb.fdeg[bp];
            if d < 0 {
                neg.push((bp, b, d));
            } else {
                by_d.entry(d).or_default().push((bp, b));
            }
        }
    }
    let monos: std::collections::BTreeMap<i32, Vec<u128>> = by_d.keys().map(|&d| (d, ring.monos_of_deg(d))).collect();
    let max_monos = monos.values().map(|v| v.len()).max().unwrap_or(0);
    let tt = 2 + max_monos;
    let npts = tt + 2;
    let off = crate::local::root_offsets(&src.labels.iter().map(|&a| a as usize).collect::<Vec<_>>());
    'redraw: for attempt in 0..=3usize {
        let pts = crate::smap::points::<F>(kw, seed.wrapping_add(attempt as u64), npts, 0);
        // values of Φ at every point (fibre coordinates, fixed-point order)
        let mut phis: Vec<Vec<Vec<F>>> = Vec::with_capacity(npts);
        let mut evals: Vec<Vec<Vec<F>>> = Vec::with_capacity(npts); // e_j of each window upright
        for (k, p) in pts.iter().enumerate() {
            let y = &p.x;
            let wl = f(y);
            if wl.nrows as usize != nt || wl.ncols() != ns {
                return Err(format!("certificate: WinLoc shape {}x{} != {}x{}", wl.nrows, wl.ncols(), nt, ns));
            }
            // block diagonality over the window top cosets (M4a review finding a): a bimodule
            // map commutes with the top operators, which are diagonal in the localized basis
            for j in 0..ns {
                for kk in wl.col_start[j] as usize..wl.col_start[j + 1] as usize {
                    let i = wl.rows[kk] as usize;
                    if !wl.vals[kk].is_zero() && tcos[i] != scos[j] {
                        return Err(format!(
                            "certificate: entry (tgt flow {}, src flow {}) = {} lies between different window top cosets ({} vs {}) at point {}: not right-linear (not a bimodule map)",
                            i, j, wl.vals[kk], tcos[i], scos[j], k
                        ));
                    }
                }
            }
            let d = wl.dense();
            let fd: Vec<Vec<F>> = (0..nt).map(|q| (0..ns).map(|pp| d[tpos[q] as usize][spos[pp] as usize].clone()).collect()).collect();
            let ls = loc_matrix(&sb, &sfp, y);
            let lt = loc_matrix(&tb, &tfp, y);
            let ltinv = inverse(&lt).ok_or("certificate: L_tgt singular")?;
            let phi = matmul(&matmul(&ltinv, &fd), &ls);
            for &(bp, b, dd) in &neg {
                if !phi[bp][b].is_zero() {
                    return Err(format!("certificate: entry ({}, {}) has degree {} < 0 but value {} at point {}", bp, b, dd, phi[bp][b], k));
                }
            }
            phis.push(phi);
            evals.push(src.labels.iter().enumerate().map(|(i, &l)| esym(crate::coset::low_mask(l as usize) << off[i], y)).collect());
        }
        let mono_val = |k: usize, m: u128| -> F {
            let mut v = F::one();
            let mut var = 0usize;
            for (i, &l) in src.labels.iter().enumerate() {
                for j in 0..l as usize {
                    for _ in 0..mono_exp(m, var) {
                        v = v.mul(&evals[k][i][j + 1]);
                    }
                    var += 1;
                }
            }
            v
        };
        for (&d, ents) in &by_d {
            let ms = &monos[&d];
            let nm = ms.len();
            let vd: Vec<Vec<F>> = (0..nm).map(|k| ms.iter().map(|&m| mono_val(k, m)).collect()).collect();
            let vinv = match inverse(&vd) {
                Some(v) => v,
                None => continue 'redraw,
            };
            // coefficients for every entry: c = V^{-1} vals
            let vals: Vec<Vec<F>> = (0..nm).map(|k| ents.iter().map(|&(bp, b)| phis[k][bp][b].clone()).collect()).collect();
            let coef = matmul(&vinv, &vals); // nm x ents
            for k in nm..npts {
                let row: Vec<F> = ms.iter().map(|&m| mono_val(k, m)).collect();
                for (e, &(bp, b)) in ents.iter().enumerate() {
                    let mut s = F::zero();
                    for (mi, r) in row.iter().enumerate() {
                        s.fma(r, &coef[mi][e]);
                    }
                    if s != phis[k][bp][b] {
                        return Err(format!(
                            "certificate: entry ({}, {}) (tgt fibre {}, src fibre {}) of degree {} is not a polynomial: residual mismatch at point {} ({} monomials)",
                            bp, b, bp, b, d, k, nm
                        ));
                    }
                }
            }
        }
        return Ok(CertReport { points: npts, max_monos, distinct_degrees: by_d.len(), entries: ns * nt, redraws: attempt, secs: t0.elapsed().as_secs_f64() });
    }
    Err("certificate: Vandermonde singular after 3 redraws".into())
}

#[cfg(test)]
mod m5c_tests {
    use super::*;
    use crate::field::F31;

    /// M5c: the templates of one hom_basis call share the localization sides of their windows
    /// (one SBim per window, not one clone per template); localizations are unchanged (A4 in
    /// tests/m4a.rs compares them with the oracle)
    #[test]
    fn hom_basis_templates_share_window_sides() {
        let a = Window::new(&[2, 2], &[Rung::e(0, 1), Rung::f(0, 1)]);
        let b = Window::new(&[2, 2], &[Rung::f(0, 1), Rung::e(0, 1)]);
        let mut lib: TmplLib<F31> = TmplLib::new(None);
        let mut ids = vec![];
        for d in 0..8 {
            ids = lib.hom_basis(&a, &b, d);
            if ids.len() >= 2 {
                break;
            }
        }
        assert!(ids.len() >= 2, "no degree with a hom space of dimension >= 2");
        for &t in &ids[1..] {
            assert!(Arc::ptr_eq(&lib.aux[t as usize].src, &lib.aux[ids[0] as usize].src));
            assert!(Arc::ptr_eq(&lib.aux[t as usize].tgt, &lib.aux[ids[0] as usize].tgt));
        }
        assert_eq!(lib.sides.len(), 2);
        // a second call on the same windows adds no side
        let _ = lib.hom_basis(&b, &a, 1);
        assert_eq!(lib.sides.len(), 2);
    }
}
