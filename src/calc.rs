//! calc: the one-command calculator's library registry, installed-library checks, the per-prime
//! box run through the localized engine with a prebuilt layer-A library (3 and 4 strands), error
//! classification, library diagnostics, prime comparison, KnotJob-style braid parsing and JSON
//! helpers (CALC-1; docs/dev/colored_design/calc/CALC_SPEC_v1.md as amended by CALC_DECISIONS.md and
//! the review's required corrections).
//!
//! Certification is exactly that of `swire::box_word` with `BoxOpts::new` production defaults
//! (determinant pass, frozen-structure replay, Freivalds D², χ = decat_chi, D14 scaling); the
//! caller runs both primes and compares ([`compare_primes`]). Every layer-A map loaded from the
//! library is re-certified on load (`SIndecOpts::layer_a_route_stamp`, sindec.rs), so the
//! library's integrity checks (sha256) are about corruption, not soundness.
use crate::colored::boxclose::PrimeField;
use crate::algebra::laurent::Poincare;
use crate::colored::sengine::EngOpts;
use crate::library::sindec::{route_header, SIndec, SIndecOpts};
use crate::colored::swire::{box_word, build_wlib, record_scan, BoxOpts, NoProvider};
use std::collections::{BTreeMap, BTreeSet};
use std::io::BufRead;
use std::path::{Path, PathBuf};
use std::time::Instant;

/// The registry (single source of truth, shared with tools/fetch_library.sh).
pub const REGISTRY_TSV: &str = include_str!("../data/libraries.tsv");

// ------------------------------------------------------------------ errors

/// A calculator failure; `code()` is the process exit code.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CalcErr {
    /// usage or unsupported input (exit 1)
    Unsupported(String),
    /// engine or certification failure (exit 2)
    Engine(String),
    /// the two primes disagree (exit 3)
    Disagree(String),
    /// library not installed, wrong version, stamp mismatch, damaged, or not covering the word (exit 4)
    Library(String),
    /// resource refusal (exit 5)
    Resource(String),
}
impl CalcErr {
    pub fn code(&self) -> i32 {
        match self {
            CalcErr::Unsupported(_) => 1,
            CalcErr::Engine(_) => 2,
            CalcErr::Disagree(_) => 3,
            CalcErr::Library(_) => 4,
            CalcErr::Resource(_) => 5,
        }
    }
    pub fn msg(&self) -> &str {
        match self {
            CalcErr::Unsupported(m) | CalcErr::Engine(m) | CalcErr::Disagree(m) | CalcErr::Library(m) | CalcErr::Resource(m) => m,
        }
    }
    pub fn kind(&self) -> &'static str {
        match self {
            CalcErr::Unsupported(_) => "unsupported",
            CalcErr::Engine(_) => "engine",
            CalcErr::Disagree(_) => "disagree",
            CalcErr::Library(_) => "library",
            CalcErr::Resource(_) => "resource",
        }
    }
}

/// Class of an engine error string (review K2: the replay-import error never reaches the caller,
/// so there is no "smapio (replay)" row; a "not in layer A" error with 0 diagnosed missing jobs
/// is reported as "present but rejected on load" by the caller).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ErrClass {
    Library,
    Retry,
    Resource,
    Engine,
}
pub fn classify(err: &str) -> ErrClass {
    const LIB: [&str; 5] = ["not in layer A", "has no closure generators", "no closure data in the library", "loaded target rejected", "idempotent required (M5i)"];
    if LIB.iter().any(|s| err.contains(s)) {
        ErrClass::Library
    } else if err.contains("redraw the nodes") {
        ErrClass::Retry
    } else if err.contains("per-point memory") || err.contains("above the budget") || err.contains("above the memory cap") {
        ErrClass::Resource
    } else {
        ErrClass::Engine
    }
}

// ------------------------------------------------------------------ registry

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RegEntry {
    pub id: String,
    /// f31 | f32b
    pub field: String,
    pub lambda: Vec<u8>,
    pub n: u8,
    /// default | m4c-lf
    pub provider: String,
    pub data_commit: String,
    /// costs table path in the data commit
    pub costs: String,
    /// "complete" or a repo path of the status file (diagnostics only)
    pub status: String,
    pub stamp_layers: Vec<usize>,
    /// (path in the data commit, sha256), restore order
    pub layers: Vec<(String, String)>,
}
impl RegEntry {
    /// `F::name()` of the entry's prime
    pub fn field_name(&self) -> &'static str {
        field_name_of(&self.field)
    }
    pub fn layers_str(&self) -> String {
        self.layers.iter().map(|(p, h)| format!("{}@{}", p, h)).collect::<Vec<_>>().join(";")
    }
}
pub fn field_name_of(f: &str) -> &'static str {
    match f {
        "f31" => "F_2147483647",
        "f32b" => "F_4294967291",
        _ => "F_unknown",
    }
}
/// Approximate download size per prime (for messages).
pub fn fetch_size(id: &str) -> &'static str {
    match id {
        "l222_n4" => "about 3.5 MB per prime",
        "l2222_n4" => "about 16-18 MB per prime",
        _ => "size unknown",
    }
}

pub fn parse_registry(s: &str) -> Result<Vec<RegEntry>, String> {
    let mut v = vec![];
    for (ln, line) in s.lines().enumerate() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() != 10 {
            return Err(format!("registry line {}: {} columns, expected 10", ln + 1, f.len()));
        }
        let bad = |what: &str| format!("registry line {}: bad {}", ln + 1, what);
        let lambda: Vec<u8> = f[2].split(',').map(|x| x.parse().map_err(|_| bad("lambda"))).collect::<Result<_, _>>()?;
        let n: u8 = f[3].parse().map_err(|_| bad("n"))?;
        let stamp_layers: Vec<usize> = f[8].split(',').filter(|x| !x.is_empty()).map(|x| x.parse().map_err(|_| bad("stamp_layers"))).collect::<Result<_, _>>()?;
        let mut layers = vec![];
        for e in f[9].split(';') {
            let (p, h) = e.rsplit_once('@').ok_or_else(|| bad("layer"))?;
            if h.len() != 64 || !h.bytes().all(|b| b.is_ascii_hexdigit()) {
                return Err(bad("layer sha256"));
            }
            layers.push((p.to_string(), h.to_string()));
        }
        if f[5].len() != 40 {
            return Err(bad("data_commit"));
        }
        if !["f31", "f32b"].contains(&f[1]) {
            return Err(bad("field"));
        }
        v.push(RegEntry { id: f[0].into(), field: f[1].into(), lambda, n, provider: f[4].into(), data_commit: f[5].into(), costs: f[6].into(), status: f[7].into(), stamp_layers, layers });
    }
    Ok(v)
}
pub fn registry() -> Vec<RegEntry> {
    parse_registry(REGISTRY_TSV).expect("data/libraries.tsv")
}
/// The registry entry for labels λ and a prime (f31 | f32b).
pub fn reg_for(lam: &[u8], field: &str) -> Option<RegEntry> {
    registry().into_iter().find(|r| r.lambda == lam && r.field == field)
}
/// Supported (λ, N max) pairs, for messages.
pub fn supported_matrix() -> String {
    let mut seen = BTreeSet::new();
    let mut v = vec![];
    for r in registry() {
        if seen.insert(r.id.clone()) {
            v.push(format!("labels ({}) on {} strands, N <= {} (library {})", r.lambda.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","), r.lambda.len(), r.n, r.id));
        }
    }
    v.join("; ")
}

/// Library root: `--lib-dir` > `$KRHOM_LIB_DIR` > `${XDG_CACHE_HOME:-$HOME/.cache}/krhom/lib`.
pub fn lib_root(cli: Option<&str>) -> PathBuf {
    if let Some(d) = cli {
        return PathBuf::from(d);
    }
    if let Some(d) = std::env::var_os("KRHOM_LIB_DIR").filter(|d| !d.is_empty()) {
        return PathBuf::from(d);
    }
    let cache = std::env::var_os("XDG_CACHE_HOME").filter(|d| !d.is_empty()).map(PathBuf::from).unwrap_or_else(|| PathBuf::from(std::env::var_os("HOME").unwrap_or_default()).join(".cache"));
    cache.join("krhom").join("lib")
}

// ------------------------------------------------------------------ installed library

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LibFile {
    pub layer: usize,
    pub sha256: String,
    pub path: String,
}
#[derive(Clone, Debug)]
pub struct InstalledLib {
    pub dir: PathBuf,
    pub kv: BTreeMap<String, String>,
    pub files: Vec<LibFile>,
    pub costs: PathBuf,
    pub status: Option<PathBuf>,
}

/// Parse LIBRARY.txt: `key: value` lines up to `files:`, then `LAYER SHA256 PATH` lines.
pub fn parse_library_txt(s: &str) -> Result<(BTreeMap<String, String>, Vec<LibFile>), String> {
    let mut kv = BTreeMap::new();
    let mut files = vec![];
    let mut in_files = false;
    for line in s.lines() {
        if in_files {
            let f: Vec<&str> = line.split_whitespace().collect();
            if f.len() != 3 {
                return Err(format!("LIBRARY.txt: bad file line {:?}", line));
            }
            let layer = f[0].parse().map_err(|_| format!("LIBRARY.txt: bad layer in {:?}", line))?;
            if f[1].len() != 64 || f[2].contains("..") || f[2].starts_with('/') {
                return Err(format!("LIBRARY.txt: bad file line {:?}", line));
            }
            files.push(LibFile { layer, sha256: f[1].into(), path: f[2].into() });
        } else if line == "files:" {
            in_files = true;
        } else if let Some((k, v)) = line.split_once(": ") {
            kv.insert(k.to_string(), v.to_string());
        } else if !line.trim().is_empty() {
            return Err(format!("LIBRARY.txt: bad line {:?}", line));
        }
    }
    if kv.get("format").map(|s| s.as_str()) != Some("krhom-library v1") {
        return Err("LIBRARY.txt: not a krhom-library v1 file".into());
    }
    if !in_files {
        return Err("LIBRARY.txt: no file list".into());
    }
    Ok((kv, files))
}

/// Check LIBRARY.txt against the registry entry; returns the first mismatch.
pub fn check_library_kv(reg: &RegEntry, kv: &BTreeMap<String, String>, nfiles: usize) -> Result<(), String> {
    let lam = reg.lambda.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",");
    let want: [(&str, String); 7] = [("id", reg.id.clone()), ("lambda", lam), ("n", reg.n.to_string()), ("field", reg.field_name().to_string()), ("provider", reg.provider.clone()), ("data_commit", reg.data_commit.clone()), ("layers", reg.layers_str())];
    for (k, v) in want.iter() {
        let got = kv.get(*k).map(|s| s.as_str()).unwrap_or("(missing)");
        if got != v {
            return Err(format!("{} is {}, this krhom expects {}", k, if *k == "data_commit" { &got[..got.len().min(12)] } else { got }, if *k == "data_commit" { &v[..12] } else { v.as_str() }));
        }
    }
    let count = |k: &str| kv.get(k).and_then(|v| v.parse::<usize>().ok());
    match (count("files_layer_a"), count("files_tmpl")) {
        (Some(a), Some(t)) if a + t == nfiles => Ok(()),
        _ => Err(format!("file counts in the header do not match the {} listed files", nfiles)),
    }
}

/// Locate and check an installed library (review K3: every listed file must be present; extra
/// files, e.g. templates written by a later run, are tolerated; hashes are checked by
/// `tools/fetch_library.sh --verify-only`).
pub fn locate(reg: &RegEntry, root: &Path) -> Result<InstalledLib, CalcErr> {
    let dir = root.join(&reg.id).join(&reg.field);
    let lt = dir.join("LIBRARY.txt");
    let fetch = format!("tools/fetch_library.sh {}", reg.id);
    let text = std::fs::read_to_string(&lt).map_err(|_| CalcErr::Library(format!("library {} ({}) is not installed under {}: run {} ({})", reg.id, reg.field_name(), root.display(), fetch, fetch_size(&reg.id))))?;
    let force = format!("tools/fetch_library.sh --force {}", reg.id);
    let (kv, files) = parse_library_txt(&text).map_err(|e| CalcErr::Library(format!("installed {}/{} is damaged ({}): run {}", reg.id, reg.field, e, force)))?;
    check_library_kv(reg, &kv, files.len()).map_err(|e| CalcErr::Library(format!("installed {}/{} does not match this krhom ({}): run {}", reg.id, reg.field, e, force)))?;
    let missing: Vec<&LibFile> = files.iter().filter(|f| !dir.join(&f.path).is_file()).collect();
    if !missing.is_empty() {
        return Err(CalcErr::Library(format!("damaged install: {} of the {} files listed in {} are missing (e.g. {}): run {}", missing.len(), files.len(), lt.display(), missing[0].path, force)));
    }
    let costs = dir.join("costs.txt");
    let head = std::fs::File::open(&costs).ok().and_then(|f| std::io::BufReader::new(f).lines().next()).and_then(|l| l.ok()).unwrap_or_default();
    if !head.starts_with("krhom-costs v2 alg=hs-kl1") {
        return Err(CalcErr::Library(format!("damaged install: {} is missing or has the wrong format: run {}", costs.display(), force)));
    }
    let status = match kv.get("status").map(|s| s.as_str()) {
        Some("complete") | None => None,
        Some(s) => Some(dir.join(s)),
    };
    Ok(InstalledLib { dir, kv, files, costs, status })
}

/// First line of a file (without reading the rest).
fn first_line(p: &Path) -> Option<String> {
    let f = std::fs::File::open(p).ok()?;
    let mut s = String::new();
    std::io::BufReader::new(f).read_line(&mut s).ok()?;
    Some(s.trim_end_matches('\n').to_string())
}

/// Stamp preflight (review K4): every layer-A file whose last writer is a stamped layer must
/// match this binary's route stamp `want` (= `route_header(sd.layer_a_stamps())`); legacy files
/// of other layers are counted only. Returns (matching stamped files, non-matching legacy files).
pub fn stamp_preflight(reg: &RegEntry, inst: &InstalledLib, want: &str) -> Result<(usize, usize), CalcErr> {
    let want = route_header(want);
    let (mut ok, mut legacy, mut bad) = (0, 0, vec![]);
    for f in inst.files.iter().filter(|f| f.path.starts_with("layer_a/")) {
        let h = first_line(&inst.dir.join(&f.path)).unwrap_or_default();
        let prefix = h.split_once(" key=").map(|(a, _)| a).unwrap_or(&h);
        let matches = route_header(prefix) == want;
        if reg.stamp_layers.contains(&f.layer) {
            if matches {
                ok += 1;
            } else {
                bad.push(f.path.clone());
            }
        } else if !matches {
            legacy += 1;
        }
    }
    if !bad.is_empty() {
        return Err(CalcErr::Library(format!(
            "{} of the {} layer-A files of the stamped layers of {}/{} do not match this binary's route stamp (provider {}; e.g. {}): rebuild krhom from a branch that contains the c27 route_header filter (cca2db4), or reinstall with tools/fetch_library.sh --force {}",
            bad.len(),
            bad.len() + ok,
            reg.id,
            reg.field,
            reg.provider,
            bad[0],
            reg.id
        )));
    }
    if ok == 0 {
        return Err(CalcErr::Library(format!("no layer-A file of {}/{} matches this binary's route stamp (provider {})", reg.id, reg.field, reg.provider)));
    }
    Ok((ok, legacy))
}

// ------------------------------------------------------------------ the per-prime box run

/// Per-run parameters of the box route.
#[derive(Clone, Debug)]
pub struct BoxParams {
    pub word: Vec<i32>,
    /// N' values to close (each in [a, reg.n])
    pub ns: Vec<usize>,
    pub reduced: bool,
    pub seed: u64,
    pub retries: u64,
    pub mem_mb: u64,
    pub check_only: bool,
    /// smapio replay-only mode (Q3: on by default; a template recompute is a retry)
    pub replay_only: bool,
}

#[derive(Clone, Debug, Default)]
pub struct PrimeOut {
    pub field: String,
    pub seed: u64,
    pub closed: BTreeMap<usize, Poincare>,
    pub secs: f64,
    pub lib_build_secs: f64,
    pub points: usize,
    pub record_secs: f64,
    pub scan_secs: f64,
    pub close_secs: f64,
    pub chi_checked: usize,
    pub closures: usize,
    pub scaling: (usize, usize, usize),
    pub det_pass: (usize, usize),
    pub layer_a: (usize, usize, usize),
    pub provider_calls: usize,
    pub stamp_ok: usize,
    pub stamp_legacy: usize,
    pub vmhwm_kb: usize,
    pub template_recompute: bool,
    pub check_only: bool,
    /// --check-library: whether a record scan was needed (some closable type had no closure data)
    pub record_scan: bool,
    pub lib_types: usize,
}

fn map_engine(field: &str, e: String) -> CalcErr {
    match classify(&e) {
        ErrClass::Library => CalcErr::Library(e),
        ErrClass::Resource => CalcErr::Resource(e),
        _ => CalcErr::Engine(format!("localized engine failed on {}: {}. This is a certification failure; please report it with --json output.", field, e)),
    }
}

/// One prime of the box route: stamp, replay-only library build from layer A, preflight, then
/// the record scan only (`check_only`) or the full `box_word` with its certification, and the
/// post-checks of review K1 (χ checked for every closure; the requested (N', mode) present).
pub fn run_prime<F: PrimeField>(reg: &RegEntry, inst: &InstalledLib, p: &BoxParams) -> Result<PrimeOut, CalcErr> {
    let t0 = Instant::now();
    let lam = reg.lambda.clone();
    // --checkpoint (src/calc_ckpt.rs): check or create this prime's manifest before any work
    let ckpt = match crate::calc_ckpt::dir() {
        Some(root) if !p.check_only => Some(crate::calc_ckpt::prepare(&root, reg, &F::name(), p)?.0),
        _ => None,
    };
    let (_, stamp) = match reg.provider.as_str() {
        "default" => crate::library::sindec::provider_for_stamped::<F>(&lam, Some(&inst.costs), 1.0),
        "m4c-lf" => crate::library::leaves::provider_for_stamped_m4c::<F>(&lam, Some(&inst.costs), 1.0, true),
        pv => return Err(CalcErr::Library(format!("unknown provider {} in the registry", pv))),
    };
    crate::library::smapio::set_replay_only(p.replay_only);
    let sd: SIndec<F> = SIndec::with_stamp(&lam, Box::new(NoProvider), stamp, SIndecOpts { cache_dir: Some(inst.dir.clone()), prune_n: Some(reg.n), layer_a_route_stamp: true, ..Default::default() });
    let (stamp_ok, stamp_legacy) = stamp_preflight(reg, inst, &sd.layer_a_stamps())?;
    let lib = build_wlib(sd, &[p.word.clone()], true).map_err(|e| map_engine(&F::name(), format!("library: {}", e)))?;
    let mut out = PrimeOut { field: F::name().to_string(), seed: p.seed, lib_build_secs: lib.build_secs, layer_a: (lib.stats.layer_a_loaded, lib.stats.layer_a_route, lib.stats.layer_a_rejected), provider_calls: lib.stats.provider_calls, stamp_ok, stamp_legacy, template_recompute: !p.replay_only, check_only: p.check_only, lib_types: lib.shapes.len(), ..Default::default() };
    let eng = EngOpts { check_d2: false, mem_point_bytes: p.mem_mb << 20, ..Default::default() };
    let mode = if p.reduced { Some(0) } else { None };
    let mut seed = p.seed;
    if p.check_only {
        // review K10: closure data for every closable type proves coverage; otherwise the record
        // scan finds out whether a final type lacks it
        if !lib.no_closure.is_empty() {
            out.record_scan = true;
            let o = BoxOpts::new(reg.n as usize, seed, p.ns.clone(), vec![mode], vec![p.word.len()], eng);
            record_scan(&lib, &p.word, &o).map_err(|e| map_engine(&F::name(), e))?;
        }
        out.secs = t0.elapsed().as_secs_f64();
        out.vmhwm_kb = crate::colored::xengine::rss_kb().1;
        return Ok(out);
    }
    let bo = loop {
        let o = BoxOpts::new(reg.n as usize, seed, p.ns.clone(), vec![mode], vec![p.word.len()], eng.clone());
        let r = match &ckpt {
            Some(d) => crate::calc_ckpt::box_word_ckpt(&lib, &p.word, &o, d),
            None => box_word(&lib, &p.word, &o),
        };
        match r {
            Ok(bo) => break bo,
            Err(e) if classify(&e) == ErrClass::Retry && seed < p.seed + p.retries => {
                eprintln!("[{}] {}; retrying with seed {}", F::name(), e, seed + 1);
                seed += 1;
            }
            Err(e) => return Err(map_engine(&F::name(), e)),
        }
    };
    // post-checks (review K1)
    let post = |m: String| CalcErr::Engine(format!("localized engine on {}: post-check failed: {}", F::name(), m));
    if bo.stages.len() != 1 || bo.stages[0].stage != p.word.len() {
        return Err(post(format!("expected one stage {}, got {:?}", p.word.len(), bo.stages.iter().map(|s| s.stage).collect::<Vec<_>>())));
    }
    let sc = bo.scaling.as_ref().ok_or_else(|| post("no scaling check".into()))?;
    let nclosed: usize = bo.stages.iter().map(|s| s.closed.len()).sum();
    if bo.chi_checked != nclosed {
        return Err(post(format!("chi checked for {} of {} closures", bo.chi_checked, nclosed)));
    }
    for &nn in &p.ns {
        let pc = bo.stages[0].closed.get(&(nn, mode)).ok_or_else(|| post(format!("no closure for N' = {} mode {:?}", nn, mode)))?;
        out.closed.insert(nn, pc.clone());
    }
    out.seed = seed;
    out.points = bo.points;
    out.record_secs = bo.record_secs;
    out.scan_secs = bo.scan_secs;
    out.close_secs = bo.close_secs;
    out.chi_checked = bo.chi_checked;
    out.closures = nclosed;
    out.scaling = (sc.report.checked, sc.report.zero_both, sc.report.declared_zero);
    out.det_pass = bo.det_pass;
    out.secs = t0.elapsed().as_secs_f64();
    out.vmhwm_kb = crate::colored::xengine::rss_kb().1;
    Ok(out)
}

/// `--box-live` (tests only, review F2): the box route on at most 2 strands with the live
/// provider and an in-memory, unpruned library (`check_d2` on), box N = max(ns), as
/// tests/swire.rs does. Never replay-only.
pub fn run_prime_live<F: PrimeField>(lam: &[u8], costs: &Path, word: &[i32], ns: &[usize], reduced: bool, seed: u64) -> Result<PrimeOut, CalcErr> {
    let t0 = Instant::now();
    if lam.len() > 2 {
        return Err(CalcErr::Unsupported("--box-live handles at most 2 strands".into()));
    }
    let n_box = *ns.iter().max().ok_or_else(|| CalcErr::Unsupported("--box-live: no N".into()))?;
    if lam.iter().any(|&a| a as usize > n_box) {
        return Err(CalcErr::Unsupported("--box-live: N below a color".into()));
    }
    crate::library::smapio::set_replay_only(false);
    let (prov, stamp) = crate::library::sindec::provider_for_stamped::<F>(lam, Some(costs), 1.0);
    let sd: SIndec<F> = SIndec::with_stamp(lam, prov, stamp, SIndecOpts { prune_n: None, ..Default::default() });
    let lib = build_wlib(sd, &[word.to_vec()], true).map_err(|e| map_engine(&F::name(), format!("library: {}", e)))?;
    let mode = if reduced { Some(0) } else { None };
    let o = BoxOpts::new(n_box, seed, ns.to_vec(), vec![mode], vec![word.len()], EngOpts { check_d2: true, ..Default::default() });
    let bo = box_word(&lib, word, &o).map_err(|e| map_engine(&F::name(), e))?;
    let nclosed: usize = bo.stages.iter().map(|s| s.closed.len()).sum();
    if bo.chi_checked != nclosed || bo.scaling.is_none() || bo.stages.len() != 1 {
        return Err(CalcErr::Engine("--box-live: post-check failed".into()));
    }
    let mut out = PrimeOut { field: F::name().to_string(), seed, lib_build_secs: lib.build_secs, points: bo.points, record_secs: bo.record_secs, scan_secs: bo.scan_secs, close_secs: bo.close_secs, chi_checked: bo.chi_checked, closures: nclosed, det_pass: bo.det_pass, provider_calls: lib.stats.provider_calls, lib_types: lib.shapes.len(), ..Default::default() };
    let sc = bo.scaling.as_ref().unwrap();
    out.scaling = (sc.report.checked, sc.report.zero_both, sc.report.declared_zero);
    for &nn in ns {
        out.closed.insert(nn, bo.stages[0].closed.get(&(nn, mode)).cloned().ok_or_else(|| CalcErr::Engine(format!("--box-live: no closure for N' = {}", nn)))?);
    }
    out.secs = t0.elapsed().as_secs_f64();
    out.vmhwm_kb = crate::colored::xengine::rss_kb().1;
    Ok(out)
}

/// Library jobs the word needs that are not `ok` in both primes' status files (spec C4
/// `diagnose_missing`; logic of examples/m5i_cover.rs). Keys `pair K (p=TYPE, rung=R)` /
/// `closure TABLE`.
pub fn diagnose_missing(lam: &[u8], n: u8, word: &[i32], status: &[PathBuf]) -> Result<Vec<String>, String> {
    let mut oks: Vec<BTreeSet<String>> = vec![];
    for s in status {
        let t = std::fs::read_to_string(s).map_err(|e| format!("{}: {}", s.display(), e))?;
        oks.push(t.lines().filter_map(|l| l.strip_suffix(" ok")).map(|l| l.to_string()).collect());
    }
    let ok = |k: &str| oks.iter().all(|o| o.contains(k));
    let mut s = crate::combinat::schur::Schur::new(lam);
    let plan = crate::library::sindec::pruned_plan(&mut s, Some(n));
    let kj = crate::library::sindec::knot_jobs(&mut s, &plan, word)?;
    let mut miss = vec![];
    for &k in &kj.pairs {
        let key = format!("pair {}", k);
        if !ok(&key) {
            let pp = &plan.pairs[k];
            miss.push(format!("{} (p={}, rung={})", key, s.table(pp.p).fmt(), crate::library::smapio::rungs_str(&[pp.rung])));
        }
    }
    for &q in &kj.closure {
        let w = crate::library::sindec::canonical_ladder(&mut s, q);
        if w.is_empty() {
            continue;
        }
        let key = format!("closure {}", s.table(q).fmt());
        if !ok(&key) {
            miss.push(key);
        }
    }
    Ok(miss)
}

/// Compare the two primes' results exactly (shared by both singular routes and the Soergel and
/// cable routes); Err lists the differing N.
pub fn compare_primes(ns: &[usize], a: &BTreeMap<usize, Poincare>, b: &BTreeMap<usize, Poincare>) -> Result<(), String> {
    let mut diffs = vec![];
    for &nn in ns {
        if a.get(&nn) != b.get(&nn) {
            diffs.push(format!("N = {}: F_2147483647 gives {}, F_4294967291 gives {}", nn, a.get(&nn).map_or("nothing".into(), |p| p.to_string()), b.get(&nn).map_or("nothing".into(), |p| p.to_string())));
        }
    }
    if diffs.is_empty() {
        Ok(())
    } else {
        Err(format!("the two primes disagree (torsion or a bug); no result printed:\n  {}", diffs.join("\n  ")))
    }
}

// ------------------------------------------------------------------ input helpers

/// Parse a braid word: signed integers ("1 -2 1 -2", commas allowed) or KnotJob letters
/// ("aBaB": a..i = σ1..σ9, A..I their inverses, letters only).
pub fn parse_braid(s: &str) -> Result<Vec<i64>, String> {
    let t = s.trim();
    if t.is_empty() {
        return Ok(vec![]);
    }
    let has_letter = t.chars().any(|c| c.is_ascii_alphabetic());
    let has_digit = t.chars().any(|c| c.is_ascii_digit());
    if has_letter && has_digit {
        return Err(format!("braid {:?} mixes letters and integers: use either signed integers (\"1 -2 1 -2\") or KnotJob letters (\"aBaB\")", s));
    }
    if has_letter {
        let mut w = vec![];
        for c in t.chars() {
            match c {
                'a'..='i' => w.push((c as i64) - ('a' as i64) + 1),
                'A'..='I' => w.push(-((c as i64) - ('A' as i64) + 1)),
                _ if c.is_whitespace() => return Err(format!("KnotJob letter braid {:?}: letters only, no spaces", s)),
                _ => return Err(format!("braid {:?}: letter {:?} is not a KnotJob generator (a..i = sigma_1..sigma_9, A..I their inverses)", s, c)),
            }
        }
        return Ok(w);
    }
    let mut w = vec![];
    for tok in t.split(|c: char| c.is_whitespace() || c == ',').filter(|x| !x.is_empty()) {
        let k: i64 = tok.parse().map_err(|_| format!("braid {:?}: {:?} is not a signed integer generator", s, tok))?;
        if k == 0 {
            return Err(format!("braid {:?}: generator 0 does not exist (generators are 1, -1, 2, -2, ...)", s));
        }
        if k.unsigned_abs() > 1000 {
            return Err(format!("braid {:?}: generator {} is too large", s, k));
        }
        w.push(k);
    }
    Ok(w)
}

/// The simplify guard of the singular dispatch (spec B step 3): use the simplified braid only if
/// it has fewer strands and, for reduced homology of a link, the marked component (the one
/// through bottom position 0) is unchanged.
pub fn use_simplified(orig_strands: usize, orig_comp0: usize, s_strands: usize, s_comp0: usize, reduced: bool, ncomp: usize) -> bool {
    s_strands < orig_strands && (!reduced || ncomp == 1 || s_comp0 == orig_comp0)
}

/// JSON string literal with escaping (review K13).
pub fn jstr(s: &str) -> String {
    let mut o = String::with_capacity(s.len() + 2);
    o.push('"');
    for c in s.chars() {
        match c {
            '"' => o.push_str("\\\""),
            '\\' => o.push_str("\\\\"),
            '\n' => o.push_str("\\n"),
            '\r' => o.push_str("\\r"),
            '\t' => o.push_str("\\t"),
            c if (c as u32) < 0x20 => o.push_str(&format!("\\u{:04x}", c as u32)),
            c => o.push(c),
        }
    }
    o.push('"');
    o
}

/// A Poincaré polynomial as JSON: string, rank and terms [t, q, k].
pub fn poly_json(p: &Poincare) -> String {
    let terms: Vec<String> = p.c.iter().map(|(&(t, q), &k)| format!("[{},{},{}]", t, q, k)).collect();
    format!("{{\"str\":{},\"rank\":{},\"terms\":[{}]}}", jstr(&p.to_string()), p.total_rank(), terms.join(","))
}

/// Per-prime provenance as JSON.
pub fn prime_json(p: &PrimeOut) -> String {
    format!(
        "{{\"field\":{},\"seed\":{},\"secs\":{:.2},\"lib_build_secs\":{:.2},\"points\":{},\"record_secs\":{:.2},\"scan_secs\":{:.2},\"close_secs\":{:.2},\"chi_checked\":{},\"closures\":{},\"scaling\":{{\"checked\":{},\"zero_both\":{},\"declared_zero\":{}}},\"det_pass\":[{},{}],\"layer_a\":{{\"loaded\":{},\"route\":{},\"rejected\":{}}},\"provider_calls\":{},\"stamp_preflight\":{{\"stamped_ok\":{},\"legacy_other\":{}}},\"template_recompute\":{},\"check_only\":{},\"record_scan\":{},\"lib_types\":{},\"vmhwm_kb\":{}}}",
        jstr(&p.field),
        p.seed,
        p.secs,
        p.lib_build_secs,
        p.points,
        p.record_secs,
        p.scan_secs,
        p.close_secs,
        p.chi_checked,
        p.closures,
        p.scaling.0,
        p.scaling.1,
        p.scaling.2,
        p.det_pass.0,
        p.det_pass.1,
        p.layer_a.0,
        p.layer_a.1,
        p.layer_a.2,
        p.provider_calls,
        p.stamp_ok,
        p.stamp_legacy,
        p.template_recompute,
        p.check_only,
        p.record_scan,
        p.lib_types,
        p.vmhwm_kb
    )
}

/// Environment switches that change the computation on the box route (review K12): refused.
pub const REFUSED_ENV: [&str; 7] = ["KRHOM_FORCE_LADDER", "KRHOM_NO_M5I", "KRHOM_DCID_REF", "KRHOM_M5I_EXACT_UNK", "KRHOM_M5I_NO_PREFIX", "KRHOM_M5I_NO_REROUTE", "KRHOM_NO_REROUTE_FALLBACK"];
/// Refused switches that are set (KRHOM_LEAVES_C27 also, when set to anything but empty).
pub fn refused_env_set() -> Vec<String> {
    let mut v: Vec<String> = REFUSED_ENV.iter().filter(|k| std::env::var_os(k).is_some()).map(|k| k.to_string()).collect();
    if std::env::var_os("KRHOM_LEAVES_C27").is_some_and(|v| !v.is_empty()) {
        v.push("KRHOM_LEAVES_C27".into());
    }
    v
}
/// All KRHOM_* variables (recorded in --json).
pub fn krhom_env() -> Vec<(String, String)> {
    let mut v: Vec<(String, String)> = std::env::vars().filter(|(k, _)| k.starts_with("KRHOM_")).collect();
    v.sort();
    v
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn registry_parses() {
        let r = registry();
        assert_eq!(r.len(), 4);
        assert_eq!(r.iter().map(|e| format!("{}/{}", e.id, e.field)).collect::<Vec<_>>(), ["l222_n4/f31", "l222_n4/f32b", "l2222_n4/f31", "l2222_n4/f32b"]);
        for e in &r {
            assert_eq!(e.data_commit, "ee548d2bea8724f88660897996678c76eab1a43d");
            assert!(e.layers.iter().all(|(_, h)| h.len() == 64));
            assert_eq!(e.n, 4);
        }
        let l4 = reg_for(&[2, 2, 2, 2], "f31").unwrap();
        assert_eq!(l4.provider, "m4c-lf");
        assert_eq!(l4.stamp_layers, vec![1, 2, 3]);
        let order: Vec<&str> = l4.layers.iter().map(|(p, _)| p.split('/').nth(1).unwrap()).collect();
        assert_eq!(order, ["l2222_n4", "l2222_n4_c25", "l2222_n4_c24", "l2222_n4_c27"]);
        assert!(l4.layers[3].1.starts_with("dcb4dda2"));
        assert!(reg_for(&[3, 3, 3], "f31").is_none());
        assert!(parse_registry("a\tb\n").is_err());
    }

    fn sample_lib(reg: &RegEntry) -> String {
        format!(
            "format: krhom-library v1\nid: {}\nlambda: 2,2,2\nn: 4\nfield: F_2147483647\nprovider: default\ndata_commit: {}\nlayers: {}\nstamp_layers: 0\ncosts: costs.txt {}\nstatus: complete\nfiles_layer_a: 1\nfiles_tmpl: 1\nfiles:\n0 {} layer_a/a-0123456789abcdef.txt\n0 {} tmpl/tmpl-0123456789abcdef.txt\n",
            reg.id,
            reg.data_commit,
            reg.layers_str(),
            "0".repeat(64),
            "a".repeat(64),
            "b".repeat(64)
        )
    }

    #[test]
    fn library_txt_parse_and_mismatch() {
        let reg = reg_for(&[2, 2, 2], "f31").unwrap();
        let s = sample_lib(&reg);
        let (kv, files) = parse_library_txt(&s).unwrap();
        assert_eq!(files.len(), 2);
        assert_eq!(files[1], LibFile { layer: 0, sha256: "b".repeat(64), path: "tmpl/tmpl-0123456789abcdef.txt".into() });
        check_library_kv(&reg, &kv, files.len()).unwrap();
        // wrong prime, wrong data commit, wrong count
        let reg32 = reg_for(&[2, 2, 2], "f32b").unwrap();
        assert!(check_library_kv(&reg32, &kv, 2).unwrap_err().contains("field"));
        let s2 = s.replace(&reg.data_commit, &"1".repeat(40));
        let (kv2, _) = parse_library_txt(&s2).unwrap();
        assert!(check_library_kv(&reg, &kv2, 2).unwrap_err().contains("data_commit is 111111111111, this krhom expects ee548d2bea87"));
        assert!(check_library_kv(&reg, &kv, 3).is_err());
        assert!(parse_library_txt("format: other\nfiles:\n").is_err());
        assert!(parse_library_txt(&s.replace("layer_a/a-", "../a-")).is_err());
    }

    #[test]
    fn locate_reports_missing_install_and_damage() {
        let reg = reg_for(&[2, 2, 2], "f31").unwrap();
        let root = std::env::temp_dir().join(format!("krhom_calc_locate_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        let e = locate(&reg, &root).unwrap_err();
        assert_eq!(e.code(), 4);
        assert!(e.msg().contains("tools/fetch_library.sh l222_n4"), "{}", e.msg());
        let dir = root.join("l222_n4/f31");
        std::fs::create_dir_all(dir.join("layer_a")).unwrap();
        std::fs::create_dir_all(dir.join("tmpl")).unwrap();
        std::fs::write(dir.join("LIBRARY.txt"), sample_lib(&reg)).unwrap();
        std::fs::write(dir.join("costs.txt"), "krhom-costs v2 alg=hs-kl1 host=x\n").unwrap();
        std::fs::write(dir.join("layer_a/a-0123456789abcdef.txt"), "x\n").unwrap();
        let e = locate(&reg, &root).unwrap_err();
        assert!(e.msg().contains("damaged install") && e.msg().contains("--force l222_n4"), "{}", e.msg());
        std::fs::write(dir.join("tmpl/tmpl-0123456789abcdef.txt"), "y\n").unwrap();
        std::fs::write(dir.join("tmpl/tmpl-ffffffffffffffff.txt"), "extra\n").unwrap();
        let inst = locate(&reg, &root).unwrap();
        assert!(inst.status.is_none());
        let _ = std::fs::remove_dir_all(&root);
    }

    #[test]
    fn classify_table() {
        assert_eq!(classify("library: dec_word(..): pair ...: replay: pairs not in layer A"), ErrClass::Library);
        assert_eq!(classify("final type 0,2/2,0 has no closure generators: x"), ErrClass::Library);
        assert_eq!(classify("no closure data in the library"), ErrClass::Library);
        assert_eq!(classify("cache: loaded target rejected: bad"), ErrClass::Library);
        assert_eq!(classify("type 1,1/1,1: canonical ladder decomposable, idempotent required (M5i)"), ErrClass::Library);
        assert_eq!(classify("generator determinant fails at box point [..] for type t (ins: 3 vectors of length 3): redraw the nodes (another --seed) before any scan"), ErrClass::Retry);
        assert_eq!(classify("per-point memory 12.0 MB above the budget 10.0 MB at stage 3: run with streaming or on the cloud"), ErrClass::Resource);
        assert_eq!(classify("staged rank stage for [2, 2] at N = 4 needs about 9 MB, above the memory cap of 1 MB (set KRHOM_MEM_CAP_MB to raise it)"), ErrClass::Resource);
        assert_eq!(classify("CHI FAIL stage 4 N=3 marked None"), ErrClass::Engine);
        assert_eq!(classify("scaling check: entry 3 weight 2"), ErrClass::Engine);
    }

    #[test]
    fn compare_primes_agree_and_disagree() {
        let mut a = BTreeMap::new();
        let mut p = Poincare::default();
        p.add_term(0, 0, 1);
        a.insert(2, p.clone());
        let b = a.clone();
        assert!(compare_primes(&[2], &a, &b).is_ok());
        let mut c = a.clone();
        let mut q = p.clone();
        q.add_term(1, 2, 1);
        c.insert(2, q);
        let e = compare_primes(&[2], &a, &c).unwrap_err();
        assert!(e.contains("N = 2") && e.contains("disagree"), "{}", e);
    }

    #[test]
    fn braid_parsing() {
        assert_eq!(parse_braid("1 -2 1 -2").unwrap(), vec![1, -2, 1, -2]);
        assert_eq!(parse_braid("1,-2, 1").unwrap(), vec![1, -2, 1]);
        assert_eq!(parse_braid("aBaB").unwrap(), vec![1, -2, 1, -2]);
        assert_eq!(parse_braid("aaa").unwrap(), vec![1, 1, 1]);
        assert_eq!(parse_braid("iI").unwrap(), vec![9, -9]);
        assert!(parse_braid("a1").unwrap_err().contains("mixes"));
        assert!(parse_braid("a b").unwrap_err().contains("no spaces"));
        assert!(parse_braid("aj").unwrap_err().contains("not a KnotJob generator"));
        assert!(parse_braid("1 0").unwrap_err().contains("generator 0"));
        assert!(parse_braid("1 x").is_err());
        assert_eq!(parse_braid("  ").unwrap(), Vec::<i64>::new());
    }

    #[test]
    fn reduced_link_simplify_guard() {
        assert!(use_simplified(3, 0, 2, 1, false, 2));
        assert!(use_simplified(3, 0, 2, 1, true, 1));
        assert!(!use_simplified(3, 0, 2, 1, true, 2));
        assert!(use_simplified(3, 0, 2, 0, true, 2));
        assert!(!use_simplified(3, 0, 3, 0, false, 1));
    }

    #[test]
    fn json_escaping() {
        assert_eq!(jstr("a\"b\\c\nd"), "\"a\\\"b\\\\c\\nd\"");
        let mut p = Poincare::default();
        p.add_term(-1, 2, 3);
        assert_eq!(poly_json(&p), "{\"str\":\"3t^-1q^2\",\"rank\":3,\"terms\":[[-1,2,3]]}");
    }
}
