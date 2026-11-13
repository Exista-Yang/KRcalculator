//! M6a-2a: the checkpointed, resumable, shardable box driver (SPEC_3C_GL6 rev1 §2 "M6a-2a", D14).
//!
//! A *box job* is (λ, words, closable stages per word, box N, node seed, prune_n, modes, closure
//! N's) over one prime. [`box_run`] scans the shard's box points of every word (record scan at
//! the generic point first, frozen structure replayed at every box point) and streams the closure
//! entries into one `BoxAcc` per (word, stage, mode), checkpointed atomically in a run directory:
//!
//! ```text
//! DIR/run.hdr                      governing header (job, library hash, engine stamp, record
//!                                  traces and bindings, shard k/K); a resume must match it
//! DIR/w{i}_s{stage}_{u|r<m>}.acc   BoxAcc checkpoint (header with the binding, done-point
//!                                  bitmap, accumulators), each written atomically
//! DIR/commits.log                  the commit and stamps of every binary that ran in DIR
//! ```
//!
//! [`box_finish`] loads the checkpoints of every shard directory, merges them (refusing
//! mismatched headers, mismatched run headers and overlapping bitmaps), checks the complete
//! down-set, assembles `BoxClosureData`, records Freivalds `D2Cert`s, closes, checks χ against
//! `decat_chi` (fatal) and runs the D14 scaling check (mandatory; test-only opt-out).
//!
//! D14 items: atomic checkpoints with a header binding prime, nodes, library (layer-A stamps and
//! content through the record values), word, prune_n, mode, N, stage, skeleton, engine version
//! and record trace; a done-point bitmap; merge refusal; the down-set check in `finish`; the
//! generator-only determinant pass at every point of the shard before its first scan (and at
//! every needed point in `box_finish`); frozen structure (the record trace, pivot values
//! included, replayed at every point). Points are scanned in batches whose values are
//! accumulated as soon as the batch finishes (memory bounded by the batch size).
use crate::boxclose::*;
use crate::field::Field;
use crate::sengine::par_map;
use crate::swire::*;
use std::collections::BTreeMap;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

pub const BOXRUN_STAMP: &str = "boxrun-v1";

/// The commit of this binary (embedded at build time by build.rs; "-dirty" if the work tree
/// had changes then).
pub fn build_commit() -> &'static str {
    option_env!("KRHOM_BUILD_COMMIT").unwrap_or("unknown")
}

/// A box job (one prime: the field is the type parameter of the drivers).
#[derive(Clone, Debug)]
pub struct BoxJob {
    pub lam: Vec<u8>,
    pub words: Vec<Vec<i32>>,
    /// closed stages per word
    pub stages: Vec<Vec<usize>>,
    /// n_box, seed, n_close, modes, engine options, Freivalds vectors, scaling, batch (the
    /// `stages` field is ignored: per word above)
    pub o: BoxOpts,
}

impl BoxJob {
    fn word_opts(&self, wi: usize) -> BoxOpts {
        let mut o = self.o.clone();
        o.stages = self.stages[wi].clone();
        o
    }
}

/// Options of one `box_run` invocation (not part of the governing header).
#[derive(Clone, Debug)]
pub struct RunOpts {
    /// shard k of K (points of the needed list with index ≡ k mod K)
    pub shard: (usize, usize),
    /// save every checkpoint at least this often (seconds; 0 = after every batch)
    pub ckpt_secs: f64,
    /// TEST / demo: stop after scanning this many points in this invocation, without saving
    /// what was scanned since the last checkpoint (a simulated kill)
    pub stop_after: Option<usize>,
}

impl Default for RunOpts {
    fn default() -> Self {
        RunOpts { shard: (0, 1), ckpt_secs: 60.0, stop_after: None }
    }
}

/// What one `box_run` invocation did.
#[derive(Clone, Debug, Default)]
pub struct RunStatus {
    /// points of this shard (union over the words)
    pub shard_points: usize,
    /// points still pending at the start of this invocation
    pub pending_at_start: usize,
    /// points scanned in this invocation
    pub scanned: usize,
    /// points of the shard complete in every accumulator at the end
    pub complete: bool,
    /// determinant pass: (points, checks)
    pub det_pass: (usize, usize),
    pub resumed: bool,
    pub checkpoints_written: usize,
    pub secs: f64,
    /// per point cost of the scanned points
    pub costs: Vec<PointCost>,
}

/// Record scans of every word, with their stage bindings.
pub struct Records<F: Field> {
    pub recs: Vec<ScanOut<F>>,
    pub secs: f64,
}

pub fn records<F: PrimeField>(lib: &WLib<F>, job: &BoxJob) -> Result<Records<F>, String> {
    let t = Instant::now();
    let mut recs = vec![];
    for (wi, w) in job.words.iter().enumerate() {
        let o = job.word_opts(wi);
        let r = record_scan(lib, w, &o)?;
        if r.stages.iter().map(|s| s.stage).collect::<Vec<_>>() != job.stages[wi] {
            return Err(format!("word {:?}: record scan closed stages {:?}, job asks {:?}", w, r.stages.iter().map(|s| s.stage).collect::<Vec<_>>(), job.stages[wi]));
        }
        recs.push(r);
    }
    Ok(Records { recs, secs: t.elapsed().as_secs_f64() })
}

fn mode_str(m: Option<usize>) -> String {
    m.map_or("u".into(), |k| format!("r{}", k))
}

/// The governing header of a job (everything a checkpoint's meaning depends on; the shard line is
/// added by `box_run`). Library hash, engine stamp and the per-stage bindings come from the
/// library and the record scans.
pub fn governing_header<F: PrimeField>(lib: &WLib<F>, job: &BoxJob, recs: &Records<F>) -> String {
    let mut s = String::new();
    let w: Vec<String> = job.words.iter().map(|w| w.iter().map(|g| g.to_string()).collect::<Vec<_>>().join(" ")).collect();
    let st: Vec<String> = job.stages.iter().map(|v| v.iter().map(|g| g.to_string()).collect::<Vec<_>>().join(",")).collect();
    let jl = |v: &[usize]| v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",");
    s += &format!("krhom-boxrun v1\nstamp={}\nfield={}\nprime={}\nlam={}\nwords={}\nstages={}\n", BOXRUN_STAMP, F::name(), F::characteristic(), jl(&job.lam.iter().map(|&x| x as usize).collect::<Vec<_>>()), w.join(";"), st.join(";"));
    s += &format!("n_box={}\nseed={}\nprune_n={}\nmodes={}\nn_close={}\nfreivalds_k={}\nd2_all_points={}\ncheck_d2={}\n", job.o.n_box, job.o.seed, lib.prune_n.map_or("none".into(), |n| n.to_string()), job.o.modes.iter().map(|&m| mode_str(m)).collect::<Vec<_>>().join(","), jl(&job.o.n_close), job.o.freivalds_k, job.o.d2_all_points, job.o.eng.check_d2);
    s += &format!("library_hash={:x}\nengine_stamp={}\nbox_code_stamp={}\n", lib.hash, engine_stamp(&job.o.eng), BOX_CODE_STAMP);
    for (wi, r) in recs.recs.iter().enumerate() {
        for so in &r.stages {
            let b = stage_binding(so, &r.ps.frozen, &job.o.eng);
            s += &format!("bind w{} s{}: skeleton={:x} trace={:x} record={:x} entries={}\n", wi, so.stage, b.skeleton_hash, b.trace_hash, b.record_hash, so.kinds.len());
        }
    }
    s
}

fn write_atomic(path: &Path, text: &str) -> Result<(), String> {
    let tmp = path.with_extension(format!("tmp{}", std::process::id()));
    let r = (|| -> std::io::Result<()> {
        let mut f = std::fs::File::create(&tmp)?;
        f.write_all(text.as_bytes())?;
        f.sync_all()?;
        drop(f);
        std::fs::rename(&tmp, path)?;
        if let Some(d) = path.parent() {
            if let Ok(d) = std::fs::File::open(d) {
                let _ = d.sync_all();
            }
        }
        Ok(())
    })();
    r.map_err(|e| {
        let _ = std::fs::remove_file(&tmp);
        format!("{}: {}", path.display(), e)
    })
}

/// Parse "shard=k/K" of a run header; the governing part is everything else.
fn split_run_header(t: &str) -> Result<(String, (usize, usize)), String> {
    let mut gov = String::new();
    let mut shard = None;
    for l in t.lines() {
        if let Some(v) = l.strip_prefix("shard=") {
            let (k, n) = v.split_once('/').ok_or("run header: bad shard line")?;
            shard = Some((k.parse().map_err(|_| "run header: bad shard")?, n.parse().map_err(|_| "run header: bad shard")?));
        } else {
            gov += l;
            gov.push('\n');
        }
    }
    Ok((gov, shard.ok_or("run header: no shard line")?))
}

/// Per word: grids (per mode), the unreduced grid, the nodes.
struct Grids<F> {
    nodes: BoxNodes<F>,
    ug: BoxGrid<F>,
    grids: Vec<(Option<usize>, Arc<BoxGrid<F>>)>,
}

fn grids<F: PrimeField>(lam: &[u8], o: &BoxOpts) -> Grids<F> {
    let nodes = BoxNodes::<F>::new(lam, o.n_box, o.seed);
    let grids = o.modes.iter().map(|&m| (m, Arc::new(BoxGrid::new(nodes.clone(), m)))).collect();
    Grids { ug: BoxGrid::new(nodes.clone(), None), nodes, grids }
}

/// The expected checkpoint header of (word, stage, mode).
fn want_header<F: PrimeField>(lib: &WLib<F>, word: &[i32], so: &StageOut<F>, rec: &ScanOut<F>, g: &Arc<BoxGrid<F>>, eng: &crate::sengine::EngOpts) -> BoxAcc<F> {
    BoxAcc::with_binding(g.clone(), so.weights.clone(), lib.hash, word_hash(&[word.to_vec()]), lib.prune_n.map(|n| n as usize), stage_binding(so, &rec.ps.frozen, eng))
}

fn acc_name(wi: usize, stage: usize, m: Option<usize>) -> String {
    format!("w{}_s{}_{}.acc", wi, stage, mode_str(m))
}

/// Empty accumulators per word, stage and mode (in memory).
fn fresh_accs<F: PrimeField>(lib: &WLib<F>, job: &BoxJob, recs: &Records<F>, gs: &Grids<F>) -> Vec<Vec<Vec<BoxAcc<F>>>> {
    job.words.iter().zip(recs.recs.iter()).map(|(w, r)| r.stages.iter().map(|so| gs.grids.iter().map(|(_, g)| want_header(lib, w, so, r, g, &job.o.eng)).collect()).collect()).collect()
}

/// The needed points of every word, their union (sorted) and the shard's part of it.
fn shard_points<F: PrimeField>(accs: &[Vec<Vec<BoxAcc<F>>>], gs: &Grids<F>, modes: &[Option<usize>], shard: (usize, usize)) -> (Vec<u32>, Vec<u32>) {
    let mut all: Vec<u32> = vec![];
    for wa in accs {
        let maxsize = wa.iter().flatten().map(|a| a.needed_size().unwrap_or(0)).max().unwrap_or(0);
        all.extend(needed_points(&gs.ug, maxsize, modes));
    }
    all.sort();
    all.dedup();
    let mine = all.iter().enumerate().filter(|(i, _)| i % shard.1 == shard.0).map(|(_, &p)| p).collect();
    (all, mine)
}

/// Scan this shard's pending box points of every word, checkpointing atomically into `dir`
/// (created if missing). Resumes from the checkpoints in `dir` after verifying the governing
/// header and every checkpoint header. The generator-only determinant pass runs over the shard's
/// points before any scan.
pub fn box_run<F: PrimeField>(lib: &WLib<F>, job: &BoxJob, dir: &Path, ro: &RunOpts) -> Result<RunStatus, String> {
    let t0 = Instant::now();
    if ro.shard.1 == 0 || ro.shard.0 >= ro.shard.1 {
        return Err(format!("bad shard {}/{}", ro.shard.0, ro.shard.1));
    }
    let recs = records(lib, job)?;
    let gov = governing_header(lib, job, &recs);
    std::fs::create_dir_all(dir).map_err(|e| format!("{}: {}", dir.display(), e))?;
    let hp = dir.join("run.hdr");
    let mut st = RunStatus::default();
    if hp.exists() {
        let t = std::fs::read_to_string(&hp).map_err(|e| format!("{}: {}", hp.display(), e))?;
        let (g, sh) = split_run_header(&t)?;
        if g != gov {
            return Err(format!("{}: the run header does not match this job (resume refused):\n{}", hp.display(), line_diff(&g, &gov)));
        }
        if sh != ro.shard {
            return Err(format!("{}: directory belongs to shard {}/{}, not {}/{}", hp.display(), sh.0, sh.1, ro.shard.0, ro.shard.1));
        }
        st.resumed = true;
    } else {
        write_atomic(&hp, &format!("{}shard={}/{}\n", gov, ro.shard.0, ro.shard.1))?;
    }
    {
        let mut f = std::fs::OpenOptions::new().create(true).append(true).open(dir.join("commits.log")).map_err(|e| e.to_string())?;
        let _ = writeln!(f, "run start: commit={} engine={} shard={}/{}", build_commit(), engine_stamp(&job.o.eng), ro.shard.0, ro.shard.1);
    }
    let gs = grids::<F>(&lib.lam, &job.o);
    let mut accs = fresh_accs(lib, job, &recs, &gs);
    // load existing checkpoints (each against its full expected header)
    for (wi, wa) in accs.iter_mut().enumerate() {
        for (si, sa) in wa.iter_mut().enumerate() {
            for (gi, a) in sa.iter_mut().enumerate() {
                let p = dir.join(acc_name(wi, recs.recs[wi].stages[si].stage, gs.grids[gi].0));
                if p.exists() {
                    let want = a.header.clone();
                    *a = BoxAcc::load_expect(&p, gs.grids[gi].1.clone(), &want)?;
                    st.resumed = true;
                }
            }
        }
    }
    let (_, mine) = shard_points(&accs, &gs, &job.o.modes, ro.shard);
    st.shard_points = mine.len();
    // generator-only determinant pass at every point of the shard, before any scan
    let det = generator_pass(lib, &gs.nodes, &gs.ug, &mine)?;
    st.det_pass = (mine.len(), det);
    let pending: Vec<u32> = mine.iter().copied().filter(|&id| accs.iter().any(|wa| point_pending(wa, &gs.grids, &gs.ug, id))).collect();
    st.pending_at_start = pending.len();
    let popts: Vec<_> = (0..job.words.len()).map(|wi| replay_opts(lib, &job.word_opts(wi))).collect();
    let save_all = |accs: &[Vec<Vec<BoxAcc<F>>>]| -> Result<usize, String> {
        let mut n = 0;
        for (wi, wa) in accs.iter().enumerate() {
            for (si, sa) in wa.iter().enumerate() {
                for (gi, a) in sa.iter().enumerate() {
                    a.save(&dir.join(acc_name(wi, recs.recs[wi].stages[si].stage, gs.grids[gi].0)))?;
                    n += 1;
                }
            }
        }
        Ok(n)
    };
    if !st.resumed {
        st.checkpoints_written += save_all(&accs)?; // empty checkpoints: the run exists
    }
    let mut last_save = Instant::now();
    let batch = job.o.batch.max(1);
    let mut dirty = false;
    for chunk in pending.chunks(batch) {
        let chunk: Vec<u32> = match ro.stop_after {
            Some(m) if st.scanned + chunk.len() > m => chunk[..m - st.scanned].to_vec(),
            _ => chunk.to_vec(),
        };
        if chunk.is_empty() {
            break;
        }
        let res: Vec<Result<(u32, Vec<(Vec<Vec<F>>, PointCost)>), String>> = par_map(&chunk, |&id| {
            let mut per = vec![];
            for (wi, w) in job.words.iter().enumerate() {
                per.push(scan_box_point(lib, w, &gs.nodes, &gs.ug, id, &recs.recs[wi], &popts[wi])?);
            }
            Ok((id, per))
        });
        for r in res {
            let (id, per) = r?;
            for (wi, (vals, cost)) in per.into_iter().enumerate() {
                accumulate_point(&mut accs[wi], &gs.grids, &gs.ug, id, &vals, true)?;
                if wi == 0 {
                    st.costs.push(cost);
                }
            }
        }
        st.scanned += chunk.len();
        dirty = true;
        if ro.stop_after.is_some_and(|m| st.scanned >= m) {
            // simulated kill: what was scanned since the last checkpoint is lost
            st.secs = t0.elapsed().as_secs_f64();
            return Ok(st);
        }
        if last_save.elapsed().as_secs_f64() >= ro.ckpt_secs {
            st.checkpoints_written += save_all(&accs)?;
            last_save = Instant::now();
            dirty = false;
        }
    }
    if dirty {
        st.checkpoints_written += save_all(&accs)?;
    }
    st.complete = mine.iter().all(|&id| !accs.iter().any(|wa| point_pending(wa, &gs.grids, &gs.ug, id)));
    st.secs = t0.elapsed().as_secs_f64();
    Ok(st)
}

fn line_diff(a: &str, b: &str) -> String {
    let (la, lb): (Vec<&str>, Vec<&str>) = (a.lines().collect(), b.lines().collect());
    let mut d = vec![];
    for i in 0..la.len().max(lb.len()) {
        let (x, y) = (la.get(i).copied().unwrap_or(""), lb.get(i).copied().unwrap_or(""));
        if x != y {
            d.push(format!("  file: {}\n  job:  {}", x, y));
        }
    }
    d.join("\n")
}

/// One word's result of `box_finish`.
#[derive(Clone, Debug)]
pub struct WordResult {
    pub word: Vec<i32>,
    pub stages: Vec<StageResult>,
    pub chi_checked: usize,
    pub scaling: Option<ScalingOut>,
    /// digest of every merged accumulator (byte-identity of runs): (stage, mode, digest, points)
    pub digests: Vec<(usize, Option<usize>, u64, usize)>,
}

#[derive(Clone, Debug)]
pub struct FinishOut {
    pub words: Vec<WordResult>,
    pub shards: Vec<(usize, usize)>,
    pub det_pass: (usize, usize),
    pub record_secs: f64,
    pub close_secs: f64,
    pub governing: String,
}

/// Merge the shard directories, check completeness, close every stage of every word.
pub fn box_finish<F: PrimeField>(lib: &WLib<F>, job: &BoxJob, dirs: &[PathBuf]) -> Result<FinishOut, String> {
    if dirs.is_empty() {
        return Err("box_finish: no shard directory".into());
    }
    let recs = records(lib, job)?;
    let gov = governing_header(lib, job, &recs);
    let mut shards = vec![];
    for d in dirs {
        let hp = d.join("run.hdr");
        let t = std::fs::read_to_string(&hp).map_err(|e| format!("{}: {}", hp.display(), e))?;
        let (g, sh) = split_run_header(&t)?;
        if g != gov {
            return Err(format!("{}: run header does not match this job (merge refused):\n{}", hp.display(), line_diff(&g, &gov)));
        }
        if shards.iter().any(|&(k, n): &(usize, usize)| n != sh.1 || k == sh.0) {
            return Err(format!("{}: shard {}/{} duplicates or does not fit the shards {:?}", d.display(), sh.0, sh.1, shards));
        }
        shards.push(sh);
    }
    let gs = grids::<F>(&lib.lam, &job.o);
    let mut accs = fresh_accs(lib, job, &recs, &gs);
    for (wi, wa) in accs.iter_mut().enumerate() {
        for (si, sa) in wa.iter_mut().enumerate() {
            for (gi, a) in sa.iter_mut().enumerate() {
                let name = acc_name(wi, recs.recs[wi].stages[si].stage, gs.grids[gi].0);
                let want = a.header.clone();
                for d in dirs {
                    let s = BoxAcc::load_expect(&d.join(&name), gs.grids[gi].1.clone(), &want)?;
                    a.merge(&s).map_err(|e| format!("{} from {}: {}", name, d.display(), e))?;
                }
            }
        }
    }
    // determinant pass at every needed point (cheap; the shards ran it on their own points)
    let (all, _) = shard_points(&accs, &gs, &job.o.modes, (0, 1));
    let det = generator_pass(lib, &gs.nodes, &gs.ug, &all)?;
    let tc = Instant::now();
    let mut words = vec![];
    for (wi, w) in job.words.iter().enumerate() {
        let o = job.word_opts(wi);
        let rec = &recs.recs[wi];
        let mut wr = WordResult { word: w.clone(), stages: vec![], chi_checked: 0, scaling: None, digests: vec![] };
        for (si, so) in rec.stages.iter().enumerate() {
            for (gi, a) in accs[wi][si].iter().enumerate() {
                wr.digests.push((so.stage, gs.grids[gi].0, a.digest(), a.points_done()));
            }
            let (sr, chi) = close_stage_box(lib, w, so, &accs[wi][si], &gs.grids, &o)?;
            wr.chi_checked += chi;
            wr.stages.push(sr);
        }
        if o.scaling {
            let ts = Instant::now();
            let ids = scaling_points(&gs.ug);
            let rep = scaling_check_frozen(lib, w, &o, &F::from_u64(SCALING_LAMBDA), &ids, rec)?;
            wr.scaling = Some(ScalingOut { points: ids, lambda: SCALING_LAMBDA, report: rep, secs_milli: ts.elapsed().as_millis() as u64 });
        }
        words.push(wr);
    }
    Ok(FinishOut { words, shards, det_pass: (all.len(), det), record_secs: recs.secs, close_secs: tc.elapsed().as_secs_f64(), governing: gov })
}

/// Scan one box point (by unreduced grid id) of every word: per word and stage, the number of
/// entries and the hash of their values (debugging / spot checks).
pub fn box_point<F: PrimeField>(lib: &WLib<F>, job: &BoxJob, id: u32) -> Result<Vec<Vec<(usize, usize, u64)>>, String> {
    let recs = records(lib, job)?;
    let gs = grids::<F>(&lib.lam, &job.o);
    if id >= gs.ug.npoints {
        return Err(format!("point {} outside the box ({} points)", id, gs.ug.npoints));
    }
    generator_pass(lib, &gs.nodes, &gs.ug, &[id])?;
    let mut out = vec![];
    for (wi, w) in job.words.iter().enumerate() {
        let (vals, _) = scan_box_point(lib, w, &gs.nodes, &gs.ug, id, &recs.recs[wi], &replay_opts(lib, &job.word_opts(wi)))?;
        out.push(vals.iter().zip(recs.recs[wi].stages.iter()).map(|(v, so)| (so.stage, v.len(), values_hash(v))).collect());
    }
    Ok(out)
}

/// JSON of the closures of a finish (one record per (word, stage, N', mode)), with every
/// certification switch.
pub fn finish_json<F: PrimeField>(lib: &WLib<F>, job: &BoxJob, out: &FinishOut) -> String {
    let poly = |p: &crate::laurent::Poincare| {
        let terms: Vec<String> = p.c.iter().map(|(&(t, q), &k)| format!("[{},{},{}]", t, q, k)).collect();
        format!("{{\"str\":\"{}\",\"rank\":{},\"terms\":[{}]}}", p, p.total_rank(), terms.join(","))
    };
    let mut res = vec![];
    let mut scal = vec![];
    let mut dig = vec![];
    for wr in &out.words {
        for sr in &wr.stages {
            for (&(nn, marked), p) in &sr.closed {
                res.push(format!("{{\"word\":{:?},\"stage\":{},\"reduced\":{},\"N\":{},\"P\":{},\"entries\":{},\"objects\":{},\"d2_blocks\":{},\"chi_ok\":true}}", wr.word, sr.stage, marked.is_some(), nn, poly(p), sr.entries, sr.objects, sr.d2_blocks[&(nn, marked)]));
            }
        }
        scal.push(match &wr.scaling {
            Some(s) => format!("{{\"word\":{:?},\"points\":{:?},\"checked\":{},\"zero_both\":{},\"declared_zero\":{},\"pass\":true}}", wr.word, s.points, s.report.checked, s.report.zero_both, s.report.declared_zero),
            None => format!("{{\"word\":{:?},\"pass\":null,\"skipped\":\"test opt-out\"}}", wr.word),
        });
        for (st, m, d, n) in &wr.digests {
            dig.push(format!("{{\"word\":{:?},\"stage\":{},\"mode\":\"{}\",\"digest\":\"{:x}\",\"points\":{}}}", wr.word, st, mode_str(*m), d, n));
        }
    }
    format!(
        "{{\"tool\":\"sbox finish\",\"commit\":\"{}\",\"field\":\"{}\",\"lam\":{:?},\"n_box\":{},\"seed\":{},\"prune_n\":{},\"library_hash\":\"{:x}\",\"engine\":\"{}\",\"check_d2\":{},\"freivalds_k\":{},\"shards\":[{}],\"det_pass\":{{\"points\":{},\"checks\":{}}},\"chi_fatal\":true,\"scaling\":[{}],\"digests\":[{}],\"results\":[\n{}\n]}}\n",
        build_commit(), F::name(), job.lam, job.o.n_box, job.o.seed, lib.prune_n.map_or("null".into(), |n| n.to_string()), lib.hash, engine_stamp(&job.o.eng), job.o.eng.check_d2, job.o.freivalds_k, out.shards.iter().map(|(k, n)| format!("[{},{}]", k, n)).collect::<Vec<_>>().join(","), out.det_pass.0, out.det_pass.1, scal.join(","), dig.join(","), res.join(",\n")
    )
}

/// Closable stages of a word (labels back to λ).
pub fn closable_stages(lam: &[u8], word: &[i32]) -> Vec<usize> {
    let mut labels = lam.to_vec();
    let mut v = vec![];
    for (m, &g) in word.iter().enumerate() {
        let (i, _) = letter(g);
        labels.swap(i, i + 1);
        if labels == lam {
            v.push(m + 1);
        }
    }
    v
}

/// Summary map (word index, stage, N', marked) → Poincaré polynomial string (comparisons).
pub fn result_map(out: &FinishOut) -> BTreeMap<(usize, usize, usize, Option<usize>), String> {
    let mut m = BTreeMap::new();
    for (wi, wr) in out.words.iter().enumerate() {
        for sr in &wr.stages {
            for (&(nn, mk), p) in &sr.closed {
                m.insert((wi, sr.stage, nn, mk), p.to_string());
            }
        }
    }
    m
}
