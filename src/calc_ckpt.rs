//! calc_ckpt: `krhom --checkpoint DIR`, the resumable box route (F-CKPT).
//!
//! The 3-4 strand box route of [`crate::calc::run_prime`] normally scans every box point in
//! memory (`swire::box_word`): an interrupted run loses everything. With a checkpoint directory
//! the same computation goes through the checkpointed driver of `colored::boxrun` instead
//! ([`box_run_with`] then [`box_finish_with`], the machinery behind `sbox box|finish`; no new
//! checkpoint format):
//!
//! ```text
//! DIR/<f31|f32b>/manifest.txt      the full input of this prime (see `manifest`); a rerun with
//!                                  a different input is refused (exit 1)
//! DIR/<f31|f32b>/seed<S>/run.hdr   boxrun's governing header (library hash, engine stamp,
//!                                  record traces, ...); boxrun refuses a mismatch itself
//! DIR/<f31|f32b>/seed<S>/*.acc     boxrun's atomic BoxAcc checkpoints (done-point bitmaps)
//! DIR/<f31|f32b>/seed<S>/commits.log
//! ```
//!
//! One subdirectory per node seed: a "redraw the nodes" retry (seed + 1) starts its own.
//! Certification is that of `box_word`: the generator-only determinant pass (at the points of the
//! run before any scan and again at every needed point in `box_finish`), frozen-structure replay,
//! Freivalds D² (`close_stage_box`), χ = decat_chi and the D14 scaling check, with the same
//! `BoxOpts` (`BoxOpts::new` production defaults); the caller compares the two primes as before.
//! Point values do not depend on how the points are batched or checkpointed, and the closures
//! read only the accumulators, so the results equal those of the uninterrupted run.
//!
//! Environment: `KRHOM_CKPT_SECS` (minimum seconds between checkpoint saves, default 60; 0 =
//! after every batch of 16 points) and the test hook `KRHOM_CKPT_STOP_AFTER=M` (stop after M
//! points of this invocation without saving the points scanned since the last checkpoint, i.e.
//! a simulated kill; the process exits with code 75). Parallelism within a run is that of the
//! in-memory route (`KRHOM_THREADS` points at a time); several machines: `sbox box --shard`.
use crate::calc::{BoxParams, CalcErr, RegEntry};
use crate::colored::boxclose::PrimeField;
use crate::colored::boxrun::{box_finish_with, box_run_with, build_commit, records, BoxJob, RunOpts, BOXRUN_STAMP};
use crate::colored::swire::{BoxOpts, BoxWordOut, WLib};
use std::path::{Path, PathBuf};
use std::sync::Mutex;

pub const CKPT_FORMAT: &str = "krhom-checkpoint v1";

static DIR: Mutex<Option<PathBuf>> = Mutex::new(None);

/// `--checkpoint DIR` of the command line (krhom `parse_args`): one braid per run (a `--file`
/// batch would put several computations in one directory) and not with `--check-library`.
pub fn configure(dir: Option<&str>, items: usize, check_library: bool) -> Result<(), CalcErr> {
    let Some(d) = dir else {
        *DIR.lock().unwrap() = None;
        return Ok(());
    };
    if d.trim().is_empty() {
        return Err(CalcErr::Unsupported("--checkpoint needs a directory".into()));
    }
    if items != 1 {
        return Err(CalcErr::Unsupported("--checkpoint takes one braid (--braid or --pd), not a --file batch: give each braid its own directory".into()));
    }
    if check_library {
        return Err(CalcErr::Unsupported("--checkpoint with --check-library: the library check runs no box scan; drop --checkpoint".into()));
    }
    *DIR.lock().unwrap() = Some(PathBuf::from(d));
    Ok(())
}

/// The configured checkpoint root, if any.
pub fn dir() -> Option<PathBuf> {
    DIR.lock().unwrap().clone()
}

fn join<T: ToString>(v: &[T], sep: &str) -> String {
    v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(sep)
}

/// The manifest of one prime: everything the computation depends on that the caller chose
/// (boxrun's run.hdr adds the library hash, engine stamp and record traces).
pub fn manifest(reg: &RegEntry, field: &str, p: &BoxParams) -> String {
    format!(
        "{}\nfield={}\nbraid={}\nstrands={}\nlabels={}\nN={}\nmode={}\nseed={}\nbox_n={}\nlibrary={}\nlibrary_field={}\ndata_commit={}\nprovider={}\nkrhom_commit={}\nboxrun={}\n",
        CKPT_FORMAT,
        field,
        join(&p.word, " "),
        reg.lambda.len(),
        join(&reg.lambda, ","),
        join(&p.ns, ","),
        if p.reduced { "reduced (strand at bottom position 0 marked)" } else { "unreduced" },
        p.seed,
        reg.n,
        reg.id,
        reg.field,
        reg.data_commit,
        reg.provider,
        build_commit(),
        BOXRUN_STAMP
    )
}

/// Lines that differ between a stored and a wanted manifest (`key: directory has X, this run has Y`).
pub fn manifest_diff(have: &str, want: &str) -> Vec<String> {
    let kv = |s: &str| -> Vec<(String, String)> { s.lines().map(|l| l.split_once('=').map_or((l.to_string(), String::new()), |(k, v)| (k.to_string(), v.to_string()))).collect() };
    let (h, w) = (kv(have), kv(want));
    let mut d = vec![];
    for (k, v) in &w {
        match h.iter().find(|(hk, _)| hk == k) {
            Some((_, hv)) if hv == v => {}
            Some((_, hv)) => d.push(format!("{}: directory has {:?}, this run has {:?}", k, hv, v)),
            None => d.push(format!("{}: missing in the directory, this run has {:?}", k, v)),
        }
    }
    for (k, v) in &h {
        if !w.iter().any(|(wk, _)| wk == k) {
            d.push(format!("{}: directory has {:?}, unknown to this krhom", k, v));
        }
    }
    d
}

/// The prime's checkpoint directory: create it with its manifest, or check the stored one
/// (mismatch: exit 1, nothing touched). Returns the directory and whether it existed.
pub fn prepare(root: &Path, reg: &RegEntry, field: &str, p: &BoxParams) -> Result<(PathBuf, bool), CalcErr> {
    let un = |m: String| CalcErr::Unsupported(m);
    let d = root.join(&reg.field);
    let mp = d.join("manifest.txt");
    let want = manifest(reg, field, p);
    if mp.exists() {
        let have = std::fs::read_to_string(&mp).map_err(|e| un(format!("--checkpoint {}: {}", mp.display(), e)))?;
        let diff = manifest_diff(&have, &want);
        if !diff.is_empty() {
            return Err(un(format!(
                "--checkpoint {}: the directory holds a different computation (resume refused; use another directory, or delete this one to start over):\n  {}",
                d.display(),
                diff.join("\n  ")
            )));
        }
        return Ok((d, true));
    }
    if d.exists() && std::fs::read_dir(&d).map(|mut r| r.next().is_some()).unwrap_or(true) {
        return Err(un(format!("--checkpoint {}: the directory is not empty and has no manifest.txt (not a krhom checkpoint); use another directory", d.display())));
    }
    std::fs::create_dir_all(&d).map_err(|e| un(format!("--checkpoint {}: {}", d.display(), e)))?;
    let tmp = d.join(format!("manifest.txt.tmp{}", std::process::id()));
    std::fs::write(&tmp, &want).and_then(|_| std::fs::rename(&tmp, &mp)).map_err(|e| un(format!("--checkpoint {}: {}", mp.display(), e)))?;
    Ok((d, false))
}

fn env_num<T: std::str::FromStr>(k: &str) -> Option<T> {
    std::env::var(k).ok().and_then(|v| v.trim().parse().ok())
}

/// `box_word` through the checkpointed driver: scan the pending points of `dir/seed<S>` (resuming
/// from its checkpoints), then merge, check and close exactly as `box_finish` does. Returns the
/// `BoxWordOut` fields the calculator reads (stages, scaling, χ count, points, determinant pass,
/// timings).
pub fn box_word_ckpt<F: PrimeField>(lib: &WLib<F>, word: &[i32], o: &BoxOpts, dir: &Path) -> Result<BoxWordOut, String> {
    let job = BoxJob { lam: lib.lam.clone(), words: vec![word.to_vec()], stages: vec![o.stages.clone()], o: o.clone() };
    let recs = records(lib, &job)?;
    let sd = dir.join(format!("seed{}", o.seed));
    let ro = RunOpts { shard: (0, 1), ckpt_secs: env_num("KRHOM_CKPT_SECS").unwrap_or(60.0), stop_after: env_num("KRHOM_CKPT_STOP_AFTER") };
    let st = box_run_with(lib, &job, &recs, &sd, &ro)?;
    eprintln!(
        "[{}] checkpoint {}: {} box points, {} already done{}, {} scanned now, {} checkpoint files written",
        F::name(),
        sd.display(),
        st.shard_points,
        st.shard_points - st.pending_at_start,
        if st.resumed { " (resumed)" } else { "" },
        st.scanned,
        st.checkpoints_written
    );
    if !st.complete {
        if let Some(m) = ro.stop_after {
            eprintln!("[{}] KRHOM_CKPT_STOP_AFTER={}: simulated kill (test hook); rerun the same command to resume", F::name(), m);
            std::process::exit(75);
        }
        return Err(format!("checkpoint {}: points still pending after the run", sd.display()));
    }
    let fo = box_finish_with(lib, &job, &recs, &[sd])?;
    let rec = &recs.recs[0];
    let wr = fo.words.into_iter().next().ok_or("box_finish: no word")?;
    let final_types = rec.stages.last().map(|s| s.skel.objs.iter().map(|o| lib.names[&o.0].clone()).collect()).unwrap_or_default();
    Ok(BoxWordOut {
        stages: wr.stages,
        record: rec.ps.stages.clone(),
        record_secs: recs.secs,
        lb: rec.lb.clone(),
        points: fo.det_pass.0,
        costs: st.costs,
        scan_secs: st.secs,
        close_secs: fo.close_secs,
        final_types,
        scaling: wr.scaling,
        chi_checked: wr.chi_checked,
        det_pass: fo.det_pass,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params(word: Vec<i32>, reduced: bool, seed: u64) -> BoxParams {
        BoxParams { word, ns: vec![2, 3, 4], reduced, seed, retries: 2, mem_mb: 3200, check_only: false, replay_only: true }
    }

    #[test]
    fn manifest_mismatch_is_refused_with_exit_1() {
        let reg = crate::calc::reg_for(&[2, 2, 2], "f31").unwrap();
        let root = std::env::temp_dir().join(format!("krhom_ckpt_manifest_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        let p = params(vec![1, -2, 1, -2], true, 201);
        let (d, existed) = prepare(&root, &reg, "F_2147483647", &p).unwrap();
        assert!(!existed && d.join("manifest.txt").is_file());
        // the same input resumes
        assert!(prepare(&root, &reg, "F_2147483647", &p).unwrap().1);
        // any difference is refused (exit 1) and names the field
        for (q, key) in [(params(vec![1, -2, 1, 2], true, 201), "braid"), (params(vec![1, -2, 1, -2], false, 201), "mode"), (params(vec![1, -2, 1, -2], true, 202), "seed"), (BoxParams { ns: vec![2, 3], ..p.clone() }, "N")] {
            let e = prepare(&root, &reg, "F_2147483647", &q).unwrap_err();
            assert_eq!(e.code(), 1, "{}", e.msg());
            assert!(e.msg().contains("resume refused") && e.msg().contains(&format!("\n  {}: directory has", key)), "{}", e.msg());
        }
        // another library / data commit
        let mut r2 = reg.clone();
        r2.data_commit = "1".repeat(40);
        assert!(prepare(&root, &r2, "F_2147483647", &p).unwrap_err().msg().contains("data_commit: directory has"));
        // a tampered krhom commit line
        let mp = d.join("manifest.txt");
        let t = std::fs::read_to_string(&mp).unwrap();
        std::fs::write(&mp, t.replace(&format!("krhom_commit={}", build_commit()), "krhom_commit=0000000")).unwrap();
        assert!(prepare(&root, &reg, "F_2147483647", &p).unwrap_err().msg().contains("krhom_commit: directory has \"0000000\""));
        // a non-empty foreign directory is refused
        let reg32 = crate::calc::reg_for(&[2, 2, 2], "f32b").unwrap();
        std::fs::create_dir_all(root.join("f32b")).unwrap();
        std::fs::write(root.join("f32b/x"), "x").unwrap();
        assert!(prepare(&root, &reg32, "F_4294967291", &p).unwrap_err().msg().contains("not a krhom checkpoint"));
        let _ = std::fs::remove_dir_all(&root);
    }

    #[test]
    fn configure_refusals() {
        assert_eq!(configure(Some("d"), 2, false).unwrap_err().code(), 1);
        assert_eq!(configure(Some("d"), 1, true).unwrap_err().code(), 1);
        assert_eq!(configure(Some(" "), 1, false).unwrap_err().code(), 1);
        configure(None, 1, false).unwrap();
        assert!(dir().is_none());
    }
}
