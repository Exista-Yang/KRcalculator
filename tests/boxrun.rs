//! M6a-2a: the checkpointed, resumable, shardable box driver (src/boxrun.rs) and the hardening
//! of the localized path (SPEC_3C_GL6 rev1 D14; review wave 3 box F1/F2, sengine F1/F2/F4/F6).
//!
//! * a run killed (simulated: unsaved work lost) and resumed = an uninterrupted run, byte for
//!   byte (accumulator digests and the finish JSON);
//! * two shards merged = one run; a missing shard is an incomplete down-set (Err);
//! * mismatched shards are refused: another seed, another engine option, a duplicated shard, a
//!   foreign checkpoint file; a resume into the directory of another job is refused;
//! * BoxAcc binding: equal weight lists with different stage / skeleton / trace / record /
//!   engine bindings never merge, and `load_expect` refuses a stale checkpoint;
//! * the D14 scaling check on the wired path catches a weight declared one too high;
//! * χ = decat_chi is fatal (a perturbed closure is an Err).
use krhom::boxclose::*;
use krhom::boxrun::*;
use krhom::field::{Field, F31};
use krhom::sengine::{EngOpts, PivotOrder};
use krhom::sindec::*;
use krhom::swire::*;
use std::path::PathBuf;
use std::sync::Arc;

type F = F31;

fn costs() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/homprov/costs_cloud.txt")
}
fn wlib(lam: &[u8], words: &[Vec<i32>]) -> WLib<F> {
    let (prov, stamp) = provider_for_stamped::<F>(lam, Some(&costs()), 1.0);
    let sd = SIndec::with_stamp(lam, prov, stamp, SIndecOpts::default());
    build_wlib(sd, words, true).unwrap()
}
fn tmp(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("boxrun_test_{}_{}", std::process::id(), name));
    let _ = std::fs::remove_dir_all(&d);
    d
}
/// (2,2) σ^3 (stages 2, 3) and σ^2 in one job, box N = 5, closures at N' = 3..5, both modes
fn job(seed: u64) -> BoxJob {
    let words = vec![vec![1, 1, 1], vec![1, 1]];
    let stages = words.iter().map(|w| closable_stages(&[2, 2], w)).collect();
    let mut o = BoxOpts::new(5, seed, vec![3, 4, 5], vec![None, Some(0)], vec![], EngOpts { check_d2: true, ..Default::default() });
    o.batch = 4;
    BoxJob { lam: vec![2, 2], words, stages, o }
}

#[test]
fn kill_resume_shards_and_refusals() {
    let j = job(31);
    let lib = wlib(&j.lam, &j.words);
    // (a) uninterrupted single run
    let da = tmp("a");
    let st = box_run(&lib, &j, &da, &RunOpts { ckpt_secs: 1e9, ..Default::default() }).unwrap();
    assert!(st.complete && !st.resumed && st.scanned == st.shard_points && st.det_pass.1 > 0, "{:?}", st);
    let fa = box_finish(&lib, &j, &[da.clone()]).unwrap();
    let ja = finish_json(&lib, &j, &fa);
    // the closures equal the in-memory loop (box_word) word by word
    for (wi, w) in j.words.iter().enumerate() {
        let mut o = j.o.clone();
        o.stages = j.stages[wi].clone();
        let bw = box_word(&lib, w, &o).unwrap();
        for (a, b) in bw.stages.iter().zip(fa.words[wi].stages.iter()) {
            assert_eq!(a.closed, b.closed, "word {:?} stage {}", w, a.stage);
        }
        assert!(fa.words[wi].scaling.as_ref().unwrap().report.checked > 0);
        assert_eq!(fa.words[wi].chi_checked, bw.chi_checked);
    }
    eprintln!("[boxrun] uninterrupted: {} points, {} closures", st.shard_points, result_map(&fa).len());

    // (b) killed twice (once with every batch checkpointed, once with nothing saved since the
    // last checkpoint), then resumed to completion: byte-identical
    let db = tmp("b");
    let s1 = box_run(&lib, &j, &db, &RunOpts { ckpt_secs: 0.0, stop_after: Some(7), ..Default::default() }).unwrap();
    assert!(!s1.complete && s1.scanned == 7);
    let s2 = box_run(&lib, &j, &db, &RunOpts { ckpt_secs: 1e9, stop_after: Some(9), ..Default::default() }).unwrap();
    assert!(s2.resumed && s2.pending_at_start == st.shard_points - 4, "the 3 unsaved points of the first kill are redone ({:?})", s2);
    let s3 = box_run(&lib, &j, &db, &RunOpts::default()).unwrap();
    assert!(s3.resumed && s3.complete && s3.pending_at_start == st.shard_points - 4, "{:?}", s3);
    let s4 = box_run(&lib, &j, &db, &RunOpts::default()).unwrap();
    assert!(s4.complete && s4.scanned == 0, "a finished shard rescans nothing");
    let fb = box_finish(&lib, &j, &[db.clone()]).unwrap();
    assert_eq!(fa.words.iter().map(|w| w.digests.clone()).collect::<Vec<_>>(), fb.words.iter().map(|w| w.digests.clone()).collect::<Vec<_>>(), "accumulator digests");
    assert_eq!(ja, finish_json(&lib, &j, &fb), "finish JSON byte-identical");
    assert!(ja.contains("\"shards\":[[0,1]]"), "finish JSON is JSON: {}", &ja[..ja.len().min(400)]);

    // (c) two shards = one run; one shard alone is incomplete
    let (d0, d1) = (tmp("s0"), tmp("s1"));
    let r0 = box_run(&lib, &j, &d0, &RunOpts { shard: (0, 2), ..Default::default() }).unwrap();
    let r1 = box_run(&lib, &j, &d1, &RunOpts { shard: (1, 2), ..Default::default() }).unwrap();
    assert_eq!(r0.shard_points + r1.shard_points, st.shard_points);
    let fc = box_finish(&lib, &j, &[d0.clone(), d1.clone()]).unwrap();
    assert_eq!(result_map(&fa), result_map(&fc));
    assert_eq!(fa.words.iter().map(|w| w.digests.clone()).collect::<Vec<_>>(), fc.words.iter().map(|w| w.digests.clone()).collect::<Vec<_>>(), "merged shards = one run, byte for byte");
    let e = box_finish(&lib, &j, &[d0.clone()]).unwrap_err();
    assert!(e.contains("incomplete down-set"), "{}", e);

    // (d) refusals
    // a duplicated shard (the same directory twice)
    let e = box_finish(&lib, &j, &[d0.clone(), d0.clone()]).unwrap_err();
    assert!(e.contains("duplicates"), "{}", e);
    // a shard of another node seed
    let j2 = job(32);
    let dx = tmp("seed32");
    box_run(&lib, &j2, &dx, &RunOpts { shard: (1, 2), ..Default::default() }).unwrap();
    let e = box_finish(&lib, &j, &[d0.clone(), dx.clone()]).unwrap_err();
    assert!(e.contains("merge refused") && e.contains("seed"), "{}", e);
    // a shard run with another pivot order (engine stamp)
    let mut j3 = job(31);
    j3.o.eng.order = PivotOrder::Rescan;
    let dy = tmp("rescan");
    box_run(&lib, &j3, &dy, &RunOpts { shard: (1, 2), ..Default::default() }).unwrap();
    let e = box_finish(&lib, &j, &[d0.clone(), dy.clone()]).unwrap_err();
    assert!(e.contains("merge refused") && e.contains("engine_stamp"), "{}", e);
    // a resume into the directory of another job, or as another shard
    let e = box_run(&lib, &j2, &d0, &RunOpts { shard: (0, 2), ..Default::default() }).unwrap_err();
    assert!(e.contains("resume refused"), "{}", e);
    let e = box_run(&lib, &j, &d0, &RunOpts { shard: (1, 2), ..Default::default() }).unwrap_err();
    assert!(e.contains("belongs to shard"), "{}", e);
    // a foreign checkpoint file under a matching run header (seed-32 accumulators copied into a
    // copy of shard 1): the checkpoint header is refused
    let dz = tmp("forged");
    std::fs::create_dir_all(&dz).unwrap();
    for f in std::fs::read_dir(&d1).unwrap() {
        let f = f.unwrap().path();
        let src = if f.extension().is_some_and(|e| e == "acc") { dx.join(f.file_name().unwrap()) } else { f.clone() };
        std::fs::copy(&src, dz.join(f.file_name().unwrap())).unwrap();
    }
    let e = box_finish(&lib, &j, &[d0.clone(), dz.clone()]).unwrap_err();
    assert!(e.contains("header does not match") || e.contains("prime, nodes"), "{}", e);
    // a truncated checkpoint is refused on resume
    let acc = std::fs::read_dir(&d0).unwrap().map(|f| f.unwrap().path()).find(|p| p.extension().is_some_and(|e| e == "acc")).unwrap();
    let raw = std::fs::read(&acc).unwrap();
    std::fs::write(&acc, &raw[..raw.len() - 9]).unwrap();
    let e = box_run(&lib, &j, &d0, &RunOpts { shard: (0, 2), ..Default::default() }).unwrap_err();
    assert!(e.contains("corrupt") || e.contains("truncated"), "{}", e);
    for d in [da, db, d0, d1, dx, dy, dz] {
        let _ = std::fs::remove_dir_all(d);
    }
}

#[test]
fn box_acc_binding_is_part_of_the_header() {
    let nodes = BoxNodes::<F>::new(&[2, 2], 4, 7);
    let grid = Arc::new(BoxGrid::new(nodes, None));
    let w = vec![Some(1u32), Some(2), None];
    let b0 = BoxBinding { stage: 3, skeleton_hash: 11, engine_stamp: "e1".into(), trace_hash: 21, record_hash: 31 };
    let mk = |b: BoxBinding| BoxAcc::<F>::with_binding(grid.clone(), w.clone(), 1, 2, None, b);
    let mut base = mk(b0.clone());
    // every binding field separates shards with equal weight lists
    for b in [
        BoxBinding { stage: 2, ..b0.clone() },
        BoxBinding { skeleton_hash: 12, ..b0.clone() },
        BoxBinding { engine_stamp: "e2".into(), ..b0.clone() },
        BoxBinding { trace_hash: 22, ..b0.clone() },
        BoxBinding { record_hash: 32, ..b0.clone() },
    ] {
        let e = base.merge(&mk(b)).unwrap_err();
        assert!(e.contains("mismatched checkpoint header") && e.contains("bind."), "{}", e);
    }
    base.merge(&mk(b0.clone())).unwrap();
    // load_expect: the full header round-trips and a stale binding is refused
    let d = tmp("bind");
    std::fs::create_dir_all(&d).unwrap();
    let p = d.join("x.acc");
    mk(b0.clone()).save(&p).unwrap();
    let want = mk(b0.clone()).header;
    assert_eq!(BoxAcc::load(&p, grid.clone()).unwrap().header, want);
    BoxAcc::load_expect(&p, grid.clone(), &want).unwrap();
    let stale = mk(BoxBinding { trace_hash: 99, ..b0 }).header;
    assert!(BoxAcc::load_expect(&p, grid.clone(), &stale).err().unwrap().contains("trace_hash"));
    let _ = std::fs::remove_dir_all(d);
}

/// Review wave 3 sengine F6: the wired path (scaled nodes → box_x → replay → check_scaling_all)
/// catches an entry declared one weight too high (and one too low).
#[test]
fn scaling_check_catches_a_misdeclared_weight() {
    let (lam, word) = (vec![2u8, 2], vec![1, 1, 1]);
    let lib = wlib(&lam, &[word.clone()]);
    let o = BoxOpts::new(5, 41, vec![3, 4], vec![None], vec![2, 3], EngOpts::default());
    let rec = record_scan(&lib, &word, &o).unwrap();
    let ug = BoxGrid::new(BoxNodes::<F>::new(&lam, 5, 41), None);
    let ids = scaling_points(&ug);
    let lambda = F::from_u64(SCALING_LAMBDA);
    let ok = scaling_check_frozen(&lib, &word, &o, &lambda, &ids, &rec).unwrap();
    assert!(ok.checked > 0);
    // a pairing entry with a nonzero value at the record point, declared at w + 1 and at w − 1
    let si = rec.stages.len() - 1;
    let e = (0..rec.stages[si].kinds.len()).rev().find(|&e| matches!(rec.stages[si].kinds[e], EntryKind::Pair { .. }) && rec.stages[si].weights[e].is_some_and(|w| w > 0) && !rec.stages[si].vals[e].is_zero()).expect("a nonzero pairing entry");
    for d in [1i64, -1] {
        let mut bad = rec.clone();
        let w = bad.stages[si].weights[e].unwrap();
        bad.stages[si].weights[e] = Some((w as i64 + d) as u32);
        let err = scaling_check_frozen(&lib, &word, &o, &lambda, &ids, &bad).unwrap_err();
        assert!(err.contains("scaling law violated"), "{}", err);
    }
}

/// Review wave 3 sengine F2: χ ≠ decat_chi is an Err (fatal in sbox and box_finish).
#[test]
fn chi_mismatch_is_an_error() {
    let (lam, word) = (vec![2u8, 2], vec![1, 1, 1]);
    let lib = wlib(&lam, &[word.clone()]);
    let mut o = BoxOpts::new(5, 51, vec![3, 4, 5], vec![None, Some(0)], vec![3], EngOpts::default());
    o.scaling = false; // TEST-ONLY opt-out (this test is about χ)
    let out = box_word(&lib, &word, &o).unwrap();
    assert!(out.scaling.is_none());
    let closed = out.stages[0].closed.clone();
    assert_eq!(chi_check(&lam, &word, 3, &closed).unwrap(), closed.len());
    for k in closed.keys() {
        let mut bad = closed.clone();
        let p = bad.get_mut(k).unwrap();
        // one more generator in homological degree 0, q^0: χ changes
        *p.c.entry((0, 0)).or_insert(0) += 1;
        let e = chi_check(&lam, &word, 3, &bad).unwrap_err();
        assert!(e.contains("CHI FAIL"), "{}", e);
    }
}

/// Review 2s_close F2/F4, 2s_core m2: the xscan JSON records the certification switches and the
/// commit of the binary, and its peak RSS is the process peak.
#[test]
fn xscan_json_provenance() {
    let d = tmp("xscan");
    std::fs::create_dir_all(&d).unwrap();
    let j = d.join("x.json");
    let o = std::process::Command::new(env!("CARGO_BIN_EXE_xscan")).args(["2,2", "1 1 1", "3,4", "--both", "--json", j.to_str().unwrap()]).env("KRHOM_THREADS", "1").output().unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let t = std::fs::read_to_string(&j).unwrap();
    for k in ["\"check_d2\":true", "\"certified\":true", "\"rung_transport\":", "\"rungwise\":", "\"seed\":", "\"threads\":1", "\"commit\":\"", "\"basepoint_ok\":true", "\"peak_rss_scope\":\"process\""] {
        assert!(t.contains(k), "missing {} in {}", k, &t[..t.len().min(600)]);
    }
    assert!(!t.contains("\"commit\":\"unknown\""));
    let o = std::process::Command::new(env!("CARGO_BIN_EXE_xscan")).args(["2,2", "1 1 1", "3", "--no-d2", "--json", j.to_str().unwrap()]).env("KRHOM_THREADS", "1").output().unwrap();
    assert!(o.status.success());
    let t = std::fs::read_to_string(&j).unwrap();
    assert!(t.contains("\"check_d2\":false") && t.contains("\"certified\":false"), "{}", &t[..t.len().min(600)]);
    let _ = std::fs::remove_dir_all(d);
}
