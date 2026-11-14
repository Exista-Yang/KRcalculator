//! sbox: the box loop through the localized path (src/swire.rs, src/boxrun.rs): sindec library →
//! sengine scans at the box points → BoxAcc → close_box. Gate G2B of M6a-1b, the first 3-strand
//! runs, and (M6a-2a) the checkpointed, resumable, shardable driver.
//!
//! usage: sbox [box|finish|point] --lam 3,3 --word "1 1 1" | --words "1 1 1;1 1 1 1 1" --n 6
//!             [--ns 3,4,5,6] [--field 31|32b] [--seed S] [--stages 2,3] [--prune N] [--cache DIR]
//!             [--costs FILE] [--provider default|m4c|m4c-lf] [--replay] [--modes u,r] [--d2-all] [--no-d2] [--mem-mb M]
//!             [--batch B] [--json OUT]
//!   (no subcommand): the one-word in-memory loop (box_word)
//!   box    --dir D [--shard k/K] [--ckpt-secs S] [--stop-after M]: scan the shard's points,
//!          checkpointing into D; rerun the same command to resume (after a kill or a restart)
//!   finish --dirs D1,D2,...: merge the shards, check completeness, close, χ, scaling, JSON
//!   point  --id I: scan one box point (value digests per stage)
//!   --stages: closable stages to close (default: every stage whose labels are λ; one word only)
//!   --replay: the library must come from layer A (--cache) with a provider that refuses calls
//!   --modes: u = unreduced extraction (serves N' unreduced), r = reduced (block 0 marked)
//! The D14 scaling check (scans at z and λz at 2 box points) and χ = decat_chi are mandatory: a
//! failure of either exits nonzero and writes no results. (`--scaling` is accepted and ignored.)
use krhom::boxclose::*;
use krhom::boxrun::*;
use krhom::field::{F31, F32b};
use krhom::calc::poly_json;
use krhom::sengine::EngOpts;
use krhom::sindec::*;
use krhom::swire::*;
use std::path::PathBuf;
use std::time::Instant;

fn opt(a: &[String], k: &str) -> Option<String> {
    a.iter().position(|x| x == k).and_then(|i| a.get(i + 1).cloned())
}
fn flag(a: &[String], k: &str) -> bool {
    a.iter().any(|x| x == k)
}
fn list<T: std::str::FromStr>(s: &str) -> Vec<T> {
    s.split([',', ' ']).filter(|x| !x.is_empty()).map(|x| x.parse().ok().expect("bad list")).collect()
}
fn hwm_kb() -> usize {
    krhom::xengine::rss_kb().1
}
fn die(msg: &str) -> ! {
    eprintln!("FATAL: {}", msg);
    std::process::exit(2)
}

/// library and job from the command line
fn setup<F: PrimeField>(a: &[String]) -> (WLib<F>, BoxJob) {
    let lam: Vec<u8> = list(&opt(a, "--lam").expect("--lam"));
    let words: Vec<Vec<i32>> = match opt(a, "--words") {
        Some(w) => w.split(';').map(list).collect(),
        None => vec![list(&opt(a, "--word").expect("--word or --words"))],
    };
    let n: usize = opt(a, "--n").expect("--n").parse().unwrap();
    let ns: Vec<usize> = opt(a, "--ns").map(|s| list(&s)).unwrap_or(vec![n]);
    let seed: u64 = opt(a, "--seed").map(|s| s.parse().unwrap()).unwrap_or(1);
    let prune: Option<u8> = opt(a, "--prune").map(|s| s.parse().unwrap());
    let cache = opt(a, "--cache").map(PathBuf::from);
    let costs = opt(a, "--costs").map(PathBuf::from).unwrap_or(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/homprov/costs_cloud.txt"));
    let modes: Vec<Option<usize>> = opt(a, "--modes").unwrap_or("u,r".into()).split(',').map(|m| if m == "u" { None } else { Some(0) }).collect();
    // C2-5: `--provider m4c | m4c-lf` builds the M4c chain (production chain + Stage 0 + Leaves
    // closure + Leaves pairs) and, above all, its layer-A stamp, so that sbox can load maps written
    // by `libbuild --provider m4c|m4c-lf` (with --replay the provider itself is never called)
    let (prov, stamp) = match opt(a, "--provider").as_deref() {
        None | Some("default") => provider_for_stamped::<F>(&lam, Some(&costs), 1.0),
        Some(pv @ ("m4c" | "m4c-lf")) => krhom::leaves::provider_for_stamped_m4c::<F>(&lam, Some(&costs), 1.0, pv == "m4c-lf"),
        Some(pv) => die(&format!("--provider {}: expected default, m4c or m4c-lf", pv)),
    };
    let prov: Box<dyn krhom::homprov::HomProvider<F>> = if flag(a, "--replay") { Box::new(NoProvider) } else { prov };
    let sd = SIndec::with_stamp(&lam, prov, stamp, SIndecOpts { cache_dir: cache, prune_n: prune, layer_a_route_stamp: flag(a, "--stamp-route"), ..Default::default() });
    let lib = build_wlib(sd, &words, true).unwrap_or_else(|e| die(&format!("library: {}", e)));
    eprintln!("library: {} types, {} closable, build {:.1}s ({} provider calls, {} layer-A loads), hash {:x}, commit {}", lib.shapes.len(), lib.closure.len(), lib.build_secs, lib.stats.provider_calls, lib.stats.layer_a_loaded, lib.hash, build_commit());
    let stages: Vec<Vec<usize>> = match opt(a, "--stages") {
        Some(s) => {
            assert!(words.len() == 1, "--stages needs a single word");
            vec![list(&s)]
        }
        None => words.iter().map(|w| closable_stages(&lam, w)).collect(),
    };
    let eng = EngOpts { check_d2: prune.is_none() && !flag(a, "--no-d2"), order: if flag(a, "--lazy-markowitz") { krhom::sengine::PivotOrder::Lazy } else { krhom::sengine::PivotOrder::Heap }, mem_point_bytes: opt(a, "--mem-mb").map_or(2048, |s| s.parse::<u64>().unwrap()) << 20, ..Default::default() };
    let mut o = BoxOpts::new(n, seed, ns, modes, stages[0].clone(), eng);
    o.d2_all_points = flag(a, "--d2-all");
    if let Some(b) = opt(a, "--batch") {
        o.batch = b.parse().unwrap();
    }
    (lib, BoxJob { lam, words, stages, o })
}

fn sub<F: PrimeField>(cmd: &str, a: &[String]) {
    let t0 = Instant::now();
    let (lib, job) = setup::<F>(a);
    match cmd {
        "box" => {
            let dir = PathBuf::from(opt(a, "--dir").expect("--dir"));
            let shard = opt(a, "--shard").map_or((0, 1), |s| {
                let (k, n) = s.split_once('/').expect("--shard k/K");
                (k.parse().unwrap(), n.parse().unwrap())
            });
            let ro = RunOpts { shard, ckpt_secs: opt(a, "--ckpt-secs").map_or(60.0, |s| s.parse().unwrap()), stop_after: opt(a, "--stop-after").map(|s| s.parse().unwrap()) };
            let st = box_run(&lib, &job, &dir, &ro).unwrap_or_else(|e| die(&format!("box run: {}", e)));
            let np = st.costs.len().max(1);
            let mean = st.costs.iter().map(|c| c.secs).sum::<f64>() / np as f64;
            println!(
                "box run shard {}/{}: {} points, {} pending at start{}, {} scanned (mean {:.3}s/point), complete {}, determinant pass {} checks at {} points, {} checkpoint files written, {:.1}s, VmHWM {} MB",
                shard.0, shard.1, st.shard_points, st.pending_at_start, if st.resumed { " (resumed)" } else { "" }, st.scanned, mean, st.complete, st.det_pass.1, st.det_pass.0, st.checkpoints_written, st.secs, hwm_kb() / 1024
            );
        }
        "finish" => {
            let dirs: Vec<PathBuf> = opt(a, "--dirs").expect("--dirs").split(',').map(PathBuf::from).collect();
            let out = box_finish(&lib, &job, &dirs).unwrap_or_else(|e| die(&format!("box finish: {}", e)));
            for wr in &out.words {
                for sr in &wr.stages {
                    for (&(nn, marked), p) in &sr.closed {
                        println!("word {:?} stage {} N={} {}: {}  (rank {})", wr.word, sr.stage, nn, if marked.is_some() { "reduced" } else { "unreduced" }, p, p.total_rank());
                    }
                }
                let s = wr.scaling.as_ref().unwrap();
                eprintln!("word {:?}: chi = decat_chi for {} closures; scaling: {} entries checked, {} zero at both points, {} declared 0", wr.word, wr.chi_checked, s.report.checked, s.report.zero_both, s.report.declared_zero);
                for (st, m, d, n) in &wr.digests {
                    eprintln!("  accumulator stage {} mode {:?}: {} points, digest {:x}", st, m, n, d);
                }
            }
            eprintln!("finish: shards {:?}, determinant pass {} checks at {} points, close {:.1}s, total {:.1}s", out.shards, out.det_pass.1, out.det_pass.0, out.close_secs, t0.elapsed().as_secs_f64());
            if let Some(path) = opt(a, "--json") {
                std::fs::write(&path, finish_json(&lib, &job, &out)).unwrap();
            }
        }
        "point" => {
            let id: u32 = opt(a, "--id").expect("--id").parse().unwrap();
            let r = box_point(&lib, &job, id).unwrap_or_else(|e| die(&format!("point: {}", e)));
            for (w, st) in job.words.iter().zip(r) {
                for (stage, n, h) in st {
                    println!("point {} word {:?} stage {}: {} entries, values hash {:x}", id, w, stage, n, h);
                }
            }
        }
        _ => unreachable!(),
    }
}

fn run<F: PrimeField>(a: &[String]) {
    let t0 = Instant::now();
    let (lib, job) = setup::<F>(a);
    assert!(job.words.len() == 1, "the one-word loop takes --word (use `box`/`finish` for several words)");
    let (lam, word, n, seed, prune) = (job.lam.clone(), job.words[0].clone(), job.o.n_box, job.o.seed, lib.prune_n);
    let o = job.o.clone();
    if flag(a, "--record-only") {
        let k: usize = lam.iter().map(|&x| x as usize).sum();
        let mut so = ScanOpts { eng: o.eng.clone(), close_stages: o.stages.clone(), k0_check: !flag(a, "--no-k0") };
        so.eng.id_coset = lib.id_coset;
        let t = Instant::now();
        let r = scan_point(&lib, &word, record_x::<F>(k, seed), 0, krhom::sengine::Mode::Record, &so).unwrap();
        for s in &r.ps.stages {
            eprintln!("stage {}: inserted {} objects {} components {} pivots {} rounds {} entries {} c_e {:.2} peak {:.1} MB secs {:.2}", s.stage + 1, s.inserted, s.objects, s.components, s.pivots, s.rounds, s.stored_entries, s.c_e, s.peak_bytes as f64 / 1e6, s.secs);
        }
        eprintln!("build {:.2}s engine {:.2}s", r.t_build, r.t_engine);
        eprintln!("record scan {:.2}s, eval memo {:.1} MB, VmHWM {} MB", t.elapsed().as_secs_f64(), r.eval_bytes as f64 / 1e6, hwm_kb() / 1024);
        return;
    }
    // box_word runs the determinant pass, χ = decat_chi (Err on mismatch) and the D14 scaling
    // check (mandatory here): any failure is FATAL (exit 2, no results)
    let out = box_word(&lib, &word, &o).unwrap_or_else(|e| die(&format!("box: {}", e)));
    let sc = out.scaling.clone().unwrap_or_else(|| die("no scaling check"));
    eprintln!("scaling check: {} entries checked, {} zero at both points, {} declared 0 ({:.1}s)", sc.report.checked, sc.report.zero_both, sc.report.declared_zero, sc.secs_milli as f64 / 1e3);
    eprintln!("chi = decat_chi: {} of {} closures; determinant pass {} checks at {} points", out.chi_checked, out.stages.iter().map(|s| s.closed.len()).sum::<usize>(), out.det_pass.1, out.det_pass.0);
    for sr in &out.stages {
        for (&(nn, marked), p) in &sr.closed {
            println!("stage {} N={} {}: {}  (rank {})", sr.stage, nn, if marked.is_some() { "reduced" } else { "unreduced" }, p, p.total_rank());
        }
    }
    // per-point costs
    let np = out.costs.len().max(1);
    let mean = out.costs.iter().map(|c| c.secs).sum::<f64>() / np as f64;
    let mx = out.costs.iter().map(|c| c.secs).fold(0.0, f64::max);
    let peak = out.costs.iter().map(|c| c.peak_bytes).max().unwrap_or(0);
    let evb = out.costs.iter().map(|c| c.eval_bytes).max().unwrap_or(0);
    let nst = out.costs.first().map_or(0, |c| c.stage_secs.len());
    let st_mean: Vec<f64> = (0..nst).map(|s| out.costs.iter().map(|c| c.stage_secs[s]).sum::<f64>() / np as f64).collect();
    let st_peak: Vec<u64> = (0..nst).map(|s| out.costs.iter().map(|c| c.stage_peak[s]).max().unwrap_or(0)).collect();
    eprintln!("{} points: scan {:.1}s wall, per point mean {:.3}s max {:.3}s, engine peak {:.1} MB, eval memo {:.1} MB; closures {:.1}s; record {:.1}s; VmHWM {} MB", out.points, out.scan_secs, mean, mx, peak as f64 / 1e6, evb as f64 / 1e6, out.close_secs, out.record_secs, hwm_kb() / 1024);
    if let Some(path) = opt(a, "--json") {
        let mut res = vec![];
        for sr in &out.stages {
            for (&(nn, marked), p) in &sr.closed {
                res.push(format!("{{\"stage\":{},\"reduced\":{},\"N\":{},\"P\":{},\"entries\":{},\"objects\":{},\"d2_blocks\":{},\"chi_ok\":true}}", sr.stage, marked.is_some(), nn, poly_json(p), sr.entries, sr.objects, sr.d2_blocks[&(nn, marked)]));
            }
        }
        let rec: Vec<String> = out.record.iter().map(|s| format!("{{\"stage\":{},\"inserted\":{},\"objects\":{},\"components\":{},\"pivots\":{},\"c_e\":{:.3},\"peak_bytes\":{},\"secs\":{:.3}}}", s.stage + 1, s.inserted, s.objects, s.components, s.pivots, s.c_e, s.peak_bytes, s.secs)).collect();
        let lb: Vec<String> = out.lb.iter().map(|(o, l)| format!("[{},{}]", o, l)).collect();
        let ft: Vec<String> = out.final_types.iter().map(|t| format!("\"{}\"", t)).collect();
        let scal = format!("{{\"points\":{:?},\"checked\":{},\"zero_both\":{},\"declared_zero\":{},\"pass\":true}}", sc.points, sc.report.checked, sc.report.zero_both, sc.report.declared_zero);
        let j = format!(
            "{{\"commit\":\"{}\",\"engine\":\"{}\",\"check_d2\":{},\"d2_all_points\":{},\"freivalds_k\":{},\"chi_fatal\":true,\"det_pass\":{{\"points\":{},\"checks\":{}}},\"field\":\"{}\",\"lam\":{:?},\"word\":{:?},\"n_box\":{},\"seed\":{},\"prune_n\":{},\"library_hash\":\"{:x}\",\"lib_types\":{},\"lib_closable\":{},\"lib_build_secs\":{:.2},\"points\":{},\"scan_secs\":{:.2},\"point_secs_mean\":{:.4},\"point_secs_max\":{:.4},\"point_engine_peak_bytes\":{},\"point_eval_bytes\":{},\"stage_secs_mean\":{:?},\"stage_peak_bytes\":{:?},\"record_secs\":{:.2},\"close_secs\":{:.2},\"total_secs\":{:.2},\"vmhwm_kb\":{},\"record\":[{}],\"lb\":[{}],\"final_types\":[{}],\"scaling\":{},\"results\":[\n{}\n]}}\n",
            build_commit(), engine_stamp(&o.eng), o.eng.check_d2, o.d2_all_points, o.freivalds_k, out.det_pass.0, out.det_pass.1, F::name(), lam, word, n, seed, prune.map_or("null".into(), |p| p.to_string()), lib.hash, lib.shapes.len(), lib.closure.len(), lib.build_secs, out.points, out.scan_secs, mean, mx, peak, evb, st_mean, st_peak, out.record_secs, out.close_secs, t0.elapsed().as_secs_f64(), hwm_kb(), rec.join(","), lb.join(","), ft.join(","), scal, res.join(",\n")
        );
        std::fs::write(&path, j).unwrap();
    }
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let cmd = a.get(1).map(|s| s.as_str()).filter(|s| ["box", "finish", "point"].contains(s));
    match (cmd, opt(&a, "--field").as_deref()) {
        (Some(c), Some("32b")) => sub::<F32b>(c, &a),
        (Some(c), _) => sub::<F31>(c, &a),
        (None, Some("32b")) => run::<F32b>(&a),
        (None, _) => run::<F31>(&a),
    }
}
