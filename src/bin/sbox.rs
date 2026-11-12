//! sbox: one-word box loop through the localized path (src/swire.rs): sindec library → sengine
//! scans at the box points → BoxAcc → close_box. Gate G2B of M6a-1b and the first 3-strand runs.
//!
//! usage: sbox --lam 3,3 --word "1 1 1" --n 6 [--ns 3,4,5,6] [--field 31|32b] [--seed S]
//!             [--stages 2,3] [--prune N] [--cache DIR] [--costs FILE] [--replay]
//!             [--modes u,r] [--d2-all] [--no-d2] [--scaling] [--mem-mb M] [--json OUT]
//!   --stages: closable stages to close (default: every stage whose labels are λ)
//!   --replay: the library must come from layer A (--cache) with a provider that refuses calls
//!   --modes: u = unreduced extraction (serves N' unreduced), r = reduced (block 0 marked)
//!   --scaling: D14 scaling check (scans at z and λz) at 2 box points
use krhom::boxclose::*;
use krhom::field::{Field, F31, F32b};
use krhom::laurent::Poincare;
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
fn poly_json(p: &Poincare) -> String {
    let terms: Vec<String> = p.c.iter().map(|(&(t, q), &k)| format!("[{},{},{}]", t, q, k)).collect();
    format!("{{\"str\":\"{}\",\"rank\":{},\"terms\":[{}]}}", p, p.total_rank(), terms.join(","))
}
fn hwm_kb() -> usize {
    krhom::xengine::rss_kb().1
}

fn run<F: PrimeField>(a: &[String]) {
    let t0 = Instant::now();
    let lam: Vec<u8> = list(&opt(a, "--lam").expect("--lam"));
    let word: Vec<i32> = list(&opt(a, "--word").expect("--word"));
    let n: usize = opt(a, "--n").expect("--n").parse().unwrap();
    let ns: Vec<usize> = opt(a, "--ns").map(|s| list(&s)).unwrap_or(vec![n]);
    let seed: u64 = opt(a, "--seed").map(|s| s.parse().unwrap()).unwrap_or(1);
    let prune: Option<u8> = opt(a, "--prune").map(|s| s.parse().unwrap());
    let cache = opt(a, "--cache").map(PathBuf::from);
    let costs = opt(a, "--costs").map(PathBuf::from).unwrap_or(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/homprov/costs_cloud.txt"));
    let modes: Vec<Option<usize>> = opt(a, "--modes").unwrap_or("u,r".into()).split(',').map(|m| if m == "u" { None } else { Some(0) }).collect();
    let (prov, stamp) = provider_for_stamped::<F>(&lam, Some(&costs), 1.0);
    let prov: Box<dyn krhom::homprov::HomProvider<F>> = if flag(a, "--replay") { Box::new(NoProvider) } else { prov };
    let sd = SIndec::with_stamp(&lam, prov, stamp, SIndecOpts { cache_dir: cache, prune_n: prune, ..Default::default() });
    let lib = build_wlib(sd, &[word.clone()], true).unwrap_or_else(|e| panic!("library: {}", e));
    eprintln!("library: {} types, {} closable, build {:.1}s ({} provider calls, {} layer-A loads), hash {:x}", lib.shapes.len(), lib.closure.len(), lib.build_secs, lib.stats.provider_calls, lib.stats.layer_a_loaded, lib.hash);
    let mut labels = lam.clone();
    let mut closable = vec![];
    for (m, &g) in word.iter().enumerate() {
        let (i, _) = letter(g);
        labels.swap(i, i + 1);
        if labels == lam {
            closable.push(m + 1);
        }
    }
    let stages: Vec<usize> = opt(a, "--stages").map(|s| list(&s)).unwrap_or(closable);
    let eng = EngOpts { check_d2: prune.is_none() && !flag(a, "--no-d2"), mem_point_bytes: opt(a, "--mem-mb").map_or(2048, |s| s.parse::<u64>().unwrap()) << 20, ..Default::default() };
    if flag(a, "--record-only") {
        let k: usize = lam.iter().map(|&x| x as usize).sum();
        let mut so = ScanOpts { eng: eng.clone(), close_stages: stages.clone(), k0_check: !flag(a, "--no-k0") };
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
    let o = BoxOpts { n_box: n, seed, n_close: ns.clone(), modes: modes.clone(), stages: stages.clone(), eng, d2_all_points: flag(a, "--d2-all"), freivalds_k: 3 };
    let out = box_word(&lib, &word, &o).unwrap_or_else(|e| panic!("box: {}", e));
    let mut scal = String::from("null");
    if flag(a, "--scaling") {
        let ug = BoxGrid::new(BoxNodes::<F>::new(&lam, n, seed), None);
        let lambda = F::from_u64(0x1234_5677);
        let ts = Instant::now();
        let (c, z) = scaling_check(&lib, &word, &o, &lambda, &[ug.npoints / 3, ug.npoints - 1]).unwrap_or_else(|e| panic!("scaling check FAILED: {}", e));
        eprintln!("scaling check: {} entries checked, {} zero at both points ({:.1}s)", c, z, ts.elapsed().as_secs_f64());
        scal = format!("{{\"checked\":{},\"zero_both\":{},\"pass\":true}}", c, z);
    }
    // χ = schur::decat_chi (unreduced), and [N choose λ_0]·χ(reduced) = χ(unreduced)
    let mut chi: std::collections::BTreeMap<(usize, usize, bool), bool> = Default::default();
    for sr in &out.stages {
        let mut sch = krhom::schur::Schur::new(&lam);
        let sopts = krhom::schur::ScanOpts { peak: false, every_stage: false, marked_block: None };
        let res = krhom::schur::decat_scan(&mut sch, &word[..sr.stage], &sopts);
        for (&(nn, marked), p) in &sr.closed {
            let want = krhom::schur::decat_chi(&mut sch, &res, &word[..sr.stage], nn);
            let got = if marked.is_some() { p.euler().mul(&krhom::laurent::LPoly::qbinom(nn as i64, lam[0] as i64)) } else { p.euler() };
            let ok = got == want;
            if !ok {
                eprintln!("CHI FAIL stage {} N={} marked {:?}: {:?} vs decat_chi {:?}", sr.stage, nn, marked, got, want);
            }
            chi.insert((sr.stage, nn, marked.is_some()), ok);
        }
    }
    eprintln!("chi = decat_chi: {} of {} closures", chi.values().filter(|&&b| b).count(), chi.len());
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
                res.push(format!("{{\"stage\":{},\"reduced\":{},\"N\":{},\"P\":{},\"entries\":{},\"objects\":{},\"d2_blocks\":{},\"chi_ok\":{}}}", sr.stage, marked.is_some(), nn, poly_json(p), sr.entries, sr.objects, sr.d2_blocks[&(nn, marked)], chi[&(sr.stage, nn, marked.is_some())]));
            }
        }
        let rec: Vec<String> = out.record.iter().map(|s| format!("{{\"stage\":{},\"inserted\":{},\"objects\":{},\"components\":{},\"pivots\":{},\"c_e\":{:.3},\"peak_bytes\":{},\"secs\":{:.3}}}", s.stage + 1, s.inserted, s.objects, s.components, s.pivots, s.c_e, s.peak_bytes, s.secs)).collect();
        let lb: Vec<String> = out.lb.iter().map(|(o, l)| format!("[{},{}]", o, l)).collect();
        let ft: Vec<String> = out.final_types.iter().map(|t| format!("\"{}\"", t)).collect();
        let j = format!(
            "{{\"field\":\"{}\",\"lam\":{:?},\"word\":{:?},\"n_box\":{},\"seed\":{},\"prune_n\":{},\"library_hash\":\"{:x}\",\"lib_types\":{},\"lib_closable\":{},\"lib_build_secs\":{:.2},\"points\":{},\"scan_secs\":{:.2},\"point_secs_mean\":{:.4},\"point_secs_max\":{:.4},\"point_engine_peak_bytes\":{},\"point_eval_bytes\":{},\"stage_secs_mean\":{:?},\"stage_peak_bytes\":{:?},\"record_secs\":{:.2},\"close_secs\":{:.2},\"total_secs\":{:.2},\"vmhwm_kb\":{},\"record\":[{}],\"lb\":[{}],\"final_types\":[{}],\"scaling\":{},\"results\":[\n{}\n]}}\n",
            F::name(), lam, word, n, seed, prune.map_or("null".into(), |p| p.to_string()), lib.hash, lib.shapes.len(), lib.closure.len(), lib.build_secs, out.points, out.scan_secs, mean, mx, peak, evb, st_mean, st_peak, out.record_secs, out.close_secs, t0.elapsed().as_secs_f64(), hwm_kb(), rec.join(","), lb.join(","), ft.join(","), scal, res.join(",\n")
        );
        std::fs::write(&path, j).unwrap();
    }
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    match opt(&a, "--field").as_deref() {
        Some("32b") => run::<F32b>(&a),
        _ => run::<F31>(&a),
    }
}
