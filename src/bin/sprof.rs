//! sprof: per-point profile of the localized scan (M6a-1c): one record point and K replay
//! points of the N-box, each single-threaded, with the time split by stage of the scan.
//!
//! usage: sprof --lam 3,3 --word "1 1 1 1" --n 6 [--field 31|32b] [--seed S] [--cache DIR]
//!              [--costs FILE] [--prune N] [--replay] [--points K] [--no-k0]
use krhom::boxclose::*;
use krhom::field::{F31, F32b};
use krhom::sengine::{EngOpts, Mode};
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

fn show<F: krhom::field::Field>(what: &str, r: &ScanOut<F>, secs: f64) {
    let p = &r.prof;
    let e = &r.eprof;
    let other = secs - p.iter().sum::<f64>() - e.0 - e.1;
    eprintln!(
        "{}: total {:.3}s | dec+obj {:.3} | dC⊗id {:.3} | rickard {:.3} | extract {:.3} | k0 {:.3} | pivots {:.3} ({} updates, {:.3e} mult-adds) | pattern hash {:.3} | rest {:.3} | engine peak {:.1} MB, eval memo {:.1} MB",
        what, secs, p[0], p[1], p[2], p[3], p[4], e.0, e.2, e.3 as f64, e.1, other, r.ps.mem.peak_bytes as f64 / 1e6, r.eval_bytes as f64 / 1e6
    );
    let m = &r.ps.mem;
    eprintln!("   memory peak {:.1} MB at {}: {} objects, {} components, complex {:.1} MB (peak complex alone {:.1} MB), previous complex {:.1} MB", m.peak_bytes as f64 / 1e6, m.peak_at, m.peak_objs, m.peak_edges, (m.peak_bytes - m.peak_prev_bytes) as f64 / 1e6, m.peak_cx_bytes as f64 / 1e6, m.peak_prev_bytes as f64 / 1e6);
    let mut acc = 0.0;
    for s in &r.ps.stages {
        acc += s.secs;
        eprintln!("   stage {}: inserted {} objects {} comps {} pivots {} entries {} peak {:.1} MB cum {:.3}s", s.stage + 1, s.inserted, s.objects, s.components, s.pivots, s.stored_entries, s.peak_bytes as f64 / 1e6, acc);
    }
}

fn run<F: PrimeField>(a: &[String]) {
    let lam: Vec<u8> = list(&opt(a, "--lam").expect("--lam"));
    let word: Vec<i32> = list(&opt(a, "--word").expect("--word"));
    let n: usize = opt(a, "--n").expect("--n").parse().unwrap();
    let seed: u64 = opt(a, "--seed").map(|s| s.parse().unwrap()).unwrap_or(1);
    let prune: Option<u8> = opt(a, "--prune").map(|s| s.parse().unwrap());
    let npts: usize = opt(a, "--points").map(|s| s.parse().unwrap()).unwrap_or(2);
    let cache = opt(a, "--cache").map(PathBuf::from);
    let costs = opt(a, "--costs").map(PathBuf::from).unwrap_or(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/homprov/costs_cloud.txt"));
    let (prov, stamp) = provider_for_stamped::<F>(&lam, Some(&costs), 1.0);
    let prov: Box<dyn krhom::homprov::HomProvider<F>> = if flag(a, "--replay") { Box::new(NoProvider) } else { prov };
    let sd = SIndec::with_stamp(&lam, prov, stamp, SIndecOpts { cache_dir: cache, prune_n: prune, layer_a_route_stamp: flag(a, "--stamp-route"), ..Default::default() });
    let lib = build_wlib(sd, &[word.clone()], true).unwrap_or_else(|e| panic!("library: {}", e));
    eprintln!("library: {} types, build {:.1}s, hash {:x}", lib.shapes.len(), lib.build_secs, lib.hash);
    let mut labels = lam.clone();
    let mut stages = vec![];
    for (m, &g) in word.iter().enumerate() {
        let (i, _) = letter(g);
        labels.swap(i, i + 1);
        if labels == lam {
            stages.push(m + 1);
        }
    }
    let k: usize = lam.iter().map(|&x| x as usize).sum();
    let mut so = ScanOpts { eng: EngOpts { check_d2: false, mem_point_bytes: 12 << 30, ..Default::default() }, close_stages: stages, k0_check: !flag(a, "--no-k0") };
    so.eng.id_coset = lib.id_coset;
    let t = Instant::now();
    let rec = scan_point(&lib, &word, record_x::<F>(k, seed ^ 0x7265_636f_7264), u32::MAX, Mode::Record, &so).unwrap();
    show("record", &rec, t.elapsed().as_secs_f64());
    {
        let mut h = krhom::sengine::Fnv::new();
        for &x in &rec.ps.trace {
            h.word(x);
        }
        for (a, b) in rec.ps.frozen.pivot_sequence() {
            h.word(a as u64);
            h.word(b as u64);
        }
        eprintln!("   record trace digest {:016x}", h.finish());
    }
    let nodes = BoxNodes::<F>::new(&lam, n, seed);
    let ug = BoxGrid::new(nodes.clone(), None);
    so.k0_check = false;
    let frozen = rec.ps.frozen.clone();
    let mut tot = 0.0;
    for j in 0..npts {
        let id = ((j + 1) * ug.npoints as usize / (npts + 1)) as u32;
        let t = Instant::now();
        let r = scan_point(&lib, &word, box_x(&nodes, &ug.mu(id)), id, Mode::Replay(frozen.clone()), &so).unwrap();
        let s = t.elapsed().as_secs_f64();
        tot += s;
        r.ps.same_structure(&rec.ps).unwrap();
        show(&format!("replay point {}", id), &r, s);
        // values digest (result preservation check across code versions)
        let mut h = krhom::sengine::Fnv::new();
        for st in &r.stages {
            for v in &st.vals {
                h.word(v.to_u64());
            }
        }
        eprintln!("   values digest {:016x}", h.finish());
    }
    eprintln!("mean replay {:.3}s; VmHWM {} MB", tot / npts.max(1) as f64, krhom::xengine::rss_kb().1 / 1024);
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    match opt(&a, "--field").as_deref() {
        Some("32b") => run::<F32b>(&a),
        _ => run::<F31>(&a),
    }
}
