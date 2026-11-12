//! xscan: run the exact singular colored engine (src/xengine.rs) on a 2-strand braid and print
//! statistics and Poincaré polynomials.
//!
//! usage: xscan LABELS WORD NS [--reduced] [--both] [--prime 31|32b] [--no-d2] [--prune N]
//!              [--full] [--plain] [--close-stages] [--json OUT]
//!   e.g. xscan 2,2 "1 1 1" 2,3,4
//!   --prune N: Foam_N pruning (`XOpts::prune_n`); then every N in NS must be <= N, and
//!     --no-d2 is refused (the pruned closure is certified through the d∘d components)
//!   --full: full build of C ⊗ X before elimination (default: streaming, `XOpts::stream`)
//!   --plain: PMat component store (default: compact `CMat`)
//!   --close-stages: also close every stage m >= 2 (the braid word[..m]) at every N
//!   --both: unreduced and reduced closures (default: unreduced, or reduced with --reduced)
//!   --json OUT: write per-stage measurements and all Poincaré polynomials to OUT
//!   --chi-only: no scan; print schur::decat_chi of every prefix word[..m] (m >= 2) at every N
use krhom::crossing::Gen;
use krhom::field::{Field, F31, F32b};
use krhom::laurent::Poincare;
use krhom::xengine::{rss_kb, xscan, XOpts};
use std::collections::BTreeMap;
use std::fmt::Write as _;

fn poly_json(p: &Poincare) -> String {
    let terms: Vec<String> = p.c.iter().map(|(&(t, q), &k)| format!("[{},{},{}]", t, q, k)).collect();
    format!("{{\"str\":\"{}\",\"rank\":{},\"terms\":[{}]}}", p, p.total_rank(), terms.join(","))
}

struct Args {
    ns: Vec<usize>,
    modes: Vec<bool>,
    close_stages: bool,
    json: Option<String>,
}

/// Returns false if a base-point check failed (the process then exits nonzero).
fn run<F: Field>(labels: &[usize], word: &[Gen], a: &Args, opts: XOpts) -> bool {
    if let Some(n0) = opts.prune_n {
        if let Some(&n) = a.ns.iter().find(|&&n| n > n0) {
            eprintln!("--prune {}: N = {} > {} is not computed by a pruned scan", n0, n, n0);
            std::process::exit(1)
        }
    }
    let t0 = std::time::Instant::now();
    let sc = match xscan::<F>(labels, word, &opts) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("xscan failed: {}", e);
            std::process::exit(2)
        }
    };
    println!("field {}", F::name());
    println!("scan {:.2}s, closure data {:.2}s", sc.timings.0, sc.timings.1);
    println!("hom_space calls {} total {:.2}s max {:.2}s; decompositions by Hom spaces {} (identity {}); max decomposed fibre {}", sc.hom_calls, sc.hom_total_secs, sc.hom_max_secs, sc.decompositions, sc.identity_decompositions, sc.max_x_fsize);
    for st in &sc.stages {
        println!("  crossing {} labels {:?}: {} -> {} objects (LB {}), {} components, {:.2}s", st.step, st.labels, st.unsimplified, st.objects, st.lower_bound, st.components, st.secs);
        if opts.prune_n.is_some() {
            println!("    pruned: {} summands dropped, {} d∘d components in the ideal of killed types", st.dropped, st.d2_ideal);
        }
    }
    if let Some(n0) = opts.prune_n {
        println!("pruned at N = {}: {} summands dropped in decompositions, {} decompositions entirely killed, {} final d∘d components (closure checked)", n0, sc.dropped_summands, sc.all_killed_decompositions, sc.final_d2.len());
    }
    println!("types: {:?}", sc.ladders);
    println!("final objects (table, h, t): {:?}", sc.final_objects);
    // (stage m, reduced) -> N -> P
    let mut res: BTreeMap<(usize, bool), BTreeMap<usize, Poincare>> = BTreeMap::new();
    let mut basepoint: BTreeMap<usize, bool> = BTreeMap::new();
    let tc = std::time::Instant::now();
    let mut stages: Vec<usize> = if a.close_stages { sc.closed_stages() } else { vec![] };
    stages.push(word.len());
    for &m in &stages {
        for &reduced in &a.modes {
            let r = sc.close_stage(m, &a.ns, reduced).unwrap_or_else(|e| panic!("{}", e));
            for (nn, p) in &r {
                println!("stage {} N = {}{}", m, nn, if reduced { " (reduced, block 0 marked)" } else { "" });
                println!("  Poincare polynomial: {}", p);
                println!("  total rank: {}", p.total_rank());
                println!("  Euler characteristic: {}", p.euler());
            }
            res.insert((m, reduced), r);
        }
        // base-point check for knots (odd 2-strand stages close to one component)
        if a.modes.contains(&true) && labels.len() == 2 && m % 2 == 1 {
            let r1 = sc.close_stage_marked(m, &a.ns, Some(1)).unwrap_or_else(|e| panic!("{}", e));
            let ok = r1 == res[&(m, true)];
            println!("stage {} base-point check (marked block 1 = marked block 0): {}", m, if ok { "PASS" } else { "FAIL" });
            basepoint.insert(m, ok);
        }
    }
    let t_close = tc.elapsed().as_secs_f64();
    println!("closures {:.2}s", t_close);
    // the process peak: xengine resets VmHWM at every crossing, so VmHWM alone is the peak since
    // the last crossing; take the maximum with every stage's peak (2s_close review F4)
    let peak_kb = sc.stages.iter().map(|s| s.mem.rss_peak_kb).max().unwrap_or(0).max(rss_kb().1);
    println!("total {:.2}s, peak RSS {} kB", t0.elapsed().as_secs_f64(), peak_kb);
    let bp_ok = basepoint.values().all(|&b| b);
    if let Some(path) = &a.json {
        let mut j = String::new();
        let w: Vec<String> = word.iter().map(|g| format!("{}", if g.positive { g.i as i64 + 1 } else { -(g.i as i64 + 1) })).collect();
        write!(j, "{{\"field\":\"{}\",\"labels\":{:?},\"word\":[{}],\"prune_n\":{},\"stream\":{},\"compact\":{},", F::name(), labels, w.join(","), opts.prune_n.map_or("null".into(), |n| n.to_string()), opts.stream, opts.compact).unwrap();
        // certification switches and provenance (2s_close F2, 2s_core m2): check_d2 also gates the
        // Π∘Ι = id composite check; a record with check_d2 = false is not certified
        write!(j, "\"check_d2\":{},\"certified\":{},\"rung_transport\":{},\"rungwise\":{},\"seed\":{},\"threads\":{},\"commit\":\"{}\",\"basepoint_ok\":{},", opts.check_d2, opts.check_d2 && bp_ok, opts.rung_transport, opts.rungwise, opts.seed, krhom::num_threads(), krhom::boxrun::build_commit(), bp_ok).unwrap();
        write!(j, "\"scan_secs\":{:.3},\"closure_data_secs\":{:.3},\"closure_secs\":{:.3},\"peak_rss_kb\":{},\"peak_rss_scope\":\"process\",\"stages\":[", sc.timings.0, sc.timings.1, t_close, peak_kb).unwrap();
        let st: Vec<String> = sc
            .stages
            .iter()
            .map(|s| {
                let m = &s.mem;
                format!(
                    "{{\"step\":{},\"unsimplified\":{},\"objects\":{},\"components\":{},\"lb\":{},\"secs\":{:.3},\"comps_before\":{},\"entries_before\":{},\"terms_before\":{},\"bytes_before\":{},\"entries_after\":{},\"terms_after\":{},\"bytes_after\":{},\"bytes_caches\":{},\"rss_tensor_kb\":{},\"rss_end_kb\":{},\"rss_peak_kb\":{},\"secs_decomp\":{:.3},\"secs_transport\":{:.3},\"secs_elim\":{:.3},\"secs_d2\":{:.3},\"stream_passes\":{},\"peak_objects\":{}}}",
                    s.step, s.unsimplified, s.objects, s.components, s.lower_bound, s.secs, m.comps_before, m.entries_before, m.terms_before, m.bytes_before, m.entries_after, m.terms_after, m.bytes_after, m.bytes_caches, m.rss_tensor_kb, m.rss_end_kb, m.rss_peak_kb, m.secs_decomp, m.secs_transport, m.secs_elim, m.secs_d2, m.stream_passes, m.peak_objects
                )
            })
            .collect();
        write!(j, "{}],\"results\":[", st.join(",")).unwrap();
        let rs: Vec<String> = res.iter().flat_map(|((m, red), r)| r.iter().map(move |(nn, p)| format!("{{\"stage\":{},\"reduced\":{},\"N\":{},\"P\":{}}}", m, red, nn, poly_json(p)))).collect();
        let bp: Vec<String> = basepoint.iter().map(|(m, ok)| format!("\"{}\":{}", m, ok)).collect();
        write!(j, "{}],\"basepoint_check\":{{{}}}}}\n", rs.join(","), bp.join(",")).unwrap();
        std::fs::write(path, j).expect("write --json");
    }
    if !bp_ok {
        eprintln!("FATAL: base-point check FAILED");
    }
    bp_ok
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 4 {
        eprintln!("usage: xscan LABELS WORD NS [--reduced] [--both] [--prime 31|32b] [--no-d2] [--prune N] [--full] [--plain] [--close-stages] [--json OUT]");
        std::process::exit(1);
    }
    let labels: Vec<usize> = args[1].split(',').map(|x| x.parse().expect("label")).collect();
    let word: Vec<Gen> = args[2]
        .split_whitespace()
        .map(|t| {
            let g: i32 = t.parse().expect("generator");
            Gen { i: (g.unsigned_abs() - 1) as usize, positive: g > 0 }
        })
        .collect();
    let ns: Vec<usize> = args[3].split(',').map(|x| x.parse().expect("N")).collect();
    let flag = |k: &str| args.iter().any(|a| a == k);
    let prime = args.iter().position(|a| a == "--prime").and_then(|i| args.get(i + 1)).cloned().unwrap_or_else(|| "31".into());
    let prune: Option<usize> = args.iter().position(|a| a == "--prune").map(|i| args.get(i + 1).and_then(|x| x.parse().ok()).expect("--prune N"));
    let json = args.iter().position(|a| a == "--json").map(|i| args.get(i + 1).cloned().expect("--json OUT"));
    let modes = if flag("--both") { vec![false, true] } else { vec![flag("--reduced")] };
    if flag("--chi-only") {
        // decategorified χ (schur::decat_chi, unreduced) of every prefix word[..m], m >= 2, at
        // every N: "chi m N [[e, c], ...]", for checking the closures of a scan
        let lam: Vec<u8> = labels.iter().map(|&a| a as u8).collect();
        let w: Vec<i32> = word.iter().map(|g| if g.positive { g.i as i32 + 1 } else { -(g.i as i32 + 1) }).collect();
        for m in 2..=w.len() {
            let mut sch = krhom::schur::Schur::new(&lam);
            let opts = krhom::schur::ScanOpts { peak: false, every_stage: false, marked_block: None };
            let res = match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| krhom::schur::decat_scan(&mut sch, &w[..m], &opts))) {
                Ok(r) => r,
                Err(_) => continue, // colors do not close up
            };
            for &nn in &ns {
                let c = krhom::schur::decat_chi(&mut sch, &res, &w[..m], nn);
                let t: Vec<String> = c.c.iter().map(|(e, k)| format!("[{},{}]", e, k)).collect();
                println!("chi {} {} [{}]", m, nn, t.join(","));
            }
        }
        return;
    }
    let a = Args { ns, modes, close_stages: flag("--close-stages"), json };
    let opts = XOpts { check_d2: !flag("--no-d2"), verbose: true, prune_n: prune, stream: !flag("--full"), compact: !flag("--plain"), close_stages: flag("--close-stages"), stage_close_ns: if flag("--close-stages") { a.ns.clone() } else { vec![] }, ..Default::default() };
    let ok = match prime.as_str() {
        "31" => run::<F31>(&labels, &word, &a, opts),
        "32b" => run::<F32b>(&labels, &word, &a, opts),
        p => {
            eprintln!("--prime must be 31 or 32b, not {}", p);
            std::process::exit(1)
        }
    };
    if !ok {
        std::process::exit(3); // a base-point FAIL is a failed run (2s_core m2)
    }
}
