//! xscan: run the exact singular colored engine (src/xengine.rs) on a 2-strand braid and print
//! statistics and Poincaré polynomials.
//!
//! usage: xscan LABELS WORD NS [--reduced] [--prime 31|32b] [--no-d2] [--prune N]
//!   e.g. xscan 2,2 "1 1 1" 2,3,4
//!   --prune N: Foam_N pruning (`XOpts::prune_n`); then every N in NS must be <= N, and
//!     --no-d2 is refused (the pruned closure is certified through the d∘d components)
use krhom::crossing::Gen;
use krhom::field::{Field, F31, F32b};
use krhom::xengine::{xscan, XOpts};

fn run<F: Field>(labels: &[usize], word: &[Gen], ns: &[usize], reduced: bool, d2: bool, prune_n: Option<usize>, stream: bool) {
    if let Some(n0) = prune_n {
        if let Some(&n) = ns.iter().find(|&&n| n > n0) {
            eprintln!("--prune {}: N = {} > {} is not computed by a pruned scan", n0, n, n0);
            std::process::exit(1)
        }
    }
    let opts = XOpts { check_d2: d2, verbose: true, prune_n, stream, ..Default::default() };
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
        if prune_n.is_some() {
            println!("    pruned: {} summands dropped, {} d∘d components in the ideal of killed types", st.dropped, st.d2_ideal);
        }
    }
    if let Some(n0) = prune_n {
        println!("pruned at N = {}: {} summands dropped in decompositions, {} decompositions entirely killed, {} final d∘d components (closure checked)", n0, sc.dropped_summands, sc.all_killed_decompositions, sc.final_d2.len());
    }
    println!("types: {:?}", sc.ladders);
    println!("final objects (table, h, t): {:?}", sc.final_objects);
    for (nn, p) in sc.close(ns, reduced) {
        println!("N = {}{}", nn, if reduced { " (reduced, block 0 marked)" } else { "" });
        println!("  Poincare polynomial: {}", p);
        println!("  total rank: {}", p.total_rank());
        println!("  Euler characteristic: {}", p.euler());
    }
    println!("total {:.2}s", t0.elapsed().as_secs_f64());
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 4 {
        eprintln!("usage: xscan LABELS WORD NS [--reduced] [--prime 31|32b] [--no-d2] [--prune N]");
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
    match prime.as_str() {
        "31" => run::<F31>(&labels, &word, &ns, flag("--reduced"), !flag("--no-d2"), prune, !flag("--full")),
        "32b" => run::<F32b>(&labels, &word, &ns, flag("--reduced"), !flag("--no-d2"), prune, !flag("--full")),
        p => {
            eprintln!("--prime must be 31 or 32b, not {}", p);
            std::process::exit(1)
        }
    }
}
