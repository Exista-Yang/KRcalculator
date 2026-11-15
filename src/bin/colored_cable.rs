//! colored_cable: exterior-colored gl_N homology of a braid closure by cabling.
//!
//! usage:
//!   colored_cable N[,N,..] (--braid WORD [--strands n] | --pd PD) [--colors c1,c2,..]
//!                 [--closure tower|exact] [--end right|left|alt] [--seed s]
//!                 [--raw] [--check-euler | --no-check] [--cube]
//!
//! Every a-colored strand becomes a parallel strands with one thick-edge insertion
//! B_{w0(S_a)}; the cable runs through the uncolored Soergel engine (one tangle-level scan
//! for all N, then the gl_N closure), each homological degree is divided by Π[a]! (one
//! factor per colored strand position) and the per-crossing shifts t^{±T} q^{∓(NT+m(m-1))}
//! (m = min(a,b), T = ab - m) are undone, so that the output is in krhom's conventions
//! (equal to the full cube `krhom --colored-method cube`). Computed over F_p, p = 2^31 - 1.
//! Colors are given per link component (default all 1); N < max color gives 0.
//!
//! --raw          also print the cable homology before division and shift
//! --check-euler  compare the Euler characteristic with the MOY state sum over all
//!                resolutions (independent oracle; done by default when there are at most
//!                2000 resolutions, forced up to 10^6 with this flag); --no-check skips it
//! --cube         also run the full cube (exponential; tiny cases only) and compare
use krhom::reference::cabling::{cable_default, colored_euler_moy, khr_colored_cable};
use krhom::input::crossing::{Gen, QConv};
use krhom::reference::cube::{khr_full_cube, BraidLink};
use krhom::uncolored::engine2::Closure;
use krhom::algebra::field::F31;
use krhom::input::pd::{parse_pd, pd_to_braid, BraidWord};
use krhom::uncolored::tower::EndPolicy;

fn usage() -> ! {
    eprintln!("usage: colored_cable N[,N,..] (--braid WORD [--strands n] | --pd PD) [--colors c1,c2,..] [--closure tower|exact] [--end right|left|alt] [--seed s] [--raw] [--check-euler | --no-check] [--cube]");
    std::process::exit(1)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 || args[1].starts_with("--") {
        usage();
    }
    let ns: Vec<usize> = args[1].split(',').map(|x| x.trim().parse().unwrap_or_else(|_| usage())).collect();
    if ns.iter().any(|&n| n < 1) {
        usage();
    }
    let get = |k: &str| args.iter().position(|a| a == k).and_then(|i| args.get(i + 1)).cloned();
    let flag = |k: &str| args.iter().any(|a| a == k);
    let bw: BraidWord = if let Some(pd) = get("--pd") {
        let pd = parse_pd(&pd).unwrap_or_else(|e| {
            eprintln!("invalid PD code: {}", e);
            std::process::exit(1)
        });
        pd_to_braid(&pd)
    } else if let Some(w) = get("--braid") {
        let w: Vec<i64> = w.split_whitespace().map(|t| t.parse().unwrap_or_else(|_| usage())).collect();
        let strands: usize = get("--strands").map(|s| s.parse().unwrap_or_else(|_| usage())).unwrap_or_else(|| w.iter().map(|x| x.unsigned_abs() as usize).max().unwrap_or(0) + 1);
        BraidWord::from_signed(strands, &w)
    } else {
        usage()
    };
    let ncomp = bw.num_components();
    let colors: Vec<usize> = match get("--colors") {
        Some(c) => c.split(',').map(|x| x.trim().parse().unwrap_or_else(|_| usage())).collect(),
        None => vec![1; ncomp],
    };
    if colors.len() != ncomp || colors.iter().any(|&c| c == 0) {
        eprintln!("need one color >= 1 per component ({} components)", ncomp);
        std::process::exit(1);
    }
    let labels: Vec<usize> = bw.components.iter().map(|&c| colors[c]).collect();
    let word: Vec<Gen> = bw.word.clone();
    let end = match get("--end") {
        Some(e) => EndPolicy::parse(&e).unwrap_or_else(|| usage()),
        None => EndPolicy::Right,
    };
    let closure = match get("--closure").as_deref() {
        None | Some("tower") => Closure::Tower(end),
        Some("exact") => Closure::Exact,
        Some(c) => Closure::parse(c).unwrap_or_else(|| usage()),
    };
    let seed: u64 = get("--seed").map(|s| s.parse().unwrap_or_else(|_| usage())).unwrap_or(7);

    println!("gl_N for N in {:?}  strands {}  crossings {}  components {}  colors {:?}", ns, bw.strands, word.len(), ncomp, colors);
    println!("braid: {}", bw.signed().iter().map(|x| x.to_string()).collect::<Vec<_>>().join(" "));
    println!("strand labels (bottom positions): {:?}", labels);
    if !cable_default(&labels) {
        eprintln!("warning: a cable of this size ({} strands, colors {:?}) used 0.4-6 GB within seconds in tests", labels.iter().sum::<usize>(), colors);
    }
    let t0 = std::time::Instant::now();
    let r = khr_colored_cable::<F31>(&labels, &word, &ns, seed, closure).unwrap_or_else(|e| {
        eprintln!("cabling failed: {}", e);
        std::process::exit(2)
    });
    let cw: Vec<String> = r.cabled.word.iter().map(|g| if g.positive { (g.i + 1).to_string() } else { format!("-{}", g.i + 1) }).collect();
    println!("cable: {} strands, {} crossings, start object B_w with w = {:?} (0-based simple reflections)", r.cabled.strands, r.cabled.word.len(), r.cabled.start);
    if flag("--raw") {
        println!("cable word: {}", cw.join(" "));
    }
    println!(
        "closure: {}; time {:.2}s (scan {:.2}s, data {:.2}s, closures {:.2}s; {} final objects)",
        match closure {
            Closure::Exact => "exact".to_string(),
            Closure::Tower(e) => format!("tower ({:?} end first)", e),
        },
        t0.elapsed().as_secs_f64(),
        r.timings.0,
        r.timings.1,
        r.timings.2,
        r.final_objects
    );
    let mut failed = false;
    for &nn in &ns {
        let p = &r.colored[&nn];
        println!("N = {}", nn);
        if flag("--raw") {
            match r.raw.get(&nn) {
                Some(raw) => println!("  cable homology (raw): {}", raw),
                None => println!("  cable homology (raw): not computed (N < max color)"),
            }
        }
        println!("  Poincare polynomial: {}", p);
        println!("  total rank: {}", p.total_rank());
        println!("  Euler characteristic: {}", p.euler());
        if !flag("--no-check") {
            let max_states = if flag("--check-euler") { 1_000_000 } else { 2000 };
            match colored_euler_moy(&labels, &word, nn, max_states) {
                Ok(Some(chi)) if chi == p.euler() => println!("  Euler characteristic vs MOY state sum: ok"),
                Ok(Some(chi)) => {
                    println!("  Euler characteristic vs MOY state sum: MISMATCH (MOY {})", chi);
                    failed = true;
                }
                Ok(None) => println!("  Euler characteristic vs MOY state sum: skipped (more than {} resolutions)", max_states),
                Err(e) => {
                    println!("  Euler characteristic vs MOY state sum: error {}", e);
                    failed = true;
                }
            }
        }
        if flag("--cube") {
            let tc = std::time::Instant::now();
            let c = khr_full_cube::<F31>(&BraidLink { n: nn, labels: labels.clone(), word: word.clone() }, QConv::A).poincare;
            let ok = &c == p;
            println!("  full cube ({:.2}s): {}", tc.elapsed().as_secs_f64(), if ok { "equal".to_string() } else { format!("DIFFERENT: {}", c) });
            failed |= !ok;
        }
    }
    if failed {
        std::process::exit(2);
    }
}
