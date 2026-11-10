//! krhom: colored Khovanov–Rozansky gl_N homology of links.
//!
//! usage:
//!   krhom N --pd "[[1,4,2,5],[3,6,4,1],[5,2,6,3]]" [--colors 1,2,...] [--cube] [--field q]
//!   krhom N --braid "1 1 2 -3" --strands 4 [--colors ...] [--cube] [--field q]
//!   uncolored only: [--closure tower|exact] [--end right|left|alt] [--reduced]
//! N may be a comma-separated list (e.g. 2,3,4,10): uncolored links are then computed
//! from a single tangle-level scan.
//! Colors are given per link component (default all 1). Uncolored links use the
//! localized Soergel engine; colored links currently use the full cube.
//! --closure: gl_N closure of the Soergel engine: the restriction tower (default) or the
//!   exact pairing-matrix closure (the oracle). --end: strand end the tower closes first.
//! --reduced: reduced homology, marking the strand at bottom position 0 of the braid
//!   (normalized so that the reduced unknot is q^0).
use krhom::crossing::{Gen, QConv};
use krhom::cube::{khr_full_cube, BraidLink};
use krhom::engine2::Closure;
use krhom::field::{F31, Q};
use krhom::pd::{parse_pd, pd_to_braid, BraidWord};
use krhom::tower::EndPolicy;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: krhom N (--pd PD | --braid WORD --strands n) [--colors c1,c2,..] [--cube] [--field q] [--closure tower|exact] [--end right|left|alt] [--reduced]");
        std::process::exit(1);
    }
    let ns: Vec<usize> = args[1].split(',').map(|x| x.parse().expect("N")).collect();
    let nn = ns[0];
    let get = |k: &str| args.iter().position(|a| a == k).and_then(|i| args.get(i + 1)).cloned();
    let flag = |k: &str| args.iter().any(|a| a == k);
    let bw: BraidWord = if let Some(pd) = get("--pd") {
        let pd = parse_pd(&pd).expect("invalid PD code");
        pd_to_braid(&pd)
    } else {
        let w: Vec<i64> = get("--braid").expect("--braid").split_whitespace().map(|t| t.parse().unwrap()).collect();
        let strands: usize = get("--strands").map(|s| s.parse().unwrap()).unwrap_or_else(|| w.iter().map(|x| x.unsigned_abs() as usize).max().unwrap_or(0) + 1);
        BraidWord::from_signed(strands, &w)
    };
    let ncomp = bw.num_components();
    let colors: Vec<usize> = match get("--colors") {
        Some(c) => c.split(',').map(|x| x.parse().unwrap()).collect(),
        None => vec![1; ncomp],
    };
    assert_eq!(colors.len(), ncomp, "need one color per component ({} components)", ncomp);
    let labels: Vec<usize> = bw.components.iter().map(|&c| colors[c]).collect();
    let word: Vec<Gen> = bw.word.clone();
    println!("gl_N for N in {:?}  strands {}  crossings {}  components {}  colors {:?}", ns, bw.strands, word.len(), ncomp, colors);
    println!("braid: {}", bw.signed().iter().map(|x| x.to_string()).collect::<Vec<_>>().join(" "));
    let t0 = std::time::Instant::now();
    let uncolored = labels.iter().all(|&l| l == 1);
    let field_q = get("--field").map(|f| f == "q").unwrap_or(false);
    let reduced = flag("--reduced");
    let end = match get("--end") {
        Some(e) => EndPolicy::parse(&e).unwrap_or_else(|| {
            eprintln!("--end must be right, left or alt");
            std::process::exit(1)
        }),
        None => EndPolicy::Right,
    };
    let closure = match get("--closure").as_deref() {
        None | Some("tower") => Closure::Tower(end),
        Some("exact") => Closure::Exact,
        Some(c) => Closure::parse(c).unwrap_or_else(|| {
            eprintln!("--closure must be tower or exact");
            std::process::exit(1)
        }),
    };
    let soergel = uncolored && !flag("--cube");
    if reduced && !soergel {
        eprintln!("--reduced is only supported for uncolored links with the Soergel engine");
        std::process::exit(1);
    }
    if soergel {
        println!("closure: {}{}", match closure { Closure::Exact => "exact".to_string(), Closure::Tower(e) => format!("tower ({:?} end first)", e) }, if reduced { ", reduced (strand at bottom position 0 marked)" } else { "" });
    }
    let results: Vec<(usize, krhom::laurent::Poincare)> = if soergel {
        if field_q {
            eprintln!("note: the Soergel engine runs over F_p (p = 2^31-1); use --cube for exact Q");
        }
        let r = krhom::engine2::khr_engine2_closure::<F31>(bw.strands, &word, &ns, &[], 7, reduced, closure);
        r.poincare.into_iter().collect()
    } else {
        ns.iter()
            .map(|&nn| {
                let link = BraidLink { n: nn, labels: labels.clone(), word: word.clone() };
                let p = if field_q { khr_full_cube::<Q>(&link, QConv::A).poincare } else { khr_full_cube::<F31>(&link, QConv::A).poincare };
                (nn, p)
            })
            .collect()
    };
    println!("time {:.2}s", t0.elapsed().as_secs_f64());
    for (nn, p) in results {
        println!("N = {}", nn);
        println!("  Poincare polynomial: {}", p);
        println!("  total rank: {}", p.total_rank());
        println!("  Euler characteristic: {}", p.euler());
    }
}
