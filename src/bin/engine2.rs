// M1 engine: one scan, then the gl_N closure for several N.
// usage: engine2 strands "word" Nlist [--closure tower|exact|both] [--end right|left|alt]
//                [--reduced] [--unit-pivots] [full] [compare]
//   e.g. engine2 4 "1 1 2 -3" 2,3,4 --closure both
// --closure: tower (default; restriction tower), exact (close_exact, the oracle), or both
//            (run both and report whether they agree).
// --end:     end variable closed first by the tower (default right).
// --reduced: reduced homology (strand at bottom position 0 marked); "reduced" also accepted.
// --unit-pivots: tower stress mode (only constant pivots, so blocks stay non-split).
// compare:   also run the old engine (unreduced only).
// --save-data FILE: write the prepared closure data (scan, bases, Gram, F) to FILE.
// --load-data FILE: read the prepared data from FILE instead of scanning (strands and
//            word must match the file).
use krhom::input::crossing::Gen;
use krhom::uncolored::engine::khr_engine;
use krhom::uncolored::engine2::{close_tower_with, close_with, load_closure_data, prepare_closure_data, save_closure_data, Closure};
use krhom::algebra::field::F31;
use krhom::uncolored::tower::{EndPolicy, TowerOpts};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 4 {
        eprintln!("usage: engine2 strands \"word\" Nlist [--closure tower|exact|both] [--end right|left|alt] [--reduced] [--save-data FILE | --load-data FILE] [full] [compare]");
        std::process::exit(1);
    }
    let n: usize = args[1].parse().unwrap();
    let word: Vec<Gen> = args[2].split_whitespace().map(|t| { let k: i32 = t.parse().unwrap(); Gen { i: (k.abs() - 1) as usize, positive: k > 0 } }).collect();
    let ns: Vec<usize> = args[3].split(',').map(|x| x.parse().unwrap()).collect();
    let get = |k: &str| args.iter().position(|a| a == k).and_then(|i| args.get(i + 1)).cloned();
    let flag = |k: &str| args.iter().any(|a| a == k);
    let reduced = flag("reduced") || flag("--reduced");
    let full = flag("full");
    let end = match get("--end") {
        Some(e) => EndPolicy::parse(&e).unwrap_or_else(|| panic!("--end: right|left|alt, got {}", e)),
        None => EndPolicy::Right,
    };
    let mode = get("--closure").unwrap_or_else(|| "tower".to_string());
    let closures: Vec<Closure> = match mode.as_str() {
        "tower" => vec![Closure::Tower(end)],
        "exact" => vec![Closure::Exact],
        "both" => vec![Closure::Exact, Closure::Tower(end)],
        m => vec![Closure::parse(m).unwrap_or_else(|| panic!("--closure: tower|exact|both, got {}", m))],
    };
    let t0 = std::time::Instant::now();
    let data = match get("--load-data") {
        Some(path) => {
            let d = load_closure_data::<F31>(&path).unwrap_or_else(|e| panic!("--load-data: {}", e));
            assert!(d.n == n && d.word == word, "--load-data: {} holds a different braid", path);
            d
        }
        None => prepare_closure_data::<F31>(n, &word, &[], 7),
    };
    if let Some(path) = get("--save-data") {
        save_closure_data(&data, &path).unwrap_or_else(|e| panic!("--save-data: {}", e));
    }
    let t_prep = t0.elapsed().as_secs_f64();
    println!(
        "engine2: data {:.2}s (scan {:.2}s, bases+Gram {:.2}s, F {:.2}s) points {} max objects {} final {}{}",
        t_prep, data.timings.0, data.timings.1, data.timings.2, data.points, data.max_objects, data.objs.len(), if reduced { "  [reduced]" } else { "" }
    );
    let mut all_agree = true;
    for &nn in &ns {
        let mut res = vec![];
        for &c in &closures {
            let t1 = std::time::Instant::now();
            let p = match c {
                Closure::Tower(policy) if flag("--unit-pivots") => close_tower_with(&data, nn, reduced, TowerOpts { policy, unit_pivots_only: true, eliminate_top: false }).poincare,
                _ => close_with(&data, nn, reduced, c),
            };
            res.push((c, p, t1.elapsed().as_secs_f64()));
        }
        for (c, p, secs) in &res {
            let s = p.to_string();
            println!("N={} {:?}: {:.2}s rank {}  {}", nn, c, secs, p.total_rank(), if full || s.len() < 400 { s } else { "(long)".into() });
        }
        if res.len() == 2 {
            let ok = res[0].1 == res[1].1;
            all_agree &= ok;
            println!("   N={} exact and tower agree: {}", nn, ok);
        }
        if flag("compare") && !reduced {
            let old = khr_engine::<F31>(n, nn, &word, 1);
            println!("   old engine agrees: {}", old.poincare == res[0].1);
        }
    }
    println!("total {:.2}s", t0.elapsed().as_secs_f64());
    if !all_agree {
        println!("DISAGREEMENT between exact and tower");
        std::process::exit(2);
    }
}
