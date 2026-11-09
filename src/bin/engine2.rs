// M1 engine: one scan, exact closure for several N; optional comparison with the old engine.
// usage: engine2 strands "word" Nlist [compare]   e.g. engine2 4 "1 1 2 -3" 2,3,4 compare
use krhom::crossing::Gen;
use krhom::engine::khr_engine;
use krhom::engine2::khr_engine2;
use krhom::field::F31;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let n: usize = args[1].parse().unwrap();
    let word: Vec<Gen> = args[2].split_whitespace().map(|t| { let k: i32 = t.parse().unwrap(); Gen { i: (k.abs() - 1) as usize, positive: k > 0 } }).collect();
    let ns: Vec<usize> = args[3].split(',').map(|x| x.parse().unwrap()).collect();
    let t0 = std::time::Instant::now();
    let r = khr_engine2::<F31>(n, &word, &ns, &[], 7);
    println!("engine2: total {:.2}s (scan {:.2}s, interp {:.2}s, close {:.2}s) points {} max objects {} final {}", t0.elapsed().as_secs_f64(), r.timings.0, r.timings.1, r.timings.2, r.points, r.max_objects, r.final_objects);
    for (nn, p) in &r.poincare {
        println!("N={}: rank {}  {}", nn, p.total_rank(), if p.to_string().len() < 400 { p.to_string() } else { "(long)".into() });
        if args.get(4).map(|s| s == "compare").unwrap_or(false) {
            let old = khr_engine::<F31>(n, *nn, &word, 1);
            println!("   old engine agrees: {}", old.poincare == *p);
        }
    }
}
