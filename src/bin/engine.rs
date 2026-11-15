// Localized Soergel engine for uncolored braid closures.
// usage: engine N strands "word" [cube]   (word: signed 1-based generators)
use krhom::input::crossing::{Gen, QConv};
use krhom::reference::cube::{khr_full_cube, BraidLink};
use krhom::uncolored::engine::khr_engine;
use krhom::algebra::field::F31;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let nn: usize = args[1].parse().unwrap();
    let n: usize = args[2].parse().unwrap();
    let word: Vec<Gen> = args[3].split_whitespace().map(|t| { let k: i32 = t.parse().unwrap(); Gen { i: (k.abs() - 1) as usize, positive: k > 0 } }).collect();
    let t0 = std::time::Instant::now();
    let r = khr_engine::<F31>(n, nn, &word, 1);
    println!("engine: {:.2}s  max objects {}  final objects {}", t0.elapsed().as_secs_f64(), r.max_objects, r.final_objects);
    println!("  {}", r.poincare);
    if args.get(4).map(|s| s == "cube").unwrap_or(false) {
        let t1 = std::time::Instant::now();
        let c = khr_full_cube::<F31>(&BraidLink { n: nn, labels: vec![1; n], word: word.clone() }, QConv::A);
        println!("cube:   {:.2}s", t1.elapsed().as_secs_f64());
        println!("  {}", c.poincare);
        println!("MATCH: {}", c.poincare == r.poincare);
    }
}
