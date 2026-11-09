// Two-prime agreement check for the M1 engine.
use krhom::crossing::Gen;
use krhom::engine2::khr_engine2;
use krhom::field::{F31, F32b};
fn main() {
    let args: Vec<String> = std::env::args().collect();
    let n: usize = args[1].parse().unwrap();
    let word: Vec<Gen> = args[2].split_whitespace().map(|t| { let k: i32 = t.parse().unwrap(); Gen { i: (k.abs() - 1) as usize, positive: k > 0 } }).collect();
    let ns: Vec<usize> = args[3].split(',').map(|x| x.parse().unwrap()).collect();
    let a = khr_engine2::<F31>(n, &word, &ns, &[], 7);
    let b = khr_engine2::<F32b>(n, &word, &ns, &[], 11);
    for nn in &ns {
        println!("N={}: p1 rank {}  p2 rank {}  equal {}", nn, a.poincare[nn].total_rank(), b.poincare[nn].total_rank(), a.poincare[nn] == b.poincare[nn]);
    }
}
