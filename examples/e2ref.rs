//! engine2 uncolored references: e2ref "WORD" N [reverse] — prints unreduced and reduced gl_N.
use krhom::crossing::Gen;
use krhom::field::F31;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let mut w: Vec<i32> = a[1].split_whitespace().map(|x| x.parse().unwrap()).collect();
    let n: usize = a[2].parse().unwrap();
    if a.get(3).map(|s| s == "reverse").unwrap_or(false) {
        w.reverse();
    }
    let strands = w.iter().map(|g| g.unsigned_abs() as usize).max().unwrap() + 1;
    let gens: Vec<Gen> = w.iter().map(|&g| Gen { i: (g.unsigned_abs() - 1) as usize, positive: g > 0 }).collect();
    for red in [false, true] {
        let r = krhom::engine2::khr_engine2_mode::<F31>(strands, &gens, &[n], &[], 7, red);
        let p = &r.poincare[&n];
        println!("word {:?} gl_{} {}: {}  (rank {})", w, n, if red { "reduced" } else { "unreduced" }, p, p.total_rank());
    }
}
