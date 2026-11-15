// Full-cube colored gl_N homology of a braid closure.
// usage: cube N labels word   e.g.  cube 2 1,1 "1 1 1"   (generators i>0 positive, -i negative, 1-based)
use krhom::input::crossing::{Gen, QConv};
use krhom::reference::cube::{khr_full_cube, BraidLink};
use krhom::algebra::field::{F31, Q};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let n: usize = args[1].parse().unwrap();
    let labels: Vec<usize> = args[2].split(',').map(|x| x.parse().unwrap()).collect();
    let word: Vec<Gen> = args[3]
        .split_whitespace()
        .map(|t| {
            let k: i32 = t.parse().unwrap();
            Gen { i: (k.abs() - 1) as usize, positive: k > 0 }
        })
        .collect();
    let link = BraidLink { n, labels, word };
    let field = args.get(4).map(|s| s.as_str()).unwrap_or("p");
    let t0 = std::time::Instant::now();
    let r = if field == "q" { khr_full_cube::<Q>(&link, QConv::A) } else { khr_full_cube::<F31>(&link, QConv::A) };
    println!("vertices={} chain_dim={} time={:.2}s", r.vertices, r.chain_dim, t0.elapsed().as_secs_f64());
    println!("Poincare: {}", r.poincare);
    println!("Euler:    {}", r.poincare.euler());
    println!("rank:     {}", r.poincare.total_rank());
}
