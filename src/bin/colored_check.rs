// Compare cabled-engine colored homology with the full cube (colored Rickard complexes).
use krhom::colored::cable;
use krhom::crossing::{Gen, QConv};
use krhom::cube::{khr_full_cube, BraidLink};
use krhom::engine::khr_engine_start;
use krhom::field::F31;
use krhom::laurent::{LPoly, Poincare};

fn divide(p: &Poincare, d: &LPoly) -> Option<Poincare> {
    let mut out = Poincare::default();
    let ts: std::collections::BTreeSet<i32> = p.c.keys().map(|k| k.0).collect();
    for t in ts {
        let mut poly = LPoly::zero();
        for (&(tt, q), &c) in &p.c { if tt == t { poly.add_term(q, c); } }
        // exact division
        let q = std::panic::catch_unwind(|| poly.div_exact(d)).ok()?;
        for (&e, &c) in &q.c { out.add_term(t, e, c); }
    }
    Some(out)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let nn: usize = args[1].parse().unwrap();
    let labels: Vec<usize> = args[2].split(',').map(|x| x.parse().unwrap()).collect();
    let word: Vec<Gen> = args[3].split_whitespace().map(|t| { let k: i32 = t.parse().unwrap(); Gen { i: (k.abs() - 1) as usize, positive: k > 0 } }).collect();
    let cb = cable(&labels, &word);
    let t0 = std::time::Instant::now();
    let r = khr_engine_start::<F31>(cb.strands, nn, &cb.word, 1, &cb.start);
    let mut fact = LPoly::one();
    for &a in &labels { fact = fact.mul(&LPoly::qfact(a as i64)); }
    // only one factor per component: divide by product over blocks of the closure's components
    println!("cable: strands {} crossings {} start {:?} time {:.2}s", cb.strands, cb.word.len(), cb.start, t0.elapsed().as_secs_f64());
    println!("  raw: {}", r.poincare);
    if args.get(4).map(|s| s == "cube").unwrap_or(false) {
        let c = khr_full_cube::<F31>(&BraidLink { n: nn, labels: labels.clone(), word: word.clone() }, QConv::A);
        println!("  cube: {}", c.poincare);
    }
    let _ = divide(&r.poincare, &fact).map(|d| println!("  raw / Π[a]!: {}", d));
}
