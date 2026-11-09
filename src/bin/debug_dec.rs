use krhom::field::{Field, F31};
use krhom::hecke::Hecke;
use krhom::indec::*;
use krhom::soergel::Ctx;
use std::collections::HashMap;
fn main() {
    let n = 3;
    let x: Vec<F31> = vec![F31::from_i64(17), F31::from_i64(1001), F31::from_i64(77777)];
    let mut ctx = Ctx::<F31>::new(n, x);
    let h = Hecke::new(n);
    let mut indecs: HashMap<usize, Indec<F31>> = HashMap::new();
    for xe in 0..ctx.g.size() { let i = indecomposable(&mut ctx, xe); indecs.insert(xe, i); }
    for (xe, s) in [(0usize, 0usize), (1, 0), (1, 1)] {
        let ix = &indecs[&xe];
        let mut word = ix.word.clone(); word.push(s);
        let e = ix.e.tensor_letter(s);
        let mult: Vec<(usize, i32, usize)> = kl_mult(&h, xe, s).into_iter().map(|(y, k, m)| (y, -k, m)).collect();
        println!("x={} word={:?} s={} mult={:?} e={:?}", xe, ix.word, s, mult, e.cols);
        let sums = decompose(&mut ctx, &word, &e, &mult, &indecs);
        for sm in &sums {
            println!("  y={} deg={} iota src {:?} tgt {:?} cols {:?}", sm.y, sm.deg, sm.iota.src, sm.iota.tgt, sm.iota.cols);
            println!("     pi cols {:?}", sm.pi.cols);
            println!("     pi∘iota {:?}", sm.pi.compose(&sm.iota).cols);
        }
    }
}
