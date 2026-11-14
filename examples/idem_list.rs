//! C2-2: decomposable canonical ladders (M5i idempotents) of λ: idem_list 2,2,2,2 [N]
use krhom::schur::Schur;
use krhom::sindec::{decomposable_ladders, idem_needed};
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let lam: Vec<u8> = a[1].split(',').map(|x| x.parse().unwrap()).collect();
    let n: Option<u8> = a.get(2).map(|x| x.parse().unwrap());
    let mut s = Schur::new(&lam);
    let all = decomposable_ladders(&mut s);
    for c in &all {
        println!("decomposable: block {:?} type {} ladder {} flows {} kill_level {} char {:?}", c.mu, c.table, c.ladder, c.flows, c.kill_level, c.char_expansion);
    }
    if let Some(n) = n {
        let v = idem_needed(&mut s, Some(n));
        println!("idem_needed({:?}; N = {}) = {} types: {:?}", lam, n, v.len(), v.iter().map(|c| c.table.clone()).collect::<Vec<_>>());
    }
    println!("decomposable total {}", all.len());
}
