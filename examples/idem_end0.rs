//! C2-2: size of End^0(L_p) (Exact hom unknowns) for the decomposable kept types (M5i scoping):
//! idem_end0 2,2,2,2 4
use krhom::library::homprov::hom_unknowns;
use krhom::combinat::schur::{flows, Schur};
use krhom::library::sindec::idem_needed;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let lam: Vec<u8> = a[1].split(',').map(|x| x.parse().unwrap()).collect();
    let n: u8 = a[2].parse().unwrap();
    let mut s = Schur::new(&lam);
    for c in idem_needed(&mut s, Some(n)) {
        let w = krhom::library::smapio::parse_rungs(&c.ladder).unwrap();
        println!("{} block {:?} ladder {} flows {} End^0 unknowns {} kill_level {}", c.table, c.mu, c.ladder, flows(&lam, &w), hom_unknowns(&lam, &w, &w, 0), c.kill_level);
    }
}
