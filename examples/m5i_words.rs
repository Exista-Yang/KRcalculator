//! C2-4 debugging: canonical ladders of the given tables and their flows; with --x TABLE RUNG the
//! word L_TABLE·RUNG, its character, and for each summand q whether L_q is commutation-equivalent.
use krhom::schur::{flows, Schur};
use krhom::sindec::canonical_ladder;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let lam: Vec<u8> = a[1].split(',').map(|x| x.parse().unwrap()).collect();
    let mut s = Schur::new(&lam);
    let mut find = |s: &mut Schur, t: &str| {
        let rows: Vec<Vec<u8>> = t.split('/').map(|r| r.split(',').map(|x| x.parse().unwrap()).collect()).collect();
        let mu: Vec<u8> = (0..rows[0].len()).map(|j| rows.iter().map(|r| r[j]).sum()).collect();
        for mu in [mu] {
            for q in s.canonical_basis(&mu).order.clone() {
                if s.table(q).fmt() == t {
                    return q;
                }
            }
        }
        panic!("no type {}", t)
    };
    let p = find(&mut s, &a[2]);
    let mut w = canonical_ladder(&mut s, p);
    println!("L_p = {} ({} flows)", krhom::smapio::rungs_str(&w), flows(&lam, &w));
    if let Some(r) = a.get(3) {
        w.extend(krhom::smapio::parse_rungs(r).unwrap());
        let cls = krhom::homprov::commutation_class(&w, 4096);
        let ch = s.ladder_char(&w);
        for (q, m) in s.expand(&ch) {
            let lq = canonical_ladder(&mut s, q);
            let ce = cls.contains(&lq);
            println!("  {}{{{}}} L_q {} ({} flows) comm-equivalent to X: {}", s.table(q).fmt(), m.fmt_ec(), krhom::smapio::rungs_str(&lq), flows(&lam, &lq), ce);
        }
        println!("  X = {} ({} flows), class size {}", krhom::smapio::rungs_str(&w), flows(&lam, &w), cls.len());
    }
}
