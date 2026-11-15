//! C2-4 (M5i scoping): for each kept idempotent type p and rung R, X = L_p·R, and per kept target
//! (q, s, cnt) of term_mult(p, R) the Exact unknowns of Hom(L_q, X) at s.
//! usage: m5i_srcjobs LAM N [TYPE]
use krhom::homprov::hom_unknowns;
use krhom::schur::{flows, Schur};
use krhom::sindec::{canonical_ladder, idem_needed};
use krhom::web::{Dir, Rung};
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let lam: Vec<u8> = a[1].split(',').map(|x| x.parse().unwrap()).collect();
    let n: u8 = a[2].parse().unwrap();
    let mut s = Schur::new(&lam);
    for c in idem_needed(&mut s, Some(n)) {
        if a.get(3).is_some_and(|t| *t != c.table) {
            continue;
        }
        let w = krhom::smapio::parse_rungs(&c.ladder).unwrap();
        let mu = c.mu.clone();
        for i in 0..lam.len() - 1 {
            for dir in [Dir::E, Dir::F] {
                let giver = if dir == Dir::E { mu[i + 1] } else { mu[i] } as usize;
                for r in 1..=giver {
                    let rg = Rung { i, dir, r };
                    let mut x = w.clone();
                    x.push(rg);
                    let mut parts = vec![];
                    let mut umax = 0u128;
                    for (q, m) in s.term_mult(c.q, &[rg]) {
                        if s.killed_at(q, n as usize) {
                            continue;
                        }
                        let lq = canonical_ladder(&mut s, q);
                        for (k, &cc) in m.c.iter().enumerate() {
                            if cc != 0 {
                                let d = m.lo + k as i32;
                                let u = hom_unknowns(&lam, &lq, &x, d);
                                umax = umax.max(u);
                                parts.push(format!("{}{{v^{}:{}}} {}", s.table(q).fmt(), d, cc, u));
                            }
                        }
                    }
                    println!("{} ⊗ {}: X {} flows, max unknowns {} | {}", c.table, krhom::smapio::rung_str(&rg), flows(&lam, &x), umax, parts.join("; "));
                }
            }
        }
    }
}
