//! C2-4 (M5i scoping): for each kept decomposable canonical ladder, the prefixes' characters
//! (indecomposable? canonical?) and the Exact unknowns of the extra-summand pairs.
//! usage: m5i_prefix LAM N
use krhom::homprov::hom_unknowns;
use krhom::sindec::{canonical_ladder, idem_needed, is_indecomposable_of};
use krhom::schur::{flows, Schur};
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let lam: Vec<u8> = a[1].split(',').map(|x| x.parse().unwrap()).collect();
    let n: u8 = a[2].parse().unwrap();
    let mut s = Schur::new(&lam);
    for c in idem_needed(&mut s, Some(n)) {
        let w = krhom::smapio::parse_rungs(&c.ladder).unwrap();
        println!("{} {} flows {}", c.table, c.ladder, c.flows);
        for k in 1..=w.len() {
            let ch = s.ladder_char(&w[..k]);
            let ex = s.expand(&ch);
            let parts: Vec<String> = ex.iter().map(|(t, m)| format!("{}{{{}}}{}", s.table(*t).fmt(), m.fmt_ec(), if canonical_ladder(&mut s, *t) == w[..k] { "*" } else { "" })).collect();
            println!("  k={} flows {} {}", k, flows(&lam, &w[..k]), parts.join(" + "));
        }
        let chw = s.ladder_char(&w);
        for (t, m) in s.expand(&chw) {
            if t == c.q { continue; }
            let lr = canonical_ladder(&mut s, t);
            let ind = is_indecomposable_of(&mut s, &lr, t);
            println!("  extra {} {{{}}} L_r {} flows {} indec {} unknowns(L_r->L_p, 0) {} (L_p->L_r) {}", s.table(t).fmt(), m.fmt_ec(), krhom::smapio::rungs_str(&lr), flows(&lam, &lr), ind, hom_unknowns(&lam, &lr, &w, 0), hom_unknowns(&lam, &w, &lr, 0));
        }
    }
}
