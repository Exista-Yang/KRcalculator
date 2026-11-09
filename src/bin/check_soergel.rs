// Compare Soergel's Hom formula with hom_space on small Bott–Samelson bimodules.
use krhom::field::F31;
use krhom::hecke::Hecke;
use krhom::mpoly::MonoCache;
use krhom::sbim::{hom_space, SBim};
use krhom::web::Rung;

fn bs(n: usize, word: &[usize]) -> SBim<F31> {
    let mut b = SBim::<F31>::identity(n, &vec![1; n]);
    for &i in word {
        b = b.add_rung(&Rung::e(i, 1)).add_rung(&Rung::f(i, 1));
    }
    b
}

fn dim_r(n: usize, d: i32) -> i128 {
    // dim of degree-d part of Q[x1..xn], deg x = 2
    if d < 0 || d % 2 != 0 {
        return 0;
    }
    let k = (d / 2) as i128;
    // C(k+n-1, n-1)
    let mut r: i128 = 1;
    for i in 0..(n as i128 - 1) {
        r = r * (k + 1 + i) / (i + 1);
    }
    r
}

fn main() {
    let n = 3;
    let h = Hecke::new(n);
    let words: Vec<Vec<usize>> = vec![vec![], vec![0], vec![1], vec![0, 1], vec![1, 0], vec![0, 0], vec![0, 1, 0], vec![1, 0, 1]];
    let mut bad = 0;
    for a in &words {
        for b in &words {
            let ba = bs(n, a);
            let bb = bs(n, b);
            let grk = h.hom_grk(&h.bs_char(a), &h.bs_char(b));
            let mut mc = MonoCache::new(ba.ring.clone());
            for d in -4..=6 {
                // normalized degree d <-> unnormalized u = d + shift(a) - shift(b)
                let u = d + ba.shift - bb.shift;
                let got = hom_space(&ba, &bb, u, &mut mc).len() as i128;
                let exp: i128 = grk.c.iter().map(|(&k, &c)| c * dim_r(n, d - k)).sum();
                if got != exp {
                    bad += 1;
                    println!("MISMATCH {:?} -> {:?} d={} got {} exp {} (grk {})", a, b, d, got, exp, grk);
                }
            }
        }
    }
    println!("mismatches: {}", bad);
}
