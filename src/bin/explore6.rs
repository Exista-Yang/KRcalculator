// Euler characteristic checks of crossing conventions via the MOY oracle.
use krhom::input::crossing::{crossing_cx, qshift, Gen, QConv};
use krhom::algebra::laurent::LPoly;
use krhom::combinat::moy::closure_value;
use krhom::combinat::web::Ladder;

/// Euler characteristic of the closure of a braid word with colors.
fn chi(n: usize, labels: &[usize], word: &[Gen], conv: QConv) -> LPoly {
    // enumerate all states
    let mut states: Vec<(Ladder, i32, i32)> = vec![(Ladder::new(n, labels.to_vec()), 0, 0)];
    for g in word {
        let mut next = vec![];
        for (lad, h, q) in &states {
            let top = lad.top();
            let cx = crossing_cx(*g, top[g.i], top[g.i + 1], n);
            for t in &cx.terms {
                let mut l2 = lad.clone();
                for r in &t.rungs {
                    l2.push(*r);
                }
                next.push((l2, h + t.hdeg, q + qshift(&cx, t.k, conv)));
            }
        }
        states = next;
    }
    let mut tot = LPoly::zero();
    for (lad, h, q) in states {
        assert_eq!(lad.top(), lad.bottom, "closure needs matching labels");
        let v = closure_value(&lad).shift(q);
        tot = if h.rem_euclid(2) == 0 { tot.add(&v) } else { tot.sub(&v) };
    }
    tot
}

fn main() {
    let p = |i| Gen { i, positive: true };
    let m = |i| Gen { i, positive: false };
    for conv in [QConv::A, QConv::EtwRev] {
        println!("== convention {:?}", conv);
        for n in 2..=4usize {
            for a in 1..n {
                let unknot = LPoly::qbinom(n as i64, a as i64);
                let r1p = chi(n, &[a, a], &[p(0)], conv);
                let r1m = chi(n, &[a, a], &[m(0)], conv);
                println!("N={} a={}: closure(s1)/unknot ok={}  closure(s1^-1) ok={}", n, a, r1p == unknot.mul(&LPoly::qbinom(1, 1)) && false || r1p == chi(n, &[a], &[], conv).mul(&LPoly::one()), r1m == chi(n, &[a], &[], conv));
                if r1p != unknot {
                    println!("    s1: {}   unknot: {}", r1p, unknot);
                }
                for b in 1..n {
                    let id = chi(n, &[a, b], &[], conv);
                    let r2 = chi(n, &[a, b], &[p(0), m(0)], conv);
                    let r2b = chi(n, &[a, b], &[m(0), p(0)], conv);
                    if r2 != id || r2b != id {
                        println!("    R2 FAIL a={} b={}", a, b);
                    }
                }
            }
        }
        // R3 uncolored and colored on 3 strands
        for n in 2..=4usize {
            for cols in [[1, 1, 1], [1, 2, 1], [2, 1, 1], [1, 1, 2]] {
                if cols.iter().any(|&c| c >= n) { continue; }
                // closure needs a permutation-compatible word; use w * w^{-1}-type check:
                // s1 s2 s1 (s2 s1 s2)^{-1} closure = identity closure
                let l = chi(n, &cols, &[p(0), p(1), p(0), m(1), m(0), m(1)], conv);
                let id = chi(n, &cols, &[], conv);
                if l != id { println!("  R3 FAIL N={} cols={:?}", n, cols); }
            }
        }
        // trefoil and figure eight, uncolored N=2
        let t = chi(2, &[1, 1], &[p(0), p(0), p(0)], conv);
        println!("trefoil N=2: {}", t);
        let f8 = chi(2, &[1, 1, 1], &[p(0), m(1), p(0), m(1)], conv);
        println!("figure-eight N=2: {}", f8);
        let t3 = chi(3, &[1, 1], &[p(0), p(0), p(0)], conv);
        println!("trefoil N=3: {}", t3);
    }
}
