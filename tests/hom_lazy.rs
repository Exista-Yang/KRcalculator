//! The lazy-verified kernel solver of `sbim::hom_space` (linsys::kernel_lazy) returns exactly the
//! reduced kernel basis of the reference solver (incremental Echelon over every equation).
use krhom::field::F31;
use krhom::mpoly::MonoCache;
use krhom::sbim::{hom_space, hom_space_reference, SBim};
use krhom::sexpr::Expr;
use krhom::web::{Dir, Rung};

fn w(s: &str) -> Vec<Rung> {
    if s.is_empty() {
        return vec![];
    }
    s.split(',')
        .map(|t| {
            let dir = if t.starts_with('E') { Dir::E } else { Dir::F };
            let (i, r) = t[1..].split_once('.').unwrap();
            Rung { i: i.parse().unwrap(), dir, r: r.parse().unwrap() }
        })
        .collect()
}

fn check(labels: &[u8], a: &str, b: &str, degs: std::ops::RangeInclusive<i32>) -> usize {
    let (ea, eb) = (Expr::new(labels, &w(a)).unwrap(), Expr::new(labels, &w(b)).unwrap());
    let (sa, sb): (SBim<F31>, SBim<F31>) = (ea.sbim(), eb.sbim());
    let mut n = 0;
    for d in degs {
        let u = d - sb.shift + sa.shift;
        let mut mc = MonoCache::new(sa.ring.clone());
        let fast = hom_space(&sa, &sb, u, &mut mc);
        let mut mc = MonoCache::new(sa.ring.clone());
        let slow = hom_space_reference(&sa, &sb, u, &mut mc);
        assert_eq!(fast.len(), slow.len(), "{:?} [{}] -> [{}] degree {}: dimensions differ", labels, a, b, d);
        for (x, y) in fast.iter().zip(&slow) {
            assert_eq!(format!("{:?}", x), format!("{:?}", y), "{:?} [{}] -> [{}] degree {}: bases differ", labels, a, b, d);
        }
        n += fast.len();
    }
    n
}

#[test]
fn lazy_kernel_equals_reference() {
    let mut total = 0;
    total += check(&[1, 3], "E0.2", "E0.1,E0.1", -3..=5);
    total += check(&[1, 3], "E0.1,E0.1", "E0.2", -3..=5);
    total += check(&[2, 2], "E0.1,F0.1", "F0.1,E0.1", -2..=4);
    total += check(&[2, 2], "E0.1,F0.1,E0.1,F0.1", "E0.1,F0.1", -2..=4);
    total += check(&[2, 2], "", "E0.1,F0.1", -2..=6);
    total += check(&[2, 2], "E0.2,F0.2", "E0.1,F0.1,E0.1,F0.1", 0..=4);
    total += check(&[2, 1, 2], "E0.1,E1.1", "E1.1,E0.1", -2..=4);
    total += check(&[1, 5, 0], "E0.1,F1.2", "F1.2,E0.1", -2..=2);
    total += check(&[2, 2, 2], "E0.1,F1.1", "F1.1,E0.1", 0..=4);
    total += check(&[3, 3], "E0.1,F0.1", "F0.1,E0.1", 0..=4);
    assert!(total > 20, "too few maps compared: {}", total);
}
