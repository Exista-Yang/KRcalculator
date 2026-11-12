//! The lazy-verified kernel solver of `sbim::hom_space` (linsys::kernel_lazy) returns exactly the
//! reduced kernel basis of the reference solver (incremental Echelon over every equation).
use krhom::field::{Field, F31};
use krhom::linsys::{kernel_lazy_counters, kernel_lazy_traced, Echelon};
use rand::{rngs::StdRng, Rng, SeedableRng};
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
    // both branches of the early stop are exercised by these hom_space calls (review solver F2)
    let (calls, early, reins) = kernel_lazy_counters();
    eprintln!("[hom_lazy] kernel_lazy: {} calls, {} early stops, {} with violated-row reinsertion", calls, early, reins);
    assert!(early >= 1, "no hom_space call took the early stop");
    assert!(reins >= 1, "no hom_space call reinserted violated rows after a check");
}

type V = Vec<(u32, F31)>;

fn rand_row(rng: &mut StdRng, ncols: usize, nnz: usize) -> V {
    let mut v: V = vec![];
    for _ in 0..nnz {
        let j = rng.gen_range(0..ncols) as u32;
        let x = F31::new(rng.gen_range(1..2_147_483_647u64));
        v.push((j, x));
    }
    v.sort_by_key(|e| e.0);
    v.dedup_by_key(|e| e.0);
    v
}
/// a random combination of rows (redundant)
fn combo(rng: &mut StdRng, rows: &[V], ncols: usize) -> V {
    let mut acc = vec![F31::zero(); ncols];
    for _ in 0..rng.gen_range(1..4) {
        let r = &rows[rng.gen_range(0..rows.len())];
        let c = F31::new(rng.gen_range(1..1000u64));
        for (j, x) in r {
            acc[*j as usize].fma(&c, x);
        }
    }
    acc.into_iter().enumerate().filter(|(_, x)| !x.is_zero()).map(|(j, x)| (j as u32, x)).collect()
}
fn reference(ncols: usize, eqs: &[V]) -> Vec<V> {
    let mut e: Echelon<F31> = Echelon::new(ncols);
    for r in eqs {
        e.insert(r);
    }
    e.kernel()
}

/// kernel_lazy == Echelon::kernel on random adversarial row orders: independent rows first, a long
/// run of redundant rows (triggers the check), then independent rows hidden near the end
/// (violated-row reinsertion); plus no-equation, full-rank and kernel dimension > 128 cases
#[test]
fn kernel_lazy_property() {
    let mut rng = StdRng::seed_from_u64(2026_1111);
    let (mut early, mut reins, mut n) = (0, 0, 0);
    for case in 0..160 {
        let ncols = rng.gen_range(1..160usize);
        let mut eqs: Vec<V> = vec![];
        let kind = case % 8;
        if kind == 0 {
            // no equations
        } else if kind == 1 {
            // full rank: a triangular system plus redundant noise
            for j in 0..ncols {
                let mut r = rand_row(&mut rng, ncols, 3);
                r.retain(|e| e.0 as usize > j);
                r.insert(0, (j as u32, F31::new(rng.gen_range(1..1000u64))));
                eqs.push(r);
            }
            let base = eqs.clone();
            for _ in 0..rng.gen_range(0..150) {
                eqs.push(combo(&mut rng, &base, ncols));
            }
            let len = eqs.len();
            for k in 0..len {
                let j = rng.gen_range(k..len);
                eqs.swap(k, j);
            }
        } else {
            let nind = rng.gen_range(0..ncols.max(1));
            let mut ind: Vec<V> = (0..nind).map(|_| { let z = rng.gen_range(1..6); rand_row(&mut rng, ncols, z) }).collect();
            let late = if kind >= 4 { rng.gen_range(1..5).min(ind.len()) } else { 0 };
            let tail: Vec<V> = ind.split_off(ind.len() - late);
            eqs.extend(ind.iter().cloned());
            if !ind.is_empty() {
                for _ in 0..rng.gen_range(60..240) {
                    eqs.push(combo(&mut rng, &ind, ncols));
                }
            }
            for t in tail {
                let at = rng.gen_range(eqs.len().saturating_sub(5)..=eqs.len());
                eqs.insert(at, t);
            }
            if kind == 7 && !eqs.is_empty() {
                // shuffled order
                let len = eqs.len();
                for k in 0..len {
                    let j = rng.gen_range(k..len);
                    eqs.swap(k, j);
                }
            }
        }
        eqs.retain(|r| !r.is_empty());
        let (got, tr) = kernel_lazy_traced(ncols, &eqs);
        let want = reference(ncols, &eqs);
        assert_eq!(got.len(), want.len(), "case {}: kernel dimension (ncols {}, {} rows)", case, ncols, eqs.len());
        assert_eq!(format!("{:?}", got), format!("{:?}", want), "case {}: kernel bases differ", case);
        if kind == 1 {
            assert!(got.is_empty(), "case {}: full rank but kernel nonempty", case);
        }
        if kind == 0 {
            assert_eq!(got.len(), ncols);
        }
        early += tr.early_stop as usize;
        reins += (tr.reinserted > 0) as usize;
        n += 1;
    }
    // a kernel of dimension > 128 (the check is off): many columns, few equations
    let ncols = 400;
    let eqs: Vec<V> = (0..300).map(|k| if k < 100 { rand_row(&mut rng, ncols, 4) } else { vec![] }).filter(|r: &V| !r.is_empty()).collect();
    let mut eqs2 = eqs.clone();
    for _ in 0..200 {
        eqs2.push(combo(&mut rng, &eqs, ncols));
    }
    let (got, tr) = kernel_lazy_traced(ncols, &eqs2);
    assert!(!tr.early_stop && tr.checks == 0, "check ran with kernel dim > 128: {:?}", tr);
    assert_eq!(format!("{:?}", got), format!("{:?}", reference(ncols, &eqs2)));
    eprintln!("[kernel_lazy] {} random systems equal Echelon::kernel; {} early stops, {} with reinsertion", n, early, reins);
    assert!(early >= 10, "too few early stops: {}", early);
    assert!(reins >= 5, "too few reinsertions: {}", reins);
}

/// M5c: the deduplicated equations (rows normalized to leading coefficient 1, duplicates
/// dropped) and the heap-driven sparse sweep of `kernel_lazy` give the reference bases on
/// windows with a zero label and with label 4 (the shapes of the (3,3,3) intermediate blocks,
/// where most rows are duplicates and rows span many columns), and on a both-sided 3-strand
/// window
#[test]
fn m5c_dedup_sparse_sweep_equals_reference() {
    let mut total = 0;
    total += check(&[2, 3, 0], "E0.1,F1.1", "F1.1,E0.1", -1..=3);
    total += check(&[4, 2, 0], "E0.1,F1.1", "F1.1,E0.1", -1..=2);
    total += check(&[1, 2, 2], "E1.1,E0.2,E1.1", "E0.1,E1.2,E0.1", -1..=3);
    total += check(&[1, 2, 2], "E0.1,E1.2,E0.1", "E1.1,E0.2,E1.1", -1..=3);
    assert!(total >= 6, "too few maps compared: {}", total);
}
