// Test hypothesis H: state space of a closed ladder = truncated closed edge ring.
use krhom::algebra::{quotient_by_span, Alg};
use krhom::field::F31;
use krhom::hom::closure_ops;
use krhom::laurent::LPoly;
use krhom::moy::closure_value;
use krhom::web::random_closed_ladder;
use rand::{rngs::StdRng, Rng, SeedableRng};

fn coinv(a: &Alg<F31>) -> LPoly {
    let ops = closure_ops(a);
    let mut gens = vec![];
    for (op, _) in &ops {
        for c in &op.cols {
            if !c.is_empty() {
                gens.push(c.clone());
            }
        }
    }
    let (keep, _) = quotient_by_span(&a.deg, gens);
    let mut co = LPoly::zero();
    for &i in &keep {
        co.add_term(a.deg[i], 1);
    }
    co
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let n: usize = args[1].parse().unwrap();
    let m: usize = args[2].parse().unwrap();
    let steps: usize = args[3].parse().unwrap();
    let count: usize = args[4].parse().unwrap();
    let mut rng = StdRng::seed_from_u64(args.get(5).map(|s| s.parse().unwrap()).unwrap_or(1));
    let (mut ok, mut bad) = (0, 0);
    for _ in 0..count {
        let bottom: Vec<usize> = (0..m).map(|_| rng.gen_range(0..=n)).collect();
        let lad = random_closed_ladder(&mut rng, n, bottom, steps);
        let a = Alg::<F31>::ladder(&lad);
        if a.dim > 30000 {
            continue;
        }
        let got = coinv(&a);
        let exp = closure_value(&lad);
        let same_shape = !exp.is_zero() && !got.is_zero() && got.shift(exp.min_deg().unwrap() - got.min_deg().unwrap()) == exp;
        if same_shape || (exp.is_zero() && got.is_zero()) {
            ok += 1;
        } else {
            bad += 1;
            if bad <= 5 {
                println!("BAD {:?} {:?}\n  got {}\n  exp {}", lad.bottom, lad.rungs.iter().map(|r| format!("{:?}{}({})", r.dir, r.i, r.r)).collect::<Vec<_>>(), got, exp);
            }
        }
    }
    println!("N={} m={} ok={} bad={}", n, m, ok, bad);
}
