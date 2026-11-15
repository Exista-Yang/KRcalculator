// Test: universal construction (pairing) state space vs MOY.
use krhom::algebra::field::F31;
use krhom::combinat::moy::closure_value;
use krhom::legacy::untrunc::{universal_state_space, UAlg};
use krhom::combinat::web::{random_closed_ladder, Dir, Ladder, Rung};
use rand::{rngs::StdRng, Rng, SeedableRng};

fn parse(n: usize, bottom: Vec<usize>, s: &str) -> Ladder {
    let mut l = Ladder::new(n, bottom);
    for tok in s.split_whitespace() {
        let dir = if tok.starts_with('E') { Dir::E } else { Dir::F };
        let i: usize = tok[1..2].parse().unwrap();
        let r: usize = tok[3..tok.len() - 1].parse().unwrap();
        l.push(Rung { i, dir, r });
    }
    l
}

fn fib_top(lad: &Ladder) -> i32 {
    let mut t = 0;
    for (lv, rg) in lad.levels().iter().zip(lad.rungs.iter()) {
        let (_, giv) = rg.ends();
        let b = lv[giv];
        t += 2 * (rg.r * (b - rg.r)) as i32;
    }
    t
}

fn test(lad: &Ladder, verbose: bool) -> bool {
    let exp = closure_value(lad);
    let n = lad.n;
    let top: i32 = 2 * lad.bottom.iter().map(|&l| (l * (n - l)) as i32).sum::<i32>();
    let cap = top + 2 * n as i32 + fib_top(lad) + 2;
    let u = UAlg::<F31>::ladder(lad, cap);
    let got = universal_state_space(&u);
    let ok = if exp.is_zero() { got.is_zero() } else { !got.is_zero() && got.shift(exp.min_deg().unwrap() - got.min_deg().unwrap()) == exp };
    if verbose || !ok {
        println!("{} {:?} {:?} dim={} cap={}\n   got {}\n   exp {}", if ok { "OK " } else { "BAD" }, lad.bottom, lad.rungs.iter().map(|r| format!("{:?}{}({})", r.dir, r.i, r.r)).collect::<Vec<_>>(), u.alg.dim, cap, got, exp);
    }
    ok
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args[1] == "one" {
        let n: usize = args[2].parse().unwrap();
        let bottom: Vec<usize> = args[3].split(',').map(|x| x.parse().unwrap()).collect();
        test(&parse(n, bottom, &args[4]), true);
        return;
    }
    let n: usize = args[1].parse().unwrap();
    let m: usize = args[2].parse().unwrap();
    let steps: usize = args[3].parse().unwrap();
    let count: usize = args[4].parse().unwrap();
    let mut rng = StdRng::seed_from_u64(args.get(5).map(|s| s.parse().unwrap()).unwrap_or(1));
    let (mut ok, mut bad) = (0, 0);
    for _ in 0..count {
        let bottom: Vec<usize> = (0..m).map(|_| rng.gen_range(0..=n)).collect();
        let lad = random_closed_ladder(&mut rng, n, bottom, steps);
        if test(&lad, false) { ok += 1 } else { bad += 1 }
    }
    println!("N={} m={} ok={} bad={}", n, m, ok, bad);
}
