// Universal construction with the polynomial-matrix SBim implementation vs MOY.
use krhom::field::F31;
use krhom::moy::closure_value;
use krhom::mpoly::MonoCache;
use krhom::sbim::{SBim, StateSpace};
use krhom::laurent::LPoly;
use krhom::web::{random_closed_ladder, Dir, Ladder, Rung};
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

fn test(lad: &Ladder, verbose: bool) -> bool {
    let exp = closure_value(lad);
    let b = SBim::<F31>::identity(lad.n, &lad.bottom).add_rungs(&lad.rungs);
    let mut mc = MonoCache::new(b.ring.clone());
    let ss = StateSpace::compute(&b, &mut mc);
    let sh = b.shift - lad.bottom.iter().map(|&l| (l * (lad.n - l)) as i32).sum::<i32>();
    let mut got = LPoly::zero();
    for (d, p) in &ss.pieces {
        got.add_term(d + sh, p.reps.len() as i128);
    }
    let ok = got == exp;
    if verbose || !ok {
        println!("{} {:?} {:?} fsize={}\n   got {}\n   exp {}", if ok { "OK " } else { "BAD" }, lad.bottom, lad.rungs.iter().map(|r| format!("{:?}{}({})", r.dir, r.i, r.r)).collect::<Vec<_>>(), b.fsize(), got, exp);
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
        if SBim::<F31>::identity(lad.n, &lad.bottom).add_rungs(&lad.rungs).fsize() > 64 { continue; }
        if test(&lad, false) { ok += 1 } else { bad += 1 }
    }
    println!("N={} m={} ok={} bad={}", n, m, ok, bad);
}
