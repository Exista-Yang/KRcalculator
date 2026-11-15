// Test model H2: F_N(closure W) = Hom_SSBim(R,B_W) / J Hom_SSBim(R,B_W).
use krhom::algebra::field::F31;
use krhom::combinat::moy::closure_value;
use krhom::legacy::untrunc::UAlg;
use krhom::combinat::web::{random_closed_ladder, Ladder, Rung, Dir};
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
    let width = exp.max_deg().unwrap_or(0) - exp.min_deg().unwrap_or(0);
    // try increasing caps until the candidate's top part is reliable
    let maxlab = lad.levels().iter().flatten().copied().max().unwrap_or(0) as i32;
    let mut cap = width + 2 * maxlab + 6;
    loop {
        let u = UAlg::<F31>::ladder(lad, cap);
        let (got, rel) = u.candidate_state_space();
        // got should be exp shifted; check degrees up to rel
        let lo = got.min_deg();
        let ok = match lo {
            None => exp.is_zero(),
            Some(lo) => {
                let sh = exp.min_deg().unwrap() - lo;
                let g2 = got.shift(sh);
                let top_got = got.max_deg().unwrap();
                // reliable if rel exceeds top of expected shifted back by 2
                let need = exp.max_deg().unwrap() - sh + 2;
                if rel < need {
                    cap += 6;
                    if cap > 80 { println!("cap too large"); return false; }
                    continue;
                }
                let _ = top_got;
                g2 == exp
            }
        };
        if verbose || !ok {
            println!("{} {:?} {:?} dim={} cap={}\n   got {}\n   exp {}", if ok {"OK "} else {"BAD"}, lad.bottom, lad.rungs.iter().map(|r| format!("{:?}{}({})", r.dir, r.i, r.r)).collect::<Vec<_>>(), u.alg.dim, cap, got, exp);
        }
        return ok;
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() > 1 && args[1] == "one" {
        let n: usize = args[2].parse().unwrap();
        let bottom: Vec<usize> = args[3].split(',').map(|x| x.parse().unwrap()).collect();
        let lad = parse(n, bottom, &args[4]);
        test(&lad, true);
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
    let _ = Rung::e(0,1);
}
