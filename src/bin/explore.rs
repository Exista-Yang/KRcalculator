use krhom::algebra::Alg;
use krhom::field::F31;
use krhom::hom::closure_graded_dim;
use krhom::moy::closure_value;
use krhom::web::{Dir, Ladder, Rung};

fn all_ladders(n: usize, m: usize, len: usize) -> Vec<Ladder> {
    let mut out = vec![];
    let mut bottoms = vec![vec![]];
    for _ in 0..m {
        let mut nb = vec![];
        for b in &bottoms {
            for l in 0..=n {
                let mut x: Vec<usize> = b.clone();
                x.push(l);
                nb.push(x);
            }
        }
        bottoms = nb;
    }
    for b in bottoms {
        let mut stack = vec![Ladder::new(n, b.clone())];
        while let Some(l) = stack.pop() {
            if l.top() == l.bottom && !l.rungs.is_empty() {
                out.push(l.clone());
            }
            if l.rungs.len() == len {
                continue;
            }
            for i in 0..m - 1 {
                for dir in [Dir::E, Dir::F] {
                    for r in 1..=n {
                        let rg = Rung { i, dir, r };
                        if rg.apply(&l.top(), n).is_some() {
                            let mut l2 = l.clone();
                            l2.push(rg);
                            stack.push(l2);
                        }
                    }
                }
            }
        }
    }
    out
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let n: usize = args.get(1).map(|s| s.parse().unwrap()).unwrap_or(2);
    let m: usize = args.get(2).map(|s| s.parse().unwrap()).unwrap_or(2);
    let len: usize = args.get(3).map(|s| s.parse().unwrap()).unwrap_or(2);
    let mut shown = 0;
    let mut ok = 0;
    let mut lads = all_ladders(n, m, len);
    lads.sort_by_key(|l| l.rungs.len());
    for l in lads {
        let a = Alg::<F31>::ladder(&l);
        let got = {
            let ops = krhom::hom::closure_ops(&a);
            let mut gens = vec![];
            for (op, _) in &ops { for c in &op.cols { if !c.is_empty() { gens.push(c.clone()); } } }
            let (keep, _) = krhom::algebra::quotient_by_span(&a.deg, gens);
            let mut co = krhom::laurent::LPoly::zero();
            for &i in &keep { co.add_term(a.deg[i], 1); }
            let exp0 = closure_value(&l);
            if !co.is_zero() && !exp0.is_zero() { co.shift(exp0.min_deg().unwrap() - co.min_deg().unwrap()) } else { co }
        };
        let exp = closure_value(&l);
        if got != exp {
            if shown < 12 {
                let kind = if got.eval_at_one() == exp.eval_at_one() {
                    if exp.min_deg().is_some() && got.shift(exp.min_deg().unwrap() - got.min_deg().unwrap()) == exp { "SHIFT" } else { "SHAPE" }
                } else {
                    "DIM"
                };
                println!("{} bottom={:?} rungs={:?}\n   got {}\n   exp {}", kind, l.bottom, l.rungs.iter().map(|r| format!("{:?}{}({})", r.dir, r.i, r.r)).collect::<Vec<_>>(), got, exp);
                shown += 1;
            }
        } else {
            ok += 1;
        }
    }
    println!("ok = {}", ok);
}
