// 6E acceptance: idempotents and decompositions for S_3, S_4.
use krhom::algebra::field::{Field, F31};
use krhom::combinat::hecke::Hecke;
use krhom::soergel::indec::*;
use krhom::soergel::soergel::Ctx;
use rand::{rngs::StdRng, Rng, SeedableRng};
use std::collections::HashMap;

fn main() {
    for n in [3usize, 4] {
        let mut rng = StdRng::seed_from_u64(11 + n as u64);
        let x: Vec<F31> = (0..n).map(|_| F31::from_i64(rng.gen_range(1..1_000_000))).collect();
        let mut ctx = Ctx::<F31>::new(n, x);
        let h = Hecke::new(n);
        let t0 = std::time::Instant::now();
        let mut indecs: HashMap<usize, Indec<F31>> = HashMap::new();
        let mut bad = 0;
        for xe in 0..ctx.g.size() {
            let ind = indecomposable(&mut ctx, xe);
            let e2 = ind.e.compose(&ind.e);
            let mut diff = e2.clone();
            diff.add_scaled(&ind.e, &F31::one().neg());
            if !diff.is_zero() { bad += 1; println!("x={} not idempotent", xe); }
            for w in 0..ctx.g.size() {
                let r = block_rank(&ctx, &ind.e, w);
                let expect = h.h(w, xe).eval_at_one() as usize;
                if r != expect { bad += 1; println!("x={} w={} block rank {} expected {}", xe, w, r, expect); }
            }
            indecs.insert(xe, ind);
        }
        println!("S_{}: {} indecomposables, errors {}, time {:.2}s", n, indecs.len(), bad, t0.elapsed().as_secs_f64());
        let t1 = std::time::Instant::now();
        let mut dbad = 0;
        let mut count = 0;
        for xe in 0..ctx.g.size() {
            for s in 0..n - 1 {
                let ix = &indecs[&xe];
                let mut word = ix.word.clone();
                word.push(s);
                let e = ix.e.tensor_letter(s);
                let mult: Vec<(usize, i32, usize)> = kl_mult(&h, xe, s).into_iter().map(|(y, k, m)| (y, -k, m)).collect();
                let sums = decompose(&mut ctx, &word, &e, &mult, &indecs);
                // Σ ι π = e
                let mut tot = LMat::zero(&word, &word);
                for sm in &sums { tot.add_scaled(&sm.iota.compose(&sm.pi), &F31::one()); }
                tot.add_scaled(&e, &F31::one().neg());
                if !tot.is_zero() { dbad += 1; println!("decomposition incomplete x={} s={} ({} summands)", xe, s, sums.len()); }
                // orthogonality
                for (a, sa) in sums.iter().enumerate() {
                    for (b, sb) in sums.iter().enumerate() {
                        let c = sa.pi.compose(&sb.iota);
                        let ok = if a == b { let mut d = c.clone(); d.add_scaled(&indecs[&sa.y].e, &F31::one().neg()); d.is_zero() } else { c.is_zero() };
                        if !ok { dbad += 1; println!("orthogonality fails x={} s={} a={} b={}", xe, s, a, b); }
                    }
                }
                count += 1;
            }
        }
        println!("S_{}: {} decompositions checked, errors {}, time {:.2}s", n, count, dbad, t1.elapsed().as_secs_f64());
    }
}
