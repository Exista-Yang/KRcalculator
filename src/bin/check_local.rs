// 6B acceptance: localization matrices invertible, top action diagonal, maps block-diagonal.
use krhom::crossing::{crossing_cx, Gen};
use krhom::field::{Field, F31};
use krhom::local::*;
use krhom::mpoly::MonoCache;
use krhom::sbim::{hom_space, SBim};
use krhom::web::{random_closed_ladder, Rung};
use rand::{rngs::StdRng, Rng, SeedableRng};

fn rand_pt(rng: &mut StdRng, k: usize) -> Vec<F31> {
    (0..k).map(|_| F31::from_i64(rng.gen_range(1..1_000_000))).collect()
}

fn check_diag(b: &SBim<F31>, x: &[F31]) -> bool {
    let pts = fixed_points(b);
    let l = loc_matrix(b, &pts, x);
    let Some(linv) = inverse(&l) else { return false };
    for (i, ops) in b.top.iter().enumerate() {
        for (j, t) in ops.iter().enumerate() {
            let tm = eval_pmat(b, t, x);
            let d = matmul(&matmul(&l, &tm), &linv);
            for (p, fp) in pts.iter().enumerate() {
                let e = esym(fp.top[i], x);
                for q in 0..pts.len() {
                    let expect = if p == q { e[j + 1].clone() } else { F31::zero() };
                    if d[p][q] != expect {
                        return false;
                    }
                }
            }
        }
    }
    true
}

fn main() {
    let mut rng = StdRng::seed_from_u64(7);
    let mut ok = 0;
    let mut bad = 0;
    for &(n, m) in &[(2usize, 3usize), (3, 3), (4, 3), (4, 4), (5, 3)] {
        for _ in 0..15 {
            let bottom: Vec<usize> = (0..m).map(|_| rng.gen_range(0..=n)).collect();
            let lad = random_closed_ladder(&mut rng, n, bottom, 4);
            let b = SBim::<F31>::identity(n, &lad.bottom).add_rungs(&lad.rungs);
            if b.fsize() > 400 { continue; }
            let k: usize = lad.bottom.iter().sum();
            let x = rand_pt(&mut rng, k);
            if check_diag(&b, &x) { ok += 1 } else { bad += 1; println!("BAD {:?}", lad); }
        }
    }
    println!("diagonal top action: ok={} bad={}", ok, bad);
    // localized Rickard differentials are block diagonal (same top configuration)
    let mut blk_bad = 0;
    for n in 2..=4usize {
        for c in 1..=n {
            for d in 1..=n {
                let cx = crossing_cx(Gen { i: 0, positive: true }, c, d, n);
                let id = SBim::<F31>::identity(n, &[c, d]);
                let webs: Vec<_> = cx.terms.iter().map(|t| id.add_rungs(&t.rungs)).collect();
                let mut mc = MonoCache::new(id.ring.clone());
                for j in 0..webs.len().saturating_sub(1) {
                    let h = hom_space(&webs[j], &webs[j + 1], 1 + webs[j].shift - webs[j + 1].shift, &mut mc);
                    let x = rand_pt(&mut rng, c + d);
                    let (p, pp, f) = localize_map(&webs[j], &webs[j + 1], &h[0], &x);
                    for a in 0..pp.len() {
                        for bq in 0..p.len() {
                            if !f[a][bq].is_zero() && pp[a].top != p[bq].top {
                                blk_bad += 1;
                            }
                        }
                    }
                }
            }
        }
    }
    println!("block-diagonality violations: {}", blk_bad);
    let _ = Rung::e(0, 1);
}
