// 6D acceptance: localized double leaves are genuine bimodule maps (same Q-combination of
// exact Hom basis at several points) and linearly independent.
use krhom::field::{Field, F31};
use krhom::local::localize_map;
use krhom::mpoly::MonoCache;
use krhom::sbim::{hom_space, SBim};
use krhom::soergel::{double_leaves, Ctx, LVec};
use krhom::web::Rung;
use rand::{rngs::StdRng, Rng, SeedableRng};
use std::collections::HashMap;

fn bs(n: usize, word: &[usize]) -> SBim<F31> {
    let mut b = SBim::<F31>::identity(n, &vec![1; n]);
    for &i in word {
        b = b.add_rung(&Rung::e(i, 1)).add_rung(&Rung::f(i, 1));
    }
    b
}

fn flatten(cols: &[LVec<F31>], tlen: usize) -> Vec<F31> {
    let nt = 1usize << tlen;
    let mut v = vec![F31::zero(); cols.len() * nt];
    for (c, col) in cols.iter().enumerate() {
        for (&b, x) in col {
            v[c * nt + b as usize] = x.clone();
        }
    }
    v
}

fn main() {
    let mut rng = StdRng::seed_from_u64(3);
    let mut pts: Vec<Vec<F31>> = vec![];
    let n = 4;
    for _ in 0..8 {
        pts.push((0..n).map(|_| F31::from_i64(rng.gen_range(1..1_000_000))).collect());
    }
    let mut ctxs: Vec<Ctx<F31>> = pts.iter().map(|x| Ctx::new(n, x.clone())).collect();
    let words: Vec<Vec<usize>> = vec![vec![], vec![0], vec![1], vec![0, 1], vec![1, 0], vec![0, 0], vec![0, 2], vec![2, 0], vec![0, 1, 0], vec![1, 0, 1], vec![0, 1, 2], vec![1, 2, 1], vec![0, 2, 1]];
    let (mut ok, mut bad, mut dep) = (0, 0, 0);
    for x in &words {
        for y in &words {
            let bx = bs(n, x);
            let by = bs(n, y);
            let mut mc = MonoCache::new(bx.ring.clone());
            for d in -3..=4 {
                let dls = double_leaves(&ctxs[0], x, y, d);
                if dls.is_empty() {
                    continue;
                }
                let u = d + bx.shift - by.shift;
                let basis = hom_space(&bx, &by, u, &mut mc);
                // localized concrete basis at each point, flattened in bit coordinates
                let mut flat_basis: Vec<Vec<Vec<F31>>> = vec![];
                for (pi, xpt) in pts.iter().enumerate() {
                    let _ = pi;
                    let mut fb = vec![];
                    for phi in &basis {
                        let (p, pp, f) = localize_map(&bx, &by, phi, xpt);
                        // convert to bits
                        let sb: Vec<u64> = p.iter().map(|fp| bits_of(fp, x, n)).collect();
                        let tb: Vec<u64> = pp.iter().map(|fp| bits_of(fp, y, n)).collect();
                        let mut cols: Vec<LVec<F31>> = vec![HashMap::new(); 1 << x.len()];
                        for (a, row) in f.iter().enumerate() {
                            for (b, v) in row.iter().enumerate() {
                                if !v.is_zero() {
                                    cols[sb[b] as usize].insert(tb[a], v.clone());
                                }
                            }
                        }
                        fb.push(flatten(&cols, y.len()));
                    }
                    flat_basis.push(fb);
                }
                let mut dl_flat0 = vec![];
                for (_, _, ops) in &dls {
                    let mut vals = vec![];
                    for c in ctxs.iter_mut() {
                        let (tw, cols) = c.map_columns(x, ops);
                        assert_eq!(&tw, y);
                        vals.push(flatten(&cols, y.len()));
                    }
                    dl_flat0.push(vals[0].clone());
                    // solve jointly over all points
                    let m = basis.len();
                    let dim = vals[0].len();
                    let mut rows: Vec<Vec<F31>> = vec![];
                    let mut rhs: Vec<F31> = vec![];
                    for k in 0..pts.len() {
                        for i in 0..dim {
                            rows.push((0..m).map(|j| flat_basis[k][j][i].clone()).collect());
                            rhs.push(vals[k][i].clone());
                        }
                    }
                    match krhom::dense::solve(&rows, m, &rhs) {
                        None => { bad += 1; println!("DL not genuine: {:?}->{:?} d={}", x, y, d); }
                        Some(_) => { ok += 1; }
                    }
                }
                if krhom::dense::rank(&dl_flat0, dl_flat0[0].len()) < dl_flat0.len() {
                    dep += 1;
                    println!("dependent DLs {:?}->{:?} d={}", x, y, d);
                }
            }
        }
    }
    println!("genuine={} not={} dependent-sets={}", ok, bad, dep);
}

fn bits_of(fp: &krhom::local::FixedPoint, word: &[usize], strands: usize) -> u64 {
    let mut pi: Vec<u64> = (0..strands).map(|u| 1u64 << u).collect();
    let mut bits = 0u64;
    for (j, &i) in word.iter().enumerate() {
        let a = pi[i];
        if fp.moved[2 * j + 1] == a {
            bits |= 1 << j;
            pi.swap(i, i + 1);
        }
    }
    bits
}
