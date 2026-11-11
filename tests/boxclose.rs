//! MBOX: the box closure layer at fixed N (src/boxclose.rs, the staged rank stage of
//! src/cclose.rs). SPEC_3C_GL6 rev1 §2 MBOX tests (a)-(d); (e) is tests/boxclose_mem.rs.
//!
//! * (a) Synthetic sparse homogeneous R^λ polynomials for (3,3), (2,2,2), (3,3,3), weights 0..top+1,
//!   N = 4..6, both modes, both primes: box extraction = `GrAlg::nf`, also restricted to every
//!   N' < N (nesting).
//! * (b) End to end on exact xengine closure data: extraction of Gram and P entries at box points,
//!   Freivalds certificate, scaling check, T0 against MOY, `close_box` (uncapped and capped) =
//!   `close_exact_colored` for N = 3..6, both modes, both extraction modes: (2,2) σ^3 and (3,2) σ^2;
//!   (3,3) σ^1..σ^3 `#[ignore]`.
//! * (c) Negative tests: each must fail as intended.
//! * (d) Shuffled point order and shard splits (through checkpoints) give identical results.
//! * The staged rank stage on P equals the exact oracle on F (and its type dims) on exact data.

use krhom::boxclose::*;
use krhom::cclose::*;
use krhom::crossing::Gen;
use krhom::field::{Field, F31, F32b};
use krhom::laurent::{LPoly, Poincare};
use krhom::xengine::*;
use rand::{rngs::StdRng, seq::SliceRandom, Rng, SeedableRng};
use std::collections::BTreeMap;
use std::sync::Arc;
use std::time::Instant;

// ------------------------------------------------------------------ helpers

fn gens(w: &[i32]) -> Vec<Gen> {
    w.iter().map(|&g| Gen { i: (g.unsigned_abs() - 1) as usize, positive: g > 0 }).collect()
}

fn rand_poly<F: PrimeField>(ring: &ERing, w: u32, nterms: usize, rng: &mut StdRng) -> EPoly<F> {
    let monos = ring.monos_of_weight(w);
    let mut t = vec![];
    for _ in 0..nterms {
        let m = monos[rng.gen_range(0..monos.len())];
        t.push((m, F::from_u64(rng.gen_range(1..F::characteristic()))));
    }
    EPoly::from_terms(t)
}

/// Extract `polys` (declared weights) at every needed point of `grid`, in the given order.
fn extract<F: PrimeField>(grid: &Arc<BoxGrid<F>>, polys: &[EPoly<F>], weights: &[Option<u32>], order: Option<&[u32]>) -> Result<BoxResult<F>, String> {
    let mut acc = BoxAcc::new(grid.clone(), weights.to_vec(), 1, 2, None);
    let pts: Vec<u32> = match order {
        Some(o) => o.to_vec(),
        None => grid.points_up_to(acc.needed_size().unwrap_or(0)),
    };
    for id in pts {
        let x = grid.coords(id);
        let vals: Vec<F> = polys.iter().map(|p| eval_epoly(&grid.ring, p, &x)).collect();
        acc.add_point(id, &vals)?;
    }
    acc.finish()
}

// ------------------------------------------------------------------ (a) synthetic

fn synthetic_case<F: PrimeField>(lam: &[u8], n: usize, marked: Option<usize>, seed: u64) -> usize {
    let ring = ERing::new(lam);
    let mut rng = StdRng::seed_from_u64(seed);
    let grid = Arc::new(BoxGrid::new(BoxNodes::<F>::new(lam, n, seed.wrapping_mul(31) + 7), marked));
    let top = grid.top;
    let mut polys = vec![];
    let mut weights = vec![];
    for w in 0..=top + 1 {
        for k in 0..2 {
            polys.push(rand_poly::<F>(&ring, w, 2 + 3 * k, &mut rng));
            weights.push(Some(w));
        }
    }
    let res = extract(&grid, &polys, &weights, None).unwrap_or_else(|e| panic!("{:?} N={} marked {:?}: {}", lam, n, marked, e));
    let mut cmp = 0;
    for np in lam.iter().map(|&x| x as usize).max().unwrap()..=n {
        let alg = GrAlg::<F>::new(lam, np, marked);
        for (e, p) in polys.iter().enumerate() {
            let got = schur_to_epoly(&grid, &res.coords[e], np, marked);
            assert_eq!(alg.nf(&got), alg.nf(p), "{:?} N={} (extracted at {}), marked {:?}, entry {} of weight {}", lam, np, n, marked, e, weights[e].unwrap());
            cmp += 1;
        }
    }
    cmp
}

#[test]
fn a_synthetic_extraction_equals_nf() {
    let t = Instant::now();
    let mut cmp = 0;
    for (lam, nmin) in [(vec![3u8, 3], 4usize), (vec![2, 2, 2], 4), (vec![3, 3, 3], 4)] {
        for n in nmin..=6 {
            for marked in [None, Some(0)] {
                let seed = (n * 10 + marked.map_or(0, |m| m + 1)) as u64 + 1000 * lam.len() as u64 + lam[0] as u64;
                cmp += synthetic_case::<F31>(&lam, n, marked, seed);
                cmp += synthetic_case::<F32b>(&lam, n, marked, seed + 1);
            }
        }
    }
    eprintln!("(a) {} extracted entries compared with nf, {:.1}s", cmp, t.elapsed().as_secs_f64());
}

#[test]
fn a_grid_basics() {
    // point counts: L(top) = product of binomials; reduced: block 0 only ∅
    let g = BoxGrid::new(BoxNodes::<F31>::new(&[3, 3, 3], 6, 5), None);
    assert_eq!(g.npoints, 8000);
    assert_eq!(g.top, 27);
    assert_eq!((0..=27).map(|s| g.count_of_size(s)).max(), Some(777));
    let r = BoxGrid::new(BoxNodes::<F31>::new(&[3, 3, 3], 6, 5), Some(0));
    assert_eq!(r.npoints, 400);
    assert_eq!(r.top, 18);
    assert_eq!((0..=18).map(|s| r.count_of_size(s)).max(), Some(48));
    // points with |μ| <= 10 at N = 6 unreduced (R4 §2: 1 834)
    assert_eq!(g.points_up_to(10).len(), 1834);
    assert_eq!(r.points_up_to(10).len(), 269);
    assert_eq!(box_points(&r, 18).len(), 400);
    // point coordinates: e_j of the roots a_{μ_j + a + 1 − j}
    let nodes = BoxNodes::<F31>::new(&[2], 4, 9);
    let x = point_coords(&nodes, &[vec![2, 1]]);
    let (r1, r2) = (nodes.a[0][3], nodes.a[0][1]);
    assert_eq!(x, vec![r1.add(&r2), r1.mul(&r2)]);
    // nodes are distinct across blocks; scaled nodes scale points
    let nb = BoxNodes::<F31>::new(&[3, 3, 3], 6, 1);
    let mut all: Vec<u64> = nb.a.iter().flatten().map(|x| x.0).collect();
    all.sort();
    all.dedup();
    assert_eq!(all.len(), 18);
    assert_ne!(BoxNodes::<F31>::new(&[3, 3], 6, 1).hash(), BoxNodes::<F31>::new(&[3, 3], 6, 2).hash());
    assert_ne!(BoxNodes::<F31>::new(&[3, 3], 6, 1).hash(), BoxNodes::<F32b>::new(&[3, 3], 6, 1).hash());
}

// ------------------------------------------------------------------ (b) end to end on xengine data

struct BoxRun<F> {
    data: BoxClosureData<F>,
    kinds: Vec<EntryKind>,
    weights: Vec<Option<u32>>,
    points: usize,
    entries: usize,
    coefs: usize,
    secs: f64,
}

/// Box extraction of the Gram and P entries of exact closure data (the stand-in for scans at the
/// box points), assembled into BoxClosureData.
fn box_from_exact<F: PrimeField>(data: &ClosureDataC<F>, pairing: &Pairing<F>, n_box: usize, mode: Option<usize>, prune_n: Option<usize>, seed: u64, weight_shift: Option<usize>) -> Result<BoxRun<F>, String> {
    let t = Instant::now();
    let edges: Vec<(usize, usize)> = pairing.iter().map(|(s, t, _)| (*s, *t)).collect();
    let (kinds, mut weights) = closure_entries(data, &edges);
    if let Some(e) = weight_shift {
        weights[e] = weights[e].map(|w| w + 1);
    }
    let grid = Arc::new(BoxGrid::new(BoxNodes::<F>::new(&data.lam, n_box, seed), mode));
    let mut acc = BoxAcc::new(grid.clone(), weights.clone(), 0, 0, prune_n);
    let pts = grid.points_up_to(acc.needed_size().unwrap_or(0));
    for &id in &pts {
        let vals = eval_entries(data, pairing, &kinds, &grid.coords(id));
        acc.add_point(id, &vals)?;
    }
    let res = acc.finish()?;
    let coefs = acc.acc_len();
    let stamp = BoxStamp { prime: F::characteristic(), node_seed_hash: grid.nodes.hash(), library_hash: 0, word_hash: 0, code_stamp: BOX_CODE_STAMP.into() };
    let bd = BoxClosureData::from_box(data, &edges, &kinds, &res, prune_n, stamp)?;
    Ok(BoxRun { data: bd, entries: kinds.len(), kinds, weights, points: pts.len(), coefs, secs: t.elapsed().as_secs_f64() })
}

/// MOY values of the final types of a scan (types matched to tables through their ladders).
fn moy_of_scan<F: Field>(sc: &XScan<F>, n: usize) -> BTreeMap<u32, LPoly> {
    let lam: Vec<u8> = sc.labels.iter().map(|&a| a as u8).collect();
    sc.type_ladders
        .iter()
        .map(|(&p, rungs)| {
            let fmt = sc.scan_ladders.iter().find(|(_, r)| *r == rungs).map(|(f, _)| f.clone()).expect("type ladder not in the scan");
            (p, moy_type_trace(&lam, &fmt, n))
        })
        .collect()
}

fn end_to_end<F: PrimeField>(labels: &[usize], word: &[i32], ns: std::ops::RangeInclusive<usize>, seed: u64) {
    let t0 = Instant::now();
    let sc = xscan::<F>(labels, &gens(word), &XOpts::default()).unwrap();
    let ts = t0.elapsed().as_secs_f64();
    let data = &sc.data;
    let pairing = pairing_from_f(data);
    let n_box = *ns.end();
    let mut lines = vec![];
    for mode in [None, Some(0)] {
        let mut run = box_from_exact(data, &pairing, n_box, mode, None, seed, None).unwrap();
        // scaling check at two points (D14, C_sound M1)
        let nodes = BoxNodes::<F>::new(&data.lam, n_box, seed);
        let lambda = F::from_u64(0x1234_5677);
        let sn = nodes.scaled(&lambda);
        let grid = BoxGrid::new(nodes.clone(), mode);
        let mut rep = ScalingReport::default();
        for id in [grid.npoints / 3, grid.npoints - 1] {
            let mu = grid.mu(id);
            let vz = eval_entries(data, &pairing, &run.kinds, &point_coords(&nodes, &mu));
            let vlz = eval_entries(data, &pairing, &run.kinds, &point_coords(&sn, &mu));
            let r = check_scaling_all(&vz, &vlz, &run.weights, &lambda).unwrap();
            rep.checked += r.checked;
            rep.zero_both += r.zero_both;
        }
        for nn in ns.clone() {
            for marked in [None, Some(0)] {
                if mode.is_some() && marked != mode {
                    assert!(close_box(&run.data, nn, marked, None, None).is_err(), "reduced box data must refuse the unreduced closure");
                    continue;
                }
                let want = close_exact_colored(data, nn, marked);
                let cert = closed_d2_freivalds(&mut run.data, nn, marked, 3, seed + nn as u64).unwrap();
                let got = close_box(&run.data, nn, marked, None, None).unwrap();
                assert_eq!(got, want, "{:?} {:?} N={} marked {:?} (box mode {:?})", labels, word, nn, marked, mode);
                // T0 against MOY, then the capped closure
                if marked.is_none() || mode.is_some() || nn == n_box {
                    let moy = moy_of_scan(&sc, nn);
                    let moy = if marked.is_some() {
                        // reduced: the MOY value of the type divided by [N choose λ_0] (block 0 marked)
                        moy.into_iter().map(|(p, v)| (p, v.div_exact(&LPoly::qbinom(nn as i64, labels[0] as i64)))).collect()
                    } else {
                        moy
                    };
                    let t0r = box_moy_check(&run.data, nn, marked, &moy, None).unwrap_or_else(|e| panic!("{:?} {:?} N={} marked {:?}: {}", labels, word, nn, marked, e));
                    assert_eq!(close_box(&run.data, nn, marked, Some(&t0r), None).unwrap(), want, "capped closure");
                }
                lines.push(format!("N={} marked {:?} mode {:?}: {} (rank {}), D2 blocks {}", nn, marked, mode, got, got.total_rank(), cert.blocks));
            }
        }
        eprintln!(
            "(b) {:?} {:?} over {} box mode {:?}: {} entries, {} points, {} accumulator coefficients, extraction {:.2}s; scaling: {} checked, {} zero at both points",
            labels,
            word,
            F::name(),
            mode,
            run.entries,
            run.points,
            run.coefs,
            run.secs,
            rep.checked,
            rep.zero_both
        );
    }
    for l in lines {
        eprintln!("    {}", l);
    }
    eprintln!("(b) {:?} {:?} over {}: scan {:.1}s, total {:.1}s", labels, word, F::name(), ts, t0.elapsed().as_secs_f64());
}

#[test]
fn b_end_to_end_22_s3() {
    end_to_end::<F31>(&[2, 2], &[1, 1, 1], 3..=6, 11);
    end_to_end::<F32b>(&[2, 2], &[1, 1, 1], 3..=6, 12);
}

#[test]
fn b_end_to_end_32_s2() {
    end_to_end::<F31>(&[3, 2], &[1, 1], 3..=6, 13);
}

#[test]
#[ignore = "MBOX (b): (3,2) sigma^2 over the second prime"]
fn b_end_to_end_32_s2_f32b() {
    end_to_end::<F32b>(&[3, 2], &[1, 1], 3..=6, 14);
}

#[test]
#[ignore = "MBOX (b): (3,3) sigma^1..sigma^3, N = 3..6 (cloud, minutes)"]
fn b_end_to_end_33() {
    for w in [vec![1], vec![1, 1], vec![1, 1, 1]] {
        end_to_end::<F31>(&[3, 3], &w, 3..=6, 15);
    }
    end_to_end::<F32b>(&[3, 3], &[1, 1, 1], 3..=6, 16);
}

#[test]
fn staged_equals_oracle_on_exact_data() {
    for (labels, word) in [(vec![2usize, 2], vec![1, -1, 1]), (vec![2, 1], vec![1, 1]), (vec![1, 2], vec![1, 1, 1, 1])] {
        let sc = xscan::<F31>(&labels, &gens(&word), &XOpts::default()).unwrap();
        let p = pairing_from_f(&sc.data);
        for nn in 2..=5 {
            for marked in [None, Some(0)] {
                let st = try_close_staged(&sc.data, &p, nn, marked, &StagedOpts::default()).unwrap();
                assert_eq!(st, close_exact_colored(&sc.data, nn, marked), "{:?} {:?} N={} marked {:?}", labels, word, nn, marked);
                let dims = try_staged_type_ranks(&sc.data, nn, marked, &StagedOpts::default()).unwrap();
                let ex = exact_type_dims_colored(&sc.data, nn, marked);
                for (x, d) in &dims {
                    assert_eq!(d.values().sum::<usize>(), ex[x]);
                }
                assert!(closed_d2_freivalds_raw(&sc.data, &p, nn, marked, &StagedOpts::default(), 2, 1).is_ok());
            }
        }
    }
}

// ------------------------------------------------------------------ (c) negative tests

fn small_run(mode: Option<usize>, prune: Option<usize>) -> (XScan<F31>, Pairing<F31>, BoxRun<F31>) {
    let opts = XOpts { prune_n: prune, ..Default::default() };
    let sc = xscan::<F31>(&[2, 2], &gens(&[1, 1, 1]), &opts).unwrap();
    let p = pairing_from_f(&sc.data);
    let run = box_from_exact(&sc.data, &p, 6, mode, prune, 21, None).unwrap();
    (sc, p, run)
}

#[test]
fn c_freivalds_detects_perturbed_p_entry() {
    let (sc, _p, run) = small_run(None, None);
    let ring = sc.data.ring();
    let nn = 4;
    let mut base = run.data.clone();
    assert!(closed_d2_freivalds(&mut base, nn, None, 3, 5).is_ok());
    // perturb one P entry of an edge whose target has outgoing edges: add c·e_{1,1}^w (nonzero in A_λ)
    let tgt_out: Vec<usize> = run.data.pairing.iter().map(|(s, _, _)| *s).collect();
    let mut failed = 0;
    let mut tried = 0;
    for ei in 0..run.data.pairing.len() {
        let (s, t, _) = run.data.pairing[ei];
        let _ = s;
        if !tgt_out.contains(&t) {
            continue;
        }
        let (os, ot) = (run.data.base.objs[s], run.data.base.objs[t]);
        let (ts, tt) = (&run.data.base.types[&os.0], &run.data.base.types[&ot.0]);
        let d = os.2 - ot.2 + ts.ins[0] + tt.outs[0];
        if d < 0 || d % 2 != 0 || d / 2 > 2 * (nn as i32 - 2) {
            continue;
        }
        let mut ex = vec![0u32; ring.nvars];
        ex[ring.var(1, 1)] = (d / 2) as u32;
        let mut bad = run.data.clone();
        bad.pairing[ei].2[0][0].add_assign(&EPoly::mono(ring.pack(&ex), F31::from_i64(7)));
        tried += 1;
        let r = closed_d2_freivalds(&mut bad, nn, None, 3, 5);
        if r.is_err() {
            failed += 1;
            // and the old certificate does not carry over to the perturbed pairing
            let mut stale = bad.clone();
            stale.d2_cert = base.d2_cert.clone();
            stale.prune_n = Some(6);
            assert!(close_box(&stale, nn, None, None, None).is_err());
        }
    }
    eprintln!("(c) perturbed P entries: {} of {} detected by Freivalds", failed, tried);
    assert!(tried > 0 && failed == tried, "Freivalds missed a perturbed P entry ({} of {})", failed, tried);
}

#[test]
fn c_freivalds_detects_reinserted_d2_component() {
    // The (3,3) scan of σσσσ^{-1}σ pruned at N = 4 leaves a d∘d component f: x -> w (h -> h+2)
    // in the ideal of types killed at 4 (`final_d2`); like every final component of the (3,3)
    // prune-4 scans of length <= 8 it closes to 0 at N' = 4, 5, 6 (probe in MBOX_report.md).
    // Re-inserted as a 3-object complex x' -> w' -> w'' (edges f and id_w) appended to the pruned
    // closure data, the closed D² is exactly Φ_{N'}(f): Freivalds must pass. Re-inserted with a
    // term outside the ideal (f + c·e_{1,1}^k, a pruning bug), Φ(f') != 0 and Freivalds must fail,
    // in agreement with `closure_kills_map` both times; through the box, the same holds for
    // `closed_d2_freivalds`, and `close_box` then refuses the uncertified pruned data.
    let (labels, word, prune) = (vec![3usize, 3], vec![1, 1, 1, -1, 1], 4usize);
    let sc = xscan::<F31>(&labels, &gens(&word), &XOpts { prune_n: Some(prune), ..Default::default() }).unwrap();
    assert!(!sc.final_d2.is_empty());
    let ring = sc.data.ring();
    let (x, w, f) = sc.final_d2[0].clone();
    let (ox, ow) = (sc.data.objs[x], sc.data.objs[w]);
    let (tx, tw) = (sc.data.types[&ox.0].clone(), sc.data.types[&ow.0].clone());
    // f' = f + c·e_{1,1}^k on one entry of degree 2k, k <= 3 (e_1^k != 0 in H*(Gr(3,4)))
    let mut fbad = f.clone();
    let mut done = false;
    'o: for c in 0..tw.ins.len() {
        for a in 0..tx.ins.len() {
            let d = ox.2 - ow.2 + tx.ins[a] - tw.ins[c];
            if d >= 0 && d % 2 == 0 && d / 2 <= 3 {
                let mut ex = vec![0u32; ring.nvars];
                ex[ring.var(1, 1)] = (d / 2) as u32;
                fbad[c][a].add_assign(&EPoly::mono(ring.pack(&ex), F31::from_i64(5)));
                done = true;
                break 'o;
            }
        }
    }
    assert!(done);
    let reinsert = |g: &Vec<Vec<EPoly<F31>>>| -> ClosureDataC<F31> {
        let mut data = sc.data.clone();
        let hb = data.objs.iter().map(|o| o.1).max().unwrap() + 10;
        let base = data.objs.len();
        data.objs.push((ox.0, hb, ox.2));
        data.objs.push((ow.0, hb + 1, ow.2));
        data.objs.push((ow.0, hb + 2, ow.2));
        let nw = tw.ins.len();
        let id: Vec<Vec<EPoly<F31>>> = (0..nw).map(|c| (0..nw).map(|a| if a == c { EPoly::one() } else { EPoly::zero() }).collect()).collect();
        data.edges.push((base, base + 1, g.clone()));
        data.edges.push((base + 1, base + 2, id));
        data.validate().unwrap();
        data
    };
    for nn in prune..=6 {
        assert!(closure_kills_map(&sc.data, x, w, &f, nn, None), "final d∘d component closes at N' = {}", nn);
    }
    for (g, closes) in [(&f, true), (&fbad, false)] {
        assert_eq!(closure_kills_map(&sc.data, x, w, g, prune, None), closes);
        let data = reinsert(g);
        let p = pairing_from_f(&data);
        let fr = closed_d2_freivalds_raw(&data, &p, prune, None, &StagedOpts::default(), 3, 9);
        assert_eq!(fr.is_ok(), closes, "raw Freivalds {:?}", fr.err());
        let run = box_from_exact(&data, &p, 6, None, Some(prune), 4, None).unwrap();
        let mut bd = run.data.clone();
        assert_eq!(closed_d2_freivalds(&mut bd, prune, None, 3, 2).is_ok(), closes);
        assert_eq!(close_box(&bd, prune, None, None, None).is_ok(), closes);
        if closes {
            assert_eq!(close_box(&bd, prune, None, None, None).unwrap(), close_exact_colored(&data, prune, None));
        }
    }
    eprintln!("(c) re-inserted final d∘d component of (3,3) {:?} (prune 4): Freivalds passes; with a term outside the ideal it fails (= closure_kills_map)", word);
}

#[test]
fn c_point_added_twice_overlap_header_downset() {
    let lam = [2u8, 2];
    let ring = ERing::new(&lam);
    let mut rng = StdRng::seed_from_u64(3);
    let grid = Arc::new(BoxGrid::new(BoxNodes::<F31>::new(&lam, 5, 3), None));
    let polys: Vec<EPoly<F31>> = (0..6).map(|w| rand_poly(&ring, w, 3, &mut rng)).collect();
    let weights: Vec<Option<u32>> = (0..6).map(Some).collect();
    let vals = |id: u32| -> Vec<F31> { let x = grid.coords(id); polys.iter().map(|p| eval_epoly(&ring, p, &x)).collect() };
    // a point added twice
    let mut acc = BoxAcc::new(grid.clone(), weights.clone(), 1, 2, None);
    acc.add_point(0, &vals(0)).unwrap();
    let e = acc.add_point(0, &vals(0)).unwrap_err();
    assert!(e.contains("already added"), "{}", e);
    // overlapping shard merge
    let mut s1 = BoxAcc::new(grid.clone(), weights.clone(), 1, 2, None);
    let mut s2 = BoxAcc::new(grid.clone(), weights.clone(), 1, 2, None);
    for id in 0..10 {
        s1.add_point(id, &vals(id)).unwrap();
    }
    for id in 9..20 {
        s2.add_point(id, &vals(id)).unwrap();
    }
    let e = s1.merge(&s2).unwrap_err();
    assert!(e.contains("overlapping"), "{}", e);
    // mismatched headers: library hash, word hash, prune_n, prime/nodes (another node seed)
    for other in [BoxAcc::new(grid.clone(), weights.clone(), 9, 2, None), BoxAcc::new(grid.clone(), weights.clone(), 1, 3, None), BoxAcc::new(grid.clone(), weights.clone(), 1, 2, Some(5))] {
        let e = s1.merge(&other).unwrap_err();
        assert!(e.contains("mismatched"), "{}", e);
    }
    let g2 = Arc::new(BoxGrid::new(BoxNodes::<F31>::new(&lam, 5, 4), None));
    assert!(s1.merge(&BoxAcc::new(g2.clone(), weights.clone(), 1, 2, None)).unwrap_err().contains("mismatched"));
    let g3 = Arc::new(BoxGrid::new(BoxNodes::<F31>::new(&lam, 5, 3), Some(0)));
    assert!(s1.merge(&BoxAcc::new(g3, weights.clone(), 1, 2, None)).unwrap_err().contains("mismatched"));
    // and a checkpoint read against another node sequence is refused
    let dir = std::env::temp_dir().join(format!("krhom_boxclose_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("s1.box");
    s1.save(&path).unwrap();
    assert!(BoxAcc::load(&path, g2).is_err());
    // truncated checkpoint
    let raw = std::fs::read(&path).unwrap();
    std::fs::write(dir.join("t.box"), &raw[..raw.len() - 8]).unwrap();
    assert!(BoxAcc::load(&dir.join("t.box"), grid.clone()).err().unwrap().contains("corrupt"));
    // incomplete down-set: all points but one of size <= needed
    let mut inc = BoxAcc::new(grid.clone(), weights.clone(), 1, 2, None);
    let need = inc.needed_size().unwrap();
    let pts = grid.points_up_to(need);
    for &id in &pts[..pts.len() - 1] {
        inc.add_point(id, &vals(id)).unwrap();
    }
    let e = inc.finish().err().unwrap();
    assert!(e.contains("incomplete down-set"), "{}", e);
    inc.add_point(*pts.last().unwrap(), &vals(*pts.last().unwrap())).unwrap();
    assert!(inc.finish().is_ok());
    // underestimated weight: the residual at w + 1 is nonzero
    let mut wlow = weights.clone();
    wlow[3] = Some(2);
    let e = extract(&grid, &polys, &wlow, None).err().unwrap();
    assert!(e.contains("residual"), "{}", e);
    // a nonzero value on an entry declared 0
    let mut wn = weights.clone();
    wn[2] = None;
    assert!(extract(&grid, &polys, &wn, None).is_err());
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn c_entry_declared_at_w_plus_1_fails_scaling() {
    let (sc, p, run) = small_run(None, None);
    let nodes = BoxNodes::<F31>::new(&sc.data.lam, 6, 21);
    let lambda = F31::from_i64(987_654_321);
    let sn = nodes.scaled(&lambda);
    let grid = BoxGrid::new(nodes.clone(), None);
    let mu = grid.mu(grid.npoints - 2);
    let vz = eval_entries(&sc.data, &p, &run.kinds, &point_coords(&nodes, &mu));
    let vlz = eval_entries(&sc.data, &p, &run.kinds, &point_coords(&sn, &mu));
    let rep = check_scaling_all(&vz, &vlz, &run.weights, &lambda).unwrap();
    assert!(rep.checked > 0);
    // declare one nonzero entry at w + 1 (an overestimated weight: extraction alone gives zeros
    // at w + 1 and the residual check cannot see it)
    let e = (0..run.weights.len()).find(|&e| run.weights[e].is_some() && !vz[e].is_zero()).unwrap();
    let mut w1 = run.weights.clone();
    w1[e] = w1[e].map(|w| w + 1);
    let err = check_scaling_all(&vz, &vlz, &w1, &lambda).unwrap_err();
    assert!(err.contains("scaling"), "{}", err);
    // the extraction with that declaration still finishes (silent zeros: no residual sees an
    // overestimated weight), which is why the scaling check is required
    assert!(box_from_exact(&sc.data, &p, 6, None, None, 21, Some(e)).is_ok());
}

#[test]
fn c_valid_n_and_d2_cert_guards() {
    // a prune-4 scan closed at N' = 5: Err (close_box and Freivalds)
    let (sc4, _, run4) = small_run(None, Some(4));
    let mut d4 = run4.data.clone();
    assert_eq!(d4.valid_n, 4);
    assert!(close_box(&d4, 5, None, None, None).unwrap_err().contains("valid_n"));
    assert!(closed_d2_freivalds(&mut d4, 5, None, 3, 1).is_err());
    // pruned data without a D2Cert: Err; with the certificate: the right answer
    let e = close_box(&d4, 4, None, None, None).unwrap_err();
    assert!(e.contains("D² certificate"), "{}", e);
    closed_d2_freivalds(&mut d4, 4, None, 3, 1).unwrap();
    assert_eq!(close_box(&d4, 4, None, None, None).unwrap(), close_exact_colored(&sc4.data, 4, None));
    // a certificate for another (N', mode) does not count
    assert!(close_box(&d4, 3, None, None, None).is_err());
    assert!(close_box(&d4, 4, Some(0), None, None).is_err());
    // the box N bounds valid_n too
    let (scu, pu, _) = small_run(None, None);
    let r = box_from_exact(&scu.data, &pu, 4, None, None, 3, None).unwrap();
    assert_eq!(r.data.valid_n, 4);
    assert!(close_box(&r.data, 5, None, None, None).is_err());
    // a T0 record of another library or N is refused
    let mut moy = BTreeMap::new();
    for (&x, _) in &r.data.base.types {
        moy.insert(x, LPoly::zero());
    }
    assert!(box_moy_check(&r.data, 4, None, &moy, None).is_err(), "T0 must fail against wrong MOY values");
    let fake = T0Record { n: 4, marked: None, library_hash: 1, ranks: BTreeMap::new() };
    assert!(close_box(&r.data, 4, None, Some(&fake), None).unwrap_err().contains("T0 record"));
}

// ------------------------------------------------------------------ (d) order and shards

#[test]
fn d_shuffled_order_and_shards_agree() {
    let mut cases = 0;
    for (lam, n, marked) in [(vec![3u8, 3], 6usize, None), (vec![2, 2, 2], 5, Some(0)), (vec![3, 3, 3], 5, None)] {
        let ring = ERing::new(&lam);
        let mut rng = StdRng::seed_from_u64(77 + n as u64);
        let grid = Arc::new(BoxGrid::new(BoxNodes::<F32b>::new(&lam, n, 5), marked));
        let polys: Vec<EPoly<F32b>> = (0..=grid.top).map(|w| rand_poly(&ring, w, 4, &mut rng)).collect();
        let weights: Vec<Option<u32>> = (0..=grid.top).map(Some).collect();
        let base = extract(&grid, &polys, &weights, None).unwrap();
        // shuffled order
        let mut pts = grid.points_up_to(grid.top);
        pts.shuffle(&mut rng);
        let sh = extract(&grid, &polys, &weights, Some(&pts)).unwrap();
        assert_eq!(base.coords, sh.coords);
        // three shards through checkpoints, merged
        let dir = std::env::temp_dir().join(format!("krhom_boxshard_{}_{}", std::process::id(), n));
        std::fs::create_dir_all(&dir).unwrap();
        let k = pts.len();
        let cuts = [0, k / 3, (2 * k) / 3, k];
        for s in 0..3 {
            let mut acc = BoxAcc::new(grid.clone(), weights.clone(), 5, 6, None);
            for &id in &pts[cuts[s]..cuts[s + 1]] {
                let x = grid.coords(id);
                acc.add_point(id, &polys.iter().map(|p| eval_epoly(&ring, p, &x)).collect::<Vec<_>>()).unwrap();
            }
            acc.save(&dir.join(format!("shard{}.box", s))).unwrap();
        }
        let mut m = BoxAcc::load(&dir.join("shard0.box"), grid.clone()).unwrap();
        assert!(m.finish().is_err(), "a single shard must not finish");
        for s in [2, 1] {
            let o = BoxAcc::load(&dir.join(format!("shard{}.box", s)), grid.clone()).unwrap();
            m.merge(&o).unwrap();
        }
        assert_eq!(m.points_done(), k);
        assert_eq!(m.finish().unwrap().coords, base.coords);
        // merging a shard twice is refused
        let o = BoxAcc::load(&dir.join("shard1.box"), grid.clone()).unwrap();
        assert!(m.merge(&o).is_err());
        std::fs::remove_dir_all(&dir).ok();
        cases += 1;
    }
    eprintln!("(d) {} cases: shuffled and sharded extraction identical", cases);
}

#[test]
fn moy_convention_matches_single_type_closure() {
    // the single-object closure of each final type of (2,2) σ^3 is its MOY value
    let sc = xscan::<F31>(&[2, 2], &gens(&[1, 1, 1]), &XOpts::default()).unwrap();
    for nn in 2..=5 {
        let moy = moy_of_scan(&sc, nn);
        for (&x, tb) in &sc.data.types {
            let d = ClosureDataC { lam: sc.data.lam.clone(), objs: vec![(x, 0, 0)], types: BTreeMap::from([(x, tb.clone())]), edges: vec![], crossings: vec![] };
            let pc: Poincare = close_exact_colored(&d, nn, None);
            assert_eq!(pc.euler(), moy[&x], "type {} N={}", x, nn);
        }
    }
}

#[test]
fn extraction_cost_per_entry() {
    // Mult-adds per entry of the streaming accumulators (outputs |ρ| = w and residuals
    // |ρ| = w + 1: Σ_μ #{ρ ⊇ μ} with nonzero (S^{-1})_{ρμ}) and memory per entry, by weight; and
    // the measured time per entry of a full extraction (one random entry per weight). Compare
    // R3 §3 (streaming: ≤ 4.2e3 for 2 blocks, ≤ 6.1e5 for 3 blocks; ≤ 384 B / 6.2 KB per entry).
    for (lam, marked) in [(vec![3u8, 3], None), (vec![3, 3], Some(0)), (vec![3, 3, 3], Some(0)), (vec![3, 3, 3], None)] {
        let grid = Arc::new(BoxGrid::new(BoxNodes::<F31>::new(&lam, 6, 3), marked));
        let mut sup = vec![];
        let mut best = (0u64, 0u32, 0usize);
        let mut tot_ops = 0u64;
        for w in 0..=grid.top {
            let mut ops = 0u64;
            for id in grid.points_up_to((w + 1).min(grid.top)) {
                for s in [w, w + 1] {
                    grid.supersets_of_size(id, s, &mut sup);
                    ops += sup.len() as u64;
                }
            }
            tot_ops += ops;
            let mem = 8 * (grid.count_of_size(w) + grid.count_of_size(w + 1));
            if ops > best.0 {
                best = (ops, w, mem);
            }
        }
        let ring = ERing::new(&lam);
        let mut rng = StdRng::seed_from_u64(1);
        let polys: Vec<EPoly<F31>> = (0..=grid.top).map(|w| rand_poly(&ring, w, 3, &mut rng)).collect();
        let weights: Vec<Option<u32>> = (0..=grid.top).map(Some).collect();
        let t = Instant::now();
        extract(&grid, &polys, &weights, None).unwrap();
        let secs = t.elapsed().as_secs_f64();
        let maxmem = (0..=grid.top).map(|w| 8 * (grid.count_of_size(w) + grid.count_of_size(w + 1))).max().unwrap();
        eprintln!(
            "cost {:?} N=6 marked {:?}: {} points; max mult-adds per entry {} (w = {}), mean over weights {:.0}; max accumulator memory per entry {} B; extraction of {} entries (one per weight, incl. evaluation) {:.3}s = {:.2} ms/entry",
            lam,
            marked,
            grid.npoints,
            best.0,
            best.1,
            tot_ops as f64 / (grid.top + 1) as f64,
            maxmem,
            polys.len(),
            secs,
            1000.0 * secs / polys.len() as f64
        );
    }
}
