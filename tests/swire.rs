//! M6a-1b: sengine wired to sindec and boxclose (src/colored/swire.rs). SPEC_3C_GL6 rev1 §2 "M6a-1b".
//!
//! * sengine against xengine at 2 points: final object multisets, K_0 = decategorified scan,
//!   d∘d = 0, identical structure at independent records and under replay; (2,2) and (3,2) by
//!   default, (3,3) `#[ignore]`.
//! * (1,1,1) final object counts against engine2.
//! * The one-word box loop against xengine's exact closure: (2,2) σ^3 (stages 2, 3) and (3,2) σ^2,
//!   N' = 3..6, unreduced and reduced; the scaling check (D14). Its negative test on the wired
//!   path is `scaling_check_catches_a_misdeclared_weight` in tests/boxrun.rs.
use krhom::colored::boxclose::*;
use krhom::input::crossing::Gen;
use krhom::algebra::field::{Field, F31, F32b};
use krhom::colored::sengine::*;
use krhom::library::sindec::*;
use krhom::colored::swire::*;
use krhom::colored::xengine::{xscan, XOpts};
use std::path::PathBuf;

fn costs() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/homprov/costs_cloud.txt")
}
fn wlib<F: Field>(lam: &[u8], prune: Option<u8>, words: &[Vec<i32>], closure: bool) -> WLib<F> {
    let (prov, stamp) = provider_for_stamped::<F>(lam, Some(&costs()), 1.0);
    let sd = SIndec::with_stamp(lam, prov, stamp, SIndecOpts { prune_n: prune, ..Default::default() });
    build_wlib(sd, words, closure).unwrap()
}
fn gens(w: &[i32]) -> Vec<Gen> {
    w.iter().map(|&g| Gen { i: (g.unsigned_abs() - 1) as usize, positive: g > 0 }).collect()
}
fn objs_of<F: Field>(lib: &WLib<F>, ps: &PointScan<F>) -> Vec<(String, i32, i32)> {
    let mut v: Vec<(String, i32, i32)> = ps.objs.iter().map(|(_, o)| (lib.names[&o.p].clone(), o.h, o.t)).collect();
    v.sort();
    v
}

/// Record at two independent generic points, replay the first trace at the second: identical
/// structure, final objects = xengine, K_0 and LB asserted every crossing, d∘d = 0 (unpruned).
fn against_xengine<F: Field>(lam: &[u8], word: &[i32], prune: Option<u8>) {
    let lib = wlib::<F>(lam, prune, &[word.to_vec()], false);
    let k: usize = lam.iter().map(|&a| a as usize).sum();
    let opts = ScanOpts { eng: EngOpts { check_d2: prune.is_none(), ..Default::default() }, close_stages: vec![], k0_check: true };
    let a = scan_point(&lib, word, record_x::<F>(k, 1), 0, Mode::Record, &opts).unwrap();
    let b = scan_point(&lib, word, record_x::<F>(k, 2), 1, Mode::Record, &opts).unwrap();
    a.ps.same_structure(&b.ps).unwrap();
    let c = scan_point(&lib, word, record_x::<F>(k, 2), 1, Mode::Replay(a.ps.frozen.clone()), &opts).unwrap();
    a.ps.same_structure(&c.ps).unwrap();
    let l: Vec<usize> = lam.iter().map(|&x| x as usize).collect();
    let xs = xscan::<F>(&l, &gens(word), &XOpts { prune_n: prune.map(|n| n as usize), ..Default::default() }).unwrap();
    assert_eq!(objs_of(&lib, &a.ps), xs.final_objects, "{:?} {:?} prune {:?}: sengine != xengine", lam, word, prune);
    eprintln!("{:?} {:?} prune {:?}: {} final objects = xengine; LB per crossing {:?}; inserted {:?}", lam, word, prune, xs.final_objects.len(), a.lb, a.ps.stages.iter().map(|s| s.inserted).collect::<Vec<_>>());
}

#[test]
fn sengine_against_xengine_22_32() {
    against_xengine::<F31>(&[2, 2], &[1, 1, 1], None);
    against_xengine::<F31>(&[2, 2], &[1, -1, 1, 1], None);
    against_xengine::<F31>(&[2, 2], &[1, 1, 1], Some(3));
    against_xengine::<F32b>(&[2, 2], &[-1, -1], None);
    against_xengine::<F31>(&[3, 2], &[1, 1], None);
    against_xengine::<F31>(&[3, 2], &[1, -1], None);
    against_xengine::<F31>(&[3, 2], &[1, 1, 1, 1], Some(4));
}

#[test]
#[ignore = "M6a-1b: (3,3) against xengine (cloud, minutes)"]
fn sengine_against_xengine_33() {
    against_xengine::<F31>(&[3, 3], &[1, 1], None);
    against_xengine::<F31>(&[3, 3], &[1, 1, 1], None);
    against_xengine::<F31>(&[3, 3], &[1, 1, 1], Some(4));
    against_xengine::<F31>(&[3, 3], &[-1, -1], Some(5));
}

#[test]
fn objects_111_against_engine2() {
    let words: Vec<Vec<i32>> = vec![vec![1, -2, 1, -2], vec![1, 1, 1, 2], vec![1, 2, -1, 2], vec![1, 2, 1, 2, 1, 2]];
    let lib = wlib::<F31>(&[1, 1, 1], None, &words, false);
    let opts = ScanOpts { eng: EngOpts { check_d2: true, ..Default::default() }, close_stages: vec![], k0_check: true };
    for w in &words {
        let a = scan_point(&lib, w, record_x::<F31>(3, 5), 0, Mode::Record, &opts).unwrap();
        let e2 = krhom::uncolored::engine2::khr_engine2::<F31>(3, &gens(w), &[2], &[], 7);
        assert_eq!(a.ps.objs.len(), e2.final_objects, "(1,1,1) {:?}: sengine {} objects, engine2 {}", w, a.ps.objs.len(), e2.final_objects);
        eprintln!("(1,1,1) {:?}: {} final objects = engine2", w, e2.final_objects);
    }
}

// ------------------------------------------------------------------ box loop against the exact closure

fn box_against_exact<F: PrimeField>(lam: &[u8], word: &[i32], stages: &[usize], seed: u64) {
    let lib = wlib::<F>(lam, None, &[word.to_vec()], true);
    let n_box = 6;
    let o = BoxOpts::new(n_box, seed, vec![3, 4, 5, 6], vec![None, Some(0)], stages.to_vec(), EngOpts { check_d2: true, ..Default::default() });
    let out = box_word(&lib, word, &o).unwrap();
    // the D14 scaling check, χ = decat_chi and the determinant pass ran inside box_word
    let sc = out.scaling.as_ref().expect("box_word runs the scaling check by default");
    assert!(sc.report.checked > 0);
    assert_eq!(out.chi_checked, out.stages.iter().map(|s| s.closed.len()).sum::<usize>());
    assert!(out.det_pass.1 > 0);
    let l: Vec<usize> = lam.iter().map(|&x| x as usize).collect();
    let xs = xscan::<F>(&l, &gens(word), &XOpts { close_stages: true, ..Default::default() }).unwrap();
    for sr in &out.stages {
        for (&(nn, marked), p) in &sr.closed {
            let want = xs.close_stage_marked(sr.stage, &[nn], marked).unwrap()[&nn].clone();
            assert_eq!(p, &want, "{:?} {:?} stage {} N' {} marked {:?}", lam, word, sr.stage, nn, marked);
        }
        eprintln!("{:?} {:?} stage {}: {} closures = xengine ({} entries, {} objects)", lam, word, sr.stage, sr.closed.len(), sr.entries, sr.objects);
    }
    let lambda = F::from_u64(0x1234_5677);
    let ug = BoxGrid::new(BoxNodes::<F>::new(lam, n_box, seed), None);
    let (checked, zero) = scaling_check(&lib, word, &o, &lambda, &[ug.npoints / 3, ug.npoints - 1]).unwrap();
    eprintln!("{:?} {:?}: {} points, scan {:.2}s; scaling {} checked, {} zero at both", lam, word, out.points, out.scan_secs, checked, zero);
    assert!(checked > 0);
}

#[test]
fn box_loop_against_exact_22_32() {
    box_against_exact::<F31>(&[2, 2], &[1, 1, 1], &[2, 3], 21);
    box_against_exact::<F32b>(&[2, 2], &[1, 1, 1], &[2, 3], 22);
    box_against_exact::<F31>(&[3, 2], &[1, 1], &[2], 23);
}

/// 4_1 = σ1 σ2^-1 σ1 σ2^-1 on (3,3,3): types of the LB complex (support of the decategorified
/// K_0 class) per prune level, final and over all stages (C_feas m9: 34 vs 38).
#[test]
fn type_counts_41_333() {
    let w = [1, -2, 1, -2];
    for prune in [None, Some(6u8), Some(5), Some(4)] {
        let mut d = Decat::new(&[3, 3, 3]);
        let mut labels = vec![3u8, 3, 3];
        let mut all = std::collections::BTreeSet::new();
        let mut fin = std::collections::BTreeSet::new();
        let mut lb = 0;
        for &g in &w {
            let (i, pos) = letter(g);
            d.step(labels[i], labels[i + 1], pos, i);
            labels.swap(i, i + 1);
            let k = d.k0(prune);
            fin = k.keys().map(|(n, _)| n.clone()).collect();
            all.extend(fin.iter().cloned());
            lb = k.values().map(|c| c.unsigned_abs()).sum::<u64>();
        }
        eprintln!("4_1 (3,3,3) prune {:?}: LB {} objects, {} types in the final LB complex, {} types over all stages", prune, lb, fin.len(), all.len());
    }
}
