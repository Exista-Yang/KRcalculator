//! C2-4 (M5i): build the idempotent e_p of decomposable canonical ladders with `SIndec::ty`, and
//! check it against the exact End^0(L_p) where affordable: when char L_p = b_p + Σ b_r (distinct
//! r, all in degree 0), End^0(L_p) ≅ k^{1+#r} (Hom^0 between distinct indecomposables is 0), so
//! dim End^0 = 1 + #r and the idempotent with top scalar 1 and the ranks of B_p is unique: e_p
//! must lie in the span of the End^0 basis (stacked blocks at 2 points, all cosets). Optionally
//! decomposes L_p ⊗ R for every rung R (dec_rung, source idempotent) and builds the closure
//! generators through e_p for types of (λ, λ).
//! usage: m5i_check LAM N|none [--field F31|F32b] [--costs T] [--cache DIR] [--types T1;T2 | --all] [--end0 MAXUNK] [--rungs] [--closure]
use krhom::algebra::field::{Field, F31, F32b};
use krhom::library::homprov::{hom_unknowns, Budget, Exact, HomProvider};
use krhom::combinat::schur::{flows, Schur, Tid};
use krhom::library::sindec::{decomposable_ladders, provider_for, SIndec, SIndecOpts};
use krhom::colored::smap::{points, Eval};
use krhom::combinat::web::{Dir, Rung};
use std::time::Instant;

fn arg(a: &[String], k: &str) -> Option<String> {
    a.iter().position(|x| x == k).and_then(|i| a.get(i + 1).cloned())
}

fn run<F: Field>(a: &[String]) {
    let lam: Vec<u8> = a[1].split(',').map(|x| x.parse().unwrap()).collect();
    let n: Option<u8> = a[2].parse().ok();
    let costs = arg(a, "--costs").map(std::path::PathBuf::from);
    let cache = arg(a, "--cache").map(std::path::PathBuf::from);
    let end0max: Option<u64> = arg(a, "--end0").map(|x| x.parse().unwrap());
    let opts = SIndecOpts { cache_dir: cache, prune_n: n, layer_a_route_stamp: true, ..SIndecOpts::default() };
    let prov = provider_for::<F>(&lam, costs.as_deref(), 1.0);
    let mut sd: SIndec<F> = SIndec::new(&lam, prov, opts);
    let cases = decomposable_ladders(&mut sd.schur);
    let want: Option<Vec<String>> = arg(a, "--types").map(|x| x.split(';').map(|t| t.to_string()).collect());
    let all = a.iter().any(|x| x == "--all");
    let mut sel: Vec<(Tid, String, Vec<u8>)> = vec![];
    for c in &cases {
        if n.is_some_and(|nn| c.kill_level > nn as usize) && !all {
            continue;
        }
        if let Some(w) = &want {
            if !w.contains(&c.table) {
                continue;
            }
        }
        sel.push((c.q, c.table.clone(), c.mu.clone()));
    }
    if let Ok(force) = std::env::var("KRHOM_FORCE_LADDER") {
        // forced (non-canonical) ladders: the types named there
        for item in force.split(';') {
            let Some((t, _)) = item.split_once('=') else { continue };
            for mu in krhom::library::sindec::visited_blocks(&lam) {
                for q in sd.schur.canonical_basis(&mu).order.clone() {
                    if sd.schur.table(q).fmt() == t.trim() && !sel.iter().any(|x| x.0 == q) {
                        sel.push((q, t.trim().to_string(), mu.clone()));
                    }
                }
            }
        }
    }
    println!("{} types selected ({} decomposable canonical ladders in the visited blocks)", sel.len(), cases.len());
    let (mut ok, mut bad) = (0, 0);
    for (q, name, mu) in sel {
        let t0 = Instant::now();
        let r = sd.ty(q).map(|t| (t.lad, t.e, t.rank));
        let (lad, e, rank) = match r {
            Ok(x) => x,
            Err(err) => {
                println!("FAIL {} block {:?}: {}", name, mu, String::from(err));
                bad += 1;
                continue;
            }
        };
        let w = sd.ar.expr_of(lad).rungs.clone();
        let log = sd.idem_log.iter().find(|x| x.0 == name).cloned().unwrap_or_default();
        println!("IDEM {} block {:?} flows {} rank B_p {} e_p {} in {:.1}s: extra {}", name, mu, flows(&lam, &w), rank, e.is_some(), t0.elapsed().as_secs_f64(), log.1);
        let Some(e) = e else { continue };
        if let Some(mx) = end0max {
            let unk = hom_unknowns(&lam, &w, &w, 0);
            if unk <= mx as u128 {
                let t1 = Instant::now();
                let mut ex = Exact { max_label: 9, max_fibre_product: u64::MAX, max_unknowns: u64::MAX, max_call_s: 1e9, ..Exact::default() };
                match HomProvider::<F>::end0(&mut ex, &mut sd.ar, &mut sd.schur, lad, &Budget::default()) {
                    Ok((basis, _)) => {
                        let nextra = { let ch = sd.schur.ladder_char(&w); sd.schur.expand(&ch).len() - 1 };
                        let cos = sd.ar.flows(lad).cosets.clone();
                        let pts = points::<F>(sd.ar.k(), 0xe0e0, 2, 0);
                        let st = Eval::stacked(&sd.ar, &basis, &pts, &cos);
                        let mut st2 = st.clone();
                        st2.extend(Eval::stacked(&sd.ar, &[e], &pts, &cos));
                        let len = st[0].len();
                        let (r1, r2) = (krhom::algebra::dense::rank(&st, len), krhom::algebra::dense::rank(&st2, len));
                        let pass = basis.len() == 1 + nextra && r1 == basis.len() && r2 == r1;
                        println!("END0 {}: unknowns {} dim End^0 {} (expected {}), e_p in span: {} -> {} ({:.1}s)", name, unk, basis.len(), 1 + nextra, r2 == r1, if pass { "PASS" } else { "MISMATCH" }, t1.elapsed().as_secs_f64());
                        if !pass {
                            bad += 1;
                        }
                    }
                    Err(m) => println!("END0 {}: not computed ({:?})", name, m.first().map(|x| x.what.clone())),
                }
            } else {
                println!("END0 {}: skipped ({} unknowns > {})", name, unk, mx);
            }
        }
        if a.iter().any(|x| x == "--rungs") {
            let nn = lam.len();
            for i in 0..nn - 1 {
                for dir in [Dir::E, Dir::F] {
                    let giver = if dir == Dir::E { mu[i + 1] } else { mu[i] } as usize;
                    for r in 1..=giver {
                        let rg = Rung { i, dir, r };
                        let t1 = Instant::now();
                        match sd.dec_rung(q, rg) {
                            Ok(d) => println!("  DEC {} ⊗ {}: {} summands {:?} killed {} ({:.1}s)", name, krhom::library::smapio::rung_str(&rg), d.summands.len(), d.summands.iter().map(|s| s.s).collect::<Vec<_>>(), d.killed.len(), t1.elapsed().as_secs_f64()),
                            Err(err) => {
                                println!("  DECFAIL {} ⊗ {}: {}", name, krhom::library::smapio::rung_str(&rg), String::from(err));
                                bad += 1;
                            }
                        }
                    }
                }
            }
        }
        if a.iter().any(|x| x == "--closure") && mu == lam {
            let t1 = Instant::now();
            match sd.closure_gens(q) {
                Ok(g) => println!("  CLOSURE {}: {} ins {} outs ({:.1}s)", name, g.ins.len(), g.outs.len(), t1.elapsed().as_secs_f64()),
                Err(err) => {
                    println!("  CLOSUREFAIL {}: {}", name, String::from(err));
                    bad += 1;
                }
            }
        }
        ok += 1;
    }
    let st = sd.stats();
    println!("SUMMARY built {} failed {}; idempotents {} provider calls {} ({:.1}s) layer A loaded {} saved {} check {:.1}s", ok, bad, st.idempotents, st.provider_calls, st.provider_secs, st.layer_a_loaded, st.layer_a_saved, st.check_secs);
    if bad > 0 {
        std::process::exit(1);
    }
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    match arg(&a, "--field").as_deref() {
        Some("F32b") => run::<F32b>(&a),
        _ => run::<F31>(&a),
    }
}
