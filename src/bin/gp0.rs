//! M4c B2 census G-P0 (M4C_PLAN §4 B2): per remaining (3,3,3) pair target, the Lemma A plan
//! (junction words, step unknowns), the exit word W_q, and the rex chain W_q → L_q
//! (src/rexplan.rs), with the budget verdict. Pure combinatorics.
//!
//! usage: gp0 [--lam 3,3,3] --n 5|6 [--pairs 313,322,...] [--out FILE.jsonl] [--costs PATH]
//!            [--visit N] [--maxlen N] [--exact-budget N] [--no-pad] [--alt-len N]
use krhom::dcoset::{junction_types, reduced_for, StepCache};
use krhom::homprov::CostTable;
use krhom::leaves::{step_score, Leaves};
use krhom::rexplan::*;
use krhom::schur::{flows, Schur, Table, Tid};
use krhom::sindec::{canonical_ladder, pruned_plan};
use krhom::smapio::rungs_str;
use krhom::web::Rung;
use std::io::Write;
use std::time::Instant;

fn opt(args: &[String], key: &str) -> Option<String> {
    args.iter().position(|a| a == key).and_then(|i| args.get(i + 1).cloned())
}

/// N = 5 and N = 6 lists after M5c (M5c_report "Remaining failure list by cause"; 24 is F32b only)
const N5: &[usize] = &[313, 322, 326, 410];
const N6: &[usize] = &[127, 212, 227, 45, 145, 458, 319, 335, 411, 383, 278, 293, 273, 78, 313, 326, 285, 385, 292, 37, 137, 332, 410, 82, 283, 418, 336, 381, 376, 281, 159, 323, 322, 71, 334, 419, 219, 467, 24];

/// predicted seconds of an exact call: the exact-shape record, else a crude model [E] from the
/// calibration (10⁴-10⁵ unknowns ≤ 48 s, 10⁵-10⁶ uncensored ≤ 419 s)
fn predict_exact(costs: &Option<CostTable>, lam: &[u8], a: &[Rung], b: &[Rung], u: u128) -> (f64, bool) {
    if u == 0 {
        return (0.0, true);
    }
    if let Some(c) = costs {
        let (fs, ft) = (flows(lam, a), flows(lam, b));
        if let Some(t) = c.predict_exact(lam, fs, ft, u, 0) {
            return (t, true);
        }
    }
    (50.0 * (u as f64 / 1e5).powf(1.3), false)
}

fn pad_table(t: &Table) -> Table {
    let (n, m) = (t.nr as usize, t.nc as usize);
    let mut e = vec![0u8; (n + 1) * (m + 1)];
    for i in 0..n {
        for j in 0..m {
            e[i * (m + 1) + j] = t.get(i, j);
        }
    }
    Table { nr: (n + 1) as u8, nc: (m + 1) as u8, e }
}

fn jstr(s: &str) -> String {
    format!("\"{}\"", s.replace('\\', "\\\\").replace('"', "\\\""))
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let lam: Vec<u8> = opt(&args, "--lam").unwrap_or_else(|| "3,3,3".into()).split(',').map(|x| x.parse().unwrap()).collect();
    let n: u8 = opt(&args, "--n").expect("--n").parse().unwrap();
    let pairs: Vec<usize> = match opt(&args, "--pairs") {
        Some(p) => p.split(',').map(|x| x.parse().unwrap()).collect(),
        None => if n == 5 { N5.to_vec() } else { N6.to_vec() },
    };
    let costs = opt(&args, "--costs").map(|p| CostTable::load(std::path::Path::new(&p)).expect("cost table"));
    let visit: usize = opt(&args, "--visit").map(|v| v.parse().unwrap()).unwrap_or(60_000);
    let maxlen: usize = opt(&args, "--maxlen").map(|v| v.parse().unwrap()).unwrap_or(12);
    let exact_budget: u128 = opt(&args, "--exact-budget").map(|v| v.parse().unwrap()).unwrap_or(100_000);
    let no_pad = args.iter().any(|a| a == "--no-pad");
    let alt_len: usize = opt(&args, "--alt-len").map(|v| v.parse().unwrap()).unwrap_or(5);
    let mut out = opt(&args, "--out").map(|p| std::fs::File::create(p).expect("out"));
    let mut s = Schur::new(&lam);
    let plan = pruned_plan(&mut s, Some(n));
    let id = s.intern(Table::identity(&lam));
    let mut leaves = Leaves::default();
    leaves.alt_len = alt_len;
    let mut pl = Planner::new();
    let mut lam4 = lam.clone();
    lam4.push(0);
    let mut s4 = Schur::new(&lam4);
    let mut pl4 = Planner::new();
    println!("G-P0 census lam {:?} N = {}: {} pairs", lam, n, pairs.len());
    for &k in &pairs {
        let pp = plan.pairs[k].clone();
        let mut x = canonical_ladder(&mut s, pp.p);
        x.push(pp.rung);
        let xf = flows(&lam, &x);
        for tg in &pp.kept {
            let t0 = Instant::now();
            let q = tg.q;
            let lq = canonical_ladder(&mut s, q);
            let rex_class = reduced_for(&mut s, &x, q).ok();
            // Lemma A plan to q
            let mut sc = StepCache::new();
            let js = junction_types(&mut s, &mut sc, &x, &[q]);
            let m = x.len();
            let wp = leaves.plan_words(&mut s, &mut sc, &x, &js, id);
            let mut max_step = 0u128;
            let mut steps = 0usize;
            let mut nonid = 0usize;
            let mut exact_secs = 0.0f64;
            let mut unrecorded = 0usize;
            let mut max_junction_flows = xf;
            let mut worst = String::new();
            for lv in 1..=m {
                for (t, wt) in &wp[lv] {
                    max_junction_flows = max_junction_flows.max(flows(&lam, wt));
                    for (tp, wtp) in &wp[lv - 1] {
                        let Some(g) = sc.step(&mut s, *tp, &x[lv - 1]).into_iter().find(|e| e.0 == *t).map(|e| e.1) else { continue };
                        steps += 1;
                        let u = step_score(&lam, wtp, &x[lv - 1], wt, &g);
                        if u > 0 {
                            nonid += 1;
                            let mut a = wtp.clone();
                            a.push(x[lv - 1]);
                            let (secs, rec) = predict_exact(&costs, &lam, &a, wt, u);
                            exact_secs += secs * (g.c.iter().filter(|&&c| c != 0).count() as f64).max(1.0);
                            if !rec {
                                unrecorded += 1;
                            }
                        }
                        if u > max_step {
                            max_step = u;
                            worst = format!("level {} {} -> {} [{}]·{} -> [{}]", lv, s.table(*tp).fmt(), s.table(*t).fmt(), rungs_str(wtp), rungs_str(&[x[lv - 1]]), rungs_str(wt));
                        }
                    }
                }
            }
            let wq: Vec<Rung> = wp[m].iter().find(|e| e.0 == q).map(|e| e.1.clone()).unwrap_or_default();
            let wq_reduced = reduced_for(&mut s, &wq, q).ok();
            // option A: final step directly into L_q
            let mut final_lq = 0u128;
            if m > 0 {
                for (tp, wtp) in &wp[m - 1] {
                    if let Some(g) = sc.step(&mut s, *tp, &x[m - 1]).into_iter().find(|e| e.0 == q).map(|e| e.1) {
                        final_lq = final_lq.max(step_score(&lam, wtp, &x[m - 1], &lq, &g));
                    }
                }
            }
            // option B: rex chain W_q -> L_q (closed forms first, then with Exact moves, then padded)
            let base = PlanOpts { max_len: maxlen.max(wq.len().max(lq.len()) + 2), max_visit: visit, ..PlanOpts::default() };
            let mut tried = vec![];
            let mut chain: Option<(Chain, &str)> = None;
            if wq_reduced {
                for (label, ex) in [("closed", 0u128), ("closed+exact", exact_budget)] {
                    let o = PlanOpts { exact_budget: ex, ..base.clone() };
                    match pl.plan(&mut s, q, &wq, &lq, &o) {
                        Ok(c) => {
                            chain = Some((c, label));
                            break;
                        }
                        Err(e) => tried.push(format!("{}: {}", label, e)),
                    }
                }
                if chain.is_none() && !no_pad {
                    let q4 = s4.intern(pad_table(s.table(q)));
                    let o = PlanOpts { exact_budget, ..base.clone() };
                    match pl4.plan(&mut s4, q4, &wq, &lq, &o) {
                        Ok(c) => chain = Some((c, "pad4")),
                        Err(e) => tried.push(format!("pad4: {}", e)),
                    }
                }
            } else {
                tried.push("W_q not reduced for q".into());
            }
            let secs = t0.elapsed().as_secs_f64();
            // verdict
            let leaves_ok = max_step <= 100_000 || (max_step <= 400_000 && unrecorded == 0);
            let (chain_ok, chain_desc, pad, clen, tflows, aflows, thm, unc, cexact, cwords, csecs, moves) = match &chain {
                Some((c, label)) => {
                    let lam_c = if *label == "pad4" { lam4.clone() } else { lam.clone() };
                    let unc = c.moves.iter().filter(|m| m.uncertifiable()).count();
                    let ok = c.max_template_flows() <= krhom::leafgen::MAX_INTERP_FLOWS as u128 && unc == 0 && c.max_exact_unknowns() <= exact_budget && c.max_word_flows(&lam_c) <= 1_000_000;
                    let mv: Vec<String> = c.moves.iter().map(|m| m.fmt()).collect();
                    (ok, label.to_string(), *label == "pad4", c.len(), c.max_template_flows(), c.max_assoc_flows(), c.theorem_moves(), unc, c.max_exact_unknowns(), c.max_word_flows(&lam_c), c.predicted_secs(), mv)
                }
                None => (false, "none".into(), false, 0, 0, 0, 0, 0, 0, 0, 0.0, vec![]),
            };
            let words_ok = max_junction_flows.max(cwords) <= 1_000_000;
            let feasible = leaves_ok && chain_ok && words_ok;
            let reason = if feasible {
                "ok".to_string()
            } else {
                let mut r = vec![];
                if !leaves_ok {
                    r.push(format!("Lemma A step {} unknowns{}", max_step, if unrecorded > 0 { " (no record)" } else { "" }));
                }
                if chain.is_none() {
                    r.push("no rex chain found".into());
                } else if !chain_ok {
                    r.push(format!("chain: template {} flows / uncertifiable assoc {} / exact {}", tflows, unc, cexact));
                }
                if !words_ok {
                    r.push(format!("word flows {}", max_junction_flows.max(cwords)));
                }
                r.join("; ")
            };
            println!(
                "pair {} q {} s {} mult {} X {} flows{}: steps {} (non-Id {}), max step {} unknowns, W_q [{}] {} flows, L_q [{}]; final-into-L_q {} unknowns; chain {} len {} (template ≤ {} flows, assoc ≤ {} flows, theorem {}, exact ≤ {}, words ≤ {}), pad {}; predicted {:.0} s exact + {:.0} s templates; {} ({:.1} s)",
                k,
                s.table(q).fmt(),
                tg.s,
                tg.mult,
                xf,
                if rex_class { " [rex class]" } else { "" },
                steps,
                nonid,
                max_step,
                rungs_str(&wq),
                flows(&lam, &wq),
                rungs_str(&lq),
                final_lq,
                chain_desc,
                clen,
                tflows,
                aflows,
                thm,
                cexact,
                cwords,
                pad,
                exact_secs,
                csecs,
                if feasible { "FEASIBLE".to_string() } else { format!("INFEASIBLE ({})", reason) },
                secs
            );
            if !worst.is_empty() {
                println!("    worst step: {}", worst);
            }
            for mv in &moves {
                println!("    {}", mv);
            }
            for t in &tried {
                println!("    tried {}", t);
            }
            if let Some(f) = out.as_mut() {
                writeln!(
                    f,
                    "{{\"n\":{},\"pair\":{},\"q\":{},\"s\":{},\"mult\":{},\"x_flows\":{},\"rex_class\":{},\"steps\":{},\"nonid_steps\":{},\"max_step_unknowns\":{},\"unrecorded_steps\":{},\"worst_step\":{},\"wq\":{},\"wq_flows\":{},\"lq\":{},\"final_into_lq_unknowns\":{},\"chain\":{},\"chain_len\":{},\"pad\":{},\"max_template_flows\":{},\"max_assoc_flows\":{},\"theorem_moves\":{},\"uncertifiable_assoc\":{},\"max_exact_unknowns\":{},\"max_word_flows\":{},\"pred_exact_s\":{:.1},\"pred_template_s\":{:.1},\"feasible\":{},\"reason\":{},\"moves\":[{}],\"secs\":{:.2}}}",
                    n,
                    k,
                    jstr(&s.table(q).fmt()),
                    tg.s,
                    tg.mult,
                    xf,
                    rex_class,
                    steps,
                    nonid,
                    max_step,
                    unrecorded,
                    jstr(&worst),
                    jstr(&rungs_str(&wq)),
                    flows(&lam, &wq),
                    jstr(&rungs_str(&lq)),
                    final_lq,
                    jstr(&chain_desc),
                    clen,
                    pad,
                    tflows,
                    aflows,
                    thm,
                    unc,
                    cexact,
                    max_junction_flows.max(cwords),
                    exact_secs,
                    csecs,
                    feasible,
                    jstr(&reason),
                    moves.iter().map(|m| jstr(m)).collect::<Vec<_>>().join(","),
                    secs
                )
                .unwrap();
                f.flush().unwrap();
            }
        }
    }
    eprintln!("planner stats {:?} / pad {:?}", pl.stats, pl4.stats);
}
