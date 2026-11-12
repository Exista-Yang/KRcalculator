//! m4c: M4c milestone A1 tools (combinatorics census G-L0 and the Stage 0 product closure run).
//!
//! m4c census [--lam 3,3,3] [--budget 100000] [--out FILE]
//!     G-L0: for every closure job of the M5b failure list (N = 6: 12 types; N = 5:
//!     1,1,1/1,1,1/1,1,1), the Lemma A plan on 3 uprights (D1 §1.3, M4C_PLAN §3.2): junction types
//!     per level on the identity paths, their junction words (prefix extension, C_m4c m3),
//!     word flows, and the exact step calls Hom^d(W_{t'}·ρ, W_t) for d in supp G_{t'->t} with
//!     their unknown counts (`homprov::hom_unknowns`, the replica of hu.py), plus the M5b
//!     whole-ladder calls for comparison. Also checks the four hu.py numbers of C_m4c.
//! m4c stage0 --field 31|32b --cache DIR [--costs PATH] [--replay] [--types all|product] [--product-first]
//!     Stage 0 (M4C_PLAN §3.2): closure generators of the product types through
//!     `leaves::ProductClosure` around the production chain, in a SIndec with layer A in DIR;
//!     `--replay` uses a provider that refuses every call (everything must load from layer A).
use krhom::dcoset::*;
use krhom::field::{Field, F31, F32b};
use krhom::homprov::*;
use krhom::schur::{flows, Lp, Schur, Table, Tid};
use krhom::sindec::*;
use krhom::smap::{Arena, ExprId, MapId};
use krhom::web::{Dir, Rung};
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::time::Instant;

fn opt(args: &[String], key: &str) -> Option<String> {
    args.iter().position(|a| a == key).and_then(|i| args.get(i + 1).cloned())
}
fn flag(args: &[String], key: &str) -> bool {
    args.iter().any(|a| a == key)
}
fn parse_lam(s: &str) -> Vec<u8> {
    s.split(',').map(|x| x.trim().parse().unwrap()).collect()
}
fn wstr(w: &[Rung]) -> String {
    if w.is_empty() {
        return "∅".into();
    }
    w.iter().map(|r| format!("{}{}.{}", if r.dir == Dir::E { "E" } else { "F" }, r.i, r.r)).collect::<Vec<_>>().join(" ")
}
fn lpstr(p: &Lp) -> String {
    let mut v = vec![];
    for (i, &c) in p.c.iter().enumerate() {
        if c != 0 {
            let e = p.lo + i as i32;
            v.push(if c == 1 { format!("v^{}", e) } else { format!("{}v^{}", c, e) });
        }
    }
    if v.is_empty() { "0".into() } else { v.join("+") }
}
fn sep(n: u128) -> String {
    let s = n.to_string();
    let mut out = String::new();
    for (i, ch) in s.chars().enumerate() {
        if i > 0 && (s.len() - i) % 3 == 0 {
            out.push(' ');
        }
        out.push(ch);
    }
    out
}

/// The M5b closure failure list (M5b_report; F31 and F32b identical): (N, type)
pub const CLOSURE_FAILS: [(u8, &str); 13] = [
    (5, "1,1,1/1,1,1/1,1,1"),
    (6, "0,0,3/3,0,0/0,3,0"),
    (6, "0,3,0/0,0,3/3,0,0"),
    (6, "1,0,2/2,0,1/0,3,0"),
    (6, "1,2,0/0,0,3/2,1,0"),
    (6, "0,3,0/1,0,2/2,0,1"),
    (6, "1,0,2/1,2,0/1,1,1"),
    (6, "1,1,1/0,2,1/2,0,1"),
    (6, "0,2,1/2,0,1/1,1,1"),
    (6, "1,1,1/1,0,2/1,2,0"),
    (6, "0,1,2/2,0,1/1,2,0"),
    (6, "0,2,1/1,0,2/2,1,0"),
    (6, "1,1,1/1,1,1/1,1,1"),
];

fn unk2(lam: &[u8], a: &[Rung], b: &[Rung], d: i32) -> u128 {
    hom_unknowns(lam, a, b, d).max(hom_unknowns(lam, b, a, d))
}

#[derive(Default)]
struct JobSummary {
    max_junctions: usize,
    max_word_flows: u128,
    max_word_len: usize,
    steps: usize,
    id_steps: usize,
    exact_steps: usize,
    max_unknowns: u128,
    over_1e4: usize,
    over_budget: usize,
    max_window_unknowns: u128,
    max_k_unknowns: u128,
    k_pad: usize,
    not_reduced: Vec<String>,
    kinds: BTreeMap<String, usize>,
}

/// Unknowns of the step t' -> t of ρ as a window on the K uprights only (D1 §1.7, SLL Thm 2.6 +
/// Lemma 3.9): the local tables m', m'' restricted to the K cells of Q and the two columns of ρ,
/// local words = canonical ladders of the local tables on the local bottom (the K cells). None if
/// the K column of Q has 3 nonzero cells (the local bottom needs a third, padded upright).
fn kwindow_unknowns(s: &mut Schur, tp: Tid, t: Tid, rho: &Rung, degs: &[i32]) -> Option<u128> {
    let (a, b) = (s.table(tp).clone(), s.table(t).clone());
    let g = grass(&a, &b, rho);
    let i = rho.i;
    // rows of L that come from column i of Q: L lists the cells of Q column by column
    let mut start = 0usize;
    for j in 0..i {
        start += (0..g.q.nr as usize).filter(|&r| g.q.get(r, j) > 0).count();
    }
    let nloc = g.local_lam.len();
    if nloc > 2 {
        return None;
    }
    if nloc < 2 {
        return Some(0);
    }
    let restrict = |m: &Table| -> Table {
        let mut e = vec![];
        for r in start..start + nloc {
            e.push(m.get(r, i));
            e.push(m.get(r, i + 1));
        }
        Table { nr: 2, nc: 2, e }
    };
    let (ma, mb) = (restrict(&g.m_from), restrict(&g.m_to));
    let lam = g.local_lam.clone();
    let mut ls = Schur::new(&lam);
    let (ia, ib) = (ls.intern(ma), ls.intern(mb));
    let mut src = canonical_ladder(&mut ls, ia);
    src.push(Rung { i: 0, dir: rho.dir, r: rho.r });
    let tgt = canonical_ladder(&mut ls, ib);
    Some(degs.iter().map(|&d| unk2(&lam, &src, &tgt, d)).max().unwrap_or(0))
}

fn census_job(s: &mut Schur, sc: &mut StepCache, p: Tid, budget: u128, out: &mut Vec<String>) -> JobSummary {
    let lam = s.lam.clone();
    let id = s.intern(Table::identity(&lam));
    let w = canonical_ladder(s, p);
    let js = junction_types(s, sc, &w, &[id]);
    let mut sum = JobSummary::default();
    // junction words, level by level
    let mut words: Vec<Vec<(Tid, Vec<Rung>, WordKind)>> = vec![vec![(id, vec![], WordKind::Canonical)]];
    for k in 1..=w.len() {
        let prev: Vec<(Tid, Vec<Rung>)> = words[k - 1].iter().map(|x| (x.0, x.1.clone())).collect();
        let mut lv = vec![];
        for j in &js[k] {
            let (ww, kind) = if j.t == id && k == w.len() { (vec![], WordKind::Canonical) } else { junction_word(s, j.t, &prev, Some(&w[k - 1]), 6, 20_000) };
            if kind == WordKind::NotReduced {
                sum.not_reduced.push(s.table(j.t).fmt());
            }
            lv.push((j.t, ww, kind));
        }
        words.push(lv);
    }
    out.push(format!("  L_p = {} ({} flows), junction types per level: {:?}", wstr(&w), flows(&lam, &w), js.iter().map(|l| l.len()).collect::<Vec<_>>()));
    for k in 1..=w.len() {
        sum.max_junctions = sum.max_junctions.max(js[k].len());
        let rho = w[k - 1];
        for (jt, (t, wt, kind)) in js[k].iter().zip(words[k].iter()) {
            let f = flows(&lam, wt);
            sum.max_word_flows = sum.max_word_flows.max(f);
            sum.max_word_len = sum.max_word_len.max(wt.len());
            let ks = match kind {
                WordKind::Prefix(_) => "prefix",
                WordKind::Canonical => "canonical",
                WordKind::Bfs => "bfs",
                WordKind::NotReduced => "NOT-REDUCED",
            };
            *sum.kinds.entry(ks.to_string()).or_default() += 1;
            out.push(format!("  level {} ({}): t = {} W_t = {} [{}; {} flows; window deg {}..{}]", k, wstr(&[rho]), s.table(*t).fmt(), wstr(wt), ks, f, jt.dmin, jt.dmax));
            for (tp, wp, _) in &words[k - 1] {
                let g = sc.step(s, *tp, &rho).into_iter().find(|x| x.0 == *t).map(|x| x.1);
                let Some(g) = g else { continue };
                // only predecessors that are junctions (on identity paths) matter
                if !js[k - 1].iter().any(|j| j.t == *tp) {
                    continue;
                }
                sum.steps += 1;
                let mut src = wp.clone();
                src.push(rho);
                if src == *wt {
                    sum.id_steps += 1;
                    out.push(format!("    step {} -> {}: G = {}: Id", s.table(*tp).fmt(), s.table(*t).fmt(), lpstr(&g)));
                    continue;
                }
                sum.exact_steps += 1;
                let mut mx = 0u128;
                let mut wmx = 0u128;
                let degs: Vec<i32> = (0..g.c.len()).filter(|&i| g.c[i] != 0).map(|i| g.lo + i as i32).collect();
                for &d in &degs {
                    mx = mx.max(unk2(&lam, &src, wt, d));
                    let ws = win_split(&lam, &src, wt);
                    if !ws.labels.is_empty() {
                        wmx = wmx.max(unk2(&ws.labels, &ws.a, &ws.b, d));
                    }
                }
                let kw = kwindow_unknowns(s, *tp, *t, &rho, &degs);
                sum.max_k_unknowns = sum.max_k_unknowns.max(kw.unwrap_or(0));
                if kw.is_none() {
                    sum.k_pad += 1;
                }
                sum.max_unknowns = sum.max_unknowns.max(mx);
                sum.max_window_unknowns = sum.max_window_unknowns.max(wmx);
                if mx > 10_000 {
                    sum.over_1e4 += 1;
                }
                if mx > budget {
                    sum.over_budget += 1;
                }
                if mx > budget {
                    // alternative reduced words of t (≤ 5 rungs): the cheapest word for this step
                    let alts = reduced_words(s, *t, 5, 400_000);
                    let mut best: Option<(u128, Vec<Rung>)> = None;
                    for aw in &alts {
                        let c = degs.iter().map(|&d| unk2(&lam, &src, aw, d)).max().unwrap_or(0);
                        if best.as_ref().map_or(true, |b| c < b.0) {
                            best = Some((c, aw.clone()));
                        }
                    }
                    if let Some((c, aw)) = best {
                        out.push(format!("    (alternatives: {} reduced words of {} with ≤ 5 rungs; cheapest for this step: {} at {} unknowns)", alts.len(), s.table(*t).fmt(), wstr(&aw), sep(c)));
                    }
                }
                out.push(format!(
                    "    step {} -> {}: G = {}: exact Hom^d({} , {}) d in {:?}: {} unknowns (window {}; K-window {}){}",
                    s.table(*tp).fmt(),
                    s.table(*t).fmt(),
                    lpstr(&g),
                    wstr(&src),
                    wstr(wt),
                    degs,
                    sep(mx),
                    sep(wmx),
                    kw.map_or("needs pad".to_string(), sep),
                    if mx > budget { "  ** OVER BUDGET **" } else { "" }
                ));
            }
        }
    }
    sum
}

fn cmd_census(args: &[String]) {
    let lam = parse_lam(&opt(args, "--lam").unwrap_or("3,3,3".into()));
    let budget: u128 = opt(args, "--budget").map(|v| v.parse().unwrap()).unwrap_or(100_000);
    let t0 = Instant::now();
    let mut out: Vec<String> = vec![];
    // the replica check (C_m4c header): 852 011, 392 995, 224 249, 406 634
    let pw = |s: &str| -> Vec<Rung> { krhom::smapio::parse_rungs(s).unwrap() };
    let checks = [
        ("E0.3,F0.3,E1.3,F1.3", 18, 852_011u128),
        ("E0.2,F0.2,E1.3,F1.3", 13, 392_995),
        ("E0.1,F0.1,E1.2,E0.1,F1.2,F0.1", 4, 224_249),
        ("E0.1,F0.1,E1.2,E0.1,F1.2,F0.1", 6, 406_634),
    ];
    out.push("# G-L0 replica check (hu.py numbers of C_m4c)".into());
    for (w, d, want) in checks {
        let w = pw(w);
        let a = hom_unknowns(&lam, &[], &w, d);
        let b = hom_unknowns(&lam, &w, &[], d);
        out.push(format!("Hom^{}(∅, {}) = {} unknowns, reverse {} (hu.py: {}) {}", d, wstr(&w), sep(a), sep(b), sep(want), if a == want && b == want { "OK" } else { "MISMATCH" }));
        assert_eq!((a, b), (want, want), "census replica mismatch");
    }
    let mut s = Schur::new(&lam);
    let mut sc = StepCache::new();
    let mut table = vec![
        "| N | type p | L_p flows | n | product (Stage 0) | whole-ladder call at d0 (M5b) | junctions/level (max) | max W_t flows (len) | steps: Id / exact | max step unknowns | > 10⁴ | > budget | max window unknowns | max K-window unknowns (needs pad) | words: prefix/canonical/bfs | not reduced |".to_string(),
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|".to_string(),
    ];
    for (n, ts) in CLOSURE_FAILS {
        let p = s.intern(Table::parse(ts));
        let w = canonical_ladder(&mut s, p);
        let cd = s.c_diag(p);
        let d0 = cd.min_exp().unwrap();
        let whole = unk2(&lam, &[], &w, d0);
        let prod = krhom::leaves::product_factorization(&mut s, &w, p, 64);
        let prods = match &prod {
            Some(f) => format!("{} · {}{}", s.table(f.a).fmt(), s.table(f.b).fmt(), if f.swaps.is_empty() { String::new() } else { format!(" ({} swap)", f.swaps.len()) }),
            None => "no".into(),
        };
        out.push(format!("\n## N = {} closure {} (c_diag = {})", n, ts, lpstr(&cd)));
        let sm = census_job(&mut s, &mut sc, p, budget, &mut out);
        table.push(format!(
            "| {} | {} | {} | {} | {} | {} | {} | {} ({}) | {} / {} | {} | {} | {} | {} | {} ({}) | {}/{}/{} | {} |",
            n,
            ts,
            sep(flows(&lam, &w)),
            cd.at1(),
            prods,
            sep(whole),
            sm.max_junctions,
            sep(sm.max_word_flows),
            sm.max_word_len,
            sm.id_steps,
            sm.exact_steps,
            sep(sm.max_unknowns),
            sm.over_1e4,
            sm.over_budget,
            sep(sm.max_window_unknowns),
            sep(sm.max_k_unknowns),
            sm.k_pad,
            sm.kinds.get("prefix").unwrap_or(&0),
            sm.kinds.get("canonical").unwrap_or(&0),
            sm.kinds.get("bfs").unwrap_or(&0),
            if sm.not_reduced.is_empty() { "-".into() } else { sm.not_reduced.join(", ") }
        ));
    }
    let mut text = vec![format!("# m4c census G-L0 (λ = {:?}, budget {} unknowns), {:.1} s", lam, budget, t0.elapsed().as_secs_f64()), String::new()];
    text.extend(table);
    text.push(String::new());
    text.extend(out);
    let body = text.join("\n");
    match opt(args, "--out") {
        Some(f) => {
            std::fs::write(&f, &body).unwrap();
            println!("{}", text[..15.min(text.len())].join("\n"));
            println!("written to {}", f);
        }
        None => println!("{}", body),
    }
}

// ------------------------------------------------------------------ Stage 0

struct NoProvider;
fn refuse(what: &str) -> Vec<Missing> {
    vec![Missing { what: format!("replay: {} not in layer A", what), degree: 0, stratum: None, got: 0, want: 1, kind: MissKind::Refused }]
}
impl<F: Field> HomProvider<F> for NoProvider {
    fn name(&self) -> &'static str {
        "NoProvider"
    }
    fn full(&mut self, _: &mut Arena<F>, _: &mut Schur, _: ExprId, _: ExprId, _: i32, _: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        Err(refuse("full"))
    }
    fn pairs(&mut self, _: &mut Arena<F>, _: &mut Schur, _: ExprId, _: Tid, _: ExprId, _: Option<MapId>, _: i32, _: usize, _: &Budget) -> Result<(Pairs<F>, SearchStats), Vec<Missing>> {
        Err(refuse("pairs"))
    }
    fn closure(&mut self, _: &mut Arena<F>, _: &mut Schur, _: ExprId, _: Tid, _: Option<MapId>, _: bool, _: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        Err(refuse("closure"))
    }
    fn end0(&mut self, _: &mut Arena<F>, _: &mut Schur, _: ExprId, _: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>> {
        Err(refuse("end0"))
    }
}

fn stage0<F: Field>(args: &[String]) {
    let lam = parse_lam(&opt(args, "--lam").unwrap_or("3,3,3".into()));
    let n: u8 = opt(args, "--n").map(|v| v.parse().unwrap()).unwrap_or(6);
    let cache = PathBuf::from(opt(args, "--cache").expect("--cache"));
    let costs = opt(args, "--costs").map(PathBuf::from);
    let replay = flag(args, "--replay");
    std::fs::create_dir_all(cache.join("layer_a")).unwrap();
    let (chain, stamp) = provider_for_stamped::<F>(&lam, costs.as_deref(), 1.0);
    let pf = flag(args, "--product-first");
    let stamp = format!("{} {}", stamp, if pf { krhom::leaves::STAGE0_STAMP_PF } else { krhom::leaves::STAGE0_STAMP });
    let prov: Box<dyn HomProvider<F>> = if replay {
        Box::new(NoProvider)
    } else {
        let mut pc = krhom::leaves::ProductClosure::new(chain);
        pc.product_first = pf;
        Box::new(pc)
    };
    let opts = SIndecOpts { cache_dir: Some(cache.clone()), prune_n: Some(n), ..SIndecOpts::default() };
    let mut sd: SIndec<F> = SIndec::with_stamp(&lam, prov, stamp, opts);
    let only_prod = opt(args, "--types").map_or(true, |v| v == "product");
    let mut ok = 0;
    let mut tried = 0;
    for (nn, ts) in CLOSURE_FAILS {
        if nn != n && !(n == 6 && nn == 5) {
            continue;
        }
        let p = sd.schur.intern(Table::parse(ts));
        let w = canonical_ladder(&mut sd.schur, p);
        let fac = krhom::leaves::product_factorization(&mut sd.schur, &w, p, 64);
        if only_prod && fac.is_none() {
            continue;
        }
        tried += 1;
        let t0 = Instant::now();
        match sd.closure_gens(p) {
            Ok(g) => {
                ok += 1;
                let st = sd.stats();
                println!("closure {} OK: ins {} outs {} ({:.2} s; layer A loaded {} saved {} rejected {})", ts, g.ins.len(), g.outs.len(), t0.elapsed().as_secs_f64(), st.layer_a_loaded, st.layer_a_saved, st.layer_a_rejected);
            }
            Err(e) => println!("closure {} FAILED ({:.2} s): {}", ts, t0.elapsed().as_secs_f64(), e),
        }
    }
    let _ = sd.ar.lib_ref().save();
    println!("STAGE0 field={} replay={} certified {}/{}", F::name(), replay, ok, tried);
}

/// G-L5 (downstream cost of a family): DAG nodes, distinct expressions and their resident
/// FlowIndex bytes, and the wall time of a full localization (`Eval::all`) of every map at one
/// fresh point (one Eval memo for the family, as M6a evaluates a library per point).
fn gl5<F: Field>(ar: &Arena<F>, ms: &[MapId]) -> (usize, usize, usize, f64, usize) {
    use krhom::smap::MapNode;
    let mut seen = std::collections::HashSet::new();
    let mut exprs = std::collections::HashSet::new();
    let mut stack: Vec<MapId> = ms.to_vec();
    while let Some(m) = stack.pop() {
        if !seen.insert(m) {
            continue;
        }
        let i = ar.info(m);
        exprs.insert(i.src);
        exprs.insert(i.tgt);
        match ar.node(m) {
            MapNode::Comp(g, f) => stack.extend([*g, *f]),
            MapNode::Lin(t) => stack.extend(t.iter().map(|x| x.1)),
            MapNode::Kron(f, _) => stack.push(*f),
            _ => {}
        }
    }
    let bytes: usize = exprs.iter().filter(|&&e| ar.has_index(e)).map(|&e| ar.flows(e).bytes()).sum();
    let pt = krhom::smap::points::<F>(ar.k(), 0x6151, 1, 7).remove(0);
    let t0 = Instant::now();
    let mut ev = krhom::smap::Eval::new(ar, pt);
    let mut entries = 0usize;
    for &m in ms {
        entries += ev.all(m).nonzero_blocks();
    }
    let secs = t0.elapsed().as_secs_f64();
    (seen.len(), exprs.len(), bytes, secs, ev.memo_bytes() + entries * 0)
}

/// `m4c leaves`: closure generators through the Leaves backend alone (A2/A3).
/// --types fails (the 13 closure jobs) | all (every type of (λ, λ)) | product | T1;T2;...
/// --direct: the provider on a fresh arena (step records and costs printed per type);
/// otherwise a SIndec with layer A in --cache (stamp: Leaves alone), `--replay` = refuse all.
fn leaves_cmd<F: Field>(args: &[String]) {
    use krhom::leaves::*;
    let lam = parse_lam(&opt(args, "--lam").unwrap_or("3,3,3".into()));
    let costs = opt(args, "--costs").map(PathBuf::from);
    let replay = flag(args, "--replay");
    let direct = flag(args, "--direct");
    let types = opt(args, "--types").unwrap_or("fails".into());
    let mut s0 = Schur::new(&lam);
    let list: Vec<String> = match types.as_str() {
        "fails" => CLOSURE_FAILS.iter().filter(|x| x.1 != "1,1,1/1,1,1/1,1,1" || x.0 == 6).map(|x| x.1.to_string()).collect(),
        "all" | "product" => {
            let order = s0.canonical_basis(&lam).order.clone();
            let mut v = vec![];
            for p in order {
                let w = canonical_ladder(&mut s0, p);
                if types == "product" && krhom::leaves::product_factorization(&mut s0, &w, p, 64).is_none() {
                    continue;
                }
                v.push(s0.table(p).fmt());
            }
            v
        }
        l => l.split(';').map(|x| x.to_string()).collect(),
    };
    let mut lv = Leaves::new();
    lv.costs = costs.as_ref().and_then(|p| CostTable::load(p).ok());
    if let Some(b) = opt(args, "--budget") {
        lv.budget = b.parse().unwrap();
    }
    if flag(args, "--no-window") {
        lv.use_window = false;
    }
    let (mut ok, mut tried) = (0, 0);
    let t_all = Instant::now();
    if direct {
        let cache = opt(args, "--cache").map(PathBuf::from);
        let mut ar = Arena::<F>::new(&lam, krhom::tmpl::TmplLib::new(cache.map(|c| c.join("tmpl"))));
        let b = Budget::default();
        let from: usize = opt(args, "--from").map(|v| v.parse().unwrap()).unwrap_or(0);
        for (ix, ts) in list.iter().enumerate().skip(from) {
            let p = s0.intern(Table::parse(ts));
            let w = canonical_ladder(&mut s0, p);
            let lp = ar.expr(krhom::sexpr::Expr::new(&lam, &w).unwrap());
            eprintln!("START #{} {} L_p = {} ({} flows)", ix, ts, wstr(&w), flows(&lam, &w));
            for outs in [false, true] {
                tried += 1;
                let t0 = Instant::now();
                let r = HomProvider::closure(&mut lv, &mut ar, &mut s0, lp, p, None, outs, &b);
                let st = &lv.last;
                let mut kinds: BTreeMap<String, usize> = BTreeMap::new();
                for r in &st.steps {
                    *kinds.entry(format!("{:?}", r.kind)).or_default() += 1;
                }
                let worst = st.steps.iter().max_by(|a, b| a.secs.partial_cmp(&b.secs).unwrap());
                match r {
                    Ok((ms, _)) => {
                        ok += 1;
                        let (nodes, nex, ibytes, esecs, memo) = gl5(&ar, &ms);
                        println!(
                            "LEAVES {} {} OK n={} wall {:.2}s G-L5 dag_nodes {} exprs {} index_bytes {} eval_all {:.4}s memo_bytes {} steps {:?} exact computed {} ({:.2}s) max unknowns {} max W flows {} leaves built {} worst step {}",
                            ts,
                            if outs { "outs" } else { "ins" },
                            ms.len(),
                            t0.elapsed().as_secs_f64(),
                            nodes,
                            nex,
                            ibytes,
                            esecs,
                            memo,
                            kinds,
                            st.exact_calls,
                            st.exact_secs,
                            sep(st.max_unknowns),
                            st.max_word_flows,
                            st.maps_built,
                            worst.map_or("-".into(), |r| format!("l{} {}->{} {:?} {:?} {}u {:.2}s", r.level, s0.table(r.from).fmt(), s0.table(r.to).fmt(), r.kind, r.degs, r.unknowns, r.secs))
                        );
                    }
                    Err(e) => println!("LEAVES {} {} FAILED ({:.2}s): {:?} {}", ts, if outs { "outs" } else { "ins" }, t0.elapsed().as_secs_f64(), e[0].kind, e[0].what),
                }
            }
        }
        let _ = ar.lib_ref().save();
        let by = ar.bytes();
        println!("LEAVES direct field={} certified {}/{} sides, wall {:.1}s, arena bytes {}", F::name(), ok, tried, t_all.elapsed().as_secs_f64(), by.total());
        return;
    }
    let n: u8 = opt(args, "--n").map(|v| v.parse().unwrap()).unwrap_or(6);
    let cache = PathBuf::from(opt(args, "--cache").expect("--cache"));
    std::fs::create_dir_all(cache.join("layer_a")).unwrap();
    let stamp = format!("homprov[{}] maxl={} {}", COST_ALG, lam.iter().max().unwrap(), LEAVES_STAMP);
    let prov: Box<dyn HomProvider<F>> = if replay { Box::new(NoProvider) } else { Box::new(lv) };
    let opts = SIndecOpts { cache_dir: Some(cache.clone()), prune_n: Some(n), ..SIndecOpts::default() };
    let mut sd: SIndec<F> = SIndec::with_stamp(&lam, prov, stamp, opts);
    for ts in &list {
        let p = sd.schur.intern(Table::parse(ts));
        tried += 1;
        let t0 = Instant::now();
        match sd.closure_gens(p) {
            Ok(g) => {
                ok += 1;
                let st = sd.stats();
                let mut din: Vec<i32> = g.ins.iter().map(|&m| sd.ar.info(m).deg).collect();
                din.sort_unstable();
                println!("closure {} OK: ins {} outs {} degrees {:?} ({:.2} s; layer A loaded {} saved {} rejected {})", ts, g.ins.len(), g.outs.len(), din, t0.elapsed().as_secs_f64(), st.layer_a_loaded, st.layer_a_saved, st.layer_a_rejected);
            }
            Err(e) => println!("closure {} FAILED ({:.2} s): {}", ts, t0.elapsed().as_secs_f64(), e),
        }
    }
    let _ = sd.ar.lib_ref().save();
    println!("LEAVES sindec field={} replay={} certified {}/{} types, wall {:.1}s, arena bytes {}", F::name(), replay, ok, tried, t_all.elapsed().as_secs_f64(), sd.ar.bytes().total());
}

/// peak resident set of this process (VmHWM, MB)
fn peak_rss_mb() -> f64 {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|t| t.lines().find(|l| l.starts_with("VmHWM:")).and_then(|l| l.split_whitespace().nth(1)).and_then(|v| v.parse::<f64>().ok()))
        .map_or(0.0, |kb| kb / 1024.0)
}

/// `m4c measure --base IN --out OUT [--cache DIR] CALL...` with CALL = "labels;src;tgt;deg"
/// (words as E0.2,E1.3; "-" = empty): one exact `hom_space` call each, in this process, recorded
/// as an uncensored exact-shape cost record (wall seconds, peak RSS) in OUT (a copy of IN plus the
/// new records). Run it under `timeout` / `ulimit -v`: a killed call writes nothing (and so stays
/// unadmitted above the 10⁵ budget).
fn measure_cmd(args: &[String]) {
    let base = opt(args, "--base").expect("--base");
    let out = opt(args, "--out").expect("--out");
    let cache = opt(args, "--cache").map(PathBuf::from);
    let mut table = CostTable::load(std::path::Path::new(&out)).or_else(|_| CostTable::load(std::path::Path::new(&base))).expect("cost table");
    let calls: Vec<&String> = args.iter().filter(|a| a.matches(';').count() == 3).collect();
    let mut lib: krhom::tmpl::TmplLib<F31> = krhom::tmpl::TmplLib::new(cache);
    for c in calls {
        let f: Vec<&str> = c.split(';').collect();
        let lab = parse_lam(f[0]);
        let pw = |x: &str| -> Vec<Rung> { if x == "-" { vec![] } else { krhom::smapio::parse_rungs(x).unwrap() } };
        let (a, b) = (pw(f[1]), pw(f[2]));
        let d: i32 = f[3].parse().unwrap();
        let (fs, ft, unk) = (flows(&lab, &a), flows(&lab, &b), hom_unknowns(&lab, &a, &b, d));
        let t0 = Instant::now();
        let before = lib.stats().computed;
        let ids = lib.try_hom_basis_cap(&krhom::tmpl::Window::new(&lab, &a), &krhom::tmpl::Window::new(&lab, &b), d, u128::MAX).unwrap();
        let secs = t0.elapsed().as_secs_f64();
        if lib.stats().computed == before {
            println!("MEASURE {} : loaded from cache, not a measurement", c);
            continue;
        }
        let rss = peak_rss_mb();
        table.record(CostRec { labels: lab.clone(), fs, ft, unknowns: unk, deg: d, secs, rss_mb: rss, censored: false });
        table.save(std::path::Path::new(&out)).unwrap();
        println!("MEASURE {} : fibres {}x{} unknowns {} dim {} {:.2}s peak RSS {:.0} MB", c, fs, ft, unk, ids.len(), secs, rss);
    }
    let _ = lib.save();
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.first().map(|s| s.as_str()) {
        Some("census") => cmd_census(&args),
        Some("stage0") => match opt(&args, "--field").as_deref() {
            Some("31") | None => stage0::<F31>(&args),
            Some("32b") => stage0::<F32b>(&args),
            Some(f) => panic!("unknown field {}", f),
        },
        Some("measure") => measure_cmd(&args),
        Some("leaves") => match opt(&args, "--field").as_deref() {
            Some("31") | None => leaves_cmd::<F31>(&args),
            Some("32b") => leaves_cmd::<F32b>(&args),
            Some(f) => panic!("unknown field {}", f),
        },
        _ => eprintln!("usage: m4c census [--lam 3,3,3] [--budget N] [--out FILE] | m4c stage0 --field 31|32b --cache DIR [--costs PATH] [--replay] [--types product|all]"),
    }
}
