//! G1 harness (SPEC_M4_M6a §4.7).
//!
//! g1 a [--only LAMBDA] [--json PATH] [--cache DIR] [--costs PATH] [--wall S] [--rss MB]
//! g1 b [--pairs N] [--seed S] [--json PATH] [--cache DIR] [--costs PATH]   (cloud)
//! g1 one LAMBDA PAIR_INDEX [--cache DIR] [--costs PATH] [--lmax L] [--cmax C]   (child process)
//! g1 closure LAMBDA TYPE_INDEX [...]                                         (child process)
//! g1 calib [--cache DIR] [--out PATH] [--per-class N] [--limit S]                (cost table)
//! exploration: g1 list LAMBDA | g1 plan LAMBDA WMAX LMAX | g1 time LABELS SRC TGT DEG
//!
//! The parent runs one child process per pair, sequentially, and polls its RSS (`ps -o rss=`)
//! and wall time; a child above the RSS cap or the wall limit is killed and recorded as
//! MissKind::Memory or MissKind::Time for its pair.
use krhom::homprov::*;
use krhom::schur::{flows, Schur};
use krhom::web::{Dir, Rung};

fn rstr(r: &Rung) -> String {
    format!("{}{}({})", if r.dir == Dir::E { "E" } else { "F" }, r.i, r.r)
}
fn wstr(rs: &[Rung]) -> String {
    rs.iter().map(rstr).collect::<Vec<_>>().join(" ")
}

/// common prefix and suffix of two words with equal bottoms and tops; the window between them
/// (on the uprights the differing rungs touch): labels, rungs (relative), fibres
fn window_info(lam: &[u8], a: &[Rung], b: &[Rung]) -> (usize, usize, Vec<u8>, Vec<Rung>, Vec<Rung>, u128, u128) {
    let mut pre = 0;
    while pre < a.len() && pre < b.len() && a[pre] == b[pre] {
        pre += 1;
    }
    let mut suf = 0;
    while suf < a.len() - pre && suf < b.len() - pre && a[a.len() - 1 - suf] == b[b.len() - 1 - suf] {
        suf += 1;
    }
    let wa = &a[pre..a.len() - suf];
    let wb = &b[pre..b.len() - suf];
    let ex = krhom::sexpr::Expr::new(lam, &a[..pre]).unwrap();
    let lv = ex.top();
    if wa.is_empty() && wb.is_empty() {
        return (pre, suf, vec![], vec![], vec![], 1, 1);
    }
    let lo = wa.iter().chain(wb).map(|r| r.i).min().unwrap();
    let hi = wa.iter().chain(wb).map(|r| r.i + 1).max().unwrap();
    let wl = lv[lo..=hi].to_vec();
    let sh = |rs: &[Rung]| -> Vec<Rung> { rs.iter().map(|r| Rung { i: r.i - lo, dir: r.dir, r: r.r }).collect() };
    let (ra, rb) = (sh(wa), sh(wb));
    let fa = flows(&wl, &ra);
    let fb = flows(&wl, &rb);
    (pre, suf, wl, ra, rb, fa, fb)
}

/// windows of a word: contiguous rungs at..at+len on uprights first..first+w (w = 2 or 3)
fn windows(w: &[Rung], wmax: usize, lmax: usize) -> Vec<WinPos> {
    let mut out = vec![];
    for ww in 2..=wmax {
        for at in 0..w.len() {
            for len in 2..=lmax {
                if at + len > w.len() {
                    break;
                }
                let lo = w[at..at + len].iter().map(|r| r.i).min().unwrap();
                let hi = w[at..at + len].iter().map(|r| r.i + 1).max().unwrap();
                if hi - lo + 1 != ww {
                    continue;
                }
                out.push(WinPos { at, len, first: lo, w: ww });
            }
        }
    }
    out
}

fn reduce(lam: &[u8], word: &[Rung], sc: &mut SchurCache, wmax: usize, lmax: usize, stats: &mut [usize; 4]) -> Option<Vec<(Vec<Rung>, i32)>> {
    let e = sc.expand_word(lam, word);
    if e.len() == 1 && e[0].1 == 0 && e[0].2 == 1 {
        return Some(vec![(word.to_vec(), 0)]);
    }
    // the cheapest decomposable window over the whole commutation class
    let mut best: Option<(u128, Vec<Rung>, WinPos, Vec<u8>, Vec<Rung>)> = None;
    for wd in commutation_class(word, 200) {
        for p in windows(&wd, wmax, lmax) {
            let lv = level_of(lam, &wd, p.at);
            let lab = lv[p.first..p.first + p.w].to_vec();
            let rel = window_word(&wd, &p);
            let ex = sc.expand_word(&lab, &rel);
            if ex.len() == 1 && ex[0].1 == 0 && ex[0].2 == 1 {
                continue;
            }
            let mut cost = 0u128;
            for (lw, sh, _, _) in &ex {
                cost = cost.max(hom_unknowns(&lab, lw, &rel, *sh)).max(hom_unknowns(&lab, &rel, lw, -*sh));
            }
            if best.as_ref().map_or(true, |b| cost < b.0) {
                best = Some((cost, wd.clone(), p.clone(), lab, rel));
            }
        }
    }
    if let Some((cost, wd, p, lab, rel)) = best {
        let ex = sc.expand_word(&lab, &rel);
        stats[p.w - 2] += 1;
        println!("  SPLIT w{} len {} [{}] cost {} lab {:?} parts {:?}", p.w, p.len, wstr(&rel), cost, lab, ex.iter().map(|e| (wstr(&e.0), e.1)).collect::<Vec<_>>());
        let mut out = vec![];
        for (lw, sh, c, _) in ex {
            let mut child: Vec<Rung> = wd[..p.at].to_vec();
            child.extend(lw.iter().map(|r| Rung { i: r.i + p.first, dir: r.dir, r: r.r }));
            child.extend_from_slice(&wd[p.at + p.len..]);
            let sub = reduce(lam, &child, sc, wmax, lmax, stats)?;
            for _ in 0..c {
                for (lf, s2) in &sub {
                    out.push((lf.clone(), sh + s2));
                }
            }
        }
        return Some(out);
    }
    stats[3] += 1;
    None
}

fn plan(lam: &[u8], wmax: usize, lmax: usize) {
    let mut s = Schur::new(lam);
    let (pairs, _) = g1_targets(&mut s);
    let mut sc = SchurCache::default();
    let (mut ok, mut stuck, mut words_eq, mut words_comm, mut words_other) = (0, 0, 0, 0, 0);
    let mut stats = [0usize; 4];
    for (k, pr) in pairs.iter().enumerate() {
        if pr.identity {
            continue;
        }
        let lp = s.canonical_basis(&pr.mu).ladder_minflow[&pr.p].clone();
        let mut x = lp.clone();
        x.push(pr.rung);
        match reduce(lam, &x, &mut sc, wmax, lmax, &mut stats) {
            None => {
                stuck += 1;
                println!("#{} STUCK X [{}]", k, wstr(&x));
            }
            Some(leaves) => {
                ok += 1;
                for (lf, sh) in &leaves {
                    // identify the leaf type
                    let e = sc.expand_word(lam, lf);
                    let lq = e[0].0.clone();
                    if *lf == lq {
                        words_eq += 1;
                    } else if commutation_class(&lq, 500).contains(lf) {
                        words_comm += 1;
                    } else {
                        words_other += 1;
                        let w = win_split(lam, lf, &lq);
                        let mut best = u128::MAX;
                        for a in commutation_class(lf, 300) { for b in commutation_class(&lq, 300) { let ww = win_split(lam, &a, &b); best = best.min(ww.fa * ww.fb); } }
                        println!("#{} leaf [{}] s {} vs Lq [{}]: window {}x{} best-after-commutation {}", k, wstr(lf), sh, wstr(&lq), w.fa, w.fb, best);
                    }
                }
            }
        }
    }
    println!("plan {:?} wmax {} lmax {}: decomposed {} stuck {}; leaves equal {} commutation {} other {}; window splits w2 {} w3 {}", lam, wmax, lmax, ok, stuck, words_eq, words_comm, words_other, stats[0], stats[1]);
}

fn parse_word(s: &str) -> Vec<Rung> {
    // "E0.1,F1.2"
    if s.is_empty() || s == "-" {
        return vec![];
    }
    s.split(',').map(|t| {
        let dir = if t.starts_with('E') { Dir::E } else { Dir::F };
        let (i, r) = t[1..].split_once('.').unwrap();
        Rung { i: i.parse().unwrap(), dir, r: r.parse().unwrap() }
    }).collect()
}

fn explore(args: &[String]) {
    if args[1] == "time" {
        use krhom::field::F31;
        let lab: Vec<u8> = args[2].split(',').map(|x| x.parse().unwrap()).collect();
        let a = parse_word(&args[3]);
        let b = parse_word(&args[4]);
        let d: i32 = args[5].parse().unwrap();
        let unk = hom_unknowns(&lab, &a, &b, d);
        let t0 = std::time::Instant::now();
        let mut lib: krhom::tmpl::TmplLib<F31> = krhom::tmpl::TmplLib::new(None);
        let ids = lib.try_hom_basis_cap(&krhom::tmpl::Window::new(&lab, &a), &krhom::tmpl::Window::new(&lab, &b), d, u128::MAX).unwrap();
        println!("time {:?} [{}] -> [{}] d {} fibres {}x{} unknowns {} dim {} secs {:.3}", lab, wstr(&a), wstr(&b), d, flows(&lab, &a), flows(&lab, &b), unk, ids.len(), t0.elapsed().as_secs_f64());
        return;
    }
    if args[1] == "planword" {
        let lam: Vec<u8> = args[2].split(',').map(|x| x.parse().unwrap()).collect();
        let w = parse_word(&args[3]);
        let mut sc = SchurCache::default();
        let mut stats = [0usize; 4];
        let r = reduce(&lam, &w, &mut sc, 3, 4, &mut stats);
        println!("leaves {:?}", r.map(|v| v.iter().map(|(a, b)| (wstr(a), *b)).collect::<Vec<_>>()));
        return;
    }
    if args[1] == "plan" {
        let lam: Vec<u8> = args[2].split(',').map(|x| x.parse().unwrap()).collect();
        let wmax: usize = args[3].parse().unwrap();
        let lmax: usize = args[4].parse().unwrap();
        plan(&lam, wmax, lmax);
        return;
    }
    let lam: Vec<u8> = args[2].split(',').map(|x| x.parse().unwrap()).collect();
    let mut s = Schur::new(&lam);
    let (pairs, ctypes) = g1_targets(&mut s);
    let mut ni = 0;
    let mut nt = 0;
    let mut sm = 0;
    for (k, pr) in pairs.iter().enumerate() {
        let lp = s.canonical_basis(&pr.mu).ladder_minflow[&pr.p].clone();
        let mut x = lp.clone();
        x.push(pr.rung);
        if pr.identity {
            continue;
        }
        ni += 1;
        println!("#{} mu {:?} p {} | Lp [{}] R {} | X flows {}", k, pr.mu, s.table(pr.p).fmt(), wstr(&lp), rstr(&pr.rung), flows(&lam, &x));
        for &(q, sd, c) in &pr.targets {
            let mq = s.table(q).cols();
            let lq = s.canonical_basis(&mq).ladder_minflow[&q].clone();
            nt += 1;
            sm += c;
            let (pre, suf, wl, wa, wb, fa, fb) = window_info(&lam, &lq, &x);
            let ca = commutation_class(&lq, 500);
            let cb = commutation_class(&x, 500);
            let mut best = u128::MAX;
            for a in &ca { for b in &cb { let w = win_split(&lam, a, b); best = best.min(w.fa * w.fb); } }
            println!("    BEST {} classes {}x{}", best, ca.len(), cb.len());
            println!("    q {} s {} x{} | Lq [{}] flows {} | prefix {} suffix {} window {:?} [{}] -> [{}] fibres {}x{}={}", s.table(q).fmt(), sd, c, wstr(&lq), flows(&lam, &lq), pre, suf, wl, wstr(&wa), wstr(&wb), fa, fb, fa*fb);
        }
    }
    println!("non-identity pairs {} targets {} sum mult {} closure types {}", ni, nt, sm, ctypes.len());
}

// ------------------------------------------------------------------ harness

use krhom::field::F31;
use krhom::sexpr::Expr;
use krhom::smap::Arena;
use krhom::tmpl::TmplLib;
use std::io::Write;
use std::path::PathBuf;
use std::time::Instant;

type F = F31;

fn opt(args: &[String], key: &str) -> Option<String> {
    args.iter().position(|a| a == key).and_then(|i| args.get(i + 1).cloned())
}
fn parse_lam(s: &str) -> Vec<u8> {
    s.split(',').map(|x| x.parse().unwrap()).collect()
}
fn lam_str(l: &[u8]) -> String {
    l.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")
}
fn jstr(s: &str) -> String {
    let mut o = String::from("\"");
    for c in s.chars() {
        match c {
            '"' => o.push_str("\\\""),
            '\\' => o.push_str("\\\\"),
            '\n' => o.push_str("\\n"),
            c if (c as u32) < 0x20 => o.push(' '),
            c => o.push(c),
        }
    }
    o.push('"');
    o
}
fn default_cache() -> PathBuf {
    std::env::var("KRHOM_TMPL_CACHE").map(PathBuf::from).unwrap_or_else(|_| std::env::temp_dir().join("krhom_tmpl_cache"))
}
fn exact_from(args: &[String]) -> Exact {
    let mut ex = Exact::default();
    if let Some(p) = opt(args, "--costs") {
        match CostTable::load(std::path::Path::new(&p)) {
            Ok(t) => ex.costs = Some(t),
            Err(e) => eprintln!("[g1] cost table {} not loaded: {}", p, e),
        }
    }
    ex
}
fn budget_from(args: &[String]) -> Budget {
    let mut b = Budget::default();
    if let Some(v) = opt(args, "--lmax") {
        b.l_max = v.parse().unwrap();
    }
    if let Some(v) = opt(args, "--cmax") {
        b.c_max = v.parse().unwrap();
    }
    b
}
fn kind_str(k: MissKind) -> &'static str {
    match k {
        MissKind::Exhausted => "Exhausted",
        MissKind::Time => "Time",
        MissKind::Memory => "Memory",
        MissKind::Refused => "Refused",
    }
}
fn tiers_str(st: &SearchStats) -> String {
    let t: Vec<String> = st.tier_hits.iter().enumerate().filter(|(_, c)| **c > 0).map(|(i, c)| format!("S{}x{}", i, c)).collect();
    if t.is_empty() { "Exact".to_string() } else { t.join("+") }
}
fn stats_json(st: &SearchStats) -> String {
    format!(
        "{{\"tiers\":{},\"candidates\":{},\"evaluated\":{},\"secs\":{:.3},\"max_flows\":{},\"points\":{},\"visited_words\":{},\"indexed_words\":{},\"peak_index_bytes\":{},\"peak_arena_bytes\":{},\"exact_calls\":{},\"exact_secs\":{:.3},\"detect\":{}}}",
        jstr(&tiers_str(st)), st.candidates, st.evaluated, st.secs, st.max_flows, st.points, st.visited_words, st.indexed_words, st.peak_index_bytes, st.peak_arena_bytes, st.exact_calls, st.exact_secs, jstr(st.detect)
    )
}
fn missing_json(ms: &[Missing]) -> String {
    let v: Vec<String> = ms.iter().map(|m| format!("{{\"kind\":\"{}\",\"degree\":{},\"stratum\":{},\"got\":{},\"want\":{},\"what\":{}}}", kind_str(m.kind), m.degree, jstr(m.stratum.as_deref().unwrap_or("")), m.got, m.want, jstr(&m.what))).collect();
    format!("[{}]", v.join(","))
}
fn max_label(l: &[u8]) -> u8 {
    l.iter().copied().max().unwrap_or(0)
}

/// child: every summand target of one pair
fn run_one(args: &[String]) {
    let lam = parse_lam(&args[2]);
    let idx: usize = args[3].parse().unwrap();
    let cache = opt(args, "--cache").map(PathBuf::from).unwrap_or_else(default_cache);
    let mut s = Schur::new(&lam);
    let (pairs, _) = g1_targets(&mut s);
    let pr = &pairs[idx];
    let lp = s.canonical_basis(&pr.mu).ladder_minflow[&pr.p].clone();
    let mut xw = lp.clone();
    xw.push(pr.rung);
    let mut ar: Arena<F> = Arena::new(&lam, TmplLib::new(Some(cache)));
    let x = ar.expr(Expr::new(&lam, &xw).unwrap());
    let mut prov: Chain<F> = Chain::for_labels(max_label(&lam), exact_from(args));
    let b = budget_from(args);
    let t0 = Instant::now();
    let mut lines = vec![];
    for &(q, sd, mult) in &pr.targets {
        let mq = s.table(q).cols();
        let lqw = s.canonical_basis(&mq).ladder_minflow[&q].clone();
        let lq = ar.expr(Expr::new(&lam, &lqw).unwrap());
        let t1 = Instant::now();
        let r = prov.pairs(&mut ar, &mut s, x, q, lq, None, sd, mult, &b);
        let line = match r {
            Ok((p, st)) => format!("{{\"q\":{},\"s\":{},\"mult\":{},\"ok\":true,\"got\":{},\"secs\":{:.3},\"stats\":{}}}", jstr(&s.table(q).fmt()), sd, mult, p.iotas.len(), t1.elapsed().as_secs_f64(), stats_json(&st)),
            Err(ms) => format!("{{\"q\":{},\"s\":{},\"mult\":{},\"ok\":false,\"secs\":{:.3},\"missing\":{}}}", jstr(&s.table(q).fmt()), sd, mult, t1.elapsed().as_secs_f64(), missing_json(&ms)),
        };
        eprintln!("[g1 one] {} {} {}: {}", lam_str(&lam), idx, s.table(q).fmt(), &line[..line.len().min(300)]);
        lines.push(line);
    }
    let by = ar.bytes();
    let _ = ar.lib_ref().save();
    println!(
        "RESULT {{\"lambda\":\"{}\",\"pair\":{},\"p\":{},\"rung\":{},\"x\":{},\"x_flows\":{},\"wall\":{:.3},\"arena_bytes\":{},\"targets\":[{}]}}",
        lam_str(&lam), idx, jstr(&s.table(pr.p).fmt()), jstr(&rstr(&pr.rung)), jstr(&wstr(&xw)), flows(&lam, &xw), t0.elapsed().as_secs_f64(), by.total(), lines.join(",")
    );
}

/// child: closure ins and outs of one type of (lambda, lambda)
fn run_closure(args: &[String]) {
    let lam = parse_lam(&args[2]);
    let idx: usize = args[3].parse().unwrap();
    let cache = opt(args, "--cache").map(PathBuf::from).unwrap_or_else(default_cache);
    let mut s = Schur::new(&lam);
    let (_, ctypes) = g1_targets(&mut s);
    let p = ctypes[idx];
    let lpw = s.canonical_basis(&lam).ladder_minflow[&p].clone();
    let mut ar: Arena<F> = Arena::new(&lam, TmplLib::new(Some(cache)));
    let lp = ar.expr(Expr::new(&lam, &lpw).unwrap());
    let mut prov: Chain<F> = Chain::for_labels(max_label(&lam), exact_from(args));
    let b = budget_from(args);
    let t0 = Instant::now();
    let mut lines = vec![];
    for outs in [false, true] {
        let t1 = Instant::now();
        let r = prov.closure(&mut ar, &mut s, lp, p, None, outs, &b);
        let n = { let cd = s.c_diag(p); cd.at1() };
        let line = match r {
            Ok((ms, st)) => format!("{{\"side\":\"{}\",\"n\":{},\"ok\":true,\"got\":{},\"secs\":{:.3},\"stats\":{}}}", if outs { "outs" } else { "ins" }, n, ms.len(), t1.elapsed().as_secs_f64(), stats_json(&st)),
            Err(ms) => format!("{{\"side\":\"{}\",\"n\":{},\"ok\":false,\"secs\":{:.3},\"missing\":{}}}", if outs { "outs" } else { "ins" }, n, t1.elapsed().as_secs_f64(), missing_json(&ms)),
        };
        eprintln!("[g1 closure] {} {}: {}", lam_str(&lam), s.table(p).fmt(), &line[..line.len().min(300)]);
        lines.push(line);
    }
    let _ = ar.lib_ref().save();
    println!("RESULT {{\"lambda\":\"{}\",\"closure_type\":{},\"p\":{},\"word\":{},\"flows\":{},\"wall\":{:.3},\"targets\":[{}]}}", lam_str(&lam), idx, jstr(&s.table(p).fmt()), jstr(&wstr(&lpw)), flows(&lam, &lpw), t0.elapsed().as_secs_f64(), lines.join(","));
}

/// run a child under the watchdog: (stdout RESULT line or None, wall, peak RSS MB, killed reason)
fn supervised(cmd: &[String], wall: f64, rss_mb: f64) -> (Option<String>, f64, f64, Option<&'static str>) {
    let exe = std::env::current_exe().unwrap();
    let mut child = std::process::Command::new(exe).args(cmd).stdout(std::process::Stdio::piped()).stderr(std::process::Stdio::inherit()).spawn().expect("spawn child");
    let pid = child.id();
    let t0 = Instant::now();
    let mut peak = 0.0f64;
    let mut killed = None;
    loop {
        match child.try_wait() {
            Ok(Some(_)) => break,
            Ok(None) => {}
            Err(_) => break,
        }
        if let Ok(o) = std::process::Command::new("ps").args(["-o", "rss=", "-p", &pid.to_string()]).output() {
            if let Ok(v) = String::from_utf8_lossy(&o.stdout).trim().parse::<f64>() {
                peak = peak.max(v / 1024.0);
            }
        }
        if peak > rss_mb {
            let _ = child.kill();
            killed = Some("Memory");
            break;
        }
        if t0.elapsed().as_secs_f64() > wall {
            let _ = child.kill();
            killed = Some("Time");
            break;
        }
        std::thread::sleep(std::time::Duration::from_millis(250));
    }
    let out = child.wait_with_output().ok();
    let line = out.and_then(|o| String::from_utf8_lossy(&o.stdout).lines().find(|l| l.starts_with("RESULT ")).map(|l| l[7..].to_string()));
    (if killed.is_some() { None } else { line }, t0.elapsed().as_secs_f64(), peak, killed)
}

fn run_gate(args: &[String], which: &str) {
    let wall: f64 = opt(args, "--wall").map(|v| v.parse().unwrap()).unwrap_or(600.0);
    let rss: f64 = opt(args, "--rss").map(|v| v.parse().unwrap()).unwrap_or(1024.0);
    let json = opt(args, "--json");
    let mut pass_args: Vec<String> = vec![];
    for k in ["--cache", "--costs", "--lmax", "--cmax"] {
        if let Some(v) = opt(args, k) {
            pass_args.push(k.into());
            pass_args.push(v);
        }
    }
    let lams: Vec<Vec<u8>> = if which == "a" {
        match opt(args, "--only") {
            Some(l) => vec![parse_lam(&l)],
            None => vec![vec![3, 3], vec![2, 2, 2]],
        }
    } else {
        vec![vec![2, 2, 2, 2], vec![3, 3, 3]]
    };
    let seed: u64 = opt(args, "--seed").map(|v| v.parse().unwrap()).unwrap_or(2026);
    let t_all = Instant::now();
    let mut records: Vec<String> = vec![];
    let (mut n_t, mut n_ok, mut n_pairs, mut n_pairs_ok) = (0usize, 0usize, 0usize, 0usize);
    let max_pairs: Option<usize> = opt(args, "--max-pairs").map(|v| v.parse().unwrap());
    for lam in &lams {
        let mut s = Schur::new(lam);
        let (pairs, ctypes) = g1_targets(&mut s);
        let mut idxs: Vec<usize> = (0..pairs.len()).filter(|&i| !pairs[i].identity).collect();
        let mut cidx: Vec<usize> = (0..ctypes.len()).collect();
        if which == "b" {
            use rand::{rngs::StdRng, seq::SliceRandom, SeedableRng};
            let (np, nc) = if lam.len() == 4 { (50, 10) } else { (20, 5) };
            let np = opt(args, "--pairs").map(|v| v.parse::<usize>().unwrap().min(np)).unwrap_or(np);
            let mut rng = StdRng::seed_from_u64(seed);
            idxs.shuffle(&mut rng);
            idxs.truncate(np);
            cidx.shuffle(&mut rng);
            cidx.truncate(nc);
        }
        if let Some(m) = max_pairs {
            idxs.truncate(m);
        }
        let mut jobs: Vec<Vec<String>> = idxs.iter().map(|&i| vec!["one".to_string(), lam_str(lam), i.to_string()]).collect();
        jobs.extend(cidx.iter().map(|&i| vec!["closure".to_string(), lam_str(lam), i.to_string()]));
        for job in jobs {
            let mut cmd = job.clone();
            cmd.extend(pass_args.iter().cloned());
            let (line, w, peak, killed) = supervised(&cmd, wall, rss);
            n_pairs += 1;
            let rec = match (line, killed) {
                (Some(l), _) => {
                    let ok_t = l.matches("\"ok\":true").count();
                    let bad_t = l.matches("\"ok\":false").count();
                    n_t += ok_t + bad_t;
                    n_ok += ok_t;
                    if bad_t == 0 {
                        n_pairs_ok += 1;
                    }
                    format!("{{\"job\":{},\"wall\":{:.2},\"peak_rss_mb\":{:.1},\"killed\":null,\"result\":{}}}", jstr(&job.join(" ")), w, peak, l)
                }
                (None, k) => {
                    // count the job's targets as failed
                    let nt = if job[0] == "one" { pairs[job[2].parse::<usize>().unwrap()].targets.len() } else { 2 };
                    n_t += nt;
                    format!("{{\"job\":{},\"wall\":{:.2},\"peak_rss_mb\":{:.1},\"killed\":{},\"result\":null}}", jstr(&job.join(" ")), w, peak, jstr(k.unwrap_or("crash")))
                }
            };
            eprintln!("[g1 {}] {} done in {:.1}s (peak {:.0} MB); targets ok {}/{} so far", which, job.join(" "), w, peak, n_ok, n_t);
            records.push(rec);
            if let Some(p) = &json {
                let _ = write_json(p, which, &records, n_t, n_ok, n_pairs, n_pairs_ok, t_all.elapsed().as_secs_f64(), false);
            }
        }
    }
    let total = t_all.elapsed().as_secs_f64();
    println!("G1({}) targets certified {}/{} ({:.1}%), jobs fully certified {}/{}, wall {:.1}s", which, n_ok, n_t, 100.0 * n_ok as f64 / n_t.max(1) as f64, n_pairs_ok, n_pairs, total);
    if let Some(p) = &json {
        write_json(p, which, &records, n_t, n_ok, n_pairs, n_pairs_ok, total, true).unwrap();
    }
}

#[allow(clippy::too_many_arguments)]
fn write_json(path: &str, which: &str, recs: &[String], n_t: usize, n_ok: usize, n_p: usize, n_p_ok: usize, wall: f64, done: bool) -> std::io::Result<()> {
    let tmp = format!("{}.{}.tmp", path, std::process::id());
    {
        let mut f = std::io::BufWriter::new(std::fs::File::create(&tmp)?);
        writeln!(f, "{{\"gate\":\"G1({})\",\"complete\":{},\"targets\":{},\"certified\":{},\"jobs\":{},\"jobs_certified\":{},\"wall_s\":{:.1},\"records\":[", which, done, n_t, n_ok, n_p, n_p_ok, wall)?;
        for (i, r) in recs.iter().enumerate() {
            writeln!(f, "{}{}", r, if i + 1 < recs.len() { "," } else { "" })?;
        }
        writeln!(f, "]}}")?;
    }
    std::fs::rename(tmp, path)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    match args.get(1).map(|s| s.as_str()) {
        Some("one") => run_one(&args),
        Some("closure") => run_closure(&args),
        Some("a") => run_gate(&args, "a"),
        Some("b") => run_gate(&args, "b"),
        Some("calib") => calib(&args),
        Some("calib-one") => calib_one(&args),
        _ => explore(&args),
    }
}

/// child: one exact window call, prints "COST <record>"
fn calib_one(args: &[String]) {
    let lab = parse_lam(&args[2]);
    let a = parse_word(&args[3]);
    let b = parse_word(&args[4]);
    let d: i32 = args[5].parse().unwrap();
    let cache = opt(args, "--cache").map(PathBuf::from);
    let mut lib: TmplLib<F> = TmplLib::new(cache);
    let t0 = Instant::now();
    let ids = lib.try_hom_basis_cap(&krhom::tmpl::Window::new(&lab, &a), &krhom::tmpl::Window::new(&lab, &b), d, u128::MAX).unwrap();
    let secs = t0.elapsed().as_secs_f64();
    let computed = lib.stats().computed;
    println!("RESULT {} {} {} {} {} {:.6} {} {}", lam_str(&lab), flows(&lab, &a), flows(&lab, &b), hom_unknowns(&lab, &a, &b, d), d, secs, ids.len(), computed);
}

fn wcode(rs: &[Rung]) -> String {
    if rs.is_empty() {
        return "-".into();
    }
    rs.iter().map(|r| format!("{}{}.{}", if r.dir == Dir::E { "E" } else { "F" }, r.i, r.r)).collect::<Vec<_>>().join(",")
}

/// Cost table calibration: >= N measured calls per label class. Classes: the whole-ladder G1(a)
/// pair calls of (3,3) and (2,2,2) (Hom^s(L_q, X) and Hom^{-s}(X, L_q), smallest fibre products
/// first) and the window calls of the local decomposition planner, each in a child process under
/// the watchdog (a kill is a censored record at the time limit).
fn calib(args: &[String]) {
    let per: usize = opt(args, "--per-class").map(|v| v.parse().unwrap()).unwrap_or(20);
    let limit: f64 = opt(args, "--limit").map(|v| v.parse().unwrap()).unwrap_or(60.0);
    let out = opt(args, "--out").unwrap_or_else(|| "costs.txt".into());
    let budget_s: f64 = opt(args, "--budget").map(|v| v.parse().unwrap()).unwrap_or(1500.0);
    let cache = opt(args, "--cache");
    let mut table = CostTable::load(std::path::Path::new(&out)).unwrap_or_default();
    if let Some(base) = opt(args, "--base") {
        if table.is_empty() {
            table = CostTable::load(std::path::Path::new(&base)).expect("base cost table");
        }
    }
    let done: std::collections::HashSet<(Vec<u8>, u128, u128, u128, i32)> = table.recs.iter().map(|r| (r.labels.clone(), r.fs, r.ft, r.unknowns, r.deg)).collect();
    let mut calls: Vec<(Vec<u8>, Vec<Rung>, Vec<Rung>, i32)> = vec![];
    for lam in [vec![3u8, 3], vec![2, 2, 2]] {
        let mut s = Schur::new(&lam);
        // closure calls: [] <-> L_p at every c_diag degree of every type of (lambda, lambda)
        let order = s.canonical_basis(&lam).order.clone();
        for &p in &order {
            let lp = s.canonical_basis(&lam).ladder_minflow[&p].clone();
            if lp.is_empty() {
                continue;
            }
            let cd = s.c_diag(p);
            for (i, &c) in cd.c.iter().enumerate() {
                if c != 0 {
                    let d = cd.lo + i as i32;
                    calls.push((lam.clone(), vec![], lp.clone(), d));
                    calls.push((lam.clone(), lp.clone(), vec![], d));
                }
            }
        }
        // pair calls: Hom^s(L_q, X) and Hom^{-s}(X, L_q), the smallest 20 by unknowns and 6 more
        // spread over the rest
        let (pairs, _) = g1_targets(&mut s);
        let mut cls: Vec<(u128, Vec<Rung>, Vec<Rung>, i32)> = vec![];
        for pr in pairs.iter().filter(|p| !p.identity) {
            let lp = s.canonical_basis(&pr.mu).ladder_minflow[&pr.p].clone();
            let mut xw = lp.clone();
            xw.push(pr.rung);
            for &(q, sd, _) in &pr.targets {
                let mq = s.table(q).cols();
                let lq = s.canonical_basis(&mq).ladder_minflow[&q].clone();
                cls.push((hom_unknowns(&lam, &lq, &xw, sd), lq.clone(), xw.clone(), sd));
                cls.push((hom_unknowns(&lam, &xw, &lq, -sd), xw.clone(), lq.clone(), -sd));
            }
        }
        cls.sort_by(|a, b| a.0.cmp(&b.0));
        cls.dedup_by(|a, b| a.1 == b.1 && a.2 == b.2 && a.3 == b.3);
        let n = cls.len();
        let mut pick: Vec<usize> = (0..per.min(n)).collect();
        for k in 1..=6 {
            let i = per + (n.saturating_sub(per)) * k / 7;
            if i < n && !pick.contains(&i) {
                pick.push(i);
            }
        }
        for i in pick {
            calls.push((lam.clone(), cls[i].1.clone(), cls[i].2.clone(), cls[i].3));
        }
    }
    calls.retain(|(l, a, b, d)| !done.contains(&(l.clone(), flows(l, a), flows(l, b), hom_unknowns(l, a, b, *d), *d)));
    eprintln!("[calib] {} calls to measure", calls.len());
    // window classes from the local decomposition planner (labels of the window)
    if let Some(p) = opt(args, "--windows") {
        if let Ok(text) = std::fs::read_to_string(&p) {
            let mut per_cls: std::collections::BTreeMap<Vec<u8>, usize> = Default::default();
            for l in text.lines() {
                let f: Vec<&str> = l.split_whitespace().collect();
                if f.len() != 4 {
                    continue;
                }
                let mut key = parse_lam(f[0]);
                let lab = key.clone();
                key.sort_unstable();
                let c = per_cls.entry(key).or_insert(0);
                if *c < per {
                    *c += 1;
                    calls.push((lab, parse_word(f[1]), parse_word(f[2]), f[3].parse().unwrap()));
                }
            }
        }
    }
    let t0 = Instant::now();
    let mut skipped = 0;
    for (lab, a, b, d) in calls {
        if t0.elapsed().as_secs_f64() > budget_s {
            skipped += 1;
            continue;
        }
        let mut cmd = vec!["calib-one".to_string(), lam_str(&lab), wcode(&a), wcode(&b), d.to_string()];
        if let Some(c) = &cache {
            cmd.push("--cache".into());
            cmd.push(c.clone());
        }
        let (line, w, peak, killed) = supervised_raw(&cmd, limit, 1024.0);
        let fs = flows(&lab, &a);
        let ft = flows(&lab, &b);
        let unk = hom_unknowns(&lab, &a, &b, d);
        match (line, killed) {
            (Some(l), None) => {
                let f: Vec<&str> = l.split_whitespace().collect();
                let secs: f64 = f[5].parse().unwrap();
                let computed: usize = f[7].parse().unwrap();
                if computed == 0 {
                    eprintln!("[calib] {:?} {} -> {} d {}: loaded from cache (not a measurement)", lab, wcode(&a), wcode(&b), d);
                    continue;
                }
                table.record(CostRec { labels: lab.clone(), fs, ft, unknowns: unk, deg: d, secs, rss_mb: peak, censored: false });
                eprintln!("[calib] {:?} {}x{} unknowns {} d {}: {:.3}s {:.0} MB", lab, fs, ft, unk, d, secs, peak);
            }
            _ => {
                table.record(CostRec { labels: lab.clone(), fs, ft, unknowns: unk, deg: d, secs: w, rss_mb: peak, censored: true });
                eprintln!("[calib] {:?} {}x{} unknowns {} d {}: KILLED ({:?}) after {:.1}s", lab, fs, ft, unk, d, killed, w);
            }
        }
        table.save(std::path::Path::new(&out)).unwrap();
    }
    println!("calibration: {} records in {} ({:?}); {} calls skipped by the time budget", table.len(), out, table.classes(), skipped);
}

/// supervised child whose RESULT line is returned raw
fn supervised_raw(cmd: &[String], wall: f64, rss_mb: f64) -> (Option<String>, f64, f64, Option<&'static str>) {
    supervised(cmd, wall, rss_mb)
}
