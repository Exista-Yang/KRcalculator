//! libbuild: the (λ) decomposition library of sindec, built one child process per job under a
//! watchdog (SPEC_3C_GL6 rev1 §2 M5b, D8).
//!
//! libbuild plan   --lam 3,3,3 --n 4|5|6 [--corner markov-ext] [--sample largex [--count 20]] [--json OUT]
//! libbuild run    --lam 3,3,3 --n N [set flags] --field 31|32b --cache DIR [--costs PATH] [--jobs 4]
//!                 [--wall 600] [--rss 2048] [--scale 1] [--json OUT] [--redo] [--max-jobs K] [--no-retry]
//! libbuild replay --lam ... --n N [set flags] --field F --cache DIR [--costs PATH] [--jobs 4] [--json OUT]
//! libbuild calib  --lam 3,3,3 --n N --field 31 --cache DIR [--base PATH] --out PATH [--per-class 12]
//!                 [--limit 120] [--budget 1800] [--rss 4096]
//! libbuild one    --lam ... --n N --field F --cache DIR [--costs PATH] [--scale S] [--replay] JOB...  (child)
//!
//! Sets (D8): without set flags the level set of N (every kept-source pair of `pruned_plan(N)`
//! plus the closure generators of the kept types of (λ, λ)); `--corner markov-ext` the extended
//! Markov corner (σ_1 terms and Rickard components, σ_2 terms both signs, from every kept type
//! p ⊗ id of (λ, λ) whose last strand is split off; closure generators of every kept type
//! reached after σ_2); `--sample largex` the `--count` kept-source pairs of `pruned_plan(N)` with
//! the largest X.
//!
//! Jobs: `pair K` (plan index), `cross P I pos|neg` (crossing σ_{I+1} on the type P), `closure P`.
//! Every job is one child process; the parent polls its RSS (`ps -o rss=`) and wall time and kills
//! it above `--rss` MB or `--wall` s (outcome Memory / Time). Results go to OUT as JSON lines (a
//! header line, one line per job attempt, a summary line); `run` skips jobs already certified in
//! OUT (resume) and every child resumes from layer A (certified targets are never recomputed).
//! Retry policy (G1(c)): a Time / Memory / Exhausted outcome is rerun once at `2 × scale`
//! (watchdog and provider Budget doubled); Time / Memory at 2x stay failures, Exhausted at 2x goes to
//! the M4c list. A replay runs the same jobs with a provider that refuses every call: everything
//! must load from layer A and pass the same checks.
use krhom::field::{Field, F31, F32b};
use krhom::homprov::*;
use krhom::schur::{crossing_terms, flows, Schur, Tid};
use krhom::sindec::*;
use krhom::smap::{Arena, ExprId, MapId};
use krhom::smapio::{labels_str, rung_str, rungs_str};
use krhom::web::Rung;
use std::collections::{BTreeMap, BTreeSet};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::time::Instant;

// ------------------------------------------------------------------ small helpers

fn opt(args: &[String], key: &str) -> Option<String> {
    args.iter().position(|a| a == key).and_then(|i| args.get(i + 1).cloned())
}
fn flag(args: &[String], key: &str) -> bool {
    args.iter().any(|a| a == key)
}
fn parse_lam(s: &str) -> Vec<u8> {
    s.split(',').map(|x| x.parse().unwrap()).collect()
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
fn git_commit() -> String {
    let run = |a: &[&str]| std::process::Command::new("git").args(a).current_dir(env!("CARGO_MANIFEST_DIR")).output().ok().filter(|o| o.status.success()).map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string());
    match run(&["rev-parse", "HEAD"]) {
        Some(h) => {
            let dirty = run(&["status", "--porcelain", "--untracked-files=no"]).map(|x| !x.is_empty()).unwrap_or(false);
            if dirty {
                format!("{}-dirty", h)
            } else {
                h
            }
        }
        None => "unknown".into(),
    }
}

/// a type of the Schur algebra by its table string (searched in the block of its column sums)
fn tid_of(s: &mut Schur, table: &str) -> Tid {
    let rows: Vec<Vec<u8>> = table.split('/').map(|r| r.split(',').map(|x| x.parse().unwrap()).collect()).collect();
    let nc = rows[0].len();
    let mu: Vec<u8> = (0..nc).map(|j| rows.iter().map(|r| r[j]).sum()).collect();
    let order = s.canonical_basis(&mu).order.clone();
    for q in order {
        if s.table(q).fmt() == table {
            return q;
        }
    }
    panic!("type {} not found in block {:?}", table, mu)
}

/// p = p' ⊗ id: the last strand of p is split off (entry (n-1, n-1) of the table = λ_{n-1})
fn last_strand_split(s: &Schur, p: Tid) -> bool {
    let t = s.table(p);
    let (nr, nc) = (t.nr as usize, t.nc as usize);
    let rows = t.rows();
    t.e[(nr - 1) * nc + (nc - 1)] == rows[nr - 1] && t.cols()[nc - 1] == rows[nr - 1]
}

// ------------------------------------------------------------------ job lists

#[derive(Clone, Debug)]
struct Job {
    /// the child's job words (pair K | cross P I pos|neg | closure P)
    spec: Vec<String>,
    /// expected kept provider targets and Σ mult (pair jobs; closure = 2 sides; cross: unknown)
    want_targets: Option<usize>,
    want_mult: Option<usize>,
    x_flows: u128,
}
impl Job {
    fn key(&self) -> String {
        self.spec.join(" ")
    }
}

struct SetDef {
    name: String,
    jobs: Vec<Job>,
}

fn set_name(args: &[String], n: u8) -> String {
    if opt(args, "--corner").is_some() {
        format!("corner-{}-n{}", opt(args, "--corner").unwrap(), n)
    } else if opt(args, "--sample").is_some() {
        format!("sample-{}{}-n{}", opt(args, "--sample").unwrap(), opt(args, "--count").unwrap_or_else(|| "20".into()), n)
    } else {
        format!("level-n{}", n)
    }
}

fn pair_job(s: &mut Schur, lam: &[u8], plan: &Plan, k: usize) -> Job {
    let pp = &plan.pairs[k];
    let mut xw = canonical_ladder(s, pp.p);
    xw.push(pp.rung);
    let trivial = pp.identity || pp.all_killed;
    Job {
        spec: vec!["pair".into(), k.to_string()],
        want_targets: Some(if trivial { 0 } else { pp.kept.len() }),
        want_mult: Some(if trivial { 0 } else { pp.kept.iter().map(|t| t.mult).sum() }),
        x_flows: flows(lam, &xw),
    }
}
fn closure_job(s: &mut Schur, lam: &[u8], p: Tid) -> Job {
    let lp = canonical_ladder(s, p);
    Job { spec: vec!["closure".into(), s.table(p).fmt()], want_targets: Some(2), want_mult: None, x_flows: flows(lam, &lp) }
}

fn build_set(args: &[String], s: &mut Schur, lam: &[u8], n: u8) -> (SetDef, Plan) {
    let plan = pruned_plan(s, Some(n));
    let mut jobs = vec![];
    if let Some(c) = opt(args, "--corner") {
        assert_eq!(c, "markov-ext", "unknown corner {}", c);
        let order = s.canonical_basis(lam).order.clone();
        let srcs: Vec<Tid> = order.into_iter().filter(|&p| last_strand_split(s, p) && !s.killed_at(p, n as usize)).collect();
        let mut reached: BTreeSet<Tid> = BTreeSet::new();
        for i in 0..lam.len() - 1 {
            for &p in &srcs {
                for pos in [true, false] {
                    let lp = canonical_ladder(s, p);
                    let mut xmax = 0;
                    for t in crossing_terms(lam[i], lam[i + 1], pos, i) {
                        let mut w = lp.clone();
                        w.extend(&t.rungs);
                        xmax = xmax.max(flows(lam, &w));
                        if i == lam.len() - 2 {
                            for (q, _) in s.term_mult(p, &t.rungs) {
                                if s.table(q).cols() == lam && !s.killed_at(q, n as usize) {
                                    reached.insert(q);
                                }
                            }
                        }
                    }
                    jobs.push(Job { spec: vec!["cross".into(), s.table(p).fmt(), i.to_string(), if pos { "pos" } else { "neg" }.into()], want_targets: None, want_mult: None, x_flows: xmax });
                }
            }
        }
        for q in reached {
            jobs.push(closure_job(s, lam, q));
        }
    } else if let Some(sm) = opt(args, "--sample") {
        assert_eq!(sm, "largex", "unknown sample {}", sm);
        let count: usize = opt(args, "--count").map(|v| v.parse().unwrap()).unwrap_or(20);
        let mut cand: Vec<Job> = (0..plan.pairs.len()).filter(|&k| !plan.pairs[k].p_killed && !plan.pairs[k].identity && !plan.pairs[k].all_killed).map(|k| pair_job(s, lam, &plan, k)).collect();
        cand.sort_by(|a, b| b.x_flows.cmp(&a.x_flows).then(a.spec[1].parse::<usize>().unwrap().cmp(&b.spec[1].parse::<usize>().unwrap())));
        cand.truncate(count);
        jobs = cand;
    } else {
        for k in 0..plan.pairs.len() {
            if !plan.pairs[k].p_killed {
                jobs.push(pair_job(s, lam, &plan, k));
            }
        }
        let order = s.canonical_basis(lam).order.clone();
        for p in order {
            if !s.killed_at(p, n as usize) && !canonical_ladder(s, p).is_empty() {
                jobs.push(closure_job(s, lam, p));
            }
        }
    }
    (SetDef { name: set_name(args, n), jobs }, plan)
}

fn cmd_plan(args: &[String]) {
    let lam = parse_lam(&opt(args, "--lam").expect("--lam"));
    let n: u8 = opt(args, "--n").expect("--n").parse().unwrap();
    let mut s = Schur::new(&lam);
    let (set, plan) = build_set(args, &mut s, &lam, n);
    println!("{}", plan.summary());
    let wt: usize = set.jobs.iter().filter_map(|j| j.want_targets).sum();
    let wm: usize = set.jobs.iter().filter_map(|j| j.want_mult).sum();
    let mx = set.jobs.iter().map(|j| j.x_flows).max().unwrap_or(0);
    let kinds: BTreeMap<String, usize> = set.jobs.iter().fold(BTreeMap::new(), |mut m, j| {
        *m.entry(j.spec[0].clone()).or_insert(0) += 1;
        m
    });
    println!("set {}: {} jobs {:?}; known targets {} (Σ mult {}); max X {} flows", set.name, set.jobs.len(), kinds, wt, wm, mx);
    if flag(args, "--list") {
        for j in &set.jobs {
            println!("  {:<40} X {:>6} targets {:?} mult {:?}", j.key(), j.x_flows, j.want_targets, j.want_mult);
        }
    }
    if let Some(p) = opt(args, "--json") {
        let recs: Vec<String> = set.jobs.iter().map(|j| format!("{{\"job\":{},\"x_flows\":{},\"want_targets\":{},\"want_mult\":{}}}", jstr(&j.key()), j.x_flows, j.want_targets.map_or("null".into(), |v| v.to_string()), j.want_mult.map_or("null".into(), |v| v.to_string()))).collect();
        std::fs::write(p, format!("{{\"set\":{},\"plan\":{},\"jobs\":[\n{}\n]}}\n", jstr(&set.name), jstr(&plan.summary()), recs.join(",\n"))).unwrap();
    }
}

// ------------------------------------------------------------------ child

/// a provider that refuses every call (replay: everything must come from layer A)
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

fn scaled_budget(scale: f64) -> Budget {
    let mut b = Budget::default();
    b.t_max_s *= scale;
    b.m_max_bytes = (b.m_max_bytes as f64 * scale) as u64;
    b.index_max_bytes = (b.index_max_bytes as f64 * scale) as u64;
    b.c_max = (b.c_max as f64 * scale) as usize;
    b
}

/// outcome class of a sindec error: Exhausted | Time | Memory | Refused | Check | Other
fn err_kind(e: &DecErr) -> &'static str {
    match e {
        DecErr::Provider { engineering, text } => {
            if !engineering {
                "Exhausted"
            } else if text.contains("Time") {
                "Time"
            } else if text.contains("Memory") {
                "Memory"
            } else {
                "Refused"
            }
        }
        DecErr::Refused(_) => "Refused",
        DecErr::Check(_) | DecErr::Cache(_) => "Check",
        DecErr::NeedsIdempotent { .. } => "NeedsIdempotent",
        _ => "Other",
    }
}

fn provider_and_stamp<F: Field>(lam: &[u8], args: &[String]) -> (Box<dyn HomProvider<F>>, String) {
    let costs = opt(args, "--costs").map(PathBuf::from);
    provider_for_stamped::<F>(lam, costs.as_deref(), 1.0)
}

fn child<F: Field>(args: &[String], job: &[String]) {
    let lam = parse_lam(&opt(args, "--lam").expect("--lam"));
    let n: u8 = opt(args, "--n").expect("--n").parse().unwrap();
    let cache = PathBuf::from(opt(args, "--cache").expect("--cache"));
    let scale: f64 = opt(args, "--scale").map(|v| v.parse().unwrap()).unwrap_or(1.0);
    let replay = flag(args, "--replay");
    let (prov, stamp) = provider_and_stamp::<F>(&lam, args);
    let prov: Box<dyn HomProvider<F>> = if replay { Box::new(NoProvider) } else { prov };
    let opts = SIndecOpts { cache_dir: Some(cache.clone()), budget: scaled_budget(scale), prune_n: Some(n), ..SIndecOpts::default() };
    let mut sd: SIndec<F> = SIndec::with_stamp(&lam, prov, stamp, opts);
    sd.target_log_file = opt(args, "--tlog").map(PathBuf::from);
    let t0 = Instant::now();
    let mut errs: Vec<(String, &'static str, String)> = vec![];
    let mut info = String::new();
    match job[0].as_str() {
        "pair" => {
            let k: usize = job[1].parse().unwrap();
            let plan = pruned_plan(&mut sd.schur, Some(n));
            let pp = plan.pairs[k].clone();
            info = format!("\"p\":{},\"rung\":{},", jstr(&sd.schur.table(pp.p).fmt()), jstr(&rung_str(&pp.rung)));
            if let Err(e) = sd.dec_rung(pp.p, pp.rung) {
                errs.push((format!("dec {}", job.join(" ")), err_kind(&e), e.to_string()));
            }
        }
        "cross" => {
            let p = tid_of(&mut sd.schur, &job[1]);
            let i: usize = job[2].parse().unwrap();
            let pos = job[3] == "pos";
            let terms = crossing_terms(lam[i], lam[i + 1], pos, i);
            let mut nsum = vec![];
            for (k, t) in terms.iter().enumerate() {
                match sd.dec_word(p, &t.rungs) {
                    Ok(d) => nsum.push(format!("{}", d.summands.len())),
                    Err(e) => {
                        nsum.push("err".into());
                        errs.push((format!("dec W_{} {}", k, rungs_str(&t.rungs)), err_kind(&e), e.to_string()));
                    }
                }
            }
            let mut nrick = 0;
            if i == 0 && errs.is_empty() {
                for j in 0..terms.len() - 1 {
                    match sd.rickard_comp_at(p, i, pos, j) {
                        Ok(r) => nrick += r.len(),
                        Err(DecErr::Killed { .. }) => {}
                        Err(e) => errs.push((format!("rickard d_{}", j), err_kind(&e), e.to_string())),
                    }
                }
            }
            info = format!("\"summands_per_term\":[{}],\"rickard_components\":{},", nsum.join(","), nrick);
        }
        "closure" => {
            let p = tid_of(&mut sd.schur, &job[1]);
            match sd.closure_gens(p) {
                Ok(g) => info = format!("\"ins\":{},\"outs\":{},", g.ins.len(), g.outs.len()),
                Err(e) => errs.push(("closure".into(), err_kind(&e), e.to_string())),
            }
        }
        other => panic!("unknown job {}", other),
    }
    let st = sd.stats();
    let by = sd.ar.bytes();
    let tl: Vec<String> = sd.target_log.iter().map(|e| format!("{{\"p\":{},\"q\":{},\"s\":{},\"mult\":{},\"secs\":{:.3},\"summary\":{}}}", jstr(&e.0), jstr(&e.1), e.2, e.3, e.4, jstr(&e.5))).collect();
    let el: Vec<String> = errs.iter().map(|(w, k, t)| format!("{{\"what\":{},\"kind\":{},\"text\":{}}}", jstr(w), jstr(k), jstr(&t[..t.len().min(2000)]))).collect();
    let _ = sd.ar.lib_ref().save();
    println!(
        "RESULT {{\"ok\":{},{}\"wall\":{:.3},\"errors\":[{}],\"stats\":{{\"decomps\":{},\"identity\":{},\"all_killed\":{},\"targets\":{},\"sum_mult\":{},\"provider_calls\":{},\"provider_secs\":{:.3},\"max_x_flows\":{},\"layer_a_loaded\":{},\"layer_a_saved\":{},\"layer_a_rejected\":{},\"rickard\":{},\"closure\":{},\"scaling_entries\":{},\"scaling_zero\":{},\"check_secs\":{:.3}}},\"arena_bytes\":{},\"fresh\":[{}]}}",
        errs.is_empty(), info, t0.elapsed().as_secs_f64(), el.join(","), st.decomps, st.identity, st.all_killed, st.targets, st.sum_mult, st.provider_calls, st.provider_secs, st.max_x_flows, st.layer_a_loaded, st.layer_a_saved, st.layer_a_rejected, st.rickard, st.closure, st.scaling_entries, st.scaling_zero, st.check_secs, by.total(), tl.join(",")
    );
}

// ------------------------------------------------------------------ watchdog and parent

/// run a child under the watchdog: (stdout RESULT line or None, wall, peak RSS MB, killed reason)
fn supervised(cmd: &[String], wall: f64, rss_mb: f64, threads: usize) -> (Option<String>, f64, f64, Option<&'static str>) {
    let exe = std::env::current_exe().unwrap();
    let mut child = std::process::Command::new(exe).args(cmd).env("KRHOM_THREADS", threads.to_string()).stdout(std::process::Stdio::piped()).stderr(std::process::Stdio::inherit()).spawn().expect("spawn child");
    let pid = child.id();
    let t0 = Instant::now();
    let mut peak = 0.0f64;
    let mut killed = None;
    // drain stdout in a thread so that a chatty child never blocks on a full pipe
    let mut so = child.stdout.take().unwrap();
    let reader = std::thread::spawn(move || {
        let mut s = String::new();
        let _ = std::io::Read::read_to_string(&mut so, &mut s);
        s
    });
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
    let _ = child.wait();
    let out = reader.join().unwrap_or_default();
    let line = out.lines().find(|l| l.starts_with("RESULT ")).map(|l| l[7..].to_string());
    (if killed.is_some() { None } else { line }, t0.elapsed().as_secs_f64(), peak, killed)
}

/// minimal field extraction from our own flat JSON (`"key":value` with value up to , or })
fn jfield<'a>(l: &'a str, key: &str) -> Option<&'a str> {
    let pat = format!("\"{}\":", key);
    let i = l.find(&pat)? + pat.len();
    let rest = &l[i..];
    if let Some(r) = rest.strip_prefix('"') {
        let mut esc = false;
        for (k, c) in r.char_indices() {
            if esc {
                esc = false;
            } else if c == '\\' {
                esc = true;
            } else if c == '"' {
                return Some(&r[..k]);
            }
        }
        None
    } else {
        let end = rest.find([',', '}', ']']).unwrap_or(rest.len());
        Some(rest[..end].trim())
    }
}

#[derive(Clone, Debug)]
struct Attempt {
    job: String,
    scale: f64,
    outcome: String,
    rec: String,
}

/// outcome of one attempt: ok | Time | Memory | Exhausted | Refused | Check | crash | ...
fn outcome_of(line: &Option<String>, killed: Option<&str>) -> String {
    match (line, killed) {
        (_, Some(k)) => k.to_string(),
        (None, None) => "crash".into(),
        (Some(l), None) => {
            if l.contains("\"ok\":true") {
                "ok".into()
            } else {
                // worst kind among the errors: Exhausted dominates engineering outcomes
                let kinds: Vec<&str> = l.match_indices("\"kind\":\"").map(|(i, _)| &l[i + 8..]).map(|r| &r[..r.find('"').unwrap_or(0)]).collect();
                for k in ["Check", "Other", "NeedsIdempotent", "Exhausted", "Refused", "Memory", "Time"] {
                    if kinds.contains(&k) {
                        return k.to_string();
                    }
                }
                "failed".into()
            }
        }
    }
}

fn field_of(args: &[String]) -> String {
    opt(args, "--field").unwrap_or_else(|| "31".into())
}

fn cmd_run(args: &[String], replay: bool) {
    let lam = parse_lam(&opt(args, "--lam").expect("--lam"));
    let n: u8 = opt(args, "--n").expect("--n").parse().unwrap();
    let field = field_of(args);
    assert!(field == "31" || field == "32b", "--field 31|32b");
    let cache = opt(args, "--cache").expect("--cache");
    let jobs_n: usize = opt(args, "--jobs").map(|v| v.parse().unwrap()).unwrap_or(4);
    let wall: f64 = opt(args, "--wall").map(|v| v.parse().unwrap()).unwrap_or(600.0);
    let rss: f64 = opt(args, "--rss").map(|v| v.parse().unwrap()).unwrap_or(2048.0);
    let scale0: f64 = opt(args, "--scale").map(|v| v.parse().unwrap()).unwrap_or(1.0);
    let retry = !flag(args, "--no-retry") && !replay;
    let mut s = Schur::new(&lam);
    let (set, _plan) = build_set(args, &mut s, &lam, n);
    let out = opt(args, "--json").unwrap_or_else(|| format!("libbuild_{}_{}_f{}.jsonl", labels_str(&lam).replace(',', ""), set.name, field));
    // resume: jobs already certified in OUT (same mode) are skipped
    let mut done: BTreeSet<String> = BTreeSet::new();
    let mut prev: Vec<String> = vec![];
    if !flag(args, "--redo") {
        if let Ok(t) = std::fs::read_to_string(&out) {
            for l in t.lines() {
                if l.starts_with("{\"attempt\"") {
                    if l.contains("\"outcome\":\"ok\"") {
                        if let Some(j) = jfield(l, "job") {
                            done.insert(j.to_string());
                        }
                    }
                    prev.push(l.to_string());
                }
            }
        }
    }
    let mut todo: Vec<Job> = set.jobs.iter().filter(|j| !done.contains(&j.key())).cloned().collect();
    if let Some(m) = opt(args, "--max-jobs") {
        todo.truncate(m.parse().unwrap());
    }
    let restricted = opt(args, "--max-jobs").is_some();
    let stamp = match field.as_str() {
        "31" => provider_and_stamp::<F31>(&lam, args).1,
        _ => provider_and_stamp::<F32b>(&lam, args).1,
    };
    let argv: Vec<String> = args.iter().map(|a| jstr(a)).collect();
    let header = format!(
        "{{\"header\":{{\"tool\":\"libbuild\",\"mode\":{},\"args\":[{}],\"set\":{},\"lam\":{},\"n\":{},\"field\":{},\"commit\":{},\"host\":{},\"stamps\":{},\"watchdog\":{{\"wall_s\":{},\"rss_mb\":{},\"jobs\":{}}},\"scale\":{},\"jobs_total\":{},\"jobs_skipped_done\":{},\"restricted\":{}}}}}",
        jstr(if replay { "replay" } else { "run" }), argv.join(","), jstr(&set.name), jstr(&labels_str(&lam)), n, jstr(&field), jstr(&git_commit()), jstr(&host_stamp()),
        jstr(&format!("{} ladfn={} tmpl_alg={} prov={}", SINDEC_VERSION, LADDER_STAMP, krhom::tmpl::TMPL_ALG, stamp)), wall, rss, jobs_n, scale0, set.jobs.len(), done.len(), restricted
    );
    eprintln!("[libbuild] {} {} F{}: {} jobs, {} already certified, {} to run, {} workers", if replay { "replay" } else { "run" }, set.name, field, set.jobs.len(), done.len(), todo.len(), jobs_n);
    let fh = Arc::new(Mutex::new({
        let mut f = std::fs::OpenOptions::new().create(true).append(true).open(&out).unwrap();
        writeln!(f, "{}", header).unwrap();
        f
    }));
    let threads = (4 / jobs_n).max(1);
    let queue: Arc<Mutex<Vec<(Job, f64)>>> = Arc::new(Mutex::new(todo.iter().rev().map(|j| (j.clone(), scale0)).collect()));
    let attempts: Arc<Mutex<Vec<Attempt>>> = Arc::new(Mutex::new(vec![]));
    let t_all = Instant::now();
    let mut base: Vec<String> = vec!["one".into()];
    for k in ["--lam", "--n", "--field", "--cache", "--costs", "--corner", "--sample", "--count"] {
        if let Some(v) = opt(args, k) {
            base.push(k.into());
            base.push(v);
        }
    }
    if replay {
        base.push("--replay".into());
    }
    let tlog = format!("{}.targets", out);
    let mut hs = vec![];
    for _w in 0..jobs_n {
        let (queue, attempts, fh, base, tlog) = (queue.clone(), attempts.clone(), fh.clone(), base.clone(), tlog.clone());
        hs.push(std::thread::spawn(move || loop {
            let item = queue.lock().unwrap().pop();
            let Some((job, sc)) = item else { break };
            let mut cmd = base.clone();
            cmd.extend(["--scale".to_string(), sc.to_string(), "--tlog".to_string(), tlog.clone(), "--".to_string()]);
            cmd.extend(job.spec.iter().cloned());
            let (line, w, peak, killed) = supervised(&cmd, wall * sc, rss * sc, threads);
            let outcome = outcome_of(&line, killed);
            // a pair job certifies exactly its plan targets (counted fresh + loaded)
            let mut outcome = outcome;
            if outcome == "ok" {
                if let (Some(l), Some(wt)) = (&line, job.want_targets) {
                    if job.spec[0] == "pair" {
                        let got: usize = jfield(l, "targets").and_then(|v| v.parse().ok()).unwrap_or(usize::MAX);
                        let gm: usize = jfield(l, "sum_mult").and_then(|v| v.parse().ok()).unwrap_or(usize::MAX);
                        if got != wt || Some(gm) != job.want_mult {
                            outcome = "Count".into();
                        }
                    }
                }
            }
            let rec = format!(
                "{{\"attempt\":1,\"job\":{},\"scale\":{},\"outcome\":{},\"x_flows\":{},\"want_targets\":{},\"want_mult\":{},\"wall\":{:.2},\"peak_rss_mb\":{:.1},\"result\":{}}}",
                jstr(&job.key()), sc, jstr(&outcome), job.x_flows, job.want_targets.map_or("null".into(), |v| v.to_string()), job.want_mult.map_or("null".into(), |v| v.to_string()), w, peak, line.clone().unwrap_or_else(|| "null".into())
            );
            writeln!(fh.lock().unwrap(), "{}", rec).unwrap();
            eprintln!("[libbuild] {} x{} -> {} in {:.1}s (peak {:.0} MB, X {})", job.key(), sc, outcome, w, peak, job.x_flows);
            let again = retry && sc < 2.0 * scale0 && ["Time", "Memory", "Exhausted"].contains(&outcome.as_str());
            attempts.lock().unwrap().push(Attempt { job: job.key(), scale: sc, outcome, rec });
            if again {
                queue.lock().unwrap().insert(0, (job, 2.0 * sc));
            }
        }));
    }
    for h in hs {
        h.join().unwrap();
    }
    // summary over this run and the certified attempts of earlier runs
    let att = attempts.lock().unwrap().clone();
    let mut last: BTreeMap<String, Attempt> = BTreeMap::new();
    for a in att {
        last.insert(a.job.clone(), a);
    }
    let mut n_ok = done.len();
    let mut fails: Vec<String> = vec![];
    let mut m4c: Vec<String> = vec![];
    for j in &todo {
        match last.get(&j.key()) {
            Some(a) if a.outcome == "ok" => n_ok += 1,
            Some(a) => {
                fails.push(format!("{} ({} at x{})", a.job, a.outcome, a.scale));
                if a.outcome == "Exhausted" && a.scale >= 2.0 * scale0 {
                    m4c.push(a.job.clone());
                }
            }
            None => fails.push(format!("{} (not run)", j.key())),
        }
    }
    let total = t_all.elapsed().as_secs_f64();
    let summary = format!(
        "{{\"summary\":{{\"set\":{},\"field\":{},\"mode\":{},\"jobs\":{},\"certified\":{},\"failed\":[{}],\"m4c\":[{}],\"wall\":{:.1},\"complete\":{}}}}}",
        jstr(&set.name), jstr(&field), jstr(if replay { "replay" } else { "run" }), set.jobs.len(), n_ok, fails.iter().map(|f| jstr(f)).collect::<Vec<_>>().join(","), m4c.iter().map(|f| jstr(f)).collect::<Vec<_>>().join(","), total, !restricted && n_ok == set.jobs.len()
    );
    writeln!(fh.lock().unwrap(), "{}", summary).unwrap();
    println!("{} F{} {}: jobs certified {}/{} ({:.1}%), failures {}, M4c {}, wall {:.1}s", set.name, field, if replay { "replay" } else { "run" }, n_ok, set.jobs.len(), 100.0 * n_ok as f64 / set.jobs.len().max(1) as f64, fails.len(), m4c.len(), total);
    for f in &fails {
        println!("  FAIL {}", f);
    }
    let _ = prev;
}

// ------------------------------------------------------------------ calibration

/// child: one Exact whole-ladder hom_space call, prints "RESULT labels fs ft unknowns d secs dim computed"
fn calib_one<F: Field>(args: &[String], job: &[String]) {
    let lam = parse_lam(&opt(args, "--lam").expect("--lam"));
    let a = krhom::smapio::parse_rungs(&job[0]).unwrap();
    let b = krhom::smapio::parse_rungs(&job[1]).unwrap();
    let d: i32 = job[2].parse().unwrap();
    let cache = opt(args, "--cache").map(|c| PathBuf::from(c).join("tmpl"));
    let mut lib: krhom::tmpl::TmplLib<F> = krhom::tmpl::TmplLib::new(cache);
    let t0 = Instant::now();
    let ids = lib.try_hom_basis_cap(&krhom::tmpl::Window::new(&lam, &a), &krhom::tmpl::Window::new(&lam, &b), d, u128::MAX).unwrap();
    let secs = t0.elapsed().as_secs_f64();
    let computed = lib.stats().computed;
    let _ = lib.save();
    println!("RESULT {} {} {} {} {} {:.6} {} {}", labels_str(&lam), flows(&lam, &a), flows(&lam, &b), hom_unknowns(&lam, &a, &b, d), d, secs, ids.len(), computed);
}

/// Cost records of the whole-ladder (λ) Exact calls of the plan targets that the label rule does
/// not admit (fibre product > cap or unknowns > cap): Hom^s(L_q, X) and Hom^{-s}(X, L_q), smallest
/// unknowns first, `--per-class` per sorted degree class spread over the range; each call in a child
/// under the watchdog (a kill is a censored record at the time limit). Output: `--base` records plus
/// the new ones, stamped v2 for this host.
fn cmd_calib(args: &[String]) {
    let lam = parse_lam(&opt(args, "--lam").expect("--lam"));
    let n: u8 = opt(args, "--n").expect("--n").parse().unwrap();
    let per: usize = opt(args, "--per-class").map(|v| v.parse().unwrap()).unwrap_or(12);
    let limit: f64 = opt(args, "--limit").map(|v| v.parse().unwrap()).unwrap_or(120.0);
    let budget_s: f64 = opt(args, "--budget").map(|v| v.parse().unwrap()).unwrap_or(1800.0);
    let rss: f64 = opt(args, "--rss").map(|v| v.parse().unwrap()).unwrap_or(4096.0);
    let out = opt(args, "--out").expect("--out");
    let cache = opt(args, "--cache").expect("--cache");
    let mut table = CostTable::new();
    if let Ok(t) = CostTable::load(Path::new(&out)) {
        for r in t.recs {
            table.record(r);
        }
    } else if let Some(b) = opt(args, "--base") {
        unsafe { std::env::set_var("KRHOM_COSTS_ANY_HOST", "1") };
        let t = CostTable::load(Path::new(&b)).expect("base cost table");
        for r in t.recs {
            table.record(r);
        }
    }
    let ex = Exact { max_label: lam.iter().copied().max().unwrap(), ..Exact::default() };
    let done: BTreeSet<(u128, u128, u128, i32)> = table.recs.iter().filter(|r| r.labels == lam).map(|r| (r.fs, r.ft, r.unknowns, r.deg)).collect();
    let mut s = Schur::new(&lam);
    let plan = pruned_plan(&mut s, Some(n));
    let mut cls: Vec<(u128, Vec<Rung>, Vec<Rung>, i32)> = vec![];
    for pp in plan.pairs.iter().filter(|p| !p.p_killed && !p.identity && !p.all_killed) {
        let mut xw = canonical_ladder(&mut s, pp.p);
        xw.push(pp.rung);
        for t in &pp.kept {
            let lq = canonical_ladder(&mut s, t.q);
            for (a, b, d) in [(lq.clone(), xw.clone(), t.s), (xw.clone(), lq.clone(), -t.s)] {
                let (fs, ft, unk) = (flows(&lam, &a), flows(&lam, &b), hom_unknowns(&lam, &a, &b, d));
                if ex.admit(&lam, fs, ft, unk, d).is_err() && !done.contains(&(fs, ft, unk, d)) {
                    cls.push((unk, a, b, d));
                }
            }
        }
    }
    cls.sort_by(|a, b| a.0.cmp(&b.0));
    cls.dedup_by(|a, b| a.1 == b.1 && a.2 == b.2 && a.3 == b.3);
    let m = cls.len();
    let mut pick: Vec<usize> = (0..per.min(m)).collect();
    for k in 1..=per {
        let i = per + m.saturating_sub(per) * k / (per + 1);
        if i < m && !pick.contains(&i) {
            pick.push(i);
        }
    }
    eprintln!("[calib] {} unadmitted (λ) calls at N = {}; measuring {}", m, n, pick.len());
    let t0 = Instant::now();
    let mut skipped = 0;
    let field = field_of(args);
    for i in pick {
        if t0.elapsed().as_secs_f64() > budget_s {
            skipped += 1;
            continue;
        }
        let (_, a, b, d) = &cls[i];
        let cmd = vec!["calib-one".to_string(), "--lam".into(), labels_str(&lam), "--field".into(), field.clone(), "--cache".into(), cache.clone(), "--".into(), rungs_str(a), rungs_str(b), d.to_string()];
        let (line, w, peak, killed) = supervised(&cmd, limit, rss, 1);
        let (fs, ft, unk) = (flows(&lam, a), flows(&lam, b), hom_unknowns(&lam, a, b, *d));
        match (line, killed) {
            (Some(l), None) => {
                let f: Vec<&str> = l.split_whitespace().collect();
                let secs: f64 = f[5].parse().unwrap();
                let computed: usize = f[7].parse().unwrap();
                if computed == 0 {
                    eprintln!("[calib] {}x{} d {}: loaded from cache (not a measurement)", fs, ft, d);
                    continue;
                }
                table.record(CostRec { labels: lam.clone(), fs, ft, unknowns: unk, deg: *d, secs, rss_mb: peak, censored: false });
                eprintln!("[calib] {}x{} unknowns {} d {}: {:.3}s {:.0} MB (dim {})", fs, ft, unk, d, secs, peak, f[6]);
            }
            (_, k) => {
                table.record(CostRec { labels: lam.clone(), fs, ft, unknowns: unk, deg: *d, secs: w, rss_mb: peak, censored: true });
                eprintln!("[calib] {}x{} unknowns {} d {}: KILLED ({:?}) after {:.1}s, {:.0} MB", fs, ft, unk, d, k, w, peak);
            }
        }
        table.save(Path::new(&out)).unwrap();
    }
    println!("calibration: {} records in {} ({:?}); {} calls skipped by the time budget", table.len(), out, table.classes(), skipped);
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let cmd = args.get(1).cloned().unwrap_or_default();
    // child job words follow "--"
    let job: Vec<String> = args.iter().position(|a| a == "--").map(|i| args[i + 1..].to_vec()).unwrap_or_default();
    let field = field_of(&args);
    match cmd.as_str() {
        "plan" => cmd_plan(&args),
        "run" => cmd_run(&args, false),
        "replay" => cmd_run(&args, true),
        "calib" => cmd_calib(&args),
        "one" if field == "31" => child::<F31>(&args, &job),
        "one" => child::<F32b>(&args, &job),
        "calib-one" if field == "31" => calib_one::<F31>(&args, &job),
        "calib-one" => calib_one::<F32b>(&args, &job),
        _ => {
            eprintln!("usage: libbuild plan|run|replay|calib --lam 3,3,3 --n 4|5|6 [--corner markov-ext] [--sample largex] --field 31|32b --cache DIR --jobs 4 --wall 600 --rss 2048 --json OUT");
            std::process::exit(2);
        }
    }
}
