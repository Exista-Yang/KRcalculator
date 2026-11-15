//! krhom: colored Khovanov–Rozansky gl_N homology of links, as a one-command calculator.
//!
//! usage (see `krhom --help` and docs/CALCULATOR.md):
//!   krhom N[,N...] --braid WORD [--strands n] [--colors c | c1,c2,...] [--reduced]
//!   krhom N[,N...] --pd PD [--colors ...] [--reduced]
//!   krhom N[,N...] --file BRAIDS [--colors ...] [--reduced] [--csv OUT]
//! WORD is signed integers ("1 -2 1 -2") or KnotJob letters ("aBaB": a..i = σ1..σ9, A..I their
//! inverses); --strands defaults to max |generator| + 1. A single color applies to every
//! component. Every result is over Q: the randomized engines run over F_p for two large primes
//! and must agree (certified over F_p for both primes; equal over Q unless both are torsion
//! primes); the full cube is exact over Q.
//!
//! Routes (chosen automatically; `--colored-method cable|cube|singular` overrides):
//! * all colors 1: the localized Soergel engine (engine2), any number of strands, two primes;
//! * colored, at most 2 strands after simplification: the exact singular Soergel engine
//!   (xengine), two primes;
//! * colored, 3-4 strands, uniform color 2, N <= 4: the localized box engine with the prebuilt
//!   layer-A library (src/calc.rs; tools/fetch_library.sh), two primes;
//! * anything else colored: exit 1 with the support matrix (`--colored-method cube`: exact over Q,
//!   tiny inputs only).
//! Exit codes: 0 ok, 1 usage/unsupported, 2 engine/certification failure, 3 primes disagree,
//! 4 library missing/damaged/not covering the word, 5 resource refusal.
use krhom::calc::*;
use krhom::reference::cabling::{cable_default, khr_colored_cable};
use krhom::input::crossing::{Gen, QConv};
use krhom::reference::cube::{khr_full_cube, BraidLink};
use krhom::uncolored::engine2::Closure;
use krhom::algebra::field::{F31, F32b, Q};
use krhom::algebra::laurent::Poincare;
use krhom::input::pd::{parse_pd, pd_to_braid, BraidWord};
use krhom::uncolored::tower::EndPolicy;
use krhom::colored::xengine::khr_singular_exact;
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::time::Instant;

const USAGE: &str = "usage: krhom N[,N...] (--braid WORD [--strands n] | --pd PD | --file BRAIDS) [--colors c | c1,c2,...] [--reduced]
             [--table] [--json OUT] [--csv OUT] [--no-simplify] [--check-library] [--lib-dir DIR]
             [--seed S] [--seed-retries K] [--mem-mb M] [--mem-cap-mb M]
             [--colored-method cable|cube|singular] [--cube] [--closure tower|exact] [--end right|left|alt]
  N        gl_N ranks, a comma-separated list (e.g. 2,3,4)
  --braid  signed integers \"1 -2 1 -2\" or KnotJob letters \"aBaB\" (a..i = sigma_1..sigma_9, A..I inverses)
  --strands  default max |generator| + 1
  --colors c applies to every component; c1,c2,... is one color per link component (default 1)
  --reduced  reduced homology (the strand at bottom position 0 is marked)
  --file   one braid per line, optional \"NAME;\" prefix, '#' comments; same N/colors/--reduced for all
  --table  also print the bigraded table of each N (rows t, columns q, blank = 0)
  --json OUT  machine-readable results (an array with --file); --csv OUT rows name,N,reduced,t,q,rank
  colored 3-4 strands (uniform color 2, N <= 4) need the prebuilt library: tools/fetch_library.sh l222_n4 l2222_n4
  results are over Q (two large primes must agree; the full cube is exact over Q); docs/CALCULATOR.md";

const VALUE_OPTS: [&str; 18] = ["--pd", "--braid", "--strands", "--colors", "--colored-method", "--field", "--closure", "--end", "--json", "--csv", "--file", "--lib-dir", "--seed", "--seed-f32b", "--seed-retries", "--mem-mb", "--mem-cap-mb", "--name"];
const FLAG_OPTS: [&str; 9] = ["--table", "--cube", "--reduced", "--single-prime", "--no-simplify", "--check-library", "--box-live", "--help", "-h"];

#[derive(Clone, Debug)]
struct Opts {
    ns: Vec<usize>,
    colors: Option<Vec<usize>>,
    strands: Option<usize>,
    method: Option<String>,
    cube: bool,
    reduced: bool,
    single_prime: bool,
    no_simplify: bool,
    check_library: bool,
    box_live: bool,
    closure: Closure,
    lib_dir: Option<String>,
    seed: u64,
    seed2: u64,
    retries: u64,
    mem_mb: Option<u64>,
    mem_cap_mb: Option<u64>,
    /// KRHOM_MEM_CAP_MB as given by the user's environment at start
    env_mem_cap: bool,
    env_any_host: bool,
    /// --table: print the bigraded table after each N block
    table: bool,
}

enum Input {
    Braid(String),
    Pd(String),
}
struct Item {
    name: Option<String>,
    input: Input,
}

/// What one item produced.
#[derive(Default)]
struct Outcome {
    route: &'static str,
    used: Option<BraidWord>,
    orig: Option<BraidWord>,
    colors: Vec<usize>,
    labels: Vec<usize>,
    results: BTreeMap<usize, Poincare>,
    zero_ns: Vec<usize>,
    primes: Vec<PrimeOut>,
    seeds: Vec<u64>,
    agree: Option<bool>,
    library: Option<String>,
    coefficients: String,
    check_library: Option<String>,
    mem_cap_mb: Option<u64>,
    mem_mb: Option<u64>,
}

fn fail(code: i32, msg: &str) -> ! {
    eprintln!("krhom: {}", msg);
    std::process::exit(code)
}

fn parse_args(args: &[String]) -> Result<(Opts, Vec<Item>, Option<String>, Option<String>), CalcErr> {
    let un = |m: String| CalcErr::Unsupported(m);
    // validate the option names first (a typo must not silently change the route)
    let mut i = 2;
    let mut seen: Vec<&str> = vec![];
    while i < args.len() {
        let a = args[i].as_str();
        if VALUE_OPTS.contains(&a) {
            if i + 1 >= args.len() {
                return Err(un(format!("{} needs a value", a)));
            }
            if seen.contains(&a) {
                return Err(un(format!("{} given twice", a)));
            }
            seen.push(a);
            i += 2;
        } else if FLAG_OPTS.contains(&a) {
            i += 1;
        } else {
            return Err(un(format!("unknown argument {:?}\n{}", a, USAGE)));
        }
    }
    let get = |k: &str| args.iter().skip(2).position(|a| a == k).and_then(|i| args.get(i + 3)).cloned();
    let flag = |k: &str| args.iter().skip(2).any(|a| a == k);
    let ns: Vec<usize> = args[1].split(',').map(|x| x.trim().parse::<usize>().map_err(|_| un(format!("N list {:?}: expected comma-separated integers >= 1 (e.g. 2,3,4)", args[1])))).collect::<Result<_, _>>()?;
    if ns.is_empty() || ns.iter().any(|&n| n == 0) {
        return Err(un(format!("N list {:?}: every N must be at least 1", args[1])));
    }
    let mut ns_sorted = ns.clone();
    ns_sorted.sort();
    ns_sorted.dedup();
    let colors = match get("--colors") {
        Some(c) => Some(c.split(',').map(|x| x.trim().parse::<usize>().ok().filter(|&v| v >= 1).ok_or_else(|| un(format!("--colors {:?}: colors are integers >= 1", c)))).collect::<Result<Vec<_>, _>>()?),
        None => None,
    };
    let strands = match get("--strands") {
        Some(s) => Some(s.trim().parse::<usize>().ok().filter(|&v| (1..=64).contains(&v)).ok_or_else(|| un(format!("--strands {:?}: expected an integer in 1..=64", s)))?),
        None => None,
    };
    let method = get("--colored-method");
    if let Some(m) = &method {
        if !["cable", "cube", "singular"].contains(&m.as_str()) {
            return Err(un(format!("--colored-method must be cable, cube or singular, not {}", m)));
        }
    }
    let end = match get("--end") {
        Some(e) => EndPolicy::parse(&e).ok_or_else(|| un("--end must be right, left or alt".into()))?,
        None => EndPolicy::Right,
    };
    let closure = match get("--closure").as_deref() {
        None | Some("tower") => Closure::Tower(end),
        Some("exact") => Closure::Exact,
        Some(c) => Closure::parse(c).ok_or_else(|| un("--closure must be tower or exact".into()))?,
    };
    let num = |k: &str, d: Option<u64>| -> Result<Option<u64>, CalcErr> {
        match get(k) {
            Some(v) => v.trim().parse::<u64>().map(Some).map_err(|_| un(format!("{} {:?}: expected a non-negative integer", k, v))),
            None => Ok(d),
        }
    };
    let seed = num("--seed", Some(201))?.unwrap();
    let seed2 = num("--seed-f32b", Some(seed))?.unwrap();
    let retries = num("--seed-retries", Some(2))?.unwrap();
    let mem_mb = num("--mem-mb", None)?;
    let mem_cap_mb = num("--mem-cap-mb", None)?;
    if get("--field").is_some() {
        eprintln!("krhom: note: --field is a developer option and has no effect; results are over Q");
    }
    let opts = Opts {
        ns: ns_sorted,
        colors,
        strands,
        method,
        cube: flag("--cube"),
        reduced: flag("--reduced"),
        single_prime: flag("--single-prime"),
        no_simplify: flag("--no-simplify"),
        check_library: flag("--check-library"),
        box_live: flag("--box-live"),
        closure,
        lib_dir: get("--lib-dir"),
        seed,
        seed2,
        retries,
        mem_mb,
        mem_cap_mb,
        env_mem_cap: std::env::var_os("KRHOM_MEM_CAP_MB").is_some(),
        env_any_host: std::env::var_os("KRHOM_COSTS_ANY_HOST").is_some(),
        table: flag("--table"),
    };
    let nin = ["--braid", "--pd", "--file"].iter().filter(|k| get(k).is_some()).count();
    if nin != 1 {
        return Err(un(format!("give exactly one of --braid, --pd, --file\n{}", USAGE)));
    }
    let mut items = vec![];
    if let Some(b) = get("--braid") {
        items.push(Item { name: get("--name"), input: Input::Braid(b) });
    } else if let Some(p) = get("--pd") {
        items.push(Item { name: get("--name"), input: Input::Pd(p) });
    } else if let Some(f) = get("--file") {
        let text = std::fs::read_to_string(&f).map_err(|e| un(format!("--file {}: {}", f, e)))?;
        for (ln, line) in text.lines().enumerate() {
            let l = line.split('#').next().unwrap_or("").trim();
            if l.is_empty() {
                continue;
            }
            let (name, w) = match l.split_once(';') {
                Some((n, w)) => (n.trim().to_string(), w.trim().to_string()),
                None => (format!("line{}", ln + 1), l.to_string()),
            };
            items.push(Item { name: Some(name), input: Input::Braid(w) });
        }
        if items.is_empty() {
            return Err(un(format!("--file {}: no braids", f)));
        }
    }
    if opts.check_library && opts.method.as_deref().is_some_and(|m| m != "singular") {
        return Err(un("--check-library implies the singular route; drop --colored-method".into()));
    }
    Ok((opts, items, get("--json"), get("--csv")))
}

fn sigstr(w: &[i64]) -> String {
    if w.is_empty() {
        return "(empty)".into();
    }
    w.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(" ")
}

/// Run an engine closure, turning a panic into an engine error (exit 2, the batch continues).
fn guard<T>(what: &str, f: impl FnOnce() -> Result<T, CalcErr>) -> Result<T, CalcErr> {
    match std::panic::catch_unwind(std::panic::AssertUnwindSafe(f)) {
        Ok(r) => r,
        Err(e) => {
            let m = e.downcast_ref::<String>().cloned().or_else(|| e.downcast_ref::<&str>().map(|s| s.to_string())).unwrap_or_else(|| "panic".into());
            Err(CalcErr::Engine(format!("{} panicked: {} (please report this)", what, m)))
        }
    }
}

/// Two-prime certification of an F_p computation: run F31 (seed a) and, unless --single-prime,
/// F32b (seed b), and compare.
fn two_primes(o: &Opts, out: &mut Outcome, ns: &[usize], f31: impl FnOnce() -> Result<BTreeMap<usize, Poincare>, CalcErr>, f32b: impl FnOnce() -> Result<BTreeMap<usize, Poincare>, CalcErr>) -> Result<BTreeMap<usize, Poincare>, CalcErr> {
    let t = Instant::now();
    let a = guard("engine (F_2147483647)", f31)?;
    eprintln!("[F_2147483647] {:.2}s", t.elapsed().as_secs_f64());
    if o.single_prime {
        out.coefficients = "F_2147483647 (developer option --single-prime: not certified over Q)".into();
        return Ok(a);
    }
    let t = Instant::now();
    let b = guard("engine (F_4294967291)", f32b)?;
    eprintln!("[F_4294967291] {:.2}s", t.elapsed().as_secs_f64());
    compare_primes(ns, &a, &b).map_err(|e| {
        out.agree = Some(false);
        CalcErr::Disagree(e)
    })?;
    out.agree = Some(true);
    eprintln!("krhom: F_2147483647 and F_4294967291 agree: certified over both primes (equal over Q unless both are torsion primes)");
    Ok(a)
}

fn set_env(k: &str, v: &str) {
    // SAFETY: called from the main thread between items, when no engine thread is running
    unsafe { std::env::set_var(k, v) };
}

fn run_item(item: &Item, o: &Opts) -> Result<Outcome, CalcErr> {
    let un = |m: String| CalcErr::Unsupported(m);
    let mut out = Outcome { coefficients: "Q".into(), ..Default::default() };
    // the name heads the block before any early return (a failed batch line is identifiable)
    if let Some(n) = &item.name {
        println!("name: {}", n);
    }
    // ---- input
    let (mut bw, is_braid) = match &item.input {
        Input::Braid(s) => {
            let w = parse_braid(s).map_err(un)?;
            let need = w.iter().map(|x| x.unsigned_abs() as usize).max().unwrap_or(0) + 1;
            let strands = match o.strands {
                Some(n) if n < need => return Err(un(format!("--strands {} is too small for the braid (generator {} needs {} strands)", n, need - 1, need))),
                Some(n) => n,
                None => need,
            };
            (BraidWord::from_signed(strands, &w), true)
        }
        Input::Pd(s) => {
            let pd = parse_pd(s).map_err(|e| un(format!("invalid PD code: {}", e)))?;
            let r = std::panic::catch_unwind(|| pd_to_braid(&pd)).map_err(|_| un("invalid PD code (cannot convert it to a braid)".into()))?;
            (r, false)
        }
    };
    let ncomp = bw.num_components();
    let colors: Vec<usize> = match &o.colors {
        None => vec![1; ncomp],
        Some(c) if c.len() == 1 => vec![c[0]; ncomp],
        Some(c) if c.len() == ncomp => c.clone(),
        Some(c) => return Err(un(format!("{} colors given, but the link has {} components (give one color for all, or one per component)", c.len(), ncomp))),
    };
    let lab = |b: &BraidWord| -> Vec<usize> { b.components.iter().map(|&c| colors[c]).collect() };
    let mut labels = lab(&bw);
    out.colors = colors.clone();
    out.orig = Some(bw.clone());
    let uncolored = labels.iter().all(|&l| l == 1);
    let method = if o.cube { Some("cube") } else { o.method.as_deref() };
    let route_kind = match (method, uncolored) {
        (Some("cube"), _) => "cube",
        (Some("cable"), false) => "cable",
        (Some("cable"), true) => "soergel",
        (Some("singular"), _) => "singular",
        (None, true) => "soergel",
        (None, false) => "singular",
        _ => unreachable!(),
    };
    if o.check_library && route_kind != "singular" {
        return Err(un("--check-library applies to colored braids on 3-4 strands".into()));
    }
    if o.reduced && (route_kind == "cube" || route_kind == "cable") {
        return Err(un("--reduced is only supported for uncolored links with the Soergel engine and with the singular route".into()));
    }
    println!("gl_N for N in {:?}  strands {}  crossings {}  components {}  colors {:?}", o.ns, bw.strands, bw.word.len(), ncomp, colors);
    println!("braid: {}", sigstr(&bw.signed()));
    // ---- N below a color: wedge^a C^N = 0 (review K14 / Q5), on every route
    let amax = *labels.iter().max().unwrap_or(&1);
    out.zero_ns = o.ns.iter().copied().filter(|&n| n < amax).collect();
    let run_ns: Vec<usize> = o.ns.iter().copied().filter(|&n| n >= amax).collect();
    if !out.zero_ns.is_empty() {
        eprintln!("krhom: note: N = {:?} is below the color {}: the homology is 0 there (wedge^a C^N = 0 for a > N); not computed", out.zero_ns, amax);
    }
    if run_ns.is_empty() && !o.check_library {
        // every N is below the color: the answer is 0 on every route, before any route refusal
        // (no library, colors > 3) and without running an engine (D4 Q5, review cli #3)
        out.route = "zero";
        println!("method: none (every N is below the color {}: the homology is 0, nothing to compute)", amax);
        out.used = Some(bw);
        out.labels = labels;
        for &n in &out.zero_ns.clone() {
            out.results.insert(n, Poincare::default());
        }
        return Ok(out);
    }
    let word: Vec<Gen> = bw.word.clone();
    let results: BTreeMap<usize, Poincare> = match route_kind {
        "soergel" => {
            out.route = "soergel";
            println!("closure: {}{}", match o.closure { Closure::Exact => "exact".to_string(), Closure::Tower(e) => format!("tower ({:?} end first)", e) }, if o.reduced { ", reduced (strand at bottom position 0 marked)" } else { "" });
            out.seeds = vec![7, 11];
            let (n, cl, red) = (bw.strands, o.closure, o.reduced);
            let (w1, w2, ns1, ns2) = (word.clone(), word.clone(), run_ns.clone(), run_ns.clone());
            two_primes(o, &mut out, &run_ns, move || Ok(krhom::uncolored::engine2::khr_engine2_closure::<F31>(n, &w1, &ns1, &[], 7, red, cl).poincare.into_iter().collect()), move || Ok(krhom::uncolored::engine2::khr_engine2_closure::<F32b>(n, &w2, &ns2, &[], 11, red, cl).poincare.into_iter().collect()))?
        }
        "cable" => {
            out.route = "cable";
            let cable_strands: usize = labels.iter().sum();
            println!("colored method: cable ({} strands)", cable_strands);
            if !cable_default(&labels) {
                eprintln!("warning: a cable of this size ({} strands, colors {:?}) used 0.4-6 GB within seconds in tests", cable_strands, colors);
            }
            out.seeds = vec![7, 11];
            let (l1, l2, w1, w2, ns1, ns2, cl) = (labels.clone(), labels.clone(), word.clone(), word.clone(), run_ns.clone(), run_ns.clone(), o.closure);
            two_primes(
                o,
                &mut out,
                &run_ns,
                move || khr_colored_cable::<F31>(&l1, &w1, &ns1, 7, cl).map(|r| r.colored.into_iter().collect()).map_err(|e| CalcErr::Engine(format!("cable method failed: {}", e))),
                move || khr_colored_cable::<F32b>(&l2, &w2, &ns2, 11, cl).map(|r| r.colored.into_iter().collect()).map_err(|e| CalcErr::Engine(format!("cable method failed: {}", e))),
            )?
        }
        "cube" => {
            out.route = "cube";
            println!("{}method: full cube (exact over Q)", if uncolored { "" } else { "colored " });
            let mut r = BTreeMap::new();
            for &nn in &run_ns {
                let link = BraidLink { n: nn, labels: labels.clone(), word: word.clone() };
                r.insert(nn, guard("full cube", || Ok(khr_full_cube::<Q>(&link, QConv::A).poincare))?);
            }
            r
        }
        _ => {
            // ---- the singular route: simplification (spec B step 3)
            if is_braid && !o.no_simplify {
                let s = bw.simplify();
                if use_simplified(bw.strands, bw.components[0], s.strands, s.components.first().copied().unwrap_or(0), o.reduced, ncomp) {
                    println!("simplified braid: {}  (strands {})", sigstr(&s.signed()), s.strands);
                    bw = s;
                    labels = lab(&bw);
                }
            }
            let word: Vec<Gen> = bw.word.clone();
            let wi: Vec<i32> = bw.signed().iter().map(|&x| x as i32).collect();
            if bw.strands <= 2 && !o.box_live {
                out.route = "exact2";
                if o.check_library {
                    return Err(un("--check-library: this braid is on at most 2 strands (after simplification); the exact 2-strand engine needs no library".into()));
                }
                if labels.iter().any(|&a| a > 3) {
                    return Err(un(format!("colors {:?}: the exact 2-strand engine handles colors 1..=3; --colored-method cube is exact over Q (tiny inputs only)", colors)));
                }
                exact2_shape(&labels, &run_ns, o.reduced).map_err(un)?;
                out.mem_cap_mb = Some(std::env::var("KRHOM_MEM_CAP_MB").ok().and_then(|v| v.trim().parse().ok()).unwrap_or(krhom::colored::cclose::DEFAULT_MEM_CAP_MB));
                println!("method: singular (exact 2-strand singular Soergel engine){}", if o.reduced { ", reduced (strand at bottom position 0 marked)" } else { "" });
                let red = o.reduced;
                let (l1, l2, w1, w2, ns1, ns2) = (labels.clone(), labels.clone(), word.clone(), word.clone(), run_ns.clone(), run_ns.clone());
                let em = exact2_err;
                two_primes(o, &mut out, &run_ns, move || khr_singular_exact::<F31>(&l1, &w1, &ns1, red).map_err(em), move || khr_singular_exact::<F32b>(&l2, &w2, &ns2, red).map_err(em))?
            } else if o.box_live {
                out.route = "box-live";
                if bw.strands != 2 || wi.is_empty() {
                    return Err(un("--box-live (tests only) handles nonempty braids on exactly 2 strands".into()));
                }
                let lam: Vec<u8> = labels.iter().map(|&a| a as u8).collect();
                println!("method: singular (localized box engine, live provider, in-memory unpruned library, box N {}, seed {}){}", run_ns.iter().max().copied().unwrap_or(0), o.seed, if o.reduced { ", reduced (strand at bottom position 0 marked)" } else { "" });
                let costs = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/homprov/costs_cloud.txt");
                if run_ns.is_empty() {
                    BTreeMap::new()
                } else {
                    let pr = std::cell::RefCell::new(vec![]);
                    let (red, s1, s2) = (o.reduced, o.seed, o.seed2);
                    let r = two_primes(
                        o,
                        &mut out,
                        &run_ns,
                        || {
                            let p = run_prime_live::<F31>(&lam, &costs, &wi, &run_ns, red, s1)?;
                            let c = p.closed.clone();
                            pr.borrow_mut().push(p);
                            Ok(c)
                        },
                        || {
                            let p = run_prime_live::<F32b>(&lam, &costs, &wi, &run_ns, red, s2)?;
                            let c = p.closed.clone();
                            pr.borrow_mut().push(p);
                            Ok(c)
                        },
                    );
                    out.primes = pr.into_inner();
                    r?
                }
            } else {
                box_route(o, &mut out, &bw, &labels, &wi, &run_ns)?
            }
        }
    };
    out.used = Some(bw);
    out.labels = labels;
    out.results = results;
    for &n in &out.zero_ns.clone() {
        out.results.insert(n, Poincare::default());
    }
    Ok(out)
}

/// The exact 2-strand closure's shape limit (cclose::check_alg_shape): every unmarked block of
/// color a needs a(N - a) <= 255 (N <= 256 for color 1, N <= 129 for color 2, N <= 88 for
/// color 3). Checked before any engine work; exit 1.
fn exact2_shape(labels: &[usize], run_ns: &[usize], reduced: bool) -> Result<(), String> {
    for &nn in run_ns {
        for (i, &a) in labels.iter().enumerate() {
            if reduced && i == 0 || a > nn {
                continue;
            }
            if a * (nn - a) > 255 {
                return Err(format!(
                    "N = {} with color {}: the exact 2-strand engine supports a(N - a) <= 255 for every color a (here {} * {} = {}), i.e. N <= 256 for color 1, N <= 129 for color 2, N <= 88 for color 3",
                    nn,
                    a,
                    a,
                    nn - a,
                    a * (nn - a)
                ));
            }
        }
    }
    Ok(())
}

/// Classify an error of the exact 2-strand engine: a memory-cap refusal is a resource refusal
/// (exit 5, with the flag to raise), an unsupported shape is exit 1, anything else exit 2.
fn exact2_err(e: String) -> CalcErr {
    if e.contains("above the memory cap") {
        CalcErr::Resource(format!("{}; with krhom use --mem-cap-mb M to raise it (the run then needs about that much memory)", e))
    } else if e.contains("top weight") || e.contains("labels above 16") {
        CalcErr::Unsupported(format!("singular method: {}", e))
    } else {
        CalcErr::Engine(format!("singular method failed: {}", e))
    }
}

/// The bigraded table: rows t ascending, columns q ascending (step 2 when every q has the same
/// parity), entries = dimensions, blank for 0.
fn table_lines(p: &Poincare) -> Vec<String> {
    if p.c.is_empty() {
        return vec!["  table: 0".into()];
    }
    let (tmin, tmax) = (p.c.keys().map(|k| k.0).min().unwrap(), p.c.keys().map(|k| k.0).max().unwrap());
    let (qmin, qmax) = (p.c.keys().map(|k| k.1).min().unwrap(), p.c.keys().map(|k| k.1).max().unwrap());
    let step = if p.c.keys().all(|k| (k.1 - qmin).rem_euclid(2) == 0) { 2 } else { 1 };
    let qs: Vec<i32> = (qmin..=qmax).step_by(step).collect();
    let w = qs.iter().map(|q| q.to_string().len()).chain(p.c.values().map(|k| k.to_string().len())).max().unwrap() + 2;
    let lw = (tmin..=tmax).map(|t| t.to_string().len()).max().unwrap().max(3);
    let mut out = vec!["  table (rows t, columns q):".to_string()];
    let mut head = format!("  {:>lw$}", "t\\q", lw = lw);
    for q in &qs {
        head.push_str(&format!("{:>w$}", q, w = w));
    }
    out.push(head);
    for t in tmin..=tmax {
        let mut row = format!("  {:>lw$}", t, lw = lw);
        for q in &qs {
            match p.c.get(&(t, *q)) {
                Some(k) => row.push_str(&format!("{:>w$}", k, w = w)),
                None => row.push_str(&" ".repeat(w)),
            }
        }
        out.push(row.trim_end().to_string());
    }
    out
}

/// RFC 4180 field: quoted (with doubled quotes) when it contains a comma, quote, CR or LF.
fn csv_field(s: &str) -> String {
    if s.contains([',', '"', '\n', '\r']) {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
}

/// Write `body` to `path` atomically when possible (temporary file in the same directory, then
/// rename), so an abort never leaves a half-written file; falls back to a direct write.
fn write_atomic(path: &str, body: &str) -> std::io::Result<()> {
    let tmp = format!("{}.tmp{}", path, std::process::id());
    if std::fs::write(&tmp, body).is_ok() && std::fs::rename(&tmp, path).is_ok() {
        return Ok(());
    }
    let _ = std::fs::remove_file(&tmp);
    std::fs::write(path, body)
}

/// The 3-4 strand box route (spec B step 5, C4).
fn box_route(o: &Opts, out: &mut Outcome, bw: &BraidWord, labels: &[usize], wi: &[i32], run_ns: &[usize]) -> Result<BTreeMap<usize, Poincare>, CalcErr> {
    let un = |m: String| CalcErr::Unsupported(m);
    out.route = "box";
    let cube_hint = "--colored-method cube is exact over Q but exponential (tiny inputs only)";
    if labels.iter().all(|&l| l == 1) {
        return Err(un("uncolored input: drop --colored-method singular; the Soergel engine handles any number of strands".into()));
    }
    if labels.iter().any(|&l| l != labels[0]) {
        return Err(un(format!("labels {:?} on {} strands are not uniform; supported on 3-4 strands: {}. {}", labels, bw.strands, supported_matrix(), cube_hint)));
    }
    let lam: Vec<u8> = labels.iter().map(|&a| a.min(255) as u8).collect();
    let Some(reg31) = reg_for(&lam, "f31") else {
        return Err(un(format!("no library for labels {:?} on {} strands; supported on 3-4 strands: {}. {}", labels, bw.strands, supported_matrix(), cube_hint)));
    };
    let reg32 = reg_for(&lam, "f32b").expect("registry has both primes");
    if let Some(&bad) = o.ns.iter().find(|&&n| n > reg31.n as usize) {
        return Err(un(format!("N = {} is above the library's N <= {} for labels {:?} (library {}); supported: {}", bad, reg31.n, labels, reg31.id, supported_matrix())));
    }
    let used: std::collections::BTreeSet<usize> = bw.word.iter().map(|g| g.i).collect();
    if let Some(i) = (0..bw.strands - 1).find(|i| !used.contains(i)) {
        return Err(un(format!("split link (sigma_{} never occurs): compute the pieces separately", i + 1)));
    }
    let refused = refused_env_set();
    if !refused.is_empty() {
        return Err(un(format!("validation-only environment switches change the computation on the box route: unset {}", refused.join(", "))));
    }
    // environment before any engine thread (spec C6; review K7: thread-aware closure cap)
    if !o.env_any_host {
        set_env("KRHOM_COSTS_ANY_HOST", "1");
    }
    let threads = krhom::num_threads() as u64;
    let cap = match (o.mem_cap_mb, o.env_mem_cap) {
        (Some(c), _) => {
            set_env("KRHOM_MEM_CAP_MB", &c.to_string());
            c
        }
        (None, true) => std::env::var("KRHOM_MEM_CAP_MB").ok().and_then(|v| v.trim().parse().ok()).unwrap_or(0),
        (None, false) => {
            let c = if bw.strands == 3 { 4500 * (threads / 2).max(1) } else { 200_000 };
            set_env("KRHOM_MEM_CAP_MB", &c.to_string());
            c
        }
    };
    out.mem_cap_mb = Some(cap);
    let mem_mb = o.mem_mb.unwrap_or(if bw.strands == 3 { 3200 } else { 3000 });
    out.mem_mb = Some(mem_mb);
    if run_ns.is_empty() && !o.check_library {
        // every N is below the color: the homology is 0, no library or engine needed (Q5)
        println!("method: singular (every N is below the color: zero, nothing to compute)");
        return Ok(BTreeMap::new());
    }
    let root = lib_root(o.lib_dir.as_deref());
    let inst31 = locate(&reg31, &root)?;
    let inst32 = if o.single_prime { None } else { Some(locate(&reg32, &root)?) };
    out.library = Some(format!("{{\"id\":{},\"fields\":[{}],\"data_commit\":{},\"dir\":{}}}", jstr(&reg31.id), if o.single_prime { "\"F_2147483647\"".to_string() } else { "\"F_2147483647\",\"F_4294967291\"".to_string() }, jstr(&reg31.data_commit), jstr(&root.join(&reg31.id).display().to_string())));
    let seeds = if o.seed == o.seed2 || o.single_prime { format!("seed {}", o.seed) } else { format!("seeds F31 {} / F32b {}", o.seed, o.seed2) };
    println!(
        "method: singular (localized box engine; library {} @ data {}, labels ({}), prune {}, box N {}, {}){}",
        reg31.id,
        &reg31.data_commit[..7],
        lam.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","),
        reg31.n,
        reg31.n,
        seeds,
        if o.reduced { ", reduced (strand at bottom position 0 marked)" } else { "" }
    );
    out.seeds = vec![o.seed, o.seed2];
    let mk = |seed: u64, replay_only: bool| BoxParams { word: wi.to_vec(), ns: run_ns.to_vec(), reduced: o.reduced, seed, retries: o.retries, mem_mb, check_only: o.check_library, replay_only };
    // one prime with the library diagnosis of review K2 and the template-recompute retry (Q3)
    let status_files: Vec<PathBuf> = [Some(&inst31), inst32.as_ref()].into_iter().flatten().filter_map(|i| i.status.clone()).collect();
    let complete = inst31.status.is_none();
    let diag = |e: CalcErr| -> CalcErr {
        if e.code() != 4 || e.msg().contains("this binary's route stamp") || e.msg().contains("does not cover this word") {
            return e; // not a library error, or the stamp preflight's own message
        }
        let first = e.msg().lines().next().unwrap_or("").chars().take(400).collect::<String>();
        if !complete {
            eprintln!("krhom: library error ({}); diagnosing which jobs the word needs ...", first);
            match diagnose_missing(&lam, reg31.n, wi, &status_files) {
                Ok(m) if !m.is_empty() => {
                    return CalcErr::Library(format!("library {} does not cover this word: {} jobs not certified for both primes: {}", reg31.id, m.len(), m.join(", ")));
                }
                Ok(_) => {}
                Err(d) => eprintln!("krhom: diagnosis failed: {}", d),
            }
            CalcErr::Library(format!("library {} files are present but were rejected on load (see the [sindec] layer A ... (rejected) lines above; template cache or stamp): {}; run tools/fetch_library.sh --force {}", reg31.id, first, reg31.id))
        } else {
            CalcErr::Library(format!("damaged install: library {} is complete but a map is missing or rejected ({}): run tools/fetch_library.sh --force {}", reg31.id, first, reg31.id))
        }
    };
    // a word that needs uncertified jobs is reported before any template-recompute retry (which
    // would only rebuild the library to fail again)
    let missing = || -> Option<CalcErr> {
        if complete {
            return None;
        }
        eprintln!("krhom: a layer-A map is missing; diagnosing which jobs the word needs ...");
        match diagnose_missing(&lam, reg31.n, wi, &status_files) {
            Ok(m) if !m.is_empty() => Some(CalcErr::Library(format!("library {} does not cover this word: {} jobs not certified for both primes: {}", reg31.id, m.len(), m.join(", ")))),
            Ok(_) => None,
            Err(d) => {
                eprintln!("krhom: diagnosis failed: {}", d);
                None
            }
        }
    };
    fn one<F: krhom::colored::boxclose::PrimeField>(reg: &RegEntry, inst: &InstalledLib, p: BoxParams, retry: BoxParams, missing: &dyn Fn() -> Option<CalcErr>) -> Result<PrimeOut, CalcErr> {
        match guard(&format!("box engine ({})", F::name()), || run_prime::<F>(reg, inst, &p)) {
            Err(e) if e.code() == 4 && e.msg().contains("not in layer A") && p.replay_only => {
                if let Some(m) = missing() {
                    return Err(m);
                }
                eprintln!("[{}] a layer-A map was not found in replay-only mode; retrying once with template recomputation enabled (deterministic; loaded maps are re-certified)", F::name());
                guard(&format!("box engine ({})", F::name()), || run_prime::<F>(reg, inst, &retry))
            }
            r => r,
        }
    }
    let report = |p: &PrimeOut| {
        if p.check_only {
            eprintln!("[{}] library: {} layer-A loads ({} rejected, {} by route stamp; preflight {} stamped files ok, {} legacy), {} types, {:.1}s; record scan {}; {:.1}s; VmHWM {} MB", p.field, p.layer_a.0, p.layer_a.2, p.layer_a.1, p.stamp_ok, p.stamp_legacy, p.lib_types, p.lib_build_secs, if p.record_scan { "run" } else { "not needed (closure data for every closable type)" }, p.secs, p.vmhwm_kb / 1024);
        } else {
            eprintln!(
                "[{}] library: {} layer-A loads ({} rejected, {} by route stamp; preflight {} stamped files ok, {} legacy), {:.1}s; {} box points, record {:.1}s, scan {:.1}s, close {:.1}s; chi = decat_chi for {} closures; scaling: {} entries checked ({} zero at both, {} declared 0); determinant pass {} checks at {} points; seed {}{}; total {:.1}s; VmHWM {} MB",
                p.field, p.layer_a.0, p.layer_a.2, p.layer_a.1, p.stamp_ok, p.stamp_legacy, p.lib_build_secs, p.points, p.record_secs, p.scan_secs, p.close_secs, p.chi_checked, p.scaling.0, p.scaling.1, p.scaling.2, p.det_pass.1, p.det_pass.0, p.seed, if p.template_recompute { "; templates recomputed" } else { "" }, p.secs, p.vmhwm_kb / 1024
            );
        }
    };
    let a = one::<F31>(&reg31, &inst31, mk(o.seed, true), mk(o.seed, false), &missing).map_err(diag)?;
    report(&a);
    out.primes.push(a.clone());
    if let Some(inst32) = &inst32 {
        let b = one::<F32b>(&reg32, inst32, mk(o.seed2, true), mk(o.seed2, false), &missing).map_err(diag)?;
        report(&b);
        out.primes.push(b.clone());
        if !o.check_library {
            if let Err(e) = compare_primes(run_ns, &a.closed, &b.closed) {
                out.agree = Some(false);
                if o.seed != o.seed2 {
                    return Err(CalcErr::Disagree(format!("{}\n  (seeds differ: {} / {}; rerun with equal seeds to separate prime dependence from seed dependence)", e, o.seed, o.seed2)));
                }
                return Err(CalcErr::Disagree(e));
            }
            out.agree = Some(true);
            eprintln!("krhom: F_2147483647 and F_4294967291 agree: certified over both primes (equal over Q unless both are torsion primes)");
        }
    } else {
        out.coefficients = "F_2147483647 (developer option --single-prime: not certified over Q)".into();
    }
    if o.check_library {
        out.check_library = Some(format!("library {} covers this word ({} prime{}; no box scan run)", reg31.id, out.primes.len(), if out.primes.len() == 1 { "" } else { "s" }));
    }
    Ok(a.closed)
}

fn item_json(item: &Item, o: &Opts, r: &Result<Outcome, CalcErr>, secs: f64) -> String {
    let name = item.name.as_deref().map_or("null".to_string(), jstr);
    let (inp, inkind) = match &item.input {
        Input::Braid(s) => (s.clone(), "braid"),
        Input::Pd(s) => (s.clone(), "pd"),
    };
    let env: Vec<String> = krhom_env().iter().map(|(k, v)| format!("{}:{}", jstr(k), jstr(v))).collect();
    let common = format!("\"tool\":\"krhom\",\"commit\":{},\"name\":{},\"input\":{{\"kind\":\"{}\",\"text\":{},\"N\":{:?},\"reduced\":{}", jstr(krhom::colored::boxrun::build_commit()), name, inkind, jstr(&inp), o.ns, o.reduced);
    match r {
        Ok(out) => {
            let bw = |b: &Option<BraidWord>| b.as_ref().map_or("null".to_string(), |b| format!("{:?}", b.signed()));
            let st = |b: &Option<BraidWord>| b.as_ref().map_or(0, |b| b.strands);
            let res: Vec<String> = out.results.iter().map(|(n, p)| format!("{{\"N\":{},\"P\":{},\"euler\":{}}}", n, poly_json(p), jstr(&p.euler().to_string()))).collect();
            format!(
                "{{{},\"braid\":{},\"strands\":{},\"colors\":{:?},\"labels\":{:?},\"used_braid\":{},\"used_strands\":{}}},\"route\":{},\"coefficients\":{},\"library\":{},\"seeds\":{:?},\"mem_mb\":{},\"mem_cap_mb\":{},\"primes\":[{}],\"agree\":{},\"zero_N\":{:?},\"check_library\":{},\"results\":[{}],\"secs\":{:.2},\"env\":{{{}}},\"error\":null}}",
                common,
                bw(&out.orig),
                st(&out.orig),
                out.colors,
                out.labels,
                bw(&out.used),
                st(&out.used),
                jstr(out.route),
                jstr(&out.coefficients),
                out.library.clone().unwrap_or("null".into()),
                out.seeds,
                out.mem_mb.map_or("null".into(), |v| v.to_string()),
                out.mem_cap_mb.map_or("null".into(), |v| v.to_string()),
                out.primes.iter().map(prime_json).collect::<Vec<_>>().join(","),
                out.agree.map_or("null".into(), |a| a.to_string()),
                out.zero_ns,
                out.check_library.as_deref().map_or("null".into(), jstr),
                res.join(","),
                secs,
                env.join(","),
            )
        }
        Err(e) => format!("{{{}}},\"results\":[],\"secs\":{:.2},\"env\":{{{}}},\"error\":{{\"code\":{},\"kind\":{},\"message\":{}}}}}", common, secs, env.join(","), e.code(), jstr(e.kind()), jstr(e.msg())),
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.iter().any(|a| a == "--help" || a == "-h") {
        println!("{}", USAGE);
        return;
    }
    if args.len() < 3 || args[1].starts_with('-') {
        fail(1, USAGE);
    }
    let (o, items, json, csv) = parse_args(&args).unwrap_or_else(|e| fail(e.code(), e.msg()));
    let batch = items.len() > 1 || matches!(items.first(), Some(Item { name: Some(_), .. })) && args.iter().any(|a| a == "--file");
    // --mem-cap-mb applies on every route: set it before any engine work (review cli #1)
    if let Some(c) = o.mem_cap_mb {
        set_env("KRHOM_MEM_CAP_MB", &c.to_string());
    }
    // outputs are written incrementally, one record per finished item, so an abort (OOM under
    // ulimit -v, SIGKILL, timeout) keeps the finished results (review cli #2); a bad path fails
    // before any work
    let mut csv_file = csv.as_ref().map(|p| std::fs::File::create(p).unwrap_or_else(|e| fail(1, &format!("--csv {}: {}", p, e))));
    let csv_write = |f: &mut Option<std::fs::File>, rows: &[String]| {
        use std::io::Write;
        if rows.is_empty() {
            return;
        }
        if let Some(f) = f.as_mut() {
            let mut body = rows.join("\n");
            body.push('\n');
            if let Err(e) = f.write_all(body.as_bytes()).and_then(|_| f.flush()) {
                fail(1, &format!("--csv {}: {}", csv.as_deref().unwrap_or(""), e));
            }
        }
    };
    csv_write(&mut csv_file, &["name,N,reduced,t,q,rank".to_string()]);
    let json_write = |jsons: &[String]| {
        if let Some(p) = &json {
            let body = if batch { if jsons.is_empty() { "[\n]\n".to_string() } else { format!("[\n{}\n]\n", jsons.join(",\n")) } } else { format!("{}\n", jsons.first().map_or("", |s| s.as_str())) };
            if let Err(e) = write_atomic(p, &body) {
                fail(1, &format!("--json {}: {}", p, e));
            }
        }
    };
    if batch {
        json_write(&[]);
    } else if let Some(p) = &json {
        // check the path is writable now (the record itself is written when the item finishes)
        if let Err(e) = std::fs::OpenOptions::new().create(true).append(true).open(p) {
            fail(1, &format!("--json {}: {}", p, e));
        }
    }
    let mut code = 0;
    let mut jsons = vec![];
    for item in &items {
        let t0 = Instant::now();
        let r = run_item(item, &o);
        let secs = t0.elapsed().as_secs_f64();
        let mut csv_rows = vec![];
        match &r {
            Ok(out) => {
                if let Some(c) = &out.check_library {
                    println!("check-library: {}", c);
                } else {
                    println!("coefficients: {}", out.coefficients);
                    println!("time {:.2}s", secs);
                    let nm = csv_field(item.name.as_deref().unwrap_or_default());
                    for (nn, p) in &out.results {
                        println!("N = {}", nn);
                        println!("  Poincare polynomial: {}", p);
                        println!("  total rank: {}", p.total_rank());
                        println!("  Euler characteristic: {}", p.euler());
                        if o.table {
                            for l in table_lines(p) {
                                println!("{}", l);
                            }
                        }
                        if p.c.is_empty() {
                            csv_rows.push(format!("{},{},{},,,0", nm, nn, o.reduced));
                        }
                        for (&(t, q), &k) in &p.c {
                            csv_rows.push(format!("{},{},{},{},{},{}", nm, nn, o.reduced, t, q, k));
                        }
                    }
                }
            }
            Err(e) => {
                let who = item.name.as_deref().map_or(String::new(), |n| format!("{}: ", n));
                eprintln!("krhom: {}{}", who, e.msg());
                if batch {
                    println!("error: {}exit {} ({})", who, e.code(), e.kind());
                }
                if code == 0 {
                    code = e.code();
                }
            }
        }
        if batch {
            println!();
        }
        {
            use std::io::Write;
            let _ = std::io::stdout().flush();
        }
        csv_write(&mut csv_file, &csv_rows);
        jsons.push(item_json(item, &o, &r, secs));
        json_write(&jsons);
    }
    std::process::exit(code);
}
