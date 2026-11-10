//! schur_scan: command-line front end of the decategorified oracle `krhom::schur` (M3a).
//!
//! schur_scan basis   --lam 2,2,2,2 [--mu ...] [--save DIR] [--json]
//!     prints: types, sum of ranks, max rank (+table), rank(p_max), kappa, sum flows (bfs, minflow), time, peak RSS
//! schur_scan scan    --labels 2,2,2,2 --word "1 1 2 -3 2 1 -3 -2 -2 -3 -3" [--peak] [--every-stage] [--marked 0] [--load DIR] [--json]
//!     prints per stage and final: unsimplified objects/ldim, LB objects/types/ldim/sum S, max rank present,
//!     kappa present, peak; Gram histogram (weight: #objects, unknowns, unknowns reduced)
//! schur_scan chi     --labels 2,2 --word "1 1 1" --n 3[,4,...] [--json]
//! --json emits one JSON object per line (no JSON dependency).
//!
//! Basis files written by --save / read by --load are named cb_<lam>_<mu>.txt (labels joined by '-').

use krhom::schur::*;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::time::Instant;

#[repr(C)]
struct RUsage {
    utime: [i64; 2],
    stime: [i64; 2],
    maxrss: i64,
    rest: [i64; 13],
}
unsafe extern "C" {
    fn getrusage(who: i32, usage: *mut RUsage) -> i32;
}
/// peak resident set size of this process in bytes (getrusage; macOS reports bytes, Linux KiB)
fn peak_rss() -> u64 {
    let mut ru = RUsage { utime: [0; 2], stime: [0; 2], maxrss: 0, rest: [0; 13] };
    let rc = unsafe { getrusage(0, &mut ru) };
    if rc != 0 {
        return 0;
    }
    if cfg!(target_os = "macos") { ru.maxrss as u64 } else { ru.maxrss as u64 * 1024 }
}

fn usage() -> ! {
    eprintln!(
        "usage:\n  schur_scan basis --lam 2,2,2,2 [--mu ...] [--save DIR] [--json]\n  schur_scan scan --labels 2,2,2,2 --word \"1 1 2 -3\" [--peak] [--every-stage] [--marked 0] [--load DIR] [--json]\n  schur_scan chi --labels 2,2 --word \"1 1 1\" --n 3[,4,..] [--json]"
    );
    std::process::exit(2)
}

struct Args {
    cmd: String,
    kv: BTreeMap<String, String>,
    flags: Vec<String>,
}
fn parse_args() -> Args {
    let a: Vec<String> = std::env::args().skip(1).collect();
    if a.is_empty() {
        usage();
    }
    let mut kv = BTreeMap::new();
    let mut flags = vec![];
    let mut i = 1;
    while i < a.len() {
        let k = a[i].as_str();
        match k {
            "--peak" | "--every-stage" | "--json" => flags.push(k.to_string()),
            "--lam" | "--mu" | "--save" | "--labels" | "--word" | "--marked" | "--load" | "--n" => {
                if i + 1 >= a.len() {
                    usage();
                }
                kv.insert(k.to_string(), a[i + 1].clone());
                i += 1;
            }
            _ => {
                eprintln!("unknown argument {}", k);
                usage();
            }
        }
        i += 1;
    }
    Args { cmd: a[0].clone(), kv, flags }
}
fn comp(s: &str) -> Comp {
    s.split(',').map(|x| x.trim().parse().unwrap_or_else(|_| usage())).collect()
}
fn comp_tag(c: &[u8]) -> String {
    c.iter().map(|x| x.to_string()).collect::<Vec<_>>().join("-")
}
fn json_arr<T: std::fmt::Display>(v: &[T]) -> String {
    format!("[{}]", v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","))
}
fn basis_file(dir: &Path, lam: &[u8], mu: &[u8]) -> PathBuf {
    dir.join(format!("cb_{}_{}.txt", comp_tag(lam), comp_tag(mu)))
}
fn mb(b: u64) -> f64 {
    b as f64 / (1024.0 * 1024.0)
}

fn cmd_basis(args: &Args) {
    let lam = comp(args.kv.get("--lam").unwrap_or_else(|| usage()));
    let mu = args.kv.get("--mu").map(|s| comp(s)).unwrap_or_else(|| lam.clone());
    let json = args.flags.iter().any(|f| f == "--json");
    let t0 = Instant::now();
    let mut s = Schur::new(&lam);
    let cb = s.canonical_basis(&mu);
    let (order, lb, lm) = (cb.order.clone(), cb.ladder_bfs.clone(), cb.ladder_minflow.clone());
    let entries: usize = cb.b.values().map(|e| e.t.len()).sum();
    let t_basis = t0.elapsed().as_secs_f64();
    let mut sum_rank: u128 = 0;
    let (mut maxr, mut arg) = (0u128, order[0]);
    for &p in &order {
        let r = s.rank(&mu, p);
        sum_rank += r;
        if r > maxr {
            maxr = r;
            arg = p;
        }
    }
    let pmax = *order.last().unwrap();
    let r_pmax = s.rank(&mu, pmax);
    let kappa = s.kappa(&mu);
    let fb: u128 = order.iter().map(|p| flows(&lam, &lb[p])).sum();
    let fm: u128 = order.iter().map(|p| flows(&lam, &lm[p])).sum();
    let longest = order.iter().map(|p| lb[p].len()).max().unwrap_or(0);
    let t_all = t0.elapsed().as_secs_f64();
    if let Some(dir) = args.kv.get("--save") {
        let dir = PathBuf::from(dir);
        std::fs::create_dir_all(&dir).expect("cannot create --save directory");
        let path = basis_file(&dir, &lam, &mu);
        s.save_basis(&mu, &path).expect("save_basis failed");
        if !json {
            println!("saved {}", path.display());
        }
    }
    let rss = peak_rss();
    if json {
        println!(
            "{{\"kind\":\"basis\",\"lam\":{},\"mu\":{},\"types\":{},\"entries\":{},\"sum_ranks\":{},\"max_rank\":{},\"max_rank_table\":\"{}\",\"pmax_table\":\"{}\",\"rank_pmax\":{},\"kappa\":{},\"sum_flows_bfs\":{},\"sum_flows_minflow\":{},\"longest_bfs_ladder\":{},\"time_basis_s\":{:.3},\"time_s\":{:.3},\"peak_rss\":{}}}",
            json_arr(&lam),
            json_arr(&mu),
            order.len(),
            entries,
            sum_rank,
            maxr,
            s.table(arg).fmt(),
            s.table(pmax).fmt(),
            r_pmax,
            kappa,
            fb,
            fm,
            longest,
            t_basis,
            t_all,
            rss
        );
    } else {
        println!("block {:?} -> {:?}", lam, mu);
        println!("  types {}  (entries in all b_p: {})", order.len(), entries);
        println!("  sum of ranks {}", sum_rank);
        println!("  max rank {} at {}", maxr, s.table(arg).fmt());
        println!("  rank(p_max) {} (p_max = {})", r_pmax, s.table(pmax).fmt());
        println!("  kappa {}", kappa);
        println!("  sum flows: bfs ladders {}, minimal-flow ladders {}; longest bfs ladder {}", fb, fm, longest);
        println!("  time: basis {:.2}s, total {:.2}s; peak RSS {:.1} MB", t_basis, t_all, mb(rss));
    }
}

fn load_dir(s: &mut Schur, dir: &Path) {
    let prefix = format!("cb_{}_", comp_tag(&s.lam));
    let rd = std::fs::read_dir(dir).expect("cannot read --load directory");
    for e in rd.flatten() {
        let name = e.file_name().to_string_lossy().to_string();
        if name.starts_with(&prefix) && name.ends_with(".txt") {
            s.load_basis(&e.path()).unwrap_or_else(|err| panic!("load_basis {}: {}", name, err));
            eprintln!("loaded {}", name);
        }
    }
}

fn cmd_scan(args: &Args) {
    let labels = comp(args.kv.get("--labels").unwrap_or_else(|| usage()));
    let wstr = args.kv.get("--word").unwrap_or_else(|| usage()).clone();
    let word: Vec<i32> = wstr.split_whitespace().map(|x| x.parse().unwrap_or_else(|_| usage())).collect();
    let json = args.flags.iter().any(|f| f == "--json");
    let opts = ScanOpts {
        peak: args.flags.iter().any(|f| f == "--peak"),
        every_stage: args.flags.iter().any(|f| f == "--every-stage"),
        marked_block: args.kv.get("--marked").map(|m| m.parse().unwrap_or_else(|_| usage())),
    };
    let t0 = Instant::now();
    let mut s = Schur::new(&labels);
    if let Some(dir) = args.kv.get("--load") {
        load_dir(&mut s, Path::new(dir));
    }
    let res = decat_scan(&mut s, &word, &opts);
    let dt = t0.elapsed().as_secs_f64();
    let rss = peak_rss();
    let wjson = format!("\"{}\"", wstr.split_whitespace().collect::<Vec<_>>().join(" "));
    for st in &res.stages {
        let peak = match st.peak_before_next {
            Some((a, b, c)) => format!("[{},{},{}]", a, b, c),
            None => "null".to_string(),
        };
        if json {
            println!(
                "{{\"kind\":\"stage\",\"labels\":{},\"word\":{},\"step\":{},\"stage_labels\":{},\"unsimplified_objects\":{},\"unsimplified_ldim\":{},\"lb_objects\":{},\"lb_types\":{},\"lb_ldim\":{},\"lb_sum_s\":{},\"max_rank_present\":{},\"kappa_present\":{},\"peak_before_next\":{}}}",
                json_arr(&labels),
                wjson,
                st.step,
                json_arr(&st.labels),
                st.unsimplified_objects,
                st.unsimplified_ldim,
                st.lb_objects,
                st.lb_types,
                st.lb_ldim,
                st.lb_sum_s,
                st.max_rank_present,
                st.kappa_present,
                peak
            );
        } else {
            println!(
                "step {:>3} labels {:?}: unsimplified {} (ldim {}); LB objects {} types {} ldim {} sum S {}; max rank {} kappa {}; peak before next {}",
                st.step,
                st.labels,
                st.unsimplified_objects,
                st.unsimplified_ldim,
                st.lb_objects,
                st.lb_types,
                st.lb_ldim,
                st.lb_sum_s,
                st.max_rank_present,
                st.kappa_present,
                peak
            );
        }
    }
    // worst peak over the scan
    let worst = res.stages.iter().filter_map(|st| st.peak_before_next.map(|p| (st.step + 1, p))).max_by(|x, y| x.1.0.cmp(&y.1.0).then(y.0.cmp(&x.0)));
    // Gram histogram by weight: (#types, #objects, unknowns, unknowns reduced)
    let mut objs: BTreeMap<Tid, u128> = BTreeMap::new();
    for &(p, _, c) in &res.k0 {
        *objs.entry(p).or_default() += c.unsigned_abs();
    }
    let mut hist: BTreeMap<i32, (usize, u128, u128, u128)> = BTreeMap::new();
    for g in &res.gram {
        let e = hist.entry(g.weight).or_insert((0, 0, g.unknowns, g.unknowns_reduced));
        e.0 += 1;
        e.1 += objs[&g.p];
    }
    if json {
        for (w, (nt, no, u, ur)) in &hist {
            println!("{{\"kind\":\"gram\",\"labels\":{},\"word\":{},\"weight\":{},\"types\":{},\"objects\":{},\"unknowns\":{},\"unknowns_reduced\":{}}}", json_arr(&labels), wjson, w, nt, no, u, ur);
        }
        let wk = match worst {
            Some((j, (a, b, c))) => format!("{{\"crossing\":{},\"objects\":{},\"ldim\":{},\"sum_s\":{}}}", j, a, b, c),
            None => "null".to_string(),
        };
        println!("{{\"kind\":\"summary\",\"labels\":{},\"word\":{},\"worst_peak\":{},\"time_s\":{:.3},\"peak_rss\":{}}}", json_arr(&labels), wjson, wk, dt, rss);
    } else {
        if let Some((j, (a, b, c))) = worst {
            println!("worst peak: crossing {}: objects {}, localized dim {}, sum S {}", j, a, b, c);
        }
        if !hist.is_empty() {
            println!("Gram histogram (weight: #types, #objects, unknowns, unknowns reduced):");
            for (w, (nt, no, u, ur)) in &hist {
                println!("  {:>3}: {:>5} {:>8} {:>10} {:>10}", w, nt, no, u, ur);
            }
        }
        println!("time {:.2}s; peak RSS {:.1} MB", dt, mb(rss));
    }
}

fn cmd_chi(args: &Args) {
    let labels = comp(args.kv.get("--labels").unwrap_or_else(|| usage()));
    let wstr = args.kv.get("--word").unwrap_or_else(|| usage()).clone();
    let word: Vec<i32> = wstr.split_whitespace().map(|x| x.parse().unwrap_or_else(|_| usage())).collect();
    let ns: Vec<usize> = args.kv.get("--n").unwrap_or_else(|| usage()).split(',').map(|x| x.trim().parse().unwrap_or_else(|_| usage())).collect();
    let json = args.flags.iter().any(|f| f == "--json");
    let t0 = Instant::now();
    let mut s = Schur::new(&labels);
    let res = decat_scan(&mut s, &word, &ScanOpts { peak: false, every_stage: false, marked_block: None });
    for n in ns {
        let chi = decat_chi(&mut s, &res, &word, n);
        if json {
            let terms: Vec<String> = chi.c.iter().map(|(e, k)| format!("[{},{}]", e, k)).collect();
            println!("{{\"kind\":\"chi\",\"labels\":{},\"word\":\"{}\",\"n\":{},\"chi\":[{}],\"chi_at_1\":{}}}", json_arr(&labels), wstr, n, terms.join(","), chi.eval_at_one());
        } else {
            println!("N = {}: chi = {}", n, chi);
            println!("        chi(1) = {}", chi.eval_at_one());
        }
    }
    if !json {
        println!("time {:.2}s; peak RSS {:.1} MB", t0.elapsed().as_secs_f64(), mb(peak_rss()));
    }
}

fn main() {
    let args = parse_args();
    match args.cmd.as_str() {
        "basis" => cmd_basis(&args),
        "scan" => cmd_scan(&args),
        "chi" => cmd_chi(&args),
        _ => usage(),
    }
}
