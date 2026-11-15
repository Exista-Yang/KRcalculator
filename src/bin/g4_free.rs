//! Gate G4 (SYNTHESIS §3.1; C2 M1(b), U3): freeness of Δ_a R_N over A = H*(Gr(a, N)), the
//! hypothesis of the (Δ_λ G, F) adapter in src/cclose.rs.
//!
//! Command-line front end of the library check `krhom::colored::cclose::g4::check`, whose module doc
//! states the certificate (A), (B) (exact, for every (a, N)) and the cross-checks (C), (D)
//! (C2's direct method). tests/cclose.rs re-runs the certificate on every pair of
//! `cclose::G4_CERTIFIED`.
//!
//! Usage: g4_free [--direct] [--gens] [--max-a A] [--max-n N] [a:N ...]
//! Default: all 2 ≤ a ≤ 5, a ≤ N ≤ 10, certificate only. Output: one line per pair, and the
//! certified list in the form of `cclose::G4_CERTIFIED`. Supported pairs: 1 ≤ a ≤ min(N, 8),
//! N ≤ 20 (`g4::MAX_N`); any other pair is a usage error (exit status 2).

use krhom::colored::cclose::g4::{binom, check, fact, pi_a};
use std::collections::BTreeMap;
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let direct = args.iter().any(|a| a == "--direct");
    let gens = args.iter().any(|a| a == "--gens");
    let get = |name: &str, def: usize| -> usize { args.iter().position(|a| a == name).and_then(|i| args.get(i + 1)).and_then(|v| v.parse().ok()).unwrap_or(def) };
    let (max_a, max_n) = (get("--max-a", 5), get("--max-n", 10));
    let mut pairs: Vec<(usize, usize)> = args.iter().filter_map(|s| s.split_once(':')).filter_map(|(x, y)| Some((x.parse().ok()?, y.parse().ok()?))).collect();
    if pairs.is_empty() {
        for a in 2..=max_a {
            for nn in a..=max_n {
                pairs.push((a, nn));
            }
        }
    }
    let mut certified = vec![];
    for (a, nn) in pairs {
        let t = Instant::now();
        let r = check(a, nn, direct, gens).unwrap_or_else(|e| {
            eprintln!("g4_free: {}", e);
            std::process::exit(2)
        });
        let expect = fact(nn) / fact(nn - a);
        let ok = r.certified();
        let mut line = format!(
            "a={} N={}: I_N*Delta=0 {} | staircase lead terms distinct ({} = N!/(N-a)!) {} | {}",
            a,
            nn,
            if r.annihilated { "yes" } else { "NO" },
            expect,
            if r.staircase_ok { "yes" } else { "NO" },
            if ok { "CERTIFIED free, rank a! = ".to_string() + &fact(a).to_string() } else { "NOT certified".to_string() }
        );
        if let Some(dim) = r.direct_dim {
            line += &format!(" | direct dim {} {}", dim, if dim == expect { "ok" } else { "MISMATCH" });
        }
        if let Some(g) = &r.gens {
            let l = a * (a - 1) / 2;
            let want: BTreeMap<usize, usize> = pi_a(a).into_iter().enumerate().filter(|&(_, c)| c > 0).map(|(k, c)| (k + l, c)).collect();
            let ngen: usize = g.values().sum();
            line += &format!(
                " | generators {} in x-degrees {:?} {} (#gens*C(N,a) = {})",
                ngen,
                g,
                if *g == want { "= t^l(w) pi_a ok" } else { "MISMATCH" },
                ngen * binom(nn, a)
            );
        }
        line += &format!(" | {:.2}s", t.elapsed().as_secs_f64());
        println!("{}", line);
        if ok {
            certified.push((a, nn));
        }
    }
    let list: Vec<String> = certified.iter().map(|(a, n)| format!("({}, {})", a, n)).collect();
    println!("certified: [{}]", list.join(", "));
}
