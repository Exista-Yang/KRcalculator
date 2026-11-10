//! Gate G4 (SYNTHESIS §3.1; C2 M1(b), U3): freeness of Δ_a R_N over A = H*(Gr(a, N)), the
//! hypothesis of the (Δ_λ G, F) adapter in src/cclose.rs.
//!
//! R_N = F[x_1..x_a]/(x_i^N), Δ = Π_{i<j} (x_i − x_j), F = F_p with p = 2^31 − 1,
//! A = Sym_a / I_N, I_N = (h_{N−a+1}, …, h_N).
//!
//! Certificate, exact, for every (a, N):
//!  (A) h_k · Δ ≡ 0 mod (x_1^N, …, x_a^N) for k = N−a+1..N. Then I_N F[x] · Δ ⊂ (x^N), so
//!      f ↦ Δf induces a surjection Q := F[x]/I_N F[x] → Δ R_N, and Δ R_N is an A-module.
//!      Q = F[x] ⊗_{Sym_a} A has dimension a!·C(N, a) = N!/(N−a)! because F[x] is free over
//!      Sym_a of rank a! (Artin basis; any field). Hence dim Δ R_N ≤ N!/(N−a)!.
//!  (B) For μ in the staircase S = {μ : μ_i ≤ N − a + i − 1} (|S| = N!/(N−a)!) the truncated
//!      products Δ·x^μ have pairwise distinct lex-leading monomials (namely x^{μ+δ},
//!      δ = (a−1, …, 0), coefficient 1, all exponents ≤ N − 1). Hence they are linearly
//!      independent and dim Δ R_N ≥ N!/(N−a)!.
//!  With (A) and (B), Q → Δ R_N is an isomorphism: Δ R_N ≅ F[x] ⊗_{Sym_a} A is free over A of
//!  rank a!, with generator degrees deg Δ + (Artin degrees), i.e. graded rank t^{l(w_a)} π_a(t)
//!  (t of x-degree 1). The checker computes (A) and (B) explicitly (it does not assume the
//!  leading-term claim, it evaluates every product).
//!
//! Cross-checks (C2's direct method, `--direct` / `--gens`, small cases):
//!  (C) dim Δ R_N by rank computations degree by degree = N!/(N−a)!;
//!  (D) minimal number of A-module generators per degree, dim M_D − dim (A_+ M)_D, equals the
//!      coefficients of t^{l(w_a)} π_a(t), and #generators · C(N, a) = dim Δ R_N.
//!
//! Usage: g4_free [--direct] [--gens] [--max-a A] [--max-n N] [a:N ...]
//! Default: all 2 ≤ a ≤ 5, a ≤ N ≤ 10, certificate only. Output: one line per pair, and the
//! certified list in the form of `cclose::G4_CERTIFIED`.

use krhom::field::{Field, F31};
use krhom::linsys::RankEchelon;
use krhom::sparse::SpVec;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::time::Instant;

type Mono = Vec<u8>;
type Poly = HashMap<Mono, F31>;

fn pmul(a: &Poly, b: &Poly, cap: Option<u8>) -> Poly {
    let mut r: Poly = HashMap::new();
    for (ma, ca) in a {
        for (mb, cb) in b {
            let m: Mono = ma.iter().zip(mb.iter()).map(|(x, y)| x + y).collect();
            if let Some(c) = cap {
                if m.iter().any(|&e| e > c) {
                    continue;
                }
            }
            r.entry(m).or_insert_with(F31::zero).fma(ca, cb);
        }
    }
    r.retain(|_, c| !c.is_zero());
    r
}

fn delta(a: usize) -> Poly {
    let mut d: Poly = HashMap::from([(vec![0u8; a], F31::one())]);
    for i in 0..a {
        for j in (i + 1)..a {
            let mut ei = vec![0u8; a];
            ei[i] = 1;
            let mut ej = vec![0u8; a];
            ej[j] = 1;
            let f: Poly = HashMap::from([(ei, F31::one()), (ej, F31::one().neg())]);
            d = pmul(&d, &f, None);
        }
    }
    d
}

/// All exponent vectors of length a with total degree k and entries ≤ cap.
fn monos(a: usize, k: usize, cap: usize) -> Vec<Mono> {
    fn rec(i: usize, left: usize, cap: usize, cur: &mut Mono, out: &mut Vec<Mono>) {
        if i == cur.len() {
            if left == 0 {
                out.push(cur.clone());
            }
            return;
        }
        for e in 0..=left.min(cap) {
            cur[i] = e as u8;
            rec(i + 1, left - e, cap, cur, out);
        }
        cur[i] = 0;
    }
    let mut out = vec![];
    rec(0, k, cap, &mut vec![0; a], &mut out);
    out
}

fn h_k(a: usize, k: usize) -> Poly {
    monos(a, k, k).into_iter().map(|m| (m, F31::one())).collect()
}

fn e_j(a: usize, j: usize) -> Poly {
    monos(a, j, 1).into_iter().map(|m| (m, F31::one())).collect()
}

fn binom(n: usize, k: usize) -> usize {
    if k > n {
        return 0;
    }
    (0..k).fold(1usize, |r, i| r * (n - i) / (i + 1))
}

fn fact(n: usize) -> usize {
    (1..=n).product()
}

/// Coefficients of π_a(t) = Π_{i=1}^{a} (1 + t + … + t^{i−1}).
fn pi_a(a: usize) -> Vec<usize> {
    let mut p = vec![1usize];
    for i in 1..=a {
        let mut q = vec![0usize; p.len() + i - 1];
        for (k, &c) in p.iter().enumerate() {
            for j in 0..i {
                q[k + j] += c;
            }
        }
        p = q;
    }
    p
}

struct Report {
    annihilated: bool,
    staircase_ok: bool,
    direct_dim: Option<usize>,
    gens: Option<BTreeMap<usize, usize>>,
}

fn check(a: usize, nn: usize, direct: bool, gens: bool) -> Report {
    let cap = (nn - 1) as u8;
    let d = delta(a);
    let l = a * (a - 1) / 2;
    // (A) h_k Δ ≡ 0 mod x^N for k = N−a+1..N
    let annihilated = ((nn - a + 1)..=nn).all(|k| pmul(&h_k(a, k), &d, Some(cap)).is_empty());
    // (B) staircase leading terms
    let mut leads: HashSet<Mono> = HashSet::new();
    let mut staircase_ok = true;
    let mut count = 0usize;
    let mut mu = vec![0u8; a];
    let lim: Vec<u8> = (0..a).map(|i| (nn - a + i) as u8).collect(); // μ_i ≤ N − a + i (0-based i)
    'outer: loop {
        let xm: Poly = HashMap::from([(mu.clone(), F31::one())]);
        let p = pmul(&d, &xm, Some(cap));
        let lead = p.iter().max_by(|x, y| x.0.cmp(y.0));
        let want: Mono = (0..a).map(|i| mu[i] + (a - 1 - i) as u8).collect();
        match lead {
            Some((m, c)) if *m == want && c.is_one() => {}
            _ => staircase_ok = false,
        }
        if !leads.insert(want) {
            staircase_ok = false;
        }
        count += 1;
        let mut i = 0;
        loop {
            if i == a {
                break 'outer;
            }
            if mu[i] < lim[i] {
                mu[i] += 1;
                break;
            }
            mu[i] = 0;
            i += 1;
        }
    }
    if count != fact(nn) / fact(nn - a) {
        staircase_ok = false;
    }
    let mut direct_dim = None;
    let mut gen_out = None;
    if direct || gens {
        // (C) rank of multiplication by Δ, degree by degree; keep a basis of M_D
        let maxdeg = a * (nn - 1);
        let mut basis: Vec<Vec<Poly>> = vec![vec![]; maxdeg + 1];
        let mut colmap: Vec<HashMap<Mono, u32>> = vec![HashMap::new(); maxdeg + 1];
        for (dd, cm) in colmap.iter_mut().enumerate() {
            let mut ms = monos(a, dd, nn - 1);
            ms.sort_by(|x, y| y.cmp(x)); // lex descending: pivots at leading terms
            for (i, m) in ms.into_iter().enumerate() {
                cm.insert(m, i as u32);
            }
        }
        let tovec = |p: &Poly, cm: &HashMap<Mono, u32>| -> SpVec<F31> {
            let mut v: SpVec<F31> = p.iter().map(|(m, c)| (cm[m], *c)).collect();
            v.sort_unstable_by_key(|e| e.0);
            v
        };
        let mut total = 0usize;
        for k in 0..=maxdeg {
            let dd = k + l;
            if dd > maxdeg {
                break;
            }
            let mut ech = RankEchelon::<F31>::new(colmap[dd].len().max(1));
            for m in monos(a, k, nn - 1) {
                let p = pmul(&d, &HashMap::from([(m, F31::one())]), Some(cap));
                if p.is_empty() {
                    continue;
                }
                if ech.insert(&tovec(&p, &colmap[dd])) {
                    basis[dd].push(p);
                }
            }
            total += ech.rank();
        }
        direct_dim = Some(total);
        if gens {
            // (D) generators per degree: dim M_D − rank of Σ_j e_j M_{D−j}
            let es: Vec<Poly> = (1..=a).map(|j| e_j(a, j)).collect();
            let mut g: BTreeMap<usize, usize> = BTreeMap::new();
            for dd in 0..=maxdeg {
                let dim = basis[dd].len();
                if dim == 0 {
                    continue;
                }
                let mut ech = RankEchelon::<F31>::new(colmap[dd].len().max(1));
                for j in 1..=a {
                    if dd < j {
                        continue;
                    }
                    for b in &basis[dd - j] {
                        let p = pmul(b, &es[j - 1], Some(cap));
                        if !p.is_empty() {
                            ech.insert(&tovec(&p, &colmap[dd]));
                        }
                    }
                }
                if dim > ech.rank() {
                    g.insert(dd, dim - ech.rank());
                }
            }
            gen_out = Some(g);
        }
    }
    Report { annihilated, staircase_ok, direct_dim, gens: gen_out }
}

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
        assert!(a >= 1 && a <= nn && a <= 8 && nn <= 64, "unsupported (a, N) = ({}, {})", a, nn);
        let t = Instant::now();
        let r = check(a, nn, direct, gens);
        let expect = fact(nn) / fact(nn - a);
        let ok = r.annihilated && r.staircase_ok;
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
