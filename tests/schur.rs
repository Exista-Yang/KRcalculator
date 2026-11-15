//! Tests of the decategorified oracle `schur.rs` (milestone M3a, SYNTHESIS §8.5).
//! Heavy cases are `#[ignore]`; run them one at a time with `--ignored`.

use krhom::input::crossing::{crossing_cx, qshift, Gen, QConv};
use krhom::combinat::hecke::Hecke;
use krhom::algebra::laurent::LPoly;
use krhom::combinat::schur::*;
use krhom::combinat::web::{Dir, Rung};
use std::collections::{BTreeMap, HashMap};
use std::path::PathBuf;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/schur").join(name)
}
fn tid(s: &mut Schur, txt: &str) -> Tid {
    s.intern(Table::parse(txt))
}
fn lp(terms: &[(i32, i128)]) -> Lp {
    let mut r = Lp::zero();
    for &(e, k) in terms {
        r.add_scaled(&Lp::one(), k, e);
    }
    r
}
fn qint(n: i64) -> Lp {
    Lp::from_lpoly(&LPoly::qint(n))
}
fn comp(s: &str) -> Vec<u8> {
    s.split(',').map(|x| x.trim().parse().unwrap()).collect()
}
fn word(s: &str) -> Vec<i32> {
    s.split_whitespace().map(|x| x.parse().unwrap()).collect()
}
/// element as a map table-string -> coefficient (independent of interning order)
fn emap(s: &Schur, e: &Elt) -> BTreeMap<String, Lp> {
    e.t.iter().map(|(p, c)| (s.table(*p).fmt(), c.clone())).collect()
}
fn multmap(s: &Schur, v: &[(Tid, Lp)]) -> BTreeMap<String, Lp> {
    v.iter().map(|(p, c)| (s.table(*p).fmt(), c.clone())).collect()
}
/// "E0:1 F0:1" (homcases.txt format)
fn parse_hom_rungs(s: &str) -> Vec<Rung> {
    s.split_whitespace()
        .map(|t| {
            let dir = if t.starts_with('E') { Dir::E } else { Dir::F };
            let (i, r) = t[1..].split_once(':').unwrap();
            Rung { i: i.parse().unwrap(), dir, r: r.parse().unwrap() }
        })
        .collect()
}
/// "1E1 0F2" (LAD format)
fn parse_lad(s: &str) -> Vec<Rung> {
    s.split_whitespace()
        .map(|t| {
            let p = t.find(['E', 'F']).unwrap();
            let dir = if &t[p..p + 1] == "E" { Dir::E } else { Dir::F };
            Rung { i: t[..p].parse().unwrap(), dir, r: t[p + 1..].parse().unwrap() }
        })
        .collect()
}

// ---------------------------------------------------------------------------------------------
// T-a

#[test]
fn t_a_table_counts() {
    let cases: &[(&[u8], usize)] = &[
        (&[1, 1, 1], 6),
        (&[1, 1, 1, 1], 24),
        (&[2, 2], 3),
        (&[3, 3], 4),
        (&[4, 4], 5),
        (&[5, 5], 6),
        (&[2, 1, 2], 11),
        (&[2, 2, 2], 21),
        (&[3, 3, 3], 55),
        (&[2, 2, 2, 2], 282),
        (&[5, 5, 5], 231),
        (&[3, 3, 3, 3], 2008),
        (&[2, 2, 2, 2, 2], 6210),
        (&[4, 4, 4, 4], 10147),
        (&[5, 5, 5, 5], 40176),
    ];
    for &(l, n) in cases {
        let ts = tables(l, l);
        assert_eq!(ts.len(), n, "number of tables of {:?}", l);
        if n <= 2008 {
            for t in &ts {
                assert_eq!(t.rows(), l.to_vec());
                assert_eq!(t.cols(), l.to_vec());
            }
            assert!(ts.windows(2).all(|w| (w[0].lmax(), &w[0].e) < (w[1].lmax(), &w[1].e)), "order of {:?}", l);
        }
    }
    // mixed margins
    assert_eq!(tables(&[2, 1], &[1, 2]).len(), 2);
    assert_eq!(tables(&[2, 1], &[2, 1]).len(), 2);
}

// ---------------------------------------------------------------------------------------------
// T-b

#[test]
fn t_b_worked_examples() {
    let mut s = Schur::new(&[2, 2]);
    let p0 = tid(&mut s, "2,0/0,2");
    let p1 = tid(&mut s, "1,1/1,1");
    let p2 = tid(&mut s, "0,2/2,0");
    let tb = |s: &Schur, p: Tid| (s.table(p).lmin(), s.table(p).lmax(), s.table(p).ncosets());
    assert_eq!(tb(&s, p0), (0, 2, 1));
    assert_eq!(tb(&s, p1), (1, 5, 4));
    assert_eq!(tb(&s, p2), (4, 6, 1));
    let b = s.canonical_basis(&[2, 2]).b.clone();
    assert_eq!(b[&p0].t, vec![(p0, Lp::one())]);
    assert_eq!(b[&p1].coef(p1), Lp::one());
    assert_eq!(b[&p1].coef(p0), lp(&[(1, 1), (3, 1)]));
    assert_eq!(b[&p1].t.len(), 2);
    assert_eq!(b[&p2].coef(p2), Lp::one());
    assert_eq!(b[&p2].coef(p1), lp(&[(1, 1)]));
    assert_eq!(b[&p2].coef(p0), lp(&[(4, 1)]));
    assert_eq!(b[&p2].t.len(), 3);
    assert_eq!([s.rank(&[2, 2], p0), s.rank(&[2, 2], p1), s.rank(&[2, 2], p2)], [1, 6, 6]);
    // Hom dimensions of (2,2), degrees -2..4 (R1 §1.6)
    let want: [[[i128; 7]; 3]; 3] = [
        [[0, 0, 1, 0, 2, 0, 5], [0, 0, 0, 1, 0, 3, 0], [0, 0, 0, 0, 0, 0, 1]],
        [[0, 0, 0, 1, 0, 3, 0], [0, 0, 1, 0, 5, 0, 14], [0, 0, 0, 1, 0, 4, 0]],
        [[0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 1, 0, 4, 0], [0, 0, 1, 0, 3, 0, 9]],
    ];
    let ps = [p0, p1, p2];
    for a in 0..3 {
        for c in 0..3 {
            assert_eq!(s.hom_gd(&b[&ps[a]], &b[&ps[c]], -2, 4), want[a][c].to_vec(), "Hom(B_p{}, B_p{})", a, c);
        }
    }

    // (2,1) -> (1,2)
    let mut s = Schur::new(&[2, 1]);
    let c0 = tid(&mut s, "1,1/0,1");
    let c1 = tid(&mut s, "0,2/1,0");
    let bc = s.canonical_basis(&[1, 2]).b.clone();
    assert_eq!(bc[&c0].t, vec![(c0, Lp::one())]);
    assert_eq!(bc[&c1].coef(c1), Lp::one());
    assert_eq!(bc[&c1].coef(c0), lp(&[(1, 1)]));
    assert_eq!(bc[&c1].t.len(), 2);
    assert_eq!([s.rank(&[1, 2], c0), s.rank(&[1, 2], c1)], [2, 3]);
    assert_eq!(s.hom_gd(&bc[&c0], &bc[&c0], -2, 4), vec![0, 0, 1, 0, 3, 0, 6]);
    assert_eq!(s.hom_gd(&bc[&c0], &bc[&c1], -2, 4), vec![0, 0, 0, 1, 0, 3, 0]);
    assert_eq!(s.hom_gd(&bc[&c1], &bc[&c1], -2, 4), vec![0, 0, 1, 0, 3, 0, 7]);
    // (2,1) -> (2,1)
    let d0 = tid(&mut s, "2,0/0,1");
    let d1 = tid(&mut s, "1,1/1,0");
    let bd = s.canonical_basis(&[2, 1]).b.clone();
    assert_eq!(bd[&d1].coef(d1), Lp::one());
    assert_eq!(bd[&d1].coef(d0), lp(&[(2, 1)]));
    assert_eq!(bd[&d1].t.len(), 2);
    assert_eq!(s.rank(&[2, 1], d1), 3);
    assert_eq!(s.rank(&[2, 1], d0), 1);
}

#[test]
fn split_merge_rung_consistency() {
    // a rung equals split followed by merge, as in schur.py
    let mut s = Schur::new(&[2, 1, 2]);
    let x = s.ladder_char(&[Rung::f(0, 1), Rung::e(1, 2)]);
    let a = s.rung(&x, &Rung::e(0, 2));
    let y = s.split_col(&x, 1, &[2, x.mu[1] - 2]);
    let y = s.merge_cols(&y, 0);
    assert_eq!(emap(&s, &a), emap(&s, &y));
    let a = s.rung(&x, &Rung::f(1, 1));
    let y = s.split_col(&x, 1, &[x.mu[1] - 1, 1]);
    let y = s.merge_cols(&y, 2);
    assert_eq!(emap(&s, &a), emap(&s, &y));
    // Lp basics
    let p = lp(&[(-3, 2), (0, 1), (2, -1)]);
    assert_eq!(p.bar(), lp(&[(3, 2), (0, 1), (-2, -1)]));
    assert_eq!(p.sym_part(), lp(&[(-3, 2), (3, 2), (0, 1)]));
    assert_eq!(Lp::from_lpoly(&p.to_lpoly()), p);
    assert_eq!(Lp::parse_ec(&p.fmt_ec()), p);
    assert_eq!(p.at1(), 2);
    assert_eq!((p.min_exp(), p.max_exp()), (Some(-3), Some(2)));
    assert_eq!(qint(2).mul(&qint(3)), lp(&[(-3, 1), (-1, 2), (1, 2), (3, 1)]));
    // pi_LR
    assert_eq!(Table::parse("2,1/0,3").pi_lr(), LPoly::qfact(2).mul(&LPoly::qfact(3)));
}

// ---------------------------------------------------------------------------------------------
// T-c, T-c2

struct CbRef {
    b: BTreeMap<String, BTreeMap<String, BTreeMap<String, Lp>>>, // lam -> p -> q -> coefficient
    lad: BTreeMap<String, BTreeMap<String, Vec<Rung>>>,
}
fn load_cb_ref() -> Option<CbRef> {
    let text = std::fs::read_to_string(data("cb_ref.txt")).ok()?;
    let mut r = CbRef { b: BTreeMap::new(), lad: BTreeMap::new() };
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let (tag, rest) = line.split_once(' ').unwrap();
        let f: Vec<&str> = rest.split('|').collect();
        match tag {
            "CB" => {
                let mut m = BTreeMap::new();
                for ent in &f[2..] {
                    let (q, poly) = ent.trim().split_once(' ').unwrap();
                    m.insert(q.to_string(), Lp::parse_ec(poly));
                }
                r.b.entry(f[0].to_string()).or_default().insert(f[1].to_string(), m);
            }
            "LAD" => {
                r.lad.entry(f[0].to_string()).or_default().insert(f[1].to_string(), parse_lad(f[2]));
            }
            _ => panic!("unknown fixture tag {}", tag),
        }
    }
    Some(r)
}

/// (sum of ranks, kappa, sum flows bfs, sum flows minflow)
fn check_reference_basis(r: &CbRef, lam: &[u8]) -> (u128, i128, u128, u128) {
    let key = lam.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",");
    let refb = &r.b[&key];
    let mut s = Schur::new(lam);
    let t0 = std::time::Instant::now();
    let cb = s.canonical_basis(lam);
    let order = cb.order.clone();
    let b = cb.b.clone();
    let lb = cb.ladder_bfs.clone();
    let lm = cb.ladder_minflow.clone();
    let dt = t0.elapsed().as_secs_f64();
    assert_eq!(order.len(), refb.len(), "number of b_p for {:?}", lam);
    let mut entries = 0;
    for &p in &order {
        let pname = s.table(p).fmt();
        assert_eq!(emap(&s, &b[&p]), refb[&pname], "b_p differs from cb.pkl for {:?} p = {}", lam, pname);
        entries += b[&p].t.len();
    }
    // the BFS ladders are those of schur.py
    let mut lad_same = 0;
    for &p in &order {
        if lb[&p] == r.lad[&key][&s.table(p).fmt()] {
            lad_same += 1;
        }
    }
    assert_eq!(lad_same, order.len(), "BFS ladders differ from schur.py for {:?}", lam);
    // T-c2: the canonical basis does not depend on the ladder that is peeled
    let bm = s.basis_from_ladders(lam, &lm);
    for &p in &order {
        assert_eq!(bm[&p], b[&p], "T-c2: minimal-flow ladder gives a different b_p for {:?}", lam);
    }
    let ranks: u128 = order.iter().map(|&p| s.rank(lam, p)).sum();
    let kappa = s.kappa(lam);
    let fb: u128 = order.iter().map(|p| flows(lam, &lb[p])).sum();
    let fm: u128 = order.iter().map(|p| flows(lam, &lm[p])).sum();
    let longest = order.iter().map(|p| lb[p].len()).max().unwrap();
    println!(
        "T-c {:?}: {} types, {} entries, sum ranks {}, kappa {}, flows bfs {} minflow {}, longest bfs ladder {}, basis time {:.2}s",
        lam,
        order.len(),
        entries,
        ranks,
        kappa,
        fb,
        fm,
        longest,
        dt
    );
    (ranks, kappa, fb, fm)
}

#[test]
fn t_c_reference_bases() {
    let r = load_cb_ref().expect("tests/data/schur/cb_ref.txt missing");
    let cases: &[(&[u8], u128, i128)] = &[
        (&[1, 1, 1], 19, 1),
        (&[1, 1, 1, 1], 219, 2),
        (&[2, 2], 13, 2),
        (&[3, 3], 63, 3),
        (&[2, 2, 2], 1045, 8),
        (&[3, 3, 3], 76399, 54),
        (&[2, 1, 2], 199, 4),
    ];
    for &(lam, ranks, kappa) in cases {
        let (rk, ka, _fb, fm) = check_reference_basis(&r, lam);
        assert_eq!((rk, ka), (ranks, kappa), "aggregates of {:?}", lam);
        if lam == [2, 2, 2] && fm > 1171 {
            println!("NOTE: minflow sum for (2,2,2) is {} > 1171 (R1's ladders)", fm);
        }
        if lam == [3, 3, 3] && fm > 105787 {
            println!("NOTE: minflow sum for (3,3,3) is {} > 105787 (R1's ladders)", fm);
        }
    }
}

#[test]
#[ignore]
fn t_c_reference_basis_2222() {
    let r = match load_cb_ref() {
        Some(r) if r.b.contains_key("2,2,2,2") => r,
        _ => {
            println!("cb_ref.txt has no (2,2,2,2) entry: skipped");
            return;
        }
    };
    let (rk, ka, _fb, fm) = check_reference_basis(&r, &[2, 2, 2, 2]);
    assert_eq!((rk, ka), (397687, 63));
    if fm > 519337 {
        println!("NOTE: minflow sum for (2,2,2,2) is {} > 519337 (R1's ladders)", fm);
    }
    // T-j: max rank 25 920 and its argmax; rank of the minimal-flow ladder of that table
    let mut s = Schur::new(&[2, 2, 2, 2]);
    let order = s.canonical_basis(&[2, 2, 2, 2]).order.clone();
    let (mut best, mut arg) = (0, 0);
    for &p in &order {
        let rk = s.rank(&[2, 2, 2, 2], p);
        if rk > best {
            best = rk;
            arg = p;
        }
    }
    assert_eq!(best, 25920);
    assert_eq!(s.table(arg).fmt(), "0,0,1,1/0,1,0,1/1,0,1,0/1,1,0,0");
    let pmax = *order.last().unwrap();
    assert_eq!(s.rank(&[2, 2, 2, 2], pmax), 2520);
    let lm = s.canonical_basis(&[2, 2, 2, 2]).ladder_minflow[&arg].clone();
    let lb = s.canonical_basis(&[2, 2, 2, 2]).ladder_bfs[&arg].clone();
    println!(
        "T-c (2,2,2,2): rank-25920 table: flows(minflow ladder) = {} ({:?}), flows(bfs ladder) = {}",
        flows(&[2, 2, 2, 2], &lm),
        lm.iter().map(|r| format!("{}{:?}{}", r.i, r.dir, r.r)).collect::<Vec<_>>(),
        flows(&[2, 2, 2, 2], &lb)
    );
}

// ---------------------------------------------------------------------------------------------
// T-d

#[test]
fn t_d_hom_formula() {
    let text = std::fs::read_to_string(data("hom60.txt")).expect("hom60.txt missing");
    let mut n = 0;
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let rest = line.strip_prefix("HOM ").unwrap();
        let f: Vec<&str> = rest.split('|').collect();
        let bottom = comp(f[0]);
        let src = parse_hom_rungs(f[1]);
        let tgt = parse_hom_rungs(f[2]);
        let dlo: i32 = f[3].trim().parse().unwrap();
        let dhi: i32 = f[4].trim().parse().unwrap();
        let mut want = vec![0i128; (dhi - dlo + 1) as usize];
        for t in f[5].split_whitespace() {
            let (d, k) = t.split_once(':').unwrap();
            want[(d.parse::<i32>().unwrap() - dlo) as usize] = k.parse().unwrap();
        }
        let mut s = Schur::new(&bottom);
        let xa = s.ladder_char(&src);
        let xb = s.ladder_char(&tgt);
        assert_eq!(s.hom_gd(&xa, &xb, dlo, dhi), want, "Hom formula differs for {}", line);
        n += 1;
    }
    assert_eq!(n, 60);
}

// ---------------------------------------------------------------------------------------------
// T-e, T-e2, T-f

/// one-line permutation of a permutation table: w[k] = row of the 1 in column k
fn perm_of(t: &Table) -> Vec<u8> {
    (0..t.nc as usize).map(|k| (0..t.nr as usize).find(|&i| t.get(i, k) == 1).unwrap() as u8).collect()
}

#[test]
fn t_e_uncolored_is_kl() {
    for n in [3usize, 4] {
        let h = Hecke::new(n);
        let lam = vec![1u8; n];
        let mut s = Schur::new(&lam);
        let cb = s.canonical_basis(&lam);
        let order = cb.order.clone();
        let b = cb.b.clone();
        for &p in &order {
            let x = h.g.index[&perm_of(s.table(p))];
            for &q in &order {
                let y = h.g.index[&perm_of(s.table(q))];
                assert_eq!(b[&p].coef(q).to_lpoly(), h.h(y, x), "c_q(b_p) != h_(y,x) for n = {}", n);
            }
            let r: i128 = (0..h.g.size()).map(|y| h.h(y, x).eval_at_one()).sum();
            assert_eq!(s.rank(&lam, p) as i128, r);
        }
    }
}

fn blocks_of(c: &[u8]) -> Vec<usize> {
    let mut b = vec![];
    for (i, &a) in c.iter().enumerate() {
        b.extend(std::iter::repeat(i).take(a as usize));
    }
    b
}
/// transfer.py table_of
fn table_of(w: &[u8], lam: &[u8], mu: &[u8]) -> Table {
    let (bl, bm) = (blocks_of(lam), blocks_of(mu));
    let mut e = vec![0u8; lam.len() * mu.len()];
    for k in 0..w.len() {
        e[bl[w[k] as usize] * mu.len() + bm[k]] += 1;
    }
    Table { nr: lam.len() as u8, nc: mu.len() as u8, e }
}
fn pmax_map(h: &Hecke, lam: &[u8], mu: &[u8]) -> HashMap<Table, usize> {
    let mut best: HashMap<Table, usize> = HashMap::new();
    for (x, w) in h.g.elems.iter().enumerate() {
        let t = table_of(w, lam, mu);
        let e = best.entry(t).or_insert(x);
        if h.g.len[x] > h.g.len[*e] {
            *e = x;
        }
    }
    best
}
fn fact(n: u8) -> u128 {
    (1..=n as u128).product()
}

#[test]
fn t_e2_singular_is_kl_of_pmax() {
    let mut hs: HashMap<usize, Hecke> = HashMap::new();
    let cases: &[(&[u8], &[u8])] = &[(&[2, 2], &[2, 2]), (&[2, 1, 2], &[2, 1, 2]), (&[3, 3], &[3, 3]), (&[2, 2, 2], &[2, 2, 2]), (&[2, 1], &[1, 2])];
    for &(lam, mu) in cases {
        let k: usize = lam.iter().map(|&a| a as usize).sum();
        let h = hs.entry(k).or_insert_with(|| Hecke::new(k));
        let pm = pmax_map(h, lam, mu);
        let mut s = Schur::new(lam);
        let cb = s.canonical_basis(mu);
        let order = cb.order.clone();
        let b = cb.b.clone();
        assert_eq!(pm.len(), order.len());
        let wmu: u128 = mu.iter().map(|&a| fact(a)).product();
        for &p in &order {
            let xmax = pm[s.table(p)];
            assert_eq!(h.g.len[xmax] as i64, s.table(p).lmax(), "l(p_max) for {}", s.table(p).fmt());
            for &r in &order {
                let rmax = pm[s.table(r)];
                assert_eq!(b[&p].coef(r).to_lpoly(), h.h(rmax, xmax), "c_r(b_p) != h_(rmax,pmax) in {:?}->{:?}", lam, mu);
            }
            let tot: i128 = (0..h.g.size()).map(|y| h.h(y, xmax).eval_at_one()).sum();
            assert_eq!(s.rank(mu, p) * wmu, tot as u128, "rank * |W_mu| for {}", s.table(p).fmt());
        }
    }
}

/// T-f (optional): Hom transfer F1, gd_k Hom(B_pmax, B_qmax) = pi_lam pi_mu gd_k Hom(B_p, B_q), k <= 10
#[test]
fn t_f_hom_transfer() {
    const D: usize = 10;
    fn series(p: &LPoly) -> Vec<i128> {
        let mut s = vec![0i128; D + 1];
        for (&e, &k) in &p.c {
            assert!(e >= 0);
            if (e as usize) <= D {
                s[e as usize] += k;
            }
        }
        s
    }
    fn smul(a: &[i128], b: &[i128]) -> Vec<i128> {
        let mut r = vec![0i128; D + 1];
        for i in 0..=D {
            for j in 0..=D - i {
                r[i + j] += a[i] * b[j];
            }
        }
        r
    }
    fn pi_comp(c: &[u8]) -> LPoly {
        let mut r = LPoly::one();
        for &a in c {
            for j in 1..=a as i32 {
                let mut f = LPoly::zero();
                for t in 0..j {
                    f.add_term(2 * t, 1);
                }
                r = r.mul(&f);
            }
        }
        r
    }
    let mut hs: HashMap<usize, Hecke> = HashMap::new();
    let cases: &[(&[u8], &[u8])] = &[
        (&[1, 1, 1], &[1, 1, 1]),
        (&[2, 2], &[2, 2]),
        (&[1, 1, 1, 1], &[1, 1, 1, 1]),
        (&[2, 1, 2], &[2, 1, 2]),
        (&[3, 3], &[3, 3]),
        (&[2, 2, 2], &[2, 2, 2]),
        (&[2, 1], &[1, 2]),
        (&[3, 1], &[1, 3]),
        (&[2, 1, 1], &[1, 2, 1]),
        (&[3, 2], &[2, 3]),
        (&[2, 2], &[1, 3]),
        (&[2, 1, 2], &[1, 2, 2]),
    ];
    let mut npairs = 0;
    for &(lam, mu) in cases {
        let k: usize = lam.iter().map(|&a| a as usize).sum();
        let h = hs.entry(k).or_insert_with(|| Hecke::new(k));
        let pm = pmax_map(h, lam, mu);
        let mut s = Schur::new(lam);
        let cb = s.canonical_basis(mu);
        let order = cb.order.clone();
        let b = cb.b.clone();
        let pil = series(&pi_comp(lam).mul(&pi_comp(mu)));
        let mut g_r = vec![0i128; D + 1];
        g_r[0] = 1;
        for _ in 0..k {
            for d in 2..=D {
                g_r[d] += g_r[d - 2];
            }
        }
        for &p in &order {
            for &q in &order {
                let (x, y) = (pm[s.table(p)], pm[s.table(q)]);
                let ord = smul(&series(&h.hom_grk(&h.kl[x], &h.kl[y])), &g_r);
                let neg = s.hom_gd(&b[&p], &b[&q], -12, -1);
                assert!(neg.iter().all(|&a| a == 0), "Hom^{{<0}} != 0");
                let sing = s.hom_gd(&b[&p], &b[&q], 0, D as i32);
                assert_eq!(ord, smul(&sing, &pil), "Hom transfer fails for {:?}->{:?}", lam, mu);
                npairs += 1;
            }
        }
    }
    println!("T-f: Hom transfer verified on {} pairs", npairs);
}

/// a * b_s in the standard basis (Lp coefficients)
fn mul_bs_lp(h: &Hecke, a: &[Lp], s: usize) -> Vec<Lp> {
    let mut out = vec![Lp::zero(); a.len()];
    for (y, c) in a.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        let ys = h.g.rmul[y][s];
        out[ys].add_scaled(c, 1, 0);
        out[y].add_scaled(c, 1, if h.g.len[ys] > h.g.len[y] { 1 } else { -1 });
    }
    out
}
/// a * b where b = sum_z b_z H_z; a H_z built along right descents (H_s = b_s - v), only for z
/// in the support of b, which is closed under taking prefixes when b is a KL element
fn hecke_mul_lp(h: &Hecke, a: &[Lp], b: &[Lp]) -> Vec<Lp> {
    fn ensure(h: &Hecke, z: usize, memo: &mut Vec<Option<Vec<Lp>>>) {
        if memo[z].is_some() {
            return;
        }
        let s = (0..h.g.n - 1).find(|&i| h.g.right_descent(z, i)).unwrap();
        let zp = h.g.rmul[z][s];
        ensure(h, zp, memo);
        let prev = memo[zp].as_ref().unwrap();
        let mut r = mul_bs_lp(h, prev, s);
        for (y, c) in prev.iter().enumerate() {
            r[y].add_scaled(c, -1, 1);
        }
        memo[z] = Some(r);
    }
    let n = a.len();
    let mut memo: Vec<Option<Vec<Lp>>> = vec![None; n];
    memo[0] = Some(a.to_vec());
    let mut out = vec![Lp::zero(); n];
    for z in 0..n {
        if b[z].is_zero() {
            continue;
        }
        ensure(h, z, &mut memo);
        for (w, c) in memo[z].as_ref().unwrap().iter().enumerate() {
            if !c.is_zero() {
                out[w] = out[w].add(&c.mul(&b[z]));
            }
        }
    }
    out
}
fn to_kl_lp(kl: &[Vec<Lp>], mut a: Vec<Lp>) -> Vec<Lp> {
    let n = a.len();
    let mut out = vec![Lp::zero(); n];
    for x in (0..n).rev() {
        if a[x].is_zero() {
            continue;
        }
        let c = a[x].clone();
        for z in 0..n {
            if !kl[x][z].is_zero() {
                a[z] = a[z].sub(&kl[x][z].mul(&c));
            }
        }
        out[x] = c;
    }
    out
}

/// T-g (optional): the 71 products of F2: b_pmax b_Wmax = [lambda]! * sum_r m_r b_rmax,
/// where b_p * [W] = sum_r m_r b_r in the Schur algebroid (transfer.py check_products).
#[test]
fn t_g_products() {
    let cases: Vec<(Vec<u8>, Vec<Vec<Rung>>)> = vec![
        (vec![2, 2], vec![vec![Rung::e(0, 1), Rung::f(0, 1)], vec![Rung::e(0, 2), Rung::f(0, 2)]]),
        (vec![3, 3], (1..=3).map(|k| vec![Rung::e(0, k), Rung::f(0, k)]).collect()),
        (vec![2, 2, 2], vec![vec![Rung::e(0, 1), Rung::f(0, 1)], vec![Rung::e(1, 2), Rung::f(1, 2)]]),
        (vec![2, 1, 2], vec![vec![Rung::e(0, 1), Rung::f(0, 1)]]),
    ];
    let mut count = 0;
    let mut hs: HashMap<usize, (Hecke, Vec<Vec<Lp>>)> = HashMap::new();
    for (lam, crossings) in cases {
        let k: usize = lam.iter().map(|&a| a as usize).sum();
        let (h, kl) = hs.entry(k).or_insert_with(|| {
            let h = Hecke::new(k);
            let kl = h.kl.iter().map(|e| e.iter().map(Lp::from_lpoly).collect()).collect();
            (h, kl)
        });
        let pm = pmax_map(h, &lam, &lam);
        let mut s = Schur::new(&lam);
        let order = s.canonical_basis(&lam).order.clone();
        let mut qf = Lp::one();
        for &a in &lam {
            qf = qf.mul(&Lp::from_lpoly(&LPoly::qfact(a as i64)));
        }
        for rs in &crossings {
            let w = s.ladder_char(rs);
            let ew = s.expand(&w);
            assert!(ew.len() == 1 && ew[0].1 == Lp::one(), "crossing term not indecomposable");
            let y = pm[s.table(ew[0].0)];
            for &p in &order {
                let sing = s.term_mult(p, rs);
                let x = pm[s.table(p)];
                let got = to_kl_lp(kl, hecke_mul_lp(h, &kl[x], &kl[y]));
                let mut want = vec![Lp::zero(); got.len()];
                for (r, m) in &sing {
                    want[pm[s.table(*r)]] = qf.mul(m);
                }
                assert_eq!(got, want, "product b_pmax * b_Wmax for {:?}, p = {}", lam, s.table(p).fmt());
                count += 1;
            }
        }
    }
    assert_eq!(count, 71);
}

// ---------------------------------------------------------------------------------------------
// T-h, T-h2

fn binom(n: u128, k: u128) -> u128 {
    if k > n {
        return 0;
    }
    (0..k).fold(1u128, |r, i| r * (n - i) / (i + 1))
}

#[test]
fn t_h_rickard_rigidity() {
    for c in 1..=5u8 {
        for d in 1..=5u8 {
            for pos in [true, false] {
                let mut s = Schur::new(&[c, d]);
                let terms = crossing_terms(c, d, pos, 0);
                let xs: Vec<Elt> = terms.iter().map(|t| s.ladder_char(&t.rungs)).collect();
                let mut ranks = vec![];
                for (k, x) in xs.iter().enumerate() {
                    let e = s.expand(x);
                    assert!(e.len() == 1 && e[0].1 == Lp::one(), "term {} of ({},{}) {} not indecomposable", k, c, d, pos);
                    let end = s.hom_gd(x, x, -60, 0);
                    assert!(end[..60].iter().all(|&a| a == 0) && end[60] == 1, "End of term {} of ({},{})", k, c, d);
                    let p = e[0].0;
                    ranks.push(s.rank(&x.mu, p));
                    if pos {
                        let k8 = k as u8;
                        let want = if c >= d { format!("{},{}/{},{}", d - k8, c - d + k8, k8, d - k8) } else { format!("{},{}/{},{}", c - k8, k8, d - c + k8, c - k8) };
                        assert_eq!(s.table(p).fmt(), want, "term table of ({},{}) k = {}", c, d, k);
                    }
                }
                for k in 0..xs.len().saturating_sub(1) {
                    let (a, b) = if pos { (&xs[k], &xs[k + 1]) } else { (&xs[k + 1], &xs[k]) };
                    let hm = s.hom_gd(a, b, -60, 1);
                    assert!(hm[..61].iter().all(|&v| v == 0) && hm[61] == 1, "differential of ({},{}) {} at k = {}: {:?}", c, d, pos, k, &hm[55..]);
                }
                for k in 0..xs.len().saturating_sub(2) {
                    let (a, b) = if pos { (&xs[k], &xs[k + 2]) } else { (&xs[k + 2], &xs[k]) };
                    assert_eq!(s.hom_gd(a, b, 2, 2), vec![0], "Hom^2(W_k, W_k+2) for ({},{})", c, d);
                }
                if c == d {
                    let a = c as u128;
                    let want: Vec<u128> = (0..=a).map(|k| binom(a, k) * binom(a + k, k)).collect();
                    assert_eq!(ranks, want, "ranks for a = {}", c);
                }
            }
        }
    }
    let mut s = Schur::new(&[5, 5]);
    let rk: Vec<u128> = crossing_terms(5, 5, true, 0)
        .iter()
        .map(|t| {
            let x = s.ladder_char(&t.rungs);
            let p = s.expand(&x)[0].0;
            s.rank(&[5, 5], p)
        })
        .collect();
    assert_eq!(rk, vec![1, 30, 210, 560, 630, 252]);
}

#[test]
fn t_h2_two_rung_ladders() {
    for c in 1..=5u8 {
        for d in 1..=5u8 {
            let mut s = Schur::new(&[c, d]);
            for mu in [vec![c, d], vec![d, c]] {
                let b = s.canonical_basis(&mu).b.clone();
                // all ladders with at most two rungs ending at mu
                let mut chars: Vec<Elt> = vec![];
                let lab0 = [c as usize, d as usize];
                let mut ladders: Vec<Vec<Rung>> = vec![vec![]];
                for d1 in [Dir::E, Dir::F] {
                    for r1 in 1..=(c + d) as usize {
                        let a = Rung { i: 0, dir: d1, r: r1 };
                        let Some(l1) = a.apply(&lab0, 64) else { continue };
                        ladders.push(vec![a]);
                        for d2 in [Dir::E, Dir::F] {
                            for r2 in 1..=(c + d) as usize {
                                let bb = Rung { i: 0, dir: d2, r: r2 };
                                if bb.apply(&l1, 64).is_some() {
                                    ladders.push(vec![a, bb]);
                                }
                            }
                        }
                    }
                }
                for l in &ladders {
                    let x = s.ladder_char(l);
                    if x.mu == mu {
                        chars.push(x);
                    }
                }
                for (p, bp) in &b {
                    assert!(chars.iter().any(|x| x == bp), "b_{} of ({},{}) -> {:?} is not a 2-rung ladder", s.table(*p).fmt(), c, d, mu);
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------------------------
// T-i

#[test]
fn t_i_decompositions() {
    let mut s = Schur::new(&[2, 2]);
    let p1 = tid(&mut s, "1,1/1,1");
    let p2 = tid(&mut s, "0,2/2,0");
    let w1 = [Rung::e(0, 1), Rung::f(0, 1)];
    let w2 = [Rung::e(0, 2), Rung::f(0, 2)];
    let two_sq = qint(2).mul(&qint(2));
    let mut want = BTreeMap::new();
    want.insert("0,2/2,0".to_string(), two_sq);
    want.insert("1,1/1,1".to_string(), qint(2));
    let got = s.term_mult(p1, &w1);
    assert_eq!(multmap(&s, &got), want);
    let got = s.term_mult(p1, &w2);
    assert_eq!(multmap(&s, &got), BTreeMap::from([("0,2/2,0".to_string(), lp(&[(-3, 1), (-1, 2), (1, 2), (3, 1)]))]));
    let got = s.term_mult(p2, &w2);
    assert_eq!(multmap(&s, &got), BTreeMap::from([("0,2/2,0".to_string(), Lp::from_lpoly(&LPoly::qbinom(4, 2)))]));

    // (2,1) with W'_k the terms of the positive (1,2) crossing
    let mut s = Schur::new(&[2, 1]);
    let c0 = tid(&mut s, "1,1/0,1");
    let c1 = tid(&mut s, "0,2/1,0");
    let wp: Vec<Vec<Rung>> = crossing_terms(1, 2, true, 0).into_iter().map(|t| t.rungs).collect();
    assert_eq!(wp[0], vec![Rung::e(0, 1)]);
    assert_eq!(wp[1], vec![Rung::f(0, 1), Rung::e(0, 2)]);
    let d0 = "2,0/0,1".to_string();
    let d1 = "1,1/1,0".to_string();
    let g = s.term_mult(c0, &wp[0]);
    assert_eq!(multmap(&s, &g), BTreeMap::from([(d1.clone(), Lp::one()), (d0, Lp::one())]));
    let g = s.term_mult(c0, &wp[1]);
    assert_eq!(multmap(&s, &g), BTreeMap::from([(d1.clone(), qint(2))]));
    let g = s.term_mult(c1, &wp[0]);
    assert_eq!(multmap(&s, &g), BTreeMap::from([(d1.clone(), qint(2))]));
    let g = s.term_mult(c1, &wp[1]);
    assert_eq!(multmap(&s, &g), BTreeMap::from([(d1, qint(3))]));

    // split: all 21 types of (2,2,2,0) times F_2^(1) stay indecomposable (t10.py, Williamson
    // Prop 7.4.3): a single canonical element times v^e with coefficient exactly 1 (a
    // multiplicity k > 1 would be k copies, i.e. decomposable)
    let lam = [2u8, 2, 2, 0];
    let mut s = Schur::new(&lam);
    let order = s.canonical_basis(&lam).order.clone();
    assert_eq!(order.len(), 21);
    for p in order {
        let e = s.term_mult(p, &[Rung::f(2, 1)]);
        assert!(e.len() == 1 && e[0].1.c == vec![1], "split of {} is not v^e times a single canonical element: {:?}", s.table(p).fmt(), multmap(&s, &e));
    }
}

// ---------------------------------------------------------------------------------------------
// T-j

#[test]
fn t_j_ranks_and_kappa() {
    let cases: &[(&[u8], u128, &str, u128)] = &[
        (&[2, 2, 2], 216, "0,1,1/1,0,1/1,1,0", 90),
        (&[3, 3, 3], 8820, "0,1,2/1,1,1/2,1,0", 1680),
        (&[2, 1, 2], 36, "", 30),
    ];
    for &(lam, maxr, arg, rpmax) in cases {
        let mut s = Schur::new(lam);
        let order = s.canonical_basis(lam).order.clone();
        let (mut best, mut am) = (0, 0);
        for &p in &order {
            let r = s.rank(lam, p);
            if r > best {
                best = r;
                am = p;
            }
        }
        assert_eq!(best, maxr, "max rank of {:?}", lam);
        if !arg.is_empty() {
            assert_eq!(s.table(am).fmt(), arg);
        }
        let pmax = *order.last().unwrap();
        assert_eq!(s.rank(lam, pmax), rpmax, "rank(B_pmax) for {:?}", lam);
        let k: u128 = lam.iter().map(|&a| a as u128).sum();
        assert_eq!(rpmax, (1..=k).product::<u128>() / lam.iter().map(|&a| fact(a)).product::<u128>());
    }
    let mut s = Schur::new(&[4, 4]);
    assert_eq!(s.kappa(&[4, 4]), 6);
    let order = s.canonical_basis(&[4, 4]).order.clone();
    let mut r: Vec<u128> = order.iter().map(|&p| s.rank(&[4, 4], p)).collect();
    r.sort();
    assert_eq!(r, vec![1, 20, 70, 90, 140]);
    let mut s = Schur::new(&[5, 5]);
    assert_eq!(s.kappa(&[5, 5]), 10);
}

// ---------------------------------------------------------------------------------------------
// T-k, T-l

struct ScanWant {
    unsimpl: u128,
    lb_obj: u128,
    lb_types: usize,
    lb_ldim: u128,
    lb_s: u128,
    max_rank: u128,
    peak: Option<(usize, u128, u128, u128)>, // (crossing, objects, ldim, sum S)
}

fn run_scan(lam: &[u8], w: &str, want: &ScanWant, marked: Option<usize>) -> (Schur, ScanResult) {
    let t0 = std::time::Instant::now();
    let mut s = Schur::new(lam);
    let opts = ScanOpts { peak: want.peak.is_some(), every_stage: false, marked_block: marked };
    let res = decat_scan(&mut s, &word(w), &opts);
    let f = res.stages.last().unwrap();
    println!(
        "T-k {:?} [{}]: unsimplified {} (ldim {}), LB objects {} types {} ldim {} sum S {} max rank {} kappa {} ({:.2}s)",
        lam,
        w,
        f.unsimplified_objects,
        f.unsimplified_ldim,
        f.lb_objects,
        f.lb_types,
        f.lb_ldim,
        f.lb_sum_s,
        f.max_rank_present,
        f.kappa_present,
        t0.elapsed().as_secs_f64()
    );
    assert_eq!(f.step, word(w).len());
    assert_eq!(
        (f.unsimplified_objects, f.lb_objects, f.lb_types, f.lb_ldim, f.lb_sum_s, f.max_rank_present),
        (want.unsimpl, want.lb_obj, want.lb_types, want.lb_ldim, want.lb_s, want.max_rank),
        "scan of {:?} [{}]",
        lam,
        w
    );
    assert_eq!(f.lb_objects, res.k0.iter().map(|e| e.2.unsigned_abs()).sum::<u128>());
    if let Some((cr, a, b, c)) = want.peak {
        // the stage after j crossings carries the peak of crossing j+1
        let (j, pk) = res
            .stages
            .iter()
            .filter_map(|st| st.peak_before_next.map(|p| (st.step + 1, p)))
            .max_by(|x, y| x.1.0.cmp(&y.1.0).then(y.0.cmp(&x.0)))
            .unwrap();
        println!("     peak at crossing {}: {:?}", j, pk);
        assert_eq!((j, pk), (cr, (a, b, c)), "peak of {:?} [{}]", lam, w);
    }
    (s, res)
}

#[test]
fn t_k_scans_default() {
    run_scan(&[2, 2], "1 1 1", &ScanWant { unsimpl: 367, lb_obj: 11, lb_types: 3, lb_ldim: 61, lb_s: 67, max_rank: 6, peak: Some((3, 81, 481, 495)) }, None);
    run_scan(&[3, 3], "1 1 1", &ScanWant { unsimpl: 12332, lb_obj: 24, lb_types: 4, lb_ldim: 507, lb_s: 693, max_rank: 30, peak: Some((3, 802, 16695, 18489)) }, None);
    run_scan(&[2, 2, 2], "1 -2 1 -2", &ScanWant { unsimpl: 533, lb_obj: 77, lb_types: 16, lb_ldim: 3049, lb_s: 6073, max_rank: 216, peak: Some((4, 325, 17641, 26733)) }, None);
    run_scan(&[1, 1, 1], "1 -2 1 -2", &ScanWant { unsimpl: 24, lb_obj: 12, lb_types: 5, lb_ldim: 33, lb_s: 33, max_rank: 4, peak: None }, None);
    run_scan(&[1, 1, 1, 1], "1 1 2 -3 2 1 -3 -2 -2 -3 -3", &ScanWant { unsimpl: 17248, lb_obj: 192, lb_types: 19, lb_ldim: 1439, lb_s: 1471, max_rank: 16, peak: Some((11, 434, 3591, 3623)) }, None);
    // extra rows of A.2 (decat_scan.py)
    run_scan(&[2, 2], "1 1 1 1 1", &ScanWant { unsimpl: 61883, lb_obj: 19, lb_types: 3, lb_ldim: 109, lb_s: 119, max_rank: 6, peak: None }, None);
    run_scan(&[3, 3], "1 1 1 1 1", &ScanWant { unsimpl: 49590996, lb_obj: 42, lb_types: 4, lb_ldim: 911, lb_s: 1253, max_rank: 30, peak: None }, None);
    // 4- and 5-colored rows of A.2 (exact values: schur.py + t8.py logic, rerun for M3a)
    run_scan(&[4, 4], "1 1 1", &ScanWant { unsimpl: 460_461, lb_obj: 47, lb_types: 5, lb_ldim: 4121, lb_s: 7379, max_rank: 140, peak: None }, None);
    run_scan(&[5, 5], "1 1 1", &ScanWant { unsimpl: 18_230_234, lb_obj: 82, lb_types: 6, lb_ldim: 32151, lb_s: 79421, max_rank: 630, peak: None }, None);
    run_scan(&[5, 5], "1 1 1 1 1", &ScanWant { unsimpl: 53_509_193_375_758, lb_obj: 196, lb_types: 6, lb_ldim: 74603, lb_s: 172_253, max_rank: 630, peak: None }, None);
    run_scan(&[4, 4], "1 1", &ScanWant { unsimpl: 1309, lb_obj: 21, lb_types: 5, lb_ldim: 1801, lb_s: 3409, max_rank: 140, peak: None }, None);
    run_scan(&[5, 5], "1 1", &ScanWant { unsimpl: 9373, lb_obj: 33, lb_types: 6, lb_ldim: 12325, lb_s: 31525, max_rank: 630, peak: None }, None);
    // cheap in Rust although listed as ignored in SYNTHESIS §8.5 (Python cost): run by default
    run_scan(&[3, 3, 3], "1 -2 1 -2", &ScanWant { unsimpl: 11748, lb_obj: 328, lb_types: 40, lb_ldim: 348153, lb_s: 1_984_617, max_rank: 7200, peak: Some((4, 5036, 6_663_069, 20_378_865)) }, None);
    run_scan(&[2, 2, 2], "1 2 1 2 1 2 1 2", &ScanWant { unsimpl: 9_151_101, lb_obj: 187, lb_types: 21, lb_ldim: 13933, lb_s: 27753, max_rank: 216, peak: None }, None);
    run_scan(&[3, 3, 3], "1 2 1 2 1 2 1 2", &ScanWant { unsimpl: 123_881_877_708, lb_obj: 1084, lb_types: 55, lb_ldim: 2_566_793, lb_s: 13_023_797, max_rank: 8820, peak: None }, None);
    run_scan(&[2, 2, 2], "1 1 -2 1 -2 1 1 -2 -2 1", &ScanWant { unsimpl: 1_620_828_641, lb_obj: 8835, lb_types: 21, lb_ldim: 546_553, lb_s: 969_741, max_rank: 216, peak: Some((10, 33369, 1_901_497, 2_897_749)) }, None);
}

#[test]
fn t_k_mixed_colors() {
    // mixed colors: stage blocks (lambda, sigma lambda) use their own canonical bases.
    // Exact values from an independent Python reference (schur.py canonical bases + t8.py
    // scan logic: review_M3a/scan_ref.py for the stages, fix_m3a_cl1/gram_ref.py for K_0 and
    // the Gram rows), so a wrong (lambda, sigma lambda) basis changes them.
    let mut s = Schur::new(&[2, 1, 1]);
    let res = decat_scan(&mut s, &word("1 2 2 1"), &ScanOpts { peak: true, every_stage: true, marked_block: Some(0) });
    type Row = (usize, Vec<u8>, [u128; 7], i128, Option<(u128, u128, u128)>);
    // (step, labels, [unsimplified objects, unsimplified ldim, LB objects, LB types, LB ldim,
    //  LB sum S, max rank present], kappa present, peak before next)
    let want: Vec<Row> = vec![
        (0, vec![2, 1, 1], [1, 1, 1, 1, 1, 1, 1], 1, Some((2, 5, 5))),
        (1, vec![1, 2, 1], [2, 5, 2, 2, 5, 5, 3], 1, Some((4, 25, 25))),
        (2, vec![1, 1, 2], [4, 25, 4, 4, 25, 25, 9], 1, Some((18, 125, 125))),
        (3, vec![1, 2, 1], [18, 125, 6, 4, 35, 35, 9], 1, Some((25, 175, 183))),
        (4, vec![2, 1, 1], [73, 625, 13, 6, 103, 111, 12], 2, None),
    ];
    let got: Vec<Row> = res
        .stages
        .iter()
        .map(|st| {
            (
                st.step,
                st.labels.clone(),
                [st.unsimplified_objects, st.unsimplified_ldim, st.lb_objects, st.lb_types as u128, st.lb_ldim, st.lb_sum_s, st.max_rank_present],
                st.kappa_present,
                st.peak_before_next,
            )
        })
        .collect();
    assert_eq!(got, want);
    // the K_0 class: (type, q-shift) -> signed multiplicity
    let k0: BTreeMap<(String, i32), i128> = res.k0.iter().map(|&(p, k, c)| ((s.table(p).fmt(), k), c)).collect();
    let k0_want: BTreeMap<(String, i32), i128> = [
        ("2,0,0/0,1,0/0,0,1", 0, 1),
        ("1,1,0/1,0,0/0,0,1", -6, 1),
        ("1,1,0/1,0,0/0,0,1", -2, -1),
        ("1,0,1/1,0,0/0,1,0", -5, -1),
        ("1,0,1/1,0,0/0,1,0", -3, 1),
        ("1,1,0/0,0,1/1,0,0", -5, -1),
        ("1,1,0/0,0,1/1,0,0", -3, 1),
        ("1,0,1/0,1,0/1,0,0", -4, 1),
        ("1,0,1/0,1,0/1,0,0", -2, -1),
        ("0,1,1/1,0,0/1,0,0", -7, 1),
        ("0,1,1/1,0,0/1,0,0", -5, -2),
        ("0,1,1/1,0,0/1,0,0", -3, 1),
    ]
    .iter()
    .map(|&(p, k, c)| ((p.to_string(), k), c))
    .collect();
    assert_eq!(k0, k0_want);
    // Gram rows in order (lmax, table): (p, weight, unknowns, unknowns reduced (block 0 marked))
    let gram: Vec<(String, i32, u128, u128)> = res.gram.iter().map(|g| (s.table(g.p).fmt(), g.weight, g.unknowns, g.unknowns_reduced)).collect();
    let gram_want: Vec<(String, i32, u128, u128)> = [
        ("2,0,0/0,1,0/0,0,1", 0, 1, 1),
        ("1,1,0/1,0,0/0,0,1", 2, 7, 3),
        ("1,0,1/1,0,0/0,1,0", 3, 13, 4),
        ("1,1,0/0,0,1/1,0,0", 3, 13, 4),
        ("1,0,1/0,1,0/1,0,0", 4, 22, 5),
        ("0,1,1/1,0,0/1,0,0", 5, 34, 6),
    ]
    .iter()
    .map(|&(p, w, u, r)| (p.to_string(), w, u, r))
    .collect();
    assert_eq!(gram, gram_want);
    for g in &res.gram {
        assert_eq!(g.unknowns_reduced, emono_count(&[1, 1], g.weight as usize));
    }
}

#[test]
fn t_k_every_stage_and_k0() {
    // every_stage records steps 0..=len; the K_0 class is the signed multiset of the final stage
    let mut s = Schur::new(&[2, 1]);
    let res = decat_scan(&mut s, &word("1 1"), &ScanOpts { peak: false, every_stage: true, marked_block: None });
    assert_eq!(res.stages.iter().map(|st| st.step).collect::<Vec<_>>(), vec![0, 1, 2]);
    assert_eq!(res.stages[1].labels, vec![1, 2]);
    assert_eq!(res.stages[0].lb_objects, 1);
    assert!(res.k0.iter().all(|e| e.2 != 0));
}

#[test]
#[ignore]
fn t_k_scan_conway_2colored() {
    let t0 = std::time::Instant::now();
    let (mut s, res) = run_scan(
        &[2, 2, 2, 2],
        "1 1 2 -3 2 1 -3 -2 -2 -3 -3",
        // exact values: decat_scan.py conway2 and peak.py conway (rerun for M3a)
        &ScanWant { unsimpl: 1_468_708_995, lb_obj: 11029, lb_types: 235, lb_ldim: 13_282_393, lb_s: 66_632_927, max_rank: 25920, peak: Some((11, 103_573, 103_035_829, 339_702_323)) },
        Some(0),
    );
    // R1 §6.4: lower bound after 4 / 6 / 8 / 10 crossings
    let lbs: Vec<u128> = [4usize, 6, 8, 10].iter().map(|&j| res.stages.iter().find(|st| st.step == j).unwrap().lb_objects).collect();
    assert_eq!(lbs, vec![63, 707, 2695, 7167]);
    // T-l: Gram rows
    assert_eq!(res.gram.len(), 235);
    let maxw = res.gram.iter().map(|g| g.weight).max().unwrap();
    let maxu = res.gram.iter().map(|g| g.unknowns).max().unwrap();
    let maxr = res.gram.iter().map(|g| g.unknowns_reduced).max().unwrap();
    assert_eq!((maxw, maxu, maxr), (22, 159_952, 13_377));
    let mut hist: BTreeMap<i32, usize> = BTreeMap::new();
    for g in &res.gram {
        *hist.entry(g.weight).or_default() += 1;
    }
    println!("Conway Gram histogram (weight: #types): {:?}", hist);
    println!("2-colored Conway scan with peaks: {:.1}s", t0.elapsed().as_secs_f64());
    let _ = s.kappa(&[2, 2, 2, 2]);
}

#[test]
fn t_l_closure_weights() {
    let cases: &[(&[u8], usize, u128)] = &[
        (&[2, 2], 4, 14),
        (&[3, 3], 9, 158),
        (&[2, 2, 2], 12, 1218),
        (&[2, 2, 2, 2], 24, 260_780),
        (&[3, 3, 3], 27, 299_958),
        (&[5, 5], 25, 40_466),
        (&[2, 2, 2, 2], 22, 159_952),
        (&[2, 2, 2], 22, 13_377),
    ];
    for &(l, w, n) in cases {
        assert_eq!(emono_count(l, w), n, "emono_count({:?}, {})", l, w);
    }
    // A.4 rows
    assert_eq!((1..=4).map(|w| emono_count(&[2, 2], w)).collect::<Vec<_>>(), vec![2, 5, 8, 14]);
    assert_eq!((1..=8).map(|w| emono_count(&[2, 2, 2, 2], w)).collect::<Vec<_>>(), vec![4, 14, 36, 85, 176, 344, 624, 1086]);
    assert_eq!(emono_count(&[2, 1, 1, 1], 9), 420);
    assert_eq!(emono_count(&[2, 2, 1, 1], 13), 2856);
    for (lam, want) in [(vec![2u8, 2], vec![0, 3, 4]), (vec![3u8, 3], vec![0, 5, 8, 9])] {
        let mut s = Schur::new(&lam);
        let res = decat_scan(&mut s, &word("1 1 1"), &ScanOpts { peak: false, every_stage: false, marked_block: Some(0) });
        let mut ws: Vec<i32> = res.gram.iter().map(|g| g.weight).collect();
        ws.sort();
        ws.dedup();
        assert_eq!(ws, want, "closure weights of the {:?} trefoil", lam);
        for g in &res.gram {
            assert_eq!(g.unknowns, emono_count(&lam, g.weight as usize));
            assert_eq!(g.unknowns_reduced, emono_count(&lam[1..], g.weight as usize));
        }
    }
}

// ---------------------------------------------------------------------------------------------
// T-m

#[test]
fn t_m_agrees_with_crossing_rs() {
    for c in 1..=5u8 {
        for d in 1..=5u8 {
            for pos in [true, false] {
                for i in [0usize, 2] {
                    let n = (c + d) as usize;
                    let cx = crossing_cx(Gen { i, positive: pos }, c as usize, d as usize, n);
                    let ours = crossing_terms(c, d, pos, i);
                    assert_eq!(ours.len(), cx.terms.len(), "number of terms ({},{}) {}", c, d, pos);
                    let eps = if pos { 1 } else { -1 };
                    for (t, ct) in ours.iter().zip(&cx.terms) {
                        assert_eq!(t.rungs, ct.rungs);
                        assert_eq!(t.hdeg, ct.hdeg);
                        assert_eq!(t.dq, qshift(&cx, ct.k, QConv::A) + eps * cx.x, "local q part ({},{}) {} k = {}", c, d, pos, ct.k);
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------------------------
// T-n

fn parse_chi(s: &str) -> LPoly {
    let mut p = LPoly::zero();
    for t in s.split_whitespace() {
        let (e, k) = t.split_once(':').unwrap();
        p.add_term(e.parse().unwrap(), k.parse().unwrap());
    }
    p
}

#[test]
fn t_n_decategorified_chi() {
    let text = std::fs::read_to_string(data("chi_ref.txt")).expect("chi_ref.txt missing");
    let (mut n_ok, mut n_alt_differs, mut n_a3) = (0, 0, 0);
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let (body, is_extra) = match line.split_once('#') {
            Some((b, _)) => (b.trim(), true),
            None => (line, false),
        };
        let f: Vec<&str> = body.strip_prefix("CHI ").unwrap().split('|').collect();
        let n: usize = f[0].trim().parse().unwrap();
        let labels = comp(f[1]);
        let w = word(f[2]);
        let want = parse_chi(f[3]);
        let mut s = Schur::new(&labels);
        let res = decat_scan(&mut s, &w, &ScanOpts { peak: false, every_stage: false, marked_block: None });
        let got = decat_chi(&mut s, &res, &w, n);
        let alt = decat_chi_conv(&mut s, &res, &w, n, true);
        println!("T-n N={} labels={:?} word={:?}: v->q {} ; v->q^-1 {}", n, labels, w, if got == want { "PASS" } else { "FAIL" }, if alt == want { "PASS" } else { "FAIL" });
        assert_eq!(got, want, "chi for N = {}, labels {:?}, word {:?} (v -> q^-1 gives {})", n, labels, w, alt);
        n_ok += 1;
        if alt != want {
            n_alt_differs += 1;
        }
        if !is_extra {
            n_a3 += 1;
        }
    }
    assert!(n_a3 >= 8, "chi_ref.txt must contain the 8 cases of A.3");
    assert!(n_alt_differs > 0, "T-n does not discriminate v -> q from v -> q^-1");
    println!("T-n: {} cases pass with v -> q; v -> q^-1 fails on {} of them", n_ok, n_alt_differs);
}

// ---------------------------------------------------------------------------------------------
// T-o

#[test]
fn t_o_cache_round_trip() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR"));
    for (lam, mu) in [(vec![2u8, 2, 2], vec![2u8, 2, 2]), (vec![2u8, 1], vec![1u8, 2])] {
        let path = dir.join(format!("schur_cb_{}_{}.txt", lam.len(), mu[0]));
        let mut s = Schur::new(&lam);
        let _ = s.canonical_basis(&mu);
        s.save_basis(&mu, &path).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        assert!(text.starts_with("schur-cb v1 lam="), "missing version header");
        let mut t = Schur::new(&lam);
        let _ = t.intern(Table::parse(if lam.len() == 3 { "0,1,1/1,0,1/1,1,0" } else { "0,2/1,0" })); // different interning order
        t.load_basis(&path).unwrap();
        type Snap = (Vec<String>, BTreeMap<String, BTreeMap<String, Lp>>, BTreeMap<String, (Vec<Rung>, Vec<Rung>)>);
        fn snapshot(s: &mut Schur, mu: &[u8]) -> Snap {
            let cb = s.canonical_basis(mu);
            let (order, b, lb, lm) = (cb.order.clone(), cb.b.clone(), cb.ladder_bfs.clone(), cb.ladder_minflow.clone());
            let o = order.iter().map(|&p| s.table(p).fmt()).collect();
            let bb = b.iter().map(|(p, e)| (s.table(*p).fmt(), emap(s, e))).collect();
            let ll = order.iter().map(|p| (s.table(*p).fmt(), (lb[p].clone(), lm[p].clone()))).collect();
            (o, bb, ll)
        }
        assert_eq!(snapshot(&mut s, &mu), snapshot(&mut t, &mu));
        // a wrong lambda is rejected
        let mut u = Schur::new(&[1, 1]);
        assert!(u.load_basis(&path).is_err());
    }
}

// ---------------------------------------------------------------------------------------------
// input validation, cache files and JSON records (review of M3a)

fn schur_scan(args: &[&str]) -> std::process::Output {
    std::process::Command::new(env!("CARGO_BIN_EXE_schur_scan")).args(args).output().expect("cannot run schur_scan")
}

#[test]
fn marked_block_out_of_range_is_an_error() {
    let lam = [2u8, 2];
    let w = word("1 1 1");
    let opts = |m: usize| ScanOpts { peak: false, every_stage: false, marked_block: Some(m) };
    assert!(check_scan(&lam, &w, &opts(0)).is_ok());
    assert!(check_scan(&lam, &w, &opts(1)).is_ok());
    for m in [2, 5] {
        let e = check_scan(&lam, &w, &opts(m)).unwrap_err();
        assert!(e.contains(&format!("marked block {} out of range", m)), "{}", e);
    }
    assert!(check_scan(&[2, 0], &[], &opts(1)).unwrap_err().contains("label 0"));
    assert!(check_scan(&lam, &word("1 -2"), &opts(0)).unwrap_err().contains("needs 3 strands"));
    assert!(check_scan(&lam, &[0], &opts(0)).is_err());
    // the CLI reports a usage error instead of unknowns_reduced = unknowns
    let out = schur_scan(&["scan", "--labels", "2,2", "--word", "1 1 1", "--marked", "5"]);
    assert_eq!(out.status.code(), Some(2), "stdout: {}", String::from_utf8_lossy(&out.stdout));
    assert!(String::from_utf8_lossy(&out.stderr).contains("marked block 5 out of range"));
    assert!(out.stdout.is_empty(), "no scan output expected");
}

#[test]
#[should_panic(expected = "decat_scan: marked block 5 out of range")]
fn decat_scan_refuses_marked_block_out_of_range() {
    let mut s = Schur::new(&[2, 2]);
    let _ = decat_scan(&mut s, &word("1 1 1"), &ScanOpts { peak: false, every_stage: false, marked_block: Some(5) });
}

#[test]
fn unreachable_blocks_are_refused_early() {
    assert!(check_block(&[2, 1], &[1, 2]).is_ok());
    assert!(check_block(&[2, 2, 2, 0], &[2, 2, 1, 1]).is_ok());
    assert!(check_block(&[2, 1], &[1, 1, 1]).unwrap_err().contains("as many entries"));
    assert!(check_block(&[2, 1], &[1, 1]).unwrap_err().contains("same sum"));
    // the CLI fails before any search
    let out = schur_scan(&["basis", "--lam", "2,1", "--mu", "1,1,1"]);
    assert_eq!(out.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&out.stderr).contains("as many entries"));
}

#[test]
#[should_panic(expected = "canonical_basis: block [2, 1] -> [1, 1, 1]: mu must have as many entries")]
fn canonical_basis_refuses_mu_of_other_length_before_searching() {
    // (the BFS would otherwise explore every reachable state and then fail its own assert)
    let mut s = Schur::new(&[2, 1]);
    let _ = s.canonical_basis(&[1, 1, 1]);
}

#[test]
fn load_basis_rejects_malformed_or_corrupted_files() {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR"));
    let lam = [2u8, 2];
    let src = dir.join("schur_cb_corrupt_src.txt");
    let mut s = Schur::new(&lam);
    let _ = s.canonical_basis(&lam);
    s.save_basis(&lam, &src).unwrap();
    let text = std::fs::read_to_string(&src).unwrap();
    Schur::new(&lam).load_basis(&src).expect("the untouched file loads");
    let line_of = |prefix: &str| text.lines().find(|l| l.starts_with(prefix)).unwrap_or_else(|| panic!("no line {:?}", prefix)).to_string();
    let lad_p1 = line_of("LADMF 2,2|1,1/1,1|").rsplit('|').next().unwrap().to_string();
    let i128_max = i128::MAX.to_string();
    // (what, fixture text, its replacement, a substring of the expected error): each corruption
    // must be rejected for its own reason, not by an earlier unrelated check
    let cases: Vec<(&str, String, String, &str)> = vec![
        // semantic corruptions that keep every coefficient in vZ[v] (accepted before)
        ("b_p1: (v + v^3) -> (v + 7 v^3)", "2,0/0,2 1:1 3:1".into(), "2,0/0,2 1:1 3:7".into(), "b_1,1/1,1 is not the canonical basis element"),
        ("b_p2: v^4 h_p0 -> v^2 h_p0 (expansions stay non-negative)", "2,0/0,2 4:1".into(), "2,0/0,2 2:1".into(), "b_0,2/2,0 is not the canonical basis element"),
        ("coefficient of h_p is 2", "|1,1/1,1 0:1".into(), "|1,1/1,1 0:2".into(), "b_1,1/1,1: coefficient of h_p is 0:2, not 1"),
        ("lower coefficient not in vZ[v]", "2,0/0,2 1:1 3:1".into(), "2,0/0,2 0:1 3:1".into(), "of h_2,0/0,2 is not in vZ[v]"),
        ("table outside the block", "|2,0/0,2 4:1".into(), "|3,0/0,1 4:1".into(), "table 3,0/0,1 is outside the block"),
        ("LADMF ladder with another leading term", line_of("LADMF 2,2|0,2/2,0|"), format!("LADMF 2,2|0,2/2,0|{}", lad_p1), "LADMF ladder of 0,2/2,0 does not have the leading term h_p"),
        ("LAD ladder ending outside mu", line_of("LAD 2,2|0,2/2,0|"), "LAD 2,2|0,2/2,0|0E1".into(), "ladder ends at [3, 1], not at mu"),
        ("missing CB line", format!("{}\n", line_of("CB 2,2|0,2/2,0|")), String::new(), "basis file misses table 0,2/2,0"),
        ("duplicate CB line", line_of("CB 2,2|2,0/0,2|"), format!("{}\n{}", line_of("CB 2,2|2,0/0,2|"), line_of("CB 2,2|2,0/0,2|")), "two CB lines for table 2,0/0,2"),
        ("mu of another length", "mu=2,2".into(), "mu=2,2,0".into(), "mu must have as many entries as lambda"),
        ("bad version", "schur-cb v1".into(), "schur-cb v2".into(), "bad header"),
        // malformed tokens (panicked before)
        ("table cell 2,x", "|2,0/0,2 4:1".into(), "|2,x/0,2 4:1".into(), "bad table \"2,x/0,2\""),
        ("ragged table", "|2,0/0,2 4:1".into(), "|2,0,0/0,2 4:1".into(), "ragged table"),
        ("bad coefficient", "2,0/0,2 4:1".into(), "2,0/0,2 4:z".into(), "bad coefficient in \"4:z\""),
        ("bad exponent", "2,0/0,2 4:1".into(), "2,0/0,2 q:1".into(), "bad exponent in \"q:1\""),
        ("exponent too large for dense storage", "2,0/0,2 4:1".into(), "2,0/0,2 5000000:1".into(), "exponent out of range [0, 6] in \"5000000:1\""),
        ("coefficient overflow", "2,0/0,2 4:1".into(), format!("2,0/0,2 4:{} 4:1", i128_max), "coefficient overflow at \"4:1\""),
        ("invalid rung", line_of("LAD 2,2|0,2/2,0|"), "LAD 2,2|0,2/2,0|0E9".into(), "invalid rung 0E9"),
        // crafted inputs of the review of the M3a fixes (panicked before)
        ("u8 row and column sums wrap around (128 + 130 = 2 mod 256)", "|2,0/0,2 4:1".into(), "|2,0/0,2 4:1|128,130/130,128 5:1".into(), "table 128,130/130,128 is outside the block"),
        ("repeated entries whose sum overflows i128", "2,0/0,2 1:1 3:1".into(), format!("2,0/0,2 1:1 3:1|2,0/0,2 1:{}", i128_max), "b_1,1/1,1: coefficient overflow in the repeated entries of h_2,0/0,2"),
        ("rung index 2^64 - 1", line_of("LAD 2,2|0,2/2,0|"), format!("LAD 2,2|0,2/2,0|{}E1", u64::MAX), "invalid rung"),
        ("exponent above the largest lmax of the block (6)", "2,0/0,2 4:1".into(), "2,0/0,2 7:1".into(), "exponent out of range [0, 6] in \"7:1\""),
        ("negative exponent", "2,0/0,2 4:1".into(), "2,0/0,2 -1:1".into(), "exponent out of range [0, 6] in \"-1:1\""),
        // a valid but long ladder: its character is computed with an i128 overflow pre-check
        ("long ladder whose character overflows i128", line_of("LADMF 2,2|1,1/1,1|"), format!("LADMF 2,2|1,1/1,1|{}", vec!["0E1 0F1"; 400].join(" ")), "LADMF ladder of 1,1/1,1: its character overflows i128 coefficients"),
    ];
    for (i, (what, from, to, want)) in cases.into_iter().enumerate() {
        assert!(text.contains(&from), "{}: fixture text {:?} not found", what, from);
        let p = dir.join(format!("schur_cb_corrupt_{}.txt", i));
        std::fs::write(&p, text.replacen(&from, &to, 1)).unwrap();
        let mut t = Schur::new(&lam);
        match std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| t.load_basis(&p))) {
            Ok(Err(e)) => {
                assert_eq!(e.kind(), std::io::ErrorKind::InvalidData, "{}", what);
                assert!(e.to_string().contains(want), "{}: rejected for another reason: {}", what, e);
            }
            Ok(Ok(())) => panic!("{}: corrupted basis accepted", what),
            Err(_) => panic!("{}: load_basis panicked instead of returning Err", what),
        }
    }
    // the CLI exits with an error (status 1) instead of a panic (101) on the crafted files
    let cli = dir.join("schur_cb_crafted_cli");
    let _ = std::fs::remove_dir_all(&cli);
    std::fs::create_dir_all(&cli).unwrap();
    for (from, to) in [("|2,0/0,2 4:1".to_string(), "|2,0/0,2 4:1|128,130/130,128 5:1".to_string()), ("2,0/0,2 1:1 3:1".into(), format!("2,0/0,2 1:1 3:1|2,0/0,2 1:{}", i128_max))] {
        std::fs::write(cli.join("cb_2-2_2-2.txt"), text.replacen(&from, &to, 1)).unwrap();
        let out = schur_scan(&["scan", "--labels", "2,2", "--word", "1 1 1", "--load", cli.to_str().unwrap()]);
        assert_eq!(out.status.code(), Some(1), "{:?}: {}", to, String::from_utf8_lossy(&out.stderr));
        assert!(String::from_utf8_lossy(&out.stderr).contains("schur_scan: error: load_basis"));
    }
}

/// value of `"key":` in a one-line JSON object written by schur_scan (string or number)
fn json_field<'a>(line: &'a str, key: &str) -> &'a str {
    let pat = format!("\"{}\":", key);
    let i = line.find(&pat).unwrap_or_else(|| panic!("no field {} in {}", key, line)) + pat.len();
    let rest = &line[i..];
    match rest.strip_prefix('"') {
        Some(r) => r.split('"').next().unwrap(),
        None => rest.split([',', '}']).next().unwrap(),
    }
}

#[test]
fn cli_json_has_gram_rows() {
    let out = schur_scan(&["scan", "--labels", "2,2", "--word", "1 1 1", "--marked", "0", "--json"]);
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let text = String::from_utf8(out.stdout).unwrap();
    let rows: Vec<(String, i32, u128, u128, u128)> = text
        .lines()
        .filter(|l| l.contains("\"kind\":\"gram_row\""))
        .map(|l| {
            assert!(l.contains("\"labels\":[2,2]"), "{}", l);
            assert_eq!(json_field(l, "word"), "1 1 1");
            let num = |k: &str| json_field(l, k).parse::<u128>().unwrap();
            (json_field(l, "p").to_string(), json_field(l, "weight").parse().unwrap(), num("objects"), num("unknowns"), num("unknowns_reduced"))
        })
        .collect();
    let mut s = Schur::new(&[2, 2]);
    let res = decat_scan(&mut s, &word("1 1 1"), &ScanOpts { peak: false, every_stage: false, marked_block: Some(0) });
    let mut objs: BTreeMap<Tid, u128> = BTreeMap::new();
    for &(p, _, c) in &res.k0 {
        *objs.entry(p).or_default() += c.unsigned_abs();
    }
    let want: Vec<(String, i32, u128, u128, u128)> = res.gram.iter().map(|g| (s.table(g.p).fmt(), g.weight, objs[&g.p], g.unknowns, g.unknowns_reduced)).collect();
    assert_eq!(rows, want);
    // the 2-colored trefoil: 3 LB types, 11 LB objects, weights {0, 3, 4} (T-k, T-l)
    assert_eq!(rows.len(), 3);
    assert_eq!(rows.iter().map(|r| r.2).sum::<u128>(), 11);
    assert_eq!(rows.iter().map(|r| r.1).collect::<Vec<_>>(), vec![0, 3, 4]);
    for r in &rows {
        assert_eq!((r.3, r.4), (emono_count(&[2, 2], r.1 as usize), emono_count(&[2], r.1 as usize)));
    }
}
