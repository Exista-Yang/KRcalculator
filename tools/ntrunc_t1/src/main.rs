//! ntrunc_t1: N-truncation study for exterior-colored gl_N KR homology (scratch crate, read-only use of krhom).
//!
//! t1 kill  LAM NMIN NMAX [moy]   agreement table: trace criterion vs RSK(col-reversed) vs RS(pmax)^T vs naive
//! t1 scan  LAM WORD N|inf        truncated decategorified scan (canonical coordinates), stats
//! t1 cmp   LAM WORD N            truncated scan + untruncated decat_scan + chi check
//! t1 rbfs  LAM N [targetonly]    restricted (labels <= N) BFS + peel, compared with the full basis
//! t1 count LAM N                 surviving type counts per block (RSK rule only, no basis)

use krhom::algebra::laurent::LPoly;
use krhom::combinat::schur::*;
use krhom::combinat::web::{Dir, Ladder, Rung};
use std::collections::{BTreeMap, BinaryHeap, HashMap, HashSet, VecDeque};
use std::time::Instant;

// ------------------------------------------------------------------------------------------------
// resources

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
fn peak_rss_mb() -> f64 {
    let mut ru = RUsage { utime: [0; 2], stime: [0; 2], maxrss: 0, rest: [0; 13] };
    let rc = unsafe { getrusage(0, &mut ru) };
    if rc != 0 {
        return 0.0;
    }
    // ru_maxrss is in bytes on macOS and in kilobytes on Linux
    if cfg!(target_os = "macos") { ru.maxrss as f64 / 1048576.0 } else { ru.maxrss as f64 / 1024.0 }
}

// ------------------------------------------------------------------------------------------------
// combinatorics

fn comps(k: usize, n: usize) -> Vec<Vec<u8>> {
    let mut out = vec![];
    fn rec(k: usize, n: usize, cur: &mut Vec<u8>, out: &mut Vec<Vec<u8>>) {
        if n == 1 {
            cur.push(k as u8);
            out.push(cur.clone());
            cur.pop();
            return;
        }
        for a in 0..=k {
            cur.push(a as u8);
            rec(k - a, n - 1, cur, out);
            cur.pop();
        }
    }
    rec(k, n, &mut vec![], &mut out);
    out
}

/// RSK row insertion of the second coordinates of the lexicographically sorted pairs; returns the shape
fn rsk_shape(mut pairs: Vec<(usize, usize)>) -> Vec<usize> {
    pairs.sort();
    let mut rows: Vec<Vec<usize>> = vec![];
    for &(_, x0) in &pairs {
        let mut x = x0;
        let mut placed = false;
        for row in rows.iter_mut() {
            if let Some(pos) = row.iter().position(|&y| y > x) {
                std::mem::swap(&mut row[pos], &mut x);
            } else {
                row.push(x);
                placed = true;
                break;
            }
        }
        if !placed {
            rows.push(vec![x]);
        }
    }
    rows.iter().map(|r| r.len()).collect()
}

fn transpose_shape(sh: &[usize]) -> Vec<usize> {
    if sh.is_empty() {
        return vec![];
    }
    (0..sh[0]).map(|c| sh.iter().filter(|&&r| r > c).count()).collect()
}

/// mu(p) = RSK shape of the table with its columns reversed
fn mu_colrev(t: &Table) -> Vec<usize> {
    let (n, m) = (t.nr as usize, t.nc as usize);
    let mut pairs = vec![];
    for i in 0..n {
        for j in 0..m {
            for _ in 0..t.get(i, j) {
                pairs.push((i, m - 1 - j));
            }
        }
    }
    rsk_shape(pairs)
}
/// RSK shape of the table itself (no reversal), for comparison
fn mu_plain(t: &Table) -> Vec<usize> {
    let (n, m) = (t.nr as usize, t.nc as usize);
    let mut pairs = vec![];
    for i in 0..n {
        for j in 0..m {
            for _ in 0..t.get(i, j) {
                pairs.push((i, j));
            }
        }
    }
    rsk_shape(pairs)
}

/// longest element of the double coset of table t, as a permutation (positions = bottom blocks,
/// values = top blocks); asserts that its length is lmax(t)
fn pmax_perm(t: &Table) -> Vec<usize> {
    let (n, m) = (t.nr as usize, t.nc as usize);
    let cols = t.cols();
    let mut toff = vec![0usize; m];
    for j in 1..m {
        toff[j] = toff[j - 1] + cols[j - 1] as usize;
    }
    let mut w = vec![];
    for i in 0..n {
        let mut vals = vec![];
        for j in 0..m {
            let below: usize = (i + 1..n).map(|i2| t.get(i2, j) as usize).sum();
            let st = toff[j] + below;
            for k in 0..t.get(i, j) as usize {
                vals.push(st + k);
            }
        }
        vals.sort_unstable_by(|a, b| b.cmp(a));
        w.extend(vals);
    }
    let mut inv = 0i64;
    for a in 0..w.len() {
        for b in a + 1..w.len() {
            if w[a] > w[b] {
                inv += 1;
            }
        }
    }
    assert_eq!(inv, t.lmax(), "pmax length != lmax for {}", t.fmt());
    w
}
fn mu_rs_pmax_t(t: &Table) -> Vec<usize> {
    let w = pmax_perm(t);
    let pairs: Vec<(usize, usize)> = w.iter().enumerate().map(|(i, &x)| (i, x)).collect();
    transpose_shape(&rsk_shape(pairs))
}

fn transpose_table(t: &Table) -> Table {
    let (n, m) = (t.nr as usize, t.nc as usize);
    let mut e = vec![0u8; n * m];
    for i in 0..n {
        for j in 0..m {
            e[j * n + i] = t.get(i, j);
        }
    }
    Table { nr: m as u8, nc: n as u8, e }
}

fn lw(c: &[u8]) -> i64 {
    c.iter().map(|&a| (a as i64) * (a as i64 - 1) / 2).sum()
}

/// coefficient list of prod_i [N choose a_i]_t (Gaussian binomials in t), index = weight
/// (= Hilbert series of (x)_i H*(Gr(a_i, N)) in units of degree 2)
fn gauss(n: usize, a: usize) -> Vec<u128> {
    // [n choose a]_t = [n-1 choose a-1]_t + t^a [n-1 choose a]_t
    let mut g: Vec<Vec<Vec<u128>>> = vec![vec![vec![1u128]]];
    for nn in 1..=n {
        let mut row = vec![];
        for aa in 0..=nn {
            let mut v = vec![0u128; aa * (nn - aa) + 1];
            if aa >= 1 {
                for (i, &x) in g[nn - 1][aa - 1].iter().enumerate() {
                    v[i] += x;
                }
            }
            if aa <= nn - 1 {
                for (i, &x) in g[nn - 1][aa].iter().enumerate() {
                    v[i + aa] += x;
                }
            }
            row.push(v);
        }
        g.push(row);
    }
    g[n][a].clone()
}
fn box_series(lam: &[u8], n: usize) -> Vec<u128> {
    let mut ser = vec![1u128];
    for &a in lam {
        if a as usize > n {
            return vec![0];
        }
        ser = conv(&ser, &gauss(n, a as usize));
    }
    ser
}
fn conv(a: &[u128], b: &[u128]) -> Vec<u128> {
    let mut r = vec![0u128; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        for (j, &y) in b.iter().enumerate() {
            r[i + j] += x * y;
        }
    }
    r
}

// ------------------------------------------------------------------------------------------------
// traces

fn closure_checked(lam: &[u8], rs: &[Rung], n: usize) -> LPoly {
    let mut lab: Vec<usize> = lam.iter().map(|&a| a as usize).collect();
    if lab.iter().any(|&a| a > n) {
        return LPoly::zero();
    }
    for r in rs {
        match r.apply(&lab, n) {
            Some(l) => lab = l,
            None => return LPoly::zero(),
        }
    }
    let lad = Ladder { n, bottom: lam.iter().map(|&a| a as usize).collect(), rungs: rs.to_vec() };
    krhom::combinat::moy::closure_value(&lad)
}

/// tr_N(b_p), p in (lam, lam): copy of schur::trace_b (v -> q)
fn trace_b(s: &mut Schur, p: Tid, n: usize, memo: &mut HashMap<Tid, LPoly>) -> LPoly {
    if let Some(v) = memo.get(&p) {
        return v.clone();
    }
    let lam = s.lam.clone();
    let lad = s.canonical_basis(&lam).ladder_minflow[&p].clone();
    let ch = s.ladder_char(&lad);
    let mult = s.expand(&ch);
    let mut v = closure_checked(&lam, &lad, n);
    for (q, m) in mult {
        if q == p {
            assert_eq!(m, Lp::one());
            continue;
        }
        let tq = trace_b(s, q, n, memo);
        v = v.sub(&m.to_lpoly().mul(&tq));
    }
    memo.insert(p, v.clone());
    v
}

fn ladder_max_label(lam: &[u8], rs: &[Rung]) -> usize {
    let mut lab: Vec<usize> = lam.iter().map(|&a| a as usize).collect();
    let mut mx = *lab.iter().max().unwrap();
    for r in rs {
        let (rec, giv) = match r.dir {
            Dir::E => (r.i, r.i + 1),
            Dir::F => (r.i + 1, r.i),
        };
        lab[giv] -= r.r;
        lab[rec] += r.r;
        mx = mx.max(lab[rec]);
    }
    mx
}

fn lp_mul(a: &Lp, b: &Lp) -> Lp {
    a.mul(b)
}

/// b_y as a Z[v,v^-1]-combination of the minimal-flow ladders L_s of its block: map s -> coeff
fn ladder_combos(s: &mut Schur, mu: &[u8]) -> HashMap<Tid, HashMap<Tid, Lp>> {
    let cb = s.canonical_basis(mu);
    let order = cb.order.clone();
    let lads = cb.ladder_minflow.clone();
    let mut out: HashMap<Tid, HashMap<Tid, Lp>> = HashMap::new();
    for &y in &order {
        let ch = s.ladder_char(&lads[&y]);
        let e = s.expand(&ch);
        let mut c: HashMap<Tid, Lp> = HashMap::new();
        c.insert(y, Lp::one());
        for (q, m) in e {
            if q == y {
                continue;
            }
            for (sid, a) in &out[&q] {
                let ent = c.entry(*sid).or_insert_with(Lp::zero);
                *ent = ent.sub(&lp_mul(&m, a));
            }
        }
        c.retain(|_, v| !v.is_zero());
        out.insert(y, c);
    }
    out
}

fn elt_scale(x: &Elt, a: &Lp) -> Elt {
    let mut r = Elt::zero(&x.mu);
    for (i, &k) in a.c.iter().enumerate() {
        if k != 0 {
            r = r.add_scaled(x, k, a.lo + i as i32);
        }
    }
    r
}

fn lpoly_nonneg(p: &LPoly) -> bool {
    p.c.values().all(|&k| k > 0)
}

// ------------------------------------------------------------------------------------------------
// kill table

fn cmd_kill(lam: Vec<u8>, nmin: usize, nmax: usize, moy_check: bool) {
    let t0 = Instant::now();
    let k: usize = lam.iter().map(|&a| a as usize).sum();
    let nparts = lam.len();
    let mus = comps(k, nparts);
    let mut s = Schur::new(&lam);
    s.canonical_basis(&lam);
    // rows: N -> counters
    #[derive(Default, Clone)]
    struct Cnt {
        types: usize,
        killed_trace: usize,
        killed_rsk: usize,
        mism_rsk: usize,
        mism_rspmax: usize,
        mism_plain: usize,
        killed_naive: usize,
        naive_false_neg: usize,
        naive_false_pos: usize,
        diag_mism: usize,
        moy_checked: usize,
        moy_mism: usize,
        tr_nonpos: usize,
        types_nt: usize,
        killed_nt: usize,
    }
    let mut cnt: BTreeMap<usize, Cnt> = BTreeMap::new();
    let mut examples: Vec<String> = vec![];
    let mut rsk_vs_rs = 0usize;
    let mut maxlab_ne_mu1 = 0usize;
    let mut maxlab_ne_examples: Vec<String> = vec![];
    let mut ntypes_total = 0usize;
    // per-N trace memo for (lam, lam)
    let mut memos: BTreeMap<usize, HashMap<Tid, LPoly>> = (nmin..=nmax).map(|n| (n, HashMap::new())).collect();
    // check positivity of tr_N(b_r)
    let lam_order = s.canonical_basis(&lam).order.clone();
    for n in nmin..=nmax {
        let memo = memos.get_mut(&n).unwrap();
        for &r in &lam_order {
            let t = trace_b(&mut s, r, n, memo);
            if !lpoly_nonneg(&t) {
                cnt.entry(n).or_default().tr_nonpos += 1;
            }
        }
    }
    eprintln!("traces of (lam,lam) types done: {:.1}s", t0.elapsed().as_secs_f64());
    for mu in &mus {
        let cb = s.canonical_basis(mu);
        let order = cb.order.clone();
        let bmap: HashMap<Tid, Elt> = order.iter().map(|p| (*p, cb.b[p].clone())).collect();
        let mf: HashMap<Tid, Vec<Rung>> = cb.ladder_minflow.clone();
        let mut sm = Schur::new(mu);
        let cbm = sm.canonical_basis(&lam);
        let order_m = cbm.order.clone();
        let lads_m: HashMap<Tid, Vec<Rung>> = cbm.ladder_minflow.clone();
        let combos_m = ladder_combos(&mut sm, &lam);
        let combos_l = if moy_check { Some(ladder_combos(&mut s, mu)) } else { None };
        // expansions E[p][s] of b_p * L_s in (lam, lam)
        let mut expn: HashMap<(Tid, Tid), Vec<(Tid, Lp)>> = HashMap::new();
        for &p in &order {
            for &sid in &order_m {
                let prod = s.mul_rungs(&bmap[&p], &lads_m[&sid]);
                let e = s.expand(&prod);
                expn.insert((p, sid), e);
            }
        }
        ntypes_total += order.len();
        for &p in &order {
            let t = s.table(p).clone();
            let mu_c = mu_colrev(&t);
            let mu_r = mu_rs_pmax_t(&t);
            let mu_p = mu_plain(&t);
            if mu_c != mu_r {
                rsk_vs_rs += 1;
            }
            let pt = transpose_table(&t);
            let pt_id = sm.intern(pt);
            let maxlab = ladder_max_label(&lam, &mf[&p]);
            if maxlab != mu_c[0] {
                maxlab_ne_mu1 += 1;
                if maxlab_ne_examples.len() < 10 {
                    maxlab_ne_examples.push(format!("{} mu={:?} maxlab={}", t.fmt(), mu_c, maxlab));
                }
            }
            for n in nmin..=nmax {
                let memo = memos.get_mut(&n).unwrap();
                let mut tvals: HashMap<Tid, LPoly> = HashMap::new();
                let mut any_nonzero = false;
                for &sid in &order_m {
                    let mut tot = LPoly::zero();
                    for (r, m) in &expn[&(p, sid)] {
                        let tr = trace_b(&mut s, *r, n, memo);
                        tot = tot.add(&m.to_lpoly().mul(&tr));
                    }
                    if !tot.is_zero() {
                        any_nonzero = true;
                    }
                    tvals.insert(sid, tot);
                }
                let killed_trace = !any_nonzero;
                // <b_p, b_p> = tr(b_p b_{p^T})
                let mut diag = LPoly::zero();
                for (sid, a) in &combos_m[&pt_id] {
                    diag = diag.add(&a.to_lpoly().mul(&tvals[sid]));
                }
                let killed_diag = diag.is_zero();
                let killed_rsk = mu_c[0] > n;
                let killed_rs = mu_r[0] > n;
                let killed_plain = mu_p[0] > n;
                let killed_naive = maxlab > n;
                let c = cnt.entry(n).or_default();
                c.types += 1;
                if mu.iter().all(|&a| a as usize <= n) {
                    c.types_nt += 1;
                    if killed_trace {
                        c.killed_nt += 1;
                    }
                }
                if killed_trace {
                    c.killed_trace += 1;
                }
                if killed_rsk {
                    c.killed_rsk += 1;
                }
                if killed_rsk != killed_trace {
                    c.mism_rsk += 1;
                    if examples.len() < 20 {
                        examples.push(format!("RSK MISMATCH N={} p={} mu_colrev={:?} trace_killed={}", n, t.fmt(), mu_c, killed_trace));
                    }
                }
                if killed_rs != killed_trace {
                    c.mism_rspmax += 1;
                }
                if killed_plain != killed_trace {
                    c.mism_plain += 1;
                }
                if killed_diag != killed_trace {
                    c.diag_mism += 1;
                }
                if killed_naive {
                    c.killed_naive += 1;
                }
                if killed_trace && !killed_naive {
                    c.naive_false_neg += 1;
                    if examples.len() < 40 && mu.iter().all(|&a| a as usize <= n) {
                        examples.push(format!("naive misses: N={} p={} mu={:?} minflow maxlabel={}", n, t.fmt(), mu_c, maxlab));
                    }
                }
                if !killed_trace && killed_naive {
                    c.naive_false_pos += 1;
                }
                // fully MOY composite-ladder check: sum_r c_pr closure(L_r L_s) == T[p][s]
                if let Some(cl) = &combos_l {
                    c.moy_checked += 1;
                    let mut bad = false;
                    for &sid in &order_m {
                        let mut tot = LPoly::zero();
                        for (r, a) in &cl[&p] {
                            let mut rs = mf[r].clone();
                            rs.extend(lads_m[&sid].iter().copied());
                            let v = closure_checked(&lam, &rs, n);
                            tot = tot.add(&a.to_lpoly().mul(&v));
                        }
                        if tot != tvals[&sid] {
                            bad = true;
                        }
                    }
                    if bad {
                        cnt.get_mut(&n).unwrap().moy_mism += 1;
                    }
                }
            }
        }
    }
    println!("# kill table lam={:?} K={} blocks={} types(all mu)={} time={:.1}s peakRSS={:.0}MB", lam, k, mus.len(), ntypes_total, t0.elapsed().as_secs_f64(), peak_rss_mb());
    println!("# RSK(colrev) != RS(pmax)^T for {} types", rsk_vs_rs);
    println!("# max label of minimal-flow ladder != mu_1 for {} types {:?}", maxlab_ne_mu1, maxlab_ne_examples);
    println!("N\ttypes\tkilled_trace\tkilled_rsk\tmism_rsk\tmism_RSpmaxT\tmism_RSKplain\tmism_diag\tkilled_naive\tnaive_false_neg\tnaive_false_pos\tmoy_checked\tmoy_mism\ttr_b_nonpositive\ttypes_mu<=N\tkilled_mu<=N");
    for (n, c) in &cnt {
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            n, c.types, c.killed_trace, c.killed_rsk, c.mism_rsk, c.mism_rspmax, c.mism_plain, c.diag_mism, c.killed_naive, c.naive_false_neg, c.naive_false_pos, c.moy_checked, c.moy_mism, c.tr_nonpos, c.types_nt, c.killed_nt
        );
    }
    for e in examples {
        println!("# {}", e);
    }
}

// ------------------------------------------------------------------------------------------------
// truncated scan (canonical coordinates)

#[derive(Default, Debug, Clone)]
struct TStats {
    max_unsimpl: u128,
    max_lb_obj: u128,
    max_lb_types: usize,
    max_lb_ldim: u128,
    max_peak_obj: u128,
    max_peak_ldim: u128,
    final_lb_obj: u128,
    final_types: usize,
    gram_max_w: i64,
    gram_max_full: u128,
    gram_max_full_red: u128,
    gram_max_box: u128,
    gram_max_box_red: u128,
    gram_sum_box: u128, // sum over final K0 entries of |c| * box unknowns
    killed_seen: usize,
}

enum Bas {
    Full,
    Restricted { b: HashMap<Tid, Elt>, cache: HashMap<(Tid, Vec<(usize, u8, usize)>), Vec<(Tid, Lp)>>, rank: HashMap<Tid, u128> },
}
fn rkey(rs: &[Rung]) -> Vec<(usize, u8, usize)> {
    rs.iter().map(|r| (r.i, if r.dir == Dir::E { 0u8 } else { 1u8 }, r.r)).collect()
}
fn rexpand(s: &Schur, x: &Elt, b: &HashMap<Tid, Elt>) -> Result<Vec<(Tid, Lp)>, Tid> {
    let mut acc: HashMap<Tid, Lp> = x.t.iter().cloned().collect();
    let mut heap: BinaryHeap<(i64, Tid)> = x.t.iter().map(|e| (s.lmax_of(e.0), e.0)).collect();
    let mut out = vec![];
    while let Some((_, t)) = heap.pop() {
        let m = match acc.remove(&t) {
            Some(c) => c,
            None => continue,
        };
        if m.is_zero() {
            continue;
        }
        assert!(m.is_nonneg(), "rexpand: negative multiplicity");
        let bt = b.get(&t).ok_or(t)?;
        for (q, c) in &bt.t {
            if *q == t {
                continue;
            }
            let present = acc.contains_key(q);
            let ent = acc.entry(*q).or_insert_with(Lp::zero);
            *ent = ent.sub(&m.mul(c));
            if !present {
                heap.push((s.lmax_of(*q), *q));
            }
        }
        out.push((t, m));
    }
    out.sort_unstable_by_key(|e| e.0);
    Ok(out)
}
impl Bas {
    fn tm(&mut self, s: &mut Schur, p: Tid, rs: &[Rung]) -> Vec<(Tid, Lp)> {
        match self {
            Bas::Full => s.term_mult(p, rs),
            Bas::Restricted { b, cache, .. } => {
                let key = (p, rkey(rs));
                if let Some(v) = cache.get(&key) {
                    return v.clone();
                }
                let y = s.mul_rungs(&b[&p], rs);
                let e = match rexpand(s, &y, b) {
                    Ok(e) => e,
                    Err(q) => panic!("restricted expand: b_q missing for {} (killed={})", s.table(q).fmt(), "?"),
                };
                cache.insert(key, e.clone());
                e
            }
        }
    }
    fn rank(&mut self, s: &mut Schur, labels: &[u8], p: Tid) -> u128 {
        match self {
            Bas::Full => s.rank(labels, p),
            Bas::Restricted { b, rank, .. } => {
                if let Some(&r) = rank.get(&p) {
                    return r;
                }
                let r = b[&p].t.iter().map(|(q, c)| c.at1() as u128 * s.ncosets_of(*q)).sum();
                rank.insert(p, r);
                r
            }
        }
    }
}

fn crossing_at(labels: &[u8], g: i32) -> (u8, u8, bool, usize) {
    let i = (g.unsigned_abs() - 1) as usize;
    (labels[i], labels[i + 1], g > 0, i)
}

struct Killer {
    n: Option<usize>,
    cache: HashMap<Tid, bool>,
    killed: HashSet<Tid>,
}
impl Killer {
    fn killed(&mut self, s: &Schur, p: Tid) -> bool {
        let n = match self.n {
            None => return false,
            Some(n) => n,
        };
        if let Some(&b) = self.cache.get(&p) {
            return b;
        }
        let b = mu_colrev(s.table(p))[0] > n;
        self.cache.insert(p, b);
        if b {
            self.killed.insert(p);
        }
        b
    }
}

fn tscan(s: &mut Schur, bas: &mut Bas, word: &[i32], n: Option<usize>, peak: bool) -> (TStats, Vec<(Tid, i32, i128)>) {
    let lam = s.lam.clone();
    let mut labels = lam.clone();
    let mut kill = Killer { n, cache: HashMap::new(), killed: HashSet::new() };
    let id = s.intern(Table::identity(&lam));
    let mut x: BTreeMap<i32, HashMap<Tid, Lp>> = BTreeMap::new();
    x.insert(0, [(id, Lp::one())].into_iter().collect());
    let mut st = TStats::default();
    let mut final_signed: BTreeMap<(Tid, i32), i128> = BTreeMap::new();
    let mut pcache: HashMap<(Tid, Vec<(usize, u8, usize)>), (u128, u128)> = HashMap::new();
    for step in 0..=word.len() {
        if step > 0 {
            let (c, d, pos, i) = crossing_at(&labels, word[step - 1]);
            let terms = crossing_terms(c, d, pos, i);
            let mut y: BTreeMap<i32, HashMap<Tid, Lp>> = BTreeMap::new();
            for (h, mp) in &x {
                for t in &terms {
                    let slot = y.entry(h + t.hdeg).or_default();
                    for (p, m) in mp {
                        let e = bas.tm(s, *p, &t.rungs);
                        for (q, mq) in e {
                            if kill.killed(s, q) {
                                continue;
                            }
                            let add = m.mul(&mq).shift(t.dq);
                            let ent = slot.entry(q).or_insert_with(Lp::zero);
                            *ent = ent.add(&add);
                        }
                    }
                    slot.retain(|_, v| !v.is_zero());
                }
            }
            y.retain(|_, m| !m.is_empty());
            x = y;
            labels.swap(i, i + 1);
        }
        // stats
        let mut uobj: u128 = 0;
        let mut signed: BTreeMap<(Tid, i32), i128> = BTreeMap::new();
        for (h, mp) in &x {
            let sg: i128 = if h.rem_euclid(2) == 0 { 1 } else { -1 };
            for (p, m) in mp {
                uobj += m.at1() as u128;
                for (i, &c) in m.c.iter().enumerate() {
                    if c != 0 {
                        *signed.entry((*p, m.lo + i as i32)).or_insert(0) += sg * c;
                    }
                }
            }
        }
        signed.retain(|_, c| *c != 0);
        let mut lb_obj: u128 = 0;
        let mut lb_ldim: u128 = 0;
        let mut types: HashSet<Tid> = HashSet::new();
        for (&(p, _), &c) in &signed {
            lb_obj += c.unsigned_abs();
            lb_ldim += c.unsigned_abs() * bas.rank(s, &labels, p);
            types.insert(p);
        }
        st.max_unsimpl = st.max_unsimpl.max(uobj);
        st.max_lb_obj = st.max_lb_obj.max(lb_obj);
        st.max_lb_types = st.max_lb_types.max(types.len());
        st.max_lb_ldim = st.max_lb_ldim.max(lb_ldim);
        if peak && step < word.len() {
            let (c, d, pos, i) = crossing_at(&labels, word[step]);
            let terms = crossing_terms(c, d, pos, i);
            let mut nl = labels.clone();
            nl.swap(i, i + 1);
            let (mut a, mut b) = (0u128, 0u128);
            for (&(p, _), &cc) in &signed {
                let mult = cc.unsigned_abs();
                for t in &terms {
                    let key = (p, t.rungs.iter().map(|r| (r.i, if r.dir == Dir::E { 0u8 } else { 1u8 }, r.r)).collect::<Vec<_>>());
                    let v = if let Some(v) = pcache.get(&key) {
                        *v
                    } else {
                        let e = bas.tm(s, p, &t.rungs);
                        let (mut x1, mut x2) = (0u128, 0u128);
                        for (q, m) in &e {
                            if kill.killed(s, *q) {
                                continue;
                            }
                            let kk = m.at1() as u128;
                            x1 += kk;
                            x2 += kk * bas.rank(s, &nl, *q);
                        }
                        pcache.insert(key, (x1, x2));
                        (x1, x2)
                    };
                    a += mult * v.0;
                    b += mult * v.1;
                }
            }
            st.max_peak_obj = st.max_peak_obj.max(a);
            st.max_peak_ldim = st.max_peak_ldim.max(b);
        }
        if step == word.len() {
            st.final_lb_obj = lb_obj;
            st.final_types = types.len();
            final_signed = signed;
        }
    }
    // gram
    if labels == lam {
        let lwl = lw(&lam);
        let red: Vec<u8> = lam[1..].to_vec();
        let nb = n.unwrap_or(lam.iter().map(|&a| a as usize).sum());
        let bs = box_series(&lam, nb);
        let bsr = box_series(&red, nb);
        let mut types: Vec<Tid> = final_signed.keys().map(|k| k.0).collect();
        types.dedup();
        let mut objs: HashMap<Tid, u128> = HashMap::new();
        for (&(p, _), &c) in &final_signed {
            *objs.entry(p).or_insert(0) += c.unsigned_abs();
        }
        for p in types {
            let w = s.lmax_of(p) - lwl;
            st.gram_max_w = st.gram_max_w.max(w);
            let wu = w as usize;
            let f = emono_count(&lam, wu);
            let fr = emono_count(&red, wu);
            let b = *bs.get(wu).unwrap_or(&0);
            let br = *bsr.get(wu).unwrap_or(&0);
            st.gram_max_full = st.gram_max_full.max(f);
            st.gram_max_full_red = st.gram_max_full_red.max(fr);
            st.gram_max_box = st.gram_max_box.max(b);
            st.gram_max_box_red = st.gram_max_box_red.max(br);
            st.gram_sum_box += objs[&p] * b;
        }
    }
    st.killed_seen = kill.killed.len();
    let k0 = final_signed.iter().map(|(&(p, k), &c)| (p, k, c)).collect();
    (st, k0)
}

fn chi_from_k0(s: &mut Schur, k0: &[(Tid, i32, i128)], word: &[i32], n: usize) -> LPoly {
    let lam = s.lam.clone();
    let mut labels = lam.clone();
    let mut glob: i64 = 0;
    for &g in word {
        let (c, d, pos, i) = crossing_at(&labels, g);
        let m = c.min(d) as i64;
        let eps = if pos { 1 } else { -1 };
        glob -= eps * m * (n as i64 - m);
        labels.swap(i, i + 1);
    }
    let mut memo = HashMap::new();
    let mut total = LPoly::zero();
    for &(p, k, c) in k0 {
        let t = trace_b(s, p, n, &mut memo);
        total = total.add(&t.shift(k).scale(c));
    }
    total.shift(glob as i32)
}

fn parse_word(w: &str) -> Vec<i32> {
    w.split(|c: char| c == ' ' || c == ',').filter(|x| !x.is_empty()).map(|x| x.parse().unwrap()).collect()
}
fn parse_lam(l: &str) -> Vec<u8> {
    l.split(',').map(|x| x.trim().parse().unwrap()).collect()
}

fn print_stats(tag: &str, st: &TStats) {
    println!(
        "{}\tunsimpl={}\tLBobj={}\tLBtypes={}\tLBldim={}\tpeakobj={}\tpeakldim={}\tfinalLB={}\tfinaltypes={}\tgramW={}\tgramFull={}\tgramFullRed={}\tgramBox={}\tgramBoxRed={}\tgramSumBox={}\tkilledTypesSeen={}",
        tag, st.max_unsimpl, st.max_lb_obj, st.max_lb_types, st.max_lb_ldim, st.max_peak_obj, st.max_peak_ldim, st.final_lb_obj, st.final_types, st.gram_max_w, st.gram_max_full, st.gram_max_full_red, st.gram_max_box, st.gram_max_box_red, st.gram_sum_box, st.killed_seen
    );
}

/// cmp: truncated scans at the given Ns, untruncated own scan, and krhom::combinat::schur::decat_scan (validation)
fn cmd_cmp(s: &mut Schur, name: &str, lam: Vec<u8>, word: Vec<i32>, ns: Vec<usize>, check_chi: bool, check_decat: bool) {
    let t0 = Instant::now();
    assert_eq!(s.lam, lam);
    let (su, k0u) = tscan(s, &mut Bas::Full, &word, None, true);
    print_stats(&format!("{}\tlam={:?}\tN=inf", name, lam), &su);
    if check_decat {
        let res = decat_scan(s, &word, &ScanOpts { peak: true, every_stage: true, marked_block: Some(0) });
        let mx_lb = res.stages.iter().map(|x| x.lb_objects).max().unwrap();
        let mx_peak = res.stages.iter().filter_map(|x| x.peak_before_next.map(|p| p.0)).max().unwrap_or(0);
        let mx_types = res.stages.iter().map(|x| x.lb_types).max().unwrap();
        let mx_un = res.stages.iter().map(|x| x.unsimplified_objects).max().unwrap();
        let gmax = res.gram.iter().map(|g| g.unknowns).max().unwrap_or(0);
        let ok = mx_lb == su.max_lb_obj && mx_peak == su.max_peak_obj && mx_types == su.max_lb_types && mx_un == su.max_unsimpl && gmax == su.gram_max_full;
        println!("{}\tdecat_scan\tLBobj={}\tLBtypes={}\tpeak={}\tunsimpl={}\tgramFull={}\tagree_with_own_untrunc={}", name, mx_lb, mx_types, mx_peak, mx_un, gmax, ok);
        let mut k0d = res.k0.clone();
        k0d.sort();
        let mut k0o = k0u.clone();
        k0o.sort();
        println!("{}\tK0 untrunc own == decat_scan: {}", name, k0d == k0o);
    }
    for &n in &ns {
        let (st, k0) = tscan(s, &mut Bas::Full, &word, Some(n), true);
        print_stats(&format!("{}\tlam={:?}\tN={}", name, lam, n), &st);
        if check_chi && n <= 6 {
            let c_t = chi_from_k0(s, &k0, &word, n);
            let c_u = chi_from_k0(s, &k0u, &word, n);
            println!("{}\tN={}\tchi_trunc==chi_untrunc: {}\tchi={}", name, n, c_t == c_u, fmt_lpoly(&c_t));
        }
    }
    println!("# {} time={:.1}s peakRSS={:.0}MB", name, t0.elapsed().as_secs_f64(), peak_rss_mb());
}

fn fmt_lpoly(p: &LPoly) -> String {
    if p.is_zero() {
        return "0".into();
    }
    p.c.iter().map(|(e, k)| format!("{}q^{}", k, e)).collect::<Vec<_>>().join("+")
}

// ------------------------------------------------------------------------------------------------
// restricted BFS

fn leading(s: &Schur, x: &Elt) -> (Vec<Tid>, bool) {
    let mut best = i64::MIN;
    let mut tops: Vec<Tid> = vec![];
    let mut clean = false;
    for (p, c) in &x.t {
        let l = s.lmax_of(*p);
        if l > best {
            best = l;
            tops.clear();
            tops.push(*p);
            clean = *c == Lp::one();
        } else if l == best {
            tops.push(*p);
        }
    }
    let clean = clean && tops.len() == 1;
    tops.sort_unstable();
    (tops, clean)
}

fn peel(s: &Schur, p: Tid, x: &Elt, b: &HashMap<Tid, Elt>) -> Result<Elt, Tid> {
    let mut acc: HashMap<Tid, Lp> = x.t.iter().cloned().collect();
    let mut heap: BinaryHeap<(i64, Tid)> = x.t.iter().filter(|e| e.0 != p).map(|e| (s.lmax_of(e.0), e.0)).collect();
    while let Some((_, q)) = heap.pop() {
        let m = match acc.get(&q) {
            Some(c) => c.sym_part(),
            None => continue,
        };
        if m.is_zero() {
            continue;
        }
        let bq = match b.get(&q) {
            Some(v) => v,
            None => return Err(q),
        };
        for (r, c) in &bq.t {
            let present = acc.contains_key(r);
            let ent = acc.entry(*r).or_insert_with(Lp::zero);
            *ent = ent.sub(&m.mul(c));
            if !present && *r != q {
                heap.push((s.lmax_of(*r), *r));
            }
        }
    }
    let mut t: Vec<(Tid, Lp)> = acc.into_iter().filter(|(_, c)| !c.is_zero()).collect();
    t.sort_unstable_by_key(|e| e.0);
    for (q, c) in &t {
        if *q == p {
            assert_eq!(*c, Lp::one());
        } else {
            assert!(c.min_exp().unwrap() > 0, "peel: non-positive coefficient");
        }
    }
    Ok(Elt { mu: x.mu.clone(), t })
}

/// BFS of schur.rs (clean states, dedup on leading tables) restricted to labels <= n.
/// Records the first ladder + character for every table of every block reached (or only the
/// block (lam, lam) when target_only).
fn rbfs(s: &mut Schur, n: usize, target_only: bool) -> (HashMap<Tid, (Vec<Rung>, Elt)>, usize, usize) {
    let lam = s.lam.clone();
    let nup = lam.len();
    let mut found: HashMap<Tid, (Vec<Rung>, Elt)> = HashMap::new();
    let mut seen: HashSet<Vec<Tid>> = HashSet::new();
    let lean = std::env::var("LEAN").is_ok();
    let mut queue: VecDeque<(Vec<Rung>, Elt)> = VecDeque::new();
    let start = s.unit();
    let (tops, _) = leading(s, &start);
    seen.insert(tops.clone());
    found.insert(tops[0], (vec![], start.clone()));
    queue.push_back((vec![], start));
    let mut maxq = 0usize;
    let mut popped = 0usize;
    while let Some((w, x0)) = queue.pop_front() {
        popped += 1;
        let x = if lean { s.ladder_char(&w) } else { x0 };
        let cur = x.mu.clone();
        for i in 0..nup - 1 {
            for dir in [Dir::E, Dir::F] {
                let (giver, recv) = if dir == Dir::E { (cur[i + 1], cur[i]) } else { (cur[i], cur[i + 1]) };
                for r in 1..=giver as usize {
                    if recv as usize + r > n {
                        break;
                    }
                    let rg = Rung { i, dir, r };
                    let y = s.rung(&x, &rg);
                    if y.is_zero() {
                        continue;
                    }
                    let (tops, clean) = leading(s, &y);
                    if !seen.insert(tops.clone()) {
                        continue;
                    }
                    if !clean {
                        continue;
                    }
                    let mut w2 = w.clone();
                    w2.push(rg);
                    if !target_only || y.mu == lam {
                        found.entry(tops[0]).or_insert_with(|| (w2.clone(), y.clone()));
                    }
                    queue.push_back((w2, if lean { Elt::zero(&y.mu) } else { y }));
                }
            }
        }
        maxq = maxq.max(queue.len());
    }
    (found, maxq, popped)
}

fn cmd_rbfs(lam: Vec<u8>, n: usize, target_only: bool, compare: bool) {
    let t0 = Instant::now();
    let k: usize = lam.iter().map(|&a| a as usize).sum();
    let mut s = Schur::new(&lam);
    let (found, maxq, popped) = rbfs(&mut s, n, target_only);
    let nb = n; // label bound of the BFS
    let n: usize = std::env::var("KILLN").ok().map(|x| x.parse().unwrap()).unwrap_or(n); // survival threshold
    let t_bfs = t0.elapsed().as_secs_f64();
    let mus: Vec<Vec<u8>> = comps(k, lam.len()).into_iter().filter(|m| m.iter().all(|&a| a as usize <= nb)).filter(|m| !target_only || *m == lam).collect();
    let (mut n_surv, mut n_surv_found, mut n_killed, mut n_killed_found) = (0, 0, 0, 0);
    let (mut peel_ok, mut peel_missing_killed, mut peel_missing_surv, mut peel_eq, mut peel_ne) = (0, 0, 0, 0, 0);
    let mut surv_not_found = vec![];
    let mut lam_surv = 0usize;
    let mut surv_only_fail = 0usize;
    let (mut xchar_eq, mut xchar_ne, mut xchar_fail) = (0usize, 0usize, 0usize);
    let compare_x = true;
    for mu in &mus {
        let tabs = tables(&lam, mu);
        let mut ids: Vec<Tid> = tabs.into_iter().map(|t| s.intern(t)).collect();
        ids.sort_by_key(|&p| (s.lmax_of(p), s.table(p).e.clone()));
        let mut b: HashMap<Tid, Elt> = HashMap::new();
        for &p in &ids {
            let killed = mu_colrev(s.table(p))[0] > n;
            if killed {
                n_killed += 1;
                if found.contains_key(&p) {
                    n_killed_found += 1;
                }
            } else {
                n_surv += 1;
                if *mu == lam {
                    lam_surv += 1;
                }
                if found.contains_key(&p) {
                    n_surv_found += 1;
                } else if surv_not_found.len() < 10 {
                    surv_not_found.push(s.table(p).fmt());
                }
            }
            if let Some((_, x)) = found.get(&p) {
                match peel(&s, p, x, &b) {
                    Ok(bp) => {
                        if !killed {
                            peel_ok += 1;
                        }
                        b.insert(p, bp);
                    }
                    Err(q) => {
                        if !killed {
                            if mu_colrev(s.table(q))[0] > n {
                                peel_missing_killed += 1;
                            } else {
                                peel_missing_surv += 1;
                            }
                        }
                    }
                }
            }
        }
        // surviving-only peel: does peeling a surviving type ever need a killed b_q?
        let mut b2: HashMap<Tid, Elt> = HashMap::new();
        for &p in &ids {
            if mu_colrev(s.table(p))[0] > n {
                continue;
            }
            if let Some((_, x)) = found.get(&p) {
                match peel(&s, p, x, &b2) {
                    Ok(bp) => {
                        b2.insert(p, bp);
                    }
                    Err(_) => surv_only_fail += 1,
                }
            }
        }
        // peel with ladder characters X_q = char(L_q) in place of b_q for killed q (no killed b_q needed)
        let mut b3: HashMap<Tid, Elt> = HashMap::new();
        for &p in &ids {
            if let Some((_, x)) = found.get(&p) {
                if mu_colrev(s.table(p))[0] > n {
                    b3.insert(p, x.clone());
                    continue;
                }
                match peel(&s, p, x, &b3) {
                    Ok(bp) => {
                        if compare_x {
                            if b.get(&p) == Some(&bp) { xchar_eq += 1 } else { xchar_ne += 1 }
                        }
                        b3.insert(p, bp);
                    }
                    Err(_) => xchar_fail += 1,
                }
            }
        }
        if compare {
            let full = s.canonical_basis(mu).b.clone();
            for (p, bp) in &b {
                if mu_colrev(s.table(*p))[0] > n {
                    continue;
                }
                if full[p] == *bp {
                    peel_eq += 1;
                } else {
                    peel_ne += 1;
                }
            }
        }
    }
    println!(
        "rbfs\tlam={:?}\tbound={}\tN={}\ttarget_only={}\tblocks={}\tsurviving={}\tsurv_found={}\tkilled={}\tkilled_found={}\tpeel_ok(surv)={}\tpeel_fail_missing_killed={}\tpeel_fail_missing_surv={}\tcmp_eq={}\tcmp_ne={}\tlamlam_surviving={}\tsurv_only_peel_fail={}\tstates_popped={}\tmax_queue={}\tbfs_s={:.1}\ttotal_s={:.1}\tpeakRSS={:.0}MB",
        lam, nb, n, target_only, mus.len(), n_surv, n_surv_found, n_killed, n_killed_found, peel_ok, peel_missing_killed, peel_missing_surv, peel_eq, peel_ne, lam_surv, surv_only_fail, popped, maxq, t_bfs, t0.elapsed().as_secs_f64(), peak_rss_mb()
    );
    println!("# peel with X_q = char(L_q) for killed q: equal to exact b_p {} / differ {} / fail {}", xchar_eq, xchar_ne, xchar_fail);
    if !surv_not_found.is_empty() {
        println!("# surviving not found: {:?}", surv_not_found);
    }
}

/// restricted-basis truncated scan: basis of (lam,lam) from the labels<=N BFS only
fn cmd_rscan(lam: Vec<u8>, n: usize, words: Vec<(String, Vec<i32>)>, validate: bool) {
    let t0 = Instant::now();
    let mut s = Schur::new(&lam);
    let (found, maxq, popped) = rbfs(&mut s, n, true);
    let nb = n;
    let n: usize = std::env::var("KILLN").ok().map(|x| x.parse().unwrap()).unwrap_or(n);
    let t_bfs = t0.elapsed().as_secs_f64();
    let mut ids: Vec<Tid> = found.keys().copied().collect();
    ids.sort_by_key(|&p| (s.lmax_of(p), s.table(p).e.clone()));
    let mut b: HashMap<Tid, Elt> = HashMap::new();
    let mut fail = 0;
    for &p in &ids {
        match peel(&s, p, &found[&p].1, &b) {
            Ok(bp) => {
                b.insert(p, bp);
            }
            Err(_) => fail += 1,
        }
    }
    drop(found);
    let nsurv = ids.iter().filter(|&&p| mu_colrev(s.table(p))[0] <= n).count();
    let nsurv_ok = b.keys().filter(|&&p| mu_colrev(s.table(p))[0] <= n).count();
    let all_surv = tables(&lam, &lam).iter().filter(|t| mu_colrev(t)[0] <= n).count();
    println!(
        "rscan_basis\tlam={:?}\tbound={}\tN={}\tfound={}\tpeeled={}\tpeel_fail={}\tsurviving_total={}\tsurviving_found={}\tsurviving_peeled={}\tpopped={}\tmaxq={}\tbfs_s={:.1}\tbasis_s={:.1}\tpeakRSS={:.0}MB",
        lam, nb, n, ids.len(), b.len(), fail, all_surv, nsurv, nsurv_ok, popped, maxq, t_bfs, t0.elapsed().as_secs_f64(), peak_rss_mb()
    );
    let mut bas = Bas::Restricted { b, cache: HashMap::new(), rank: HashMap::new() };
    for (name, w) in words {
        let t1 = Instant::now();
        let (st, k0) = tscan(&mut s, &mut bas, &w, Some(n), true);
        print_stats(&format!("{}\tlam={:?}\tN={}\trestricted", name, lam, n), &st);
        if validate {
            let (st2, k02) = tscan(&mut s, &mut Bas::Full, &w, Some(n), true);
            let mut a = k0.clone();
            a.sort();
            let mut c = k02.clone();
            c.sort();
            println!("{}\tvalidate restricted==full-truncated: K0 {} stats {}", name, a == c, format!("{:?}", st) == format!("{:?}", st2));
        }
        println!("# {} scan_s={:.1} peakRSS={:.0}MB", name, t1.elapsed().as_secs_f64(), peak_rss_mb());
    }
}

fn cmd_count(lam: Vec<u8>, n: usize) {
    let k: usize = lam.iter().map(|&a| a as usize).sum();
    let all = tables(&lam, &lam);
    let surv = all.iter().filter(|t| mu_colrev(t)[0] <= n).count();
    let mut tot_all = 0usize;
    let mut tot_surv = 0usize;
    let mut blocks_le_n = 0usize;
    for mu in comps(k, lam.len()) {
        if mu.iter().any(|&a| a as usize > n) {
            continue;
        }
        blocks_le_n += 1;
        let t = tables(&lam, &mu);
        tot_all += t.len();
        tot_surv += t.iter().filter(|t| mu_colrev(t)[0] <= n).count();
    }
    // weight histogram of surviving (lam, lam) types
    let lwl = lw(&lam);
    let mut maxw_all = 0;
    let mut maxw_surv = 0;
    for t in &all {
        let w = t.lmax() - lwl;
        maxw_all = maxw_all.max(w);
        if mu_colrev(t)[0] <= n {
            maxw_surv = maxw_surv.max(w);
        }
    }
    let bs = box_series(&lam, n);
    let maxbox = bs.iter().max().copied().unwrap_or(0);
    println!(
        "count\tlam={:?}\tN={}\tlamlam_types={}\tlamlam_surviving={}\tblocks_labels<=N={}\ttables_in_those={}\tsurviving_in_those={}\tmax_weight_all={}\tmax_weight_surv={}\temono_full(maxw_all)={}\tbox_dim_total={}\tbox_max_coeff={}",
        lam, n, all.len(), surv, blocks_le_n, tot_all, tot_surv, maxw_all, maxw_surv, emono_count(&lam, maxw_all as usize), bs.iter().sum::<u128>(), maxbox
    );
}

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    match a[0].as_str() {
        "kill" => cmd_kill(parse_lam(&a[1]), a[2].parse().unwrap(), a[3].parse().unwrap(), a.get(4).map(|x| x == "moy").unwrap_or(false)),
        "cmp" => {
            // cmp NAME LAM WORD NS [nochi] [nodecat]
            let ns: Vec<usize> = a[4].split(',').map(|x| x.parse().unwrap()).collect();
            let chi = !a.iter().any(|x| x == "nochi");
            let dec = !a.iter().any(|x| x == "nodecat");
            let lam = parse_lam(&a[2]);
            let mut s = Schur::new(&lam);
            cmd_cmp(&mut s, &a[1], lam, parse_word(&a[3]), ns, chi, dec)
        }
        "rbfs" => cmd_rbfs(parse_lam(&a[1]), a[2].parse().unwrap(), a.iter().any(|x| x == "targetonly"), !a.iter().any(|x| x == "nocmp")),
        "rscan" => {
            // rscan LAM N FILE [validate]   FILE lines NAME;WORD
            let txt = std::fs::read_to_string(&a[3]).unwrap();
            let words: Vec<(String, Vec<i32>)> = txt.lines().filter(|l| !l.trim().is_empty() && !l.starts_with('#')).map(|l| {
                let f: Vec<&str> = l.split(';').collect();
                (f[0].to_string(), parse_word(f[1]))
            }).collect();
            cmd_rscan(parse_lam(&a[1]), a[2].parse().unwrap(), words, a.iter().any(|x| x == "validate"))
        }
        "count" => cmd_count(parse_lam(&a[1]), a[2].parse().unwrap()),
        "mu" => {
            let t = Table::parse(&a[1]);
            println!("colrev {:?} plain {:?} rs_pmax_T {:?} pmax {:?}", mu_colrev(&t), mu_plain(&t), mu_rs_pmax_t(&t), pmax_perm(&t));
        }
        "files" => {
            // batch: each line "NAME;LAM;WORD;NS"
            let txt = std::fs::read_to_string(&a[1]).unwrap();
            let chi = !a.iter().any(|x| x == "nochi");
            let dec = !a.iter().any(|x| x == "nodecat");
            let mut schurs: HashMap<Vec<u8>, Schur> = HashMap::new();
            for line in txt.lines() {
                if line.trim().is_empty() || line.starts_with('#') {
                    continue;
                }
                let f: Vec<&str> = line.split(';').collect();
                let ns: Vec<usize> = f[3].split(',').map(|x| x.parse().unwrap()).collect();
                let lam = parse_lam(f[1]);
                let s = schurs.entry(lam.clone()).or_insert_with(|| Schur::new(&lam));
                cmd_cmp(s, f[0], lam, parse_word(f[2]), ns, chi, dec);
            }
        }
        _ => panic!("unknown command"),
    }
}
