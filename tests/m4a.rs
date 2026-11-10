//! M4a acceptance tests A1-A16 (SPEC_M4_M6a.md §3.5): the localized substrate (`coset`, `sexpr`,
//! `tmpl`, `smap`) against the exact oracles `sbim::hom_space`, `local::localize_map` and `schur`.
//!
//! Default suite: A1-A16 (A10's "passes for every template of A4-A9" and A13's round trip run
//! inside the tests that build those templates). Ignored (run once each, timings reported):
//! A7(iv) a + b <= 9, A9 (3,3), A11 on (2,2), the release performance targets of §3.6 and the
//! (2,2,2,2) FlowIndex.
use krhom::coset::*;
use krhom::field::{Field, F31};
use krhom::local::{fixed_points, localize_map};
use krhom::mpoly::Poly;
use krhom::pmat::PMat;
use krhom::sbim::SBim;
use krhom::schur::{crossing_terms, Elt, Schur, Tid};
use krhom::sexpr::*;
use krhom::smap::*;
use krhom::tmpl::*;
use krhom::web::{Dir, Rung};
use rand::{rngs::StdRng, Rng, SeedableRng};
use std::collections::{BTreeMap, HashMap};
use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;

type F = F31;

// ------------------------------------------------------------------ helpers

fn ex(b: &[u8], r: &[Rung]) -> Expr {
    Expr::new(b, r).unwrap_or_else(|e| panic!("{}", e))
}
fn win(e: &Expr) -> Window {
    Window::new(&e.bottom, &e.rungs)
}
fn e(i: usize, r: usize) -> Rung {
    Rung::e(i, r)
}
fn f(i: usize, r: usize) -> Rung {
    Rung::f(i, r)
}
fn shifted(rs: &[Rung], by: usize) -> Vec<Rung> {
    rs.iter().map(|r| Rung { i: r.i + by, dir: r.dir, r: r.r }).collect()
}
fn apply(l: &[u8], r: &Rung) -> Option<Vec<u8>> {
    let (rec, giv) = rung_ends(r);
    if giv >= l.len() || rec >= l.len() || (l[giv] as usize) < r.r || r.r == 0 {
        return None;
    }
    let mut v = l.to_vec();
    v[giv] -= r.r as u8;
    v[rec] += r.r as u8;
    Some(v)
}
fn cache_dir(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("m4a_cache").join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// localize_map of an exact map, reordered into sorted (FlowIndex) order: rows tgt, cols src
fn oracle_sorted(src: &Expr, tgt: &Expr, phi: &PMat<F>, x: &[F]) -> Vec<Vec<F>> {
    let (sb, tb) = (src.sbim::<F>(), tgt.sbim::<F>());
    let (p, pp, fm) = localize_map(&sb, &tb, phi, x);
    let (sf, tf) = (FlowIndex::build(src), FlowIndex::build(tgt));
    let sp: Vec<usize> = p.iter().map(|q| sf.locate(&q.moved).expect("fixed point not a flow") as usize).collect();
    let tp: Vec<usize> = pp.iter().map(|q| tf.locate(&q.moved).expect("fixed point not a flow") as usize).collect();
    let mut out = vec![vec![F::zero(); sp.len()]; tp.len()];
    for (q, row) in fm.iter().enumerate() {
        for (pi, v) in row.iter().enumerate() {
            out[tp[q]][sp[pi]] = *v;
        }
    }
    out
}
fn eval_dense(ar: &Arena<F>, m: MapId, pt: &Point<F>) -> Vec<Vec<F>> {
    Eval::new(ar, pt.clone()).winloc(m).dense()
}
fn nnz(m: &[Vec<F>]) -> usize {
    m.iter().map(|r| r.iter().filter(|x| !x.is_zero()).count()).sum()
}
/// c with a = c·b entrywise (zeros must coincide), the same at every point
fn ratio(ar: &Arena<F>, a: MapId, b: MapId, pts: &[Point<F>]) -> Result<F, String> {
    let mut c: Option<F> = None;
    for p in pts {
        let (da, db) = (eval_dense(ar, a, p), eval_dense(ar, b, p));
        for (ra, rb) in da.iter().zip(&db) {
            for (x, y) in ra.iter().zip(rb) {
                if y.is_zero() {
                    if !x.is_zero() {
                        return Err(format!("closed form nonzero where the exact map vanishes (point {})", p.id));
                    }
                    continue;
                }
                let r = x.mul(&y.inv());
                match c {
                    None => c = Some(r),
                    Some(c0) if c0 == r => {}
                    Some(c0) => return Err(format!("ratio not constant: {} vs {} (point {})", c0, r, p.id)),
                }
            }
        }
    }
    c.ok_or_else(|| "exact map is zero".to_string())
}

/// requests made of a template library (to replay them against a fresh library: A13)
#[derive(Clone, Debug)]
enum Req {
    Hom(Window, Window, i32),
    Rick(u8, u8, bool),
}
fn replay(lib: &mut TmplLib<F>, r: &Req) -> Vec<TmplId> {
    match r {
        Req::Hom(s, t, d) => lib.hom_basis(s, t, *d),
        Req::Rick(c, d, p) => lib.rickard(*c, *d, *p).1,
    }
}
/// A13: save, then a fresh library on the same directory loads identical PMats without computing
fn roundtrip(lib: &TmplLib<F>, reqs: &[Req], name: &str) {
    lib.save().unwrap();
    let dir = lib.cache_dir().expect("library without cache dir").to_path_buf();
    let mut fresh: TmplLib<F> = TmplLib::new(Some(dir));
    let mut fresh_ids = vec![];
    for r in reqs {
        fresh_ids.push(replay(&mut fresh, r));
    }
    let st = fresh.stats();
    assert_eq!(st.computed, 0, "{}: fresh library recomputed {} calls", name, st.computed);
    assert_eq!(st.rejected, 0, "{}: fresh library rejected cache files", name);
    // compare with the original library's templates of the same request
    let mut n = 0;
    for (r, ids2) in reqs.iter().zip(&fresh_ids) {
        let ids1: Vec<TmplId> = (0..lib.len() as TmplId)
            .filter(|&t| {
                let tm = lib.get(t);
                match (r, &tm.source) {
                    (Req::Hom(s, tg, d), TmplSource::HomSpace { .. }) => tm.src == *s && tm.tgt == *tg && tm.deg == *d,
                    (Req::Rick(c, dd, p), TmplSource::Rickard { c: c2, d: d2, positive, .. }) => c == c2 && dd == d2 && p == positive,
                    _ => false,
                }
            })
            .collect();
        assert_eq!(ids1.len(), ids2.len(), "{}: {:?}: basis sizes differ after reload", name, r);
        for (&a, &b) in ids1.iter().zip(ids2) {
            let (pa, pb) = (&lib.get(a).phi, &fresh.get(b).phi);
            assert!(pa.sub(pb).is_zero(), "{}: {:?}: reloaded PMat differs", name, r);
            n += 1;
        }
    }
    eprintln!("[A13] {}: {} calls, {} templates reloaded identically ({} files loaded)", name, reqs.len(), n, st.loaded);
}

/// A10: certificate of a template (must pass)
fn cert_template(lib: &TmplLib<F>, t: TmplId, seed: u64) -> CertReport {
    let tm = lib.get(t);
    let (s, tg, d) = (tm.src.clone(), tm.tgt.clone(), tm.deg);
    poly_certificate::<F>(&s, &tg, d, &mut |y: &[F]| lib.localize(t, y), seed).unwrap_or_else(|e| panic!("certificate of template {:?} -> {:?} (deg {}) failed: {}", s, tg, d, e))
}

/// random rung word on `bottom` (as web::random_closed_ladder draws rungs, without closing)
fn random_rungs(rng: &mut StdRng, bottom: &[u8], nrungs: usize) -> Vec<Rung> {
    let n = bottom.len();
    let mut lab = bottom.to_vec();
    let mut out = vec![];
    for _ in 0..nrungs {
        for _try in 0..20 {
            let i = rng.gen_range(0..n - 1);
            let dir = if rng.gen_bool(0.5) { Dir::E } else { Dir::F };
            let giv = if dir == Dir::E { i + 1 } else { i };
            let g = lab[giv] as usize;
            if g == 0 {
                continue;
            }
            let r = rng.gen_range(1..=g);
            let rg = Rung { i, dir, r };
            lab = apply(&lab, &rg).unwrap();
            out.push(rg);
            break;
        }
    }
    out
}

/// the 280 random ladders of A2 (seeds 1..=280); flows capped at 20 000 (redrawn)
fn a2_ladders() -> Vec<Expr> {
    let mut out = vec![];
    for s in 1..=280u64 {
        let mut rng = StdRng::seed_from_u64(s);
        let n = rng.gen_range(2..=4usize);
        let lab: Vec<u8> = loop {
            let l: Vec<u8> = (0..n).map(|_| rng.gen_range(0..=3u8)).collect();
            let k: u32 = l.iter().map(|&x| x as u32).sum();
            if (1..=8).contains(&k) {
                break l;
            }
        };
        let ex = loop {
            let nr = rng.gen_range(1..=6usize);
            let rs = random_rungs(&mut rng, &lab, nr);
            if rs.is_empty() {
                continue;
            }
            let x = ex(&lab, &rs);
            if x.nflows() <= 20_000 {
                break x;
            }
        };
        out.push(ex);
    }
    out
}

/// (block, rung) pairs of crossing terms on a uniformly colored lambda (probe P1's work list)
fn work_pairs(lam: &[u8]) -> Vec<(Vec<u8>, Rung)> {
    let n = lam.len();
    let mut work: Vec<(Vec<u8>, Rung)> = vec![];
    for i in 0..n - 1 {
        for pos in [true, false] {
            for t in crossing_terms(lam[i], lam[i + 1], pos, i) {
                let mut l = lam.to_vec();
                for r in &t.rungs {
                    if !work.contains(&(l.clone(), *r)) {
                        work.push((l.clone(), *r));
                    }
                    l = apply(&l, r).unwrap();
                }
            }
        }
    }
    work
}

/// canonical ladders of all types of the blocks reached, and X = L_p R of every work pair
fn canonical_and_x(lam: &[u8]) -> (Vec<Expr>, Vec<Expr>) {
    let mut s = Schur::new(lam);
    let work = work_pairs(lam);
    let mut blocks: Vec<Vec<u8>> = vec![];
    for (mu, _) in &work {
        if !blocks.contains(mu) {
            blocks.push(mu.clone());
        }
    }
    let mut can = vec![];
    for mu in &blocks {
        let order = s.canonical_basis(mu).order.clone();
        for p in order {
            can.push(ex(lam, &s.canonical_basis(mu).ladder_minflow[&p].clone()));
        }
    }
    let mut xs = vec![];
    for (mu, r) in &work {
        let order = s.canonical_basis(mu).order.clone();
        for p in order {
            let mut l = s.canonical_basis(mu).ladder_minflow[&p].clone();
            l.push(*r);
            xs.push(ex(lam, &l));
        }
    }
    (can, xs)
}

/// canonical ladders of block (lam, lam) with an extra rung appended
fn block_ladders(lam: &[u8], extra: Option<Rung>) -> Vec<Expr> {
    let mut s = Schur::new(lam);
    let order = s.canonical_basis(lam).order.clone();
    order
        .iter()
        .map(|p| {
            let mut l = s.canonical_basis(lam).ladder_minflow[p].clone();
            if let Some(r) = extra {
                l.push(r);
            }
            ex(lam, &l)
        })
        .collect()
}

/// a dummy SBim carrying only bottom labels and rungs (enough for local::fixed_points)
fn fp_sbim(x: &Expr) -> SBim<F> {
    let lab: Vec<usize> = x.bottom.iter().map(|&a| a as usize).collect();
    let (ring, voff) = krhom::sbim::bottom_ring(&lab);
    SBim { n: 2 * x.k(), ring, voff, bot_labels: lab.clone(), top_labels: lab, fdeg: vec![], top: vec![], shift: 0, rungs: x.rungs.clone() }
}

fn partitions_in_box(rows: usize, cols: usize) -> Vec<Vec<u8>> {
    // partitions with at most `rows` parts, each <= cols
    let mut out = vec![];
    fn rec(rows: usize, maxp: usize, cur: &mut Vec<u8>, out: &mut Vec<Vec<u8>>) {
        out.push(cur.clone());
        if cur.len() == rows {
            return;
        }
        for p in 1..=maxp {
            cur.push(p as u8);
            rec(rows, p, cur, out);
            cur.pop();
        }
    }
    rec(rows, cols, &mut vec![], &mut out);
    out
}

/// permutations of 0..n with signs
fn perms(n: usize) -> Vec<(Vec<usize>, bool)> {
    if n == 0 {
        return vec![(vec![], false)];
    }
    let mut out = vec![];
    for (p, odd) in perms(n - 1) {
        for pos in 0..=p.len() {
            let mut q = p.clone();
            q.insert(pos, n - 1);
            // inserting n-1 at pos creates (len - pos) inversions
            out.push((q, odd ^ ((p.len() - pos) % 2 == 1)));
        }
    }
    out
}
fn conj(nu: &[u8]) -> Vec<usize> {
    let l = nu.first().copied().unwrap_or(0) as usize;
    (1..=l).map(|i| nu.iter().filter(|&&x| x as usize >= i).count()).collect()
}

// ------------------------------------------------------------------ A1

#[test]
fn a1_cosets() {
    let rows: Vec<(&[u8], &[u8], u64)> = vec![
        (&[2, 2, 2], &[2, 2, 2], 90),
        (&[2, 2, 2], &[3, 1, 2], 60),
        (&[2, 2, 2], &[4, 0, 2], 15),
        (&[2, 2, 2], &[2, 3, 1], 60),
        (&[2, 2, 2], &[2, 4, 0], 15),
        (&[3, 3], &[3, 3], 20),
        (&[3, 3], &[4, 2], 15),
        (&[3, 3], &[5, 1], 6),
        (&[3, 3], &[6, 0], 1),
        (&[2, 2], &[2, 2], 6),
        (&[2, 1, 2], &[2, 1, 2], 30),
        (&[1, 1, 1], &[1, 1, 1], 6),
        (&[3, 3, 3], &[3, 3, 3], 1680),
        (&[2, 2, 2, 2], &[2, 2, 2, 2], 2520),
    ];
    let mut rng = StdRng::seed_from_u64(11);
    for (lam, mu, want) in rows {
        let cs = CosetSpace::new(lam, mu);
        assert_eq!(cs.count(), want, "count {:?} -> {:?}", lam, mu);
        if want <= 1680 {
            for r in 0..want {
                let c = cs.unrank(r);
                assert_eq!(cs.rank(&c), r, "rank(unrank) {:?} -> {:?}", lam, mu);
                let mut union = 0u64;
                for (j, m) in c.0.iter().enumerate() {
                    assert_eq!(m.count_ones() as u8, mu[j]);
                    assert_eq!(union & m, 0, "unrank: overlapping masks");
                    union |= m;
                }
                assert_eq!(union, low_mask(cs.roots.k), "unrank: masks do not cover the roots");
            }
        }
        let tabs = krhom::schur::tables(lam, mu);
        let mut total = 0u64;
        let mut all = vec![];
        for t in &tabs {
            let cos = cs.cosets_of(t);
            assert_eq!(cos.len() as u128, t.ncosets(), "cosets_of({}) for {:?} -> {:?}", t.fmt(), lam, mu);
            assert!(cos.windows(2).all(|w| w[0] < w[1]));
            for &c in &cos {
                assert_eq!(&cs.table(&cs.unrank(c)), t, "table of a coset of {}", t.fmt());
            }
            total += cos.len() as u64;
            all.extend(cos);
            let w = cs.rep(t);
            assert_eq!(&cs.table(&w), t, "table(rep(t)) for {}", t.fmt());
        }
        assert_eq!(total, want, "sum of N(t) over tables of {:?} -> {:?}", lam, mu);
        all.sort();
        all.dedup();
        assert_eq!(all.len() as u64, want);
        if lam == [2, 2, 2] && mu == [2, 2, 2] {
            assert_eq!(tabs.len(), 21);
        }
        // act preserves tables: random block-preserving permutations
        let roots = Roots::new(lam);
        for _ in 0..20 {
            let mut sigma: Vec<usize> = (0..roots.k).collect();
            for i in 0..lam.len() {
                let o = roots.off[i];
                let l = lam[i] as usize;
                for a in (1..l).rev() {
                    let b = rng.gen_range(0..=a);
                    sigma.swap(o + a, o + b);
                }
            }
            let c = cs.unrank(rng.gen_range(0..want));
            let d = cs.act(&sigma, &c);
            assert_eq!(cs.table(&c), cs.table(&d), "act changed the table");
        }
    }
}

// ------------------------------------------------------------------ A2, A3, A16 (FlowIndex bytes)

#[test]
fn a2_a3_flows_deodhar_and_fixed_points() {
    let lads = a2_ladders();
    assert_eq!(lads.len(), 280);
    let mut schurs: HashMap<Vec<u8>, Schur> = HashMap::new();
    let mut maxf = 0;
    for x in &lads {
        let fi = FlowIndex::build(x);
        let s = schurs.entry(x.bottom.clone()).or_insert_with(|| Schur::new(&x.bottom));
        deodhar_check(&fi, s).unwrap_or_else(|e| panic!("A2: {}", e));
        // A3: multiset of moved masks = local::fixed_points; total = schur::flows
        let mut ours: Vec<Vec<u64>> = (0..fi.nflows() as u32).map(|f| fi.flow(f).to_vec()).collect();
        let mut theirs: Vec<Vec<u64>> = fixed_points(&fp_sbim(x)).into_iter().map(|p| p.moved).collect();
        ours.sort();
        theirs.sort();
        assert_eq!(ours, theirs, "A3: flows of {:?} {:?} differ from local::fixed_points", x.bottom, x.rungs);
        assert_eq!(fi.nflows() as u128, x.nflows());
        // A16: packed index size
        let bound = 8 * fi.nr * fi.nflows() + 12 * fi.support().len() + 1024;
        assert!(fi.bytes() <= bound, "A16: FlowIndex::bytes {} > {} for {:?} {:?}", fi.bytes(), bound, x.bottom, x.rungs);
        // find / locate consistency
        for (i, &c) in fi.support().iter().enumerate() {
            for g in fi.start[i]..fi.start[i + 1] {
                assert_eq!(fi.find(c, fi.flow(g)), Some(g - fi.start[i]));
                assert_eq!(fi.locate(fi.flow(g)), Some(g));
            }
        }
        maxf = maxf.max(fi.nflows());
    }
    eprintln!("[A2] 280 random ladders pass deodhar_check and match fixed_points (max flows {})", maxf);
}

#[test]
fn a2_scan_blocks_and_maxima() {
    for (lam, ncan, nx) in [(vec![2u8, 2, 2], 63usize, 126usize), (vec![3, 3], 10, 18)] {
        let (can, xs) = canonical_and_x(&lam);
        assert_eq!((can.len(), xs.len()), (ncan, nx), "number of canonical ladders and X for {:?}", lam);
        let mut s = Schur::new(&lam);
        let mut mb = (0, 0);
        for (k, x) in can.iter().chain(&xs).enumerate() {
            let fi = FlowIndex::build(x);
            deodhar_check(&fi, &mut s).unwrap_or_else(|e| panic!("A2 {:?}: {}", lam, e));
            if k < can.len() {
                mb.0 = mb.0.max(fi.max_block());
            } else {
                mb.1 = mb.1.max(fi.max_block());
            }
        }
        eprintln!("[A2] {:?}: {} canonical (max block {}), {} X (max block {}) pass", lam, can.len(), mb.0, xs.len(), mb.1);
    }
    // recorded maxima [V-P5]
    let mb = |v: &[Expr]| v.iter().map(|x| FlowIndex::build(x).max_block()).max().unwrap();
    let c222 = block_ladders(&[2, 2, 2], None);
    let x222 = block_ladders(&[2, 2, 2], Some(e(0, 1)));
    let x33 = block_ladders(&[3, 3], Some(e(0, 2)));
    assert_eq!(c222.len(), 21);
    assert_eq!(mb(&c222), 8, "A2: max_block of the (2,2,2) canonical ladders");
    assert_eq!(mb(&x222), 20, "A2: max_block of (2,2,2) X = L_p E_0^(1)");
    assert_eq!(mb(&x33), 18, "A2: max_block of (3,3) X = L_p E_0^(2)");
    // A16 / §3.6: the (2,2,2) X with 432 flows and 7 rungs
    let big = x222.iter().max_by_key(|x| x.nflows()).unwrap();
    assert_eq!((big.nflows(), big.rungs.len()), (432, 7));
    let t0 = Instant::now();
    let fi = FlowIndex::build(big);
    let dt = t0.elapsed().as_secs_f64();
    assert!(fi.bytes() <= 26_000, "A16: FlowIndex of the 432-flow X takes {} B > 26 KB", fi.bytes());
    eprintln!("[A16] (2,2,2) X 432 flows x 7 rungs: FlowIndex {} B, built in {:.2} ms (this profile)", fi.bytes(), dt * 1e3);
}

// ------------------------------------------------------------------ A4

#[test]
fn a4_whole_ladder_bridge() {
    let dir = cache_dir("a4");
    let mut lib: TmplLib<F> = TmplLib::new(Some(dir));
    let mut reqs = vec![];
    let mut pairs: Vec<(Expr, Expr)> = vec![];
    let mut seed = 4000u64;
    while pairs.len() < 20 {
        seed += 1;
        let mut rng = StdRng::seed_from_u64(seed);
        let n = rng.gen_range(2..=3usize);
        let lab: Vec<u8> = (0..n).map(|_| rng.gen_range(0..=3u8)).collect();
        let k: u32 = lab.iter().map(|&x| x as u32).sum();
        if !(2..=6).contains(&k) {
            continue;
        }
        let mut by_top: BTreeMap<Vec<u8>, Vec<Expr>> = BTreeMap::new();
        for _ in 0..12 {
            let nr = rng.gen_range(1..=4usize);
            let rs = random_rungs(&mut rng, &lab, nr);
            if rs.is_empty() {
                continue;
            }
            let x = ex(&lab, &rs);
            if x.nflows() <= 40 {
                let v = by_top.entry(x.top()).or_default();
                if !v.contains(&x) {
                    v.push(x);
                }
            }
        }
        if let Some(v) = by_top.values().find(|v| v.len() >= 2) {
            pairs.push((v[0].clone(), v[1].clone()));
        }
    }
    let mut ntm = 0;
    let mut certs = vec![];
    for (a, b) in &pairs {
        let mut s = Schur::new(&a.bottom);
        let (ca, cb) = (s.ladder_char(&a.rungs), s.ladder_char(&b.rungs));
        let gd = s.hom_gd(&ca, &cb, -40, 40);
        let i = gd.iter().position(|&x| x > 0).expect("Hom = 0 in -40..40");
        let d = -40 + i as i32;
        let ts = lib.hom_basis(&win(a), &win(b), d);
        reqs.push(Req::Hom(win(a), win(b), d));
        assert_eq!(ts.len() as i128, gd[i], "A4: hom_basis dimension != hom_gd for {:?} -> {:?} degree {}", a.rungs, b.rungs, d);
        let mut ar = Arena::new(&a.bottom, std::mem::replace(&mut lib, TmplLib::new(None)));
        let (sa, sb) = (ar.expr(a.clone()), ar.expr(b.clone()));
        let pts = points::<F>(a.k(), 40 + ntm as u64, 2, 0);
        for &t in &ts {
            let m = ar.add(MapNode::Local { src: sa, tgt: sb, at: 0, first: 0, tmpl: t }, "A4").unwrap();
            assert_eq!(ar.info(m).deg, d);
            assert!(ar.info(m).certified);
            for p in &pts {
                let want = oracle_sorted(a, b, &ar.lib_ref().get(t).phi, &p.x);
                assert_eq!(eval_dense(&ar, m, p), want, "A4: {:?} -> {:?} degree {} at point {}", a.rungs, b.rungs, d, p.id);
            }
            certs.push(cert_template(ar.lib_ref(), t, 7 + t as u64));
            ntm += 1;
        }
        lib = ar.into_lib();
    }
    roundtrip(&lib, &reqs, "A4");
    let maxm = certs.iter().map(|c| c.max_monos).max().unwrap();
    eprintln!("[A4] 20 pairs, {} templates match localize_map at 2 points; certificates pass (max |Mono_D| {})", ntm, maxm);
}

// ------------------------------------------------------------------ A5 (+ A16 WinLoc::dense)

struct WinCase {
    name: &'static str,
    labels: Vec<u8>,
    a: Vec<Rung>,
    a2: Vec<Rung>,
    rick: Option<(u8, u8)>,
    deg: i32,
    middle: bool,
}

fn a5_cases() -> Vec<WinCase> {
    vec![
        WinCase { name: "E2 (1,1) [E(1),F(1)] -> [] deg 1", labels: vec![1, 1], a: vec![e(0, 1), f(0, 1)], a2: vec![], rick: None, deg: 1, middle: true },
        WinCase { name: "E2 (1,1) [] -> [E(1),F(1)] deg 1", labels: vec![1, 1], a: vec![], a2: vec![e(0, 1), f(0, 1)], rick: None, deg: 1, middle: false },
        WinCase { name: "dots (2,2) [E(1)] -> [E(1)] deg 2", labels: vec![2, 2], a: vec![e(0, 1)], a2: vec![e(0, 1)], rick: None, deg: 2, middle: true },
        WinCase { name: "Rickard (1,1) d0 top", labels: vec![1, 1], a: vec![], a2: vec![], rick: Some((1, 1)), deg: 1, middle: false },
        WinCase { name: "Rickard (1,1) d0 middle", labels: vec![1, 1], a: vec![], a2: vec![], rick: Some((1, 1)), deg: 1, middle: true },
        WinCase { name: "Rickard (2,1) d0 top", labels: vec![2, 1], a: vec![], a2: vec![], rick: Some((2, 1)), deg: 1, middle: false },
        WinCase { name: "Rickard (2,1) d0 middle", labels: vec![2, 1], a: vec![], a2: vec![], rick: Some((2, 1)), deg: 1, middle: true },
        WinCase { name: "3-upright [E0(1),E1(1)] -> itself deg 2", labels: vec![1, 1, 1], a: vec![e(0, 1), e(1, 1)], a2: vec![e(0, 1), e(1, 1)], rick: None, deg: 2, middle: true },
    ]
}

#[test]
fn a5_local_application_and_kron() {
    let dir = cache_dir("a5");
    let mut lib: TmplLib<F> = TmplLib::new(Some(dir));
    let mut reqs = vec![];
    let mut checked = 0;
    for (ci, cs) in a5_cases().into_iter().enumerate() {
        let (ts, wa, wa2) = match cs.rick {
            Some((c, d)) => {
                let (wins, ids) = lib.rickard(c, d, true);
                reqs.push(Req::Rick(c, d, true));
                (vec![ids[0]], wins[0].clone(), wins[1].clone())
            }
            None => {
                let (w1, w2) = (Window::new(&cs.labels, &cs.a), Window::new(&cs.labels, &cs.a2));
                reqs.push(Req::Hom(w1.clone(), w2.clone(), cs.deg));
                (lib.hom_basis(&w1, &w2, cs.deg), w1, w2)
            }
        };
        assert!(!ts.is_empty(), "A5 {}: empty Hom", cs.name);
        let w = cs.labels.len();
        // random prefix P (1-3 rungs on 3 uprights, K <= 6) reaching the window labels
        let mut rng = StdRng::seed_from_u64(500 + ci as u64);
        let (bottom, pre, first) = loop {
            let lab: Vec<u8> = (0..3).map(|_| rng.gen_range(0..=3u8)).collect();
            let k: u32 = lab.iter().map(|&x| x as u32).sum();
            if k > 6 {
                continue;
            }
            let np = rng.gen_range(1..=3usize);
            let pr = random_rungs(&mut rng, &lab, np);
            if pr.is_empty() {
                continue;
            }
            let first = if w == 3 { 0 } else { rng.gen_range(0..=1usize) };
            let lv = ex(&lab, &pr).top();
            if lv[first..first + w] == cs.labels[..] {
                break (lab, pr, first);
            }
        };
        let pa = ex(&bottom, &[pre.clone(), shifted(&wa.rungs, first)].concat());
        let pa2 = ex(&bottom, &[pre.clone(), shifted(&wa2.rungs, first)].concat());
        let q = loop {
            let nq = if cs.middle { rng.gen_range(1..=2usize) } else { 0 };
            let q = random_rungs(&mut rng, &pa.top(), nq);
            if q.len() == nq {
                break q;
            }
        };
        let src = pa.concat(&q).unwrap();
        let tgt = pa2.concat(&q).unwrap();
        let fq = ex(&pa.top(), &q).nflows() as usize;
        let psb = ex(&bottom, &pre).sbim::<F>();
        let mut ar = Arena::new(&bottom, std::mem::replace(&mut lib, TmplLib::new(None)));
        let (si, ti) = (ar.expr(src.clone()), ar.expr(tgt.clone()));
        let (pai, pa2i) = (ar.expr(pa.clone()), ar.expr(pa2.clone()));
        let pts = points::<F>(src.k(), 900 + ci as u64, 2, 0);
        for &t in &ts {
            let m = ar.add(MapNode::Local { src: si, tgt: ti, at: pre.len(), first, tmpl: t }, "A5").unwrap();
            let fm = ar.add(MapNode::Local { src: pai, tgt: pa2i, at: pre.len(), first, tmpl: t }, "A5 f").unwrap();
            let km = ar.add(MapNode::Kron(fm, q.clone()), "A5 kron").unwrap();
            assert_eq!((ar.info(km).src, ar.info(km).tgt), (si, ti), "Kron ends");
            let tm = ar.lib_ref().get(t);
            let phi = tensor_exact(&psb, first, tm, fq);
            let phik = tensor_exact(&psb, first, tm, 1).kron_id(fq);
            for p in &pts {
                let want = oracle_sorted(&src, &tgt, &phi, &p.x);
                assert_eq!(eval_dense(&ar, m, p), want, "A5 {}: Local != localize_map(tensor_exact) (P {:?}, first {}, Q {:?})", cs.name, pre, first, q);
                let wantk = oracle_sorted(&src, &tgt, &phik, &p.x);
                assert_eq!(eval_dense(&ar, km, p), wantk, "A5 {}: Kron != localize_map(f.kron_id)", cs.name);
                // A16: WinLoc::dense = reordered localize_map of the window
                let y: Vec<F> = p.x.iter().take(tm.src.labels.iter().map(|&a| a as usize).sum()).cloned().collect();
                assert_eq!(ar.lib_ref().localize(t, &y).dense(), oracle_sorted(&tm.src.expr(), &tm.tgt.expr(), &tm.phi, &y), "A16: WinLoc::dense of {}", cs.name);
            }
            cert_template(ar.lib_ref(), t, 50 + t as u64);
            checked += 1;
        }
        lib = ar.into_lib();
    }
    roundtrip(&lib, &reqs, "A5");
    eprintln!("[A5] {} window templates: Local and Kron match localize_map(tensor_exact) at 2 points", checked);
}

// ------------------------------------------------------------------ A6

/// s_nu as a polynomial in the variables e_1..e_l (indices off..off+l) of one upright
fn schur_poly(nu: &[u8], off: usize, l: usize) -> Poly<F> {
    let cj = conj(nu);
    let n = cj.len();
    let ek = |k: isize| -> Poly<F> {
        if k == 0 {
            Poly::one()
        } else if k < 0 || k as usize > l {
            Poly::zero()
        } else {
            Poly::var(off + k as usize - 1)
        }
    };
    let mut s = Poly::zero();
    for (p, odd) in perms(n) {
        let mut t = Poly::one();
        for i in 0..n {
            t = t.mul(&ek(cj[i] as isize - i as isize + p[i] as isize));
        }
        s = if odd { s.sub(&t) } else { s.add(&t) };
    }
    s
}
/// s_nu of the top operators of upright u of an SBim (Jacobi-Trudi with PMat products)
fn schur_top(b: &SBim<F>, u: usize, nu: &[u8]) -> PMat<F> {
    let cj = conj(nu);
    let n = cj.len();
    let fsz = b.fsize();
    let l = b.top_labels[u];
    let ek = |k: isize| -> PMat<F> {
        if k == 0 {
            PMat::identity(fsz)
        } else if k < 0 || k as usize > l {
            PMat::zero(fsz, fsz)
        } else {
            b.top[u][k as usize - 1].clone()
        }
    };
    let mut s = PMat::zero(fsz, fsz);
    for (p, odd) in perms(n) {
        let mut t = PMat::identity(fsz);
        for i in 0..n {
            t = t.matmul(&ek(cj[i] as isize - i as isize + p[i] as isize));
        }
        s = if odd { s.sub(&t) } else { s.add(&t) };
    }
    s
}

#[test]
fn a6_dots() {
    let x = ex(&[2, 1], &[e(0, 1), f(0, 2)]); // levels (2,1) -> (3,0) -> (1,2)
    let xb = x.sbim::<F>();
    let mut ar = Arena::<F>::new(&[2, 1], TmplLib::new(None));
    let xi = ar.expr(x.clone());
    let nus: Vec<Vec<u8>> = vec![vec![1], vec![1, 1], vec![2], vec![2, 1], vec![3], vec![1, 1, 1]];
    let pts = points::<F>(3, 66, 2, 0);
    let mut n = 0;
    // top level: e_j = SBim::top[u][j-1] and Jacobi-Trudi products, both uprights
    for u in 0..2 {
        for nu in &nus {
            let m = ar.add(MapNode::Dot { e: xi, level: 2, upright: u, nu: nu.clone() }, "A6").unwrap();
            let exact = schur_top(&xb, u, nu);
            if nu.iter().all(|&a| a == 1) && nu.len() <= xb.top_labels[u] {
                assert!(exact.sub(&xb.top[u][nu.len() - 1]).is_zero());
            }
            for p in &pts {
                assert_eq!(eval_dense(&ar, m, p), oracle_sorted(&x, &x, &exact, &p.x), "A6: top dot {:?} on upright {}", nu, u);
            }
            n += 1;
        }
    }
    // middle (1) and bottom (0) levels: Local of the window identity times a dot (tensor_exact)
    for level in [0usize, 1] {
        let lv = x.level(level);
        let pre = ex(&[2, 1], &x.rungs[..level]).sbim::<F>();
        let fq = ex(&lv, &x.rungs[level..]).nflows() as usize;
        for u in 0..2 {
            if lv[u] == 0 {
                continue;
            }
            for nu in &nus {
                let m = ar.add(MapNode::Dot { e: xi, level, upright: u, nu: nu.clone() }, "A6").unwrap();
                let wdw = Window::new(&[lv[u]], &[]);
                let t = Template::new(wdw.clone(), wdw.clone(), 2 * nu.iter().map(|&a| a as i32).sum::<i32>(), PMat::scalar(1, &schur_poly(nu, 0, lv[u] as usize)), TmplSource::HomSpace { index: 0 });
                let exact = tensor_exact(&pre, u, &t, fq);
                for p in &pts {
                    assert_eq!(eval_dense(&ar, m, p), oracle_sorted(&x, &x, &exact, &p.x), "A6: dot {:?} at level {} upright {}", nu, level, u);
                }
                n += 1;
            }
        }
    }
    eprintln!("[A6] {} dots match the exact top operators / tensor_exact", n);
}

// ------------------------------------------------------------------ A7

/// E2 and E3 constants, products and digon checks for a + b <= amax (E2) and c <= amax (E3)
fn e2e3_checks(amax: usize, tag: &str, do_e2: bool, do_e3: bool) {
    let mut consts: Vec<String> = vec![];
    let mut certs = 0;
    let mut reqs_all = 0;
    let t0 = Instant::now();
    let dir = cache_dir(&format!("a7_{}", tag));
    let mut lib: TmplLib<F> = TmplLib::new(Some(dir));
    let mut reqs = vec![];
    // E2 over (a, b), both orientations
    for k in (2..=amax).filter(|_| do_e2) {
        for a in 1..k {
            let b = k - a;
            let (a8, b8) = (a as u8, b as u8);
            let lab = [a8, b8];
            let wr = ex(&lab, &[e(0, b), f(0, b)]);
            let wl = ex(&lab, &[f(0, a), e(0, a)]);
            let emp = ex(&lab, &[]);
            let d = (a * b) as i32;
            let get = |s: &Expr, t: &Expr, lib: &mut TmplLib<F>, reqs: &mut Vec<Req>| -> TmplId {
                let ts = lib.hom_basis(&win(s), &win(t), d);
                reqs.push(Req::Hom(win(s), win(t), d));
                assert_eq!(ts.len(), 1, "A7: dim Hom^{} = {} for {:?} -> {:?}", d, ts.len(), s.rungs, t.rungs);
                ts[0]
            };
            let te = [get(&wr, &emp, &mut lib, &mut reqs), get(&emp, &wr, &mut lib, &mut reqs), get(&wl, &emp, &mut lib, &mut reqs), get(&emp, &wl, &mut lib, &mut reqs)];
            let mut ar = Arena::new(&lab, std::mem::replace(&mut lib, TmplLib::new(None)));
            let (ri, li, ei) = (ar.expr(wr.clone()), ar.expr(wl.clone()), ar.expr(emp.clone()));
            let kinds = [(ElemKind::Eps, ri, ei), (ElemKind::Eta, ei, ri), (ElemKind::EpsL, li, ei), (ElemKind::EtaL, ei, li)];
            let pts = points::<F>(k, 700 + k as u64 * 10 + a as u64, 3, 0);
            let mut ms = vec![];
            for (n, &(kind, s, t)) in kinds.iter().enumerate() {
                let cf = ar.add(MapNode::Elem { kind, src: s, tgt: t, at: 0, upright: 0, scalar: F::one() }, "A7").unwrap();
                assert!(!ar.info(cf).certified, "uncertified before certify_elem");
                let xm = ar.add(MapNode::Local { src: s, tgt: t, at: 0, first: 0, tmpl: te[n] }, "A7 exact").unwrap();
                assert_eq!(ar.info(cf).deg, ar.info(xm).deg, "A7: degree of {:?}", kind);
                let c = ratio(&ar, cf, xm, &pts).unwrap_or_else(|e| panic!("A7 {:?} ({},{}): {}", kind, a, b, e));
                assert_eq!(c, F::one(), "A7: constant of {:?} at labels ({},{}) is {}, expected 1", kind, a, b, c);
                ar.certify_elem(kind, (a8, b8), c);
                consts.push(format!("{:?}({},{})={}", kind, a, b, c));
                // single nonzero entry
                for p in &pts {
                    assert_eq!(nnz(&eval_dense(&ar, cf, p)), 1, "A7: {:?} ({},{}) has more than one nonzero entry", kind, a, b);
                }
                ms.push(cf);
            }
            // products: Eps∘Eta = e(B,A), EpsL∘EtaL = e(A,B)
            let am = low_mask(a);
            let bm = low_mask(b) << a;
            let pe = ar.add(MapNode::Comp(ms[0], ms[1]), "eps eta").unwrap();
            let pl = ar.add(MapNode::Comp(ms[2], ms[3]), "epsL etaL").unwrap();
            for p in &pts {
                let lay = Layout::of(ar.flows(ei), ar.flows(ei));
                assert_eq!(lay.cosets.len(), 1);
                let mut ev = Eval::new(&ar, p.clone());
                let c = lay.cosets[0];
                assert_eq!(ev.block(pe, c).a, vec![e_pair(bm, am, &p.x)], "A7: Eps∘Eta != e(B,A) at ({},{})", a, b);
                assert_eq!(ev.block(pl, c).a, vec![e_pair(am, bm, &p.x)], "A7: EpsL∘EtaL != e(A,B) at ({},{})", a, b);
            }
            // certified now (new nodes), genuine
            let again = ar.add(MapNode::Elem { kind: ElemKind::Eps, src: ri, tgt: ei, at: 0, upright: 0, scalar: F::one() }, "A7").unwrap();
            assert!(ar.info(again).certified);
            ar.assert_genuine(pe).unwrap();
            for &t in &te {
                cert_template(ar.lib_ref(), t, 90 + t as u64);
                certs += 1;
            }
            lib = ar.into_lib();
            reqs_all += 4;
            if amax > 6 {
                eprintln!("[A7 {}] E2 ({},{}) done at {:.1}s", tag, a, b, t0.elapsed().as_secs_f64());
            }
        }
    }
    // E3 digons (c; b), both orientations
    for c in (2..=amax).filter(|_| do_e3) {
        for b in 1..c {
            let a = c - b;
            let d = -((a * b) as i32);
            for mirror in [false, true] {
                let lab: Vec<u8> = if mirror { vec![0, c as u8] } else { vec![c as u8, 0] };
                let dg = if mirror { ex(&lab, &[e(0, b), f(0, b)]) } else { ex(&lab, &[f(0, b), e(0, b)]) };
                let emp = ex(&lab, &[]);
                let ti = lib.hom_basis(&win(&emp), &win(&dg), d);
                let tt = lib.hom_basis(&win(&dg), &win(&emp), d);
                reqs.push(Req::Hom(win(&emp), win(&dg), d));
                reqs.push(Req::Hom(win(&dg), win(&emp), d));
                assert_eq!((ti.len(), tt.len()), (1, 1), "A7: digon ({};{}) Hom dims", c, b);
                let mut ar = Arena::new(&lab, std::mem::replace(&mut lib, TmplLib::new(None)));
                let (di, ei) = (ar.expr(dg.clone()), ar.expr(emp.clone()));
                let (ki, kt) = if mirror { (ElemKind::IotaL, ElemKind::TraceL) } else { (ElemKind::Iota, ElemKind::Trace) };
                let pts = points::<F>(c, 800 + c as u64 * 10 + b as u64 + if mirror { 5 } else { 0 }, 3, 0);
                let io = ar.add(MapNode::Elem { kind: ki, src: ei, tgt: di, at: 0, upright: 0, scalar: F::one() }, "A7").unwrap();
                let tr = ar.add(MapNode::Elem { kind: kt, src: di, tgt: ei, at: 0, upright: 0, scalar: F::one() }, "A7").unwrap();
                let xi = ar.add(MapNode::Local { src: ei, tgt: di, at: 0, first: 0, tmpl: ti[0] }, "A7 exact").unwrap();
                let xt = ar.add(MapNode::Local { src: di, tgt: ei, at: 0, first: 0, tmpl: tt[0] }, "A7 exact").unwrap();
                for (cf, xm, kind) in [(io, xi, ki), (tr, xt, kt)] {
                    assert_eq!(ar.info(cf).deg, d);
                    let cst = ratio(&ar, cf, xm, &pts).unwrap_or_else(|e| panic!("A7 {:?} ({};{}): {}", kind, c, b, e));
                    assert_eq!(cst, F::one(), "A7: constant of {:?} at digon ({};{}) is {}, expected 1", kind, c, b, cst);
                    ar.certify_elem(kind, (c as u8, b as u8), cst);
                    consts.push(format!("{:?}({};{})={}", kind, c, b, cst));
                }
                // (iii) Trace∘Iota = 0; Trace∘s_nu(S)∘Iota = 0 for nu strictly inside the b x a box;
                // nu = (a^b): a nonzero constant
                let ti0 = ar.add(MapNode::Comp(tr, io), "tr io").unwrap();
                let s_up = if mirror { 0 } else { 1 };
                let mut top_const = None;
                for nu in partitions_in_box(b, a) {
                    let size: usize = nu.iter().map(|&x| x as usize).sum();
                    let dot = ar.add(MapNode::Dot { e: di, level: 1, upright: s_up, nu: nu.clone() }, "A7 dot").unwrap();
                    let m1 = ar.add(MapNode::Comp(dot, io), "").unwrap();
                    let m2 = ar.add(MapNode::Comp(tr, m1), "").unwrap();
                    for p in &pts {
                        let v = eval_dense(&ar, m2, p)[0][0];
                        if size < a * b {
                            assert!(v.is_zero(), "A7: Trace∘s_{:?}∘Iota != 0 at ({};{}) {}", nu, c, b, if mirror { "L" } else { "" });
                        } else {
                            assert!(!v.is_zero());
                            match top_const {
                                None => top_const = Some(v),
                                Some(t) => assert_eq!(t, v, "A7: Trace∘s_(a^b)∘Iota not constant at ({};{})", c, b),
                            }
                        }
                    }
                }
                for p in &pts {
                    assert!(eval_dense(&ar, ti0, p)[0][0].is_zero(), "A7: Trace∘Iota != 0");
                }
                consts.push(format!("Trace∘s_({}^{})∘Iota{}({};{})={}", a, b, if mirror { "L" } else { "" }, c, b, top_const.unwrap()));
                for &t in ti.iter().chain(&tt) {
                    cert_template(ar.lib_ref(), t, 120 + t as u64);
                    certs += 1;
                }
                lib = ar.into_lib();
                reqs_all += 2;
                if amax > 6 {
                    eprintln!("[A7 {}] E3 ({};{}){} done at {:.1}s", tag, c, b, if mirror { " L" } else { "" }, t0.elapsed().as_secs_f64());
                }
            }
        }
    }
    roundtrip(&lib, &reqs, &format!("A7 {}", tag));
    eprintln!("[A7 {}] {} exact maps; {} certificates pass; {:.1}s", tag, reqs_all, certs, t0.elapsed().as_secs_f64());
    eprintln!("[A7 {}] constants: {}", tag, consts.join(" "));
}

#[test]
fn a7_e2_e3_closed_forms() {
    e2e3_checks(6, "le6", true, true);
}

#[test]
#[ignore = "A7(iv) E2: a + b <= 9 ((3,3,3) windows); run before G1(b); cloud job (> 10 min on the laptop)"]
fn a7_iv_e2_up_to_9() {
    e2e3_checks(9, "e2le9", true, false);
}

#[test]
#[ignore = "A7(iv) E3: c <= 9 ((3,3,3) windows); run before G1(b)"]
fn a7_iv_e3_up_to_9() {
    e2e3_checks(9, "e3le9", false, true);
}

// ------------------------------------------------------------------ A8

#[test]
fn a8_assoc() {
    let mut lib: TmplLib<F> = TmplLib::new(Some(cache_dir("a8")));
    let mut reqs = vec![];
    let mut consts = vec![];
    for (a, b, c) in [(1u8, 1u8, 1u8), (2, 1, 1), (1, 2, 1), (1, 1, 2), (2, 2, 2)] {
        let (bu, cu) = (b as usize, c as usize);
        for mirror in [false, true] {
            let (lab, s_r, t_r) = if !mirror {
                (vec![a, b, c], vec![e(0, bu), e(1, cu), e(0, cu)], vec![e(1, cu), e(0, bu + cu)])
            } else {
                (vec![c, b, a], vec![f(1, bu), f(0, cu), f(1, cu)], vec![f(0, cu), f(1, bu + cu)])
            };
            let (sx, tx) = (ex(&lab, &s_r), ex(&lab, &t_r));
            let mut s = Schur::new(&lab);
            let (cs, ct) = (s.ladder_char(&sx.rungs), s.ladder_char(&tx.rungs));
            assert_eq!(s.hom_gd(&cs, &ct, 0, 0), vec![1], "A8: dim Hom^0 for ({},{},{}) mirror {}", a, b, c, mirror);
            let ts = lib.hom_basis(&win(&sx), &win(&tx), 0);
            reqs.push(Req::Hom(win(&sx), win(&tx), 0));
            assert_eq!(ts.len(), 1);
            let mut ar = Arena::new(&lab, std::mem::replace(&mut lib, TmplLib::new(None)));
            let (si, ti) = (ar.expr(sx.clone()), ar.expr(tx.clone()));
            let fw = ar.add(MapNode::Elem { kind: ElemKind::Assoc, src: si, tgt: ti, at: 0, upright: 0, scalar: F::one() }, "A8").unwrap();
            let bw = ar.add(MapNode::Elem { kind: ElemKind::Assoc, src: ti, tgt: si, at: 0, upright: 0, scalar: F::one() }, "A8").unwrap();
            let xm = ar.add(MapNode::Local { src: si, tgt: ti, at: 0, first: 0, tmpl: ts[0] }, "A8 exact").unwrap();
            let pts = points::<F>(lab.iter().map(|&x| x as usize).sum(), 81, 3, 0);
            let cst = ratio(&ar, fw, xm, &pts).unwrap_or_else(|e| panic!("A8 ({},{},{}): {}", a, b, c, e));
            let key = ar.elem_key(fw).unwrap().1;
            ar.certify_elem(ElemKind::Assoc, key, cst);
            consts.push(format!("Assoc({},{},{}{})={}", a, b, c, if mirror { " F" } else { "" }, cst));
            let rt = ar.add(MapNode::Comp(bw, fw), "A8 round trip").unwrap();
            let id = ar.add(MapNode::Id(si), "id").unwrap();
            for p in &pts {
                assert_eq!(eval_dense(&ar, rt, p), eval_dense(&ar, id, p), "A8: Assoc ∘ reverse != id");
            }
            cert_template(ar.lib_ref(), ts[0], 33);
            lib = ar.into_lib();
        }
    }
    roundtrip(&lib, &reqs, "A8");
    eprintln!("[A8] constants: {}", consts.join(" "));
    assert!(consts.iter().all(|s| s.ends_with("=1")), "A8: an Assoc constant differs from 1: {:?}", consts);
}

// ------------------------------------------------------------------ A9 (+ A16 WinLoc::dense)

fn rickard_checks(cases: &[(u8, u8)], tag: &str) -> Vec<(String, f64)> {
    let mut lib: TmplLib<F> = TmplLib::new(Some(cache_dir(&format!("a9_{}", tag))));
    let mut reqs = vec![];
    let mut times = vec![];
    for &(c, d) in cases {
        for positive in [true, false] {
            let t0 = Instant::now();
            let (wins, ids) = lib.rickard(c, d, positive);
            let dt = t0.elapsed().as_secs_f64();
            times.push((format!("({},{}) {}", c, d, if positive { "+" } else { "-" }), dt));
            reqs.push(Req::Rick(c, d, positive));
            let lab = [c, d];
            let mut s = Schur::new(&lab);
            let chars: Vec<Elt> = wins.iter().map(|w| s.ladder_char(&w.rungs)).collect();
            for j in 0..ids.len() {
                let (a, b) = if positive { (j, j + 1) } else { (j + 1, j) };
                let gd = s.hom_gd(&chars[a], &chars[b], -60, 1);
                assert!(gd[..61].iter().all(|&x| x == 0) && gd[61] == 1, "A9: Hom^<=1 dims of ({},{}) map {}", c, d, j);
            }
            let mut ar = Arena::new(&lab, std::mem::replace(&mut lib, TmplLib::new(None)));
            let ws: Vec<ExprId> = wins.iter().map(|w| ar.expr(w.expr())).collect();
            let ms: Vec<MapId> = ids
                .iter()
                .enumerate()
                .map(|(j, &t)| {
                    let (a, b) = if positive { (j, j + 1) } else { (j + 1, j) };
                    ar.add(MapNode::Local { src: ws[a], tgt: ws[b], at: 0, first: 0, tmpl: t }, "A9").unwrap()
                })
                .collect();
            let k = (c + d) as usize;
            let pts = points::<F>(k, 900 + c as u64 * 7 + d as u64, 2, 0);
            // d∘d = 0
            for j in 0..ms.len().saturating_sub(1) {
                let dd = if positive { ar.add(MapNode::Comp(ms[j + 1], ms[j]), "dd").unwrap() } else { ar.add(MapNode::Comp(ms[j], ms[j + 1]), "dd").unwrap() };
                for p in &pts {
                    assert!(Eval::new(&ar, p.clone()).is_zero(dd), "A9: d∘d != 0 for ({},{}) {} at j = {}", c, d, positive, j);
                }
            }
            for (j, &t) in ids.iter().enumerate() {
                let tm = ar.lib_ref().get(t);
                let (se, te) = (tm.src.expr(), tm.tgt.expr());
                for p in &pts {
                    // block diagonality of the dense localize_map output
                    let (fp_s, fp_t, fm) = localize_map(&se.sbim::<F>(), &te.sbim::<F>(), &tm.phi, &p.x);
                    for (q, row) in fm.iter().enumerate() {
                        for (pi, v) in row.iter().enumerate() {
                            if !v.is_zero() {
                                assert_eq!(fp_t[q].top, fp_s[pi].top, "A9: nonzero entry between different top cosets");
                            }
                        }
                    }
                    // A16: WinLoc::dense = reordered localize_map
                    assert_eq!(ar.lib_ref().localize(t, &p.x).dense(), oracle_sorted(&se, &te, &tm.phi, &p.x), "A16: WinLoc::dense of Rickard ({},{}) map {}", c, d, j);
                    assert_eq!(eval_dense(&ar, ms[j], p), oracle_sorted(&se, &te, &tm.phi, &p.x), "A9: Local of Rickard map {}", j);
                }
                // equivariance under a transposition of two roots of one bottom block
                let blk = if c >= 2 { Some(0usize) } else if d >= 2 { Some(1) } else { None };
                if let Some(bi) = blk {
                    let o = if bi == 0 { 0 } else { c as usize };
                    let mut sigma: Vec<usize> = (0..k).collect();
                    sigma.swap(o, o + 1);
                    let (fs, ft) = (ar.flows(ar.info(ms[j]).src).clone(), ar.flows(ar.info(ms[j]).tgt).clone());
                    for p in &pts {
                        let sx: Vec<F> = (0..k).map(|i| p.x[sigma[i]]).collect();
                        let mut e1 = Eval::new(&ar, p.clone());
                        let mut e2 = Eval::new(&ar, Point { id: 99, x: sx });
                        for &cc in fs.support() {
                            let b1 = e1.block(ms[j], cc);
                            let sc = fs.top.rank(&fs.top.act(&sigma, &fs.top.unrank(cc)));
                            let b2 = e2.block(ms[j], sc);
                            let (s0, t0r) = (fs.range(cc).start, ft.range(cc).start);
                            for col in 0..b1.cols {
                                let sf: Vec<u64> = fs.flow(s0 + col).iter().map(|&m| permute_mask(&sigma, m)).collect();
                                let col2 = fs.find(sc, &sf).unwrap() as usize;
                                for row in 0..b1.rows {
                                    let tf: Vec<u64> = ft.flow(t0r + row).iter().map(|&m| permute_mask(&sigma, m)).collect();
                                    let row2 = ft.find(sc, &tf).unwrap() as usize;
                                    assert_eq!(b1.get(row as usize, col as usize), b2.get(row2, col2), "A9: equivariance of ({},{}) map {}", c, d, j);
                                }
                            }
                        }
                    }
                }
                cert_template(ar.lib_ref(), t, 300 + t as u64);
            }
            lib = ar.into_lib();
        }
    }
    roundtrip(&lib, &reqs, &format!("A9 {}", tag));
    times
}

#[test]
fn a9_rickard_templates() {
    let t = rickard_checks(&[(1, 1), (2, 1), (1, 2), (2, 2)], "small");
    eprintln!("[A9] local_cx times: {:?}", t);
}

#[test]
fn a9_rickard_three_labels() {
    // default if each local_cx call takes < 10 s (measured; see the report)
    let t = rickard_checks(&[(3, 1), (1, 3), (3, 2), (2, 3)], "three");
    eprintln!("[A9] local_cx times: {:?}", t);
    for (name, s) in &t {
        assert!(*s < 10.0, "A9: local_cx {} took {:.1}s >= 10 s (move it to the ignored set)", name, s);
    }
}

#[test]
#[ignore = "A9: (3,3) Rickard templates (about 2 min per sign per prime)"]
fn a9_rickard_33() {
    let t = rickard_checks(&[(3, 3)], "33");
    eprintln!("[A9] (3,3) local_cx times: {:?}", t);
}

// ------------------------------------------------------------------ A10

#[test]
fn a10_certificates() {
    // closed forms with dots, as WinLoc functions of the window roots
    let mut reports = vec![];
    for (a, b) in [(1u8, 1u8), (2, 1), (2, 2)] {
        let lab = [a, b];
        let mut ar = Arena::<F>::new(&lab, TmplLib::new(None));
        let w = ex(&lab, &[e(0, b as usize), f(0, b as usize)]);
        let emp = ex(&lab, &[]);
        let (wi, ei) = (ar.expr(w.clone()), ar.expr(emp.clone()));
        let eps = ar.add(MapNode::Elem { kind: ElemKind::Eps, src: wi, tgt: ei, at: 0, upright: 0, scalar: F::one() }, "").unwrap();
        let eta = ar.add(MapNode::Elem { kind: ElemKind::Eta, src: ei, tgt: wi, at: 0, upright: 0, scalar: F::one() }, "").unwrap();
        let dot = ar.add(MapNode::Dot { e: wi, level: 1, upright: 0, nu: vec![1] }, "").unwrap();
        let eps_dot = ar.add(MapNode::Comp(eps, dot), "").unwrap();
        let dot_eta = ar.add(MapNode::Comp(dot, eta), "").unwrap();
        for (m, s, t) in [(eps, &w, &emp), (eta, &emp, &w), (eps_dot, &w, &emp), (dot_eta, &emp, &w)] {
            let d = ar.info(m).deg;
            let r = poly_certificate::<F>(&win(s), &win(t), d, &mut |y: &[F]| Eval::new(&ar, Point { id: 0, x: y.to_vec() }).winloc(m), 5).unwrap_or_else(|e| panic!("A10: closed form ({},{}) failed: {}", a, b, e));
            reports.push(r);
        }
        // (a) Eta multiplied at its coset by 1/(x0 - x1): rejected, naming the entry
        if (a, b) == (2, 2) {
            let r = poly_certificate::<F>(
                &win(&emp),
                &win(&w),
                ar.info(eta).deg,
                &mut |y: &[F]| {
                    let mut wl = Eval::new(&ar, Point { id: 0, x: y.to_vec() }).winloc(eta);
                    let s = y[0].sub(&y[1]).inv();
                    for v in wl.vals.iter_mut() {
                        *v = v.mul(&s);
                    }
                    wl
                },
                5,
            );
            let err = r.expect_err("A10(a): a non-polynomial Eta passed");
            assert!(err.contains("entry"), "A10(a): error does not name the entry: {}", err);
            eprintln!("[A10] (a) rejected: {}", err);
        }
    }
    for (c, bq) in [(2u8, 1usize), (3, 1), (4, 2)] {
        let lab = [c, 0];
        let mut ar = Arena::<F>::new(&lab, TmplLib::new(None));
        let dg = ex(&lab, &[f(0, bq), e(0, bq)]);
        let emp = ex(&lab, &[]);
        let (di, ei) = (ar.expr(dg.clone()), ar.expr(emp.clone()));
        let io = ar.add(MapNode::Elem { kind: ElemKind::Iota, src: ei, tgt: di, at: 0, upright: 0, scalar: F::one() }, "").unwrap();
        let tr = ar.add(MapNode::Elem { kind: ElemKind::Trace, src: di, tgt: ei, at: 0, upright: 0, scalar: F::one() }, "").unwrap();
        let dot = ar.add(MapNode::Dot { e: di, level: 1, upright: 1, nu: vec![1] }, "").unwrap();
        let tr_dot = ar.add(MapNode::Comp(tr, dot), "").unwrap();
        let dot_io = ar.add(MapNode::Comp(dot, io), "").unwrap();
        for (m, s, t) in [(io, &emp, &dg), (tr, &dg, &emp), (tr_dot, &dg, &emp), (dot_io, &emp, &dg)] {
            let d = ar.info(m).deg;
            let r = poly_certificate::<F>(&win(s), &win(t), d, &mut |y: &[F]| Eval::new(&ar, Point { id: 0, x: y.to_vec() }).winloc(m), 6).unwrap_or_else(|e| panic!("A10: digon ({};{}) failed: {}", c, bq, e));
            reports.push(r);
        }
        // (c) a sign flipped on one non-identity flow (Iota, Trace: nonzero on every flow)
        for (m, s, t) in [(io, &emp, &dg), (tr, &dg, &emp)] {
            let d = ar.info(m).deg;
            let r = poly_certificate::<F>(
                &win(s),
                &win(t),
                d,
                &mut |y: &[F]| {
                    let mut wl = Eval::new(&ar, Point { id: 0, x: y.to_vec() }).winloc(m);
                    let last = wl.vals.len() - 1;
                    wl.vals[last] = wl.vals[last].neg();
                    wl
                },
                6,
            );
            let err = r.expect_err("A10(c): a sign-flipped closed form passed");
            assert!(err.contains("entry"), "A10(c): error does not name the entry: {}", err);
            eprintln!("[A10] (c) digon ({};{}) rejected: {}", c, bq, err);
        }
    }
    // (c) literally on Eps: Eps vanishes off the identity flow, so the flip is a no-op (recorded)
    {
        let lab = [1u8, 1];
        let mut ar = Arena::<F>::new(&lab, TmplLib::new(None));
        let w = ex(&lab, &[e(0, 1), f(0, 1)]);
        let emp = ex(&lab, &[]);
        let (wi, ei) = (ar.expr(w.clone()), ar.expr(emp.clone()));
        let eps = ar.add(MapNode::Elem { kind: ElemKind::Eps, src: wi, tgt: ei, at: 0, upright: 0, scalar: F::one() }, "").unwrap();
        let wl = Eval::new(&ar, Point { id: 0, x: vec![F::from_i64(3), F::from_i64(5)] }).winloc(eps);
        assert_eq!(wl.vals.len(), 1, "Eps has a single nonzero (the identity flow)");
        eprintln!("[A10] (c) literal 'Eps with a sign flipped on a non-identity flow' is a no-op: Eps has 1 nonzero of {} flows", wl.ncols());
    }
    // (b) Rickard d_0 of (2,1) rescaled at every point so that its first nonzero entry is 1
    let mut lib: TmplLib<F> = TmplLib::new(None);
    let (_, ids) = lib.rickard(2, 1, true);
    let t = ids[0];
    let tm = lib.get(t).clone();
    let r = poly_certificate::<F>(
        &tm.src,
        &tm.tgt,
        tm.deg,
        &mut |y: &[F]| {
            let mut wl = lib.localize(t, y);
            let s = wl.vals[0].inv();
            for v in wl.vals.iter_mut() {
                *v = v.mul(&s);
            }
            wl
        },
        8,
    );
    let err = r.expect_err("A10(b): a point-normalized Rickard map passed");
    assert!(err.contains("entry"), "A10(b): error does not name the entry: {}", err);
    eprintln!("[A10] (b) rejected: {}", err);
    // T and max |Mono_D| for the (2,2) and (3,2) Rickard maps
    for (c, d) in [(2u8, 2u8), (3, 2)] {
        for (j, &t) in lib.rickard(c, d, true).1.clone().iter().enumerate() {
            let r = cert_template(&lib, t, 9);
            eprintln!("[A10] Rickard ({},{}) d_{}: T = {}, max |Mono_D| = {}, distinct D = {}, entries = {}, redraws = {}, {:.3}s", c, d, j, r.points - 2, r.max_monos, r.distinct_degrees, r.entries, r.redraws, r.secs);
        }
    }
    eprintln!("[A10] {} closed-form certificates pass (max |Mono_D| {})", reports.len(), reports.iter().map(|r| r.max_monos).max().unwrap());
}

// ------------------------------------------------------------------ A11

/// idempotent e = ι π / G on X for the summand L_q{-1}; checks reduced coordinates at every coset
fn idempotent_check(lab: &[u8], lq: &Expr, x: &Expr) {
    let mut lib: TmplLib<F> = TmplLib::new(None);
    let t0 = Instant::now();
    let iota = lib.hom_basis(&win(lq), &win(x), -1);
    let pis = lib.hom_basis(&win(x), &win(lq), 1);
    eprintln!("[A11] {:?}: Hom^-1 dim {}, Hom^1 dim {}, {:.1}s", lab, iota.len(), pis.len(), t0.elapsed().as_secs_f64());
    let mut s = Schur::new(lab);
    let (cq, cx) = (s.ladder_char(&lq.rungs), s.ladder_char(&x.rungs));
    assert_eq!(s.hom_gd(&cq, &cx, -1, -1)[0], iota.len() as i128, "A11: Hom^-1 dimension");
    assert!(!iota.is_empty());
    let mut ar = Arena::new(lab, lib);
    let (qi, xi) = (ar.expr(lq.clone()), ar.expr(x.clone()));
    let io = ar.add(MapNode::Local { src: qi, tgt: xi, at: 0, first: 0, tmpl: iota[0] }, "iota").unwrap();
    let pts = points::<F>(x.k(), 1111, 2, 0);
    let mut chosen = None;
    for &t in &pis {
        let pm = ar.add(MapNode::Local { src: xi, tgt: qi, at: 0, first: 0, tmpl: t }, "pi").unwrap();
        let g = ar.add(MapNode::Comp(pm, io), "G").unwrap();
        // G = π∘ι is a scalar (End^0(L_q) = k): the same at every coset and point
        let mut gv: Option<F> = None;
        let mut scalar = true;
        for p in &pts {
            let lm = Eval::new(&ar, p.clone()).all(g);
            for i in 0..lm.lay.cosets.len() {
                let blk = lm.block(i);
                let n = lm.lay.rows[i] as usize;
                for r in 0..n {
                    for c in 0..n {
                        let v = blk[r * n + c];
                        if r == c {
                            match gv {
                                None => gv = Some(v),
                                Some(g0) => scalar &= g0 == v,
                            }
                        } else {
                            scalar &= v.is_zero();
                        }
                    }
                }
            }
        }
        assert!(scalar, "A11: π∘ι is not a scalar");
        if let Some(g0) = gv {
            if !g0.is_zero() {
                chosen = Some((pm, g0));
                break;
            }
        }
    }
    let (pm, g0) = chosen.expect("A11: no π with π∘ι != 0");
    let ip = ar.add(MapNode::Comp(io, pm), "ι π").unwrap();
    let em = ar.add(MapNode::Lin(vec![(g0.inv(), ip)]), "e = ι π / G").unwrap();
    let xf = ar.flows(xi).clone();
    let bq = s.expand(&cq);
    assert_eq!(bq.len(), 1, "L_q is indecomposable");
    let bqe = s.canonical_basis(&cq.mu).b[&bq[0].0].clone();
    for p in &pts {
        let mut ev = Eval::new(&ar, p.clone());
        for &c in xf.support() {
            let ec = ev.block(em, c);
            let red = reduced_of(&ec).unwrap_or_else(|e| panic!("A11: reduced_of at coset {}: {}", c, e));
            let tab = xf.top.table(&xf.top.unrank(c));
            let tid: Tid = s.intern(tab.clone());
            let mc = bqe.coef(tid).at1() as usize;
            assert_eq!(red.m, mc, "A11: rank E_c != m_c(B_q) at table {}", tab.fmt());
            assert_eq!(red.q.mul(&red.p), Block::identity(red.m));
            assert_eq!(red.p.mul(&red.q), *ec);
        }
    }
    eprintln!("[A11] {:?}: reduced coordinates at {} cosets, G = {}", lab, xf.support().len(), g0);
}

#[test]
fn a11_reduced_coordinates() {
    let bs = ex(&[1, 1], &[e(0, 1), f(0, 1)]);
    let bss = ex(&[1, 1], &[e(0, 1), f(0, 1), e(0, 1), f(0, 1)]);
    idempotent_check(&[1, 1], &bs, &bss);
}

/// default: hom_space here takes about 0.1 s (the spec makes it ignored only above 30 s)
#[test]
fn a11_reduced_coordinates_22() {
    let x = ex(&[2, 2], &[e(0, 1), f(0, 1), e(0, 1), f(0, 1)]);
    let mut s = Schur::new(&[2, 2]);
    let cx = s.ladder_char(&x.rungs);
    // the summand with multiplicity v + v^-1
    let summ = s.expand(&cx);
    let (q, m) = summ.iter().find(|(_, m)| m.coef(-1) == 1 && m.coef(1) == 1).expect("a summand B_q{±1}").clone();
    let mu = s.table(q).cols();
    let lq = ex(&[2, 2], &s.canonical_basis(&mu).ladder_minflow[&q].clone());
    eprintln!("[A11] (2,2): summand {} with multiplicity {}, L_q = {:?}", s.table(q).fmt(), m.fmt_ec(), lq.rungs);
    idempotent_check(&[2, 2], &lq, &x);
}

// ------------------------------------------------------------------ A12

#[test]
fn a12_point_independence() {
    let lab = [2u8, 2];
    let w = Window::new(&lab, &[e(0, 1)]);
    let mut lib: TmplLib<F> = TmplLib::new(None);
    let ts = lib.hom_basis(&w, &w, 2);
    assert!(ts.len() >= 2, "need several dots");
    let mut ar = Arena::new(&lab, lib);
    let xi = ar.expr(w.expr());
    let ms: Vec<MapId> = ts.iter().map(|&t| ar.add(MapNode::Local { src: xi, tgt: xi, at: 0, first: 0, tmpl: t }, "").unwrap()).collect();
    let coefs: Vec<F> = (0..ms.len()).map(|i| F::from_i64(3 * i as i64 - 7)).collect();
    let lin = ar.add(MapNode::Lin(coefs.iter().cloned().zip(ms.iter().cloned()).collect()), "lin").unwrap();
    let pts = points::<F>(4, 1212, 3, 0);
    for p in &pts {
        for (&m, &t) in ms.iter().zip(&ts) {
            assert_eq!(eval_dense(&ar, m, p), oracle_sorted(&w.expr(), &w.expr(), &ar.lib_ref().get(t).phi, &p.x), "A12: template at point {}", p.id);
        }
        let parts: Vec<Vec<Vec<F>>> = ms.iter().map(|&m| eval_dense(&ar, m, p)).collect();
        let (nr, nc) = (parts[0].len(), parts[0][0].len());
        let want: Vec<Vec<F>> = (0..nr).map(|i| (0..nc).map(|j| parts.iter().zip(&coefs).fold(F::zero(), |s, (pm, c)| s.add(&c.mul(&pm[i][j])))).collect()).collect();
        assert_eq!(eval_dense(&ar, lin, p), want, "A12: Lin does not commute with evaluation");
    }
    // the same map at two points differs (non-constant) but each matches the oracle
    assert_ne!(eval_dense(&ar, ms[0], &pts[0]), eval_dense(&ar, ms[0], &pts[1]));
}

// ------------------------------------------------------------------ A13

#[test]
fn a13_cache_header_and_corruption() {
    let dir = cache_dir("a13");
    let w1 = Window::new(&[1, 1], &[e(0, 1), f(0, 1)]);
    let w0 = Window::new(&[1, 1], &[]);
    let mut lib: TmplLib<F> = TmplLib::new(Some(dir.clone()));
    let ids = lib.hom_basis(&w1, &w0, 1);
    let _ = lib.rickard(2, 1, true);
    lib.save().unwrap();
    let path = lib.hom_cache_path(&w1, &w0, 1).unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    let header = text.lines().next().unwrap().to_string();
    assert!(header.starts_with("krhom-tmpl v1 "), "header {:?}", header);
    assert!(header.contains(&format!("field={}", F::name())), "header without the field: {:?}", header);
    // a fresh library loads without computing
    let mut l2: TmplLib<F> = TmplLib::new(Some(dir.clone()));
    let ids2 = l2.hom_basis(&w1, &w0, 1);
    assert_eq!(l2.stats().computed, 0);
    assert!(lib.get(ids[0]).phi.sub(&l2.get(ids2[0]).phi).is_zero());
    // corrupted header: rejected (read_cache_file Err; the library recomputes and counts it)
    let bad = text.replacen("krhom-tmpl v1", "krhom-tmpl v0", 1);
    std::fs::write(&path, bad).unwrap();
    assert!(read_cache_file::<F>(&path, &header).is_err(), "corrupted header accepted");
    let mut l3: TmplLib<F> = TmplLib::new(Some(dir.clone()));
    let ids3 = l3.hom_basis(&w1, &w0, 1);
    assert_eq!((l3.stats().rejected, l3.stats().computed), (1, 1), "corrupted header not rejected");
    assert!(lib.get(ids[0]).phi.sub(&l3.get(ids3[0]).phi).is_zero());
    // a different field does not read F31 files (the name and header carry the field)
    let mut l4: TmplLib<krhom::field::F32b> = TmplLib::new(Some(dir));
    let _ = l4.hom_basis(&w1, &w0, 1);
    assert_eq!(l4.stats().loaded, 0, "an F32b library loaded an F31 cache file");
}

// ------------------------------------------------------------------ A14

#[test]
fn a14_determinism() {
    let run = || -> (u64, Vec<u64>, Vec<Vec<u64>>) {
        use std::hash::{Hash, Hasher};
        let lab = [2u8, 1];
        let mut lib: TmplLib<F> = TmplLib::new(None);
        let (wins, ids) = lib.rickard(2, 1, true);
        let mut ar = Arena::new(&lab, lib);
        let w0 = ar.expr(wins[0].expr());
        let w1 = ar.expr(wins[1].expr());
        let m = ar.add(MapNode::Local { src: w0, tgt: w1, at: 0, first: 0, tmpl: ids[0] }, "").unwrap();
        let pts = points::<F>(3, 1414, 2, 0);
        let mut h = std::collections::hash_map::DefaultHasher::new();
        let mut raw = vec![];
        for p in &pts {
            let lm = Eval::new(&ar, p.clone()).all(m);
            for v in &lm.data {
                v.0.hash(&mut h);
            }
            raw.push(lm.data.iter().map(|v| v.0).collect());
        }
        let x222 = block_ladders(&[2, 2, 2], Some(e(0, 1)));
        let fi = FlowIndex::build(&x222[7]);
        let mut v = fi.masks.clone();
        v.extend(fi.start.iter().map(|&s| s as u64));
        v.extend(fi.cosets.iter().cloned());
        (h.finish(), v, raw)
    };
    let (h1, v1, r1) = run();
    let (h2, v2, r2) = run();
    assert_eq!(h1, h2, "A14: block hashes differ between runs");
    assert_eq!(r1, r2);
    assert_eq!(v1, v2, "A14: FlowIndex orders differ between runs");
    // points are deterministic from (seed, id)
    let p1 = points::<F>(5, 3, 4, 10);
    let p2 = points::<F>(5, 3, 2, 12);
    assert_eq!(p1[2].x, p2[0].x);
    assert_eq!(p1[3].x, p2[1].x);
}

// ------------------------------------------------------------------ A15

#[test]
fn a15_arena_checks() {
    let lab = [2u8, 2];
    let mut lib: TmplLib<F> = TmplLib::new(None);
    let w = Window::new(&lab, &[e(0, 1)]);
    let dots = lib.hom_basis(&w, &w, 2);
    let mut ar = Arena::new(&lab, lib);
    let x0 = ar.expr(ex(&lab, &[]));
    let x1 = ar.expr(ex(&lab, &[e(0, 1)]));
    let x2 = ar.expr(ex(&lab, &[e(0, 2), f(0, 2)]));
    let id0 = ar.add(MapNode::Id(x0), "id0").unwrap();
    let id1 = ar.add(MapNode::Id(x1), "id1").unwrap();
    let d1 = ar.add(MapNode::Dot { e: x1, level: 1, upright: 0, nu: vec![1] }, "dot").unwrap();
    // Comp with mismatched expressions
    let err = ar.add(MapNode::Comp(id0, id1), "bad comp").unwrap_err();
    assert!(err.contains("Comp"), "{}", err);
    // Lin with mixed degrees, mixed ends, no terms
    let err = ar.add(MapNode::Lin(vec![(F::one(), id1), (F::one(), d1)]), "bad lin").unwrap_err();
    assert!(err.contains("degree"), "{}", err);
    let err = ar.add(MapNode::Lin(vec![(F::one(), id0), (F::one(), id1)]), "bad lin").unwrap_err();
    assert!(err.contains("ends"), "{}", err);
    assert!(ar.add(MapNode::Lin(vec![]), "empty lin").is_err());
    // Local whose window does not match the words
    let err = ar.add(MapNode::Local { src: x2, tgt: x2, at: 0, first: 0, tmpl: dots[0] }, "bad local").unwrap_err();
    assert!(err.contains("window"), "{}", err);
    let err = ar.add(MapNode::Local { src: x1, tgt: x1, at: 1, first: 0, tmpl: dots[0] }, "bad local").unwrap_err();
    assert!(err.contains("window"), "{}", err);
    // Elem (2,2) Eps before certify_elem: not certified; assert_genuine names kind and labels
    let eps = ar.add(MapNode::Elem { kind: ElemKind::Eps, src: x2, tgt: x0, at: 0, upright: 0, scalar: F::one() }, "eps").unwrap();
    assert!(!ar.info(eps).certified);
    let c = ar.add(MapNode::Comp(id0, eps), "comp").unwrap();
    assert!(!ar.info(c).certified);
    let err = ar.assert_genuine(c).unwrap_err();
    assert!(err.contains("Eps") && err.contains("(2, 2)"), "{}", err);
    ar.certify_elem(ElemKind::Eps, (2, 2), F::one());
    assert_eq!(ar.elem_certified(ElemKind::Eps, (2, 2)), Some(F::one()));
    ar.assert_genuine(c).unwrap();
    let eps2 = ar.add(MapNode::Elem { kind: ElemKind::Eps, src: x2, tgt: x0, at: 0, upright: 0, scalar: F::one() }, "eps").unwrap();
    assert!(ar.info(eps2).certified);
    // wrong window for a closed form
    assert!(ar.add(MapNode::Elem { kind: ElemKind::Eps, src: x1, tgt: x0, at: 0, upright: 0, scalar: F::one() }, "").is_err());
    // drop_index / index: identical rebuild
    let h = |fi: &FlowIndex| (fi.masks.clone(), fi.start.clone(), fi.cosets.clone());
    let before = h(ar.flows(x2));
    ar.drop_index(x2);
    assert!(!ar.has_index(x2));
    let after = h(ar.index(x2));
    assert_eq!(before, after, "A15: rebuilt FlowIndex differs");

    // scratch regions: 50 random nodes over 10 new expressions, keep 3
    let old_nodes = ar.n_nodes();
    let old_infos: Vec<(ExprId, ExprId, i32, String)> = (0..old_nodes as MapId).map(|m| {
        let i = ar.info(m);
        (i.src, i.tgt, i.deg, i.provenance.clone())
    }).collect();
    let old_exprs: Vec<Expr> = (0..ar.n_exprs() as ExprId).map(|i| ar.expr_of(i).clone()).collect();
    let mark = ar.mark();
    let words: Vec<Vec<Rung>> = vec![
        vec![e(0, 2)],
        vec![f(0, 1)],
        vec![f(0, 2)],
        vec![e(0, 1), f(0, 1)],
        vec![e(0, 1), f(0, 2)],
        vec![f(0, 1), e(0, 1)],
        vec![e(0, 1), e(0, 1)],
        vec![f(0, 1), f(0, 1)],
        vec![e(0, 2), f(0, 1)],
        vec![f(0, 2), e(0, 2)],
    ];
    let new_e: Vec<ExprId> = words.iter().map(|w| ar.expr(ex(&lab, w))).collect();
    assert_eq!(ar.n_exprs(), mark.exprs as usize + 10);
    let mut rng = StdRng::seed_from_u64(1515);
    let mut new_m: Vec<MapId> = vec![];
    while new_m.len() < 50 {
        let kind = rng.gen_range(0..4);
        let r = match kind {
            0 => ar.add(MapNode::Id(new_e[rng.gen_range(0..10)]), "scratch id"),
            1 => {
                let ee = new_e[rng.gen_range(0..10)];
                let nr = ar.expr_of(ee).rungs.len();
                let nus = [vec![1u8], vec![2], vec![1, 1]];
                ar.add(MapNode::Dot { e: ee, level: rng.gen_range(0..=nr), upright: rng.gen_range(0..2), nu: nus[rng.gen_range(0..3)].clone() }, "scratch dot")
            }
            2 if !new_m.is_empty() => {
                let g = new_m[rng.gen_range(0..new_m.len())];
                let cands: Vec<MapId> = new_m.iter().cloned().filter(|&h| ar.info(h).src == ar.info(g).tgt).collect();
                let h = cands[rng.gen_range(0..cands.len())];
                ar.add(MapNode::Comp(h, g), "scratch comp")
            }
            _ if !new_m.is_empty() => {
                let g = new_m[rng.gen_range(0..new_m.len())];
                let cands: Vec<MapId> = new_m.iter().cloned().filter(|&h| ar.info(h).src == ar.info(g).src && ar.info(h).tgt == ar.info(g).tgt && ar.info(h).deg == ar.info(g).deg).collect();
                let h = cands[rng.gen_range(0..cands.len())];
                ar.add(MapNode::Lin(vec![(F::from_i64(rng.gen_range(1..9)), g), (F::from_i64(-3), h)]), "scratch lin")
            }
            _ => continue,
        };
        new_m.push(r.unwrap());
    }
    assert_eq!(ar.n_exprs(), mark.exprs as usize + 10, "scratch nodes created no further expressions");
    let keep: Vec<MapId> = new_m[47..].to_vec();
    // expected survivors
    let mut reach = std::collections::BTreeSet::new();
    let mut st = keep.clone();
    while let Some(m) = st.pop() {
        if (m as usize) < mark.nodes as usize || !reach.insert(m) {
            continue;
        }
        match ar.node(m) {
            MapNode::Comp(g, h) => st.extend([*g, *h]),
            MapNode::Lin(t) => st.extend(t.iter().map(|x| x.1)),
            MapNode::Kron(g, _) => st.push(*g),
            _ => {}
        }
    }
    let mut need_e = std::collections::BTreeSet::new();
    for &m in &reach {
        for x in [ar.info(m).src, ar.info(m).tgt] {
            if x >= mark.exprs {
                need_e.insert(x);
            }
        }
    }
    let pts = points::<F>(4, 1516, 2, 0);
    let before_vals: Vec<Vec<Vec<u64>>> = keep.iter().map(|&m| pts.iter().map(|p| Eval::new(&ar, p.clone()).all(m).data.iter().map(|v| v.0).collect()).collect()).collect();
    let before_ends: Vec<(Expr, Expr, i32)> = keep.iter().map(|&m| (ar.expr_of(ar.info(m).src).clone(), ar.expr_of(ar.info(m).tgt).clone(), ar.info(m).deg)).collect();
    let bytes_before = ar.bytes().total();
    let newids = ar.commit(mark, &keep);
    assert_eq!(ar.n_nodes(), mark.nodes as usize + reach.len(), "A15: surviving node count");
    assert_eq!(ar.n_exprs(), mark.exprs as usize + need_e.len(), "A15: surviving expression count");
    for m in 0..old_nodes as MapId {
        let i = ar.info(m);
        assert_eq!((i.src, i.tgt, i.deg, i.provenance.clone()), old_infos[m as usize], "A15: an id below the mark changed");
    }
    for (i, x) in old_exprs.iter().enumerate() {
        assert_eq!(ar.expr_of(i as ExprId), x);
    }
    let after_vals: Vec<Vec<Vec<u64>>> = newids.iter().map(|&m| pts.iter().map(|p| Eval::new(&ar, p.clone()).all(m).data.iter().map(|v| v.0).collect()).collect()).collect();
    assert_eq!(before_vals, after_vals, "A15: renumbered maps evaluate differently");
    for (k, &m) in newids.iter().enumerate() {
        let i = ar.info(m);
        assert_eq!((ar.expr_of(i.src).clone(), ar.expr_of(i.tgt).clone(), i.deg), before_ends[k]);
    }
    let bytes_after = ar.bytes().total();
    assert!(bytes_after < bytes_before, "A15: bytes did not decrease ({} -> {})", bytes_before, bytes_after);
    eprintln!("[A15] commit kept {} of 50 scratch nodes and {} of 10 expressions; bytes {} -> {}", reach.len(), need_e.len(), bytes_before, bytes_after);
    // re-interning: an expression dropped by commit gets a fresh id again
    let again = ar.expr(ex(&lab, &words[0]));
    assert!(again >= mark.exprs);
}

// ------------------------------------------------------------------ A16

fn _ss<T: Send + Sync>() {}

#[test]
fn a16_storage_and_thread_safety() {
    _ss::<Arena<F>>();
    _ss::<TmplLib<F>>();
    _ss::<FlowIndex>();
    _ss::<Layout>();
    _ss::<LocMap<F>>();
    // LocMap bytes and shared layouts
    let lab = [2u8, 1];
    let mut lib: TmplLib<F> = TmplLib::new(None);
    let (wins, ids) = lib.rickard(2, 1, true);
    let mut ar = Arena::new(&lab, lib);
    let w0 = ar.expr(wins[0].expr());
    let w1 = ar.expr(wins[1].expr());
    let m = ar.add(MapNode::Local { src: w0, tgt: w1, at: 0, first: 0, tmpl: ids[0] }, "").unwrap();
    let m2 = ar.add(MapNode::Lin(vec![(F::from_i64(2), m)]), "").unwrap();
    let mut ev = Eval::new(&ar, points::<F>(3, 16, 1, 0).remove(0));
    let (l1, l2) = (ev.all(m), ev.all(m2));
    assert!(Arc::ptr_eq(&l1.lay, &l2.lay), "A16: maps with equal ends do not share their Layout");
    assert_eq!(l1.bytes(), 8 * l1.lay.len as usize + std::mem::size_of::<LocMap<F>>());
    assert!(!l1.is_zero());
    assert!(l1.nonzero_blocks() >= 1);
    // gemm_acc against the naive product
    let mut rng = StdRng::seed_from_u64(1616);
    for _ in 0..100 {
        let (mm, kk, nn) = (rng.gen_range(0..=20usize), rng.gen_range(0..=20usize), rng.gen_range(0..=20usize));
        let rnd = |rng: &mut StdRng, n: usize| -> Vec<F> { (0..n).map(|_| if rng.gen_bool(0.3) { F::zero() } else { F::from_i64(rng.gen_range(-1000..1000)) }).collect() };
        let a = rnd(&mut rng, mm * kk);
        let b = rnd(&mut rng, kk * nn);
        let mut c = rnd(&mut rng, mm * nn);
        let alpha = F::from_i64(rng.gen_range(2..50));
        let mut want = c.clone();
        for i in 0..mm {
            for j in 0..nn {
                let mut s = F::zero();
                for t in 0..kk {
                    s = s.add(&a[i * kk + t].mul(&b[t * nn + j]));
                }
                want[i * nn + j] = want[i * nn + j].add(&alpha.mul(&s));
            }
        }
        gemm_acc(&mut c, &a, &b, mm, kk, nn, alpha);
        assert_eq!(c, want, "A16: gemm_acc differs from the naive product ({}x{}x{})", mm, kk, nn);
    }
}

// ------------------------------------------------------------------ §3.6 performance (release)

#[test]
#[ignore = "§3.6 performance targets: run with --release"]
fn perf_targets() {
    // FlowIndex of the (2,2,2) X with 432 flows
    let x222 = block_ladders(&[2, 2, 2], Some(e(0, 1)));
    let big = x222.iter().max_by_key(|x| x.nflows()).unwrap().clone();
    let t0 = Instant::now();
    let reps = 50;
    let mut fi = FlowIndex::build(&big);
    for _ in 1..reps {
        fi = FlowIndex::build(&big);
    }
    let t_fi = t0.elapsed().as_secs_f64() / reps as f64;
    eprintln!("[perf] FlowIndex (2,2,2) X 432 flows: {:.3} ms, {} B (targets < 2 ms, <= 26 KB)", t_fi * 1e3, fi.bytes());
    // Eval::all of a Local map on the 432-flow ladder: a degree-2 dot template on its last rung
    let lv = big.level(big.rungs.len() - 1);
    let last = *big.rungs.last().unwrap();
    let (lo, hi) = (last.i, last.i + 1);
    let w = Window::new(&lv[lo..=hi], &[Rung { i: 0, dir: last.dir, r: last.r }]);
    let mut lib: TmplLib<F> = TmplLib::new(None);
    let ts = lib.hom_basis(&w, &w, 2);
    let mut ar = Arena::new(&[2, 2, 2], lib);
    let xi = ar.expr(big.clone());
    let m = ar.add(MapNode::Local { src: xi, tgt: xi, at: big.rungs.len() - 1, first: lo, tmpl: ts[0] }, "").unwrap();
    let pts = points::<F>(6, 36, 5, 0);
    let t1 = Instant::now();
    for p in &pts {
        let mut ev = Eval::new(&ar, p.clone());
        let lm = ev.all(m);
        assert_eq!(lm.lay.cosets.len(), fi.support().len());
    }
    let t_all = t1.elapsed().as_secs_f64() / pts.len() as f64;
    eprintln!("[perf] Eval::all of a Local map on the 432-flow ladder: {:.2} ms (target < 20 ms)", t_all * 1e3);
    // template localization at window fibre <= 36: the largest (2,2) Rickard map (fibre 6) and the
    // A11 projection Hom^1([E(1),F(1),E(1),F(1)], L_q) over (2,2) (source fibre 36)
    let mut lib: TmplLib<F> = TmplLib::new(None);
    let (_, ids) = lib.rickard(2, 2, true);
    let (t, fib) = ids.iter().map(|&t| (t, lib.get(t).src_flows.len().max(lib.get(t).tgt_flows.len()))).max_by_key(|x| x.1).unwrap();
    let x36 = Window::new(&[2, 2], &[e(0, 1), f(0, 1), e(0, 1), f(0, 1)]);
    let lq = Window::new(&[2, 2], &[e(0, 1), f(0, 1)]);
    let t36 = lib.hom_basis(&x36, &lq, 1)[0];
    assert_eq!(lib.get(t36).src_flows.len(), 36);
    let y = points::<F>(4, 37, 1, 0).remove(0).x;
    let mut t_loc = 0.0f64;
    for (tt, fb) in [(t, fib), (t36, 36)] {
        let t2 = Instant::now();
        let n = 200;
        for _ in 0..n {
            let _ = lib.localize(tt, &y);
        }
        let dt = t2.elapsed().as_secs_f64() / n as f64;
        eprintln!("[perf] template localization at window fibre {}: {:.3} ms (target < 1 ms)", fb, dt * 1e3);
        t_loc = t_loc.max(dt);
    }
    assert!(t_fi < 2e-3 && fi.bytes() <= 26_000, "FlowIndex target missed");
    assert!(t_all < 20e-3, "Eval::all target missed");
    assert!(t_loc < 1e-3, "localization target missed");
}

#[test]
#[ignore = "§3.6: FlowIndex of the (2,2,2,2) X with 51 840 flows (release)"]
fn perf_flowindex_2222() {
    let lam = [2u8, 2, 2, 2];
    let t0 = Instant::now();
    let (_, xs) = canonical_and_x(&lam);
    let big = xs.iter().max_by_key(|x| x.nflows()).unwrap().clone();
    eprintln!("[perf] (2,2,2,2): {} X ladders from the canonical bases in {:.1}s; largest {} flows, {} rungs", xs.len(), t0.elapsed().as_secs_f64(), big.nflows(), big.rungs.len());
    let t1 = Instant::now();
    let fi = FlowIndex::build(&big);
    let dt = t1.elapsed().as_secs_f64();
    let mut s = Schur::new(&lam);
    deodhar_check(&fi, &mut s).unwrap();
    eprintln!("[perf] FlowIndex (2,2,2,2) X {} flows x {} rungs: {:.3} s, {} B (support {} cosets; targets < 0.3 s, <= 4.6 MB + support)", fi.nflows(), fi.nr, dt, fi.bytes(), fi.support().len());
    assert_eq!(big.nflows(), 51_840);
    assert!(dt < 0.3);
    assert!(fi.bytes() <= 4_600_000 + 12 * fi.support().len() + 1024);
}
