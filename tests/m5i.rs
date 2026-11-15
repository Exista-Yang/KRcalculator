//! C2-4 (M5i): idempotents of decomposable canonical ladders and their use in sindec / swire.
//!
//! Default suite: (3,1,3) (the decomposable type 2,0,1/0,1,0/1,0,2, 48 flows): e_p built by
//! peeling the extra summand (provider pairs), checked against the exact End^0(L_p) (dim 2, e_p in
//! its span: the unique idempotent with top scalar 1 and the ranks of B_p), decompositions with the
//! source idempotent and the commutation reroute, closure generators through e_p, reduced
//! coordinates (QP = id, PQ = e_p) and frozen-data reduced shapes; the (2,2,2,2) type of block
//! (2,4,0,2) (96 flows) through prefix propagation. Ignored (minutes, layer A recommended): the
//! (2,2,2,2) End^0 comparison (196 k unknowns) and a forced-ladder (2,2,2) scan of 4_1.
use krhom::algebra::field::{Field, F31};
use krhom::library::homprov::{hom_unknowns, Budget, Exact, HomProvider};
use krhom::combinat::schur::{Schur, Tid};
use krhom::library::sindec::*;
use krhom::colored::smap::{points, Eval, MapId};
use krhom::combinat::web::{Dir, Rung};
use std::path::PathBuf;

type F = F31;

fn costs() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/homprov/costs_cloud.txt")
}
fn cache_dir(name: &str) -> PathBuf {
    let d = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("m5i_cache").join(name);
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}
fn builder(lam: &[u8], prune: Option<u8>, cache: Option<PathBuf>) -> SIndec<F> {
    let (prov, stamp) = provider_for_stamped::<F>(lam, Some(&costs()), 1.0);
    SIndec::with_stamp(lam, prov, stamp, SIndecOpts { cache_dir: cache, prune_n: prune, ..Default::default() })
}
fn first_decomposable(s: &mut Schur, table: Option<&str>) -> Tid {
    decomposable_ladders(s).into_iter().find(|c| table.is_none_or(|t| c.table == t)).map(|c| c.q).expect("no decomposable canonical ladder")
}
/// e_p lies in the span of an End^0(L_p) basis (stacked blocks at 2 points, all cosets)
fn in_end0_span(sd: &mut SIndec<F>, lad: u32, e: MapId) -> (usize, bool) {
    let mut ex = Exact { max_label: 9, max_fibre_product: u64::MAX, max_unknowns: u64::MAX, max_call_s: 1e9, ..Exact::default() };
    let (basis, _) = HomProvider::<F>::end0(&mut ex, &mut sd.ar, &mut sd.schur, lad, &Budget::default()).expect("End^0");
    let cos = sd.ar.flows(lad).cosets.clone();
    let pts = points::<F>(sd.ar.k(), 0xe0e0, 2, 0);
    let st = Eval::stacked(&sd.ar, &basis, &pts, &cos);
    let mut st2 = st.clone();
    st2.extend(Eval::stacked(&sd.ar, &[e], &pts, &cos));
    let n = st[0].len();
    let (r1, r2) = (krhom::algebra::dense::rank(&st, n), krhom::algebra::dense::rank(&st2, n));
    (basis.len(), r1 == basis.len() && r2 == r1)
}

#[test]
fn m5i_313_idempotent_end0_dec_closure() {
    let lam = [3, 1, 3];
    let mut sd = builder(&lam, None, Some(cache_dir("313")));
    let p = first_decomposable(&mut sd.schur, Some("2,0,1/0,1,0/1,0,2"));
    let (lad, e, rank) = {
        let t = sd.ty(p).expect("M5i e_p");
        (t.lad, t.e.expect("e_p"), t.rank)
    };
    assert_eq!(rank, 44);
    // exact reference: dim End^0 = 2 (b_p + b_r), e_p in its span
    let (dim, inspan) = in_end0_span(&mut sd, lad, e);
    assert_eq!(dim, 2);
    assert!(inspan, "e_p not in End^0(L_p)");
    // decompositions with the source idempotent (Σ ιπ = Kron(e_p, R) and Π∘Ι = δ checked inside)
    for rg in [Rung { i: 0, dir: Dir::E, r: 1 }, Rung { i: 1, dir: Dir::E, r: 1 }, Rung { i: 0, dir: Dir::F, r: 2 }] {
        let d = sd.dec_rung(p, rg).unwrap_or_else(|e| panic!("dec_rung {:?}: {}", rg, String::from(e)));
        let want: usize = sd.schur.term_mult(p, &[rg]).iter().map(|(_, m)| m.at1() as usize).sum();
        assert_eq!(d.summands.len(), want);
    }
    // closure generators through e_p: c_diag(b_p) of them, in the image of e_p
    let g = sd.closure_gens(p).expect("closure with e_p");
    assert_eq!(g.ins.len(), sd.schur.c_diag(p).at1() as usize);
    // reduced coordinates from the frozen data: QP = id, PQ = block(e_p), Σ m_c = rank B_p
    let dd = sd.freeze();
    let pt = points::<F>(dd.ar().k(), 7, 1, 0).remove(0);
    let mut ev = Eval::new(dd.ar(), pt);
    let cos = dd.ar().flows(lad).cosets.clone();
    let mut tot = 0u128;
    for &c in &cos {
        let r = dd.reduced(&mut ev, p, c).expect("reduced");
        let b = ev.block(e, c);
        assert_eq!(r.q.mul(&r.p), krhom::colored::smap::Block::identity(r.m));
        assert_eq!(r.p.mul(&r.q), *b);
        tot += r.m as u128;
    }
    assert_eq!(tot, rank);
}

#[test]
fn m5i_2402_prefix_propagation() {
    // (2,2,2,2), block (2,4,0,2): L_p = L_{p'}·F2.1 with p' decomposable (prefix propagation)
    let lam = [2, 2, 2, 2];
    let mut sd = builder(&lam, Some(4), Some(cache_dir("2402")));
    let p = first_decomposable(&mut sd.schur, Some("1,1,0,0/0,1,0,1/1,1,0,0/0,1,0,1"));
    let t = sd.ty(p).unwrap_or_else(|e| panic!("{}", String::from(e)));
    assert!(t.e.is_some());
    assert_eq!(t.rank, 86);
    assert!(sd.idem_log.iter().any(|l| l.1.starts_with("Kron(")), "{:?}", sd.idem_log);
}

#[test]
#[ignore]
fn m5i_2402_end0_reference() {
    // 196 k Exact unknowns (≈ 75 s): e_p (prefix propagation) is the unique End^0 idempotent
    let lam = [2, 2, 2, 2];
    let mut sd = builder(&lam, Some(4), Some(cache_dir("2402e")));
    let p = first_decomposable(&mut sd.schur, Some("1,1,0,0/0,1,0,1/1,1,0,0/0,1,0,1"));
    let (lad, e) = {
        let t = sd.ty(p).unwrap();
        (t.lad, t.e.unwrap())
    };
    let w = sd.ar.expr_of(lad).rungs.clone();
    assert_eq!(hom_unknowns(&lam, &w, &w, 0), 196_112);
    let (dim, inspan) = in_end0_span(&mut sd, lad, e);
    assert_eq!(dim, 2);
    assert!(inspan);
}
