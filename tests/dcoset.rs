//! M4c A1, gate G-L1 (combinatorics of src/combinat/dcoset.rs): singular Deodhar sums of the independent
//! rung-step generating functions against `Schur::ladder_char` on every canonical ladder of the
//! blocks (λ, μ) (every μ) of (3,3), (2,1,2), (2,2,2), (3,3,3) and on every one-rung extension; flows per
//! coset (`sexpr::deodhar_check`); Lemma 8.7 sums (identity paths = c_diag, extensions =
//! Σ term_mult·b_q); minimal representatives, ψ and Bruhat order; rung-level Grassmannian
//! reduction (G_local = G_global); degree sign (E2 Eta degree = identity-path defect).
use krhom::combinat::dcoset::*;
use krhom::combinat::schur::{tables, Elt, Lp, Schur, Table, Tid};
use krhom::library::sindec::canonical_ladder;
use krhom::combinat::web::{Dir, Rung};
use std::collections::BTreeMap;

const LAMS: [&[u8]; 4] = [&[3, 3], &[2, 1, 2], &[2, 2, 2], &[3, 3, 3]];

/// every block μ of λ: compositions of K with as many parts as λ (zeros allowed)
fn visited_blocks(lam: &[u8]) -> Vec<Vec<u8>> {
    let k: u8 = lam.iter().sum();
    let n = lam.len();
    let mut out = vec![];
    fn rec(n: usize, left: u8, cur: &mut Vec<u8>, out: &mut Vec<Vec<u8>>) {
        if cur.len() == n - 1 {
            cur.push(left);
            out.push(cur.clone());
            cur.pop();
            return;
        }
        for a in 0..=left {
            cur.push(a);
            rec(n, left - a, cur, out);
            cur.pop();
        }
    }
    rec(n, k, &mut vec![], &mut out);
    // λ first (closure types), then the rest
    out.sort_by_key(|m| m.as_slice() != lam);
    out
}

fn as_map(e: &Elt) -> BTreeMap<Tid, Lp> {
    e.t.iter().cloned().collect()
}

fn rungs_on(top: &[u8]) -> Vec<Rung> {
    let mut out = vec![];
    for i in 0..top.len() - 1 {
        for dir in [Dir::E, Dir::F] {
            let giver = if dir == Dir::E { top[i + 1] } else { top[i] } as usize;
            for r in 1..=giver {
                out.push(Rung { i, dir, r });
            }
        }
    }
    out
}

/// every (block, type, ladder) of λ: minimal-flow and BFS canonical ladders
fn ladders(s: &mut Schur) -> Vec<(Vec<u8>, Tid, Vec<Rung>)> {
    let lam = s.lam.clone();
    let mut out = vec![];
    for mu in visited_blocks(&lam) {
        let cb = s.canonical_basis(&mu);
        let mut ord = cb.order.clone();
        ord.sort();
        for t in ord {
            let cb = s.canonical_basis(&mu);
            let a = cb.ladder_minflow[&t].clone();
            let b = cb.ladder_bfs[&t].clone();
            out.push((mu.clone(), t, a.clone()));
            if b != a {
                out.push((mu.clone(), t, b));
            }
        }
    }
    out
}

#[test]
fn gl1_deodhar_sums_equal_ladder_char() {
    let mut total = 0usize;
    for lam in LAMS {
        let mut s = Schur::new(lam);
        let mut sc = StepCache::new();
        let lads = ladders(&mut s);
        let mut n = 0usize;
        for (_, _, w) in &lads {
            let mut words = vec![w.clone()];
            for r in rungs_on(&top_labels(lam, w)) {
                let mut w2 = w.clone();
                w2.push(r);
                words.push(w2);
            }
            for x in words {
                let fw = forward(&mut s, &mut sc, &x);
                let ch = s.ladder_char(&x);
                assert_eq!(fw.last().unwrap(), &as_map(&ch), "G-L1: Deodhar sum != ladder_char for {:?} {:?}", lam, x);
                n += 1;
            }
        }
        eprintln!("G-L1 Deodhar: λ = {:?}: {} canonical ladders, {} words (with one-rung extensions) agree", lam, lads.len(), n);
        total += n;
    }
    assert!(total > 1000);
}

#[test]
fn gl1_flows_per_coset() {
    // flows of the FlowIndex at every coset = c_t(char)(1) (sexpr::deodhar_check), on the
    // canonical ladders up to 20 000 flows
    for lam in LAMS {
        let mut s = Schur::new(lam);
        let mut n = 0;
        for (_, _, w) in ladders(&mut s) {
            if krhom::combinat::schur::flows(lam, &w) > 20_000 {
                continue;
            }
            let e = krhom::colored::sexpr::Expr::new(lam, &w).unwrap();
            let fi = krhom::colored::sexpr::FlowIndex::build(&e);
            krhom::colored::sexpr::deodhar_check(&fi, &mut s).unwrap();
            n += 1;
        }
        eprintln!("G-L1 flows: λ = {:?}: {} ladders checked", lam, n);
    }
}

#[test]
fn gl1_lemma87_sums() {
    for lam in LAMS {
        let mut s = Schur::new(lam);
        let mut sc = StepCache::new();
        let id = s.intern(Table::identity(lam));
        // identity paths of indecomposable L_p = c_diag(b_p); the junction products are constant
        let mut nd = 0;
        let ord = s.canonical_basis(lam).order.clone();
        for &p in &ord {
            let w = canonical_ladder(&mut s, p);
            if !krhom::library::sindec::is_indecomposable_of(&mut s, &w, p) {
                continue;
            }
            let js = junction_types(&mut s, &mut sc, &w, &[id]);
            let cd = s.c_diag(p);
            for lv in &js {
                let mut acc = Lp::zero();
                for j in lv {
                    acc = acc.add(&j.fwd.mul(&j.back));
                }
                assert_eq!(acc, cd, "Lemma 8.7: Σ_t fwd·back at level {} != c_diag for {}", lv.first().map_or(0, |j| j.level), s.table(p).fmt());
            }
            nd += 1;
        }
        // one-rung extensions of every canonical ladder: Σ_q term_mult_q · b_q
        let mut ne = 0;
        for mu in visited_blocks(lam) {
            let ordm = s.canonical_basis(&mu).order.clone();
            for &p in &ordm {
                let w = canonical_ladder(&mut s, p);
                if !krhom::library::sindec::is_indecomposable_of(&mut s, &w, p) {
                    continue;
                }
                for r in rungs_on(&mu) {
                    let mut w2 = w.clone();
                    w2.push(r);
                    let fw = forward(&mut s, &mut sc, &w2);
                    let tm = s.term_mult(p, &[r]);
                    let mu2 = top_labels(lam, &w2);
                    let mut acc = Elt::zero(&mu2);
                    for (q, m) in &tm {
                        let bq = s.canonical_basis(&mu2).b[q].clone();
                        for (e, &k) in m.c.iter().enumerate() {
                            if k != 0 {
                                acc = acc.add_scaled(&bq, k, m.lo + e as i32);
                            }
                        }
                    }
                    assert_eq!(fw.last().unwrap(), &as_map(&acc), "Lemma 8.7: extension sums != Σ term_mult·b_q for {} · {:?}", s.table(p).fmt(), r);
                    ne += 1;
                }
            }
        }
        eprintln!("G-L1 Lemma 8.7: λ = {:?}: {} closure types (identity paths = c_diag), {} extensions (= Σ term_mult·b_q)", lam, nd, ne);
    }
}

#[test]
fn gl1_minperm_psi_bruhat() {
    for lam in LAMS {
        let mut s = Schur::new(lam);
        let mut nt = 0;
        let mut bad_reduced = vec![];
        for mu in visited_blocks(lam) {
            let ts = tables(lam, &mu);
            let id_like = ts.iter().min_by_key(|t| t.lmax()).unwrap().clone();
            for t in &ts {
                let w = min_perm(t);
                assert_eq!(perm_len(&w), t.lmin(), "ℓ(p_) != lmin for {}", t.fmt());
                assert_eq!(table_of_perm(&w, lam, &mu), *t, "table_of_perm(min_perm) for {}", t.fmt());
                assert_eq!(psi_weighted_len(t), t.lmin(), "weighted ℓ(ψ) != lmin for {}", t.fmt());
                assert!(bruhat_le(&id_like, t), "the minimal table is not below {}", t.fmt());
                // lmax(p) < lmax(q) is necessary for p < q
                nt += 1;
            }
            for t in &ts {
                for u in &ts {
                    if u != t && bruhat_le(u, t) {
                        assert!(u.lmax() < t.lmax() && u.lmin() < t.lmin(), "Bruhat {} < {} but lengths do not decrease", u.fmt(), t.fmt());
                    }
                }
            }
            // every canonical ladder is a reduced word of its type (support Bruhat-below)
            let ord = s.canonical_basis(&mu).order.clone();
            for &p in &ord {
                let w = canonical_ladder(&mut s, p);
                let r = reduced_for(&mut s, &w, p);
                if !r.ok() {
                    bad_reduced.push((s.table(p).fmt(), r));
                }
            }
        }
        eprintln!("G-L1 tables: λ = {:?}: {} tables; canonical ladders not reduced: {:?}", lam, nt, bad_reduced);
        assert!(bad_reduced.is_empty(), "canonical ladders that are not reduced words: {:?}", bad_reduced);
    }
}

#[test]
fn gl1_grassmannian_reduction() {
    for lam in LAMS {
        let mut n = 0;
        for mu in visited_blocks(lam) {
            for t in tables(lam, &mu) {
                for r in rungs_on(&mu) {
                    for st in step_terms(&t, &r) {
                        let g = grass(&t, &st.to, &r);
                        let qm = min_perm(&g.q);
                        assert_eq!(min_perm(&g.z), qm, "z_ != Q_");
                        for (m, tt) in [(&g.m_from, &t), (&g.m_to, &st.to)] {
                            let mp = min_perm(m);
                            assert_eq!(perm_comp(&qm, &mp), min_perm(tt), "t_ != Q_ m_ ({} via {:?})", tt.fmt(), r);
                            assert_eq!(perm_len(&qm) + perm_len(&mp), tt.lmin(), "lengths do not add");
                        }
                        // the local step m' -> m'' of the same rung has the same G (Lemma 3.9)
                        let loc = step_gf_table(&g.m_from, &r);
                        let gl = loc.iter().find(|x| x.0 == g.m_to).map(|x| x.1.clone()).unwrap_or_else(Lp::zero);
                        assert_eq!(gl, st.gf, "G_local != G_global for {} -> {} via {:?}", t.fmt(), st.to.fmt(), r);
                        n += 1;
                    }
                }
            }
        }
        eprintln!("G-L1 Grassmannian reduction: λ = {:?}: {} steps, G_local = G_global", lam, n);
    }
}

#[test]
fn gl1_degree_sign() {
    use krhom::algebra::field::{Field, F31};
    use krhom::colored::sexpr::Expr;
    use krhom::colored::smap::{Arena, ElemKind, MapNode};
    use krhom::colored::tmpl::TmplLib;
    // the E2 unit Eta: R -> [E0.b F0.b] over (a, b) has krhom degree +ab; the identity path of
    // E0.b F0.b has defect +ab: EKLP defect = krhom normalized degree (no sign flip)
    for (a, b) in [(1u8, 1u8), (2, 1), (1, 2), (2, 2), (3, 3), (3, 6)] {
        let lab = [a, b];
        let w = vec![Rung { i: 0, dir: Dir::E, r: b as usize }, Rung { i: 0, dir: Dir::F, r: b as usize }];
        let mut ar = Arena::<F31>::new(&lab, TmplLib::new(None));
        let wi = ar.expr(Expr::new(&lab, &w).unwrap());
        let ei = ar.expr(Expr::new(&lab, &[]).unwrap());
        let eta = ar.add(MapNode::Elem { kind: ElemKind::Eta, src: ei, tgt: wi, at: 0, upright: 0, scalar: F31::one() }, "").unwrap();
        let mut s = Schur::new(&lab);
        let mut sc = StepCache::new();
        let id = s.intern(Table::identity(&lab));
        let fw = forward(&mut s, &mut sc, &w);
        let c = fw.last().unwrap().get(&id).cloned().unwrap_or_else(Lp::zero);
        assert_eq!(c, Lp::mono((a as i32) * (b as i32), 1), "identity path of E0.{} F0.{} on ({},{})", b, b, a, b);
        assert_eq!(ar.info(eta).deg, (a as i32) * (b as i32), "Eta degree");
    }
}
