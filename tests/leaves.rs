//! M4c A1, Stage 0 (src/leaves.rs): the 5 product closure types of (3,3,3) at N = 6 through
//! `ProductClosure` (product rule first, inner = the production chain), certified by
//! `certify_closure` inside the provider and again here; non-product types are refused by the
//! product rule.
use krhom::field::F31;
use krhom::homprov::*;
use krhom::leaves::*;
use krhom::schur::{Schur, Table};
use krhom::sexpr::Expr;
use krhom::sindec::{canonical_ladder, provider_for};
use krhom::smap::Arena;
use krhom::tmpl::TmplLib;

const PRODUCTS: [&str; 5] = ["0,0,3/3,0,0/0,3,0", "0,3,0/0,0,3/3,0,0", "1,0,2/2,0,1/0,3,0", "1,2,0/0,0,3/2,1,0", "0,3,0/1,0,2/2,0,1"];

#[test]
fn stage0_factorizations() {
    let lam = [3u8, 3, 3];
    let mut s = Schur::new(&lam);
    for ts in PRODUCTS {
        let p = s.intern(Table::parse(ts));
        let w = canonical_ladder(&mut s, p);
        let f = product_factorization(&mut s, &w, p, 256).unwrap_or_else(|| panic!("{} is not found as a product", ts));
        assert!(f.canonical, "{}: factors are not canonical ladders", ts);
        assert!(f.swaps.len() <= 1, "{}: {} swaps", ts, f.swaps.len());
        assert_eq!(s.c_diag(f.a).mul(&s.c_diag(f.b)), s.c_diag(p));
    }
    // the 7 non-product failing types (L3 §3.2) are not products
    for ts in ["1,1,1/1,1,1/1,1,1", "1,0,2/1,2,0/1,1,1", "1,1,1/0,2,1/2,0,1", "0,2,1/2,0,1/1,1,1", "1,1,1/1,0,2/1,2,0", "0,1,2/2,0,1/1,2,0", "0,2,1/1,0,2/2,1,0"] {
        let p = s.intern(Table::parse(ts));
        let w = canonical_ladder(&mut s, p);
        assert!(product_factorization(&mut s, &w, p, 256).is_none(), "{} unexpectedly a product", ts);
    }
}

#[test]
fn stage0_product_closure_certified() {
    let lam = [3u8, 3, 3];
    let mut s = Schur::new(&lam);
    let mut ar = Arena::<F31>::new(&lam, TmplLib::new(None));
    let inner = provider_for::<F31>(&lam, Some(std::path::Path::new("tests/data/homprov/costs_333_m5c.txt")), 1.0);
    let mut pc = ProductClosure::new(inner);
    pc.product_first = true;
    let b = Budget::default();
    let r0 = ar.expr(Expr::new(&lam, &[]).unwrap());
    for ts in PRODUCTS {
        let p = s.intern(Table::parse(ts));
        let w = canonical_ladder(&mut s, p);
        let lp = ar.expr(Expr::new(&lam, &w).unwrap());
        let want = cdiag_degrees(&mut s, p);
        for outs in [false, true] {
            let (ms, _) = HomProvider::closure(&mut pc, &mut ar, &mut s, lp, p, None, outs, &b).unwrap_or_else(|e| panic!("stage0 {} outs={}: {:?}", ts, outs, e));
            certify_closure(&ar, &lam, r0, lp, None, &ms, &want, outs, 0xa1).unwrap();
        }
    }
    // the product rule refuses a non-product type
    let p = s.intern(Table::parse("1,1,1/1,1,1/1,1,1"));
    let w = canonical_ladder(&mut s, p);
    let lp = ar.expr(Expr::new(&lam, &w).unwrap());
    let e = pc.product_closure(&mut ar, &mut s, lp, p, true, &b).expect_err("non-product accepted");
    assert_eq!(e[0].kind, MissKind::Refused);
}
