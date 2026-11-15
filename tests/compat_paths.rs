//! The pre-restructure module paths (every module top-level: `krhom::sbim`, `krhom::colored`, …)
//! keep resolving through the re-exports in src/lib.rs and the group modules.
#![allow(unused_imports)]

use krhom::algebra::Alg; // former top-level `algebra` (now legacy::algebra)
use krhom::colored::{cable, khr_colored_cable, Cabled}; // former top-level `colored` (now reference::cabling)
use krhom::soergel::{Ctx, Op}; // former top-level `soergel` (now soergel::soergel)
use krhom::{
    boxclose, boxrun, cclose, coset, crossing, cube, dcoset, dense, engine, engine2, field, hecke, hom, homprov, indec, laurent, leafgen,
    leaves, linsys, local, moy, mpoly, pd, pmat, poly, rexplan, sbim, schur, sengine, sexpr, sindec, smap, smapio, sparse, swire, symfun,
    tmpl, tower, untrunc, web, xengine,
};

#[test]
fn old_and_new_paths_are_the_same_items() {
    // a function reachable through both paths is one item
    let a: fn(&[usize], &[krhom::crossing::Gen]) -> krhom::reference::cabling::Cabled = krhom::colored::cable;
    let b: fn(&[usize], &[krhom::input::crossing::Gen]) -> Cabled = krhom::reference::cabling::cable;
    assert_eq!(a as usize, b as usize);
    assert_eq!(krhom::sindec::SINDEC_VERSION, krhom::library::sindec::SINDEC_VERSION);
    assert_eq!(krhom::leaves::LEAVES_STAMP, krhom::library::leaves::LEAVES_STAMP);
}

/// The stamp and key strings stored in published libraries (branch data-gl6) are literals;
/// moving modules must not change them.
#[test]
fn stamps_are_unchanged() {
    use krhom::colored::{boxclose::BOX_CODE_STAMP, boxrun::BOXRUN_STAMP, sengine::SENGINE_STAMP, swire::SWIRE_STAMP, tmpl::TMPL_ALG};
    use krhom::library::{homprov::COST_ALG, leaves, sindec};
    assert_eq!(sindec::SINDEC_VERSION, "krhom-sindec v2");
    assert_eq!(sindec::LADDER_STAMP, "schur.ladder_minflow.v1");
    assert_eq!(sindec::M5I_STAMP, "m5i-peel1");
    assert_eq!(TMPL_ALG, "hs1.lcx1");
    assert_eq!(COST_ALG, "hs-kl1");
    assert_eq!(SENGINE_STAMP, "sengine-v2");
    assert_eq!(SWIRE_STAMP, "swire-v1");
    assert_eq!(BOX_CODE_STAMP, "boxclose-v2");
    assert_eq!(BOXRUN_STAMP, "boxrun-v1");
    assert_eq!(leaves::STAGE0_STAMP, "stage0=ProductClosure(v1,inner-first)");
    assert_eq!(leaves::LEAVES_STAMP, "leaves=Leaves(v1,closure,threeup,prefix-ext+costalt5,window-whole,budget=1e5/4e5/os2e4l4)");
    assert_eq!(leaves::C27_STAMP, "c27=joint(sfx4,keep40,w12)+eqpairs+detours");
    assert_eq!(leaves::LEAF_PAIRS_STAMP, "leafpairs=LeafPairs(v1,last,rex-chain|lemmaA+exit-entry,plan=exact1e5/len12/visit6e4,assoc-as-template,persistable-only)");
}
