//! Singular Soergel bimodules: singular Bott–Samelson bimodules of ladder webs, equivariant
//! localization, the localized Soergel calculus for S_n and indecomposables B_x.
//!
//! The items of [`soergel`](self::soergel) are also re-exported here, so the old path
//! `krhom::soergel::Ctx` keeps working.

/// singular Bott–Samelson bimodules (untruncated), Hom spaces
pub mod sbim;
#[allow(clippy::module_inception)]
/// localized Soergel calculus for uncolored braids
pub mod soergel;
/// indecomposables B_x: idempotents and decompositions
pub mod indec;
/// equivariant localization at torus fixed points
pub mod local;

// Old path `krhom::soergel::*` (the former top-level module `soergel`).
pub use self::soergel::*;
