//! Obsolete models kept only for the experiments documented in docs/dev/PLAN.md: truncated
//! edge algebras and their bimodules (`algebra`, `hom`), and the cut-off untruncated bimodules
//! (`untrunc`, superseded by [`crate::soergel::sbim`]). Not used by the calculator.
//!
//! The old path `krhom::algebra::Alg` (of the former top-level module `algebra`) is kept by a
//! re-export in [`crate::algebra`].

/// truncated edge algebras (obsolete)
pub mod algebra;
/// Hom spaces of edge-algebra bimodules (obsolete)
pub mod hom;
/// cut-off untruncated singular Bott–Samelson bimodules (superseded)
pub mod untrunc;
