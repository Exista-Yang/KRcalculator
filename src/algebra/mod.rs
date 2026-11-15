//! Coefficient fields and linear algebra: prime fields and Q, Laurent and multivariate
//! polynomials, polynomial matrices, sparse and dense matrices, sparse linear systems, and the
//! universal structure constants of relative Grassmannians (`symfun`).

/// prime fields F_p and Q
pub mod field;
/// small i64 polynomials (structure-constant precomputation)
pub mod poly;
/// multivariate polynomials with packed monomials
pub mod mpoly;
/// Laurent polynomials in q and Poincaré polynomials in (t, q)
pub mod laurent;
/// sparse matrices with polynomial entries
pub mod pmat;
/// sparse vectors and matrices
pub mod sparse;
/// dense row reduction, rank, kernel
pub mod dense;
/// sparse linear systems: kernel and rank
pub mod linsys;
/// relative Grassmannians Gr(r, b): structure constants, GrTable cache
pub mod symfun;

// Old path `krhom::algebra::*` (the former top-level module `algebra`, now `legacy::algebra`).
pub use crate::legacy::algebra::*;
