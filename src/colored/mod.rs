//! The colored engines: the exact 2-strand singular engine (`xengine`) with the colored closure
//! (`cclose`); the localized multi-strand engine (`sengine`, `swire`) on point-independent map DAGs
//! (`smap`, `sexpr`) and exact window maps (`tmpl`); the fixed-N box closure (`boxclose`) and the
//! checkpointed, sharded box driver (`boxrun`).
//!
//! The items of the cabling reference method ([`crate::reference::cabling`], the former top-level
//! module `colored`) are re-exported here, so old paths such as `krhom::colored::khr_colored_cable`
//! keep working.

/// exact small (2-strand) colored engine
pub mod xengine;
/// colored closure layer
pub mod cclose;
/// localized multi-strand scanning engine (core)
pub mod sengine;
/// sengine wired to the library and the box closure
pub mod swire;
/// global maps between ladder expressions (DAG)
pub mod smap;
/// ladder expressions with zero labels and their flows
pub mod sexpr;
/// exact window maps, localized per window root tuple
pub mod tmpl;
/// closure at a fixed N (box closure)
pub mod boxclose;
/// checkpointed, resumable, sharded box driver
pub mod boxrun;

// Old path `krhom::colored::*` (the former top-level module `colored`, now `reference::cabling`).
pub use crate::reference::cabling::*;
