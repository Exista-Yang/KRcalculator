//! The uncolored engine: Bar-Natan style scan in the localized Soergel calculus (`engine`),
//! the scalable gl_N closure for several N (`engine2`, used by krhom for uncolored braids) and the
//! restriction-tower closure (`tower`).

/// scanning engine for uncolored braid closures (milestone 6F)
pub mod engine;
/// one scan, closure for several N (krhom's uncolored route)
pub mod engine2;
/// restriction-tower closure
pub mod tower;
