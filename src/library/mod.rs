//! The decomposition library of the 3-4-strand route: the Hom provider and certification
//! (`homprov`), singular indecomposables and decompositions with Foam_N pruning (`sindec`), the
//! layer-A text format (`smapio`), and the light-leaves providers (`leaves`, `leafgen`, `rexplan`).

/// Hom provider, stacked-point certification, cost tables
pub mod homprov;
/// singular indecomposables, merge decompositions, layer-A cache
pub mod sindec;
/// text export/import of layer-A sub-DAGs
pub mod smapio;
/// light-leaves provider (stage 0, Leaves, LeafPairs)
pub mod leaves;
/// closed-form generators in localized form
pub mod leafgen;
/// singlestep forms and rex-chain planner
pub mod rexplan;
