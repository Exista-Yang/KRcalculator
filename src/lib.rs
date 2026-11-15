//! krhom: colored Khovanov–Rozansky gl_N link homology.
//!
//! # Module map
//!
//! | module | contents |
//! |---|---|
//! | [`calc`] | the calculator behind `krhom`: argument model, route selection, library registry, output |
//! | [`input`] | `crossing` (braid generators, colored Rickard complexes), `pd` (PD codes → braids) |
//! | [`algebra`] | `field`, `poly`, `mpoly`, `laurent`, `pmat`, `sparse`, `dense`, `linsys`, `symfun` |
//! | [`combinat`] | `web`, `moy`, `hecke`, `schur` (q-Schur oracle), `coset`, `dcoset` |
//! | [`soergel`] | `sbim` (singular Bott–Samelson bimodules), `local`, `soergel`, `indec` |
//! | [`uncolored`] | `engine`, `engine2` (uncolored route), `tower` |
//! | [`colored`] | `xengine`, `cclose` (exact 2-strand route); `sengine`, `swire`, `smap`, `sexpr`, `tmpl`, `boxclose`, `boxrun` (3-4-strand box route) |
//! | [`library`] | `homprov`, `sindec`, `smapio`, `leaves`, `leafgen`, `rexplan` (decomposition library) |
//! | [`reference`] | `cube` (full cube), `cabling` (reference methods) |
//! | [`legacy`] | `algebra`, `hom`, `untrunc` (obsolete models, not used by the calculator) |
//!
//! Main dependency direction: `algebra` → `combinat` → `soergel` → `uncolored` / `colored` ⇄
//! `library` → `calc`; `input` and `reference` sit beside them.
//!
//! # Old paths
//!
//! Until the restructure every module was top-level (`krhom::sbim`, `krhom::xengine`, …). Those
//! paths keep working through the re-exports below. `krhom::algebra`, `krhom::soergel` and
//! `krhom::colored` are now the group modules; the items of the former modules of these names are
//! re-exported inside them (from `legacy::algebra`, `soergel::soergel` and `reference::cabling`).
//! New code should use the grouped paths.
//!
//! Stamps and cache keys (`sindec`, `tmpl`, `homprov`, `leaves`, `smapio`, `boxclose`, `boxrun`,
//! `sengine`, `swire`) are literal strings and do not depend on module paths.

pub mod algebra;
pub mod calc;
pub mod colored;
pub mod combinat;
pub mod input;
pub mod legacy;
pub mod library;
pub mod reference;
pub mod soergel;
pub mod uncolored;

// Old top-level paths (before the restructure), kept for external callers.
pub use algebra::{field, poly, mpoly, laurent, pmat, sparse, dense, linsys, symfun};
pub use combinat::{hecke, schur, coset, dcoset, web, moy};
pub use soergel::{sbim, indec, local};
pub use legacy::{hom, untrunc};
pub use uncolored::{engine, engine2, tower};
pub use colored::{xengine, cclose, sengine, swire, smap, sexpr, tmpl, boxclose, boxrun};
pub use library::{homprov, sindec, smapio, leaves, leafgen, rexplan};
pub use input::{pd, crossing};
pub use reference::cube;

/// Number of worker threads: the environment variable KRHOM_THREADS if set
/// (useful to keep a shared machine responsive), otherwise all available cores.
pub fn num_threads() -> usize {
    std::env::var("KRHOM_THREADS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|&k| k > 0)
        .unwrap_or_else(|| std::thread::available_parallelism().map(|x| x.get()).unwrap_or(4))
}
