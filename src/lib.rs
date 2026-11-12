//! krhom: colored Khovanov–Rozansky gl_N link homology.

pub mod dense;
pub mod field;
pub mod laurent;
pub mod poly;
pub mod sparse;
pub mod symfun;
pub mod moy;
pub mod web;
pub mod algebra;
pub mod hom;
pub mod untrunc;
pub mod linsys;
pub mod mpoly;
pub mod pmat;
pub mod sbim;
pub mod crossing;
pub mod cube;
pub mod hecke;
pub mod local;
pub mod soergel;
pub mod indec;
pub mod engine;
pub mod pd;
pub mod colored;
pub mod engine2;
pub mod tower;
pub mod schur;
pub mod cclose;
pub mod xengine;
pub mod coset;
pub mod sexpr;
pub mod smap;
pub mod tmpl;
pub mod homprov;
pub mod smapio;
pub mod sindec;
pub mod dcoset;
pub mod leaves;
pub mod leafgen;

/// Number of worker threads: the environment variable KRHOM_THREADS if set
/// (useful to keep a shared machine responsive), otherwise all available cores.
pub fn num_threads() -> usize {
    std::env::var("KRHOM_THREADS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .filter(|&k| k > 0)
        .unwrap_or_else(|| std::thread::available_parallelism().map(|x| x.get()).unwrap_or(4))
}
