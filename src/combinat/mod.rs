//! Combinatorics of webs, Hecke algebras and cosets: ladder webs, the MOY oracle, the Hecke
//! algebra with its Kazhdan–Lusztig basis, the q-Schur algebroid (decategorified oracle), and cosets
//! and double cosets.

/// Hecke algebra of S_n, Kazhdan–Lusztig basis
pub mod hecke;
/// q-Schur algebroid: Hom formula, decat_chi, kill rule
pub mod schur;
/// cosets W/W_μ as top configurations of roots
pub mod coset;
/// double-coset combinatorics for the light-leaves backend
pub mod dcoset;
/// ladder webs
pub mod web;
/// MOY evaluation of closed ladder webs (independent oracle)
pub mod moy;
