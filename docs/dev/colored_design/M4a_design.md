# M4a design note (localized substrate)

Implements SPEC_M4_M6a.md §3 (revision 2). Worktree `.claude/worktrees/m4`, branch `colored-m4`
(from dev 7070814). New files: `src/coset.rs`, `src/sexpr.rs`, `src/tmpl.rs`, `src/smap.rs`,
`tests/m4a.rs`; `src/lib.rs` gets 4 `pub mod` lines. No other repository file changes
(`schur.rs`, `cclose.rs` untouched).

## Module layout

- **coset.rs**: `Roots`, `Conf`, `CosetSpace`. The rank is a mixed radix over uprights, with
  upright 0 most significant. Each digit is the colex rank of mask_j among the roots that masks
  0..j-1 do not use, with radix C(remaining, mu_j). The arithmetic is u128 and the count is
  asserted to fit u64. `rep` follows the spec rule (column by column, the first unused roots of
  each bottom block). `cosets_of` enumerates the cosets by subset choice and sorts them.
- **sexpr.rs**: `Expr` (validated; zero labels allowed), `Flow`, and the packed `FlowIndex`:
  - `build` runs a forward DFS (as `local::fixed_points`), ranks each flow's top configuration,
    sorts the flows by (rank, masks), and writes `cosets`/`start`/`masks` with exact capacity;
  - `find` binary-searches the coset, then the masks inside it;
  - additive helpers: `locate(moved)` (global id from moved masks, by replay), `nflows`,
    `top_conf`.

  `deodhar_check` uses `Schur::ladder_char` and `CosetSpace::cosets_of`.
- **tmpl.rs**: `Window`, `Template` (all fields pub, as in the spec), `WinLoc` (CSC), and
  `TmplLib`. The lib keeps per-template auxiliary data in a private parallel vector:
  - the window SBims, fixed points, and the fixed-point-to-sorted permutation;
  - the window FlowIndexes.

  `localize` = `local::loc_matrix`/`inverse`/`eval_pmat` on the stored window SBims, reordered
  into CSC over the sorted window flows. The lib has no interior mutability (Send + Sync).

  The disk cache is one text file per `hom_basis`/`rickard` call, named
  `tmpl-<fnv64(key)>.txt`. The header is
  `krhom-tmpl v1 field=<F::name()> kind=<hom|rickard> key=<canonical key>`. Coefficients are
  written with Display and read back with `from_i64`. A bad header is counted in
  `stats().rejected` and the maps are recomputed.

  `tensor_exact` substitutes the prefix's top operators for the window variables (the
  w-upright generalization of `cube::apply_local_poly`), then applies `kron_id`.

  `poly_certificate` follows §3.4.5:
  - one point set for all entries;
  - one Vandermonde inverse per distinct degree D;
  - a redraw of the whole set with seed+1 when a Vandermonde is singular, at most 3 times;
  - the residual is checked on all points that the solve did not use.
- **smap.rs**: `Point`/`points`, `Block`, `gemm_acc`, `Layout`, `LocMap`, `ElemKind`, `MapNode`,
  `MapInfo`, `Arena`, `ArenaMark`, `ArenaBytes`, `Eval`, `Reduced`/`reduced_of`.
  - The Arena keeps:
    - interned exprs, with `Option<FlowIndex>` per expr and a generation stamp per built index;
    - nodes, infos, and a per-node window record (`at`, `first`, `w`, `|A|`, `|A'|`, window ids);
    - a window table (window `Expr` to window FlowIndex, never dropped, like templates);
    - the `TmplLib` and the elem constants.
  - `Local` and `Elem` are evaluated by one window routine (§3.4.2): it remaps the source flow's
    window masks to window-local roots, locates the window column, walks the CSC column, and
    maps each target window flow back into the global target flow (`find`, asserted).
  - Closed forms produce the same CSC from the window flows at window roots y (§3.4.3, every
    constant 1). Eval memoizes:
    - blocks per (MapId, coset);
    - window maps per (template or elem key, window masks);
    - layouts per (src, tgt).
  - `commit` keeps the new nodes reachable from `keep`, renumbered in creation order (children
    precede parents). It keeps their expressions and indices, also renumbered. It drops indices
    built after the mark for old expressions that no surviving node references. Template and
    window tables are kept.

## Decisions where the spec is silent or self-contradictory (also in the report)

1. **`hom_basis` signature.** `hom_basis` returns `Vec<TmplId>` (binding), but §3.4.4 says that
   windows above fsize 1 600 return Err. Added `try_hom_basis -> Result<Vec<TmplId>, String>`;
   `hom_basis` panics on that Err.
2. **Assoc window.** Assoc (E4) windows are not 2-upright (A8 uses 3-upright windows).
   - The window is everything after the common prefix `at`, minus the longest common rung
     suffix of src and tgt. It spans uprights `upright..=max touched upright`.
   - Its labels key for `certify_elem` is (K_w, w): the window total label and the width.
   - The bijection with top configurations and the equality of supports are asserted in `add`.
3. **Digon labels key.** For Iota/Trace/IotaL/TraceL the key is (c, b) (the digon (c; b)). For
   E2 it is (a, b).
4. **Bottom labels.** `Arena::expr` requires `bottom == lam`, since points have length K and a
   map needs equal bottoms. Tests that need other labels build separate arenas.
   `Arena::into_lib` (additive) moves the template lib between them.
5. **Additive helpers.** `Arena::assert_genuine(m) -> Result<(), String>` is the M4a form of
   M4b's `assert_genuine`: it walks the DAG and names the first unrecorded (kind, labels).
   Other additive helpers: `Arena::lib_ref`, `Eval::winloc(m)` (whole map as CSC over sorted
   global flows; used for certificates of composites), `TmplLib::add_template`,
   `TmplLib::stats`, `Template::new`, and `tmpl::read_cache_file`.
6. **Point seeding.** Point id `i` uses
   `StdRng::seed_from_u64(seed ^ (i * 0x9E3779B97F4A7C15))`. Draws are uniform in 1..2e9 and
   redrawn on a collision, so a point depends only on (seed, id).
7. **A10(c).** "Eps with a sign flipped on one non-identity flow" is a no-op, because Eps
   vanishes off the identity flow. The test applies the flip to Iota and to Trace, whose values
   are nonzero on every flow. It also keeps the literal Eps case and records that it passes.
8. **A15 "non-constant use".** Lin coefficients are field constants by type. The test covers an
   empty Lin, mixed degrees, and mixed ends.
9. **A10 and A13 placement.** A10's "passes for every template of A4, A5, A7, A8, A9" and A13's
   "fresh lib loads identical PMats for every template of A4-A9" run inside those tests, on the
   templates they build: `cert_template` and `roundtrip`. This avoids recomputing the templates.
   The `a10_certificates` and `a13_cache_header_and_corruption` tests hold the remaining parts:
   - A10: closed forms with dots, the failure cases (a)-(c), and the T / |Mono_D| report;
   - A13: the header check, corruption rejection, and field separation.
10. **Which tests are ignored, from measured costs** (2026-11-10):
    - A11 on (2,2) takes 0.12 s, so it is in the default suite. The spec makes it ignored only
      above 30 s.
    - A9 (3,1), (1,3), (3,2), (2,3): each `local_cx` call takes < 0.08 s, so they are default.
    - A9 (3,3): `local_cx` takes 108 s per sign, so it is ignored. It passed: 215 s, 38 MB.
    - A7(iv) is split into E2 and E3 halves, which print per-case progress. Both exceed 10 min on
      the laptop, so they are cloud jobs.
