# M4b design note (Hom provider, certification, gate G1)

Implements SPEC_M4_M6a.md §4 (revision 2) on top of M4a. Worktree `.claude/worktrees/m4b`,
branch `colored-m4b` (from dev c77434c). New files: `src/homprov.rs`, `src/bin/g1.rs`,
`tests/homprov.rs`; `src/lib.rs` gets `pub mod homprov;`. M4a review fixes touch `src/tmpl.rs`,
`src/smap.rs` and `tests/m4a.rs` only.

## 0. M4a review fixes (done first; each has a test in tests/m4a.rs, `review_*`)

| finding | fix |
|---|---|
| (a) certificate accepts non-bimodule maps | `poly_certificate` checks block diagonality over the window top cosets at every point (Err names the entry) |
| (b) `add_template` yields certified nodes | new `TmplSource::Added { certified }`; `add_template` forces `certified: false`; `TmplLib::certify_added` = homogeneous degree + exact commutation with every top operator + `poly_certificate`; `Arena::certify_template` wraps it; `Local` nodes and `assert_genuine` consult `TmplLib::is_certified` |
| (c) Rickard uniqueness not asserted | `rickard` asserts dim Hom^1 = 1, Hom^{<1} = 0 (`schur::hom_gd` on term characters) and d_j != 0, for computed and cached maps |
| (d) stale `certified` flag | `Arena::refresh_certified` recomputes every node after `certify_elem`, `certify_assoc`, `certify_template` (children precede parents) |
| (e) Assoc keys | `add` asserts equal window shifts; Assoc constants are recorded per window pair (`certify_assoc`, `assoc_certified`, `elem_windows`); `certify_elem(Assoc, ..)` panics |
| (f) cache | header `krhom-tmpl v2 alg=<TMPL_ALG> field=..`; other stamps are rejected; tmp names carry pid and a process counter |
| (g) point seeds | `StdRng::from_seed(point_seed(seed, id))` with (seed, id) in disjoint bytes (injective) |

## 1. Measured facts that shape the provider (this milestone)

- `g1_targets` reproduces the spec counts: (3,3) 12 non-identity pairs, 72 targets, Σ mult 121,
  4 closure types; (2,2,2) 98 / 358 / 404 / 21.
- Whole-ladder Exact is far out of reach for (2,2,2) pairs (fibre products up to 52 000).
- The maximal non-common window of L_q and X = L_p R (common prefix and suffix, restricted to
  the uprights the differing rungs touch) is <= 1 600 for 280 of 358 (2,2,2) targets; after
  canonicalizing both words by rung commutation (degree-0 isomorphisms, below) for 304.
- hom_space time is not a function of fibre or unknowns alone: a 3-upright 4-rung window,
  36 x 15 flows, 2 580 unknowns took 11.7 s, while 18 x 15 / 872 unknowns took 0.33 s.

## 2. Search: local decomposition tree (tiers S1-S3) plus exact windows (S6)

**Rung commutation.** Adjacent rungs on disjoint upright pairs, or on pairs sharing one upright
when both merge into it or both split from it, commute up to a degree-0 isomorphism (associativity
of a three-part merge/split). Realized as an exact window template (3 or 4 uprights, 2 rungs,
`hom_basis` in degree 0, dimension 1 asserted).

**Local decomposition tree (the S3 "paths", moves M4/M5 made systematic).** For a word Z:
1. if `expand(ladder_char(Z))` is a single b_t at shift 0 with coefficient 1, Z is a leaf;
2. otherwise pick, over the commutation class of Z, the decomposable window W (2 or 3 adjacent
   uprights, 2-4 contiguous rungs) with the smallest predicted exact cost; decompose W in its
   own block (Hecke algebroid on the window labels: W = ⊕ L_t{s}^c with L_t the window block's
   minimal-flow canonical ladders); exact ι, π on the window (`hom_basis`), intersection form,
   rref selection and Neumann orthogonalization exactly as `xengine::decompose_direct`, but on
   the window; each local summand gives a child word with the window replaced by L_t;
3. recurse (only into children whose character still contains the requested (q, s)).

The composite ι = (commutation isos) ∘ (id ⊗ ι_loc ⊗ id) ∘ ... and π likewise are Arena DAGs of
certified leaves (window templates from hom_space). Along one root-to-leaf path π∘ι = id; across
different leaves π_a ∘ ι_b = 0 at their first divergence. A leaf word Z with b_t = b_q but
Z != L_q is connected to L_q by degree-0 maps φ: Z -> L_q, ψ: L_q -> Z found as exact windows on
the non-common window of the best commutation-equivalent pair (any nonzero degree-0 map between
indecomposables of the same type is an isomorphism, C2 M6); the form entry becomes the top scalar
of φ∘ψ (nonzero, checked).

**Tiers as implemented.**
- S0: identity (equal words, s = 0, mult 1; `Id` for full/end0 when dim Hom^d = 1 and d = 0).
- S1: extremes: tree leaves of (q, s_lo) (ι) and (q, s_hi) (π); certified by rule 5 when the
  space is 1-dimensional. Closure: the lowest in/out generator (exact window [] <-> L_p).
- S2: dot multiples Dot_f ∘ ι_lo, ι_lo ∘ Dot_g, π_hi ∘ Dot, Dot ∘ π_hi (products of Schur dots on
  segments, order: total degree, segment, partition), greedy form rank; closure g = Dot ∘ g_0.
- S3: tree leaves at the target (q, s) directly (pairs); for full/end0 the S4 shape below.
- S4: for full/end0: leaves of the trees of src and tgt with equal type r give
  ι_tgt ∘ Dots(L_r) ∘ π_src; plus S6 candidates with dots.
- S5: **not implemented** (no candidates). Its purpose (re-using decompositions of the reverse
  rung) is subsumed by the window tree; recorded as a spec deviation.
- S6: exact window on the maximal non-common window of the best commutation-equivalent pair of
  words, `Local`/`Kron` with commutation isos, all basis elements at the remaining degree, with
  dots.

All window exact calls inside Search are admitted by a predicted cost (cost table when it has the
key, otherwise an unknown-count cap); refused windows are skipped, never waited for.

## 3. Certification (§4.3) and providers

- `hom_gens`, `m_points`, `detecting_cosets` (witnesses), the counting bound with fallbacks A
  (orbits) and B (more points), recorded in `SearchStats::detect`; `stacked_rank` by a dense
  incremental echelon; `certify_full` (|S| = hom_gd, stacked rank, 2 fresh points, genuineness).
  The "residual at 2 fresh points" adds columns, so it cannot lower the rank; it is kept as
  specified (cheap) and documented as such.
- `assert_genuine` = `Arena::assert_genuine` (Elem constants per (kind, labels) / Assoc per window,
  added templates certified) plus a check that every node's degree bookkeeping is consistent.
- Pairs: form at w_q = rep(q) in q's block, rank = mult, the same matrix at a second point; the
  provider returns exactly `mult` ι and `mult` π with an invertible form.
- Closure: greedy degree-by-degree determinant selection on the identity coset; degree multiset =
  c_diag(b_p); nonsingular at two points. Exact closure uses `hom_basis` per needed degree
  ([] <-> L_p windows) instead of `rlambda_bases` (same admission rule; avoids Nakayama over all
  degrees, which is what made P9 slow).
- Exact admission (§4.5): label rule (max λ_i <= 2, fibre product <= 1 600, unknowns <= 5 000)
  or cost table <= 10 s; a cost-table prediction above 10 s vetoes the label rule (the per-call
  ceiling). Unknowns are computed combinatorially (`hom_unknowns`, Gaussian binomials of the fibre
  degrees).
- Closed forms used by Search are certified on demand per (kind, labels) by the A7 comparison
  with the exact unique window map at 3 points (constant must be 1) before `certify_elem`.
- Scratch regions: every provider call is `mark` ... `commit(keep = returned maps)`. Nothing is
  cached across calls in the Arena; window templates are cached in the TmplLib (memory and disk).

## 4. Cost table and G1 harness

- `CostTable`: records (labels, src/tgt fibres, unknowns, degree, seconds, RSS) and keys
  (sorted labels, fibre product rounded up to a power of 2, degree) -> max seconds; text format
  with a header; censored records (watchdog kills) store the kill time as a lower bound.
  Built by `g1 calib` from >= 20 measured calls per label class, each call in a child process
  under the watchdog.
- `g1 a|b [--pairs N] [--seed S] [--json PATH] [--cache DIR] [--costs PATH]`: one child process per
  pair (`g1 one ...`), sequential, nice 19, watchdog polling `ps -o rss=` (kill above 1 GB) and
  wall time (600 s); a kill is recorded as `MissKind::Time` or `Memory` for the running pair.
  Per target: tiers used, candidates against target rank, time, max flows, points; per pair:
  peak index and arena bytes, visited words. G1(b) reruns `Exhausted` pairs with l_max + 2 and
  4·c_max.

## 5. Changes made while implementing (measured reasons)

- **Window splits use "extremes plus dots" locally.** A local family L_t{s_lo..s_hi} of a window
  needs ι only at s_lo and π only at s_hi exactly; the other degrees are exact dot products
  (Schur dots on interior segments, `dot_pmat` = `tensor_exact` of a one-upright scalar template)
  applied to them. Reason: Hom^{+3}([F0(2)], [F0(3) E0(1)]) over (6,0) has 3 139 unknowns, dim 17,
  17-25 s, while its extreme partner has 849 unknowns, 0.3 s.
- **Rex moves for leaf isomorphisms.** A leaf of type q whose word is not commutation-equivalent to
  L_q is connected to it by rewriting indecomposable windows to their canonical local word (the
  unique degree-0 window maps, normalized to be mutually inverse), searched breadth first from both
  ends until the words meet; only then an exact degree-0 window between the best pair.
- **Closure tier S4 (double-leaf shape).** Generators of L_p are π_leaf ∘ Dot_Y ∘ (id ⊗ unit of a
  closed window [R, R^∨]) ∘ g' with g' a generator of a smaller type p' (recursively, depth <= 4)
  and π_leaf from the decomposition tree of Y = L_{p'}·[R, R^∨] (outs symmetric). Reason: the
  lowest generators of the large (2,2,2) types lie in degrees 6-12, where whole-ladder exact calls
  take 30 s to over 60 s (and mirror-twin shapes differ by 20x).
- **Search window admission** (not the Exact backend): an exact-shape cost record decides when the
  table has one; otherwise the bucket maximum may admit (<= 30 s) but not veto; otherwise the
  unknown cap (2 500). Commutation windows have their own cap (20 000 unknowns; measured 2 936
  unknowns in 0.09 s). The Exact backend keeps the spec rule (label rule or bucket <= 10 s; a
  bucket above 10 s vetoes the label rule).
- **Closure of the empty canonical ladder** (the identity type) is the identity map (S0).
- **Cloud session 2026-11-11** (see HANDOFF_M4b.md): `sbim::hom_space` solves its equations with
  `linsys::kernel_lazy` (dense-accumulator elimination, verified early stop; same reduced kernel
  basis as before, so cached templates stay valid), 100-1000x faster on the G1 windows. Caps
  re-measured with it: Exact unknown cap 5 000 -> 40 000; Search window caps 1 600 fibres /
  2 500 unknowns -> 6 400 / 40 000 (`Search::max_window_fibre`, `max_window_unknowns`).
  `certify_full` retries at the full W_lambda-orbits (M+2, then 2(M+2) points) when the planned
  detecting set gives rank < dim. Cost table recalibrated on the cloud box:
  tests/data/homprov/costs_cloud.txt (costs_m4b.txt holds old-solver timings, e.g. censored
  60 s records that veto windows now taking 1 s).
