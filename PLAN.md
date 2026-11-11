# KRcalculator — plan (v2, 2026-11-08)

Goal: compute bigraded, exterior-colored Khovanov–Rozansky gl_N link homology
(ETW-type conventions, see below) over Q and F_p, plus a one-parameter deformation
over Q[t] sufficient for concordance invariants. Priority: complete and correct
first, performance second. Performance target: uncolored gl_4 homology of the Conway
knot (11n34) in under one hour on a personal laptop (Apple M5, 10 cores, 16 GB).

## What was learned in milestones 1–2 (important)

* The "truncated edge algebra" model (bimodules over R_λ = ⊗ H*(Gr(λ_i,N)) given by
  edge rings with facet truncations) is WRONG. Counterexample (N=3): the square web
  E F 1_(2,1) gives a non-self-dual bimodule (graded dim 1,4,6,5,2) violating
  EF 1 ≅ FE 1 ⊕ 1; Hom(R, A_W) has wrong grading. There is no 2-functor
  Foam_N → (R_λ-bimodules) of this kind (R_λ local ⇒ K_0 too small).
* Even closed webs: the truncated edge ring of a closed ladder is not always the
  state space (N=3, bottom (1,2,0), 6 rungs: 1,6,13,15,10,3 vs MOY 2,8,14,14,8,2).
  gl_N state spaces are not cyclic over the edge ring in general.
* "Untruncated invariants modulo facet relations of the web itself" is also too big
  (26 vs 24 in the same example).
* CORRECT MODEL (validated on hundreds of random ladders, N ≤ 6, all labels):
  - Tangle level: untruncated singular Soergel (Bott–Samelson) bimodules B_W of
    ladder webs, a genuine 2-functor (Khovanov–Robert–Wagner 2026, arXiv:2603.22957).
    Represented as free modules over R_λ = ⊗ Sym(λ_i) with polynomial-matrix top actions.
  - gl_N closure = universal construction:
      F_N(closure W) = Hom_SSBim(R, B_W) / radical of (f,g) ↦ τ_N(g∘f),
    g ∈ Hom_SSBim(B_W, R), τ_N = top coefficient of ⊗ H*(Gr(λ_i,N)).
    Graded dims agree exactly with the MOY oracle with normalization
    shift(W) = -Σ_i λ_i(N-λ_i) - Σ_rungs r(b-r).
* Colored Rickard complexes: between consecutive terms the space of degree-1
  bimodule maps is exactly 1-dimensional (checked N ≤ 5, all colors, both signs),
  so differentials are computed, not hand-coded; d² = 0 automatically.

## Conventions (fixed)

* Positive crossing (left strand over, labels (a,b) at the bottom), ETW Fig. 5/(3-2):
  term k = F^{(a-b+k)} E^{(k)} 1_(a,b) (a ≥ b; mirrored if a < b), k = 0..min(a,b),
  homological degree +k, q-shift q^{-(x+k)}, x = min(a,b)(N - min(a,b)).
  Negative crossing: vertical flip, degree -k, q^{x+k}.
  (Read in ETW's own shift convention, ETW's q^{k-x} differs from this only by a
  per-crossing framing factor; this choice is the R1-invariant one.)
* Direction: positive knots in negative q-degree (as KR, Lewark–Lobb), e.g. N=2
  trefoil = q^-1 + q^-3 + t^2 q^-5 + t^3 q^-9 (Khovanov's with q ↦ q^-1).

## Milestones

1. [done] Infrastructure: fields Q / F_p, Laurent polys, sparse & dense linear algebra.
2. [done] MOY oracle (skew Howe quantum trace); relative Grassmannian structure
   constants; singular BS bimodules over polynomial rings; Hom spaces; universal
   construction state spaces. Acceptance: state space = MOY on random webs.
3. [done, replaced] Foams are not hand-coded: Rickard differentials = unique
   degree-1 bimodule maps.
4a. [done] Colored Rickard complexes, d² = 0, uniqueness.
4b. [done for braid closures] Full cube, uncolored and colored. Verified: N=2
   Khovanov (trefoil, figure-eight); trefoil = KR formula for N=2,3,4; Markov
   stabilization invariance; colored unknots.
   TODO: compare with foamho sl_3 tables, colored Hopf link (Hogancamp–Rose–Wedrich),
   colored T(2,n) (J. Wang), knots up to 6 crossings (Carqueville–Murfet).
5. Input: PD code → braid (Vogel's algorithm) + colors per component; braid words. CLI.
6. Efficient engine (audited 2026-11-09). Audit finding: the generic Hom solver
   (unknown |F|×|F| polynomial matrix) does not scale: End^0 of a Bott–Samelson
   bimodule in S_4 takes 0.01 s (length 4), 0.43 s (length 5), > 10 min (length 6);
   w_0 has length 6 and tensoring with a crossing gives length 7. So the engine must
   not solve Hom systems at scale. Design: localized Soergel calculus + light leaves,
   evaluated numerically at M generic points over F_p.
   6A. Hecke algebra oracle for S_n: KL polynomials, KL basis, μ, b_x·b_s, Soergel
       graded ranks of Hom. Accept: known S_3/S_4 values (P_{e,3412} = 1+q etc.),
       b_x b_s formula, agreement with hom_space dims on small Bott–Samelsons.
   6B. Localization of uncolored Bott–Samelson bimodules: fixed points =
       subexpressions, loc maps from the concrete model, block structure by w.
       Accept: L invertible, top action diagonal, localized concrete maps are
       block-diagonal.
   6C. Elementary maps (start/end dots, trivalents, 4- and 6-valent vertices) taken
       from exact small concrete Homs (length ≤ 3), stored as localized rational
       formulas. Accept: Soergel-calculus relations hold numerically (barbell,
       needle, isotopy, two-colour relations).
   6D. Rex moves and light leaves / double leaves in localized numeric form.
       Accept: number of double leaves in each degree = Soergel graded rank (6A),
       linearly independent at generic points.
   6E. Indecomposables B_x (idempotents via local intersection forms) and
       decomposition data for B_x B_s. Accept: localized block ranks = P_{w,x}(1),
       End^0(B_x) one-dimensional, idempotent relations.
   6F. Scanning engine: tensor crossings one at a time, decompose, Gaussian
       elimination of degree-0 isomorphisms; morphisms = small block matrices at M
       points. Closure Φ_N from e-blocks (interpolated) and precomputed
       Hom(R,B_x), Hom(B_x,R). Accept: equals the full cube (4b) on all small
       braids for N = 2..5.
   6G. Colored strands via clasp idempotents (∧^a strand = a parallel strands with
       the thick-edge idempotent) or singular light leaves. Accept: equals 4b on
       colored examples.
7. Performance: Conway knot N=4 benchmark; profiling; parallelism.
8. Deformation over Q[t] (equivariant parameters E_j(t)): τ_N and facet relations
   with E ≠ 0; Smith normal form; reduced theory; s_N / ℷ_N; checks vs Lewark–Lobb.

## Risks

* 6a cost: Hom spaces between BS bimodules of length ~8 over Q[x1..x4] are large
  linear systems; mitigate by computing only degree-0/low-degree pieces, using
  adjunction, caching to disk.
* Small characteristic: Soergel theory changes (p-canonical bases); restrict to Q and
  large p unless validated.

## Status (2026-11-09)

* 6A–6F done and validated: Hecke/KL oracle; localization; elementary maps and
  light/double leaves (548 double leaves checked genuine against exact Homs);
  idempotents for all of S_3, S_4 (block ranks = P_{w,x}(1)); all 84 B_x⊗B_s
  decompositions; scanning engine equals the full cube on 9 test braids (N=2..4).
* 5 done (subagent): PD -> braid (Vogel + simplifier), Jones checks on 18 links.
  CLI: `krhom N --pd ... | --braid ... --strands n [--colors ...]`.
* Conway knot 11n34 (braid 1 1 2 -3 2 1 -3 -2 -2 -3 -3): uncolored gl_4 homology in
  31 s on the laptop (N=2: 0.9 s, N=3: 4.5 s). N=2 Euler characteristic equals the
  Knot Atlas Jones polynomial; full-cube cross-check at N=2 running.
* 6G (colored) via cabling + one B_{w_0(S_a)} insertion per component, divided by
  Π[a]!: agrees in shape with the ETW full cube on (2,1),(2,2),(3,1) Hopf links and
  2-colored kinks for N=3,4; per positive crossing shift t^T q^{-(NT+m(m-1))},
  m = min(a,b), T = ab - m (m(m-1) term only tested for m ≤ 2).
* Next: (a) scalable closure: R-bases of Hom(R,B_x), interpolate the e-block
  matrices of final maps from few points, exact τ_N (the current τ_N weights need
  C(n(N-1)+n-1, n-1) points, infeasible for ≥ 6 strands); lazy indecomposables;
  (b) verify the colored shift for m = 3 and derive it; (c) milestone 8
  (deformation); (d) exact Q results (rational reconstruction or two primes).

## Plan v4 (2026-11-09): links ≤ 12 crossings, arbitrary colors, N ≤ 10

Deformation / concordance postponed. The target is deliberately loose: the goal is
an engine whose cost degrades gracefully in (strands b, colors, N), plus a measured
feasibility envelope, not a promise to finish every (link, coloring, N).

Design consequences
* The tangle-level (Soergel) computation is N-independent: one minimal complex per
  (braid, coloring) serves every N ≤ 10; N only enters the closure (and optional
  pruning of objects that vanish in Foam_N for small N).
* Colors a and N-a are exchanged by reversing a component; use min(a, N-a)
  (re-braid the reoriented diagram), so effective colors are ≤ N/2 ≤ 5.
* Ordinary Bott–Samelson representatives do not scale in b (2^{l(w0)} fixed
  points: b=6 gives 32768, b=7 gives 2^21) and cabling multiplies b by the color.
  The unified core must be the singular (thick) calculus: an indecomposable is
  represented through a small singular ladder (e.g. B_{w0} via full merge/split,
  rank n! instead of 2^{l(w0)}), uncolored braids included.

Milestones (each with automated acceptance tests)
M1  Scalable closure and engine hygiene (current uncolored engine).
    R-bases of Hom(R,B_x) by graded Nakayama on multi-point vectors; interpolate only
    the small e-block matrices of final maps; exact τ_N by polynomial arithmetic.
    Vectorize all localized data over evaluation points; lazy indecomposables; disk
    cache of per-(b, colors) data; two-prime agreement for Q-ranks.
    Accept: identical results on all existing tests; Conway N = 2..10; a 5-strand
    12-crossing knot; point count independent of strand number.
M2  Batch infrastructure. Ingest PD codes of all knots/links ≤ 12 crossings
    (KnotInfo / LinkInfo), braid conversion with per-component orientation choice,
    cost model (b, crossings, colors), resumable job runner, JSON outputs, automatic
    checks (Euler characteristic vs oracle when cheap, invariance under random Markov
    moves on samples). Produce uncolored tables for N ≤ 10 early.
    Accept: uncolored N=2 table agrees with Khovanov homology data; N=3 with foamho.
M3  Decategorified oracle for colors. q-Schur algebra: canonical basis on double
    cosets (singular KL), singular Soergel Hom formula, decomposition multiplicities,
    N-truncation (which canonical basis elements vanish for gl_N). A colored
    MOY / R-matrix evaluator fast enough for Euler-characteristic checks of table
    entries. Accept: agrees with ordinary KL (maximal double coset elements) and with
    hom_space dims and the skew Howe oracle on small cases.
M4  Singular localized calculus. Fixed points = flows (already in local.rs);
    elementary thick maps from exact small Homs; singular light leaves after
    Elias–Ko–Libedinsky–Patimo (JEMS 2025, arXiv:2401.03053) in localized numeric
    form. Accept: double-leaf counts = singular Hom formula; genuineness against
    exact Homs at many points; linear independence. (Main risk: the companion
    "singular Soergel calculus" paper is not out; fallback = spanning Hom^0 by
    random compositions of elementary maps certified by M3 dimensions.)
M5  Singular indecomposables and decomposition data (idempotents via local
    intersection forms; B ⊗ crossing-term webs), lazily generated and cached.
    Accept: localized block ranks = singular KL values; multiplicities = M3.
M6  Singular scanning engine for colored braids (any colors; uncolored large b as a
    special case), optional Foam_N pruning, colored closure.
    Accept: equals the ETW full cube and the cabling engine on small colored links,
    equals the uncolored engine; colored unknots, Hopf links (Hogancamp–Rose–Wedrich),
    T(2,n) (J. Wang); ∧^2-colored gl_4 Euler characteristic vs the so(6)
    specialization of the Kauffman polynomial (normalization to be fixed).
M7  Performance and envelope: profiling, parallel scheduling, run the tables in
    increasing cost, publish what was computed and where it stops.
Later: deformation over Q[t] and concordance invariants; colored HOMFLY via the
    Hochschild closure of the same minimal complexes.

## M1 status (2026-11-09)

Done: exact closure from R-bases of Hom(R,B_x), Hom(B_x,R) (graded Nakayama),
interpolated Gram G_x and edge matrices F (adaptive point count, residual checks),
τ_N by coefficient extraction, one scan for all N; lazy indecomposables; rank via
pairing matrices (no inverse Gram), packed monomials, term-ordered columns
(near fill-free elimination), parallel (h,q) blocks; two-prime agreement; CLI takes
an N list; regression test vs the full cube.
Measurements (Apple M5): Conway knot all N = 2..10 in 105 s (N=4 < 1 s). Conway
N=2 full cube (independent model, 2048 vertices, 4844 s) gives the identical
Poincaré polynomial. 5-strand 12-crossing braid: scan 14 s; closure N=4 1.1 s,
N=5 12.6 s, N=6 93 s.
Not done: vectorization over points and disk cache (scanning is no longer the
bottleneck).
New bottleneck: the closed complex has dimension ~ N^(strands) while the homology
is linear in N. Remedy (proposed M1b): close strands one at a time with the Foam_N
partial trace (Markov trace) and Gaussian-eliminate after each step, so complexes
stay small; this is the Bar-Natan principle of closing as early as possible.

## M1b decision (2026-11-09, design panel: 4 designs + 4 adversarial critiques + synthesis)

* Chosen: algebraic "restriction tower" closure (close one end variable at a time by
  restriction of scalars A_k -> A_{k-1}, Hankel/Toeplitz matrices, Gaussian elimination
  over A_{k-1}; exact at stage 0). Proofs (P1-P6) checked by two reviewers; correctness
  does not depend on splitting quality. Not the refuted relative trace (B_s gives N-1).
* Reserve: Foam_N partial-trace Smith data (for non-split types), SAC extrapolation
  (knots only, labelled "predicted", never mixed into computed tables), drop sets.
* Synthesis: scratchpad closure_synthesis.md; full designs/critiques: design_panel.json.
* M0 done: removed dead inverse Gram / index map, scan data freed after interpolation,
  per-variable cap vector, reduced mode (marked strand cap 0, reduced unknot = q^0),
  phase timings, KRHOM_THREADS. Checks: reduced trefoil = q^{-2(N-1)} + t^2 q^{-2(N+1)} +
  t^3 q^{-4N}; Conway reduced ranks 33, then 49 for N>=3; Lewark-type relation
  [N]·P_red - P_unred divisible by (1+t) with nonnegative quotient (N=2..6).
  Conway N=2..10 with 2 threads: unreduced 85 s / 995 MB peak; reduced 2.5 s / 86 MB.
* New findings: closure-data preparation (bases + Gram) is slow for 5 strands (~160 s);
  the ordinary Soergel scan of a 6-strand 15-crossing braid reached ~77 GB, so 6-7
  strands need the singular calculus (M4-M6) regardless of the closure.

## M1b.1 status (2026-11-09): restriction tower implemented (src/tower.rs)

* `src/tower.rs`: packed polynomials (6 bits per variable, n ≤ 10, untruncated products
  asserted to degree ≤ 63, truncated products tested by one add + mask), Hankel
  restriction of one end variable (order Right / Left / Alternate; `Middle` is a
  diagnostic interior-first order), 1x1 exact-division pivots (unit pivots first), junk
  coordinates, connected components of the remainder, Markowitz-ordered elimination with
  lazy cost updates; isomorphism test = matching coordinate degrees + invertible F_0 +
  equal dimension (same interned generator matrix, else `gdim` by recursive
  decomposition); inverse over A by the Neumann series; stage 0 asserts that no map is
  left and every object is one-dimensional. Decompositions and dimensions are cached per
  interned generator matrix; N values run sequentially, caches freed per N.
* `src/engine2.rs`: `prepare_closure_data` (scan, bases, Gram, F; scan data freed) +
  `close_with(data, N, reduced, Closure::{Exact, Tower(end)})`; `close_exact` kept as the
  oracle; `khr_engine2` / `khr_engine2_mode` unchanged (exact); `khr_engine2_closure`
  chooses. CLIs: `engine2 ... --closure tower|exact|both --end right|left|alt --reduced
  [--unit-pivots]`, `krhom ... --closure tower|exact --end ... --reduced` (default tower).
* Tests (tests/acceptance.rs, all against close_exact or independent data): tower == exact
  on the full-cube regression braids, trefoil, figure-eight, Hopf, T(2,-4), a 3-component
  chain and the 2-component unlink, N=2..6, unreduced and reduced, right/left/alternate
  ends; tower == full cube (N=2,3) inside soergel_engine_matches_full_cube; Conway
  N=2..10 ranks 34..426, ends agree, tower == exact for N ≤ 8 (also checked once at
  N=9,10); reduced trefoil formula N=2..10; Conway reduced 33, then 49, == exact N=2..10;
  Lewark check ([N]P_red − P_unred divisible by 1+t, quotient ≥ 0) for Conway and three
  small knots, N=2..6; chain level: tower gdim per type == close_exact state-space dims,
  dim_unreduced = N·dim_reduced per type (Conway N=2..5); stress mode with unit pivots
  only (non-split blocks) == exact; synthetic non-minimal complex with a non-cyclic
  isomorphism between different generator matrices (Neumann inverse, gdim certificate,
  elimination update) == exact; second prime F_4294967291. Unit tests: division, caps
  test, end orders, inverse over A, dim Φ_N(B_s) = N(N-1); 6-strand braid with a
  non-splitting right end (below). Full suite 43 tests, 35 s wall with KRHOM_THREADS=2
  (peak RSS 0.9 GB with the tests running concurrently).
* Measurements (Apple M5, KRHOM_THREADS=2, nice 19, one process per row, Conway N=2..10):

  | closure | mode | N=10 closure | whole run N=2..10 | peak RSS |
  |---|---|---|---|---|
  | tower | unreduced | 0.11 s | 1.55 s (scan 1.1 s) | 54 MB |
  | tower | reduced | 0.11 s | 1.71 s | 62 MB |
  | exact | unreduced | 79 s | 130 s | 951 MB |
  | exact | reduced | 1.06 s | 3.97 s | 83 MB |

  Stage sizes at N=10 (right end): 1890 objects / 31672 maps / 87k terms → 120 → 68 → 49
  → 426 (as predicted by O4). 5-strand braid "-3 -4 2 -3 2 1 3 4 -1 4 -1 -3": tower ==
  exact for N=2..5, ranks 10 + 9(N−2) up to N=10, N=10 closure 0.06 s (2356 → 56 → 36 →
  18 → 9), but data preparation 19 s (scan 8 s, bases + Gram 6–7 s): the end-to-end
  bottleneck is now M1c.
* 6-strand braid "-2 -3 -1 -5 2 4 -5 -2 1 -2 1 -2 5" (92 final objects): data 30 s (scan
  12 s, bases + Gram 18 s at 800 points), 156 MB. Tower == exact for N=2,3,4 (exact N=4:
  3.0 s); ranks 3N−2; right and left ends identical for N=4..10. Right end: non-split
  blocks as predicted by the critique (156 → 876 non-split coordinates for N=4 → 10,
  peak 6662 objects before elimination), N=10 closure 2.37 s, 405 MB; left end: no
  non-split coordinates, N=10 0.01 s, 179 MB. So the adaptive end choice of M1b.2 matters
  for cost, not for correctness. (Correction after verification: the right end grows like
  N^7–7.5 here for N=7..10, not ~N^4.)
* Observed: with exact-division pivots every type of the 4- and 5-strand examples split
  into cyclic children at either end; only one real isomorphism between non-cyclic
  objects occurred (6-strand braid, right end), so the synthetic test covers that path. The interior-first order (`Middle`) and
  the unit-pivot stress mode leave hundreds to thousands of non-split coordinates and
  still agree with exact (stress mode at N=6 peaks at 1.3 GB: correctness does not depend
  on splitting, cost does).

## M1b.1 fixes after independent verification (2026-11-09)

* Oracle crash (wrong-result risk): `close_exact` packed monomials with 5 bits per variable;
  the sum of a representative and a P_ab term (each ≤ N-1 per variable) carried into the
  next field for N ≥ 17 (bogus columns, RankEchelon index panic, possibly wrong ranks).
  Now u128 packing with w = bitlen(2·max cap) bits per variable (≥ 5), asserted n·w ≤ 128;
  dual index keyed by (object, b, monomial). Regression tests: monomial-map complexes on
  2 and 3 variables at N = 16, 17, 18, 20, 32 (closed-form rank, tower == exact) and N = 64
  (oracle only, beyond the tower's N ≤ 32); figure-eight and "1 1 1 2" at N = 17, 20, 32.
* Tower memory: `norm` kept the capacity of the uncombined product list (up to |a|·|b|)
  in every stored polynomial (5-strand N=10: 4.2 GB RSS, swapping). It now compacts in
  place and shrinks.
* Tower speed: long homogeneous products and long accumulations are combined by direct
  addressing in the exponent box (lowest varying variable implied by the degree; index
  order = packed order, so no sort), falling back to sorting when the box is large; the
  untruncated Schur update of `decompose_raw` is fused (a − q·b by shifted merges for a
  short factor). `gdim` is keyed by the set of closed variables and no longer follows the
  tower's order through non-split blocks: it takes the next tower variable if that splits
  completely, else an end of the open interval that does, else the fewest non-split
  coordinates; constant matrices and 1x1 monomials are counted directly.
* Data preparation: light-leaf matching uses only (end element, degree), computed without
  braid-move searches; op lists are built once per occurring leaf; `rex_path` is memoized
  per context; `closure_gens` enumerates all degrees in one pass (it used to redo all 2^ℓ
  light leaves for each of 4ℓ+5 degrees). Gens-only interpolation points are evaluated
  only for the types whose own Gram weight needs them (5-strand: one type needs 723
  points, most need ≤ 218). Closure data are unchanged (byte-identical files).
* `engine2 --save-data FILE / --load-data FILE` (text format, prime fields): re-close a
  prepared braid without rescanning (`save_closure_data` / `load_closure_data`).
* Measurements (Apple M5, KRHOM_THREADS=2, nice 19, one process per row):

  | case | before | after |
  |---|---|---|
  | 5-strand "1 -2 3 -4 2 1 -3 4 -2 3 1 -4", tower right N=7/8/9/10 | 0.95/2.56/6.42/15.1 s, 4.2 GB RSS | 0.62/1.44/3.23/6.91 s, 0.96 GB RSS |
  | same braid, data preparation | 262–410 s | 114–153 s (bases+Gram 89–120 s; varies with machine load) |
  | 6-strand "-2 -3 -1 -5 2 4 -5 -2 1 -2 1 -2 5", right end N=10 | 2.87 s (x6 stage 2.63 s) | 1.92 s (x6 stage 0.012 s) |
  | same braid, data preparation | 30–37 s | 10.5–12.7 s (bases+Gram 2.9 s) |
  | Conway N=2..10 tower unreduced / reduced (whole run) | 1.55 / 1.71 s | 1.24 / 1.35 s, 50 / 60 MB |
  | Conway N=2..10 exact unreduced / reduced (whole run) | 130 / 3.97 s | 110 / 2.95 s, 1.04 GB / 96 MB |

* Tests: 48 (17 unit, 20 acceptance, 11 pd), all pass; new: the two oracle regression
  tests, closure-data save/load round trip, and unit tests of `norm` (dense and sorted
  paths), fused and dense products (homogeneous or not, truncated or not).
* Still open (not clear fixes): when neither end splits (5-strand braid above) a non-split
  block at one stage becomes thousands of cyclic children with dense block-Toeplitz maps at
  the next (x4 stage at N=10: 6010 objects, 541k maps, 12.7M terms), so the closure still
  grows ~N^7; this needs block pivots (or another reduction of non-split blocks), not the
  adaptive end choice. For the 6-strand right end the cost moved from `gdim` to the x4
  stage (decomposing the x5 non-split blocks, which the tower needs anyway); the left end
  avoids it. The remaining preparation cost is the per-point idempotent of the largest
  types (`indecomposable`: all degree-0 double leaves as 2^ℓ-column localized matrices),
  e.g. B_{w0} on 5 strands as a start object; that is M1c. `Middle` stays a diagnostic order.

## M1b.1 wrap-up (2026-11-09 evening)

- The verification workflow was stopped during its second fix round because its local experiments
  (up to 2.5 GB, earlier 6.7 GB) made the laptop lag. Sources were restored to the state that
  passed re-verification 1. The unfinished round-2 work (2x2 block pivots in tower.rs, which made
  the 6-strand right end split; a 5-strand attempt in progress with debug code) is archived, not
  merged, under the scratchpad folder `fix2_partial/` (tower.rs.round2, tower_round2.diff,
  acceptance.rs.round2). It is the starting point for M1b.2 block pivots.
- Fixed by hand after re-verification 1, each with a regression test:
  - close_exact's state-space degree loop now runs over [min(ins), topdeg - min(outs)]; the old
    upper bound dropped top degrees when min(outs) < min(ins) (never in braid data).
  - The tower adds parallel edges (same source and target) like close_exact; it used to keep only
    the last one (never in braid data).
  - prepare_closure_data counts the extra scan points needed for F as scan time and prints them.
- cargo test --release: 50 passed (17 lib, 22 acceptance, 11 pd), peak RSS 0.8 GB.
- Open: 5-strand braids whose first end does not split grow ~N^6-7 in the tower (block pivots,
  M1b.2); data preparation is slow and erratic for long 5-strand words (M1c); B_{w0} start
  objects on 5 strands do not get through preparation.

## Priorities after M1b.1 (2026-11-09, user)

- Primary goal: compute colored (wedge^a per component) gl_N link homology itself. Engineering and
  scale are secondary (decided by the implementer); concordance invariants are a side goal.
- Consequence: the critical path is the singular (colored) engine (old M3-M6). M1b.2 block pivots
  and M1c data-preparation speed for uncolored 5-strand braids are demoted to secondary tasks.
- Scale reality check for the showcase target (5-colored gl_10 of the Conway knot): Wedrich
  (arXiv:1602.02769, Thm 4) gives reduced wedge^5 HOMFLY-PT rank >= 49^5 for large N; 2-colored
  measurements (cabling engine: unreduced rank / C(N,2) at N=2k is about 30-40% of r^k) suggest
  about 1e7 reduced and 1e9-1e10 unreduced generators at N=10. Stretch goal, not the next milestone.
- Localized singular calculus: rank of an indecomposable = number of cosets in its double coset,
  at most prod (a_i)! (2-colored 4 strands: 16; 3-colored: 1296; 5-colored 4 strands: about 1.3e7).
- Design workflow (research, three competing designs, adversarial critiques, synthesis) started;
  outputs in the scratchpad folder colored_design/.

## Colored engine: design decision and first milestones (2026-11-10)

- Design workflow (3 research notes, 3 competing designs, 3 adversarial critiques, synthesis; kept
  outside the repo in the session scratchpad, colored_design/SYNTHESIS.md). Chosen: D1, localized
  singular Soergel calculus (objects B_p for contingency tables p, coset-wise localized maps, scalar
  elimination, Hom bases certified against the EKLP Hom formula), with D2's induction bridge as an
  oracle layer and D3's Hom-basis coordinates as a later scale-out layer. Order:
  M3a (q-Schur oracle) || CL1 (colored closure) -> M6.0 (exact small colored engine) ->
  M4a/M4b/M5 (localized substrate, Hom provider, decompositions) -> M6a (2-colored 3-strand,
  3-colored 2-strand) -> M6b/M6c -> CL2/CL3/M7a (scale-out, gated by measurements).
- Corrected scale: the localized rank of B_p is sum_{q<=p} h(1) N(q), not the coset count
  (2-colored 4 strands: max 25 920; 5-colored 4 strands: single objects >= 1.2e10). 2-colored Conway:
  >= 11 029 objects in any minimal complex, peak about 1e5 objects / 3.4e8 sum S per point,
  1.6e5 closure unknowns per Gram entry: the hardest near-term target. 5-colored on >= 3 strands is
  out of reach for every design (showcase stays a stretch goal).
- Multi-strand oracle: duality at N = a+1 (wedge^a gl_{a+1} = uncolored gl_{a+1} on the reversed
  braid; expected from Queffelec-Rose, labelled "via duality (expected)").
- Cabling reference tool: src/bin/colored_cable.rs, krhom --colored-method cable|cube; equal to the
  full cube for colors (2,1), (2,2), (3,1); (3,2) and (3,3) cables do not fit in memory.
- M3a done: src/schur.rs (contingency tables, canonical bases, minimal-flow ladders, EKLP Hom formula,
  expansions, Rickard rigidity, decategorified scans with lower bound and peak, Gram weights,
  decategorified Euler characteristic), schur_scan CLI, tests/schur.rs. Identical to the Python
  reference on all fixtures; (2,2,2,2) basis 1.7 s / 100 MB; 2-colored Conway scan with peak 2.2 s.
  T-n confirms v -> q. 5-strand bases grow fast ((2,2,2,2,1) > 1.5 GB): cloud.
- CL1 done: src/cclose.rs (ring R^lambda, contract v2, exact colored closure oracle over
  (x)_i H*(Gr(lambda_i, N)), reduced mode, global_q, (Delta G) adapter to the restriction tower,
  exact closure data of small ladders), g4_free (gate G4 certified for 2 <= a <= 5, a <= N <= 10),
  tests/cclose.rs. Checked against MOY, StateSpace, Wang's Hopf formula, the full cube and engine2.

## Target narrowed: 3-colored gl_6 (2026-11-10, user)

- For now the showcase is wedge^3-colored gl_6 homology (knots, up to Conway size), instead of
  5-colored gl_10. First 2-strand T(2,n) (Wang: 3-colored trefoil gl_6 rank 136), then 3-strand knots.
- Fixed N = 6 allows (to be validated): (1) N-truncation: drop every type B_p whose Lusztig cell weight
  mu(p) (conjectured: RSK shape of the column-reversed table) has mu_1 > N, at every stage of the scan;
  (2) closure only at N = 6, with unknowns restricted to the Gr(3,6) box (<= 20 monomials per block).
- 2-strand (3,3) at N = 6 prunes nothing (max internal label 6), so M4b (fast certified Hom spaces) stays
  on the critical path.

### N-truncation validated (2026-11-10, branch trunc-n)

- Rule: B_p maps to 0 in Foam_N iff mu(p)_1 > N, mu(p) = RSK shape of the column-reversed table
  (= RS(p_max)^T; Greene: mu_1 = largest entry sum along a down-left chain, an O(nm) DP).
  Decategorified: 0 mismatches against the MOY trace kernel on 9 833 (type, N) pairs (51 (lambda, N) rows,
  up to (3,3,3) and (2,2,2,2)); literature: CKM Thm 4.4.1, Lusztig 29.2 / Doty Thm 4.2, Cui-Luo-Wang
  Prop 5.6(3). The categorified statement is assumed via the engine's Hom formula; dropping killed
  types mid-scan only needs the decategorified kernel + positivity + the closure contract.
- Categorified: pruned = unpruned 2-strand scans on every feasible case (tests/ntrunc.rs); new:
  3-colored Hopf/trefoil/T(2,5) at N = 3, 4 (Wang, Hopf formula, chi, duality with uncolored gl_4).
- Do NOT use "the minimal-flow ladder has a label > N" as a kill criterion (false for (3,3,3,1)).
- Envelope at N = 6: 3-strand (3,3,3) knots LB ~0.49x, kept-summand peak ~0.36x of untruncated
  (pre-decomposition tensors still ~0.62x); 4_1 LB 276, 8_19 LB 776, 12a_1248 LB 7.4e5.
  (3,3,3,3): 790 of 2 008 types survive; closure in the Gr(3,6)^4 box instead of 1.9e9 e-monomials
  (weight vs box-degree convention still to be checked for weights 37-40).
- N = 5 on (3,3) still needs a 5401 s Hom between kept types; N = 6 prunes nothing on 2 strands:
  M4b remains the critical path. Pruning must be built into M5/M6a (decompositions drop killed summands).
