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
