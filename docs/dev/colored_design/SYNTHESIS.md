# SYNTHESIS: colored gl_N engine, chosen design, plan, first task

Date: 2026-11-10. Author: synthesis agent S.

Inputs, all read in full:
- `PLAN.md` (What was learned, Conventions, Plan v4, M1b decision, M1b.1 status and wrap-up, Priorities after M1b.1);
- research notes R1 (singular theory), R2 (code map), R3 (closure and validation);
- designs D1, D2, D3 and their critiques C1, C2, C3;
- the reference scripts that the first task ports: `scripts/schur.py`, `t8.py`, `t3.py`, `t10.py`, `cmp_hom.py`, `proto_D1/peak.py`, `proto_D1/decat_scan.py`, `proto_D2/transfer.py`, and `r3proto/src/bin/chi.rs`.

New numbers in this note:
- e-monomial counts and Grassmannian fixed-point counts (`colored_design/synth/emono.py`);
- Euler characteristics from the `chi` binary (Appendix A.3).

No repository file was edited.

Tags (as in R1-R3 and the designs):

| tag | meaning |
|---|---|
| [L] | literature, theorem number given |
| [V] | verified numerically (by whom is said) |
| [D] | derived, not independently checked |
| [?] | open |
| [E] | estimate |

Abbreviations: B_p is the singular indecomposable for a contingency table p; K = sum of labels; R^lambda = (x)_i Sym(lambda_i); A_lambda = (x)_i H*(Gr(lambda_i, N)); "LB" = Grothendieck lower bound.

---------------------------------------------------------------------------------------------------

## 0. Decision in one page

**Chosen design.** The backbone is D1: localized singular Soergel calculus.
- Objects: B_p, represented by a canonical *minimal-flow* reduced ladder L_p and an idempotent e_p.
- Maps: evaluated one coset at a time; scalar Gaussian elimination.
- Crossings: the half-step pipeline (split = re-indexing, merge = decomposition by adjunction and local intersection forms).
- Hom bases: certified against the Hom formula.
- Closure: the universal-construction closure over R^lambda.

D1 is adopted with the corrections the critiques required, and with grafts from D2 and D3.

| from | what is grafted |
|---|---|
| **D2** | The induction bridge, used as an *oracle layer*: Hom transfer F1, the [mu]! law F2, closure transfer F5. The existing, validated `tower.rs` run on (Delta_lambda G, F) serves as an interim production closure for K <= 10 (C2 M1(c), gated on a freeness check G4). Genuineness and point-independence discipline. Per-coset B-coordinates (U_c, V_c). A Hom provider with an exact backend and a certified-composition backend. |
| **D3** | Coordinates in certified Hom^e bases as the *scale-out* representation ("R-coord"). Detecting sets (Lemma L3) with residual checks. The "intrinsic basis words only" rule. Closure data from identity-coset samples, cached per type and per basis word. The chain-level comparison "multiset of (type, shift, h) equals engine2's". |
| **C1-C3** | Stacked multi-point evaluation everywhere linear independence over a Hom space is decided (one point is not faithful on Hom^d, C1 I1). Closure-first ordering with an exact small engine (C3 M1). Polynomiality certificates for every template and closed-form map (C3 H3). The decategorified peak model (C2 H2). The streaming-elimination rule (C1 I5). Duality at N = a+1 as the primary rank-level oracle (C2 M3, C3 M2). Grassmannian-box-restricted fixed-point interpolation (C1 I2(b)). Explicit measured gates. |

**Rejected:**
- one-point coordinate extraction and one-point certification (C1 I1: false inference, verified numerically);
- D3 as the primary architecture (C3: it front-loads certified bases for *every* receiving Hom space before any result);
- D2 as a separate engine (C2: its folded data model *is* D1's R-full representation);
- learning rex moves by `hom_space` beyond 2 colors (C2 H3(b));
- dense Vandermonde interpolation of closure data beyond about 2e3 unknowns (C2 H1).

**Order of work.**
1. **M3a** decategorified oracle `schur.rs`: the FIRST TASK, section 8. In parallel, **CL1** colored closure layer.
2. **M6.0** exact small colored engine: first colored homology through the new path, bigraded, checked against Wang, the Hopf formula, the full cube and cabling.
3. **M4a/M4b/M5** localized substrate, Hom provider, decompositions; gate G1 on the certified-composition ("Search") provider.
4. **M6a** R-full localized engine on <= 3 strands with multi-point closure: 3-colored 2-strand, 2-colored 3-strand knots and links. This is the first genuinely new data.
5. **M6b/M6c** validation sweep, uncolored regression, M-point scans with final-component coordinates.
6. Scale-out, each step gated by measurement: CL2 structured closure interpolation, CL3 colored tower, M7a R-coord.
7. **M4c** EKLP singular light leaves only if the Search provider fails its gates.

**Realistic deliverables of this plan** [E; section 1.3 has the table]:
- all 2-strand colored braids with colors <= 3, at all N <= 10 (laptop);
- 2-colored knots and links on 3 strands up to 12 crossings, all N <= 10 (laptop, table runs on cloud). N >= 7 needs the (Delta G) tower after G4, or CL3;
- 3-colored 3-strand knots: small N on the laptop, N = 10 on the cloud after CL2 and CL3;
- links with one or two thick components and the rest label 1, on 4-5 strands: probably cheap, measured in M3b;
- 2-colored 4-strand knots (Conway), reduced theory at N <= 6: plausible after R-coord (M7a).

**Not promised:**
- 2-colored 4-strand at N = 10 unreduced;
- 2-colored 5-strand and 3-colored 4-strand at N = 10. Closure data alone reach 1e7-1e9 unknowns per Gram entry (C1 I2, C3 E4) and need research, not engineering.
- 5-colored on >= 3 strands, including the 5-colored gl_10 Conway knot. Single objects have localized rank >= 1.2e10, and the reduced answer is >= 49^5 for large N by Wedrich arXiv:1602.02769 Thm 4.

The showcase stays a stretch goal. All three designs and all three critiques agree on this.

**Correction to the user-level estimates.** In the previous message, 2-colored Conway was called "about 1e5, easy once the engine is ready". The *homology* is about 1e5 unreduced at N = 10. But the tangle-level minimal complex has at least 11 029 objects of localized dimension >= 1.33e7. The scan peaks at about 1.0e5 objects with sum S(q,q) = 3.4e8, which is about 7-10 GB of localized map storage per evaluation point (C2 H2). One closure type needs 1.6e5 interpolation unknowns per Gram entry. It is the hardest near-term target, not an easy one (R1 §6.4, C1 I2, C2 H1-H2).

---------------------------------------------------------------------------------------------------

## 1. Why this design

### 1.1 What the three designs share (and is therefore settled)

Every design keeps the following:
- the validated model: untruncated SSBim plus the KRW universal-construction closure;
- the N-independent tangle-level scan;
- the colored Rickard complexes of `crossing.rs` (convention A) with unique differentials (R1 §5.1 [V] for c, d <= 5, both signs);
- the closure contract of `engine2` with the ring changed to R^lambda (R3 A3, proved);
- the Hecke algebroid (q-Schur) oracle for types, characters, Hom formula and multiplicities.

They differ only in three places:
1. how morphisms are represented (fully localized blocks; folded blocks, which are the same thing; or coordinates);
2. where Hom bases come from;
3. how closure data are obtained.

So the synthesis is mostly about choosing among these three axes and fixing their flaws.

### 1.2 Why D1 as the backbone

- **Simplest correctness argument.** Every number is a localized value of a genuine bimodule map built from validated pieces (`local.rs`, `sbim::hom_space` templates, closed-form elementary maps with certificates). The scan needs Hom bases only for decompositions and the closure, not for every stage-level Hom space (C3 H1(a)).
- **Independently checkable.** A fully localized engine can be compared object by object and map by map with the exact small engine M6.0, with engine2 (uncolored), and with the cabling engine (K <= 5).
- **The critiques converge on it.**
  - C1: "adopt with changes".
  - C2: D2's folded model *is* D1's R-full, and D2 is best kept as the bridge and oracle layer.
  - C3: D3 should be "a variant of D1's R-coord layer". D3's components K1-K5 fit there after the needed Hom spaces have been enumerated from real scans.
- **R-coord is the only representation that might reach 2-colored 4-strand braids** (C2 H4, C3 K1). It is planned explicitly (M7a) but gated: C3 measured that type-level cache reuse collapses to 1.3-2.4x per braid once there are about 100 types (C3 E3). Its time advantage over R-full is unproven (C3 H4); its memory advantage is clear (3.5-43x on uncolored scans, C3 E1).

### 1.3 Envelope (corrected; every colored row is an estimate unless marked)

Numbers are from R1 §6.4, D1 §8.3, D2 §6.2, C1 I2, C2 H1-H2 and C3 E4, recomputed or cross-checked there.

| case | composition | LB objects / localized dim | peak objects / sum S per point | max Gram unknowns (unreduced / reduced) | verdict |
|---|---|---|---|---|---|
| 2-colored 2-strand (Hopf, T(2,n)) | (2,2), (2,1) | 11 / 61 (trefoil) | 81 / 495 | 14 | laptop, M6.0 (exact) |
| 3-colored 2-strand | (3,3) | 24 / 507 (trefoil) | 802 / 18 489 | 158 / 12 | laptop, M6a |
| 4-5-colored 2-strand | (5,5) | 82 / 32 151 (5-colored trefoil) | - | up to 40 466 | later (M8): (5,5) Rickard words need certificates; closure needs CL2 + CL3 |
| 2-colored 3-strand, <= 12 crossings | (2,2,2) | 77 / 3 049 (4_1); 8 835 / 5.5e5 (10 crossings) | 33 369 / 2.9e6 (10 crossings) | <= 1 218 | laptop per knot (about 0.1 GB per point); table runs on cloud |
| 3-colored 3-strand | (3,3,3) | 328 / 3.5e5 (4_1); 1 084 / 2.6e6 (8_19) | 5 036 / 2.0e7 (4_1) | 88 233 / 3 866 (4_1) | small N laptop (box-restricted: 64 / 1 000 / 8 000 points at N = 4 / 5 / 6); N = 10 cloud |
| one thick component, rest label 1, 4-5 strands | e.g. (2,1,1,1), (2,2,1,1) | not computed (M3b) | - | W = 9 gives about 4e2; W = 13 gives about 3e3 | probably laptop; measure in M3b |
| 2-colored 4-strand (Conway) | (2,2,2,2) | 11 029 / 1.33e7 | 103 573 / 3.4e8 (7-10 GB per point R-full) | 159 952 / 13 377 | needs R-coord (M7a) and CL2. Reduced, N <= 6: plausible on cloud. N = 10 unreduced: open |
| 3-colored 4-strand; 2-colored 5-strand | (3,3,3,3); (2,2,2,2,2) | not computed (M3b, cloud) | - | up to 1.9e9; up to 1.0e8 | out of plan unless G0 shows that high-weight types do not occur |
| 5-colored, >= 3 strands | (5,5,5), (5,5,5,5) | - | single objects >= 7.6e5 / 1.2e10 | - | out of reach for every design |

**Two cheap facts change what is computable soon.**

1. **Duality at N = a+1** (R3 V3; categorified version expected from Queffelec-Rose, not proven here [?]). For a knot K = cl(beta) colored a:

       KhR_{a+1}(K; wedge^a) = KhR_{a+1}(cl(rev beta); wedge^1)    (bigraded, no shift)

   The right-hand side is uncolored and is computed today by engine2 on the reversed word, for every knot in `data/knots_upto12.jsonl`. This is the primary rank-level oracle for multi-strand colored braids. Results obtained this way must be labelled "via duality (expected)", never "computed".
2. **Mixed colorings are much cheaper.** The closure weight bound is W = C(K,2) - sum_i C(lambda_i, 2), and K is small when most labels are 1. Links with one thick component are common in the link table.

### 1.4 Disputed points and their resolution

| topic | D1 / D2 / D3 position | critique | resolution in this plan |
|---|---|---|---|
| faithfulness at one point | D1: orbit-reduced witnesses at one point decide Hom^d | C1 I1: false for Hom^d. ker(ev_x) = (+)_j (R^lambda_{d-d_j} ∩ m_x) b_j. Verified 9/9 cases | All rank, coordinate, certification and structure-map decisions use **stacked evaluation at M(p,q,d) + 2 points**, where M = max_j dim R^lambda_{(d-d_j)/2} with generator degrees d_j from the oracle. Witness cosets (one per double coset) remain a constant-factor saving within each point. Single-map statements (zero test of one map, top-stratum scalars, pivots, minimal polynomial of one element) may use one point. |
| idempotents | D1: re-solve e in a basis at one point | C1 I1: silently non-idempotent elsewhere | e_p = P(a), a polynomial in one fixed random combination a of certified End^0 words (seeded). Never re-solved. Assert e^2 = e at every new point where e is evaluated. |
| scans per braid | D1: one point + coordinate extraction; D2: multi-point or type-local tables; D3: one scan, coordinates | C1 I1(c), C3 H4 | M6a: multi-point scans (engine2 style), enough for 2-strand braids and (2,2,2). M6c: an M-point scan (M = tens to hundreds) plus coordinates of the *final* components only, against certified bases of the final (p, q, e) spaces, with closure tables per basis word. M7a: full R-coord, gated by G3. |
| Hom generator | D1: numerical EKLP leaves in 5-6 days; D2: B1 exact, B3 search, B4 EKLP; D3: Tier 1 search, Tier 2 EKLP | all three: unproven spanning, budget 2-3x, no gate | Backends in order: exact `hom_space` (small fibres), then **Search** (certified compositions, specified as an algorithm in §2.4, with budgets), then EKLP (M4c) only if gates G1/G2 fail. Certification = count from the Hom formula plus stacked rank plus a fresh-point residual plus genuineness certificates of the generators. |
| rex moves | D2: learned by hom_space on <= 3 uprights | C2 H3(b): infeasible from 3 colors (fibres up to 8 000) | Closed forms: associativity is the identity in value coordinates (E4 [D]); square switches are composites of E2/E3; singular braid relations follow the explicit diagrams of EKLP (102)-(106) (as cited by C2; numbers to be re-checked when implemented). `hom_space` only *verifies* them on small labels. Any degree-0 map with nonzero top scalar is an isomorphism on B_q (C2 M6), so rex paths are an optimization. |
| genuineness | D3: rank certification suffices | C3 H3: rank does not certify genuineness; a Rickard word with d = 1 passes for any nonzero candidate | **Polynomiality certificate** for every window template, closed-form op and composite Rickard differential: Phi = L'^{-1} f L must interpolate to polynomials of the predicted degrees in the window's bottom e-variables, with residuals at fresh points (window fibres <= 630). Also d^2 = 0 on samples for every Rickard word. |
| point-independent scalars | D1: "normalized by uniqueness" | C1 I7 | Every template carries a constant scalar fixed once (by comparison with the exact PMat, or a closed form). "Make this entry 1 at the current point" is forbidden. |
| closure algorithm | all three: a new colored tower | C2 M1(c): tower.rs on (Delta_lambda G, F) divided by gr(V) is provably correct given freeness | Reference: the colored exact oracle (Schur/GrTable coordinates, CL1). Interim production for K <= 10: tower.rs on (Delta G, F), only after G4 (freeness of Delta_a R_N over H*(Gr(a,N)) checked for all a <= 5, N <= 10; C2 verified it in some ranges). Colored tower (CL3) is a performance milestone validated against both. |
| closure generators | Nakayama via `select_basis` at many points | C2 M2: one generic point suffices (determinant argument) | Choose n = c_diag(b_p)(1) genuine maps with the right degree multiset whose n x n identity-coset value matrix is nonsingular at one generic point; re-verify at a second point. Regression test against `engine2::select_basis` on the uncolored suite. |
| closure interpolation | dense Interp; "thousands of points"; factorial Schur | C1 I2, C2 H1, C3 H2: 1e5-1e9 unknowns, O(m^3) | Dense Interp only up to about 2e3 unknowns. Beyond that, CL2 (proposed in C1 I2(b) and C2 H1(a)): interpolation at Grassmannian torus fixed points z_{mu+delta} (factorial Schur vanishing), restricted per N to the lambda_i x (N - lambda_i) boxes and weight <= w, solved per block by triangular transforms on the lower set [D]. A spike measures it before any target needs it. Reduced theory first. |
| first deliverable | D1: M6a after about 15 days; D2/D3: M5 | C3 M1: closure first, exact objects | M6.0 (exact 2-strand engine, exact closure data, no localization) right after M3a and CL1. |
| streaming elimination | D1: eliminate between consecutive terms | C1 I5: wrong if a pivot source lies in the newest term | With terms 0..k+1 present (positive crossing), only pivots whose *source* lies in terms <= k are eligible; mirrored for negative crossings. Unit test against full-then-eliminate on random complexes. |
| peak and memory model | uncolored calibration (2.3x) | C2 H2: 4-33x for colored braids | The decategorified peak model from M3a (B_p (x) W_k expansions of the LB complex). Colored c_e and LB/actual ratios are measured in M6a, and every colored envelope row is labelled [E] until then. |
| V15 (ordinary closure of folded data) | D2: "catches closure data" | C2 M1(a): both sides consume the same (G, F) | Relabelled: tests the closure *algorithm* only. Closure data are tested by single-type dim Phi_N(B_p) = MOY, generator degrees = c_diag, and end-to-end oracles. |
| V13 (cabling object lists = {p_max}) | D2: acceptance at M4 | C2 M3: rests on A1 [?] and omits the homological shift | Diagnostic only, until A1's precondition passes (every final cabling type is a maximal double-coset representative). Compare (p_max, h + H, q + Q). |
| canonical ladder | D1: min-flow; D2: first found by BFS | C2 L5 | Min-flow, with a deterministic tie-break. The canonical basis must not depend on the choice (test). |
| uncolored w0-type objects | D3: BS(rex x) | C3 L3: B_{w0} then needs the 6.7 GB idempotent | Thick (merge-split) ladders, as PLAN v4 says. |
| pivot coset mid-scan | D3: identity coset | C3 L4: none exists when mu != lambda | Any coset of the top stratum of q (one flow there). |
| component degree | D1: d = s - s'; R1 P2: s - s' + 1 | C1 I10 | e = t_src - t_tgt (D3 §2.3; checked: R{0} -> B_s{-1} by the start dot has e = 1). Test that every stored component has total q-degree 0. |
| object model | n uprights | C3 H1(d): EKLP Y_r passes through LR(r) with up to n·m parts | Expressions are ladders on any number of uprights with zero labels (`web::Ladder` with bound = K, so no truncation). `local.rs` already supports this. The Search provider may pass through more uprights. |

---------------------------------------------------------------------------------------------------

## 2. Architecture

### 2.1 Modules (new files under `src/`; existing validated modules unchanged except for making a few helpers `pub`)

| module | milestone | content |
|---|---|---|
| `schur.rs` | M3a | Hecke algebroid in the table model: tables, l_min, l_max, N(p), split/merge/rung rules, canonical basis, min-flow canonical ladders, Hom formula, multiplicities, c_diag, kappa, K_0 class, decategorified scan with LB, peak model and closure-weight histogram, decategorified closure (Euler characteristic) through `moy.rs` |
| `cclose.rs` | CL1 | ring R^lambda (weighted e-variables); `ClosureData` contract v2; per-block H*(Gr(a,N)) arithmetic (GrTable first, Schur/LR tables later); colored exact closure oracle; reduced mode; global_q; the (Delta_lambda G, F) adapter to `tower::close_tower` |
| `xengine.rs` | M6.0 | exact colored engine for small fibres: objects are 2-rung ladders, maps are exact `PMat` over R^lambda, decompositions by `hom_space` plus scalar intersection forms, exact closure data |
| `coset.rs`, `sexpr.rs`, `smap.rs`, `tmpl.rs` | M4a | cosets and flows at a coset; expressions with zero labels; elementary maps E1-E4 (closed forms, constant scalars); window templates (exact, certified); map arena with coset-local evaluation at stacked point vectors |
| `homprov.rs` | M4b | `HomProvider` (Exact, Search; Leaves later); stacked-point certification; genuineness certificates |
| `sindec.rs` | M5 | idempotents e_p = P(a); merge decompositions by adjunction and local intersection forms; per-coset U/V coordinates; caches |
| `sengine.rs` | M6a | typed complexes, colored Rickard insertion (half-steps), Markowitz scalar elimination with the streaming rule, multi-point closure-data preparation |
| `scoord.rs` | M6c / M7a | coordinates in certified Hom^e bases (final components first, then all), detecting sets, per-basis closure tables |
| `cinterp.rs` | CL2 | structured closure-data interpolation (fixed points, per-N boxes) |
| `ctower.rs` | CL3 | colored restriction tower (whole thick blocks, LR-Toeplitz maps, unit then block pivots) |
| `leaves.rs` | M4c (gated) | EKLP singular light leaves (type A specialization) |

### 2.2 Data flow

```
braid word + colors per bottom position (CLI, pd.rs)
  -> schur.rs: types, b_p, multiplicities, c_diag, K_0 checks              (N-independent, cached per lambda)
  -> scan (xengine for tiny cases; sengine otherwise): minimal complex of B_p{t} at h   (N-independent)
       per stage: K_0 class == decategorified scan; #objects >= LB
  -> closure data over R^lambda: generators (one-point determinant), Gram G_p, edge matrices F
       (exact in xengine; interpolated in sengine; per-basis tables in scoord)
  -> per N: colored exact oracle | tower on (Delta G, F) / gr(V) | ctower; reduced mode; global_q(N)
  -> Poincare polynomial (two primes for final results)
```

### 2.3 Fixed conventions for all new code

- **Labels and orientation.** Convention A of `crossing.rs` (PLAN "Conventions (fixed)").
  - Object q-degrees carry only N-independent shifts: -eps*k of the term (eps = +1 positive, -1 negative), plus decomposition shifts in EKLP normalization. Krhom ladders are EKLP singular Bott-Samelson bimodules with the same grading (R1 §0.1 [V]).
  - global_q(N) = - sum_crossings eps_c m_c (N - m_c) - sum_{unmarked blocks} lambda_i (N - lambda_i), with m_c = min of the two labels at c (R3 A4). For labels 1 this is `engine2::closure_q_shift`.
- **Component degree.** A component B_p{t} -> B_q{t'} lies in Hom^e(B_p, B_q) with e = t - t'.
- **N-independence.** Rickard complexes are built untruncated (`crossing_cx` with n >= c + d). Phi_N is assumed to kill webs with internal labels > N; this is tested (V2) before relying on it, and pruning stays an optional optimization. Colors > N give 0; N < max color is skipped.
- **Arithmetic.** F_p with p = 2^31-1. Final results also over F_4294967291. Every stacked-point decision is re-checked at fresh points.

### 2.4 The Search provider (C3 H1(c): specified as an algorithm)

**Input.** Source and target expressions X, Y (each with its idempotent), a degree d, and a target: either "full Hom^d" (dimension from the Hom formula) or "modulo lower terms at stratum q" (rank from the multiplicity [v^s] m_q). A budget: maximal path length L_max, candidate count C_max, wall time T_max, memory M_max.

**Search space.** Words of the double-leaf shape (EKLP Def 7.44):

    X --path--> Z --(dot in a segment ring)--> Z --path--> Y

- Z ranges over expressions (any number of uprights, zero labels allowed) with c_r(Z) != 0 for some stratum r shared by X and Y.
- A path is a sequence of local moves of known degree:
  - associativity and commutation, degree 0;
  - square switches EF <-> FE (composites of E2/E3), degree 0 on the top summand;
  - digon creation and annihilation (E2 unit/counit, degree +ab; E3 inclusion/trace, degree -ab);
  - fusion and fission of consecutive rungs, E^(r)E^(s) <-> E^(r+s), with Schur decorations;
  - dots (E1, degree 2j).
- Words are enumerated in a fixed order (path length, then dot degree, then lexicographic), with fixed seeds.

**Algorithm.** For each candidate:
1. Evaluate it at the stacked detecting set (one coset per stratum in S_d, D3 L3; M points from the formula of C1 I1).
2. Append it to an incremental echelon.
3. Stop on reaching the target rank.

On failure, report the missing (degree, stratum) pairs and return an error. Never accept a partial basis.

**Termination.** The budget guarantees it. **Completeness** is not guaranteed [?]. EKLP Thm 7.49 guarantees that double leaves exist in the singlestep calculus. Whether the moves above, restricted to bounded expressions, generate all Hom spaces is the open question of C3 §5.1 (a fullness result of the Mackaay-Stosic-Vaz type would settle part of it; no reference has been verified).

**Gate.** G1 (M4b), section 3.

---------------------------------------------------------------------------------------------------

## 3. Milestone plan

Estimates are agent-days [E]. The critiques found every design's schedule optimistic by 1.5-3x; these figures already include that correction. "Where" follows the resource guidelines: laptop means < 1 GB RSS and minutes, KRHOM_THREADS=2, nice 19, one process at a time; cloud means a 4-core, 15 GB box.

Mapping to earlier names: M3 = D1 M3, D2-M0, D3-M0. CL1 = D1's C1 (closure layer) and the closure half of D2-M5 and D3-M5; CL2 = D1's C2 (structured interpolation); CL3 = the colored tower of R3 A6. (CL is used for closure milestones so they are not confused with the critiques C1-C3.) M4a/M4b/M5 = D1 M4a/M4b/M5a/M5b, D2-M1..M3, D3-M1..M3. M6a = D1 M6a/M6b, D2-M4/M5, D3-M4. M7a = D1 M7a, D3's coordinate layer.

### 3.1 The milestones

**M3a. Decategorified oracle `schur.rs`** (FIRST TASK, §8). Laptop, 2-3 days.
- Acceptance:
  - canonical bases identical to `cb.pkl` (8 compositions);
  - Hom formula = `hom_space` on the 60 R1 cases;
  - (1^n) canonical basis = `hecke.rs` KL basis; c_r(b_p) = h_{r_max,p_max} for K <= 6;
  - Rickard rigidity for c, d <= 5;
  - K_0 lower bounds and peaks as tabulated (§8.5);
  - decategorified Euler characteristic = `chi` on 8 cases.

**M3b. Envelope measurements.** Cloud, 1 day of runs.
- Runs:
  - canonical bases of (3,3,3,3), (2,2,2,2,2), (5,5,5), (4,4,4) and the mixed compositions;
  - decategorified scans of the 3-colored Conway knot, of 2-colored 5-strand 12-crossing knots, of 2- and 3-colored scans of all knots in `data/knots_upto12.jsonl` with braid index <= 4, and of mixed colorings of links in `data/links_upto11.jsonl`.
- Output: a JSON table of LB, peak, kappa, max Gram weight, unknowns (unreduced and reduced) and types per case.
- Gate **G0**: fixes which rows of §1.3 enter the plan.

**CL1. Colored closure layer `cclose.rs`.** Laptop, 3-4 days. Runs in parallel with M3a; §9.1.
- Acceptance:
  - tau = `sbim::Tau`;
  - dim Phi_N of random small closed ladders from exact R^lambda-bases = `StateSpace` = MOY;
  - colored exact oracle = (Delta G) tower divided by gr(V);
  - on uncolored engine2 data, the new closure = `close_exact` (Conway N = 2..6);
  - reduced mode on single types.
- Gate **G4**: freeness of Delta_a R_N for all a <= 5, N <= 10 (the cases with a·N large run on cloud).

**M6.0. Exact small colored engine `xengine.rs`.** Laptop, 2-3 days. Needs M3a and CL1. §9.2.
- First colored homology through the new path: 2-strand braids with colors in {1,2}² (3 tried within a time budget), N <= 10, unreduced and reduced.
- Acceptance:
  - V1 unknots; V2 labels > N and wedge^N;
  - V4 Hopf (2,1), (2,2) bigraded for N = 2..6;
  - V5 Wang: 2-colored trefoil ranks 7, 22, 46, 79, 121, 172, 301 (N = 3..8, 10), T(2,5) 46 and 106 (N = 4, 5), T(2,7) 78 (N = 4), bigraded against `wang.py`;
  - V9 full cube (Hopf, kinks);
  - V8 cabling (2-colored T(2,n), n <= 5);
  - V3: (2,2)@N=3 trefoil = uncolored gl_3 trefoil;
  - χ = M3a's decategorified χ;
  - two primes.

**M4a. Localized substrate.** Laptop, 3 days.
- Content:
  - cosets and flows at a coset (backward search);
  - expressions with zero labels;
  - E1-E4 closed forms with constant scalars;
  - window templates from `hom_space`, with polynomiality certificates;
  - Rickard templates (colors <= 3 exact via `cube::local_cx` made `pub`);
  - map arena with stacked-point coset-local evaluation;
  - per-coset U/V coordinates.
- Acceptance:
  - flows per coset = singular Deodhar coefficient (280 R1 ladders);
  - coset blocks = `local::localize_map`;
  - E2/E3 identities (eps∘eta = Euler class, d∘iota = 0);
  - Rickard d^2 = 0, block-diagonality, uniqueness, equivariance f_{sigma c}(x) = f_c(sigma^{-1} x);
  - polynomiality certificate passes for all templates and fails on a deliberately perturbed one.

**M4b. Hom provider and certification.** Laptop, 2 days, plus a G1 spike of 3-5 days.
- Content: `Exact` and `Search` backends (§2.4); stacked certification (M from the oracle); fresh-point residuals.
- Acceptance:
  - certified dimension = Hom formula on the 60 R1 cases and on all type pairs of (2,2), (3,3), (2,1,2), (2,2,2) in degrees <= 4;
  - the cases of critique C1's faith table need exactly the predicted M.
- Gate **G1**: Search certifies (a) every merge-decomposition and closure-generator Hom needed for the (3,3) 2-strand blocks and all 21 types of (2,2,2), and (b) a fixed random sample of 50 (2,2,2,2) and 20 (3,3,3) pairs, within T_max = 10 min and 1 GB per pair. Record candidates tried versus target. If G1 fails: 2-strand colors <= 2 continue on the exact backend; M4c starts before the 3-strand part of M6a.

**M5. Singular idempotents and merge decompositions.** Laptop, 3 days.
- Content: half-step pipeline (split = re-index, W Prop 7.4.3(1) [L]); merge decomposition by adjunction Hom(B_q, B (x) M) = Hom(B_{q+}, B) (checked in critique C1) and local intersection forms at a top-stratum coset; e_p = P(a); orthogonalization; disk cache.
- Acceptance, for all types of (2,2), (3,3), (2,1,2), (2,2,2) and the half-step blocks used:
  - block ranks of e_p at one coset per stratum = h_{r,p}(1);
  - End^0 = 1;
  - pi∘iota = id, iota∘pi orthogonal idempotents;
  - ranks add up;
  - multiplicities = schur;
  - e^2 = e at fresh points;
  - (1,1,1,1): all 84 B_x B_s decompositions agree with `indec.rs` multiplicities.

**M6a. R-full localized engine on <= 3 strands.** Laptop, 3-4 days.
- Content:
  - typed complex; colored Rickard insertion per source object (C1 I9);
  - Markowitz scalar elimination with the streaming rule;
  - multi-point closure-data preparation (engine2 style; generators by the one-point determinant criterion; Gram from generator-only points; F by rescans; dense Interp <= 2e3 unknowns).
- Acceptance:
  - equals M6.0 on every overlapping case (Poincaré polynomials, and final multisets of (type, t, h));
  - K_0 class per stage = M3a;
  - V4 for all 13 R3 Hopf cases; V5 for k = 3 (N = 6, 7: 136, 307);
  - 2-colored 4_1, 5_2, 6_2, 6_3 and 8_19 at N = 2..6. N = 2 is the wedge^N rule (V2). Every result must also pass χ = `decat_chi`, V10 and two primes, and the cloud cube reference for 2-colored 4_1 at N = 3, 4 once it exists;
  - V3 at N = 3 for every 3-strand knot ≤ 10 crossings (2-colored gl_3 = uncolored gl_3 of the reversed word);
  - V10 Markov, conjugation and mirror on 3-strand words;
  - two primes; two point sets.
- Measured: colored c_e, the LB/actual ratio, the degree histogram of the final components, and peak RSS per point (all recorded in the envelope).
- Gate **G2**, at the end of M6a: the Search provider certified every Hom the M6a runs needed (no exact-backend fallbacks beyond fibre about 40), and it certifies a fixed sample of the (2,2,2,2) merge-decomposition and closure-generator Homs (the types of the 2-colored Conway LB complex) within budget. If G2 fails, start M4c before any 4-strand work.

**M6b. Validation sweep and integration.** Laptop, with table runs on cloud; 2-3 days.
- Content:
  - uncolored regression through the singular engine (V11);
  - CLI: colors per component, mixed colorings, reduced mode, a ↔ N-a option labelled "via duality";
  - batch runs on the knot and link tables for 2-colored 3-strand inputs and mixed colorings.
- Acceptance:
  - engine2's ranks: Conway 34..426, reduced 33 then 49, the regression set;
  - V6 Euler characteristics (transfer-matrix MOY; so(6) Kauffman for wedge^2 gl_4 on knots);
  - V7 (marked-block independence; [N choose a] P_red - P_unred = (1+t) Q with Q >= 0).

**M6c. M-point scan with final-component coordinates.** Laptop, 3 days.
- Content: certified bases of the *final* (p, q, e) spaces (Search); coordinates by stacked solve with residual; per-basis closure tables from identity-coset samples (no rescans).
- Acceptance: identical results to M6a on all of M6a's cases; coordinates consistent across M + 2 points; scan count per braid drops from #F-monomials to M.
- Gate **G3**: the list of (p, q, e) spaces that real scans touch (C3 H1(b)) and the per-component cost of R-full versus a coordinate prototype on the 2-colored trefoil and figure-eight decide M7a.

**CL2. Structured closure interpolation.** Spike 2 days (start before M6c needs weights > 12); full 3-5 days.
- Content: Grassmannian fixed-point nodes per block with disjoint node sequences (all roots distinct); per-N box restriction; per-block triangular transforms on the lower set; oversampling with residuals; identity-block invertibility checked at every point (C1 I8).
- Acceptance: random weighted-homogeneous symmetric polynomials of weight <= 22 on (2,2,2,2) recovered exactly; equality with dense Interp up to 2e3 unknowns; measured cost per weight-17 and weight-22 entry.

**CL3. Colored tower `ctower.rs`.** 5-7 days. Laptop development, cloud for large N.
- Content: R3 A5-A6 (one thick block per stage, Schur/LR coordinates, LR-Toeplitz maps); unit pivots, then M1b.2 block pivots (archived `fix2_partial/`).
- Acceptance: equals the colored exact oracle and the (Delta G) tower on every CL1/M6 case; stage statistics recorded (open question O4 for colors).

**M7a. R-coord layer.** 5-8 days, only if G3 passes.
- Content: D3 K1-K5 (all components as coordinates; transport, diff and composition tables at stacked detecting sets; intrinsic words only).
- Acceptance: identical to R-full on all M6 cases; memory and time measured.

**M4c. EKLP singular light leaves** (only if G1 or G2 fails). 10-20 days.
- Content: type A dictionary (D1 §2.10), checked on EKLP Examples 2.5 and 2.8; cores; rex-move search; sinister ELLs evaluated through closed forms. A gate at the end of the combinatorics part: path counts = characters (EKLP Cor 3.7) and Hom counts (Cor 4.17).
- Acceptance: double leaves certified (count, stacked rank, polynomiality) on random expression pairs on 3-4 uprights.

**M7b. Envelope runs.** Cloud, ongoing.
- Runs: 2-colored tables (<= 3 strands; 4 strands where G0/G3 allow); 3-colored 3-strand; 4-5-colored 2-strand (after (5,5) Rickard certificates, CL2 and CL3); 2-colored Conway reduced at N <= 6, then higher.
- Output: results with timings, peak RSS and the "computed / via duality / predicted" status per entry.

**M8 (later). Deformation and concordance invariants.** Rose-Wedrich / Wu generic deformations: the same minimal complexes; only the closure ring changes (equivariant A_lambda; `symfun::GrTable` with c != 0).
- Uncolored s_N first, against the 2977 s values in `data/knots_upto12.jsonl`; then colored.
- Note [D]: evaluating closure data at Grassmannian fixed points with *generic* equivariant parameters computes the generically deformed homology directly (rank at the fixed points). The CL2 machinery therefore serves this goal too.

### 3.2 Critical path and first results [E]

```
M3a ─┬─> M3b (cloud)                                   ─> G0
     └─> M6.0 <── CL1 (parallel with M3a)               ─> first colored results via the new path (~1.5 weeks)
M3a + CL1 ─> M4a ─> M4b (+G1 spike) ─> M5 ─> M6a        ─> first new colored data: 2-colored 3-strand, 3-colored 2-strand (~5-6 weeks)
                                         ─> M6b ─> M6c (+CL2 spike) ─> G3 ─> M7a, CL3 ─> 2-colored Conway reduced, N <= 6 (>= 9-11 weeks, uncertain)
G1 or G2 fails ─> M4c (adds 2-4 weeks)
```

---------------------------------------------------------------------------------------------------

## 4. Validation strategy

### 4.1 Layers (each layer runs automatically from the milestone that introduces it)

| layer | checks | catches | first milestone |
|---|---|---|---|
| L0 decategorified | schur vs the Python references (`cb.pkl`, `homres.txt`) and `hecke.rs` (K <= 6); Rickard rigidity; K_0 class of *every* stage complex equals the decategorified scan; #objects >= LB; decategorified Euler characteristic vs `chi` | KL and normalization errors; lost, duplicated or mis-shifted summands; v ↔ q orientation | M3a |
| L1 type-level certificates | Hom dimensions by stacked rank = Hom formula; block ranks of e_p = KL values at 1; End^0 = 1; pi∘iota = id; e^2 = e at fresh points; polynomiality of templates; d^2 = 0; fresh-point residuals of every solve | wrong bases, non-genuine maps, point-dependent normalizations | M4a-M5 |
| L2 closure | tau = `Tau`; single-type dim Phi_N(B_p) = MOY of the closure of b_p (and `StateSpace` for small ladders); generator degrees = c_diag(b_p); exact oracle = (Delta G) tower divided by gr(V) = ctower; uncolored special case = `close_exact` | Grassmannian arithmetic, normalization, Nakayama over R^lambda, tower bugs | CL1 |
| L3 end-to-end small | V1, V2, V4 (closed Hopf formula, 13 cases), V5 (Wang T(2,2m+1), bigraded), V9 full cube (Hopf, kinks), V8 cabling (where the shift is verified: m <= 2 and (3,1)), M6.0 = M6a on the overlap | conventions, differentials for mixed labels, h-grading | M6.0 |
| L4 multi-strand, rank level | **V3 duality at N = a+1** against engine2 on the reversed word (primary); V10 Markov, braid relations, conjugation, mirror; V7 reduced checks; V11 uncolored regression through the singular engine; one-off cloud references (cube for 2-colored 4_1 at N = 3, 4; cabling for 2-colored 4_1, 6 cable strands, if memory allows) | engine bugs invisible on 2 strands; reduced-theory errors | M6a |
| L5 sweeps | V6 Euler characteristics over the tables (transfer-matrix MOY; so(6) Kauffman for wedge^2 gl_4 on knots; Knot Atlas A3/A4); Lewark-type relation | missing terms, grading errors at scale (blind to differentials, C3 H3) | M6b |

### 4.2 Standing rules

- Every final result is computed over two primes. Every stacked-point decision is re-verified at points not used to make it.
- Every scan asserts the K_0 class per stage. It is cheap, and it is the strongest scan-side invariant.
- Euler-characteristic agreement is never accepted as evidence about differentials. Rank-level oracles (duality, cube, cabling, Wang, Hopf) are required before a new color pattern counts as validated.
- Duality-derived numbers are reference data, labelled "via duality (expected)". A disagreement is investigated on both sides.
- Cloud reference runs record the exact braid word, colors, N, mirror status and commit hash. R3 found that the stored cube result was the mirror input: `trefoil_2colored_gl3_cube.txt` is σ1^{-3}.

---------------------------------------------------------------------------------------------------

## 5. Risks and fallbacks

| # | risk | effect | early signal | mitigation and fallback |
|---|---|---|---|---|
| K1 | Search provider fails to span (heuristic moves miss maps), or is too slow on large types | no decompositions or closure generators beyond small fibres | G1 failure rate; candidates per certified dimension | exact backend for small fibres; adjunction keeps Homs inside one block and low degree; M4c EKLP light leaves (Thm 7.49 guarantees existence); the envelope shrinks to what certifies |
| K2 | Closure data too large (1.6e5 unknowns per entry for 2-colored Conway; 1e7-1e9 for (2^5), (3^4)) | larger targets blocked regardless of the scan | M3b Gram-weight histograms; CL2 spike timings | reduced theory first; per-N box restriction (small N cheap); per-basis tables cached across braids (M6c/M7a); research question: closed forms for Gram matrices of large types (D1 §9.3) |
| K3 | R-full memory (7-10 GB per point for 2-colored Conway) and rescans | 4-strand 2-colored out of reach for R-full | M6a measured c_e and peaks | M-point scans (M6c), R-coord (M7a), streaming per source object, u32 storage |
| K4 | Non-genuine closed-form maps (scalars [D], E4 identity [D], composite Rickard words) | silently wrong homology | polynomiality certificate failures; d^2 != 0 | certificates are mandatory (C3 H3); compare with exact `hom_space` on small windows |
| K5 | Point dependence (per-point normalizations, seed-dependent choices) | inconsistent F, silent errors | stacked-point and fresh-point residuals; e^2 = e checks | constant scalars only; seeded deterministic choices; abort on disagreement |
| K6 | Phi_N does not kill internal labels > N for some colored web | wrong results in the N-independent scan | V2 tests; duality at N = a+1 (internal labels 2a > a+1) | prune by `crossing_cx(n = N)` per N (costs one scan per N) |
| K7 | Categorified duality false in some case | oracle mismatch | V3 disagreement while all other oracles agree | treat as a research finding; keep the cube and cabling references |
| K8 | Characteristic p effects (p-canonical basis at p = 2^31-1) | wrong multiplicities | multiplicity asserts; two-prime disagreement | two primes; the K_0 asserts catch multiplicity changes |
| K9 | Colored tower blow-up at large N (non-split blocks times C(N,a)) | closure for N >= 8 slow | ctower stage statistics | reduced mode; block pivots; exact oracle and (Delta G) tower for small cases |
| K10 | Schedule (all designs underestimated) | late results | milestone slippage | M6.0 delivers validated colored results early; gates prevent sinking weeks into a dead end |
| K11 | Python BFS for canonical bases does not scale to (3^4), (2^5), (5^3) | no envelope numbers for stretch targets | M3a/M3b timings | run on cloud; these runs are measurements, not gates (C1 I11). If the BFS is too large, a later variant builds b_p block by block from canonical parents (b_p' * rung, then peel), keeping only the blocks a scan needs |

---------------------------------------------------------------------------------------------------

## 6. Laptop versus cloud

The rule: anything predicted to exceed about 1 GB RSS or about 10 minutes of CPU goes to a cloud session (4 cores, 15 GB), never to the laptop. Cloud sessions cannot reply by message, so each cloud job writes its results, exact inputs, commit hash, timings and `/usr/bin/time -v` output to files (or a branch) that the next local session reads.

| task | where | why |
|---|---|---|
| All implementation, unit tests, default acceptance suite (kept under about 60 s extra) | laptop | small, interactive |
| M3a heavy tests: (2,2,2,2) and (3,3,3) canonical bases, 2-colored Conway decategorified scan with peaks | laptop, `#[ignore]` tests, one at a time | expected < 1 min and < 0.5 GB in Rust (Python: 366 s, 0.97 GB) |
| **M3b** canonical bases (3,3,3,3), (2,2,2,2,2), (5,5,5), (4,4,4); scans of 3-colored Conway, 2-colored 5-strand knots, all knots with braid index <= 4 colored 2 and 3, and mixed-color links | **cloud** | (3^4) has about 1.6e5 clean BFS states and (2^5) about 7.6e5 (estimate: several GB); table sweeps are embarrassingly parallel |
| G4 freeness check of Delta_a R_N over H*(Gr(a,N)) for a·N large (e.g. (4,9), (5,10): dim up to 3.0e4 times A-structure) | cloud (small cases on laptop) | linear algebra on spaces up to 1e5 |
| Duality reference tables: engine2 on reversed words of all knots in `data/knots_upto12.jsonl` at N = 3..6 (wedge^{N-1} of a knot = uncolored gl_N of the reversed knot; the same identity also relates other pairs, e.g. 3-colored gl_5 = 2-colored gl_5, used later) | cloud, added to the existing batch | thousands of runs; 5-strand words take minutes each |
| Full-cube references: 2-colored 4_1 at N = 3 and 4 (81 vertices); Hopf (3,1), (3,2) at N <= 5; (2,1,1) links at N = 3, 4 | cloud, once, with a memory watchdog | the 2-colored trefoil at gl_3 took 42 min on the laptop |
| Cabling references: 2-colored T(2,n), n <= 7, N = 2..6; 2-colored 4_1 (6 cable strands, 16 cable crossings) | cloud | 6-strand ordinary scans once reached 77 GB; needs a cap |
| M6a/M6b table runs (2-colored 3-strand, all N <= 10; mixed colorings) | cloud | hundreds of knots times multi-point scans |
| CL2 spikes on weight 17-22 entries; closure-data libraries for (2,2,2,2) and (3,3,3) | cloud | 1e5-1e6 point evaluations |
| CL3 / M7b large-N closures (N >= 8 unreduced), 3-colored N = 10, 2-colored Conway | cloud | GB-scale towers |

---------------------------------------------------------------------------------------------------

## 7. Items the critiques left unresolved (explicit)

| # | item | status | where it is settled |
|---|---|---|---|
| U1 | Spanning of the Search provider; whether bounded ladder moves (on any number of uprights) generate all Hom spaces (C3 §5.1); completeness of light leaves with arbitrary degree-0 maps in place of rex moves (D1 Q6, C1 I3) | open [?] | G1, G2 (empirical); M4c if needed |
| U2 | A fast closure-data solve at weights 17-24 (Newton-type transforms on lower sets of multipartitions); invertibility of the identity-coset block at structured points (C1 I8) | open; per-block triangular transform proposed here [D] | CL2 spike |
| U3 | Freeness of Delta_a R_N over H*(Gr(a,N)) in general (expected Delta_a R_N = H*(Fl(1..a; N))); verified only for (2, 2..6), (2,10), (3, 3..6), (3,8), (4, 5..7), (5,6) (C2) | open [?] | G4 (all a <= 5, N <= 10) before production use |
| U4 | Phi_N kills webs with internal labels > N, for colored labels (3 probes only) | open [?] | V2 in M6.0; duality at N = a+1 |
| U5 | Categorified duality a ↔ N-a (Queffelec-Rose expected; one rank-level case, the trefoil at (2,2), N = 3) | expected [?] | used as an oracle; mismatches investigated |
| U6 | Reduced colored theory: marked-point compatibility of the KRW identification (R3 C1); marking different blocks of one component gives the same result | expected [?] | V7 in M6b |
| U7 | Projector absorption A1 (C2): affects only the diagnostic V13 | open [?] | V13 stays diagnostic |
| U8 | Braid and Markov invariance of the SSBim colored Rickard complexes: the Webster-Williamson theorem number is unchecked (ETW Thm 3.5 covers Foam_N) | literature check pending [?] | V10 in every milestone from M6a |
| U9 | Lemma L3's [D] step (restriction to a stratum is injective on Hom modulo lower terms); supported on 40 two-upright cases (C3 E5) | [D] | detecting sets are certified anyway; T4 test in M4b |
| U10 | Colored degree spreads of differential components (they set M and the R-coord sizes); colored c_e; LB/actual ratio; peak reduction from eager elimination | unmeasured | M6a measurements |
| U11 | Cache reuse across braids for 282 and 2 008 types (C3 Q4) | unmeasured | G3 |
| U12 | O4 for colors: are colored tower stage sizes N-independent after the first restriction? | unmeasured | CL3 |
| U13 | Is every 2-upright canonical basis element exactly the character of a 2-rung ladder (so e_p = id on 2 strands)? (D1 Q4) | expected from canonical bases of U-dot(sl_2) | M3a test T-h2 |
| U14 | Bruhat order on tables via the rank criterion (R1 [?]) | avoided: only l_max is used as a linear extension | - |
| U15 | Characteristic p = 2^31-1 versus the p-canonical basis | assumed harmless for K <= 20 [?] | two primes |
| U16 | Cabling shift t^T q^{-(NT+m(m-1))} verified only for m <= 2 and (3,1) | partially verified | cabling comparisons restricted to those cases |
| U17 | Wang arXiv:2211.08409 Prop 1.6 / Thm 3.5 exponent: the extracted text is 2k too high; the corrected form matches MOY (R3) | journal version unchecked | the Hopf formula of R3 is used, not the printed one |
| U18 | Kauffman/so(6) recipe for links (linking-number factor) | open | V6 on knots only until calibrated |
| U19 | Whether high-Gram-weight types actually occur in minimal complexes of 4-5-strand knots (if not, K2 softens) | measurable | M3b (decategorified LB types) |
| U20 | D2's claim that the singular engine also fixes 6-7-strand uncolored data preparation (thick w0 ladders) | unsupported (C2 L2) | measure in M6b (V11 on the 6-strand braids) |
| U21 | Research side question from the previous message: do colored s-type invariants carry information beyond the uncolored ones? | open | after M8, on 2- and 3-strand knots |

---------------------------------------------------------------------------------------------------

## 8. FIRST IMPLEMENTATION TASK: M3a, the decategorified oracle `schur.rs`

### 8.1 Goal and scope

**Goal.** Port the validated Python Hecke-algebroid reference to Rust and extend it with what every later milestone needs. The reference consists of:
- `scripts/schur.py`: the table model, split/merge/rung rules, canonical basis, Hom formula;
- `scripts/t8.py`: the decategorified scan and the Grothendieck lower bound;
- `proto_D1/peak.py`: the peak model;
- parts of `proto_D2/transfer.py`: the cross-checks against ordinary KL polynomials.

**Why first.**
- Every design puts it first (D1 M3, D2-M0, D3-M0).
- Everything later consumes it: types, canonical bases, multiplicities for decompositions, Hom dimensions for certification, closure generator degrees, the K_0 class checked after every crossing, and the envelope measurements behind gate G0.
- It is pure combinatorics with an existing validated reference, so it is low risk.
- It pins down the v ↔ q and shift conventions early, through the decategorified Euler characteristic test T-n against `chi`.

**Non-goals.** No bimodules, no localization, no homology. The Bruhat order on tables via a rank criterion is not implemented (R1 marks it [?]): use l_max only. No Foam_N pruning.

**Repository permissions for this task.**
- Add the new files listed in §8.2, plus one line `pub mod schur;` in `src/lib.rs`. Change nothing else.
- The working tree has local, uncommitted edits (`PLAN.md`, `src/bin/krhom.rs`, `src/colored.rs`, `tests/acceptance.rs`, `src/bin/colored_cable.rs`). Leave them untouched: no stash, revert or reformat.
- Work on a new branch from `dev` (e.g. `colored-m3a`); commit only as the orchestrator instructs.
- Resource rules:
  - `export PATH="$HOME/.cargo/bin:$PATH"`;
  - `KRHOM_THREADS=2 CARGO_BUILD_JOBS=2 nice -n 19`, one process at a time, < 1 GB RSS;
  - `CARGO_TARGET_DIR=/private/tmp/claude-501/-Users-existair-Documents-KRcalculator/bf6f8870-08ad-436b-9ff7-0414655cc59a/scratchpad/colored_design/target` (never the repo's `target/`).

### 8.2 Files

| file | content |
|---|---|
| `src/schur.rs` | the library (§8.3, §8.4) |
| `src/lib.rs` | add `pub mod schur;` |
| `src/bin/schur_scan.rs` | CLI (§8.6) |
| `tests/schur.rs` | integration tests (§8.5) |
| `tests/data/schur/cb_ref.txt` | canonical bases of (1,1,1), (1,1,1,1), (2,2), (3,3), (2,2,2), (3,3,3), (2,1,2), (2,2,2,2), exported from `cb.pkl`: 19 088 entries for (2,2,2,2), about 1 MB in total. If the size is unwelcome in the repo, keep (2,2,2,2) in the scratchpad and skip that test when the file is absent |
| `tests/data/schur/hom60.txt` | `scripts/homcases.txt` and `scripts/homres.txt` merged (60 ladder pairs, dims per normalized degree) |
| `tests/data/schur/chi_ref.txt` | the Euler characteristics of Appendix A.3 |
| `colored_design/scripts/export_schur_fixtures.py` (scratchpad, not repo) | writes `cb_ref.txt` from `cb.pkl` and `hom60.txt` from `homres.txt`, whose field 6 already holds the normalized-degree dimensions that `cmp_hom.py` compares with `hom_gd`. `chi_ref.txt` is transcribed from Appendix A.3, or regenerated with the `chi` binary. Run with `venv/bin/python -I` |

Fixture line formats (all whitespace-separated; tables written row-major as `r0c0,r0c1/r1c0,r1c1`; Laurent polynomials as `exp:coef` lists):

```
CB <lam>|<p>|<q> <e>:<c> <e>:<c> ...|<q> ...           one line per canonical basis element b_p
LAD <lam>|<p>|<i><E|F><r> ...                            the ladder found by schur.py's BFS (informational)
HOM <bottom>|<src rungs>|<tgt rungs>|<dlo>|<dhi>|<d>:<dim> ...   rungs as in homcases.txt ("E0:1 F0:1")
CHI <N>|<labels>|<word>|<e>:<c> ...
```

### 8.3 Types and functions (signatures are binding; internals free)

```rust
//! src/schur.rs — Hecke algebroid (q-Schur algebroid) of type A in the contingency-table model.
//! Conventions: EKLP arXiv:2401.03053 §3 (eqs. (44), (45)), Williamson arXiv:1010.1283v3 §2.3,
//! R1_singular_theory.md §1 and Appendix B. Rows = bottom blocks (lambda), columns = top blocks (mu).
use crate::laurent::LPoly;
use crate::web::{Dir, Rung};

pub type Comp = Vec<u8>;                     // composition, zeros allowed, K = sum <= 64

#[derive(Clone, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub struct Table { pub nr: u8, pub nc: u8, pub e: Vec<u8> }     // row-major
impl Table {
    pub fn identity(lam: &[u8]) -> Table;    // diag(lambda)
    pub fn get(&self, i: usize, j: usize) -> u8;
    pub fn rows(&self) -> Comp;              // = lambda
    pub fn cols(&self) -> Comp;              // = mu
    pub fn lmin(&self) -> i64;               // sum_{i>i', j<j'} M_ij M_i'j'
    pub fn lmax(&self) -> i64;               // lmin + l(w_lambda) + l(w_mu) - sum C(M_ij, 2)
    pub fn ncosets(&self) -> u128;           // N(p) = prod_i lambda_i! / prod_j M_ij!
    pub fn pi_lr(&self) -> LPoly;            // prod_ij [M_ij]!  (symmetric quantum factorials)
    pub fn parse(s: &str) -> Table;          // "2,0/0,2"
    pub fn fmt(&self) -> String;
}
pub fn tables(lam: &[u8], mu: &[u8]) -> Vec<Table>;      // all tables with these margins, sorted by (lmax, e)

/// Dense Laurent polynomial, i128 coefficients, checked arithmetic (unsimplified counts reach 1e11+).
#[derive(Clone, PartialEq, Eq, Hash, Debug, Default)]
pub struct Lp { pub lo: i32, pub c: Vec<i128> }
impl Lp {
    pub fn zero() -> Lp;  pub fn one() -> Lp;  pub fn mono(e: i32, k: i128) -> Lp;
    pub fn is_zero(&self) -> bool;  pub fn at1(&self) -> i128;  pub fn bar(&self) -> Lp;
    pub fn shift(&self, s: i32) -> Lp;  pub fn mul(&self, o: &Lp) -> Lp;
    pub fn add_scaled(&mut self, o: &Lp, k: i128, shift: i32);   // self += k * v^shift * o
    pub fn sym_part(&self) -> Lp;            // unique bar-invariant m with self - m in vZ[v] (schur.py sym_part)
    pub fn max_exp(&self) -> Option<i32>;  pub fn min_exp(&self) -> Option<i32>;
    pub fn to_lpoly(&self) -> LPoly;  pub fn from_lpoly(p: &LPoly) -> Lp;
}

pub type Tid = u32;                          // interned table id (per Schur instance)

/// Element of the block (lambda, mu) in the standard basis h_p.
#[derive(Clone, Debug, PartialEq)]
pub struct Elt { pub mu: Comp, pub t: Vec<(Tid, Lp)> }      // sorted by Tid, no zero coefficients

pub struct CanonicalBasis {
    pub mu: Comp,
    pub order: Vec<Tid>,                      // all tables of (lambda, mu) by ascending lmax
    pub b: std::collections::HashMap<Tid, Elt>,              // b_p, with coefficient 1 at p
    pub ladder_bfs: std::collections::HashMap<Tid, Vec<Rung>>,     // first found (schur.py order)
    pub ladder_minflow: std::collections::HashMap<Tid, Vec<Rung>>, // min flows, then length, then lexicographic
}

pub struct Schur { pub lam: Comp, /* interner, rung-map cache, bases cache: private */ }
impl Schur {
    pub fn new(lam: &[u8]) -> Schur;
    pub fn intern(&mut self, t: Table) -> Tid;
    pub fn table(&self, id: Tid) -> &Table;
    pub fn unit(&mut self) -> Elt;                                       // h_{diag(lambda)}
    pub fn split_col(&mut self, x: &Elt, j: usize, parts: &[u8]) -> Elt; // EKLP (44)
    pub fn merge_cols(&mut self, x: &Elt, j: usize) -> Elt;              // EKLP (45)
    pub fn rung(&mut self, x: &Elt, r: &Rung) -> Elt;                    // E_i^(r): split col i+1 as (r, b-r), merge cols i, i+1;
                                                                         // F_i^(r): split col i as (a-r, r), merge cols i+1, i+2
    pub fn mul_rungs(&mut self, x: &Elt, rs: &[Rung]) -> Elt;
    pub fn ladder_char(&mut self, rs: &[Rung]) -> Elt;                   // unit * rungs
    pub fn canonical_basis(&mut self, mu: &[u8]) -> &CanonicalBasis;     // lazy, cached
    pub fn expand(&mut self, x: &Elt) -> Vec<(Tid, Lp)>;                 // x = sum m_q b_q; panics unless all m_q in N[v,v^-1]
    pub fn hom_gd(&self, x: &Elt, y: &Elt, dlo: i32, dhi: i32) -> Vec<i128>; // EKLP Thm 4.16 (y self-dual), degrees dlo..=dhi
    pub fn rank(&mut self, mu: &[u8], p: Tid) -> u128;                   // rank over R^lambda = sum_q c_q(b_p)(1) N(q)
    pub fn c_diag(&mut self, p: Tid) -> Lp;                              // coefficient of h_diag in b_p (mu = lambda)
    pub fn kappa(&mut self, mu: &[u8]) -> i128;                          // max_{p,q} c_q(b_p)(1)
    pub fn s_storage(&mut self, mu: &[u8], p: Tid, q: Tid) -> u128;      // S(p,q) = sum_r c_r(b_p)(1) c_r(b_q)(1) N(r)
    pub fn term_mult(&mut self, p: Tid, rs: &[Rung]) -> Vec<(Tid, Lp)>;  // expand(b_p * [rungs])
    pub fn save_basis(&self, mu: &[u8], path: &std::path::Path) -> std::io::Result<()>; // header "schur-cb v1 lam=.. mu=.."
    pub fn load_basis(&mut self, path: &std::path::Path) -> std::io::Result<()>;
}

pub struct Term { pub rungs: Vec<Rung>, pub hdeg: i32, pub dq: i32 }
pub fn crossing_terms(c: u8, d: u8, positive: bool, i: usize) -> Vec<Term>; // = t8.py cross_terms; zero rungs dropped
pub fn flows(bottom: &[u8], rs: &[Rung]) -> u128;                       // prod over rungs of C(giver, r)
pub fn emono_count(lam: &[u8], w: usize) -> u128;                       // e-monomials of weight w in R^lambda (weights j)

pub struct ScanOpts { pub peak: bool, pub every_stage: bool, pub marked_block: Option<usize> }
pub struct StageStats {
    pub step: usize, pub labels: Comp,
    pub unsimplified_objects: u128, pub unsimplified_ldim: u128,
    pub lb_objects: u128, pub lb_types: usize, pub lb_ldim: u128, pub lb_sum_s: u128,
    pub max_rank_present: u128, pub kappa_present: i128,
    pub peak_before_next: Option<(u128, u128, u128)>,   // (objects, ldim, sum S) of LB (x) next crossing, expanded
}
pub struct GramRow { pub p: Tid, pub weight: i32, pub unknowns: u128, pub unknowns_reduced: u128 }
pub struct ScanResult {
    pub stages: Vec<StageStats>,
    pub k0: Vec<(Tid, i32, i128)>,     // signed multiset (type, q-shift, sum_h (-1)^h mult): the K_0 class
    pub gram: Vec<GramRow>,            // for the final LB types (block (lambda, lambda))
}
pub fn decat_scan(s: &mut Schur, word: &[i32], opts: &ScanOpts) -> ScanResult;
/// Decategorified closure: chi(KhR_N) = q^{-sum_c eps_c m_c (N - m_c)} * sum_{(p,k)} signed * q^k * tr_N(b_p),
/// tr_N(b_p) = moy::closure_value(L_p) - sum_{q != p} m_q(q) tr_N(b_q), char(L_p) = sum m_q b_q,
/// with closure_value := 0 when L_p has an intermediate label > N.
pub fn decat_chi(s: &mut Schur, res: &ScanResult, word: &[i32], n_gl: usize) -> LPoly;
```

`Rung` and `Dir` are the existing `crate::web` types: `Rung { i, dir, r }`, where `E` moves r from upright i+1 to i and `F` moves r from i to i+1. This is the same convention as `crossing.rs` and `schur.py`.

### 8.4 Algorithms (follow the Python reference unless stated)

1. **Tables, l_min, l_max, N(p), pi_LR.** As in `schur.py`.
2. **Split and merge.** As `split_col` and `merge_cols` in `schur.py`.
   - Merge coefficient: v^{lmin(new) - lmin(old)} * pi_LR(new) / pi_LR(old) (exact division, a product of quantum binomials).
   - Split coefficient: v^{lmax(old) - lmax(new)}.
   - Cache the sparse "rung map" (source table to a list of (target table, Lp)) per (source mu, rung). Scans apply the same rungs many times.
3. **Canonical basis.** Port `canonical_basis(lam, mu)`:
   - BFS from 1_lambda over rungs;
   - extend only *clean* states (a unique l_max-maximal table with coefficient exactly 1);
   - deduplicate on (top composition, sorted set of leading tables);
   - record the first ladder found for each table of the target block;
   - peel in ascending l_max, subtracting sym_part(x[q]) * b_q for q in descending l_max;
   - assert that every remaining coefficient lies in vZ[v].

   Peeling uses the canonical bases of the *target* block only. Store characters compactly (interned tables, `Lp`). Free BFS states that can no longer be extended.
4. **Minimal-flow ladders.** A uniform-cost search over the same clean states, ordered by (flows, length, lexicographic rung list); the first pop of each table gives `ladder_minflow`.
   - The canonical basis must come out identical whichever ladder is peeled (test T-c2). b_p is unique: bar-invariant, h_p + vZ[v]-span of lower terms.
5. **Hom formula.** gd Hom(x, y) = sum_M x[M] y[M] gd(R_M), with gd(R_M) = prod_ij prod_{k=1..M_ij} 1/(1 - v^{2k}), truncated to [dlo, dhi]. Port `hom_gd`.
6. **Expansion into the canonical basis.** Port `expand_cb` (t8.py): repeatedly take the l_max-maximal table and subtract its coefficient times b_top. Assert non-negative integer coefficients.
7. **Decategorified scan.** Port `scan` (t8.py) with `crossing_terms`:
   - per homological degree h an `Elt`; apply each term's rungs and the shift dq;
   - labels swap after each crossing;
   - stage statistics in the block (lambda, current labels), using that block's canonical basis (lazy; mixed colors need blocks (lambda, sigma lambda));
   - LB = sum over (p, q-shift) of |sum_h (-1)^h coefficient|;
   - the peak as in `peak.py`, with the term multiplicities cached per (p, rungs).
8. **Gram rows.** For the final LB types:
   - weight w_p = max exponent of c_diag(b_p). It must equal l(p_max) - l(w_lambda) (assert; D2 uses the latter);
   - `unknowns` = emono_count(lambda, w_p);
   - `unknowns_reduced` = emono_count(lambda with the marked block removed, w_p).
9. **decat_chi.** As in the doc comment of `decat_chi`.
   - `moy::closure_value` is evaluated on `web::Ladder { n: N, bottom: lambda, rungs: L_p }`, and only when every level satisfies `Rung::apply(.., N)` (labels <= N). Otherwise the value is 0 (skew Howe: that weight space is zero).
   - The trace is linear and v-semilinear with v ↦ q: tr(v^n x) = q^n tr(x).
   - **If T-n fails**, first test the alternative identification v ↦ q^{-1} and report which one passes. Do not silently change conventions.
10. **Performance targets.** Laptop, one thread:
    - (2,2,2,2) canonical basis < 60 s and < 500 MB (Python: 366 s, 0.97 GB);
    - (3,3,3) < 30 s;
    - 2-colored Conway scan with peaks < 60 s and < 500 MB (Python: 28-36 s for the scan alone);
    - all non-ignored tests < 60 s in total.

### 8.5 Tests and expected outputs (`tests/schur.rs`; heavy ones `#[ignore]`, run once with `--ignored` and reported)

**T-a. Table counts.**

| lambda = mu | count |
|---|---|
| (1,1,1) | 6 |
| (1,1,1,1) | 24 |
| (2,2) | 3 |
| (3,3) | 4 |
| (4,4) | 5 |
| (5,5) | 6 |
| (2,1,2) | 11 |
| (2,2,2) | 21 |
| (3,3,3) | 55 |
| (2,2,2,2) | 282 |
| (5,5,5) | 231 |
| (3,3,3,3) | 2 008 |
| (2,2,2,2,2) | 6 210 |
| (4,4,4,4) | 10 147 |
| (5,5,5,5) | 40 176 (enumeration only) |

**T-b. Worked examples (R1 §1.6).**
- (2,2), with p0 = [[2,0],[0,2]], p1 = [[1,1],[1,1]], p2 = [[0,2],[2,0]]:
  - (l_min, l_max) = (0,2), (1,5), (4,6); N(p) = 1, 4, 1;
  - b_p1 = h_p1 + (v + v^3) h_p0;
  - b_p2 = h_p2 + v h_p1 + v^4 h_p0;
  - ranks 1, 6, 6.
- (2,1) → (1,2): b_c1 = h_c1 + v h_c0, with c0 = [[1,1],[0,1]] and c1 = [[0,2],[1,0]]; ranks 2, 3.
- (2,1) → (2,1): b_d1 = h_d1 + v^2 h_d0, with d1 = [[1,1],[1,0]]; rank 3.
- Hom dimensions of (2,2), degrees -2..4: the 3x3 table of R1 §1.6 (e.g. Hom(B_p1, B_p1) = 0 0 1 0 5 0 14).

**T-c. Reference canonical bases.** Every b_p of the 8 compositions in `cb_ref.txt` is identical.
- Aggregates to report: sum of ranks 19, 219, 13, 63, 1 045, 76 399, 199, 397 687; kappa 1, 2, 2, 3, 8, 54, 4, 63 (order as in §8.2).
- **T-c2.** Canonical basis from `ladder_minflow` = from `ladder_bfs`.
- Report the sum over p of flows(L_p^minflow). Expected values are at most 1 171 for (2,2,2), 105 787 for (3,3,3) and 519 337 for (2,2,2,2), the sums for R1's ladders. The search deduplicates states, so it is not provably optimal: a larger value is reported, not counted as a failure.
- Report flows(L_p^minflow) for the rank-25 920 table (R1's ladder has 34 560).

**T-d. Hom formula.** All 60 `HOM` lines of `hom60.txt`, degrees -12..8, identical. The normalized-degree convention is that of `scripts/cmp_hom.py` (R1 §0.1).

**T-e. Uncolored = KL.**
- For lambda = mu = (1^n), n = 3, 4: c_q(b_p) = `Hecke::h(y, x)`, where x and y are the permutations of the permutation tables p and q. KL polynomials are invariant under simultaneous inversion, so the one-line convention cannot make this test pass wrongly. rank(B_x) = sum_y h_{y,x}(1).
- **T-e2.** For (2,2), (2,1,2), (3,3), (2,2,2) and the block (2,1) → (1,2) (K <= 6), with p_max found by brute force over `SymGroup` using `table_of` from `transfer.py`:
  - l(p_max) = lmax(p);
  - c_r(b_p) = h_{r_max, p_max} (D2 §1.1);
  - rank(B_p) * |W_mu| = sum_y h_{y,p_max}(1).
- **T-f (optional).** Hom transfer F1 up to degree 10 on the 12 blocks of `transfer.py`.
- **T-g (optional).** The 71 products of F2. Needs general Hecke multiplication in the test file.

**T-h. Rickard rigidity (port of t3.py).** For all c, d <= 5 and both signs:
- every term ladder is indecomposable: ladder character = a single b_p with coefficient 1; End^0 = 1; End^{<0} = 0;
- Hom^{<1}(W_k, W_{k+1}) = 0 and dim Hom^1 = 1, in the direction of the differential;
- Hom^2(W_k, W_{k+2}) = 0;
- for c >= d (positive), the term tables are p_k = [[d-k, c-d+k], [k, d-k]];
- for c = d = a, the ranks are C(a,k) C(a+k,k): a = 5 gives 1, 30, 210, 560, 630, 252.

**T-h2 (U13).** For every 2-upright block (c,d) → (c,d) and (c,d) → (d,c) with c, d <= 5, every canonical basis element is the character of a 2-rung ladder (so e_p = id on 2 strands).

**T-i. Decompositions (R1 §1.6).**
- (2,2):
  - B_p1 (x) W_1 = [2]^2 B_p2 + [2] B_p1;
  - B_p1 (x) W_2 = [2][3] B_p2 = (v^-3 + 2v^-1 + 2v + v^3) B_p2;
  - B_p2 (x) W_2 = [4 choose 2] B_p2.
- (2,1), with W'_k the terms of the positive (1,2) crossing:
  - C0 (x) W'_0 = D1 + D0;
  - C0 (x) W'_1 = C1 (x) W'_0 = [2] D1;
  - C1 (x) W'_1 = [3] D1.
- Split: all 21 types of (2,2,2,0) times F_2^(1) are a single canonical element times a monomial (t10.py).

**T-j. Ranks and kappa.**
- Max rank and its argmax:
  - (2,2,2): 216 at [[0,1,1],[1,0,1],[1,1,0]];
  - (3,3,3): 8 820 at [[0,1,2],[1,1,1],[2,1,0]];
  - (2,2,2,2): 25 920 at [[0,0,1,1],[0,1,0,1],[1,0,1,0],[1,1,0,0]], ignored test;
  - (2,1,2): 36.
- rank(B_{p_max}) = K!/prod lambda_i!: 90, 1 680, 2 520.
- kappa of (4,4) = 6 and of (5,5) = 10; the (4,4) ranks are {1, 20, 90, 140, 70}.

**T-k. Decategorified scans.** Expected values are in Appendix A.2, from `t8.py`, `decat_scan.py` and `peak.py`. In the default suite: 2-colored trefoil, 3-colored trefoil, 2-colored figure-eight, uncolored figure-eight, uncolored Conway. As ignored tests: 2-colored Conway, 3-colored figure-eight, 8_19 cases, the 10-crossing 3-braid.

**T-l. Closure weights.**
- e-monomial counts: emono_count((2,2), 4) = 14; ((3,3), 9) = 158; ((2,2,2), 12) = 1 218; ((2,2,2,2), 24) = 260 780; ((3,3,3), 27) = 299 958; ((5,5), 25) = 40 466; ((2,2,2,2), 22) = 159 952; ((2,2,2), 22) = 13 377.
- Final types of the 2-colored trefoil have weights {0, 3, 4}; of the 3-colored trefoil {0, 5, 8, 9}.
- 2-colored Conway (ignored test): 235 LB types, max weight 22, max unknowns 159 952 (reduced 13 377).

**T-m. Agreement with `crossing.rs`.** For c, d <= 5 and both signs, `crossing_terms(c, d, pos, i)` has the same rungs, hdeg and local q part as `crossing::crossing_cx(Gen{i, positive}, c, d, n = c + d)`. The local q part is `qshift(.., QConv::A)` + eps*m(N-m), which is N-independent.

**T-n. Decategorified Euler characteristic = `chi`** (Appendix A.3, all 8 cases). This is the earliest end-to-end check of the scan's shift bookkeeping and of global_q(N).

**T-o. Cache round trip.** `save_basis` then `load_basis` gives identical bases. The file starts with a version header.

### 8.6 CLI `schur_scan`

```
schur_scan basis   --lam 2,2,2,2 [--mu ...] [--save DIR]
    prints: types, sum of ranks, max rank (+table), rank(p_max), kappa, sum flows (bfs, minflow), time, peak RSS
schur_scan scan    --labels 2,2,2,2 --word "1 1 2 -3 2 1 -3 -2 -2 -3 -3" [--peak] [--every-stage] [--marked 0] [--load DIR]
    prints per stage and final: unsimplified objects/ldim, LB objects/types/ldim/sum S, max rank present,
    kappa present, peak; Gram histogram (weight: #objects, unknowns, unknowns reduced)
schur_scan chi     --labels 2,2 --word "1 1 1" --n 3
--json emits one JSON object per line (fields as in StageStats/GramRow plus labels, word, timings, peak_rss),
written with format! (no new crate dependency).
```

Table sweeps (M3b, cloud) read `data/knots_upto12.jsonl` and `data/links_upto11.jsonl` with a small Python driver kept in the scratchpad (`colored_design/scripts/schur_sweep.py`, run with `venv/bin/python -I`). It calls `schur_scan scan --json` once per record and color pattern, so the crate gains no JSON dependency. Colors per bottom position follow the record's component permutation. Uniform colors are trivially closable.

### 8.7 Done criteria and report

- `cargo test --release --test schur` passes in < 60 s. The ignored tests pass when run one at a time (report time and peak RSS of each).
- `cargo test --release` for the whole crate still passes, with no change to existing tests.
- The report to the orchestrator contains:
  - timings and peak RSS for (2,2,2,2), (3,3,3) and the 2-colored Conway scan;
  - the T-n outcome and which v ↔ q identification passed;
  - any mismatch with the Python reference, with a minimal example;
  - whether the BFS needed more memory than predicted (input for M3b on cloud).
- Afterwards, as a separate step and only when the user asks: a PLAN.md entry for M3a.

---------------------------------------------------------------------------------------------------

## 9. The next two tasks (enough detail to start; full specs after M3a)

### 9.1 CL1, colored closure layer `src/cclose.rs` (parallel with M3a; needs nothing from it)

**Ring.** R^lambda with variables e_{i,j}, weight j (1 <= j <= lambda_i). Monomials packed in a u128, 6 bits per variable, at most 20 variables, with overflow asserts.
- `EPoly<F>`: sparse sorted list of (packed monomial, coefficient).
- Conversion to root polynomials for the (Delta G) adapter.

**Contract v2.**

```rust
pub struct TypeBasisC<F> { pub ins: Vec<i32>, pub outs: Vec<i32>, pub gram: Vec<Vec<EPoly<F>>> } // degrees, G[b][a]
pub struct ClosureDataC<F> {
    pub lam: Comp,                    // closable: final labels = lam
    pub objs: Vec<(u32 /*type*/, i32 /*h*/, i32 /*q, N-independent*/)>,
    pub types: BTreeMap<u32, TypeBasisC<F>>,
    pub edges: Vec<(usize, usize, Vec<Vec<EPoly<F>>>)>,   // F[c][a]: rows = target ins, cols = source ins
    pub crossings: Vec<(bool /*positive*/, u8 /*m = min label*/)>,
}
```

**Closure algorithms.**
- **Exact oracle** `close_exact_colored(data, N, marked: Option<usize>) -> Poincare`. A_lambda is (x)_i H*(Gr(lambda_i, N)) in the `symfun::GrTable(lambda_i, N)` e-monomial box basis at c = 0, with tau from `sbim::Tau` (R3 A2: the same functional as the Schur top class). The chain group of object o is A_lambda^{ins} shifted by deg; Phi_N(B) = Im(G (x) A_lambda); ranks per (h, q) by `RankEchelon`.
  - Reduced mode: substitute e_{m,*} = 0 in G and F and drop the marked block's term of global_q (R3 A7).
- **(Delta G) adapter.** Expand G and F in root variables, multiply G by Delta_lambda = prod over blocks of Vandermondes, call `tower::close_tower` with caps N-1 per root (K <= 10), then divide the Poincaré polynomial by prod_i gr(V_{lambda_i}).
  - The normalization of gr(V_a) (t^{l(w_a)} pi_a(t), t of q-degree 2, plus a shift) is fixed on single types against the exact oracle.
  - Reduced: drop the marked block's roots and its Delta factor (C2 M1(c)).
  - **Usable in production only after G4.**

**Tests.**
1. tau = `sbim::Tau` on random elements.
2. For 30 random small closed ladders (2-3 uprights, K <= 6), R^lambda-bases of Hom(R^lambda, L) from exact `hom_space` (graded Nakayama, exact), Gram by composition, closure by the exact oracle: graded dims = `StateSpace::compute` = `moy::closure_value`.
3. Exact oracle = (Delta G) adapter on those cases.
4. On engine2 `ClosureData` (lambda = 1^n), the new exact closure = `close_exact` for the Conway knot at N = 2..6, unreduced and reduced.
5. G4 script: freeness of Delta_a R_N, extending C2's `vandermonde_free.py` to all a <= 5, N <= 10. Large cases go to the cloud.

### 9.2 M6.0, exact small colored engine `src/xengine.rs` (after M3a and CL1)

**Objects.** (type p from schur, h, t). The ladder is the 2-rung canonical ladder (T-h2), and e_p = id.

**Maps.** Exact `pmat::PMat` over R^lambda (from `sbim`):
- composition = `pmat::matmul`;
- f (x) id_W = `PMat::kron_id(fsize(W))`, valid because W's fibre basis is free over the top ring of the source;
- id (x) d_k = `cube::apply_local_poly` (make it `pub`, or copy it) with d_k from `cube::local_cx` (made `pub`), built untruncated with n >= c + d.

**Tensoring and decomposition.**
- X = L_p (x) W_k (`SBim::add_rungs`).
- For each summand type q with multiplicity m_q(v) from `schur::term_mult`, top-down by l_max:
  - iota_i = `hom_space(L_q, X, s)` and pi_j = `hom_space(X, L_q, -s)`;
  - the intersection form G_ij = scalar of pi_j ∘ iota_i, read exactly, because End^0(L_q) = k * id;
  - require rank G = [v^s] m_q; orthogonalize as in `indec.rs`.
- Assert that the summand ranks add up to fsize(X).

**Elimination.** Scalar pivots (End^0 = k). Markowitz order. Koszul sign (-1)^h on id (x) d.

**Closure data, exact.**
- Generators: `hom_space(1_lambda, L_q, d)` in the degrees of c_diag(b_q). Choose R^lambda-generators by the one-point determinant criterion (C2 M2), plus an exact rank check.
- Gram = h_b ∘ g_a, exact. F by an exact solve over R^lambda, degree by degree.
- Closure: CL1's exact oracle.

**Scope.**
- 2-strand braids with labels (1,1), (2,1), (1,2), (2,2).
- (3,1), (3,2), (3,3) within a 10-minute budget per `hom_space` call. Record which ones fit: X fibres reach 900 for (3,3), C2 H3(c).

**Acceptance.** As in §3.1 M6.0:
- V1, V2 (including internal labels > N: N = 2, 3 for (2,2));
- V4: the R3 V4 Hopf polynomials for (2,1), (2,2) at N = 3, 4, 5, 6, bigraded;
- V5: `wang.py` bigraded for k = 2, T(2,3), T(2,5), T(2,7) at N up to 10 where cheap;
- V9: `cube::khr_full_cube` on Hopf (2,1), (2,2) at N <= 4 and on (2,2) kinks;
- V8: cabling;
- reduced unknot = 1, and the Lewark-type relation;
- χ = `schur::decat_chi`;
- two primes.

---------------------------------------------------------------------------------------------------

## Appendix A. Reference numbers (expected outputs)

### A.1 Algebroid data

These are in tests T-a, T-b, T-c, T-j (§8.5). Sources: R1 §1.6 and §6.2, D3 §10.2-10.3, `cb.pkl` aggregates (recomputed here with `synth/cbsize.py`).

| lambda | types | entries in all b_p | Laurent terms | sum of ranks | kappa | longest BFS ladder |
|---|---|---|---|---|---|---|
| (1,1,1) | 6 | 19 | 19 | 19 | 1 | 4 |
| (1,1,1,1) | 24 | 213 | 219 | 219 | 2 | 8 |
| (2,2) | 3 | 6 | 7 | 13 | 2 | 2 |
| (3,3) | 4 | 10 | 15 | 63 | 3 | 2 |
| (2,2,2) | 21 | 175 | 248 | 1 045 | 8 | 6 |
| (3,3,3) | 55 | 1 016 | 2 284 | 76 399 | 54 | 6 |
| (2,1,2) | 11 | 62 | 78 | 199 | 4 | 5 |
| (2,2,2,2) | 282 | 19 088 | 35 970 | 397 687 | 63 | 10 |

### A.2 Decategorified scans

Sources: D1 §8.3, D2 §6.2, D3 §10.4, C2 H2. All were produced by `t8.py`, `decat_scan.py` and `peak.py`.
- Values written as 9.15e6 and similar are rounded. The Rust tests should use exact values obtained by rerunning the Python reference once (cheap, except the two Conway rows and the 3-colored 8_19 row, which take seconds to minutes).
- "Peak" is the expanded LB complex tensored with the next crossing, before elimination, at the worst crossing: objects / localized dim / sum S.

| case | lambda | word | unsimplified objects | LB objects | LB types | LB ldim | LB sum S | max rank | peak |
|---|---|---|---|---|---|---|---|---|---|
| 2-colored trefoil | (2,2) | 1 1 1 | 367 | 11 | 3 | 61 | 67 | 6 | 81 / 481 / 495 |
| 3-colored trefoil | (3,3) | 1 1 1 | 12 332 | 24 | 4 | 507 | 693 | 30 | 802 / 16 695 / 18 489 |
| 4-colored trefoil | (4,4) | 1 1 1 | 460 461 | 47 | - | 4 121 | 7 379 | 140 | - |
| 5-colored trefoil | (5,5) | 1 1 1 | 1.8e7 | 82 | - | 32 151 | 79 421 | 630 | - |
| 2-colored T(2,5) | (2,2) | 1 1 1 1 1 | 61 883 | 19 | 3 | 109 | 119 | 6 | - |
| 3-colored T(2,5) | (3,3) | 1 1 1 1 1 | 4.96e7 | 42 | 4 | 911 | 1 253 | 30 | - |
| 5-colored T(2,5) | (5,5) | 1 1 1 1 1 | 5.4e13 | 196 | - | 74 603 | 172 253 | 630 | - |
| 4- / 5-colored Hopf | (4,4) / (5,5) | 1 1 | - | 21 / 33 | - | 1 801 / 12 325 | - | - | - |
| uncolored figure-eight | (1,1,1) | 1 -2 1 -2 | 24 | 12 | 5 | 33 | 33 | 4 | - |
| 2-colored figure-eight | (2,2,2) | 1 -2 1 -2 | 533 | 77 | 16 | 3 049 | 6 073 | 216 | 325 / 17 641 / 26 733 |
| 3-colored figure-eight | (3,3,3) | 1 -2 1 -2 | 11 748 | 328 | 40 | 348 153 | 1 984 617 | 7 200 | 5 036 / 6.66e6 / 2.04e7 |
| 2-colored T(3,4) = 8_19 | (2,2,2) | 1 2 1 2 1 2 1 2 | 9.15e6 | 187 | 21 | 13 933 | 27 753 | 216 | - |
| 3-colored 8_19 | (3,3,3) | 1 2 1 2 1 2 1 2 | 1.24e11 | 1 084 | 55 | 2 566 793 | 13 023 797 | 8 820 | - |
| 2-colored 3-braid, 10 crossings | (2,2,2) | 1 1 -2 1 -2 1 1 -2 -2 1 | 1.62e9 | 8 835 | 21 | 546 553 | 969 741 | 216 | 33 369 / 1.90e6 / 2.90e6 |
| uncolored Conway | (1,1,1,1) | 1 1 2 -3 2 1 -3 -2 -2 -3 -3 | 17 248 | 192 | 19 | 1 439 | 1 471 | 16 | 434 / 3 591 / 3 623 |
| 2-colored Conway | (2,2,2,2) | same | 1.47e9 | 11 029 | 235 | 13 282 393 | 66 632 927 | 25 920 | 103 573 / 1.03e8 / 3.40e8 |

Closure weights of the final LB types (D1 §6.3, C3 E4):
- 2-colored trefoil: {0, 3, 4};
- 3-colored trefoil: {0, 5, 8, 9};
- 2-colored Conway: weights 0..22, bulk 11-17, max unknowns 159 952 (reduced 13 377);
- 3-colored figure-eight: max 88 233 (reduced 3 866).

### A.3 Euler characteristics from `chi`

Produced with `colored_design/target/release/chi N labels word` (cube of resolutions, `moy::closure_value`, convention A). These are the expected values of test T-n.

| N | labels | word | chi |
|---|---|---|---|
| 3 | 2,2 | 1 1 1 | q^-2 + q^-4 + 2q^-6 + q^-8 - q^-12 - q^-14 |
| 4 | 2,2 | 1 1 1 | q^-4 + q^-6 + 3q^-8 + 3q^-10 + 3q^-12 + q^-14 - 2q^-18 - 2q^-20 - 2q^-22 - q^-24 + q^-30 |
| 3 | 2,1 | 1 1 | 1 + 2q^-2 + 2q^-4 + 2q^-6 + q^-8 + q^-10 |
| 4 | 2,2 | 1 1 | 1 + q^-2 + 3q^-4 + 4q^-6 + 6q^-8 + 6q^-10 + 6q^-12 + 4q^-14 + 3q^-16 + q^-18 + q^-20 |
| 3 | 2,2,2 | 1 -2 1 -2 | q^8 + q^6 - 1 + q^-6 + q^-8 |
| 4 | 2,2,2 | 1 -2 1 -2 | q^16 + q^12 + q^10 + q^-10 + q^-12 + q^-16 |
| 6 | 3,3 | 1 1 1 | q^-9 + q^-11 + 3q^-13 + 5q^-15 + 7q^-17 + 8q^-19 + 10q^-21 + 9q^-23 + 7q^-25 + 4q^-27 - q^-29 - 5q^-31 - 8q^-33 - 10q^-35 - 10q^-37 - 7q^-39 - 5q^-41 - q^-43 + q^-45 + 3q^-47 + 3q^-49 + 3q^-51 + 2q^-53 + q^-55 - q^-63 |
| 4 | 2,1,1 | 1 2 2 1 | q^-2 + 3q^-4 + 5q^-6 + 8q^-8 + 11q^-10 + 13q^-12 + 14q^-14 + 13q^-16 + 11q^-18 + 8q^-20 + 5q^-22 + 3q^-24 + q^-26 |

Consistency checks already visible in the table:
- (2,2) at N = 3 on "1 1 1" is the uncolored gl_3 trefoil of R3 V5: t^0(q^-2+q^-4+q^-6) + t^2(q^-6+q^-8) + t^3(q^-12+q^-14).
- (2,2,2) at N = 3 on the figure-eight equals the uncolored gl_3 figure-eight (`chi 3 1,1,1 "1 -2 1 -2"` gives the same polynomial). This is duality at N = a+1, decategorified.
- (2,2) at N = 4 on "1 1" equals the R3 V4 Hopf polynomial at t = -1.

### A.4 Closure-data and fixed-point counts (`synth/emono.py`, `synth/emono2.py`)

**e-monomials of R^lambda by weight 1..8:**

| lambda | weights 1..8 |
|---|---|
| (2,2) | 2, 5, 8, 14 (up to W = 4) |
| (3,3) | 2, 5, 10, 18, 30, 49, 74, 110 |
| (2,2,2) | 3, 9, 19, 39, 69, 119, 189, 294 |
| (3,3,3) | 3, 9, 22, 48, 96, 182, 324, 555 |
| (2,2,2,2) | 4, 14, 36, 85, 176, 344, 624, 1 086 |
| (2,1,1,1) | 4, 11, 24, 46, 80, 130, 200, 295; weight 9: 420 |
| (2,2,1,1) | 4, 12, 28, 58, 108, 188, 308, 483; weight 13: 2 856 |
| (2,2,2,2,2) | at W = 40: 1.03e8 |

These give the stacked-point counts M of C1 I1. For example, (2,2,2,2) needs M = 4, 14, 36 for (d - d_j)/2 = 1, 2, 3.

**Grassmannian fixed points per N** (box-restricted closure-data interpolation, upper bound before the weight cut):

| lambda | N = 3 | N = 4 | N = 5 | N = 6 | N = 8 | N = 10 |
|---|---|---|---|---|---|---|
| (2,2) | 9 | 36 | 100 | 225 | 784 | 2 025 |
| (3,3) | 1 | 16 | 100 | 400 | 3 136 | 14 400 |
| (2,2,2) | 27 | 216 | 1 000 | 3 375 | 21 952 | 91 125 |
| (3,3,3) | 1 | 64 | 1 000 | 8 000 | 175 616 | 1 728 000 |
| (2,2,2,2) | 81 | 1 296 | 10 000 | 50 625 | 614 656 | 4 100 625 |

## Appendix B. Provenance

All paths are under `/private/tmp/claude-501/-Users-existair-Documents-KRcalculator/bf6f8870-08ad-436b-9ff7-0414655cc59a/scratchpad/colored_design/`.

| source | what it provides |
|---|---|
| `R1_singular_theory.md`, `R2_code_map.md`, `R3_closure_and_validation.md` | research notes |
| `D1.md`, `D2.md`, `D3.md`; `C1.md`, `C2.md`, `C3.md` | designs and critiques |
| `scripts/schur.py`, `scripts/t3.py`, `scripts/t5.py`, `scripts/t8.py`, `scripts/t10.py`, `scripts/cmp_hom.py`, `scripts/cb.pkl`, `scripts/homcases.txt`, `scripts/homres.txt` | Python reference implementation and fixtures for M3a |
| `proto_D1/peak.py`, `proto_D1/decat_scan.py` | decategorified scan reference |
| `proto_D2/transfer.py` | Hom transfer and ordinary-KL cross-checks |
| `proto_D2/closurecost.py` | top Gram weight = l(p_max) - l(w_lambda) |
| `c1scripts/*`, `proto_C2/*`, `proto_C3/*` | critique evidence |
| `r3proto/src/bin/chi.rs` (binary `target/release/chi`) | Euler characteristics, A.3 |
| `r3scripts/wang.py`, `r3scripts/hopf_ours.py` | bigraded references for M6.0 |
| `synth/emono.py`, `synth/emono2.py`, `synth/cbsize.py` | numbers in A.1 and A.4 (this note) |

Literature as cited in the notes:
- EKLP arXiv:2401.03053: Thm 4.13, Thm 4.16, Cor 4.17, Prop 4.19, Prop 4.21, Thm 6.11, Cor 6.12, Thm 7.49, Thm 7.51, eqs. (44), (45), (77a/b); (102)-(106) as cited by C2, not re-checked.
- Williamson arXiv:1010.1283v3: Thm 3, Prop 7.4.3, Thm 7.4.1, Thm 7.5.1, erratum §8.
- Elias, Thick Soergel calculus in type A, arXiv:1009.2120.
- ETW arXiv:1703.06691: Def 3.3, Thm 3.5, Conv 3.1.
- KRW arXiv:2603.22957.
- Wedrich arXiv:1602.02769: Def 2.17, Thm 4.
- J. Wang arXiv:2211.08409 (Def 2.23, Prop 1.5, Lemma 3.4) and arXiv:2604.25174 (Thm 1.2).
- Hogancamp-Rose-Wedrich arXiv:2107.09590, Prop 9.4.
- Rose-Wedrich arXiv:1501.02567; Wu arXiv:1011.2254 (deformation side goal).
