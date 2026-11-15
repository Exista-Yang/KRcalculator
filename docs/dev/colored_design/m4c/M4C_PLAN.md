# M4C_PLAN: decision and plan for M4c (complete Hom bases for the (3,3,3) library)

Date 2026-11-11. Read-only synthesis; nothing in the repository was compiled or run.
Inputs (this folder): L1_eklp.md (paper extraction), L2_alt.md (alternatives), L3_code.md (code map), D1_eklp.md
(design: Lemma A, rung-level), D2_light.md (design: padded cell-local EKLP), C_m4c.md (adversarial critique).
Code: dev 277a21c. M5b data: origin/colored-m5 docs/dev/colored_design/gl6/M5b_report.md, m4c/data/costs_333.txt.
Paper keys as in L1: SLL = EKLP arXiv 2401.03053 (JEMS 2025); CORE = arXiv 2402.08673v2; EKo = arXiv 2105.08563;
QR = arXiv 1405.5920; W = Williamson IMRN 2011.
Tags: [V] verified for this synthesis (code read, paper text via L1/C, or the laptop-safe replica m4c/c_tools/hu.py,
< 1 s per call); [K] inference, to be checked by a gate; [E] estimate.

## 0. Decision

1. **Pursue the hybrid "M4c-H": localized light leaves built by D1's rung-level recursion (Lemma A).** Concretely:
   - **Construction:** D1's Lemma A recursion. Its completeness argument rests on SLL Prop 8.32 and Cor 8.33.
   - **Stage 0 first:** D2's product rule for closure.
   - **Generators and moves:** D2's generator library, as closed forms only. Braid and rex moves are composed from
     closed forms after CORE Lemmas 3.12 and 3.13 and SLL Thm 6.7. They are never solved for as exact windows.
   - **Padding:** D1's trailing-zero `Pad`.
   - **Hard census budget** (C §3), enforced before any generator code.
2. **Rejected parts:**
   - D2's per-singlestep recursion through 4-block levels (C M2);
   - D2's `Repad` (C M7) and its lazy FlowIndex (C M3);
   - D1's exact braid and rex windows (C C1);
   - D1's "Leaves before Exact" production order (C M4).

   Section 2 gives the reasons.
3. **Budget rule (gates G-L0 and G-P0)**, by construction:
   - no exact `hom_space` call above 10⁵ unknowns; a call up to 4·10⁵ is allowed only with an uncensored exact-shape
     cost record;
   - no interned expression above 10⁶ flows;
   - otherwise the step is built from closed forms, or the target goes to a fallback.

   Calibration from costs_333.txt [V]:

   | unknowns | records | censored | max uncensored | max RSS |
   |---|---|---|---|---|
   | 10⁴-10⁵ | 74 | 1 | 48 s | 0.5 GB |
   | 10⁵-10⁶ | 136 | 92 | 419 s | 2.6 GB |
   | ≥ 10⁶ | 32 | all | - | - |
4. **Order: closure, then small-λ completeness, then (3,3,3) pairs.**
   - **Closure** needs no frozen file and no new mathematics beyond Lemma A. It closes the 12 + 1 closure jobs (N = 6 and
     N = 5) in about 3 sessions.
   - **Small-λ completeness:** G1(a) and B2 with Leaves alone, using exact steps, which are affordable there. About 2
     sessions.
   - **(3,3,3) pairs** carry the research risk: localized closed forms for rex moves, plus memory. They go through a
     decision gate (G-P0) after 2 sessions of generator and planner work.
5. **Outlook** (details in §8):
   - The closure list closes with about 85 % probability.
   - Every N = 6 gap closes, with M5c's help, with about 35-40 % probability in the expected 15-16 sessions, and about
     50-55 % with the reserve (up to about 22 sessions).

## 1. What was verified for this synthesis

- **Replica reproduces M5b.** hu.py reproduces measured M5b numbers, e.g. 852 011 unknowns for Hom(∅, E0.3 F0.3 E1.3
  F1.3) at d = 18 [V]. For pair 322, the whole-word exact call X → L_q at d = 0 has **919 495 095** unknowns [V].
- **Braid windows (C1) [V], braid.py.** Unknowns at d = 0 for a braid move on cells (x,y,z):

  | cells | (2,1,2) | (2,2,2) | (3,1,3) | (2,3,2) | (3,2,3) | (3,3,3) |
  |---|---|---|---|---|---|---|
  | unknowns | 18 019 | 711 580 | 2 966 475 | 12 609 774 | 3.6·10⁸ | 1.6·10¹⁰ |

  A far commutation on cells (3,3,2,1) has 192 484.
- **New finding: even single generators exceed the cap as exact windows at large blocks** [V, hu.py].
  - The sideways crossing [Is+t−s] ⇄ [Is−s+t] on blocks (a,b,c) is, in krhom, the 2-upright pair F0.b ⇄ E0.c F0.(b+c)
    over labels (a+b, c). Its degree is ac (SLL; checked by hand from ℓ(I)+ℓ(Ist)−ℓ(Is)−ℓ(It)).
  - Exact unknowns at |d| = ac:

    | blocks (a,b,c) | (2,2,2) | (2,3,2) | (3,3,3) |
    |---|---|---|---|
    | at d = +ac | 1 036 | 3 479 | **233 282** |
    | at d = −ac | 85 | 380 | 3 537 |

  - So D2's claim "≤ 10⁴ unknowns per generator window" fails at blocks (3,3,3). Large sideways crossings must be
    closed-form templates. The sign of krhom's normalized degree against SLL's is checked by G-L1.
- **Elem validation windows are affordable** [V].
  - E2 Eta/Eps over (a,b): at most 112 390 unknowns, at (4,5); 42 872 at (3,6).
  - E3 Iota/Trace over (c,0): 1 unknown.
  - Production validation of Elem constants is therefore cheap, either by exact comparison or by `tmpl::poly_certificate`
    (tmpl.rs:677) with a nonvanishing check, since these Hom spaces are 1-dimensional (L3 §4.4).
- **Cell-form flows (C M2) [V], yall.py.**
  - 3-column levels: at most 944 784 flows (the decomposable type 0,1,2/1,1,1/2,1,0), otherwise at most 419 904;
    43 of 811 tables exceed 10⁵.
  - 4-column levels: up to 13 436 928; 990 of 4 380 tables exceed 10⁵.
  - `FlowIndex::build` (sexpr.rs:130-175) holds `flat` and `masks` at once, so the build peak is twice the final index.
- **The trailing-zero `Pad` is the identity** [V]. `CosetSpace::rank` (coset.rs:116) multiplies by radix C(0,0) = 1 and
  adds 0 for an empty trailing upright. D2's `Repad` re-indexes cosets and is a trusted translation layer.
- **Code facts [V]:**
  - `Arena::expr` asserts bottom == λ (smap.rs:580). `Arena::flows` panics without an index (smap.rs:633). The Arena is
    documented `Send + Sync`, with no interior mutability (smap.rs:1-6, :499). So C M3 holds.
  - `certify_elem` and `certify_assoc` store any scalar given (smap.rs:821, :831). So C M5 holds.
  - D11 freezes homprov.rs, g1.rs, sbim.rs, linsys.rs, smap.rs and tmpl.rs, but not smapio.rs (SPEC_3C_GL6_rev1:165-168).
    So C m8 holds.
- **N = 5 failures** (M5b_report, both primes):
  - Refused pairs 313, 326, 406 and 410, with X = 600-3 600; they are cheap pairs whose windows moved;
  - pair 322, Exhausted, rex class;
  - closure 1,1,1/1,1,1/1,1,1, Exhausted.

## 2. Why this design: findings adopted or rejected

| finding (source) | verdict | consequence here |
|---|---|---|
| Both designs rest on proved results (Prop 8.32, Cor 8.33, Thm 7.49); no generators-and-relations paper is needed (L1 §4, C §0) | adopt | SPEC's "companion paper" risk is dropped; certification stays the only acceptance criterion |
| D1 Lemma A is correct; closure via the identity coset needs no lower terms (C, verified strengths) | adopt | closure = Lemma A on 3 uprights; invertibility is generic (rational functions) and checked only at random points, so it is a diagnostic and never accepts anything |
| C1: exact windows are unbounded in both designs | adopt, and extend (sideways crossings, §1) | budget rule §0.3; every braid, rex and large generator is a closed form; census before code |
| m3: choose W_t := W_{t*}·ρ (prefix extension) when it is reduced | adopt | top steps become `Id`, which removes the rex-type closure steps of C1 |
| m2: types with no 3-upright word that has exactly one flow at t | adopt | G-L0 lists them; they route through closed-form steps (B-phase) |
| m4: greedy selection must include all lower degrees of supp G | adopt | prevents spurious rank shortfalls |
| M2: 4-block levels of up to 1.3·10⁷ flows (D2) | adopt, reject D2's recursion | the recursion stays at rung level; junction words live at 3-column levels; 4-upright words appear only inside budgeted rex chains |
| M3: a lazy FlowIndex is a redesign of the frozen smap.rs | adopt | no laziness; budget on flows instead (≤ 10⁶ per interned expression ≈ 80 MB of masks at 10 rungs, peak 2×) |
| M7: `Repad` enlarges the trusted base; trailing-zero `Pad` is the identity | adopt D1 `Pad` | one additive smap.rs change (D11 exception), only if G-P0 needs a 4th upright |
| M5: Elem/Assoc constants are caller-trusted | adopt | `validate_elem` in leaves.rs: exact comparison when ≤ 1.2·10⁵ unknowns, else `poly_certificate` plus nonvanishing. Assoc only when both windows have the same multistep EKLP expression and ≤ 500 flows; never an unvalidated constant |
| M4: downstream M6a cost of deep Leaves DAGs | adopt | production chain puts Leaves **last**, so only failing targets get Leaves maps; gate G-L5 measures t_lib and resident index bytes before M6a-3 uses them |
| M1, M6: D2's normalizer is incomplete; D1's factored and chart layer is undesigned | adopt | the rex planner and chart layer are explicit B-phase work, with a decision gate (G-P0) |
| M8: effort about 2× the design estimates | adopt | §8: about 15-16 sessions expected, not 8-11 |
| m5: double leaves need left R^λ monomials (Dot at level 0) | adopt | in the `full`/`end0` milestone |
| L2: KL/foam thick spanning set "S7" as the primary route | reject | it needs the same padding; its streams are redundant (10⁴-10⁸ for large X); thick crossings on 9 roots are also beyond the cap as exact windows. Kept only as a source of closed-form crossing formulas (KL3 §5-6 polynomial action) |
| (new) Williamson bending for pairs, Hom(L_q, X) ≅ Hom(R, X·L̄_q) | reject | X·L̄_q has up to 21 600 × 7 200 ≈ 1.6·10⁸ flows, so the global FlowIndex alone is ≈ 15 GB (at about 12 rungs, 8 B per mask) |
| m11: L2 says the EKLP risk stands | reject | CORE (arXiv 2402.08673v2) is in papers/ and is used constructively only; reducedness is gated by `ladder_char` (m9) |

## 3. Construction summary

### 3.1 Lemma A (localized light leaves; D1 §1.3, C verified)
Let X = X'·ρ, with fixed reduced junction words W_{t'} at the top level of X' and W_t at the top level of X. Suppose:
- (i) each family LL_{t'} ⊂ Hom(X', W_{t'}) has an invertible block at rep(t');
- (ii) each step family S_{t'→t} ⊂ Hom(W_{t'}·ρ, W_t) has an invertible top-to-top matrix at rep(t).

Then LL_t = {S ∘ Kron(ℓ, [ρ])} has an invertible block at rep(t), and the degrees add.

- **Completeness** (C verified, with C m1's corrected citation): SLL Prop 8.32 at the expression W_{t'}·ρ, together with
  the injectivity argument of the proof of Thm 8.31 (Prop 8.26 + Prop 4.20). The step families with degree multiset
  G_{t'→t} exist (Lemma 8.7), so greedy selection never runs dry.
- **Dual side:** Cor 8.33 (columns, for maps W_t → X).
- **Step generating functions** come from `Schur::split_col`/`merge_cols` (schur.rs:913/:931). Their Deodhar sums are
  checked against `ladder_char` (:997) and `sexpr::deodhar_check` (:267) (gate G-L1).

### 3.2 Closure (outs: L_p → R; ins: R → L_p)
- **Stage 0 (D2 §2.1, C "do it first").** Covers the 5 product types:
  0,0,3/3,0,0/0,3,0; 0,3,0/0,0,3/3,0,0; 1,0,2/2,0,1/0,3,0; 1,2,0/0,0,3/2,1,0; 0,3,0/1,0,2/2,0,1.
  - L_p ≃ A·B, after at most one commutation. Their c_diag(p) = c_diag(A)·c_diag(B) [L3 §3.2, D].
  - Generators: Comp(b_j, Kron(a_i, B)) and the flips, from the factors' certified closure maps (layer A, level ≤ 4).
  - The commutation iso is an Assoc. Both windows are the same multistep expression; e.g. E0.2/F1.3 at labels (3,6,0)
    has 60 flows, in bijection with the top configurations [V]. Its validation window is 34 189 unknowns [V].
  - `certify_closure` (homprov.rs:945) decides.
- **Lemma A on 3 uprights** for the 7 other types (2 700, 2 880, 4 320 and 7 200 flows):
  - the target is the identity coset, the junction words are short 3-upright words chosen by the prefix-extension rule
    (m3), and the last junction word is ∅;
  - step families, in order: `Id` → E2/E3 `Elem` with validated constants → exact window under the budget → closed-form
    steps (only if G-L0 flags the step; B-phase generators);
  - n = c_diag(b_p)(1) leaves with degree multiset c_diag(b_p). `certify_closure` checks rank n at the identity coset at
    2 points;
  - the 1,1,1/1,1,1/1,1,1 type has n = 54 (degrees 4-18); it is the N = 5 job too.

### 3.3 Pairs (π: X → L_q of degree −s, ι: L_q → X of degree s)
- **Small λ: (3,3), (2,2,2), (2,1,2), (2,2).**
  - Lemma A with W_q := L_q. The last step is an exact window S ⊂ Hom(W_{t'}·ρ, L_q). That is affordable: the Exact
    backend already certifies 264 G1(a) targets with such calls (HANDOFF_M4b).
  - No rex planner and no `Pad` are needed.
  - This validates selection, flips, degree windows and `certify_pairs` (homprov.rs:790) on the full G1(a) set.
- **(3,3,3):**
  - X is consumed rung by rung. Junction words are prefix extensions where they are reduced, so the top-path steps are
    `Id`. Lower steps are single-step leaves after the Grassmannian reduction (SLL Thm 2.6, Lemma 3.9).
  - Each lower step is realized as id_Z ⊗ s:
    - with s an exact window on the K uprights **only if** the third column is non-empty and G-P0 shows ≤ 10⁵ unknowns;
    - otherwise s is composed from closed forms (§3.4).
  - The final junction word W_q is a prefix extension of X's top path. The exit and entry maps W_q ⇄ L_q are **rex
    chains of closed-form moves**:
    - up-up and down-down (E4 Assoc);
    - type-A switchbacks with d = 2 (EKo eq. (108)) as φ/ψ composites (SLL Thm 6.7: right/left-facing crossings, ccw
      cups and caps);
    - atomic braid moves = switchback ∘ (up-up, down-down) ∘ switchback through the merged expression (CORE Lemma 3.12);
    - far commutation (CORE Lemma 3.13).
  - Chains are planned combinatorially (EKo Thm 5.30 is constructive; CORE Thm 3.14 for cores) on singlestep forms,
    then realized as krhom words, with a trailing-zero `Pad` when a 4th block needs its own upright.
  - Every intermediate word is checked reduced (`ladder_char`: c_t(1) = 1, support ≤ t).
  - Any degree-0 map with top scalar ≠ 0 between reduced words of q is admissible (Prop 4.19, Thm 6.11; L1 §5.2).
  - π = exit ∘ LL and ι = LL^∨ ∘ entry. One exit/entry pair per q is shared by all its targets.
  - `certify_pairs` selects a rank-mult minor at rep(q). The form rank = mult rests on Thm 7.52 plus the local
    intersection form [K]. A shortfall is a visible failure, never a wrong answer.
- **The rex class** (71, 159, 322, 323, 325, 334, 419: X ≅ B_q, s = 0, mult 1) is exactly "exit/entry for X itself": no
  leaves are needed beyond the identity. These pairs are therefore the cleanest test of the closed-form rex machinery.
  Pair 322 is also the N = 5 blocker.

### 3.4 Generators and how each is trusted
| generator | realization | trust |
|---|---|---|
| polynomial | `Dot` (E1) | none needed |
| cw cup/cap, whole block | `Elem` Eps/Eta(+L) | `validate_elem`: exact 2-upright comparison (≤ 112 390 unknowns [V]) or `poly_certificate` + nonzero |
| ccw cup/cap (digon) | `Elem` Iota/Trace(+L) | same (1 unknown [V]) |
| general cw/ccw cup/cap E^(r)F^(r) ↔ 1, r < block | Added template from the Frobenius formulas (SLL §5.2) | `TmplLib::certify_added` (tmpl.rs:392: degree, exact commutation, `poly_certificate`) |
| up-up/down-down (same multistep expression) | `Elem` Assoc | combinatorial check that the multistep expressions are equal, plus `poly_certificate` when ≤ 500 flows; refuse otherwise |
| sideways crossings | Added template from SLL eqs. (91)-(92) (right-facing: multiplication; left-facing: f ↦ ∂Δ(f)), or the eq. (84) composite | `certify_added`; windows of 20 x 84 flows at blocks (3,3,3), so cheap to certify although the exact call is 2.3·10⁵ unknowns |
| switchback φ/ψ, braid, far commutation | `Comp`/`Local` DAGs of the rows above (SLL Thm 6.7; CORE Lemmas 3.12, 3.13) | genuine because the children are; ψ∘φ = id and top scalar ≠ 0 checked as diagnostics |
| small residual windows | `TmplLib::try_hom_basis_cap` (tmpl.rs:441) under the budget | HomSpace, genuine |

Where a generator's Hom space has dimension > 1, use the closed form, never "any element" (D2 R5, C m6).

### 3.5 Double leaves (B2 gate; later `full`/`end0`)
- Candidates: LL^∨_r(W_r → B) ∘ Dot_g ∘ LL_r(A → W_r), times left R^λ monomials (C m5). g runs over the Schur basis of
  R^{LR(r)} over R^λ.
- Both sides use the same junction word W_r, so no exit rex is needed.
- `certify_full_plan` (homprov.rs:635) selects at the Hom-formula count (Thm 7.49 says the span is complete).
- At small λ the steps are exact windows. This is the strongest small-case test of the Lemma A machinery in all degrees.

### 3.6 Certification and failure semantics (unchanged)
- **Acceptance** happens only through the existing checks:
  - `certify_closure`, `certify_pairs` (form at rep(q) at 2 points, genuineness), `certify_full_plan`;
  - sindec `dechom_postcheck` (sindec.rs:330, D14 scaling), `closure_gens` re-certification on load (:1007), replay.
- **Internal rank tests** are diagnostics.
- **Failure kinds:**
  - `Refused`: a step violates the budget, or a constant is unvalidated;
  - `Time` / `Memory`: the budget was exceeded at run time;
  - `Exhausted`: a rank shortfall after every window was computed. This is a red-flag bug report, since §3.1 says it
    cannot happen with correct code.

## 4. Milestones (cloud sessions of about 4-5 h)

| id | branch | content (files) | tests and gate (exit criterion) | sessions |
|---|---|---|---|---|
| A1 | m4c-census | new src/dcoset.rs: `conf_type`, `step_gf`, `reach_back`, `junction_word` (prefix-extension BFS + `ladder_char` reducedness), `grass`, `psi`; census in new src/bin/m4c.rs; Stage 0 in new src/leaves.rs (`product_closure`, inner = current chain) | tests/dcoset.rs: **G-L1** (Deodhar sums = `ladder_char` for every canonical ladder of (3,3), (2,1,2), (2,2,2), (3,3,3) and their one-rung extensions; Lemma 8.7 sums = `c_diag`/`term_mult`; degree sign of §1). **G-L0**: per closure job, step-call unknowns and junction-word flows, and types without a one-flow word. **Stage 0**: 5/5 product types at N = 6, F31 + F32b, replay | 1 |
| A2 | m4c-closure | leaves.rs: `build_leaves` (Out/In), `step_family` (Id → Elem → exact under budget), greedy selection over all lower degrees, `closure()`, `validate_elem`; re-implemented private helpers (`cdiag_degrees` :989, `identity_coset` :538, `rep_rank` :529); libbuild `--provider leaves` | tests/leaves.rs: **V0** vs `sbim::hom_space` at λ ∈ {(2), (1,1), (2,1), (1,2,1)}; SLL Ex 7.39-7.41 (G-L2 closure part); **G1(a) closure part, Leaves alone**: (3,3) 8/8 and (2,2,2) 42/42 targets, both primes | 1 |
| A3 | m4c-closure | (3,3,3) closure: level-4/5 closure regression (5 + 19 jobs) with Leaves alone; then the 13 failing closure jobs; wiring in `provider_for_stamped` (sindec.rs:421) **after** S3-S6, closure calls only; new stamp; merge to dev | **closure list 13/13** (N = 5: 1,1,1/1,1,1/1,1,1; N = 6: the 12 types), both primes; replay; **G-L5** (t_lib and resident index bytes of the new maps against the M6a budget, R4_m6a §4) | 1 (+1 reserve if G-L0 flags steps) |
| A4 | m4c-small | `pairs()` via Lemma A with W_q = L_q and exact final steps; `end0`/`full` stubs refuse | **G1(a) with Leaves alone: 480/480** ((3,3) 72 + 8, (2,2,2) 358 + 42), both primes | 1 |
| A5 | m4c-dll | double leaves `full`/`end0` (§3.5) | **B2 with Leaves alone**, degrees 0..4: (2,2) 45/45, (2,1,2) 605/605 (Search alone: 581), (2,2,2) 2 205/2 205 (Chain: 2 064), (3,3) 80/80 | 1 |
| B1 | m4c-gen | new src/leafgen.rs: general cups/caps and sideways crossings as Added templates (SLL §5.2, eqs. (84), (91)-(92)); switchback φ/ψ (Thm 6.7); braid and far moves per CORE Lemmas 3.12/3.13; Assoc validator | unit tests at small blocks and cells ((1,1,1), (1,2,1), (2,1,2), (2,2,2)) against exact `hom_space`: degree, top scalar ≠ 0, ψ∘φ = id, braid composite = a rex move modulo lower (rank of the 1×1 form at the top coset); `certify_added` on every template up to blocks (3,3,3) | 1-2 |
| B2 | m4c-rexplan | dcoset.rs: singlestep forms, rex-chain planner (EKo Thm 5.30 constructive; CORE Thm 3.14), chart words; census G-P0 for all (3,3,3) pair targets on the N = 5/6 lists | **G-P0** (decision): per target, chain length, pad z, max flows of interned words, max residual exact unknowns, predicted wall. Go for every target within budget; the others go to §7 fallbacks. Also one `KRHOM_TRACE` run per rex-class pair, confirming X ≅ B_q (L3 §3.3) | 1 |
| B3 | m4c-pad | smap.rs: `Arena::expr` accepts bottom λ·0^z (z ≤ `MAX_PAD`); `MapNode::Pad` (identity block); arms in `add`, `children`, `assert_genuine`, `commit`, exhaustive matches (smap.rs :809, :862, :1103); smapio.rs `P` record. **D11 exception, user approval first; reviewed alone** | G-R: Pad block = identity at all cosets of random words of (3,3), (2,2,2) and (3,3,3); Pad vs Id on equal words; negative test (non-trailing zero); smapio round trip. Only if G-P0 needs a 4th upright | 0-1 |
| B4 | m4c-pairs | leaves.rs: exit/entry from rex chains; rex-class pairs | **pair 322 at N = 5**, then 71, 159, 322, 323, 325, 334, 419 at N = 6, both primes; replay | 1-2 |
| B5 | m4c-pairs | general (3,3,3) pairs: Lemma A on X with closed-form or budgeted lower steps | split class 268, 336, 361, 368; Memory 195, 196, 199, 444, 445, 460, 475; then the Refused pairs that M5c leaves (N = 5: 313, 326, 406, 410; N = 6: up to 38), both primes | 2-3 |
| B6 | m4c-run | production chain (Leaves last), full level-5 and level-6 reruns, replay, G-L5 on every Leaves map, persist to data-gl6 (D12) | **level 5 177/177 and level 6 377/377**, both primes; replay = run | 1 |

**Totals:**
- A1-A5: 5 sessions, plus 1 reserve.
- B1-B6: 7-10 sessions, plus 2-3 reserve.
- Expected: about 15-16 sessions.

**Stopping points:**
- after A3: closure-only M4c, self-contained, no frozen file;
- after A5: the mathematics is validated in all degrees at small λ;
- after B2: G-P0 decides how much of the (3,3,3) pair list M4c takes.

## 5. Build and validate on small cases first
1. **Pure combinatorics (A1, seconds):**
   - G-L1 on every canonical ladder;
   - the census with hu.py as an independent replica; the Rust census must reproduce 852 011, 392 995, 224 249 and
     406 634 (C header);
   - the degree-sign check of the sideways crossing.
2. **Exact-oracle checks at tiny λ (A2):** leaves against `sbim::hom_space` bases (stacked rank), SLL Ex 7.39 (two
   leaves, degrees ±1), Ex 7.40/7.41 (4 and 5 leaves), Ex 7.46 (10 double leaves, in A5).
3. **Known-complete sets with Leaves alone:**
   - G1(a) closure part (A2), then G1(a) 480/480 (A4);
   - B2, including (2,1,2) 605 where Search alone fails 24 and (2,2,2) 2 205 where the Chain fails 141 (A5).

   These test Lemma A, flips, degree windows, sprinkle-free selection and both certifiers before any (3,3,3) scaling
   work.
4. **(3,3,3) regression before the failures (A3):** level-4 and level-5 closure jobs with Leaves alone; the certified
   multisets must match layer A.
5. **Closed-form generators at small cells first (B1).** Every template is compared with the exact window wherever that
   window is ≤ 10⁵ unknowns, e.g. sideways crossings at blocks (2,2,2) (1 036) and (2,3,2) (3 479). Only then are
   (3,3,3)-size templates admitted, through `certify_added` alone.
6. **Rex machinery on the smallest rex-class pair first.** Pair 159 (X = 240 flows), then 71 (600), then 322 (4 320).

## 6. Interaction with M5c and M6a

**M5c (parallel engineering: calibration, memory, admission).**
- **Division of labour.** M5c owns the Search and Exact fixes listed in L3 §2.3:
  - (a) S4 unreachable in `Search::tiers(3,6)`;
  - (b) `leaf_iso` reports admission as Exhausted (homprov.rs:2371);
  - (c) the 2x retry does not raise l_max;
  - (d) Exact buckets mix closure and pair calls.

  Most small-X Refused pairs (N = 5: 313, 326, 406, 410) are M5c's. M4c touches no homprov.rs code.
- **Acceptance list.** M4c's B-phase list = the M5b list minus what M5c certifies, refreshed from M5c's run before B4.
  Fix (a) may also certify some product closure types. Stage 0 stays as an independent second route.
- **Shared files.**
  - sindec.rs `provider_for_stamped` (one insertion, A3) and libbuild.rs flags. Whoever merges second rebases.
  - Both change the provider stamp. Schedule **one** layer-A rebuild after both merges (≈ 26 min per prime warm,
    M5b_report).
  - M5c's exact-shape cost records may admit M4c step windows of 10⁵-4·10⁵ unknowns (budget rule §0.3).
- **If M5c changes Exact admission semantics,** the census cap (10⁵) stays M4c's own rule, applied in leaves.rs before
  calling `try_hom_basis_cap`.

**M6a.**
- **Not blocked:** M6a-1b (wiring, G2B) and M6a-2a/2b at N = 4 need nothing from M4c.
- **Blocked: M6a-3** (prune 6; 4_1 at N = 4, 5, 6) needs the level-5/6 library.
  - Recommended cheap step (cloud, combinatorics): the "touched set" of the M6a-3/M6a-4 braid words (4_1, 5_2, 6_2) at
    N = 5, 6, i.e. which library jobs their pruned plans reach. M4c then works on the touched failures first.
  - If 4_1 avoids pair 322 and the 12 closure types, M6a-3 can start after A3.
- **G-L5 must pass** before Leaves maps enter the library that M6a-3 evaluates.
  - The fallback is a compaction pass: replace a Leaves map by a cheaper certified map spanning the same form when one
    exists.
  - Leaves stays last in the chain, so all other targets keep their Search/Exact maps.
- **Merge order.** B3's smap.rs change must land after M6a-1b and after any M6a branch that touches smap.rs. It must
  keep A16 (`Send + Sync`, no interior mutability).
- **Cache keys.** The library hash in M6a's cache keys changes with the new stamp. Rebuild once, together with M5c.

## 7. Risks and fallbacks

| # | risk | likelihood [E] | detection | fallback |
|---|---|---|---|---|
| R1 | closed-form rex moves cannot be written in krhom's localized form, or chains are too long or too large (C1, M1) | medium-high | B1 unit tests; G-P0 | (a) prefer cell-aligned chains (CORE Thm 3.14) and accept longer chains; (b) **word-agnostic types**: sindec keyed by (q, word) so a rex-class X serves as its own representative (L3 §5; 2-3 sessions; invalidates layer A; needs the M6a contract); (c) leave the residue to M5c calibration at ≤ 4·10⁵ unknowns |
| R2 | G-L0 flags closure steps above the cap (C1 for ThreeUp) | medium-low | A1 census | closed-form steps from B1, i.e. move B1 before A3; or M5c measures those exact shapes |
| R3 | intermediate words too large (M2); FlowIndex is global (M3) | medium (pairs), low (closure) | G-L0/G-P0 flow budget; RSS in runs | prefix-extension words; minimal-flow chart choice (Z is free, SLL §7.10); `drop_index` between targets; run the worst jobs alone at 6-8 GB |
| R4 | downstream M6a cost of Leaves DAGs (M4) | medium | G-L5 | Leaves last; compaction pass; memoized per-point LocMaps in M6a |
| R5 | an Elem/Assoc constant or an Added template is wrong (M5) | low | `validate_elem`, `certify_added`, exact comparison at small labels | never record an unvalidated constant; refuse the step |
| R6 | D11 exception for `Pad` refused or delayed | low-medium | — | A-phase needs none; B4 tries 3-upright chains first; otherwise wait for the merge window |
| R7 | theory gaps [K]: rung-level combination of Thm 2.6, form rank via Thm 7.52, degree normalization | low | G-L1, V0, A4, A5 | they affect existence only; certification rejects anything wrong |
| R8 | CORE is a preprint (v2) | low | `ladder_char` reducedness on every chart and chain word | constructive use only |
| R9 | effort overrun (C M8) | high | session exits | the stopping points of §4; no B5 before G-P0 is green for the rex class |
| R10 | M5c and M4c both change stamps and sindec | medium | merge review | a single coordinated rebuild; Leaves last keeps the other maps byte-identical in route |

## 8. Honest assessment

**Probabilities** [E]:

| outcome | probability | main reason |
|---|---|---|
| closure list 13/13 (N = 5 and N = 6) after A1-A3 (+1 reserve) | **≈ 85 %** | Lemma A is verified; G-L0 is the only open quantity, and closure steps are between ≤ 4-rung words (2-7·10⁴ unknowns in C's sample) |
| G1(a) 480/480 with Leaves alone (A4) | ≈ 85 % | exact steps are affordable at small λ |
| B2 (2,2)/(2,1,2) with Leaves alone | ≈ 85 % | |
| B2 (2,2,2) 2 205 with Leaves alone | ≈ 65-70 % | large ladders at d ≥ 1 may still need closed-form steps |
| rex class at (3,3,3) (pair 322 + 7 pairs at N = 6) | **≈ 55 %** | needs closed-form localized rex moves (new derivation) and chain sizes within budget; pair 322 differs from L_q on 5 of 6 rungs on all 3 uprights |
| remaining N = 6 pairs (split 4, Memory 7, Refused residue), given the rex machinery works | ≈ 70 % | the same machinery; the risks are memory and downstream cost |

- **All N = 6 gaps closed (M4c + M5c)** ≈ 0.85 × 0.55 × 0.75 ≈ **35-40 %** in the expected 15-16 sessions, and about
  **50-55 %** with the reserve (up to about 22 sessions).
- **N = 5 complete** (closure in A3, pair 322 in B4, the 4 Refused by M5c) ≈ 50 %. Pair 322 is the single bottleneck.
- **Expected sessions:**
  - closure: 3 (range 3-4);
  - small-λ completeness: 2;
  - (3,3,3) pairs: 9-12;
  - total about 15-16 (range 12-22).

  This is above D1's 8 and D2's 9-11, in line with C M8. C's 15-25 for pairs assumed a full EKLP pair route; this plan
  narrows it to prefix-extension leaves plus exit/entry chains.

**What M4c guarantees regardless.**
- Mathematical completeness of the candidate stream (Prop 8.32, Cor 8.33, Thm 7.49) whenever every step is computable.
- Unchanged, sound certification.
- A census that says per target, before any long run, whether it fits.

It does not guarantee feasibility at (3,3,3) for pairs. That depends on the B1/B2 results (G-P0).

## 9. First actions
1. Launch A1 (m4c-census) on the cloud now. Its census and Stage 0 are useful whatever happens to the pair route.
2. Ask the user now whether a D11 exception for an additive, reviewed `MapNode::Pad` in smap.rs is acceptable after the
   M6a-1b merge (B3). Nothing before B3 needs it.
3. Ask M5c to publish its refreshed failure list and any new exact-shape cost records. They define the B-phase
   acceptance set and the 10⁵-4·10⁵ admission band.
4. Request the M6a "touched set" for 4_1 / 5_2 / 6_2 at N = 5, 6, to order the B-phase targets.
5. Before B1, derive by hand one sideways crossing (blocks (1,1,1) and (2,1,2)) and one switchback (EKo eq. (108), SLL
   Ex 6.1) in krhom's localized form, and check them against exact windows on the cloud (≤ 10³ unknowns). This is the
   first test of R1.
