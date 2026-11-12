# C_m4c: adversarial critique of the M4c designs D1 (EKLP, Lemma A) and D2 (light, padded EKLP)

Date 2026-11-11. Read-only study. Nothing in the repository was compiled or run; no cargo, no binaries.

**Read in full:** D1_eklp.md, D2_light.md, L1_eklp.md, L2_alt.md, L3_code.md.

**Checked against:**
- the SLL TeX (papers/eklp_2401.03053/src): Thm 2.6, Prop 4.19-4.21, Lemma 8.2, Lemma 8.7, Prop 8.26, Thm 8.31, Prop 8.32, Cor 8.33, Thm 8.34, eqs. (84) and (91)-(92), the §7.7 footnote;
- CORE Lemmas 3.12 and 3.13 (papers/eklp_2402.08673/src/typeA.tex);
- QR Prop 3.15 (papers/qr_foam_1405.5920/paper.flat);
- dev 277a21c: smap.rs, sexpr.rs, coset.rs, tmpl.rs, homprov.rs, tests/m4a.rs;
- M5b data (m4c/data: costs_333.txt, run_n6.jsonl).

**Numerical tool (new, laptop-safe, < 1 s each):** m4c/c_tools/hu.py is a line-by-line replica of `homprov::hom_unknowns` (homprov.rs:344), `fdeg_hist` (:280), `word_shift` and `schur::emono_count` (schur.rs:1684). It reproduces the measured M5b numbers exactly:
- 852 011 (Hom(∅, E0.3 F0.3 E1.3 F1.3), d = 18);
- 392 995 (L_p of 1,0,2/2,0,1/0,3,0 → ∅, d = 13);
- 224 249 and 406 634 (L_p of 1,1,1/1,1,1/1,1,1, d = 4 and 6).

Helpers: braid.py (atomic braid windows), far.py (far commutations), yall.py (cell-form flows of every (λ,μ) table, 3 and 4 columns) and anti.py (types whose braid window is infeasible).

Tags: [V] verified here (paper text, code or computation); [K] inference; [E] estimate.

## 0. Verdict

1. **Mathematics.**
   - Both designs rest on proved results: SLL Prop 8.32 and Cor 8.33 for closure; Thm 7.52 plus the standard local-intersection-form argument for pairs. Neither needs the unpublished generators-and-relations presentation: krhom evaluates bimodule maps, so SPEC's "companion paper" risk does not apply [V].
   - D1's Lemma A (§1.3) is correct as stated [V, §2.1 below]. D2 follows SLL Def 7.28 almost verbatim, and its rex substitutes are admissible (Prop 8.26 uses only "degree 0, 1⊗ ↦ 1⊗") [V].
   - Certification is unchanged, sound and sufficient: wrong candidates fail visibly.
2. **The central engineering claim is not established for either design.** The claim is that every exact call becomes a small local window, so that Refused and Memory disappear and every target can be certified.
   - Both designs keep exact `hom_space` windows whose size has no bound by construction.
   - The replica shows that representative windows of both designs are 10⁵-10¹⁰ unknowns. The measured wall is ≈ 4-5·10⁵ unknowns: (3,3,3) windows of 108×120 flows at d = 0 have 437 311 unknowns and time out at 60 s with 2.7 GB (costs_333.txt).
   - The completeness argument (D1 §1.4, D2 §3) therefore reduces to "complete if every local window is admitted". That is the same failure class (Refused) as today.
3. **Smaller issues.**
   - Downstream (M6a) evaluation cost and memory of the padded intermediate words are not budgeted.
   - Genuineness certificates for E2/E3/E4 leaves are only test code today.
   - Effort is underestimated by about 2×.
4. **Recommended shape** (§3):
   1. D2's Stage 0 now.
   2. A census run with this replica.
   3. Closure through D1's Lemma A on 3 uprights, with prefix-extension junction words.
   4. Pairs through D1's rung-level recursion, built from closed-form generators only (D2's ELL; switchback-composed rex moves after CORE Lemma 3.12). No exact window above about 10⁴-10⁵ unknowns, by construction.
   5. D1's trailing-zero `Pad`, not D2's `Repad`.
   6. A downstream gate before production.

## 1. Ranked findings

### Critical

**C1 (both). Exact windows of unbounded size remain, so "no window spans 9 roots" / "Refused and Memory disappear" is false.**

- **D1, closure on 3 uprights (`Chart::ThreeUp`).** Every step S_{t'→t} is an exact call Hom^d(W_{t'}·ρ, W_t) between whole 3-upright words on all 9 roots (D1 §2.3 `step_family` item 3).
  - It is cheap when W_t is short: Hom(E0.2 F0.2 E1.3 E0.1, E1.3 E0.1) has 2.1-7.4·10⁴ unknowns for d = 0..4.
  - It is not cheap when the step is rex-type, that is, when W_{t'}·ρ is itself reduced for t but differs from the fixed W_t: two 90-flow reduced 4-rung words at d = 0 give 141 831 unknowns, and two 150-flow words give 1 119 356 [V, hu.py].
  - D1 §6 estimates "typically 10³-10⁴ unknowns, 0.5-20 s" without a census.
- **D1, `Chart::Factored`.** D1 §1.7 and §6 say the K windows are "the measured class of (2,2,2) and (3,3) windows, 10²-10⁴ unknowns". This is false whenever the third column is empty: K then holds all 9 roots and the local bottom is (3,3,3).
  - It happens in the failing words [V, by levels]: E0.1 at level (3,6,0) in L_p of 0,1,2/2,0,1/1,2,0 (E0.2 F0.2 E1.3 E0.1 F1.3 F0.1), and E0.2 at level (4,5,0) in pair 334.
  - Those local windows belong to the (3,3,3) class, where 10²-flow windows already time out (costs_333.txt).
- **Both: atomic braid moves as exact degree-0 windows.** These are D1 §1.8(c) and D2 §2.4 and §2.6. A braid move is H₀H₁H₀ ↔ H₁H₀H₁ on cells (x,y,z), with H = E_u^(b)F_u^(a). Unknowns at d = 0 [V, braid.py]:

  | cells | unknowns at d = 0 |
  |---|---|
  | (1,1,1) | 111 |
  | (1,2,1) | 946 |
  | (2,1,2) | 18 019 |
  | (2,2,2) | 711 580 |
  | (3,1,3) | 2 966 475 |
  | (2,3,2) | 12 609 774 |
  | (3,3,3) | 1.6·10¹⁰ |

  - In a 3-column table the only pairwise-inverted triple is the anti-diagonal, so the braid window is fixed by the type. 6 of the 55 (λ,λ) types (for example 0,1,2/1,2,0/2,0,1 and 0,0,3/0,3,0/3,0,0) and 36 of the 811 3-column (λ,μ) tables have an infeasible braid window (> 4·10⁵ unknowns) [V, anti.py]. 4-column levels have more such triples.
  - Far commutations are also non-trivial: cells (3,3,2,1) give 192 484 unknowns [V, far.py].
- **D2, rex fallback.** The fallback "one exact degree-0 window on the differing cells (≤ 6 cell uprights)" (D2 §2.6, R3) has the same problem: its worst case is 9 roots on 6 uprights, not "a regular S₆ piece of 10⁴-10⁵".
- **Fix.**
  - Never solve for rex moves or braid moves exactly above about 10⁴ unknowns. Compose them from closed forms instead:
    - CORE Lemma 3.12: switchback ∘ (up-up, down-down) ∘ switchback, passing through the merged expression [I+s_a+s_c−s_d−s_e];
    - CORE Lemma 3.13 for far commutation;
    - switchback maps φ and ψ from SLL Thm 6.7 (§6.1), built from right- or left-facing crossings and ccw cups and caps. Both have closed forms: eq. (84) defines the left crossing by a twist, and eqs. (91)-(92) give ∂Δ.
  - D2's generator library already contains these pieces; D1's has none.
  - Then make the census (§3 step 2) a hard gate that bounds every exact call.

### Major

**M1 (D2). The cell rex normalizer is incomplete, so "no general rex planner is needed" (D2 §0.1, §2.6) does not hold.**
- The previous target in Def 7.28 is Z'∘X_{n'}, of the form split → cross → merge.
- The next source is Z∘X_m (SLL Notation 7.1 / 7.22), of the form [[λ⊃LR(z)]] ∘ atomic(ψ(z)) ∘ [[L⊃LR(m)]] ∘ core(m) ∘ [[RR(m)⊂J]]: split, cross, split, cross, merge.
- Turning one into the other moves splits through atomic crossings. These are fork slides and type-A switchbacks (EKo eq. 108), not up-up/down-down (E4) or Coxeter moves on ψ.
- D2 handles these only by the exact fallback, which hits C1. D1 at least lists them as rule (b).

**M2 (D2). Intermediate targets at −s levels can be 10⁶-10⁷ flows.**
- D2 recurses per EKLP single step. After every −s step the target is a reduced expression of a (λ,ν)-coset with ν of 4 blocks.
- The cell forms of such cosets reach 13 436 928 flows (ν = (3,3,2,1)), and 990 of the 4 380 tables exceed 10⁵ [V, yall.py].
- D2 §7 budgets at most 104 976 flows (944 784 only for the decomposable type). At about 40 padded rungs, the FlowIndex costs 8·40 B per flow, so 3+ GB per expression, and `FlowIndex::build` (sexpr.rs:130-175) holds `flat` and `masks` at once, so the build peak is about 2× that.
- D1 avoids this: its steps are whole rungs (G_{t'→t} summed over the middle u), so its junction words live only at 3-column levels, where the cell forms are at most 419 904 flows.

**M3 (D2). Mitigation R2(d), a lazy per-coset FlowIndex "in sexpr.rs (not frozen) behind the same methods", is not a local change.**
- `Arena` is documented as "Send + Sync … no Rc, RefCell or Cell" (smap.rs:496).
- `Eval` takes `&Arena` and calls `ar.flows(id)`, which panics unless the index is built (smap.rs:633).
- `Eval::block`, `kron_block` and `window_block` (:1092-1250) use global flow ids (`range`, `flow`, `find`).
- So laziness needs interior mutability, which is a redesign of smap.rs. Without it, every expression in a certified DAG must have its full index resident while `certify_*` evaluates it, for all cosets, even though only one coset is read.

**M4 (both). Downstream M6a cost is not budgeted.**
- M6a evaluates every library ι/π at all cosets, at every point (R4_m6a.md §4: t_lib ≈ 10 s and M_ev ≈ 0.8 GB per point at prune 6).
- A Leaves map routes through intermediate words of 10⁵ flows (D1 charts) up to 10⁶-10⁷ flows (D2), against X ≤ 21 600 flows (run_n6.jsonl). All their indices must stay resident.
- D1 states 5-20× (G-L5) without a basis; D2 says nothing.
- D1 also puts Leaves before Exact in production (§2.5). Targets that Exact certifies today with 5-20-node maps would then get 50-300-node Leaves maps.
- Keep Leaves last, or use it only for targets the chain fails. Make G-L5 (t_lib and resident index bytes on the first certified closure types) a gate before the pair stage.

**M5 (both). Genuineness of Elem/Assoc leaves is a caller-trusted record, validated only in tests.**
- `certify_elem` and `certify_assoc` store whatever scalar they are given (smap.rs:821-833). The only validation is exact 1-dimensional `hom_space` comparisons in tests/m4a.rs: A7 at :760-790 and A8 at :918-960. No provider uses Elem today (L3 §1.3).
- Leaves would be the first production caller, so the validation must become reviewed production code with size caps:
  - E2 (4,5) costs 112 390 unknowns, which is acceptable;
  - an Assoc that reorders a 9-root split into (3,3,3), for example F0.6 F1.3 ↔ F0.3 F1.3 F0.3, costs 2.7·10⁷, which is infeasible [V, hu.py].
- Restrict Assoc to merge-only windows (1 unknown) and small splits, or add a structural certificate. Never record an unvalidated constant.

**M6 (D1). The factored step "S = id_Z ⊗ s" (D1 §1.7) is under-specified with trailing-zero padding.**
- Comp, Kron and certification need equal tops. So the chart word C_K(t) must end at μ_k·0^z, with its 4-block level L_K⁺ merged back and blocks moved back from the pad uprights.
- That leaves trivial "unspread" rungs between T'_K and the next rung ρ, so the window is not contiguous.
- Each step is then conjugated by spread/unspread rungs and far commutations: Assoc or commutation windows, each needing its own certified constant (M5). D1 does not design this layer.

**M7 (D2). `Repad` enlarges the trusted base.**
- `Repad` is about 150 lines of coset translation and flow reordering in the frozen smap.rs (D2 §5.2). It runs in every downstream evaluation.
- Certification tests ranks and degrees at rep(q) or the identity coset only. A translation error at other cosets would be a certified wrong map, invisible to the D14 scaling check and to two-prime agreement.
- G-R (Repad∘Pad f = f at random cosets and points) must be a hard gate over all cosets. Comparisons with layer-A (3,3) and (2,2,2) maps are also needed.
- D1's trailing-zero `Pad` is the identity: `CosetSpace::rank` is invariant under trailing zeros, with radix C(0,0) = 1 (coset.rs:116-140) [V]. QR Prop 3.15 also pads with trailing zeros [V]. Prefer D1's `Pad`.

**M8 (both). The effort is underestimated.**
- D1 plans 8 sessions, but S5 alone holds factored steps, exit and entry leaves, the rex planner and `pairs()`.
- D2 plans 9-11, but sessions 4-5 hold:
  - Added sideways-crossing templates written as PMats in sbim coordinates from the Sweedler formulas (91)-(92);
  - the ELL compiler (Lemma 8.2), including zero-upright placement;
  - sprinkles, the LL recursion and the normalizer.
- SYNTHESIS gave 10-20 days for M4c. A realistic estimate is 3-5 sessions for closure and 15-25 for pairs, conditional on C1 being fixed.

### Minor

- **m1 (D1 §1.4).** "Localization at rep(t) is injective on Hom_{≮t}" is true, but its citations are wrong.
  - Prop 4.19 covers only degree-0 maps between reduced expressions of the same coset; W Thm 7.3.5 is a stalk statement.
  - The correct argument is the proof of Thm 8.31: Prop 8.26 together with Prop 4.20, im f ∩ R^{LR(p)}·1⊗ = 0. The localization of 1⊗ at the single top flow is nonzero.
- **m2 (D1, 3-upright closure).** Lemma A needs a 3-upright word with exactly one flow at type t for every junction type.
  - Types without one (likely the decomposable 0,1,2/1,1,1/2,1,0 and its (3,4,2) analogue [K]) break Lemma A under the "fallback canonical_ladder".
  - G-L0 must list them and route them through Factored.
- **m3 (D1).** Choose W_t := W_{t*}·ρ whenever that word is reduced for t. Lemma A allows any word, so the top steps become `Id` and the rex-type step calls of C1 disappear.
- **m4 (D1).** Pruning step families by degree window must still run greedy selection over all lower degrees of supp G. Otherwise an R_t⁺-multiple of a lower step can be selected. This would not be a wrong answer: it shows up as a rank shortfall.
- **m5 (D1 §1.9).** The double-leaf candidates omit left R^λ monomials (Dot at level 0). `certify_full_plan` counts k-dimensions at a fixed degree (SLL Rem 7.48/7.50); D2 §2.10 includes them.
- **m6 (D2 §2.7).** "A flip is D(leaf) up to one scalar" is false when a window's Hom space has dimension > 1. It is the flip of another admissible choice (Thm 8.34, Cor 8.33), so ranks are unaffected; this is a wording fix.
- **m7 (D2 §2.2).** The fixed layout (row i on uprights 3i..3i+2, column j on 3j..3j+2) does not cover 4-block levels. A pure split needs an empty neighbouring upright, and none is free when 9 blocks sit on 9 uprights. A zero-placement normalization (Assoc windows, see M5) is needed.
- **m8 (D2 §5.2).** smapio.rs is not frozen: D11 lists homprov, g1, sbim, linsys, smap and tmpl (SPEC_3C_GL6_rev1 §D11). Only smap.rs needs the exception. Adding a MapNode variant touches the exhaustive matches in smap.rs (:809, :862, :1103) and smapio.rs:215. homprov.rs:953 has a wildcard arm.
- **m9 (both).** The reducedness of Y_q and of the chart words comes from CORE Thm 3.11. CORE is an arXiv v2 preprint, used here only constructively. Gate it cheaply with `ladder_char` (c_t(1) = 1 and support ≤ t), as D1 §1.1 already does for junction words; D2 has no such gate for Y_q.
- **m10 (D2 §0.2).** QR Prop 3.15 says fullness holds after adding zero uprights for some m′ ≥ m [V]. That is existence, not necessity, so "padding is unavoidable" is a judgement, not a theorem.
- **m11 (notes).** L2 §0 and §2.8 say no type-A cores companion has appeared and that the "EKLP risk stands". That contradicts L1 §0 and §4: CORE 2402.08673v2 is in papers/ and both designs use it.
- **m12 (both).** The number of leaves is not the cost driver:
  - closure: n = c_diag(b_p)(1) ≤ 54, which equals the Hom rank;
  - D2 pairs: [v^{∓s}] c_q(char X), small and ≥ mult;
  - D1 pairs: D1 builds full junction families, Σ_t c_t(X_{≤k})(1) maps per level (hundreds to thousands for X = 21 600 [E]). It should report this in G-L0. The Eval memo holds a block per (map, coset).

### Verified strengths (keep)

- **D1 Lemma A.**
  - W_{t'} has flows only at types ≤ t', so the matrix is block lower-triangular in a linear extension of the Bruhat order.
  - On the diagonal, after multiplying by blockdiag(ℓ(c'))⁻¹, it is T ⊗ I.
  - Squareness follows from flows at a coset = c_t(1). This is checked in code by `deodhar_check` (sexpr.rs:267).
  - The existence argument (§1.4) holds: the top-path subfamily of the Prop 8.32 basis has an invertible top-to-top block, because leaves with lower prefixes vanish there.
- **Closure.** The leaves ending at the identity coset form an R^λ-basis of Hom(L_p, R^λ). The identity coset is Bruhat-minimal, so there are no lower terms (Prop 8.32 and the §7.7 footnote) [V]. This matches the determinant argument of `certify_closure` (homprov.rs:945-985).
- **D2 Stage 0 (product rule) is sound.**
  - The identity-coset matrix is the Kronecker product of two nonsingular matrices when n_p = n_A·n_B.
  - It is cheap: a 20×20 commutation window.
  - It needs neither new mathematics nor a D11 exception, and it certifies 5 of the 12 closure types. **Do it first.**
- `certify_pairs` and `certify_closure` (homprov.rs:790, :945) are unchanged and remain the only acceptance criterion. Both designs respect the arena mark/commit discipline.

## 2. Answers to the review questions

| question | D1 | D2 |
|---|---|---|
| Basis of the right Hom space, ring, degrees, singular? | yes: R^λ-basis for closure in all degrees; R_q-basis modulo lower for pairs; Lemma A is new but proven [V] | yes, SLL verbatim (Prop 8.32, Cor 8.33; Thm 7.52 for pairs [K, standard]) |
| Unpublished results? | none essential; uses CORE (preprint) for charts; [K] rung-level Thm 2.6 | none essential; CORE Thm 3.11 for Y-forms; [K] padding and flip |
| Implementable in krhom? | closure yes (no new nodes); pairs need `Pad` (D11) plus an undesigned spread/unspread layer (M6) | needs `Repad` (D11), sideways-crossing PMats, an ELL compiler, zero placement (m7), a planner (M1) |
| Hidden costs | 3-upright and K = 9 windows, braid windows (C1); chart indices; M6a cost (M4) | braid templates and fallback (C1); 4-block −s levels up to 1.3·10⁷ flows (M2); laziness impossible as stated (M3); M6a cost (M4) |
| Certification soundness | sound; new trust: `Pad` (trivial), Elem/Assoc constants (M5) | sound at the family level; new trust: `Repad` translation (M7), Elem/Assoc constants (M5) |
| Effort | 8 sessions → realistically 3-5 (closure), 15-25 (pairs) | 9-11 → 1 (Stage 0) + 15-25 |

## 3. Recommended M4c path

1. **Stage 0 now** (D2 §2.1): the 5 product closure types, both primes. One session, no D11 exception.
2. **Census G-L0, extended, before any generator code.** Combinatorics only, using hu.py or its Rust twin. For each of the 68 failures, list:
   - junction types per level, with G and word lengths;
   - `hom_unknowns` of every 3-upright step call, using prefix-extension junction words (m3);
   - every K window, whether the third column is empty, and the cells of every braid triple;
   - the flows of every chart word.

   Hard rule: no exact call above about 10⁵ unknowns (the measured wall is 4-5·10⁵).
3. **Closure:** D1's Lemma A on 3 uprights for the 7 non-product types, if the census passes. Otherwise use step 4's generators.
4. **Pairs:**
   - recursion at rung level (D1), so no 4-block word is ever materialized;
   - trailing-zero `Pad` (D1);
   - step maps from closed forms: E2/E3 (existing), E4 restricted to merges, sideways crossings (eq. 84 twist, or (91)-(92) as Added templates) and ELL sinister diagrams (D2);
   - rex and braid moves composed from switchbacks (SLL Thm 6.7, CORE Lemmas 3.12 and 3.13), not solved for;
   - exact windows only where the census shows ≤ 10⁴ unknowns.
5. **Gates:**
   - G-L5 (downstream t_lib and resident index bytes) after the first certified closure types;
   - production validation of every Elem/Assoc constant (M5);
   - Leaves last in the chain.
