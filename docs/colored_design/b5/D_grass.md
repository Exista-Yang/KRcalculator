# D_grass: B5 through the rung-local Grassmannian factorization of the Lemma A steps

Date 2026-11-12. This is a read-only design against `m4c-rexplan` (4e886f6, worktree `revrex`). Nothing was run in Rust.

**Inputs:**
- `B5_census_analysis.md` (this folder). Its Python replica reproduced `step_score` on all 111 non-rex targets.
- Design and review notes: D1 §1.3-1.8 and §2.4; L1 §1-3; A1 (G-L1); B1; B2 (open problems 1 and 2); C_m4c (C1, M4, M5, M6, m1, m3).
- Code: `src/{leaves,dcoset,leafgen,rexplan,smap,homprov}.rs`, `src/bin/gp0.rs`.

**New numbers in this file** come from `tools/dgrass_win.py` (same replica, < 1 s).

**Tags:** [P] = paper, [V] = computed with the replica, [K] = inference that a gate checks, [E] = estimate.

## 0. Verdict

**Go, behind gates, with no `Pad` on the primary route.**

- **What to take from the Grassmannian factorization.** It tells us how to build the steps, but its literal chart Z·T_K is the wrong place to lay out words.
  - 40 of the 42 distinct worst steps have 3 K cells and a non-empty other column, so Z·T_K needs a 4th upright.
  - The detours W ↔ Z·T_K on λ·0 then cost 10⁶-10⁹ unknowns as single calls. Only 313, 212, 281 and 381 fit under 10⁵ [V].
- **Projected chart.** B5 keeps the factorized form S = id_P ⊗ s on 3 uprights. P is Z with its pad rungs dropped, and the joint word search finds it.
- **Window certificate (Prop. W, §1.3).** This replaces the missing theorem. It is four combinatorial checks at one configuration. When they pass, the window step family is complete a priori, by SLL Prop 8.32 applied locally. That settles B2 open problem 2 for every window that passes, before any exact call is made.
- **Predicted coverage, worst step only:**

  | groups | targets | how they are covered |
  |---|---|---|
  | A, B, D1, D2 | 74 of 101 | windows of at most 78 371 unknowns, Pull/Push templates, or two records |
  | E | 8 | rex chains, no step call |
  | C | 19 | need a dry rank test (windows ≤ 6 320) |

- **Sizes [V].** Every chosen step window has ≤ 720 flows per side and a maximum label ≤ 6, so no GrTable(·, 9) build is needed (that build was the B1 blocker). Every junction word has ≤ 38 880 flows, against the 10⁶ budget.
- **`Pad` is a narrow fallback.** The only padded step that is fully within budget and not trivial is 313 (local window 484, detours 51 496 / 19 764). Build `Pad` only if 313 or 410 is still the last N = 5 blocker after B5c.
- **Honest outlook [E]:**
  - all of N = 5 complete: 45-60 %;
  - N = 6: about 19-25 of the 32 non-feasible pairs;
  - 3-4 cloud sessions, plus 1-2 if `Pad` is needed;
  - stop at the end of B5c.

## 1. Construction

### 1.1 Rung-level reduction (what is proved, and where)

Take one Lemma A step t' →ρ t at level k, with ρ on uprights K = (i, i+1) (M4C_PLAN §3.1).

1. **Reduction.** `dcoset::grass(t', t, ρ)` returns:
   - Q = merge_columns(t', i), which also contains t;
   - L = RR(Q);
   - z with z_ = Q_;
   - the local tables m' and m'' on `local_lam`, the K cells (c₀, c₁, c₂).

   Source: SLL Thm 2.6, Prop 2.9, eq. (11) [P].
2. **Same generating function.** G_{t'→t} = G_{m'→m''}, with the right side from `step_gf` on the local tables. This is SLL Lemma 3.9 for one single step [P]. At the rung level (a split followed by a merge, with the same z) it is [K], but A1 G-L1 verified it for all 28 000 (3,3,3) steps, together with "t_ = Q_·m_, lengths add".
3. **Factorized light leaves.** A single-step light leaf is SSLL = id_Z ⊗ ELL([m,n]), with X_p = Z∘X_m (SLL Notation 7.22, Def 7.23, Lemma 7.24).
   - ELL is the sinister diagram of Def 7.3, built by Lemma 8.2; its flip is eq. (133).
   - Polynomials enter only as sprinkles (Def 7.33).
   - So every step family can be taken of the form id_Z ⊗ s, with s ⊂ Hom(T(m')·ρ, T(m'')) living on the K cells.

### 1.2 Three charts for the prefix, and the joint words

| chart | applies to | prefix | uprights | worst steps |
|---|---|---|---|---|
| true, 3 uprights | K column of Q has ≤ 2 cells (D1 class (a)) | Z, a row swap E1.3·F1.3 | 3 | 2 (281, 381) |
| projected, 3 uprights | 3 K cells and a non-empty other column (class (a-pad)) | P = Z minus its rungs on the pad upright; the other column rides with the last K cell (labels (c₀, c₁, c₂+o)) | 3 | 40 |
| padded, 4 uprights | fallback only (§1.6) | Z·T(m)·U on λ·0 | 4 | – |

Class (b), where K holds all 9 roots (C_m4c C1), does not occur among the worst steps [V].

**Joint words (A1 open problem 1).** Write W'_{t'} = P·A' and W'_t = P·B', where A' and B' are words on the labels ν after P.
- `leaves::joint_pair` searches pairs from:
  - reduced words of t' (≤ 6 rungs) × reduced words of t (≤ 7 rungs), indexed by a prefix trie so that only pairs sharing a prefix are scored;
  - seeded words P·A, with suffixes of ≤ 4 or 5 rungs.
- Ranking: W-certificate first, then `step_score`, then flows.
- For a final step, W'_q is free, because the exit chain absorbs it.
- With "G_loc = G" as a filter, this reproduces the census-analysis windows: 30 steps have such a window ≤ 10⁵ [V].

### 1.3 Proposition W: completeness of a window step family (the completeness argument)

**Notation.**
- c = rep(t).
- f* = the single flow of W'_t at c. π* is its P-part and c_loc its window part, which is a configuration of the window problem on ν.
- Φ = the top-to-top columns: flows of W'_{t'}·ρ at c whose configuration after |W'_{t'}| rungs has type t'. `step_family` already asserts |Φ| = G(1).

**Checks** (all combinatorial):
- **(W1)** A' is reduced for a local type a', and B' for a local type b. Test with `dcoset::reduced_for` on `Schur::new(ν)`.
- **(W2)** b occurs in step_gf(a', ρ_loc), and G_{a'→b} = G_{t'→t} as Laurent polynomials.
- **(W3)** The local type of c_loc is b.
- **(W4)** Every f ∈ Φ has P-part π*. The map f ↦ (window part of f) is a bijection from Φ onto Φ_loc, the local flows of A'·ρ at c_loc whose configuration after |A'| rungs has type a'.

**Prop W.** Assume (W1)-(W4). Take k-bases of Hom^d_ν(A'·ρ, B') for d ∈ supp G, embed them as `Local` nodes id_P ⊗ s, and run greedy selection by ascending degree. The result has degree multiset G and an invertible top-to-top matrix at rep(t). That is exactly hypothesis (ii) of Lemma A (M4C_PLAN §3.1).

*Proof.*
1. **Local completeness.** Apply D1 §1.4 to the window problem: expression A'·ρ, target b, reduced word B'.
   - Prop 8.32 gives an R^{LR(b)}-basis rex∘LL of Hom_{≮b}(A'·ρ, B').
   - Localization at a type-b configuration is injective on this basis, by the proof of Thm 8.31 (Prop 8.26 + Prop 4.20; this is the corrected citation of C_m4c m1).
   - The leaves through the top path of A' are unique by (W1). They have degrees G_{a'→b} (Lemma 8.7) and an invertible top-to-top block. All other leaves vanish on these columns (Lemma A triangularity).
   - By (W3), c_loc = w·rep_loc(b) for some w ∈ W_ν. Genuine maps are R^ν-bilinear, so the block at c_loc is the w-translate of the block at rep_loc(b) and stays invertible over k(y).
   - So greedy selection finds n = G_loc(1) maps.
2. **Transfer to the full step.** At c, a `Local` node's block is block-diagonal in prefix flows (smap.rs `Local`: id_P ⊗ tmpl ⊗ id_Q).
   - The target has only the flow f*, so only columns with P-part π* can be nonzero.
   - By (W4), these columns are all of Φ, and they match Φ_loc one to one.
   - Hence T_S(c) = T_s(c_loc), evaluated at the point permuted by π*.
3. **Degrees.** The degree of id_P ⊗ s is the degree of s, and G_loc = G by (W2). □

**Remarks.**
- **True chart.** Take P = Z, A' = T(m'), B' = T(m''). Then (W1)-(W4) are the rung-level content of Thm 2.6 + Lemma 3.9: z_ = Q_, and ρ preserves the Q-merged configuration, so every column in Φ passes through the top flow of Z. This is [P] per single step and [K] at the rung level, where A1 verified it. B5 never relies on this: it runs the checks for every step it uses.
- **Projected chart.** Here the checks *are* the proof. (W1)+(W2) ("G_loc = G") hold for 30 of the 42 chosen windows [V]. (W3) and (W4) are new and have not been checked [K]; B5a computes them.
- **Cost.** One `FlowIndex` of W'_{t'}·ρ (≤ 38 880 flows [V]) plus local Deodhar on ν: milliseconds.
- **A window that fails W is not wrong**, only not known to be complete. It goes to the dry rank test: `step_family` with Whole disabled, one exact call ≤ 10⁵ at one point.
- **The B2 caveat.** Pair 313, level 5: a 137 964-unknown window fell back to a 1.96 M whole-word call. Prop W is designed to flag that kind of failure before any exact call; V1c checks that it does.

### 1.4 Connecting to the planned words, and rex-type steps

- **Adopt.** Use W'_{t'} and W'_t as the junction words whenever every step touching t' or t (at both levels, all incoming t'') stays within budget. This is the spirit of C_m4c m3. A final W'_q is always adopted, exactly as `alt_final` does now.
- **Detour.** Otherwise:
  - Out: S = R_out ∘ (id_P ⊗ s) ∘ Kron(R_in, [ρ]).
  - In: S^∨ = Kron(R'_in, [ρ]) ∘ (id_P ⊗ s^∨) ∘ R'_out.
  - Here R_in: W_{t'} → W'_{t'} and R_out: W'_t → W_t. The primed R' are chains planned in the opposite direction, at types t' and t (`rexplan::Planner::plan`).
  - **Why Lemma A still holds.** Each R has degree 0 and sends 1⊗ ↦ 1⊗ between reduced words of one type (SLL Thm 6.11, Prop 4.19). So at every configuration of that type its block is the scalar 1 on the unique top flows. `Kron` evaluates R_in at each column's prefix configuration. Hence T_S(c) = T_{id_P⊗s}(c).
- **Rex (group E).** When G = 1 and W_{t'}·ρ is reduced for t, s is the planner chain W_{t'}·ρ → W_t. Its 1×1 top-to-top matrix is 1 (Thm 6.11). In SLL terms this is the forward-step identity "ELL([n,m]) ∘ rex = id" (Ex 7.18, Lemma 8.14).
- **Detour costs [V].**
  - Only 10 of the 42 detour pairs fit a single Exact call ≤ 10⁵; 21 have an in-side ≤ 10⁵ and 14 an out-side ≤ 10⁵.
  - The rest need planner chains. B2 found all 118 exit chains on 3 uprights, with ≤ 6 moves, Exact ≤ 99 175 and templates ≤ 90 flows.

### 1.5 Realizing s (try in this order; the first that passes is used)

1. **`Id` or `Rex`** (§1.4).
2. **Closed form**, when G = v^d with coefficient 1:
   - `CfKind::Pull` for Out with a flag-type source window, or `Push` with a flag-type target (B1);
   - degree check d = shift difference (`ClosedForm::new`);
   - the top-to-top entry must be nonzero: 1 for Pull, e_X/e_Y for Push;
   - genuineness is decided by `TmplLib::certify_added`.

   Candidates: group A (2 steps), D1 (3 steps) and 8 B steps. For D1 this replaces a 1.9-3.0·10⁵-unknown exact call by templates of 80×32, 90×90 and 100×120 flows, with maximum label 6 [V]. Their time is estimated with `rexplan::predict_template_secs`.
3. **Exact window** ≤ 10⁵, or ≤ 4·10⁵ with an uncensored exact-shape record (`Leaves::admit`), through `try_hom_basis_cap` and the existing greedy selection.
4. **Not on the B5 path:** closed forms for multi-term G. These would be ELLP: Pull/Push ∘ `Dot` by the sprinkles f_j = ∂_{L_j}∂_{y_j}(P_{R_{j+1}}) (Def 7.33, Lemma 7.34) [K]. Every multi-term window is ≤ 78 371 exact anyway, except D2, which uses records.

### 1.6 Padded chart with trivial rungs (fallback; fixes C_m4c M6)

- **Setting.** Bottom λ·0 = (3,3,3,0). The chart word is C(t) = Z·T(m)·U.
  - Z is reduced for z·0: `bfs_reduced` on gp0's `pad_table`.
  - T(m) is the local word on the 3 K uprights.
- **K = {1,2}.** Layout (o | c₀ c₁ c₂) on uprights 0 | 1 2 3. T(m) empties upright 3, so C(t) already ends at μ·0 and **no unspread rung is needed**.
- **K = {0,1}.** Layout (c₀ c₁ c₂ | o).
  - T(m) empties upright 2, and U = E2.o moves o back. U is a trivial rung: one flow per configuration, degree 0.
  - Then W·ρ = Z·T'·U·ρ. One far-commutation `Swap` (B1 `CfKind::Swap`: degree 0, permutation localization) turns it into Z·T'·ρ·U.
  - s is then a single contiguous `Local` node at at = |Z|, with U as its suffix (smap `Local`: src = P·A·Q).
  - This is exactly M6's objection to D1 §1.7 (the window was not contiguous), resolved.
- **Step.** S = Unpad ∘ R4_out ∘ Local(s) ∘ Local(Swap) ∘ Kron(R4_in, [ρ]) ∘ Pad. The R4 chains live on λ·0 (Planner, gp0 `pl4`).
  - Prop W holds, by the true-chart remark in §1.3.
  - Swap and Pad have top scalar 1.
  - It needs `MapNode::Pad` as in D1 §2.4 (identity block, since `CosetSpace::rank` ignores trailing zeros). This is the approved edit "only if strictly needed".
- **Costs [V].**
  - Local windows: 45-57 666 unknowns, median 1 274.
  - Single-call 4-upright detours ≤ 10⁵ only for 313 (51 496 / 19 764), 212, 281 and 381.
  - All others are 10⁶-10⁹ and would need a planner on λ·0, which has never been run.

## 2. Coverage by census group (worst step; window unknowns are the max over supp G and both directions)

| group | targets (pairs) | chart and s | max exact for s | templates | connection | completeness | fallback |
|---|---|---|---|---|---|---|---|
| A | 2 (281, 381) | true chart, 3 uprights; exact window (Pull/Push optional) | 4 603 | Pull/Push on [3,3], 20×30 flows | single Exact detours 3 184 / 34 189 | Prop W (true chart) | – |
| B | 57 (24, 45, 78, 82, 127, 137, 145, 212, 227, 273, 278, 283, 285, 292, 293, 326, 335, 336, 376, 383, 385, 411, 418) | projected; exact window | 78 371 (326, level 5); ≤ 57 736 for the others | Pull/Push on 8 steps (optional) | adopt; otherwise planner chains (single-call bounds up to 9.2·10⁸) | Prop W if W3 and W4 pass | dry rank test; next candidate; record |
| C | 19 (313, 410 at N = 5 and 6; 137, 145, 283, 293, 383) | projected; window not W-certifiable (G_loc ≠ G) | 6 320 | – | as B | dry rank test in B5a | W-certifiable windows with a record: 283 842 (313), 318 975 (137/145); padded chart for 313 |
| D1 | 6 (82, 285, 385, 458) | projected; Pull/Push | 0 (if `certify_added` fails: record of 194 000-295 452) | 80×32, 90×90, 100×120 flows, max label 6 | as B | Prop W, top entry ≠ 0 | record |
| D2 | 9 (326 at N = 5 and 6; 385, 467) | projected; exact window with a record | 113 288 and 389 029 | – | as B | Prop W | nothing cheap: padded detours are 4·10⁷ for 326 (final), 242 953 / 292 969 for 385/467 |
| E | 8 (219, 319, 332) | degenerate (G = 1); `Rex` chain | 0 | planner templates (≤ 90 flows in B2) | single-call bounds 0.70-12.7 M, so planner chains | Thm 6.11 | adopt W_t := W_{t'}·ρ |

**Outlook per pair [E].**
- **N = 5:**
  - 313 ≈ 85 %: three independent routes (the 484 window with a rank test, the 283 842 record, the in-budget padded chart).
  - 326 ≈ 60 %: three over-budget steps (the B step with out-detour 300 838, the D2 record 113 288, and a hidden final step of 147 091).
  - 410 ≈ 60 %: one cheap route only (the 3 583 window with a rank test, G = v⁰); the padded fallback needs a 4-upright exit chain.
  - All three together: 45-60 %, which is lower than the census analysis's 60-70 %.
- **N = 6:**
  - the 22 clean pairs: 15-18 certified;
  - the 7 group-C pairs: 3-5;
  - the 3 group-D2 pairs: 1-2;
  - total 19-25 of 32.

## 3. Code changes (no frozen file on the primary route)

| file | function | change |
|---|---|---|
| src/dcoset.rs | `chart(s, tp, rho, t) -> Chart { g: Grass, z3: Option<Vec<Rung>>, z4: Vec<Rung>, p3: Vec<Rung>, mirrored, needs_pad }` | new: `grass`, then Z by `bfs_reduced` on z (3 uprights) or z·0 (4 uprights), then projected prefix p3; helpers `mirror_rung` and `mirror_table` for K = {1,2} |
| src/dcoset.rs | `local_type(sw, w) -> Option<Tid>`, `window_gloc(lam, at, a, rho, b) -> Option<(Tid, Tid, Lp)>` | new: checks W1 and W2 on `Schur::new(ν)` |
| src/leaves.rs | `window_cert(s, lam, wfrom, rho, wto, tp, t, g) -> WinCert` | new: W1-W4; W4 uses the `FlowIndex`/`CosetSpace` code of `step_family`'s `cols`, without an Arena |
| src/leaves.rs | `Leaves::joint_pair(..) -> Vec<JointCand>` | new: prefix-trie search over `reduced_words(t', 6)` × `reduced_words(t, 7)` plus seeds P·A; keeps the top k |
| src/leaves.rs | `Leaves::plan_b5(..) -> Plan { words, steps: Vec<Vec<StepPlan>> }` | new, run after `plan_words`: chooses adopt or detour per step; dry `Planner::plan` for R_in, R_out; closed-form shape test |
| src/leaves.rs | `StepKind` | add `Rex`, `Detour`, `ClosedForm`, `Chart4` (Chart4 only in B5d) |
| src/leaves.rs | `step_family` | dispatch on `StepPlan`; the Window/Whole code stays the default path; `no_whole` flag for the dry rank test |
| src/leaves.rs | `build_leaves`, `pairs_leaves` | consume the plan; `pairs_leaves` adds the B4 check (form⁻¹·π)∘ι = id blockwise at every coset of L_q |
| src/leaves.rs | `Leaves::b5` flag, `LEAVES_STAMP` | flag off by default, so A closure and B4 are byte-identical; suffix `b5=grass(v1,…)` only when it is on |
| src/leafgen.rs | `step_closed_form(labels, a, b, dir, d) -> Option<ClosedForm>` | new: Pull/Push via `ClosedForm::new`, then `build`/`build_in` |
| src/bin/gp0.rs | `--all-steps --b5 --wcert --rank-dry` | census of every step (not just `worst_step`), with route, W flags and rank verdict per step in the jsonl |
| src/bin/m4c.rs | `measure` | records for D2 (113 288, 389 029) and for any 1-4·10⁵ fallback |
| tests/b5grass.rs | new | validation V0 (§5) |
| src/smap.rs, src/smapio.rs | `MapNode::Pad`, a `P` record | **B5d only**, as in D1 §2.4 (about 120 lines, reviewed on its own) |

Untouched: homprov.rs (`certify_pairs`, `win_split`, `hom_unknowns` are used as they are), rexplan.rs, sbim, linsys, tmpl, sindec.

## 4. Correctness and trust (certification unchanged)

- **Final acceptance.**
  - `certify_pairs`: count, form rank = mult at rep(q) at two points, genuineness.
  - Then the π∘ι check against id_{L_q}.
  - Intermediate levels keep `build_leaves`' diagnostic (square, invertible at rep(t)).
- **Map nodes.** Only Comp, Kron, Local, Id and Elem-Assoc, as in B2/B4. No new `MapNode` on the primary route.
  - Exact windows are `hom_space` templates.
  - Closed forms are `Added` templates, re-certified when imported.
  - Chain Assocs are numerically validated. The theorem-level Assoc stays available to the planner, but B2 never needed it (largest Assoc window: 60 flows).
- **What Prop W is used for.** It only predicts success: planning and budget. It never admits a map. If the prediction is wrong, the result is `Exhausted` or a rank shortfall, never a wrong basis.
- **`Pad` (B5d only).** It adds one trusted node, an identity: C_m4c M7 rates D1's trailing-zero `Pad` as trivial new trust (`CosetSpace::rank` ignores trailing zeros).

## 5. Validation plan

- **V0: small λ, against exact ground truth** ((3,3), (2,1,2), (2,2,1), (2,2,2); minutes).
  - (a) For every Lemma A step of every closure plan and every G1(a) pair plan, and for every pair of reduced words (≤ 4 rungs) of t', t that share a prefix: compare the W verdict with the actual rank of the window family (exact window + greedy at one point).
    - Requirement: 0 cases of "W passes, rank short".
    - Also report the share of "W fails, rank full", which measures how conservative W is.
  - (b) Every Pull/Push candidate is genuine and lies in the exact Hom^d window (≤ 2·10⁴ unknowns), as in B1.
  - (c) `pairs_leaves` with `b5` on, over all (2,1,2) and (2,2,2) G1(a) targets: the certified set must equal the Exact backend's (264 targets, HANDOFF_M4b).
- **V1: reproduce the (3,3,3) results.**
  - (a) A closure, 12/12, with `b5` on: certification identical, maps may differ.
  - (b) The 7 B4 rex-class targets, both primes (no Lemma A step, so no change).
  - (c) The non-rex target of pair 313: its level-5 step (the B2 fallback to 1.96 M) must fail W in B5a and be rerouted.
- **V2: the list.** Order:
  1. N = 5 (313, 326, 410), both primes.
  2. N = 6, groups A/B/D1/E.
  3. Group C.
  4. Group D2.

  For each target:
  - `certify_pairs` and the π∘ι check;
  - replay from layer A in a fresh arena (Added templates re-certified);
  - G-L5 counters (DAG nodes, resident index bytes) for M6a (C_m4c M4).

  The downstream checks (sindec wiring, d∘d, K₀, χ) follow B2 §6.3 in the M6a integration, not in B5.

## 6. Sessions, gates, stopping point

| session | content | gate |
|---|---|---|
| B5a (1) | All-steps census, charts, joint words, W flags, closed-form shapes, dry planner runs for detours. Dry rank test on every chosen window: at most about 60 exact calls ≤ 78 371 unknowns, under 1 h [E]. | **G-B5a:** if at least 70 of the 101 targets have *every* step on an in-budget route, go. With 40-69, do N = 5 and the clean pairs only. Below 40, no-go: document and stop. |
| B5b (1-1.5) | leaves.rs step kinds, `plan_b5`, `step_closed_form`, V0, V1. | V0(a) has no counterexample; V1 unchanged. |
| B5c (1-1.5) | V2 runs, both primes. `m4c measure` for D2 and for any fallback record. | **Stopping point:** report certified / refused per target. |
| B5d (1-2, optional) | `Pad`, `Chart4`, Planner on λ·0. | Trigger: 313 or 410 is the last N = 5 blocker after B5c. |

Total: 3-4 sessions, plus 1-2 if B5d is needed. Any pair still failing after B5c stays a documented gap and is reported as `Refused`.

## 7. Risks

| # | risk | likelihood [E] | detection | consequence |
|---|---|---|---|---|
| 1 | Hidden over-budget steps: the census lists only the worst step; 12 targets are known to have a second step > 10⁵ | medium | B5a all-steps census | the same pipeline per step; targets that still fail go to the fallbacks |
| 2 | W3/W4 fail for projected prefixes (P is not a Thm 2.6 chart: L' = (c₀, c₁, c₂+o) ⊄ Js) | medium | B5a W flags | no a-priori completeness; the dry rank test still decides cheaply |
| 3 | Rank test fails on group-C windows (G_loc ≠ G) | medium | B5a | records of 283 842 / 318 975, or the padded chart (313) |
| 4 | Detour chains are not found within budget (32 of 42 exceed a single 10⁵ call) | low-medium | dry `Planner::plan` in B5a | adopt the joint words instead, or the next candidate |
| 5 | Pull/Push is not genuine on the D1 windows | low-medium | `certify_added` | records of 194 000-295 452 unknowns (A2 measured 434 787 in 179-272 s, 2.3-3.2 GB) |
| 6 | D2 records are censored (92 of the 136 M5b records in 10⁵-10⁶ were censored) | medium | `m4c measure` | 326 has no cheap fallback (padded detours 4·10⁷), so N = 5 stays open |
| 7 | Downstream M6a cost of deeper DAGs (2 chains per moved step) | low-medium | G-L5 counters | compaction; Leaves stays last in the chain |
| 8 | B5d: exact windows with a zero label on λ·0 and the Swap(U, ρ) templates are untested (B1's (3,3,2,1) far commutation took 820 s) | medium | B5d V0 on (2,2,2)·0 | Pad route limited to 313 |
