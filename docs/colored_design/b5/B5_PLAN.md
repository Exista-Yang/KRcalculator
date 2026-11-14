# B5 plan: Lemma A lower steps for the remaining (3,3,3) pair targets, N = 5 and N = 6

Date 2026-11-12. Synthesis of four read-only documents in this folder:
- `B5_census_analysis.md`: the 42 worst steps, 3-upright joint words, groups A-E;
- `D_grass.md`: the rung-local Grassmannian route, Prop W, Detour/Rex steps;
- `D_sinister.md`: localized sinister leaves, PullTw/PushTw;
- `C_B5.md`: the adversarial critique, findings F1-F10.

It is checked against branch `m4c-rexplan` 4e886f6 (worktree `revrex`): `b2/gp0_n5b.jsonl`, `gp0_n6.jsonl`, B2 §4-6,
M4C_PLAN §3.1/§3.3/§8, and `docs/colored_design/gl6/M5c_coverage.md`. No cargo was run.

**Tags:**
- [V] computed (Python replica or jsonl);
- [C] read in the code;
- [P] paper (SLL = EKLP 2401.03053);
- [K] inference that a gate checks;
- [E] estimate.

## 0. Decision

**Hybrid. GO for B5a only (1 cloud session). Commit to N = 5 behind a hard gate. N = 6 is a conditional extension, not
the goal.**

1. **Primary route: D_grass's.**
   - Exact windows on 3-upright joint junction words.
   - B1 `Pull`/`Push` closed forms where the shape fits.
   - `Detour`/`Rex` steps through B2 planner chains.
   - Prop W used as a planning certificate.
   - No new `MapNode` and no frozen file.
2. **C_B5's B5a corrections are mandatory:**
   - an all-steps census taken from the executor's own plan;
   - per-predecessor final words with their own exits;
   - a label guard in the planner;
   - records counted per (direction, degree).
3. **D_sinister's `PullTw`/`PushTw`** is an optional add-on (B5d), used only on named steps where it replaces records
   that would otherwise be censored.
4. **`Pad`/`Chart4`** stays deferred (B5d) and is used only for pair 383 or for a last N = 5 blocker.
5. **The literal sinister-chart route is rejected.** It would need 4-5-upright charts, `Pad` twice and planners that
   have never run. D_sinister §0.1 and C_m4c M2 both say no-go.
6. **Why not STOP now.**
   - B5a is cheap and decides the question.
   - Pair 313 probably closes with records alone (§4.1).
   - The N = 5 odds (≈ 30 %) are too uncertain to call either way without the all-steps census.
7. **Why not promise more.**
   - Partial completion of a level yields **no knot values** (§1, fact K).
   - N = 6 needs all 25 universal pairs, including the non-local 383.

## 1. What is established, and what moved the odds

| # | fact | source | consequence |
|---|---|---|---|
| A | No worst step needs a 9-root K block; rung-local K windows are 45-57 666 unknowns with G_loc = G | census §0, [V] | a small factored step exists mathematically for every target |
| B | 3-upright joint words (common prefix = Grassmannian prefix with pad rungs dropped) give windows ≤ 10⁵ for 37/42 and ≤ 4·10⁵ for 42/42 worst steps; labels ≤ 6; junction words ≤ 38 880 flows | census §0.3, C_B5 §2 [V] | no `Pad` and no GrTable(·, 9) on the primary route |
| C | Locality (W3, W4 with the P-part **and** the Q-part) holds for 39/40 chosen non-Id windows; it fails only for 383 level 7 (0/6 columns local) | C_B5 F4 [V] | Prop W (D_grass §1.3) is complete a priori for the 30 windows with G_loc = G; 383 needs `Pad` or stays a gap |
| D | Corrected template cost: every candidate template ≤ 204 s (B2 fit; about 3× more for n = 720) | C_B5 F5 [V] | closed forms are affordable; 458 is not over the cap |
| E | The census `worst_step` is stale after the alt-final swap (gp0.rs:174-177 updates `max_step`, not `worst`) | C_B5 F1a [C] | the 326, 410, 383, 283, 293, 385 and 82 strings describe pre-swap steps |
| F | Several final steps share W_q, so hidden final predecessors exist in 9/12 checkable groups, including **every N = 5 target** | C_B5 F1b, F2 [V] | final words and exits must be chosen per predecessor (`pairs_leaves` today composes a single exit, leaves.rs:1179) |
| G | 74/101 targets have unenumerated non-final steps (up to 29 non-Id steps per target) | C_B5 F1c | coverage counts in D_grass (74/101) and D_sinister (67 closed) are upper bounds |
| H | B2's "window not complete" evidence (313: 1.96 M fallback) was a **budget refusal** of a hidden final step's 122 792 window, not a rank shortfall | C_B5 F3, `b2/b4_313_la_f31.log` l.3 | no (3,3,3) window family has been seen to lack rank; D_grass V1(c) is replaced (§3) |
| I | Records are keyed by (labels, fs, ft, unknowns, deg), so each (direction, degree) above 10⁵ needs its own record | C_B5 F7, homprov.rs:1104 [C] | about 14-20 record runs at 1.6-3.9·10⁵ unknowns if everything goes by records |
| J | `PlanOpts` (rexplan.rs:308) has no label guard; a Collapse can create a label-9 level, which needs GrTable(r, 9) (> 90 min) | C_B5 F6 [C] | add `max_label ≤ 8` |
| **K** | **Knot coverage is all-or-nothing.** At N = 5 every level-5 failure blocks every 3-strand knot. At N = 6, 57 of 68 failures do. 4_1 at gl₆ needs exactly the 57 universal ones. Of the 32 non-feasible N = 6 pairs, **25 are universal**: 82, 127, 212, 273, 278, 281, 283, 285, 292, 293, 313, 319, 326, 332, 335, 336, 376, 381, **383**, 385, 410, 411, 418, 458, 467. The other 7 (24, 45, 78, 137, 145, 219, 227) are needed by 104-161 of 185 knots | M5c_coverage.md (LB), intersected with gp0_n6.jsonl [V] | a certified pair has library value only. Knot values at gl₅ need 313, 326 and 410 all closed; at gl₆ they need all 25 universal pairs |

**Why the odds fell.**
- D_grass (45-60 %) and D_sinister (70-80 %) priced only the stale worst steps.
- E-G add, for N = 5:
  - a 295 452 hidden final step for 410;
  - a 122 792 hidden final step for 313;
  - 6 unenumerated non-Id steps for 326 q₁;
  - a post-swap final step of 15.7 M for 326 q₂, which is the same t' → q as the analysed D2 step but a different word.

## 2. Design (what B5 builds)

### 2.1 Per step, attempt in this order (`Leaves::step_family`, new `StepPlan` from `plan_b5`)

1. **`Id` / `Rex`.** Use when G = 1 and W_{t'}·ρ is reduced for t. Realize the planner chain W_{t'}·ρ → W_t with
   `leafgen::realize_chain`. The 1×1 top block is 1 (SLL Thm 6.11, Ex 7.18). This covers group E: 219, 319, 332, and
   the rex-type hidden predecessors of 283, 383, 293.
2. **Closed form.** Use when G = v^d with coefficient 1, a side is flag-type and the degree matches.
   - Build B1 `CfKind::Pull` (Out) or `Push` (In) on the joint window with `leafgen::build`/`build_in`.
   - `TmplLib::certify_added` decides.
   - Candidates: group A, D1, 8 B steps, and 410's hidden 295 452 step (S2 shape, 100×120 flows).
3. **Exact window.** ≤ 10⁵ unknowns, or ≤ 4·10⁵ with an uncensored record for that (direction, degree), via
   `Leaves::admit` and `try_hom_basis_cap`, then the existing greedy selection.
   - Prop W (W1-W4, Q-part included) marks the window complete a priori.
   - Otherwise the B5a dry rank test must have passed (group C: 410, 283, 293, 313-alt, 137/145).
4. **Connection.**
   - **Adopt** the joint words W'_{t'}, W'_t as junction words when every step touching t' or t stays in budget.
   - **Otherwise detour:** S = R_out ∘ Local(id_P ⊗ s) ∘ Kron(R_in, [ρ]), with R_in, R_out degree-0 planner chains
     between reduced words of one type. Their top scalar is 1 (SLL Prop 4.19, Thm 6.11), so the top-to-top block of
     S equals that of s (D_grass §1.4).
5. **Final steps.** Use one final word W'_{q,i} per predecessor i, each with its own exit chain R_i : W'_{q,i} → L_q
   (B2 planner on 3 uprights, as for the 118 exit chains). This is a `Detour` whose R_out is an exit. It is mandatory
   for 313, 410, 283, 383 and 293 (C_B5 F2).
6. **No route.** The target is reported `Refused`, never accepted.

### 2.2 Completeness and certification

- **Prop W** is used only for planning and budget (D_grass §1.3; local Prop 8.32 + Thm 8.31 proof + Lemma 8.7).
  It never admits a map. A wrong prediction gives `Exhausted` or a rank shortfall, not a wrong basis.
- **Acceptance is unchanged:**
  - `homprov::certify_pairs`: count, form rank = mult at rep(q) at 2 points, genuineness;
  - then π∘ι = id_{L_q} blockwise at every coset of L_q, at 2 points (B4 procedure);
  - both primes;
  - replay from layer A, with Added templates re-certified on import.
- **Map nodes** are only Comp, Kron, Local, Id and Elem-Assoc, as in B2/B4.
  - The theorem-level Assoc certificate (approved) stays available to the planner, but no B2 chain needed it.
  - If a B5 chain does need it, the approved cross-checks and labelling apply.

### 2.3 Optional add-ons (B5d, only on trigger)

- **`CfKind::PullTw`/`PushTw`** (D_sinister §1.2): a Pull/Push twisted by one incidence Euler factor ε(M_j, R), with
  sprinkles as `Dot`/`Lin` nodes. Templates are `Added` and re-certified, so there is no new trust.
  - Named uses: 313's G_loc = G window (283 842; saves 4 records) and 326 q₂'s D2 window (113 288, ε = ε(M₁, M₂),
    hand-checked; saves 2 records).
  - It does **not** help 385/467 or 137/145, whose cover degree is too high.
- **`MapNode::Pad` + `Chart4` + Planner on λ·0** (D_grass §1.6, approved "only if strictly needed"; about 120 lines in
  smap.rs/smapio.rs, reviewed separately). It is the only known route for 383, whose local window is 690, and a
  fallback for 410/313.

## 3. Milestones (cloud sessions of about 4-5 h; never end a turn while jobs run)

The base is `m4c-rexplan` merged with the result of the running B4w session (pairs wiring, levels 5/6 rebuild,
leafgen OOM gate). B5 starts after B4w is merged into dev.

| session | content (file:function) | gate | targets expected to certify |
|---|---|---|---|
| **B5a** (1) | see B5a below | **G-B5a**, §3.1 | 313 non-rex target, N = 5 and N = 6 (records only) |
| **B5b** (1-1.5) | see B5b below | V0(a) has no "W passes, rank short" case; V0(c) equals the Exact backend; V1 is byte-identical with `b5` off and certification-identical with it on | none new (validation) |
| **B5c** (1) | see B5c below | N = 5 library 177/177, both primes, then 4_1 at gl₅ | **N = 5: 326 (both q), 410**, plus 313 if B5a has not already closed it. N = 6 by-product: the pairs B5a found fully routed with no record. Expected 8-15, chiefly A/B/E: 281, 381, 212, 376, 411, 24, 127, 273, 278, 336, 418, 219, 319, 332 |
| **B5d** (1-2, conditional) | PullTw/PushTw (leafgen.rs, new `wsl.rs`: `shape`, `find_family`, `sprinkle_dag`) and/or `MapNode::Pad` + `StepKind::Chart4` + gp0 `pl4` chains | trigger: the N = 6 gate passes, **or** exactly one N = 5 step is blocked and one of these fixes it | 383 (Pad); 313/326 without records |
| **B5e** (2-3, conditional) | N = 6 runs of everything routed; record campaign (about 14-20 `m4c measure` runs at 1.6-3.9·10⁵: 137/145, 385/467, D1 fallbacks); group-C rank-tested windows; level-6 rebuild; 4_1 at gl₆ | the 25 universal pairs certified, both primes | N = 6 universal set, then the 7 non-universal pairs |

**B5a content:**
1. `gp0 --all-steps --b5` takes the steps from `Leaves::build_leaves` → `plan_words` after the alt-final swap. Fix the
   stale `worst` string (gp0.rs:174-177).
2. Per step, record:
   - the `leaves::joint_pair` candidates;
   - `dcoset::window_gloc`;
   - `leaves::window_cert` (W1-W4 with the Q-part);
   - the Pull/Push shape;
   - the records needed per (direction, degree).
3. Per final predecessor: W'_{q,i} and a dry exit `Planner::plan`.
4. Dry `Planner::plan` for every detour, with `PlanOpts.max_label = 8` (rexplan.rs:308; the same guard in
   `joint_pair`).
5. Dry rank test (`step_family` with Whole disabled, one point) for every chosen window ≤ 10⁵: about 60-100 calls,
   1-2 h.
6. Dry λ·0 planner (`pl4`) for 383's chart detours.
7. **Quick win:** `m4c measure` 313's two final windows (137 964 and the hidden 122 792), every (direction, degree)
   above 10⁵. Then `gp0 --b4 31|32b --lemma-a` on 313's non-rex target.

**B5b content:**
- leaves.rs:
  - `StepKind::{Rex, Detour, ClosedForm}`;
  - `Leaves::plan_b5` (adopt or detour);
  - `joint_pair` called from `plan_words` (leaves.rs:798);
  - per-predecessor exits in `pairs_leaves` (leaves.rs:1179);
  - the `b5` flag, off by default, with the stamp suffix only when it is on.
- leafgen.rs: `step_closed_form`.
- Validation:
  - **V0** on (3,3), (2,1,2), (2,2,1), (2,2,2), all three parts against exact ground truth:
    - (a) W verdict against actual rank;
    - (b) Pull/Push lie in the exact window;
    - (c) `pairs_leaves` with `b5` on, over the 264 G1(a) targets.
  - **V1:** A closure 12/12 and the 7 B4 rex-class targets.
  - **V1(c)** is replaced: 313's hidden final step must be planned per predecessor and must not fall back to Whole.

**B5c content:**
- N = 5 runs, both primes: 313, then 326 q₁ and q₂, then 410.
  - Each target gets `certify_pairs` + π∘ι, and replay.
  - The fallback records are measured in-session.
- Level-5 `libbuild` rebuild.
- 4_1, then the 185 knots at gl₅ via the M6a boxrun driver.
- G-L5 counters.

### 3.1 G-B5a, the go/no-go gate

"Routed" means every step of the target, taken from the executor's plan, is one of:
- (i) ≤ 10⁵ with W-cert or a passed dry rank test;
- (ii) a closed-form shape with a predicted template ≤ 300 s;
- (iii) ≤ 4·10⁵ with an uncensored record measured in B5a;
- (iv) Id/Rex/detour/exit, with planner chains found within the B2 budget (templates ≤ 1 200 flows, Exact ≤ 10⁵,
  label ≤ 8).

| outcome | action |
|---|---|
| all of 313, 326 q₁, 326 q₂, 410 routed (a closed-form shape still to be certified counts as routed) | **GO N = 5**: run B5b and B5c |
| exactly one N = 5 step unrouted, with a named B5d fix | GO, with B5d moved before B5c |
| two or more N = 5 steps unrouted, or one with no named fix | **NO-GO**: stop B5 (§6) |
| N = 6: at least 23 of the 25 universal pairs routed, and 383's λ·0 detours found by the dry `pl4` run | N = 6 extension (B5d + B5e) becomes eligible after B5c |
| otherwise | N = 6 3-strand stays a documented gap. Routed pairs are certified as a B5c by-product only, with no extra sessions |

## 4. Target expectations

### 4.1 N = 5 (4 targets in 3 pairs; all must close for any gl₅ knot value)

| target (q) | non-Id steps | known binding steps [V] | route | P [E] |
|---|---|---|---|---|
| 313 non-rex (1,2,0/2,0,1/0,1,2) | 3 (max non-final 53) | final 1,2,0/2,0,1/0,2,1 → q, G = v⁰+v², 137 964 (one record measured in B2, `b2/costs_b2.txt`); **hidden** final 1,2,0/2,1,0/0,1,2 → q, G = v², window [5,1,3] 120×40 flows, 122 792 at d = 2, refused for lack of a record | B5a records + rerun. Otherwise per-predecessor words: 18 461 local S2 window; group C 484 window with rank test; PullTw | **≈ 85 %** |
| 326 q₂ (1,2,0/1,0,2/1,1,1) | 3 (max non-final 15 375) | final 1,2,0/1,1,1/1,1,1 → q, F1.1, G = v¹+v³: 15.7 M post-swap; joint window 113 288 [5,3,1], G_loc = G, detours 0 / 6 425 | 2 records (one degree, both directions), or PullTw ε(M₁, M₂) | **≈ 70 %** |
| 326 q₁ (2,1,0/1,0,2/0,2,1) | **8** | level 5 3,0,0/2,0,1/0,2,1 → 2,1,0/1,1,1/0,2,1, G = v²: 956 320 → joint 78 371, Pull/Push shape; out-detour single bound 300 838 at an **intermediate** type (never planned there); alt final 147 091 (unanalysed); 6 steps unenumerated | adopt or detour + Pull/Push; final per predecessor | **≈ 55 %** |
| 410 (1,1,1/2,0,1/0,2,1) | 3 (max non-final 24) | two finals that no shared W_q (≤ 6 rungs) keeps under 1.9 M: 1,2,0/2,0,1/0,3,0 → q, G = v⁰, 3 583 window with G_loc = v⁰+v² (rank test; local columns ⊃ Φ, so triangularity suggests it passes [K]); **hidden** 1,2,0/2,1,0/0,2,1 → q, G = v¹, 295 452 [5,3,1] 100×120, S2 shape | per-predecessor finals + rank test + Pull/Push (template ≈ 150-200 s); fallbacks: 2 records at ≈ 3·10⁵ (censoring risk), Pad local 3 583 | **≈ 55 %** |

### 4.2 N = 6 (97 non-feasible targets in 32 pairs; 66 targets in the 25 universal pairs)

| class | pairs | expectation |
|---|---|---|
| clean worst step (A/B/D1/E) | 24, 45, 78, 82, 127, 212, 219, 227, 273, 278, 281, 285, 292, 319, 332, 335, 336, 376, 381, 411, 418, 458 | about 12-17 certify if B5 works. Hidden steps: 82 (1.04 M, 14.6 M), 285 (10.4 M, 151 774), 292 (356 983), 411 (135 902); 292/335 need 7-rung q words |
| group C | 137, 145, 283, 293, 313, 383, 410 | 313 and 410 as at N = 5. 283/293: rank tests on 58-unknown windows, plus a rex-type hidden predecessor (per-step exit). 137/145: 4 records or a rank test. **383 is non-local: Pad or a gap** |
| group D2 | 326, 385, 467 | 385/467: 4 records each at 2.3-3.9·10⁵, plus a hidden alt final of 40.3 M and 1.17 M, and non-final 1.69 M. PullTw does not apply |

## 5. Probabilities (honest; from today, before B5a)

| outcome | P [E] | main driver |
|---|---|---|
| 313 certified (N = 5 and N = 6) | ≈ 85 % | records alone may suffice |
| **N = 5 complete** (313, 326 ×2, 410), giving gl₅ for all 185 3-strand knots | **≈ 30 % (20-40 %)** | 326 q₁'s 6 unenumerated steps and intermediate-type detours; 410's two finals |
| G-B5a says GO for N = 5 | ≈ 40 % | if it does, P(N = 5 complete \| GO) ≈ 70 % |
| N = 6 pairs certified by the end of B5c | 10-18 of 32 (library value only) | hidden steps, records |
| **N = 6 knot values** (all 25 universal pairs, so 4_1 at gl₆) | **≈ 7 % (3-15 %)** | 383 non-local (Pad + never-run λ·0 planner), 385/467 record chains with hidden 40 M finals, 326, 410 |
| all 32 N = 6 pairs | ≤ 5 % | as above, plus the non-universal C and D2 pairs |

**Effort.**
- N = 5 path: B5a + B5b + B5c = 3-3.5 sessions.
- N = 6 extension: a further 3-5 sessions (B5d + B5e).
- So 6-8.5 in total, consistent with C_B5's 5-7 to the stopping point. M4C_PLAN R9 (overrun) remains high.

## 6. Recommended stopping point

1. **Default: stop after B5c**, whatever happens there, unless the N = 6 eligibility row of G-B5a holds. The
   expected cost of N = 6 is ≥ 3 more sessions for a ≈ 7 % chance of a gl₆ 3-strand value.
2. **If G-B5a says NO-GO: stop after B5a** (1 session spent) and write up what exists:
   - 2-strand 3-colored gl₆: T(2,3/5/7/9) at N = 3-8 equal Wang; T(2,2k) links are new;
   - 3-strand: 185/185 knots at N = 3, 4. These equal uncolored gl₄ of the reversed braid by duality, so they
     validate the singular engine but are not new invariants;
   - the certified (3,3,3) library: closure 12/12, 7 rex-class pairs, 313 if closed, the level 5/6 status;
   - the B5a census as a precise map of what blocks gl₅/gl₆, plus the method (localized singular light leaves, Lemma A,
     planner).

   Be explicit that "partial N = 5" gives **no** gl₅ knot values (fact K). The publishable unit is the library and the
   method, not partial tables.
3. **If B5c completes N = 5:** add 3-strand 3-colored gl₅ for all 185 knots. These are the first genuinely new 3-strand
   values. By ∧^a ↔ ∧^{N−a} duality they equal 2-colored gl₅ of the reversed knot, which gives an independent
   cross-check if the (2,2,2) route is cheap [K]. N = 6 3-strand is then a documented gap with its census.
4. **Continue to N = 6 only if the gate holds**, capped at +4 sessions. Stop at the cap or at the first universal pair
   that has no route.

## 7. Main risks

| # | risk | P [E] | detection | consequence |
|---|---|---|---|---|
| 1 | unenumerated steps (326 q₁: 6; 74/101 targets) are over budget with no route | medium-high | B5a all-steps census | N = 5 NO-GO, or that target is Refused |
| 2 | detours at intermediate types are not found (planner only ever ran at final types) | medium | B5a dry `Planner::plan` | adopt the joint words; next candidate |
| 3 | records censored (92/136 M5b records in 10⁵-10⁶ were censored) | medium | `m4c measure` in B5a/B5c | PullTw (313, 326 q₂) or Pad; otherwise a gap |
| 4 | Pull/Push not genuine on a window | low-medium | `certify_added` | exact window with a record |
| 5 | group-C rank tests fail (410, 283/293) | low-medium (F4) | B5a dry rank test | Pad local chart (B5d) or a gap |
| 6 | M6a downstream cost of deeper DAGs (2 chains per moved step) | low-medium | G-L5 counters | compaction; Leaves last |

## 8. Corrections recorded against the inputs

- **Census §5.3**, "12 targets known to hide a second step": this is a lower bound (C_B5 F1c).
- **D_grass:**
  - §1.3, "137 964 → 1.96 M": the 1.96 M was the hidden final step's whole-word fallback, after its 122 792 window was
    refused by budget.
  - V1(c) is replaced (§3, B5b).
  - W4 must include the Q-suffix.
- **D_sinister:**
  - template costs per direction are wrong; all are ≤ 204 s by the B2 fit (C_B5 F5);
  - "313, 326 closed" ignored their hidden final steps;
  - "458 over cap" is withdrawn.
- **Both designs:**
  - "a final W'_q is always adopted / one exit" is wrong when q has several predecessors (C_B5 F2);
  - their N = 5 odds (45-60 % and 70-80 %) are superseded by §5.
