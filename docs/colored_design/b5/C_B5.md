# C_B5: adversarial critique of D_grass.md and D_sinister.md (B5 design)

Date 2026-11-12. This review is read-only, against `m4c-rexplan` 4e886f6 (worktree `revrex`). No cargo was run.

**Inputs:**
- D_grass.md, D_sinister.md, B5_census_analysis.md (this folder);
- B1/B2 reports, M4C_PLAN §3.1/§3.3, L1;
- `src/{leaves,leafgen,rexplan,homprov,tmpl,sbim,local}.rs` and `src/bin/gp0.rs`;
- `b2/gp0_n5b.jsonl`, `gp0_n6.jsonl`, `b2/b4_313_la_f31.log`.

**New checks.** The scripts are in `b5/tools_c/`. They reuse the census replica `tools/kr.py` and `hu.py`, and every run took under 20 s.

| script | what it checks |
|---|---|
| `mtrue.py` | template monomials and cost |
| `w4.py`, `w4c.py` | locality of the global columns |
| `w3.py`, `gloc.py` | local types and G_loc |
| `finals.py`, `preds.py`, `shared.py`, `hid410.py`, `hidden313.py` | hidden final steps |

**Tags:**
- [V] computed here;
- [C] read in the code;
- [K] inference.

## 0. Verdict

**Conditional go for B5a only (one session, combinatorics plus dry rank tests). Both designs' coverage and N = 5 odds are too optimistic.**

- **The math holds up.** Locality of the chosen windows, the label caps and the template costs are better than either design could show.
- **The work list is incomplete and partly stale (F1, F2):**
  - the census `worst_step` string is stale after the alt-final swap;
  - several final steps share one W_q;
  - 74 of 101 targets have unbounded hidden steps.
- **In every N = 5 target where X could be inferred, there is a second over-budget final step that neither design analysed.**
- **Revised outlook [E]:**
  - N = 5 complete: **20-35 %**. D_grass said 45-60 %, D_sinister 70-80 %.
  - N = 6: about 10-18 of 32 pairs.
  - About 5-7 sessions to the stopping point, not 3-4.
- **Recommended shape.** Use D_grass's route as primary: no new kinds, exact windows plus B1 Pull/Push, and Detour/Rex steps. Add D_sinister's PullTw only for the named steps where it saves record runs. Make B5a's gate count **all** steps, from the executor's own plan.

## 1. Findings, most severe first

### F1 (critical): the census basis misses binding steps. Coverage counts are worst-step upper bounds.

**(a) Stale worst step [C].**
- In `gp0.rs` (lines 174-177), the alt-final branch sets `max_step = max_nonfinal.max(best.0)` but never updates `worst`.
- So in 12 (N, pair, q) groups with an `alt_final`, the `worst_step` string describes the step *before* the W_q swap:
  - 326 ×4, 410 ×2, 383, 283, 293, 385 ×2, 82.
- The census analysis, and both designs, analysed those strings. The binding post-swap step is often a different final step. [V]:
  - for 410, with the alt W_q `E0.2,E1.3,F0.1,F1.3,F0.1`, the census step costs 1 011 648 and the gp0 maximum of 1 935 144 is a different predecessor;
  - for 383 and 283, the gp0 maximum of 27 245 686 is a hidden predecessor.

**(b) Several final steps share W_q [V] (`finals.py`).**
- I inferred X = W_t'·ρ from the worst-step string. It matched `x_flows` for 12 groups.
- In 9 of those 12, q has a second predecessor that nobody analysed:

  | pair | hidden predecessor → q | G | best window alone | with a shared W_q |
  |---|---|---|---|---|
  | 313 | 1,2,0/2,1,0/0,1,2 | v² | 18 461 (local, S2) | min over 12 q words: 122 792 |
  | 410 | 1,2,0/2,1,0/0,2,1 | v¹ | 295 452 on [5,3,1], 100×120, S2 shape, local | no W_q of ≤ 6 rungs gets both finals under 1.9 M |
  | 283, 383, 293 | 0,3,0/2,0,1/1,0,2 and 0,3,0/1,0,2/2,0,1 | v⁰ | rex-type: all 19 reduced words W satisfy "W·ρ reduced for q", so Id with its own W'_q | – |
  | 335, 292 | – | – | size unknown (q needs 7-rung words) | – |

**(c) Hidden non-final steps.**
- 74 of 101 targets (34 of 50 (N, pair, q) groups) have their worst step, or `max_nonfinal`, above 10⁵. The other steps of those targets are unknown: up to 29 non-Id steps per target (227).
- The census analysis says "12 targets known". That is a lower bound, not an estimate.

**Consequences:**
- D_grass's "74/101" and D_sinister's "67 closed" are worst-step-only upper bounds.
- Per-pair odds, especially for N = 5, must be recomputed after an all-steps census.
- The census must be emitted from the **same plan the executor uses** (`Leaves::build_leaves` → `plan_words`, after the alt-final swap), listing every step.

### F2 (high): the policy for final words is wrong for several predecessors

- D_grass §1.4 says "a final W'_q is always adopted". D_sinister §3.6 says "exit chain W'_q → L_q (as `alt_final`)". Both assume one free W'_q per step.
- But every final step into q shares W_q, and `pairs_leaves` composes **one** exit chain with one family word (leaves.rs:1188-1200) [C].
- 313, 410, 283, 383 and 293 need per-predecessor final words: R_out,i : W'_q,i → L_q, or into a common W_q. That is D_grass's `Detour` with an exit per step.
- **Code changes needed:**
  - per-step exits in `pairs_leaves`;
  - a dry run of `Planner::plan` for each W'_q,i in B5a.

  Neither design lists them.

### F3 (high): the B2 "open problem 2" evidence is misread. D_grass §1.3 and V1(c) are wrong.

- **The failing call.** `b4_313_la_f31.log` line 3 shows the call that failed: `[E0.2,E1.1,F0.2,F1.1] -> [E0.2,E1.3,F1.3,F0.2] d=2: 1963540`.
  - This is the **hidden** final step 1,2,0/2,1,0/0,1,2 → q.
  - It is not the census step, whose whole-word calls are 2 515 760 and 4 561 122 [V].
- **What happened.** Its window is [5,1,3], 120×40 flows, 122 792 unknowns at d = 2 [V].
  - `Leaves::admit` refused it: over budget, no record.
  - `step_family` swallowed the window error, the greedy selection then ran short, and it fell back to Whole, which was refused.
- **So** this was a budget refusal, not a rank shortfall. There is **no observed (3,3,3) case of a local window family lacking rank**.
- **What to fix:**
  - V1(c), "must fail W in B5a", tests the wrong thing; replace it.
  - The "137 964 → 1.96 M" sentence in D_grass §1.3 should be corrected.

### F4 (medium; good news that changes the risk ranking): locality (W3, W4) holds for 39 of the 40 chosen non-Id windows [V]

**W4 check (`w4.py`).** Φ = flows of W'_t'·ρ at a type-t coset whose configuration after |W'_t'| rungs has type t'; this is exactly `step_family`'s `cols`. Results:
- every f ∈ Φ has the target's P-part **and** Q-part;
- |Φ| = G(1);
- this holds for every `pick` window except 383 level 7: 0 of 6 columns are local. Pair 37's window also fails, but 37 is census-feasible.

**W3 check (`w3.py`).** The window parts have the local maximal types, both at c_loc and before ρ.

**Consequences:**
- D_grass Prop W now has all of W1-W4 for the 30 windows with G_loc = G. They are a-priori complete.
- D_grass risk 2 ("W3/W4 fail, medium") is retired.
- D_sinister asserted "global columns = window columns" with no argument. It is now [V].
- D_grass's W4 omits the **Q-suffix** condition: `Local` = id_P ⊗ s ⊗ id_Q. Every non-local variant I found has a suffix (45, 335 and 285 joint, 383, 37). Add the Q-part to W4.

**Group C ([V], `gloc.py`).** The group-C windows are local too:
- 410 (3 583), 283/293 (58), 313 (484), 137/145 (6 320).
- In each, G_loc = G + higher terms, and the local columns are a strict superset of Φ. Examples:
  - 410: v⁰ against v⁰ + v²;
  - 313: v⁰ + v² against v⁰ + v² + v⁴.
- By light-leaf triangularity, the lowest leaves should restrict invertibly [K]. So the group-C rank tests are more likely to pass than "medium risk" suggests.

### F5 (medium): D_sinister's template-cost formula is wrong per direction

- **The error.** `interpolate` (leafgen.rs:323-334) and `poly_certificate` (tmpl.rs:728) use M = monomials at the **maximal** entry degree d − sh(src) − sh(tgt). That is symmetric in the direction.
- `tcost.py` uses d − 2·sh(b) for In and d − 2·sh(aρ) for Out.
- **Effect.** In is underestimated and Out overestimated [V] (`mtrue.py`):

  | step | In, D_sinister | Out, D_sinister | correct (both directions) |
  |---|---|---|---|
  | 458 | 6.5 s | 1 900 s | 204 s (M = 3 029) |
  | 385/467 | 16 s | 2 785 s | 154 s |
  | 137/145 | 21 s | 744 s | 183 s |
  | 313 | 17.5 s | 235 s | 54 s |

- **Net effect:**
  - all 39 templates are ≤ 204 s by the B2 fit;
  - 458 is *not* over the cap;
  - the maximum label over all levels of both window words is ≤ 6, so GrTable(r, ≤ 6) suffices. The claim holds [V].
- **Caveats:**
  - the fit's 5e-7·M·n² term ignores the n³ inverse of L_tgt at each of the M + 2 points;
  - for n = 720 (137/145, In direction) expect about 3×, roughly 10 min;
  - the fit was calibrated on n ≤ 225.

### F6 (medium): the GrTable(·, 9) risk moved into the chains

- `PlanOpts` (rexplan.rs:308) has no label guard. A Collapse into `canonical_max` of a 3-upright, 9-root window creates a level of label 9.
- The SBim build that `certify_added` needs then waits on GrTable(r, 9), which ran > 90 min unfinished. B2's chains avoided this only because their templates were ≤ 90 flows.
- B5 adds many new chains (detours at levels 3-6, per-step exits). **Add `max_label ≤ 8`** to `PlanOpts` and to `joint_pair`. GrTable(r, 8) is 3-7 min once and persisted.

### F7 (medium): record runs are undercounted

`CostTable::predict_exact` (homprov.rs:1104) keys records by the exact (labels, fs, ft, unknowns, deg). So each direction and each degree above 10⁵ needs its own record [V]:

| window | degrees above 10⁵ | records |
|---|---|---|
| 313, 283 842 | d = 0 (157 488), d = 2 | 4 |
| 137/145 | d = −1 (167 402), d = 1 (318 975) | 4 |
| 385/467 | d = 4 (226 585), d = 6 (389 029) | 4 |
| 326 D2 | one | 2 |
| each D1 step, if the closed form fails | one | 2 |

- That is about 14-20 runs at 1.6-3.9·10⁵ unknowns, not "two records".
- For scale: A2's 434 787-unknown run took 179-272 s and 2.3-3.2 GB, and 92 of 136 M5b records in 10⁵-10⁶ were censored.

### F8 (medium): detours and adoption at intermediate types are untested

- Only 10 of 42 detours fit a single Exact call ≤ 10⁵. The planner has only run at final types q (118 exit chains).
- Adopting W'_t changes every step into and out of t. With up to 29 non-Id steps per target, that coupling is the second-largest unknown after F1.
- B5a must dry-run `Planner::plan` for every detour, with the label guard.

### F9 (low-medium): D_sinister specifics

- **`msearch.py` columns.** Its T is a local proxy: the maximal local type. It matches the global Φ for every ρ-last window except the one-sided 283/293 windows (4 local columns against 3 global) [V]. Its 13 PullTw successes are local-rank facts.
- **Flag-target lemma.** The ZMT argument is fine: finite fibres over fixed points give quasi-finiteness by semicontinuity, so the map is an isomorphism onto a normal Schubert variety, and restriction from G/P is surjective. Genuineness is decided by `certify_added` anyway.
- **What PullTw buys.**
  - 313 at 283 842: avoids 4 records;
  - 326 D2: avoids 2 records;
  - 458, 82, 285/385: S2, which B1 already has.

  So PullTw is optional, not the primary route.
- **Pair 383.** Its non-locality is confirmed and stronger than stated: 0 of 6 global columns are local. It needs `Pad` (local window 690) or it stays a gap.

### F10 (effort): 3-4 sessions is not realistic

B5a in Rust must do all of the following:
- re-derive the joint search and the W checks;
- run an all-steps census from the executor's plan;
- give per-step final words;
- dry-run the detours;
- dry-run the rank tests.

After that come about 20 record runs, V0 on small λ, both primes, and, for D_sinister, `wsl.rs` plus a new `CfKind`.

B2 did a comparable amount in one session, but with no hidden-step surprises. Estimate 5-7 sessions to the stopping point.

## 2. What holds [V]

| claim | status |
|---|---|
| The 42 picked windows reproduce their unknowns. | Holds: 0 mismatches against `hu.py`. |
| The picked W'_t' and W'_t are reduced for t' and t. | Holds for all 42. |
| Junction words stay small. | Holds: ≤ 38 880 flows. |
| Window labels stay small. | Holds: max label ≤ 6 at every level. |
| Group E steps are Id. | Holds. |
| 383 is non-local. | Holds. |
| Template costs are within the cap. | Within the cap once corrected (F5). |
| Certification is unchanged in both designs. | Holds: `certify_pairs` and π∘ι; PullTw templates are `Added` and re-certified. |

## 3. Required changes before B5b (the B5a gate)

1. **`gp0 --all-steps`.** Take the steps from `build_leaves`' own plan, after the alt-final swap. Print every step (level, t' → t, words, G, window, whole) and every final predecessor. Fix the stale `worst` string (gp0.rs:174-177).
2. **Per step, record:**
   - the best local window (W1-W4 *with* the Q-part);
   - G_loc;
   - the closed-form shape;
   - the records needed (direction × degree);
   - the final word per predecessor.
3. **Dry runs.** Run the rank tests for windows ≤ 10⁵, and `Planner::plan` (with `max_label ≤ 8`) for every detour and per-step exit.
4. **Gate on all steps.**
   - **Go for N = 5** if 313, 410 and 326 each have every step on an in-budget route.
   - **N = 6:** take only the (pair, q) groups that are fully in budget.
   - **Stop and report** otherwise. Use `Pad` only for 383 or as the last N = 5 blocker.

## 4. Reproduction

`T=/private/tmp/claude-501/-Users-existair-Documents-KRcalculator/bf6f8870-08ad-436b-9ff7-0414655cc59a/scratchpad/colored_design/b5/tools_c`

| finding | command |
|---|---|
| F5 | `python3 $T/mtrue.py` |
| F4 | `python3 $T/w4.py`; `python3 $T/w3.py`; `python3 $T/gloc.py` |
| F1 | `python3 $T/finals.py 0 60`; `python3 $T/shared.py 1,1,1/2,0,1/0,2,1 F1.2 "1,2,0/2,0,1/0,3,0;1,2,0/2,1,0/0,2,1" 5 6` |
| F1, F3 | `python3 $T/hid410.py 1,2,0/2,0,1/0,1,2 F1.1 1,2,0/2,1,0/0,1,2 5 6` |
