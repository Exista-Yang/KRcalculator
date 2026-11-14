# B5 census analysis: the worst Lemma A step of every non-rex (3,3,3) pair target, N = 5 and N = 6

Date 2026-11-12. Read-only analysis of branch `m4c-rexplan` (4e886f6) in `.claude/worktrees/revrex`. No cargo, no
binaries. Inputs: `docs/colored_design/m4c/b2/gp0_n5b.jsonl`, `gp0_n6.jsonl` (field `worst_step`), `gp0_n6.log`,
M4C_PLAN §3.1/§3.3, D1 §1.7-1.8, A1/A2/B1/B2 reports, `src/leaves.rs` (`step_score`, `step_family`, `plan_words`),
`src/dcoset.rs` (`step_terms`, `grass`, `reduced_for`, `reduced_words`), `src/homprov.rs` (`win_split`,
`hom_unknowns`), `src/bin/gp0.rs`, `tools/c_tools/hu.py`.
Tags: [V] computed here with the Python replica; [K] inference to be confirmed by a gate; [E] estimate.

## 0. Verdict (go / no-go)

**Go for B5, without `Pad` in the primary route, behind one cheap census gate (B5a).**

1. **No worst step needs a full 9-root block.** Every worst step is a rung step whose merged K column has 2 or 3 row
   cells; in 40 of 42 distinct steps it has 3 cells and the other column is non-empty, so the D1 §1.7 chart
   Z·T_K needs a 4-block level (Pad). The (b) case "K holds all 9 roots" never occurs. [V]
2. **The rung-local K-windows are tiny**: 45-57 666 unknowns, median 1 274, 41 of 42 ≤ 10⁴, and G_local = G_global
   on all 46 (step, word) variants (SLL Lemma 3.9, as A1's G-L1). [V] But they need `Pad` plus 4-upright rex
   detours W ↔ Z·T_K whose single-Exact bounds are mostly 10⁶-10⁹ (the planner has never run on λ·0).
3. **New: the 3-upright joint word search closes the gap without `Pad`.** A1 open problem 1 ("vary source and target
   junction words jointly") was never run. Choosing W'_{t'} and W'_t together, with a common prefix equal to the
   Grassmannian prefix Z with its pad rungs dropped (the other column merged into the last K cell on upright 2),
   gives step windows ≤ 10⁵ for 37 of 42 steps and ≤ 4·10⁵ for all 42 [V]. Examples: 971 748 783 → 34 696
   (pair 293), 927 678 854 → 295 452 (285/385), 217 411 216 → 6 320 (137/145, rank test needed), 132 967 552 → 18 461 (45),
   13 109 102 → 29 841 (273/278), 1 548 184 → 13 051 (127).
4. **What remains open** (all gate-checkable, none decided here):
   - (i) the census lists only the worst step; 12 targets are known to have a second step > 10⁵ (§5.3), others may;
   - (ii) the joint words differ from the plan's words, so each moved step needs rex detours (degree-0 rex chains
     between reduced words of one type, B2's planner on 3 uprights); only 10 of 42 detour pairs fit one Exact call
     ≤ 10⁵, the rest need planner chains;
   - (iii) window completeness: 30 steps have a window ≤ 10⁵ with G_local = G (a priori complete [K]); the
     group-C windows (19 targets) need the rank test, with a `Pad` fallback whose local windows are 58-3 583;
   - (iv) group D2 (9 targets: 326, 385, 467) stays at 113 288 / 389 029 unknowns: two exact-shape records needed.
5. **Honest outlook [E].** With B5a (census) + B5b (detours, joint words, closed-form Pull) + runs:
   - N = 5 complete (313, 326, 410): ≈ 60-70 %;
   - N = 6: the 22 "clean" pairs ≈ 70 %, all 32 non-feasible pairs ≈ 35-45 %;
   - effort 3-4 sessions, +1-2 if the `Pad` fallback (B3) is needed.

   Reasonable stopping point: after B5a, certify the clean groups, measure D2, and stop. Any pair still failing
   then stays as a documented gap; `Pad` only if it is the last blocker of N = 5.

## 1. Method and replica validation

- **Replica** (`tools/kr.py`): line-by-line Python of `dcoset::step_terms` (SLL Def 3.5 defects, Gaussian-binomial
  poly factor), `forward`, `reduced_for` (c_t = 1, support ≤ t by `min_perm` + tableau Bruhat), `grass` (SLL Thm 2.6
  at rung level), `reduced_words` (BFS, cap), `homprov::win_split`, and `hu.py` for `hom_unknowns`.
- **Check [V]: the census `step_score` was reproduced for all 111 non-rex targets (0 mismatches).** That covers the
  worst string, supp G, both directions, window vs whole, and the alt-final substitution.
- **Per worst step:**
  - the level, t' → t, the words, ρ and G;
  - the per-degree (whole, window) unknowns and the window labels;
  - Q = merged table, K cells (`grass` `local_lam`), other column;
  - the local tables m' → m'' and G_local (always = G);
  - the best local K-window over local reduced words of ≤ 5 rungs.
- **4-upright chart** (`tools/chart2.py`):
  - the chart word is C_K(t') = Z·T'_K(·E2.o), with Z found by BFS on λ·0 = (3,3,3,0) for the z table of `grass`;
  - K = {1,2} steps were analysed in the mirror frame (uprights reversed, E↔F, tables reversed). Unknowns are
    mirror-invariant; in the real frame Z is 1-2 rungs longer;
  - the detours W_{t'}·0 → C_K(t') and C_K(t) → W_t·0 (or L_q) are bounded by one Exact Hom⁰.
- **3-upright joint search** (`tools/joint_all.py`, `seeded_all.py`, `seeded_ok.py`, `seeded_hard.py`, `jfilter.py`):
  - candidates: all reduced words of t' (≤ 5 rungs) × t (≤ 6 rungs), BFS cap 2.5·10⁵;
  - for types whose words have 6-7 rungs (BFS too shallow in Python), seeded words P·A, with P = the projected
    Grassmannian prefix (Z minus the rungs on upright 3), suffixes ≤ 4/5 rungs;
  - score = `step_score`. **G_local test:** the window after the common prefix must be a single local step m' → m''
    (local types on the window labels) with `step_gf` = G.
- **Closed form (B1 kinds)** for the chosen window: Pull (Out) / Push (In) when G is one monomial v^d with
  coefficient 1, one window side is flag-type (all standard coefficients ≤ 1 at v = 1), and d = shift(Y) − shift(X).
  This is a candidate: `leafgen::build` (interpolate + `certify_added`) decides.
- **Run times.** Two runs took 47 s and 58 s (one core, < 200 MB), above the 20 s guideline; every other run was split
  into ≤ 18 s batches with pickled caches.

## 2. What the worst steps look like

- **Coverage.** 101 non-feasible non-rex targets in 50 (N, pair, q) have **40 distinct worst steps** (t', t, ρ). The
  plan is per (pair, q), so all s of a q share one step. Two more steps occur only on census-feasible targets
  (37, 285-q3) because the alt-final word already fixed them: 42 rows in §4.
- **Levels.**
  - 22 of the 42 steps are **final steps** (t = q, level 5-7). They are freely re-choosable: the exit chain
    W_q → L_q absorbs the target word.
  - 20 are **non-final** (levels 3-6).
- **Windows.** All 42 have census windows on all 9 roots:
  - 19 have a window on (3,3,3), i.e. no usable common prefix or suffix;
  - the others sit on (3,6,0), (4,5,0), (5,1,3), (1,5,3), (4,2,3), (2,4,3), (2,6,1) or (6,3,0).
  - Whole-word sizes are 9.6·10⁵-1.4·10¹⁰.
- **D1 §1.7 class:**

  | class | meaning | steps | targets |
  |---|---|---|---|
  | (a) | K column has 2 cells: Z on 3 uprights (a row swap E1.3·F1.3), no pad | 2 (281, 381; + 127's feasible level-2 step) | 2 |
  | (a-pad) | K column has 3 cells and the other column is non-empty: Z level (c0, c1, c2, o) has 4 blocks | 40 | 99 |
  | (b) | K column holds all 9 roots (local = global) | **0** | 0 |
  | (c) | rex-type: G = 1 and W_{t'}·ρ reduced for t | 3 (219, 319, 332), a subset of (a-pad) | 8 |
- **Why the census windows are big.** The plan's junction words are chosen one type at a time: prefix extension,
  canonical ladder, then a ≤ 5-rung descent with the neighbour fixed. For a 3-cell K column with a non-empty other
  column, no single-sided change can create a common prefix. A **joint** change can: the common prefix is the
  Grassmannian prefix with its pad rung dropped. Below, F2.x is the rung that moves the other column to the pad
  upright.

  | pair | t' → t | 4-upright Z | 3-upright prefix P | window |
  |---|---|---|---|---|
  | 127 | 1,1,1/2,1,0/0,1,2 → 2,0,1/3,0,0/0,1,2 | F0.1,F1.1,F2.3 | F0.1,F1.1 | 13 051, labels (2,3,4) |
  | 137/145 | 1,1,1/1,1,1/1,1,1 → 2,0,1/2,0,1/1,1,1 | F0.1,F1.2,F2.3 | F0.1,F1.2 | 6 320 |
  | 293 | 2,0,1/2,0,1/1,1,1 → 1,1,1/1,1,1/1,1,1 | F0.1,F1.2,F2.3 | F0.1,F1.2 | 34 696 |

  This is a partial Grassmannian reduction on 3 uprights: upright 2 holds c₂ ∪ other.

## 3. Group table (by step shape and remedy)

Unknowns are max over supp G of min(whole, window), both directions (the census `step_score`).

| group | step shape | distinct steps | (N, pair, q) / targets | pairs | census worst step | best 3-upright window | 4-up chart local K-window | B1 closed form | B5 remedy |
|---|---|---|---|---|---|---|---|---|---|
| **A** | K = 2 cells (3,3), Z = row swap E1.3·F1.3 on 3 uprights | 2 | 2 / 2 | 281, 381 | 4.2-6.3 M | 3 310-4 603 (G_loc = G) | same (3 310-4 603) | **Pull/Push** (G = v¹, v⁴; E1.3·F1.3·E0.3 flag-type) | chart words; detours 3 184 / 34 189 (single Exact) |
| **B** | 3 K cells + other column; projected-Z prefix gives a complete-shaped window | 24 | 28 / 57 | 326(N5,N6), 24, 45, 78, 82, 127, 137, 145, 212, 227, 273, 278, 283, 285, 292, 293, 335, 336, 376, 383, 385, 411, 418 | 0.2 M-972 M | **238-78 371, G_loc = G** | 45-57 666 | Pull/Push on 8 steps (G one monomial); otherwise exact window | joint words + detours (B2 planner) |
| **C** | as B, but the ≤ 10⁵ window has G_loc ≠ G (local types coarser) or is not of the form A·ρ → B | 6 | 9 / 19 | 313(N5,N6), 410(N5,N6), 137, 145, 283, 293, 383 | 0.14 M-772 M | 58-6 320 (rank test needed); G_loc = G alternatives 283 842 (313), 318 975 (137/145), none ≤ 4·10⁵ (410, 283, 293, 383) | 58-3 583 | none (G has 2-4 terms) except 410 (G = v⁰, not flag-type) | rank test in the census; fallback: G_loc = G window with a record, or Pad chart |
| **D1** | 3-upright best in 1-4·10⁵, G one monomial | 3 | 4 / 6 | 82, 285, 385, 458 | 2.5 M-928 M | 194 000-295 452 | 510-5 265 | **Pull/Push candidate** (G = v¹, v¹, v⁵) | closed-form template instead of the exact call |
| **D2** | 3-upright best in 1-4·10⁵, G two monomials | 2 | 4 / 9 | 326(N5,N6), 385, 467 | 19 M-22.7 M | 113 288 (326 final), 389 029 (385/467) | 1 274-1 397 | none (G = v¹+v³, v⁴+v⁶) | measured exact-shape records (two shapes); Pad chart fallback |
| **E** | rex-type step: G = 1, W_{t'}·ρ is reduced for t | 3 | 3 / 8 | 219, 319, 332 | 0.70-12.7 M | 0 (Id) after W_t := W_{t'}·ρ | 0 | rex chain (Collapse / Expand / Assoc / FlagMatch / Swap) | `StepKind::Rex`: realize the B2 planner chain W_{t'}·ρ → W_t, no exact call |

**Totals.** 101 targets / 50 (N, pair, q) / 40 steps:
- A + B + D1 + E (closed or certifiable with no new trust): 73 targets;
- C (rank test): 19 targets;
- D2 (records): 9 targets.

**Per pair** (worst-step level only):
- **N = 6.** 22 of the 32 non-feasible pairs are "clean" (A/B/D1/E only): 24, 45, 78, 82, 127, 212, 219, 227, 273,
  278, 281, 285, 292, 319, 332, 335, 336, 376, 381, 411, 418, 458. Seven need C: 137, 145, 283, 293, 313, 383, 410.
  Three need D2: 326, 385, 467.
- **N = 5.** 313 → C; 410 → C; 326 → B + D2 (+ an unanalysed second final step of 147 091, §5.3).

Local K-window sizes (D1 §1.7) are listed for comparison: they are what the `Pad` route would cost per step. They
show that a small factored step exists mathematically for every target.

## 4. Per-step table (all 42 distinct worst steps)

- lv = level of the step; F = final step.
- "census whole / step" = max-degree whole-word unknowns / census `step_score`.
- K cells are the merged K column of Q = merge_columns(t', i) (mirror frame where marked).
- Detour bounds are one Exact Hom⁰ (min of whole and window, both directions): in = W_{t'} → W'_{t'};
  out = W'_t → W_t, or W'_q → L_q for final steps. A planner chain is usually much cheaper than this bound.

| grp | pairs (targets) | lv | t' → t | W_t'·ρ → W_t (census words) | G | census whole / step (window labels) | Q: K column cells, other col | 4-up chart local K-window | 3-up joint best (labels) | 3-up best with G_loc = G | detour Exact bound in / out | B1 closed form for the chosen window |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | 281 (1) | 4 | 3,0,0/0,0,3/3,0,0 → 1,2,0/0,0,3/2,1,0 | [E1.3,E0.3,F1.3]·F0.3 → [E1.3,E0.2,F1.3,F0.2] | v^1 | 25 327 405 / 6 342 787 [3, 6, 0] | K={0,1}, cells [3, 3], other 3 | 4 603 | 4 603 [3, 3] | 4 603 | 3 184 / 34 189 | Pull(Out)/Push(In) |
| A | 381 (1) | 4 | 3,0,0/0,0,3/3,0,0 → 2,1,0/0,0,3/1,2,0 | [E1.3,E0.3,F1.3]·F0.3 → [E1.3,E0.1,F1.3,F0.1] | v^4 | 17 553 206 / 4 221 776 [3, 6, 0] | K={0,1}, cells [3, 3], other 3 | 3 310 | 3 310 [3, 3] | 3 310 | 3 184 / 34 189 | Pull(Out)/Push(In) |
| B | 293 (2) | 6 | 2,0,1/2,0,1/1,1,1 → 1,1,1/1,1,1/1,1,1 | [E0.2,E1.3,E0.1,F0.1,F1.3]·F0.2 → [E0.1,F0.1,E1.2,E0.1,F1.2,F0.1] | v^1 | 971 748 783 / 971 748 783 [3, 3, 3] | K={0,1}, cells [2, 2, 2], other 3 | 4 695 | 34 696 [2, 3, 4] | 34 696 | 3 074 782 / 919 495 095 | Pull(Out)/Push(In) |
| B | 336 (2) | 7F | 0,2,1/2,0,1/1,2,0 → 0,1,2/2,0,1/1,2,0 | [E0.2,E1.3,E0.1,F1.1,F0.3,F1.1]·F1.1 → [E0.2,F0.2,E1.3,E0.1,F1.3,F0.1] | v^-1 + v^1 | 14 342 746 161 / 853 377 920 [5, 1, 3] | K={1,2}, cells [2, 1, 3], other 3 (mirror) | 752 | 752 [1, 5] | 752 | 884 / 477 989 990 | none (G has >1 term: Pull∘Dot sprinkles [K]) |
| B | 292 (3) | 7F | 0,1,2/2,0,1/1,2,0 → 0,2,1/2,0,1/1,2,0 | [E0.2,F0.2,E1.3,E0.1,F1.3,F0.1]·E1.1 → [E0.2,E1.3,E0.1,F1.1,F0.3,F1.1] | v^-2 + v^0 | 7 973 688 874 / 477 989 990 [5, 1, 3] | K={1,2}, cells [2, 1, 3], other 3 (mirror) | 394 | 394 [1, 5] | 6 073 | 477 989 990 / 2 413 242 | none (G has >1 term: Pull∘Dot sprinkles [K]) |
| B | 82 (1) | 6F | 1,0,2/2,1,0/1,1,1 → 1,0,2/1,2,0/1,1,1 | [F0.2,E1.1,E0.3,E1.1,F1.2]·F0.1 → [F0.2,E0.1,E1.2,E0.1,F1.2] | v^0 + v^2 | 3 051 105 496 / 183 104 704 [1, 5, 3] | K={0,1}, cells [1, 3, 2], other 3 | 3 968 | 27 192 [1, 4, 4] | 27 192 | 1 425 066 / 38 334 726 | none (G has >1 term: Pull∘Dot sprinkles [K]) |
| B | 418 (1) | 6F | 1,2,0/1,1,1/1,2,0 → 1,1,1/1,0,2/1,2,0 | [E0.1,E1.3,F1.1,F0.2,E0.1]·F1.2 → [E0.1,F0.1,E1.3,E0.1,F1.3,F0.1] | v^0 + v^2 | 2 104 672 327 / 182 643 095 [4, 2, 3] | K={1,2}, cells [2, 2, 2], other 3 (mirror) | 6 239 | 6 239 [4, 2] | 6 239 | 16 028 343 / 708 247 147 | none (G has >1 term: Pull∘Dot sprinkles [K]) |
| B | 45 (2) | 6 | 1,0,2/2,0,1/1,2,0 → 1,0,2/1,1,1/1,2,0 | [F0.2,E0.2,E1.3,E0.1,F1.3]·F0.1 → [F0.2,E0.1,E1.3,E0.1,F1.3] | v^1 | 2 239 829 737 / 132 967 552 [1, 5, 3] | K={0,1}, cells [1, 2, 3], other 3 | 539 | 18 461 [1, 3, 5] | 18 461 | 15 819 429 / 74 745 807 | Pull(Out)/Push(In) |
| B | 45 (2) | 6 | 1,1,1/2,0,1/1,1,1 → 1,1,1/1,1,1/1,1,1 | [E0.2,E1.3,E0.1,F0.1,F1.3,F0.1]·F0.1 → [E0.1,F0.1,E1.2,E0.1,F1.2,F0.1] | v^0 | 1 540 351 073 / 68 972 978 [3, 3, 3] | K={0,1}, cells [2, 2, 2], other 3 | 6 963 | 11 [2, 2] | 51 816 | 137 424 990 / 919 495 095 | none (flag 00, d=0, shifts -10,-8) |
| B | 292 (3) | 6 | 1,0,2/3,0,0/0,2,1 → 1,0,2/2,1,0/0,2,1 | [E0.3,F0.2,E1.3,F1.3]·F0.1 → [E0.2,E1.2,F0.2,F1.2] | v^1 | 40 165 062 / 40 165 062 [3, 3, 3] | K={0,1}, cells [1, 3, 2], other 3 | 45 | 18 461 [1, 3, 5] | 18 461 | 1 971 122 / 37 736 389 | Pull(Out)/Push(In) |
| B | 78 (1) | 5F | 2,0,1/1,1,1/1,1,1 → 2,0,1/0,2,1/1,1,1 | [F0.1,E0.1,E1.2,E0.1,F1.2]·F0.1 → [F0.1,E1.2,E0.1,F1.2] | v^1 + v^3 | 246 281 307 / 21 134 393 [2, 4, 3] | K={0,1}, cells [2, 2, 2], other 3 | 1 274 | 13 051 [2, 3, 4] | 13 051 | 3 319 170 / 1 915 652 | none (G has >1 term: Pull∘Dot sprinkles [K]) |
| B | 273, 278 (4) | 4 | 3,0,0/1,0,2/1,1,1 → 2,1,0/0,1,2/1,1,1 | [E0.1,E1.3,E0.1,F1.3]·F0.2 → [E1.2,E0.1,F1.2,F0.1] | v^2 | 13 109 102 / 13 109 102 [3, 3, 3] | K={0,1}, cells [3, 1, 2], other 3 | 2 694 | 29 841 [3, 2, 4] | 29 841 | 59 769 / 5 291 812 | Pull(Out)/Push(In) |
| B | 273, 285 (3) | 5F | 1,2,0/0,1,2/2,0,1 → 1,2,0/0,2,1/2,0,1 | [E1.2,E0.2,F1.2,F0.2]·E1.1 → [E1.2,F1.1,F0.2,E0.2] | v^0 + v^2 | 171 284 288 / 10 430 174 [3, 5, 1] | K={1,2}, cells [1, 3, 2], other 3 (mirror) | 159 | 57 736 [5, 3, 1] | 57 736 | 2 331 249 / 1 891 814 | none (G has >1 term: Pull∘Dot sprinkles [K]) |
| B | 137 (3) | 6F | 1,1,1/1,2,0/1,0,2 → 2,0,1/2,1,0/1,0,2 | [F0.2,E1.1,E0.2,F1.1]·E0.2 → [E1.1,E0.3,F0.1,F1.1] | v^-4 + 2v^-2 + v^0 | 8 435 824 / 8 435 824 [3, 3, 3] | K={0,1}, cells [2, 3, 1], other 3 | 1 365 | 752 [2, 4] | 51 816 | 20 863 646 / 564 256 | none (G has >1 term: Pull∘Dot sprinkles [K]) |
| B | 283 (2) | 7F | 0,2,1/2,0,1/1,1,1 → 0,2,1/2,0,1/1,2,0 | [E0.2,E1.3,E0.1,F0.1,F1.3,F0.2]·E1.1 → [E0.2,E1.3,E0.1,F1.1,F0.3,F1.1] | v^-1 + v^1 | 3 387 028 370 / 4 089 813 [6, 3, 0] | K={1,2}, cells [2, 1, 3], other 3 (mirror) | 2 524 | 238 [4, 2] | 238 | 52 562 544 / 112 288 049 | none (G has >1 term: Pull∘Dot sprinkles [K]) |
| B | 145, 227 (9) | 6 | 1,2,0/2,0,1/1,0,2 → 1,2,0/2,0,1/0,1,2 | [E1.1,E0.3,F1.1,F0.2]·F0.1 → [E0.2,E1.3,F1.3,F0.2] | v^1 | 4 050 681 / 4 050 681 [3, 3, 3] | K={0,1}, cells [3, 2, 1], other 3 | 8 652 | 29 [3, 3] | 3 002 | 773 / 0 | none (flag 01, d=1, shifts -11,-6) |
| B | 383 (3) | 6 | 2,0,1/3,0,0/0,1,2 → 1,1,1/2,1,0/0,1,2 | [E0.3,E1.3,F0.1,F1.3]·F0.2 → [E0.2,F0.2,E1.1,F1.1] | v^1 | 4 050 681 / 4 050 681 [3, 3, 3] | K={0,1}, cells [2, 3, 1], other 3 | 922 | 34 696 [2, 3, 4] | 34 696 | 15 039 / 5 291 812 | Pull(Out)/Push(In) |
| B | 385 (3) | 6F | 2,1,0/0,1,2/1,1,1 → 2,1,0/0,3,0/1,1,1 | [E1.2,E0.1,F1.2,F0.1]·E1.2 → [F0.1,E1.2,E0.1] | v^-2 + v^0 + v^2 | 2 633 881 / 2 633 881 [3, 3, 3] | K={1,2}, cells [2, 3, 1], other 3 (mirror) | 517 | 3 571 [4, 4, 1] | 19 371 | 2 579 / 73 262 | none (G has >1 term: Pull∘Dot sprinkles [K]) |
| B | 127 (4) | 5F | 1,1,1/2,1,0/0,1,2 → 2,0,1/3,0,0/0,1,2 | [E0.2,F0.2,E1.1,F1.1]·E0.2 → [E0.3,F0.1,E1.3,F1.3] | v^-3 + 2v^-1 + 2v^1 + v^3 | 1 548 184 / 1 548 184 [3, 3, 3] | K={0,1}, cells [2, 3, 1], other 3 | 353 | 13 051 [2, 3, 4] | 13 051 | 5 291 812 / 15 039 | none (G has >1 term: Pull∘Dot sprinkles [K]) |
| B | 82 (1) | 3 | 1,2,0/0,3,0/0,1,2 → 2,1,0/1,2,0/1,0,2 | [F0.2,E1.1]·E0.3 → [F0.1,E1.1,E0.2] | v^1 + v^3 | 1 036 000 / 1 036 000 [3, 3, 3] | K={0,1}, cells [3, 3, 1], other 2 | 57 666 | 27 286 [3, 4] | 27 286 | 293 / 141 831 | none (G has >1 term: Pull∘Dot sprinkles [K]) |
| B | 326(N5), 326 (2) | 5 | 3,0,0/2,0,1/0,2,1 → 2,1,0/1,1,1/0,2,1 | [E0.2,E1.3,F1.2]·F0.2 → [E0.1,E1.2,F1.1,F0.1] | v^2 | 956 320 / 956 320 [3, 3, 3] | K={0,1}, cells [3, 2, 2], other 2 | 997 | 78 371 [3, 2, 4] | 78 371 | 2 244 / 300 838 | Pull(Out)/Push(In) |
| B | 376 (1) | 4 | 3,0,0/1,0,2/2,0,1 → 1,2,0/1,0,2/1,1,1 | [E0.1,E1.3,E0.2,F1.3]·F0.3 → [E0.1,E1.3,E0.1,F1.3,F0.2] | v^2 | 33 075 182 / 928 861 [4, 5, 0] | K={0,1}, cells [3, 1, 2], other 3 | 7 230 | 3 471 [4, 2] | 3 471 | 415 / 1 825 | Pull(Out)/Push(In) |
| B | 285 (1) | 6F | 2,0,1/0,2,1/1,1,1 → 2,0,1/0,3,0/1,1,1 | [F0.1,E1.2,F1.2,E0.1]·E1.1 → [F0.1,E1.1,E0.1,F1.1,E1.1] | v^-2 + v^0 + v^2 | 57 917 453 / 447 163 [2, 4, 3] | K={1,2}, cells [2, 3, 1], other 3 (mirror) | 690 | 210 [3, 3] | 4 588 | 1 915 652 / 5 498 | none (G has >1 term: Pull∘Dot sprinkles [K]) |
| B | 24 (2) | 5F | 2,0,1/0,2,1/1,1,1 → 2,0,1/0,2,1/2,0,1 | [F0.1,E1.2,E0.1,F1.2]·E0.1 → [F0.1,E1.2,E0.2,F1.2] | v^-1 + v^1 | 89 776 752 / 439 729 [2, 6, 1] | K={0,1}, cells [2, 2, 2], other 3 | 238 | 238 [2, 4] | 238 | 64 502 / 162 489 | none (G has >1 term: Pull∘Dot sprinkles [K]) |
| B | 212 (1) | 4 | 1,2,0/2,1,0/0,3,0 → 1,2,0/2,0,1/0,1,2 | [E0.2,E1.3,F0.2]·F1.3 → [E0.2,E1.3,F1.3,F0.2] | v^6 | 14 550 258 / 392 995 [5, 4, 0] | K={1,2}, cells [3, 1, 2], other 3 (mirror) | 3 310 | 1 516 [2, 4] | 1 516 | 749 / 19 764 | Pull(Out)/Push(In) |
| B | 335 (2) | 7F | 1,1,1/1,1,1/1,2,0 → 1,1,1/1,0,2/1,2,0 | [E0.1,F0.1,E1.3,F1.2,F0.1,E0.1]·F1.1 → [E0.1,E1.3,F1.2,F0.2,E0.1,F1.1] | v^-1 + v^1 | 4 029 903 339 / 372 384 [4, 2, 3] | K={1,2}, cells [2, 2, 2], other 3 (mirror) | 9 822 | 497 [4, 2] | 9 822 | 676 166 158 / 708 247 147 | none (G has >1 term: Pull∘Dot sprinkles [K]) |
| B | 411 (1) | 5F | 2,1,0/0,2,1/1,2,0 → 2,1,0/0,0,3/1,2,0 | [E1.3,E0.1,F1.1,F0.1]·F1.2 → [E1.3,E0.1,F1.3,F0.1] | v^0 + v^2 + v^4 | 46 175 082 / 200 667 [4, 5, 0] | K={1,2}, cells [2, 3, 1], other 3 (mirror) | 4 668 | 4 668 [1, 5] | 4 668 | 1 825 / 10 260 | none (G has >1 term: Pull∘Dot sprinkles [K]) |
| C | 293 (2) | 7F | 0,2,1/1,0,2/2,1,0 → 0,3,0/1,0,2/2,1,0 | [E0.1,F0.1,E1.3,E0.2,F1.3,F0.2]·E1.1 → [E1.2,E0.3,F1.2,F0.3,E1.1] | v^-1 + v^1 + v^3 | 7 906 160 352 / 771 786 049 [3, 3, 3] | K={1,2}, cells [1, 2, 3], other 3 (mirror) | 58 | 58 [2, 4] | none ≤ 4e5 | 801 113 158 / 17 645 | none (G has >1 term: Pull∘Dot sprinkles [K]) |
| C | 137, 145 (7) | 7F | 1,1,1/1,1,1/1,1,1 → 2,0,1/2,0,1/1,1,1 | [E0.1,F0.1,E1.2,E0.1,F1.2,F0.1]·E0.2 → [E0.2,E1.3,E0.1,F0.1,F1.3] | v^-3 + 2v^-1 + v^1 | 217 411 216 / 217 411 216 [3, 3, 3] | K={0,1}, cells [2, 2, 2], other 3 | 838 | 6 320 [2, 3, 4] | 318 975 | 919 495 095 / 3 074 782 | none (G has >1 term: Pull∘Dot sprinkles [K]) |
| C | 283 (2) | 7F | 0,2,1/2,0,1/1,1,1 → 0,3,0/2,0,1/1,1,1 | [E0.2,E1.3,E0.1,F0.1,F1.3,F0.2]·E1.1 → [E1.1,E0.3,F1.1,F0.3,E1.1] | v^-1 + v^1 + v^3 | 895 177 792 / 85 872 528 [3, 3, 3] | K={1,2}, cells [2, 1, 3], other 3 (mirror) | 58 | 58 [2, 4] | none ≤ 4e5 | 52 562 544 / 884 | none (G has >1 term: Pull∘Dot sprinkles [K]) |
| C | 383 (4) | 7F | 0,2,1/2,0,1/1,1,1 → 0,3,0/2,0,1/1,2,0 | [E0.2,E1.3,E0.1,F0.1,F1.3,F0.2]·E1.2 → [E1.1,E0.3,F1.1,F0.3,E1.2] | v^-3 + 2v^-1 + 2v^1 + v^3 | 895 177 792 / 85 872 528 [3, 3, 3] | K={1,2}, cells [2, 1, 3], other 3 (mirror) | 690 | 61 [3, 3] | none ≤ 4e5 | 52 562 544 / 2 413 242 | none (G has >1 term: Pull∘Dot sprinkles [K]) |
| C | 410(N5), 410 (2) | 5F | 1,2,0/2,0,1/0,3,0 → 1,1,1/2,0,1/0,2,1 | [E0.2,E1.3,F1.1,F0.2]·F1.2 → [E0.2,F0.1,E1.3,F1.3,F0.1] | v^0 | 34 282 850 / 2 116 515 [5, 1, 3] | K={1,2}, cells [3, 1, 2], other 3 (mirror) | 3 583 | 3 583 [2, 4] | none ≤ 4e5 | 19 764 / 923 496 | none (flag 00, d=0, shifts -9,-7) |
| C | 37 (2) | 6F | 1,0,2/1,2,0/1,1,1 → 1,0,2/2,1,0/1,1,1 | [F0.2,E0.1,E1.2,E0.1,F1.2]·E0.1 → [E0.2,F0.2,E1.2,F1.2,E0.1] | v^-1 + v^1 | 1 990 717 305 / 816 723 [3, 3, 3] | K={0,1}, cells [1, 3, 2], other 3 | 1 881 | 2 749 [3, 3, 3] | none ≤ 4e5 | 2 951 542 / 344 933 968 | none (G has >1 term: Pull∘Dot sprinkles [K]) |
| C | 313(N5), 313 (2) | 5F | 1,2,0/2,0,1/0,2,1 → 1,2,0/2,0,1/0,1,2 | [E0.2,E1.3,F1.2,F0.2]·F1.1 → [E0.2,E1.3,F1.3,F0.2] | v^0 + v^2 | 4 561 122 / 137 964 [5, 4, 0] | K={1,2}, cells [3, 1, 2], other 3 (mirror) | 484 | 484 [2, 4] | 283 842 | 0 / 19 764 | none (G has >1 term: Pull∘Dot sprinkles [K]) |
| D1 | 285, 385 (4) | 5 | 1,1,1/0,2,1/0,2,1 → 1,1,1/1,1,1/1,1,1 | [F0.1,E1.2,F1.2,F0.1]·E0.2 → [E0.1,F0.1,E1.2,E0.1,F1.2,F0.1] | v^1 | 927 678 854 / 927 678 854 [3, 3, 3] | K={0,1}, cells [2, 2, 2], other 3 | 4 695 | 295 452 [1, 3, 5] | 295 452 | 242 953 / 919 495 095 | Pull(Out)/Push(In) |
| D1 | 82 (1) | 5 | 1,2,0/2,1,0/1,1,1 → 1,1,1/2,0,1/1,1,1 | [E0.2,E1.2,E0.1,F0.2]·F1.2 → [E0.2,E1.3,E0.1,F0.1,F1.3,F0.1] | v^1 | 265 625 333 / 16 221 263 [5, 1, 3] | K={1,2}, cells [2, 1, 2], other 4 (mirror) | 510 | 194 000 [1, 2, 6] | 194 000 | 3 902 / 167 576 156 | Pull(Out)/Push(In) |
| D1 | 458 (1) | 4 | 3,0,0/2,0,1/1,0,2 → 2,1,0/1,1,1/0,1,2 | [E1.1,E0.3,F1.1]·F0.3 → [E0.1,F1.1,F0.1,E1.1] | v^5 | 2 468 710 / 2 468 710 [3, 3, 3] | K={0,1}, cells [3, 2, 1], other 3 | 5 265 | 210 002 [3, 2, 4] | 210 002 | 773 / 92 | Pull(Out)/Push(In) |
| D2 | 326(N5), 326 (2) | 6F | 1,2,0/1,1,1/1,1,1 → 1,2,0/1,0,2/1,1,1 | [E0.1,E1.2,E0.1,F1.1,F0.2]·F1.1 → [E0.1,E1.3,E0.1,F1.3,F0.2] | v^1 + v^3 | 264 300 190 / 22 654 063 [4, 2, 3] | K={1,2}, cells [2, 2, 2], other 3 (mirror) | 1 274 | 113 288 [5, 3, 1] | 113 288 | 0 / 6 425 | none (G has >1 term: Pull∘Dot sprinkles [K]) |
| D2 | 385, 467 (7) | 5 | 1,1,1/0,2,1/0,2,1 → 2,0,1/1,1,1/0,2,1 | [F0.1,E1.2,F1.2,F0.1]·E0.2 → [F0.1,E0.1,E1.2,F1.2] | v^4 + v^6 | 218 487 888 / 19 008 170 [2, 4, 3] | K={0,1}, cells [2, 2, 2], other 3 | 1 397 | 389 029 [2, 2, 5] | 389 029 | 242 953 / 292 969 | none (G has >1 term: Pull∘Dot sprinkles [K]) |
| E | 332 (1) | 4 | 1,2,0/0,1,2/0,3,0 → 1,2,0/0,1,2/2,1,0 | [E1.3,F1.2,F0.2]·E0.2 → [E1.3,E0.2,F1.2,F0.2] | v^0 | 49 493 015 / 12 738 410 [3, 6, 0] | K={0,1}, cells [3, 1, 3], other 2 | 0 | 0 (Id) | 0 | 0 / 12 738 410 | Id |
| E | 219 (5) | 4 | 1,2,0/1,2,0/0,2,1 → 1,2,0/1,2,0/1,1,1 | [F0.2,E0.1,E1.2]·E0.1 → [E1.1,E0.2,F0.2,E1.1] | v^0 | 2 951 542 / 2 951 542 [3, 3, 3] | K={0,1}, cells [3, 3, 2], other 1 | 0 | 0 (Id) | 0 | 0 / 2 951 542 | Id |
| E | 319 (2) | 4 | 2,1,0/0,1,2/0,3,0 → 2,1,0/0,1,2/1,2,0 | [E1.3,F1.2,F0.1]·E0.1 → [E1.3,E0.1,F1.2,F0.1] | v^0 | 2 405 483 / 701 771 [3, 6, 0] | K={0,1}, cells [3, 1, 3], other 2 | 0 | 0 (Id) | 0 | 0 / 701 771 | Id |

## 5. Implications for the B5 design (concrete)

### 5.1 B5a: census gate, combinatorics only (Rust, minutes; ~0.5-1 session)

- **`src/bin/gp0.rs`.**
  - Emit **every** step with `step_score` > budget, not only `worst`: (level, t', t, ρ, words, supp G, whole, window,
    labels). Today 12 targets are known to hide a second over-budget step (§5.3).
  - Add the joint search below and the rank test, and print the new verdict per target.
- **Joint junction-word search** (new `leaves::joint_pair(s, t', ρ, t, opts)`, called from `Leaves::plan_words` after
  the descent, for each step still > budget).
  - **Candidates** are `dcoset::reduced_words(t', ≤ 6) × reduced_words(t, ≤ 7)` (cap 4·10⁵; cheap in Rust), plus
    seeded words P·A.
  - **Seed P** is the projected Grassmannian prefix: `dcoset::grass` → z table on λ·0 → `bfs_reduced` → drop the
    rungs on upright 3. For K = {1,2}, mirror the step: on 3 uprights the mirror is exact.
  - **Score:** `step_score`. Prefer windows with **G_local = G**. Implement this as new
    `dcoset::window_gloc(lam, a, ρ, b) -> Option<Lp>`: the local types of the window words on the window labels,
    then `step_gf`.
  - **Final steps:** the chosen W'_q replaces W_q, and the exit chain W'_q → L_q is re-planned (already done for
    `alt_final`).
- **Rank test (B2 open problem 2)** for each chosen window ≤ 10⁵ unknowns. Run `Leaves::step_family` (Window kind) in
  a dry arena, one point: an exact window basis plus greedy selection on the top-to-top columns. This decides group C
  and confirms groups A and B before any production run.
- **Detours.** For every moved step, call `rexplan::Planner::plan` on
  - (W_{t'} → W'_{t'}) at type t', and
  - (W'_t → W_t) at type t, non-final steps only (and on L_q for final steps, already done);

  with B2's budget: templates ≤ 1 200 flows, Exact ≤ 10⁵, theorem Assoc only for up-down shapes. Record the chain
  length and costs.
- **Rex-type steps (group E).** When G = 1 and W_{t'}·ρ is reduced for t, plan W_{t'}·ρ → W_t; no exact call.
- **Exit criterion (G-P5).** Per target, every step either:
  - (1) ≤ 10⁵ with a passing rank test, or a certifiable closed form; or
  - (2) is a 1-4·10⁵ exact call with a measured record; or
  - (3) is a planner chain within budget.

  Targets that fail go to §5.4.

### 5.2 B5b: implementation in leaves.rs (no frozen file; ~1-1.5 sessions + runs)

- **`Leaves::step_family`: new `StepKind`s.**
  - **`Rex`** (group E): `leafgen::realize_chain` of the planner chain W_{t'}·ρ → W_t (degree 0, 1⊗ ↦ 1⊗). Its
    top-to-top matrix is the 1×1 scalar 1, which is valid since G = 1.
  - **`Detour`** (groups A-D):
    - Out: S = R_out ∘ S_win ∘ Kron(R_in, [ρ]);
    - In: S^∨ = Kron(R_in⁻¹, [ρ]) ∘ S_win^∨ ∘ R_out⁻¹, with R⁻¹ the reversed chain (Collapse ↔ Expand; Assoc,
      FlagMatch and Swap are symmetric; reverse Exact moves are their own windows).

    **Lemma A validity:** R_in and R_out are degree-0 maps between reduced words of one type with top scalar 1, and
    Kron is block-diagonal over prefix configurations. So top-to-top(S) = T(S_win)·diag(R_in(c′)), which is
    invertible iff T(S_win) is (same argument as M4C_PLAN §3.1; SLL Prop 4.19, Thm 6.11).
  - **`ClosedForm`** (groups A, D1 and 8 B steps): when G = v^d with coefficient 1 and the window satisfies the Pull
    or Push shape, use `leafgen::build` (Pull for Out, Push for In) on the window as a `Local` template.
    `certify_added` decides; if it fails, fall back to the exact window (≤ 4·10⁵ with a record in D1).
- **`Leaves::plan_words`:** call `joint_pair` for steps still above budget. Optionally adopt W'_{t'} / W'_t as the
  junction words themselves when the neighbouring steps stay within budget (no detour); otherwise use `Detour`.
- **Wiring:** `leaves::pairs_leaves` (exists) with the new step kinds. π = exit ∘ LL and ι = LL^∨ ∘ entry, then
  `homprov::certify_pairs` unchanged (count, form rank at rep(q) at 2 points, genuineness), then π∘ι = id blockwise
  at every coset of L_q, as in B4.
- **Nothing new is trusted.**
  - No new `MapNode`: only Comp, Kron, Local, Id and Elem-Assoc, as in B2/B4.
  - Exact windows are `hom_space` templates. Closed forms are Added templates, re-certified on import.
- **Records (group D2).** `m4c measure` the 113 288 shape (326 final, d ∈ {1,3}) and the 389 029 shape (385/467,
  d ∈ {4,6}). They are comparable to A2's 111 497 (47-66 s, 0.7 GB) and 434 787 (179-272 s, 2.3-3.2 GB) records.
- **Order.**
  - N = 5 first: 326 (B + D2), 313 (C), 410 (C).
  - Then N = 6: the clean 22 pairs, then C, then D2.
  - Both primes. Downstream checks (sindec wiring, d∘d, K₀, χ, replay) follow B2 §6 item 3.

### 5.3 Known hidden steps (the census reports only the worst step)

**Targets with a second step above 10⁵ are known** from `max_nonfinal` vs the final step and from `alt_final`:

| target | second step |
|---|---|
| 326 q 2,1,0/1,0,2/0,2,1 (N = 5 and 6) | final step 147 091 after the alt-final word |
| 145 q 2,0,1/2,0,1/1,1,1 | non-final 68 972 978: **the same step shape as pair 45's** (1,1,1/2,0,1/1,1,1 → 1,1,1/1,1,1/1,1,1), resolved here at 51 816 (G_loc = G) |
| 137 q 2,0,1/2,0,1/1,1,1 | non-final 91 724 472 |
| 137 q 2,0,1/2,1,0/1,0,2 | non-final 1 483 163 |
| 293 q 0,3,0/1,0,2/2,1,0 | non-final 1 192 238 |
| 385 q 2,1,0/0,3,0/1,1,1 | non-final 1 685 822 |
| 385 q 1,2,0/1,1,1/1,2,0 | alt final 40 312 945 |
| 385 q 2,0,1/1,2,0/0,3,0 | alt final 1 169 288 |
| 292 q 0,2,1/2,0,1/1,2,0 | non-final 356 983 |
| 285 q 1,2,0/0,2,1/2,0,1 | non-final 151 774 |
| 411 | non-final 135 902 |

Other targets may hide more (only the maximum is in the jsonl). The B5a census must list all of them. These steps
have the same structure (3-column levels, 3-cell K columns), so the joint search is expected to apply [E]; this is
the largest source of uncertainty in the per-pair estimates.

### 5.4 Fallbacks, in order

1. **Another candidate from the joint search** (keep the top k per step), or a deeper search (≤ 7/8 rungs).
2. **A record** for a 1-4·10⁵ window with G_loc = G:
   - 313: 283 842;
   - 137/145: 318 975.
3. **`Pad` (B3, approved "if strictly needed") with the 4-upright chart C_K(t) = Z·T_K(·E2.o).**
   - Local windows are 45-57 666, with G_local = G by Lemma 3.9, so they are complete a priori [K, rung-level
     transfer as D1 §1.7].
   - The cost moves to 4-upright rex detours: their single-Exact bounds are 10⁶-10⁹, and `rexplan` has never been
     run on λ·0 (gp0's `pl4`).
   - Only worth it if it is the last blocker of N = 5 (313 or 410), or of several N = 6 pairs.
4. **Leave the target uncertified.** It is reported as Refused; certification semantics are unchanged.

## 6. Risks and how each is detected

| # | risk | likelihood [E] | detection | consequence |
|---|---|---|---|---|
| 1 | hidden second steps do not yield to the joint search | medium | B5a census (all steps listed) | those targets → fallbacks §5.4 |
| 2 | detour chains (3 uprights) not found within budget; 32 of 42 detours exceed one 10⁵ Exact | low-medium (B2 found all 118 exit chains on 3 uprights, ≤ 6 moves, Exact ≤ 99 175) | B5a planner runs | adopt the joint words in the plan instead, or the next candidate |
| 3 | G_loc ≠ G windows (group C) fail the rank test | medium | B5a rank test | record of the G_loc = G alternative, or Pad |
| 4 | G_loc = G windows still fail the rank test (the [K] transfer for non-Grassmannian prefixes) | low | same | as 3 |
| 5 | Pull/Push candidates not genuine for the chosen windows | low-medium | `certify_added` (`interpolate` returns Err) | exact window: D1 then needs records of 194 000-295 452 unknowns |
| 6 | D2 records censored (1-4·10⁵ calls: 92 of 136 M5b records in 10⁵-10⁶ were censored) | medium | `m4c measure` | Pad chart for 326 (local 1 274) and 385/467 (local 1 397) |
| 7 | downstream (M6a) cost of deeper DAGs (detours add 2 chains per moved step) | low-medium | G-L5 counters | compaction; Leaves last |

## 7. Appendix: chosen 3-upright junction-word pairs (real frame)

"G_loc = G: no" means the window needs the rank test (group C). The 4-up chart Z is in the mirror frame where marked.

| t' → t (ρ) | chosen W'_t' | chosen W'_t | step unknowns | G_loc = G | 4-up chart Z (frame) |
|---|---|---|---|---|---|
| 2,0,1/2,0,1/1,1,1 → 1,1,1/1,1,1/1,1,1 (F0.2) | F0.1,F1.2,E1.1,E0.3,E1.1 | F0.1,F1.2,E1.1,E0.2,F0.1,E1.1 | 34 696 | yes | F0.1,F1.2,F2.3 |
| 1,1,1/0,2,1/0,2,1 → 1,1,1/1,1,1/1,1,1 (E0.2) | F0.1,F1.2,F0.1,E1.2 | F0.1,F1.2,F0.1,E0.1,E1.2,E0.1 | 295 452 | yes | F0.1,F1.2,F2.3 |
| 0,2,1/2,0,1/1,2,0 → 0,1,2/2,0,1/1,2,0 (F1.1) | E0.2,E1.1,E0.1,F1.1,E1.2,F0.3,F1.1 | E0.2,E1.1,E0.1,F1.1,E1.2,F0.3,F1.2 | 752 | yes | F1.2,F0.1,F1.1,F2.3 (mirror) |
| 0,2,1/1,0,2/2,1,0 → 0,3,0/1,0,2/2,1,0 (E1.1) | E0.1,E1.2,E0.2,F1.2,E1.1,F0.3,F1.1 | E0.1,E1.2,E0.2,F1.2,E1.1,F0.3 | 58 | no (rank test) | F1.1,F0.2,F1.2,F2.3 (mirror) |
| 0,1,2/2,0,1/1,2,0 → 0,2,1/2,0,1/1,2,0 (E1.1) | E0.2,E1.1,E0.1,F1.1,F0.3,F1.2,E1.2 | E0.2,E1.1,E0.1,F1.1,F0.3,E1.2,F1.1 | 6 073 | yes | F1.2,F0.1,F1.1,F2.3 (mirror) |
| 1,1,1/1,1,1/1,1,1 → 2,0,1/2,0,1/1,1,1 (E0.2) | F0.1,F1.2,E0.1,F0.1,E1.2,E0.1 | F0.1,F1.2,E0.2,E1.2,E0.1 | 318 975 | yes | F0.1,F1.2,F2.3 |
| 1,0,2/2,1,0/1,1,1 → 1,0,2/1,2,0/1,1,1 (F0.1) | F0.2,E1.1,F1.2,E0.3,E1.1 | F0.2,E1.1,F1.2,E0.2,E1.1 | 27 192 | yes | F0.2,F1.2,F2.3 |
| 1,2,0/1,1,1/1,2,0 → 1,1,1/1,0,2/1,2,0 (F1.2) | E1.1,E0.2,F1.2,F0.2,E1.3 | E1.1,E0.2,F1.2,F0.2,F1.1,E1.2 | 6 239 | yes | F0.1,F1.2,F2.3 (mirror) |
| 1,0,2/2,0,1/1,2,0 → 1,0,2/1,1,1/1,2,0 (F0.1) | F0.2,E1.1,F1.3,E0.3,E1.2 | F0.2,E1.1,F1.3,E0.2,E1.2 | 18 461 | yes | F0.2,F1.3,F2.3 |
| 0,2,1/2,0,1/1,1,1 → 0,3,0/2,0,1/1,2,0 (E1.2) | E0.2,E1.1,E0.1,F1.1,F0.3,F1.1,E1.1 | E0.2,E1.1,E0.1,F1.1,F0.3,E1.2 | 61 | no (rank test) | F1.2,F0.1,F1.1,F2.3 (mirror) |
| 0,2,1/2,0,1/1,1,1 → 0,3,0/2,0,1/1,1,1 (E1.1) | E0.2,E1.1,E0.1,F1.1,E1.1,F0.3,F1.1 | E0.2,E1.1,E0.1,F1.1,E1.1,F0.3 | 58 | no (rank test) | F1.2,F0.1,F1.1,F2.3 (mirror) |
| 1,1,1/2,0,1/1,1,1 → 1,1,1/1,1,1/1,1,1 (F0.1) | F0.1,F1.2,E1.1,E0.3,F0.1,E1.1 | F0.1,F1.2,E1.1,E0.2,F0.1,E1.1 | 51 816 | yes | F0.1,F1.2,F2.3 |
| 1,0,2/3,0,0/0,2,1 → 1,0,2/2,1,0/0,2,1 (F0.1) | F0.2,F1.2,E0.3,E1.2 | F0.2,F1.2,E0.2,E1.2 | 18 461 | yes | F0.2,F1.2,F2.3 |
| 1,2,0/1,1,1/1,1,1 → 1,2,0/1,0,2/1,1,1 (F1.1) | E0.1,E1.2,E0.1,F1.1,F0.2 | E0.1,E1.2,E0.1,F1.2,F0.2 | 113 288 | yes | F0.1,F1.2,F2.3 (mirror) |
| 2,0,1/1,1,1/1,1,1 → 2,0,1/0,2,1/1,1,1 (F0.1) | F0.1,E1.1,F1.2,E0.2,E1.1 | F0.1,E1.1,F1.2,E0.1,E1.1 | 13 051 | yes | F0.1,F1.2,F2.3 |
| 1,1,1/0,2,1/0,2,1 → 2,0,1/1,1,1/0,2,1 (E0.2) | F0.1,F1.2,F0.1,E1.2 | F0.1,F1.2,E0.1,E1.2 | 389 029 | yes | F0.1,F1.2,F2.3 |
| 1,2,0/2,1,0/1,1,1 → 1,1,1/2,0,1/1,1,1 (F1.2) | E1.1,E0.3,E1.1,F0.2 | E1.1,E0.3,E1.1,F1.1,F0.2,F1.1 | 194 000 | yes | F0.1,F1.3,F2.4 (mirror) |
| 3,0,0/1,0,2/1,1,1 → 2,1,0/0,1,2/1,1,1 (F0.2) | E1.1,F1.2,E0.2,E1.1 | E1.1,F1.2,E0.1,F0.1,E1.1 | 29 841 | yes | F1.2,F2.3 |
| 1,2,0/0,1,2/0,3,0 → 1,2,0/0,1,2/2,1,0 (E0.2) | E1.3,F1.2,F0.2 | E1.3,F1.2,F0.2,E0.2 | 0 | Id | F1.2,F2.2 |
| 1,2,0/0,1,2/2,0,1 → 1,2,0/0,2,1/2,0,1 (E1.1) | E1.2,E0.2,F1.2,F0.2 | E1.2,E0.2,F1.1,F0.2 | 57 736 | yes | F0.2,F1.2,F2.3 (mirror) |
| 1,1,1/1,2,0/1,0,2 → 2,0,1/2,1,0/1,0,2 (E0.2) | F0.1,F1.1,F0.1,E1.1,E0.2 | F0.1,F1.1,E1.1,E0.3 | 51 816 | yes | F0.1,F1.1,F2.3 |
| 3,0,0/0,0,3/3,0,0 → 1,2,0/0,0,3/2,1,0 (F0.3) | E1.3,F1.3,E0.3 | E1.3,F1.3,E0.2,F0.2 | 4 603 | yes | E1.3,F1.3 |
| 3,0,0/0,0,3/3,0,0 → 2,1,0/0,0,3/1,2,0 (F0.3) | E1.3,F1.3,E0.3 | E1.3,F1.3,E0.1,F0.1 | 3 310 | yes | E1.3,F1.3 |
| 0,2,1/2,0,1/1,1,1 → 0,2,1/2,0,1/1,2,0 (E1.1) | E0.2,E1.1,E0.1,F1.1,F0.3,F1.1,E1.1 | E0.2,E1.1,E0.1,F1.1,F0.3,F1.1,E1.2 | 238 | yes | F1.2,F0.1,F1.1,F2.3 (mirror) |
| 1,2,0/2,0,1/1,0,2 → 1,2,0/2,0,1/0,1,2 (F0.1) | E0.2,E1.3,F1.3,E0.1,F0.2 | E0.2,E1.3,F1.3,F0.2 | 3 002 | yes | F1.1,F2.3 |
| 2,0,1/3,0,0/0,1,2 → 1,1,1/2,1,0/0,1,2 (F0.2) | F0.1,F1.1,E0.3,E1.1 | F0.1,F1.1,E0.2,F0.1,E1.1 | 34 696 | yes | F0.1,F1.1,F2.3 |
| 1,2,0/1,2,0/0,2,1 → 1,2,0/1,2,0/1,1,1 (E0.1) | F0.2,E0.1,E1.2 | F0.2,E0.1,E1.2,E0.1 | 0 | Id | F2.1 |
| 2,1,0/0,1,2/1,1,1 → 2,1,0/0,3,0/1,1,1 (E1.2) | E1.1,E0.1,E1.1,F1.2,F0.1 | E1.1,E0.1,F0.1,E1.1 | 19 371 | yes | F0.1,F1.1,F2.3 (mirror) |
| 3,0,0/2,0,1/1,0,2 → 2,1,0/1,1,1/0,1,2 (F0.3) | F1.1,E1.1,E0.3 | F1.1,E0.1,F0.1,E1.1 | 210 002 | yes | F1.1,F2.3 |
| 1,2,0/2,0,1/0,3,0 → 1,1,1/2,0,1/0,2,1 (F1.2) | E0.2,F1.1,F0.2,E1.3 | E0.2,F1.1,F0.2,F1.1,E1.2 | 3 583 | no (rank test) | F1.2,F2.3 (mirror) |
| 1,1,1/2,1,0/0,1,2 → 2,0,1/3,0,0/0,1,2 (E0.2) | F0.1,F1.1,E0.2,F0.1,E1.1 | F0.1,F1.1,E0.3,E1.1 | 13 051 | yes | F0.1,F1.1,F2.3 |
| 1,2,0/0,3,0/0,1,2 → 2,1,0/1,2,0/1,0,2 (E0.3) | E1.1,F0.2 | E1.1,E0.2,F0.1 | 27 286 | yes | F2.2 |
| 3,0,0/2,0,1/0,2,1 → 2,1,0/1,1,1/0,2,1 (F0.2) | F1.1,E0.2,E1.2 | F1.1,E0.1,F0.1,E1.2 | 78 371 | yes | F1.1,F2.2 |
| 3,0,0/1,0,2/2,0,1 → 1,2,0/1,0,2/1,1,1 (F0.3) | E0.1,E1.3,F1.3,E0.2 | E0.1,E1.3,F1.3,E0.1,F0.2 | 3 471 | yes | F1.2,F2.3 |
| 1,0,2/1,2,0/1,1,1 → 1,0,2/2,1,0/1,1,1 (E0.1) | E0.1,E1.2,E0.1,F0.2,F1.2 | E0.2,E1.2,F0.2,F1.2,E0.1 | 2 749 | no (rank test) | F0.2,F1.2,F2.3 |
| 2,1,0/0,1,2/0,3,0 → 2,1,0/0,1,2/1,2,0 (E0.1) | E1.3,F1.2,F0.1 | E1.3,F1.2,F0.1,E0.1 | 0 | Id | F1.2,F2.2 |
| 2,0,1/0,2,1/1,1,1 → 2,0,1/0,3,0/1,1,1 (E1.1) | F0.1,E1.1,E0.1,F1.2,E1.1 | F0.1,E1.1,E0.1,E1.1,F1.1 | 4 588 | yes | F0.1,F1.1,F2.3 (mirror) |
| 2,0,1/0,2,1/1,1,1 → 2,0,1/0,2,1/2,0,1 (E0.1) | F0.1,E1.2,F1.2,E0.1 | F0.1,E1.2,F1.2,E0.2 | 238 | yes | F0.1,F1.2,F2.3 |
| 1,2,0/2,1,0/0,3,0 → 1,2,0/2,0,1/0,1,2 (F1.3) | E0.2,F1.1,F0.2,E1.4 | E0.2,F1.1,F0.2,E1.1 | 1 516 | yes | F1.2,F2.3 (mirror) |
| 1,1,1/1,1,1/1,2,0 → 1,1,1/1,0,2/1,2,0 (F1.1) | E1.1,E0.2,F1.2,F0.2,E1.3,F1.1 | E1.1,E0.2,F1.2,F0.2,F1.1,E1.2 | 9 822 | yes | F0.1,F1.2,F2.3 (mirror) |
| 2,1,0/0,2,1/1,2,0 → 2,1,0/0,0,3/1,2,0 (F1.2) | E1.3,E0.1,F1.5,F0.1,E1.4 | E1.3,E0.1,F1.5,F0.1,E1.2 | 4 668 | yes | F0.1,F1.1,F2.3 (mirror) |
| 1,2,0/2,0,1/0,2,1 → 1,2,0/2,0,1/0,1,2 (F1.1) | E0.2,E1.3,F1.2,F0.2 | E0.2,F1.1,F0.2,E1.1 | 283 842 | yes | F1.2,F2.3 (mirror) |

## 8. Reproduction (scratchpad tools; Python 3, no cargo)

`TOOLS=/private/tmp/claude-501/-Users-existair-Documents-KRcalculator/bf6f8870-08ad-436b-9ff7-0414655cc59a/scratchpad/colored_design/b5/tools`

| step | command | output |
|---|---|---|
| 1 | `python3 $TOOLS/analyze.py` | parse, validate `step_score` (0 mismatches), `grass`, local K-windows → `res.pkl` |
| 2 | `python3 $TOOLS/zpre.py; python3 $TOOLS/chart2.py` | 4-upright charts and detour bounds → `chart2.pkl` |
| 3 | `python3 $TOOLS/joint_all.py`, repeated until "done 42/42" | 3-upright joint search, cached (`wc.pkl`, `joint.pkl`) |
| 4 | `seeded_all.py`, `seeded_ok.py`, `seeded_hard.py` | seeded search for 6-7-rung types; `jfilter.py` (G_loc filter) |
| 5 | `consolidate.py; detours.py; groups.py; targets.py; mdtable.py` | tables (`steps2.json`, `targets.json`, `table.md`, `appendix.md`) |
