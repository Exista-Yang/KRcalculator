VERIFICATION (area: leaves). Worktree /Users/existair/Documents/KRcalculator/.claude/worktrees/calc at 1a9c8de. I only read files; nothing was run.

== Finding 1: per-predecessor final words accepted above the budget, and the census and test hide it. Verdict: REAL. Severity: minor (diagnostics and test coverage only; results are never wrong) ==

Confirmed by reading the code:
- src/leaves.rs:1185-1210 (`final_pred_words`):
  - The function only starts looking when s0 > self.budget (1185).
  - It keeps any candidate with `sc1 < s0` (1194). There is no check that `sc1 <= self.budget`.
  - It then takes the first of the 8 cheapest candidates that has a degree-0 rex chain both ways (`out.insert((q, i), c)`, 1210).
  - So if no candidate at or under 1e5 has a chain, the function returns a word above budget, for example sc1 = 2 M when s0 = 3.6 M.
- src/bin/gp0.rs:472-488 (`steps_census`, `--steps`):
  - For every replaced final step whose original `score_w` is above 1e5, the code runs `over -= 1` unconditionally (480-485).
  - `mx` scores those steps as 0 (488).
  - The replacement score `sw` (476) is only printed.
  - `over` cannot underflow: plan_steps (leaves.rs:592-613) emits one step per (level, tp, t), and each decremented step was counted earlier (gp0.rs:430-432).
- tests/c27_joint.rs:22 (`plan_max`): the replaced final steps are filtered out completely, so `per_predecessor_final_pair_1444` (40-45) checks only `nf >= 1` and the other steps. A regression where W'_{q,i} scores above 1e5 would still pass.

What happens at run time:
- build_leaves (leaves.rs:1256-1258) calls `step_family` on W'_{q,i}.
- Every exact call goes through `exact_local` -> `admit` (leaves.rs:786 -> 751-779). A call is admitted at ≤ 1e5 unknowns. Between 1e5 and 4e5 it needs an uncensored cost record. Above 4e5 it is refused.
- So a 2 M step produces MissKind::Refused. It can never produce a wrong result: certify_pairs, dechom_postcheck and check_pi_iota are still downstream.
- One caveat on the stated scenario: `step_score_w` uses whole-word cost whenever W1/W2 fails (leaves.rs:543-564), and step_family still tries the window first. So 2 M is an upper estimate, and the real call could sometimes be cheaper and succeed. The census defect stands either way: it reports "steps > 1e5: 0" without ever checking `sw`.

Nothing else in the code reads these numbers. plan_steps and final_pred_words are used outside leaves.rs only by gp0.rs:396/426/472 and tests/c27_joint.rs. No production route uses the census to choose a path, so the impact is a misreported feasibility figure in C2-7 §2 and a test that cannot catch this regression.

Suggested fix: in gp0, count `sw` in `mx` and `over` instead of 0 and the decrement. In `plan_max`, take the max with `step_score_w` of the replacement word. Optionally, require `sc1 <= self.budget` (or ≤ hard_cap) at leaves.rs:1194, or prefer such candidates.

== Other claims in the review that are not findings (checked only briefly) ==
- Item 6: the column-restricted localization is ON by default, confirmed at src/smap.rs:1215 (`src_fixed_points(t) >= LOC_COLS_MIN (48) && KRHOM_LOC_FULL unset`). So the C2-7 claim that the output is "byte-identical when off" does not cover this path. That is a scope note for whoever reviews smap/tmpl/local. It is not a defect with a failure scenario in the leaves area.
- I did not re-check the review's soundness answers 1-5 line by line. The code I read for this finding (build_leaves 1222-1312, admit 751-779, step_score_w and window_w12 525-564) agrees with them: joint words and per-predecessor final words are planning only, the W'_{q,i} leaves are composed with a realized rex chain back to W_q (1259-1277), and build_leaves still checks at rep(t) that the leaf block is square and invertible (1295-1306).

Summary: 1 finding. REAL: 1 (minor). NOT-REAL: 0. UNCERTAIN: 0. No blocker in the leaves area.