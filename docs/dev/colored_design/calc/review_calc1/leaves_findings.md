AREA REVIEW: B4w leaf-pair wiring + C2-7 joint words / per-predecessor finals / Prop W (read-only; worktree /Users/existair/Documents/KRcalculator/.claude/worktrees/calc at 1a9c8de)

VERDICT: I found no soundness defect, data-loss path or crash in this area. Joint words, per-predecessor final words and e_q Leaves pairs never admit a map on their own. Every result still goes through certify_pairs and then sindec's full checks. One low-severity defect, in the diagnostics and a test, is listed below.

== Finding 1 (low; diagnostic / test coverage) ==
- Where: src/bin/gp0.rs:482-488, tests/c27_joint.rs:22, together with src/leaves.rs:1194 and src/leaves.rs:1210.
- What is wrong: `final_pred_words` accepts a per-predecessor final word W' as soon as its W-aware score sc1 is lower than s0 (`if sc1 < s0`). It does not require sc1 ≤ the 1e5 budget. Both the census and the test then treat the replaced final step as solved:
  - gp0 decrements `over` without condition and scores the step as 0 in `mx` (gp0.rs:484, 488).
  - `plan_max` in the test filters the step out (c27_joint.rs:22).
  - Neither looks at the replacement step's own score `sw`. gp0 prints `sw` but does not use it.
- Scenario: the original final step scores s0 = 3.6 M, the best W' scores 2 M, and a rex chain W' ⇄ W_q exists.
  - The census prints "steps > 1e5: 0" and a max that comes only from the other steps.
  - `per_predecessor_final_pair_1444` would still pass.
  - At run time, `build_leaves` (leaves.rs:1258) makes the step into W' an exact call of 2 M unknowns. `admit` (leaves.rs:786) refuses it because it is above the 4e5 hard cap, so the job is Refused.
- Impact: the planner is misreported as feasible, and the test cannot catch a regression. Results are never wrong.
- Published numbers: the C2-7 report §2 column "joint + per-predecessor finals" comes from this metric. The actual alt final windows for 1444/1470 were 10 530 unknowns, so the published figures hold, but nothing checks them.
- Fix: in gp0, put `sw` into `mx`/`over` instead of 0 / the decrement, and assert `sw ≤ 1e5` in the test.

== Answers to the review questions ==

1. Is Lemma A certification still sound with joint words?
- Yes. `joint_descent` (leaves.rs:1051-1167) and `final_pred_words` (1169-1215) are pure planning. They only pick words that pass `reduced_for` (`suffix_words`, 469-507), and `build_leaves` re-checks each one (1243-1245).
- Every step family still has to pass `step_family`'s test (809-952): degree multiset exactly G, and full rank on the top-to-top columns at rep(t).
- Each level still gets `build_leaves`' square/invertible diagnostic at rep(t) (1295-1306).
- Final acceptance is unchanged:
  - certify_pairs (homprov.rs:790-814): form rank = mult, a second-point equality check, genuineness;
  - sindec `dechom_postcheck` (sindec.rs:378);
  - `check_pi_iota` (sindec.rs:1411): Π∘Ι = δ with e_q on the diagonal, and Σιπ = e_X when nothing is killed.
- Per-predecessor finals are composed with realized degree-0 rex chains (Out: R∘S, In: S∘R; leaves.rs:1254-1281). Every leaf still lands in W_q, and the rank diagnostic plus certify_pairs re-verify the result.
- The predecessor index fi lines up with plan[m−1]: both are built in js[k] order (plan_words 965-972, build_leaves 1242/1310).

2. Prop W: is the window family checked complete before an exact call, and what happens when a check fails?
- W1/W2 run in planning only: `window_w12` (525-541) and `step_score_w` (543-564).
- W3/W4 are not implemented anywhere.
- `step_family` never consults W. It always tries the window family first, then the whole words.
- This departs from D_grass §1.3, which wants W to flag a bad window "before any exact call" and to use a dry rank test with Whole disabled. The cost is time only, because every exact call is still capped:
  - window and whole calls both go through `admit` (786): ≤ 1e5 unknowns, or ≤ 4e5 with an uncensored cost record;
  - a refused whole call returns Refused (908-910);
  - a rank shortfall returns Exhausted (952).
- A failed check therefore always ends in refusal, never silent acceptance. This matches D_grass §4: "Prop W ... never admits a map".

3. e_q pairs (c27):
- The e_q contract (`check_idempotent`) and mult are enforced by `dechom_precheck` (sindec.rs:346-375) before any provider runs.
- After wrapping (ι∘e_q, e_q∘π; leaves.rs:1752-1755), π'∘ι' lies in e_q·End⁰(L_q)·e_q, which is one-dimensional with top scalar 1 at rep(q). So the form check at rep(q) is enough.
- The rex shortcut is excluded when e_q is present (leaves.rs:1497).
- Maps that use Elem constants layer A cannot store are refused before commit (1447-1476, 1514-1517).

4. Provider stamps (KRHOM_LEAVES_C27):
- The provider flags and the stamp token come from the same env read (leaves.rs:1411-1418), so a provider cannot carry the wrong stamp.
- Without `--stamp-route`, the full header is compared exactly, so a c27 file is rejected by a non-c27 binary and vice versa (sindec.rs:1749-1757). Mismatches are detected.
- With `--stamp-route` (as in the C2 runs and the calculator's run_prime, `layer_a_route_stamp: true`, calc.rs:398-401), `route_header` drops the `c27=` token (sindec.rs:458-468). Cross-loading in either direction is then silent by design; it only shows up in the `layer_a_route` counter.
- That is safe:
  - layer-A keys do not depend on the provider (sindec.rs:1700-1738);
  - loaded pairs go through `dechom_postcheck` (1241) and then `check_pi_iota`;
  - closures go through `certify_closure` (1655), with loaded generators re-wrapped with the current e_p;
  - smapio re-certifies Added templates and theorem-level Assoc nodes on import.
- If an embedded e_q differs from the session's e_q, the result is a hard DecErr::Check, not acceptance.
- The calculator refuses the env var outright (calc.rs:651-657).
- `b4wstamp` only restamps files whose route header equals the chain stamp's. Everything else is skipped and reported, and the write is atomic (tmp + rename).

5. Default-off guarantee (c27 unset):
- `joint_descent` runs only when `joint` is set (1033), and `final_pred_words` returns an empty map (1172), so `alt` is None and steps are built exactly as before.
- e_q targets are returned before `leaf_pairs` runs (1556) and are refused inside it too (1493).
- The stamp is byte-identical to B4w.
- The wrappers `pairs_leaves`, `rex_class_pairs` and `realize_chain` pass eq = None / persistable = false, so they are unchanged.
- The only visible difference from B4w for e_q targets is the failure text: the "[Leaves] not implemented" entry is gone. The engineering classification is unchanged.

6. Other changes in the area that alter behaviour by default:
- `interp_gate` (leafgen.rs:80-104, called at 358, 697, 733-741, 947; commit 1879eba) is on by default but sits only on the Leaves/leafgen paths. No production-chain code calls `interpolate`/`build_in`.
  - It can turn a Leaves interpolation that used to succeed into a refusal.
  - It cannot reject persisted templates: `certify_added`/`poly_certificate` on import are not gated.
- Outside my area, but relevant to "default-off": the C2-7 column-restricted localization (smap.rs window_block with LOC_COLS_MIN = 48, tmpl.rs `localize_cols`, local.rs `inverse_cols`) is ON by default; `KRHOM_LOC_FULL=1` turns it off. So C2-7 §3's "byte-identical when off" claim does not cover it. It is exact arithmetic, and one values-hash comparison is reported. Whoever reviews smap/tmpl/local should check it.

7. Non-defects noted:
- The Leaves pair route ignores the `Budget` time limit, and there is no time cap on the joint planning loops. Their sizes are bounded by `joint_cap` / `joint_keep` / planner `max_visit`.
- `realize_chain_opts` caches persistable-Assoc Pull templates under the FlagMatch key. That is consistent, because FlagMatch builds the same CfKind::Pull closed form (leafgen.rs:1281-1298 vs 1347).
- Templates are rebuilt per `build_leaves` call (separate caches for the Out and In builds). This duplicates templates and costs time, but is otherwise harmless.