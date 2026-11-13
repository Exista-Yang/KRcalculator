# Review wave 5: theorem-level Assoc certificate (revrex, m4c-rexplan 4e886f6)

Scope: src/leafgen.rs (lines 703-1281), src/smapio.rs (export 288-300, import 439-459), tests/thmassoc.rs,
plus the callers: rexplan.rs (assoc_move, needs_theorem, uncertifiable) and leaves.rs (rex_class_pairs, pairs_leaves).
I read the code only. I ran no cargo, binaries or computations.

## Verdict

No critical or major defect. I found nothing that can produce a wrong Poincare polynomial, a false certification, silent
corruption or a crash. There are two minor findings, both test gaps (section 6).

## 1. Acceptance condition: matches the approval, nothing broader

`theorem_assoc_check` (leafgen.rs:767-790) checks:
- equal bottoms and equal tops;
- equal `multistep` expressions;
- `is_updown`;
- both words flag-type (`is_flag_type`, at most 1 flow per top coset);
- equal `support()`. The coset ranks share one CosetSpace because bottoms and tops are equal, so this gives a flow
  bijection with equal top configurations;
- equal `shift()`.

`is_updown` (758-764) works on wall sets. More walls means a finer composition, which means a smaller parabolic.
- With 3 terms, the middle wall set must be contained in both ends. So the middle is the coarser parabolic W ⊃ I, J,
  which is [[I ⊂ W ⊃ J]] (merge, then split).
- A down-up shape has its middle ⊇ both ends and is refused. Check: a single rung F0.1 on (2,1) gives
  [{2},{1,2},{1}] and is refused. canonical_max gives [{2},{},{1}] and is accepted. The (3,3,3) log shows
  [[3,6],[],[3,6]].
- With ≤ 2 terms the shape is monotone: [I ⊂ J] or [I ⊃ J]. These are the degenerate up-down cases (W = J or W = I).
  The reduction keeps adjacent terms comparable (transitivity, and simultaneous removals in one pass are consistent
  after dedup), so a 2-term result is always monotone.

`multistep` (575-617) is a valid contraction of monotone runs:
- the raw sequence is L_t ⊆ R_{t+1} ⊇ L_{t+1} in walls;
- the refinement puts the moved part next to the receiver for both E and F rungs;
- zero labels are absorbed by `walls`.

The approval does not ask for a flow cap, and the code has none. realize_chain (1210-1221) routes ≤ 500-flow windows to
`validate_assoc`.

Only `certify_assoc_theorem` (795-800) registers a pair, and it always runs the check first.

`add_thm_assoc` (806-817) skips the theorem in one case only: the arena already holds a constant for the pair and the
pair is not registered. That constant can only be numeric.

The planner's `uncertifiable` (rexplan.rs:248) is advisory only. realize_chain re-runs the full theorem check on the
arena's own trimmed window (`elem_windows`). assoc_move strips the window exactly as `smap::elem_spec` does
(at = a0, first = lo).

## 2. Localized map: identity matching, scalar 1, no sign (correct)

Two facts:
- In krhom's convention 1⊗1 localizes to the all-ones vector: column 0 of L is the weight-0 fibre monomial
  (symfun.rs:468 sorts by weight).
- A bimodule map preserves top cosets.

So any bimodule map between two flag-type words with equal support is diagonal, with one entry d_c per coset. Requiring
1⊗1 ↦ 1⊗1 forces d_c = 1 for every c.

So the SLL (88) isomorphism (up-up/down-down contraction, 1⊗1 ↦ 1⊗1) must localize to exactly `elem_winloc(Assoc)`
(smap.rs:476-482): match by equal top configuration, value 1, no sign. Eval applies the window map with equal
configurations outside the window (smap.rs:1189-1250, unchanged). The node scalar is F::one() (leafgen.rs:807).

An edited scalar in a layer-A file still gives a nonzero multiple of a genuine map. That is not a false certification.

## 3. Cross-checks: what they can and cannot detect

**(i) `theorem_vs_numeric` (872-889).** It runs the same `poly_certificate` that `validate_assoc` runs, on the same pair
and the same closed form. The small-window test asserts `thm only == 0`, so every theorem acceptance was checked
numerically. Braid middles have no common prefix or suffix (M starts E0.y and ends F0.y; M' starts E1.z and ends
F1.x), so the pair cross-checked is the pair `braid_dag_thm` certifies.

Limitation: the report records **no numeric failure at all**. There are 4 866 theorem refusals, and every one passes
numerically. So (i) has no negative instances. It would also pass for a check that ignored the multistep and up-down
conditions. It confirms that identity matchings are polynomial, but it does not test the acceptance condition.

**(ii) `spot_check` (1028-1140).** The ring and degree logic is correct:
- polynomial in y along a random line, with degree ≤ D/2;
- same u and D as `poly_certificate`;
- LightLoc fdeg and mixed-radix order are consistent (and tested against `loc_matrix`);
- Newton interpolation and Horner evaluation are correct;
- `m + 1 ≤ npts` always holds.

Equivariance under S_λ makes polynomial-in-y equivalent to membership in R^λ for these closed forms.

For the theorem Assoc, however, two parts hold by construction:
- part (a) checks only that cosets are preserved, at 16 cosets drawn with replacement;
- column 0: a 0/1 matching sends ones to ones, so Φ[:,0] = e_0 identically.

So at (3,3,3) all the information comes from the 5 sampled columns [29, 201, 2, 43, 10]. No test shows that this part
rejects anything (finding 2).

## 4. Marking and import re-check

- `export` (smapio.rs:290-300) writes `thm-assoc k` for every exported Assoc node whose window pair is registered.
- `import` (439-459) re-runs `certify_assoc_theorem` on exactly those nodes. It refuses a record on a node that is not
  an Assoc, and refuses a pair that fails the check.
- Without the record, a fresh arena has no constant, and `assert_genuine` fails (smap.rs:853-880).

A record-less load is genuine only if the arena already holds the constant from an earlier certification (theorem or
numeric) in the same arena. That is legitimate. Nothing persists Assoc constants outside smapio; I checked revm2a's
libbuild, boxclose, boxrun and swire as well.

Over-marking (no correctness effect):
- The registry is process-wide and keyed only by window pair, so `relies_on_theorem` and `export` can mark a node whose
  constant in that arena is numeric. The importer then re-certifies it by the theorem, which still passes.
- The provenance tag is set even on the numeric early return (807, 812). `braid_dag_thm` tags the composite whether or
  not the theorem was used (1171).
- After a smapio round trip, provenance reads "smapio import" and the tag is lost. Marking survives through the record
  and the registry. No code reads provenance.

## 5. Theorem use outside the approved rules?

No path accepts a theorem Assoc outside the combinatorial rules.

One policy observation, which is not a defect: `add_thm_assoc` and `braid_dag_thm` apply the theorem to windows the
numeric certificate could handle (≤ 500 flows; for example `braid_dag_thm` at (1,2,1) in the marking test) without
running it. realize_chain does prefer numeric there. Cross-check (i) covers those braid cells offline.

## 6. Findings

### F1 (minor, test does not test what it claims): import refusal of a forged record on a non-theorem pair is untested
tests/thmassoc.rs:237-246. The comment says such a pair "is not expressible (add refuses)". That is false:
- `Arena::add` refuses an Assoc only when there is no bijection or the shifts differ (smap.rs:701-712);
- the B2 report itself counts 4 866 bijective pairs, with equal shifts, that the theorem refuses (down-up, or a
  different multistep).

The substitute check moves the record to node 0. The test then asserts genuineness only if the import succeeds. It
never asserts `r4.is_err()`.

The comment at line 232 is also wrong: it uses the same registered pair. The assertion there still holds, because
genuineness is per arena.

The importer code is correct as written. This is a test gap: it would not catch a regression in the most important
negative case of "re-checked on import".

Suggested fix: numerically validate a bijective, equal-shift, down-up Assoc (≤ 500 flows), export it, inject
`thm-assoc k`, import into a fresh arena, and assert the "thm-assoc re-check failed" error. Also assert `r4.is_err()`
when node 0 is not an Assoc.

### F2 (minor, weak negative control): the spot-check corruption test only exercises column 0
tests/thmassoc.rs:151-177. Scaling one coset by 2 is always caught deterministically by column 0:
- Φ[:,0] = e_0 + L_tgt⁻¹ e_c;
- every row bp ≥ 1 has fdeg ≥ 2, so its degree D is negative and the entry must be zero;
- L_tgt⁻¹ e_c ∝ e_0 is impossible when there is more than one coset.

Yet the test only asserts that at least half of the corruptions are detected. More importantly, no corruption keeps
f(1⊗1) = 1⊗1 and coset preservation while breaking polynomiality in columns ≥ 1. Those columns are the only part of the
spot check that tells us anything about the theorem Assoc (section 3).

Effect: a regression in the columns ≥ 1 branch would make the (3,3,3) spot check vacuous, and this test would still
pass. Examples: a wrong D, LightLoc column order wrong for a shape the 3 LightLoc test words do not cover, or a skip
condition that discards entries.

Suggested fix:
- assert that every corruption is detected (`bad == support.len() - 1`);
- add a corruption g with g·ones = 0 inside a coset block with more than one flow (for example on a Rex window,
  +δ / −δ on two target flows of the same coset), check with `cols` excluding 0, and assert Err.
