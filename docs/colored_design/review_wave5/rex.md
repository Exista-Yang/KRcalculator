# Review wave 5: rex planner and rex-class pairs (revrex, m4c-rexplan 4e886f6)

Lens: src/rexplan.rs, leafgen::realize_chain, leaves::rex_class_pairs / pairs_leaves, src/bin/gp0.rs,
symfun.rs GrTable persistence (+ grwarm), tests/rexplan.rs. Read-only; nothing compiled or run.
Line numbers refer to the revrex worktree.

## Verdict

No defect in this lens can produce a wrong Poincare polynomial or a false certification. Every
realized chain is a DAG of certified leaves, and the pair certificates (certify_pairs, plus pi.iota = id
for rex-class pairs) do not depend on the planner being right: a planner or realization mistake
can only lead to an Err, never to a wrong certified map. I found three minor issues: one latent
stale-cache risk, one gap between a code comment and the approval, and one possible
multi-hour stall that the current N = 5/6 lists do not reach.

## What was verified (claims vs code)

- **Templates (Collapse/Expand/FlagMatch/Swap)**: realize_chain builds the ClosedForm (Rex / Pull / Pull / Swap{na:1}),
  refuses cf.deg != 0 (leafgen.rs:1264), certifies with build_in -> interpolate (fit checked at extra points)
  -> certify_added (degree check, commutation with top operators, poly_certificate) (leafgen.rs:1267, tmpl.rs:410-440).
  Rex normalizes with sum 1/e_Y and Pull/FlagMatch are identity matchings, so the all-ones vector maps to the
  all-ones vector (1(x) -> 1(x)) by construction. ClosedForm::new enforces equal shifts for Rex, a flag-type source
  for Pull (canonical_max and both FlagMatch ends are flag-type), and disjoint uprights for Swap. Reversed moves
  from the backward BFS half (rexplan.rs:676-704) keep every precondition (Collapse<->Expand swap; Assoc/FlagMatch/Swap symmetric).
- **Assoc**: Elem node; the arena window from elem_spec (prefix from `at`, common suffix stripped, smap.rs:301-328)
  is the same as the planner's assoc_move window (a0 = first differing rung, rexplan.rs:427-460). The node is
  certified by validate_assoc when <= 500 flows, otherwise by add_thm_assoc, which refuses non-up-down shapes.
- **Exact**: try_hom_basis_cap(ws, wt, 0, ...) (leafgen.rs:1228). The first basis element whose row sums are
  c * ones with c != 0 is chosen and scaled by c^-1 through Lin (leafgen.rs:1246). A single evaluation
  point is enough: both words are reduced for the same window type, so the shifts are equal and u = 0. Column 0 of Phi
  (fibre degree 0) can then only hit target fibre degree 0, so Phi e_0 = c e_0 exactly.
- **Exact budget direction**: neighbours checks hom_unknowns(src -> z) (rexplan.rs:579-580). Moves found from
  the `to` side are used in reverse, but at d = 0 with equal shifts (equal top fibre degree) and palindromic
  fibre-degree histograms, hom_unknowns is symmetric. No hidden over-budget reverse call.
- **rex_class_pairs** (leaves.rs:1123-1170): it checks that X is reduced for q (so B_q occurs once, at s = 0: sdeg 0 and mult 1
  are correct), plans X->L_q and L_q->X (plan() rejects ends that are not reduced for q), and calls
  certify_pairs(x, lq, 0, [io], [pi], rep_rank_of(q), 1) (check_info ends/degrees, form rank, equal form at a 2nd
  point, assert_genuine). pi.iota = id is then checked blockwise on every coset of layout(lq, lq), at 2 points
  (leaves.rs:1139-1153). Composition order is right: Comp(pi, io) = pi o io : L_q -> L_q.
- **pairs_leaves**: exit o Out-leaf and In-leaf o entry (leaves.rs:1194, 1199). Arena::add refuses mismatched ends and
  certify_pairs re-checks the ends, so a wrong orientation would fail and could not certify anything.
- **B4 logs** match the report table: 7 rex-class targets, both primes, form [[1]], theorem-Assoc nodes [].
- **GrTable load** (symfun.rs:398-424): magic, version, (r, b), length, FNV-1a hash of the payload, a structural parse
  (bi < nb, row counts, exact end of buffer, e_mul.len() == r, f_mul.len() == b - r), and basis/weights equal to
  fibre_basis(r, b). Writes are atomic (pid temp file + rename). Bit flips, truncation, a wrong (r, b) and a
  layout change (with a version bump) are all rejected.
- **tests/rexplan.rs**: realized_chains_small really checks genuineness, degree 0, ones -> ones at every coset,
  scalar 1 at rep(q), and membership in the exact Hom^0 span at 2 points. rex_class_pairs_small asserts form [[1]]
  and no theorem nodes. No vacuous assertion found.

## Findings

### m1 (minor): GrTable cache content is not bound to the code that produced it
symfun.rs:236-249, 398-424, 429. load validates only integrity (hash) and the cheap basis/weights. e_mul and
f_mul, the actual structure constants, are never checked against anything. If GrTable::compute or a helper it
uses (poly, solve_mod) changes without a manual GR_FORMAT bump, or if KRHOM_GRCACHE points to a directory shared
with another branch or worktree, an old table is loaded silently. sbim/local/algebra then use wrong
multiplication tables, and the M6a-2a binding of checkpoints to commit/library hash does not cover this file.
compute has not changed since the initial commit, so the risk is latent, not live.
Fix: on load, run a cheap algebraic self-check, for a few random basis vectors a and j:
sum_i e_i f_{j-i} a = c_j a, and e_i e_k a = e_k e_i a, computed with the loaded tables. Alternatively,
put a build.rs hash of symfun.rs/poly.rs into the header.

### m2 (minor): no per-use spot check for theorem-certified Assoc, despite the comment
leafgen.rs:722-725 says spot_check_assoc "tests the bimodule property of each (3,3,3) use", but
certify_assoc_theorem / add_thm_assoc (leafgen.rs:795-818), the realize_chain theorem branch (leafgen.rs:1218)
and the smapio import re-check (smapio.rs:458) never call it. Approval condition (ii) was met once, on the
(3,3,3) braid middle step (spot333.log). A later chain with a different > 500-flow up-down Assoc would be
certified on the theorem alone. No current map uses the theorem (census max Assoc 60 flows; every B4
theorem_assoc = []), so nothing is affected today.
Fix: call spot_check_assoc inside certify_assoc_theorem whenever flows > ASSOC_CERT_FLOWS (17 s at
1 680 flows), or correct the comment and record that (ii) is a one-time cross-check.

### m3 (minor, not reached by the N = 5/6 lists): planner budgets ignore the GrTable(r, b >= 8) build cost
Template windows are capped only by flows (<= 1 200, rexplan.rs:532) and Exact windows only by unknowns
(rexplan.rs:579). A window in which some rung has a giver label of 9 is admissible: for example canonical_max or a
merge-all/split word on a (3,3,3) sub-window with a top such as (6,3,0) has 84 flows. Building it calls
GrTable::get(r, 9) through sbim. The report says GrTable(3,9) and (6,9) did not finish in 90 min (2-3 GB).
GrTable::get holds the global table mutex during the build (symfun.rs:194-205), so every other GrTable user in the
process blocks too. From the logs I computed the maximum giver label over all census and B4 moves: it is 6,
so nothing hits this today. B5's general Lemma A words could.
Fix: in neighbours, reject (or cost) windows whose words have a giver label >= 8 unless
GrTable::cache_path(r, b) exists.

## Non-defect notes (no action needed for correctness)
- exact_candidates / assoc_candidates are memoized on (labels, top) only (rexplan.rs:414, 494), but maxlen
  depends on end - at (rexplan.rs:568). Whichever window size comes first fixes the candidate set, so the census
  reach can depend on visit order. Every chain is still checked when it is realized.
- all_words says "prune: ... L1 distance" (rexplan.rs:722), but no pruning is implemented. This costs speed only.
- rex_class_pairs hard-codes s = 0, mult = 1 and gp0 --b4 does not compare them with tg.s / tg.mult. This is correct
  because reduced_for implies c_q = 1, but an assert would protect the provider wiring planned in report §6.3.
