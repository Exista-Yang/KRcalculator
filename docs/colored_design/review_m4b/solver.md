# M4b review: linsys::kernel_lazy and its use in sbim::hom_space

Branch colored-m4b, HEAD 8ad4ac6, diff 3c3323b..HEAD. Read-only review. No cargo was run. I wrote a
faithful Python port of `kernel_lazy` to cross-check the reasoning (well under 20 s; details in §6).

## Verdict

The solver is **exact**. It returns **the same reduced kernel basis, in the same order and
representation** as the old `Echelon::kernel` path. I found no defect that can produce a wrong
Hom basis, a false certification, or a hang. The template stamp `TMPL_ALG = "hs1.lcx1"`
(src/tmpl.rs:183-186) was correctly left unchanged.

There are two minor findings:
- a memory-scaling regression in the final dense kernel;
- a test-coverage gap: the new early-stop branch is never shown to be exercised.

## 1. Equation generation (sbim.rs:202-322 vs hom_space_reference sbim.rs:326-444)

- `hom_space_reference` is byte-identical to dev's `hom_space`. I checked this by extracting
  both function bodies and diffing them.
- The new `hom_space` differs from the reference in exactly three places:
  1. `Echelon::new` is replaced by `all_rows` (sbim.rs:228);
  2. `ech.insert(&merged)` is replaced by `all_rows.push(merged)` (sbim.rs:290);
  3. a sort plus `kernel_lazy` replaces `ech.kernel()` (sbim.rs:296-297).
- The unknown order (`unk`, sbim.rs:208-222) is the same in both. The rows are the same merged,
  zero-free sparse rows. The PMat assembly (sbim.rs:298-321) is identical.
- All rows of all operator pairs exist before the solver runs, so the early stop always sees
  the complete system. No rows are generated lazily or in chunks.
- Memory: all rows are now held at once. The reference instead freed each operator's rows after
  inserting them. Measured RSS is acceptable: 478 MB at 112 536 unknowns, the censored record in
  tests/data/homprov/costs_cloud.txt.
- The comment at sbim.rs:295 says the order is "deterministic". That is not quite true. Rows
  with equal length and equal index pattern but different coefficients keep the per-process
  randomized HashMap order, because the sort is unstable. Only timing is affected, not the
  result, since the basis is unique.

## 2. Field arithmetic (field.rs)

- `Fp<P>(u64)` keeps canonical values in [0, P). `add` (field.rs:98-101) sums to less than 2P,
  which is below 2^33.
- `fma` (field.rs:129-131) computes `(self + a*b) % P`. The largest value is
  (P-1) + (P-1)^2 = P^2 - P.
  - For F32b, P = 4294967291 = 2^32 - 5, so P^2 - P = 2^64 - 11*2^32 + 30, which is below 2^64.
  - So there is no overflow for F31 or F32b, in either `acc` or the back substitution.
- No code outside field.rs builds a raw `Fp(x)`; I grepped src/. So unreduced values cannot
  reach the solver.
- Over Q the code is also correct. The dense buffers are just much larger (see F1).

## 3. Elimination (`insert`, linsys.rs:204-254)

- **Precondition.** `acc` is zero on entry. Every write goes to [lo, hi]:
  - the row's own indices;
  - pivot-row entries j > col >= lo, with `hi` widened as needed (linsys.rs:226-228).
  
  [lo, hi] is cleared at the end (linsys.rs:250-252), so the precondition holds for the next
  call. Duplicate indices would be summed (`add_assign`, linsys.rs:212); hom_space merges them
  in any case.
- **Pivoting.** The scan runs over increasing columns:
  - every nonzero column that has a pivot is eliminated with the implicit leading 1
    (`c = -acc[col]`, `acc[col] = 0`, `acc[j] += c*y`);
  - the first nonzero column without a pivot becomes `lead`.
  
  Pivot rows store only entries strictly to the right of their pivot. Elimination at `col`
  therefore only touches columns above `col`, so a single left-to-right pass reduces the new row
  against every existing pivot. The stored row (linsys.rs:239-246) is normalized: it has its
  lead at l and entries only at columns that were free at insertion time.
- Older pivot rows are not back-reduced, so the stored form is echelon but not reduced. That is
  harmless:
  - the pivot set equals {lead(v) : v in rowspace}, an invariant of the row space;
  - the old `Echelon` has the same pivot set, since its pivot is also the smallest column
    (linsys.rs:57-65).

## 4. Kernel and normalization (`kernel_dense`, linsys.rs:256-279)

- For each free column f, the code sets x[f] = 1 and x = 0 elsewhere, then back-substitutes for
  pivots p < f in descending order: x[p] = -Σ row_p[j] x[j].
- Pivots p > f stay 0, which is correct. By downward induction, x[p] for p > f depends only on
  coordinates j > p > f, and those are all 0.
- The pivot rows for p < f may contain later pivot columns. Those coordinates are already
  computed when they are needed, because the loop runs in descending order.
- Result: the unique kernel vector whose free coordinates are e_f. This is exactly what
  `Echelon::kernel` (linsys.rs:84-101) computes from the RREF, as e_f - Σ_p rref_p[f] e_p.
- Output representation matches:
  - free columns come in increasing order;
  - entries are sorted with zeros filtered (linsys.rs:329-331);
  - the old version drops zeros in `Acc::take` (sparse.rs:48-61).
  
  The bases are therefore identical element by element, and cached templates remain valid.

## 5. Early stop (linsys.rs:289-326)

- **Trigger** (linsys.rs:305): at least `run_max` (64) consecutive redundant rows, kernel
  dimension at most 128, and at least one row left.
- **Check** (linsys.rs:308-309): every row j >= i with `!done[j]` is tested against every
  kernel vector (`violated` uses `any` over all of `ker`, linsys.rs:280-288).
  - Rows below i are all done: the main loop sets `done[i]` before `i += 1`.
  - Rows marked done were inserted, so they lie in the span.
  - So when the check reports no violated row, every row of `eqs` vanishes on
    K = ker(inserted). Then ker(all) = ker(inserted), the row spaces are equal, the pivot sets
    are equal, and the bases are equal (§4).
- **Termination and cost.** A violated row is not in the current row space, so the first insert
  in `bad` raises the rank.
  - In fact, after one round every remaining row was either violated, and is now inserted, or
    not violated, and already lay in the old row space.
  - So the row space equals span(all rows) after one round, and the second check (linsys.rs:308)
    always comes back empty.
  - The inner loop therefore runs at most twice, with no hang. Each round costs one
    `kernel_dense` (at most 128 × ncols × 8 B, about 41 MB at 40 000 unknowns) and one scan.
- Once triggered, the scan always finishes. `run_max *= 2` (linsys.rs:324) is unreachable, and
  the doc comment's "the check repeats" refers only to the inner loop. This is not a defect.
- **Degenerate cases.**
  - No equations: identity basis, same as the old solver (but see F1).
  - Empty kernel: `kernel_dense` returns [], so the Hom is empty.
  - All rows redundant or empty: the empty-row guard at linsys.rs:205-207 handles them, and
    hom_space never passes empty rows anyway (sbim.rs:288-291).
  - ncols = 0: hom_space returns before calling the solver (sbim.rs:223-225). `kernel_lazy(0, [])`
    returns [].
  - `ncols - rank` cannot underflow, because rank <= ncols.

## 6. Cross-check by simulation

The script is
/private/tmp/claude-501/-Users-existair-Documents-KRcalculator/bf6f8870-08ad-436b-9ff7-0414655cc59a/scratchpad/colored_design/m4b_review/sim/port.py.
It is a line-by-line Python port of `kernel_lazy` (same constants 64 and 128, p = 2^31 - 1),
compared against an independent dense RREF kernel.

- I used 120 random adversarial systems: 20-140 columns, most independent rows first, then
  60-230 redundant rows, then 0-4 independent rows placed near the end.
- The two kernels were identical in all 120 cases.
- The early stop fired in 118 cases, and violated rows were inserted after the trigger in 94.
- The inner loop never ran more than 2 iterations.
- The degenerate cases (0 columns; no equations; full rank) also agree.

## 7. Dense accumulator memory

- `acc` is ncols × 8 B: about 320 KB at 40 000 unknowns and 0.9 MB at 112 536. This is
  negligible.
- The check-time kernel is bounded by 128 × ncols × 8 B.
- The **final** `kernel_dense` is not bounded; see F1.

---

## Findings

### F1 (minor): the final kernel is materialized densely, all vectors at once

**Where:** linsys.rs:256-279, called at linsys.rs:328.

When the early stop does not fire (true kernel dimension above 128, or never 64 redundant rows
in a row), `kernel_dense` allocates `dim` dense vectors of length ncols before converting any of
them to sparse.

- **Peak memory:** dim × ncols × 8 B over Fp. Over Q it is about 64 B per entry (BigRational).
- **Old solver:** O(nnz of the RREF and the sparse kernel).
- **No-equation case:** ncols² × 8 B. Example: `hom_space(id, id, d)` from
  `StateSpace::compute` (sbim.rs:577-581), which the default `krhom` cube path reaches via
  cube.rs:142. Its old cost was O(ncols).

**Failure scenario.**
- A Hom with dim ≈ 3 300 at 40 000 unknowns (the label-rule cap, homprov.rs:997) needs 1.06 GB
  of transient memory. That is above the laptop G1 watchdog `--rss 1024` (HANDOFF_M4b.md:19).
  The call is killed and recorded as a Memory miss, not a wrong result.
- `cclose::nakayama` (cclose.rs:2128-2149) calls `hom_space(id, b, d)` in high degrees. There
  Hom(R^λ, B)_d ≈ rank × monomial count, so dim/ncols ≈ rank/fibres. For M5/M6a-size ladders
  this can reach many GB.
- Identity state in the cube at n = 5 strands, N = 6, d = top: ncols = C(29,4) = 23 751, so
  4.5 GB.
- The cost-table branch of `Exact::admit` (homprov.rs:1005-1007) admits calls whatever their
  unknown count, so ncols above 40 000 is reachable.

**Current exposure:** small. G1 peak RSS was 856 MB, and the cost-table maximum is 117 MB for
calls of at most 27 494 unknowns.

**Fix:** reuse one dense buffer per free column and push its sparse form immediately. Keep the
dense form only for the 128-bounded check, or evaluate the check from sparse kernel vectors.

### F2 (minor): tests/hom_lazy.rs never shows that the new early-stop path ran; there is no direct unit test

**Where:** tests/hom_lazy.rs:22-55.

The test does compare map by map against `hom_space_reference`:
- dimension check, then full derived-`Debug` equality of every PMat (pmat.rs:9-10, mpoly.rs:100);
- 10 web pairs over several degrees, with a `total > 20` guard.

What it does not do:
- It never asserts or counts that `kernel_lazy` took the early stop (linsys.rs:305-323). It never
  checks that rows were inserted after a trigger (linsys.rs:314-319).
- There is no unit test of `kernel_lazy` on synthetic systems: the adversarial order (an
  independent row after more than 64 redundant ones), no equations, full rank, or
  kernel dimension above 128.
- Only F31 is used.

The larger cases (for example [2,2] E0.2,F0.2 -> E0.1,F0.1,E0.1,F0.1, about 10^3 unknowns and
several thousand rows) very likely trigger the early stop, but the test cannot show it. If the
early stop or the violated-row insertion broke, the test could still pass on inputs that never
reach those lines.

My port (§6) shows the current code is correct, so this is a coverage gap only.

**Fix:**
- add a small counter or return flag from `kernel_lazy`, and assert that at least one case
  triggers the early stop with violated rows;
- add a direct `kernel_lazy` vs `Echelon::kernel` property test covering the adversarial order
  and the degenerate cases.

## Non-findings (checked)

- Overflow in F32b `fma` and in the accumulators: not possible (§2).
- Early stop skipping rows: not possible; every unprocessed row is checked against every kernel
  vector (§5).
- Basis normalization differing from the old solver: not possible (§4). The cache stamp is
  correctly unchanged.
- Hang in the inner loop: not possible; at most 2 rounds (§5).
- Duplicate or zero entries in rows: both are handled.
