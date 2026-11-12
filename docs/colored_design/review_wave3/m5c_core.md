# Review wave 3: M5c changes to the frozen core (revm5, colored-m5 f2fcb1d vs dev 370f5b1)

Scope: src/linsys.rs (kernel_lazy sparse sweep), src/sbim.rs (hom_space row dedup), src/tmpl.rs
(Arc<Side> sharing), src/homprov.rs (Exact admission by call kind, leaf_iso Refused), and the tests
that cover them (tests/hom_lazy.rs, tmpl::m5c_tests, tests/homprov.rs m5c test). Read-only: no cargo run.

## Verdict

No defect in the core changes that can produce a wrong Poincare polynomial, a false certification,
silent corruption or a crash/hang. The solver changes preserve results exactly, the side sharing
changes no template data, and the homprov changes affect only which calls are attempted and how a
failure is classified. Two minor test-quality findings follow.

## Verified (with reasoning)

### linsys.rs kernel_lazy_traced sparse sweep (lines 241-308)
- Columns enter the min-heap once (`mark`), at first touch. A pivot row stored at lead l holds only
  columns > l (row built from `rest`, which contains only columns popped after the lead, line 286).
  So elimination at column `col` only touches columns > col, and pops are monotone increasing. A
  column's value is final when it is popped, and a popped column can never be touched again. The
  `continue` on zero (line 268) therefore cannot lose a later nonzero.
- The lead is the first nonzero non-pivot column, and every later nonzero pivot column is eliminated,
  as in the dense scan. `rest` gets pushes only after the lead is set (if/else-if/else chain). The
  stored row is (rest nonzero) * inv(lead), which equals the old `(l+1)..=hi` scan, because every
  column the old scan saw as nonzero was touched. Arithmetic is exact over F_p, so the pivot rows
  are identical.
- The heap is drained fully in every insert. `touched` resets both `acc` and `mark` (line 304),
  including on the redundant-row path. There is no state leak between inserts or into the
  reinsertion loop.
- Early stop (lines 384-407) is unchanged. Rows before i are all inserted; every not-done row
  j >= i is checked against every kernel basis vector before `early_stop`. The final kernel is
  kernel_sparse(piv): one vector per free column, in increasing order. That is the unique reduced
  basis given the row space and the column order.
- `kernel_lazy` has exactly one caller, sbim::hom_space (sbim.rs:328).

### sbim.rs hom_space row dedup (lines 298-325)
- Rows are `merged` after `retain(nonzero)`, and empty rows are never pushed (sbim.rs ~287-291). So
  r[0] exists and r[0].1 != 0. `Fp::inv` asserts on zero (field.rs:114-117). Scaling by a nonzero
  scalar keeps every entry nonzero.
- `merged` is sorted by column, so r[0] is the leading column. After normalization, two rows are
  equal iff they are scalar multiples, and Fp/Q equality is canonical (Fp values always reduced;
  BigRational normalized).
- The sort key (len, column sequence) is a total order on patterns, so equal patterns are
  contiguous. The `start` run logic compares each row against every kept row of its run, so all
  duplicates are dropped and the first occurrence is kept.
- The row space is unchanged, so the kernel is unchanged (see above). Rows of equal pattern still
  come in HashMap order, as before; this affects only the trace, not the result.
- Run-local O(run^2) comparison: rows of length 1 collapse to a single kept row. Distinct
  same-pattern rows with integer-structured coefficients are few. This is not a hang risk.

### tmpl.rs Arc<Side> sharing (lines 131-147, 370-389, 415-440, 560-593, 621-624)
- `sides` is keyed by `Window {labels, rungs}`, which fully determines `Expr::new(labels, rungs)` and
  hence `sbim()`. The only non-None SBim passed in (try_hom_basis_cap, lines 477-496) is
  `Expr::new(&src.labels, &src.rungs)?.sbim()`, the same object that `side()` builds for None
  (`w.expr().sbim()`). The determinism of `sbim()` across calls is already relied on by the disk cache
  (phi loaded from disk is applied to a freshly built SBim).
- TmplLib<F> is per field, so sides of different fields cannot mix. `aux` is only appended, never
  truncated or rolled back. certify_added and localize read `a.src`/`a.tgt` exactly as they read the
  old per-template fields.
- `memory_bytes` now counts each side once, which matches the real memory.

### homprov.rs admission / classification
- `call_kind`, the kind-keyed buckets, `predict_kind` and `predict_call` (lines 1040-1137) only feed
  Exact::admit (1269) and Search::win_ok (1768). An admitted Exact call still returns the exact
  hom_space basis and goes through the same postchecks (sindec dechom_postcheck, certify_closure).
- leaf_iso (2383-2412): `best` is chosen by exactly the same condition as before. The only change is
  that when `best` is None, the error kind is Refused instead of Exhausted, and the cheapest refused
  window goes into the trace. MissKind is consumed by Ctx::end_kind/stuck_by (1620-1640) and
  sindec::prov_fail `engineering` (sindec.rs:372-377), and then only in libbuild outcome/retry
  classification. No acceptance path depends on it.
- The S4 closure leaf_iso caller (2979) only logs the error, so nothing changes there.
- The layer-A route stamp (sindec.rs `route_header`, `layer_a_load`) is outside this lens but was
  checked. It keeps the default (exact header) and, when enabled, ignores only the `costs=` token.
  Loaded pairs still go through dechom_postcheck (sindec.rs:721) and closures through certify_closure
  (sindec.rs:1080).

### Tests
- `kernel_lazy_property` (pre-existing, unchanged) now runs the new sparse sweep against
  `Echelon::kernel`, an independent HashMap/RREF implementation (linsys.rs:9-100). It covers
  adversarial orders, full rank, no equations and kernel dimension > 128. This is a real oracle test
  of the heap sweep.
- `m5c_dedup_sparse_sweep_equals_reference` compares hom_space (dedup + heap) with
  hom_space_reference (no dedup, Echelon) by Debug string on 4 window pairs, so it really tests
  equality. It does not assert that any row was actually deduplicated ("most rows are duplicates"
  is not checked), so the claim that the removal branch is exercised is not verified by the test.
  This is a note, not a finding: equality is the property that matters.
- `tmpl::m5c_tests::hom_basis_templates_share_window_sides` is a real test (Arc::ptr_eq, side count
  stable across a second call). Localization equality is covered by m4a A4 (unchanged), which
  compares localizations with the oracle.

## Findings

### F1 (minor) tests/hom_lazy.rs:58-61: the hom_space early-stop coverage assertion can be satisfied by the new M5c test
`lazy_kernel_equals_reference` asserts `early >= 1` and `reins >= 1` from `kernel_lazy_counters()`.
These are process-global atomics (linsys.rs:226-235), and the test's own comment says the branches are
"exercised by these hom_space calls". Before M5c it was the only test in the binary that called
`kernel_lazy` (the property test uses `kernel_lazy_traced`, which does not touch the counters). M5c
adds `m5c_dedup_sparse_sweep_equals_reference`, which also calls hom_space, and libtest runs the two
tests in parallel threads (no test-thread config in the repo).
Failure scenario: row dedup removes the scalar-duplicate rows, which are the cheapest redundant rows.
If that leaves fewer than 64 consecutive redundant rows in the reference test's windows, its own
calls no longer take the early stop or reinsertion branches. The assertion still passes whenever the
M5c test's calls have incremented the counters first, so a loss of hom_space-level early-stop coverage
goes undetected. The result depends on scheduling. Branch coverage itself is still guaranteed
separately by kernel_lazy_property (early >= 10, reins >= 5), so this is test-quality only.
Fix: compute counter deltas inside the test while holding a static Mutex shared by both tests, or
make lazy_kernel_equals_reference call kernel_lazy_traced on the hom_space equations and count its
own traces.

### F2 (minor) tests/homprov.rs:1007-1012: the "censored exact shape refuses even if the kind bucket were fast" check is vacuous
`t3` records a fast 64x64 call (fibre 4096, bucket 12) and a censored 60x70 call (fibre 4200,
bucket 13). The censored record sets its own kind bucket (kind 1, bucket 13) to +inf
(CostTable::record, homprov.rs:1083). The fast record is in a different bucket. So both the
exact-shape path and the bucket path give +inf, and the assertion cannot tell exact-first precedence
from bucket-only logic. More generally, a censored exact shape always makes its own kind bucket +inf,
so the situation the comment describes cannot occur. The precedence that matters (a fast exact shape
inside a +inf bucket admits) is tested by `ex2`. No product impact.
Fix: drop the sub-case, or reword it as "a censored shape also closes its kind bucket".

## Not findings (checked, out of scope or harmless)
- `COST_ALG` stays `hs-kl1` although the solver is faster. Old tables over-predict, which is
  conservative. The new `costs_333_m5c.txt` carries the same stamp, so a pre-M5c binary would
  under-predict with it. This is admission only.
- tests/sindec.rs:790 sets `KRHOM_COSTS_ANY_HOST` process-wide in a parallel test binary. No other
  test in that binary asserts host rejection, so no claim is falsified.
