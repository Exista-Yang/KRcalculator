# Review wave 4: M6a-1c engine speedups (revm6c, colored-m6a 11921fa)

Scope: `git diff c95517f..HEAD` restricted to the M6a-1c commits (bf19c46, d285ba6 and the run commits);
code changes in src/sengine.rs, src/swire.rs, src/bin/sprof.rs, src/bin/sbox.rs (one flag). The merged
colored-m5 content is excluded because it was already reviewed. No tests changed in this range. Read-only; nothing was run.

## Verdict

No defect found that could change a Poincare polynomial, give a false certification, or silently
corrupt a component. All four speedups are exactly result-preserving: same values, and same pivot order in Record mode.
Replay PointDependent detection is intact. There are two **minor** verification/test gaps (below).

## What was checked (with the equivalence argument)

### 1. d_C (x) id: `dcid_step` / `assemble` (swire.rs:648-855) vs the per-pair reference (KRHOM_DCID_REF, swire.rs:444-571)
- Cache `asm` keyed by (Tid, coset). It is cleared at the start of every term k (swire.rs:436), and `decs` is keyed by Tid
  for the same k, so each (Tid, coset) has exactly one Decomp. The cache is local to `scan_point`, so its values
  are per point and no point-dependent value leaks to another point.
- Row and column counts match. `assemble` uses kr = |flows(dec.x) at cc|. That equals `kron_idx`'s cols
  (a side; colinfo length) and rows (b side) (swire.rs:360), and it equals `Eval::block`'s rows/cols for
  iota/pi (smap.rs:1099). The reference made the same assumptions (`ib.a[g*m_s..]`, `pb_.a[i*kr+r]`).
- The arithmetic is the same. KI_all = (f(x)id)[iota_1|...|iota_n]. Its column block s is exactly the reference KI_s. The
  stacked pi_all*KI_all has block (it, ia) = pi_it (f(x)id) iota_ia. Per entry, the accumulation order (r ascending, then
  colinfo/rp order) is the same as the reference, and the fields are exact anyway (Fp, Q). The `last_f` index cache
  (swire.rs:761-770) only skips a repeated `f.lay.index(cp)`.
- Zero and skip decisions are the same. The reference inserts a coset block iff pi_it*KI_ia != 0. The new code copies
  the block iff it is nonzero (swire.rs:822). A LocMap is created lazily on the first nonzero block, and `add_map` is called
  in the same (a, b, ia, it) order (swire.rs:842-850). The extra `any = true` for all-zero iota rows only yields
  an empty `nzr` and a zero `out`, so nothing is inserted.
- Scratch buffers do not leak stale data. `ki`, `out` and `nzr` are `clear()`ed and `resize`d with zeros per coset
  (swire.rs:758-759, 797-798). Each slot block is overwritten once per coset with `clone_from_slice`, and every coset is visited once
  (sorted support).

### 2. Elimination: `CompPlan` and the per-pivot plan cache (sengine.rs:782-827, 1540-1611)
- The plan key is (x type, pivot type, y type). g = d_ay has layout reg(p_a, p_y), and f = d_xb has layout reg(p_x, p_b) with
  p_b = p_a, which `apply_pivot_inner` enforces (sengine.rs:1520-1522). Every stored component has a layout equal in content
  to the registry layout of its type pair: `PEngine::add_map` rejects others (sengine.rs:1029-1031), and fill-ins
  use `reg.layout`. So the offsets cached in the first plan are valid for every later use. Shapes are registered
  once per Tid (`ShapeReg::register` refuses a second shape), and the engine is per point.
- `CompPlan::apply` follows `compose_acc`/`gemm_acc` loop for loop (same coset set = dst cosets where g and f
  have blocks; alpha != 0 always). It has the same zero skips and the same `fma`.
- Deferring the removal of zeroed edges to after the y loop (sengine.rs:1601-1608) cannot change anything:
  each (x, y) is visited once per pivot, and x != y by the h grading.
- `touched` as a sorted, deduplicated Vec with a, b removed gives the same set and the same order as the old BTreeSet.

### 3. Buffer pool (sengine.rs:1580-1590, 1627-1635)
- A reused buffer is `clear()`ed and then `resize(n, 0)`d (n = lay.len = the pool key), so no stale values remain. Buffers move by
  ownership (`take_map` -> `recycle`). They are never shared with `prev` (Arc<SCx>) or with live components, so no
  aliasing is possible. `cx.bytes` uses capacity both at insertion and at removal, so the counter stays consistent.
  The pool is driven by structure only, so Record and Replay behave identically. Its 16 MB counts toward the budget, which
  can only cause a refusal, never a wrong value.

### 4. `markowitz_heap_fast` with compaction (sengine.rs:1304-1410) vs `markowitz_rescan`
- Candidates: edges with key[x] == key[y] (type, t) and elig[x]. This is the same as `is_candidate`. Eligibility depends
  only on `steps` and `cur_step`, which are fixed during a round. The cost formula matches `cost()`, and `max(1)` matters only for
  entries whose edge no longer exists, which are dropped by `cx.get(a,b).is_none()` (sengine.rs:1357).
- Invariant: every current candidate has an entry with current versions and current cost. Every object
  whose degree changes is in `touched`. Fill-in and zeroing affect only x in in(b) and y in out(a), and the removal
  of a and b affects only their neighbours. All touched versions and degrees are updated before any re-push
  (sengine.rs:1380-1399). `c_out`/`c_in` contain every key-equal edge (initial edges plus `new_edges`). An edge can
  only be (re)created by a pivot that touches both ends, so an entry with current versions has the current cost.
- Compaction (sengine.rs:1343-1352) keeps exactly the entries that are alive and version-current (and so valid),
  then sorts and deduplicates them (duplicates are identical tuples), so the pop order of the valid set is unchanged. The popped entry
  is outside the heap while compaction runs. Ties are broken by (cost, a, b), the same as rescan.
- Conclusion: the pivot sequence equals the rescan order (and HeapV1).

### 5. Replay PointDependent detection
- This code is unchanged: the recorded pivot must be alive, have an edge, and be a candidate (sengine.rs:1212-1215); Pivot errors map to
  PointDependent; the pattern hash is checked after every round (1228) and at every checkpoint (914-935); the stage
  determinism hash is checked; and `finish` requires rec == frozen. Zeroing in the new code is value-dependent exactly as
  before, so it is still visible in the pattern.

## Findings (minor)

### m1. No test reaches heap compaction
`compact_at = max(2*len, 1<<16)` (sengine.rs:1339, 1351). The synthetic complexes used by
`markowitz_heap_equals_rescan_order` (tests/sengine.rs:144, c_params/w_params: <= 4 types, a few dozen
objects per factor) never come near 65 536 heap entries. No test anywhere reaches 65 536 either: the swire tests use
short words, and the full suite runs with `--release`, so debug_asserts are off. The test passes on 400 scans, but it says
nothing about the compaction branch, which is the one place where order preservation depends on a predicate
(`retain` at 1347). Its only evidence is the sprof sigma^5/sigma^6 record runs.
Failure scenario: a later edit to the retain predicate (for example also dropping entries with no edge, or a
version off-by-one) drops a valid candidate during the final `all` round. The pivot is then never taken, the
complex is left non-minimal, and every unit test still passes.
Fix: make the threshold a test-visible parameter (EngOpts or a `#[cfg(test)]` constant) and run the
heap=rescan test with threshold 0 (compaction after every pop).

### m2. The before/after "record trace digest" comparison cited in M6a1c_report.md section 2 does not exist
The report says the value digests "and the record trace digest are identical before and after at sigma^4, sigma^5,
sigma^6 (runs_m6a1c/prof_s*.base.txt vs *.after.txt)". The base files contain no `record trace digest` line,
because sprof printed it only from d285ba6 on. So the before/after identity of the record pivot sequence
was never compared directly. The conclusion still holds on the evidence that does exist: the per-stage
inserted/objects/components/pivots/entries are identical, the update and multiply-add counts are identical
(302 089 / 1 534 101 / 5 950 811), and the replay value digests are identical. Those replay values depend on the frozen
pivot sequence, and the closure entries cover every closable stage, so the digests are meaningful. The
`KRHOM_DCID_REF` "(same digests)" claim likewise has no artifact, and no test runs the reference path, so
the reference will drift unnoticed. The new path itself is covered end to end by tests/swire.rs
(`box_loop_against_exact_22_32` and the xengine comparisons).
Fix: record base/after record-trace digests (or remove the claim), and add a small test that runs
`scan_point` with and without KRHOM_DCID_REF and compares `stages[*].vals` and `ps.frozen`.
