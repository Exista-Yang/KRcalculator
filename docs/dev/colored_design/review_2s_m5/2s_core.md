# Review: colored-2s (M2S′) engine core, xengine.rs equivalence to the full build

Branch colored-2s at 0fb5dc4 (worktree .claude/worktrees/rev2s), diff ce91f0e..HEAD. Read-only review:
no cargo, no binaries. Line numbers refer to src/xengine.rs at 0fb5dc4 unless a file is named.

## Verdict

**No defect found that can give a wrong Poincaré polynomial, a false certification or silent corruption.**
The streaming build, the compact store, rung-by-rung transport, threaded transport and fill-in, the lazy
heap and the sort-based product are all exactly equivalent to the previous full build: the same
operations, or the same eliminations in a different order. In the second case the minimal complex is the
same up to isomorphism, so the closure is the same. There are two minor findings, both about the test
and the provenance record, not about the engine.

As an independent check, I recomputed the published results from the JSON files with a Python script of
under 1 s. Every Wang block present in tests/data/wang_ref.txt matches the stage closures exactly
(σ^9: 21 comparisons, T(2,3/5/7/9) at every N that has a block). σ^9 F31 equals σ^9 F32b on all 96
closures. σ^7 F32b and all three σ^5 files, including the threaded one, are subsets of σ^9 F31. In
every file the "str" and "terms" fields agree.

## Verified, no defect

1. **Eligibility rule, streaming (1353-1462).**
   - The slices are by h(o) of the C object. Slice s brings:
     - all objects (o,k,·) with h(o) = s;
     - id⊗d inside each o⊗X;
     - d_C⊗id from slice s−1.

     So at the end of pass s:
     - every object of a slice ≤ s has all its incoming raw components;
     - every object of a slice < s has all its outgoing raw components.
   - out_done = slice < s and in_done = slice ≤ s (1453-1454). Note that in_done is always true by
     construction.
   - The raw components still to come all satisfy slice(u) ≥ s and slice(v) ≥ s+1. Let a′ → b′ be a
     pivot eliminated in pass s′ ≤ s. Then:
     - slice(a′) < s′, so a′ is never a later source;
     - slice(b′) ≤ s′, so b′ is never a later target.

     A later raw component therefore never feeds d_{x b′} or d_{a′ y}. It enters an elimination only as
     the additive term d_{xy}, so linearity gives the same complex as the full build with that pivot
     order.
   - The only later raw components that touch an eliminated object start at some b′ of slice s−1.
     These are exactly the ones the `contains_key(&srcs[sx])` check skips (1421). The full build drops
     them too, since components out of b vanish.
   - A later component into an eliminated a′ is impossible. If it happened, `add_map` would panic
     (fail closed).
   - The last pass runs `simplify(None)`, so the result is fully minimal. Minimal complexes are unique
     up to isomorphism (Krull-Schmidt; End^0(B_p) = k).
   - Gaps in the h-slices (`hs[si-1] != hs_cur-1`) are handled.
   - Signs (−1)^{h(o)} and the (ks, kt) orientation match the full build (1329-1338 vs the stream body).
2. **Rung-by-rung transport (`transport_view`, 711-756).**
   - The identity π_t(f⊗id)ι_s = π_2((π_1(f⊗id)ι_1)⊗id)ι_2 holds exactly, because `decompose_word`
     stores ι = (ι_1⊗id_fr)·ι_2 and π = π_2·(π_1⊗id_fr) without re-orthogonalizing (1053-1056).
   - kron_id is multiplicative and associative in the row·k+b convention.
   - Summand order [first × rest] matches `offs`.
   - Killed or empty sub-decompositions are handled (the early return; offsets use rest lengths).
   - The kron factor is asserted (720-721).
   - `rw` falls back to direct transport when `rungwise` is false (1411), because the composite is
     then a direct decomposition.
3. **Threading.**
   - `par_map` (681-704) returns its results in job order, and transport results are added in job
     order (1417-1425).
   - Fill-in pairs (x,y) are unique per pivot. They are read before and written after each chunk, in
     pair order (1522-1543). The results are deterministic, and the field operations are exact and
     commutative.
   - No shared mutable state is touched in the workers except the monomial interner (an RwLock;
     read → drop → write with a re-check, 341-376; no nested lock acquisition, so no deadlock).
4. **Compact store (`CMat`, 320-418).**
   - Fp is always canonical in [0, P): every constructor in field.rs reduces, and no direct `Fp(x)`
     construction exists outside field.rs. P(F32b) = 4294967291 < 2^32, so `x.0 as u32` is lossless
     and `from_i64` restores the value exactly.
   - The tptr and colptr casts happen before the final length assert. That assert covers the maximum,
     so a truncation would fail closed.
   - Interning is keyed by the full u128 monomial, so there are no collisions, and indices are never
     reused. The table is ring-agnostic, so sharing it across scans and fields is safe.
   - A u32 index wrap needs 2^32 distinct monomials (≫ 100 GB), which is not reachable.
   - Other fields stay `Plain`.
5. **Sort-based product `mm` (622-674).**
   - It sorts by (row, mono) and sums equal keys. The constant scale is applied after the sum (by
     distributivity), and zero sums are skipped.
   - Row grouping after a cancelled group is correct (the `cur` handling).
   - The output is in canonical form, the same as `PMat::matmul`. The fill-in formula
     `mm(d_ay, d_xb, −c^{-1})` equals the old `scale_const(h.matmul(g), f)`.
6. **Lazy heap (1483-1565).**
   - Every edge whose cost (|in(v)|−1)(|out(u)|−1) can change has u ∈ in(a)∪in(b) or
     v ∈ out(a)∪out(b). Both sets are captured before fill-in (1518-1519), and every current edge on
     them is re-pushed (1550-1563), including new fill-in edges.
   - A popped entry is used only if its edge exists and its cost is current (1506-1511).
   - So the pivot sequence equals the full rescan's: min (cost, a, b) over eligible edges, with the
     same tie-break.
7. **Certification checks in every mode.**
   - `check_degrees` runs on every raw component before any elimination: per pass in streaming (1446),
     before `simplify` in the full build.
   - d∘d is checked on the final complex in both modes (1797-1808).
   - The K_0 check and the lower bound run in both modes.
   - Stage closures (1864-1879) run only after these checks, use `crossings[..m]`, and keep the pruned
     guards (`close_stage_marked` 1956-1972, `close_guarded` 1977-1989).
8. **Tests.**
   - `streaming_equals_full_build` compares final object multisets, K_0, and reduced and unreduced
     polynomials. It covers stream + compact + rung transport, stream + plain + direct, and full +
     compact against the full plain baseline, on mixed-sign and pruned words.
   - The helpers moved to tests/common/mod.rs are byte-identical to the originals.

## Findings (minor)

### m1. The G2X test can pass without comparing against Wang (tests/xstream.rs:226-227)
`if let Some(want) = wang_ref().get(&(3, nn, mm))` skips silently when a block is missing.
- If the block 3 6 2 were absent or mis-keyed in tests/data/wang_ref.txt, the gate's headline value
  (T(2,5) = 396 at gl_6) would not be checked, and the test would still pass.
- The SPEC gate item "rank 1 at N = 3" for T(2,5) is not asserted either: there is no block 3 3 2, and
  no fallback rank check.
- T(2,7) and T(2,9) against Wang exist only as a manual check of the JSON. I re-verified them: they
  pass.

Fix:
- assert that the blocks exist for N = 4..8 (and stage 3 at N = 3..8);
- assert total rank 1 at N = 3;
- for M2S-run, add a test or script that checks the gl6_t2n data against wang_ref, with presence
  required.

### m2. The JSON record omits the switches that certification depends on (src/bin/xscan.rs:94, 156, 84-85)
- The record has `stream`, `compact` and `prune_n`, but not `check_d2`, `rung_transport`, `rungwise`
  or the engine commit.
- `--no-d2` (allowed for unpruned scans, line 156) also turns off the composite Π∘Ι check
  (`eng.check_orth = opts.check_d2`, 1751). A JSON from such a run looks the same as a certified one.
- A base-point FAIL is printed and stored as `false`, but the exit status is still 0.
- C_sound m5 makes `check_d2` the only exact guard against transport and cache bugs. The JSON files
  are the inputs for gl6_t2n.txt (M2S-run).

Fix:
- write `check_d2`, `rung_transport`, `rungwise`, `KRHOM_THREADS` and `git rev-parse HEAD` into the
  JSON;
- exit nonzero on a base-point FAIL;
- have the M2S-run ingestion refuse records with `check_d2 = false`.

## Notes (not defects)
- The lazy heap keeps stale duplicates. Its size is bounded by the total number of re-pushes, which the
  σ^9 run absorbed: 2.78 GB peak RSS in all.
- `check_degrees` decodes the whole partial complex on every pass. This only costs time.
- in_done is always true by construction. If a future change inserted d_C⊗id after a pass's `simplify`,
  the rule would still look satisfied. An assertion that every C component of slice s−1 is transported
  before pass s eliminates would make this structural invariant explicit.
