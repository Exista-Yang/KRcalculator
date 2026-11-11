# Review colored-2s (0fb5dc4 vs ce91f0e): closures and outputs

Scope: `close_stages`, `close_stage` / `close_stage_marked`, `close_guarded`, the xscan outputs (JSON, `--chi-only`), the
docs M2S_G2X.md / M2S_step0.md / PLAN.md, the Wang blocks in tests/data/wang_ref.txt, and the tests/ diff.
Read-only review; no cargo. Checks run: Python on the committed JSONs, plus tools/wang.py for blocks that finish in
<= 15 s.

## Verdict
No defect in the closure path can produce a wrong Poincaré polynomial. Four minor defects are in the test and
reporting layer.

## Verified (no defect)
- **Stage closures.** In src/xengine.rs:1856-1873 a stage closure is built from `ncx` after that stage's d∘d, K_0
  and LB checks (1808-1840). It uses `crossings.clone()`, which holds exactly the first m crossings (pushed at 1765),
  and `last_d2` of the same stage (1819). It runs only when `labs == labels` and 2 <= m < len; the final stage
  goes through `data`.
- **Guards.** `close_stage_marked` (1956-1974) refuses N' > prune_n for both final and intermediate stages.
  `close_guarded` (1977-1995) asserts N' <= prune_n and runs `closure_kills_map` on every stage-m d∘d component with
  the same `marked`. These are the same semantics as `close` at ce91f0e. The stage polynomials printed during the
  scan come after that stage's checks, and they also go through `close_guarded`.
- **global_q and grading.** A wrong crossing list would shift q by multiples of 9. Stages 3/5/7/9 equal Wang
  bigraded, so this is confirmed.
- **JSON (6 files).**
  - All 96 closures agree across sigma5 F31 / F31_parallel / F32b, sigma7 F32b and sigma9 F31 / F32b.
  - Stage 2 equals `hopf_formula` at N=3..8.
  - Odd stages equal Wang 3 N m wherever a block exists. At N = 3 they give t^0q^0, which is what tools/wang.py
    gives for 3 3 2/3/4.
  - All rank and per-stage tables in M2S_G2X.md match the JSON. Base point is true at 3, 5, 7 and 9.
  - χ(unred) = [N,3]·χ(red) holds for knots and links.
- **Independent χ check (Rosso–Jones).** The χ of all stages 2..9 at N=4..8, the links T(2,4/6/8) included, fits
  Σ_j d_j (ε_j q^{c_j} f)^m, where c_j are the contents of the ∧^3⊗∧^3 summands. The fit uses 4 stages; the other
  4 are predicted, with alternating ε and f = ±q^{-3(N-4)}.
- **Large-N and N = 4 behaviour.**
  - Reduced knot ranks stabilize at 3^3, 5^3, 7^3 and 9^3 for N ≥ 6.
  - At N = 4 (∧^3 ≅ dual), the ranks match the uncolored N+2k(N−1) for knots and N^2+2(k−1)(N−1) for links.
- **Wang blocks.**
  - Format: 39 blocks, 4 lines each, rank = Σ coefficients, chi line = Euler characteristic, terms sorted.
  - Regenerated with tools/wang.py and byte-identical: 3 4 3, 3 4 4, 3 5 3, 3 5 4, 3 7 1, 3 8 1.
  - Not regenerated (each takes over 20 s): 3 8 2, 3 7 4, 3 8 3, 3 8 4. They equal the two-prime engine output.
- **Tests.**
  - tests/xengine.rs only moved `parse_poincare`, `wang_ref` and `hopf_formula` to tests/common/mod.rs. No test
    was deleted, weakened or newly `#[ignore]`d.
  - XOpts defaults are now stream/compact/rung_transport, so the existing tests exercise the new path.
  - `close_stage_equals_separate_scan` and `close_stage_refuses_n_above_prune_n` test what they claim.

## Findings

### F1 (minor, test) g2x_sigma5_33 silently skips missing Wang blocks; stage 5 at N = 3 is never compared
- **Where.** tests/xstream.rs:227 `if let Some(want) = wang_ref().get(&(3, nn, mm))`. The docstring at 199-200
  claims stage 5 = Wang 3 N 2 "at N = 3..6".
- **Gap.** wang_ref.txt has no block 3 3 2, so N = 3 is skipped with no message. Deleting or mis-keying any block
  also turns that check into a no-op. The same pattern would apply to T(2,7) and T(2,9) at N = 3 (no 3 3 3 or
  3 3 4 blocks).
- **Fix.**
  - Append the blocks 3 3 2, 3 3 3 and 3 3 4. In review, tools/wang.py gives t^0q^0 / rank 1 / chi +1q^0 for each.
  - Make the test `expect` a block for every N in the range it claims.

### F2 (minor, output provenance) The xscan JSON does not record check_d2, so an uncertified run looks certified
- **Where.** src/bin/xscan.rs:94-95 records field, labels, word, prune_n, stream and compact. It does not record
  check_d2, rungwise/rung_transport, seed, threads or the git rev.
- **Why it matters.** `--no-d2` disables the exact d∘d = 0 assertion. It also disables the Π∘Ι = id composite
  check, because xengine.rs:1751 sets `check_orth = check_d2`.
  - C_sound m5 says this is "the only exact guard against transport and cache bugs".
  - The JSONs are the deliverable from which M2S-run writes tests/data/gl6_t2n.txt.
- **Fix.** Write `"check_d2"` and the remaining options to the JSON. Refuse `--json` with `--no-d2`, or mark that
  output uncertified.

### F3 (minor, crash / silent omission in a check tool) `--chi-only` catches the wrong panic
- **Where.** src/bin/xscan.rs:143-148.
- **Problem.** The code wraps `decat_scan` in `catch_unwind` with the comment "colors do not close up". But
  `decat_scan` does not panic on colors that fail to close (schur.rs:1840-1858, 1917). It is `decat_chi` that
  asserts `labels == lam` (schur.rs:1997), and that call is outside the catch.
- **Consequences.**
  - Any (a,b) word with a ≠ b aborts the tool at m = 3.
  - A genuine panic in `decat_scan` is swallowed. The χ lines for that m are then silently missing, and a
    comparison over the printed lines passes.
  - The doc's "χ = decat_chi for all 48 unreduced closures" for sigma9 depends on this path. The count of 48 is
    right there, and an independent χ check confirms it.
- **Fix.** Test closability from the label permutation before calling either function, and do not catch panics.

### F4 (minor, gate measurement) "peak RSS" is the peak since the start of the last crossing, not the process peak
- **Where.** xengine.rs:1769 calls `reset_hwm()` at every crossing. xscan.rs:90 and :95 report `rss_kb().1` as
  the run's peak RSS. The memory gates (≤ 6 GB, ≤ 13 GB) read that value.
- **Effect.**
  - For the committed JSONs the per-stage peaks increase from stage to stage, so the reported value is the true
    peak (checked).
  - A word whose largest crossing is not the last one is under-reported.
  - On macOS /proc is missing: 0 is reported and the reset fails silently.
- **Fix.** Report the maximum of `stages[].mem.rss_peak_kb` and the post-closure HWM, or skip the reset.

## Notes (not defects)
- **Base-point check is weak.** On σ_1^m the check (marked 1 = marked 0) follows already from the left-right
  flip symmetry of the braid. It tests that cclose handles different marked blocks consistently, not base-point
  independence on a nontrivial braid.
- **Stale source path.** The wang_ref.txt header and tests/common/mod.rs:29 name `colored_design/r3scripts/wang.py`,
  which does not exist. The generator is tools/wang.py.
