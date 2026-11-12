# Review wave 3: box closure (revm6a c95517f)

Scope: `src/boxclose.rs` (the whole file), the additive MBOX section of `src/cclose.rs` (lines 2363-3044), `tests/boxclose.rs`, `tests/boxclose_mem.rs`, and the `boxclose` call sites in `src/swire.rs` and `src/bin/sbox.rs`. The review was read-only: nothing was compiled or run. The only computation was a check, under one second, that the scaling λ has large multiplicative order in both fields.

**Verdict.** I found no critical or major defect in this lens.

These parts are mathematically correct and are guarded at runtime:
- the factorial Schur extraction;
- restriction to N' and to reduced mode;
- the Freivalds test;
- the staged rank stage.

There are two minor findings. Both concern guards that are weaker than the report and the docs claim.

## Verified correct (with evidence)

### Factorial Schur extraction
- **Point roots.** At z_μ the roots of each block are a_{μ_j+a+1−j}. This is consistent in three places: `block_e_values` (boxclose.rs:243), `BoxGrid::new` (:432) and `swire::box_x` (swire.rs:662).
- **Factorial Schur.** The factorial powers are (x|a)^k = Π_{t≤k}(x−a_t), so the vanishing points are x_j = a_{μ_j+a−j+1}. That is the Macdonald/Molev convention.
- **Vanishing and inverse.** The constructor checks, for every block:
  - vanishing (s_ρ(z_μ|a) = 0 unless ρ ⊆ μ) and a nonzero diagonal, over the whole box (:436-444);
  - that the forward-substitution inverse is triangular for containment (:465-473).

  A wrong convention would therefore panic; it could not corrupt data silently.
- **Ordering.** Sorting by (size, lex) (`box_parts` :235) is a linear extension of ⊆, so the forward substitution in :447-462 is valid.
- **Components outside the box.** σ ⊄ box(N) never contributes at box points: σ_1 > N−a ≥ μ_1. So the in-box coefficients come out correct whatever the nodes past a_N are.
- **Top-degree coefficients.** For homogeneous f of weight w, the factorial coefficients satisfy c'_ρ = c_ρ (the ordinary Schur coefficient) for |ρ| = w, and c'_ρ = 0 for |ρ| > w.
- **Tensor structure.** S = ⊗S^(i). `supersets_of_size` (:534) enumerates (S^{-1})_{ρμ} = Π_i (S^(i))^{-1}[ρ^i][μ^i] correctly, and the `rem` pruning is sound.
- **Down-set completeness.** Two facts make the completeness check exact:
  - `add_point` skips μ with |μ| > w+1 (:718). That is correct because ρ ⊇ μ forces |ρ| ≥ |μ|.
  - `finish` (:752) requires every point with |μ| ≤ max_w min(w+1, top). That set is exactly the union of the down-sets of all output and residual ρ: the box poset is graded, so every μ below top extends to some ρ.
- **Restriction to N'.** `schur_to_epoly` (:569) keeps ρ_1 ≤ N'−a per block, and returns 0 if a label exceeds N'. This is the class in H*(Gr(a,N')) because c_ρ does not depend on N.
- **Reduced mode.** The marked block has only μ = ∅, and its roots a_1..a_λm are nonzero. Only ρ^m = ∅ survives, which matches `GrAlg`'s trivial marked block (it kills every monomial in e_m, cclose.rs:949).

### Dual Jacobi-Trudi
`jacobi_trudi` (:322) computes det[e_{ρ'_r − r + c}] with the correct cofactor signs. Test (a) also checks it against `GrAlg::nf` at every N' (nesting), for both primes and both modes.

### BoxAcc checkpoints
- **What the header covers.** The header has prime, node hash, library hash, word hash, prune_n, mode, N, λ, code stamp and entries hash. The payload hash and payload length cover the bitmap, the weights and the accumulators.
- **`load` (:864) refuses** each of the following:
  - a truncated or corrupt payload, or a missing header end;
  - a mismatch with the grid in prime, nodes, N, λ, mode or npoints;
  - a different code stamp;
  - a mismatched entries hash or payload size;
  - a done count that differs from the bitmap popcount.
- **`merge` (:782) refuses** any header difference and any overlap between the done bitmaps.
- **`save`** is atomic (temp file, fsync, rename, directory fsync).

### Scaling check
- `check_scaling` (:945) tests v(λz) = λ^w v(z).
- The scaled nodes λa give roots λx, and each e_j scales by λ^j. This is consistent with `scaling_check` (swire.rs:831).
- λ = 0x12345677 has no order ≤ 200 in F_{2^31−1} or in F_4294967291 (checked), so weight differences of practical size are visible.
- The two scaling points are box points made from random nodes. They are effectively random points, so "zero at both points" means the polynomial is identically zero, with high probability.

### close_box guards
- **valid_n.** `BoxClosureData::new` (:1093) sets valid_n = min(prune_n, n_box).
- **check_request** (:1128) refuses N' > valid_n, and refuses reduced-mode data asked to serve another mode. `close_box`, `closed_d2_freivalds` and `box_moy_check` all call it.
- **Pruned data** (:1158-1163) needs a D2Cert stored under the key (N', marked), whose pairing hash and library hash are recomputed from the current data.
- **T0** must match n, marked and the recomputed library hash (:1165-1170).
- **The rank cap is exact.** It equals the uncapped ranks of the same library, and greedy selection stops at the same rows, so capped = uncapped.

### Freivalds test (cclose.rs:2921)
- u is uniform over F_p and drawn per (h, q) block. `k.max(1)` vectors are used, so k = 0 does not pass vacuously, and sbox passes k = 3.
- y = D_h u is computed on every A_λ coordinate.
- Then x = Π^{-T} y|duals, followed by two exact checks:
  - Σ x_i rep_i = y on every coordinate. This is membership in Im(G_t), so the projection the rank stage uses is injective on the image.
  - D_{h+1} x = 0 on every coordinate.
- A false pass has probability ≤ p^{-k} per block for each of the two linear conditions.
- Blocks with no (h+1, q) group are skipped. This is correct for D²: the chain group there is 0, and the rank stage assigns rank 0.

### Staged rank stage = oracle
- **Row selection.** `type_state_space_staged` (:2520) runs the same rows (a, m1), in the same order, through a greedy `RankEchelon` insertion. So when uncapped, the representatives are identical to `type_state_space` (:1079).
- **Duals.** These come from a different column order (order_key instead of discovery order). Any r independent columns give an injective projection of Φ(B_t)_d, so the ranks are unchanged.
- **Rest of the stage.** `closed_cx` builds the same groups, dual ordering and global shift as `try_close_exact_colored_capped`. The rank loop (:2822-2846) is the oracle's loop (:1360-1385), with `pairing[ei].2[b][a]` in place of `pab[ei][a][b]`.
- **Missing type in a rank cap.** It yields zero reps (`EMPTY_DEG_RANKS`). This is harmless through `close_box`, because the T0 ranks list every type of the hashed library.

### Negative tests are real
For each test I checked that the code path it targets is the only possible source of the failure (where the assertion does not check the message):
- the double add;
- each kind of header mismatch;
- a checkpoint loaded against other nodes;
- truncation;
- an incomplete down-set;
- an underestimated weight caught by the residual;
- a nonzero value on an entry declared 0: removing the `None` guard would let `finish` succeed, so the test would fail;
- overlapping shards and merging a shard twice;
- the stale certificate after a P perturbation: only the hash mismatch can make it Err;
- valid_n with prune 4 and N' = 5;
- a missing certificate at another (N', mode);
- T0 against wrong MOY values;
- a T0 record from a foreign library.

The two Freivalds tests cover the two different branches:
- `c_freivalds_detects_perturbed_p_entry` mainly exercises the Im(G) consistency branch.
- `c_freivalds_detects_reinserted_d2_component` exercises the genuine D² ≠ 0 branch: there P = G·F' is consistent and D² = Φ(f').

## Findings

### F1 (minor). The scaling check is opt-in, and the residual check misses underestimated weights in reduced mode
- **Docs vs. code.** The module doc (boxclose.rs:19-20) and MBOX_report §1 say the w+1 residual "catches an underestimated weight", and that only an overestimate needs the scaling check.
- **Where the residual is blind.** In reduced mode the residual sees only ρ^m = ∅. Counterexample: λ = (2,2), block 0 marked, and f = e_{0,1}·g with g of weight w in block 1, declared at weight w (true weight w+1).
  - At reduced points, e_{0,1} = a_1 + a_2 =: c ≠ 0, so the extraction returns c·g at weight w.
  - The residual at w+1 with ρ^0 = ∅ is 0, because every Schur component of f has ρ^0 = (1). So `finish` succeeds.
  - The true reduced value is f|_{e_0=0} = 0.

  The same blind spot appears in any mode when f's weight-(w+1) Schur components all lie outside box(N_box), and for w ≥ top (already noted as O4).
- **Who catches it.** Only `check_scaling` does.
- **Why that is a gap.** Nothing enforces the scaling check:
  - `BoxClosureData` carries no scaling record;
  - `close_box` does not ask for one;
  - `swire::box_word` (swire.rs:730-829) never runs it;
  - sbox runs it only with `--scaling` (sbox.rs:80-86).

  The G2B and stretch runs did pass `--scaling` (their JSON files have `"scaling":{...,"pass":true}`), so the published numbers are covered. A future driver, or a run without the flag, could still produce wrong reduced Poincaré polynomials from a bookkeeping error, and no guard would fail.
- **Fix.** Add a `ScalingCert` to `BoxClosureData`: node hash, λ, point ids, word and library hash, checked count. Have `close_box` refuse data without one, as it does for D2Cert. Correct the doc so it says the residual is sufficient only in unreduced mode, for in-box components with w+1 ≤ top.

### F2 (minor). The checkpoint header binds the entry list only by its weights
- `entries_hash` (boxclose.rs:616) hashes the weight sequence only. It does not hash the `EntryKind` list (type, b, a / edge, b, a), the skeleton (objs, edges), the stage index, or the engine version. `BOX_CODE_STAMP` (:54) versions boxclose only.
- **Failure scenario.** Two shards with equal weight lists but different entry meanings merge without complaint (`merge` :782). Examples:
  - two closable stages of the same word;
  - replays under a different `PivotOrder` or a changed sengine. The frozen pivot trace changes the minimal-complex basis, and with it P, while every header field stays the same.
- **Consequences.**
  - For entries with w+1 ≤ top, the mixed interpolation usually leaves a nonzero residual, so `finish` fails loudly.
  - For entries with w ≥ top (no residual) the top coefficient is silently a mixture.
  - In reduced mode, the F1 blind spot applies as well.
- `load(path, grid)` also cannot compare library, word or prune_n against what the caller expects. A single-shard resume from a stale file of another run with the same weight list is accepted.
- **Impact today: none.** `box_word` does not checkpoint (M6a1b report §6.4). It matters for the M6a-2a box_run/box_finish drivers.
- **Fix.** Include in the header:
  - a hash of (kinds, skel.objs, edges, stage);
  - an engine or library code stamp (sengine/swire version, pivot order);
  - the record-point seed.

  Give `load` an expected header and compare all of it.

## Notes (not defects in this lens)
- **Memory in `box_word`.** It buffers every point's values for every stage (swire.rs:764-798) before calling `add_point`, which defeats the streaming design of `BoxAcc`. Memory is about points × entries × 8 B:
  - about 0.8 GB at the R3 4_1 estimate, reduced (400 points × 2.6e5 entries);
  - about 17 GB unreduced (8 000 points), though that closure is refused by O1 in any case.

  This is worth fixing in M6a-2a.
- **D2Cert binding.** `D2Cert.vectors` records k rather than `k.max(1)`, and `close_box` does not check `vectors`. The pairing and library hashes do not cover `base.objs`. Both are harmless for the current drivers.
- **Skipped Freivalds blocks.** For blocks whose (h+1, q) group is empty, P is not checked for membership in Im(G_t). That does not affect the closure, but P corrupted only there would go unnoticed.
