# MBOX report: box closure layer at fixed N (branch `colored-box`)

Date 2026-11-11. Branch `colored-box` from dev ce91f0e. Spec: SPEC_3C_GL6_rev1 §2 MBOX and §6 "Brief MBOX"
(D5, D7, D14), R3_closure §2-6, R4_m6a §2, C_sound M1, M2, M5. Files touched: new `src/boxclose.rs`,
`tests/boxclose.rs`, `tests/boxclose_mem.rs`, `examples/probe_final_d2.rs`; additive section at the end of
`src/cclose.rs` (the `ClosureDataC` struct is unchanged, no existing function changed); `src/lib.rs` (module line).
xengine, homprov, smap, tmpl, linsys and sindec are not edited.

## 1. API

### src/boxclose.rs
- **Nodes and points.**
  - `BoxNodes::new(lam, n, seed)`: n nodes per block, all distinct and nonzero across blocks.
  - `scaled(λ)` gives the node sequence λ·a, whose points are λ·z_μ. `hash()` is the checkpoints' node-seed hash.
  - `box_parts(a, n)`: partitions in a × (n−a), ordered by size.
  - `point_coords(nodes, μ)`: the values of e_{i,j} at z_μ. The roots of block i are a_{μ_j + λ_i + 1 − j}.
  - `eval_epoly`.
- **`BoxGrid::new(nodes, marked)`.** For each block it builds the factorial Schur matrix S^(i)[μ][ρ] = s_ρ(z_μ | a)
  (determinant formula) and its inverse.
  - The constructor asserts the vanishing property (S^(i) is triangular for containment, with a nonzero diagonal).
    It also asserts that the inverse is triangular.
  - Reduced mode: the marked block has only μ = ∅.
  - Point ids use a mixed radix. The grid lists points by size, and `supersets_of_size(μ, s)` gives the
    (S^{-1})_{ρμ} with ρ ⊇ μ and |ρ| = s.
  - `box_points(grid, max_size)` lists the down-set L(max_size).
- **`schur_to_epoly(grid, coords, N', marked)`.** Turns Schur coordinates into Σ c_ρ Π_i s_{ρ^i}, using dual
  Jacobi-Trudi in the e-variables of each block. It restricts to ρ^i ⊆ λ_i × (N' − λ_i) and, in reduced mode, to
  ρ^m = ∅. The image in A_λ at N' is the entry mod I_{N'}.
- **`BoxAcc`** (streaming, D14 / C_sound M2):
  - `new(grid, weights, library_hash, word_hash, prune_n)`. Weights are declared per entry. `None` declares an
    entry that must be 0 (negative or odd degree).
  - Per entry it holds the coefficients |ρ| = w (outputs) and |ρ| = w + 1 (residuals).
  - `add_point(id, values)` refuses a point that was already added, a wrong number of values, and a nonzero value
    on an entry declared 0.
  - `finish()` refuses an incomplete down-set. The check covers all points with |μ| ≤ min(w_max + 1, top), which
    is exactly the union of the down-sets of every output and residual ρ. It also refuses any nonzero residual.
  - `merge` refuses a mismatched header (every field) and overlapping done-bitmaps.
  - `save` is atomic: it writes a temp file in the same directory, fsyncs it, renames it, then fsyncs the
    directory. The text header has prime, node_seed_hash, library_hash, word_hash, prune_n, marked, n, lam,
    code_stamp, entries_hash, entries, npoints, ndone, payload_len and payload_hash. The binary payload holds the
    bitmap, the weights and the accumulators.
  - `load(path, grid)` refuses a truncated or corrupt payload (length or hash), and a prime, nodes, N, λ, mode or
    code stamp that differs from the grid.
- **`check_scaling(v_z, v_λz, λ, w)` and `check_scaling_all`** (C_sound M1). Each returns Err on v(λz) ≠ λ^w v(z).
  The report counts the entries that are 0 at both points.
- **Closure entries.** `closure_entries(skeleton, edges)` lists every Gram entry, then every P entry, with the
  declared weights from the bookkeeping. Gram: 2w = ins + outs. P: 2w = q_s − q_t + ins_s + outs_t.
  `eval_entries` evaluates exact data; it stands in for a scan in tests and the 2-strand regression.
- **`BoxClosureData`** (D5):
  - fields: `{ base: ClosureDataC (types and objects; Gram from the box; edges empty), prune_n, n_box, valid_n =
    min(prune_n, n_box), mode, pairing, d2_cert: BTreeMap<(N', marked), D2Cert>, stamp }`;
  - `mode` is an addition to the spec's list. Data extracted in reduced mode (points with μ^0 = ∅ only) serve only
    that reduced closure; unreduced data serve both modes;
  - constructors: `from_box(...)` assembles the data from a `BoxResult`; `new(...)`.
- **`close_box(data, N', marked, t0, cap)`.** It refuses:
  - N' > valid_n;
  - a mode the data cannot serve;
  - pruned data without a `D2Cert` for (N', marked) that matches the current pairing hash and library hash;
  - a T0 record of another (N', mode, library hash).

  It then runs `cclose::try_close_staged` on P, capped by the MOY ranks when a T0 record is given.
- **`closed_d2_freivalds(&mut data, N', marked, k, seed)`.** Refuses N' > valid_n. It runs the Freivalds test and,
  on success, records `D2Cert { n, marked, vectors, blocks, seed, pairing_hash, library_hash }`.
- **`box_moy_check(data, N', marked, moy, cap)`** (T0). It computes the uncapped staged graded dims of every type
  and compares them with the MOY values, then returns `T0Record { n, marked, library_hash, ranks }`. That record is
  the only route to the rank cap in `close_box`.
  - `moy_type_trace(lam, table, N)` gives tr_N(b_p) through `schur::decat_chi`.
  - The reduced MOY value is tr_N(b_p) / [N choose λ_0]; test (b) checks this.
- **Hashes.** `library_hash`, `pairing_hash`, `word_hash` and `entries_hash` use FNV-1a, which is stable across
  builds. `BOX_CODE_STAMP = "boxclose-v1"`. `PrimeField` (u64 serialization) is implemented for `Fp<P>`.

### src/cclose.rs (additive)
- **Pairing types and helpers.**
  - `Pairing<F>`: the closed pairing per edge, P[b][a].
  - `pairing_from_f(data)`: P = G_t·F, the oracle's `pab`.
  - `validate_pairing`.
- **`try_close_staged(data, pairing, N, marked, &StagedOpts { cap, rank_cap })`.** The rank stage on P:
  - State spaces are built per type, only for types that carry objects, one type per thread.
  - Rows m·G_{·a} are streamed into a `RankEchelon` and only the independent rows are kept.
  - With `rank_cap`, selection stops at r_d, and it is an Err if r_d is not reached.
  - Columns are indexed per degree without a hash map of discovered columns, in `order_key` order.
  - Without a cap, the representatives are the oracle's (greedy selection in the same row order).
- **`try_staged_type_ranks`**: graded dims per type and degree.
- **`staged_mem_estimate`.** This is the refit of `exact_mem_estimate` for the staged path. It sums:
  - the A_λ tables;
  - the k largest per-type builds (k = `KRHOM_THREADS`; `staged_type_mem_estimate`);
  - 64 B per bound on chain coordinates;
  - the k largest rank-stage echelons;
  - 48 B per pairing term;
  - 16 MB fixed.

  Echelons and kept rows are charged as dense: 64 B per rank × column entry, the worst case of fill-in. The
  state-space bound per type and degree is min(#candidate rows, #columns), or the cap. `exact_mem_estimate` itself
  is unchanged, so existing behaviour does not change.
- **`closed_d2_freivalds_raw(data, pairing, N, marked, opts, k, seed)`.** For k random u per (h, q) block it
  computes y = D_h u on every A_λ coordinate and recovers x = Π^{-T} y. Π is the rep × dual matrix per type and
  degree (dense LU, blocks up to `FREIVALDS_MAX_BLOCK` = 4096).
  - It asserts that Σ x_i·rep_i = y on every coordinate, i.e. the image lies in Im(G_t). This catches P entries
    that are inconsistent with the Gram data.
  - It asserts D_{h+1} x = 0 on every coordinate.
  - The false-pass probability is ≤ p^-k per block.

## 2. Tests

`KRHOM_THREADS=4 cargo test --release --test boxclose` takes 2 s for 13 default tests. `--test boxclose_mem` takes
8 s.

| test | what | result |
|---|---|---|
| (a) `a_synthetic_extraction_equals_nf` | sparse random homogeneous polynomials (2 per weight, 2-5 terms), weights 0..top+1; (3,3), (2,2,2), (3,3,3); extraction at N = 4, 5, 6; both modes; F31 and F32b; compared with `GrAlg::nf` at every N' from max λ to N (nesting) | 3 920 entries equal, 2.0 s |
| (a) `a_grid_basics` | point counts (8 000 / 400, L(10) = 1 834 / 269, max box 777 / 48, as in R4 §2 and R3 §2), node distinctness, point coordinates, hashes | ok |
| (b) `b_end_to_end_22_s3` | (2,2) σ^3, F31 and F32b, N = 3..6, unreduced and reduced, unreduced and reduced extraction: Freivalds certificate, scaling check at 2 points, T0, `close_box` uncapped and capped = `close_exact_colored` | ok (ranks 7 / 22 / 46 / 79 unreduced, = Wang T(2,3)) |
| (b) `b_end_to_end_32_s2` | (3,2) σ^2, same checks, F31 | ok (rank 300 at N = 6) |
| (b) `b_end_to_end_32_s2_f32b` `#[ignore]` | the same over F32b | ok (run in this session) |
| (b) `b_end_to_end_33` `#[ignore]` | (3,3) σ^1, σ^2, σ^3 F31; σ^3 F32b | ok, 184 s: N = 6 ranks 20 / 400 (Hopf) / 136 (trefoil, both primes) |
| `staged_equals_oracle_on_exact_data` | `try_close_staged(pairing_from_f)` = `close_exact_colored`, staged dims = `exact_type_dims_colored`, Freivalds passes: (2,2) σσ^{-1}σ, (2,1) σ^2, (1,2) σ^4, N = 2..5, both modes | ok |
| `moy_convention_matches_single_type_closure` | the single-object closure of each final type equals `moy_type_trace` | ok |
| (c) `c_freivalds_detects_perturbed_p_entry` | P entry + 7·e_{1,1}^w on every eligible edge | 10 / 10 detected; the stale certificate is refused by `close_box` |
| (c) `c_freivalds_detects_reinserted_d2_component` | final d∘d component of the (3,3) σσσσ^{-1}σ prune-4 scan, re-inserted (see below) | passes on f (closes to 0); fails on f + 5·e_{1,1}^k; equal to `closure_kills_map`, raw and through the box |
| (c) `c_point_added_twice_overlap_header_downset` | double add; overlapping shards; mismatched header (library, word, prune_n, node seed, mode); checkpoint read against another node sequence; truncated checkpoint; incomplete down-set; an entry declared at w − 1 (nonzero residual); a nonzero value on an entry declared 0 | each Err as intended |
| (c) `c_entry_declared_at_w_plus_1_fails_scaling` | one nonzero entry declared at w + 1 | `check_scaling_all` Err. The extraction alone finishes (silent zeros), which shows that the scaling check is needed. |
| (c) `c_valid_n_and_d2_cert_guards` | prune-4 scan closed at N' = 5; pruned data without a D2Cert; a certificate for another (N', mode); box N < N'; T0 against wrong MOY values; a T0 record of another library | each Err as intended; with the certificate, the result = oracle |
| (d) `d_shuffled_order_and_shards_agree` | (3,3) N = 6, (2,2,2) N = 5 reduced, (3,3,3) N = 5: shuffled order; 3 shards saved, loaded and merged; merging a shard twice | identical; Err |
| `extraction_cost_per_entry` | §3 | prints |
| (e) `e_staged_rank_stage_41_size_reduced_within_estimate` | §4 | peak within estimate |
| (e) `..._unreduced_refused_before_allocating` `#[ignore]` | §4 | refuses with Err, peak 0 MB |

**Negative test "pruned scan with a final_d2 edge re-inserted" (interpretation).**
- A probe (`examples/probe_final_d2.rs`) found no 2-strand scan with labels ≤ 2, nor (3,1) or short (3,2), that
  has final d∘d components.
- The (3,3) scans pruned at 4 do have them. Every final component of every (3,3) word of length ≤ 8 at prune 4
  closes to 0 at N' = 4, 5 *and* 6 (`closure_kills_map`). So a genuine re-inserted component can never make
  Freivalds fail.
- The test therefore re-inserts the component f as an appended complex x' → w' → w'' (edges f and id). The closed
  D² of that complex is exactly Φ(f). The test checks:
  - Freivalds passes on f;
  - it fails on f + c·e_{1,1}^k, which simulates a pruning bug that leaves a component outside the killed ideal;
  - both verdicts equal `closure_kills_map`, for the raw call and for `closed_d2_freivalds` + `close_box` on box
    data.

**Full suite.** See §6.

## 3. Extraction cost per entry, against R3 §3

Measured by `extraction_cost_per_entry`: N = 6, F31, one thread. Mult-adds count outputs (|ρ| = w) **and**
residuals (|ρ| = w + 1). R3's streaming column counts outputs only, so the comparable figure is about half of mine.

| case | points | max mult-adds per entry (w) | mean over weights | R3 §3 streaming (outputs only) | accumulator memory per entry, max (outputs + residuals, 8 B) | R3 (outputs only) | measured time per entry (incl. evaluation at all points) |
|---|---|---|---|---|---|---|---|
| 2 unmarked blocks: (3,3) unreduced, (3,3,3) reduced | 400 | 8 268 (12) | 3 224 | ≤ 4.2e3 | 744 B | ≤ 384 B | 0.20-0.24 ms |
| (3,3,3) unreduced | 8 000 | 1 195 897 (18) | 382 812 | ≤ 6.1e5 | 12.4 KB | ≤ 6.2 KB | 16.8 ms |
| (3,3) reduced | 20 | 66 (6) | 35 | — | 48 B | — | 0.01 ms |

So the streaming cost matches R3 within the factor of about 2 that the residuals add.
- At R3's 4_1 estimate of 2.6e5 P entries, unreduced extraction would take about 2.6e5 × 17 ms ≈ 1.2 h of CPU per
  prime. Reduced would take about 1 min. Both exclude the scans.
- The per-block sweep (R3: ≤ 2.1e5 per entry) is not implemented. It would cut the unreduced cost by about 3-6×,
  at 10× the memory.

In test (b), the (3,3) σ^3 extraction (470 entries, 392 points, 35 394 coefficients) takes 0.2 s unreduced and
0.01 s reduced.

## 4. Memory of the staged rank stage at 4_1 size

Synthetic data (`tests/boxclose_mem.rs`):
- types: λ = (3,3,3), the 4_1 = σ_1σ_2^{-1}σ_1σ_2^{-1} types of `decat_scan` kept at N = 6. That is **34 types**,
  agreeing with R3's 34;
- ins: the exponents of c_diag(b_p), 281 in all, max n_p = 54;
- objects: from the K_0 class, **276**, agreeing with R3/R4's 276 final objects;
- Gram and P: random sparse homogeneous entries (3, resp. 2 terms) of the bookkeeping weights; 2 045 pairing edges.

The peak is measured with a counting global allocator (bytes requested, all threads), with `KRHOM_THREADS` = 4.

| case | chain dim | measured peak | `staged_mem_estimate` | oracle `exact_mem_estimate` (same data) | time |
|---|---|---|---|---|---|
| reduced N = 6 (dim A = 400), uncapped | 62 534 | **133.0 MB** | 3 399 MB (ratio 0.04) | 2 095 MB | 8 s with the capped rerun |
| reduced N = 6, rank cap = uncapped ranks | 62 534 | 133.4 MB | 2 492 MB | — | |
| unreduced N = 6 (dim A = 8 000), uncapped | — | 0 MB (refused) | 2.6e6 MB | — | Err before any allocation |

- The estimate is a safe bound: every measured case is within it. It is loose, because echelon fill is charged as
  dense; real fill is 4% of that.
- Random Gram data give generic ranks, so the rank cap does not shrink the synthetic case. On real data, r_d =
  dim Φ(B_p)_d is far below min(rows, cols).
- Unreduced 3-strand N = 6 needs the cap from a real T0 plus a sparse fill model before it can run (§7).

## 5. Commits

See the final answer of this session (git log on `colored-box`).

## 6. Full suite

`KRHOM_THREADS=4 cargo test --release`: result recorded at the end of the session (§8).

## 7. Open problems

1. **The unreduced staged estimate is unusable at 3 strands, N = 6** (2.6 TB with dense charges).
   - What is needed: MOY caps from T0 on the real (3,3,3) library, plus a measured fill model for `RankEchelon` on
     m·G rows to replace the dense charge.
   - Calibration data so far: reduced at 4_1 size, the peak was 4% of the model.
   - The duals stage also keeps r selected rows of ≤ C columns.
2. **The Freivalds Π^{-1} is a dense LU per (type, degree)**, refused above r = 4096. That covers reduced 3-strand
   cases (≤ 54 × 48 = 2 592 rows) but not unreduced ones (r up to about 42 000). A sparse LU, or the capped Π, is
   needed there.
3. **The deterministic P2 variant (R3 §5.3) is not implemented.** The generator-only determinant pass, frozen
   structure and the `box_run` / `box_finish` drivers belong to M6a-2a.
4. **The residual check is in place, but it does not cover entries of weight ≥ top** (O2 of R3): unreduced
   entries of weight ≥ 27 and reduced ones of weight ≥ 18. There, only the scaling check, two node sequences and
   two primes apply.
   - The scaling check here runs on exact data at 2 points. A production run must apply it to scan values at z and
     λz.
5. **The checkpoint header's library hash and word hash are supplied by the caller.**
   - `BoxClosureData::from_box` recomputes `library_hash` from the extracted Gram data. T0 records and D² certificates
     are bound to that recomputed hash and to the pairing hash.
   - The checkpoint's library hash should be the M5 layer-A hash once that exists.
6. **Re-inserted d∘d test.** No genuine final component fails to close at N' ≤ 6 among the (3,3) prune-4 words of
   length ≤ 8, so the failing branch of the test uses a synthetic non-ideal term (§2).
7. **The R3 sweep alternative** (per-block S^{-1}) for unreduced extraction is not implemented (§3).

## 8. Full-suite result
(filled below)
