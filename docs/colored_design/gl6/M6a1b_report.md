# M6a-1b report: library wiring and the 2-strand box gate (branch colored-m6a)

Date 2026-11-12. Machine: cloud container, 4 cores (Xeon @ 2.10 GHz), 15 GB. Spec: SPEC_3C_GL6 rev1 §0,
§1 (D5, D7, D14), §2 "M6a-1b", §6 "Brief M6a-1b". Results over F_p for p = 2^31−1 and 4294967291; equal
to the result over Q unless both primes are torsion primes.

## 0. Verdict

| item | result |
|---|---|
| wiring (`dec_word`, `rickard_comp_at`, d_C ⊗ id per coset, K_0 / LB / degree asserts, `check_d2`, λz scans, box loop) | **done** (src/swire.rs, src/bin/sbox.rs) |
| tests: sengine vs xengine at 2 points ((2,2), (3,2) default; (3,3) ignored), (1,1,1) vs engine2, box loop vs exact closure | **pass** |
| **G2B as specified (σ^1..σ^9 at N = 6 by box)** | **not reached: σ^1..σ^5 PASS; σ^6..σ^9 not run unreduced** (cost, §4) |
| G2B σ^1..σ^5 unreduced + reduced, N' = 3..6, both primes × two node sequences, scaling | **PASS** (Hopf 400, T(2,3) 136, T(2,4) 864, T(2,5) 396 and all reduced, full polynomials = m2s / Wang) |
| σ^7 reduced (20 box points) | __S7__ |
| 4_1 type counts (C_feas m9) | **34** types in the final LB complex at N = 6 (276 objects); 40 over all stages; "38" is none of the counts (§5) |
| stretch: (3,3,3) prune 4 at N = 4 through the localized path | **PASS**: 4_1 = engine2 gl_4 via duality (expected); Markov words σ_1^m σ_2 (m = 1..5) = 2-strand σ^m exactly; both primes × two node sequences; χ; scaling |

The localized path is correct wherever it was run (every closure equal to the exact references,
term by term). It is not yet fast enough for the full gate: per-point cost grows about 4× per
crossing on (3,3) (σ^5: 9.8 s and 130 MB per point; σ^6 record 65-87 s and 0.4-1 GB; σ^8 record
> 2 GB), so σ^9 at 400 box points would take on the order of 400 × 1 h per prime and node sequence.
§4 has the measurements and §7 the open problems.

## 1. What was built

**src/swire.rs** (new; no frozen file touched; `src/sengine.rs` got performance work, §4):
- `WLib` / `build_wlib(SIndec, words, closure)`: breadth-first walk over every type of every
  unsimplified stage of the words (every kept term at the prune level, every Rickard component
  of adjacent kept terms, closure generators of every type reached at a stage with labels λ),
  then `freeze`. Library-side release asserts, once: Rickard components run between the canonical
  ladders with degree t_src − t_tgt; every decomposition summand has ι: L_q → X of degree s and
  π: X → L_q of degree −s (`MapInfo`); closure generator degrees = exponents of c_diag(b_p)
  (as xengine). Shapes from `lib_shape` with rep(p) taken on the builder's Schur; table names;
  the identity type and coset; the library hash. A types without closure generators is an error
  only if it is final.
- `scan_point(lib, word, x, id, Mode, ScanOpts)`: one scan at one point (record or replay).
  Per crossing: terms ascending (positive) or descending (negative), one engine step per term;
  objects (q_j, h + hdeg_k, t + dq_k + s_j) from `dec_word`; **d_C ⊗ id** per coset as
  π_τ ∘ (f ⊗ id_W) ∘ ι_σ with the Kron rule applied to f's localized blocks (a cached,
  point-independent flow map L_a ⊗ W → L_b ⊗ W per coset; sparse in f); **(−1)^h id ⊗ d** from
  `rickard_comp_at`, evaluated by `Eval` on the engine's layouts. Asserts after every crossing
  (record point): **K_0 = the decategorified scan** (own `Decat` in its own Schur, killed types
  removed at prune N, matched by table) and **#objects ≥ LB**; `check_d2` (unpruned record point).
  Closure entries of the requested stages at the identity coset: Gram G[b][a] = h_b · g_a and
  pairing P[b][a] = H^y_b · d_id · G^x_a, in `closure_entries` order (weights 2w = ins + outs,
  resp. q_s − q_t + ins_s + outs_t).
- `box_word(lib, word, BoxOpts)`: the one-word box loop (≈ 100 lines): record at a generic point,
  replay (`Mode::Replay`, frozen structure, D14) at every needed box point of the unreduced grid in
  parallel (`par_map`), `same_structure` and identical closure skeletons asserted per point,
  `BoxAcc` per (stage, mode) — the reduced grid's points are the μ^0 = ∅ subset of the same scans;
  `finish` (down-set and residual checks), `BoxClosureData::from_box`, `closed_d2_freivalds`
  (D2Cert, k = 3; also on unpruned data), `close_box` at every N' ≤ N.
- `scaling_check`: scans at z_μ and λ·z_μ (scaled node sequence) at 2 box points, `check_scaling_all`
  on every entry of every stage (D14).
- `box_x`: Arena coordinates of z_μ = the roots a_{μ_j + λ_i + 1 − j} block by block.
- `NoProvider`: a provider refusing every call (library from layer A only).

**src/bin/sbox.rs**: `sbox --lam --word --n --ns --field 31|32b --seed --stages --prune --cache
--costs --replay --modes u,r --scaling --json OUT [--record-only] [--lazy-markowitz]`; prints and
stores the closures, χ = `schur::decat_chi` for every closure (reduced: [N choose λ_0]·χ),
per-point costs (per stage), record statistics, LB, final types.

**tools/g2b_compare.py** (closures vs m2s JSON and Wang blocks), **tools/stretch_compare.py**,
**examples/e2ref.rs** (engine2 uncolored references), run scripts runs_g2b/*.sh, runs_stretch/run333.sh.

**Restoring the (3,3,3) library from data-gl6/m5b** (sha256 of both chunks = MANIFEST): the cost
table is host-stamped (M5b ran on a 2.30 GHz Xeon) and is rejected on this machine, which changes
the provider stamp and hence every layer-A key (libbuild's own replay fails the same way). With
`KRHOM_COSTS_ANY_HOST=1` the table is admitted, the stamp matches, and every map loads from layer A
(re-certified on load; the provider refuses every call). Recorded as an open problem (§7).

## 2. Tests (tests/swire.rs; default 5 tests, 4.2 s release)

| test | what | result |
|---|---|---|
| `sengine_against_xengine_22_32` | (2,2) σ^3, σσ^{-1}σσ, σ^3 prune 3, σ^{-2} (F32b); (3,2) σ^2, σσ^{-1}, σ^4 prune 4: independent records at 2 points = identical structure; replay of point 0 at point 1 identical; final object multisets = xengine; K_0 and LB every crossing; d∘d = 0 (unpruned) | pass |
| `sengine_against_xengine_33` (ignored) | (3,3) σ^2, σ^3, σ^3 prune 4, σ^{-2} prune 5 | pass (195 s) |
| `objects_111_against_engine2` | (1,1,1) 4_1, σ_1^3σ_2, σ_1σ_2σ_1^{-1}σ_2, (σ_1σ_2)^3: final object counts = engine2 (12, 8, 10, 15) | pass |
| `box_loop_against_exact_22_32` | box loop (record + replays at the box points of N = 6) vs xengine `close_stage_marked`: (2,2) σ^3 stages 2, 3 (F31 and F32b), (3,2) σ^2; N' = 3..6, unreduced and reduced; scaling check at 2 points | pass (every closure equal) |
| `type_counts_41_333` | §5 | prints |

sengine's own suite still passes after the engine changes (11 + 1 ignored), including
`markowitz_heap_equals_rescan_order` (the exact heap keeps the rescan order).

## 3. Gate G2B

Runs: `runs_g2b/runw.sh M FIELD SEED` = `sbox --lam 3,3 --word σ^M --n 6 --ns 3,4,5,6 --modes u,r
--scaling` (unreduced extraction at all 399-400 points of L(18), reduced from the μ^0 = ∅ subset);
comparison `tools/g2b_compare.py docs/colored_design/gl6/m2s/sigma9_33_F31.json runs_g2b/*.json`.
Node sequences: seeds 101/103 (F31), 102/104 (F32b).

| stage | link | unreduced N'=3/4/5/6 | reduced N'=3/4/5/6 | reference | result |
|---|---|---|---|---|---|
| 1 | σ (unknot) | 1/4/10/20 | 1/1/1/1 | [N choose 3]_q (χ) | PASS |
| 2 | σ^2 Hopf | 1/16/100/**400** | 1/4/10/20 | m2s; Hopf formula | PASS |
| 3 | σ^3 T(2,3) | 1/10/46/**136** | 1/3/9/27 | m2s; Wang 3 N 1 | PASS |
| 4 | σ^4 T(2,4) link | 1/22/184/**864** | 1/6/26/88 | m2s | PASS |
| 5 | σ^5 T(2,5) | 1/16/106/**396** | 1/5/25/125 | m2s; Wang 3 N 2 | PASS |

Every run (F31 seeds 101, 103; F32b seeds 102, 104) gives all 32 stage-2..5 closures equal to the
m2s JSON of its prime as **full polynomials** (`g2b_compare.py`: 0 differences), the Wang blocks
3 N 1 and 3 N 2 at N = 3..6 as full polynomials, χ = `decat_chi` on all 40 closures, and the
scaling check (3 724 entries checked at z and λz at 2 box points, all stages; 3 808 zero at both).
Freivalds D2Cert (k = 3) recorded for every closure. Scan wall 940-1 174 s per run (399 points,
4 threads), record point 21-22 s, closures ≈ 1 s, VmHWM 0.9-1.1 GB.

Not run: σ^6..σ^9 unreduced (T(2,6) 1568, T(2,7) 848, T(2,8) 2560, T(2,9) 1540): the σ^9 attempt
(F31 seed 101) was refused by the per-point memory budget (2 GB) at crossing 8 of the record point
after ≈ 25 min (runs_g2b/s9_f31_seed101.err) — see §4.

## 4. Per-point costs

(3,3) σ^M at N = 6, F31, replay at box points with 4 threads (cumulative engine seconds per point
after each crossing, mean over 399 points; engine peak bytes = the per-point budget counter).

| M | s per point (replay) | engine peak per point | record point s (exact Markowitz, K_0 + d∘d) | inserted → minimal | pivots |
|---|---|---|---|---|---|
| 1 | 0.0004 | 18 KB | 0.0 | 4 → 4 | 0 |
| 2 | 0.08 | 12 MB | 0.09 | 189 → 15 | 87 |
| 3 | 0.45 | 19 MB | 0.43 | 928 → 40 | 444 |
| 4 | 2.45 | 48 MB | 3.4 | 2 599 → 85 | 1 257 |
| 5 | 9.75 | 130 MB | 17.5 | 5 580 → 156 | 2 712 |
| 6 | __S6__ | 364 MB (record) | 48-66 | 10 249 → 259 | 4 995 |
| 7 | __S7C__ | | | | |
| 8 | — | > 2 GB at the record point (refused by the budget) | > 20 min | | |

σ^5 box run: 977 s scan wall for 399 points (4 threads), 21 s record, 1 s for all 40 closures
(including Freivalds), VmHWM 0.9 GB; total 18 min per prime and node sequence.

**Where the time goes** (profiling counters `KRHOM_SPROF`, gdb stack samples; σ^5/σ^6 record):
- the arithmetic is negligible (crossing 5: 1.2 M pivot updates, 2.6e7 multiply-adds, 1.8 s);
- the exact Markowitz heap's bookkeeping dominated the record point (crossing 5: 2.9 M re-pushes in
  one round). A per-object index of equal-(type, t) neighbours (the only possible candidates) keeps
  the order identical and helps little, because the candidate graphs of (3,3) are dense (4 types).
  `PivotOrder::Lazy` (lazy-update Markowitz) removes the overhead but raises fill-in 6×
  (4.4 M → 28 M updates at crossing 6), so it is off by default (`--lazy-markowitz`);
- replay points skip the heap; their cost is the construction of C ⊗ W (π K ι per summand pair and
  coset, now sparse in f and over the nonzero rows of K ι) and per-update allocation in the
  elimination (one LocMap and BTreeMap entry per fill-in edge);
- growth is ≈ 4× per crossing in time and ≈ 2.5-3× in memory, like xengine's (its crossing 9 took
  2 460 s symbolically). Projected σ^9 replay: of the order of 1 h and several GB per point.

## 5. 4_1 type counts (C_feas m9)

From the decategorified scan of 4_1 = σ_1σ_2^{-1}σ_1σ_2^{-1} on (3,3,3) (`type_counts_41_333`):

| prune | LB objects | types in the final LB complex | types over all stages |
|---|---|---|---|
| none | 328 | 40 | 54 |
| 6 | **276** | **34** | 40 |
| 5 | 77 | 16 | 20 |
| 4 | 12 | 5 | 6 |

Closable types used by the actual scan (prune 4, both primes, both node sequences): **5**, with 12
final objects = LB. At N = 6, R3's "34 of the 40 types" is right (34 types in the final LB complex
at prune 6; 40 = all types of the unpruned final complex, or all types touched at prune 6 over the
stages). The "38" of SPEC/R4 §9 matches none of these counts and should be replaced by 34 (T0 must
cover the 34 final types; 40 if every stage is closed).

## 6. Stretch: first 3-strand 3-colored runs through the localized path (N = 4, prune 4)

`runs_stretch/run333.sh`: λ = (3,3,3), prune 4, library restored from data-gl6/m5b (F31 and F32b;
27 layer-A loads, 0 provider calls, 6 types, all closable), box N = 4 (64 points unreduced, 16 reduced),
closures N' = 3, 4, both modes, every stage; seeds 201, 202 per prime.

| word | N = 4 unreduced | N = 4 reduced | checks |
|---|---|---|---|
| 4_1 = [1,−2,1,−2] | rank 16: t^-2q^7 + t^-2q^9 + t^-2q^11 + t^-1q^-1 + t^-1q^1 + t^-1q^3 + t^0q^-3 + t^0q^-1 + t^0q^1 + t^0q^3 + t^1q^-3 + t^1q^-1 + t^1q^1 + t^2q^-11 + t^2q^-9 + t^2q^-7 | rank 5: t^-2q^8 + t^-1q^2 + t^0q^0 + t^1q^-2 + t^2q^-8 | = engine2 uncolored gl_4 of the reversed braid [−2,1,−2,1] (both modes): **via duality (expected)**; χ = decat_chi; N = 3 rank 1 |
| [1,2] | 4 | 1 | = unknot [4 choose 3] |
| [1,1,2] | 16 | 4 | = 2-strand σ^2 (Hopf) of m2s at N = 3, 4, both modes |
| [1,1,1,2] | 10 | 3 | = 2-strand σ^3 (T(2,3), = Wang 3 4 1) |
| [1,1,1,1,2] | 22 | 6 | = 2-strand σ^4 (T(2,4) link) |
| [1,1,1,1,1,2] | 16 | 5 | = 2-strand σ^5 (T(2,5), = Wang 3 4 2) |

The Markov comparisons hold as full polynomials with zero shift (the global q-shift absorbs the
stabilization). All 24 runs (6 words × 2 primes × 2 node sequences) are identical per word; χ =
`schur::decat_chi` on every closure (192); scaling check on 346-3 146 entries per run; Freivalds
D2Cert recorded before every pruned closure (`close_box` refuses without it). Per point 0.04 s,
2.7 MB. The intermediate stages are closable too (labels (3,3,3) throughout): σ_1 alone gives the
2-component split link (16 = 4·4), σ_1σ_2^{-1} the unknot (4).

## 7. Open problems

1. **Engine speed and memory at depth** (blocks G2B σ^6..σ^9 unreduced and every larger 3-strand
   run): per-update allocation (LocMap + BTreeMap per fill-in edge), the Markowitz bookkeeping on
   dense candidate graphs, and per-point construction of C ⊗ W. Options: an arena/slab store for
   components; batched Gaussian elimination of each term's isomorphism block (the d_C ⊗ id part of
   a term is a direct sum of copies of C's differential: eliminate it once per type pattern);
   reusing point 0's construction structure (which π K ι blocks are nonzero) at replay points.
   σ^8 needs > 2 GB per point at the record point.
2. The cost table's host stamp enters the layer-A key, so a library built on one machine is not
   loadable on another without `KRHOM_COSTS_ANY_HOST=1` (M5b open problem 4: stamp the provider
   route, not the table).
3. `check_d2` and the K_0/LB asserts run at the record point only (they are point independent given
   the frozen structure); `--d2-all` turns d∘d on at every box point.
4. The checkpoint/shard machinery of `BoxAcc` is not used by the one-word loop (all points in one
   process); M6a-2a's `box_run` / `box_finish` add it.
5. The "38 types" figure in SPEC §2 / R4 §9 should be corrected to 34 (§5).
