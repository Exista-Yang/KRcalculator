# M6a-1c report: a faster per-point scan, G2B at depth, and the 3-strand check at N = 3, 4

Date 2026-11-12. Branch colored-m6a (merged colored-m5 f2fcb1d first). Machine: cloud container, 4 cores
(Xeon @ 2.10 GHz), 15 GB. Results over F_p, p = 2^31−1 and 4294967291.

## 0. Verdict

| item | result |
|---|---|
| profile of the per-point scan (σ^4..σ^7 on (3,3)) | done: docs/colored_design/gl6/M6a1c_profile.md |
| result-preserving speedups | replay point **3.6-3.9×**, record point **5-8×** (σ^4..σ^7); identical closure values (value digests), identical pivot traces, tests pass |
| G2B σ^6 at N = 6 by box, unreduced + reduced | **PASS** both primes (T(2,6) link unreduced 1568, reduced 252 at N = 6; all 48 closures of stages 2..6 = m2s, Wang, χ, scaling) |
| G2B σ^8 reduced | **PASS** both primes (T(2,8) reduced 1/10/82/**560** at N' = 3..6; stages 2..8 = m2s) |
| G2B σ^7 unreduced + reduced | **PASS** (T(2,7) unreduced 1/22/190/**848**, reduced 343 at N = 6; stages 2..7 = m2s, Wang 3 N 3) — both primes |
| G2B σ^9 reduced | **PASS** both primes (T(2,9) reduced 1/9/81/**729**; all reduced stages 2..9 = m2s) |
| G2B σ^8, σ^9 unreduced | **not run** (projection §4) |
| (3,3,3) 185 knots at N = 3, 4 (prune 4), unreduced + reduced, both primes | **185/185 PASS**: N = 4 = uncolored gl_4 of the reversed braid (engine2) **exactly** (no shift), "via duality (expected)"; N = 3 = 1 (t^0 q^0) for every knot, both modes; χ = `decat_chi`; primes identical |
| pruned vs unpruned (3,3,3) | **not feasible**: prune 5 / 6 types are not in layer A and the provider refuses them (§5) |

## 1. Profile (summary; details in M6a1c_profile.md)

Before (c9d97b9, one thread, F31): replay σ^4/5/6 = 2.6 / 10.5 / 33.6 s; 72-78 % in **d_C ⊗ id**
(π K ι per summand pair: ≈ 400 pairs per (component, coset), one hash lookup + one allocation per pair,
1.6e9 multiply-adds at σ^6), 22-27 % in **elimination** (0.7-1.5 µs per update for ≈ 21 multiply-adds:
overhead-bound), < 3 % everything else (dec_word, Rickard maps, closure extraction, pattern hashes, K_0 —
library evaluation is memoized per point by `Eval`, 11.8 MB). Record point: 56-65 % in the **Markowitz
heap** (σ^5: 5.6 M pushes and pops, heap of 5.4 M stale entries = 130 MB uncounted). Memory: at the peak
(σ^5) 4 332 objects and 404 370 components during the insertion of the last step; 292 bytes per
component (≈ 164 data, ≈ 128 fixed overhead).

## 2. Optimizations (all in src/sengine.rs, src/swire.rs; no frozen file touched)

1. **d_C ⊗ id** (`swire::dcid_step`, `Asm`, `assemble`): ι_all = [ι_1 | ι_2 | …] and π_all = [π_1; π_2; …]
   (CSC by columns of X) assembled once per (type, coset) and step at the point; per (component, coset)
   one K ι_all (sparse in f) and one π_all · (K ι_all) over the nonzero rows, then cut into the summand
   pairs straight into the LocMaps. Exact arithmetic, same blocks, same insertion order. The old per-pair
   path is kept behind `KRHOM_DCID_REF=1` (same digests).
2. **Elimination** (`CompPlan`, `apply_pivot_inner`): block plans per (x type, pivot type, y type)
   (offsets and m, k, n; no binary searches), per-pivot local cache of (layout, plan), one BTreeMap
   `entry` per update, touched set as a sorted Vec, a bounded pool (16 MB, counted in the budget) of freed
   component buffers reused for fill-in edges.
3. **Record point** (`markowitz_heap_fast`, now `PivotOrder::Heap`; the M6a-1b version is
   `PivotOrder::HeapV1`): the same order (minimum of (cost, a, b) = the rescan order;
   `markowitz_heap_equals_rescan_order` passes on 400 scans) with dense per-object arrays (degrees,
   versions, keys, eligibility fixed within a round, candidate lists) and **compaction** of the heap
   (drop stale entries and duplicates when it exceeds 2× its last valid size; pop order of valid
   entries unchanged).
4. Memory accounting: component count kept incrementally (`SCx::components` O(1)); the peak's
   objects/components/previous-complex bytes recorded (`MemStats`); `EngOpts::peak_detail`.
5. Tools: `src/bin/sprof.rs` (per-part profile, value and trace digests), `--stamp-route` for sbox/sprof
   (load the mixed-stamp m5c layer A), `tools/k3_compare.py`.

**Verification.** Value digests of all closure entries at two replay points and the record trace digest
are identical before and after at σ^4, σ^5, σ^6 (runs_m6a1c/prof_s*.base.txt vs *.after.txt; σ^7 digests
unchanged between the intermediate and final code); G2B σ^5 box with the new code = the M6a-1b run term by
term (all 40 closures) and = m2s; sengine (11) and swire (4) suites; full suite §6.

### Speedups (one thread, F31, mean of 2 replay points; after = final code under load of 3 other jobs)

| M | replay before | replay after | × | d_C⊗id before → after | pivots before → after | record before → after | × |
|---|---|---|---|---|---|---|---|
| 4 | 2.63 s | 0.71 s | 3.7 | 2.06 → 0.44 | 0.47 → 0.18 | 4.45 → 0.89 s | 5.0 |
| 5 | 10.5 s | 3.06 s (2.4 unloaded) | 3.4-4.4 | 7.9 → 1.6 | 2.4 → 1.1 | 24.4 → 3.4 s | 7.2 |
| 6 | 33.6 s | 9.1 s | 3.7 | 24.3 → 4.5 | 8.8 → 4.2 | 104.6 → 12.6 s | 8.3 |
| 7 | 87.7 s (M6a-1b, 3 thr) | 24.8 s | 3.5 | — → 10.1 | — → 13.8 | 231.9 → 36.3 s | 6.4 |

Box runs (3 threads + one 1-thread job alongside): σ^5 399 points 977 s → 383 s wall; σ^6 400 points
1 190-1 613 s wall (M6a-1b: not run); σ^8 record point 108 s (M6a-1b: refused at > 2 GB after 25 min, budget
now 3.5 GB). Memory is unchanged by these changes (the per-point peak is the dense unsimplified step):
48 / 131 / 323 / 758 / 1 682 MB at σ^4..σ^8.

## 3. G2B at N = 6 (runs_m6a1c/g2b, `tools/g2b_compare.py` against docs/colored_design/gl6/m2s/sigma9_33_*.json)

| run | points | per point | record | result |
|---|---|---|---|---|
| σ^5 u+r F31 seed 101 | 399 | 2.9 s | 3.5 s | 40 closures = m2s, Wang; = M6a-1b run |
| σ^6 u+r F31 seed 101 | 400 | 12.1 s (contended) | 11.3 s | 48 closures = m2s (stage 6 unreduced 1/28/292/**1568**, reduced 1/8/50/252), Wang 3 N 1, 3 N 2; χ 48/48; scaling 8 746 entries |
| σ^6 u+r F32b seed 102 | 400 | 8.9 s | 14.9 s | same, 0 differences |
| σ^8 r F31 seed 101 | 20 | 72 s | 108 s | 32 reduced closures of stages 2..8 = m2s (σ^7: 1/7/49/343; σ^8: 1/10/82/**560**); χ 32/32; scaling 43 822 entries |
| σ^8 r F32b seed 102 | 20 | 70 s | 106 s | same, 0 differences |

| σ^9 r F31 seed 101 | 20 | 204 s | 286 s | 36 reduced closures of stages 2..9 = m2s (σ^9: 1/9/81/**729**); χ 36/36; scaling 85 856 entries; engine peak 3.28 GB per point, VmHWM 9.0 GB (2 threads) |
| σ^9 r F32b seed 102 | 20 | 189 s | 268 s | same, 0 differences |
| σ^7 u+r F31 seed 101 | 400 | 26.7 s (2 threads, contended) | 36 s | 56 closures = m2s (stage 7 unreduced 1/22/190/**848**, = Wang 3 N 3); χ 56/56; scaling 20 726 entries |
| σ^7 u+r F32b seed 102 | 400 | 28.5 s (3 threads, contended) | 40 s | same, 0 differences (848 = Wang 3 N 3) |

One node sequence per prime (F31 seed 101, F32b seed 102); the second node sequences of M6a-1b
(seeds 103/104) were not rerun at depth (σ^5..σ^7 reduced already passed them in M6a-1b). Not run:
σ^8 and σ^9 unreduced (§4). Every closure of every run equals the m2s JSON of its prime as a full
polynomial (`g2b_compare.py`: 0 differences).

## 4. Projection for the rest of G2B

Replay cost now grows ≈ 2.7-3× per crossing (σ^5 → σ^8: 2.4, 9.1, 24.8, ≈ 70 s), memory ≈ 2.2-2.5×
(131, 323, 758, 1 682 MB). Projected σ^9: ≈ 200-220 s and ≈ 4 GB per point. σ^8 unreduced: 400 × 70 s /
3 threads ≈ 2.6 h per prime and node sequence; σ^9 unreduced: 400 × 210 s / 2-3 threads (memory) ≈ 8-12 h
per prime — beyond one session; reduced σ^9 (20 points) ≈ 1 h. At σ^7, elimination is 55 % (23.6 M
updates, memory-latency bound: ≈ 0.6 µs each on a 2.5 M-component complex) and d_C ⊗ id 40 %
(5e9 multiply-adds class, mod-p bound).

## 5. (3,3,3) at N = 3, 4 (runs_m6a1c/k3, docs/colored_design/gl6/gl4_333.json, runs_m6a1c/k3_summary.md)

Library: data-gl6/m5c layer A (sha256 of both chunks = MANIFEST), `--replay --stamp-route`,
`KRHOM_COSTS_ANY_HOST=1`, costs_333_m5c.txt; 0 provider calls. Prune 4, box N = 4 (64 points unreduced,
16 reduced), closures N' = 3, 4, both modes, the final stage, scaling check at 2 points; F31 and F32b,
node seed 201.

- **185 of 185 knots** of tools/dev/ntrunc_t1/data/knots3_all.txt (braid lengths 4 … 12): every closure
  computed for both primes, the primes identical, χ = `schur::decat_chi` on every closure (740).
- **N = 4: equal, as full polynomials and with zero shift, to the uncolored gl_4 homology of the
  reversed braid (engine2, examples/e2ref) — unreduced and reduced, all 185 knots** ("via duality
  (expected)": ∧^3 C^4 ≅ (C^4)^* ⊗ det). Unreduced ranks 16 … 796, reduced 5 … 265; unreduced = 3·reduced + 1 for 163 of the 185
  knots (the thin ones).
- **N = 3: t^0 q^0 (rank 1) for every knot, unreduced and reduced** (∧^3 C^3 is the trivial rep, so the
  colored invariant of a knot is that of the unknot = 1 at this normalization).
- Second node sequence (seed 202) on the 58 knots of ≤ 10 letters: all 116 runs (58 knots × 2 primes) identical to
  seed 201, closure by closure.
- Cost: per point 0.04 s (4 letters) … 0.28 s mean at 12 letters; per knot ≤ 53 s; VmHWM ≤ 302 MB.
- **Pruned vs unpruned: not feasible.** 4_1 at prune 5 and 6 needs dec_word(1,2,0/2,0,1/0,1,2 ⊗ E1.1,F1.1)
  (resp. 2,0,1/1,2,0/0,1,2 ⊗ E0.2,F0.2) which is not in m5c layer A; with the provider allowed it refuses
  (Exact not admitted: labels (3,3,3) caps, fibre product 432 000, no cost-table entry; Search tiers refused
  by admission) (runs_m6a1c/k3u). The engine2 comparison is an independent unpruned reference at N = 4.

## 6. Tests

`KRHOM_THREADS=4 cargo test --release` (final code, while production jobs ran): **240 passed, 0 failed,
28 ignored** (runs_m6a1c/full_suite.log). sengine 11 passed (incl. `markowitz_heap_equals_rescan_order`,
`streaming_equals_full_and_known_model`, `memory_accounting_and_refusal`), swire 4 passed
(`box_loop_against_exact_22_32`, `sengine_against_xengine_22_32`).

## 7. Open problems

1. **Elimination at depth** is now the main cost (55 % at σ^7): ≈ 0.6 µs per update on a 2.5 M-component
   complex, memory-latency bound (BTreeMap out/inn, one heap allocation per component). Options: u32
   component storage (both primes < 2^32; −28 % bytes per component), a slab/arena of components with
   compact per-object adjacency arrays, and the hash of the pattern computed incrementally.
2. **The dense unsimplified step** (≈ 93 components per object at the peak; σ^8: 1.7 GB per point)
   bounds the threads per machine at σ^9 (≈ 4 GB per point). Delayed elimination (apply pivots inside a
   step and propagate later arrows through the stored Schur complement) would cut the peak, at the cost of
   a more complex replay.
3. d_C ⊗ id arithmetic: delayed modular reduction (u128 accumulators) in π_all · (K ι_all); the assembled
   ι/π could be shared across crossings with identical (type, term) (σ^M on (3,3): every crossing).
4. G2B σ^8 / σ^9 unreduced (≈ 2.6 h and ≈ 8-12 h per prime at 400 points): needs (1)-(2) or the sharded
   box driver of M6a-2a (`box_run` / `box_finish`) across sessions.
5. Unpruned (3,3,3) types at levels 5/6 for the pruned-vs-unpruned check (library, M5c open items).
