# R4: M6a retargeted to 3-colored 3-strand knots at gl_6 (Foam_N pruning)

Design note, 2026-11-10, revision 2. It replaces the first draft at this path; that draft is kept as `gl6/.R4_m6a_draft1.md`.

Inputs:
- SPEC_M4_M6a.md rev 2 (§6, Appendices M and R), SYNTHESIS.md, NTRUNC_STUDY.md, HANDOFF_M4b.md;
- M4a/M4b code (worktree m4b, 59eaea2: smap.rs, homprov.rs) and, on dev, xengine.rs (prune_n) and cclose.rs;
- the sibling notes gl6/R1_2strand.md, gl6/R2_m5prune.md and gl6/R3_closure.md;
- the truncated scans ntrunc/t1/out/scan_knots3*.txt and scan_t2n.txt.

Tags: [M] measured earlier; [C] computed for this note (gl6/py/*.py, < 1 s each); [E] estimate. Every cost row is [E] until the calibration run S3 (§7).

What changed from the first draft:
1. **Markov words are much cheaper than the draft said** (§5.2): their peak P is 20-500x smaller. The draft multiplied the 2-strand rank by C(9,3), which is wrong. As a result, unreduced N = 6 checks against Wang become affordable.
2. **Wang's references are unreduced.** Reduced N = 6 results are therefore checked against R1's exact 2-strand closure, not against Wang.
3. The library cost now agrees with R2.
4. u32 storage is added.
5. The median-LB knot is 12a_859, named explicitly.
6. All of M6a uses one closure path (box points).
7. The memory taken by the evaluated library is now counted per thread.
8. Sections are renumbered. R3 cites the draft's §1.2 (box points), which is §2 here, and the draft's §2-3 (cost model), which are §4-5 here.

## 0. Recommendation

1. **Closure data at fixed N, from Grassmannian box points** (SYNTHESIS CL2 restricted to N ≤ 6; R3 owns the details).
   - SPEC §6.5 interpolates in e-monomials of R^λ. At the 3-strand Gram weight 22 that is 88 233 unknowns per entry (3 866 reduced), so §6.5 refuses even 4_1.
   - Box extraction needs at most 48 unknowns per entry reduced (777 unreduced). The scan count per prime is 400 for reduced N = 6 and 8 000 for unreduced.
   - The N = 6 point set contains the N = 4 and N = 5 point sets.
2. **Realistic M6a targets on the cloud box** (4 cores, 15 GB, sessions of 5 h; two primes; central model):
   - **Markov words** [1,1,2], [1,1,1,2], [1,1,1,1,2], [1,1,1,1,1,2] at N = 4, 5, 6. P ≤ 6e4; the unreduced N = 6 run takes about 1-2 h wall per word (4-7 h at the high end of τ).
     - Unreduced, they must equal Wang: ranks 136 and 396 at N = 6.
     - Reduced, they must equal R1's exact results.
     - Only a small corner of the (3,3,3) library is needed, so these runs can start before G1(c) finishes.
   - **Class A, reduced N = 6:** 4_1 2.4 h wall, 5_2 7.9 h, 6_2 9.7 h.
   - **Class B, reduced N = 6** (P 1e7-2e7: 7_3, 8_2, 8_19, 9_3, 8_20, 10_2, 12a_722): 17-29 h wall each with u32 storage, i.e. 4-6 sessions each. A stretch goal, attempted only after S3.
   - **N = 3 and 4, all 185 knots:** about 5 h CPU (one session), checked against uncolored gl_4 by duality.
   - **N = 5:** reduced for the 86 knots with P ≤ 1.3e7. Unreduced for 4_1 through 8_19, at ≤ 3 h wall each.
3. **Unreduced N = 6 for genuine 3-strand knots needs M6c.**
   - Box extraction costs 20x the reduced run: 4_1 would take 47 h wall.
   - M6c needs emono(w_c) + 8 scans per prime, valid at every N: 1 469 scans if the coordinate weight w_c is ≤ 10, i.e. 4_1 in about 9 h wall. w_c is unmeasured for (3,3,3).
4. **Beyond class B, both memory and time per point are too large.**
   - The median-LB knot 12a_859 has P = 6.8e8: 110-210 GB and 3.8 h per point. 12a_1248 needs 2.6-4.9 TB and 91 h per point.
   - 175 of the 185 knots have P > 2e7 at N = 6.
   - These need M7a (R-coord) or a new idea. Time × 800 points stays the binding limit even if per-source streaming solves the memory.
5. **Minimal M6a scope** (§7):
   - the sengine core with prune_n, term streaming, K_0/LB checks and the determinism trace;
   - the box closure layer;
   - a resumable multi-point driver that runs a batch of words at each point;
   - a staged exact closure;
   - u32 storage, only if S3 shows class B is memory-bound.

   §8 lists what is dropped or deferred.
6. **The dominant unknown is the one-time (3,3,3) library** (M5 plus G1(c)).
   - Size: 42 / 158 / 338 pairs with a kept source at N = 4 / 5 / 6 (R2 §1).
   - [E] 2-10 h CPU per prime at R2's 20-100 s per pair, rising to 20-200 h if the cost per target grows with the X flows (up to 31 500).
   - G1(c) level 4 measures it first.
7. **G1(a) does not block this target.** Its failures are all (2,2,2). M6a's precondition becomes G1(a)-(3,3) (passed, 80/80) plus G1(c) plus T0 on the kept (3,3,3) types.
8. **First certified 3-colored gl_6 3-strand results**, in order:
   - (a) the Markov words at N = 6, unreduced (= Wang) and reduced (= R1). They run on the 3-strand engine but their closures are 2-strand knots.
   - (b) 4_1 reduced at N = 6, the first new value, together with 4_1 at N = 4 (duality) and at N = 5. Two primes on two node sequences.

## 1. (3,3,3) object sizes [C]

SPEC Appendix R holds only 2-colored reference values (χ, duality at N = 3, Wang, Hopf). The (3,3,3) numbers below come from `scripts/cb.pkl` (`gl6/py/ranks333.py`, `cdiag333.py`).

Definitions. Localization runs over the 1 680 cosets of S_9/S_3^3. For a type p:
- rank (number of flows): r(p) = Σ_{q≤p} N(q) h_qp(1);
- End block size: S(p) = Σ_q N(q) h_qp(1)²;
- N(q) = Π_i 3!/Π_j q_ij! is the number of cosets in stratum q.

| prune N | types kept | Σ r | max r | max S | max h_qp(1) | Σ c_diag(1) = closure gens | Gram entries Σ c_diag² | max Gram weight (scans) | box top degree, unreduced / reduced |
|---|---|---|---|---|---|---|---|---|---|
| 4 | 6 | 409 | 144 | 324 | 9 | 31 | 217 | 10 | 9 / 6 |
| 5 | 20 | 9 649 | 4 320 | 52 920 | 54 | 172 | 3 808 | 18 | 18 / 12 |
| 6 | 40 | 40 489 | 7 200 | 68 688 | 54 | 344 | 6 800 | 22 | 27 / 18 |
| none | 55 | 76 399 | 8 820 | 68 688 | 54 | 420 | 7 422 | 27 | - |

- **Largest kept types at N = 6.**
  - 0,1,2/2,0,1/1,2,0 and its transpose: r = 7 200, S = 68 688 (0.55 MB per End block at 8 B), 27 closure generators in degrees 10-22.
  - (1,1,1)^3: r = 4 320, h = 54 at the identity stratum, 54 generators in degrees 4-18; kept from N = 5 on (μ_1 = 5).
- **Ranks multiply under tensor product:** rank(L_p·W) = r(p)·r(W).
  - The four terms of a (3,3) crossing have r = 1, 12, 30, 20 (Σ 63, R1 §1). So the peak of one crossing step is at most 63 times the LB rank before it.
  - A (3,3) type p ⊗ id_3 keeps its 2-strand rank, r = r_2(p): the third block of its table is diagonal, so N(q) gains no factor.
- **Idempotents.** The only decomposable (λ,λ) ladder, 0,1,2/1,1,1/2,1,0 (r = 8 820, μ_1 = 7), is killed at N ≤ 6, so e_p = id for every kept (λ,λ) type. Two intermediate-block cases are still open (R2 §2.4).
- **Ranks, not object counts, set the sizes.** At N = 6, over all 185 knots, the localized rank per object is 310-1 020 at the peak and 960-2 370 on the LB complex. ΣS/Σr = κ_S is 3.1 at the untruncated 4_1 peak [M, C2 H2] and 6.7 over all kept types.
- **Gram weight against box top degree.** Already at 3 strands the maximal Gram weight exceeds the box top degree: 10 > 9 at N = 4 (unreduced), and 22 > 18 at N = 6 (reduced). Such entries are 0 modulo I_N.
  - This is the same phenomenon as the (3,3,3,3) question (40 > 36), which R3 §1 resolves: it is neither a bug nor a convention mismatch.
  - T0 at N = 4 on the 6 kept types exercises it cheaply.

## 2. Closure data at a fixed N (box points; only what M6a needs)

**Construction.**
- Node sequences: a^(i) for i = 0, 1, 2, six random values per block, all 18 distinct.
- Points: x_μ = (a^(0)_{μ0+δ}, a^(1)_{μ1+δ}, a^(2)_{μ2+δ}), with μ_i ⊆ 3 × (N−3).
- Vanishing of factorial Schur functions makes the solve triangular on the lower set.
- For an entry homogeneous of weight w, the top coefficients are the Schur coefficients, so f mod I_N = Σ_{ρ⊆box, |ρ|=w} c_ρ s_ρ.
- Converted to GrAlg's e-monomial box basis, this gives a `ClosureDataC` valid at that N only (new field `valid_n`). `close_exact_colored` accepts it unchanged.

**Variants.**
- Reduced theory (block 0 marked, augmentation, `drop_block`): only μ0 = ∅ is needed.
- Nesting: with the same nodes, the N' < N boxes are lower subsets of the N box.

**Scans per prime** (lower set |μ| ≤ w) [C, `counts.py`]: 64 / 16 points at N = 4, 1 000 / 100 at N = 5, and 8 000 / 400 at N = 6 (unreduced / reduced). At N = 6, w ≤ 22 needs 7 920 / 400 points and w ≤ 10 needs 1 834 / 269.

**How each part is obtained.**
- **Gram:** generator-only evaluation, no scan.
- **Differential data:** one scan per point. R3 §0 refines this:
  - ship the closed pairing P = H_id·d_id·G_id instead of F, so no solve with G'_id is needed at structured points;
  - accumulate per point in a streaming way, with output in Schur coordinates;
  - certify the pruned d∘d by Freivalds on the closed complex, with `closure_kills_map` kept as the oracle.

**Checks.**
- Residuals exist only where points with |μ| > w exist; there are none for w ≥ 18 (reduced) or w ≥ 27 (unreduced). There, two primes on two node sequences take their place.
- The closure-generator determinant is asserted at every point; this guards against degenerate box points.

**Reduced or unreduced.** Wang's blocks in `tests/data/wang_ref.txt` (and Appendix R.3) are **unreduced**: k = 2, N = 3, T(2,3) has rank 7, which equals R.2's unreduced 7 (reduced 3). Unreduced costs 20x reduced at N = 6, 10x at N = 5 and 4x at N = 4.

**2 strands.** R1's exact xengine route produces the gl_6 T(2,n) results. M6a's own 2-strand runs (box, 20 / 400 points at N = 6) are a regression against them.

## 3. The (3,3,3) library (M5 pruned, gate G1(c)), per prime [E]

- **Size.**
  - 42 / 158 / 338 pairs with a kept source at N = 4 / 5 / 6, of 494 (R2 §1). There are 3 151 unpruned targets [M, P11]; `pruned_plan` gives the kept-target counts.
  - Markov words need only the pairs with source p ⊗ id_3 under one crossing on strands (2,3): about 16 composite decompositions with X ≤ 900 flows.
- **What G1(a) measured (laptop, cold template cache)** [M, `g1a_run1_cold_partial.json`]:
  - (3,3), 12 pairs: mean 5.6 s, max 22.6 s.
  - (2,2,2), 46 pairs that returned (including the failed pairs 18, 25 and 32): median 0.3 s, mean 24 s, max 496 s (180 flows). Pair 20 hit the 600 s limit.
  - Wall time does not correlate with X flows: template construction dominates.
  - No (3,3,3) Hom has been certified yet, and (3,3,3) X reach 31 500 flows.
- **Estimate.**
  - N = 4: under 1 h.
  - N = 6: 2-10 h CPU at 20-100 s per pair (R2), up to 20-200 h at 50-500 s per target on about 1 700 targets. Multiply by 1.5 for the cloud box. One process per pair over 4 cores divides the wall time by 4.
  - The second prime replays the certified DAG: about 0.3-0.5x of the first.
- **Persistence.** Layer A goes on disk and is pushed as an archive (R2 §4.1), because sessions last 5 h.

## 4. Resource model per evaluation point (cloud Xeon, 4 cores, 15 GB)

**Input.** P is the decategorified peak localized dimension at prune N: Σ r over the kept summands of LB ⊗ next crossing, from the scan files. The pre-decomposition tensors exist only inside the sandwich π(f⊗id)ι, so they cost time, not memory.

**Memory per point:** M_pt = μ·P + M_ev.
- μ = c_e·κ_S·σ·b, where:
  - c_e = stored entries / ΣS: 5-30, unmeasured for colors;
  - κ_S: 3-7;
  - σ, the streaming factor: 0.5-1;
  - b = 8.5 B for a u64 `Fp` plus shared offsets, or 4.5 B with u32 entries (both primes are < 2^32).
- Central μ = 300 B (u64) or 160 B (u32); the band is 64-1 785 B.
- M_ev holds the library's ι/π blocks evaluated at the point. That is about 2 Σ_pairs S(X) ≈ 1e8 numbers at prune 6, i.e. 0.8 GB if cached for the whole scan (0.4 GB with u32). Evicting them per crossing trades this memory for time.

**Time per point:** t_pt = τ·P + t_lib.
- τ: central 2e-5 s, band 5e-6 to 1e-4. Upper bound: engine2's proxy times κ_S. Lower bound: about 100 modular multiply-adds per stored entry at 3e8 per second.
- t_lib: about 10 s at prune 6 (3-30), 3 s at prune 5, 0.1 s at prune 4. It is paid once per point for a whole batch of words.

**Threads:** min(4, ⌊12 GB / M_pt⌋). The other 3 GB hold 2 GB for the library DAG and template cache, plus 1 GB of headroom.

**Scans** (two primes on two node sequences): at N = 6, 800 reduced or 16 000 unreduced; at N = 5, 200 / 2 000; at N = 4, 32 / 128. Each point writes its samples to disk, so jobs resume across sessions.

## 5. Estimates per knot class (central; [C] `gl6/py/r4model.py` from `n4/n5/n6.txt`)

### 5.1 N = 6, prune 6

| knot | LB objects / ldim | peak objects / P | M_pt u64 / u32 | t_pt | threads u64 / u32 | reduced: CPU → wall (u64 / u32) | unreduced |
|---|---|---|---|---|---|---|---|
| 4_1 | 276 / 2.7e5 | 2 602 / 1.64e6 | 1.3 / 0.66 GB | 43 s | 4 / 4 | 9.5 h → **2.4 h** | 190 h CPU, 47 h wall (M6c) |
| 5_2 | 680 / 7.4e5 | 9 550 / 6.6e6 | 2.8 / 1.5 GB | 142 s | 4 / 4 | 32 h → **7.9 h** | 157 h wall (M6c) |
| 6_2 | 1 022 / 1.0e6 | 11 604 / 8.2e6 | 3.3 / 1.7 GB | 174 s | 3 / 4 | 39 h → 12.9 / **9.7 h** | M6c |
| 8_19 | 776 / 1.3e6 | 18 144 / 1.5e7 | 5.3 / 2.8 GB | 312 s | 2 / 4 | 69 h → 35 / **17 h** | M6c |
| 12a_722 (largest in B) | 2 662 / 3.3e6 | 60 932 / 1.9e7 | 6.5 / 3.4 GB | 391 s | 1 / 3 | 87 h → 87 / **29 h** | no |
| 12a_859 (median LB) | 22 452 / 2.7e7 | 1.2e6 / 6.8e8 | 210 / 110 GB | 3.8 h | 0 | 3 000 h CPU; does not fit | no |
| 12a_1248 (largest) | 7.4e5 / 1.0e9 | 2.3e7 / 1.6e10 | 4.9 / 2.6 TB | 91 h | 0 | no | no |

Bands with all constants at their low or high end:
- 4_1: M 0.1-3.7 GB, t 11-194 s.
- 5_2: M 0.4-13 GB, t 36-690 s.
- 8_19: M 1-28 GB, t 78-1 500 s.

So the class boundaries can move by a factor of 2-4 once S3 has measured the constants.

### 5.2 Markov words (3-strand words whose closures are 2-strand links)

[C] from `scan_t2n.txt` and rank multiplicativity (§1). This is an upper bound, since it ignores pruning.

| word | closure | LB ldim before σ_2 | P (σ_2 step) | peak of the σ_1 phase |
|---|---|---|---|---|
| [1,1,2] | Hopf | 265 | ≤ 1.7e4 | 4.0e3 |
| [1,1,1,2] | T(2,3) | 507 | ≤ 3.2e4 | 1.7e4 |
| [1,1,1,1,2] | T(2,4) | 709 | ≤ 4.5e4 | 3.2e4 |
| [1,1,1,1,1,2] | T(2,5) | 911 | ≤ 5.7e4 | 4.5e4 |

Cost:
- t_pt ≈ 0.4-1.2 s plus t_lib, and t_lib is small: only the p ⊗ id_3 corner of the library is used ([E] 0.1-1 s).
- Reduced N = 6: under 0.5 h CPU per word.
- Unreduced N = 6 (16 000 scans): 5-8 h CPU per word, i.e. **1.2-2 h wall** (4-7 h at τ = 1e-4). The same scans also give N = 4 and 5 by nesting.
- The unknot words [1,2], [1,-2] are smaller still.

### 5.3 N ≤ 5

- **N = 4, prune 4, all 185 knots.**
  - P ≤ 9.2e4 (median 2.8e4), M_pt ≤ 25 MB, t_pt ≤ 2 s.
  - 128 scans per knot; the reduced points are a subset of these.
  - In total about 5 h CPU, i.e. 1.3 h wall.
  - N = 3 is a single object (wedge^3 = det).
- **N = 5, prune 5.**

  | knot | P | reduced, wall | unreduced, wall |
  |---|---|---|---|
  | 4_1 | 1.4e5 | 0.1 h | 0.8 h |
  | 8_19 | 8.7e5 | 0.3 h | 2.8 h |
  | 12a_859 | 1.9e7 | 7 h (u32, 3 threads) | 70 h (M6c) |
  | 12a_1248 | 1.5e8 | 24 GB per point: does not fit | does not fit |

  86 knots have P(5) ≤ 1.3e7.

### 5.4 Closure step and one-time state spaces

- **Per-type closure state spaces at N = 6** (T0 data, independent of the knot, cached):
  - reduced: at most 54 × 48 rows per q-degree, seconds;
  - unreduced: up to 54 × 777 ≈ 42 000 rows, so it needs sparse or per-degree echelon.
- **`cclose::exact_mem_estimate`** multiplies dim A by all Gram terms of all types. For 3 colors it extrapolates to about 40 GB reduced and TB unreduced. It must stage per type and be recalibrated (risk R4).
- **Knot-level closure** for 4_1: total dimension about 1e5-1e6 reduced (minutes), 2e6-2e7 unreduced (hours).

## 6. Classification of the 3-strand knots at N = 6 [C counts from n6.txt]

| class | P | count | knots | layer |
|---|---|---|---|---|
| Markov | ≤ 6e4 | words | [1,1,2] … [1,1,1,1,1,2], unknot words | **M6a**, unreduced included |
| A | ≤ 1e7 | 3 | 4_1, 5_2, 6_2 | **M6a**, reduced |
| B | 1e7-2e7 | 7 | 7_3, 8_2, 8_19, 9_3, 8_20, 10_2, 12a_722 | M6a + u32, after S3 (17-29 h wall each) |
| C | 2e7-1e8 | 28 | 6_3, 8_5, 11a_364, 10_124, … | M7a; or per-source streaming (cuts memory towards 2-4x the LB rank) together with τ at the low end. At central τ: about 56 h wall per knot |
| D | > 1e8 | 147 | median P of all 185 knots: 5.2e8 | not planned |
| all | N = 3, 4 | 185 | | **M6a**, one session |
| unreduced N = 6, classes A-B | | | | **M6c** |

## 7. Minimal M6a scope for the first certified 3-strand gl_6 results

**S0. Prerequisites from other notes.**
- M5 with pruning, layers A and B, `restrict(N')` and `pruned_plan` (R2).
- G1(c) at levels 4, 5 and 6.
- The box closure spike (R3).

**S1. sengine core (§6.3-6.4)** with `SOpts.prune_n`:
- `crossing_cx` term truncation, asserting the cell rule;
- streaming on for (3,3,3);
- per stage:
  - K_0 equals the truncated decategorified class;
  - #objects ≥ LB, with actual/LB recorded;
  - degree asserts;
  - the determinism trace;
- `check_d2` only on unpruned runs; the final d∘d is certified on the closed complex (R3: Freivalds, with `closure_kills_map` as the oracle).

**S2. Box driver.**
- Point sets per (prime, node seed). Gram at generator-only points.
- One scan per point, for a batch of words. Each point's outputs go to disk.
- Triangular extraction, `ClosureDataC.valid_n`.
- Staged exact closure, with `close_plan` using the recalibrated estimate.

**S3. Calibration, the first cloud job.**
- Single-point scans of 4_1 at prune 4, 5 and 6, and of [1,1,1,2] at prune 6.
- Record time, RSS, c_e, κ_S, σ, t_lib, M_ev and actual/LB per stage.
- These replace the constants of §4; §5 and §6 are then re-issued.

**S4. u32 `LocMap` storage**, with arithmetic in u64. Only if S3 shows class B is memory-bound.

**Runs, in order (each one gates the next):**
- (a) Regression.
  - (3,3) Hopf, T(2,3) and T(2,5) at N ≤ 6, unreduced and reduced, must equal R1's exact xengine results.
  - Uncolored 3-strand words (T2) must equal engine2; this tests the box path against dense interpolation.
- (b) Prune 4: all 185 knots and the Markov words at N = 3, 4 (duality; Wang at N = 4: ranks 10 and 16).
- (c) The Markov corner of the prune-6 library, then the Markov words at N = 4, 5, 6, unreduced and reduced. This is the **first gl_6 result on 3 strands**.
- (d) G1(c) levels 5 and 6 complete.
- (e) 4_1:
  - N = 6 reduced: the **first new value**;
  - N = 5, unreduced and reduced;
  - N = 4 from the prune-6 scans by nesting;
  - mirror and rotation.
- (f) 5_2 and 6_2 reduced at N = 6, then class B.

**Cost once the library exists:** (b) 1.3 h, (c) 5-8 h wall for the 4 words, (e) about 3.5 h. In total about 2-3 sessions.

## 8. SPEC §6 items for this target

| item | decision | reason |
|---|---|---|
| §6.1 M6a.1 (2 strands), T1 | keep; add (3,3) against R1's exact results | regression of the localized path |
| §6.1 M6a.2 2-colored targets, T8 | **defer** | off the target |
| T10 (cube, cabling) | **drop** | 3-colored 3-strand cables are 9-strand uncolored: they do not fit |
| precondition "G1(a) passed; T0 for every (2,2,2) type" | replace by G1(a)-(3,3), G1(c) and T0 on the kept (3,3,3) types | the (2,2,2) failures are off the path |
| §6.2 interface | keep. Add `prune_n`, `CloseMethod::Box{n}`, node seeds and `valid_n`. Drop `points0`, `interp_cap`, `residual`, `f_mem_bytes`, `time_budget_s` | box path |
| §6.3 PEngine, trace, K_0/LB, degree and memory asserts | keep. K_0 is the truncated class; the memory budget is the per-thread cloud budget, not 512 MiB | |
| §6.3 `check_d2` after every crossing | **change**: unpruned runs only; pruned runs certify the final d∘d on the closed complex at each N' (Freivalds; `closure_kills_map`) | d∘d = 0 only modulo the killed ideal |
| §6.4 Markowitz heap, streaming rule, T5 | keep; streaming on for (3,3,3) | memory is binding |
| §6.5 step 1 (builder, freeze, parallel points) | keep; batch words per point; library from disk | t_lib amortized; 5-h sessions |
| §6.5 step 2 (generators, determinant at every point) | keep, at the box points | guards against degenerate points |
| §6.5 steps 3-4 (dense e-monomial interpolation, `interp_cap`, "needs CL2", rescans, time and F-memory rules) | **drop for all of M6a** | 88 233 unknowns at w = 22; T2 against engine2 tests the box path |
| §6.5 step 5: `DeltaTower`, `Auto`, G4 | **drop** | box data are mod I_N; the tower needs R^λ polynomials |
| §6.5 laptop policy, laptop rows of §6.6, the 768 MB cap | **drop** | the laptop does no compute |
| §6.5 step 6, two primes | keep, on distinct node sequences | |
| §6.6 table | replaced by §5 | |
| T0 | **retarget**: 6 / 20 / 40 kept (3,3,3) types at N = 4 / 5 / 6, killed types → 0, (3,3) at N = 3..8; drop the (2,2,2) parts | Gram, box extraction, weight convention |
| T2, T3, T4, T6 | keep | cheap |
| T6b refusals | **replace**: per-point memory above the budget; predicted wall above one session (split into resumable chunks); `valid_n` mismatch | new scheduler |
| T7 | keep as a regression against R1, N = 3..7 | |
| T9 | **retarget** to the (3,3,3) Markov words at N = 4, 5, 6 | first gl_6 3-strand check |
| §6.8 calibration (10-crossing (2,2,2) word) | replace by S3 | |
| §6.8 G2 ((2,2,2,2) sample) | **defer** | 4 strands are out of scope |
| §6.9 "duality N = 3, weight ≤ 6" | replace by duality at N = 4 and Markov words at N = 6 | |
| unpruned (3,3,3) scans | **defer** | need e_p for the killed 8 820-rank type; strictly larger than prune 6 |

## 9. Acceptance tests and oracles at 3 strands, 3 colors

| oracle | N | what it checks | cost | label |
|---|---|---|---|---|
| wedge^3 gl_3 = det: every knot gives t^0 q^0 | 3 | `global_q` and Rickard-term pruning (only W_0 survives). Trivial for the engine: the pruned complex is 1 object in all 185 scans | none | computed |
| duality: uncolored gl_4 (engine2) of the reversed word, bigraded, unreduced and reduced | 4 | 3-strand engine on the 6 prune-4 types; box closure up to weight 9 (unreduced) or 6 (reduced) | 1 session for all 185 | via duality (expected) |
| Markov words (§5.2), their rotations, [1,1,1,−2] | 4, 5, 6 | prune-6 3-strand scans and the box closure at all weights. Unreduced = Wang (N = 4: 10, 16; N = 5: 46, 106; N = 6: 136, 396); [1,1,2] = `hopf_formula`(6;3,3); reduced = R1 exact | 1-2 h wall per word (N = 6 unreduced) | computed = reference |
| unknot words [1,2], [1,−2] | 4-6 | unreduced = [N choose 3]_q; reduced = 1 | small | computed = known |
| T0, per type | 4, 5, 6 | box Gram gives Φ_N(B_p) = MOY; killed types give 0; includes Gram weight > box top | minutes | |
| nested pruning: prune-6 scans closed at N' = 4, 5 equal the prune-4 and prune-5 runs | 4, 5 | the killed ideal through elimination; box nesting | the N' subsets of the points (16 / 100 reduced) | |
| pruned vs unpruned | 3, 4 (2 strands) | the kill rule categorified (xengine [V]). 3 strands deferred: needs the idempotent | | |
| `decat_chi`; per-stage K_0 and LB equal to `scan_knots3` (4_1 N = 6: LB 276, 38 types) | all | missing or shifted summands | free | |
| two primes on two node sequences, identical bigraded results | all | prime and point dependence; replaces the high-weight residuals | counted in §4 | |
| final d∘d closes to 0; F∘F = 0 in A_N | all | differential modulo the killed ideal | small | |
| symmetries: 4_1 amphichiral, P(t,q) = P(1/t,1/q); mirror word; rotation [−2,1,−2,1]; reversed word (wedge^3 is self-dual at N = 6) | 4-6 | grading, conjugation invariance | free / 2x | |
| Lewark-type relation [N choose 3] P_red − P_unred = (1+t)Q with Q ≥ 0 (as M6.0) | 4, 5; 6 for the Markov words | reduced against unreduced | when both exist | |
| duality: 2-colored gl_5 of the reversed word | 5 | cross-λ check | needs (2,2,2) M6a; optional | via duality (expected) |

## 10. Risks and open checks

- **R1. (3,3,3) Hom provider: coverage and cost unmeasured.** X has up to 31 500 flows, and Search left stuck branches on (2,2,2). G1(c) level 4 comes first, then the Markov corner. M4c is the fallback.
- **R2. Per-point constants.** The μ and τ bands are 20-28x wide; S3 decides. Class B can move up or down by a factor of 2-4.
- **R3. actual/LB.** If the minimal complexes stay above LB, P grows. Measured per stage in S3.
- **R4. Closure at unreduced N = 6.** State spaces reach about 42 000 rows per type. Stage per type and recalibrate `exact_mem_estimate`.
- **R5. Box-point degeneracy.** Mitigation: random nodes, the determinant assert, redraw.
- **R6. Polynomiality of F at high weight** is not tested by the box extraction. It rests on the certified generators, the two node sequences, F∘F = 0 and χ.
- **R7. Idempotents** for the 2 open intermediate-block ladders (R2 §2.4).
- **R8. The kill rule on 3 strands** rests on assumption A1. Its first tests are duality at N = 4, the Markov words at N = 6 and nested pruning.
- **R9. Type coverage.** The Markov words visit only the p ⊗ id corner. 4_1 and T0 cover the other kept types; record each oracle's type coverage.
- **R10. Memory headroom.** The library DAG, the template cache and M_ev per thread are unmeasured. Budget 2 GB plus M_ev.

## 11. Provenance

- `gl6/py/ranks333.py`, `cdiag333.py`: §1, from `scripts/cb.pkl`.
- `counts.py`: emono, boxes and lower sets, §2.
- `r4model.py`: §4-§6, from `n4.txt`, `n5.txt` and `n6.txt`. These hold the N = 4, 5, 6 rows of `ntrunc/t1/out/scan_knots3*.txt` (185 knots; 12n_751 is listed twice).
- `scan_t2n.txt`: the Markov sizes of §5.2.
- `g1a_run1_cold_partial.json`: G1(a) timings.
- C2 H2: κ_S at the 4_1 peak.
- R2 §1: pair counts.

Before S3, run the decategorified scan tool (ntrunc t1, seconds) on the Markov words to replace the §5.2 bounds.
