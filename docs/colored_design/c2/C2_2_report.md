# C2-2 report: (2,2,2,2) library and the first 4-strand 2-colored gl_4 values (branch colored-2c4)

Date 2026-11-13. Machine: 4 cores, 15 GB. Base: colored-2c4 + origin/dev 202b371 (merge 572255b). Tags: [M]
measured, [E] estimate. Outputs: runs_c2/k4/, docs/colored_design/c2/C2_2_*, docs/colored_design/c2/gl4_2222.json;
library data on data-gl6 under c2/l2222_n4/.

## Summary

- **(2,2,2,2) library, prune 4, measured for both primes.** Level set (1 032 pairs + 125 closure types = 1 157 jobs):
  **1 087 / 1 157 certified (93.9 %) for F31 and for F32b, with the identical 70 failures**; Markov corner
  **144 / 144** for both. Failures by cause (both primes): **43 NeedsIdempotent** (40 pairs + 3 closure types),
  19 pair failures by admission / search (9 Refused, 10 Exhausted at 2x), 8 closure types Refused. Max X certified
  15 552 flows (the maximum of the level set); job wall mean 18-20 s, p90 25 s; peak child RSS 0.35 / 2.1 GB.
  Replay (provider refusing every call, all from layer A): every certified job reproduces, both primes (level
  1 087 / 1 157 = all certified, corner 144 / 144).
- **New blocker: M5i (idempotents).** 7 kept types of (2,2,2,2) at N = 4 (25 over all levels) have a decomposable
  minimal-flow canonical ladder (char L_p = b_p + b_r, one with two extra summands), so `SIndec::ty` returns
  NeedsIdempotent and the engine has no reduced coordinates for them (sengine `lib_shape`). C2-1's census did not
  test this (`idem_needed` is not part of `plan`). **Every knot of knots4_cloud.txt visits them**: class A needs
  18 M5i jobs, class B 30, Conway 37 (`libbuild cover`, section 2). So **no 4-strand knot is computable until M5i
  exists**; library calibration alone cannot fix it.
- **Cover: 0 of 22 knots** fully covered (both primes). Class A is blocked by 29 jobs: 18 M5i + 8 Refused + 3
  Exhausted.
- **First 4-strand (2,2,2,2) gl_4 values** (the words that avoid the M5i types): the Markov words
  σ_1^m σ_2^{±1} σ_3^{±1} (m = 1…5), σ_1σ_2σ_1σ_3, and the stabilized 3-strand knots 4_1·σ_3 and 5_2·σ_3:
  **11 words, all PASS** (χ = decat_chi on all 6 closures, D14 scaling, N' = 3 = engine2 gl_3 of the reversed
  braid, N' = 2 rank 1, F31 = F32b closure by closure). The Markov words equal the 2-strand (2,2) σ^m values exactly at N' = 2, 3, 4,
  reduced and unreduced; 4_1·σ_3 and 5_2·σ_3 equal C2-1's 3-strand (2,2,2) values exactly (N' = 3, 4, both
  modes); σ_1σ_2σ_3 is the colored unknot [4 choose 2]_q. These are the first runs of the 4-strand localized
  pipeline on a (2,2,2,2) library: they validate the library (layer A, replay), the sharded box size (1 296
  points) and the 4-strand closure at N = 4.
- **Measured 4-strand point cost** (step 2, on the largest computable words, not 6_1 which needs M5i): 0.05-0.65 s
  per box point (1 thread) at envelope peak Hom dimension 2.2 k-71 k; VmHWM ≤ 0.41 GB; closure ≤ 26 s. The staged-
  rank memory *estimate* is far too pessimistic at 4 strands (4_1·σ_3: estimate 13.9 GB, real VmHWM 0.29 GB; 5_2·σ_3:
  estimate 55 GB, real 0.41 GB): runs need KRHOM_MEM_CAP_MB raised (ulimit kept as the real guard).
- **Projection** (section 5): once M5i exists, 6_1 / 7_6 / 7_7 cost 1.2-8 core-h unreduced (0.2-1.3 reduced-only)
  per prime; Conway 11-150 core-h unreduced, 2-25 reduced. The spread is the two fits (4-strand words vs 3-strand
  knots) extrapolated 10-200x.
- **Code change**: src/bin/libbuild.rs — a child `pair K` job locates its plan pair without building the whole plan
  (`plan_pair_at`, same enumeration as `pruned_plan`; 23 s → 8 s per child at (2,2,2,2); KRHOM_LIBBUILD_FULLPLAN=1
  restores the old path; identical (p, rung) checked on 4 pairs). New examples: idem_list, idem_end0, ladder_search
  (M5i scoping). Full suite `KRHOM_THREADS=4 cargo test --release`: **282 passed, 0 failed, 36 ignored**
  (runs_c2/k4/full_suite.log).

## 1. (2,2,2,2) library (runs_c2/k4/lib, runs_c2/k4/calib) [M]

Commands (runs_c2/k4/lib/build2222.sh): `KRHOM_COSTS_ANY_HOST=1 libbuild run --lam 2,2,2,2 --n 4 --field F
--cache C --costs T --jobs 4 --wall 900 --rss 3000 [--corner markov-ext]`, then `libbuild replay` with the same set
flags. F31 ran with T = tests/data/homprov/costs_cloud.txt (no (2,2,2,2) records); the Refused closure calls of
that run were then calibrated (`libbuild calib --windows` on the `plan --closure-calls` list of the failing closure
types, limit 125 s, 4 jobs: 88 calls, 34 measured, 54 censored = killed at 125 s or 3 GB) into
runs_c2/k4/calib/costs_2222.txt, and the failed F31 jobs were rerun with it (`runs_c2/k4/lib/fix.sh`, route
stamp). F32b ran directly with costs_2222.txt. The default production provider chain was used throughout (Exact,
Search S0-S6); the m4c / m4c-lf providers were not used (their layer-A stamps are not readable by `sbox`, which
has no `--provider`; and the dominant failure, NeedsIdempotent, is not a provider failure).

| set | prime | jobs | run | after fix (calibrated table) | replay | wall (run) | max X | failures |
|---|---|---|---|---|---|---|---|---|
| level-n4 (1 032 pairs + 125 closure) | F31 | 1 157 | 1 066 (92.1 %) | **1 087 (93.9 %)** | 1 087 (all certified) | 4 992 s | 15 552 | 70 |
| level-n4 | F32b | 1 157 | **1 087 (93.9 %)** | (calibrated from the start) | 1 087 (all certified) | 4 763 s* | 15 552 | 70 |
| corner-markov-ext-n4 (96 cross + 48 closure) | F31 | 144 | 141 | **144** | 144 | 362 s | 1 296 | 0 |
| corner-markov-ext-n4 | F32b | 144 | **144** | – | 144 | 281 s | 1 296 | 0 |

\* F32b wall includes a restart (3 → 4 workers, resumed) and contention with the 4-strand scans.

**G1(c).** Not passed for either prime (70 level failures each), but the two primes agree job by job (the same 70).
Failures by cause and job kind (both primes):

| cause | pairs | closure | note |
|---|---|---|---|
| NeedsIdempotent | 40 | 3 | source or target type has a decomposable canonical ladder (M5i); 22 with a decomposable source, 21 with a decomposable target |
| Refused (admission) | 9 | 8 | Exact over caps (fibre product / unknowns) with no cost record, Search windows refused; pair windows were not calibrated this session |
| Exhausted at 2x | 10 | 0 | Search tiers S0-S6 exhausted (also at 2x budget); M4c list |

Provider tiers (fresh calls, F31 level): S1 1 527, Exact 531, S2 346, mixed 179; seconds per call mean 1.4,
p90 4.2, max 14 (F32b: max 186 s, one 2 GB child). Before the calibration 32 closure types failed (16 Refused, 13
Exhausted, 3 NeedsIdempotent); the 88 calibrated closure calls recovered 21 of them (and the 3 corner ones).

**Cover** (`libbuild cover`, both primes, runs_c2/k4/cover/): **0 of 22 knots** of knots4_cloud.txt. Blocking jobs
after the fix pass:

| class | knots | blocking jobs | of which M5i | Refused | Exhausted |
|---|---|---|---|---|---|
| A | 6_1, 7_2, 7_4, 7_6, 7_7, 8_4, 8_6, 8_8, 8_11, 8_13, 8_14, 9_10, 9_11, 9_13, 9_17, 9_18, 9_20, 9_23 | 29 | 18 | 8 | 3 |
| B | 8_15, 9_22, 11n_42 | 44 | 30 | 9 | 5 |
| C | Conway | 61 | 37 | 14 | 10 |

Also blocked: (σ_1σ_2σ_3)^k (T(4,2), T(4,3) = 8_19, T(4,5)) and other "interleaved" words (runs_c2/k4/cover/
cover_torus_f31_pass1.md). Covered: the Markov / stabilized words of section 3.

**Persistence.** data-gl6 c2/l2222_n4/{f31,f32b} (layer A + templates, ≤ 45 MB chunks, MANIFEST with stamps,
sha256, restore), c2/costs_2222.txt, c2/runs2222/ (all run / fix / replay JSONL and logs). Layer-A files carry
the costs_cloud or the costs_2222 stamp: read them with `--costs costs_2222.txt --stamp-route`.

## 2. M5i: decomposable canonical ladders at (2,2,2,2) [M]

`examples/idem_list 2,2,2,2 4` (runs_c2/k4/idem_2222.txt): 25 types of the visited blocks have a decomposable
minimal-flow ladder; 7 are kept at N = 4 (`idem_needed`):

| type | block | flows of L_p | kill level | char L_p | Exact End^0(L_p) unknowns |
|---|---|---|---|---|---|
| 1,0,1,0/0,1,0,1/1,0,1,0/0,1,0,1 | (2,2,2,2) | 1 296 | **3** | b_p + b_{2,0,0,0/0,0,2,0/0,2,0,0/0,0,0,2} | 41.3 M |
| 1,0,1,0/0,0,1,1/1,1,0,0/0,1,0,1 | (2,2,2,2) | 576 | 4 | b_p + b_r | 14.4 M |
| 0,1,0,1/1,0,1,0/0,1,0,1/1,0,1,0 | (2,2,2,2) | 7 776 | 4 | b_p + b_r | 2.93 G |
| 1,0,1,0/0,1,0,1/1,1,0,0/0,1,0,1 | (2,3,1,2) | 384 | 4 | b_p + b_r | 5.45 M |
| 1,1,0,0/0,1,0,1/1,1,0,0/0,1,0,1 | (2,4,0,2) | 96 | 4 | b_p + b_r | 0.20 M |
| 1,0,1,0/1,0,0,1/0,1,1,0/0,1,1,0 | (2,2,3,1) | 576 | 4 | b_p + b_r | 10.8 M |
| 1,0,1,0/0,0,1,1/1,1,0,0/0,1,1,0 | (2,2,3,1) | 1 152 | 4 | b_p + b_r1 + b_r2 | 63.4 M |

All extra summands sit in degree 0 (multiplicity v^0), the singular analogue of B_{s2}B_{s1}B_{s3}B_{s2} =
B_{s2s1s3s2} ⊕ B_{s2} in S_4. A breadth-first search over ladder words (examples/ladder_search, 40 000 distinct
characters, depth ≤ 10) found no indecomposable alternative ladder for any of the 7 (inconclusive beyond that
depth, but the degree-0 defect suggests none exists). The first type is kept already at N = 3.

What M5i needs (SPEC_M4_M6a §5.2 has the design; status in the code):
1. **e_p.** The spec's route E = provider.end0(L_p) + Krylov P(a) is infeasible here: the whole-ladder End^0 is
   0.2 M-2.9 G Exact unknowns (calibrated Exact limit ≈ 0.2-0.3 M in 125 s). Use instead the decomposition
   machinery that already works: (ι_r, π_r) = provider.pairs(L_p, r, L_r, s = 0, cnt) for each extra summand r
   (the same call `dec_rung` makes for its targets; L_r is indecomposable), normalize by G, orthogonalize, and set
   e_p = id − Σ ι_r π_r; checks e² = e at 3 points, block ranks = KL values, top scalar 1 (spec §5.2 step 6).
2. **sindec.** `ty` stores e_p; `dec_rung` applies the source idempotent (ι ↦ Kron(e_p, R)∘ι, π ↦ π∘Kron(e_p, R))
   and passes e_q for decomposable targets (the provider's `pairs` / `closure` already take `Option<MapId>`
   idempotents and `check_idempotent`); closure keys carry ep (now "ep=none"); layer-A keys change only for these
   types.
3. **engine.** sengine `lib_shape` refuses e_p; swire evaluates map blocks in flow coordinates (`ev.block`,
   Kron indices of L_p ⊗ R). It needs the reduced coordinates P (flows × m), Q (m × flows) per coset
   (`smap::reduced_of`, `DecData::reduced`, currently an identity stub) applied to every block into / out of an
   idempotent type, Kron(P_p, id_R) in the crossing step, and the closure generators composed with e_p.
   This is the substantial part (swire 1.4 k lines, sengine 2.4 k).
Estimated as one session for 1-2 plus one for 3 and its validation (forced non-canonical ladders on (2,2) and
(1,1,1,1) as in spec D2, then 4_1·σ_3 and the Markov words again, then 6_1).

## 3. 4-strand (2,2,2,2) words at N = 4 (runs_c2/k4/markov, tools/c2/k4_markov.py) [M]

Command per (word, prime) (runs_c2/k4/markov/one.sh): `sbox --lam 2,2,2,2 --word W --n 4 --ns 2,3,4 --field F
--seed 201 --prune 4 --replay --cache LIB --costs costs_2222.txt --stamp-route --stages len(W) --mem-mb 3000`
(one-word in-memory box loop, 1 296 points, u and r; χ and the D14 scaling check mandatory; KRHOM_MEM_CAP_MB raised,
see section 4). References: 2-strand (2,2) σ^m xscan of C2-1 (runs_c2/xscan/s{m}_f{F}_p4.json), C2-1's 3-strand
(2,2,2) values (gl4_222.json), engine2 gl_3 of the reversed braid (runs_c2/k4/e2ref), N' = 2 rank 1.

| word | braid | N'=4 u / r | N'=3 u / r | N'=2 | N'=3 vs gl_3 reversed (u / r) | reference | comparison | χ | scaling | primes | status |
|---|---|---|---|---|---|---|---|---|---|---|---|
| m1 | 1 2 3 | 6 / 1 | 3 / 1 | 1 / 1 | eq / eq | colored unknot | u4:eq r4:eq u3:eq | yes | yes | F31 = F32b | PASS |
| m2 | 1 1 2 3 | 36 / 6 | 9 / 3 | 1 / 1 | eq / eq | 2-strand sigma^2 (xscan) | u2:eq r2:eq u3:eq r3:eq u4:eq r4:eq | yes | yes | F31 = F32b | PASS |
| m3 | 1 1 1 2 3 | 22 / 9 | 7 / 3 | 1 / 1 | eq / eq | 2-strand sigma^3 (xscan) | u2:eq r2:eq u3:eq r3:eq u4:eq r4:eq | yes | yes | F31 = F32b | PASS |
| m4 | 1 1 1 1 2 3 | 64 / 18 | 13 / 5 | 1 / 1 | eq / eq | 2-strand sigma^4 (xscan) | u2:eq r2:eq u3:eq r3:eq u4:eq r4:eq | yes | yes | F31 = F32b | PASS |
| m5 | 1 1 1 1 1 2 3 | 46 / 25 | 11 / 5 | 1 / 1 | eq / eq | 2-strand sigma^5 (xscan) | u2:eq r2:eq u3:eq r3:eq u4:eq r4:eq | yes | yes | F31 = F32b | PASS |
| m3n | 1 1 1 -2 3 | 22 / 9 | 7 / 3 | 1 / 1 | eq / eq | 2-strand sigma^3 (xscan) | u2:eq r2:eq u3:eq r3:eq u4:eq r4:eq | yes | yes | F31 = F32b | PASS |
| m3nn | 1 1 1 2 -3 | 22 / 9 | 7 / 3 | 1 / 1 | eq / eq | 2-strand sigma^3 (xscan) | u2:eq r2:eq u3:eq r3:eq u4:eq r4:eq | yes | yes | F31 = F32b | PASS |
| m3x | 1 1 1 -2 -3 | 22 / 9 | 7 / 3 | 1 / 1 | eq / eq | 2-strand sigma^3 (xscan) | u2:eq r2:eq u3:eq r3:eq u4:eq r4:eq | yes | yes | F31 = F32b | PASS |
| u_s3 | 1 2 1 3 | 36 / 6 | 9 / 3 | 1 / 1 | eq / eq | 2-strand sigma^2 (xscan) | u2:eq r2:eq u3:eq r3:eq u4:eq r4:eq | yes | yes | F31 = F32b | PASS |
| w41s | 1 -2 1 -2 3 | 46 / 25 | 11 / 5 | 1 / 1 | eq / eq | 3-strand 4_1 (C2-1) | u3:eq u4:eq r3:eq r4:eq | yes | yes | F31 = F32b | PASS |
| w52s | 1 1 1 2 -1 2 3 | 78 / 49 | 15 / 7 | 1 / 1 | eq / eq | 3-strand 5_2 (C2-1) | u3:eq u4:eq r3:eq r4:eq | yes | yes | F31 = F32b | PASS |

**11 words, 11 PASS, both primes.** (u_s3 = σ_1σ_2σ_1σ_3 is the Hopf link, = 2-strand σ^2.) F32b runs are
3-10x slower per point than F31 here (field arithmetic plus contention with the replays).

Unreduced N' = 4 Poincaré polynomials (full set in docs/colored_design/c2/gl4_2222.json):
- σ_1σ_2σ_3 (unknot): q^-4 + q^-2 + 2 + q^2 + q^4 = [4 choose 2]_q; reduced 1.
- σ_1^3 σ_2 σ_3 (trefoil): rank 22, reduced t^0q^-8 + t^2q^-14 + t^2q^-12 + t^3q^-18 + t^3q^-16 + t^4q^-20 + t^5q^-24 +
  t^5q^-22 + t^6q^-26 (rank 9) — equal to the 2-strand value, for every sign choice of σ_2, σ_3.
- 4_1·σ_3: rank 46 / reduced 25, 5_2·σ_3: 78 / 49 — equal to the 3-strand values (C2-1).

## 4. Measured cost of a 4-strand box point (step 2) [M]

6_1 itself cannot be measured (its library jobs include 18 M5i ones; `sbox` stops while building the library).
Measured instead on the computable 4-strand words, F31, `point_secs_mean` = single-thread seconds per box point:

| word | letters | lib types | envelope peak Hom dim (N=4) | s / point | scan (1 296 pts) | close | VmHWM |
|---|---|---|---|---|---|---|---|
| σ_1σ_2σ_3 | 3 | 27 | 2 197 | 0.051 | 18 s (4 thr) | 5 s | 100 MB |
| σ_1^3σ_2σ_3 | 5 | 27 | 10 309 | 0.064 | 23 s (4 thr) | 7 s | 109 MB |
| σ_1^5σ_2σ_3 | 7 | 27 | 18 421 | 0.170† | 60 s | 8 s | 122 MB |
| σ_1σ_2σ_1σ_3 | 4 | 48 | 7 735 | 0.365† | 127 s | 9 s | 226 MB |
| 4_1·σ_3 | 5 | 48 | 33 397 | 0.499† | 334 s (2 thr) | 14 s | 294 MB |
| 5_2·σ_3 | 7 | 48 | 71 461 | 0.653† | 438 s (2 thr) | 25 s | 408 MB |

† run alongside the F32b library build (contended). Library load from layer A: 21-28 s per run.

**Staged-rank memory.** `cclose::staged_mem_estimate` (dense upper model, scaled by KRHOM_THREADS) refused
σ_1^5σ_2σ_3 (4.7 GB), 4_1·σ_3 (13.9 GB) and 5_2·σ_3 (54.9 GB) under the 4 GB cap; with the cap lifted (KRHOM_MEM_CAP_MB
= 200 000, ulimit -v 14 GB as the real guard) the whole runs peaked at 122-408 MB. At 4 strands the estimate is
30-130x too high; runs_c2/k4/markov/one.sh therefore sets MEMCAP high. (C2-1's 3-strand note — "staged rank needs
1-3.4 GB" — was the estimate too, not the use.) A recalibration of the estimate for 4 strands is worth doing
before large scans.

## 5. Projection for class A / B / C (once M5i and the remaining library jobs exist) [E]

Per-point fit on the 11 measured 4-strand words: point_secs ≈ e^-9.64 · peakdim^0.80 (noisy: 27 vs 48 library
types, contended runs); C2-1's 3-strand fit: e^-11.06 · peakdim^1.04. Both extrapolated to the knots' envelope
(runs_c2/census/env4_2222.txt; 10-200x beyond the measured range):

| knot | peak Hom dim | s / point (4-strand fit … 3-strand fit) | unreduced 1 296 pts, core-h | reduced-only 216 pts, core-h |
|---|---|---|---|---|
| 6_1 | 0.78 M | 3.4 … 22 | 1.2 … 8.0 | 0.2 … 1.3 |
| 7_6 | 1.02 M | 4.2 … 30 | 1.5 … 10.6 | 0.25 … 1.8 |
| 7_7 | 1.31 M | 5.1 … 38 | 1.8 … 13.8 | 0.3 … 2.3 |
| 9_23 (largest A) | 9.1 M | 24 … 291 | 8.7 … 105 | 1.5 … 17 |
| Conway | 12.8 M | 32 … 414 | 11 … 149 | 1.9 … 25 |

Per prime; both primes double it. Memory per point: ≤ 0.4 GB measured at peak dim 71 k; the box driver's 4 shards
fit in 15 GB unless the engine peak grows faster than linearly (unmeasured above 71 k).

**Order of work.** (1) M5i (section 2): two sessions. (2) The 8 Refused + 3 Exhausted non-M5i jobs of class A:
trace their Search windows (`KRHOM_TRACE`, tools/trace_windows.py) and calibrate (`calib --windows`); pairs 958, 1235
(Refused) and 1444, 1470, 1894 (Exhausted) plus 6 closure types (Refused at the calibrated table: censored calls of
degree 8-9, 1 296 x 1 fibres). (3) Scans: 6_1, 7_6, 7_7 unreduced + reduced, both primes (≈ 1-14 core-h each,
`sbox box --shard k/4 --ckpt-secs 600`, runs_c2/k4/scan/box.sh is ready), then reduced-only for the rest of A.

**Conway.** Library: needs the 37 M5i jobs and 24 others (14 Refused, 10 Exhausted — the Exhausted ones also at 2x,
i.e. M4c territory). Scan: 2-25 core-h reduced-only per prime by the projection, 11-149 unreduced. So Conway is
**not reachable this round** (blocked by M5i and the library), but if the 4-strand fit's lower end holds it becomes
a 1-2 session scan after M5i, not the 148 core-h of C2-1's projection; the next measurement (6_1, 0.78 M) decides.

## 6. Checks (mandatory list) [M]

- χ = decat_chi: all 6 closures of every run (sbox exits nonzero otherwise).
- D14 scaling: passed in every run (2 points; 3 340-63 082 entries checked).
- Both primes identical: all 11 words, closure by closure (6 closures each).
- N' = 3 duality: N' = 3 = uncolored gl_3 of the reversed braid (engine2), equal as polynomials, unreduced and
  reduced, all 11 words.
- N' = 2: rank 1 (unreduced and reduced), all 11 words.
- Markov: σ_1^m σ_2^{±1} σ_3^{±1} = 2-strand σ^m at N' = 2, 3, 4 (u and r), m = 1…5 (m = 1: colored unknot).
- Library: replay (provider refusing every call) certifies every job certified by the run, both primes (level 1 087, corner 144).
- The 4-strand knots themselves (6_1, 7_6, 7_7, …): not computed (M5i).

## 7. Open problems

1. **M5i** (section 2) is the blocker for every 4-strand knot at (2,2,2,2), N = 4. The spec's end0 route is
   infeasible at these sizes; the pairs-peeling route reuses certified machinery. The engine part (reduced
   coordinates in swire/sengine) is the larger half.
2. **Library completion outside M5i**: 9 Refused + 10 Exhausted pairs, 8 Refused closure types (identical for both
   primes). Pair windows were not traced/calibrated this session; 54 of the 88 calibrated closure calls were
   censored (≥ 0.25 M unknowns or 3 GB). M4c (Leaves) could help closures, but its layer-A stamp is not readable by
   `sbox` (no `--provider` there): either give `sbox` the provider flag or let `--stamp-route` accept the m4c chain.
3. **Staged-rank memory estimate** is 30-130x too high at 4 strands (section 4); runs must set KRHOM_MEM_CAP_MB.
4. **Projection spread** (×7 for 6_1): one measured 6_1 point after M5i settles it.
5. **No independent N' = 4 reference** at 4 strands (as at 3): the Markov and stabilization identities above are
   the cross-checks; they exercise the (2,2,2,2) library, the 4-strand box (1 296 points) and the closure, but not
   the M5i types that real 4-strand knots need.
