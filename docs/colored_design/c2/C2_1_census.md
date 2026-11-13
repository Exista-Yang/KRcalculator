# C2-1 census: (2,2,2,2) at N = 4 (2-colored gl_4, 4-strand knots)

Session C2-1, branch colored-2c4, 2026-11-13. Machine: 4 cores, 15 GB. No provider calls were made for
(2,2,2,2); everything here is plan counting, the truncated decategorified envelope, and projections from
the 3-strand (2,2,2) measurements of this session and the (3,3,3) level-6 build (M5c). Tags: [M] measured,
[E] estimate / projection.

## 1. Library plan (libbuild plan, `pruned_plan`) [M]

| λ | N | pairs | kept-source pairs | provider pairs (non-identity, with targets) | kept targets (Σ mult) | killed targets | closure types | max X | plan time |
|---|---|---|---|---|---|---|---|---|---|
| (2,2,2) | 4 | 126 | 90 | 67 | 201 (229) | 56 | 15 | 432 | < 1 s |
| (2,2,2) | 6 (unpruned) | 126 | 126 | 98 | 358 (404) | 0 | 20 | 432 | < 1 s |
| (2,2,2,2) | 4 | 2 520 | 1 032 | 869 | 2 675 (2 990) | 1 514 | 125 | 15 552 | 23 s |
| (2,2,2,2) | 8 (unpruned) | 2 520 | 2 520 | 2 191 | 9 255 (10 322) | 0 | 281 | 51 840 | 23 s |

Pruning at N = 4 keeps 41 % of the (2,2,2,2) pairs (1 032 / 2 520), 29 % of the targets by count, 44 % of the
closure types (125 / 281) and cuts max X from 51 840 to 15 552 flows (runs_c2/census/plan_2222_n{4,8}.{out,json}).

X distribution (provider pairs / closure jobs, flows; tools/c2/xdist.py):

| set | (0,100] | (100,500] | (500,1k] | (1k,2k] | (2k,5k] | (5k,10k] | (10k,20k] | max |
|---|---|---|---|---|---|---|---|---|
| (2,2,2) N=4 pairs | 61 | 6 | 0 | 0 | 0 | 0 | 0 | 432 |
| (2,2,2,2) N=4 pairs | 176 | 377 | 77 | 112 | 94 | 24 | 9 | 15 552 |
| (2,2,2,2) N=4 closure | 31 | 59 | 5 | 25 | 2 | 3 | 0 | 7 776 |
| (3,3,3) N=6 pairs (M5c, for comparison) | 57 | 70 | 62 | 74 | 30 | 29 | 8 (+8 > 20k) | 21 600 |

## 2. Coverage per knot (libbuild cover, `knot_jobs`) [M]

The library jobs a pruned scan of each knot visits (upper bound; lower bound = types with nonzero K_0 class),
runs_c2/census/cover4_2222.txt, with the X of the needed jobs (runs_c2/census/coverx4.txt):

| knot class | knots | jobs UB pairs + closure (LB) | provider pairs needed by X bucket (as above) | closure by X bucket |
|---|---|---|---|---|
| A: σ_1σ_2σ_3 words with 7-11 letters | 6_1, 7_2, 7_4, 7_6, 7_7, 8_4, 8_6, 8_8, 8_11, 8_13, 8_14, 9_10, 9_11, 9_13, 9_17, 9_18, 9_20, 9_23 | 450 + 93 (414-450 + 93) | 107 / 132 / 20 / 30 / 25 / 3 / 1 (318) | 30 / 41 / 3 / 18 / 0 / 1 (93) |
| B | 8_15, 9_22, 11n_42 | 555 + 93 (507-555 + 93) | 121 / 170 / 31 / 50 / 45 / 4 / 2 (423) | same 93 |
| C | Conway (11n_34) | 846 + 125 (828 + 117) | 161 / 301 / 58 / 78 / 66 / 18 / 6 (688) | 31 / 59 / 5 / 25 / 2 / 3 (125) |

As for (3,3,3) (M5c), a knot visits almost the whole level set of its class after a few crossings, so
coverage is all-or-nothing in the library's failure list: one knot of a class needs essentially the same
library as every other knot of that class. Even 6_1 needs the X = 15 552 pairs.

## 3. Truncated decategorified envelope (tools/ntrunc_t1 `files ... nodecat`, labels 2,2,2,2) [M]

runs_c2/census/env4_2222.txt (all 22 knots of knots4_cloud.txt; 1.6-6.5 s each). Per knot, the maximum over
the stages of the K_0 lower bound on objects (LB obj) and types, the peak object count, the peak Hom
dimension of the peak complex, and the Gram-entry sizes: "box" = unknowns of one Gram entry restricted to the
Gr(2,4)^4 box (unreduced / reduced), against the unboxed full entry of the untruncated scan. χ of the truncated
K_0 class = untruncated χ at N = 4 for all 22 knots (and for all 185 three-strand knots, env3_222.txt).

| knot | letters | LB obj N=4 (N=inf) | LB types N=4 (N=inf) | peak obj N=4 (N=inf) | peak Hom dim N=4 | Gram box max unknowns u / r (N=inf full u / r) | Σ Gram box |
|---|---|---|---|---|---|---|---|
| Conway | 11 | 4225 (11029) | 120 (241) | 29553 (103573) | 12828709 | 196 / 38 (159952 / 13377) | 403619 |
| 6_1 | 7 | 789 (1298) | 94 (152) | 2567 (5724) | 778699 | 196 / 38 (94523 / 9009) | 94387 |
| 7_2 | 9 | 1617 (3120) | 94 (153) | 7123 (16972) | 2287375 | 196 / 38 (94523 / 9009) | 191935 |
| 7_4 | 9 | 1771 (3354) | 94 (153) | 6979 (17420) | 2220679 | 196 / 38 (94523 / 9009) | 212173 |
| 7_6 | 7 | 1203 (1908) | 94 (152) | 3833 (7826) | 1022431 | 196 / 38 (94523 / 9009) | 146883 |
| 7_7 | 7 | 1351 (1860) | 94 (151) | 4417 (8754) | 1314271 | 196 / 38 (94523 / 9009) | 163489 |
| 8_11 | 9 | 2433 (4046) | 94 (153) | 8351 (20130) | 2526439 | 196 / 38 (94523 / 9009) | 295387 |
| 8_13 | 9 | 2693 (4300) | 94 (152) | 10633 (22298) | 3371851 | 196 / 38 (94523 / 9009) | 316533 |
| 8_14 | 9 | 2565 (4152) | 94 (152) | 9007 (20806) | 2774887 | 196 / 38 (94523 / 9009) | 302981 |
| 8_15 | 9 | 2399 (4042) | 94 (152) | 8835 (20908) | 2576329 | 196 / 38 (94523 / 9009) | 267965 |
| 8_4 | 9 | 2115 (3696) | 94 (153) | 9367 (20370) | 2929807 | 196 / 38 (94523 / 9009) | 251551 |
| 8_6 | 9 | 1889 (3268) | 94 (152) | 6345 (14484) | 1787287 | 196 / 38 (94523 / 9009) | 228249 |
| 8_8 | 9 | 2389 (4168) | 94 (153) | 14295 (28452) | 3959023 | 196 / 38 (94523 / 9009) | 286089 |
| 9_10 | 11 | 3341 (6164) | 94 (153) | 12935 (32100) | 3954379 | 196 / 38 (94523 / 9009) | 406303 |
| 9_11 | 9 | 2261 (3822) | 94 (152) | 7599 (16568) | 2029147 | 196 / 38 (94523 / 9009) | 274831 |
| 9_13 | 11 | 4039 (7412) | 94 (153) | 16723 (39878) | 5368459 | 196 / 38 (94523 / 9009) | 469729 |
| 9_17 | 9 | 3057 (4848) | 94 (153) | 10777 (24242) | 3152083 | 196 / 38 (94523 / 9009) | 374925 |
| 9_18 | 11 | 4319 (7808) | 94 (153) | 16689 (40358) | 5261707 | 196 / 38 (94523 / 9009) | 508305 |
| 9_20 | 9 | 3103 (5386) | 94 (153) | 18711 (36852) | 4892983 | 196 / 38 (94523 / 9009) | 376257 |
| 9_22 | 9 | 3663 (5600) | 94 (151) | 12731 (26128) | 3257419 | 196 / 38 (94523 / 9009) | 429727 |
| 9_23 | 11 | 4713 (8670) | 94 (153) | 30921 (67416) | 9148747 | 196 / 38 (94523 / 9009) | 537389 |
| 11n_42 | 13 | 9981 (30347) | 94 (153) | 42585 (184303) | 17787529 | 196 / 38 (94523 / 9009) | 1028725 |

For comparison, (2,2,2) 3-strand at N = 4 (env3_222.txt): 4_1 LB obj 69, 16 types, peak 205 obj, peak Hom
dimension 6 721, Gram box 38 / 8 unknowns; the largest 3-strand knots reach a peak Hom dimension of 6.5 M.
Box points: 6^4 = 1 296 (unreduced) and 6^3 = 216 (reduced) for 4 strands, against 216 / 36 for 3 strands.

The Gram systems are tiny (≤ 196 unknowns per entry): the cost is the localized scan (complex size), not the
Gram solve.

## 4. Projections [E]

**Library build (2,2,2,2) level 4, per prime.** With the (3,3,3) level-6 per-X-bucket outcomes and walls
(M5c F31 run, tools/c2/xfail.py: 100 / 99 / 92 / 91 / 73 / 62 / 75 % certified in the buckets above; mean
certified wall 0.6-22 s), the 869 provider pairs cost about 5 000 s of job wall, i.e. **about 25 min with 4
workers**, plus the G1(c) 2x retries of the failures (up to 20 min each) → **0.5-5 h per prime**. The same
rates predict **about 57 pair failures** (A: 14, B: 23, C: 41 of the needed ones) and closure failures in the
X ≥ 2 000 buckets (the (3,3,3) level-6 closure types with X ≥ 2 000 all failed; the 25 (2,2,2,2) closure types in
1 000-2 000 have no (3,3,3) analogue). These rates are a pessimistic transfer: the (2,2,2) builds of this
session were 100 % (level 4: 105 / 105, unpruned level 6: 146 / 146, both primes, about 30-150 s), labels of
(2,2,2,2) ladders at N = 4 are ≤ 4 against ≤ 6 for (3,3,3) at N = 6, and M4c (Leaves, B4w) has since lifted
(3,3,3) level 6 to 345 / 377. Only a measured build decides.

**Scan cost per knot.** Power-law fit of the measured (2,2,2) per-point scan time (1 thread, 4 jobs on the
4 cores) against the envelope's peak Hom dimension: point_secs ≈ e^-7.17 · peakdim^0.66 over the first 32 runs
(residuals ×0.4-4; tools/c2/project.py, runs_c2/census/projection.txt). Extrapolated (two orders of magnitude,
and without the 94 vs 16 types per object, so treat as a lower estimate within ×10):

| class | per point | unreduced (1 296 points) | reduced only (216 points) |
|---|---|---|---|
| A (6_1 … 9_23) | 6-30 s | 2-11 core-h | 0.4-1.8 core-h |
| B (8_15, 9_22, 11n_42) | 13-47 s | 5-17 core-h | 0.8-2.8 core-h |
| C (Conway) | 38 s | 14 core-h | 2.3 core-h |

Memory: the 3-strand runs peak at ≤ 0.3 GB (VmHWM) at 12 letters; the 4-strand peak complexes are 3-30× larger
in objects, so 1-4 GB per point is the expected range [E]; 2-4 parallel points fit in 15 GB.

## 5. Go / no-go for 4 strands at N = 4

| class | library | scan (if library complete) | verdict |
|---|---|---|---|
| A: 6_1, 7_x, 8_x (8_4, 8_6, 8_8, 8_11, 8_13, 8_14), 9_x (9_10, 9_11, 9_13, 9_17, 9_18, 9_20, 9_23) | 450 + 93 jobs (subset of level 4), max X 15 552; projected 0-14 failures | 2-11 core-h unreduced per knot and prime | **GO for a measured library build** (next step); scans GO if the build certifies the class (1-3 knots per session per prime unreduced, all reduced) |
| B: 8_15, 9_22, 11n_42 | 555 + 93 jobs | 5-17 core-h | conditional GO (after A's library; 105 more pairs) |
| C: Conway | 846 + 125 jobs (nearly the full level 4) | 14 core-h unreduced, 2.3 reduced | **no-go this round**; reduced-only after the full level-4 build |
| any knot, unpruned (N = 8 library) | 2 191 pairs, max X 51 840 | - | no-go (pruning is required; pruned = unpruned verified on 3 strands) |

The decisive unknown is the (2,2,2,2) library certification rate; it is cheap to measure (≤ 1 h per prime for a
first pass). Recommended first step of the next session: `libbuild run --lam 2,2,2,2 --n 4` for F31 (then F32b),
`libbuild cover` per class, then 6_1 through sbox (box driver, sharded, `--ckpt-secs 600`).
