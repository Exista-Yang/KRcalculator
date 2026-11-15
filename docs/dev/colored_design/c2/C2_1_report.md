# C2-1 report: 2-colored (∧² on every component) gl_4 homology, 2-, 3- and 4-strand (branch colored-2c4)

Date 2026-11-13. Machine: 4 cores, 15 GB. Base origin/dev 7c64a78; **no source file of the crate was changed**
(frozen files untouched); new: run scripts (runs_c2/), analysis tools (tools/c2/), docs (docs/dev/colored_design/c2/),
library data (data-gl6/c2/). Tags: [M] measured, [E] estimate.

## Summary

- **2-strand sanity** (step 1): (2,2) σ^3 and σ^5 at N' = 2, 3, 4, both primes, pruned (prune 4) and unpruned:
  all 24 unreduced closures equal tests/data/wang_ref.txt (blocks 2 N m, N = 2, 3, 4, m = 1, 2); every run
  certified with a passing base-point check. σ^5 at N = 4: rank 46 unreduced, 25 reduced.
- **(2,2,2) library, prune 4** (step 2): `libbuild run` level-n4 **105 / 105** jobs (90 pairs, 15 closure types;
  201 targets, Σ mult 229) and the Markov corner (identity closure included) **21 / 21**, for **F31 and F32b**;
  replay (provider refusing every call) 100 % for both. **G1(c) passed for both primes.** Wall 29 s (F31) / 35 s
  (F32b) with 4 workers, max X 432, peak child RSS 125 MB, no retries. The unpruned (level 6) library was built
  too: 146 / 146 + 21 / 21, both primes (104 s / 148 s). Persisted: data-gl6/c2/ (98e79b8).
- **3-strand 2-colored gl_4** (step 3): see section 3 (numbers filled at the end of the session).
- **4-strand census** (step 4): docs/dev/colored_design/c2/C2_1_census.md. The (2,2,2,2) level-4 library is
  1 032 pairs (869 provider pairs) + 125 closure types, max X 15 552; every 4-strand knot needs 450-846 pairs and
  93-125 closure types of it. The truncated envelope keeps Gram entries tiny (≤ 196 unknowns), so the costs
  are the library build (unknown certification rate) and the scan (projected 8-210 core-h unreduced per knot).
  Verdict: **GO for measuring the (2,2,2,2) level-4 library next**; unreduced scans of the smallest knots
  (6_1, 7_6, 7_7) are one session per prime with the box driver, reduced-only for the rest; Conway no-go this round.
- Full suite `KRHOM_THREADS=4 cargo test --release`: **282 passed, 0 failed, 36 ignored** (runs_c2/full_suite.log).

## 1. 2-strand sanity (runs_c2/xscan/, tools/dev/xscan_vs_wang.py) [M]

`xscan 2,2 "σ^m" 2,3,4 --both --prime P [--prune 4]` for m = 3, 5 (P = 31, 32b), plus m = 2, 4 at prune 4 for the
Markov check. All 24 (m = 3, 5) unreduced closures equal Wang (k = 2): ranks 1 / 7 / 22 (σ^3) and 1 / 11 / 46 (σ^5)
at N' = 2 / 3 / 4; reduced at N' = 4: 9 (σ^3), 25 (σ^5). Pruned = unpruned, F31 = F32b. Each run < 1 s.

## 2. (2,2,2) library (runs_c2/lib/, data-gl6/c2/) [M]

Commands (runs_c2/lib/build222.sh): `KRHOM_COSTS_ANY_HOST=1 libbuild run|replay --lam 2,2,2 --n 4 --field F
--cache C --costs tests/data/homprov/costs_cloud.txt --jobs 4 --wall 900 --rss 3000`, then the same with
`--corner markov-ext`. Provider: the default production chain (the M4b-calibrated cost table
costs_cloud.txt, which has (2,2,2) records); m4c / m4c-lf were not needed (no failure to fix).

| set | prime | jobs | certified | replay | wall (run) | max X | retries |
|---|---|---|---|---|---|---|---|
| level-n4 (90 pairs + 15 closure) | F31 | 105 | **105** | 105 | 29.1 s | 432 | 0 |
| level-n4 | F32b | 105 | **105** | 105 | 35.1 s | 432 | 0 |
| corner-markov-ext-n4 (12 cross + 9 closure) | F31 / F32b | 21 | **21 / 21** | 21 / 21 | 2.9 / 2.8 s | 36 | 0 |
| level-n6 = unpruned (126 pairs + 20 closure) | F31 / F32b | 146 | **146 / 146** | 146 / 146 | 104 / 148 s | 432 | 0 |
| corner-markov-ext-n6 | F31 / F32b | 21 | 21 / 21 | 21 / 21 | 5.5 s | 36 | 0 |

Note: the level set contains closure generators of the kept *non-identity* types only; the identity type
2,0,0/0,2,0/0,0,2 (the final type of 4_1, for example) comes from the Markov corner set. Without it sbox stops with
"final type … has no closure generators: replay: closure not in layer A". So the library for knots is
level ∪ corner.

`libbuild cover --lam 2,2,2 --n 4` with both run files: **185 of 185 knots of knots3_all.txt fully covered** (upper
bound; runs_c2/census/cover3_222.md).

Persistence: data-gl6/c2/l222_n4/{f31,f32b} (prune 4) and l222_n6_unpruned/{f31,f32b}, one chunk each (≤ 18 MB),
MANIFEST.txt with stamps, commit, sha256 and the restore command; runs/ holds all run/replay JSONL and logs.

## 3. 3-strand (2,2,2) at N = 4 (runs_c2/k3/, results/colored_gl4/gl4_222.json) [M]

Command per (knot, prime) (runs_c2/k3/one.sh): `sbox --lam 2,2,2 --word W --n 4 --ns 2,3,4 --field F --seed 201
--prune 4 --mem-mb 3200 --replay --cache LIB --costs costs_cloud.txt --stages len(W)` (box N = 4: 216 points
unreduced, 36 reduced; one-word in-memory box loop; scaling check and χ mandatory), 4 jobs × 1 thread, queue
ordered by the truncated envelope's peak Hom dimension (runs_c2/k3/order.txt).

**34 knots computed (both primes), 34 PASS, 0 failures.** Braid lengths 4-12; 68 runs (34 knots × 2 primes).
Every run: F31 = F32b closure by closure; χ = `schur::decat_chi` on all 6 closures (chi_ok); the D14 scaling
record present and passing (2 points); base-point check in the engine. All checks held for all 34 knots.

**N' = 3 duality check.** wedge^2 C^3 ≅ (C^3)^* ⊗ det, so the (2,2,2) gl_3 invariant should equal the uncolored
gl_3 homology of the *reversed* braid (engine2, examples/e2ref), up to the global shift the engine already applies.
It does, as full polynomials (zero residual shift) for all 34 knots, unreduced and reduced. N' = 2: wedge^2 C^2 is
the determinant (1-dimensional), so the colored invariant is that of the unknot = rank 1; this holds for all 34,
unreduced and reduced.

**N' = 4 is the new value.** wedge^2 C^4 is self-dual, so there is no duality shortcut; these N' = 4 polynomials
(results/colored_gl4/gl4_222.json) are the genuinely new 2-colored gl_4 invariants. Unreduced ranks 46 … 2 206,
reduced 25 … 2 025; the first 15 knots (by length) below, the full table in gl4_222_summary.md:

| knot | word | N=4 unred | N=4 red | N=3 unred | N=3 red | N=3 vs gl_3 reversed (u / r) | χ | scaling | primes | status |
|---|---|---|---|---|---|---|---|---|---|---|
| 4_1 | 1 -2 1 -2 | 46 | 25 | 11 | 5 | eq / eq | yes | yes | yes | PASS |
| 5_2 | 1 1 1 2 -1 2 | 78 | 49 | 15 | 7 | eq / eq | yes | yes | yes | PASS |
| 6_2 | 1 1 1 -2 1 -2 | 166 | 121 | 23 | 11 | eq / eq | yes | yes | yes | PASS |
| 6_3 | 1 1 -2 1 -2 -2 | 222 | 169 | 27 | 13 | eq / eq | yes | yes | yes | PASS |
| 7_3 | 1 1 1 1 1 2 -1 2 | 222 | 169 | 27 | 13 | eq / eq | yes | yes | yes | PASS |
| 7_5 | 1 1 1 1 2 -1 2 2 | 358 | 289 | 35 | 17 | eq / eq | yes | yes | yes | PASS |
| 8_2 | 1 1 1 1 1 -2 1 -2 | 358 | 289 | 35 | 17 | eq / eq | yes | yes | yes | PASS |
| 8_5 | 1 1 1 -2 1 1 1 -2 | 526 | 441 | 43 | 21 | eq / eq | yes | yes | yes | PASS |
| 8_7 | -1 -1 -1 -1 2 -1 2 2 | 622 | 529 | 47 | 23 | eq / eq | yes | yes | yes | PASS |
| 8_9 | 1 1 1 -2 1 -2 -2 -2 | 726 | 625 | 51 | 25 | eq / eq | yes | yes | yes | PASS |
| 8_10 | -1 -1 -1 2 -1 -1 2 2 | 838 | 729 | 55 | 27 | eq / eq | yes | yes | yes | PASS |
| 8_19 | 1 1 1 2 1 1 1 2 | 90 | 29 | 19 | 11 | eq / eq | yes | yes | yes | PASS |
| 8_20 | 1 1 1 -2 -1 -1 -1 -2 | 118 | 81 | 19 | 9 | eq / eq | yes | yes | yes | PASS |
| 8_21 | 1 1 1 2 -1 -1 2 2 | 286 | 225 | 31 | 15 | eq / eq | yes | yes | yes | PASS |
| 9_3 | 1 1 1 1 1 1 1 2 -1 2 | 438 | 361 | 39 | 19 | eq / eq | yes | yes | yes | PASS |

**Structure.** reduced rank(N'=4) = reduced rank(N'=3)^2 for 27 of the 34 knots (the thin ones); the 7 exceptions
are 8_19 (29 vs 121), 10_124 (61 vs 289), 10_125 (125 vs 169), 10_139 (133 vs 441), 10_161 (173 vs 361),
12n_235 (181 vs 361), 12n_242 (101 vs 529) — informational, not a check. The unreduced N' = 4 Poincaré polynomial
of every amphichiral knot in the set is mirror-symmetric P(t,q) = P(1/t,1/q) (e.g. 4_1).

**Checks against the library.** Pruned (prune 4) = unpruned (level-6 library) on 4_1, 5_2, 6_2, 6_3, 7_3, both
closure-by-closure, all six closures each (runs_c2/k3/unpruned/, tools/dev/same.py). Markov words σ_1^m σ_2^{±1}
(m = 2, 3, 4, 5) on (2,2,2) equal the 2-strand σ^m xscan of (2,2) at N' = 2, 3, 4 (reduced and unreduced), both
primes, all 12 words (tools/c2/markov.py): this is the destabilization invariance of the σ_2 strand.

**Cube / cable reference (4_1, colors 2).** Both infeasible here within 5 GB: `krhom --colored-method cube` ran
> 90 min without finishing (it had not left the cube-assembly stage); `--colored-method cable` (6 strands) ran out
of memory (OOM at ~3 GB). The N' = 3 duality reference (engine2) and the pruned-vs-unpruned agreement are the
independent cross-checks that stand in for it.

**Cost (prune 4, box N = 4, 216 points unreduced + 36 reduced, 1 thread, 4 jobs on 4 cores).** 45 s (4_1) to
37 min (12a_1214, unreduced rank 2 206) per (knot, prime); VmHWM ≤ 0.3 GB except a handful of 10-12 letter knots
that need the staged-rank memory cap raised (KRHOM_MEM_CAP_MB, --mem-mb; set to 4500 / 3200 in runs_c2/k3/one.sh).
The queue was ordered by the truncated envelope's peak Hom dimension (runs_c2/k3/order.txt), a good cost proxy
(section 4 fit). The remaining 151 of 185 knots are the longer braids; the projection is 150 core-h for all of them.

## 4. 4-strand census (docs/dev/colored_design/c2/C2_1_census.md) [M + E]

See the census for the tables. Measured: `libbuild plan` (2,2,2,2) at N = 4 and 8; `libbuild cover` for the 22
knots of knots4_cloud.txt (labels 2,2,2,2); the truncated decategorified envelope (tools/dev/ntrunc_t1) for the 22
knots (χ truncated = untruncated for all 22, and for all 185 3-strand knots); per-X-bucket outcomes of the (3,3,3)
level-6 build for the library projection; the 3-strand per-point scan cost fit for the scan projection.

## 5. Open problems

1. **3-strand coverage is partial by time, not by capability.** 34 of 185 knots computed this session (all braid
   lengths up to 12 among them), 0 failures. The (2,2,2) library is complete (level 4 + corner, both primes, replay
   100 %); the remaining 151 knots are only a matter of compute (≈ 150 core-h, projected), and are resumable per
   (knot, prime) with runs_c2/k3/queue.sh (order.txt). No mathematical or engine obstruction was seen.

2. **No cube/cable reference at colors 2 here.** Both the full-cube and the 6-strand cable methods exceed 5 GB /
   90 min on 4_1. The N' = 3 duality (engine2 reversed braid) and pruned = unpruned are the cross-checks for N' = 4;
   an independent N' = 4 reference would need either a bigger machine for the cube or a dedicated colored reference
   implementation. Worth adding before trusting 4-strand N' = 4 values that have no duality shortcut.

3. **Staged-rank memory at N = 4.** A few 10-12 letter 3-strand knots need the staged-rank step's memory cap raised
   (1-3.4 GB observed) above the 768 MB default; set with KRHOM_MEM_CAP_MB / --mem-mb. At 4 strands this step will be
   larger (the box is Gr(2,4)^4, 1 296 points, peak complexes 3-30× the 3-strand ones): expect 1-4 GB per point, so
   2-3 parallel points per 15 GB machine, not 4.

4. **(2,2,2,2) library certification rate is unmeasured.** The 4-strand go/no-go rests on a projection from the
   (3,3,3) level-6 per-X-bucket rates, which is pessimistic (the (2,2,2) builds here were 100 %, labels ≤ 4 vs ≤ 6).
   The cheap decisive next step is `libbuild run --lam 2,2,2,2 --n 4` for one prime (≤ ~1 h projected), then `cover`.

5. **4-strand scan projection is the softest number.** The fit (point_secs ∝ peakdim^1.04 over 45 3-strand runs,
   residual ×0.3-3) is extrapolated one-to-two orders of magnitude; a single measured 4-strand box point (sbox
   `point`, once the library exists) would replace it. Unreduced Conway (148 core-h projected) is out of reach this
   way; reduced-only (216 points) is the realistic 4-strand unreduced-free target for the hard knots.

6. **Duality was used only as a check, not a shortcut.** wedge^2 C^4 self-dual means the N' = 4 values are new; the
   wedge^2 C^N for N = 2, 3 are dual/trivial and only serve as references. There is no analogue of the (3,3,3)
   "N = 3 = trivial" collapse at N = 4 for (2,2,2).
