# C2-1 report: 2-colored (∧² on every component) gl_4 homology, 2-, 3- and 4-strand (branch colored-2c4)

Date 2026-11-13. Machine: 4 cores, 15 GB. Base origin/dev 7c64a78; **no source file of the crate was changed**
(frozen files untouched); new: run scripts (runs_c2/), analysis tools (tools/c2/), docs (docs/colored_design/c2/),
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
- **4-strand census** (step 4): docs/colored_design/c2/C2_1_census.md. The (2,2,2,2) level-4 library is
  1 032 pairs (869 provider pairs) + 125 closure types, max X 15 552; every 4-strand knot needs 450-846 pairs and
  93-125 closure types of it. The truncated envelope keeps Gram entries tiny (≤ 196 unknowns), so the costs
  are the library build (unknown certification rate) and the scan (projected 8-210 core-h unreduced per knot).
  Verdict: **GO for measuring the (2,2,2,2) level-4 library next**; unreduced scans of the smallest knots
  (6_1, 7_6, 7_7) are one session per prime with the box driver, reduced-only for the rest; Conway no-go this round.
- Full suite `KRHOM_THREADS=4 cargo test --release`: **282 passed, 0 failed, 36 ignored** (runs_c2/full_suite.log).

## 1. 2-strand sanity (runs_c2/xscan/, tools/c2/xscan_vs_wang.py) [M]

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

## 3. 3-strand (2,2,2) at N = 4 (runs_c2/k3/, docs/colored_design/c2/gl4_222.json) [M]

Command per (knot, prime) (runs_c2/k3/one.sh): `sbox --lam 2,2,2 --word W --n 4 --ns 2,3,4 --field F --seed 201
--prune 4 --mem-mb 3200 --replay --cache LIB --costs costs_cloud.txt --stages len(W)` (box N = 4: 216 points
unreduced, 36 reduced; one-word in-memory box loop; scaling check and χ mandatory), 4 jobs × 1 thread, queue
ordered by the truncated envelope's peak Hom dimension (runs_c2/k3/order.txt).

RESULTS_PLACEHOLDER

## 4. 4-strand census (docs/colored_design/c2/C2_1_census.md) [M + E]

See the census for the tables. Measured: `libbuild plan` (2,2,2,2) at N = 4 and 8; `libbuild cover` for the 22
knots of knots4_cloud.txt (labels 2,2,2,2); the truncated decategorified envelope (tools/ntrunc_t1) for the 22
knots (χ truncated = untruncated for all 22, and for all 185 3-strand knots); per-X-bucket outcomes of the (3,3,3)
level-6 build for the library projection; the 3-strand per-point scan cost fit for the scan projection.

## 5. Open problems

OPEN_PLACEHOLDER
