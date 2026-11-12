# M5c report: shrinking the (3,3,3) failure list; knot coverage at N = 5, 6 (branch colored-m5)

Date 2026-11-12/13. Machine: 4 cores (Intel Xeon @ **2.10 GHz**; M5b ran on a 2.30 GHz box of the same class), 15 GiB.
Base 9fd8850 (M5b) merged with origin/dev as d2f910b. Warm start: layer A + templates of data-gl6/m5b/{f31,f32b}
restored per MANIFEST. Tags: [M] measured here, [E] estimate.

## Summary

- **Level 5: 172 / 177 jobs for both primes (M5b: 171).** **Level 6: 330 / 377 (F31) and 329 / 377 (F32b) (M5b: 309
  each).** Replay (provider refusing every call) certifies exactly the same jobs for every level and prime. No check
  was weakened; every change is either result-preserving (solver, memory, sharing), a classification fix, an admission
  rule (which calls are tried) or calibration data.
- Fixed: all 7 Memory pairs (a template-library memory bug), 4 of the 11 "Exhausted" pairs (they were leaf-isomorphism
  *admission* refusals, now reported as such and calibrated), 7 Refused pairs at level 6 and pair 406 at levels 5 and 6
  (calibration rounds), 3 level-6 closure types (resume after their templates were cached); the exact solver is about
  2x faster with identical output.
- **Knot coverage: no 3-strand knot is fully computable at N = 5 or N = 6, before or after.** Every knot needs (lower
  bound) every level-5 failure, and at least 38 of the 47-48 level-6 failures. N = 5 for all 185 knots needs exactly
  5 more jobs: pairs 313, 322, 326, 410 and closure 1,1,1/1,1,1/1,1,1.

## Commits (colored-m5)

| commit | step | what |
|---|---|---|
| d2f910b | - | merge origin/dev |
| c203964 | 1 | `libbuild cover`, `sindec::knot_jobs` (UB / K_0 LB), M5c_coverage.md, tests/knotcov.rs |
| 35db8e2 | 2 | Exact cost buckets by call kind, exact-shape record first (`predict_call`), test |
| c5da2d5 | 3 | optional layer-A route stamp (`--stamp-route`), test |
| a31c6f8 | 4 | window SBims shared between templates (Memory pairs), test |
| 431e48a | 5 | leaf-iso admission refusals are Refused, `[leaf_iso]` trace, tools/trace_windows.py, `KRHOM_HS_TIMING`, test |
| 0312b5f | 5 | solver: equation dedup + heap sparse sweep in `kernel_lazy` (identical results), hom_lazy test |
| f21e6e6 | 3 | cost table tests/data/homprov/costs_333_m5c.txt (351 records) |
| (this) | 6 | M5c_report.md, coverage after M5c |

data-gl6: 57f0471 (F31 snapshot after the first level-6 run), fd418f2 (final F31 + F32b, logs, README).

## 1. Knot coverage (M5c_coverage.md) [M]

`libbuild cover` / `sindec::knot_jobs` list the library jobs a pruned localized scan of a braid word needs: plan pairs
visited rung by rung through every Rickard term kept by the cell rule, and the closure types of the final chain
objects; as an upper bound (no cancellation) and a lower bound (types with nonzero K_0 class, which survive any
Gaussian elimination; checked against `decat_scan`). For the 185 three-strand knots of knots3_all.txt:

- **No knot is covered at N = 5 or N = 6 with the M5b library, not even by the lower bound.** 4_1 (4 crossings)
  already visits every kept type of (λ, λ) after 4 crossings and 89 of the 158 level-5 provider pairs.
- N = 5: each of the 6 level-5 failures is needed (LB) by **every** knot. N = 6: 57 of the 68 failures are needed by
  every knot, the other 11 by 104-181 knots.
- Consequence: knot coverage is all-or-nothing in the failure list; the engineering target is the list itself.

After M5c (section 6): with the M5c run files of both primes: still 0 of 185 at N = 5 and at N = 6 (UB and LB). N = 5: the 5 remaining
level-5 failures are each needed by every knot. N = 6: 38-48 blocking failures per knot (LB), 38 of them universal.

## 2. Exact admission by call kind (M5b open problem 3)

`CostTable` buckets are keyed by (labels, call kind, fibre-product bucket, degree); kind 0 = one-sided (closure calls
Hom(R, L_p), Hom(L_p, R)), kind 1 = two-sided (pair and window calls). `Exact::admit` uses `predict_call`: the record
of exactly this call shape first, else the bucket of its kind. Test `homprov::m5c_exact_admission_by_call_kind`.
Effect on the failure list: none by itself [M] — the whole-ladder pair calls of the refused pairs are 2-900 M unknowns
(e.g. pair 406: 2.2 M), far beyond any record; the fix removes a wrong veto, it does not make those calls feasible.

## 3. Refused windows: calibration rounds, route stamp

**Host stamp.** This container is a 2.10 GHz host; the M5b table (2.30 GHz stamp) is rejected here unless
`KRHOM_COSTS_ANY_HOST=1` (the first traces of this session ran without a table by accident; noticed from the stamp
`costs=rejected`, 4 layer-A files that had been overwritten were restored from the tarball). All runs below use
`KRHOM_COSTS_ANY_HOST=1`.

**Route stamp** (M5b open problem 4): `SIndecOpts::layer_a_route_stamp` / `libbuild --stamp-route` accepts a layer-A
file whose header equals the current one up to the cost-table token (`costs=ALG:RECS:HASH` compared as `costs=ALG`).
Soundness does not depend on it: every loaded map is re-certified on load exactly as a fresh one (count, ends,
degrees, D14 scaling; Π∘Ι = id, form and genuineness in `dec_one`; `certify_closure`), which replay already relies
on; the cost table only steers which calls a provider attempts. Solver generation, chain, caps, slowdown, ladder and
template stamps and the key must still match; files are written with the full stamp; default off. Test
`sindec::m5c_layer_a_route_stamp`. This let calibration and building interleave without invalidating layer A.

**Search order.** `choose_window` already proposes the cheapest admissible window (minimum unknowns over all
admissible decomposable windows of the commutation class); refusals happen only when *no* window is admissible. The
lever is therefore the table (and the solver speed), not the order.

**Rounds** (`tools/trace_windows.py` turns `KRHOM_TRACE` refusals — `choose_window` and the new `leaf_iso` lines — into
`calib --windows` input; limit 125 s ≥ `max_window_s` = 120 s, 3 jobs, `--retry-censored`):

| round | input | calls | completed | censored (125 s) | pairs certified afterwards (trace) |
|---|---|---|---|---|---|
| 1 | windows ≤ 260 k unknowns of the 4 + 38 refused pairs | 66 | 32 | 34 | 406 (N = 5, 6), 24, 27, 66 |
| 2 | windows ≤ 330 k of the 49 still failing (incl. leaf-iso windows), after the solver speed-up | 60 | 41 | 19 | 163, 268, 325, 361, 368, 392, 455 |
| 3 | windows ≤ 200 k of the 45 still failing | 7 | 2 | 5 | (see G1(c) below) |
| closure | Hom^6(R, L_p), Hom^6(L_p, R), p = 1,1,1/1,1,1/1,1,1 (406 634 unknowns) | 2 | 2 (262 s, 342 s) | 0 | closure reaches degree 8 (710 k: no) |

Result table `tests/data/homprov/costs_333_m5c.txt` (351 records; costs_333.txt unchanged for M5b reproduction).
Smallest blocking window of every still-refused pair after round 3: 227 (158 774), 212 (172 997), 127 (180 200), 45 / 145 (309 704), 458 (321 709), 319 / 335
(355 749), 411 (366 280), 383 (407 587), 278 / 293 (424 126), 273 (476 880), 78 (482 275), 219 (483 355; Time at 2x
for F31), 313 / 326 (516 446), 285 (601 752), 385 (622 068), 292 (701 771), 37 (828 843), 137 (910 974), 332
(958 629), 410 (1.26 M), 82 (1.43 M), 283 (1.56 M), 418 (1.90 M), 336 (2.02 M), 381 (2.32 M), 376 (8.78 M), 281
(19.9 M); leaf isos 159 / 323 (735 960), 322 (5.46 M), 71 (6.58 M), 334 / 419 (146.6 M). Pair 467 (107 034) runs out
of time at 2x (1 200 s). Windows above about 250 k unknowns do not complete in 125 s on this machine even with the
faster solver.

## 4. Memory failures (pairs 195, 196, 199, 444, 445, 460, 475) [M]

Diagnosis (gdb on an allocation failure under `ulimit -v`): the memory went into `TmplLib::push`, which stored a
**clone of both whole-ladder SBims (and their fixed points) in every template** of a hom-space call; the whole-ladder
Exact pair calls of these pairs return hom spaces of dimension up to 42, i.e. 42 copies of two SBims of 240-7 200
flows. Fix: one `Arc<Side>` (SBim, fixed points, positions) per window in `TmplLib::sides`, shared by all its
templates; data and localizations unchanged. Test `tmpl::m5c_tests` (sharing) and m4a A4 (localizations vs oracle).

| pair | X | M5b | M5c (F31, 1 thread) |
|---|---|---|---|
| 195 | 360 | Memory at 2x (> 4 GB; 5.5-6.8 GB at --rss 6500) | ok, 163 MB, 194 s |
| 196 | 360 | Memory | ok, 60 MB, 5 s |
| 199 | 240 | Memory | ok, 106 MB, 70 s |
| 444 | 360 | Memory | ok, 185 MB, 9 s |
| 445 | 360 | Memory | ok, 85 MB, 100 s |
| 460 | 600 | Memory | ok, 193 MB, 35 s |
| 475 | 7 200 | Memory | ok, 174 MB, 125 s |

The hom_space equations and `kernel_lazy` were not the problem; nevertheless the solver got two result-preserving
speed-ups while measuring (section 5b).

## 5. Exhausted pairs and closure generators

**5a. The Exhausted pairs were admission failures.** `Search::leaf_iso` returned Exhausted when every candidate
isomorphism window between a leaf word and the canonical ladder of q was refused by admission. All 11 level-6
"Exhausted at 2x" pairs (71, 159, 268, 322, 323, 325, 334, 336, 361, 368, 419; 322 also at level 5) fail exactly
there. Now reported as Refused (with the cheapest refused window, traced for calibration). Test
`sindec::m5c_leaf_iso_refusal_is_engineering`. Window sizes of the iso: 151 774 unknowns (268, 325, 361, 368 —
certified after round 2), 735 960 (159, 323), 5.4 M (322: (4,4,1) E0 F1 F0(2) F1 ↔ F0 E0 F1(2) F0), 6.6 M (71),
146.6 M (334, 336, 419: (4,5,0)/(5,4,0) words of 2 400 flows).

**5b. Faster exact solver, identical results.** `hom_space` drops scalar-multiple equation rows (normalize, dedup:
about 1/3 of the rows); `kernel_lazy` visits the accumulator's nonzero columns through a min-heap instead of a dense
scan between the row's first and last column (same pivots, order and arithmetic). The reduced kernel basis depends
only on the row space and the column order, so every result is unchanged; tests hom_lazy (property test vs
`Echelon::kernel`, `hom_space` vs `hom_space_reference`, new m5c case). (4,5,0) 68 857 unknowns: 14.3 s → 7.5 s
(kernel 11.8 → 4.6 s). The closure call Hom^6(R, L_{1,1,1/1,1,1/1,1,1}) went from > 600 s (censored) to 262 s.
`COST_ALG` stays `hs-kl1`: old records are upper bounds for the faster solver.

**5c. Closure 1,1,1/1,1,1/1,1,1 (level 5) and the top-degree closure generators.** Search finds 3 of the 4 degree-6
generators of Hom^6(R, L_p); the Exact closure call is per side and all-or-nothing over the c_diag degrees 4..18:
degree 6 is now feasible (262 / 342 s, 3.2-3.7 GB), degree 8 has 710 441 unknowns and degrees 10-18 1.2-7.5 M. So
the job still fails, now at degree 8. Not implemented (for M4c): (i) a hybrid closure provider that takes Search's
generators where it finds them and asks Exact only for the missing degrees; (ii) the duality of closure generators
(ins of p ↔ outs of the transposed type via the anti-involution of the Schur algebroid; p = p^T here) — this needs a
"flip" map node and its certification in smap, which does not exist yet; (iii) degree reduction by dots: a generator
of degree d+2 as (dot) ∘ (generator of degree d) is what Search's S4 tries; its failure at degree 6 suggests that the
missing generator is not a dot multiple.

## 6. G1(c) rebuild, levels 5 and 6, both primes [M]

`libbuild run` (route stamp, table costs_333_m5c, 4 workers, wall 600 s, RSS 2 048 MB, G1(c) retry at 2x), then
`libbuild replay` (provider refusing every call, RSS 4 096 MB).

| set | prime | jobs | M5b certified | M5c certified | M5c replay | run wall (s) | failures M5c |
|---|---|---|---|---|---|---|---|
| level-n5 | F31 | 177 | 171 | **172** | 172 | 58 | 4 Refused, 1 closure Exhausted at 2x |
| level-n5 | F32b | 177 | 171 | **172** | 172 | 285 | same |
| level-n6 | F31 | 377 | 309 | **330** | 330 | 2 234 + 2 413 (resume) | 36 Refused, 2 Time at 2x (219, 467), 9 closure Exhausted at 2x |
| level-n6 | F32b | 377 | 309 | **329** | 329 | 1 820 + 2 186 (resume) | 38 Refused (the 36 of F31 + 24, 219), 1 Time at 2x (467), 9 closure Exhausted at 2x |

New at level 5: pair 406. New at level 6 (F31): 18 pairs — the 7 Memory pairs, 4 ex-Exhausted (268, 325, 361,
368), 7 ex-Refused (24, 27, 66, 163, 392, 406, 455) — and 3 closure types; F32b the same without pair 24.

The level-6 run was done twice per prime: the first pass certified 327 (F31) / 325 (F32b); a resume pass (only
uncertified jobs) certified 3 more closure types (1,0,2/2,0,1/0,3,0, 0,3,0/1,0,2/2,0,1, 1,2,0/0,0,3/2,1,0) and, for
F32b, pair 455 (after the three F31-only calibration templates had been computed for F32b, see below). The reason for
the closures: in the first pass the closure ins side was found and saved at 2x, the outs side failed; templates
computed in that attempt are cached and cached templates are always admitted, so the resume found the outs side at
1x in 0.4 s. This is legitimate (each map is certified) but shows that results depend on the order of work.

**Prime difference (M5b open problem 5).** F31's template cache holds the calibration templates, F32b's does not:
three windows that completed in 120-125 s (above `max_window_s` = 120 s, so refused by the table) were cached for F31
only; they were computed for F32b with `calib-one --field 32b` (100-106 s each), as M5b did for the closure
templates. Pair 24 still differs (F31 ok, F32b Refused): Search's window choice prefers cached windows (cost 0),
so F31 takes a different decomposition route through a cached (2,4,3) window, while F32b's route ends in a (3,3,3)
window of 857 625 unknowns. Pairs 45, 145 and 219 fail for both primes (F31 at 2x, F32b at 1x).

Costs of the M5c pass: level 6 per prime about 37 min (first pass) + 40 min (resume, dominated by the 1 200 s 2x
retries of 219 and 467); level 5 1-5 min; replays 6-8 min.

## Remaining failure list by cause

Level 5 (both primes): pairs 313, 326 (Search windows of 516 446 unknowns), 410 (1.26 M), 322 (leaf iso, 5.46 M),
closure 1,1,1/1,1,1/1,1,1 (Exact closure needs degree 8: 710 441 unknowns; Search finds 3 of 4 generators at degree 6).

Level 6 (F31; F32b additionally pair 24):

| cause | jobs |
|---|---|
| Refused, Search window 150-500 k unknowns (calibration at a longer limit might admit some) | 127, 212, 227, 45, 145, 458, 319, 335, 411, 383, 278, 293, 273, 78 |
| Refused, Search window > 500 k unknowns | 313, 326, 285, 385, 292, 37, 137, 332, 410, 82, 283, 418, 336, 381, 376, 281 |
| Refused, leaf isomorphism window 0.7-147 M unknowns | 159, 323, 322, 71, 334, 419 |
| Time at 2x (1 200 s) | 219, 467 |
| Exhausted at 2x: closure generators at top degrees (Exact closure call censored / > 600 s) | 0,0,3/3,0,0/0,3,0, 0,3,0/0,0,3/3,0,0, 1,1,1/1,1,1/1,1,1, 0,2,1/2,0,1/1,1,1, 1,0,2/1,2,0/1,1,1, 1,1,1/0,2,1/2,0,1, 1,1,1/1,0,2/1,2,0, 0,1,2/2,0,1/1,2,0, 0,2,1/1,0,2/2,1,0 |

No Memory and no completeness (Search exhausted without admission refusal) failure remains among the pairs.

## Persistence

data-gl6/m5c/{f31,f32b}: `layerA_tmpl.tar.gz.part-000` (one chunk each, 16.7 / 13.9 MB; 1 542 / 1 540 layer-A files,
1 283 / 1 122 templates), MANIFEST.txt (stamps, commit, sha256, restore command with `--stamp-route`), runs/ (run and
replay logs of levels 5 and 6), m5c/costs_333_m5c.txt. The M5b data (data-gl6/m5b) is unchanged.

## Open problems

1. **Large windows** (≥ ~250 k unknowns at 120 s per window): the refused pairs whose cheapest window is 300 k-20 M
   unknowns, and the iso windows of 0.7-147 M. Options: a higher per-window limit with a matching job wall; a
   cheaper decomposition of those windows (they all contain a label ≥ 4 merge block such as (4,5,0), (5,1,3), (6,0,3),
   (1,5,3)); a better exact solver (the kernel is still a semi-dense elimination; the systems are 15-25x
   overdetermined; structured/blocked elimination by the bidegree of the unknowns was not tried).
2. **Leaf isomorphisms**: 5 pairs need an iso between two words of one indecomposable whose differing window is
   5-147 M unknowns. Iso through intermediate words (iso_reach) fails within 400 words; a rex-move (braid relation)
   library for divided powers of sl_3 type in the (λ) labels would give these isos without solving.
3. **Closure generators at top degrees** (section 5c): hybrid Search + Exact per degree, duality, and the 12 level-6
   closure types.
4. The cost table mixes measurements of two solver generations and two clock speeds (both conservative directions are
   documented); a clean recalibration of the whole table with the faster solver would admit more windows (most
   censored records are from the slower solver under 3-way contention).
5. Coverage: every knot needs all of the level's failures; only the complete list unlocks knots (section 1).
