# M5b report: the (3,3,3) library and gate G1(c) (branch colored-m5)

Date 2026-11-11. Machine: 4 cores (Intel Xeon @ 2.30 GHz), 15 GiB. Base bf326b7 (M5a) merged with origin/dev
ce91f0e (M4b hardening) as bb52c19. Measurements before the long runs: M5b_measure.md.
Tags: [M] measured here, [E] estimate.

## Result

**G1(c) is not passed.** Level 4 is 100 % certified for both primes; the extended Markov corner and the large-X
sample are not (same failures for F31 and F32b, nothing uncertified was accepted, no check was weakened).

| set | prime | jobs | certified | targets (Σ mult) certified | closure types | failures | G1(c) |
|---|---|---|---|---|---|---|---|
| level-n4 (42 pairs + 5 closure) | F31 | 47 | **47** | 19 (19) | 5 | - | **pass** |
| level-n4 | F32b | 47 | **47** | 19 (19) | 5 | - | **pass** |
| corner-markov-ext-n6 (16 crossing + 16 closure jobs) | F31 | 32 | 30 | 166 (264) in the 16 crossing jobs | 14 of 16 | 2 closure, Exhausted at 2x | fail |
| corner-markov-ext-n6 | F32b | 32 | 30 | 166 (264) | 14 of 16 | same 2 | fail |
| sample-largex20-n6 | F31 | 20 | 12 | 149 (207) of 188 (255) | - | 7 Refused, 1 Exhausted at 2x | fail |
| sample-largex20-n6 | F32b | 20 | 12 | 149 (207) | - | same 8 | fail |

Beyond the brief (for the projection, and as M5b-run work already done):

| set | prime | jobs | certified | targets (Σ mult) | failures |
|---|---|---|---|---|---|
| level-n5 (158 pairs + 19 closure) | F31 | 177 | 171 | 317 (357) of 324 (364) | 4 Refused, pair 322 and closure 1,1,1/1,1,1/1,1,1 Exhausted at 2x |
| level-n5 | F32b | 177 | 171 | 317 (357) | same 6 |
| sample-random30-n6 (30 of the 286 level-6 provider pairs, seed 2026) | F31 | 30 | 24 | 130 (227) | 4 Refused, pairs 336, 323 Exhausted at 2x |
| sample-random30-n6 | F32b | 30 | 24 | 130 (227) | same 6 |
| level-n6 (full) | F31 | 377 | 309 | 1 363 (2 173) of 1 592 (2 566) | 38 Refused, 23 Exhausted at 2x, 7 Memory at 2x (section "Level 6") |
| level-n6 (full) | F32b | 377 | 309 | 1 363 (2 173) | same 68 |

Every certified job was **replayed** with a provider that refuses every call (`libbuild replay`): all maps load from
layer A and pass the same checks (form at rep(q) at 2 points, genuineness, ends, degrees, D14 scaling, Π∘Ι = id,
kept-multiset assert, rank bookkeeping). Replay = run for every set and both primes (the two 2x-memory pairs 135 and
382 replay at `--rss 4096`).

## What was built

- **src/bin/libbuild.rs** `plan|run|replay|calib` (+ children `one`, `calib-one`):
  - sets: level N (every kept-source pair of `pruned_plan(N)` + closure generators of the kept non-identity types of
    (λ, λ)); `--corner markov-ext` (D8: for every kept type p' ⊗ id_3 of (3,3,3), the σ_1 terms with all Rickard
    components and the σ_2 terms, both signs; closure generators of the kept types reached after σ_2);
    `--sample largex|random --count K [--seed S]`;
  - one child process per job under a watchdog (`ps` RSS poll every 250 ms, wall), `--jobs 4` workers, children with
    KRHOM_THREADS = 4/jobs; JSON lines output (header with args, commit, host, layer-A stamps, watchdog; one record
    per attempt; summary); resume skips jobs certified in OUT, and every child resumes target by target from layer A;
  - G1(c) policy: Time / Memory / Exhausted rerun once at 2x (watchdog and provider Budget doubled); Exhausted at 2x
    goes to the M4c list; Refused is not retried (the budget does not change admission);
  - a pair job must certify exactly its plan targets and Σ mult (else outcome `Count`);
  - `calib`: cost records of (λ) whole-ladder pair calls, or of a list of calls (`--windows FILE`: Search windows parsed
    from `KRHOM_TRACE` choose_window lines, or `plan --closure-calls` lists), each in a child under the watchdog
    (killed = censored record); `--retry-censored` remeasures censored shapes at a longer limit. The table is v2-stamped
    for this host.
- **src/sindec.rs**: `target_log` carries the provider summary (tier hits) and closure calls; `target_log_file`
  appends each entry at once (timings survive a watchdog kill); `provider_for_stamped` at label 3 uses
  `Exact { max_label: 3, max_call_s: EXACT_L3_MAX_CALL_S = 600 }` (stamped). No frozen file was edited.
- **tests/data/homprov/costs_333.txt**: table v3, 273 records (126 censored) measured here: 138 [3,3,3] (whole-ladder
  closure calls and (3,3,3) windows) and 135 windows of the intermediate blocks ([0,3,6], [0,4,5], [1,3,5], [2,3,4],
  [1,4,4], [1,2,6], [4,5], [3,5]).
- **tests/libbuild.rs** (2 tests, 5.4 s): plan counts of the three (3,3,3) sets; a (3,3) N = 5 run → resume (skips
  all) → replay 100 % with no new layer-A file → replay with an emptied layer A fails.
- tools/libbuild_report.py (tables of this report), tools/persist_layer_a.sh (data-gl6 chunks).

Tests: `KRHOM_THREADS=4 cargo test --release --test sindec --test homprov --test libbuild`:
13 + 17 + 2 passed, 0 failed (6 ignored, as before).

## Per set and prime (table v3)

Seconds are provider seconds per fresh provider call (a target's ι/π pair or one closure side), from the
`.targets` logs; job wall is per certified job. F31 ran with the template cache of the calibration (the windows and
closure calls measured there are cached); F32b computed every template itself except the 19 closure templates computed
once with `calib-one --field 32b` (75-420 s each, the same calls F31 has from calibration), so **F32b is the cold cost**.

| set | prime | fresh calls | tiers | s/call mean | p90 | max | job wall sum | max X | peak RSS |
|---|---|---|---|---|---|---|---|---|---|
| level-n4 | F31 | 29 | S1 19, S1+S2 10 | 0.42 | 0.84 | 6.1 | 35 s | 576 | 138 MB |
| level-n4 | F32b | 29 | S1 19, S1+S2 10 | 0.23 | 0.90 | 1.8 | 29 s | 576 | 53 MB |
| corner | F31 | 92 | S2 38, S1 34, S1+S2 16, Exact 2, S0 2 | 0.16 | 0.49 | 1.0 | 80 s | 900 | 212 MB |
| corner | F32b | 92 | same | 0.16 | 0.48 | 0.9 | 74 s | 900 | 212 MB |
| largex | F31 | 151 | S1 96, S2 53, S3 2 | 0.31 | 0.45 | 6.3 | 143 s | 21 600 | 2 749 MB |
| largex | F32b | 151 | same | 0.27 | 0.40 | 5.4 | 122 s | 21 600 | 2 749 MB |
| level-n5 | F31 | 258 | S1 178, S2 64, S1+S2 16 | 0.83 | 0.83 | 56 | 317 s | 12 960 | 896 MB |
| level-n5 | F32b | 258 | same | 2.69 | 1.68 | 104 | 797 s | 12 960 | 896 MB |
| random30-n6 | F31 | 86 | S1 41, S2 33, S3 11, S6 1 | 2.05 | 0.71 | 150 | 249 s | 21 600 | 976 MB |
| random30-n6 | F32b | 86 | same | 7.43 | 0.87 | 559 | 709 s | 21 600 | 976 MB |

Wall per set with 4 workers (F31 / F32b): level 4 10 / 7 s; corner 21 / 19 s; large-X 64 / 56 s; level 5 162 / 283 s;
random level-6 sample 154 / 575 s. The two pairs that need 2x memory (135, 382: X = 12 960, 33 targets each) peak at
2.6-2.7 GB. No Time outcome occurred. No target was found by the Exact backend except 2 corner closure sides; Exact's
role is closure generators (whole-ladder Exact pair calls are 10^8-10^10 unknowns at level 6).

Layer A: about 660 files per prime for these sets (with templates: 2.9-4.4 MB compressed per prime).
Regeneration cost (cold, F32b figures): about 30 min per prime for these sets (calibration templates included:
the 19 closure calls take about 10 min on 4 cores).

## Failures and the M4c list

Same jobs for both primes. Search "Exhausted" means the candidate stream ended (completeness); "Refused" means every
stuck branch was an admission refusal (a window above the caps without a usable cost record).

**M4c list (Exhausted at 2x):**

| job | set | X | what fails |
|---|---|---|---|
| closure 1,0,2/2,0,1/0,3,0 | corner | 600 | top closure degrees: Search finds no generator; the Exact calls of degree >= 15 exceed 3.5 GB / 600 s (censored) |
| closure 0,0,3/3,0,0/0,3,0 | corner | 400 | same, degree 18 (852 011 unknowns) |
| closure 1,1,1/1,1,1/1,1,1 | level 5 | 4 320 | same, degrees >= 6 of the 4 320-flow ladder (406 634 unknowns: > 600 s) |
| pair 336 (0,2,1/2,0,1/1,2,0 ⊗ F1.1) | large-X, random | 14 400 | Search S0-S6 exhausted at degree -1, no admission refusal |
| pair 322 (1,2,0/1,1,1/1,1,1 ⊗ F1.1) | level 5 | 4 320 | Search exhausted at degree 0 |
| pair 323 (0,3,0/2,0,1/1,1,1 ⊗ F1.1) | random level 6 | 2 880 | Search exhausted at degree 0 |

**Refused (admission; also for M4c, since raising the budget does not help):** large-X pairs 45, 145, 292, 293, 392
(X = 21 600), 335 (10 800), 82 (8 640); level-5 pairs 313, 326, 406, 410 (X = 600-3 600); random level-6 pairs 137, 418
(and 292, 293). The remaining candidate windows of the large-X pairs are >= 495 308 unknowns and measured infeasible
(> 120 s at 3.3 GB, or > 3.5 GB in 70-95 s). The level-5 ones are cheap pairs whose windows moved: each calibration
round admits some windows and Search then proposes new ones (e.g. pair 406, X = 600: a [4,5,0] window 40 x 200 with
68 857 unknowns that no round has measured yet). Further `calib --windows` rounds (about 30-60 min each) should shrink
this list; every round changes the table stamp and so invalidates layer A (a rebuild of these sets is about 10 min warm).

A further weakness found here: Exact admission uses the bucket maximum (labels, fibre-product bucket, degree), and
closure calls (fibre product = flows(L_p)) share buckets with pair calls; one censored closure call closes the bucket
for pair calls too (pair 406's whole-ladder Exact call is vetoed as "predicts inf" for this reason). M4c should key the
Exact bucket by call kind or use the exact-shape record first, as Search does.

## Projection for levels 5 and 6 (M5b-run) [M]

- **Level 5 is already built** here, both primes: 171 / 177 jobs, 162 s (F31 warm) and 283 s (F32b cold) wall. What is
  left is the 6 failures above (M4c).
- **Level 6**, from the uniform random sample of 30 of the 286 provider pairs (F32b, cold): 709 s of job wall for 24
  certified pairs plus the failures; per provider pair about 25 s of job wall, with a heavy tail (one pair 560 s).
  286 pairs x 25 s / 4 workers ≈ **30 min**, plus 52 identity/all-killed pairs (seconds), plus 39 closure jobs: on
  level 5 the 19 closure jobs cost about 2 min, but level-6 closure types reach 4 320-10 500 flows and Exact closure
  calls of 100-600 s each [M, calibration], so **0.5-1.5 h** [E from the measured calls] → **about 1-2 h wall per
  prime**, with a failure rate near the sample's 20 % of provider pairs (Refused and Exhausted) and several closure
  types expected in the M4c list. The full level-6 runs of this session (section "Level 6") replace this estimate:
  26 min per prime, 82 % of the jobs.
- Memory: 2 GB per child suffices except for 2 of 50 measured large pairs (2.6-2.7 GB at 2x); 4 workers fit in 15 GB.

## Level 6 (full, both primes) [M]

`libbuild run --n 6` (338 pairs + 39 closure jobs), table v3, 4 workers, wall 600 s, RSS 2 048 MB, then the
G1(c) retries; afterwards one extra pass of the failures at `--rss 6500 --wall 1200` (no change).

| prime | jobs | certified | pair targets (Σ mult) | closure types | wall | failures |
|---|---|---|---|---|---|---|
| F31 | 377 | 309 (82.0 %) | 1 363 (2 173) of 1 592 (2 566) | 27 of 39 | 1 536 s (+ the killed first attempt, see below) | 38 Refused, 11 pairs + 12 closure Exhausted at 2x, 7 Memory at 2x |
| F32b | 377 | 309 (82.0 %) | 1 363 (2 173) | 27 of 39 | 1 587 s | identical list |

Tiers of the 847 fresh provider calls: S2 398, S1 378, S3 46, Exact 10, S6 6, S1+S2 8, S6+Exact 1. Provider seconds
per call: mean 1.4, p90 0.8-0.9, max 214-224 s; certified job wall mean 4.7-5.0 s, p90 8 s, max 227-235 s; peak RSS of a
certified job 4.3-4.4 GB (2x). The failure lists of F31 and F32b are identical.

Level-6 failures:
- **Exhausted at 2x (M4c):** pairs 71, 159, 268, 322, 323, 325, 334, 336, 361, 368, 419; closure types 0,3,0/1,0,2/2,0,1,
  1,0,2/2,0,1/0,3,0, 1,2,0/0,0,3/2,1,0, 0,0,3/3,0,0/0,3,0, 0,3,0/0,0,3/3,0,0, 1,1,1/1,1,1/1,1,1, 0,2,1/2,0,1/1,1,1,
  1,0,2/1,2,0/1,1,1, 1,1,1/0,2,1/2,0,1, 1,1,1/1,0,2/1,2,0, 0,1,2/2,0,1/1,2,0, 0,2,1/1,0,2/2,1,0.
- **Memory at 2x:** pairs 195, 196, 199, 444, 445, 460, 475. Several have tiny X (360-600 flows) but exceed 4 GB (and
  5.5-6.8 GB at `--rss 6500`): the memory goes into a window or whole-ladder exact call, not into X. For M4c.
- **Refused (admission):** pairs 24, 27, 37, 45, 66, 78, 82, 127, 137, 145, 163, 212, 219, 227, 273, 278, 281, 283, 285,
  292, 293, 313, 319, 326, 332, 335, 376, 381, 383, 385, 392, 406, 410, 411, 418, 455, 458, 467.

Incident: the first F31 level-6 run died after 39 jobs because the libbuild binary was rebuilt under it
(`current_exe()` became "… (deleted)", spawn failed, the parent panicked). Fixed in 1274015 (the parent strips
" (deleted)" and records a spawn failure as an outcome); the rerun resumed from layer A and the run file.

Replay of level 6 (provider refusing every call, `--rss 4096`): 309 / 377 for both primes, the same jobs (206 s / 187 s);
during it the kernel OOM-killed one child of an uncertified job (4 replay children of the Memory pairs at once).

So the projection above is replaced by measurement: **level 6 costs about 26 min wall per prime** on this machine
with a warm level-5 cache (pairs ≈ 25 min, closure jobs included), and certifies 82 % of the jobs (85.6 % of the pair
targets, 84.7 % of Σ mult). What M5b-run has left is not compute but the 68 failures per prime, i.e. M4c.

## Persistence (data-gl6)

data-gl6/m5b/{f31,f32b}: `layerA_tmpl.tar.gz.part-000` (layer A + template library through level 6, one chunk each, 14-16 MB; 1 469 layer-A files per prime) with
MANIFEST.txt (λ, field, code commit, layer-A stamps incl. the cost-table hash, file counts, sha256, restore command),
the run and replay logs (runs/*.jsonl, *.targets) and costs_333.txt. Restore:
`cat layerA_tmpl.tar.gz.part-* | tar -C CACHE -xzf -` then `libbuild run|replay ... --cache CACHE --costs
tests/data/homprov/costs_333.txt` (the table must be byte-identical, its hash is in the stamp).

## Open problems

1. The M4c list above (6 Exhausted jobs: 3 closure types whose top-degree closure generators Search does not find and
   Exact cannot compute within 3.5 GB / 600 s; 3 pairs where Search S0-S6 is exhausted).
2. 13 distinct pairs Refused by window admission. Large-X: windows >= 495 k unknowns are beyond the exact solver here.
   Small-X (level 5): more calibration rounds, or Search proposing the cheapest windows first.
3. Exact bucket admission mixes closure and pair calls (see above).
4. The cost table in the layer-A stamp makes calibration and building sequential; a stamp of the provider *route*
   (chain and caps) without the table hash would keep layer A valid across calibration rounds, since every loaded map
   is re-certified anyway. Not changed here (cache-key policy, D14).
5. F31's template cache contains every calibration template; a cached template is admitted by Search whatever the table
   says. To keep the primes on the same route, F32b got the same closure templates computed once (`calib-one`); a
   cleaner rule would be "admit only by the table" or record the template set in the stamp.
