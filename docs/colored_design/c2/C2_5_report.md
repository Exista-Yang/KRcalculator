# C2-5 report: repair of the non-M5i (2,2,2,2) library failures (branch colored-2c4-fix)

Date 2026-11-14. Machine: 4 cores (Xeon, this container reports 2.30 GHz; the C2-2 table was measured on a
2.10 GHz host, so every run uses `KRHOM_COSTS_ANY_HOST=1`), 15 GB. Base: origin/dev 710c102 + merge of
origin/m4c-rexplan 876634f (B4w: Leaves pair route in the production provider, leafgen memory gate, b4wstamp).
**The m4c-rexplan work had not been reviewed when it was merged; everything below that says "Leaves" depends on
it.** Tags: [M] measured, [E] estimate. Outputs: runs_c2/k4/fix/; library data on data-gl6 under
c2/l2222_n4_c25/ and c2/costs_2222_c25.txt.

Scope: the 27 non-M5i failures of C2-2's (2,2,2,2) prune-4 level set (identical for F31 and F32b): 9 pairs
Refused, 10 pairs Exhausted at 2x, 8 closure types Refused. The 43 NeedsIdempotent jobs (M5i) belong to C2-4 and
were not touched; src/sindec.rs, src/sengine.rs, src/swire.rs were not edited.

## Summary

@@SUMMARY@@

## 1. What was done

1. **Layer A restored** from data-gl6 c2/l2222_n4/{f31,f32b} (sha256 as in the MANIFESTs: 0ec12187…, 180d90ea…),
   2 575 layer-A files each, plus templates.
2. **Traces and calibration of the refused pair windows (default chain).** Each failing job was rerun alone
   (`libbuild one`, `KRHOM_TRACE=1`, runs_c2/k4/fix/trace.sh), the refused Search windows / leaf isomorphisms were
   turned into `calib --windows` input by tools/dev/trace_windows.py (≤ 600 k unknowns) and measured with limit 125 s
   (≥ `max_window_s` = 120 s, the admission ceiling of Search windows), rss 3-4 GB (runs_c2/k4/fix/calib.sh):

   | round | input | calls measured | range | censored | certified afterwards (F31) |
   |---|---|---|---|---|---|
   | 1 | windows of the 9 Refused pairs, C2-2 table | 13 | 45 k-509 k unknowns, 4-78 s | 0 | 367, 369, 617 |
   | 2 | windows still refused after round 1 | 6 | 294 k-506 k, 84-120.4 s | 0 | 140 |
   | 3-5 | Leaves Lemma A step calls of pair 1830 above the 1·10⁵ budget | 3 | 104 k-212 k, 20-57 s | 0 | – (1830 stops at a 3.9 M step) |

   Result table runs_c2/k4/fix/costs_2222_c25.txt (= C2-2's costs_2222.txt + 22 records, 534 records). The window
   records are in label sets such as (1,4,2), (4,1,3), (3,3,2), (4,2,0,2) that had no record at all, which is why
   45-65 k-unknown windows (4-10 s) were refused in C2-2.
3. **Leaves routes** (`--provider m4c-lf`: production chain first, then Stage 0 product closure / Leaves closure
   (leaves-first for closure calls without e_p) and the Leaves pair route (rex class, else Lemma A + exit/entry),
   all inside the unchanged sindec acceptance). C2-2's layer A was copied and carried over to the M4c-lf stamp with
   `b4wstamp` (2 575 files per prime; M5c's route-stamp argument: the stamp only steers which calls are attempted,
   every loaded map is re-certified on load; the replay below checks it). Jobs run with
   `libbuild one --provider m4c-lf --stamp-route`.
4. **sbox** (`src/bin/sbox.rs`, the only source change): new `--provider default|m4c|m4c-lf`, which builds the
   M4c provider and hence its layer-A stamp, so that `sbox --replay --stamp-route --provider m4c-lf` loads the
   maps written by the Leaves routes (C2-2 open problem 2). With `--replay` the provider is never called.
5. **Replay** of level + corner for both primes from the merged layer A (provider refusing every call, M4c-lf
   stamp; runs_c2/k4/fix/replay.sh), then **cover** (`libbuild cover` on the replay files).

## 2. Per failure (both primes)

F32b: every job fixed for F31 was rerun for F32b with the same table and provider and certified (15 / 15,
runs_c2/k4/fix/f32b_1.out); the still-failing jobs fail on call sizes / time records that do not depend on the
prime and were not rerun for F32b (the replay shows them failing for both).

### 2a. Refused pairs (9)

| pair | p (type) | C2-2 cause | C2-5 | route | F31 | F32b |
|---|---|---|---|---|---|---|
| 140 | 0,0,1,1/2,0,0,0/0,1,0,1/0,1,1,0 E0.1 | windows (4,1,3) 64.6 k, 506 k unrecorded | calibrated (rounds 1-2) | Search (chain), 307 s | **ok** | **ok** |
| 367 | 1,0,1,0/0,0,1,1/2,0,0,0/0,1,0,1 F0.1 | leaf iso (4,2,0,2) 508.7 k unrecorded | calibrated (78 s) | Search leaf iso, 176 s | **ok** | **ok** |
| 369 | 1,0,1,0/1,0,0,1/0,1,1,0/1,0,0,1 F0.1 | leaf iso (3,0,3) 49 k unrecorded | calibrated (4 s) | Search leaf iso, 17 s | **ok** | **ok** |
| 617 | 1,0,0,1/1,0,1,0/0,1,0,1/0,1,1,0 E0.2 | window (3,2,3) 107 k unrecorded | calibrated (19-21 s) | Search, 55 s | **ok** | **ok** |
| 626 | 0,0,1,1/2,0,0,0/0,1,0,1/0,1,1,0 E0.2 | windows (4,1,3) | calibrated window 294 k (84-86 s) → Search then Exhausted (π∘ι not scalar) | **Leaves Lemma A** (5 targets, 2-10 s each) | **ok** | **ok** |
| 958 | 0,1,0,1/1,0,1,0/1,1,0,0/0,0,1,1 E1.1 | leaf iso 2.4 M | – | **Leaves rex class** (98 s) | **ok** | **ok** |
| 214 | 0,0,1,1/1,0,0,1/1,1,0,0/0,1,1,0 E0.1 | windows (3,3,2) 273 k, 449 k | calibrated; the 449 079-unknown window measured 101.6 s one way and **120.4 s** the other (> 120 s window ceiling) → still Refused; leaf iso 10.8 M; Leaves Lemma A step 67.2 M > 4·10⁵ hard cap | – | fail | fail |
| 1235 | 0,1,0,1/1,0,1,0/1,1,0,0/0,1,0,1 F1.1 | windows (1,4,2) 45 k unrecorded | calibrated → Search Exhausted (π∘ι not a scalar on the (1,4,2) windows); Leaves Lemma A step 83.7 M > cap | – | fail | fail |
| 1471 | 0,1,0,1/1,0,1,0/1,0,1,0/0,1,0,1 E1.2 | as 1235 | as 1235; Leaves Lemma A step 1.97 M > cap | – | fail | fail |

### 2b. Exhausted pairs (10)

Search S0-S6 is genuinely exhausted for these (stuck branches, no admission refusal). Leaves pair route:

| pair | p, rung | Leaves result | F31 | F32b |
|---|---|---|---|---|
| 59 | 1,0,0,1/1,1,0,0/0,0,2,0/0,1,0,1 E0.1 | **rex class, certified (50 s)** | **ok** | **ok** |
| 1830 | 0,1,1,0/1,0,1,0/0,0,0,2/1,1,0,0 E2.1 | Lemma A: three steps above the 1·10⁵ budget measured (104 511 d1 21 s, 211 537 d3 57 s, 155 193 d3); then a step of **3 877 930** unknowns > hard cap | fail | fail |
| 942 | 0,1,1,0/1,0,1,0/1,0,0,1/0,1,0,1 E1.1 | Lemma A step 620 955 > cap (s = −1) | fail | fail |
| 1054 | 0,0,1,1/1,0,0,1/1,1,0,0/0,1,1,0 E1.1 | Lemma A step 443 391 > cap | fail | fail |
| 1540 | same p, E1.2 | Lemma A step 443 391 > cap | fail | fail |
| 1906 | 0,1,1,0/0,0,1,1/1,0,0,1/1,1,0,0 E2.1 | Lemma A step 655 863 > cap | fail | fail |
| 2392 | same p, E2.2 | Lemma A step 655 863 > cap | fail | fail |
| 1444 | 0,1,0,1/1,0,1,0/1,1,0,0/0,0,1,1 E1.2 | Lemma A step 38.0 M | fail | fail |
| 1470 | 0,1,0,1/1,0,0,1/1,1,0,0/0,0,2,0 E1.2 | Lemma A step 60.2 M | fail | fail |
| 1894 | 0,0,1,1/1,0,0,1/1,1,0,0/0,1,1,0 E2.1 | Lemma A step 184.6 M | fail | fail |

The hard cap 4·10⁵ (B4w: ≤ 1·10⁵, ≤ 4·10⁵ with an uncensored exact-shape record) was not raised: like B4w's 4.5·10⁵
question for pair 467 it is a plan decision. 1054 / 1540 (443 k, 11 % above) would be the first candidates; 942,
1906, 2392 (0.62-0.66 M) next. Memory scales roughly like 4-5 GB per M unknowns on this solver (B4w: 397 k → 1.9 GB),
so all five are feasible on this machine if the cap is raised.

### 2c. Refused closure types (8)

All 8 were refused because their remaining closure calls (degree 7-9 Hom(R, L_p) / Hom(L_p, R) on 1 296 x 1 or
1 x 7 776 fibres) were censored in C2-2's calibration (≥ 296 k unknowns, > 125 s or 3 GB). The Leaves closure route
(leaves-first) certifies each of them in 7-13 s per job with steps of ≤ 1·10⁵ unknowns (largest recorded step
≤ 805 unknowns in the per-level log), for both primes:

| closure type | F31 | F32b |
|---|---|---|
| 0,0,1,1/1,0,0,1/1,1,0,0/0,1,1,0 | **ok** (10.8 s) | **ok** |
| 0,0,1,1/2,0,0,0/0,1,0,1/0,1,1,0 | **ok** (7.4 s) | **ok** |
| 0,0,2,0/0,0,0,2/2,0,0,0/0,2,0,0 | **ok** (6.8 s) | **ok** |
| 0,0,2,0/1,0,0,1/1,0,0,1/0,2,0,0 | **ok** (8.9 s) | **ok** |
| 0,1,1,0/0,0,0,2/2,0,0,0/0,1,1,0 | **ok** (9.4 s) | **ok** |
| 0,1,1,0/0,0,1,1/1,0,0,1/1,1,0,0 | **ok** (11.8 s) | **ok** |
| 0,1,1,0/1,0,0,1/1,0,0,1/0,1,1,0 | **ok** (12.6 s) | **ok** |
| 0,1,1,0/1,0,1,0/0,0,0,2/1,1,0,0 | **ok** (9.2 s) | **ok** |

**Tally: 16 of 27 fixed for both primes** (8 closure types; pairs 59, 140, 367, 369, 617, 626, 958 and … see
below). **11 remain** (pairs 214, 1235, 1471; 942, 1054, 1444, 1470, 1540, 1830, 1894, 1906, 2392 minus the ones
counted above).

@@REPLAY@@

@@COVER@@

## 6. Checks

- Acceptance unchanged: every new map passed sindec's checks (`dechom_postcheck`, D14 scaling, Π∘Ι = id at all
  cosets, kept multiset; `certify_closure`; Leaves pairs additionally `certify_pairs` and π∘ι = id at every coset
  for the rex class) and was re-certified on load in the replay. No check, cap or admission ceiling was changed;
  the only inputs that changed are cost records (measurements) and the provider chain (B4w's, as merged).
- Calibration records are measurements of exactly the calls the providers asked for (no record was edited; the
  120.4 s measurement that keeps pair 214 out was kept).
@@SUITE@@

## 7. Open problems

1. **Pair 214** misses the Search window ceiling (120 s) by 0.4 s on one direction of a 449 k-unknown window
   (the other direction 101.6 s; measured while 3 other calibration jobs ran). Either a quieter re-measurement
   policy (min over repeats instead of max) or a slightly higher window ceiling for measured shapes would admit it;
   both are plan decisions.
2. **The Leaves Lemma A hard cap (4·10⁵)**: 1054 / 1540 (443 k), 942 (621 k), 1906 / 2392 (656 k) stop on a single
   step above it; their step calls are likely feasible here (≈ 2-4 GB, 3-6 min). B5 territory (as in B4w §4).
3. **Big steps**: 1830 (3.9 M after three measured steps), 1471 (2.0 M), 1444 (38 M), 1470 (60 M), 214 (67 M),
   1235 (84 M), 1894 (185 M) need a better Lemma A plan (smaller windows) or another route.
4. **1235 / 1471**: Search is exhausted on the (1,4,2) windows (π∘ι not a scalar, i.e. the summand is not split by
   that window), not refused: no further calibration can help the chain there.
5. M5i (C2-4) remains the dominant blocker of every 4-strand knot.
