# M2S′ step 0: measurements of the exact 2-strand engine on (3,3) σ^n

Date 2026-11-11, cloud (4 cores, 15 GB), branch colored-2s (dev eb6d8d5 + colored-m4b 402b06a, i.e.
with 31f9c13 sparse final kernel), F31, `xscan 3,3 "1 ..." 3,4,5,6`, KRHOM_THREADS=4.
Baseline = the unchanged engine (full build of C ⊗ X, `PMat` store) plus measurement code only
(StageMem). Bytes are estimates of heap use of the stored components (32 B per `PVec` entry,
32 B per term (u128 Mono + Fp u64), vector headers); RSS from /proc (VmHWM reset per stage).

## Baseline per stage (σ^4 run; stages 1-3 = the σ^3 run, which gave the same numbers)

| crossing | objects before → after (LB) | comps before | entries before | terms before | MB before | terms after | MB after | caches MB | RSS peak MB | s decomp | s transport | s elim | s d∘d |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 4 → 4 (4) | 3 | 219 | 300 | 0.0 | 300 | 0.0 | 0.1 | 20 | 0.0 | 0.36 | 0.00 | 0.00 |
| 2 | 189 → 15 (13) | 1 538 | 156 154 | 541 611 | 26 | 11 438 | 0.7 | 81 | 121 | 48.3 | 1.8 | 0.05 | 0.02 |
| 3 | 928 → 40 (24) | 23 974 | 2 159 218 | 7 604 737 | 480 | 96 344 | 6.6 | 81 | 606 | 0.0 | 20.2 | 3.2 | 0.38 |
| 4 | 2 599 → 85 (33) | 131 048 | 11 604 504 | 48 810 056 | 3 015 | 497 062 | 35 | 81 | 3 285 | 0.0 | 145.0 | 57.6 | 4.2 |
| 5 (σ^5) | ≈ 5 300 (est.) | — | — | — | — | — | — | 81 | **killed at 12 012 MB** | 0.0 | > 560 (still in transport) | — | — |

- Closure data: 0.58 s (σ^3), 0.90 s (σ^4); closures at N = 3..6 < 0.1 s. σ^3 total 74.6 s,
  σ^4 total 283.9 s, peak RSS 0.61 / 3.29 GB.
- σ^5: crossing 5 started at ≈ 290 s; RSS grew linearly at ≈ 21 MB/s during the transport of
  d_C ⊗ id and hit the 12 GB watchdog at ≈ 850 s, before any elimination.
- Hom spaces: 144 calls, 46.7 s, max 2.3 s, all in crossings 1-2 (all 4 (3,3) types exist after
  crossing 1). The sparse final kernel (31f9c13) changes nothing visible here: σ^3 takes 74 s /
  0.61 GB, as in run H with the dense kernel (70 s / 0.6 GB). The memory problem is not the solver.

## Where the bytes are

- **Stored components of the unsimplified C ⊗ X dominate**: 480 MB of 606 MB at crossing 3,
  3 015 MB of 3 285 MB at crossing 4; caches (decompositions incl. composite ι/π, dloc) are a
  constant 81 MB. Transients (kron_id products) are not resident: RSS = components + caches.
- **Poly terms dominate the components**: 3.5-4.2 terms per entry, about 62 B per term in total.
  The minimal complex is 1-2 % of the unsimplified one (3 015 MB → 35 MB at crossing 4).
- Growth per crossing (unsimplified bytes): ×18 (2→3), ×6.3 (3→4); objects ×4.9, ×2.8 and the
  bytes per object ×2.3 (0.52 → 1.16 MB) because the polynomial degrees grow.
- **Time: transport dominates** from crossing 3 on (70 % of crossing 4: 145 of 208 s);
  elimination (full Markowitz rescan per pivot) is 28 %. This meets the step-4 condition
  (transport > 50 % of the time), so rung-by-rung transport is built too.

## After steps 1, 2, 4 (same machine, σ^3 and σ^4 to completion)

Streaming by C-degree (step 1), compact u32 store (step 2), rung-by-rung transport (step 4);
the σ^4 row ran next to the σ^5 baseline and three wang.py jobs (CPU shared).

| run | crossing | objects inserted / held at once | comps held | terms held | MB held | RSS peak MB | s transport | s elim |
|---|---|---|---|---|---|---|---|---|
| stream only, σ^3 | 3 | 928 / 492 | 10 248 | 3 537 803 | 214 | 324 | 16.6 | 3.0 |
| stream + compact, σ^3 | 3 | 928 / 492 | 10 248 | 3 537 803 | 39 | 134 | 17.0 | 2.3 |
| stream + compact + rung transport, σ^4 | 3 | 928 / 492 | 10 248 | 3 537 803 | 39 | 147 | 8.0 | 2.3 |
| stream + compact + rung transport, σ^4 | 4 | 2 599 / 965 | 40 382 | 18 012 911 | 187 | 291 | 56.4 | 31.7 |

- Peak RSS σ^4: 3 285 MB → 291 MB (×11); held bytes at crossing 4: 3 015 MB → 187 MB (×16:
  ×2.7 from streaming, ×6 from the compact store). The residual RSS is mostly the 81 MB caches.
- Same Poincaré polynomials at N = 3..6 (σ^3, σ^4) and the same final object multisets as the
  baseline.
