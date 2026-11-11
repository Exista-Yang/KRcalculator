# M2S′ gate G2X: (3,3) σ^5 with close_stages at N = 3..8

Date 2026-11-11, cloud (4 cores, 15 GB), branch colored-2s. Command
`xscan 3,3 "1 1 1 1 1" 3,4,5,6,7,8 --close-stages --both [--prime 32b] --json OUT`
(streaming, compact store, rung-by-rung transport, check_d2 on). Results over F_p for
p = 2^31−1 and 4294967291; equal to the result over Q unless both primes are torsion primes.
JSON (all polynomials, per-stage measurements): `m2s/sigma5_33_F31.json`, `m2s/sigma5_33_F32b.json`,
`m2s/sigma5_33_F31_parallel.json` (the threaded engine of 995e8c0's parent 487aa40+, same results).

## Verdict: PASS

| item | result |
|---|---|
| peak RSS ≤ 6 GB | **PASS**: 555 MB (serial engine, F31), 525 MB (F32b); 1 066 MB with 4 threads |
| T(2,5) = Wang 3 N 2, full polynomials | **PASS** at N = 4, 5, 6, 7, 8 (16 / 106 / 396 / 1011 / 2076); N = 3 rank 1 |
| stage 2: Hopf = `hopf_formula` | **PASS** at N = 3..8 (1 / 16 / 100 / 400 / 1225 / 3136) |
| stage 3: T(2,3) = Wang 3 N 1 (= run H) | **PASS** at N = 3..8 (1 / 10 / 46 / 136 / 307 / 586) |
| stage 4: T(2,4) link | recorded (below); χ = decat_chi, both primes |
| both primes identical | **PASS**: all 48 closures (stages 2-5 × 2 modes × N = 3..8) identical |
| exact d∘d = 0 | **PASS**: asserted after every crossing (unpruned scan, check_d2 on) |
| χ = `schur::decat_chi` | **PASS** at every stage 2-5 and N = 3..8 (test `g2x_sigma5_33`) |
| reduced, base-point check | **PASS**: marked block 1 = marked block 0 at stages 3 and 5 (both primes); χ(unreduced) = [N choose 3] χ(reduced) |
| σ^9 projected ≤ 12 GB | **PASS** (projection below; the σ^9 run itself is reported in the session result) |

New Wang blocks generated for the gate with tools/wang.py and appended to tests/data/wang_ref.txt:
3 8 1 (rank 586), 3 8 2 (rank 2076), 3 7 4 (rank 4483).

## Ranks (F31 = F32b)

| stage (link) | mode | N=3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|
| 2 (Hopf) | unreduced | 1 | 16 | 100 | 400 | 1225 | 3136 |
| 2 (Hopf) | reduced | 1 | 4 | 10 | 20 | 35 | 56 |
| 3 (T(2,3)) | unreduced | 1 | 10 | 46 | 136 | 307 | 586 |
| 3 (T(2,3)) | reduced | 1 | 3 | 9 | 27 | 27 | 27 |
| 4 (T(2,4)) | unreduced | 1 | 22 | 184 | 864 | 2817 | 7326 |
| 4 (T(2,4)) | reduced | 1 | 6 | 26 | 88 | 139 | 202 |
| 5 (T(2,5)) | unreduced | 1 | 16 | 106 | 396 | 1011 | 2076 |
| 5 (T(2,5)) | reduced | 1 | 5 | 25 | 125 | 125 | 125 |

The reduced knot ranks stabilize at N ≥ 6 = 2·3 to 3^3 and 5^3 (the cube of the reduced uncolored
rank), as expected from the large-N behaviour; the T(2,4) link values are new (no independent formula).

T(2,4), gl_4, unreduced: t^0q^-12 + t^0q^-10 + t^0q^-8 + t^0q^-6 + t^2q^-14 + t^2q^-12 + t^2q^-10 +
t^3q^-22 + t^3q^-20 + t^3q^-18 + t^4q^-24 + 2t^4q^-22 + 3t^4q^-20 + 3t^4q^-18 + 2t^4q^-16 + t^4q^-14;
reduced: t^0q^-9 + t^2q^-13 + t^3q^-19 + t^4q^-21 + t^4q^-19 + t^4q^-17. (All N in the JSON.)

## Per-stage cost of σ^5 (F31, serial engine of 7e2afb9)

| crossing | objects inserted → minimal | held at once (objects / MB / terms) | RSS peak MB | s |
|---|---|---|---|---|
| 1 | 4 → 4 | 4 / 0.0 / 300 | 19 | 0.4 |
| 2 | 189 → 15 | 143 / 4.4 / 0.41 M | 132 | 51.6 (Hom spaces 49.8) |
| 3 | 928 → 40 | 492 / 39 / 3.5 M | 147 | 10.9 |
| 4 | 2 599 → 85 | 965 / 187 / 18.0 M | 296 | 94.1 |
| 5 | 5 580 → 156 | 1 634 / 394 / 38.5 M | 524 | 323.1 |

Closure data 0.4-0.75 s per stage; the 48 closures at N = 3..8 take 28 s; total 510 s.
With 4 threads (487aa40/995e8c0): crossing 5 in 192 s (transport 192 → 60 s, elimination 124 s),
RSS 1.07 GB (thread arenas and chunked plain results).

## σ^9 projection

Held bytes grew ×2.11 from crossing 4 to 5, i.e. as n^3.3 (n^5.4 from 3 to 4; the exponent falls).
With n^3.3: 720 / 1 197 / 1 860 / 2 743 MB held at crossings 6-9; with n^4: 4.1 GB at crossing 9;
even n^5 gives 7.5 GB. RSS = held × 1.33 (serial) or × 2.7 (4 threads): σ^9 at ≈ 3.6 GB serial or
≈ 7.4 GB threaded with n^3.3, ≤ 11 GB threaded with n^4. All ≤ 12 GB. Time is the binding constraint:
crossing time grew as n^5.5 (94 → 323 s serial, 51 → 192 s threaded), projecting ≈ 2.5 h threaded
for σ^9 per prime and sign.

## σ^7 and σ^9 runs (same session, after G2X)

Commands as above with words σ^7 / σ^9 (`--close-stages --both`, N = 3..8), MALLOC_ARENA_MAX=4,
KRHOM_THREADS=4; two runs shared the 4 cores most of the time.
- **σ^9, F31** (engine 995e8c0): complete. JSON `m2s/sigma9_33_F31.json` (96 closures: stages 2-9 ×
  2 modes × N = 3..8). Scan 5 686 s, closures 236 s, total 5 925 s, **peak RSS 2.78 GB**.
- **σ^7, F32b** (engine 5e27894 + chi mode): complete, peak RSS 1.66 GB; `m2s/sigma7_33_F32b.json`.
  All 72 stage closures (stages 2-7) identical to the F31 σ^9 run.
- **σ^9, F32b**: see the session report (it was still running when this section was written).

Checks on the F31 σ^9 run (all PASS): χ = `schur::decat_chi` (xscan `--chi-only`) for all 48
unreduced closures; stage 2 = Hopf formula; T(2,3), T(2,5), T(2,7), T(2,9) = Wang 3 N m as full
polynomials wherever a block exists (T(2,7) at N = 4..8, T(2,9) at N = 4..8); base-point check
(marked block 1 = 0) at stages 3, 5, 7, 9; d∘d = 0 exactly after every crossing.
New Wang blocks: 3 8 3 (5006), 3 8 4 (9856), 3 4 3, 3 5 3, 3 4 4, 3 5 4.

| stage | knot/link | unreduced N=3..8 | reduced N=3..8 |
|---|---|---|---|
| 6 | T(2,6) | 1 / 28 / 292 / 1568 / 5321 / 13916 | 1 / 8 / 50 / 252 / 363 / 492 |
| 7 | T(2,7) | 1 / 22 / 190 / **848** / 2339 / 5006 | 1 / 7 / 49 / 343 / 343 / 343 |
| 8 | T(2,8) | 1 / 34 / 424 / 2560 / 8929 / 23386 | 1 / 10 / 82 / 560 / 755 / 974 |
| 9 | T(2,9) | 1 / 28 / 298 / **1540** / 4483 / 9856 | 1 / 9 / 81 / 729 / 729 / 729 |

Per crossing of σ^9 (F31, 4 threads, CPU shared with a second run):

| crossing | inserted → minimal | held objects | held MB | RSS peak MB | s | s transport | s elim |
|---|---|---|---|---|---|---|---|
| 5 | 5 580 → 156 | 1 634 | 394 | 798 | 201 | 61 | 132 |
| 6 | 10 249 → 259 | 2 304 | 657 | 1 104 | 386 | 89 | 280 |
| 7 | 16 984 → 400 | 3 310 | 800 | 1 457 | 937 | 204 | 700 |
| 8 | 26 163 → 585 | 4 446 | 1 213 | 2 129 | 1 421 | 309 | 1 058 |
| 9 | 38 164 → 820 | 5 774 | 1 630 | 2 776 | 2 460 | 477 | 1 889 |

The σ^9 projection of the gate (≤ 12 GB; 2.7-4.1 GB held) held: 1.63 GB held, 2.78 GB RSS.
Elimination is now 77 % of the time (fill-in products, decode/encode of the compact store).
