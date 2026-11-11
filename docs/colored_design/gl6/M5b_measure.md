# M5b measurements: the (3,3,3) library before the long runs

Date 2026-11-11. Machine: 4 cores (Intel Xeon @ 2.30 GHz; host stamp `4cpu-Intel_R__Xeon_R__Processor_@_2.30GHz`,
not the 2.10 GHz class of costs_cloud.txt), 15 GiB. Branch colored-m5 (bf326b7 + merge of origin/dev ce91f0e = bb52c19,
libbuild ffa27be). Field F31, `libbuild run --jobs 4 --wall 600 --rss 2048` (each child KRHOM_THREADS=1), cache
/tmp/m5b/meas (throw-away). Provider: `provider_for_stamped((3,3,3), costs, 1.0)` = [Search S0-S2, Exact(label 3),
Search S3-S6].

## 1. First pass without a (3,3,3) cost table

| set | jobs | certified | failures | wall | max X | peak RSS |
|---|---|---|---|---|---|---|
| level-n4 (42 pairs + 5 closure) | 47 | **47** | 0 | 13.8 s | 576 | 148 MB |
| sample-largex20-n6 | 20 | 6 | 13 Refused, 1 Exhausted (pair 336, also at 2x) | 100 s | 21 600 | 2 608 MB (pair 135 at 2x) |
| corner-markov-ext-n6 (16 cross + 16 closure) | 32 | 26 | 6 closure Exhausted (also at 2x) | 24 s | 900 | 208 MB |

Level 4: 19 pair targets (Σ mult 19) + 5 closure types, every target through Search S1 (no Exact call); slowest
pair 308 (X = 576 flows) 9.2 s.

Large-X sample (the 20 kept-source level-6 pairs with the largest X: 10 at X = 21 600, 2 at 14 400, 4 at 12 960, 2 at
10 800, 4 at 8 640 [sic: ties broken by index]):
- certified: pairs 90, 282, 36, 39, 136 and 135 (135 needed 2x: 2 140 MB > 2 048 MB at 1x, 2 608 MB at 2x, 38 s);
  provider seconds per target 0.1-2 s; pair 282 (X = 12 960, 16 targets) 19.6 s, 648 MB.
- **Refused (13)**: Search's local decomposition tree finds no *admissible* decomposable window: every candidate
  window is above the window caps (fibre product 6 400, unknowns 40 000) and has no cost record. Examples from
  `KRHOM_TRACE` (pair 392): `[3,6,0] [E0(1) F1(3) F0(1) E1(2)]` vs summand `[F1(1) F0(1) E0(1)]` degree -4 (108 x 720,
  338 139 unknowns); `[3,3,3] [E1(3) E0(1) F1(3) E1(2)]` (30 x 180, 40 322 unknowns: 0.8 % above the cap). S6 also
  refused all 16 of its windows. Exact on the whole ladders is hopeless here (fibre products 1.5e7-1e8, unknowns
  4e8-9e9).
- **Exhausted**: pair 336 (X = 14 400): Search S0-S2 and S3-S6 exhausted with no admission refusal, at 1x and 2x.

Corner (D8, built at prune 6, which contains levels 4 and 5 through `restrict`): the sources are the 4 kept types
p' ⊗ id_3 of (3,3,3) (the (3,3) types with the third strand split off). All 16 crossing jobs certify: σ_1 terms with
all their Rickard components (max X 900 flows, 17 s per job), σ_2 terms both signs (every σ_2 term of p' ⊗ id_3
decomposes by the identity shortcut: one summand, multiplicity 1, X = its canonical ladder). Of the 16 closure types
reached after σ_2, 10 certify and 6 end Exhausted: Search finds no closure generator in the top degree (e.g. 0,0,3/3,0,0/0,3,0
degree 18), Exact is not admitted (fibre product 400 but 852 011 unknowns, no record).

## 2. Calibration (window and closure cost records)

Both failure modes are admission, which the D8 mechanism addresses with measured (3,3,3) cost records:
`libbuild calib --windows FILE` measures each listed exact call (labels, a, b, d) in a child under the watchdog
(limit 60 s, a kill is a censored record = never admitted) and writes a v2 table for this host.
The list: the 403 distinct windows refused in the large-X trace (both ι and π directions) and the 414 whole-ladder
closure calls Hom^d(R, L_p), Hom^d(L_p, R) of the corner and level-6 closure types.

Exact's per-call ceiling at label 3 is raised from 10 s to Search's window ceiling, 120 s (`EXACT_L3_MAX_CALL_S`,
part of the provider stamp): the first measured (3,3,3) closure calls take 16-20 s. Admission still requires a
measured, uncensored record; no cap of the label rule was changed and no check was touched.

(Calibration results and the second pass: section 3.)
