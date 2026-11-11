# M5a report: sindec core with Foam_N pruning (branch colored-m5)

Date 2026-11-11. Machine: 4 cores (Xeon 2.10 GHz), 15 GiB, no swap. Base: origin/dev eb6d8d5, merged with
origin/colored-m4b 402b06a (hardening) at the start and 22ff843 (HANDOFF update, docs only) at the end.
No frozen file (homprov, g1, sbim, linsys, smap, tmpl, xengine, cclose) was edited; lib.rs gained
`pub mod smapio; pub mod sindec;`.

## What was built

**src/sindec.rs**
- SPEC_M4_M6a §5.1 API (`SType`, `SSummand`, `Decomp`, `SIndecOpts`, `SIndec`, `DecData`, `DecStats`)
  plus `SIndecOpts.prune_n: Option<u8>` and `DecErr` (Killed, NeedsIdempotent, Refused, Provider
  {engineering}, Check, NotPrecomputed, Invalid, Cache). All methods return `Result<_, DecErr>`
  (`DecErr: Into<String>`).
- `killed(q)` = `Schur::killed_at(q, N)` only (D3). Killed types get no ladder, no ι/π, no rungs.
- `DecHom` trait with the blanket impl over every `HomProvider` (the only code touching it), with
  the D4 guards before any provider work (`dechom_precheck`): lq is the canonical ladder of q,
  L_q indecomposable of type q unless e_q (then `check_idempotent`), mult = [v^s] mult of b_q in
  char X, tops; and after the call (`dechom_postcheck`): mult maps each way, ends and degrees from
  `ar.info`, and the **D14 degree-scaling check** at the top-stratum coset rep(q).
- **Scaling law (measured here):** every localized entry of a map of normalized degree d satisfies
  v(λz) = λ^u v(z) with 2u = d − sh(tgt) + sh(src), uniformly over flows and cosets (checked on all
  (2,2) Exact pairs, then asserted on every target). `check_scaling` rejects a declared degree off by
  1 or 2 (negative test). (3,3) unpruned: 3 880 entries checked, 128 zero at both points.
- `provider_for` / `provider_for_stamped(lam, costs, slowdown)`: max label ≤ 2 → `Chain::for_labels`;
  max label 3 → own Chain [Search S0-S2, Exact { max_label: 3 } with the cost table scaled by
  `slowdown` via `CostTable::new()+record()`, Search S3-S6]. The stamp (COST_ALG, table hash,
  slowdown, chain) enters the layer-A header.
- `ty`, `dec_rung`, `dec_word` (rung by rung), all-killed shortcut (release assert killed_rank =
  flows(X), no provider call) and identity shortcut; per target M5's own re-certification (form at
  rep(q) at 2 points equal and invertible, genuineness); π' = G^{-1}π as `Lin`; Neumann
  orthogonalization among kept summands (unitriangularity asserted); Π∘Ι = id at 2 points on every
  coset of X's support; Σιπ = e_X when nothing is killed; rank bookkeeping; the **kept-multiset
  release assert** (`assert_eq!`) for composites.
- `rickard_comp_at(p, i, sign, j)` (3-strand position i; `rickard_comp(p, c, d, sign, j)` keeps the
  §5.1 form for position 0) with the cell-rule assert per crossing and Err on a killed term;
  component degrees asserted = 1 + s_σ − s_τ.
- `closure_gens` (provider closure, re-certified by `certify_closure`), `freeze`, `DecData::restrict(n)`
  (Err for n > prune_n), `prebuild(&Plan)`.
- Layer A: one file per certified target (pairs, closure) in `cache_dir/layer_a`, named by the hash
  of the key (λ, field, table strings of p and q, canonical-ladder rung words of L_p, L_q, rung, s,
  mult); the header holds the stamps (`krhom-sindec v2`, `LADDER_STAMP`, `TMPL_ALG`, provider stamp).
  A mismatching header is rejected and recomputed; loaded maps go through the same checks as fresh
  ones (form, genuineness, ends, degrees, scaling).
- `pruned_plan(schur, N)`, `visited_blocks`, `decomposable_ladders`, `idem_needed`.
- A frozen DecData returns `Err(NotPrecomputed)` on a missing pair, Rickard list, closure or type.

**src/smapio.rs**: Arena sub-DAG export/import on the public API. Templates by content: hom_space
templates by call + index, Rickard by (c, d, sign, j), both re-obtained from the library and checked
against a content hash; Search's added window templates stored in full and re-certified on import.
Prime fields only.

**Deviations.** `rickard_comp` needs a strand position on 3 strands (`rickard_comp_at`). `ty` refuses
decomposable kept types (`NeedsIdempotent`); idempotents are M5i. `reduced()` is the identity stub
(R2 §7.5).

## Tests (tests/sindec.rs)

Default suite (13 tests, 4.1 s): D1 types of (2,1), (1,2), (2,2), (3,2), (3,1), (1,1,1) and the
(3,3,3) kill counts (219 types; 12 / 58 / 138 kept at N = 4/5/6; 40, 30, 15, 4 per block at N = 6);
D3/D10 all G1 pairs of (2,1), (2,2), (3,2), (1,1,1) unpruned and pruned (targets and Σ mult = plan);
D4 all crossing terms; D5 a mini localized scan (decompositions, Kron transport, Rickard components,
Gaussian elimination at a point) **equal to xengine's final objects** on (2,1), (1,1), (2,2), (3,2),
pruned and unpruned, plus the braid relation and R2 on (1,1,1); D7 layer-A round trip (no provider
call on reload, identical blocks at 2 points), stale TMPL_ALG stamp rejected, restrict(3/4) of a
prune-5 build = fresh builds, restrict(6) Err; D8 freeze and missing pair Err; D9 Send+Sync and
2 threads; negative DecHom: decomposable L_q ((3,1,3) 2,0,1/0,1,0/1,0,2), wrong lq, wrong mult,
wrong degree (provider wrapper adding a dot), wrong declared degree in the scaling check.

Ignored, run in this session (197.8 s total): `d3_33_all_pairs_and_decdata_timing`,
`d5_against_xengine_33` ((3,3) σ, σ² at N = 3, 4 and unpruned, σ^-2 at N = 4: all equal to xengine;
σ² unpruned 15 objects, xengine 45.6 s vs sindec mini scan 79.8 s), `gate_plan_333`.

Full suite `KRHOM_THREADS=4 cargo test --release` on the final merged code (c8c676b + report): **196 passed, 0 failed, 20 ignored** (sindec 13 / 3 ignored in 3.9 s; homprov 17 / 3 in 56 s; acceptance 25 in 92 s).

## Gate output

pruned_plan(3,3,3) (494 pairs over 7 blocks):

| N | pairs with kept source | non-identity pairs with targets | all-killed | kept targets | Σ mult | killed targets |
|---|---|---|---|---|---|---|
| 4 | 42 | 11 | 24 | 19 | 19 | 150 |
| 5 | 158 | 92 | 40 | 324 | 364 | 646 |
| 6 | 338 | 286 | 0 | 1 592 | 2 566 | 646 |
| none | 494 | 429 | 0 | 3 151 (= R2) | 5 146 | 0 |

(killed targets counted over kept-source pairs; equal counts at N = 5 and 6 are a coincidence of
that convention.) (3,3): unpruned 72 targets, Σ mult 121 (= G1(a)); N = 5: 18 / 21; N = 4: 2 / 2;
N = 3: 0.

**idem_needed** over the visited (3,3,3) blocks: empty at N = 4, 5, 6. The decomposable canonical
ladders are exactly three, all with kill level 7:
- (3,3,3) 0,1,2/1,1,1/2,1,0 (10 500 flows), char = b + b_{1,0,2/0,2,1/2,1,0};
- (3,4,2) 0,2,1/1,1,1/2,1,0 (5 250 flows), char = b + b_{0,2,1/1,2,0/2,0,1};
- (3,5,1) 0,3,0/1,1,1/2,1,0 (1 050 flows), char = b + b_{0,3,0/1,2,0/2,0,1}.

This settles R2 §9's open point: **M5i is not needed for N ≤ 6 on 3 strands.**

## (3,3) measurements (F31, cold, KRHOM_THREADS=1 for the provider)

- All 12 non-identity G1 pairs, unpruned: 72 provider calls, 11.7 s; 12.7 s total.
- DecData of every crossing term (both signs) of every (3,3) type, with all Rickard components and
  closure generators, unpruned: **56.4 s** (provider 10.4 s in 80 calls, M5 checks 0.26 s; the rest is
  the Rickard nonzero tests and composite Π∘Ι checks); 10 types, 34 decompositions, 24 Rickard lists,
  4 closure types; max X 90 flows. Layer A: 80 files, **89 KB**. Arena 180 MB (flow indices).

## Open problems / what M5b needs

- `SIndec::prebuild` exists (visits kept-source pairs of a plan); M5b adds src/bin/libbuild.rs (one
  process per pair, resumable from layer A, `calib`), the (3,3,3) build order of D8 and G1(c).
- No (3,3,3) provider call was made here: Search/Exact cost on (3,3,3) X up to 31 500 flows is
  unmeasured. costs_cloud.txt has no 3,3,3 records; `provider_for` admits Exact at label 3 by the
  label caps (fibre ≤ 1 600, unknowns ≤ 40 000) only.
- Rickard nonzero tests and composite checks dominate the (3,3) DecData time (≈ 45 s of 56 s): for
  (3,3,3) consider testing nonzero at the top stratum first or caching evaluators.
- Arena memory (180 MB for (3,3)) is mostly flow indices; for (3,3,3) drop indices of intermediate
  composites or measure first.
- Layer A is per prime (templates are); F32b needs its own files (keys contain the field).
