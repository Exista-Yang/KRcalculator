# SPEC_3C_GL6 rev1: plan to certified 3-colored gl_6 homology (wedge^3 on every component, N = 6)

Date 2026-11-11. This revises SPEC_3C_GL6.md after two critiques, gl6/C_sound.md (soundness) and
gl6/C_feas.md (feasibility and ordering). Every finding was checked against the code on dev f2d8be9 (M4b
merged as d49ee31) and on origin/colored-m4b 31f9c13 (the M4b hardening, still in progress). §7 records each
finding as adopted or rejected, with the reason. Nothing was compiled or run for this note.
Tags: [M] measured, [V] verified in code or data, [E] estimate, [D] decided here.

## 0. State on 2026-11-11

- **Run H** (cloud, colored-m4b 8ad4ac6; PLAN.md, last section) [M]. Exact unpruned xengine on (3,3):
  - trefoil at gl_3..gl_6: ranks 1 / 10 / 46 / 136, bigraded equal to Wang at every N, both primes,
    about 70 s and 0.6 GB;
  - reduced trefoil at gl_6 consistent with χ; pruned = unpruned for N ≤ 5;
  - Hopf link: 1 / 16 / 100 / 400, equal to the closed formula;
  - T(2,5), T(2,7), T(2,9): **out of memory at crossing 5 (13.6 GB)**.
- **Cause of the memory failure** [V, code]. It is not the Hom computations:
  - decompositions are cached per (p, sign, k) and `dloc` per (p, sign, j) (xengine.rs:568-604, 833-861);
  - all 4 (3,3) types exist after crossing 1, so every `hom_space` call happens in crossings 1-2.

  What fills crossing 5 is the whole unsimplified C⊗X, which is built before any elimination
  (`tensor_crossing`, xengine.rs:863-913). It is stored as `PMat<Poly>` at 32 B per term (u128 `Mono` +
  `Fp(u64)`, mpoly.rs:8, :101), plus fill-in from `simplify`. T(2,9) (peak about 5 000 objects against
  2 154) extrapolates to at least 50 GB [E].
- **The M4b hardening session is running on colored-m4b.** 31f9c13 makes the final kernel of `kernel_lazy`
  sparse (review solver F1/F2). That session owns homprov.rs, linsys.rs, g1.rs, sbim.rs, smap.rs, tmpl.rs
  and tests/hom_lazy.rs.
- **Run G** (N = 5, pruned, old solver) is redundant, since H already shows pruned = unpruned for N ≤ 5.
  Stop it (D2).
- **Wang references** (tests/data/wang_ref.txt; unreduced, positive σ_1^(2m+1), krhom conventions):
  - N = 3: rank 1;
  - N = 4, m = 1, 2: 10, 16;
  - N = 5, m = 1, 2: 46, 106;
  - N = 6, m = 1..4: 136, 396, 848, 1540.

  The generator script `r3scripts/wang.py` exists only in the laptop scratchpad. The main session commits
  it to tools/wang.py before launch (§3).

## 1. Decisions

**D1. 2-strand route** [revised].
- **Primary route: exact xengine after a memory fix** (milestone M2S′: one session, time-boxed, gate G2X).
  The certificate is exact end to end (C_sound §0):
  - exact decompositions with Π∘Ι = id and Σ ranks = rank X;
  - exact d∘d = 0;
  - exact closure.
- **Parallel route: sengine plus box points at N = 6 on (3,3)** (C_feas C3). Measured or estimated costs:
  ≤ 21 MB and 0.4-1.1 s per point at T(2,5) (SPEC_M4_M6a §6.6); 309 points unreduced and 20 reduced
  (R3 §4.3). This route is gate G2B of M6a-1b + MBOX in wave 2, and it runs whatever G2X does:
  - if G2X passes, G2B is the 2-strand regression of the localized path;
  - if G2X fails, G2B delivers T(2,7) and T(2,9).
- **Dropped:**
  - M2S part B (dagexact, `HomSource::Provider`);
  - the wdecomp, dloc and Rickard disk caches;
  - R-linear transport tables.

  Run H shows that Hom computations are no longer a cost. Rung-by-rung transport is built only if the
  measurements of M2S′ step 0 call for it.

**D2. Run G stops** [revised]. H gives its values and the pruned = unpruned check for N ≤ 5. F2b stays
cancelled.

**D3. Pruning lives only in M5/M6a** [unchanged].
- `HomProvider` stays N-agnostic.
- `SIndecOpts.prune_n` is fixed per builder.
- `killed(q)` is `Schur::killed_at` and nothing else.
- `DecData::restrict(n')` serves every n' ≤ prune_n.

**D4. Certification with dropped summands** [revised]. R2 §3 applies, plus the following guards, all as
release asserts:
- **Kept multiset.** In `dec_word`, the composite kept multiset must equal the kept part of `term_mult`.
  This assertion is mandatory and never debug-only. It is what protects kills in intermediate blocks
  (C_sound m9).
- **Π∘Ι = id** among the kept summands, at 2 points on all cosets; rank bookkeeping; the Rickard cell rule.
- **DecHom adapter guards** (M4b review F1/F2/F3; C_sound m6 = C_feas m7). Per call:
  - `expand(ladder_char(word(lq))) == [(q, 1)]`, using the same ladder function as `idem_needed`;
  - lq is the canonical ladder of q;
  - `mult == [v^s] expand(char X)[q]`;
  - the ends and degree of every map: ι: lq → X of degree s, and π: X → lq of degree −s, read from
    `ar.info`.
- **Degree scaling** per target (D14): the ι and π values on the top-stratum coset scale with the declared
  degree.
- Σιπ = e_X is checked only on unpruned (3,3).

**D5. Closure** [revised].
- **2 strands:** xengine uses `close_exact_colored`. The sengine route uses the box.
- **3 strands:** Grassmannian box points at fixed N, with the closed pairing
  P = H^y_id · d_id · G^x_id. This is the `pab` product of cclose.rs:1336-1362 (C_sound: correct).
- **Wrapper type for box data** (C_feas M5). Box data travel in
  `BoxClosureData { base: ClosureDataC, prune_n, n_box, valid_n, pairing, d2_cert: BTreeMap<(N', mode), D2Cert>, stamp }`
  with valid_n = min(prune_n, n_box). `ClosureDataC` itself stays unchanged, so its struct literals stay
  valid (xengine.rs:1085, 1301; cclose.rs:498, 2358; tests/cclose.rs:410, 544, 641, 758).
- **`close_box` refuses** (C_sound M5):
  - N' > valid_n;
  - pruned data without a recorded D² certificate for that (N', mode). The rank formula
    dim − r_out − r_in (cclose.rs:1395-1399) needs D² = 0.
- **Dropped for this whole target, 2 strands included** (C_feas m10): the dense e-monomial interpolation of
  SPEC_M4_M6a §6.5, `DeltaTower`, `Auto`, G4 and the laptop policy.

**D6. Gram weight 40 > 36 is settled** [unchanged]. C_sound confirms that `validate` uses
2·weight = ins + outs (cclose.rs:452-456). T0 at N = 4 covers the two (3,3,3) types with w = 10 > 9. Off
the path: an assertion over the 790 kept (3,3,3,3) types in a later basis run.

**D7. Pruned d∘d and the oracles** [revised].
- **Freivalds test.** D_{h+1} Π^{-1} D_h u = 0 for k = 3 random u per (h, q) block, per N' and per mode;
  the false-pass probability is ≤ p^-3 per block. The verdict is stored as `D2Cert` in `BoxClosureData`.
  Without it, ranks cannot be read (D5). The deterministic P2 extraction is used on small cases.
- **New gate GP3** (C_sound M4): pruned = unpruned on 3 strands, before any prune-6 production
  (M6a-2a).
- **Labels for the oracles** (C_sound m8):
  - Certifying: exact checks, two primes × two node sequences, the Wang and Hopf formulas, pruned =
    unpruned (GP3), nesting (prune 6 closed at N' = 4, 5 equals prune 4, 5), engine2 at labels 1, and
    base-point independence (m3).
  - Blocking only (they fail closed but cannot certify): duality at N = 4 ("expected"; it sees weights
    ≤ 9 only), the Lewark-type relation (1+t)Q ≥ 0, amphichirality, and the Gukov-Stošić Λ^3
    superpolynomial (optional).

  Every reported polynomial lists the oracles it passed, by label.

**D8. The (3,3,3) library** [revised]. It is fixed, knot-independent and prebuilt. Build order:
1. level 4, complete (42 pairs);
2. the **extended Markov corner** (C_feas m5):
   - the σ_1 decompositions and Rickard components on p⊗id_3 (the (3,3) pairs embedded in (4,2,3),
     (5,1,3) and (6,0,3)) at every level;
   - the σ_2 pairs from sources p⊗id_3, both signs (for [1,1,1,−2]);
   - closure generators for every type reached after σ_2;
3. a **large-X sample**: the about 20 level-6 pairs with the largest X (C_feas M2);
4. level 5 (158 pairs);
5. level 6 (338 pairs).

**Exact must be opened explicitly for (3,3,3).** As the code stands it is closed for (3,3,3) [V]:
- `Exact::default` has max_label 2 (homprov.rs:997);
- `Chain::for_labels` sets `max_label = 0` for labels ≥ 3 (homprov.rs:2821);
- costs_cloud.txt has no 3,3,3 records.

`provider_for` therefore builds its own `Chain` (public tuple struct): [Search S0-S2, `Exact { max_label:
3, .. }` with a cost table, Search S3-S6]. The cost table is scaled by the slowdown through
`CostTable::new()` + `record()`, both public, and gets its records from `libbuild calib`.

**Costs.** The (2,2,2) figure of 1.3 s per target [M] does not transfer: there X ≤ 432 flows, while
(3,3,3) reaches 31 500. M5b-run is planned at 1-4 sessions until the large-X sample has measured it.
F32b runs full certification at 1x. A tier replay (`Search::tiers(t, t)` per pair, public) is added only
if it is measured to save time.

**D9. Idempotents: minimal** [unchanged].
- `idem_needed(λ, N)` classifies the 2 open intermediate-block ladders in M5a.
- If it is non-empty for a kept type: M5i (SPEC_M4_M6a §5.2, e = P(a)) runs before level 6.
- 1 of the 55 (3,3,3) minimal-flow ladders is decomposable (review F1).

**D10. Order of targets** [revised].
1. **2 strands (M2S′, then M2S-run and G2B):** first T(2,5) and the T(2,4) link; then T(2,7), T(2,9), the
   T(2,6) and T(2,8) links, mirrors and reduced versions. The T(2,2m) links have no independent formula,
   so they are the first new values.
2. **3 strands, in this order:**
   - GP3, then T0;
   - the plain Markov words: unreduced at N = 6 only for [1,2], [1,1,2], [1,1,1,2] and [1,1,1,1,1,2];
     reduced for all of them;
   - rotations and [1,1,1,−2];
   - 4_1 reduced at N = 6 (the first new 3-strand value);
   - 5_2 and 6_2 reduced.
- Unreduced N = 6 for real 3-strand knots is decided after S3: M6c if w_c ≤ 10, otherwise a
  multi-session box run.
- 4 strands are out of scope.

**D11. Branch base and frozen files** [revised].
- All new branches start from dev f2d8be9, which contains M4b.
- New branches must not edit these files, which belong to the M4b hardening: homprov.rs, g1.rs,
  sbim.rs, linsys.rs, smap.rs, tmpl.rs.
- M5a's Arena export and import go in a new file src/smapio.rs. It uses only the public API:
  `Arena::{node, add, expr_of, info, certify_*}`, the pub enum `MapNode`, `TmplLib::save` and
  `read_cache_file`. (`tmpl::write_cache_file` is private: tmpl.rs:203.)
- Any branch that runs large-X Search (M5b and later) must first merge the hardening (≥ 31f9c13, sparse
  final kernel), from colored-m4b or from dev after that merge.

**D12. Laptop and persistence** [revised] (C_feas M6).
- The laptop does no compute.
- Layer-A archives, template caches and box checkpoints go to a data branch `data-gl6`, as compressed
  chunks of ≤ 45 MB with a manifest (hash and stamps per chunk).
- Every report states the cost of regenerating its data.
- Nothing larger than 50 MB goes into a code branch.

**D13. Certification wording** (C_sound m1). Results are certified "over F_p for p = 2^31−1 and 4294967291;
equal to the result over Q unless both primes are torsion primes". The headers in tests/data use the same
wording, plus the code commit.

**D14. Integrity of the localized path** [new] (C_sound M1, M2, M3, m2, m4).
- **Scaling check.** Since z_μ(λa) = λ·z_μ(a), every identity-coset value satisfies P(λz) = λ^w P(z), and
  likewise for the Gram and generator entries.
  - Cost: 2 extra scans per word and prime, at z and λz, with λ of large multiplicative order.
  - Assert the scaling law for every final edge, Gram and generator entry against its declared weight.
  - Report the count of entries that are 0 at both points.
  - Negative test: declare one entry at w + 1; the check must fail.
- **Checkpoints.**
  - Written atomically (write to a temp file, then rename).
  - Header: prime, node-seed hash, library stamp, word-list hash, prune_n, mode, N, code stamp.
  - A done-point bitmap.
  - Merging refuses mismatched headers and overlapping bitmaps.
  - `finish()` asserts that the complete down-set {μ ⊆ ρ} was received for every output ρ and every
    residual ρ.
- **Determinant pass.** A generator-only pass checks the generator determinant at all points before the
  first scan. Redrawing the nodes after scanning has started invalidates the whole node sequence.
- **Frozen structure.** The sparsity pattern and the pivot sequence are frozen from point 0. Values at
  box points never decide structure.
- **Caches.**
  - Every on-disk cache key includes: field, prune_n, the canonical-ladder rung words of p and of every
    summand, an algorithm stamp, and the library hash.
  - Positional data store their summand list. On load they are re-checked (Π∘Ι = id, Π invertible).
- **T0.** T0 runs uncapped once per (type, N', mode), on the exact library file used in production (hash
  recorded). After that, the MOY rank cap is only an assertion.

## 2. Milestones

Dependencies:
- H (done) → M2S′ → M2S-run.
- M5a → M5b → M5b-run.
- M6a-1a has no dependencies. {M5a, M6a-1a} → M6a-1b. {M6a-1b, MBOX} → gate G2B.
- {M6a-1b, MBOX, M5b} → M6a-2a → M6a-2b. {M2S-run or G2B} → M6a-2b (reduced 2-strand references).
- {M6a-2b, M5b-run} → M6a-3 → M6a-4.
- Conditional: M5i if `idem_needed` is non-empty (before level 6); M4c if G1(c) ends Exhausted at 2x budget.

| id | kind | needs | gate | sessions |
|---|---|---|---|---|
| H | compute | — | done: trefoil 136, Hopf 400 | done |
| M2S′ | dev + gate | dev | G2X: σ^5 ≤ 6 GB, T(2,5) = Wang, σ^9 projected ≤ 12 GB | 1 (time-boxed) |
| M2S-run | compute | G2X | T(2,7) 848, T(2,9) 1540, links, mirrors, reduced → gl6_t2n.txt | 1 |
| M5a | dev | dev | sindec suite; `pruned_plan` and `idem_needed` printed | 1-2 |
| MBOX | dev | dev | box suite with negative tests; staged rank stage at 4_1 size | 1 |
| M6a-1a | dev | — (against the SPEC_M4_M6a §5.1 interface) | synthetic engine suite | 1 |
| M6a-1b | dev + gate | M5a, M6a-1a, MBOX | G2B: (3,3) σ^1..σ^9 at N = 6 by box = Wang | 1 |
| M5b | dev + compute | M5a, hardening merge | G1(c): level 4, extended corner, large-X sample | 1-2 |
| M5b-run | compute | M5b | G1(c): levels 5 and 6, both primes | 1-4 |
| M6a-2a | dev + gate | M6a-1b, MBOX, M5b | GP3; T0 at N = 4 (5, 6 when the library exists) | 1 |
| M6a-2b | compute | M6a-2a, M2S-run or G2B | GM: Markov words = Wang at N = 4..6; 185 knots at prune 4 | 1-2 |
| M6a-3 | compute (+ small code) | M6a-2b, M5b-run | S3 at prune 6; G41: 4_1 at N = 4, 5, 6 | 2-3 |
| M6a-4 | dev + compute | M6a-3 | 5_2 and 6_2 reduced at N = 6 | 4+ |
| M5i / M4c | dev | conditional | `idem_needed` / Exhausted set | 1 / 1+ |

Critical path to 4_1: M5a → M6a-1b → M6a-2a → M6a-2b → M6a-3, which is 6-9 sessions. It runs in parallel
with M5a → M5b → M5b-run (3-8 sessions). Together this is about **7-11 sessions, i.e. waves 4-5** (§4).

### M2S′. Memory-first exact 2-strand engine (branch `colored-2s`)
- **Step 0: measure** (≤ 45 min).
  - Runs: σ^3 and σ^4 to completion; σ^5 under an RSS watchdog at 12 GB.
  - Record per stage: objects, components and Poly terms before and after `simplify`; bytes split into
    stored components, caches and transients; RSS; seconds split into decomposition, transport,
    elimination, d∘d and closure.
- **Step 1: stream insertion by C-degree in `tensor_crossing`.**
  - Insert o⊗X slice by slice in h(o).
  - Eliminate a pivot a → b only once in(b) and out(a) are complete (the eligibility rule of
    SPEC_M4_M6a §6.4).
  - Keep the full build behind `XOpts::stream = false` for tests.
- **Step 2: compact storage** of the stored components, if step 0 shows that Poly terms dominate. The
  format: an xengine-local interned u32 monomial index plus a u32 coefficient (both primes are < 2^32),
  converted to `PMat` on use. pmat.rs and mpoly.rs are not edited, because they are shared.
- **Step 3: `XOpts::close_stages`** (C_sound m5).
  - Closure data after every crossing m ≥ 2, read by `XScan::close_stage(m, ns, reduced)`.
  - Uses `crossings[..m]` for `global_q`.
  - Keeps `close`'s pruned guards: N' ≤ prune_n, and the stage-m `final_d2` closes to 0.
  - `check_d2` stays on in production.
- **Step 4, conditional:** rung-by-rung transport (R1 §3.3). Only if step 0 shows that transient `kron_id`
  products are > 30 % of peak RSS or that transport is > 50 % of the time.
- **Housekeeping:** move `hopf_formula` (private, tests/xengine.rs:156) into tests/common/mod.rs; add
  `StageInfo` fields for RSS and bytes; add the xscan flags `--stream --close-stages --json OUT`.
- **Tests (default ≤ 2 min; (3,3) cases `#[ignore]`, run in the cloud):**
  - streaming = full build: same Poincaré polynomials at N = 3..6 and the same final object multisets,
    on (2,1), (2,2), (3,1) and (3,2) words, and on (3,3) σ^1 and σ^2 pruned at N = 3, 4;
  - the compact store round-trips to the same `PMat`;
  - `close_stage(m)` equals a separate scan of σ^m for (2,2), m = 2, 3;
  - `close_stage` refuses N' > prune_n.
- **Gate G2X**, (3,3) σ^5 with close_stages at N = 3..8, peak RSS ≤ 6 GB:
  - T(2,5) = Wang 3 N 2 as full polynomials (16 / 106 / 396 at N = 4, 5, 6), and rank 1 at N = 3;
  - stage closures: Hopf = `hopf_formula`, T(2,3) = run H, and the T(2,4) link;
  - both primes identical; exact d∘d = 0; χ = `decat_chi`;
  - reduced results recorded, with the base-point check (marked block Some(1) = Some(0)) for knots;
  - the per-stage growth projects σ^9 at ≤ 12 GB.

  If G2X passes, M2S-run follows. If it fails, stop and report the measurements; T(2,7..9) then come from
  G2B.
- Budget: 5 h, KRHOM_THREADS=4.

### M2S-run. 2-strand production (compute, branch `colored-2s`)
- **Runs:** σ^9 and σ^−9 with close_stages (σ^10 optional, for T(2,10)), N = 3..8, reduced and
  unreduced, F31 and F32b.
- **Acceptance:**
  - T(2,7) = 848 and T(2,9) = 1540 at N = 6, and every Wang block at N = 4, 5, 7, 8 that exists;
  - mirrors give P(1/t, 1/q);
  - links: χ, the two primes, and Hopf = the formula;
  - the base-point check for knots.
- **Output:** tests/data/gl6_t2n.txt, with the D13 header. It is the reduced 2-strand reference for M6a.
- Budget ≤ 5 h, RSS watchdog at 13 GB. Wall time [E] tens of minutes per prime and sign (the trefoil
  takes 70 s; growth is roughly quadratic in the peak object count).

### MBOX. Box closure layer at fixed N (branch `colored-box`)
- **New src/boxclose.rs:**
  - `BoxNodes`, `box_points`, `point_coords`, and per-block factorial Schur S^(i) with its triangular
    inverse;
  - `schur_to_epoly`, with restriction to N' and to reduced;
  - `BoxAcc` with the D14 checkpoints: atomic writes, header, bitmap, merge refusal, `finish()` down-set
    check, residuals at weight w + 1;
  - `check_scaling(v_z, v_λz, λ, w)`;
  - `BoxClosureData` (the D5 wrapper) and `close_box`, which carries the valid_n and D² guards;
  - `closed_d2_freivalds`, which records a `D2Cert`;
  - `box_moy_check` (T0), uncapped, with the library hash recorded.
- **src/cclose.rs** (C_feas M7). Additive changes only, leaving the `ClosureDataC` struct untouched:
  - a staged per-type rank stage: today `close_exact_colored` builds every `type_state_space` up front
    (cclose.rs:1299);
  - the MOY rank cap;
  - `exact_mem_estimate` refitted;
  - the rank stage takes P in place of G_y·F.
- **Tests** (tests/boxclose.rs; default < 60 s):
  - (a) Synthetic: sparse random homogeneous R^λ polynomials for (3,3), (2,2,2) and (3,3,3), at weights up
    to top + 1, N = 4..6, both modes, both primes. Box extraction = `nf`. They are sparse because dense
    weight-27 polynomials are too large (88 233 e-monomials at w = 22).
  - (b) End to end, on exact xengine data, against `close_exact_colored` for N = 3..6, both modes:
    (2,2) σ^3 and (3,2) σ^2 by default; (3,3) σ^1..σ^3 marked `#[ignore]` (cloud, about 70 s).
  - (c) Negative tests, each of which must fail:
    - Freivalds on one perturbed P entry;
    - Freivalds on a pruned scan with one `final_d2` edge re-inserted;
    - a point added twice;
    - an overlapping shard merge;
    - a mismatched header;
    - an incomplete down-set;
    - an entry declared at w + 1 (scaling check);
    - a prune-4 scan closed at N' = 5 (Err);
    - pruned data without a `D2Cert` (Err).
  - (d) Shuffled point order and a split into shards give the same result.
  - (e) The staged rank stage on synthetic data the size of 4_1 (34-40 types, the measured `ins`) stays
    within its memory estimate.
- Budget: 1 session.

### M5a. sindec core with pruning (branch `colored-m5`)
- **New src/sindec.rs:**
  - the SPEC_M4_M6a §5.1 interface plus `SIndecOpts.prune_n` and `DecErr`;
  - `DecHom`, with the D4 guards and the per-target degree-scaling check (D14);
  - `provider_for` (D8 chain with `Exact { max_label: 3 }` and a scaled cost table);
  - `ty`, `dec_rung`, `dec_word` (release kept-multiset assert), `rickard_comp`, `closure_gens`,
    `freeze`, `restrict`;
  - `pruned_plan`, `idem_needed`;
  - layer-A cache keyed by λ, field, table strings, canonical-ladder rung words, the ladder-function stamp,
    `TMPL_ALG` and the provider stamp, and re-checked on load.
- **New src/smapio.rs:** Arena export and import on top of the public API (D11).
- **A frozen DecData returns Err on a missing pair, never skips it** (C_sound m7).
- **Tests** (tests/sindec.rs; default ≤ 2 min on (2,1), (2,2), (3,2), (1,1,1); the (3,3) comparisons
  against xengine are marked `#[ignore]` and run in the cloud): D1, D3-D5, D7-D10 of the old SPEC, plus:
  - negative DecHom tests: a decomposable L_q, a wrong lq, a wrong mult, a wrong degree;
  - a stale stamp on load is rejected;
  - a missing pair gives Err.
- **Gate:** the suite passes; `pruned_plan(3,3,3; 4, 5, 6)` and `idem_needed` over the visited (3,3,3)
  blocks are printed.

### M6a-1a. sengine core without a library (branch `colored-m6a`)
- Built against the binding §5.1 interface (C_feas M8):
  - `SCx` / `LocMap` over Arc'd `Layout`s;
  - a Markowitz heap with streaming eligibility (SPEC_M4_M6a §6.4);
  - sparsity pattern and pivot trace frozen from point 0 (D14);
  - determinism trace, memory accounting and memory refusal;
  - the `PointScan` and `gram_at_point` types.
- **Tests:** synthetic complexes with known homology; streaming = full elimination; identical traces at
  2 points; the memory refusal.
- Launched only when a session slot is free. If the §5.1 interface drifts in M5a, M6a-1b absorbs the
  change.

### M6a-1b. Library wiring and the 2-strand box gate (branch `colored-m6a`)
- **Wiring:**
  - crossings via `dec_word` and `rickard_comp`;
  - d_C⊗id rung by rung, blockwise per coset;
  - K_0, LB and degree asserts; `check_d2` on unpruned scans;
  - scans at λz for the scaling check.
- **Tests:** sengine against xengine at 2 points, on (2,2) and (3,2) by default and on (3,3) with
  `#[ignore]`; (1,1,1) object counts against engine2.
- **Gate G2B**, using a one-word box loop of about 50 lines with `BoxAcc`:
  - (3,3) σ^1..σ^9 at N = 6, unreduced = Wang 136 / 396 / 848 / 1540, and Hopf = 400;
  - reduced = gl6_t2n.txt where it exists;
  - both primes, two node sequences, and the scaling check passes.
  - Cost [E] ≤ 1.5 h.
- **Report:** the 4_1 type count from the scan, as both "types in the LB complex" and "closable types
  used"; this settles the 34 vs 38 discrepancy (C_feas m9).

### M5b. (3,3,3) library and gate G1(c) (branch `colored-m5`, after M5a and the hardening merge)
- **Build:**
  - `SIndec::prebuild`;
  - src/bin/libbuild.rs with `plan|run|replay|calib --lam 3,3,3 --n 4|5|6 --corner markov-ext --sample
    largex --field 31|32b --cache DIR --jobs 4 --wall 600 --rss 2048 --json OUT`;
  - one process per pair, resumable from layer A;
  - `calib` writes (3,3,3) cost records (C_feas M4).
- **Contents, in this order:** level 4, then the extended corner (D8), then the large-X sample, plus the
  closure generators of the kept (λ,λ) types at level 4. F31 first, then F32b with full certification.
- **G1(c):** 100 % of the targets of each set certified, for both primes.
  - A Time or Memory outcome is rerun at 2x budget.
  - An Exhausted outcome at 2x budget goes to the M4c list.
  - If `idem_needed` is non-empty, M5i comes first.
- **Persistence:** layer A goes to `data-gl6` in chunks (D12).
- **Report:** per set and prime: targets, Σ mult, tier histogram, seconds (mean, p90, max), max X, peak
  RSS, failures, and a **measured** projection for levels 5 and 6.

### M5b-run. Levels 5 and 6 (compute)
- `libbuild run --n 5`, then `--n 6`, F31, then F32b, resuming from `data-gl6`. G1(c) must reach 100 %.
- 1-4 sessions, according to M5b's projection.

### M6a-2a. Box driver, GP3 and T0 (branch `colored-m6a`)
- **`box_run`** (batch of words per point) and **`box_finish`** (returns `BoxClosureData`, runs the staged
  closure from MBOX and records Freivalds `D2Cert`s), with all of D14; CLI `sscan box|finish|point`.
- **Gate GP3** (C_sound M4):
  - (a) (1,1,1) pruned at N = 2 (B_{w0} killed) and N = 3, on all 185 three-strand knots, against
    engine2's uncolored gl_2 and gl_3 results (validated on the tables);
  - (b) (2,2,2) pruned at N = 3, 4 against unpruned (2,2,2) on 4_1, 5_2 and 6_2, both modes. G1(a)
    certified all (2,2,2) targets.
  - The report lists, per oracle, the code paths, types and library pairs it touched.
- **T0**, uncapped: N = 4 on 6 types; N = 5 and 6 on 20 and 40 types once those levels exist.

### M6a-2b. First 3-strand gl_6 checks (compute)
- **Gate GM.** The plain Markov words at prune 6, using the extended corner only:
  - unreduced at N = 4, 5 for all words;
  - unreduced at N = 6 for [1,2] (= [6 choose 3]_q), [1,1,2] (Hopf 400), [1,1,1,2] (T(2,3), 136) and
    [1,1,1,1,1,2] (T(2,5), 396);
  - reduced at N = 4..6 for all words, = gl6_t2n.txt;
  - [1,1,1,1,2] (the T(2,4) link) = gl6_t2n.txt;
  - [1,1,1,−2] if the corner has the negative sign;
  - P2 = 0, Freivalds passes, and the scaling check passes.
- **The 185 knots at prune 4,** N = 3, 4, both modes:
  - N = 3 gives rank 1;
  - N = 4: duality with engine2's uncolored gl_4 (labelled blocking only);
  - χ; two primes; two node sequences.
- 1-2 sessions; the Markov words at N = 6 unreduced take [E] 1.2-2 h each.

### M6a-3. First new 3-strand value (compute + small code)
- **First hour: S3 at prune 6.** Single points of 4_1 at prune 4, 5, 6, and of [1,1,1,2]. Record c_e, κ_S,
  σ, τ, t_lib, M_ev and actual/LB, then re-issue R4 §5.
- **Then:**
  - rotations of the Markov words (these need level-5/6 σ_1 pairs from M5b-run);
  - 4_1 reduced at N = 6 (2.4 h [E], band from S3);
  - 4_1 at N = 5, both modes;
  - 4_1 at N = 4 by nesting (= the M6a-2b result);
  - base-point check (block 1 marked) at N = 4 and 5, and at N = 6 if S3 puts it under 3 h;
  - amphichirality, mirror word, rotation [−2,1,−2,1];
  - two primes × two node sequences, and the scaling check.
- **Output:** tests/data/gl6_333.txt, with the D13 header and the passed oracles by label.
- 2-3 sessions [E] 8-14 h, sharded through the checkpoints.

### M6a-4. Scale within 3 strands
- u32 `LocMap` storage, if S3 shows a memory bound.
- 5_2 (7.9 h [E]) and 6_2 (9.7 h [E]), reduced at N = 6, sharded over 2 sessions each.
- Class B if S3 allows it.
- The M6c decision from w_c.

## 3. What to launch now

**Launch now (2 sessions), both from dev f2d8be9:**
1. **M2S′** on `colored-2s`. Edits src/xengine.rs, src/bin/xscan.rs, tests/xengine.rs, and the new
   tests/common/mod.rs.
2. **M5a** on `colored-m5`. Edits the new src/sindec.rs, the new src/smapio.rs, and tests/sindec.rs.

Neither touches the M4b hardening files (homprov.rs, linsys.rs, g1.rs, sbim.rs, smap.rs, tmpl.rs). Both
consume homprov and sbim read-only.

**Main session, before launch (no compute beyond seconds):**
- stop run G;
- commit `r3scripts/wang.py` as tools/wang.py;
- append the Wang blocks 3 7 m and 3 8 m (m = 1..4), 3 5 3, 3 5 4, 3 4 3 and 3 4 4 to wang_ref.txt;
- create the branch `data-gl6`.

**Next, as slots free:**
- MBOX on `colored-box`: new src/boxclose.rs, additive changes in src/cclose.rs, tests/boxclose.rs, and
  appended tests in tests/cclose.rs;
- M6a-1a on `colored-m6a`: new src/sengine.rs and tests/sengine.rs.

**Merges:** only src/lib.rs (`pub mod` lines) and Cargo.toml (`[[bin]]` entries) are shared. Merge order
into dev: hardening, then M2S′ / M5a / MBOX / M6a-1a in any order.

## 4. Waves and deliverables

| wave | sessions | certified result at the end of the wave |
|---|---|---|
| 0 | H (done) | trefoil T(2,3) = Wang 136 and Hopf 400 at gl_6 (and N = 3..5), exact, both primes |
| 1 | M2S′, M5a (now); MBOX, M6a-1a | if G2X passes: T(2,5) = 396 and the T(2,4) link at N = 3..8, exact |
| 2 | M2S-run, M5b, M6a-1b (G2B) | T(2,7) 848, T(2,9) 1540, the T(2,6) and T(2,8) links, mirrors, reduced; G2B = Wang through the localized path |
| 3 | M6a-2a, M5b-run | GP3 (pruned = unpruned on 3 strands); T0 |
| 4 | M6a-2b, M5b-run (continued) | the Markov words = Wang through the 3-strand engine at N = 4..6; the 185 knots at N = 3, 4 |
| 5 | M6a-3 | S3; **4_1 reduced at gl_6, the first new 3-strand value** |
| 6+ | M6a-4 | 5_2, 6_2 reduced at gl_6 |

- **Already certified:** the 3-colored gl_6 trefoil (136) and Hopf link (400). Exact, with Wang and the
  closed formula as independent checks.
- **First new values:** the T(2,4) link at gl_6, in wave 1 (one session), if G2X passes. Otherwise the
  T(2,2m) links come from G2B in wave 2.
- **First new 3-strand value:** 4_1 reduced at N = 6, in wave 4-5. It slips by 1-2 sessions if M5i or M4c
  is needed.

## 5. Risks and fallbacks

| risk | signal | fallback |
|---|---|---|
| xengine memory not reducible enough | G2X step 0 or the σ^9 projection > 12 GB | G2B (box route) for T(2,7..9); xengine stays the reference for σ^≤5 |
| (3,3,3) Search incomplete or slow at large X (X ≤ 31 500; only X ≤ 432 measured) | M5b large-X sample; G1(c) Exhausted at 2x | Exact at label 3 with calib records; M4c light leaves on the failing set |
| a kept type needs an idempotent | `idem_needed` non-empty | M5i, 1 session |
| silent wrong answer in the localized path that passes two-prime agreement (degree bookkeeping, accumulators, caches, pruning code on general (3,3,3) types) | D14 scaling, checkpoint and down-set checks; GP3; GM; kept-multiset assert; D2Cert | stop production; bisect on (1,1,1) / (2,2,2) / Markov words |
| pruning implementation on 3 strands (the rule itself is backed by the `kernel_of_gl_n_quotient` unit test and CL1; C_sound §0) | GP3, nesting, GM or Freivalds fail | diagnose on small words with prune_n = None; the timeline slips |
| per-point constants for (3,3,3) at prune 6 (4_1 band 11-194 s and 0.1-3.7 GB per point) | S3 | reduced first; u32 storage (M6a-4); unreduced N = 6 only via M6c |
| 5-h sessions and sharded data | missing or mismatched checkpoints or archives | `data-gl6` chunks; headers and bitmaps (D14); resumable drivers |
| merge conflicts with the hardening | — | frozen files (D11); `BoxClosureData` wrapper; new files only |
| Gram weight above the box top (D6) | T0 fails | revisit the weight convention (R3 §1.1) |

## 6. Session briefs

**Brief M2S′** (branch colored-2s from dev f2d8be9; 5 h; KRHOM_THREADS=4)
- Read PLAN.md (last section), xengine.rs:1-80 and 553-1000, and this SPEC §0 and §2 M2S′.
- Step 0: σ^3, σ^4 and σ^5 (RSS watchdog at 12 GB). Record per-stage objects, components and Poly terms
  (before and after simplify), bytes by category, RSS and seconds by phase. Push the table first.
- Steps 1-3: streaming by C-degree with the eligibility rule; compact u32 component storage local to
  xengine.rs (if terms dominate); close_stages with `crossings[..m]` and the pruned guards; `check_d2`
  stays on.
- Step 4 (conditional): rungwise transport.
- Do not edit homprov, g1, sbim, linsys, smap, tmpl, cclose, pmat or mpoly. No disk caches, no part B.
- Tests: streaming = full build; the compact store round trip; close_stage = a separate scan;
  close_stage refuses N' > prune_n. Default ≤ 2 min; (3,3) tests `#[ignore]`.
- Gate G2X: σ^5 at N = 3..8 within 6 GB; T(2,5) = Wang 3 N 2 (396 at N = 6); Hopf; T(2,3); the T(2,4)
  link; both primes; exact d∘d = 0; χ; the base-point check; σ^9 projected ≤ 12 GB.
- Report: the step-0 table, polynomials as JSON, peak RSS per stage, and go/no-go for M2S-run.

**Brief M2S-run** (compute, colored-2s)
- σ^±9 (σ^10 optional) with `--stream --close-stages`, N = 3..8, both modes, F31 and F32b.
- Accept only if: T(2,7) = 848 and T(2,9) = 1540 at N = 6, plus every Wang block present; the mirrors
  hold; the primes agree; χ holds; the base-point check passes.
- Write tests/data/gl6_t2n.txt with the D13 header. Budget ≤ 5 h, RSS ≤ 13 GB.

**Brief M5a** (branch colored-m5 from dev f2d8be9; 1-2 sessions)
- Read gl6/R2_m5prune.md, SPEC_M4_M6a §5, review_m4b/certify.md (F1-F3), and this SPEC D3, D4, D8,
  D11, D14.
- New src/sindec.rs:
  - the §5.1 API with prune_n and `DecErr`;
  - `DecHom` with the F1/F2/F3 guards and the degree-scaling check;
  - `provider_for` (own Chain, `Exact { max_label: 3 }`, scaled cost table);
  - all-killed and identity shortcuts; M5's own form re-certification at 2 points; Neumann among the kept
    summands;
  - the kept-multiset release assert;
  - `rickard_comp` with the cell-rule assert; `closure_gens`; `freeze`; `restrict`;
  - layer A with full stamps;
  - `pruned_plan`, `idem_needed`;
  - a missing pair gives Err.
- New src/smapio.rs: Arena export and import via the public API.
- Do not edit homprov, g1, sbim, linsys, smap, tmpl, xengine or cclose.
- Tests: as in §2 M5a, plus the negative DecHom tests. Default ≤ 2 min; (3,3) against xengine
  `#[ignore]`.
- Report: kept-target counts and Σ mult per level; the `idem_needed` list; the (3,3) layer-A size;
  provider seconds; the (3,3) DecData timing (G2B needs it).

**Brief MBOX** (branch colored-box from dev f2d8be9; 1 session)
- Read gl6/R3_closure.md §2-6, R4_m6a.md §2, this SPEC D5, D7, D14, and C_sound M2/M5.
- New src/boxclose.rs: nodes, points, factorial Schur with its triangular inverse; `BoxAcc` with D14
  checkpoints; `schur_to_epoly`; `check_scaling`; `BoxClosureData` and `close_box` with the guards;
  `closed_d2_freivalds` → `D2Cert`; `box_moy_check` (uncapped, library hash recorded).
- cclose.rs, additive only: the staged per-type rank stage with the MOY cap, a refitted
  `exact_mem_estimate`, and P in the rank stage. Do not change the `ClosureDataC` struct.
- Do not edit xengine, homprov, smap, tmpl or linsys.
- Tests: §2 MBOX (a)-(e), including every negative test. Default < 60 s; (3,3) σ^1..σ^3 `#[ignore]`.
- Report: extraction cost per entry, against R3 §3; the memory of the staged rank stage at 4_1 size.

**Brief M6a-1a** (branch colored-m6a from dev f2d8be9; only if a slot is free)
- Read SPEC_M4_M6a §5.1 and §6.3-6.4, and gl6/R4_m6a.md §0-1 and §7-9.
- New src/sengine.rs: `SCx`/`LocMap`/`Layout`; a Markowitz heap with streaming; frozen sparsity and
  pivot trace from point 0; the determinism trace; memory accounting and refusal; the `PointScan` and
  `gram_at_point` types; a trait matching §5.1 for the library.
- Tests on synthetic complexes. Do not edit any existing file except lib.rs.

**Brief M6a-1b** (colored-m6a after the M5a, M6a-1a and MBOX merges)
- Wire `dec_word` and `rickard_comp`; rung-by-rung coset-blockwise transport; K_0, LB and degree
  asserts; `check_d2` on unpruned scans; λz scans.
- Tests against xengine at 2 points ((3,3) `#[ignore]`); (1,1,1) against engine2.
- Gate G2B: (3,3) σ^1..σ^9 at N = 6 via the one-word box loop = Wang 136 / 396 / 848 / 1540 and Hopf 400;
  reduced = gl6_t2n where it exists; two primes × two node sequences; the scaling check.
- Report: the per-point time and RSS for σ^1..σ^9, and both 4_1 type counts.

**Brief M5b** (colored-m5 after M5a, merged with the hardening ≥ 31f9c13)
- `SIndec::prebuild`; libbuild plan/run/replay/calib.
- Build, F31 then F32b: level 4 → extended Markov corner → large-X sample, plus the level-4 closure
  generators.
- G1(c) at 100 % per set: 2x reruns for Time and Memory outcomes; Exhausted → the M4c list. If
  `idem_needed` is non-empty, M5i first.
- Push layer A to `data-gl6` in ≤ 45 MB chunks.
- Report: per set and prime, targets, tiers, seconds (mean, p90, max), max X, peak RSS, failures; the
  measured level-5/6 projection.

**Brief M5b-run** (compute)
- `libbuild run --n 5`, then `--n 6`, F31, then F32b, resuming from `data-gl6`. G1(c) at 100 %.
- Report as M5b. 1-4 sessions.

**Brief M6a-2a** (colored-m6a; needs M6a-1b, MBOX and M5b)
- `box_run` and `box_finish` with all of D14 (checkpoint headers, bitmaps, down-set check, generator-only
  determinant pass, frozen structure); `sscan box|finish|point`.
- Gate GP3: (1,1,1) prune 2, 3 on the 185 knots (tools/dev/ntrunc_t1/data/knots3_all.txt) against engine2
  gl_2 / gl_3; (2,2,2) prune 3, 4 against unpruned on 4_1, 5_2 and 6_2.
- T0 uncapped at N = 4.
- Report: per oracle, the types and pairs it touched; the T0 table.

**Brief M6a-2b** (compute, colored-m6a; needs gl6_t2n.txt or the G2B reduced results)
- Gate GM: the Markov words as in §2 M6a-2b.
- Then the 185 knots at prune 4, N = 3, 4, both modes (duality labelled blocking only).
- Report: the polynomials, oracle verdicts by label, and per-point cost.

**Brief M6a-3** (compute, colored-m6a; needs M5b-run)
- First hour: S3 at prune 6; re-issue R4 §5.
- Then: rotations and [1,1,1,−2]; 4_1 (reduced at N = 6; both modes at N = 5; N = 4 by nesting;
  base-point check; amphichirality, mirror, rotation; two primes × two node sequences; scaling check).
- Shard over 2-3 sessions. Write tests/data/gl6_333.txt (D13 header, oracle labels).

**Brief M6a-4** (after M6a-3)
- u32 `LocMap` if S3 says so (test: u32 = u64 on 4_1 at N = 4).
- 5_2 and 6_2 reduced at N = 6, 2 sessions each.
- Propose M6c or a multi-session unreduced plan with the measured constants.

## 7. Critique log

Verdicts: **A** = adopted, **P** = partly adopted, **R** = rejected. "[V]" means checked in code or data
for this revision.

### C_sound (soundness)

| finding | verdict | reason, and where it lands |
|---|---|---|
| M1. Declared degree ≠ actual degree gives silent box zeros | A | [V] `certify_pairs` checks no ends or degrees (review certify.md F3; tests.md #1). sengine has no analogue of `check_degrees`. → D14 scaling check; D4 DecHom; gates in MBOX (c), M6a-1b and GM. |
| M2. Box accumulators can be corrupted across checkpoints, shards and node redraws | A | No residual exists at w ≥ top, and systematic corruption can agree across two runs. → D14 checkpoints; MBOX tests (c) and (d). |
| M3. Cache stamps missing; dloc is positional | A (largely moot) | [V] dloc stores summand positions (xengine.rs:833-860), placed via `place[&(o,ks)][si]` (:908). The M2S disk caches are dropped (C_feas m1). The stamping rule stays for layer A and state-space caches (D14). |
| M4. No 3-strand pruned = unpruned test | A | [V] B_{w0} is killed at N = 2; G1(a) certified all (2,2,2) targets. → gate GP3 in M6a-2a. |
| M5. valid_n and D² not enforced by the API | A | [V] xengine has these guards (xengine.rs:1140-1144, 1276-1290); the rank formula (cclose.rs:1395-1399) needs D² = 0. → D5 `BoxClosureData` / `close_box`; MBOX negative tests. |
| m1. Two primes certify F_p, not Q | A | → D13. |
| m2. The T0 rank cap is circular | A | → D14 (T0 uncapped once, hash recorded). |
| m3. No reference for reduced mode | A | → base-point checks in G2X, M2S-run and M6a-3 (at N = 6 for 4_1 only if S3 allows: it doubles the reduced points). |
| m4. Zero-dropping changes the pivot order | A | → D14 frozen structure; M6a-1a. |
| m5. close_stage guards; keep check_d2 | A | → M2S′ step 3. |
| m6. DecHom must carry F1/F2/F3; F32b must re-certify | A | → D4, D8 (F32b at full certification). |
| m7. Rotations need pairs outside the corner; Err on a missing pair | A | → D8 extended corner; rotations in M6a-3 after M5b-run; M5a Err. |
| m8. Conjectural oracles cannot certify | A | → D7 labels. |
| m9. Risk weight belongs on the implementation | A | → §5 risk rows; the kept-multiset assert is a release assert. |

### C_feas (feasibility and ordering)

| finding | verdict | reason, and where it lands |
|---|---|---|
| C1. The SPEC is stale; the failure is memory | A | [V] PLAN.md f2d8be9. → §0; wave 0 done; M2S replaced by M2S′. |
| C2. Root cause is the stored exact complex | A | [V] xengine.rs:568-604, 833-913; mpoly.rs:8, :101. → M2S′ (measure, stream, compact). The compact store stays local to xengine.rs, because pmat.rs and mpoly.rs are shared. |
| C3. sengine + box is the shortest route; invert D1 for n ≥ 5 | P | G2B adopted as a wave-2 gate and as the fallback. **The inversion is rejected:** the box route needs M5a + M6a-1a/b + MBOX (3-4 sessions), while M2S′ is one session and keeps exact certification. G2B becomes primary only if G2X fails. |
| M1. Session sizes are optimistic | A | → §2 session column. |
| M2. Calibrations come after the plans that depend on them | A | → large-X sample in M5b; S3 at prune 6 in M6a-3's first hour. |
| M3. The library cost rests on X ≤ 432 | A | → M5b-run 1-4 sessions until measured. The 1 GB dense-kernel part is already fixed on colored-m4b 31f9c13 (sparse final kernel), so M5b must contain the hardening (D11). |
| M4. Exact is closed for (3,3,3) | P | [V] homprov.rs:997, 1019-1024, 2821; no 3,3,3 cost records. → D8. **Sub-claim rejected:** "a slowdown needs a homprov edit". `CostTable::new()` and `record()` are public (homprov.rs:843, 851), so `provider_for` can build a scaled copy. |
| M5. ClosureDataC struct-literal conflict | A | [V] all 8 literals. → D5 wrapper. |
| M6. Persistence across sessions is not designed | A | → D12 `data-gl6`. |
| M7. The staged rank stage has no owner | A | [V] cclose.rs:1299 builds all state spaces up front. → MBOX. |
| M8. Missing edges; M6a-1 serialized for no reason | A | → M2S-run/G2B → M6a-2b edge; M6a-1a split out; plain Markov words in M6a-2b. |
| m1. Drop part B and the caches; stop run G | A | → D1, D2. |
| m2. Branch base; the private writer | A | [V] tmpl.rs:203. → D11, smapio.rs. |
| m3. Test budgets | A | → `#[ignore]` cloud tests; sparse synthetic polynomials. |
| m4. Paths: knot list, wang.py, hopf_formula | A | [V] all three. → §3 pre-launch actions; M2S′ housekeeping. |
| m5. The Markov corner is under-counted | A | → D8 extended corner. |
| m6. M6a-3 is 8-14 h | A | → 2-3 sessions; unreduced N = 6 limited to 4 Markov words. |
| m7. DecHom must assert that L_q is indecomposable | A | = C_sound m6 → D4. |
| m8. No replay hook | P | [V] there is none. But a replay can be built from the public `Search::tiers(lo, hi)`. The saving is unmeasured, so F32b is budgeted at 1x (D8). |
| m9. 4_1 type count, 38 vs 34 | A | T0 covers all kept types (6 / 20 / 40), so coverage does not depend on the count. M6a-1b reports both counts. |
| m10. Dense interpolation: kept for 2 strands or not | A | Dropped everywhere (D5). |
| Revised wave table | P | Adopted with GP3 and T0 added before the Markov words, so 4_1 falls in wave 4-5 rather than wave 4. |
