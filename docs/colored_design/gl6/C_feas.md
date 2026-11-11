# C_feas: adversarial critique of SPEC_3C_GL6 rev 1 (feasibility and ordering)

Date 2026-11-10/11. Read-only: nothing compiled or run (one < 1 s Python summary of
tests/data/homprov/g1a_cloud.json). Sources checked: dev = origin/dev f2d8be9 (M4b merged as d49ee31;
PLAN.md last section), src/xengine.rs, src/pmat.rs, src/mpoly.rs, src/linsys.rs, src/homprov.rs,
src/tmpl.rs, src/cclose.rs, src/cube.rs, src/bin/g1.rs, tests/xengine.rs, tests/cclose.rs,
docs/colored_design/{SPEC_M4_M6a.md, HANDOFF_M4b.md, review_m4b/*.md}, the four gl6 notes.
Tags: [V] verified in code or data, [E] estimate.

## 0. Bottom line

The SPEC was written against a state that has already moved on. Run H is done. The trefoil result
(rank 136) is in hand, and the SPEC's assumed next bottleneck (transport time) is not the one that
was measured. The measured one is memory: T(2,5), T(2,7) and T(2,9) run out of memory at crossing 5
(13.6 GB). The SPEC's 2-strand milestones (M2S, M2S-run) carry RSS budgets of ≤ 4 GB and have no
memory work in scope. On 3 strands the order is mostly right, but there are four problems:
- the sessions are sized at the optimistic end;
- the decisive calibrations (the library cost at large X, S3 at prune 6) are scheduled after the
  plans that depend on them;
- persistence across 5-h sessions is not designed;
- two pieces of critical-path work have no owner: the staged closure rank stage, and the
  ClosureDataC change that crosses file ownership.

## 1. Critical

**C1. The SPEC is stale. Run H's measured failure is memory, and M2S does not address it.** [V]
- PLAN.md (dev f2d8be9, "First 3-colored gl_6 results"):
  - trefoil at gl_3..gl_6 in ~70 s / 0.6 GB: ranks 1, 10, 46, 136 = Wang, bigraded, both primes;
  - Hopf 1/16/100/400 = the formula;
  - pruned = unpruned for N ≤ 5;
  - T(2,5), T(2,7), T(2,9): **out of memory at crossing 5 (13.6 GB)** in the unpruned exact scan.
- Consequences for the SPEC:
  - §4 wave 0 has already been delivered.
  - M2S Gate B (sigma^5) is exactly the case that fails.
  - The M2S budget "RSS ≤ 4 GB" and the M2S-run plan (sigma^10, ≤ 4 GB) contradict the
    measurement. R1 §3.4 ("memory stays under 2 GB") is falsified.
  - M2S part A (rung-by-rung transport, disk caches, close_stages, timing fields) and part B
    (homprov hook) are all about time. Run H shows time is fine: the whole trefoil took 70 s,
    including `cube::local_cx`, which also calls `hom_space` (cube.rs:62).

**C2. Root cause: the stored exact complex. Neither `kernel_lazy` nor the decompositions.** [V, code]
- Why it is not the Hom computations:
  - decompositions are cached per (p, sign, k) (`decompose`, xengine.rs:568-604) and dloc per
    (p, sign, j) (xengine.rs:833-861);
  - (3,3) has only 4 closable types, and all 4 are objects after crossing 1 (W_k = L_{t_k}, identity
    shortcuts; m6_run33k.log: "4 objects -> 4, 0 decomposition Hom calls");
  - so every `hom_space` call, and every possible dense `kernel_dense` (review_m4b/solver.md F1),
    happens in crossings 1-2;
  - crossing 5 makes no Hom calls at all. PLAN's first suspect ("dense final kernel of
    kernel_lazy") is therefore excluded by the code.
- What crossing 5 holds:
  - the whole unsimplified C⊗X, built before elimination (tensor_crossing, xengine.rs:863-913;
    module doc: "The whole tensor complex is built before eliminating");
  - each component is a `PMat` of `Poly` with terms (u128 Mono, Fp(u64)), i.e. 32 B per term plus a
    24 B Vec header per entry (mpoly.rs:8, :101-103; pmat.rs:10-14; field.rs:48);
  - plus Gaussian fill-in in `simplify` (xengine.rs:933-965).
- Scale:
  - at the T(2,5) peak, the localized End data are ΣS = 48 849 scalars (SPEC_M4_M6a App. M:
    "2 154 / 44 667 / 48 849"), against > 13.6 GB exact;
  - T(2,9), with peak about 5 000 objects (R1 §1), extrapolates to about 4-5x that, i.e. ≥ 50 GB [E].
- Rung-by-rung transport removes only the transient 900×900 `kron_id` products, not the stored
  components. Term streaming (SPEC_M4_M6a §6.4) saves ≤ 20 %, because W_2 + W_3 hold 50 of the 63
  ranks.
- Required M2S re-scope (memory first, time-boxed to one session):
  - (a) record per-stage RSS and unsimplified component and term counts (StageInfo has only seconds
    planned);
  - (b) **streaming by C-degree**: insert o⊗X slice by slice in h(o), and eliminate a pivot a→b
    only once in(b) and out(a) are complete (the §6.4 eligibility rule in the C-direction).
    This keeps about 2 slices live, about n/2 fewer [E];
  - (c) compact storage: a per-degree u32 monomial index (MonoCache) plus a u32 coefficient,
    about 4x [E];
  - (d) gate: sigma^5 at N = 3..8 within 6 GB, with a sigma^9 projection ≤ 12 GB; otherwise stop
    and use C3;
  - drop part B and the Rickard/dloc caches (m1).

**C3. The shortest route to T(2,5..9) at gl_6 is now sengine plus box, not xengine.** [V data, E cost]
- SPEC_M4_M6a §6.6 for (3,3) T(2,5): ≤ 21 MB per point, 0.4-1.1 s per point.
- R3 §4.3: N = 6 needs 309 points unreduced and 20 reduced, with ≤ 48 unknowns per entry.
- T(2,9) is then minutes per prime. Every ingredient is already on the plan:
  - M5a (3,3) DecData: the 80 (3,3) targets certified in 9.1 s on the cloud, g1a_cloud.json [V];
  - M6a-1 `scan_at_point`;
  - MBOX `BoxAcc`.
- D1's ranking ("(iii) only a regression of M6a against M2S") should be inverted for n ≥ 5. Make
  "(3,3) sigma^1..sigma^9 at N = 6, unreduced = Wang 136/396/848/1540, both primes and node
  sequences" a joint M6a-1 + MBOX gate in wave 2. It needs only a ~50-line one-word box loop, not
  M6a-2's driver.
- Wang's closed formula is the independent certificate for the knots. The exact xengine stays the
  reference wherever it fits: sigma^≤4, and links or reduced results after C2's fix.

## 2. Major

**M1. Session sizes are at the optimistic end.**
- M6a-1 in one session means all of SPEC_M4_M6a §6.3-6.4:
  - PEngine, SCx/LocMap/Layout, an Eval per thread;
  - a Markowitz heap with streaming, and T5 with 200 seeds and a sensitivity check;
  - traces, memory accounting, PointScan/gram_at_point.
  For scale, xengine alone is 1 313 lines.
- M6a-2 in one session (≤ 5 h) means the box driver, checkpoints, finish, staged closure, CLI and
  T0, **plus** S3 (≤ 1 h), 185 knots (1.3 h wall) and the 2-strand regression.
- M5a: the §5.1 API, pruning, layers A/B, `Arena` export/import, D1-D10.
- Track record: M4b took a WIP session, a finishing session, a review, and hardening that is still
  running.
- Realistic sizes: M5a 1-2 sessions, M6a-1 2, M6a-2 2 (split dev from runs). The first new value
  (4_1) then needs about 7-9 sessions, not "about 5" (§4).

**M2. The calibrations come after the decisions that depend on them.**
- (a) S3 at prune 6 cannot run in M6a-2.
  - 4_1 uses 34 of the 40 level-6 types (R3 §4.1).
  - M6a-2 depends only on M5b level 4 plus the corner, and runs in parallel with M5b-run (§2 graph).
  - So τ, μ and c_e for types with r up to 7 200 are first measured inside M6a-3, after M6a-3 and
    M6a-4 have been planned.
- (b) The library cost has no large-X sample: level 4 has r ≤ 144, and the corner has X ≤ 900.
  - Fix: M5b adds a stratified sample of the ~20 largest-X level-6 pairs.
  - S3 at prune 6 moves to the first hour of M6a-3, or M5b-run is serialized before M6a-2.

**M3. D8's library estimate rests on (2,2,2) data with X ≤ 432 flows.** [V]
- The G1 jobs run sequentially (g1.rs:468-470, `supervised`), so 523 s / 400 targets = 1.3 s per
  target is correct. But (3,3,3) X reaches 31 500 flows.
- In G1(a), per-target time does not track X: the mean is ≤ 4 s at every X.
- RSS: one job with X = 216 peaked at 856 MB. Against `--rss 2048` × 4 jobs, Memory reruns at 2x
  are likely.
- The F1 dense final kernel can reach 1 GB at dim 3 300 × 40 000 unknowns (linsys.rs:256-279). Its
  fix is in hardening scope (linsys.rs is frozen by D11).
- So "0.2-3 h wall per prime, 1-2 sessions" is not grounded. Plan 1-4 sessions until M2(b) has
  measured it.

**M4. The Exact fallback is closed by default for λ = (3,3,3).** [V]
- `Exact::default` has `max_label: 2` (homprov.rs:997), and `basis` admits on `ar.lam()`
  (homprov.rs:1019-1024). costs_cloud.txt has no (3,3,3) records. So every whole-ladder Exact call
  on (3,3,3) is Refused.
- Only Search windows (6 400 fibres / 40 000 unknowns, `win_ok`) remain.
- The recalibration tool (`g1 calib`) lives in the frozen g1.rs.
- `CostTable` has no slowdown factor, and its `map` is private (homprov.rs:829-832), so
  `provider_for(.., slowdown)` cannot apply one without editing homprov.
- Fix: `provider_for` sets `max_label: 3`, and libbuild gets its own calib mode.

**M5. A merge conflict hides in MBOX.** [V]
- MBOX adds `valid_n` and `pairing` to `ClosureDataC` (cclose.rs:404-410). But the struct is built
  by literal in several places:
  - xengine.rs:1085 and :1301, which M2S owns and MBOX "must not touch";
  - cclose.rs:498 and :2358;
  - tests/cclose.rs:410, 544, 641 and 758, so the change is not append-only.
- Fix: keep the box fields in a wrapper (`BoxClosureData { base, valid_n, pairing }`) or add a
  constructor. Otherwise both branches edit the same lines.

**M6. Persistence across 5-h sessions is not designed.**
- Template cache: G1(a) on the laptop left 4.6 MB for 233 targets with X ≤ 432 (largest file
  414 KB). Level 6 has about 10x the targets, larger windows and two primes, so about 50-200 MB [E].
- Box accumulators for 5_2 / 6_2 reduced: about 0.6-1e6 P entries × 26 × 8 B ≈ 130-200 MB, scaled
  from R3 §7.
- Both exceed the "< 50 MB in git" rule (D12), yet M5b-run and M6a-4 shard their work across
  sessions.
- Fix: compressed archives in ≤ 45 MB chunks on a data branch, plus a stated cost to regenerate per
  session.

**M7. The 3-strand closure rank stage has no owner.**
- At N = 6 it needs staged per-type state spaces, the MOY rank cap and a refitted
  `exact_mem_estimate` (R3 §6; R4 §5.4: today's estimate reaches about 40 GB reduced and TB
  unreduced). `close_exact_colored` builds every `type_state_space` up front (cclose.rs:1295-1300).
- The MBOX brief omits this. M6a-2 lists it but does not own cclose.rs.
- Fix: assign it to MBOX and test it on data the size of 4_1.

**M8. Edges are missing from the graph, and M6a-1 is serialized for no reason.**
- M2S-run → M6a-2(a) and → M6a-3: the reduced references come from gl6_t2n.txt. With C1, reduced
  T(2,5) and the T(2,4) link have no reference yet.
- M6a-1 waits for the M5a merge. Its heap, streaming, LocMap, T5 and trace parts need no library,
  so it can start in wave 1 against the binding §5.1 interface. That saves one wave.
- The plain Markov words (σ_2 last, no rotations) need only the corner, so they can move into
  M6a-2 and give the first 3-strand gl_6 check one wave earlier. Rotations and 4_1 still wait for
  M5b-run.

## 3. Minor

- **m1.** Drop M2S part B (dagexact, `HomSource::Provider`) and the Rickard/dloc disk caches: run H
  shows Hom is no longer a cost. Stop run G (N = 5, old solver, the 5 401 s call class): H already
  gives pruned = unpruned for N ≤ 5. That frees the cloud box.
- **m2.** Branch base: branch from dev f2d8be9, not "origin/colored-m4b 8ad4ac6".
  - Hardening may touch smap.rs and tmpl.rs (review_m4b/tests.md §5; `refresh_certified` is
    O(#nodes) per certification, a scaling note for the M5 Arenas).
  - Put M5a's export/import in a new src/smapio.rs on top of the public API.
  - `tmpl::write_cache_file` is private (tmpl.rs:203), so the M2S caches would need a tmpl.rs edit.
- **m3.** Test budgets:
  - MBOX (b) runs exact (3,3) sigma^1..3 (the trefoil alone takes ~70 s) inside "< 60 s";
  - M5a D5 and the M6a-1 tests compare against unpruned (3,3) xengine inside "≤ 2 min";
  - (3,3,3) synthetic polynomials of weight 27 are dense (≥ 88 233 e-monomials at w = 22).
  - Fix: `#[ignore]` or cached fixtures, and sparse random polynomials.
- **m4.** Paths:
  - the knot list is tools/ntrunc_t1/data/knots3_all.txt (185 lines), not tests/data;
  - wang.py exists only in the laptop scratchpad (r3scripts/wang.py). Commit it to tools/ so cloud
    sessions can extend the references;
  - `hopf_formula` is private in tests/xengine.rs:156 and must move to a shared test module.
- **m5.** The Markov corner is under-counted ("about 16 composites"). It also needs:
  - the σ_1 decompositions and Rickard components on p⊗id_3, i.e. the (3,3) pairs embedded in
    (4,2,3), (5,1,3) and (6,0,3);
  - negative-sign components for [1,1,1,-2];
  - closure generators for every type reached after σ_2.
- **m6.** M6a-3 budget: 4-6 words × 1.2-2 h, plus 4_1 at 2.4 h, plus N = 5 at 1 h, comes to about
  8-14 h, i.e. 2-3 sessions. At the high end of τ this triples.
- **m7.** review_m4b/certify.md F1 (major): pairs are certified when L_q is decomposable. Until
  hardening lands, `DecHom` must itself assert `expand(ladder_char(lq)) == [(q, 1)]`. `idem_needed`
  classifies the types but is not a per-call guard.
- **m8.** "replay F32b along F31's tier choices": homprov has no replay hook; only
  `Search::tiers(lo, hi)` is public. The 0.3-0.5x saving is speculative.
- **m9.** 4_1 type counts disagree: "38 types" (SPEC, R4 §9) against "34" (R3 §4.1, R4 §5.1). T0
  coverage depends on which is right.
- **m10.** R3 §2 keeps dense interpolation for 2 strands, while R4 §8 drops §6.5 "for all of M6a".
  Settle this; box points are cheap for (3,3).

## 4. Suggested re-ordering

| wave | sessions (parallel) | result |
|---|---|---|
| 0 | done (H) | trefoil 136, Hopf, N = 3..6 |
| 1 | M2S' (memory first, time-boxed, gate C2(d)); M5a; MBOX (+ rank stage, ClosureDataC wrapper); M6a-1a (engine without library) | sigma^5 exact if M2S' passes |
| 2 | M6a-1b (DecData wiring + 2-strand box loop); M5b (+ large-X sample); M2S-run if gate passed | **T(2,5..9) at gl_6 = Wang** (C3), both primes |
| 3 | M6a-2a (driver, T0, Markov words N = 4..6, prune-4 duality on 185 knots); M5b-run | first 3-strand gl_6 checks |
| 4 | M6a-2b (S3 at prune 6) → M6a-3 (4_1) | first new value |
