# M4b handoff (2026-11-11, laptop -> cloud)

Design documents in this folder were written in a local scratchpad; absolute /private/tmp/... paths in
them refer to that scratchpad and do not exist elsewhere. Read in this order: SYNTHESIS.md (design D1),
SPEC_M4_M6a.md (rev 2; section 4 = M4b, 4.7 = gate G1, section 5 = M5), M4b_design.md (implementation
notes, section 5 = changes made while implementing), NTRUNC_STUDY.md (Foam_N pruning rule).

## State of branch colored-m4b (WIP, not reviewed)

- src/homprov.rs: Hom provider (Exact + Search tiers S0-S6), certification against the EKLP Hom formula,
  template cache, cost-table admission. src/bin/g1.rs: gate G1 harness (`g1 a`, `g1 one`, `g1 closure`,
  `g1 calib`). tests/homprov.rs. M4a review minors (a)-(g) fixed in smap.rs / tmpl.rs / tests/m4a.rs.
- Not compiled since the WIP commit was made (it compiled and ran on the laptop at the time).
- Cost table calibrated on an Apple M5 laptop: tests/data/homprov/costs_m4b.txt (the cloud box is
  about 1.5x slower; admission caps may need a recalibration with `g1 calib`).

## Gate G1(a), laptop run (stopped for heat after 63 of 110 jobs, 1791 s)

Command: `g1 a --cache DIR --costs tests/data/homprov/costs_m4b.txt --wall 600 --rss 1024 --json OUT`
(cold template cache). Partial result: tests/data/homprov/g1a_laptop_partial.json.

- 243 targets reached, 233 certified; 59 of 63 jobs certified. All (3,3) jobs certified.
- Failures, all on lambda = (2,2,2): Exact refused (caps: fibre product 14400 / 1600, unknowns
  188680 / 5000, or cost-table prediction 60 s > 10 s) and Search "tiers S0-S6 exhausted (2 stuck branches)":
  - pair 18, p = 0,0,2/1,1,0/1,1,0, x = E0(2) F0(3) E1(2) E0(1) F1(2) E0(1):
    q = 0,0,2/1,1,0/2,0,0 (s = -1, 1) and q = 0,0,2/2,0,0/1,1,0 (s = -1, 1);
  - pair 25, p = 1,0,1/2,0,0/0,1,1, x = E0(2) F0(1) E1(2) F1(2) F0(1): q = 0,1,1/2,0,0/0,1,1 (s = 0);
  - pair 32, p = 0,0,2/2,0,0/1,1,0, x = E0(2) E1(2) E0(1) F0(2) F1(2) F0(1):
    q = 0,1,1/1,1,0/1,0,1 and q = 0,0,2/1,1,0/1,1,0 (s = 0);
  - pair 20: killed at the 600 s wall limit (Time).
- Not yet run: jobs 64-110 and the 4 + 21 closure types (`g1 closure`).

## Next steps

1. Build, run tests (cargo test --release --test homprov --test m4a).
2. Rerun G1(a) completely (cold cache) and diagnose the stuck Search branches on the pairs above.
3. Independent adversarial review of homprov.rs, then fixes.
4. If G1(a) passes: M5 (SPEC section 5) with Foam_N pruning (schur::Schur::killed_at) built into the
   decompositions. If it fails: list the failing targets; the fallback is M4c (EKLP light leaves).

## Cloud session 2026-11-11

Machine: 4 cores, 15 GiB (no swap). Rust stable, release builds, KRHOM_THREADS=1 for the gate.
dev merged into colored-m4b first (9ca07ad; no conflicts).

### What changed (commits 1a72074, dc58d0c, d23ccff, this one)

1. **`sbim::hom_space` solver** (root cause of pairs 18 and 20). The stuck time was in single
   exact window calls, mostly degree-0 isomorphisms between words of one indecomposable type
   (leaf isomorphisms, rex moves): (4,0,2) [E1(2) E0(1) F0(2)] -> [F0(3) E1(2) E0(2)] took 93 s,
   [F0(4) E1(2) E0(3)] -> [F0(3) E1(2) E0(2)] 270 s. `Echelon::reduce` re-sorted the whole row
   after every pivot step and every one of the (many, mostly redundant) equations was reduced.
   New `linsys::kernel_lazy`: dense-accumulator reduction (rows sorted short first); after a run
   of redundant rows (kernel dim <= 128) the current kernel is evaluated on all remaining rows,
   violated rows are inserted, and the scan stops only when none is violated. Exact (kernel(all)
   ⊆ kernel(inserted); a basis of the latter annihilating every row gives equality) and returns
   the same reduced kernel basis, so template caches stay valid. tests/hom_lazy.rs compares it
   map by map with `hom_space_reference` (the old solver, kept). Now 0.20 s and 0.34 s.
2. **Caps re-measured with the new solver** (spec §4.5 caps were sized for the old one):
   - Exact label rule: unknowns 5 000 -> 40 000 (every Hom^d, d = 0..4, between canonical
     ladders of (2,1,2) and (2,2,2) with fibre product <= 1 600: <= 7.4 s, max 39 165 unknowns).
   - Search windows: fibre product 1 600 -> 6 400 (own field `max_window_fibre`), unknowns
     2 500 -> 40 000 (measured: 32x72 / 8 762 unknowns 0.8-1.0 s; 60x90 / 21 300: 1.8 s;
     120x50 / 33 134: 5.6 s). These unblocked pairs 32, 77, 82 ("no admissible decomposable
     window", the "2 stuck branches").
3. **`certify_full` orbit fallback.** If the stacked rank at the planned detecting set is below
   the dimension, retry at the full W_lambda-orbits of the common tables with M+2 and then
   2(M+2) points (recorded as detect "orbits-fallback"/"-2x"). Sound: rank = |S| at any coset set
   proves independence. Fixed b2_hom60_through_exact ((1,3) F0(1)E0(3)F0(2) -> E0(3)F0(2)F0(1),
   d = 6: exact basis, witness rank 141 < 145 at 18 points), which failed before this session.
4. **Cost table recalibrated** (`g1 calib --windows <KRHOM_WINLOG of a G1 run> --per-class 20
   --limit 60`): tests/data/homprov/costs_cloud.txt, 424 records, 1 censored ((2,2,2) 432x24,
   112 536 unknowns, d = 1, > 60 s). Reason: costs_m4b.txt holds old-solver timings; e.g. its
   censored 60 s record for (2,4,0) 32x72 vetoed a window that now takes 1 s. Keep using
   costs_cloud.txt (or recalibrate) with the new solver.
5. Tracing: KRHOM_TRACE=1 also logs every computed exact call (time, fibres, unknowns) and the
   windows `choose_window` refuses; KRHOM_WINLOG=PATH logs every exact window call.

Certification rules are unchanged otherwise (|S| = Hom formula, stacked rank, fresh points,
genuineness; pairs: form rank = mult at two points; closure: degree multiset and nonsingular
identity-coset matrix at two points).

### Tests

- Before the fixes: homprov 8 passed, 2 failed (b2_hom60_through_exact: certify false negative;
  b2_type_pairs_22_212_through_chain: 560/605, Exact refused by the 5 000 cap, Search
  incomplete; it took 659 s).
- After: full suite (`cargo test --release`) 176 passed, 0 failed, 17 ignored; homprov, m4a,
  hom_lazy rerun on the final code: homprov 10 passed / 3 ignored, m4a 25 / 5 ignored, hom_lazy 1.
  B2 (2,1,2): 605/605 in 62 s. (The full suite ran before the last Search-cap change, which only
  touches Search admission; homprov/m4a/hom_lazy were rerun after it.)

### Gate G1(a): PASSED

Command: `g1 a --cache /tmp/tmpl_g1a (cold) --costs tests/data/homprov/costs_cloud.txt --wall 900
--rss 4096 --json tests/data/homprov/g1a_cloud.json` (`g1 a` runs the pair jobs and the closure
jobs of both lambdas).

| lambda | pair jobs | (q,s) targets | Σ mult | closure types | closure targets (ins+outs) / Σ n | wall |
|---|---|---|---|---|---|---|
| (3,3) | 12 / 12 | 72 / 72 | 121 | 4 / 4 | 8 / 8, n 16 | 9.1 s |
| (2,2,2) | 98 / 98 | 358 / 358 | 404 | 21 / 21 | 42 / 42, n 96 | 523.5 s |

Total 480 / 480 targets certified, 135 / 135 jobs, 533 s wall (spec prediction <= 30 min).
Peak RSS 856 MB (job one 2,2,2 116). Slowest jobs: one 2,2,2 56 (76 s), 53 (28.6 s), 119
(28.6 s), 116 (28.6 s), 60 (27.0 s), 82 (26.8 s). Tiers per target: Exact 264, S1 150, S2 60,
S1+S2 4, S0 2. No Time, Memory or Exhausted outcome; no rerun needed.

Robustness check: the same gate without any cost table (Exact by the label rule only, Search by
its caps) also certifies 480 / 480 (135 / 135 jobs) in 412 s.

Laptop failures, now: pair 18 certified (was Time / Exhausted), 20 (was 600 s Time) about 3 s,
25 and 32 certified by S1. Pre-session code on this box (old solver, laptop table): pairs 18
and 20 hit the time limit (first 42 jobs); with the new solver but old caps, pairs 77 and 82
failed with stuck branches.

### Remaining problems / notes

- The admission caps are machine- and solver-specific; recalibrate (`g1 calib`) after solver
  changes. costs_m4b.txt is obsolete for the new solver.
- `Search` S5 is still not implemented (spec deviation recorded in M4b_design §2).
- Not done in this session (per instructions): independent adversarial review of homprov.rs,
  G1(b) (queue on the cloud: `g1 b --cache DIR --costs tests/data/homprov/costs_cloud.txt
  --json OUT`, seed 2026), M5.

## Hardening session 2026-11-11

This session works through the independent review (docs/dev/colored_design/review_m4b/certify.md,
solver.md, tests.md). The review found no defect that gives a wrong result today. The items
below are preconditions and hardening, so that M5 cannot misuse the provider without noticing.
Machine: 4 cores (Xeon 2.10 GHz), 15 GiB.

### What changed, per item

1. **Precondition for `pairs`** (review F1). `homprov::pairs_precheck` runs in Exact and in
   Search before any other work. A failure is returned as `Refused`. It checks:
   - X and L_q have the top of q.
   - With `eq = None`: `expand(ladder_char(L_q)) == [(q, 1)]`. Otherwise the call is refused
     with "L_q decomposable: e_q required" (this also catches an L_q of another type).
   - With `eq = Some(e)`: `check_idempotent`, the e_q/e_p contract below.
2. **`mult` validation** (review F2). `mult` must equal [v^sdeg] of the multiplicity of b_q in
   char X. Otherwise the call is `Refused` and the message gives both numbers.
3. **Ends and degrees** (review F3).
   - `certify_pairs` checks every ι = (L_q -> X, s) and every π = (X -> L_q, -s).
   - `certify_closure` checks that the ends are R -> L_p (ins) or L_p -> R (outs).
   - `detect_plan` asserts that the Schur and Arena lambdas agree.
   - `TmplLib::instance()` is a nonce, unique per process. `Search` binds its splits/comms
     TmplId caches to it and clears them when it sees another library.
   - `certify_pairs` and `certify_closure` are now `pub`, for tests.
4. **Closure with e_p** (review F4, SPEC §4.3 item 7). This is implemented, not refused.
   - `closure_precheck`: without e_p, a decomposable L_p is `Refused`; with e_p,
     `check_idempotent` must pass.
   - With e_p, every selected map must be the DAG node `Comp(e_p, g)` (ins) or `Comp(g, e_p)`
     (outs), so it lies in Hom(R, B_p) or Hom(B_p, R).
   - The n x N identity-coset matrix must have rank n at two points. This is the same as
     nonsingularity in reduced coordinates of Im e_p, and the determinant argument applies to
     the free rank-n module Hom(R, B_p).
   - The S2 candidates Dot∘g are wrapped with e_p again, and so is the S0 identity.
   - The `Selector` accepts blocks of width N.
   - Test `closure_with_idempotent_313`: in (3,1,3), the type 2,0,1/0,1,0/1,0,2 has
     char L_p = b_p + b_{3,0,0/0,0,1/0,1,2} (48 flows). e_p = P(a) is computed from an Exact
     End^0 basis (the SPEC §5.2 algorithm, in the test). Both Exact and Search certify ins and
     outs. A family without the e_p wrapper is rejected.
   - A scan of the minimal-flow ladders of (λ, λ) found these decomposable types:
     - (2,3,2): 1 of 25 (640 flows);
     - (3,1,3): 1 of 16 (48 flows);
     - (2,2,2,1): 2 of 120;
     - (3,3,3): 1 of 55. This is 0,1,2/1,1,1/2,1,0, 10 500 flows, char = b_p + b_{1,0,2/0,2,1/2,1,0}.
   - None of (2,1,2), (1,2,2), (2,2,1), (2,2,2), (3,1,2), (1,3,2), (3,2,2), (3,3,1),
     (1,1,1,1), (2,1,1,2) has one.
5. **Refused vs Exhausted** (review F5, F6).
   - `Ctx` counts stuck branches whose cause was admission (`stuck_adm`) and tiers or degrees
     whose window admission refused (`tier_refused`). The cases are: `choose_window` refused a
     decomposable window, a split/commutation/leaf iso was refused, or S6 or a closure window
     was refused.
   - A call ends `Refused` iff every cause was admission. Otherwise it ends `Exhausted`.
   - Before it reports `Exhausted`, Search `full` re-ranks every offered candidate at the orbit
     cosets with the fallback point counts (`orbit_fallback`, shared with `certify_full_plan`).
   - S6 had no window at all for Hom^d(X, X) with equal words. It now uses the whole ladder
     there. This fixed B2 (2,2) through Search alone: 43 -> 45 of 45.
   - `window_maps` distinguishes "no window" (Exhausted) from "windows refused" (Refused).
6. **Cost table** (review tests finding 3).
   - A censored record sets its bucket prediction to +infinity, so it can never admit; the
     exact-shape rule is unchanged.
   - The header is now `krhom-costs v2 alg=<COST_ALG> host=<host_stamp()>`. `COST_ALG` is
     `"hs-kl1"`; `host_stamp()` is the CPU count and model name, because the host name is
     random in cloud containers.
   - `load` rejects a v1 table, another solver, or another host (`KRHOM_COSTS_ANY_HOST=1`
     accepts another host with a warning). g1 prints the reason and runs without the table.
   - tests/data/homprov/costs_cloud.txt is converted to v2 (records unchanged; it was measured
     on this machine class). costs_m4b.txt stays v1 and is now rejected (it was already obsolete).
   - `Search::window_admitted` is exposed for tests.
7. **`TmplLib::has_hom`** counts a disk file only if its header matches exactly (version,
   TMPL_ALG, field, key).
8. **`linsys::kernel_lazy`** builds the final kernel sparsely, with one reused dense buffer and
   each vector pushed sparse at once. The dense kernel is now built only for the early-stop
   check, where its dimension is at most 128. Peak memory of the result is
   O(ncols + nnz(kernel)) instead of dim · ncols; the 40 000-unknown admission no longer
   implies dim · 320 KB.
9. **Tests and stats.**
   - `kernel_lazy_traced` returns a `KernelTrace`, and `kernel_lazy_counters()` gives
     process-wide counts.
   - tests/hom_lazy.rs asserts that both the early stop and the violated-row reinsertion are
     taken.
   - `kernel_lazy_property` compares against `Echelon::kernel` element by element: 160 random
     adversarial systems, plus no equations, full rank, and kernel dimension > 128.
   - `SearchStats` gets `window_calls` / `window_secs` for Search's own exact window calls,
     which were dropped before. `exact_calls` counts only the Exact backend, so B8's assertion
     keeps its meaning.
   - The trait method `HomProvider::failed_stats` (default None) is additive. `Chain` absorbs
     the stats of failed providers, and g1 records them for failed targets.
   - B2 records counts per composition.
   - New test `b2_type_pairs_22_212_search_alone`; results under "Test results".
   - B4b forces the witness set through the new `certify_full_at` and asserts that the
     counting assertion, not the rank test, fails. It then certifies through the fallback.
   - New `b9_*` tests: preconditions; end and degree checks; cache binding; cost-table
     admission and stamps; `has_hom` header.
10. **g1 JSON header.** The header records:
    - `run`: args, commit (with -dirty), host stamp, seed, watchdog, Budget, Exact caps and
      the cost-table stamp;
    - `finished`;
    - `complete` = finished and unrestricted (no `--only`, `--max-pairs` or `--pairs`). A
      restricted run is also flagged in the summary line.

### Test results

- Full suite, `KRHOM_THREADS=4 cargo test --release`: 183 passed, 0 failed, 17 ignored.
- B2 through Search S0-S6 alone:
  - (2,2): 45 / 45.
  - (2,1,2): 581 / 605. The 24 misses are all Exhausted and none is Time, Memory or Refused.
    They are endomorphism-like spaces of types such as 1,0,1/0,1,0/1,0,1 and
    0,0,2/1,0,0/1,1,0 in degrees 2-4; dots, S4 and S6 windows span too little. This is a
    Search completeness gap (U1). The Chain covers it with Exact.
  - The test asserts this bound, which may only improve.
- B2 through the Chain, (2,2) / (2,1,2): 45 / 45 and 605 / 605.

### Gate G1(a), cold cache, rerun on the hardened code

Command: `g1 a --cache /tmp/tmpl_g1a --costs tests/data/homprov/costs_cloud.txt --wall 900
--rss 4096 --json` at commit 0408229. Result file: tests/data/homprov/g1a_cloud_hardened.json
(`complete: true`).

- 480 / 480 targets certified, 135 / 135 jobs, 534.3 s wall: (3,3) 11.6 s, (2,2,2) 522.5 s.
- Peak RSS 856 MB (job one 2,2,2 116).
- Slowest jobs (seconds): one 2,2,2 56 (77.6), 119 (30.2), 53 (27.8), 116 (27.1), 60 (26.3).
- Tiers per target: Exact 264, S1 150, S2 60, S1+S2 4, S0 2. These are identical to the
  previous run.
- 531 Exact-backend calls and 412 Search window calls were computed. Search window calls are
  now counted.
- No precondition refused any G1(a) target. Every canonical ladder of (3,3) and (2,2,2) has
  character b_p, and every `mult` from `term_mult` matched char X.

### What M5 must provide (the e_q / e_p contract)

- **When `eq` / `ep` may be None.** Only when `expand(ladder_char(L)) == [(type, 1)]`, i.e.
  `ladder_char(L_p) == b_p` (SPEC §5.2, "e = None"). Any other call is `Refused` with
  "L_q decomposable: e_q required" or "L_p decomposable: e_p required". Search for that string
  instead of treating the refusal as a completeness failure.
- **What an idempotent must satisfy.** `check_idempotent` (pub) checks, at two points:
  - it is a genuine Arena map L -> L of degree 0 (Lin/Comp of certified maps);
  - b_p occurs exactly once, in degree 0, in char L;
  - at the witness coset rep(t) of every table t of char L, block(e) is idempotent of rank
    c_t(b_p)(1);
  - the scalar at rep(p) is 1.

  This is SPEC §5.2 step 6. Build e = P(a) from `end0` as in the test helper
  `idempotent_from_end0` (tests/homprov.rs), and run `check_idempotent` yourself before
  caching an e.
- **What the provider does with e.**
  - Pairs: it wraps ι ↦ ι∘e_q and π ↦ e_q∘π.
  - Closure: it wraps every candidate as e_p∘g (ins) or g∘e_p (outs) and certifies rank n on
    the identity coset.

  M5 must not wrap again. To use the source idempotent (e_p ⊗ id_R) of SPEC §5.2 step 3, M5
  composes it on top itself, as specified there.
- **`mult` must be exact.** Pass `mult` = [v^s] of the multiplicity of b_q in char X, where X
  is the actual word. For a decomposable L_p, `term_mult(p, [R])` (expand of b_p·[R]) is
  smaller than the multiplicity in char(L_p R) and is refused. Use expand(ladder_char(X))
  together with the source idempotent, or decompose L_p first.
- **One Search may serve several Arenas.** Its caches are cleared when the library changes.
  TmplIds of one TmplLib are never valid in another.
- **Interpreting failures.**
  - `MissKind::Refused` means admission or a precondition declined the call. Raising l_max or
    c_max will not help; raise a cap or supply e or a cost record.
  - `Exhausted` is the completeness signal.
  - Search `full` / `end0` are not complete on (2,1,2), see above. Use the Chain (Exact first
    for labels <= 2). For (3,3) and (3,3,3) the Chain has Exact only through the cost table,
    so measure `end0` there before relying on it.

### Open problems

- Search `full` completeness: 24 / 605 of (2,1,2) are Exhausted through Search alone.
- Ignored B2 run through the Chain (`cargo test --release --test homprov -- --ignored`,
  914 s), test `b2_type_pairs_33_222_through_chain`: **FAILED**.
  - (3,3): 80 / 80, all through Search (dim0 52, S6 20, S0 4, S2 4).
  - (2,2,2): 2064 / 2205 (dim0 1650, Exact 316, S6 84, S0 7, S2 7), so 141 are missing.
    - In every miss, Exact was refused by the label-rule caps (fibre product above 1 600, or
      more than 40 000 unknowns) with no cost-table entry.
    - Of the 40 misses printed, Search ended 19 Exhausted and 21 Refused (window admission).
    - The misses are high-fibre type pairs, d = 1..4, between types such as 0,1,1/1,0,1/1,1,0,
      1,0,1/0,2,0/1,0,1 and 0,0,2/1,1,0/1,1,0.
  - So Search `full` / `end0` is not complete on (2,2,2) at d >= 1 for large ladders. The
    ignored test documents this and is left failing (not weakened).
  - G1(a) does not need these spaces (pairs and closure only). M5's end0 needs d = 0: none of
    the 40 printed misses is at d = 0, but the test prints only the first 40 of 141. Before M5 relies on `full` at d > 0, add cost records or raise the caps,
    or improve Search S4 / S6.
  - The other ignored tests pass: b3_b4_faith_e1f1e1f1 and b4b_detecting_set_fallback (the
    forced witness set fails the counting assertion; the basis certifies through the fallback).
- The (3,3,3) decomposable type 0,1,2/1,1,1/2,1,0 has 10 500 flows, so its End^0 through Exact
  is not admitted by the label rule. M5 needs Search `end0` (or a cost record) for its e_p.
- costs_cloud.txt is stamped for this machine class only. Recalibrate (`g1 calib`) on another.
- S5 is still not implemented (unchanged).
