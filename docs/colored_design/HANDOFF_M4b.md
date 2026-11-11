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
