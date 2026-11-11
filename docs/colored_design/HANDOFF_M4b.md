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
