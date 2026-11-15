# C2-7 report: unblocking the (2,2,2,2) class-A library jobs and the first 4-strand 2-colored gl_4 knots (branch colored-2c4-unblock)

Date 2026-11-14. Machine: 4 cores (Xeon, 2.10 GHz reported; cost table measured at 2.30 GHz, so every run uses
`KRHOM_COSTS_ANY_HOST=1`), 15 GB. Base: origin/colored-2c4-m5i 73cbd7c merged with origin/dev f525aec (B5 design
study). Tags: [M] measured, [E] estimate. Outputs: runs_c2/c27/, tests/c27_joint.rs; library delta on data-gl6
under c2/l2222_n4_c27/.

## Summary

- **All 6 class-A blockers certified for both primes, run = replay** (pairs 1235, 1444, 1470, 1894 and the M5i
  pairs 1769, 1812), without weakening certification: every Lemma A exact call ≤ 86 967 unknowns (inside the
  10⁵ budget, no cost record needed), no cap or check changed, every map re-certified on load in the replay.
  Routes (src/leaves.rs, `KRHOM_LEAVES_C27=1`, off by default): **joint junction words** (B5's common-prefix idea:
  1235's 83.7 M step → 9 584), **W1/W2 Prop-W checks** in the planner score (they expose windows the old score
  wrongly reported as small), **per-predecessor final words** with rex chains (1444, 1470: two predecessors of the
  final type, C_B5 F2), and **Leaves pairs for e_q targets** (1769, 1812: the rerouted X is itself a reduced word of
  the decomposable target; exit chain + e_q wrap). Later added: per-step detours at every level (class B).
- **All-steps census** of the executor's own plan (`gp0 --steps`), replacing gp0's stale worst-step string.
- **Cover: all 18 class-A knots of knots4_cloud.txt have 0 missing jobs** (both primes; 450 + 93 jobs each).
- **First genuine 4-strand 2-colored (∧²) gl_4 knots: 6_1, 7_2, 7_4, 7_6, 7_7, reduced, both primes, all checks PASS**
  (χ = decat_chi, D14 scaling, F31 = F32b, N' = 3 equal to engine2 gl_3 of the reversed braid, N' = 2 rank 1).
  Reduced N' = 4 ranks 81, 121, 225, 361, 441 (= det²). Unreduced runs (1 296 box points, ≈ 4 h per knot and prime here)
  were not done.
- Engine: column-restricted template localization for large windows (identical values; generator pass −23 %).
- Class B: 7 of its 11 missing jobs certified for both primes, run = replay (one with a measured cost record); the
  4 left (1412, 1435, 1472, 1483) stop on one 2.77 M step into the decomposable type T5.
- Layer-A delta on data-gl6 c2/l2222_n4_c27 (both primes); full suite: 290 passed, 0 failed, 37 ignored.

## 1. Setup

- Layer A restored from data-gl6 in the order C2-2 (c2/l2222_n4) → C2-5 (c2/l2222_n4_c25) → C2-4 delta
  (c2/l2222_n4_c24): 5 871 files for F31, 5 840 for F32b (2 861 / 2 830 layer-A + templates), cost table
  c2/costs_2222_c25.txt. All jobs: `libbuild one --lam 2,2,2,2 --n 4 --field F --cache C --costs costs_2222_c25.txt
  --provider m4c-lf --stamp-route [--replay] -- JOB` (runs_c2/c27/lib/runjobs.sh; ulimit -v 14 GB; 1 thread).
- The 6 blockers fail in these provider calls (C2-4/C2-5 logs; X = L_P·R, target q, degree s):

  | pair | failing call(s) (first call; reroute) | C2-5 / C2-4 Leaves verdict |
  |---|---|---|
  | 1235 | L(0,1,0,1/1,0,1,0/1,1,0,0/0,1,0,1)·F1.1 → 0,0,2,0/1,0,0,1/1,1,0,0/0,1,0,1, s 0 | step 83.7 M > cap (Search exhausted) |
  | 1444 | L(0,1,0,1/1,0,1,0/1,1,0,0/0,0,1,1)·E1.2 → 0,2,0,0/1,0,0,1/1,1,0,0/0,1,0,1, s −1, 0, 1; reroute L(0,1,1,0/1,0,1,0/1,1,0,0/0,0,1,1)·E1.2 | step 38.0 M; reroute 5.4 M |
  | 1470 | L(0,1,0,1/1,0,0,1/1,1,0,0/0,0,2,0)·E1.2 → 0,2,0,0/1,0,0,1/1,1,0,0/0,1,0,1; reroute L(0,1,1,0/1,0,1,0/1,1,0,0/0,0,2,0)·E1.2 | 60.2 M; reroute 1.77 M |
  | 1894 | L(0,0,1,1/1,0,0,1/1,1,0,0/0,1,1,0)·E2.1 → 0,0,2,0/1,0,0,1/1,1,0,0/0,1,1,0, s −1; reroute L(0,1,0,1/1,0,1,0/1,1,0,0/0,1,1,0)·F1.1, s 0 | 184.6 M; reroute 282.7 M |
  | 1769 (M5i) | reroute L(1,1,0,0/1,0,0,1/0,2,0,0/0,1,1,0)·F1.2 → T6 = 1,0,1,0/1,0,0,1/0,1,1,0/0,1,1,0 (e_q) | Leaves pairs refused (e_q not implemented); Exact 10.8 M |
  | 1812 (M5i) | reroute L(1,1,0,0/0,1,0,1/1,1,0,0/0,1,1,0)·F1.2 → T7 = 1,0,1,0/0,0,1,1/1,1,0,0/0,1,1,0 (e_q) | as 1769; Exact 63 M |

## 2. All-steps census of the executor's own plan (work item 1)

`gp0 --lam 2,2,2,2 --steps "P|R|Q;..." [--joint]` (new, src/bin/gp0.rs `steps_census`) prints **every** Lemma A step
of the plan that `build_leaves` executes (`Leaves::plan_words`, then the per-predecessor finals of §3.3), with the
window (common prefix / suffix, labels), whole and window unknowns per degree of supp G, the W1/W2 flag and the
W-aware score of §3.2, the cost-record status of every call above 10⁵, and the exit chain W_q → L_q (B2 planner dry
run). It replaces gp0's single stale `worst` string (C_B5 F1a: gp0.rs:174-177) by the full list; files
runs_c2/c27/steps_all.txt (default plan), steps_all--joint.txt, steps_hard_pp.txt (joint + per-predecessor finals).

Max step (unknowns, both directions, max over supp G) per failing call [M]:

| call | default plan (step_score) | default, W-aware | joint (W-aware) | joint + per-predecessor finals |
|---|---|---|---|---|
| 1235 | 83 720 537 (1 step > 10⁵) | 83 720 537 | **81 747** | 81 747 |
| 1444 first | 7 390 521 (2) | 7 390 521 | 3 674 403 | **86 967** |
| 1444 reroute | 86 967 (0!) | **11 064 924** | 1 045 467 | 1 045 467 |
| 1470 first | 11 720 095 (2) | 11 720 095 | 3 674 403 | **86 967** |
| 1470 reroute | 3 468 296 (2) | 3 468 296 | 1 771 135 | **86 967** |
| 1894 first | 731 246 505 (2) | 3 172 681 034 | 3 172 681 034 | (not needed) |
| 1894 reroute | 54 979 129 (1) | 54 979 129 | **86 967** | 86 967 |
| 1769 / 1812 | 0 (X is reduced for q: every step Id) | 0 | 0 | 0 |

Findings:
1. **The default score hides incomplete windows.** For 1444's rerouted call `step_score` reports 86 967 (a
   suffix window on (2,2,2)), but the run failed with a 5.4 M whole-word call: the window E0.1…F0.1 → E1.1…F0.2
   shares the suffix E2.1·E1.2, i.e. the step rung ρ sits in the common suffix, the window family has too small a
   rank, and `step_family` falls back to the whole words. The W-aware score (§3.2) predicts the actual run
   (11.1 M) — census = executor.
2. **Hidden final predecessors** (C_B5 F2) are what blocks 1444 and 1470: two predecessors of q at the last level,
   one Id into W_q, the other needing a different final word. A common W_q cannot serve both; per-predecessor
   final words (§3.3) take both below 10⁵.
3. For the M5i pairs the rerouted X is **itself a reduced expression of the decomposable target** (X = B_q ⊕
   B_r): the Lemma A plan is all-Id; only the exit chain X → L_q (B2 planner, 4 moves, exact ≤ 47 930 unknowns,
   templates ≤ 9 flows) and the idempotent e_q are needed.

## 3. Implementation (work item 2; src/leaves.rs, off by default)

Switch: `KRHOM_LEAVES_C27=1` in `provider_for_stamped_m4c` sets `Leaves::joint` and `Leaves::eq_pairs` and appends
the route token `c27=joint(sfx4,keep40,w12)+eqpairs` (`C27_STAMP`) to the provider stamp. Without it the provider,
its stamp and every plan are byte-identical to C2-4 (`joint` and `eq_pairs` default false; `plan_words` only calls
`joint_descent` when `joint`). `sindec::route_header` drops the `c27=` token (the only change to sindec.rs), so with
`--stamp-route` C2-7 files load in C2-4/C2-5 configurations and vice versa; every loaded map is re-certified on
load, as for the cost-table and M4c-lf route tokens.

### 3.1 Joint junction words (`Leaves::joint_descent`, B5 §2.1 "adopt")

After the existing coordinate descent, for every step t' -ρ-> t whose W-aware score exceeds the budget (10⁵):
candidate pairs (W_{t'}, W_t) = (P·A, P·B) with a common prefix P over the prefixes of the current words of t' and
t, A, B suffix words of ≤ 4 rungs (`suffix_words`: BFS from the top labels of P, ≤ 60 000 suffixes, reduced for t' /
t by `reduced_for`); the 40 pairs with the smallest score of this step are evaluated on the plan objective (max
W-aware score over every step into or out of t' and t, then flows) and the best is adopted if it improves it.
3 sweeps. The plan structure (one word per junction type per level) is unchanged, so **no detour maps** are
needed; the identity type keeps ∅; the final word of a non-identity target may change (the exit chain absorbs it).
This is the B5 census's "common prefix" idea, found here by search over prefixes of the current words instead of
the Grassmannian prefix (which for (2,2,2,2) would need a 5th, pad upright); the windows found are on 3-4 uprights
(labels (2,2,2), (3,1,2), (4,0,2,2), (4,1,1,2)).

### 3.2 Prop W checks (`window_w12`, `step_score_w`)

W1 + W2 of D_grass §1.3: the window ends in ρ (no common suffix), A' and B' are reduced for local types a', b on the
window labels ν (`local_type`: the table t of char w with w reduced for t, on `Schur::new(ν)`), and b occurs in
step_gf(a', ρ_loc) with G_{a'→b} = G. The joint planner scores a window only when it passes W1 + W2, else the whole
words (what `step_family` falls back to). W3/W4 (flow-level) are not checked a priori: `step_family`'s rank test
decides at run time, and a failure is Exhausted/Refused, never a wrong basis. Planning only — acceptance is
unchanged.

### 3.3 Per-predecessor final words (`Leaves::final_pred_words`, B5 §2.1 item 5, C_B5 F2)

With `joint`, for every predecessor t'_i of a final non-identity type q whose final step is over budget: the word
W'_{q,i} reduced for q minimizing its own W-aware final-step score (candidates P·B, P over the prefixes of
W_{t'_i}·ρ, |B| ≤ 4), kept only if strictly better and if degree-0 rex chains W'_{q,i} ⇄ W_q exist (B2 planner,
`leaf_pair_plan_opts`: exact ≤ 10⁵, ≤ 12 rungs). `build_leaves` builds predecessor i's final step into W'_{q,i} and
composes it with the realized chain (Out: R∘S, In: S∘R), so every leaf lands in W_q; the chains are degree-0 maps
between reduced words of q with top scalar 1, so the top-to-top block stays invertible (the diagnostic rank check of
`build_leaves` and then `certify_pairs` re-verify it).

### 3.3b Per-step detours at every level (`Leaves::step_words`, B5 §2.1 item 4)

Added after the 6 blockers were certified (they did not need it), when the class-B job 1054 stopped on a
1.05 M step that the pairwise joint search could not lower (a local minimum: the joint move that fixed a later
step had moved W_t away from the prefix extension W_{t'}·ρ). `step_words` generalizes §3.3 to every level: for
an over-budget step t'_i -ρ-> t at level k, a word W'_{t,i} for this step only (often W_{t'}·ρ itself: an Id step)
plus degree-0 rex chains W'_{t,i} ⇄ W_t, so every leaf at t still lands in the plan word W_t. Same acceptance
(rank check at rep(t) in `build_leaves`, then certify_pairs / certify_closure). On with `KRHOM_LEAVES_C27=1`
(`KRHOM_LEAVES_C27_DETOURS=0` restricts it to the final level); `gp0 --steps ... --joint --detours` shows it
(1054's call: max step 1 045 467 → 50 221). The chains are planned twice (dry in `step_words`, again in
`build_leaves`); a cache would remove the second planning (job 949 spent most of its 31 min in the planner).

### 3.4 Leaves pairs for targets with e_q (`pairs_leaves_eq`)

`leaf_pairs` no longer refuses e_q targets when `eq_pairs`: the Lemma A leaves and exit / entry chains are built as
for an indecomposable target, then ι ↦ ι∘e_q, π ↦ e_q∘π (as the production backends' `with_idem`) before
`certify_pairs`; e_q is the identity at rep(q) (the extra summands of L_q have types < q), so the form at rep(q) is
that of the leaves. The rex-class shortcut is not used for e_q targets. Acceptance unchanged: certify_pairs, then
sindec's M5i checks (Π∘Ι = δ with e_q on the diagonal, Σ ιπ = E, D14 scaling, layer-A save and re-certified load).

### 3.5 Engine: column-restricted localization (src/tmpl.rs, src/local.rs, src/smap.rs)

Profiling the 4-strand box run (gdb samples, `KRHOM_LPROF`, `KRHOM_GPROF`) showed template localization dominating
(e.g. 1.48 M localizations / 143 s per record scan; the generator determinant pass 8.5-11 s per box point, mostly
whole-word closure templates with 1 296 / 1 728 source fixed points). `Eval::window_block` now localizes templates
whose source window has ≥ 48 fixed points only on the source columns the coset needs
(`TmplLib::localize_cols`: L_t·Φ·(L_s⁻¹)[:, P] with `local::inverse_cols`, n³/3 + n²k instead of the full
inverse and two n³ products), one per distinct window configuration; `KRHOM_LOC_FULL=1` restores the old path.
Exact arithmetic, so the values are identical: point 7 of 6_1 gives the same values hash (c504c75f99310de4) with
both paths. The generator pass drops from 11.1 s to 8.5 s per point; the scan itself is dominated by many mid-size
localizations and did not get faster.

## 4. The 6 blockers: certified (work items 3-4) [M]

Runs (provider on, `KRHOM_LEAVES_C27=1`) and replays (`--replay`: provider refusing every call, maps only from
layer A, re-certified on load) — runs_c2/c27/lib/c27_f{31,32b}_{run,replay}.jsonl:

| pair | route (the call that failed before) | F31 run | F32b run | F31 replay | F32b replay |
|---|---|---|---|---|---|
| 1235 | Leaves Lemma A, joint words: the 83.7 M step becomes a (2,2,2) window of 9 584 (level 8) after a joint choice of the level-7/8 words; max step 81 747 (whole, level 5) | ok 46 s | ok 80 s | ok 19 s | ok 14 s |
| 1444 | first call, joint + per-predecessor final word [E0.1,E1.2,E0.1,F1.1,E2.2,F0.2,E1.1,F2.2] (final window (4,1,1,2), 10 530) + 1-move exit chain; targets at s = −1, 0, 1 (5 calls); max step 86 967 | ok 694 s | ok 645 s | ok 24 s | ok 17 s |
| 1470 | first call as 1444 (same final step shape); max step 86 967 | ok 390 s | ok 346 s | ok 27 s | ok 17 s |
| 1894 | first call refused (3.2 G step, no joint route), reroute (Leaves, joint words: 55.0 M → 86 967) | ok 51 s | ok 64 s | ok 39 s | ok 25 s |
| 1769 | M5i: reroute → Leaves pairs with e_T6: X reduced for T6, all-Id Lemma A, exit chain (4 moves) + e_q wrap; certify_pairs + M5i checks | ok 157 s | ok 167 s | ok 29 s | ok 24 s |
| 1812 | M5i: as 1769 with e_T7 | ok 150 s | ok 153 s | ok 32 s | ok 32 s |

**6 / 6 certified for both primes; replay = run.** Every Lemma A exact call was ≤ 86 967 unknowns (≤ the 10⁵
budget: no call needed a cost record, the hard cap was never approached); the largest planner exact move 47 930.
No cap, check or admission rule was changed; the Exact fallback caps of C2-4 were not raised. (The first F31 runs
of 1444 / 1470 failed with joint words alone — the failures that motivated §3.3 — and are kept in the jsonl.)

Layer-A delta: data-gl6 c2/l2222_n4_c27/{f31,f32b} (12 layer-A + 170 template files per prime, 0.27 MB;
runs_c2/c27/persist.sh; MANIFEST with sha256).

## 5. Cover of knots4_cloud.txt (work item 5) [M]

`examples/m5i_cover` (sindec::knot_jobs, upper bound and K_0 lower bound) against the union of the replay
statuses (C2-5 full replay level + corner, C2-4 M5i replay, C2-7 replay; runs_c2/c27/status.py →
runs_c2/c27/cover/status_f{31,32b}.txt, 1 236 / 1 254 jobs ok for each prime; cover_f{31,32b}.txt, identical):

| knots | jobs | missing |
|---|---|---|
| **class A (18): 6_1, 7_2, 7_4, 7_6, 7_7, 8_4, 8_6, 8_8, 8_11, 8_13, 8_14, 9_10, 9_11, 9_13, 9_17, 9_18, 9_20, 9_23** | 450 + 93 each | **0** (both primes) |
| class B: 8_15, 9_22, 11n_42 | 555 + 93 | 11 → **4** after §8 (1412, 1435, 1472, 1483) |
| Conway | 846 + 125 | 16 → 9 after §8 (942, 1412, 1435, 1472, 1483, 1830, 1906, 2392, closure T3) |

So every class-A knot of knots4_cloud.txt is computable from layer A alone (`--replay`).

## 6. 4-strand knots (work item 6) [M]

Driver: `sbox box` (4 shards, `--ckpt-secs 600`, KRHOM_THREADS=1 each) + `sbox finish` (4 threads), prune 4, box
N = 4, closures N' = 2, 3, 4, seed 201, `--replay --provider m4c-lf --stamp-route`, library = the restored layer A
plus the C2-7 delta, `KRHOM_MEM_CAP_MB=200000` (raised), `ulimit -v 14 GB` as the real guard
(runs_c2/c27/k4/box.sh, finish.sh, queue.sh; outputs runs_c2/c27/k4/out/; comparison runs_c2/c27/k4/compare.py →
compare.json). Library per run: 94 types (all closable), 2 idempotent types visited (T1, m = 1 290; T2, m = 516),
1 119 layer-A loads, 0 provider answers (replay), hash c40342369384a331 for 6_1.

**Reduced-only runs** (`--modes r`): a reduced closure needs only the box points with μ⁰ = ∅ (216 of 1 296). Cost
per run [M]: library 140-220 s per process, generator determinant pass 54 points × 8.5 s per shard, scan
26-40 s per point (F31 ≈ F32b; 30-40 % slower while class-B library jobs ran beside), finish 12-14 min
(determinant pass over all 216 points again, 4 threads, close ≈ 2 min): **≈ 45-60 min wall per (knot, prime)** for 7-letter words; 9-letter words (7_2, 7_4) need ≈ 3.4 GB per shard, so 3 shards (4 were OOM-killed by the 14 GB limit), 35-42 s per point, ≈ 1.2 h per (knot, prime). A container restart interrupted 7_2 F32b; the shards resumed from their checkpoints (box.sh is resumable).

Not run (queued, runs_c2/c27/k4/queue2.txt): 8_11, 8_14, 8_4, 8_6 and the other class-A knots (all covered by the library).

| knot | braid (knots4_cloud) | N'=4 reduced | N'=3 reduced = gl_3 reversed | N'=2 | χ | scaling | F31 = F32b | status |
|---|---|---|---|---|---|---|---|---|
| **6_1** | 1 1 2 −1 −3 2 −3 | rank 81 | rank 9, eq | rank 1 | yes | yes | yes | **PASS** |
| **7_6** | 1 1 −2 1 3 −2 3 | rank 361 | rank 19, eq | rank 1 | yes | yes | yes | **PASS** |
| **7_7** | −1 2 −1 2 −3 2 −3 | rank 441 | rank 21, eq | rank 1 | yes | yes | yes | **PASS** |
| **7_2** | 1 1 1 2 -1 2 3 -2 3 | rank 121 | rank 11, eq | rank 1 | yes | yes | yes | **PASS** |
| **7_4** | 1 1 2 -1 2 2 3 -2 3 | rank 225 | rank 15, eq | rank 1 | yes | yes | yes | **PASS** |

Reduced N' = 4 Poincaré polynomials (t = homological, q = quantum; full records in gl4_2222.json, key `c2_7`):
- 6_1 (rank 81): t^-4q^12 + t^-3q^8 + t^-3q^10 + t^-2q^4 + 3t^-2q^6 + t^-2q^8 + t^-1q^0 + 4t^-1q^2 + 3t^-1q^4 + 3t^0q^-2 + 6t^0q^0 + t^0q^2 + 2t^1q^-6 + 6t^1q^-4 + 4t^1q^-2 + t^2q^-10 + 4t^2q^-8 + 6t^2q^-6 + t^2q^-4 + 2t^3q^-12 + 5t^3q^-10 + 3t^3q^-8 + t^4q^-16 + 3t^4q^-14 + 4t^4q^-12 + t^4q^-10 + t^5q^-18 + 3t^5q^-16 + 2t^5q^-14 + t^6q^-20 + 2t^6q^-18 + t^7q^-22 + t^7q^-20 + t^8q^-24
- 7_6 (rank 361): t^-4q^4 + 2t^-3q^0 + 2t^-3q^2 + t^-2q^-4 + 6t^-2q^-2 + 3t^-2q^0 + 7t^-1q^-6 + 9t^-1q^-4 + 2t^-1q^-2 + 4t^0q^-10 + 15t^0q^-8 + 9t^0q^-6 + t^0q^-4 + t^1q^-14 + 14t^1q^-12 + 19t^1q^-10 + 6t^1q^-8 + 8t^2q^-16 + 24t^2q^-14 + 15t^2q^-12 + 2t^2q^-10 + 2t^3q^-20 + 19t^3q^-18 + 24t^3q^-16 + 7t^3q^-14 + 9t^4q^-22 + 26t^4q^-20 + 14t^4q^-18 + t^4q^-16 + 2t^5q^-26 + 17t^5q^-24 + 19t^5q^-22 + 4t^5q^-20 + 6t^6q^-28 + 17t^6q^-26 + 8t^6q^-24 + t^7q^-32 + 9t^7q^-30 + 9t^7q^-28 + t^7q^-26 + 2t^8q^-34 + 6t^8q^-32 + 2t^8q^-30 + 2t^9q^-36 + 2t^9q^-34 + t^10q^-38
- 7_7 (rank 441): t^-8q^24 + 2t^-7q^20 + 2t^-7q^22 + t^-6q^16 + 6t^-6q^18 + 3t^-6q^20 + 8t^-5q^14 + 10t^-5q^16 + 2t^-5q^18 + 5t^-4q^10 + 17t^-4q^12 + 10t^-4q^14 + t^-4q^16 + t^-3q^6 + 17t^-3q^8 + 22t^-3q^10 + 6t^-3q^12 + 10t^-2q^4 + 29t^-2q^6 + 17t^-2q^8 + 2t^-2q^10 + 3t^-1q^0 + 24t^-1q^2 + 29t^-1q^4 + 8t^-1q^6 + 12t^0q^-2 + 32t^0q^0 + 17t^0q^2 + t^0q^4 + 3t^1q^-6 + 22t^1q^-4 + 24t^1q^-2 + 5t^1q^0 + 9t^2q^-8 + 22t^2q^-6 + 10t^2q^-4 + t^3q^-12 + 12t^3q^-10 + 12t^3q^-8 + t^3q^-6 + 3t^4q^-14 + 9t^4q^-12 + 3t^4q^-10 + 3t^5q^-16 + 3t^5q^-14 + t^6q^-18
- 7_2 (rank 121): t^0q^-8 + t^1q^-12 + t^1q^-10 + t^2q^-16 + 3t^2q^-14 + t^2q^-12 + t^3q^-20 + 4t^3q^-18 + 3t^3q^-16 + t^4q^-24 + 4t^4q^-22 + 6t^4q^-20 + t^4q^-18 + 3t^5q^-26 + 7t^5q^-24 + 4t^5q^-22 + 2t^6q^-30 + 6t^6q^-28 + 7t^6q^-26 + t^6q^-24 + t^7q^-34 + 4t^7q^-32 + 7t^7q^-30 + 4t^7q^-28 + 2t^8q^-36 + 5t^8q^-34 + 6t^8q^-32 + t^8q^-30 + t^9q^-40 + 3t^9q^-38 + 5t^9q^-36 + 3t^9q^-34 + t^10q^-42 + 3t^10q^-40 + 4t^10q^-38 + t^10q^-36 + t^11q^-44 + 3t^11q^-42 + 2t^11q^-40 + t^12q^-46 + 2t^12q^-44 + t^13q^-48 + t^13q^-46 + t^14q^-50
- 7_4 (rank 225): t^0q^-8 + 2t^1q^-12 + 2t^1q^-10 + 3t^2q^-16 + 6t^2q^-14 + t^2q^-12 + 2t^3q^-20 + 8t^3q^-18 + 6t^3q^-16 + t^4q^-24 + 8t^4q^-22 + 12t^4q^-20 + 2t^4q^-18 + 6t^5q^-26 + 14t^5q^-24 + 8t^5q^-22 + 3t^6q^-30 + 12t^6q^-28 + 14t^6q^-26 + 3t^6q^-24 + t^7q^-34 + 7t^7q^-32 + 14t^7q^-30 + 8t^7q^-28 + 3t^8q^-36 + 10t^8q^-34 + 12t^8q^-32 + 2t^8q^-30 + t^9q^-40 + 5t^9q^-38 + 10t^9q^-36 + 6t^9q^-34 + t^10q^-42 + 5t^10q^-40 + 7t^10q^-38 + t^10q^-36 + 2t^11q^-44 + 5t^11q^-42 + 3t^11q^-40 + 2t^12q^-46 + 3t^12q^-44 + t^13q^-48 + t^13q^-46 + t^14q^-50

**Unreduced runs were not done**: they need all 1 296 box points, i.e. ≈ 1 296 × 27 s / 4 shards ≈ 2.4 h of
scanning plus ≈ 46 min of determinant passes per shard and another ≈ 46 min in finish — about 4 h per
(knot, prime) on this machine [E from the measured per-point costs]. The 4-strand unreduced values (6_1 first)
are the next runs (`MODES=u,r runs_c2/c27/k4/box.sh 6_1 31 "1 1 2 -1 -3 2 -3" 4`, resumable).

## 7. Checks

- Library acceptance unchanged: every new map passed sindec's checks (`dechom_postcheck`, D14 scaling, Π∘Ι = id or
  δ with e_q on the diagonal, Σιπ = E, kept multiset; Leaves pairs additionally `certify_pairs` and the
  `build_leaves` rank check at rep(t) of every level) and was re-certified on load in the replays (provider
  refusing every call). No cap, budget, admission ceiling or check was changed; one cost record pair was added by
  measurement (§8, job 986; not used by any class-A job).
- Knots: χ = decat_chi on every closure (fatal), D14 scaling at 2 box points (fatal), the generator determinant pass
  at every box point (in every shard and again in finish), F31 = F32b closure by closure, N' = 3 = engine2
  uncolored gl_3 of the reversed braid (**equal as polynomials**, no shift, for every knot), N' = 2 rank 1.
- Engine change (§3.5) checked by identical point-value hashes with `KRHOM_LOC_FULL=1` (old path) on 6_1.
- Unit tests: tests/c27_joint.rs (3 tests: 1235's plan ≤ 10⁵ only with joint words, 1444 needs a per-predecessor
  final word, W1/W2 rejects the suffix window of 1444's rerouted call); full suite below.
- Observation (not used as a check): the reduced N' = 4 ranks are det(K)² (6_1: 9², 7_2: 11², 7_4: 15², 7_6: 19², 7_7: 21²), as for
  every earlier (2,2,2) / Markov-stabilized value (trefoil 9 = 3², 4_1 25 = 5², 5_2 49 = 7²); the reduced
  N' = 3 ranks equal det(K).
- Full suite `KRHOM_THREADS=4 cargo test --release --no-fail-fast` (runs_c2/c27/full_suite.log): **290 passed, 0 failed, 37 ignored** (C2-4: 287 + the 3 new tests).

## 8. Class B (side work while the knot runs used the cores) [M]

The 11 jobs missing for 8_15, 9_22, 11n_42 were run with the C2-7 provider (one job at a time beside the box
shards; runs_c2/c27/lib/runjobs_classb*.log, classb_f32b.sh):

| job | result | route / blocker |
|---|---|---|
| 1471, 926, 949, 997, 1054, 1540 | **certified F31 + F32b, replay ok (both primes)** | Leaves Lemma A (joint words; 1054 / 1540 with a per-step detour at level 7, §3.3b); 926, 949, 997 (C2-4's T4-target jobs) through the e_q Leaves pairs; 949 / 997 spend 15-40 min in the rex planner |
| 986 | **certified F31 + F32b, replay ok** | one Lemma A step of 271 848 unknowns (2,2,2,2) [E1.2,E0.1,F1.2,F0.1,E2.1,E1.1] ⇄ [E1.2,F1.1,F0.1,E0.1,E2.1], d = 0, admitted by two measured records (`libbuild calib --windows`, 85.0 s / 1 271 MB and 48.5 s / 1 008 MB, limit 125 s, uncensored; runs_c2/c27/calib, table costs_2222_c27.txt = c25 + 2 records) — the B4w admission rule unchanged |
| 1412, 1435, 1472, 1483 | fail (both runs) | a Lemma A step into T5 = 1,1,0,0/0,1,0,1/1,1,0,0/0,1,0,1 of 2.77 M unknowns (whole words; no joint pair, no target-side detour lowers it; runs_c2/c27/steps_1412.txt) |

Cover after this (both primes, identical): class B misses 4 jobs (1412, 1435, 1472, 1483), Conway 9.

## 9. Open problems

1. **Unreduced 4-strand runs** (all 1 296 box points) cost ≈ 4 h per (knot, prime) here [E]; 6_1 unreduced is
   the first to run (box.sh is resumable). The engine time is template localization (≈ 1.5 M localizations per
   scan point, mostly mid-size windows of 36-108 fixed points); a cache of localizations keyed by the window's
   variable values, or a cheaper evaluation of the e_p-wrapped closure maps, would be the next speed-up.
2. **Class B** (8_15, 9_22, 11n_42): 7 / 11 missing jobs now certified for both primes, run = replay (1471, 926,
   949, 997, 1054, 1540, 986; §8); 1412, 1435, 1472, 1483 stop on the same Lemma A step into the decomposable
   target T5 = 1,1,0,0/0,1,0,1/1,1,0,0/0,1,0,1 ([E1.2,F1.1,F0.1,E0.1,E2.2,F2.2]·F1.1 → [E1.1,E0.1,F1.1,F0.1,E2.1,
   F2.1], 2.77 M unknowns; neither a joint pair nor a target-side detour lowers it: runs_c2/c27/steps_1412.txt).
   A source-side detour (W_{t'} → W'_{t'} by a rex chain, Kron with ρ) is the missing move.
3. The joint search is a local search (prefixes of the current words, suffixes ≤ 4 rungs, pairwise moves); a
   bottleneck-path search over the levels (candidate words per type, minimize the max step) would find plans the
   pairwise moves miss (1054 needed a detour instead).
4. W3/W4 of Prop W are not checked a priori (the runtime rank test decides; no incomplete W1/W2-certified window
   was met in these runs).
5. Conway (16 missing) and the closure of T3 (C2-4) were not attempted; the closure T3 needs a Leaves closure
   route through e_p. (Several C2-4 T4-target jobs — 926, 949, 986, 997 — are now certified for F31 by the
   e_q Leaves pairs.)
