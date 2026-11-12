# M4c B2 report: theorem-level Assoc certificate, rex-chain planner, census G-P0, B4 start

Date 2026-11-12. Branch `m4c-rexplan` = origin/m4c-gen (f5fe90c) + merge origin/m4c-closure (b266a1a) + merge
origin/dev (9ed360d); both merges were clean (no lib.rs/Cargo.toml conflicts). No frozen file (smap.rs, homprov.rs,
sbim.rs, linsys.rs, tmpl.rs, sindec.rs) was edited. Fields F31 (2³¹−1) and F32b (2³²−5).

Changed/new files: `src/symfun.rs` (persistence), `src/bin/grwarm.rs`, `src/leafgen.rs` (theorem certificate, light
localization, spot check, `braid_dag_thm`, `realize_chain`), `src/smapio.rs` (`thm-assoc` record), `src/rexplan.rs`
(new: singlestep forms, planner), `src/leaves.rs` (`rex_class_pairs`, `pairs_leaves`), `src/bin/gp0.rs` (census and B4
driver), `tests/thmassoc.rs`, `tests/rexplan.rs`, data in `docs/colored_design/m4c/b2/`, `tools/gp0_summary.py`.

## 0. Summary

- **Theorem-level Assoc certificate** implemented with all four approval conditions. Cross-check (i): 9 646 theorem
  cases on small windows and the braid middle steps at (1,1,1), (1,2,1), (2,1,2), (2,2,2), (3,1,3), (2,3,2): the
  numeric certificate agrees on every one (0 disagreements). (ii): the (3,3,3) braid middle step (1 680 flows) passes
  the bimodule spot check in 17 s; corruptions are detected. (iv): marked in provenance, a registry and layer A.
- **Rex-chain planner** (src/rexplan.rs) with moves Collapse (Rex ψ), Expand (Pull φ), Assoc (numeric or theorem),
  FlagMatch (Pull identity matching), Swap, Exact (budgeted window). At small λ it connects every pair of reduced
  words of every type (452 chains, 0 failures) and the realized chains lie in the exact Hom⁰ windows (356 checked).
- **G-P0 census:** every rex chain needed at (3,3,3) was found **on 3 uprights: no `Pad` (B3) is needed** for any
  target of the N = 5 and N = 6 lists. No chain needs the theorem-level Assoc (all Assoc windows ≤ 60 flows); the
  largest template window is 90 flows, the largest chain Exact move 99 175 unknowns.
  **The binding constraint is not the exit/entry rex chain but the Lemma A lower steps** (exact step windows of
  0.4-970 M unknowns) for most non-rex pairs.
  - N = 5 (4 pairs): 322 feasible; 313 feasible after one measured record (done: 28 s / 32 s); 326, 410 infeasible
    (Lemma A steps 15.7 M and 1.9 M unknowns).
  - N = 6 (39 pairs incl. 24): 7 feasible (159, 323, 322, 71, 334, 419, 37), 4 measurable (212, 335, 411, 313:
    one step 1·10⁵-4·10⁵ without record), 28 infeasible (Lemma A steps).
- **B4 started: all 6 leaf-iso rex-class pairs certified with Leaves rex chains, both primes**: 159, 71, 323, 322,
  334, 419 at N = 6 and 322 (+ the rex-class target of 313) at N = 5. π, ι are realized rex chains, `certify_pairs`
  passes (form [[1]] at rep(q), two points, genuine), π∘ι = id at every coset of L_q. 19-213 s per pair. These were
  the M5c "leaf-isomorphism window 0.7-147 M unknowns" refusals; 322 is the N = 5 bottleneck.
  Not yet done: the downstream pipeline checks (sindec wiring, d∘d, K₀, χ, replay) — see §6.

## 1. GrTable persistence (engineering decision)

`symfun::GrTable::get` now loads `<dir>/gr_<r>_<b>.bin` (dir = `$KRHOM_GRCACHE`, default `target/grcache`,
`KRHOM_GRCACHE=off` disables; tables with b ≥ `$KRHOM_GRCACHE_MIN_B`, default 7). File: magic `KRGRTBL\0`, format
version (1), r, b, payload length, FNV-1a-64 of the payload, payload (basis, weights, e_mul, f_mul). `load` verifies
magic, version, (r, b), length, hash, and that the basis/weights equal the recomputed fibre basis; any mismatch → None
→ recompute and rewrite (atomic temp + rename). Test `symfun::persist_roundtrip_and_corruption` (round trip, wrong
(r, b), flipped byte, truncation). `grwarm r:b ...` builds and persists tables.

GrTable(3, 9) and GrTable(6, 9) were started in the background at 18:23 (`grwarm`); at report time (> 75 min) neither
had finished (RSS 3.1 GB and 1.2 GB). They were not needed by anything in this session (no 9-root template was
certified: the planner's template windows are ≤ 90 flows). Status at the final commit: see §7.

## 2. Theorem-level Assoc certificate (src/leafgen.rs)

`theorem_assoc_check(ws, wt)`: (1) equal window bottoms and tops; (2) equal multistep expressions of up-down shape
[[I ⊂ W ⊃ J]] (`is_updown`: ≤ 3 terms, the middle wall set contained in both ends); (3) both words flag-type, equal
supports (flow bijection with equal top configurations), equal shifts. Then `certify_assoc_theorem` records the
constant 1 (`Arena::certify_assoc`) and enters the pair in a process-wide registry; `add_thm_assoc` adds the Elem node
with provenance tag `[thm-assoc SLL(88)]`. Down-up shapes ([[I ⊃ W ⊂ J]]) are **not** accepted by the theorem (the
approval names R^I ⊗_{R^W} R^J); they go to the numeric certificate only.

**(i) Theorem vs numeric** (`theorem_vs_numeric`; tests/thmassoc.rs):

| set | theorem accepts, numeric passes | theorem accepts, numeric fails | theorem refuses, numeric passes |
|---|---|---|---|
| all ordered pairs of flag-type words with equal supports, ≤ 4 rungs (≤ 3 for K ≥ 5), on 11 window shapes (labels (2,1) … (2,1,2)) | 9 646 | **0** | 4 866 |
| braid middle M → M' at cells (1,1,1), (1,2,1), (2,1,2), (2,2,2) | 4 | 0 | 0 |
| braid middle at (3,1,3) (140 flows, numeric 273.5 s), (2,3,2) (210 flows, numeric 1 604.2 s) | 2 | 0 | 0 |

The theorem is conservative: every refusal in the table is a numeric pass (identity matchings between flag-type
words of equal support were never non-polynomial in this sample). The negative control is the spot check.

**(ii) Spot check** (`spot_check`, `spot_check_assoc`): without GrTable multiplication tables (`LightLoc`: only fibre
bases; test `lightloc_matches_sbim` against `loc_matrix`): (a) right-linearity/commutation with the localized top
action at random cosets; (b) polynomiality of Φ = L_tgt⁻¹ f L_src on sampled fibre columns along a random line
(degree bound D/2, D/2 + 2 points), negative-degree entries zero. Results: accepts the theorem Assocs at (1,1,1),
(2,1,2), (2,2,2) and certified Pull templates; rejects ≥ half of the single-coset corruptions (scale one coset by 2) at
each size (test `spot_check_accepts_and_rejects`). **(3,3,3) braid middle (1 680 flows): theorem ok (multistep
[[3,6],[],[3,6]]), spot check passed: 16 cosets, 6 columns, 7 points, 10 080 entries, max degree 10, 17.0 s.**

**(iii)** Nothing downstream was changed: certification of pairs is still `certify_pairs` and the B4 results in §5
used no theorem node.

**(iv) Marking.** `relies_on_theorem(ar, m)` lists the theorem-certified Assoc nodes of a DAG; `smapio::export`
writes `thm-assoc <node>` before `end` for each; `import` re-runs the theorem check on exactly those nodes (a record on
a non-Assoc node is an error; without the record a fresh arena reports the map not genuine). Test
`marking_and_smapio_roundtrip` (export, import in a fresh arena, genuine, same localization). `braid_dag_thm` is the
B1 braid with the theorem middle step.

## 3. Singlestep forms and the rex-chain planner (src/rexplan.rs)

- `singlestep` writes a word as EKLP's ±s sequence (a rung = walls added by the giver's split, then removed by the
  merge; trivial steps dropped); `max_table`, `window_of`, `splice`, `monos_of_deg`, `predict_template_secs`
  (fit to the B1 logs: 5·10⁻⁷·M·n² + 7·10⁻⁹·M³ s, ± factor 2).
- Moves (all degree 0, 1⊗ ↦ 1⊗, so every composite between reduced words of q has top scalar 1, SLL Thm 6.11, Prop
  4.19): **Collapse** (window word reduced for the window maximum → `canonical_max`, Rex ψ template), **Expand**
  (reverse, Pull φ), **Assoc** (flag-type, equal multistep; numeric ≤ 500 flows, else theorem if up-down, else
  uncertifiable), **FlagMatch** (flag-type, equal support and shift, different multistep: Pull identity matching,
  template), **Swap** (disjoint uprights), **Exact** (two reduced expressions of the same window coset, degree-0
  `hom_space` window ≤ budget; the basis element with 1⊗ ↦ c·1⊗, c ≠ 0, scaled by 1/c).
- Search: bidirectional BFS over words reduced for q (every visited word checked with `reduced_for`, gate m9); closed
  moves first, then with Exact moves, then (if still nothing) on λ·0 (4 uprights, would mean `Pad`).
- `leafgen::realize_chain` turns a chain into a DAG (Local templates certified by `build_in`, Elem Assoc via
  `validate_assoc` or `add_thm_assoc`, Exact windows via `try_hom_basis_cap`).
- Why FlagMatch and Exact were needed (found by the small-λ gate): reduced words differing by an EKLP far
  up-down commutation (e.g. (1,2,1): E1.1 F1.1 E0.1 vs E1.1 E0.1 F1.1, singlesteps [+3 −3 −2 +1] vs [+3 −2 +1 −3])
  and by a switchback whose route passes a 4-block level (e.g. (2,2) window F1.1 E1.1 vs E1.1 F1.1 inside a 3-upright
  word whose third upright is occupied). The first is a FlagMatch (certified template), the second a small Exact
  window; neither needs a 4th upright.

**Tests (tests/rexplan.rs):** `singlestep_forms`, `max_table_and_canonical_max`, `window_splice_roundtrip`,
`planner_connects_small` (λ = (1,1,1), (1,2,1), (2,1,2), (2,2,1), (2,2,2): from L_q to every reduced word of ≤ 4
rungs of every type: 33 + 123 + 161 + 166 + 93 = 576 searches incl. identities, 0 failures), `realized_chains_small`
(356 realized chains: genuine, degree 0, 1⊗ ↦ 1⊗ at every coset, scalar 1 at rep(q), in the exact Hom⁰ window when
≤ 2·10⁴ unknowns), `rex_class_pairs_small` (all 53 rex-class pair targets of the (2,1,2) and (2,2,2) G1 plans: certified, 64 s).

## 4. Census G-P0 (src/bin/gp0.rs)

Per pair target (p, ρ, q, s, mult): X = L_p·ρ; rex class (X reduced for q); Lemma A plan (`Leaves::plan_words` with
target q, A2's cost descent); step unknowns (`step_score`: max over supp G and both directions of min(window, whole));
the final word W_q (new: if the final step exceeds 10⁵, W_q is re-chosen among the reduced words of q of ≤ 7 rungs,
since the exit chain absorbs it); the chain W_q → L_q; Pad; budget verdict (exact calls ≤ 10⁵, or ≤ 4·10⁵ with an
uncensored record; template windows ≤ 1 200 flows; no uncertifiable Assoc; words ≤ 10⁶ flows). Data:
`b2/gp0_n5b.{log,jsonl}`, `b2/gp0_n6.{log,jsonl}`, tables `b2/gp0_n5_summary.md`, `b2/gp0_n6_summary.md`.

**Pad: never needed.** Every chain (all 118 + 7 targets that have a reduced W_q) was found on 3 uprights; the 4-upright
search was never entered. **B3 (MapNode::Pad) is not needed for the N = 5/6 lists.**

**Theorem-level Assoc: never needed by a chain** (max Assoc window 60 flows, numerically certifiable).

N = 5:

| pair | targets | verdict | reason / chain |
|---|---|---|---|
| 313 | 2 | **feasible** (after record) | rex-class target: chain of 2 Exact moves (≤ 64 161 unknowns); other target: W_q = L_q, one Lemma A step of 137 964 unknowns → measured (28.0 s / 32.4 s, 542 MB; record in `b2/costs_b2.txt`). Caveat below. |
| 322 | 1 | **feasible** | rex class; chain of 3 Exact moves (≤ 50 016 unknowns), words ≤ 6 480 flows |
| 326 | 2 | infeasible | Lemma A steps 956 320 and 15.7 M unknowns (chains themselves fine: ≤ 3 moves, template ≤ 20 flows) |
| 410 | 1 | infeasible | Lemma A step 1.94 M unknowns |

N = 6 (summary; full table in `b2/gp0_n6_summary.md`): **7 feasible** — 159, 323, 322, 71, 334, 419 (rex class:
chains of 3-5 moves, templates ≤ 60 flows, Exact ≤ 99 175 unknowns) and 37 (Lemma A ≤ 92 493 unknowns, chain of 6
moves); **4 measurable** (one step in 1·10⁵-4·10⁵ without a record: 212 (392 995), 335 (372 384), 411 (200 667), 313
(137 964, now measured)); **28 infeasible**, all by Lemma A lower steps: 127 (1.5 M), 227/145/383 (4.1 M …), 458
(2.5 M), 319 (701 771), 278/273 (13.1 M), 293 (972 M), 78 (21 M), 326, 285 (928 M), 385, 292 (478 M), 137 (217 M), 332
(12.7 M), 410, 82 (16 M), 283 (27 M), 418 (183 M), 336 (853 M), 381 (4.2 M), 376 (928 861), 281 (6.3 M), 219 (3.0 M),
467 (19 M), 24 (439 729), 45 (133 M). Several of these pairs have some feasible targets (127: 3 of 7; 278: 2 of 4; 273:
3 of 5; 78, 326, 285, 410: 1 each).

Predicted wall: chain templates are negligible (≤ 90-flow windows, < 1 s); chain Exact moves 0.6-10 s each (measured
in B4: 19-213 s per rex-class pair including both chains and certification); Lemma A steps: exact-shape records where
they exist, else the crude model 50·(u/10⁵)^1.3 s [E] (in the jsonl, `pred_exact_s`).

**Caveat (found by running it).** The census takes min(window, whole) per step, but `step_family` falls back to the
whole words when the window family lacks full top-to-top rank (window steps are not a priori complete, A1 §5(a)).
Running Lemma A for pair 313's non-rex target (`gp0 --b4 31 --lemma-a`) hit exactly that: a level-5 step fell back to
a whole-word call of 1 963 540 unknowns and was Refused. So "feasible/measurable" for non-rex targets is a lower
bound on cost; rex-class verdicts are exact (they were run, §5).

## 5. B4 start: rex-class pairs (`leaves::rex_class_pairs`, `gp0 --b4 31|32b`)

π = chain X → L_q, ι = chain L_q → X, realized; `certify_pairs(x, lq, 0, [ι], [π], rep(q), 1)`; then π∘ι = id_{L_q}
checked blockwise at every coset of L_q at two points. Logs `b2/b4_n5_f31.log`, `b2/b4_n5_f32b.log`,
`b2/b4_n6_f31.log`, `b2/b4_n6_f32b.log`.

| pair | X flows | q | F31 | F32b | moves (π + ι) | max Exact unknowns | wall F31 / F32b |
|---|---|---|---|---|---|---|---|
| 159 | 240 | 0,2,1/3,0,0/0,1,2 | certified | certified | 2 templates, 2 Exact | 370 | 20.5 / 18.7 s |
| 71 | 600 | 0,1,2/3,0,0/0,2,1 | certified | certified | 2 templates, 2 Exact | 4 251 | 29.6 / 27.9 s |
| 323 | 2 880 | 0,2,1/2,0,1/1,1,1 | certified | certified | 2 templates, 2 Exact | 370 | 45.8 / 46.3 s |
| 322 | 4 320 | 1,1,1/1,1,1/1,1,1 | certified | certified | 6 Exact | 50 016 | 38.8 / 34.9 s (N = 5 run: 24.1 / 22.2 s) |
| 334 | 7 200 | 0,2,1/1,0,2/2,1,0 | certified | certified | 2 templates, 4 Exact | 99 175 | 212.8 / 193.6 s |
| 419 | 7 200 | 0,1,2/2,0,1/1,2,0 | certified | certified | 2 templates, 4 Exact | 74 312 | 164.5 / 163.8 s |
| 313 (rex target 1,1,1/2,0,1/0,2,1) | 720 | | certified | certified | 1 template, 4 Exact | 64 161 | 23.0 / 20.2 s |

Form [[1]] everywhere; no map uses the theorem-level Assoc. For comparison, M5c's leaf-isomorphism windows: 735 960
(159, 323), 5.46 M (322), 6.58 M (71), 146.6 M (334, 419) unknowns.

`leaves::pairs_leaves` (Lemma A + exit/entry chains, the general pair route) is implemented and was tried on pair
313 (caveat in §4); it is not wired into a provider.

## 6. Open problems

1. **Lemma A lower steps are the real blocker** for 28 of 39 N = 6 pairs and for 326, 410 at N = 5: whole-word or
   window step calls of 0.4 M-970 M unknowns. Needed: the rung-local Grassmannian factorization of steps (D1 §1.7,
   S = id_Z ⊗ s on the K uprights) with closed-form or small exact s — i.e. B5. The census tool reports the worst step
   per target (`worst_step` in the jsonl), which is the work list.
2. **Window steps are not complete a priori**: the census's min(window, whole) is optimistic (pair 313: window 137 964
   vs the whole-word fallback 1.96 M). B5 must either prove the window families complete or plan words so that the
   whole-word fallback stays under budget; a census refinement could run the cheap window rank test.
3. **Downstream checks for the B4 maps** (approval (iii) for anything using the theorem; for these maps: count, form
   rank, genuineness and π∘ι = id are done; d∘d, K₀, χ, sindec `dechom_postcheck`, scaling, replay are not): needs
   `Leaves::pairs` wired into `provider_for_stamped_m4c` (pairs calls for rex-class targets), a new stamp, and a
   level-5/6 `libbuild run` + `replay` with the M5c layer A restored. Exact chain moves are `hom_space` templates
   (persisted in layer A by reference); FlagMatch/Collapse templates are Added templates (re-certified on import).
4. **Measurable pairs** (212, 335, 411; and the N = 6 instance of 313): measure their 1-4·10⁵ steps with
   `m4c measure` (as done for 313) — subject to caveat 2.
5. **GrTable(3, 9)/(6, 9)**: still unbuilt in this session (see §7); the persistence is in place, so a completed build
   is paid once. Not on the critical path of any census-feasible target.
6. **Pad (B3): not needed** for the current lists; revisit only if B5's local steps introduce 4-block levels.
7. The planner's Exact moves enumerate replacement words of ≤ 4 (2 uprights) / ≤ 3 (3 uprights) rungs; the census
   found every chain within these limits.

## 7. Commits and reproduction

Commits on `m4c-rexplan` (all pushed): merges; GrTable persistence + grwarm; theorem Assoc certificate + marker +
spot check + tests; planner + realize_chain + tests; gp0 census; B4 rex-class pairs + logs; census data; this report
(see `git log origin/m4c-rexplan`).

```
cargo test --release --test thmassoc                     # (i), (ii), (iv) default cases
cargo test --release --test thmassoc -- --ignored braid_middle_313_232 spot_check_333
cargo test --release --test rexplan
gp0 --n 5|6 --costs tests/data/homprov/costs_333_m5c.txt --out FILE.jsonl
python3 docs/colored_design/m4c/tools/gp0_summary.py FILE.jsonl
gp0 --b4 31|32b --n 6 --pairs 159,71,323,322,334,419
m4c measure --base tests/data/homprov/costs_333_m5c.txt --out docs/colored_design/m4c/b2/costs_b2.txt "5,4,0;F1.2,F0.2,F1.1;F1.3,F0.2;2"
```
