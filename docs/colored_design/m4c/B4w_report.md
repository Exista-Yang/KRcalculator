# M4c B4w report: Leaves pairs and closure in the production provider; level 5/6 rebuild

Date 2026-11-13. Branch `m4c-rexplan` (4e886f6 + merge origin/dev f4c865a). Data: `data-gl6/m4c_b4w`. Machine:
4 cores (Xeon @ 2.10 GHz), 15 GiB; every run with `KRHOM_COSTS_ANY_HOST=1`. Frozen files (smap.rs, homprov.rs,
sbim.rs, linsys.rs, tmpl.rs, sindec.rs) and the colored-m6a files (sengine.rs, swire.rs, boxclose.rs, boxrun.rs,
sbox.rs, libbuild.rs, smapio.rs) were not edited. Changed: src/leafgen.rs, src/leaves.rs. New: src/bin/b4wstamp.rs,
tests/b4w.rs, docs/colored_design/m4c/b4w/.

## 0. Summary

| set | prime | jobs | M5c certified | **B4w certified** | B4w replay | new |
|---|---|---|---|---|---|---|
| level 5 | F31 | 177 | 172 | **174** | 174 | pair 322, closure 1,1,1/1,1,1/1,1,1 |
| level 5 | F32b | 177 | 172 | **174** | 174 | same |
| level 6 | F31 | 377 | 330 | **345** | 345 | pairs 71, 159, 322, 323, 334, 419; the 9 closure types |
| level 6 | F32b | 377 | 329 | **344** | 344 | same |

Nothing certified by M5c was lost. Replay (a provider that refuses every call, so every map is loaded from layer A
and re-certified) certifies exactly the same jobs for every level and prime. **All 12 M5c closure failures and all
6 rex-class (leaf-isomorphism) pair failures are now in the production library.** The 32 (F31) / 33 (F32b) remaining
level-6 failures and the 3 remaining level-5 failures are all non-rex pairs. Each needs a Lemma A step above the
4·10⁵ hard cap, which is the B5 work list (§4).

- **leafgen memory gate** (review wave 4, finding 1): fixed and tested (§1).
- **Measured pairs:** none of 212, 335, 411 or 313 (N = 6) can be certified. Their failing calls are whole-word
  fallbacks of 1.96 M to 2.4 G unknowns, so no cost record admits them (B2 caveat 2). The 1·10⁵-4·10⁵ band calls of
  127, 137, 227, 385 and 467 were measured (9 records). Each of those pairs then stops at a step above the hard
  cap (§3).
- **Coverage:** still 0 of 185 three-strand knots at N = 5 and at N = 6. N = 5 now needs only pairs 313, 326 and 410
  (M5c: 5 jobs). At N = 6, 25-33 blocking jobs per knot (M5c: 38-48), 25 of them universal. 4_1 has 258/283 jobs
  certified (M5c: 243).
- **Full suite** `KRHOM_THREADS=4 cargo test --release`: **278 passed, 0 failed.**

## 1. leafgen: interpolation / certificate memory gate (review wave 4, finding 1)

`leafgen::interp_gate(src, tgt, deg)` checks three limits from the words alone (no SBim, no GrTable) and refuses
with Err if any is exceeded:
- **Flows:** the flow counts (`Expr::nflows`) must be ≤ `MAX_INTERP_FLOWS` (1 200). This check now runs before any
  SBim is built.
- **Monomials per entry degree:** at most `MAX_INTERP_MONOS` (6 000) in any one entry degree.
  - Fibre degrees of a window word lie in [0, −2·shift].
  - So entry degrees lie in [u − maxf_tgt, u + maxf_src], with u = deg − shift_tgt + shift_src.
  - Counts use `rexplan::monos_of_deg` on the bottom labels.
- **Stored localization entries:** (monomials + 4) · ns · nt must be ≤ `MAX_INTERP_CELLS` (2.5·10⁸). The + 4 is the
  point count of the certificate.

Where the gate runs:
- `interpolate` (and so `build` and `build_in`);
- `validate_assoc` (before `poly_certificate`);
- `theorem_vs_numeric` (skips the numeric certificate when over budget);
- `braid_dag`, which decides the middle step's route **before** building ψ_A / φ_B. It returns Err when neither the
  numeric Assoc certificate nor the Pull fallback fits the budget.

The certificate function itself (`tmpl::poly_certificate`, a frozen file) is unchanged. It is only reached through
these gated callers.

Test `tests/leafgen.rs::interp_gate_refuses_braid_fallback_323`:
- At cells (3,2,3), (3,3,2) and (2,3,3), the middle step (560 flows) has 39 148 monomials, matching the review's
  generating-function count.
- `interp_gate`, `interpolate` (with a closed form that panics if evaluated), `validate_assoc` and `braid_dag` all
  refuse immediately (< 1 ms).
- Small windows still pass the gate and build.

The existing leafgen / rexplan / thmassoc suites pass unchanged, so every B1/B2 template is within the gate.
Findings 2 and 3 of the review (minor) are not addressed.

## 2. Wiring (src/leaves.rs)

`provider_for_stamped_m4c(lam, costs, slowdown, leaves_first)` is what `libbuild --provider m4c | m4c-lf` uses; no
libbuild change was needed. It is now built as follows:

- **Closure:** Stage 0 product closure first, then the production chain, then Leaves (`m4c`). With `m4c-lf`, Leaves
  runs first for closure calls without e_p, as A2 recommended. This rebuild used `m4c-lf`.
- **Pairs (new):** the production chain first. If it fails, `leaf_pairs` runs (Leaves last):
  - **rex class** (X reduced for q, s = 0, mult 1): `rex_class_pairs_opts`, where π and ι are realized rex chains;
  - **otherwise:** `pairs_leaves_opts`. Lemma A leaves at q obey the Leaves budget rule (≤ 1·10⁵ unknowns;
    ≤ 4·10⁵ with an uncensored exact-shape record; the whole-word fallback included), plus rex-chain exit/entry.
  - Planner options: Exact moves ≤ 1·10⁵ unknowns, words ≤ 12 rungs, ≤ 6·10⁴ visited words.
  - Everything runs inside a mark/commit.
- **Everything else** (`full`, `end0`, closure with e_p, pairs with e_q): the chain alone.
- **Stamp:** chain stamp + `STAGE0_STAMP_PF` + `LEAVES_STAMP(order)` + `LEAF_PAIRS_STAMP` (see the MANIFEST).
- **A3 constructor:** the closure-only wiring is kept as `provider_for_stamped_m4c_a3`.

**Persistable maps only.** A numerically validated Assoc constant lives in the arena and is not part of layer A (the
A2 finding for Elem). The first sindec-level test proved this directly: 21 of the 67 (2,1,2) jobs were refused by the
new `unpersistable_elems` check. As a fix, `leafgen::realize_chain_opts(.., persistable = true)` writes numeric Assoc
moves (≤ 500 flows) as the same identity matching in the form of a Pull closed-form template. That template is
certified by `certify_added`, and smapio re-certifies Added templates on import. Theorem-level Assoc nodes (marked,
re-checked on import) are kept. `leaf_pairs` refuses any result that still contains an unpersisted Elem node.

**Acceptance is unchanged.** A Leaves pair result is accepted only through:
- `certify_pairs` (form at rep(q) at two points, genuineness), and π∘ι = id at every coset for the rex class;
- then sindec's own checks, exactly as for any provider: `dechom_postcheck`, D14 scaling, Π∘Ι = id at all cosets,
  the kept multiset, layer-A save, and the re-certified load / replay.

**Test** `tests/b4w.rs::leaf_pairs_through_sindec_212`:
- λ = (2,1,2) with the production chain refusing **every** pair call, so all targets go through `leaf_pairs`.
- Run: 67 jobs, 169 targets (rex class and Lemma A) certified through SIndec, 169 layer-A files saved.
- Replay in a fresh SIndec with a provider that refuses every call: 169 loaded, 0 rejected, 0 provider calls.

`m4c_stamp` checks the stamps.

**Carrying M5c's layer A over (src/bin/b4wstamp.rs).**
- The M4c stamp extends the chain stamp, so M5c's 1 542 / 1 540 layer-A files would be rejected and rebuilt (about
  80 min per prime).
- `b4wstamp` appends the M4c suffix to the provider token of every layer-A file whose token equals the chain stamp up
  to the cost token. Keys and bodies are unchanged; other headers are left alone. All M5c files qualified.
- The runs then load them with `--stamp-route`. This is M5c's route-stamp argument: the stamp only steers which calls
  are attempted, and every loaded map is re-certified on load. The replays confirm it.

**Cost table** `b4w/costs_b4w.txt` = costs_333_m5c (351) + the A2 and B2 records (367 records).
`b4w/costs_b4w_meas.txt` adds the 9 band records of §3; it was used only for the §3 probes, not for the library runs.

## 3. Level 5/6 rebuild, measured pairs

Procedure:
1. Restore data-gl6/m5c/{f31,f32b} (sha256 as in its MANIFEST) and run `b4wstamp --provider m4c-lf`.
2. `libbuild run --lam 3,3,3 --n 5|6 --field 31|32b --stamp-route --provider m4c-lf --jobs 4 --rss 3300` (wall
   600 s, retry at 2x), then `libbuild replay` (rss 4096).

Logs: `b4w/runs/`.

| run | wall |
|---|---|
| level 5 F31 / F32b | 48 / 54 s (replays 39 / 34 s) |
| level 6 F31 | ≈ 80 min (first pass; container restart during the last 2x retry, see below) |
| level 6 F32b | 2 659 s (3 workers) |
| level-6 replays | 1 586 s (F31, 1 worker) / 445 s (F32b) |

New maps:
- **Closure types:** certified at levels 5 and 6. Each per-side Leaves closure took ≤ 175 s, the slowest being
  0,2,1/1,0,2/2,1,0 with its 111 497-unknown recorded step.
- **Rex-class pairs:** 18-21 s each for 71 and 159 at level 6, in line with B2.
- Pair 322 and the closure 1,1,1/1,1,1/1,1,1 were computed in the level-5 run and loaded from layer A at level 6.

**Interruption (F31 level 6).** The worker restart killed the run during the 2x retry of pair 467.
- 467 had been Time at 1x, and its F32b 2x retry was Time as well.
- The F31 jsonl therefore has no summary line. The tally (345) comes from the per-job lines; the replay summary
  (345/377) is complete.

**Measured pairs (open problem 4 of B2).** Probes `gp0 --b4 31 --n 6 --lemma-a` (log
`b4w/measurable_n6_f31.log`) and the production runs. All four "measurable" pairs fail on a whole-word fallback above
the hard cap, not on their 1-4·10⁵ window step:

| pair | failing call (whole words after the window family lacked full rank) | unknowns |
|---|---|---|
| 212 | Hom^5(E0.1 F0.1 E1.3 F1.3, E0.1 F1.1 F0.1 E1.1) | 2 468 710 |
| 335 | Hom^−1 on 7-/6-rung words | 2 423 909 807 |
| 411 | Hom^0(E1.3 E0.1 F1.1 F0.1 F1.2, E1.3 E0.1 F1.3 F0.1) | 16 446 366 |
| 313 | Hom^2(E0.2 E1.1 F0.2 F1.1, E0.2 E1.3 F1.3 F0.2) | 1 963 540 |

So recording their band windows would not change anything; this is B2 caveat 2 at full strength.

Instead, the band calls that actually blocked production runs were measured (`m4c measure`; log
`b4w/measure_band_n6_f31.log`):

| call | measured |
|---|---|
| 123 240 unknowns | 48 s, 0.8 GB |
| 180 155 (both directions) | 75 / 91 s, 1.2 GB |
| 218 515 | 69 s, 1.0 GB |
| 333 995 | 211 s, 2.2 GB |
| 347 450 | 213 s, 2.3 GB |
| 396 704 | 143 s, 1.9 GB |
| 101 507 | 26 s, 0.4 GB |
| 225 943 | 73 s, 1.1 GB |

With these records, each affected pair moves on to a step above the cap:

| pair | step above the cap | unknowns |
|---|---|---|
| 127 | d = 4 | 641 683 |
| 137 / 219 | d = 4 | 615 219 |
| 227 | next step | 1 192 238 |
| 385 | other step | 31.4 M |
| 467 | d = 4 | 412 515 |

None certifies. Pair 467 misses the hard cap by 3 %. A cap of 4.5·10⁵ for recorded shapes (A2 open problem 1) is a
plan decision and was not taken.

## 4. Remaining failures by cause (all need B5 unless stated)

Level 5 (both primes): pairs **313, 326, 410**:

| pair | Lemma A step call | unknowns |
|---|---|---|
| 313 | whole-word fallback | 1.96 M |
| 326 | step call | 154.8 M |
| 410 | step call | 31.4 M |

Level 6 (F31 32 jobs; F32b 33 = the same + pair 24), with the Leaves refusal reason from the run jsonl:

| cause | pairs |
|---|---|
| Lemma A step above the hard cap, 0.4-3 M unknowns | 82 (1.70 M), 212 (2.47 M), 313 (1.96 M), 319 (2.41 M), 381 (2.97 M), 458 (2.47 M), 467 (0.41 M after band records; Time in production: the chain uses the 600 s wall first) |
| Lemma A step 3-100 M unknowns | 376 (4.05 M), 273 (15.7 M), 411 (16.4 M), 410 (31.4 M), 281 (40.2 M), 45 / 145 / 293 (44.1 M), 332 (49.5 M), 24 (F32b only; 51.7 M), 78 (98.4 M) |
| Lemma A step > 100 M unknowns | 326 (155 M), 383 (185 M), 278 (251 M), 292 (293 M), 285 (928 M), 37 (1.18 G), 418 (1.28 G), 283 (2.05 G), 335 (2.42 G), 336 (8.85 G) |
| band steps measured here, next step above the cap | 127 (642 k), 137 (615 k), 227 (1.19 M), 385 (31.4 M) |
| Time (chain uses the wall; Leaves alone stops at 615 k) | 219 (F31 Time at 1x and 2x; F32b Refused) |

**Every remaining failure is a Lemma A step** (window step or whole-word fallback) **above the budget.** B5 must cover
all of them:
- the rung-local Grassmannian factorization S = id_Z ⊗ s of Lemma A steps (D1 §1.7), with closed-form or small exact
  s, or
- a planner that keeps the whole-word fallback under the cap.

No rex-chain, Pad, theorem-Assoc, Memory or closure failure remains. Pair 24 (F32b only) is the M5c route difference.
Its Leaves route fails for both primes (51.7 M).

## 5. Coverage (`libbuild cover`, sindec::knot_jobs; b4w/cover_n5.md, cover_n6.md)

A job counts only if it is certified for both primes.

| N | knots fully covered (UB / LB) | blocking jobs per knot | needed by every knot | M5c |
|---|---|---|---|---|
| 5 | 0 / 0 | 3 | 3 of 3: pairs 313, 326, 410 | 5 (+ 322, closure 1,1,1/1,1,1/1,1,1) |
| 6 | 0 / 0 | 25-33 (UB list) | 25 of 33 | 38-48, 38 universal |

At N = 6, 4_1, 5_2 and 6_2 each have 258/283 jobs certified (M5c 4_1: 243). Coverage is still all-or-nothing in the
failure list: N = 5 for all 185 three-strand knots is exactly pairs 313, 326 and 410 away.

## 6. Open problems

1. **B5 (rung-local factorization of Lemma A steps)** is now the only blocker for every remaining job at N = 5 and
   N = 6 (§4). Smallest first:
   - 467 (412 k, just above the cap);
   - the 0.6-3 M group: 137/219, 127, 82, 313, 212/458, 319, 381;
   - then 227, 376 and the larger ones.
2. **Hard cap for recorded shapes:** 4.5·10⁵ would admit 467's d = 4 step (expected ≈ 3 min, 2.5 GB from the 397 k
   record). Pair 467 also needs its production wall: the chain's Search uses all 600 s / 1 200 s before Leaves runs.
   Options: Leaves first for pairs with a census-feasible plan, or a per-provider time split. Plan decision.
3. **Census refinement** (B2 caveat 2): every "measurable/feasible" non-rex verdict of G-P0 was optimistic here; 37,
   212, 313, 335 and 411 all fall back to whole words of ≥ 1.2 M unknowns. A cheap window rank test in gp0 would
   make the census exact.
4. Leafgen review findings 2 (Rex on non-reduced equal-shift windows panics inside winloc) and 3 (exact-window tests
   do not assert dim = 1) are still open (minor).
5. **Order dependence:** level 6 reused the level-5 layer A (same keys), so 322 and closure 1,1,1/1,1,1/1,1,1 were
   loaded, not recomputed, at level 6. This is as intended, and sound because every load is re-certified.
6. `b4wstamp` is a migration tool. A future coordinated sindec change could make the route stamp ignore a declared
   "superset" provider suffix instead.

## 7. Commits and reproduction

m4c-rexplan: merge origin/dev; leafgen gate + test; Leaves pair wiring + persistable chains + tests/b4w.rs;
b4wstamp + costs_b4w; level-5 / level-6 logs, probes, coverage; this report. data-gl6: 5e1e299 (`m4c_b4w/`).

```
cargo test --release --test leafgen interp_gate; cargo test --release --test b4w
cat data-gl6/m5c/f31/layerA_tmpl.tar.gz.part-* | tar -C C31 -xzf -
b4wstamp --cache C31 --costs docs/colored_design/m4c/b4w/costs_b4w.txt --provider m4c-lf
KRHOM_COSTS_ANY_HOST=1 libbuild run|replay --lam 3,3,3 --n 5|6 --field 31 --cache C31 \
    --costs docs/colored_design/m4c/b4w/costs_b4w.txt --stamp-route --provider m4c-lf --jobs 4 --rss 3300 --json OUT
libbuild cover --lam 3,3,3 --n 6 --knots tools/dev/ntrunc_t1/data/knots3_all.txt --runs RUN31.jsonl,RUN32b.jsonl --md OUT
```
