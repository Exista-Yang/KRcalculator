# M4b review, lens: tests, caps and gate harness

Branch colored-m4b, HEAD 8ad4ac6, diff 3c3323b..HEAD. Read-only review: no cargo, no binaries. The only
computations were Python one-liners over tests/data/homprov/*.json and *.txt, each under 1 s.
Paths are relative to /Users/existair/Documents/KRcalculator/.claude/worktrees/m4b.

## Answers to the five questions

### (1) Weakened, deleted or ignored tests and asserts

None found.

- `git diff 3c3323b..HEAD -- tests src` removes no `assert`, `#[test]` or `should_panic`. The only
  removed lines are:
  - the v1 cache-header assertion and its corruption string in tests/m4a.rs. Both are replaced
    by stricter v2 and `alg=` checks, a check that v1 headers are rejected, and a check that no
    tmp files are left (tests/m4a.rs:1353-1376).
  - `certify_elem(Assoc, ..)` in A8. It is replaced by `certify_assoc` per window, with added
    assertions that the reverse pair is not certified and that the reverse constant is 1
    (tests/m4a.rs:944-958).
- The added `#[ignore]` tests are all new: tests/homprov.rs:165 (B2 (3,3)/(2,2,2)), :252 (B3
  E1F1E1F1) and :258 (B4b). tests/homprov.rs is unchanged between the WIP commit 59eaea2 and HEAD
  (`git diff 59eaea2 HEAD -- tests/` touches only data, hom_lazy.rs and ntrunc.rs).
- The asserts added in src only tighten checks:
  - tmpl.rs:491-505 (Rickard uniqueness and nonzero);
  - smap.rs:822 (Assoc must use `certify_assoc`);
  - smap.rs:709-711 (Assoc shift equality);
  - tmpl.rs:660-662 and :702-714 (poly_certificate: tops equal, block diagonal).
- One caveat: two failing tests were made to pass by changing the code, not the tests.
  - b2_hom60: fixed by the sound certify_full orbit fallback (homprov.rs:640-672).
  - b2_type_pairs_22_212: fixed by raising the Exact unknown cap, commit dc58d0c. This changes
    what that test exercises; see finding 2.

### (2) Caps 5 000 -> 40 000 (Exact) and 1 600/2 500 -> 6 400/40 000 (Search)

The caps are used only for admission. No truncation, overflow or correctness bound depends on them.

- Exact: `max_unknowns` is read only in `admit` (homprov.rs:1009). Search: `max_window_fibre` and
  `max_window_unknowns` are read only in `win_ok` (homprov.rs:1428, 1439). All window and ladder
  calls then go through `try_hom_basis_cap(.., u128::MAX)` (homprov.rs:955), so the
  `tmpl::HOM_WINDOW_CAP` path (tmpl.rs:172, 413) is unchanged and is used only by `hom_basis`.
- Every exact call is admitted first:
  - rex moves through `iso_ok` (homprov.rs:1979);
  - commutation windows through their own cap (homprov.rs:1479);
  - leaf isomorphisms through `iso_ok` (homprov.rs:2041).
- Overflow:
  - `hom_space` indexes unknowns as u32 (sbim.rs `uidx`), and so does `kernel_lazy`. 40 000 is
    far below 2^32. Even an uncapped cost-table admission would need more than 4·10^9 unknowns
    to overflow.
  - Products of fibre counts are u128 (`flows` is checked_mul, schur.rs:1678).
    `hom_unknowns` saturates (homprov.rs:363).
- `max_iso_unknowns` (6 000) and `max_comm_unknowns` (20 000) are now below the window cap. They
  matter only when fp > 6 400. This is harmless.
- Memory: see finding 8. The final kernel in `kernel_lazy` is dense, (dim)·ncols field elements,
  and admission does not account for it. Measured RSS stays well below 1 GB:
  - costs_cloud.txt: max 117 MB among completed calls (27 494 unknowns);
  - G1(a): peak 856 MB per job (whole Arena).
- The cap raise changes which Search path is taken (`choose_window` at homprov.rs:1826-1838,
  `local_cands` at :1756). Every path ends in certification, so this changes cost, not results.

### (3) g1.rs `g1 a`

A target counts as ok only when the provider returned Ok. A killed or crashed child never counts
as ok. Closure types are included. The JSON summary matches its records.

- `supervised` drops the RESULT line whenever the watchdog killed the child (g1.rs:421). A
  child that dies without RESULT goes to the `(None, k)` arm. All of its targets are added to
  `n_t` and none to `n_ok` (g1.rs:484-489).
- `"ok":true` is printed only on `Ok` from `Chain::pairs` or `Chain::closure` (g1.rs:342-346,
  369-376). Every provider's Ok path goes through `certify_pairs` or `certify_closure`:
  - Exact: homprov.rs:1077, 1122;
  - Search: homprov.rs:2284, 2413, 2501.
- Escaped strings cannot produce a spurious `"ok":true` substring (jstr escapes `"`).
- Closure jobs exist for every type of (λ, λ) in `g1 a` (g1.rs:452, 467).
- I recomputed tests/data/homprov/g1a_cloud.json from its records:

  | | (3,3) | (2,2,2) |
  |---|---|---|
  | pair jobs | 12 | 98 |
  | targets | 72 | 358 |
  | Σ mult | 121 | 404 |
  | closure jobs | 4 | 21 |
  | closure targets | 8 | 42 |
  | Σ n | 16 | 96 |

  - Totals: 480 targets, all ok; 135 jobs, all fully ok; 0 killed. This equals the header and
    SPEC §4.7.
  - got == mult and got == n in every record.
  - The 110 (λ, p, rung) pairs and the 25 closure types are all distinct, and each job label
    matches the index the child reports. So no pair was run twice or skipped.
  - Max job wall 76 s; sampled peak RSS 856 MB.
- Caveats: findings 1, 5 and 7. They concern what "certified" checks for pairs and how the run
  is recorded, not the counting.

### (4) Cost tables

A stale or censored record cannot cause a false certification. The cost table is read only by
admission (`Exact::admit` homprov.rs:1005-1017, `Search::win_ok` :1416-1440, `iso_ok` :1444-1453).
Every admitted call is an exact `hom_space` followed by the usual certification.

A censored record can, however, cause a false **admission**: an uninterruptible call that the
record itself shows to be slow. See finding 3. Tables carry no solver or machine stamp, so stale
tables load silently. No test covers `CostTable` or `admit`.

### (5) smap.rs / tmpl.rs (M4a review minors)

I found no unintended change in M4a test behavior. The intended changes are below.

- `points` reseeding (smap.rs:187-205) changes every random point in every suite. No golden data
  depends on point values; tests/data has only homprov, schur and wang_ref.
- `add_template` now always records `Added { certified: false }` (tmpl.rs:351-357). The only
  non-test callers are Search's window pieces. They call `certify_template` right away
  (homprov.rs:1731-1739). `dot_pmat` (homprov.rs:2220) only builds a PMat.
- `certify_elem(Assoc)` now panics (smap.rs:822). No other caller passes Assoc:
  tests/m4a.rs:789 and :851 pass Eps/Eta/EpsL/EtaL and Iota/Trace. src has no Assoc user outside
  smap.rs.
- `poly_certificate` now rejects different tops and off-block entries. Genuine bimodule maps are
  block diagonal identically, so this cannot cause a false rejection. b5 still fails on the
  rescaled copy.
- `refresh_certified` is O(#nodes) per certification (smap.rs:807-819, called at :824, 833,
  841). This is a scaling note for M5/M6a Arenas, not a defect at G1 sizes.
- The new hard assert in `tmpl::rickard` (tmpl.rs:491-505) runs by default only for
  min(c, d) <= 2 (tests/m4a.rs review_c2). For (3,3), which the 3-colored target needs, it runs
  only in the ignored a9_rickard_33. I believe it holds (Rickard rigidity), but it is unexercised.

## Findings (most severe first)

### 1. [major] certify_pairs (and the G1 "ok") never checks the ends or the degree of ι and π

- `certify_pairs` (homprov.rs:750-770) checks:
  - the form rank at coset c;
  - the same form at a second point;
  - genuineness.

  It never compares `ar.info(ι)` with (L_q -> X, degree s) or `ar.info(π)` with (X -> L_q, -s).
  `certify_full` does check this (homprov.rs:622-627).
- `certify_closure` (homprov.rs:772-797) checks the degree multiset but not the ends.
- Failure scenario:
  - A Search leaf's bookkeeping shift `lf.shift` (used to select leaves, homprov.rs:2080) differs
    from the true degree δ of its maps. Its ι and π then come from the same leaf, with degrees
    s+δ and −s−δ.
  - π∘ι has degree 0, so its top-coset entry is a constant. It is invertible and equal at both
    points.
  - Result: the pair is certified at the wrong shift s. `g1 one` prints `"ok":true`, and M5/xengine
    would place the summand at the wrong grading.
- Only B6 (tests/homprov.rs:325-363, span of the exact Hom^s for (2,2), (2,1), (1,2)) would catch
  this. G1(a) and the (3,3)/(2,2,2) targets do not.
- I found no live instance: degrees compose in the Arena and `certify_added` checks template
  degrees (tmpl.rs:376-386). But certification is meant to be the single trust boundary, and for
  pairs it omits exactly the label downstream uses.
- Fix: pass (x, lq, sdeg) to `certify_pairs` and check `info` for every ι and π. Check the ends
  in `certify_closure`.

### 2. [major] B2 "through the Chain" no longer exercises Search `full`; B2 (3,3)/(2,2,2) not run

- Before dc58d0c, `b2_type_pairs_22_212_through_chain` failed 45 of 605: Exact was refused by the
  5 000 cap and Search could not complete them (commit message; HANDOFF "Tests").
- The fix raised the Exact cap to 40 000. The commit states that every (2,1,2)/(2,2,2) Hom^d,
  d = 0..4, between canonical ladders with fp <= 1 600 is now admitted. So the 45 cases, and most
  likely all of (2,2)/(2,1,2), are certified by Exact. The test can no longer tell:
  `chain_type_pairs` discards SearchStats (tests/homprov.rs:141-142), and SPEC §4.6 asks for
  "counts per composition".
- The B2 instance for the gate λs, (3,3) and (2,2,2) (tests/homprov.rs:165-176), is `#[ignore]`.
  HANDOFF lists "homprov 10 passed / 3 ignored" with no report of running it.
- For (3,3) the Chain has Exact only through the cost table (homprov.rs:2816-2822). So Search
  `full`/`end0` is the main path there.
- M5's first call is `provider.end0(L_p)` (SPEC §5.2 item 1). G1(a) exercises only pairs and
  closure.
- Result: Search `full`/`end0` completeness on (3,3), the 3-colored case, is unmeasured. On
  (2,1,2) it was known to be incomplete (45/605) before the Search cap raise, and is now hidden.
- No wrong result follows from this, but the B2 test does not test what its name and SPEC claim.
- Fix:
  - record and assert per-provider counts in `chain_type_pairs`;
  - run B2 with `Search::tiers(0, 6)` alone for (2,2)/(2,1,2);
  - run and report the ignored (3,3)/(2,2,2) B2 before M5.

### 3. [minor] Censored cost records admit calls (inverted meaning); tables are unversioned and untested

- In `Search::win_ok`, an exact-shape record decides: censored counts as infinite (homprov.rs:876-891,
  1425-1427). But the **bucket** fallback (homprov.rs:1431-1437) uses `CostTable::predict`.
  `predict` is the bucket maximum, and censored records are included at their kill time
  (homprov.rs:851-857).
- `g1 calib` kills at `--limit` 60 s by default and records `secs = w`, about 60 s, as censored
  (g1.rs:653, 670). That is below `max_window_s = 120` (homprov.rs:1399).
- So a bucket whose maximum is a censored "slower than 60 s" record admits **every** window in
  that bucket (sorted labels, fp within 2x, same degree), whatever its unknown count. The
  `max_window_unknowns` cap is bypassed. M4b_design.md:129 says the bucket may admit only at
  <= 30 s.
- The same holds in `Exact::admit` (homprov.rs:1011) for a censored record with secs <= 10 s.
  This happens on a watchdog memory kill (calib uses RSS 1 024 MB, g1.rs:653) or with
  `--limit` <= 10.
- Measured with a scan of the shipped tables:
  - costs_m4b.txt (the laptop table in the HANDOFF command): 16 buckets with fp <= 6 400 whose
    maximum is a censored 60 s record. Examples: ((2,2,2), fp 1 025-2 048, d 0), ((3,3), fp
    1 025-2 048, d 0), ((0,1,5), fp 4 097-8 192, d ±1). With that table Search admits any
    window in these buckets.
  - costs_cloud.txt: its only censored record has fp 10 368 > 6 400, so it has no effect on the
    gate as run.
- The header `krhom-costs v1` (homprov.rs:901) has no solver or machine stamp, so an obsolete
  table loads silently (g1.rs `exact_from`).
- No test covers `predict`, `predict_exact` or `admit` (only Send/Sync, tests/homprov.rs:506).
- Effect: only uninterruptible long or large `hom_space` calls (SPEC §4.5: admission is the only
  in-process guard). Never a false result.
- Fix:
  - skip censored records in `predict`, or treat them as infinite there;
  - add a solver/host stamp to the header;
  - unit-test the admission rules.

### 4. [minor] A stale cache file makes Search admit an uncapped call and rank it as free

- `TmplLib::has_hom` (tmpl.rs:404-406) tests only that the cache file exists. The path hash
  (tmpl.rs:304-307) does not include the format version or TMPL_ALG.
- A rejected file therefore still counts as "cached":
  - old v1 header files from dev/M4a runs, at the same path in the default
    `temp_dir()/krhom_tmpl_cache` (g1.rs `default_cache`);
  - any other stamp mismatch.
- `Search::win_ok` then admits (homprov.rs:1417-1419), and `win_cost` returns 0 (homprov.rs:1455),
  so `choose_window` prefers that window.
- `exact_call` then recomputes with `u128::MAX` (homprov.rs:955), because `load` rejects the file
  (tmpl.rs:321-329).
- Effect: one unguarded call per stale file (cost only).
- Fix: `has_hom` should validate the header, or the path should include the stamp.

### 5. [minor] Exact-call counters miss Search's hom_space calls; Chain drops failed providers' stats

- `Search::exact` ignores its `_st` (homprov.rs:1461-1466), as do the commutation, rex and
  leaf-iso calls (homprov.rs:1486, 1932-1933, 2049-2050). They count into `self.window_calls`,
  never into `SearchStats::exact_calls` or `exact_secs`.
- `chain_call!` never adds a failed provider's stats to `acc` (homprov.rs:2826-2849,
  `let _ = &mut acc`), although the doc at :2804 says they are accumulated.
- So the G1 JSON per-target fields are wrong:
  - `exact_calls` and `exact_secs` are 0 for Search-certified targets even when Search ran
    hom_space windows;
  - the time and candidates spent in providers that failed earlier in the chain are missing;
  - the `"tiers":"Exact"` label (g1.rs `tiers_str`) is right only because every Search Ok path
    bumps a tier.
- B8's `assert_eq!(st.exact_calls, 0)` (tests/homprov.rs:457) cannot detect hom_space calls by
  Search. B8 still holds through its tier_hits assertion (:458).
- The summary counts are unaffected.

### 6. [minor] B4b's "forced witness set fails the counting assertion" is not tested

- tests/homprov.rs:281-285 never calls `certify_full` with a forced witness set.
- Its assertion `plan.rows_per_point * (plan.m + 2) < 69 || ...` is automatically true when
  `plan.kind == "more-points"`. In that case `rows_per_point = E_w`, and the non-witness choice
  already means (M+2)·E_w < dim (homprov.rs:592, 604-605).
- The second disjunct counts cosets, not entries.
- SPEC §4.6 B4b requires showing that the counting assertion, not the rank test, fails. The test
  is also `#[ignore]` and was not reported as run.

### 7. [minor] Gate runs record neither their parameters nor truncation; the G1(b) rules are not implemented

- `write_json` always ends with `"complete":true` (g1.rs:501, 506-517). Neither the JSON nor the
  summary line records `--only`, `--max-pairs`, `--wall`, `--rss`, `--costs` or the commit
  (g1.rs:424-499). So a truncated or relaxed run looks like the full gate.
- The cloud G1(a) ran with `--wall 900 --rss 4096`, not SPEC's 600 s / 1 GB (HANDOFF). The
  recorded max wall (76 s) and sampled peak (856 MB, polled every 250 ms, g1.rs:398-418) imply a
  pass at SPEC limits, but this is inferred, not measured.
- For the queued G1(b):
  - the SPEC §4.7 rerun rule (Exhausted -> l_max+2, 4·c_max, report both) is not implemented;
  - the summary counts closure jobs as "jobs", while the criterion is >= 95 % of sampled pairs
    and no Exhausted pair.

  These are conservative (a false fail), not false passes.

### 8. [minor] kernel_lazy builds the final kernel densely; admission does not account for it

- `kernel_dense` (linsys.rs:256-279), used for the result at :328, allocates (ncols − rank) dense
  vectors of length ncols. The old `Echelon::kernel` was sparse (linsys.rs:84-101).
- At the new 40 000-unknown cap one vector is 320 KB (F31 is a u64): 46 MB for dim 145 (hom60
  case), 640 MB for dim 2 000.
- Admission counts unknowns, not dim Hom^d. The early stop is off when dim > 128
  (linsys.rs:305), so such calls also scan every row.
- Realistic dims at d <= 4 are at most a few hundred, so the risk is low.

### 9. [minor] hom_lazy does not show that the early-stop and reinsertion paths run; it is the only guard on hom_space

- hom_space templates are genuine by construction (tmpl.rs:361-364). certify_full checks only
  count = formula and stacked rank (homprov.rs:621-678). A wrong but independent basis of the
  right size, made of non-bimodule maps, would be certified. So `kernel_lazy` must be exactly
  right.
- I checked the algorithm by reading it and found it correct:
  - elimination with the dense accumulator and `hi` extension (linsys.rs:205-253);
  - back-substitution with zero at pivot columns > f (:256-279);
  - the early stop as exact verification (:305-325).
- tests/hom_lazy.rs compares with the reference on 10 windows. But it does not show that any case
  takes the early stop or the violated-row reinsertion branch (linsys.rs:309-321).
- Fix: add a counter or debug hook and assert that both branches are taken.

## Checked, no defect

- Exact refuses before any call (homprov.rs:1019-1031). B8's (3,2) 144 -> 8 refusal comes from
  the label rule (max label 3 > 2), so the cap raise does not change it.
- Exact and Search admitted calls go through hom_space with full certification (homprov.rs:1049,
  2769).
- The Chain order follows SPEC (homprov.rs:2810-2823).
- The tmp cache names are per process and per call, and failed renames are cleaned up
  (tmpl.rs:203-229). A13 checks for leftovers.
- calib skips cached loads as non-measurements (g1.rs:661-664).
