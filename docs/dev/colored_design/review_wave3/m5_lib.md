# Review wave 3: M5b/M5c library building and layer A (revm5, colored-m5 f2fcb1d vs dev 370f5b1)

Scope: src/bin/libbuild.rs, the src/sindec.rs changes (route stamp, knot_jobs, target log), and the parts of
homprov/smapio/tmpl that layer-A load and replay depend on. Read-only; nothing was run.

## Verdict

- **Can a wrong Poincaré polynomial or a falsely certified map come out of this?** No. Every map loaded from layer A
  goes through the same certificate as a fresh one, and the M6a scan uses a provider that refuses every call.
- **The defects are in libbuild's bookkeeping:**
  - **Resume ignores the output file's header.** Its run/replay summaries can report jobs as certified (including
    `"complete": true`) that were never run in that mode, with that field or N, or against that cache.
  - **Replay can still run the exact solver.** Templates the layer-A files reference are recomputed if they are not
    already in the template cache, so replay does not prove the library is self-contained.

## Findings

### F1 (major): resume ignores the output file's header; a replay into the default output file reports the run's results as its own
- **Default output name has no mode in it.** `libbuild.rs:517` builds
  `libbuild_{lam}_{set}_f{field}.jsonl` for both `run` and `replay`.
- **Resume reads every attempt line marked ok, with no header check.** `libbuild.rs:518-534` puts every
  `{"attempt"...,"outcome":"ok"}` line of OUT into `done`. It never looks at the preceding header's mode, lam, n,
  field, `--cache`, stamps or commit. The comment "(same mode)" at line 518 is not implemented.
- **Concrete failure.**
  1. `libbuild run --lam 3,3,3 --n 6 --field 31 --cache C`, then `libbuild replay` with the same flags and no
     `--json`.
  2. The replay reads the run's file and skips every job the run certified (330). It runs only the failing jobs.
  3. Its summary says mode "replay", certified 330, `complete` true if the run was complete
     (`libbuild.rs:620`, `:638`).
  4. So "replay = run" always holds, and nothing was loaded through NoProvider.
- **Same problem with an explicit `--json` reused across settings:**
  - Another field: F32b skips F31's jobs, while the layer-A keys are per field.
  - Another `--cache`: jobs are skipped although their certificates are not in this cache.
  - Another N: `pruned_plan` pairs have the same indices at every N (`sindec.rs:1404-1430`; only kept/killed
    changes). A pair certified for its N=5 targets is skipped at N=6, whose extra targets were never built.
  - After a stamp/version change: the old layer-A files are rejected by the new code, but the run file still counts
    those jobs as done.
- **`cover` has the same blind spot.** `run_status` (`libbuild.rs:856-869`) takes any ok attempt line from the given
  files without looking at the header's n, field or mode.
- **Tests do not catch it.** tests/libbuild.rs uses separate run.jsonl and replay.jsonl files.
- **The reports cannot tell.** M5b and M5c state "replay = run" and "replay certifies exactly the same jobs", and
  cannot show which way the replays were run. One hint they really ran: M5b had to replay pairs 135 and 382 at
  `--rss 4096`.
- **Fix:**
  - Put the mode (and seed, see F2) in the default name.
  - On resume, accept an attempt line only if its governing header matches: mode, lam, n, field, set, cache, stamps.
  - Never let a replay resume from run lines.
  - Make `cover` check the header's n and field.

### F2 (minor; it does set a false `"complete": true`): the certified count includes jobs outside the current set, and the default name ignores `--seed`
- **Where.** `libbuild.rs:620`: `n_ok = done.len()` counts every ok key in OUT, not only those in `set.jobs`.
  `libbuild.rs:638`: `complete = !restricted && n_ok == set.jobs.len()`.
- **Why sample runs collide.** `set_name` (`libbuild.rs:127-135`) leaves out `--seed`, so `--sample random` runs with
  different seeds share one default OUT.
- **Concrete failure.**
  1. `run --sample random --count 20 --seed 1`: 20 jobs certified.
  2. `run --sample random --count 20 --seed 2`: 20 different jobs, and all of them fail.
  3. Result: `n_ok = 20 + 0 = 20 == 20`. It prints "certified 20/20 (100.0%)" and the summary has
     `"complete": true`, next to 20 FAIL lines.
  4. In general the count can also go above 100%.
- **Fix:** count `|done ∩ keys(set.jobs)|`, and add the seed to `set_name`.

### F3 (minor; the replay check proves less than it claims): replay can still run the exact solver, and the test's "no new files" check is weak
- **What replay blocks, and what it doesn't.** In replay only the HomProvider is NoProvider (`libbuild.rs:262-283`,
  `:327`). Loading a layer-A file runs `smapio::import`. For every `H` template line it calls
  `ar.lib().try_hom_basis_cap(&s, &t, deg, u128::MAX)` (`smapio.rs:326`). If the template is not in the tmpl cache,
  `TmplLib` solves `hom_space` with no size limit and no admission check, then writes the result to the shared cache
  (`tmpl.rs:479-490`).
- **Consequences.**
  - A replay can pass while running exact solves, and it changes the cache it is supposed to be checking.
  - So "everything loads from layer A" holds for provider calls, not for templates.
  - The content-hash check at `smapio.rs:327-329` keeps this sound. It is not a proof that the library is
    self-contained, and the solve is bounded only by the watchdog.
- **When a template can be missing.** `layer_a_save` renames the layer-A file into place before
  `lib_ref().save()` (`sindec.rs:1186-1190`). Hom calls are stored as soon as they are computed, so the gap matters
  mostly for a cache restored or snapshotted out of order: `persist_layer_a.sh` copies layer_a and tmpl one after
  the other while a run may still be writing.
- **Test weaknesses.**
  - tests/libbuild.rs:53 uses `!rt.contains("\"layer_a_saved\":1")`. That only rejects counts whose decimal starts
    with 1; 2 to 9 pass.
  - Nothing checks that the tmpl directory is unchanged or that the template library's computed count stays 0
    during replay.
- **Fix:** give replay a tmpl library that only loads and never computes, and report and check `computed == 0`
  (or the tmpl file count) per replay child.

### F4 (minor; wrong record of which code ran): after a rebuild, children run the new binary but the header names the old commit
- **Where.** `supervised` (`libbuild.rs:395-399`) starts `current_exe()` with " (deleted)" stripped. After a rebuild
  during a run (this happened in M5b, per M5b_report:164), every later child runs the new binary.
- **What goes wrong.**
  - The OUT header's `commit` and `stamps` are computed once in the parent (`:545-549`).
  - The child's RESULT line carries neither.
  - So certified results are recorded under the wrong code version. Layer-A headers carry the child's real stamps;
    the run file does not.
- **Fix:** have each child print its commit and stamps in RESULT, and have the parent reject a mismatch.

### F5 (minor; a test that may fail or pass depending on the host and on test order): `sindec::m5c_layer_a_route_stamp`
- **Why the route check fails off the cloud host.**
  - `builder()` and the test's own providers use tests/data/homprov/costs_cloud.txt, stamped for the
    2.10 GHz cloud host.
  - On any other host, `CostTable::load` rejects the table and the stamp becomes `costs=rejected`
    (`sindec.rs:455-458`).
  - The route check then compares `costs=hs-kl1` (the rewritten files) with `costs=rejected` (the current header),
    rejects every file, and `assert_eq!(sc.provider_calls, 0)` fails.
- **Why it depends on test order.** The test passes off-cloud only if `m5c_leaf_iso_refusal_is_engineering` has
  already run `std::env::set_var("KRHOM_COSTS_ANY_HOST", "1")` (tests/sindec.rs:790). That call sets the variable
  for the whole test process.
- **Fix:** set the env var (or load a table stamped for the current host) inside the route test itself.

## Lens questions: what was checked

**Can libbuild count an uncertified job as certified?**
- **Child crash or partial output: no.**
  - `supervised` returns a RESULT line only if the child was not killed (`:442-443`).
  - Without a RESULT line the outcome is "crash" (`:481`).
  - RESULT is the child's last action, printed after the save (`:384-388`).
  - The library has no Drop impls and does not print to stdout.
  - `ok` is set only when `errs` is empty (`:387`).
  - Every error path in pair, cross and closure jobs pushes to `errs`.
- **Count mismatch: no.** A wrong count turns into outcome "Count" (`:587-597`), which is not retried and is counted
  as a failure.
  - `targets` and `sum_mult` in the child count only the kept (q, s) targets of the one `dec_one`, and only after
    each has passed the form, genuineness and rank checks (`sindec.rs:689-707`).
  - `pruned_plan`'s kept list has the same (q, s) keys, so the counts match.
- **Resume from a stale or foreign OUT: yes.** See F1 and F2.

**Does replay load everything from layer A?**
- **Provider level: yes.** NoProvider refuses all four calls. The resulting `MissKind::Refused` maps to
  `DecErr::Provider{engineering: true}` and the job fails. There is no fallback path in `target_pairs`
  (`sindec.rs:715-748`) or in `closure_gens` (`:1041-1087`).
  - The identity and all-killed shortcuts need no maps.
  - The Rickard components and orthogonalization only evaluate.
- **Template level: no.** See F3.

**Is the route stamp sound?**
- **Pairs: verified that a loaded map is certified the same way as a fresh one.**
  - Loaded maps (`sindec.rs:717-726`): `dechom_postcheck`, which checks the count (2·mult), ends, degrees and D14
    scaling.
  - Fresh maps: `DecHom::pairs` (`:381-392`) runs `dechom_precheck` (target only), the provider, and
    `dechom_postcheck` at a second seed. Then `target_pairs` runs `dechom_postcheck` again (`:731`).
  - Both then go through `dec_one`: the form at rep(q) at 2 points with an inverse, `assert_genuine` on every map, the
    rank-bookkeeping assert, `orthogonalize` and `check_pi_iota` (`:666-712`).
  - What a loaded map skips:
    - `dechom_precheck`: it checks the target, not the maps, and the target comes from the same `split_mult`.
    - The duplicate postcheck at the second seed.
  - Neither of these is needed to certify the decomposition. The deciding check is an invertible form into
    End⁰(L_q) plus genuineness.
- **Closures: verified.** Loaded and fresh maps both go through `certify_closure` (`:1083`): ends, the c_diag degree
  multiset, identity-coset rank at 2 points, genuineness.
- **Template content:** `smapio::import` checks it by content hash.
- **Conclusion: soundness does not depend on the stamp.** The cost table only decides which calls are tried.
  `route_header` rewrites only the `costs=` token. Keys contain no spaces, so the token split is safe.
  `costs=none` and `costs=rejected` never match a real table.

**Does the K_0 lower bound in `knot_jobs` hold?**
- **As a lower bound: yes.**
  - A type whose signed, q-graded class is nonzero is in every complex homotopy equivalent to the stage
    (K_0(K^b(A)) ≅ K_0^⊕(A), and K_0 is free on the non-killed types).
  - So the scan has to decompose that type through every kept Rickard term, rung by rung, and close the final ones.
- **The class computation is correct.**
  - `add_scaled(cp·m, ±1, dq)` with the hdeg sign and dq shift is right.
  - Dropping killed types at every step agrees with dropping them at the end, because the killed types span a
    two-sided ideal.
  - The term-kill rule `max(c,d)+k > N` is the same as `swire::kept_terms` and is asserted against the cell rule in
    `rickard_comp_at`.
- **The upper bound is what the wired engine actually needs.** `swire::build_wlib` (revm6a src/swire.rs:106-160)
  calls `dec_word` on every upper-bound type with no cancellation and fails on any error. It also tries
  `closure_gens` for every type at every stage, but those failures are not fatal.
- **The coverage doc's "full (LB)" verdict is stated as conditional,** and the "exactly these jobs" claim for N = 5
  rests on the upper bound being contained in the level set, which `tests/knotcov.rs` asserts.
- **Gaps (not code defects):**
  - The test checks only the unpruned closure lower bound against `decat_scan`; `pairs_lb` and the pruned case are
    not checked independently.
  - The pruned lower bound assumes the non-killed canonical types stay indecomposable and pairwise non-isomorphic at
    level N, which is a design premise.
