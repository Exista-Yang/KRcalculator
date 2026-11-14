# Spec: a one-command colored calculator (`krhom … --colored-method singular` on 3–4 strands) and a random-example checker

This spec is for one cloud development session of about 6–8 h on a 4-core, 15 GB machine. It merges four read-only reports (cli, lib, oracle, docs) on worktree `colored-2c4-m5i` @ 73cbd7c, plus my own read-only spot checks. Citations are file:line on 73cbd7c unless another ref is named.

**Two facts changed after the reports were written. The spec uses the new state.**

1. **data-gl6 is now at ee548d2.** It has a fourth (2,2,2,2) layer, `c2/l2222_n4_c27/{f31,f32b}`, holding the 6 class-A blocker jobs: 12 or 14 layer-A files and 170 or 162 tmpl files per prime. The restore order is "l2222_n4, l2222_n4_c25, l2222_n4_c24, then c27" (`git show ee548d2:c2/l2222_n4_c27/f31/MANIFEST.txt`).
2. **origin/colored-2c4-unblock is at 194a731.**
   - All 6 blockers are certified on both primes (commit message; `runs_c2/c27/lib/runjobs_6_f31_replay.log`).
   - Union job status is **1236 ok / 18 fail of 1254 jobs per prime** (`runs_c2/c27/cover/status_{f31,f32b}.txt`).
   - **18 class-A 4-strand knots now miss 0 jobs on both primes**: 6_1, 7_2, 7_4, 7_6, 7_7, 8_4, 8_6, 8_8, 8_11, 8_13, 8_14, 9_10, 9_11, 9_13, 9_17, 9_18, 9_20, 9_23.
   - 8_15, 9_22 and 11n_42 still miss 11 pairs; the Conway knot misses 16 (`runs_c2/c27/cover/cover_f31.txt`).
   - No 4-strand knot value exists yet: `runs_c2/c27/k4/box_6_1_f31.log` is empty.
   - The task text's figures ("1102/1157, a parallel session is certifying the rest") are therefore out of date.
   - Reading c27 files needs the 3-line `route_header` change from cca2db4, which filters `c27=` tokens (`git diff 73cbd7c origin/colored-2c4-unblock -- src/sindec.rs`). The m5i binary rejects those files.

---

## A. Scope, non-goals, base

**Goal (the user's own words, relayed):** "I need a calculator, not printed tables; check correctness on some random examples."

The deliverable is:
- one command that takes a braid or PD code, colors and N, and prints the Poincaré polynomial for 2, 3 or 4 strands;
- a library install tool;
- a random-example checker.

**Priorities**
- **P0 (must, about 4.5 h):**
  - base branch and reader fix (A.3);
  - library registry and fetch tool (C2, C3, C10, D);
  - `src/calc.rs` (C4);
  - krhom dispatch, output, exit codes and `--json` (C6);
  - default-suite tests (F1, F2);
  - acceptance H1–H8.
- **P1 (should, about 2 h):**
  - `tools/randcheck.py` (E) with a `--quick` run and one default run (H10–H11);
  - `--check-library`;
  - ignored library tests (F3);
  - `docs/CALCULATOR.md` (G).
- **P2 (could):**
  - checker Mode C (table knots);
  - sbox reusing the calc helpers;
  - a tiny in-repo library fixture;
  - resumable `--ckpt`.

**Non-goals (hard)**
- No table sweeps. Do not queue or slice runs over knot tables. Do not run class-A 4-strand knots such as 6_1 to completion: they are projected at 0.2–1.3 core-h reduced and 1.2–8 core-h unreduced (C2_2_report.md:189-191). Only `--check-library` may touch them.
- Do not edit:
  - `src/leaves.rs`;
  - `src/bin/libbuild.rs`;
  - the planner: `sindec::pruned_plan` (sindec.rs:1999), `sindec::knot_jobs` (sindec.rs:2109), `g1_work` (homprov.rs:37-54).
  Calling them is fine.
- **Certification stays as it is.** Every certification check runs inside `swire::box_word` (swire.rs:1375-1428) with `BoxOpts::new` production defaults: scaling on, Freivalds k=3 (swire.rs:1127-1130). The checks are:
  - the generator-determinant pass;
  - frozen-structure replay;
  - Freivalds D²;
  - χ = `decat_chi`;
  - the D14 scaling check.
  On top of those, two primes must agree.
  - Never set `BoxOpts::scaling = false` (test-only, swire.rs:1117-1120).
  - Never expose `--no-k0` or `--no-d2` (sbox.rs:150-162, 83).
  - Never relax layer-A stamp checks beyond `layer_a_route_stamp`. That flag is already sound because every loaded map is re-certified (sindec.rs:103-110).
- No new crate dependencies. Cargo.toml has only num-* and rand.
- No new library builds, no mixed-label libraries, no (3,3,3) support this session.
- Do not change:
  - sbox CLI behaviour, or the production scripts that use it (runs_c2/k3/one.sh:14, runs_c2/k4m/markov/one.sh:11);
  - the 2-strand exact path's results;
  - the default `--colored-method` choice.

**A.3 Base and push.**
- Create branch `colored-calc`. Push only to it, never to main or dev. End commit messages with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Recommended base: origin/colored-2c4-unblock @ 194a731.**
  - 73cbd7c is its ancestor (`git log origin/colored-2c4-unblock` lists 73cbd7c).
  - Its src diff against 73cbd7c is only src/bin/gp0.rs, src/leaves.rs and the 3-line `route_header` filter in src/sindec.rs.
  - It is needed to read the c27 layer.
- **Fallback, if the orchestrator insists on 73cbd7c:**
  - cherry-pick only the `route_header` hunk from cca2db4 (src/sindec.rs:455-466 region) as its own commit;
  - then run acceptance H9 (replay of the 6 c27 jobs with this binary).
- This choice is open question Q1.

**Compute:** cloud only. Use `KRHOM_THREADS=4` for single runs, and 2 per process when 2 run at once.

---

## B. CLI design

```
krhom N[,N...] (--pd PD | --braid WORD [--strands n]) [--colors c | c1,c2,...]
      [--colored-method cable|cube|singular] [--reduced] [--single-prime] [--json OUT]
  singular, 3-4 strands (localized box engine + prebuilt library):
      [--lib-dir DIR] [--seed S] [--seed-f32b S] [--seed-retries K]
      [--mem-mb M] [--mem-cap-mb M] [--no-simplify] [--check-library]
```

The calculator's command is, for example:

`krhom 2,3,4 --braid "1 -2 1 -2" --strands 3 --colors 2 --colored-method singular`

**New flags.** All are box-route only except `--json`, `--colors` broadcast and `--no-simplify`.

| Flag | Default | Meaning / source |
|---|---|---|
| `--lib-dir DIR` | `$KRHOM_LIB_DIR`, else `${XDG_CACHE_HOME:-$HOME/.cache}/krhom/lib` | Library root (D) |
| `--seed S` | 201 (all C2 production runs: runs_c2/k3/one.sh:14; runs_c2/k4m/markov/one.sh:11) | Box node seed (`BoxNodes::new(lam, n_box, seed)`, swire.rs:1378) |
| `--seed-f32b S` | same as `--seed` | Seed for the second prime. The checker uses a different one (oracle O2). |
| `--seed-retries K` | 2 | Redraw seeds S+1..S+K, **only** on "generator determinant fails … redraw the nodes" (swire.rs:1220) |
| `--mem-mb M` | 3200 at 3 strands, 3000 at 4 strands (one.sh lines above) | Per-point budget, `EngOpts.mem_point_bytes` (sbox.rs:83; refused by `SErr::Memory`, sengine.rs:950-951) |
| `--mem-cap-mb M` | `$KRHOM_MEM_CAP_MB` if set, else 4500 at 3 strands and 200000 at 4 strands (same scripts) | Closure memory cap. It is read only from the env var (cclose.rs:749-752, 2725-2728). krhom sets the var before any thread starts (precedent: libbuild.rs:781 `unsafe set_var`). The 4-strand estimate is 30–130× too high (C2_2_report.md:180-186), so the real guard is `ulimit -v` (documented). |
| `--no-simplify` | off | Keep the given word (see dispatch step 3) |
| `--check-library` | off | Load the library for both primes and run the record scan, without the box scan. Exit 0 if the word is covered, 4 otherwise. |
| `--json OUT` | none | Machine-readable results for **every** route (C6) |
| `--box-live` | hidden; tests only | Box route with the live provider and no library, ≤2 strands only (F2) |

Kept with their current meaning: `--reduced` (bottom position 0 marked, krhom.rs:26-27; box mode `Some(0)`), `--single-prime`, and the `--field q` note (krhom.rs:119-121).

`--colors`: one value per component, as today. **A single value now applies to all components**; today that case panics at the assert on krhom.rs:59.

**Dispatch.** Steps run in order.
0. **Parse errors exit 1, not panic.** Replace `assert!` at krhom.rs:59-60 and the `unwrap`s on `--braid`, `--strands` and `--colors` tokens (krhom.rs:50-56).
1. `singular` = `--colored-method singular && !--cube` (krhom.rs:84), unchanged. Steps 2–5 apply only when singular. Every other route is unchanged.
2. Colors: broadcast if a single value is given; otherwise the count must equal ncomp.
3. **Simplification** (braid input, unless `--no-simplify`):
   - Compute `s = bw.simplify()` (pd.rs:316-318: greedy pass plus deterministic random search, with Markov destabilization; component ids are carried, pd.rs:567-583).
   - **Use `s` only if** `s.strands < bw.strands`, **and** either:
     - not `--reduced`, or
     - ncomp == 1, or
     - `s.components[0] == bw.components[0]` (otherwise the marked component could change: conjugation swaps components).
   - Otherwise keep `bw`. Same-strand rewrites are not used, because library jobs are per word (`knot_jobs` walks the word, sindec.rs:2109-2175) and a rewrite can leave a covered word.
   - Print `simplified braid: <word>  (strands k)` after the `braid:` line.
   - Recompute labels from the chosen word, as at krhom.rs:61.
   - PD input is already simplified by `pd_to_braid` (pd.rs:147-149); do not simplify again.
4. **strands ≤ 2:** the existing exact path (krhom.rs:117-142), unchanged. This includes 1 strand: the colored unknot is tested at tests/xengine.rs:273-279.
5. **strands ≥ 3:**
   - a. All labels 1 → exit 1: "uncolored input: drop --colored-method singular; the Soergel engine handles any number of strands."
   - b. Labels not all equal → exit 1. Name the labels and list the supported λ. Suggest `--colored-method cable` when `colored::cable_default(labels)` holds (krhom.rs:16-18, 100).
   - c. No registry entry for λ (covers ≥5 strands, colors ≠ 2, and (3,3,3)) → exit 1, listing what is supported: (2,2,2) N≤4 and (2,2,2,2) N≤4.
   - d. Every N must satisfy `a ≤ N ≤ reg.n`, where a is the color, else exit 1. The upper bound is `valid_n = min(prune, box N)` (boxclose.rs:1171-1176, error at 1208-1211). Check before loading anything. (N < a: Q5.)
   - e. Split braid (some σ_i with 1 ≤ i < strands never occurs) → exit 1: "split link; compute the pieces separately." (Q11.)
   - f. Locate and validate the installed library (C4 `locate`) → exit 4 on failure, with the exact fetch command.
   - g. Stamp preflight (C4) → exit 4 on mismatch.
   - h. Run F31, then F32b (skip F32b with `--single-prime`; skip it if F31 failed). F32b is 1.5–3× slower per point (C2_2_report.md:155-156). Classify errors as in C4.
   - i. Compare primes exactly as krhom.rs:131-139 (shared function) → exit 3.
   - j. Print the results.

**Output.** stdout matches the 2-strand path line for line (krhom.rs:63-64, 140, 167-173). tests/xengine.rs:631, 645 parse these lines.
```
gl_N for N in [2, 3, 4]  strands 3  crossings 4  components 1  colors [2]
braid: 1 -2 1 -2
[simplified braid: ...  (strands k)]
method: singular (localized box engine; library l222_n4 @ data ee548d2, labels (2,2,2), prune 4, box N 4, seed 201)[, reduced (strand at bottom position 0 marked)]
fields: F_2147483647 and F_4294967291 (results agree)        | fields: F_2147483647 only (--single-prime)
time X.XXs
N = 2
  Poincare polynomial: ...
  total rank: ...
  Euler characteristic: ...
```
- If the seeds differ, the method line says `seeds F31 201 / F32b 202`.
- Progress goes to stderr, per prime. Example: `[F_2147483647] library: 178 layer-A loads (0 rejected, k by route stamp), 3.2s; 216 box points, scan 55.0s, close 2.1s; chi = decat_chi for 3 closures; scaling: c entries checked; VmHWM m MB`. These are the same quantities sbox prints (sbox.rs:75, 167-168, 183).

**Exit codes**

| Code | Meaning |
|---|---|
| 0 | Success |
| 1 | Usage or unsupported input (existing usage: krhom.rs:41, 72, 81, 88, 103) |
| 2 | Engine or certification failure, including determinant failure after retries (existing: krhom.rs:122-125) |
| 3 | The primes disagree; no result printed (existing: krhom.rs:131-139) |
| **4 (new)** | Library not installed, wrong version, stamp mismatch, or the word needs uncertified jobs (they are named) |
| **5 (new)** | Resource refusal; the message names the flag to raise |
| 101 | Remaining panics; should no longer occur for bad input after step 0 |

**Message templates** (stderr, prefix `krhom: `)
- Exit 4, not installed: `library l222_n4 (F_2147483647) is not installed under <root>: run tools/fetch_library.sh l222_n4 (≈3.5 MB per prime)`.
- Exit 4, version: `installed l2222_n4/f31 is data <sha>, this krhom expects <sha>: run tools/fetch_library.sh --force l2222_n4`.
- Exit 4, stamp: `none of the <k> layer-A files of l2222_n4/f31 matches this binary's route stamp (provider m4c-lf); rebuild krhom from a branch that contains the c27 route_header filter (cca2db4)`.
- Exit 4, jobs: `library l2222_n4 does not cover this word: 11 jobs not certified for both primes: pair 926 (p=…, rung=…), …`.
- Exit 4, complete library: if a library whose status is `complete` reports a missing map, print `damaged install: tools/fetch_library.sh --force l222_n4`.
- Exit 5:
  - `per-point memory above --mem-mb 3000: rerun with a larger --mem-mb`;
  - `closure memory estimate above --mem-cap-mb 4500: rerun with --mem-cap-mb 200000 (the estimate is 30–130× too high at 4 strands) under ulimit -v`.
- Exit 2: `localized engine failed on F_2147483647: <error>. This is a certification failure; please report it with --json output.`

---

## C. Code changes, file by file

**C1. `src/sindec.rs`.** Only if the base is 73cbd7c: the `route_header` `c27=` filter, byte-identical to origin/colored-2c4-unblock. No other sindec edits.

**C2. `data/libraries.tsv` (new).** The single source of truth for both the shell tool and Rust (`include_str!`).
- Columns: `id field lambda n provider data_commit costs status layers`.
- `layers` lists, in restore order and separated by `;`, entries of the form `<data-gl6 path>@<sha256>`.
- Fill with the full SHA of ee548d2. These sha256 values come from the MANIFESTs at ee548d2:

| id/field | provider, costs, status | layers |
|---|---|---|
| l222_n4/f31 | default, `c2/costs_cloud.txt`, `complete` | `c2/l222_n4/f31/layerA_tmpl.tar.gz.part-000@7da63d261639995e1d828d788303ee23d1896862c01818257d5bb55adac1debd` |
| l222_n4/f32b | same | `…/f32b/…part-000@c111cc94fee71e8dc176ba8407c8f5586e59d84e6113e60f8d0a800668679407` |
| l2222_n4/f31 | m4c-lf, `c2/costs_2222_c25.txt`, `data/libstatus/l2222_n4_f31.txt` | `c2/l2222_n4/f31/…part-000@0ec12187…3ae9`; `c2/l2222_n4_c25/f31/…part-000@dd3de29a…2c49`; `c2/l2222_n4_c24/f31/layerA_tmpl_delta.tar.gz.part-000@6fdf3059…17fb`; `c2/l2222_n4_c27/f31/layerA_tmpl_delta.tar.gz.part-000@dcb4dda2…cde0` |
| l2222_n4/f32b | same, `…_f32b.txt` | `180d90ea…1de1`; `5d037a75…6e8f`; `07724ac5…b903`; `a407774a…7438` |

- Copy the full hashes from `git show ee548d2:<dir>/MANIFEST.txt`.
- Notes on provenance:
  - l222_n4 is a complete level set: 105/105 plus corner 21/21 (MANIFEST `sets:`).
  - The provider of l2222_n4 is m4c-lf, per the c25/c24/c27 MANIFEST `use:` lines.
  - Including the old `c2/l2222_n4` base follows the c27 restore line; see Q2.

**C3. `data/libstatus/l2222_n4_{f31,f32b}.txt` (new).** Copied verbatim from `origin/colored-2c4-unblock:runs_c2/c27/cover/status_{f31,f32b}.txt`, one `JOB ok|fail` line per job, 1254 lines. **Used only for diagnostics**, never as a gate. Pair indices depend on plan order (libbuild.rs:340-353).

**C4. `src/calc.rs` (new, `pub mod calc;` in src/lib.rs).** It lifts the setup in sbox.rs:51-90, minus the CLI, into Result-returning functions.
- **Types**
  - `enum CalcErr { Unsupported, Engine, Disagree, Library, Resource }(String)`, with `code()` returning 1, 2, 3, 4 or 5.
  - `struct RegEntry` and `fn registry() -> Vec<RegEntry>`, parsed from `include_str!("../data/libraries.tsv")`.
  - `fn lib_root(cli: Option<&str>) -> PathBuf`.
- **`fn locate(reg, root) -> Result<InstalledLib, CalcErr>`**
  - `<root>/<id>/<field>/LIBRARY.txt` must exist and parse (D).
  - id, lambda, n, field (`F::name()`: `F_2147483647` / `F_4294967291`, field.rs:125-127), provider and data_commit must equal the registry.
  - The `layer_a/` and `tmpl/` file counts must equal LIBRARY.txt.
  - `costs.txt` must start with `krhom-costs v2 alg=hs-kl1`.
- **`fn run_prime<F: PrimeField>(reg, inst, word: &[i32], ns, reduced, seed, retries, mem_mb, check_only) -> Result<PrimeOut, CalcErr>`**
  1. Stamp, by provider:
     - default: `(_, stamp) = sindec::provider_for_stamped::<F>(&lam, Some(&inst.costs), 1.0)`;
     - m4c-lf: `leaves::provider_for_stamped_m4c::<F>(&lam, Some(&inst.costs), 1.0, true)` (sbox.rs:67-71; leaves.rs:1030).
  2. `smapio::set_replay_only(true)` (smapio.rs:43). Precedent: libbuild.rs:360-363. sbox never sets it, so a missing tmpl would otherwise be recomputed silently (smapio.rs:350-351). See Q3.
  3. `sd = SIndec::with_stamp(&lam, Box::new(swire::NoProvider), stamp, SIndecOpts { cache_dir: Some(inst.dir), prune_n: Some(reg.n), layer_a_route_stamp: true, ..Default::default() })` (sbox.rs:72-73).
  4. **Stamp preflight.**
     - `want = route_header(&sd.layer_a_stamps())` (sindec.rs:1694, 457).
     - For every `layer_a/*.txt`, read the first line and strip the trailing ` key=…` (the header is stamps + ` key=` + key, sindec.rs:1742-1744). Count the lines where `route_header(prefix) == want`.
     - Zero matches → `Library` (stamp). Non-matching legacy files are normal; note their count on stderr.
  5. `lib = swire::build_wlib(sd, &[word.to_vec()], true)` (sbox.rs:74; swire.rs:144). Record the `lib.stats` fields `layer_a_loaded`, `layer_a_route`, `layer_a_rejected` and `provider_calls` (sindec.rs:119-135).
  6. `o = BoxOpts::new(reg.n as usize, seed, ns, vec![if reduced {Some(0)} else {None}], vec![word.len()], EngOpts { check_d2: false, mem_point_bytes: mem_mb << 20, ..Default::default() })`.
     - `check_d2` is false because the run is pruned; this is sbox's `prune.is_none()` (sbox.rs:83). D² is then certified by Freivalds, and `close_box` refuses pruned data without that certificate (boxclose.rs:1236-1240).
     - **Stages are fixed to `[word.len()]`.** sbox's default `closable_stages` is every prefix for uniform λ (boxrun.rs:485-496).
     - Leave `scaling`, `freivalds_k` and `batch` at their defaults.
  7. If `check_only`: run `swire::record_scan(&lib, word, &o)` (swire.rs:1320) and return. This catches missing closures, which only surface at stage extraction ("final type … has no closure generators", swire.rs:987).
  8. `out = swire::box_word(&lib, word, &o)`. On a determinant error, retry with seed+1 up to `retries` times, reusing `lib`.
  9. **Post-checks.** Any failure is `Engine`.
     - `out.stages.len()==1` and `out.stages[0].stage == word.len()`;
     - `out.scaling.is_some()` (sbox.rs:166);
     - `out.chi_checked == ns.len()`;
     - `out.stages[0].closed` contains `(n, mode)` for every n.
  10. Return `BTreeMap<usize, Poincare>` plus provenance: seed used, `lib_build_secs`, points, `scan_secs`, `close_secs`, `chi_checked`, scaling report, `det_pass`, lib stats, VmHWM (`xengine::rss_kb`, sbox.rs:42-44).
- **`fn classify(err: &str) -> CalcErr`.** One table with unit tests, using these exact source strings:

| Error text | Source | Class |
|---|---|---|
| `"not in layer A"` | swire.rs:1472 | Library |
| `"has no closure generators"` | swire.rs:987 | Library |
| `"no closure data in the library"` | swire.rs:976 | Library |
| `"smapio (replay)"` | smapio.rs:350-351 | Library |
| `"redraw the nodes"` | swire.rs:1220 | Retry, then Engine |
| `"per-point memory"` / `"above the budget"` | sengine.rs:73-78 | Resource |
| `"above the memory cap"` | cclose.rs:2727 | Resource |
| anything else (χ, scaling, point-dependent, Freivalds) | — | Engine |

- **`fn diagnose_missing(reg, word) -> Result<Vec<String>, String>`**
  - Runs only after a `Library` error, and only when `reg.status` is a file.
  - Logic as in examples/m5i_cover.rs:13-62: `Schur::new`, `pruned_plan(&mut s, Some(n))`, `knot_jobs`, keys `pair K` / `closure <table>`. A key counts as missing unless it is `ok` in **both** primes' status files.
  - Print `p` and `rung` next to each key (libbuild.rs:340-353 explains why).
  - Cost: about 23 s for the (2,2,2,2) plan (runs_c2/census/plan_2222_n4.out:3, per the lib report). It is spent only on failure.
- **`fn box_calc(...)`** runs `run_prime::<F31>`, then `run_prime::<F32b>`, then `compare_primes`.
- **`fn compare_primes(ns, a, b) -> Result<(), String>`**, refactored out of krhom.rs:131-139 and used by both singular routes.
- **`pub fn poly_json`**, moved from sbox.rs:38-41, plus `fn result_json(...)` (C6).

**C5. `src/lib.rs`:** add `pub mod calc;`.

**C6. `src/bin/krhom.rs`**
- Step-0 parse hardening.
- `--colors` broadcast.
- Replace krhom.rs:90-93 with dispatch steps 3–5.
- Before any engine thread on the box route, set the env vars:
  - `KRHOM_COSTS_ANY_HOST=1` unless already set. Required: a host mismatch makes the stamp `costs=rejected` (homprov.rs:1169-1173; sindec.rs:484-499), every file is then rejected and the error reads "not in layer A". Every shipped table says `4cpu-Intel…` (`git show ee548d2:c2/costs_cloud.txt | head -1`).
  - `KRHOM_MEM_CAP_MB`, from `--mem-cap-mb` or the default, unless the env var is set.
- Print the method line. Keep "fields:" and the N blocks byte-identical.
- `--json OUT` for **all** routes:
  ```
  {"tool":"krhom","commit":…,"route":"exact2|box|cable|cube|soergel",
   "input":{"braid":[…],"strands":…,"colors":[…],"labels":[…],"used_braid":[…],"used_strands":…,"reduced":…},
   "library":{id,field(s),data_commit,dir}|null,
   "primes":[{"field","seed","secs","lib_build_secs","points","scan_secs","close_secs","chi_checked","scaling":{checked,zero_both,declared_zero},"det_pass":[p,c],"layer_a":{loaded,route,rejected},"vmhwm_kb"}],
   "agree":true|false|null,
   "results":[{"N":n,"P":<poly_json>,"euler":"…"}]}
  ```
  Write it by hand with `format!`, as sbox does (sbox.rs:184-199).
- Add a hint on stderr when the default method resolves to the cube and λ is a supported uniform label: "hint: --colored-method singular". The default itself does not change (Q6).
- Update the usage strings (krhom.rs:3-27, 40).

**C7. `src/bin/sbox.rs`:** only replace the local `poly_json` with `calc::poly_json` (a pure move). No other change.

**C8. `tests/xengine.rs:650-652`:** update. See F1.

**C9. `tests/calc.rs` (new).** See F.

**C10. `tools/fetch_library.sh` (new).** See D.

**C11. `tools/randcheck.py` (new, Python 3 stdlib).** See E.

**C12. `docs/CALCULATOR.md` (new).** See G.

---

## D. Library fetch, verify and layout

**Canonical installed layout.** One full snapshot per (λ, N, prime):
```
<root>/                       # --lib-dir > $KRHOM_LIB_DIR > ${XDG_CACHE_HOME:-$HOME/.cache}/krhom/lib
  l222_n4/f31/  l222_n4/f32b/  l2222_n4/f31/  l2222_n4/f32b/
     LIBRARY.txt  costs.txt  [status.txt]  MANIFESTS/<layer>.txt  layer_a/a-<16hex>.txt  tmpl/tmpl-<16hex>.txt
```
(sindec.rs:1739-1741; tmpl.rs:334-338)

`LIBRARY.txt` uses key: value lines:
```
format: krhom-library v1
id: l2222_n4
lambda: 2,2,2,2
n: 4
field: F_2147483647
provider: m4c-lf
data_commit: <full sha>
layers: <path@sha256>; ...
costs: costs.txt <sha256>
status: status.txt | complete
files_layer_a: <count>
files_tmpl: <count>
installed: <UTC>
fetched_by: krhom <commit>
```

**Tool:** `tools/fetch_library.sh [--root DIR] [--git-dir DIR | --repo URL] [--force] [--verify-only] ID[/field] ...`
1. **Transport.**
   - Default: the current clone. Run `git fetch origin <data_commit>` (objects only; the work tree is untouched), then `git show <commit>:<path> > staging/<file>` for each chunk, the MANIFESTs and the costs table.
   - `--repo URL`: `git clone --filter=blob:none --no-checkout --depth 1 -b data-gl6 URL tmp`, then `git -C tmp show`.
   - Raw HTTPS is not used, because whether the repo is public is unknown.
2. **Verify each chunk.**
   - sha256 must equal both the registry value and the MANIFEST value. Use `sha256sum`, or `shasum -a 256` on macOS.
   - `tar -tvzf` on each layer's concatenated chunks (`cat parts | …`): every entry must be a regular file or directory matching `^(\./)?(layer_a|tmpl)/?$`, `layer_a/a-[0-9a-f]{16}\.txt` or `tmpl/tmpl-[0-9a-f]{16}\.txt`. No `..`, no absolute paths, no links.
   - The `layer_a` and `tmpl` entry counts must equal the MANIFEST `files:` line. Its word order varies ("files: layer_a 233, tmpl 471" vs "files: 249 layer_a, 233 tmpl …"), so parse both orders.
3. **Extract** into a staging dir **in registry order**: l2222_n4 → c25 → c24 → c27, as the c27 MANIFEST restore line says. Later layers may overwrite earlier names.
4. **Install.**
   - Copy `costs.txt`, `status.txt` (from data/libstatus) and the MANIFESTs.
   - Count the final files and write LIBRARY.txt.
   - `mv` staging → `<root>/<id>/<field>`. With `--force`, move the old dir aside first.
   - Idempotent: if LIBRARY.txt matches data_commit and counts, print "already installed" and exit 0.
5. **`--verify-only`:** recount, recheck LIBRARY.txt, and exit nonzero on mismatch.

**Sizes** (chunk bytes, `git ls-tree -l ee548d2`):
- l222_n4: 3.50 / 3.51 MB (f31 / f32b).
- l2222_n4: 6.97 + 7.68 + 1.24 + 0.27 = 16.2 MB (f31); about 18.2 MB (f32b).

**Integrity vs soundness:** the hashes catch corruption. Soundness comes from re-certification on load (sindec.rs:103-110), so a bad map is rejected rather than used.

**Load cost per run:**
- 3 strands: 0.16–3.2 s (runs_c2/k3/markov/s3p_f31.json `lib_build_secs` 0.16; 4_1 about 3.2 s per the lib report).
- 4 strands: 21–25 s (runs_c2/k4m/markov/m1,m3,w41s_f31.json: 21.19, 21.31, 24.88).

---

## E. Random-example checker: `tools/randcheck.py`

**Interface:** `randcheck.py --seed S [--n 24 | --quick | --full] [--tiers 2,3,4] [--only ID] [--krhom target/release/krhom] [--jobs 2] [--out runs_random] [--inject-fault]`
- `--jobs` processes run concurrently, each with `KRHOM_THREADS=2` and `RLIMIT_AS` 8 GB (`resource.setrlimit` in `preexec_fn`). Production guards were `ulimit -v` 6 GB at 3 strands and 14 GB at 4 (runs_c2/k3/one.sh:13; runs_c2/k4m/markov/one.sh:9). Measured VmHWM at 4 strands is 99–302 MB (m1/m3/w41s json).
- Results are parsed only from `krhom --json`.
- Each example gets `subseed = rng.getrandbits(32)`, written to the report.
- **Relations are fixed in code before running. There are no tolerances.**

**Harness self-tests.** Failure of any → exit 2.
1. **Convention check.** `krhom 2,3,4 --braid "1 1 1" --strands 2 --colors 2 --colored-method singular` must equal `python3 tools/wang.py 2 N 1` unreduced, for N = 2, 3, 4.
   - wang.py prints in krhom conventions: ours(t,q) = Wang(1/t,1/q) (wang.py:1-4, 175-182).
   - Read wang.py to confirm it is unreduced. The expected N=4 unreduced value is in runs_c2/k3/markov/s3p_f31.out (rank 22).
2. **Anchor.** `krhom 2,3,4 --braid "1 -2 1 -2" --strands 3 --colors 2 --colored-method singular` (with and without `--reduced`) must equal the strings in H3.
3. Every tier's library is installed; otherwise exit 2 with the fetch command, unless `--tiers` excludes that tier.
4. **`--inject-fault`** perturbs one coefficient of the first box result in memory. That example must then be FAIL and the exit code 1, which proves the comparisons are not vacuous.

**Oracles.** Only those the reports show valid.
- Labels: **C** = certifying. **B** = blocking: a mismatch is FAIL, but a pass is not counted as confirmation. **E** = expected, not proven: a mismatch is WARN and never counts as a pass.
- A "reduced" comparison between two words of a **link** is made only when the component through bottom position 0 is provably the same. That holds for mirror (same permutation), reversal (inverse permutation, same cycles), and stabilization at the top. Knots are always fine (tests/xengine.rs:600-601 rule).

| # | Check | Relation and normalization | Label | Source |
|---|---|---|---|---|
| O1 | Two primes | Built into krhom; exit 3 means FAIL | C | krhom.rs:126-141 |
| O2 | Seed independence | `--seed S1 --seed-f32b S2`, S1≠S2. On exit 3, rerun with S1=S2 to separate prime dependence from seed dependence. | C | SPEC_3C_GL6_rev1.md:111; sbox.rs:59 |
| O3 | χ = decat_chi | Built in and fatal on the box route; the checker records `chi_checked == #N` | C (missing summands only) | swire.rs:1182-1203, 1366 |
| O4 | Isotopy: rotation, conjugation σ_j^{±1}, braid relation, cancelling pair σ_jσ_j^{-1} | Equal, no shift | C | tests/xengine.rs:600-611 (2 strands) |
| O5 | Markov: insert σ_n^{±1} anywhere; the new strand has the same color | Equal, no shift (t^0q^0) | C | C2_1_report.md:103-105; C2_2_report.md:143-153; runs_c2/k4m/markov/compare.txt |
| O6 | Mirror (negate every exponent) | Terms (t,q,k) → (−t,−q,k), both modes | C | tests/xengine.rs:95-101, 499-510, 593-598 |
| O7 | Reversal rev(w), N'=3,4 | Equal | E | not run before (R4_m6a.md:322 plans the (3,3,3) analogue) |
| O8 | N'=3 duality | Box N'=3 of w equals `krhom 3 --braid "rev(w)" --strands n [--reduced] --json` (uncolored Soergel, explicit `--strands`); no shift | B | C2_1_report.md:69-73; C2_2_report.md:220-221 |
| O9 | N'=2 | Exactly `t^0q^0`, both modes, every link | C (weak) | k3_c2_compare.py:50 |
| O10b | Lewark-type relation | [N 2]_q·P_red − P_unred = (1+t)·Q with Q ≥ 0, by peeling (1+t) at each q from the lowest t (tests/xengine.rs:126-160). Knots only. | E | tests/xengine.rs:515-521 ("not a theorem") |
| O11 | 2-strand closed forms, n = exponent sum | odd \|n\|≥3: wang.py `2 N (|n|-1)/2`, mirrored when n<0; \|n\|=1: [N 2]_q unreduced, 1 reduced; n=0: [N 2]_q² unreduced; even \|n\|≥4: none (skip) | C | tests/xengine.rs:282-301, 391-419 |
| O12 | Unknot words: each σ_i exactly once, any order and signs, optionally conjugated | Unreduced [N 2]_q, i.e. N'=2 `t^0q^0`, N'=3 q^-2+1+q^2, N'=4 q^-4+q^-2+2+q^2+q^4; reduced `t^0q^0` | C | compare.txt m1 ranks 1,1,3,1,6,1 |
| O13a | 2-strand cross-engine | `--box-live` equals the exact path at N'≤4, both modes | C | tests/swire.rs:88-119 |

**Excluded:**
- an independent so(6) Kauffman χ (not implemented; SYNTHESIS.md:349, 389);
- pruned vs unpruned (O13b: 143–278 s per prime, P2);
- amphichirality from tables (P2 Mode C).

**Generator.** Colors are 2 on every component; signs are ±1 with p = 1/2. Base words contain every σ_1..σ_{n−1} and no cyclically adjacent σσ⁻¹. Box-route calls always pass `--no-simplify`; otherwise Mode A words destabilize back to 2 strands.

- **T2 (2 strands; exact path; seconds per call):**
  - length 1–9; N' = 2..6;
  - partners: the mirror word and one rotation;
  - checks: O1, O4, O6, O8, O9, O11, O13a (N'≤4), O10b.
- **T3 (3 strands; l222_n4):**
  - **Mode A (half of the examples).**
    - The core is a 2-strand word of length 1–4; its value comes from the exact path, itself checked by O11 when it applies.
    - Insert σ2^{±1} at a uniform position, then apply 1–2 moves from {rotate, conjugate, braid relation where the pattern occurs cyclically, cancelling insert}. Total length ≤ 7.
    - The result must equal the core value (O4+O5).
    - Two of these examples use the core σ1^{±1}, giving an unknot (O12).
  - **Mode B (the other half).**
    - A free word of length 4–6 with both generators.
    - One partner chosen uniformly from {rotate, conjugate, braid relation or cancelling insert (O4), mirror (O6), reverse (O7)}; partner length ≤ 8.
  - Every box output also gets O2, O3, O8 (N'=3), O9 and O10b (N'=3,4).
  - Word lengths are capped because 8-letter 3-strand knots take 100–830 s per prime (runs_c2/k3/out); 4–6-letter knots take 54–75 s per prime with both modes.
- **T4 (4 strands; l2222_n4; `--reduced` only, plus unreduced when length ≤ 3):**
  - Only knots (one component) and unknot words.
  - **Mode A:** the core is a 3-strand knot word of length ≤ 3, or a 2-strand word with σ2 inserted. Insert σ3^{±1}, optionally rotate once; total length ≤ 5.
  - The reference is the core's value: box route at 3 strands for 3-strand cores, the exact path for 2-strand cores.
  - **Coverage filter, required:**
    - Draw up to 50 candidates.
    - Run `target/release/examples/m5i_cover 2,2,2,2 4 cands.txt both.txt` once, where `both.txt` = the jobs `ok` in **both** data/libstatus files (examples/m5i_cover.rs:1-4, 12, 63).
    - Keep only `missing 0`. If none survive, the example is SKIP(library).
  - Measured 4-strand cost per prime (F31, both modes, 2 threads): σ1σ2σ3 90 s, "1 1 1 2 3" 98 s, "1 -2 1 -2 3" 558 s (runs_c2/k4m/markov/*_f31.json `total_secs`).
- **Timeouts per krhom call:** T2 120 s, T3 900 s, T4 1500 s. A timeout is TIMEOUT, never PASS.
- **Sizes:**
  - default `--n 24`: T2 8, T3 12, T4 4, plus the 2 anchors. Wall time on 4 cores with `--jobs 2` is about 40 min — **an estimate, not measured** (oracle report §6).
  - `--quick`: 4 / 6 / 0, about 15 min (estimate).
  - `--full`: 10 / 20 / 10, about 1.5 h (estimate).

**Report**
- `runs_random/<seed>/examples.jsonl`, one line per example: `{id, subseed, tier, lam, word, components, mode, core{word,value}, partner{kind,word}, runs[{word, args, exit, secs, json}], checks{name:{result:PASS|FAIL|WARN|NA|SKIP, label:C|B|E, detail}}, status:PASS|FAIL|SKIP(library)|TIMEOUT|ERROR}`.
  - `detail` on FAIL is the term diff, or the negative Q coefficient for O10b.
  - Every FAIL prints its repro commands.
- `summary.md`: one row per example (id, tier, word, partner, comps, N'=4 ranks u/r, checks, status, secs).
- Final line: `seed=… n=24 PASS a FAIL b WARN c SKIP d TIMEOUT e; checks: primes+seeds x, isotopy y, markov z, mirror …, duality-N3 …, reversal-N4 …, lewark …, unknot …, closed-form …, cross-engine …`.
- **Exit codes:** 0 if no FAIL or ERROR; 1 on any FAIL; 2 on harness failure.
- This is a validation log, not a table. Commit only the acceptance runs' `summary.md` and `examples.jsonl`.

---

## F. Tests to add

**F1. Default suite: CLI refusals and plumbing** (`tests/calc.rs`, using `CARGO_BIN_EXE_krhom` with `KRHOM_LIB_DIR` set to a fresh empty temp dir and `KRHOM_THREADS=1`)
- **Update tests/xengine.rs:650-652.**
  - `krhom 3 --braid "1 2" --strands 3 --colors 2 --colored-method singular` now **succeeds**: it destabilizes to the 1-strand unknot. Assert `simplified braid:` and `total rank: 3` ([3 2]_q).
  - Move the failure assertion to the next bullet.
- `krhom 2,3,4 --braid "1 -2 1 -2" --strands 3 --colors 2 --colored-method singular` → exit 4; stderr contains `tools/fetch_library.sh l222_n4`.
- Exit-1 refusals, each with a message substring:
  - mixed labels: `--braid "1 1 2 2" --strands 3 --colors 2,2,1` (a 3-component chain; simplify cannot drop strands because the search stops when strands ≤ ncomp, pd.rs:343);
  - colors (3,3,3): `--braid "1 -2 1 -2" --strands 3 --colors 3`;
  - 5 strands: `--braid "1 1 2 2 3 3 4 4" --strands 5 --colors 2`;
  - N above the library: `krhom 5 --braid "1 -2 1 -2" --strands 3 --colors 2`;
  - uncolored with singular on 3 strands: `--braid "1 -2 1 -2" --strands 3`;
  - wrong number of colors: `--braid "1 1" --strands 2 --colors 2,2,2` (exit 1, not 101);
  - split braid: `--braid "1 1" --strands 3 --colors 2 --no-simplify`.
- `"1 1 1 2"` on 3 strands (default simplify) gives the same N blocks as the 2-strand `"1 1 1"`.
- `--colors 2` equals `--colors 2,2` for `--braid "1 1" --strands 2`.
- **Cross-engine parity:** `krhom 2,3,4 --braid "1 1 1" --strands 2 --colors 2 --colored-method singular --box-live` has the same N blocks as the exact path, with and without `--reduced`, and prints `fields: F_2147483647 and F_4294967291 (results agree)`.
  - `--box-live` builds the library in memory, unpruned (`check_d2: true`), with box N = max(ns), as tests/swire.rs:22-25, 88-96 do.
  - It never sets replay-only and refuses more than 2 strands.
- `--json` smoke test on the exact route: the keys are present, and `results[].P.rank` matches stdout.

**F2. Unit tests in src/calc.rs**
- registry parse: 4 rows, 64-hex hashes, layer order;
- LIBRARY.txt parse and mismatch errors;
- `classify` table with the exact source strings listed in C4;
- `compare_primes` agree and disagree;
- the reduced-link simplify guard.

**F3. `#[ignore]` tests that need `KRHOM_LIB_DIR` with both libraries installed**
- 4_1 strings (H3);
- `"1 1 1 2" --no-simplify` (H4);
- `"1 2 3"` and `"1 1 1 2 3" --no-simplify --reduced` (H5);
- `--check-library` for 6_1 and 8_15 (H6).

**P2:** a tiny fixture under tests/data/lib222_tiny for the stabilized trefoil (36 layer-A loads, 1.7 s: runs_c2/k3/markov/s3p_f31.json). Which files to extract is unknown: no load log exists.

---

## G. Documentation

There is **no README.md** on this base: the worktree root has none, and neither does the origin/dev tree. The README draft is local to the orchestrator.

The session writes **`docs/CALCULATOR.md`** with these sections:
1. Usage, with examples;
2. Support matrix: λ, N, strands, library id, size;
3. Installing libraries (fetch tool, root, `--force`, `--verify-only`);
4. What runs (both primes, χ, scaling, Freivalds, determinant pass, seed 201, prune 4);
5. Exit codes and messages;
6. Cost guide with the measured numbers above, plus `--reduced` for speed and `ulimit -v`;
7. Random checker usage and report format;
8. sbox as the advanced/production tool.

For the orchestrator, these places in the README draft must change (line numbers per the docs report; R = README_draft.md):

| Lines | Change |
|---|---|
| R:6-8, 10-11, 18-19, 22, 24, 26 | Reframe from table production to calculator plus support matrix |
| R:43-46, 65-69, 71-72, 49-51 | Box-route results are pruned; krhom runs both primes on 3–4 strands; fixed seed; PD→braid may exceed 4 strands |
| R:95-121 | Quick start: fetch, then one krhom command; keep the 1/11/46 and 1/5/25 example |
| R:122 | Test counts; this branch was 287/37 (C2_4_report.md:35) |
| R:129, 131, 140 | krhom and sbox rows; sbox needs `--provider` |
| R:149-150 | KRHOM_MEM_CAP_MB and KRHOM_COSTS_ANY_HOST are now set by krhom on the box route |
| R:166-169 | χ and scaling happen in "sbox and krhom" |
| R:187-261 | Results become a short "validated examples" section |
| R:263-278 | Library layout and fetch; the (2,2,2,2) row is now l2222_n4 + c25 + c24 + c27 with m4c-lf |
| R:284-300, 307 | Limitations: M5i and C2-7 done; 18 class-A knots covered; mixed labels refused |
| R:322-323, 344, 138 | Repository map additions |

The docs report also flags stale README claims: (a) the Cargo.toml description, (h) gl6 paths, (j) the KRHOM_GRCACHE default.

---

## H. Acceptance criteria

Run on the cloud machine, release build, `KRHOM_THREADS=4` unless stated.

**H1. Test suite.** `cargo test --release`: 0 failed. Report the passed/ignored counts against the base; m5i was 287/37.

**H2. Fetch tool.**
- `tools/fetch_library.sh l222_n4 l2222_n4` installs 4 dirs.
- A second run prints "already installed".
- Deleting one `layer_a` file makes `--verify-only` fail, and krhom then exits 4 with a `--force` hint.

**H3. 4_1 on 3 strands.** `krhom 2,3,4 --braid "1 -2 1 -2" --strands 3 --colors 2 --colored-method singular`: exit 0, "results agree", ranks 1, 11, 46.
- N=3: `t^-2q^6 + t^-2q^8 + t^-1q^0 + t^-1q^2 + t^0q^-2 + t^0q^0 + t^0q^2 + t^1q^-2 + t^1q^0 + t^2q^-8 + t^2q^-6`
- N=4: identical to the unreduced N=4 string in runs_c2/k3/out/4_1_f31.json (rank 46).
- With `--reduced`: ranks 1, 5, 25. N=3 is `t^-2q^6 + t^-1q^2 + t^0q^0 + t^1q^-2 + t^2q^-6`; N=4 must equal the JSON's reduced string.
- N=2 is `t^0q^0` in both modes.
- Per-prime time is about 1 min (62 s measured with u,r).

**H4. Stabilized trefoil.** `krhom 2,3,4 --braid "1 1 1 2" --strands 3 --colors 2 --colored-method singular --no-simplify`: same N blocks as `--braid "1 1 1" --strands 2`. Ranks 1, 7, 22 (reduced 1, 3, 9).
- N=3 unreduced: `t^0q^-6 + t^0q^-4 + t^0q^-2 + t^2q^-8 + t^2q^-6 + t^3q^-14 + t^3q^-12`
- N=3 reduced: `t^0q^-4 + t^2q^-8 + t^3q^-12`
- Source: runs_c2/k3/markov/s3p_f31.out; consistent with tests/xengine.rs:647.

**H5. 4 strands, reduced and unreduced** (Markov against 2 strands):
- `--braid "1 2 3" --strands 4 --no-simplify`: ranks 1, 3, 6; reduced 1, 1, 1.
- `--braid "1 1 1 2 3" --strands 4 --no-simplify`: N blocks equal to H4's 2-strand output, both modes (compare.txt m3: 1,1,7,3,22,9).
- About 1.5–3 min for F31 and 2–6 min for F32b per call (estimate from the 98 s F31 u,r run).

**H6. `--check-library` at 4 strands:**
- 6_1 `--braid "1 1 2 -1 -3 2 -3" --strands 4 --colors 2 --check-library` → exit 0.
- 8_15 `"1 1 -2 1 3 2 2 2 3"` → exit 4, listing pairs 926, 949, 986, 997, 1054, 1412, 1435, 1471, 1472, 1483, 1540 (origin/colored-2c4-unblock runs_c2/c27/cover/cover_f31.txt).

**H7. Refusal and exit-code matrix:** all F1 cases produce the specified code and message.

**H8. Library untouched:** after H3–H6, the installed library dirs are unchanged (file counts, and `find <lib> -newer LIBRARY.txt` is empty). Otherwise report what was written (Q4).

**H9. Reader compatibility** (only for the fallback base 73cbd7c + cherry-pick):
- Replay the 6 c27 jobs on a scratch copy of the installed l2222_n4, for each prime, with `KRHOM_LEAVES_C27` unset:
  `KRHOM_COSTS_ANY_HOST=1 target/release/libbuild one --lam 2,2,2,2 --n 4 --field F --cache COPY --costs COPY/costs.txt --provider m4c-lf --stamp-route --replay -- "pair 1235"` (and likewise for 1444, 1470, 1894, 1769, 1812).
- All `"ok":true`. Command form from origin/colored-2c4-unblock runs_c2/c27/lib/runjobs.sh.

**H10. Checker quick run:** `tools/randcheck.py --seed 20261115 --quick` → exit 0. Also `--inject-fault --quick --only <first T3 id>` → exit 1.

**H11. Checker default run:** `tools/randcheck.py --seed 20261115` → exit 0, no FAIL. Commit `runs_random/20261115/{summary.md,examples.jsonl}`. **If any FAIL appears, stop and report the repro; do not tune.**

**Time budget**

| Work | Estimate |
|---|---|
| Implementation (P0+P1) | about 5–6 h |
| H3–H6 | about 25 min |
| H9 (if needed) | about 6 min |
| H10 | about 15 min |
| H11 | about 40 min |

---

## Open questions (decide or report; each has a default)

- **Q1. Base branch.** Default: colored-2c4-unblock 194a731, which carries unreviewed C2-7 commits. Fallback: 73cbd7c plus the cherry-picked `route_header` hunk plus H9.
- **Q2.** Is the old `c2/l2222_n4` base layer needed under c25? The c27 restore line includes it, while the lib report advised against it. Default: follow the MANIFEST order. Report the final counts.
- **Q3. `set_replay_only(true)`.** No production sbox run used it (only libbuild.rs:363 sets it), so whether every tested word has all its templates in the libraries is unknown. Default: on. If H3–H5 hit `smapio (replay)`, report it, and switch the default off only with the orchestrator's approval. Recomputing a template is deterministic and does not affect soundness.
- **Q4.** Does a replay write into the library dir (`create_dir_all` at sindec.rs:558-562; Rickard templates re-obtained "disk cache or recomputation", smapio.rs:20-24)? Answered by H8.
- **Q5.** N' below the color on the box route: refuse (default) or print 0 (cclose.rs:2721-2722 returns 0 for label > N)? Untested end to end in the box route.
- **Q6.** An `auto` default that picks singular when a library is installed. Default: no; the session only prints a hint. Today's default for (2,2,2) is the full cube, which ran more than 90 min on 4_1 (C2_1_report.md:106-109).
- **Q7.** (3,3,3), and l222_n6_unpruned for N'=5,6, are excluded. m5b n4 is 47/47, but n5/n6 are partial (m4c_b4w MANIFEST), and N'=5,6 at box N 6 was never run.
- **Q8. Resumable `--ckpt`** (`box_run` + `box_finish`, boxrun.rs; sharded precedent origin/colored-2c4-unblock runs_c2/c27/k4/box.sh) for multi-hour 4-strand words. Deferred to a follow-up.
- **Q9. Memory guard.** At 4 strands, `--mem-cap-mb 200000` effectively disables the estimate check. Should krhom set `RLIMIT_AS` itself? That needs libc or `extern "C"`, and no new dependencies are allowed. Default: document `ulimit -v`.
- **Q10.** Different seeds per prime by default. Default: no (201/201, as in production). The checker uses different seeds.
- **Q11.** Split braids on ≥3 strands are refused. A product formula could be added later.
- **Q12.** Labels: mirror is C here (theorem; exact on 2 strands), while the repo labels amphichirality B (SPEC_3C_GL6_rev1.md:114-116). Reversal and Lewark are E.
- **Q13.** Single-mode box grids (unreduced only, or reduced only) were never run in production: C2 runs used u,r. H3 and H5 validate them against u,r values.
- **Q14.** Exact cost of `knot_jobs` itself (plan about 23 s at (2,2,2,2)). It is used only on failure.

## Report back (final message of the session)

- branch and commits;
- test counts;
- stdout of H3–H6;
- the H7 table;
- H8 and H9 results;
- checker final lines for H10 and H11, with any WARN details;
- per-prime timings and VmHWM;
- answers to Q2–Q4 and Q13;
- anything that deviated from this spec.
