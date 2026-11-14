## Adversarial review of the calculator/randcheck spec against the code (worktree m5i @ 73cbd7c, read-only)

**Verdict.** Most citations are correct. Three defects would break the build as specified, two would give misleading library diagnostics, and several numbers are wrong. The spec matches the user's request (a calculator, no table sweeps, a few random examples).

The three defects that would break it:
- **Every unreduced box run would fail.** The `chi_checked` post-check is wrong (K1).
- **Missing templates are never classified correctly.** The "smapio (replay)" row in the classifier can never fire (K2).
- **A healthy install would be reported as damaged.** The library file-count gates break once a run writes template files (K3).

### 1. "Two facts changed": all CONFIRMED
- **origin/data-gl6 = ee548d2.**
  - c2/l2222_n4_c27/{f31,f32b} holds 12/170 and 14/162 layer-A/tmpl files (MANIFESTs).
  - Restore order is l2222_n4, c25, c24, then c27.
  - All 10 sha256 values in C2 match the MANIFESTs.
  - Sizes: 16.2 MB (f31) and 18.2 MB (f32b) (`git ls-tree -l`).
- **origin/colored-2c4-unblock = 194a731.**
  - Status files: 1254 lines each, 1236 ok / 18 fail, identical between primes.
  - The 6 blockers replay ok on both primes (runjobs_6_{f31,f32b}_replay.log).
  - cover_f31: the 18 class-A knots miss 0 jobs; 8_15, 9_22 and 11n_42 miss 11; Conway misses 16. The 8_15 pair list matches H6.
  - box_6_1_f31.log is 0 bytes.
- **route_header filter.** The `c27=` filter is at origin/colored-2c4-unblock src/sindec.rs:455-466.
  - On m5i, c27 headers keep the token, so they are rejected (sindec.rs:1749-1755).
- **Diff 73cbd7c..194a731 (WRONG in detail).**
  - Under src/ it is only gp0.rs, leaves.rs and sindec.rs.
  - It also changes Cargo.toml (description line only, from the dev merge) and adds tests/c27_joint.rs.
  - With `KRHOM_LEAVES_C27` unset, the m4c-lf stamp is unchanged: the token is appended only `if c27` (unblock leaves.rs diff, around new line 1405).

### 2. Spot checks: CONFIRMED (file:line)

**krhom.rs**
- Lines 43, 50-60, 61, 63-64, 84, 90-93, 117-142, 131-139 and 167-173.
- `exit(1)` at 41, 72, 81, 88, 91 and 103.

**sbox.rs**
- seed default 1 (59).
- `--stamp-route` (73).
- `check_d2 = prune.is_none()` (83), with mem default 2048.
- `poly_json` (38-41).

**Production scripts**
- k3/one.sh:13-14: seed 201, mem 3200, MEMCAP 4500, ulimit 6 GB, THREADS=2, no `--stamp-route`.
- k4m/markov/one.sh:9-11: 3000 / 200000 / 14 GB, `--stamp-route`.

**swire.rs**

| Item | Line(s) |
|---|---|
| `build_wlib` | 144 |
| `extract_stage` errors | 976, 987 |
| `BoxOpts` (`scaling`, `freivalds_k`) | 1117-1129 |
| `chi_check` | 1182-1199 |
| determinant error | 1220 |
| `record_scan` | 1320 |
| `box_word` | 1375-1426 (`BoxNodes::new` is at 1379, not 1378) |
| `NoProvider` "not in layer A" | 1472 |
| `check_pivot_values` forced on | 464 |

**sindec.rs**
- 103-110: `layer_a_route_stamp` is sound because loaded maps are re-certified (dechom_postcheck at 1236-1238).
- 126-134: stats.
- 457: route_header.
- 484-499: `costs=rejected`.
- 558-561: `create_dir_all`.
- 1694 and 1742-1744: header format.
- 1999 and 2109: planner.

**Other files**

| File | Line(s) | What |
|---|---|---|
| homprov.rs | 37-54 | `g1_work` |
| homprov.rs | 1169-1173 | host check |
| smapio.rs | 43, 350-351 | replay-only |
| sengine.rs | 73-78, 950-951 | per-point memory |
| cclose.rs | 749-752, 2721-2728 | memory cap |
| boxclose.rs | 1171-1176, 1208-1211, 1236-1240 | box closure checks |
| pd.rs | 147-149, 316-318, 343, 568-585 | simplification |
| libbuild.rs | 340-353, 360-363, 781 | plan order, replay-only, `set_var` |
| leaves.rs | 1030 | m4c provider (m5i base) |
| field.rs | 125-127 | field names |
| tmpl.rs | 334-338 | template file names |
| examples/m5i_cover.rs | 1-63 | cover tool |

**Runs and reports**
- c2/costs_cloud.txt is byte-identical to tests/data/homprov/costs_cloud.txt (424 records, matching the l222_n4 stamp exactly).
- l222_n4 replay is 105/105 + 21/21 on both primes (data-gl6 c2/runs/*_replay.out).
- Plan time is 23.1 s (runs_c2/census/plan_2222_n4.out:3).
- Suite was 287/37 (C2_4_report.md:35; runs_c2/k4m/full_suite.log).

**Acceptance values**
- **H3:** runs_c2/k3/out/4_1_f31.json has the N=3 strings as given, ranks 1/11/46 and reduced 1/5/25, total 62.45 s.
- **H4:** s3p_f31.out matches; 36 layer-A loads, 1.7 s.
- **H5:** m1 ranks are [1,1,3,1,6,1]; m3 matches H4.
- No README.md at root, in origin/dev or in unblock (but data/README.md exists).

### 3. Oracle validity (E)
- **O1 and O2: valid.** SPEC_3C_GL6_rev1.md:111 lists "two primes × two node sequences" as certifying.
- **O3: valid.** chi_check is fatal (swire.rs:1366); for reduced it multiplies by [N λ0] (1189-1191).
- **O4: valid.** Braid relation and cancelling pair are equal braids, so the marked point is unchanged; list them as safe for reduced links too. Rotation and conjugation can move position 0 to another component, so the spec's restriction is needed.
- **O5 (Markov, no shift): valid, but the citation is WRONG.**
  - runs_c2/k4m/markov/compare.txt compares C2-4 against C2-2 outputs (compare.py header). The Markov evidence is runs_c2/k4/markov/compare.txt and C2_2_report.md:143-153, 223.
  - Both signs give no shift: s3p = s3n (3 strands); m3 = m3n = m3nn = m3x (4 strands, identical strings).
  - Interior insertion was never run (runs_c2/k4m/stab_words.txt has no results) but is a conjugate of an end stabilization, so the theorem covers it.
- **O6 (mirror): valid.**
  - Unreduced: tests/xengine.rs:499-506. Reduced: 593-598.
  - Amphichiral N'=4 symmetry: C2_1_report.md:99-100.
- **O7: E is the right label.**
- **O8 (N'=3 duality, unreduced and reduced): valid.**
  - 34 knots: C2_1_report.md:69-73. 11 words: C2_2_report.md:220-221.
  - Links: only Hopf (tests/xengine.rs:442-449; u_s3 and m2 in runs_c2/k4/markov/compare.txt). B is fine.
- **O9 for links: supported.**
  - Knots: k3_c2_compare.py:50.
  - Hopf: tests/xengine.rs:311-316 gives q^{-2 lk a(N-a)} = q^0 when a = N.
  - u_s3 Hopf at 4 strands has rank 1.
- **O10b: E is fine.** tests/xengine.rs:127-159, 517-519.
- **O11: partially unverified.**
  - wang.py was compared against the engine only for (2,2..10,1), (2,2..6,2) and (2,3..4,3) (tests/data/wang_ref.txt keys; tests/xengine.rs:391-409). (2,2,1) = t^0q^0 is CONFIRMED.
  - T2 at length ≤ 9 and N' ≤ 6 reaches m=4 and (m=3, N=5,6), which were never cross-checked. Use label B there or restrict.
  - Add the Hopf closed form for |n|=2 (tests/xengine.rs:343-350, hopf_formula).
  - Add reduced n=0 = [N 2]_q (tests/xengine.rs:613-619).
- **O12: valid.** m1 matches.
- **O13a: valid.** tests/swire.rs:88-119; the post-check pattern is at tests/swire.rs:96.

### 4. Findings

**K1 (critical): the `chi_checked` post-check is wrong for unreduced runs.**
- With modes=[None], close_stage_box closes both (N',None) and (N',Some(0)) for every N'. The reduced closure is skipped only if `o.modes` contains Some(0) (swire.rs:1351-1363, line 1356).
- chi_check counts every closure (1187-1198), so `chi_checked` is 2·|ns| in unreduced mode and |ns| in reduced mode.
- The spec's `chi_checked == ns.len()` therefore turns every unreduced 3- and 4-strand call into exit 2.
- Fix: use `chi_checked == Σ closed.len()` (as tests/swire.rs:96 does), require (n, mode) to be present for every n, and filter `closed` by mode when printing and comparing primes.
- Note: an unreduced run also yields reduced closures from the unreduced extraction (boxclose.rs:1162-1163) at no extra scan cost. Production never used that path (all C2 runs were u,r), and no box_word test runs a single mode (tests/boxrun.rs:179 runs only record and scaling), so Q13 holds.

**K2 (high): "smapio (replay)" never reaches the caller.**
- layer_a_load swallows import errors: it counts them in `layer_a_rejected` and eprints "[sindec] layer A …: … (rejected)" (sindec.rs:1756-1771).
- NoProvider then fails with "replay: … not in layer A".
- build_wlib consumes `sd`, so `lib.stats` is unavailable on failure.
- As a result, a missing template or a per-file stamp rejection is reported as "does not cover this word", while `diagnose_missing` finds 0 missing jobs.
- Required: when a Library error has 0 diagnosed missing jobs (or status is `complete`), print "library files present but rejected on load (see [sindec] lines; template cache or stamp); fetch --force".

**K3 (high): runs write template files into the library.**
- Replay never writes layer A: saves happen only after provider success (sindec.rs:1265-1277, 1638-1645), and NoProvider disables the M5i Exact fallback (1252, 1627). This partly answers Q4.
- But TmplLib auto-stores any hom basis or Rickard template it computes (tmpl.rs:479-490, 511-529). replay_only guards only H lines in layer-A imports (smapio.rs:350).
- `locate` and `--verify-only`/idempotence use count equality, so after the first run that writes a template they would report a healthy install as damaged.
- Fix: keep a file list with hashes and allow extra files, or install read-only (store errors are ignored: `let _ =` at tmpl.rs:489, 528).

**K4 (high): the stamp preflight cannot catch the case it targets.**
- On an m5i binary only the 12/14 c27 files fail the route stamp; about 2.9k other files match, so the zero-matches check never fires.
- Fix: require every file of the c27 layer (or every m4c-lf-stamped layer) to match, using per-layer file lists in LIBRARY.txt.

**K5 (high): F32b cost factor is wrong.**
- C2_2_report.md:155-156 says 3-10× slower (with contention), not 1.5-3×.
- Measured without contention:
  - 4 strands: 1.0-3.1× total (k4m json, e.g. m3 98→263 s, w41s 558→1735 s, m3nn 183→189 s).
  - 3 strands: 0.9-1.2× (runs_c2/k3/out).
- 4-strand library load is 20-34 s (m3n 33.7, u_s3 34.3), not 21-25 s.

**K6 (high): T3 timings and timeouts.**
- 8-letter 3-strand knots take 100-1151 s per prime, not 100-830 s; 10-letter knots take 247-2564 s (runs_c2/k3/out).
- With a 900 s timeout per call covering both primes, partners of length 7-8 can time out.
- Cap partner length at 6 or raise the T3 timeout to about 2400 s.
- The 40-minute default-run estimate is optimistic (roughly 60-90 min).

**K7 (high): the closure memory cap scales with thread count.**
- `staged_mem_estimate` sums the top `num_threads()` per-type and per-stage terms (cclose.rs:2654-2658, 2703-2704).
- Production ran 2 threads with a 4500 MB cap. H runs use 4 threads, and an unset `KRHOM_THREADS` means all cores (lib.rs:48-55), so a word that passed in production can be refused with exit 5.
- Fix: scale the default cap by threads/2, or document it.

**K8 (high): test and acceptance commands missing `--colored-method singular`.**
- The F1 command `krhom 5 --braid "1 -2 1 -2" --strands 3 --colors 2` falls back to the default method. cable_default is false for (2,2,2) (krhom.rs:16-18), so it runs the full cube, which ran more than 90 min on 4_1 (C2_1_report.md:107-109).
- The H6 commands also lack the flag. `--check-library` should require or imply singular.
- The updated tests/xengine.rs:650 test must set `KRHOM_LIB_DIR` to an empty directory.

**K9 (medium): step 0 misses several panics.**
- krhom.rs:43 (N list parse) and :47 (PD parse).
- `Pd::diagram()` panics on a parseable but invalid PD (pd.rs:140-142, via 147-154).
- `from_signed` asserts on generator 0 or ≥ strands (pd.rs:240).
- `BoxNodes::new` asserts label ≤ n (boxclose.rs:179); `--box-live` hits this when max(ns) < 2.

**K10 (medium): `--check-library` cost and coverage.**
- The library build for 6_1 (450+93 jobs) has never been measured.
- Cheaper sufficient check: build_wlib plus `lib.no_closure.is_empty()` (swire.rs:211-240) proves coverage. Run the record scan only if that set is non-empty.
- 6_1 visits 5 idempotent types (cover_f31.txt), and no 4-strand box run has ever gone through e_p (C2_4_report.md:29-34, 200). H6 is therefore the first exercise of that path; a failure there is a real finding.

**K11 (medium): fetch with `--repo`.**
- `--depth 1 -b data-gl6` loses ee548d2 once data-gl6 advances. Fetch the SHA explicitly.
- Whether tar entry counts equal the MANIFEST counts is UNVERIFIABLE without listing the tarballs.

**K12 (medium): environment hygiene on the box route.**
- These validation-only switches change the computation: `KRHOM_FORCE_LADDER` (sindec.rs:208-220), `KRHOM_NO_M5I` (636), `KRHOM_DCID_REF` (swire.rs:483-486), `KRHOM_M5I_NO_*`.
- Refuse them or record them in the JSON.

**K13 (low): JSON and m5i_cover details.**
- Hand-written JSON needs string escaping (precedent: libbuild.rs:56 `jstr`).
- m5i_cover needs input lines of the form "name;labels;word;N" (m5i_cover.rs:16-20) and a build with `cargo build --release --example m5i_cover`.

**K14 (low): Q5 inconsistency.**
- The exact path prints 0 for color > N (tests/xengine.rs:307-310).
- After simplification, some 3-strand inputs reach the exact path and print 0, while others are refused. Choose one behaviour.

**K15 (Q1, base branch).**
- The fallback (73cbd7c plus the route_header cherry-pick) is functionally equivalent for replay: the stamp is unchanged with c27 off, and the library hash uses route_header (swire.rs:258-260).
- The unblock base adds about 412 unreviewed lines in leaves.rs. Given the review-and-merge policy, the fallback is the safer default.
- Q2: keeping the base l2222_n4 layer is harmless. c25 restamped the C2-2 files (c25 MANIFEST "provider:" line) and holds 2612 layer-A files against 2575 in the base.

**UNVERIFIABLE**
- README_draft line references (file not found).
- Q7 m4c_b4w MANIFEST (not checked).
- Library build time for 6_1.

### 5. Required corrections
1. C4 step 9: replace `chi_checked == ns.len()` with `== Σ stages[].closed.len()`. Print and compare only the requested mode. Decide whether to expose the reduced closures from unreduced runs, and label that path as new.
2. Classifier: drop the "smapio (replay)" row. When a "not in layer A" error has 0 diagnosed missing jobs, report "present but rejected on load / damaged" and point to the "[sindec] layer A … (rejected)" stderr lines. Add rows for "loaded target rejected" (sindec.rs:1238) and "idempotent required (M5i)", both as Library or damaged.
3. Replace count-equality checks (`locate`, fetch idempotence, `--verify-only`) with per-file hash lists that tolerate extra tmpl files, or make installs read-only. State the Q4 answer: layer A is never written; tmpl may be.
4. Strengthen the preflight so that every c27-layer file (or every m4c-lf layer file) must match the route stamp.
5. Correct the cost figures:
   - F32b factor: about 1× at 3 strands, 1-3× at 4 strands.
   - 8-letter 3-strand knots: 100-1151 s per prime.
   - 4-strand library load: 20-34 s.
   - Cap T3 partner length at 6 or raise the timeout.
   - Re-estimate H11 (about 60-90 min).
6. Make the default `--mem-cap-mb` thread-aware, or fix closure threads at 2.
7. Add `--colored-method singular` to every F1 and H6 command. Set `KRHOM_LIB_DIR` in the updated tests/xengine.rs:650 test. Make `--check-library` require or imply singular.
8. Extend step 0 to cover krhom.rs:43 and :47, pd.rs:140-142 and :240, and boxclose.rs:179.
9. `--check-library`: use `no_closure.is_empty()` first, scan only if needed, and record that 6_1 is the first e_p run at 4 strands.
10. Fix the O5 citation to runs_c2/k4/markov/compare.txt. In O4, list the braid relation and cancelling pair as safe for reduced links.
11. O11: label B, or skip, for (m=4) and (m=3, N≥5). Add the Hopf formula for |n|=2 and the reduced unlink for n=0.
12. Fetch tool: fetch by SHA rather than `--depth 1 -b`.
13. Box route: refuse or record the validation-only env switches. Escape strings in the JSON.
14. Document the m5i_cover input format and its build step. Resolve Q5 consistently across both singular routes.
15. Note in A.3 that the 73cbd7c..194a731 diff also touches Cargo.toml and tests/c27_joint.rs. Prefer the fallback base unless the orchestrator has reviewed C2-7.
