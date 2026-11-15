# M6a-2a report: hardening of the localized path and the checkpointed box driver (branch colored-m6a)

Scope: the brief of this session. It covers step (0), the merged dev; step (1), the hardening items from review
waves 3 and 4 and from review_2s_m5; and step (2), the resumable, shardable box driver of SPEC D14 / M6a-2a,
demonstrated on (3,3) σ^5 at N = 6. GP3 and T0 from the SPEC's M6a-2a section were not part of this brief and
were not run (§5).

## 0. Verdict

- **Merged dev** (f4c865a): the full suite passed before any change (§1).
- **Hardening**: every listed item is done, and each has a test that fails without the change.
  - The one partial item is pivot values. They are compared exactly where the frozen structure implies
    equality (every swire scan). Synthetic complexes keep the old structure-only replay, since their cone
    scalars depend on the point by design.
- **Box driver**: `box_run` / `box_finish` / `box_point` (src/boxrun.rs) and `sbox box|finish|point` are in place.
- **Demonstration**: (3,3) σ^5 at N = 6, F31, seed 101, u+r (399 box points). The run was split into 2 shards, and
  shard 0 was SIGKILLed after 150 s and resumed. The merged result:
  - equals the M6a-1c run (40/40 closures identical);
  - equals the m2s references and Wang (0 differences);
  - has accumulators byte-identical (digests) to an uninterrupted single-shard run.

## 1. Step (0): the merged code

`git merge --no-ff origin/dev` (f4c865a) into colored-m6a gave 8c81d4b.
`KRHOM_THREADS=4 cargo test --release` exited 0. Every test binary passed (passed / ignored):

| binary | result | binary | result |
|---|---|---|---|
| lib unit tests | 25 / 1 | tests/libbuild.rs | 2 / 0 |
| tests/acceptance.rs | 25 / 0 | tests/m4a.rs | 25 / 5 |
| tests/boxclose.rs | 13 / 2 | tests/ntrunc.rs | 12 / 1 |
| tests/boxclose_mem.rs | 1 / 1 | tests/pd.rs | 11 / 0 |
| tests/cclose.rs | 24 / 2 | tests/schur.rs | 26 / 2 |
| tests/dcoset.rs | 6 / 0 | tests/sengine.rs | 11 / 1 |
| tests/hom_lazy.rs | 3 / 0 | tests/sindec.rs | 15 / 3 |
| tests/homprov.rs | 18 / 3 | tests/swire.rs | 4 / 1 |
| tests/knotcov.rs | 2 / 0 | tests/xengine.rs | 19 / 3 |
| tests/leafgen.rs | 9 / 6 | tests/xstream.rs | 4 / 3 |
| tests/leaves.rs | 8 / 0 | bin unit tests | 0 tests each |

There were no failures, so nothing had to be fixed before the changes.

## 2. Hardening (step 1)

| item | done | where | test (fails without the change) |
|---|---|---|---|
| D14 scaling check mandatory in `box_word` / sbox; opt-out only for tests | yes | `BoxOpts::scaling` (default true via `BoxOpts::new`), `box_word` runs `scaling_check_frozen` on its own record trace; `box_finish` too; sbox has no opt-out (`--scaling` accepted, ignored); a failure is FATAL (exit 2, no results) | `boxrun::scaling_check_catches_a_misdeclared_weight` (wired path: an entry declared w+1 and w−1 → Err; review sengine F6); `swire::box_loop_against_exact_22_32` asserts `out.scaling` present with checked > 0 |
| χ (decat_chi) mismatch FATAL | yes | `swire::chi_check` (Err), called in `close_stage_box` (box_word and box_finish); sbox exits 2 | `boxrun::chi_mismatch_is_an_error` |
| compare tools use chi_ok | yes | tools/g2b_compare.py and tools/stretch_compare.py fail on a result without `chi_ok: true` or a missing / failed scaling record | rerun on all committed runs: 0 differences, 0 failures |
| "zero at both points" reported apart from entries declared 0 | yes | `ScalingReport::declared_zero` | (reported in the JSON) |
| replay compares pivot values | yes, where the structure implies equality | `Round::values` (a fingerprint per pivot scalar) is recorded at point 0. With `EngOpts::check_pivot_values`, a replay refuses a differing scalar (`PointDependent`). swire's `scan_point` forces it on: a degree-0 component between equal types of genuine maps is c·id, the same at every point, including λz. Synthetic complexes keep it off. `SENGINE_STAMP = sengine-v2`. | `sengine::replay_compares_pivot_values` (refused at another point, passes at the same point and with the option off); all swire / boxrun tests and the σ^5 demo replay with it on |
| BoxAcc / checkpoint header binds library, skeleton, stage, engine version | yes | `BoxHeader.bind: BoxBinding {stage, skeleton_hash, engine_stamp, trace_hash, record_hash}`; `BoxAcc::with_binding`, `load_expect` (full header comparison), `header_diff`; checkpoint format v2 (`BOX_CODE_STAMP = boxclose-v2`); `WLib.hash` now also includes the layer-A stamps (route form). The record hash (entry values at the record point) binds the library *content* in its basis. | `boxrun::box_acc_binding_is_part_of_the_header` (each binding field separates equal weight lists; a stale `load_expect` is refused) |
| libbuild: resume checks the governing header | yes | header field `governing` = mode, λ, N, field, set, canonical cache path, seed, provider, stamp-route, stamps. An attempt line counts for a resume only under a header with the same key. | `libbuild::resume_governing_header_and_replay_only` (a replay into the run's own file resumes nothing: `jobs_skipped_done: 0`) |
| libbuild: mode and seed in the default output name | yes | `libbuild_{λ}_{set}_f{field}_{run|replay}.jsonl`; sample set names carry `-s{seed}` | same test (the default name ends in `_run.jsonl`) |
| libbuild: certified counts only over the current set | yes | `done ⊆ set` and `n_ok` counted over `set_keys` | same test (an ok line for `pair 99999` keeps 100.0 %) |
| libbuild: replay never runs the exact solver | yes | `smapio::set_replay_only`: in replay, `import` refuses an `H` line whose Hom basis is in neither memory nor the template cache; the child also reports `tmpl_computed` and fails on > 0 | same test (template cache deleted → the replay FAILs instead of recomputing) |
| libbuild: commit of the binary actually run | yes | build.rs embeds `KRHOM_BUILD_COMMIT` (with `-dirty`). The child prints `commit` and `stamps` in RESULT, and the parent marks a mismatch with its own as outcome `Commit`. The header records the build commit, not `git rev-parse` at run time. | same test (every attempt line carries the header's commit) |
| libbuild `cover`: only attempt lines of the right λ and N | yes (extra) | `run_status(path, λ, n)` | — |
| xscan JSON: check_d2, engine options, commit | yes | `check_d2`, `certified` (= check_d2 ∧ base point), `rung_transport`, `rungwise`, `seed`, `threads`, `commit`, `basepoint_ok` | `boxrun::xscan_json_provenance` |
| xscan: a base-point FAIL changes the exit status | yes | exit 3 after writing the JSON (`basepoint_ok: false`) | not triggerable on correct code; the path is the `bp_ok` flag that `certified` also uses |
| xscan: peak RSS is the process peak | yes | max over the stages' `rss_peak_kb` and the final VmHWM; `peak_rss_scope: "process"` | same test (field present) |
| box_word memory (review box note, sengine F5) | yes (extra) | points are scanned in batches (`BoxOpts::batch`, default 16) and accumulated as each batch finishes | — |

## 3. The box driver (step 2): src/boxrun.rs, `sbox box|finish|point`

- **Job.** λ, words, closed stages per word, box N, node seed, prune_n, modes, N' list, engine options. One
  library serves the whole batch of words, and each point is scanned for every word.
- **Run directory.** `run.hdr` is the governing header: the job, the library hash, the engine stamp, and a line
  per (word, stage) with the skeleton, trace and record hashes, plus `shard=k/K`. Each (word, stage, mode) has
  its own `w{i}_s{stage}_{u|r0}.acc` (a BoxAcc v2 checkpoint, written atomically: temp file, fsync, rename, dir
  fsync). `commits.log` lists every binary that ran there.
- **`box_run`.**
  1. Record scans. The governing header must equal `run.hdr` if it exists (resume), or is written if it does
     not.
  2. Every existing checkpoint is loaded with `load_expect` against its full expected header.
  3. Shard k/K takes the needed points (union over the words) with index ≡ k mod K.
  4. The generator-only determinant pass runs at every point of the shard (identity-coset in- and
     out-generator matrices nonsingular for every closable type) before any scan.
  5. Pending points are those any accumulator still needs. They are scanned in batches (replay of the frozen
     trace, pivot values included) and accumulated per file, skipping points that file's bitmap already has.
     That makes a kill between two file writes harmless.
  6. Checkpoints are written every `--ckpt-secs` and at the end.
  7. `RunOpts::stop_after` simulates a kill: work since the last checkpoint is lost.
- **`box_finish`.**
  1. Each shard directory's `run.hdr` must equal the job's header, and the shards must be distinct with the same
     K.
  2. Each accumulator is merged from every shard, after `load_expect` and `merge`, which refuses overlap and any
     header difference.
  3. The determinant pass runs at all needed points.
  4. Per stage: `finish` (complete down-set, residuals), `BoxClosureData`, Freivalds `D2Cert`s, `close_box`, and χ
     (fatal).
  5. Per word: the D14 scaling check (mandatory).
  6. The JSON records the commit, engine, check_d2, Freivalds k, shards, determinant pass, scaling and every
     accumulator digest.
- **`box_point`.** Scans one box point and prints a value hash per stage, for spot checks.

**Tests** (`tests/boxrun.rs::kill_resume_shards_and_refusals`, (2,2) σ^3 + σ^2 in one job, box N = 5, both modes,
1.4 s):
- **Uninterrupted run.** It equals `box_word` stage by stage.
- **Kill and resume.** The run is killed twice: after 7 points with a checkpoint after every batch, then after 9
  more with nothing saved. It is then resumed to completion, and a finished shard rescans nothing. The result
  has identical accumulator digests and a byte-identical finish JSON. The 3 points lost to the first kill are
  redone.
- **Shards.** Shards 0/2 + 1/2 equal the single run, digests included. One shard alone is refused with an
  incomplete down-set.
- **Refusals:**
  - the same shard given twice;
  - a shard with another node seed;
  - a shard with another pivot order (engine stamp);
  - a resume into another job's directory;
  - a resume as another shard;
  - a forged directory: the right `run.hdr` with a foreign checkpoint, refused by the checkpoint header;
  - a truncated checkpoint.

## 4. Demonstration: (3,3) σ^5, N = 6, sharded and killed (runs_m6a2a/)

The script is runs_m6a2a/demo_s5.sh, with the M6a-1c G2B parameters: `--n 6 --ns 3,4,5,6 --field 31 --seed 101 --modes u,r`.
- **Library.** Built once into the cache: 80 provider calls, 28.7 s. Every later invocation loads layer A (80
  loads, 0 calls), hash cabf041aaa23a517.
- **Shard 1/2.** 199 points, 2.73 s per point on 2 threads, 279 s.
- **Shard 0/2.** SIGKILLed by `timeout -s KILL 150` (exit 137) after 80 points had been checkpointed. The resume
  found 120 pending, scanned them (2.70 s per point) and completed.
- **Determinant pass.** 8 generator matrices per point (4 types × in/out): 1 600 and 1 592 checks on the two
  shards; 3 192 at all 399 points in finish.
- **Finish** (40 s):
  - 40 closures (stages 1..5, N' = 3..6, unreduced and reduced), all with χ = decat_chi;
  - scaling: 3 724 entries checked, 3 808 zero at both points, 0 declared 0;
  - every closure equals runs_m6a1c/g2b/s5_f31_seed101.json (the M6a-1c result: 40/40 identical);
  - `tools/g2b_compare.py` against docs/dev/colored_design/gl6/m2s/sigma9_33_F31.json and Wang: 0 differences
    (stage 5 at N = 6: rank 396 unreduced, 125 reduced).
- **Uninterrupted control.** One shard, 4 threads, 280 s: every accumulator digest is identical to the merged
  sharded and killed run, and the finish JSON is byte-identical apart from the `shards` field (and the commit
  of the regenerated file).

The files are runs_m6a2a/s5/{shard0,shard1}/ (run headers, checkpoints and commit logs, 4.7 MB), *.out / *.err,
finish.json and finish_single.json.

## 5. Open problems

1. **GP3 and T0** (SPEC M6a-2a gates) were not part of this brief and are not run. The driver can now carry them:
   `--words` batches a knot list per prune level.
2. **The scaling check uses 2 points per word** in finish. It is one run per word at the end, not per shard.
   A shard that ran a buggy binary is caught by the engine stamp, the trace and record hashes, the commit log,
   and the scaling check of the finishing binary. It is not caught per shard.
3. **The pivot-value check is a fingerprint** (SipHash of the field element), not the value. Equal values always
   pass. A differing value is missed only on a 64-bit collision.
4. **The record hash binds the library content only through the closure entries' values at the record point.**
   Two libraries that differ only in maps that never reach a closure entry at the record point are not
   separated. Their results cannot differ either, unless the trace differs, and the trace hash binds that.
5. **`box_run` rewrites every checkpoint file of the shard on each save.** At σ^5 that is 20 files and 0.2 MB.
   At 4_1 scale (≈ 2.6e5 entries × L(top) coefficients) a save is hundreds of MB, so `--ckpt-secs` should be
   minutes. Incremental (append-only) checkpoints are possible but not implemented.
6. **build.rs reruns on every commit** (`.git/refs/heads`) and on src changes, so the library recompiles once
   after each commit.
7. **xscan's `--chi-only` panic handling** (2s_close F3) and the Wang-block presence checks (2s_close F1,
   2s_core m1) are untouched. They were not in this brief.
8. **Unchanged from M6a-1c:** the open items for elimination speed at depth (§7.1-3 there) and for σ^8/σ^9
   unreduced. The sharded driver now allows those runs to span sessions.

## 6. Commits (colored-m6a)

- 8c81d4b merge origin/dev (f4c865a)
- 5d5ec83 pivot values, BoxAcc binding, mandatory scaling, fatal χ, determinant pass, boxrun driver
- e863d6b libbuild governing-header resume, names, counts, replay-only import, binary commit; xscan provenance
- 7c638e0 tests (driver, binding, scaling negative, χ, pivot values); compare tools
- 4dbf8ba libbuild and xscan tests
- 9b35adb σ^5 sharded / killed demonstration; finish JSON fix
- 01e304b, final: this report

**Final full suite** (on 9b35adb; the later commits change only this report): `KRHOM_THREADS=4 cargo test --release`
exited 0 with 0 failures. New or changed counts against §1:
- tests/boxrun.rs: 5 passed (new);
- tests/libbuild.rs: 3 passed (was 2);
- tests/sengine.rs: 12 passed (was 11);
- every other binary: as in §1.
