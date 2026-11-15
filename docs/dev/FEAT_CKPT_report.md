# F-CKPT: `krhom --checkpoint DIR` (resumable box route)

Base: origin/main 0e2ea24. Branch: `feat-checkpoint`. User documentation:
[docs/CALCULATOR.md §9](../CALCULATOR.md#9-resumable-runs---checkpoint-dir).

## Design

- No new checkpoint format: the box route goes through `colored::boxrun`'s checkpointed driver
  (the machinery of `sbox box|finish`) instead of the in-memory `swire::box_word`.
  `boxrun` gains `box_run_with` / `box_finish_with`, which take the record scans already done
  (krhom runs the record scan once per prime for both); `box_run` / `box_finish` delegate to them
  and are otherwise unchanged.
- `src/calc_ckpt.rs` (new, isolated): `configure` (the flag; refuses `--file` batches and
  `--check-library`), `manifest` / `manifest_diff` / `prepare` (per-prime manifest, refusal with
  exit 1 listing the differing keys), `box_word_ckpt` (box_run_with + box_finish_with returning a
  `BoxWordOut`, so `run_prime`'s post-checks, retries and reporting are shared).
- `src/calc.rs::run_prime`: checks/creates the manifest first (before the library build), then
  calls `box_word_ckpt` instead of `box_word` when a checkpoint directory is configured. Seed
  retries and the template-recompute retry work as before (`seed<S>/` per node seed).
- `src/bin/krhom.rs`: the `--checkpoint` value option (usage text) and one call
  `krhom::calc_ckpt::configure(...)` at the end of `parse_args`.
- Layout: `DIR/<f31|f32b>/manifest.txt`, `DIR/<f31|f32b>/seed<S>/{run.hdr, w0_s<stage>_<mode>.acc, commits.log}`.
- Manifest keys: format, field, braid (used, after simplification), strands, labels, N (the N
  computed, i.e. those >= the color), mode, seed, box_n, library, library_field, data_commit,
  provider, krhom_commit, boxrun stamp. boxrun's `run.hdr` and acc headers add library hash,
  engine stamp, record trace / skeleton bindings (refused on mismatch by boxrun itself).
- Certification identical to `box_word`: same `BoxOpts::new` options, generator-only determinant
  pass (before the scan of the run, again over all needed points in finish), frozen replay,
  Freivalds via `close_stage_box`, χ = decat_chi (fatal), D14 scaling (mandatory), then
  `run_prime`'s post-checks and the two-prime comparison in krhom.
- Parallelism: points are scanned `KRHOM_THREADS` at a time inside boxrun (`par_map` over
  batches of 16), as in the in-memory route; no `--jobs` flag. Multi-machine: `sbox box --shard`
  (documented).
- Environment: `KRHOM_CKPT_SECS` (save interval, default 60 s), test hook
  `KRHOM_CKPT_STOP_AFTER=M` (simulated kill, exit 75).

## Validation

- `cargo build --release`: ok.
- `KRHOM_THREADS=4 cargo test --release` (default suite): **308 passed, 0 failed**, 42 ignored
  (base: 305; new: `calc_ckpt::tests::manifest_mismatch_is_refused_with_exit_1`,
  `calc_ckpt::tests::configure_refusals`, `tests/checkpoint.rs::checkpoint_refusals_exit_1`).
- `cargo test --release --test checkpoint -- --include-ignored` with l222_n4 / l2222_n4 installed
  (`tools/fetch_library.sh l222_n4 l2222_n4`): 2 passed (77 s). `lib_checkpoint_resume_equals_plain_run`:
  `"1 -2 1 -2" --colors 2 --reduced`, stopped by the hook after 20 of 36 F31 points (exit 75, 16
  saved), resumed (16 already done, 20 scanned), stdout equal to the plain run except `time`;
  rerun of the finished directory scans 0 points and prints the same; changed mode, braid, N and
  seed are each refused with exit 1.
- Manual 3-strand run (same word): plain and interrupted+resumed outputs identical (except
  `time`); `--checkpoint` on the same DIR without `--reduced` → exit 1,
  `mode: directory has "reduced (strand at bottom position 0 marked)", this run has "unreduced"`.
- 4-strand check, real SIGKILL: `krhom 2,3,4 --braid "1 2 3" --colors 2 --no-simplify --reduced
  --checkpoint D4` started with `KRHOM_THREADS=1 KRHOM_CKPT_SECS=0`, killed (`-9`) right after
  the first batch checkpoint of F31; resumed with `KRHOM_THREADS=4`:
  `[F_2147483647] checkpoint D4/f31/seed201: 216 box points, 16 already done (resumed), 200 scanned now`,
  F32b fresh (216 scanned), primes agree; stdout identical to the plain run except `time`
  (N = 2, 3, 4: `t^0q^0`). Determinant pass 11664 checks at 216 points, scaling 1080 entries,
  χ for 3 closures — the same numbers as the plain run.

## Notes / not done

- stderr differs from the plain run by the checkpoint lines and timings; `close` in the stderr
  summary includes the scaling check on the checkpoint route (boxrun times them together).
- A finished prime is not cached: a rerun after F31 finished but F32b was interrupted repeats
  F31's record scan, determinant pass, closing and checks (no point scans).
- No lock against two processes on one DIR (documented: do not).
- When feat-cli's `--mark` lands, the marked strand must be added to the manifest (`mode=`
  line in `calc_ckpt::manifest`; today reduced always marks bottom position 0).
- `--json` records `KRHOM_CKPT_*` through the existing `env` field only; there is no checkpoint
  field in the per-prime JSON.
