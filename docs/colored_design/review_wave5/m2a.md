# Review wave 5: M6a-2a hardening and box driver (revm2a, colored-m6a c242970)

Scope: `git diff f4c865a HEAD -- src tests build.rs` (f4c865a is the already-reviewed merged dev).
I read the whole diff and did not run cargo or any binary.

## Verdict

No defect found that can produce a wrong Poincaré polynomial or a false certification on the
current code. I found three minor defects: one provenance gap, one resume regression in libbuild
whose test does not catch it, and one source of false refusals across toolchains.

## Lens items I checked and found sound

- **Scaling and χ are mandatory and fatal.**
  - `BoxOpts::new` (swire.rs:1015) is the only constructor in src/, and it sets `scaling: true`.
    No CLI flag or src caller sets it to false.
  - `box_word` (swire.rs:1295-1307) and `box_finish` (boxrun.rs:422-431) call `close_stage_box`.
    That calls `chi_check` and propagates Err with `?`.
  - `scaling_check_frozen` runs before Ok. sbox `die`s (exit 2) before it prints or writes JSON.
  - `chi_ok:true` and `pass:true` are hard-coded in the JSON. That is justified, because the JSON
    is only reachable after the checks have passed.
  - The weights checked are the record's (`sr.weights`, swire.rs:1340). These are the weights the
    BoxAccs were built from.
- **Checkpoints.**
  - `BoxAcc::save` (boxclose.rs:877) and `write_atomic` (boxrun.rs:148) write a temp file, fsync,
    rename, then fsync the directory. The temp names differ per file and per pid.
  - `load` checks the payload length and hash, and the grid. `load_expect` then compares the full
    header, `bind` included.
  - `merge` refuses any header difference and any bitmap overlap.
  - Accumulation is exact F_p arithmetic, so order does not matter, and resume is byte-identical.
- **Kill between two file writes.** `pending` is "any accumulator lacks the point".
  `accumulate_point(.., skip_done=true)` skips files that already have it. No point is ever added
  twice to a file.
- **Shards.**
  - The governing part of `run.hdr` must equal the finishing job's header (boxrun.rs:387).
  - Shards must be distinct with the same K (boxrun.rs:390). A missing shard fails the down-set
    check in `finish`.
  - A foreign checkpoint is refused by `load` or `load_expect`.
  - A copied shard is refused by the overlap check.
  - If two processes write one directory, the last whole-file writer wins. Each file stays
    self-consistent, because a value depends only on the bound header and the point.
- **Pivot-value replay.**
  - `apply_pivot_inner` pushes exactly one fingerprint per applied pivot, after `scalar_of`
    (sengine.rs:1567). Every record-mode caller propagates its Err.
  - Replay compares `piv_vals.last()` with `r.values[k]` (sengine.rs:1243), and the round checks
    `values.len() == pivots.len()`.
  - The invariant holds: an evaluated component a→b of equal (p, t) that is a genuine degree-0 map
    is c·id with c in k (pruning included), so it is the same c at every point and at λz.
  - `Fp` is canonical (`x % P`), so equal values hash equally. Within one process I see no
    false refusal at valid points.
- **libbuild replay-only.**
  - `has_hom` accepts a cache file only on a valid header. A corrupt body would be recomputed,
    and the child then reports `tmpl_computed > 0`, which becomes a `Check` error. That is double
    coverage.
  - A child commit or stamps mismatch becomes outcome `Commit`. The stamp format is the same on
    both sides (`layer_a_stamps` and the parent's `format!`), so this is not a false refusal.
- **xscan.** `certified = check_d2 && bp_ok`, exit 3 on a base-point FAIL after the JSON is
  written, and the peak is the process peak. That is consistent with the documented definition.

## Findings

### F1 (minor, provenance): a shard from another binary merges silently, and finish names only its own commit
- **Where:** boxrun.rs:253-255, 376-435, 479-480; libbuild.rs:572.
- **Defect:**
  - The checkpoint header binds manual stamps (`SWIRE_STAMP`, `SENGINE_STAMP`, `BOX_CODE_STAMP`),
    the trace hash and the record-point values. It does not bind the commit.
  - `commits.log` is write-only: nothing reads it, and `box_finish` does not check it.
  - `finish_json` records only `build_commit()` of the finishing binary.
- **Consequence:** shard runs across sessions are the purpose of D14, and dev moves between
  sessions. Two kinds of shard pass every header check and merge silently:
  - a shard scanned by an older or `-dirty` binary with unchanged stamps;
  - a shard whose code change only affects box-point values (not the record point's trace or
    values).
  - The result JSON then attributes the accumulators to the finishing commit.
  - Report §5.2 says such a shard "is caught by … the commit log", but no code consults the log.
  - libbuild.rs:572 says resume requires "the same binary commit". The governing key has no
    commit, so a rebuilt binary resumes from the old binary's certified lines.
- **Scenario:** shard 0 is scanned at commit A. A later commit B changes box-point evaluation
  without a stamp bump, and shard 1 is scanned at B. `finish` at B merges both. If χ happens to
  agree, the JSON says `"commit":"B"` for a mixed result.
- **Fix:**
  - Put the build commit into `run.hdr` or the BoxAcc header, or have `box_finish` parse every
    `commits.log`.
  - Refuse a mix of commits (or any `-dirty` commit) without an explicit flag, and list every
    shard's commits in the finish JSON.
  - Correct the libbuild comment, or add the commit to the governing key.

### F2 (minor, resume / test): the libbuild governing key changes once the cache directory exists, so resume is silently disabled
- **Where:** libbuild.rs:554 and 563-564; tests/libbuild.rs:76-99.
- **Defect:**
  - `cache_abs = canonicalize(cache).unwrap_or(cache)`, and canonicalize fails on a path that
    does not exist yet.
  - Only the children create the cache (sindec.rs:513, tmpl.rs:311).
  - So the first run's header carries the raw `--cache` string. Every later run carries the
    canonical path. `cur_ok` (libbuild.rs:581) is never true for the first run's lines.
- **Consequence:**
  - With a relative `--cache` (e.g. `--cache C`, as in review_wave3/m5_lib.md:25) or a symlinked
    path on a fresh cache, a restart after a timeout or kill ignores every certified job. All of
    them are rerun.
  - This wastes compute and does not give a wrong result.
- **Test gap:** `resume_governing_header_and_replay_only` never checks for a positive resume. It
  only checks `(100.0%)`, which a full rerun also prints.
  - On macOS, `temp_dir()` is /var/folders/…, a symlink to /private/var/…. The second `run`
    resumes nothing.
  - There the "pair 99999 is not counted" check passes only because the whole header is foreign,
    not because of the set filter.
  - On the Linux cloud (/tmp, an absolute non-symlinked path) it happens to work.
- **Fix:**
  - `create_dir_all(cache)` before canonicalizing (or canonicalize the parent and join the
    name).
  - Assert `jobs_skipped_done == jobs_total` on the second run.

### F3 (minor, false refusal only): persisted hashes depend on std `DefaultHasher`
- **Where:** sengine.rs:574-578 (`value_fingerprint`), used by `trace_hash` (swire.rs:1037) and
  `values_hash` (swire.rs:1055).
- **Defect:** these feed `run.hdr` and the BoxAcc binding (`trace_hash`, `record_hash`). Rust does
  not guarantee the `DefaultHasher` algorithm across releases.
- **Consequence:**
  - A shard resumed or finished by a binary built with another rustc may have every hash differ.
    Resume and merge are then refused, and the whole box has to be rescanned.
  - Within one process, the pivot-value comparison is consistent, so valid points are never
    refused.
- **Fix:** fingerprint with the crate's own `Fnv` over `to_u64()`, as boxclose already does.
  Bump `SENGINE_STAMP`.

## Notes (not defects)

- The scaling check certifies only entries that are nonzero at its 2 points. In the σ^5 demo,
  3 808 entries were zero at both points, against 3 724 checked.
  - A weight declared too high on such an entry would not be caught. The residual check catches
    only weights declared too low.
  - This is the D14 design, and the report states it. More points would raise coverage.
- `box_run` scans every word of a job at the union of needed points. That is extra cost for short
  words in a GP3 batch, not a correctness issue.
