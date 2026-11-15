# F-DIST: distribution, smoke test, CI, changelog (branch `feat-dist`)

Base: origin/main 0e2ea24. No Rust source was changed except the new test `tests/smoke.rs`;
README.md was not edited (see "For the orchestrator" below).

## 1. `tools/fetch_library.sh` without the data-gl6 history

Transport, in order:
1. the data commit is already in `--git-dir` (default: this checkout): read there, no download;
2. **partial fetch (new default):** a temporary bare repository (gc and maintenance off) gets a
   remote `krhom` = `--repo URL` or the checkout's origin, and runs
   `git fetch --no-tags --filter=blob:none --depth 1 krhom <sha>` (trees only, about 10 KB).
   `git show <sha>:<path>` then lazily fetches exactly the blobs a registry row needs: its chunks,
   their `MANIFEST.txt` and the costs file. The user's checkout is not modified (no promisor
   config, no shallow entry, no objects). Nothing in the script asks for blob sizes
   (`git ls-tree -l` would fetch every blob of the commit; found while testing);
3. **fallback / `--full-fetch`:** the old method, unchanged (by SHA with history into the checkout,
   or into a temporary bare repository with `--repo`; last resort the `data-gl6` branch). It is used
   automatically when the partial fetch fails. If a server ignores the filter (warning
   "filtering not recognized"), the depth-1 snapshot is used and a note is printed.

All checks are unchanged and apply to every transport: chunk sha256 = registry = MANIFEST, tar
listing (regular files `layer_a/a-<16hex>.txt`, `tmpl/tmpl-<16hex>.txt` only), counts = MANIFEST,
`LIBRARY.txt` with per-file sha256, `--verify-only`, "already installed", `--force` with the old copy
moved aside. All options are kept; new: `--url BASE`, `--full-fetch`. The script prints the bytes
downloaded at the end.

### Measurements (this 4-core cloud machine, fresh empty `--root`, both libraries, both primes)

| method | downloaded | wall time (fetch + verify + install) |
|---|---|---|
| old (`--full-fetch --repo URL`, = previous behaviour) | 181 176 KiB git objects (185.5 MB; data commit + 13 ancestors) | 53 s |
| new partial fetch (`tools/fetch_library.sh --root EMPTY l222_n4 l2222_n4`) | **40 456 KiB** git objects | 54–56 s |
| `--url` release assets (local http.server) | 40 160 KiB | 45 s |

The network here is fast (the old fetch alone took 8.6 s), so the wall time is dominated by the
sha256 checks of about 12 000 extracted files (≈ 45 s); on a slow link the 4.5× smaller download
dominates. A single library: `--force l222_n4/f31` downloaded 3 480 KiB.

### Fallback tests
- SHA fetch refused (protocol v0 via `GIT_CONFIG_COUNT`, local mirror whose `data-gl6` tip is a
  child of ee548d2): partial fetch refused → full fetch by SHA refused → `data-gl6` branch fetch →
  installed, `LIBRARY.txt` identical.
- Filter ignored (`file://` mirror without `uploadpack.allowFilter`): note printed, depth-1 snapshot
  (149 820 KiB) used, installed, identical.
- `--full-fetch` against GitHub: installed, identical.
- Install equivalence: `LIBRARY.txt` (all fields and the per-file sha256 lists, except the
  `installed:`/`fetched_by:` lines) is identical for all four (library, prime) rows between the
  partial fetch, the old method, and the `--url` install.
- Second run: `already installed` ×4; `--verify-only`: `verified` ×4; `--url` with `--full-fetch`:
  exit 1 (usage).

## 2. Release assets

- `tools/dev/make_release_assets.sh [--git-dir DIR | --repo URL] [--check] OUTDIR` writes into a
  local directory (refuses a non-empty one) `krhom-lib-<id>-<field>-ee548d2.tar.gz` for the 4
  registry rows, `SHA256SUMS` and `release_assets.tsv`. Each asset holds, at their data-commit
  paths, the row's chunks, their MANIFESTs and the costs file; chunks are checked against the
  registry and MANIFEST before packing. Archives are reproducible (GNU tar `--sort=name
  --no-recursion --owner=0 --group=0 --numeric-owner --mode=0644 --mtime=@<commit time>`,
  `gzip -n -9`): built from GitHub (partial fetch) and from a local mirror, the SHA256SUMS are
  identical, and `--check` matches `data/release_assets.tsv`. Nothing is uploaded; no release or
  tag was created.
- `data/release_assets.tsv` (new; `data/libraries.tsv` and its parser in `src/calc.rs` unchanged):
  `id field data_commit asset sha256 bytes`.

| asset | bytes |
|---|---|
| krhom-lib-l222_n4-f31-ee548d2.tar.gz | 3 502 935 |
| krhom-lib-l222_n4-f32b-ee548d2.tar.gz | 3 511 424 |
| krhom-lib-l2222_n4-f31-ee548d2.tar.gz | 16 058 461 |
| krhom-lib-l2222_n4-f32b-ee548d2.tar.gz | 18 051 277 |

- `fetch_library.sh --url BASE`: curl (`--proto =https,http,file`, 3 retries) or wget; sha256 against
  `data/release_assets.tsv` **before** extraction; listing: only regular files/directories, no `..`,
  exactly the expected set of paths; extracted with `--no-same-owner` into the staging directory;
  then the normal install path.
- Tests: `python3 -m http.server` on 127.0.0.1 (all 4 rows installed, identical to git) and
  `file://` (l222_n4/f31). Rejected as expected: an asset with one appended byte (sha256 mismatch),
  an asset with an extra file `c2/extra.txt` whose sha256 was put in a patched tsv ("unexpected tar
  entry"), an unreachable server ("cannot download"), `--url ftp://...` (exit 1). No staging
  directory was left behind.
- When the repository is public and a release exists, the URL will be
  `https://github.com/Exista-Yang/KRcalculator/releases/download/<tag>` (not created here).

## 3. Smoke test

`tools/smoke.sh [--bin PATH] [--lib-dir DIR] [--no-libs | --require-libs]`; uses
`target/release/krhom` (builds only `--bin krhom` if missing). Exit 2 on any mismatch or failure.

Without libraries (`KRHOM_LIB_DIR=<empty>`), 1 s:
```
ok    4_1 uncolored, N=2,3: ranks 6 11 (0 s)
ok    trefoil color 2, N=2,3,4: ranks 1 7 22 (0 s)
ok    trefoil color 2 reduced, N=2,3,4: ranks 1 3 9 (1 s)
skip  4_1 color 2 reduced on 3 strands: library l222_n4 not installed under .../empty (tools/fetch_library.sh l222_n4 installs it)
smoke: all checks passed (1 s)
```
With libraries (`--lib-dir` = the fresh partial-fetch install), 22 s:
```
ok    4_1 color 2 reduced, 3 strands (library l222_n4): ranks 1 5 25 (22 s)
smoke: all checks passed (22 s)
```
`--require-libs` with no library: FAIL, exit 2. A fake binary printing wrong ranks: 3 × FAIL, exit 2.
The library check is `fetch_library.sh --verify-only l222_n4` under the chosen root.

`tests/smoke.rs` (`cargo test --release --test smoke`): the three library-free cases, 0.11 s.

## 4. CI

`.github/workflows/ci.yml` on push and pull_request: ubuntu-latest, `dtolnay/rust-toolchain@stable`
with a check rustc ≥ 1.85, `Swatinem/rust-cache@v2`, `KRHOM_THREADS=4`,
`CARGO_PROFILE_RELEASE_DEBUG=false`; steps `cargo build --release`, `cargo test --release` (default
features), `cargo build --release --features dev`, `tools/smoke.sh --no-libs`; job timeout 40 min.

Slow tests: the six tests that take 34–95 s alone are skipped in CI with `--skip` (no test file was
edited; they still run with a plain `cargo test --release`):
`per_predecessor_final_pair_1444` (c27_joint, 95 s), `rex_class_pairs_small` (rexplan, 69 s),
`b2_type_pairs_22_212_through_chain` (homprov, 39 s), `exact2_memory_cap_and_shape_refusals`
(calc, 35 s), `leaf_pairs_through_sindec_212` (b4w, 35 s), `m5i_2402_prefix_propagation` (m5i, 34 s).

Timing of the CI steps here (4 cores): test run with the skips 198 s (302 passed, 0 failed, 41
ignored, 6 filtered) instead of 508 s; dev build after the default build 3 min 52 s; smoke < 1 s;
first release build plus test compile about 4.5 min (with debug info). Estimated cold CI run on a
4-vCPU runner: about 12 min (less with the cache and without debug info); a 2-vCPU runner (private
repositories) may take about 20 min, hence the 40-min timeout.

YAML: `python3 yaml.safe_load` ok; `actionlint 1.7.12` (pip `actionlint-py`): no findings. The
workflow has not run on GitHub yet (it will on the first push of this branch, if Actions are
enabled for the repository).

## 5. CHANGELOG.md, version

`CHANGELOG.md` (Keep a Changelog): `[0.1.0] - 2026-11-15` (features, inputs, outputs, routes, exit
codes, libraries, tools, validation, the restructure) and `[Unreleased]` (checkpoint, --mark,
--knot/--link, --latex, --estimate, clean refusals, library download without git history, release
assets, smoke test, CI). `Cargo.toml` already has `version = "0.1.0"`. No tag.

## Validation summary

- `cargo build --release`: ok (4 pre-existing warnings in the lib).
- `KRHOM_THREADS=4 cargo test --release`: **308 passed, 0 failed, 41 ignored** (305 before + 3 smoke),
  785 s including compilation.
- `tools/fetch_library.sh l222_n4 l2222_n4` in a fresh empty `--root`: partial fetch, 40 456 KiB,
  54–56 s, all 4 rows installed and verified.
- URL mode against `python3 -m http.server` and `file://`: installed, identical; tamper cases rejected.
- `tools/smoke.sh`: passes without libraries (1 s, library case skipped) and with libraries (22 s).

## For the orchestrator (README.md, not edited here)

- README line 20: "run once, about 2 min" → about 1 min, 40 MB download.
- README lines 24–27 ("Fetching libraries"): the download is now about 40 MB (partial fetch of
  data commit `ee548d2`, only the needed blobs; the checkout is not modified); `--full-fetch` gives
  the old 185 MB method, `--url BASE` installs from release assets; requirements add curl or wget
  for `--url`.
- Mention `tools/smoke.sh`, `CHANGELOG.md` and CI in the repository map / developer tools.
