# Changelog

All notable changes to krhom are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] - 2026-11-15

First release of the calculator `krhom` (crate `krhom`, Rust edition 2024, Rust 1.85 or newer).

### Added
- `--checkpoint`: resumable long runs (the box route can be interrupted and continued).
- `--mark`: choose the marked strand of reduced homology.
- `--knot` / `--link`: inputs by table name.
- `--latex`: Poincaré polynomials and tables as LaTeX.
- `--estimate`: the cost of an input (time and memory) without running it.
- Clean refusals: unsupported inputs end with one message and the exit code of the table, never a
  panic.
- `tools/fetch_library.sh` without the `data-gl6` history: a partial fetch
  (`--filter=blob:none --depth 1` into a temporary bare repository) downloads only the blobs the
  libraries need, about 40 MB for `l222_n4` and `l2222_n4` instead of about 185 MB. The checkout is
  not touched. The old method stays as `--full-fetch` and as the automatic fallback when a server
  refuses partial or SHA fetches. Checks (registry and MANIFEST sha256, tar listing, counts) are
  unchanged.
- Release assets: `tools/dev/make_release_assets.sh` writes one reproducible tarball per (library,
  prime) and `SHA256SUMS` into a local directory; `tools/fetch_library.sh --url BASE` installs from
  such assets over HTTPS (curl or wget) or `file://`, checked against `data/release_assets.tsv`
  before extraction and then exactly like a git install.
- Smoke test: `tools/smoke.sh` (seconds without libraries; the library case runs only if
  `l222_n4` is installed) and `cargo test --release --test smoke`.
- Continuous integration: `.github/workflows/ci.yml` (build, tests, dev build, smoke test).
- **One command** `krhom N[,N...] --braid WORD [--colors ...] [--reduced]` that picks the engine,
  runs it modulo two large primes (2^31 − 1 and 4294967291) and prints, for every N, the Poincaré
  polynomial, total rank and Euler characteristic of undeformed Khovanov–Rozansky gl_N homology,
  uncolored or colored by exterior powers ∧^a C^N. A result is printed only if the two primes agree
  (exit 3 otherwise).
- **Inputs:** braid words as signed integers (`"1 -2 1 -2"`) or KnotJob letters (`aBaB`), optional
  `--strands`; PD codes (`--pd`, Vogel's algorithm, then simplification); batch files
  (`--file`, one braid per line, errors reported per line). One color for all components or one
  per component; reduced homology with `--reduced`.
- **Outputs:** text, `--table` (bigraded tables), `--json` (input, route, library, per-prime
  provenance, results; rewritten atomically after every line of a batch) and `--csv`.
- **Routes:**
  - uncolored: the localized Soergel engine with the restriction-tower closure, N ≤ 32, up to about
    5 strands in practice;
  - colored on ≤ 2 strands (after simplification): the exact singular Soergel engine, colors 1..3,
    a(N − a) ≤ 255;
  - uniform color 2 on 3 and 4 strands, N ≤ 4: the localized box engine with the prebuilt certified
    layer-A libraries `l222_n4` and `l2222_n4` (4 strands experimental: only words the library
    covers; `--check-library` answers coverage without the box scan);
  - `--colored-method cube`: the full cube, exact over Q, for tiny inputs.
- **Exit codes** 0–5 with one-line messages (usage/unsupported, certification failure, prime
  disagreement, library missing or not covering the word, resource refusal), and memory guards
  (`--mem-mb`, `--mem-cap-mb`).
- **Libraries:** `tools/fetch_library.sh` installs and verifies the libraries from data commit
  `ee548d2` (registry `data/libraries.tsv`: chunk sha256 against the registry and the MANIFEST,
  tar listing checks, per-file sha256 list in `LIBRARY.txt`, `--verify-only`); krhom re-checks the
  install and the build stamps on every run and re-certifies every loaded map.
- **User tools:** `tools/randcheck.py` (random color-2 braids against fixed relations),
  `tools/batch_check.py` (table checks), `tools/wang.py` (Wang's closed formulas).

### Validation
- Uncolored: acceptance tests (Khovanov trefoil and figure-eight, the KR trefoil formula, the
  reduced trefoil for N = 2..10, cube and exact-closure agreement, the Conway knot at N = 2..10);
  every KnotInfo/LinkInfo record with ≤ 4 strands (1212 knots, 2138 oriented links) against
  Khovanov homology at N = 2 and HOMFLY-PT χ at N = 3..6, 0 failures.
- Colored, 2 strands: full bigraded agreement with Wang's T(2, 2m+1) formulas for ∧² and ∧³;
  colored Hopf links against a closed formula, the full cube and cabling.
- ∧² gl_4: all 101 computed 3-strand knots and 6_1, 7_6, 7_7 on 4 strands pass the two-prime,
  χ, scaling, N = 3 duality and N = 2 checks; 4-strand Markov words reproduce the 2- and
  3-strand values; the random checker passed every recorded case.
- Test suite: 305 tests with default features, 313 with `--features dev`, 0 failed (tests that
  need the libraries are `#[ignore]`d by default).

### Changed
- Repository restructure (no computed result changed; krhom output byte-identical to the previous
  main; the old layout is on branch `archive-pre-restructure`):
  - `src/` grouped into submodules (old module paths re-exported; every stamp literal pinned by
    `tests/compat_paths.rs`);
  - `cargo build --release` builds only `krhom`; developer binaries and examples are behind the
    feature `dev`;
  - user tools in `tools/`, developer scripts in `tools/dev/`; results in `results/`; developer
    documents in `docs/dev/`; the `runs_*` directories removed.

[Unreleased]: https://github.com/Exista-Yang/KRcalculator/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/Exista-Yang/KRcalculator/releases/tag/v0.1.0
