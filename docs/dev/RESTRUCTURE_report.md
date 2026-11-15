# Repository restructure (branch `restructure`)

Base: origin/main 00f47df. Goal: a clean, logical layout without changing any computed result.
A snapshot of the old layout is on branch `archive-pre-restructure`.

## Commits

1. Results: colored gl_4 data -> `results/colored_gl4/` (+ `results/README.md`); sharded box driver -> `tools/dev/`
2. Remove the `runs_*` directories (kept on `archive-pre-restructure`)
3. src: group the flat modules into submodules, old paths re-exported
4. Binaries: only `krhom` by default, developer binaries/examples behind feature `dev`; delete unreferenced one-offs
5. tools: user tools at top, developer scripts in `tools/dev/`, archived-run scripts deleted
6. docs: `PLAN.md`, `docs/colored_design/` -> `docs/dev/`; README map and links
7. This report

## Stamp / key check (done first)

- No `std::any::type_name`, `module_path!`, `file!()` or `line!()` anywhere in src/, tests/, examples/.
  (`std::any::Any` appears only for downcasts in `cclose`/`xengine` and two tests; `stringify!` only
  in a diagnostic message of `boxclose` that names struct fields, not module paths.)
- Every stamp and cache-key component is a literal `const &str`: `sindec::SINDEC_VERSION` ("krhom-sindec v2"),
  `LADDER_STAMP`, `M5I_STAMP`; `tmpl::TMPL_ALG` ("hs1.lcx1"); `homprov::COST_ALG` ("hs-kl1") and
  `host_stamp()` (machine-derived, not path-derived); `leaves::STAGE0_STAMP`, `STAGE0_STAMP_PF`, `LEAVES_STAMP`,
  `C27_STAMP`, `LEAF_PAIRS_STAMP`; `sengine::SENGINE_STAMP`, `swire::SWIRE_STAMP`, `boxclose::BOX_CODE_STAMP`,
  `boxrun::BOXRUN_STAMP`. Layer-A keys and route stamps (`sindec` `stamps()`, `swire` engine stamp,
  `smapio`) are built only from these literals and data (λ, primes, hashes). `smapio` has no stamp of its own.
- `include_str!("../data/libraries.tsv")` is in `src/calc.rs`, which did not move. Runtime paths use
  `CARGO_MANIFEST_DIR` (unchanged).
- Conclusion: nothing depended on module paths; nothing needed freezing. `tests/compat_paths.rs` pins
  every stamp literal so a later edit cannot change one silently.
- The automatic `crate::m` -> `crate::group::m` rewrite was checked to touch no string literal.
- Confirmed in practice: krhom loads the published `l222_n4` / `l2222_n4` libraries (data ee548d2) and
  reproduces main's output byte for byte (validation below).

## Module layout

| group | modules |
|---|---|
| `algebra/` | field, poly, mpoly, laurent, pmat, sparse, dense, linsys, symfun |
| `combinat/` | hecke, schur, coset, dcoset, web, moy |
| `soergel/` | sbim, soergel, indec, local |
| `uncolored/` | engine, engine2, tower |
| `colored/` | xengine, cclose, sengine, swire, smap, sexpr, tmpl, boxclose, boxrun |
| `library/` | homprov, sindec, smapio, leaves, leafgen, rexplan |
| `input/` | pd, crossing |
| `reference/` | cube, cabling (= former `colored.rs`) |
| `legacy/` | algebra, hom, untrunc |
| top | lib.rs, calc.rs |

Deviations from the proposed grouping:
- `algebra.rs`, `hom.rs` and `untrunc.rs` went to a new `legacy/` group, not `soergel/`. Their headers
  mark them OBSOLETE or superseded, and only the deleted `explore*` binaries used them.
- `colored.rs` became `reference/cabling.rs`, because the name `colored` is now the group.
- Name clashes: `krhom::algebra`, `krhom::soergel` and `krhom::colored` are now group modules. Their old
  items still resolve, because each group glob re-exports them: `algebra` from `legacy::algebra::*`,
  `soergel` from `soergel::soergel::*`, `colored` from `reference::cabling::*`. Every other old top-level
  path is a `pub use` in lib.rs. Internal code, binaries, tests and examples use the new paths.

## Old path -> new path

The 225 files under `docs/colored_design/` moved to `docs/dev/colored_design/` with the same relative paths,
except the gl4 result files listed below.

| old | new |
|---|---|
| `PLAN.md` | `docs/dev/PLAN.md` |
| `docs/colored_design/c2/gl4_222.json` | `results/colored_gl4/gl4_222.json` |
| `docs/colored_design/c2/gl4_2222.json` | `results/colored_gl4/gl4_2222.json` |
| `docs/colored_design/c2/gl4_222_s0.json` | `results/colored_gl4/gl4_222_s0.json` |
| `docs/colored_design/c2/gl4_222_s1.json` | `results/colored_gl4/gl4_222_s1.json` |
| `docs/colored_design/c2/gl4_222_s2.json` | `results/colored_gl4/gl4_222_s2.json` |
| `docs/colored_design/c2/gl4_222_summary.md` | `results/colored_gl4/gl4_222_summary.md` |
| `src/algebra.rs` | `src/legacy/algebra.rs` |
| `src/boxclose.rs` | `src/colored/boxclose.rs` |
| `src/boxrun.rs` | `src/colored/boxrun.rs` |
| `src/cclose.rs` | `src/colored/cclose.rs` |
| `src/colored.rs` | `src/reference/cabling.rs` |
| `src/coset.rs` | `src/combinat/coset.rs` |
| `src/crossing.rs` | `src/input/crossing.rs` |
| `src/cube.rs` | `src/reference/cube.rs` |
| `src/dcoset.rs` | `src/combinat/dcoset.rs` |
| `src/dense.rs` | `src/algebra/dense.rs` |
| `src/engine.rs` | `src/uncolored/engine.rs` |
| `src/engine2.rs` | `src/uncolored/engine2.rs` |
| `src/field.rs` | `src/algebra/field.rs` |
| `src/hecke.rs` | `src/combinat/hecke.rs` |
| `src/hom.rs` | `src/legacy/hom.rs` |
| `src/homprov.rs` | `src/library/homprov.rs` |
| `src/indec.rs` | `src/soergel/indec.rs` |
| `src/laurent.rs` | `src/algebra/laurent.rs` |
| `src/leafgen.rs` | `src/library/leafgen.rs` |
| `src/leaves.rs` | `src/library/leaves.rs` |
| `src/linsys.rs` | `src/algebra/linsys.rs` |
| `src/local.rs` | `src/soergel/local.rs` |
| `src/moy.rs` | `src/combinat/moy.rs` |
| `src/mpoly.rs` | `src/algebra/mpoly.rs` |
| `src/pd.rs` | `src/input/pd.rs` |
| `src/pmat.rs` | `src/algebra/pmat.rs` |
| `src/poly.rs` | `src/algebra/poly.rs` |
| `src/rexplan.rs` | `src/library/rexplan.rs` |
| `src/sbim.rs` | `src/soergel/sbim.rs` |
| `src/schur.rs` | `src/combinat/schur.rs` |
| `src/sengine.rs` | `src/colored/sengine.rs` |
| `src/sexpr.rs` | `src/colored/sexpr.rs` |
| `src/sindec.rs` | `src/library/sindec.rs` |
| `src/smap.rs` | `src/colored/smap.rs` |
| `src/smapio.rs` | `src/library/smapio.rs` |
| `src/soergel.rs` | `src/soergel/soergel.rs` |
| `src/sparse.rs` | `src/algebra/sparse.rs` |
| `src/swire.rs` | `src/colored/swire.rs` |
| `src/symfun.rs` | `src/algebra/symfun.rs` |
| `src/tmpl.rs` | `src/colored/tmpl.rs` |
| `src/tower.rs` | `src/uncolored/tower.rs` |
| `src/untrunc.rs` | `src/legacy/untrunc.rs` |
| `src/web.rs` | `src/combinat/web.rs` |
| `src/xengine.rs` | `src/colored/xengine.rs` |
| `tools/c2/coverx.py` | `tools/dev/coverx.py` |
| `tools/c2/k4m_status.py` | `tools/dev/k4m_status.py` |
| `tools/c2/persist_c2.sh` | `tools/dev/persist_c2.sh` |
| `tools/c2/same.py` | `tools/dev/same.py` |
| `tools/c2/xdist.py` | `tools/dev/xdist.py` |
| `tools/c2/xfail.py` | `tools/dev/xfail.py` |
| `tools/c2/xscan_vs_wang.py` | `tools/dev/xscan_vs_wang.py` |
| `tools/libbuild_report.py` | `tools/dev/libbuild_report.py` |
| `tools/ntrunc_t1/Cargo.toml` | `tools/dev/ntrunc_t1/Cargo.toml` |
| `tools/ntrunc_t1/data/k2.txt` | `tools/dev/ntrunc_t1/data/k2.txt` |
| `tools/ntrunc_t1/data/knots3_all.txt` | `tools/dev/ntrunc_t1/data/knots3_all.txt` |
| `tools/ntrunc_t1/data/knots4_cloud.txt` | `tools/dev/ntrunc_t1/data/knots4_cloud.txt` |
| `tools/ntrunc_t1/data/t2n.txt` | `tools/dev/ntrunc_t1/data/t2n.txt` |
| `tools/ntrunc_t1/src/main.rs` | `tools/dev/ntrunc_t1/src/main.rs` |
| `tools/persist_layer_a.sh` | `tools/dev/persist_layer_a.sh` |
| `tools/trace_windows.py` | `tools/dev/trace_windows.py` |
| `runs_c2/c27/k4/compare.json` | `results/colored_gl4/gl4_2222_c27_compare.json` |
| `runs_c2/k3/out/4_1_f31.json` | `results/colored_gl4/gl4_222_4_1_f31.json` |
| `runs_c2/c27/k4/box.sh`, `finish.sh` | `tools/dev/box.sh`, `finish.sh` (rewritten: `CACHE`/`COSTS`/`OUTD`/`KRHOM_REPO` instead of hard-coded `/home/user/...` paths; `COSTS` defaults to the `costs.txt` that `fetch_library.sh` installs, identical to `costs_2222_c25.txt`) |

## Deleted files (all kept on `archive-pre-restructure`)

- `runs_m6a1c`, `runs_c2`, `runs_calc`, `runs_calc2`, `runs_calc3`, `runs_stretch`, `runs_m6a2a`, `runs_g2b`, `runs_random` (3162 files).
- Binaries that nothing referenced except old build logs: `src/bin/{explore,explore2..explore7,debug_dec,bench_hom,check_indec,check_leaves,check_local,check_soergel,colored_check}.rs`.
- Examples that nothing referenced: `examples/{m5i_force_search,m5i_prefix,m5i_srcjobs,m5i_words}.rs`.
- Scripts that only served archived runs: `tools/{k3_compare,stretch_compare,g2b_compare}.py`, `tools/c2/{k3_c2_compare,k4_markov,lib_summary,markov,project}.py`.

Kept behind `--features dev`, because tests, tools or docs use them: `sbox`, `xscan`, `libbuild`, `engine2`,
`pd2braid`, `twoprime`, `gp0`, `g1`, `g4_free`, `schur_scan`, `m4c`, `sprof`, `grwarm`, `b4wstamp`,
`colored_cable`, `cube`, `engine`, and the examples `e2ref`, `m5i_cover`, `idem_end0`, `idem_list`,
`ladder_search`, `m5i_check`, `probe_final_d2`.

Other changes:
- Tests that run developer binaries are gated: `tests/libbuild.rs` (whole file, `required-features`);
  4 schur_scan tests in `tests/schur.rs`; `xscan_json_provenance` in `tests/boxrun.rs`; the g4_free call
  inside one `tests/cclose.rs` test (cfg block).
- `randcheck.py` now builds `m5i_cover` with `--features dev`. Its default `--out` is `target/randcheck`
  instead of `runs_random`.
- `tools/ntrunc_t1` moved to `tools/dev/ntrunc_t1`; its path dependency is now `../../..`.

## Validation

Builds: `KRHOM_THREADS=4 cargo build --release` builds only `krhom` (checked with
`--message-format=json`). `cargo build --release --features dev` builds all 18 binaries. Both have no errors.
The warnings are the ones main already had.

Tests (`KRHOM_THREADS=4`, `--release`, default `#[ignore]` set):
- default features: 305 passed, 0 failed, 41 ignored. Main's suite (runs_calc3/full_suite.log) had 311.
  The difference: −8 tests that need dev binaries (libbuild 3, schur 4, boxrun 1), +2 `compat_paths`.
- `--features dev`: 313 passed, 0 failed, 41 ignored (= main's 311 + the 2 `compat_paths` tests).
- `tools/dev/ntrunc_t1`: `cargo check` passes against the new paths.

krhom checks, run with `KRHOM_THREADS=4` against the libraries installed by `tools/fetch_library.sh l222_n4 l2222_n4`
(data ee548d2). The main binary was built from origin/main 00f47df before any change.

| command | exit main / new | time main / new | stdout except `time` | stderr (timings masked) |
|---|---|---|---|---|
| `krhom 2,3 --braid aBaB`  | 0 / 0 | 0.05s / 0.04s | identical | identical |
| `krhom 2,3,4 --braid aaa --colors 2`  | 0 / 0 | 0.06s / 0.07s | identical | identical |
| `krhom 2,3,4 --braid aaa --colors 2 --reduced`  | 0 / 0 | 0.06s / 0.06s | identical | identical |
| `krhom 2,3,4 --braid "1 -2 1 -2" --colors 2 --reduced`  | 0 / 0 | 21.53s / 20.64s | identical | identical |
| `krhom 2,3,4 --braid "1 2 3" --colors 2 --no-simplify --reduced`  | 0 / 0 | 66.58s / 69.87s | identical | identical |

Diff of the stdout N blocks: **empty** for all five. The only stdout line that differs is the wall-clock `time` line.

### Outputs (restructure branch)

`krhom 2,3 --braid aBaB`

```
gl_N for N in [2, 3]  strands 3  crossings 4  components 1  colors [1]
braid: 1 -2 1 -2
closure: tower (Right end first)
coefficients: Q
time 0.04s
N = 2
  Poincare polynomial: t^-2q^5 + t^-1q^1 + t^0q^-1 + t^0q^1 + t^1q^-1 + t^2q^-5
  total rank: 6
  Euler characteristic: q^5 + q^-5
N = 3
  Poincare polynomial: t^-2q^6 + t^-2q^8 + t^-1q^0 + t^-1q^2 + t^0q^-2 + t^0q^0 + t^0q^2 + t^1q^-2 + t^1q^0 + t^2q^-8 + t^2q^-6
  total rank: 11
  Euler characteristic: q^8 + q^6 - 1 + q^-6 + q^-8
```

`krhom 2,3,4 --braid aaa --colors 2`

```
gl_N for N in [2, 3, 4]  strands 2  crossings 3  components 1  colors [2]
braid: 1 1 1
method: singular (exact 2-strand singular Soergel engine)
coefficients: Q
time 0.07s
N = 2
  Poincare polynomial: t^0q^0
  total rank: 1
  Euler characteristic: 1
N = 3
  Poincare polynomial: t^0q^-6 + t^0q^-4 + t^0q^-2 + t^2q^-8 + t^2q^-6 + t^3q^-14 + t^3q^-12
  total rank: 7
  Euler characteristic: q^-2 + q^-4 + 2q^-6 + q^-8 - q^-12 - q^-14
N = 4
  Poincare polynomial: t^0q^-12 + t^0q^-10 + 2t^0q^-8 + t^0q^-6 + t^0q^-4 + t^2q^-14 + 2t^2q^-12 + 2t^2q^-10 + t^2q^-8 + t^3q^-22 + 2t^3q^-20 + 2t^3q^-18 + t^3q^-16 + t^4q^-16 + t^5q^-24 + t^5q^-22 + t^6q^-30
  total rank: 22
  Euler characteristic: q^-4 + q^-6 + 3q^-8 + 3q^-10 + 3q^-12 + q^-14 - 2q^-18 - 2q^-20 - 2q^-22 - q^-24 + q^-30
```

`krhom 2,3,4 --braid aaa --colors 2 --reduced`

```
gl_N for N in [2, 3, 4]  strands 2  crossings 3  components 1  colors [2]
braid: 1 1 1
method: singular (exact 2-strand singular Soergel engine), reduced (strand at bottom position 0 marked)
coefficients: Q
time 0.06s
N = 2
  Poincare polynomial: t^0q^0
  total rank: 1
  Euler characteristic: 1
N = 3
  Poincare polynomial: t^0q^-4 + t^2q^-8 + t^3q^-12
  total rank: 3
  Euler characteristic: q^-4 + q^-8 - q^-12
N = 4
  Poincare polynomial: t^0q^-8 + t^2q^-14 + t^2q^-12 + t^3q^-18 + t^3q^-16 + t^4q^-20 + t^5q^-24 + t^5q^-22 + t^6q^-26
  total rank: 9
  Euler characteristic: q^-8 + q^-12 + q^-14 - q^-16 - q^-18 + q^-20 - q^-22 - q^-24 + q^-26
```

`krhom 2,3,4 --braid "1 -2 1 -2" --colors 2 --reduced`

```
gl_N for N in [2, 3, 4]  strands 3  crossings 4  components 1  colors [2]
braid: 1 -2 1 -2
method: singular (localized box engine; library l222_n4 @ data ee548d2, labels (2,2,2), prune 4, box N 4, seed 201), reduced (strand at bottom position 0 marked)
coefficients: Q
time 20.64s
N = 2
  Poincare polynomial: t^0q^0
  total rank: 1
  Euler characteristic: 1
N = 3
  Poincare polynomial: t^-2q^6 + t^-1q^2 + t^0q^0 + t^1q^-2 + t^2q^-6
  total rank: 5
  Euler characteristic: q^6 - q^2 + 1 - q^-2 + q^-6
N = 4
  Poincare polynomial: t^-4q^12 + t^-3q^8 + t^-3q^10 + 2t^-2q^6 + t^-2q^8 + 2t^-1q^2 + 2t^-1q^4 + t^0q^-2 + 3t^0q^0 + t^0q^2 + 2t^1q^-4 + 2t^1q^-2 + t^2q^-8 + 2t^2q^-6 + t^3q^-10 + t^3q^-8 + t^4q^-12
  total rank: 25
  Euler characteristic: q^12 - q^10 + 2q^6 - 2q^4 - q^2 + 3 - q^-2 - 2q^-4 + 2q^-6 - q^-10 + q^-12
```

`krhom 2,3,4 --braid "1 2 3" --colors 2 --no-simplify --reduced`

```
gl_N for N in [2, 3, 4]  strands 4  crossings 3  components 1  colors [2]
braid: 1 2 3
method: singular (localized box engine; library l2222_n4 @ data ee548d2, labels (2,2,2,2), prune 4, box N 4, seed 201), reduced (strand at bottom position 0 marked)
coefficients: Q
time 69.87s
N = 2
  Poincare polynomial: t^0q^0
  total rank: 1
  Euler characteristic: 1
N = 3
  Poincare polynomial: t^0q^0
  total rank: 1
  Euler characteristic: 1
N = 4
  Poincare polynomial: t^0q^0
  total rank: 1
  Euler characteristic: 1
```

## Not done / notes

- `docs/dev/colored_design/**` keep their historical `runs_*` and flat `src/<module>.rs` paths. `docs/dev/README.md`
  explains that they refer to `archive-pre-restructure` and to this table. Paths to moved results and
  tools were updated.
- The ∧³ σ^9 data stay in `docs/dev/colored_design/gl6/m2s/`, because only `uncolored/` and `colored_gl4/`
  were requested under `results/`.
- `tools/dev/box.sh` / `finish.sh` pass `bash -n` but were not run, because a run takes hours.

## New tree (top two levels)

```
.gitignore
CITATION.cff
Cargo.lock
Cargo.toml
LICENSE-APACHE
LICENSE-MIT
README.md
build.rs
data
  README.md
  knots_upto12.jsonl
  libraries.tsv
  libstatus
  links_upto11.jsonl
  tools
docs
  CALCULATOR.md
  dev
examples
  e2ref.rs
  idem_end0.rs
  idem_list.rs
  ladder_search.rs
  m5i_check.rs
  m5i_cover.rs
  probe_final_d2.rs
results
  README.md
  colored_gl4
  uncolored
src
  algebra
  bin
  calc.rs
  colored
  combinat
  input
  legacy
  lib.rs
  library
  reference
  soergel
  uncolored
tests
  acceptance.rs
  b4w.rs
  boxclose.rs
  boxclose_mem.rs
  boxrun.rs
  c27_joint.rs
  calc.rs
  cclose.rs
  common
  compat_paths.rs
  data
  dcoset.rs
  hom_lazy.rs
  homprov.rs
  knotcov.rs
  leafgen.rs
  leaves.rs
  libbuild.rs
  m4a.rs
  m5i.rs
  ntrunc.rs
  pd.rs
  rexplan.rs
  schur.rs
  sengine.rs
  sindec.rs
  swire.rs
  thmassoc.rs
  xengine.rs
  xstream.rs
tools
  batch_check.py
  dev
  fetch_library.sh
  randcheck.py
  wang.py
```
