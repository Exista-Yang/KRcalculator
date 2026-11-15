# krhom calculator

`krhom` computes colored Khovanov–Rozansky gl_N homology of a link closed from a braid, over Q,
in one command. You give the N values, the braid word, and the colors. The program chooses the
engine, runs it, and prints the Poincaré polynomial for every N.

```
cargo build --release
tools/fetch_library.sh l222_n4 l2222_n4          # only for colored braids on 3-4 strands (once)
target/release/krhom 2,3,4 --braid aBaB --colors 2
```

Contents:
1. [Usage](#1-usage)
2. [What is supported](#2-support-matrix)
3. [Installing the libraries](#3-installing-the-libraries)
4. [What runs, and what "over Q" means](#4-what-runs)
5. [Exit codes and messages](#5-exit-codes-and-messages)
6. [Cost guide](#6-cost-guide)
7. [The random-example checker](#7-the-random-example-checker)
8. [Advanced and developer options; sbox](#8-advanced-and-developer-options)

## 1. Usage

```
krhom N[,N...] (--braid WORD [--strands n] | --pd PD | --file BRAIDS | --knot NAME | --link NAME)
               [--colors c | c1,c2,...] [--reduced [--mark K]]
               [--table] [--latex OUT] [--estimate]
               [--json OUT] [--csv OUT] [--no-simplify] [--check-library] [--lib-dir DIR]
               [--seed S] [--seed-retries K] [--mem-mb M] [--mem-cap-mb M]
               [--colored-method cable|cube|singular] [--cube] [--closure tower|exact] [--end right|left|alt]
```

| Argument | Meaning |
|---|---|
| `N[,N...]` | The values of N (the first argument), e.g. `2,3,4`. |
| `--braid WORD` | Either signed integers, `"1 -2 1 -2"` (k = σ_k, −k = σ_k⁻¹; commas also work), or KnotJob letters, `aBaB` (a..i = σ1..σ9, A..I their inverses; letters only, no spaces). A word that mixes the two forms is refused. |
| `--strands n` | Optional. Default: max \|generator\| + 1 (the KnotJob convention). If given, it must be at least that; extra strands are unknotted, split components. |
| `--pd PD` | A PD code instead of a braid, e.g. `"[[1,4,2,5],[3,6,4,1],[5,2,6,3]]"` (converted with Vogel's algorithm and simplified). |
| `--colors` | A single color applies to every component (`--colors 2`); a list gives one color per link component, numbered by the smallest bottom strand position. Default: all 1 (uncolored). |
| `--reduced` | Reduced homology: the strand at bottom position 0 is marked; the reduced unknot is `t^0q^0`. |
| `--mark K` | With `--reduced`: mark link component K instead (0-based, krhom's numbering: components are numbered by their smallest bottom strand position, the order of `--colors`). krhom conjugates the braid so that component K passes through bottom position 0 (a cyclic rotation of the word when one does it, else (σ1⋯σp) β (σp⁻¹⋯σ1⁻¹), freely reduced); a conjugate braid has the same closure and the engines are unchanged. The line `conjugated braid: ...` shows the word used, and the method line says `reduced (component K marked ...)`. For a knot only `--mark 0` is accepted; `--mark` without `--reduced` is refused. |
| `--knot NAME` | The braid of a prime knot from `data/knots_upto12.jsonl` (`3_1` … `12n_888`; the KnotInfo braid, so KnotInfo's chirality). The header starts with `name: NAME` and `source: knot NAME from data/knots_upto12.jsonl (KnotInfo braid; KnotInfo chirality)`. `--strands` must be absent or equal the table's. |
| `--link NAME` | The braid of an oriented prime link from `data/links_upto11.jsonl`, named with LinkInfo's orientation suffix, e.g. `L2a1{1}` (σ1²) or `L2a1{0}` (σ1⁻²); an unoriented name (`L2a1`) is refused with the list of its orientations. An unknown name (knot or link) exits 1 with up to 5 close matches. The tables are read from the source tree's `data/` (override: `$KRHOM_DATA_DIR`). |
| `--file BRAIDS` | Batch: one braid per line, optional `NAME;` prefix, blank lines and `#` comments ignored. N, colors, `--strands` and `--reduced` apply to every line. One output block per line, headed `name: NAME` (printed first, also for a line that fails to parse). An error on one line is reported (stdout: `error: NAME: exit C (kind)`; stderr: `krhom: NAME: message`) and the batch continues; the exit code is that of the first failing line. |
| `--json OUT` | Machine-readable results (an array with `--file`): input (`kind`, the input string as `text`, N, mode, the parsed `braid` array, ...), route, library, per-prime provenance (seeds, timings, χ and scaling checks, layer-A loads, VmHWM), results with terms `[t,q,k]`, and the error if any. With `--file` the array is rewritten after every finished line (atomically: temporary file + rename), so an abort (out of memory under `ulimit -v`, SIGKILL, a timeout) keeps the finished records. |
| `--csv OUT` | Rows `name,N,reduced,t,q,rank` (a zero polynomial is one row with empty t, q and rank 0). The name is quoted per RFC 4180 when it contains a comma, a quote or a line break (`"L6a5{0,1}"`). The header is written at start and each line's rows are appended and flushed when it finishes. |
| `--table` | After each N block, the bigraded table: rows t ascending, columns q ascending (step 2 when all q have one parity), entries = dimensions, blank for 0 (KnotJob-style). Also in `--file` mode. |
| `--latex OUT` | A LaTeX file with one `tabular` per N (rows t, columns q, as `--table`) in a `table` float whose caption names the knot/braid, the colors, N, reduced/unreduced (and the marked component) and the coefficients (`Q`). With `--file`: one `\section*` per input (a failed line gets a "Not computed: exit C" note). The file is a complete document, `\documentclass{article}` + `\usepackage{graphicx}` (wide tables are scaled to the line width), that compiles standalone with `pdflatex OUT`; the part between `% ---- krhom tables begin` and `% ---- krhom tables end` can be pasted into a document that loads graphicx. Rewritten after every finished input. |
| `--estimate` | Do not compute; estimate. On the 3-4 strand box route krhom loads the library, runs the record scan and times the determinant pass and replay of ONE box point per prime (the single-point path of `sbox point`), then prints `estimate:` lines: the number of box points (216 reduced / 1296 unreduced for a 4-strand word), and an extrapolated wall time and peak memory for the full run at the current thread count (`KRHOM_THREADS`). It is an estimate: the closure phase is not included (minutes on 4 strands), points differ in cost, and threads contend for memory bandwidth. On other routes it prints `estimate: fast route (ROUTE), no estimate needed` (the full cube: no estimate) and computes nothing. Exit 0. Not with `--check-library`. |

Examples:

```
krhom 2,3,4 --braid aBaB --colors 2                     # 4_1, color 2, 3 strands (library l222_n4)
krhom 2,3,4 --braid aBaB --colors 2 --reduced
krhom 2,3,4 --braid "1 1 2 -1 -3 2 -3" --colors 2 --check-library   # 6_1 on 4 strands: is it covered?
krhom 2,3,4,5 --braid aaa --colors 3                    # trefoil, color 3 (2 strands, exact engine)
krhom 2,3 --braid "1 1" --colors 2,1                    # Hopf link, colors (2,1)
krhom 2,3,4 --braid aBaB                                # uncolored 4_1 (Soergel engine)
krhom 3 --pd "[[1,4,2,5],[3,6,4,1],[5,2,6,3]]" --reduced
krhom 2,3 --file knots.txt --colors 2 --csv out.csv --json out.json
krhom 2,3,4 --knot 4_1 --colors 2 --reduced               # same as --braid "1 -2 1 -2"
krhom 2,3,4 --knot 6_1 --colors 2 --reduced --estimate    # time one box point, extrapolate
krhom 2,3,4 --braid "1 1" --colors 2,1 --reduced --mark 1 # Hopf link, the color-1 component marked
krhom 2,3 --link "L6a4{0,0}" --reduced --mark 2           # Borromean rings, component 2 marked
krhom 2,3 --knot 3_1 --colors 2 --latex trefoil.tex       # pdflatex trefoil.tex
```

Output (stdout):

```
gl_N for N in [2, 3, 4]  strands 3  crossings 4  components 1  colors [2]
braid: 1 -2 1 -2
method: singular (localized box engine; library l222_n4 @ data ee548d2, labels (2,2,2), prune 4, box N 4, seed 201)
coefficients: Q
time 81.84s
N = 2
  Poincare polynomial: t^0q^0
  total rank: 1
  Euler characteristic: 1
N = 3
  Poincare polynomial: t^-2q^6 + t^-2q^8 + t^-1q^0 + t^-1q^2 + t^0q^-2 + t^0q^0 + t^0q^2 + t^1q^-2 + t^1q^0 + t^2q^-8 + t^2q^-6
  total rank: 11
  Euler characteristic: q^8 + q^6 - 1 + q^-6 + q^-8
N = 4
  ...
```

Terms are sorted by t, then q; `2t^-2q^8` means rank 2 in homological degree −2 and q-degree 8.
The convention is that of the repository (ours(t,q) = Wang(1/t,1/q), see tools/wang.py): positive
knots sit in negative q-degrees. krhom's q-grading is the inverse of Bar-Natan / KnotJob / KnotInfo
(q → q⁻¹, t unchanged): the right-handed trefoil `aaa` at N = 2 is `t^0q^-3 + t^0q^-1 + t^2q^-5 +
t^3q^-9` here and `q + q^3 + t^2q^5 + t^3q^9` there. On the
2-strand and Soergel routes the method line reads `method: singular (exact 2-strand singular
Soergel engine)` or `closure: tower (Right end first)`. If simplification removed strands, a line
`simplified braid: WORD  (strands k)` follows the braid. Progress and the per-prime certification
lines go to stderr.

## 2. Support matrix

| Input | Route (automatic) | N | Strands |
|---|---|---|---|
| all colors 1 | localized Soergel engine (engine2), two primes | 1 ≤ N ≤ 32 (tower closure; N > 32 is refused with exit 1, `--closure exact` has no such cap) | ≤ 10 (tower closure; more strands are refused with exit 1 before any engine work); about 5 in practice |
| colored, ≤ 2 strands after simplification | exact singular Soergel engine (xengine), two primes | a(N − a) ≤ 255 for every color a (N ≤ 256 for color 1, N ≤ 129 for color 2, N ≤ 88 for color 3); dim A_λ = ∏ C(N, a_i) over the unmarked strands ≤ 5·10⁷ (the exact oracle's limit; e.g. color 3 needs N ≤ 35 and color 2 N ≤ 119 unreduced; refused with exit 1 before any engine work), and the closure's memory estimate under `--mem-cap-mb` (default 768 MB; e.g. color 3 at N = 30 needs about 1.5 TB and is refused with exit 5) | 1–2, colors 1..3 |
| colored, uniform color 2, 3 strands | localized box engine + library `l222_n4` (labels 2,2,2), two primes | N ≤ 4 | 3 |
| colored, uniform color 2, 4 strands | localized box engine + library `l2222_n4` (labels 2,2,2,2), two primes | N ≤ 4 | 4 |
| anything else colored | refused (exit 1) with this matrix; `--colored-method cube` is exact over Q but exponential (tiny inputs only) | | |

- N below a color: the homology is 0 (Λ^a C^N = 0 for a > N). krhom prints the zero polynomial
  for those N on every route, without running an engine, and says so on stderr. When every
  requested N is below the largest color, the answer is printed before any route refusal (`method:
  none (every N is below the color ...)`), e.g. `krhom 1,2 --braid aBaB --colors 3` prints 0 and
  exits 0 although no (3,3,3) library exists.
- Colored braids are simplified first (cyclic reduction, Markov destabilization, a deterministic
  random search). The simplified braid is used only if it has fewer strands, and for `--reduced`
  on a link only if the marked component is unchanged. `--no-simplify` keeps the given word.
- Split braids on 3–4 strands (some σ_i never occurs) are refused: compute the pieces separately.
- Mixed colors on 3–4 strands, colors other than 2, (3,3,3), 5 or more strands, and N > 4 on 3–4
  strands are refused: no certified library exists for them.
- 4-strand coverage: the (2,2,2,2) library certifies 1236 of the 1254 jobs of its plan for both
  primes. The 18 class-A 4-strand knots of the C2-7 cover list (6_1, 7_2, 7_4, 7_6, 7_7, 8_4, 8_6,
  8_8, 8_11, 8_13, 8_14, 9_10, 9_11, 9_13, 9_17, 9_18, 9_20, 9_23 in their listed braids) miss no
  job; 8_15, 9_22 and 11n_42 miss 11 jobs and the Conway knot 16. A word that needs a missing job
  is refused with exit 4 and the list of jobs. `--check-library` answers the coverage question for
  a word without running the expensive box scan.

## 3. Installing the libraries

Colored braids on 3 and 4 strands need prebuilt, certified layer-A libraries from the data
commit `ee548d2` (branch `data-gl6`). Install them once:

```
tools/fetch_library.sh l222_n4 l2222_n4      # both primes each; about 2 min
tools/fetch_library.sh --verify-only l222_n4 l2222_n4
```

| Library | Labels | N | Download per prime | Installed per prime | Layers (restore order) |
|---|---|---|---|---|---|
| `l222_n4` | (2,2,2) | ≤ 4 | 3.5 MB | 21 MB (233 layer-A, 471 template files) | c2/l222_n4 |
| `l2222_n4` | (2,2,2,2) | ≤ 4 | 16.2 MB (F31), 18.2 MB (F32b) | 68–75 MB (2873 layer-A, about 3040 template files) | c2/l2222_n4, c2/l2222_n4_c25, c2/l2222_n4_c24, c2/l2222_n4_c27 |

- Root: `--lib-dir DIR` > `$KRHOM_LIB_DIR` > `${XDG_CACHE_HOME:-$HOME/.cache}/krhom/lib`. The tool
  takes `--root DIR` (same default).
- Layout: `<root>/<id>/<f31|f32b>/{LIBRARY.txt, costs.txt, status.txt, MANIFESTS/, layer_a/, tmpl/}`.
- Transport: the tool fetches the data commit by SHA into the current clone (`git fetch origin
  <sha>`; the work tree is untouched) or, with `--repo URL`, into a temporary bare repository. It
  falls back to fetching the `data-gl6` branch if the server refuses SHA fetches.
- Every install step is checked explicitly (`|| die`; errexit is inactive inside the install
  function), an empty or missing staging path is refused, and the old copy is kept if the new one
  cannot be moved into place.
- Checks: every chunk's sha256 must equal both the registry (`data/libraries.tsv`) and the
  MANIFEST of the data commit; the tar listing must contain only regular files `layer_a/a-<16hex>.txt`
  and `tmpl/tmpl-<16hex>.txt` (no links, no `..`), with the MANIFEST's counts.
- `LIBRARY.txt` records the registry fields, the sources of `costs.txt` (`costs_source:` the
  path in the data commit) and `status.txt` (`status_source:` the checkout path and its sha256),
  and a per-file list `LAYER SHA256 PATH`. "already installed" and `--verify-only` re-check both
  sources against the registry and the current `data/libstatus` file, so a corrected status file
  or costs path triggers a reinstall (installs made before CALC-2 lack these lines and are
  reinstalled once).
  `--verify-only` checks every listed file's sha256; extra files (templates a later run may write)
  are tolerated. A second install run prints `already installed`; `--force` reinstalls.
- krhom itself checks on every run that LIBRARY.txt matches the registry, that every listed file
  exists, and that every layer-A file of the stamped layers carries this binary's route stamp
  (exit 4 otherwise, with the command to run).
- Integrity is not soundness: every layer-A map loaded from the library is re-certified on load
  (ends, degrees, scaling, Π∘Ι, forms, closure certificates), so a bad map is rejected, not used.

## 4. What runs

Every result is over Q. Precisely: **certified over F_p for two large primes (p = 2^31 − 1 and
p = 4294967291); equal over Q unless both are torsion primes.** The full cube
(`--colored-method cube`) is exact over Q. Small characteristic is not supported. The per-prime
details go to stderr and `--json`; stdout says `coefficients: Q`.

- Uncolored (Soergel engine2): F_2147483647 with seed 7 and F_4294967291 with seed 11 (as
  `src/bin/twoprime.rs`); the results must be equal, otherwise exit 3.
- Colored, ≤ 2 strands (exact singular engine): both primes, compared, exit 3 on disagreement.
- Colored, 3–4 strands (box route), for each prime:
  - the library is loaded from layer A in replay-only mode (the provider refuses every call;
    loaded maps are re-certified);
  - `swire::box_word` with the production options: the generator-determinant pass at every box
    point (a failure redraws the nodes with seed+1, at most `--seed-retries` = 2 times), the
    record scan and frozen-structure replay at every box point, Freivalds D² (k = 3), χ =
    `decat_chi` for every closure (fatal), and the D14 scaling check (fatal);
  - box N = 4, Foam_N pruning at 4, node seed 201 for both primes (`--seed`);
  - post-checks: one stage (the whole word), the scaling report present, χ checked for every
    closure, every requested (N, mode) present.

  Then the two primes must agree exactly (exit 3 otherwise).

## 5. Exit codes and messages

| Code | Meaning | Example message (stderr, prefix `krhom: `) |
|---|---|---|
| 0 | success | |
| 1 | usage or unsupported input | `N = 130 with color 2: the exact 2-strand engine supports a(N - a) <= 255 ...`; `N = 36 with labels [3, 3]: the exact 2-strand engine's algebra A_lambda has dimension 50979600 ..., above its limit 5e7 ...`; `N = 33: the tower closure supports N <= 32 ...; use --closure exact for larger N`; `11 strands: the tower closure of the Soergel engine supports at most 10 strands ...`; `--knot 4_2: unknown knot name (data/knots_upto12.jsonl); close matches: 4_1, 5_2, 6_2, 7_2, 8_2`; `--mark 1: a knot has one component; only --mark 0 is accepted`; `labels [2, 2, 1] on 3 strands are not uniform; supported on 3-4 strands: ...`; `no library for labels [3, 3, 3] on 3 strands ...`; `N = 5 is above the library's N <= 4 ...`; `split link (sigma_2 never occurs): compute the pieces separately`; `3 colors given, but the link has 1 components ...`; `braid "a1" mixes letters and integers ...`; `uncolored input: drop --colored-method singular ...`; `validation-only environment switches change the computation on the box route: unset KRHOM_NO_M5I` |
| 2 | engine or certification failure | `localized engine failed on F_2147483647: <error>. This is a certification failure; please report it with --json output.` |
| 3 | the two primes disagree; no result printed | `the two primes disagree (torsion or a bug); no result printed: N = 3: ...` |
| 4 | library missing, wrong version, stamp mismatch, damaged, or not covering the word | `library l222_n4 (F_2147483647) is not installed under <root>: run tools/fetch_library.sh l222_n4 (about 3.5 MB per prime)`; `damaged install: 1 of the 5921 files listed in .../LIBRARY.txt are missing (...): run tools/fetch_library.sh --force l2222_n4`; `library l2222_n4 does not cover this word: 11 jobs not certified for both primes: pair 926 (p=..., rung=...), ...` |
| 5 | resource refusal | `per-point memory ... above the budget ...` (raise `--mem-mb`); `... above the memory cap of ... MB ...; with krhom use --mem-cap-mb M to raise it` (on every route, including the exact 2-strand engine) |

When a library job is missing, krhom diagnoses which jobs the word needs from the shipped status
files (`data/libstatus/`; this costs about 25 s for the (2,2,2,2) plan, only on failure). If the
files are all present but a map was rejected on load, the message says so and points to the
`[sindec] layer A ... (rejected)` lines; if a template was missing in replay-only mode, krhom
retries once with template recomputation enabled (deterministic, does not affect soundness; the
JSON records `"template_recompute": true`).

## 6. Cost guide

Measured on the 4-core cloud machine (`KRHOM_THREADS=4`, both primes, this branch):

| Input | Time (F31 + F32b) | Peak memory |
|---|---|---|
| 2-strand colored braids (exact engine) | well under 1 s | small |
| `aBaB --colors 2` (4_1, 3 strands), unreduced | 44 s + 38 s | 166 MB |
| same, `--reduced` | 17 s + 16 s | 160 MB |
| `"1 1 1 2" --no-simplify` (3 strands) | 0.9 s + 0.9 s | 16 MB |
| `"1 2 3" --no-simplify` (4 strands), unreduced | 83 s + 126 s | 119 MB |
| `"1 2 3" --no-simplify --reduced` | 46 s + 56 s | 116 MB |
| `"1 1 1 2 3" --no-simplify` (4 strands), unreduced / reduced | 96 s + 91 s / 41 s + 55 s | 128 / 121 MB |
| 6_1 `"1 1 2 -1 -3 2 -3" --check-library` (library build only) | 222 s + 217 s | 1046 MB |
| 8_15 `--check-library` (not covered: exit 4 with the 11 missing jobs) | 101 s in total | |

- Loading the 4-strand library takes 20–34 s per prime for short words and about 220 s for 6_1-sized words (1119 layer-A loads); 3 strands 0.1–6 s.
- From earlier production runs (2 threads, per prime, both modes): 8-letter 3-strand knots
  100–1151 s, 10-letter 247–2564 s. F32b is about 1× F31 at 3 strands and 1–3× at 4 strands.
- 4-strand class-A knots, measured reduced runs (C2-7, docs/dev/colored_design/c2/C2_7_report.md §6;
  results/colored_gl4/gl4_2222_c27_compare.json, per-shard logs in runs_c2/c27/k4/out/ on branch archive-pre-restructure): `sbox box` with 4 shards of
  1 thread over the 216 reduced box points (54 per shard) + `sbox finish` (4 threads), per prime:

  | knot | prime | s per point | shard wall (s) | finish (s) | wall | shard core-hours | VmHWM per shard |
  |---|---|---|---|---|---|---|---|
  | 6_1 | F31 | 26.3–26.5 | 1983–1992 | 704 | 0.75 h | 2.21 | 2681–2682 MB |
  | 6_1 | F32b | 26.8–27.1 | 1928–1946 | 758 | 0.75 h | 2.15 | 2758–2768 MB |
  | 7_6 | F31 | 36.5–36.8 | 2639–2661 | 847 | 0.97 h | 2.95 | 2811 MB |
  | 7_6 | F32b | 39.5–39.9 | 2821–2832 | 822 | 1.01 h | 3.14 | 2898 MB |
  | 7_7 | F31 | 35.9–36.3 | 2597–2626 | 706 | 0.93 h | 2.90 | 2764 MB |
  | 7_7 | F32b | 25.7–25.9 | 1825–1840 | 672 | 0.70 h | 2.04 | 2850–2854 MB |

  (wall = slowest shard + finish; shard core-hours = sum of the 4 single-thread shards, without the
  4-thread finish; library build 168–196 s per process in finish.) All PASS: reduced N = 4 ranks
  81, 361, 441, F31 = F32b, χ and scaling checks. Unreduced runs were not done: they need all 1 296
  box points instead of 216, i.e. at least 1 296 × 26.3 s ≈ 9.5 core-hours of scanning per prime
  for 6_1 (a lower bound: determinant passes and finish come on top; the C2-7 report estimates about
  4 h wall per (knot, prime) with 4 shards). Use `--reduced` for speed. Runs are not resumable yet
  (the production `sbox box/finish` driver is).
- Memory: `--mem-mb` is the per-point budget (default 3200 at 3 strands, 3000 at 4 strands); the
  closure estimate cap `KRHOM_MEM_CAP_MB` is set from `--mem-cap-mb` before any engine work, on
  every route. Without the flag: the environment variable if set; else on the box route 4500 ×
  max(1, threads/2) at 3 strands and 200000 at 4 strands (where the estimate is 30–130× too
  high), and on the exact 2-strand route the library default of 768 MB. The real guard is the shell: run long jobs under
  `ulimit -v 14000000` (or less).
- `KRHOM_THREADS` sets the thread count (default: all cores).

## 7. The random-example checker

`tools/randcheck.py` draws random colored braids (color 2) and checks fixed relations; it never
sweeps knot tables.

```
tools/randcheck.py --seed 20261115 --quick        # T2 4, T3 6 examples (about 15-25 min)
tools/randcheck.py --seed 20261115                # 24 examples: T2 8, T3 12, T4 4 (about 1-1.5 h)
tools/randcheck.py --seed 20261115 --full         # 10 / 20 / 10
tools/randcheck.py --seed 20261115 --quick --inject-fault --only T3-01   # must FAIL (exit 1)
```

Options: `--tiers 2,3,4`, `--only ID`, `--krhom PATH`, `--jobs 2` (each krhom then gets
`KRHOM_THREADS=2` and an 8 GB address-space limit), `--out DIR` (default `target/randcheck`), `--no-anchor`.

Self-tests first (exit 2 if one fails): the 2-strand exact path equals `tools/wang.py` for the
(2,2) trefoil at N = 2,3,4 and the Hopf closed form of O11 at N = 2..5; the libraries of the requested tiers verify; the 4_1 anchor (both
modes, against the known strings; skipped with `--only` or `--no-anchor`).

Tiers: T2 (2 strands, exact path, N = 2..6), T3 (3 strands; Mode A: a 2-strand core stabilized by
σ2^{±1} plus 1–2 isotopy moves, length ≤ 7, two of them unknots; Mode B: a free word of length
4–6 with one partner: an isotopy move, the mirror or the reverse, partner length ≤ 6), T4 (4
strands; a 2-strand knot core stabilized twice, optionally rotated, length ≤ 5; only words the
library covers, filtered with `examples/m5i_cover`; reduced only, plus unreduced for length ≤ 3).
Box-route calls pass `--no-simplify` and use different seeds for the two primes (201 and
202 + subseed mod 97).

Relations (C = certifying: a mismatch is FAIL; B = blocking: a mismatch is FAIL, a pass is not
counted as confirmation; E = expected, not proven: a mismatch is WARN):

| # | Check | Label |
|---|---|---|
| O1 | two primes agree (krhom exit 3 = FAIL) | C |
| O2 | seed independence (different seeds per prime; on exit 3 rerun with equal seeds) | C |
| O3 | χ = decat_chi recorded for every closure of every prime | C |
| O4 | isotopy: rotation, conjugation, braid relation, cancelling pair | C |
| O5 | Markov: stabilized word = core value (exact 2-strand path) | C |
| O6 | mirror: (t,q,k) → (−t,−q,k) | C |
| O7 | reversal | E |
| O8 | colored-2 N'=3 of w = uncolored N=3 of rev(w), same mode | B |
| O9 | N'=2 is exactly t^0q^0 | C |
| O10b | [N 2]_q·P_red − P_unred = (1+t)·Q with Q ≥ 0 (knots) | E |
| O11 | 2-strand closed forms: Wang T(2,2m+1) via tools/wang.py (B for m = 4 and for m = 3 at N ≥ 5), unknot, unlink, Hopf \|n\| = 2 unreduced (Wang Thm 3.5 = tests/common `hopf_formula`, mirrored for n = −2) | C |
| O12 | unknot words: [N 2]_q unreduced, t^0q^0 reduced | C |
| O13a | `--box-live` (in-memory unpruned library, 2 strands) = exact path, N ≤ 4 | C |

Report: `<out>/<seed>/examples.jsonl` (one JSON line per example: id, subseed, tier, word,
components, mode, core, partner, runs with arguments, exit codes, times and results, checks with
label and detail, status), `summary.md` (one row per example, details and repro commands for every
non-PASS example), and the final line
`seed=... n=... PASS a FAIL b WARN c SKIP d TIMEOUT e ERROR f; checks: primes x, seeds y, ...`.
Exit 0 if there is no FAIL or ERROR, 1 otherwise, 2 on a harness failure: a self-test, the
convention call, `tools/wang.py`, the build or run of `examples/m5i_cover`, or any other exception
of the harness itself. krhom exit 4 (library missing, damaged, stamp/registry mismatch, map
rejected on load) is FAIL: T3 runs only with the complete `l222_n4` and T4 only on words the
m5i_cover filter found covered. TIMEOUT is never a PASS. Per-call timeouts: 120 s for the 2-strand
exact-path calls (T2 base, mirror and rotation; the T3/T4 cores; the convention self-test), 300 s
for the T2 `--box-live` call, 600 s for the uncolored O8 calls, 2400 s for the 3-4 strand box
calls (and the 4_1 anchor).

## 8. Advanced and developer options

Advanced (documented, supported):
- `--colored-method cable|cube|singular` overrides the automatic route. `cube` (also `--cube`) is
  the full cube of colored Rickard complexes, exact over Q and exponential. `cable` computes a
  colored link from a cable through the uncolored Soergel engine over both primes (sizes beyond
  about 4 cable strands used 0.4–6 GB in tests). `singular` is the default for colored input.
- `--closure tower|exact`, `--end right|left|alt`: the Soergel engine's closure.
- `--no-simplify`, `--check-library`, `--lib-dir`, `--seed`, `--seed-retries`, `--mem-mb`,
  `--mem-cap-mb`: see above.
- Environment: `KRHOM_THREADS`, `KRHOM_LIB_DIR`, `KRHOM_MEM_CAP_MB`. krhom sets
  `KRHOM_COSTS_ANY_HOST=1` on the box route (the shipped cost tables name another host). The
  validation-only switches `KRHOM_FORCE_LADDER`, `KRHOM_NO_M5I`, `KRHOM_DCID_REF`, `KRHOM_M5I_*`,
  `KRHOM_NO_REROUTE_FALLBACK` and `KRHOM_LEAVES_C27` change the computation and are refused on the
  box route; all `KRHOM_*` variables are recorded in `--json`.

Developer options (not in the usage text; for tests and the checker):
- `--seed-f32b S`: a different node seed for the second prime (oracle O2).
- `--single-prime`: run only F_2147483647; stdout then says
  `coefficients: F_2147483647 (developer option --single-prime: not certified over Q)`.
- `--field ...`: accepted for old scripts and ignored (a note on stderr).
- `--box-live`: the box route on 2 strands with the live provider and an in-memory unpruned
  library (cross-engine test O13a).
- `--name NAME`: a name for a single `--braid`/`--pd` input (as in `--file`).

`sbox` (src/bin/sbox.rs) is the production tool behind the box route: checkpointed, resumable,
sharded runs (`sbox box --dir D --shard k/K`, then `sbox finish`), any stage list, both modes in
one grid, and the library build options (`--provider default|m4c|m4c-lf`, `--replay`,
`--stamp-route`, `--cache`, `--costs`). Use it for multi-hour 4-strand words; krhom runs a word
in one process.
`sbox` and the other developer binaries are built only with the `dev` feature
(`cargo build --release --features dev`); `tools/dev/box.sh` and `tools/dev/finish.sh` drive a
sharded 4-strand run.
