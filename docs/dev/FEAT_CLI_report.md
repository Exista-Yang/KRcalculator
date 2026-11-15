# F-CLI report: limit refusals, --mark, --knot/--link, --latex, --estimate

Branch `feat-cli` from `origin/main` (0e2ea24). New code lives in `src/bin/krhom/cli_ext.rs`
(module `cli_ext`, included by `src/bin/krhom.rs` via `#[path]`) and in a new section of
`src/calc.rs` (`EstOut`, `estimate_prime`; `run_prime` is unchanged). `krhom.rs` gets only
hooks: two separate option arrays (`CLI_VALUE_OPTS`, `CLI_FLAG_OPTS`; the existing
`VALUE_OPTS`/`FLAG_OPTS` lines are untouched), three `Opts` fields, the `Input::Table`
variant, an `Outcome::estimate` field, one-line calls at the route points, the usage text and
the `--latex` writer in `main`. Engines, certification checks and library stamps are unchanged.

## What was implemented

1. **Refusals before engine work (exit 1).** Soergel and cable routes with the tower closure:
   more than `tower::MAX_VARS` = 10 strands (cable: the cable's strand count), or N > 32
   (`Caps::new` truncation cap N − 1 ≤ 31); the message suggests `--closure exact`, which has
   no N cap (N = 33 trefoil runs in 0.02 s). Exact 2-strand route: dim A_λ = ∏ C(N, a_i) over
   the unmarked strands > 5·10⁷, computed as `cclose::alg_dim` does (unit test: C(36,3) =
   C(120,2) = 7140, so color 3 at N ≥ 36 and color 2 at N ≥ 120 unreduced); `exact2_err` also maps
   the oracle's "too large for the exact oracle" error to exit 1.
2. **`--mark K`** (needs `--reduced`; knots accept only 0). The braid is conjugated so that
   component K is at bottom position 0: a cyclic rotation of the word when one does it (same
   length), else (σ1⋯σp) β (σp⁻¹⋯σ1⁻¹) freely reduced (no cyclic reduction, which would undo
   it). Components keep their ids, so `--colors` stays in krhom's numbering. The output shows
   `conjugated braid: ...` and the method line `reduced (component K marked ...)`.
3. **`--knot NAME` / `--link NAME`** read `data/knots_upto12.jsonl` / `data/links_upto11.jsonl`
   (`$KRHOM_DATA_DIR` overrides the source-tree path). Oriented link names (`L2a1{1}`); an
   unoriented name lists its orientations; unknown names list up to 5 close matches
   (edit distance); a link name given to `--knot` (or vice versa) is redirected. Header line
   `source: knot 4_1 from data/knots_upto12.jsonl (KnotInfo braid; KnotInfo chirality)`.
4. **`--latex OUT`**: a complete `article` document (graphicx only), one `table`/`tabular`
   per N (rows t, columns q, blank = 0, wide tables scaled to the line width), caption with
   the knot/braid, colors, N, reduced/unreduced (marked component), coefficients; with
   `--file` one `\section*` per input (errors noted). Rewritten after every input.
5. **`--estimate`**: on the box route, after the library loads, per prime: record scan, then the
   determinant pass and replay (`swire::generator_pass` + `swire::scan_box_point`, the
   single-point path of `sbox point`) of the middle needed box point, timed. Prints the number
   of box points and an extrapolated wall time and peak memory at the current thread count;
   the closure phase is not included. Other routes: `fast route (ROUTE), no estimate needed`.
   Nothing is computed; exit 0.

## Tests

`KRHOM_THREADS=4 cargo test --release` (default suite): **314 passed, 0 failed, 42 ignored**
(305 before; new: 6 integration tests in tests/calc.rs — `fcli_limits_refused_before_engine_work`,
`fcli_mark_component`, `fcli_knot_and_link_names`, `fcli_tables_parse_completely`, `fcli_latex`,
`fcli_estimate_fast_routes` — and 3 unit tests in `cli_ext`). One new ignored library test,
`lib_fcli_knot_name_and_estimate` (not run as a test; the same commands were run by hand below).

`fcli_mark_component` checks for the Hopf link with colors (2,1), N = 2, 3, 4:
χ(unreduced) = χ(reduced, K marked) · [N choose c_K]_q for K = 0, 1; `--mark 1` equals
`--colors 1,2 --reduced`; `--mark 0` equals plain `--reduced`; N = 2 by hand (Λ²C² is a line,
so marking the color-2 component gives the unreduced q⁻¹ + q⁻³; marking the color-1 one gives
q⁻²); ranks 2/9/24 = (2·1, 3·3, 4·6). The 3-component chain `1 1 2 2` satisfies the Euler
identity for every mark, and the Borromean rings `L6a4{0,0}` give the same result for every mark.

LaTeX compilation was **not** checked: no TeX installation in the container (installing
texlive was not permitted here). The file is plain `article` + `graphicx`; the test checks
the structure (balanced environments, escaping, the expected table rows).

## Validation outputs

Libraries installed with `tools/fetch_library.sh l222_n4 l2222_n4` (both primes, exit 0).

### `krhom 2,3,4 --knot 4_1 --colors 2 --reduced` (KRHOM_THREADS=4)

```
name: 4_1
source: knot 4_1 from data/knots_upto12.jsonl (KnotInfo braid; KnotInfo chirality)
gl_N for N in [2, 3, 4]  strands 3  crossings 4  components 1  colors [2]
braid: 1 -2 1 -2
method: singular (localized box engine; library l222_n4 @ data ee548d2, labels (2,2,2), prune 4, box N 4, seed 201), reduced (strand at bottom position 0 marked)
coefficients: Q
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
exit 0
```

`krhom 2,3,4 --braid "1 -2 1 -2" --colors 2 --reduced` prints the same N blocks (diff of the
N/Poincaré/rank/Euler lines: identical; 25.0 s vs 25.3 s).

### `krhom 2,3,4 --knot 6_1 --colors 2 --reduced --estimate`

```
name: 6_1
source: knot 6_1 from data/knots_upto12.jsonl (KnotInfo braid; KnotInfo chirality)
gl_N for N in [2, 3, 4]  strands 4  crossings 7  components 1  colors [2]
braid: 1 1 2 -1 -3 2 -3
method: singular (localized box engine; library l2222_n4 @ data ee548d2, labels (2,2,2,2), prune 4, box N 4, seed 201), reduced (strand at bottom position 0 marked)
estimate: F_2147483647: 216 box points (reduced grid), library 199s + record scan 32s + 216 x 33.7s per point (replay + determinant pass) / 4 threads + scaling check -> about 34m37s wall, peak memory about 9518 MB
estimate: F_4294967291: 216 box points (reduced grid), library 190s + record scan 31s + 216 x 35.2s per point (replay + determinant pass) / 4 threads + scaling check -> about 35m49s wall, peak memory about 9627 MB
estimate: total: about 1h10m wall for 2 primes at 4 threads, peak memory about 9627 MB (an estimate from one timed point per prime; the closure phase is not included and adds minutes on 4 strands; nothing was computed)
exit 0

[F_2147483647] estimate: library 198.6s, record scan 32.3s, one box point 26.1s (determinant pass 7.65s, peak 1223 MB), VmHWM 2597 MB
[F_4294967291] estimate: library 190.3s, record scan 31.5s, one box point 26.8s (determinant pass 8.39s, peak 1223 MB), VmHWM 2709 MB
```

The measured 26.1 s / 26.8 s per point (single point, 4 threads available) agree with the
C2-7 production runs of 6_1 (26.3–27.1 s per point, docs/CALCULATOR.md §6). The memory figure
is deliberately pessimistic (VmHWM + one point's engine peak per extra thread).

### Refusals (exit 1, before any engine work)

```
$ krhom 33 --braid aaa
exit 1
krhom: N = 33: the tower closure supports N <= 32 (its truncation cap N - 1 must be <= 31 in the packed monomials); use --closure exact for larger N

$ krhom 36 --braid aaa --colors 3
exit 1
krhom: N = 36 with labels [3, 3]: the exact 2-strand engine's algebra A_lambda has dimension 50979600 (product of C(N, a) over the unmarked strands), above its limit 5e7 (the exact oracle stores dense tables of that size); e.g. color 3 needs N <= 35 and color 2 N <= 119 unreduced

$ krhom 2 --braid 1,2,3,4,5,6,7,8,9,10
exit 1
krhom: 11 strands: the tower closure of the Soergel engine supports at most 10 strands (one packed variable per strand, src/uncolored/tower.rs MAX_VARS); about 5 strands are practical

$ krhom 120 --braid aaa --colors 2
exit 1
krhom: N = 120 with labels [2, 2]: the exact 2-strand engine's algebra A_lambda has dimension 50979600 (product of C(N, a) over the unmarked strands), above its limit 5e7 (the exact oracle stores dense tables of that size); e.g. color 3 needs N <= 35 and color 2 N <= 119 unreduced

$ krhom 2 --knot 4_2
exit 1
krhom: --knot 4_2: unknown knot name (data/knots_upto12.jsonl); close matches: 4_1, 5_2, 6_2, 7_2, 8_2

$ krhom 2 --link L2a1
exit 1
krhom: --link L2a1: give an oriented name; LinkInfo lists 2 orientations: L2a1{0}, L2a1{1}

```

Before this branch `krhom 33 --braid aaa` panicked in the engine (exit 2, "truncation cap 32
too large") and `krhom 36 --braid aaa --colors 3` exited 5 (memory cap).

### `--mark` (exact 2-strand route and Soergel route)

```
$ krhom 2,3,4 --braid "1 1" --colors 2,1
method: singular (exact 2-strand singular Soergel engine)
N = 2
  Poincare polynomial: t^0q^-3 + t^0q^-1
  total rank: 2
N = 3
  Poincare polynomial: t^0q^-6 + 2t^0q^-4 + 2t^0q^-2 + t^0q^0 + t^2q^-10 + t^2q^-8 + t^2q^-6
  total rank: 9
N = 4
  Poincare polynomial: t^0q^-9 + 2t^0q^-7 + 3t^0q^-5 + 3t^0q^-3 + 2t^0q^-1 + t^0q^1 + t^2q^-15 + 2t^2q^-13 + 3t^2q^-11 + 3t^2q^-9 + 2t^2q^-7 + t^2q^-5
  total rank: 24

$ krhom 2,3,4 --braid "1 1" --colors 2,1 --reduced --mark 0
method: singular (exact 2-strand singular Soergel engine), reduced (component 0 marked: --mark 0, at bottom position 0)
N = 2
  Poincare polynomial: t^0q^-3 + t^0q^-1
  total rank: 2
N = 3
  Poincare polynomial: t^0q^-4 + t^0q^-2 + t^2q^-8
  total rank: 3
N = 4
  Poincare polynomial: t^0q^-5 + t^0q^-3 + t^2q^-11 + t^2q^-9
  total rank: 4

$ krhom 2,3,4 --braid "1 1" --colors 2,1 --reduced --mark 1
conjugated braid: 1 1  (component 1 moved from bottom position 1 to 0; same closure)
method: singular (exact 2-strand singular Soergel engine), reduced (component 1 marked: --mark 1, at bottom position 0)
N = 2
  Poincare polynomial: t^0q^-2
  total rank: 1
N = 3
  Poincare polynomial: t^0q^-4 + t^0q^-2 + t^2q^-8
  total rank: 3
N = 4
  Poincare polynomial: t^0q^-6 + t^0q^-4 + t^0q^-2 + t^2q^-12 + t^2q^-10 + t^2q^-8
  total rank: 6

$ krhom 2,3,4 --braid "1 1" --colors 1,2 --reduced
method: singular (exact 2-strand singular Soergel engine), reduced (strand at bottom position 0 marked)
N = 2
  Poincare polynomial: t^0q^-2
  total rank: 1
N = 3
  Poincare polynomial: t^0q^-4 + t^0q^-2 + t^2q^-8
  total rank: 3
N = 4
  Poincare polynomial: t^0q^-6 + t^0q^-4 + t^0q^-2 + t^2q^-12 + t^2q^-10 + t^2q^-8
  total rank: 6

$ krhom 2,3 --braid "1 1 2 2" --reduced --mark 0
closure: tower (Right end first), reduced (component 0 marked: --mark 0, at bottom position 0)
  Poincare polynomial: t^0q^-2 + 2t^2q^-6 + t^4q^-10
  Poincare polynomial: t^0q^-4 + 2t^2q^-10 + 2t^2q^-8 + t^4q^-16 + 2t^4q^-14 + t^4q^-12

$ krhom 2,3 --braid "1 1 2 2" --reduced --mark 1
conjugated braid: 1 2 2 1  (component 1 moved from bottom position 1 to 0; same closure)
closure: tower (Right end first), reduced (component 1 marked: --mark 1, at bottom position 0)
  Poincare polynomial: t^0q^-2 + 2t^2q^-6 + t^4q^-10
  Poincare polynomial: t^0q^-4 + 2t^2q^-10 + 2t^2q^-8 + t^4q^-16 + 2t^4q^-14 + t^4q^-12

$ krhom 2,3 --braid "1 1 2 2" --reduced --mark 2
conjugated braid: 1 2 1 1 2 -1  (component 2 moved from bottom position 2 to 0; same closure)
closure: tower (Right end first), reduced (component 2 marked: --mark 2, at bottom position 0)
  Poincare polynomial: t^0q^-2 + 2t^2q^-6 + t^4q^-10
  Poincare polynomial: t^0q^-4 + 2t^2q^-10 + 2t^2q^-8 + t^4q^-16 + 2t^4q^-14 + t^4q^-12

$ krhom 2,3 --braid "1 1 2 2"
  Poincare polynomial: t^0q^-3 + t^0q^-1 + 2t^2q^-7 + 2t^2q^-5 + t^4q^-11 + t^4q^-9
  Poincare polynomial: t^0q^-6 + t^0q^-4 + t^0q^-2 + 2t^2q^-12 + 4t^2q^-10 + 4t^2q^-8 + 2t^2q^-6 + t^4q^-18 + 3t^4q^-16 + 4t^4q^-14 + 3t^4q^-12 + t^4q^-10

$ krhom 2,3 --link "L6a4{0,0}" --reduced --mark 0
source: link L6a4{0,0} from data/links_upto11.jsonl (LinkInfo braid; LinkInfo orientation as named)
  Poincare polynomial: t^-3q^6 + 3t^-2q^4 + 2t^-1q^2 + 4t^0q^0 + 2t^1q^-2 + 3t^2q^-4 + t^3q^-6
  Poincare polynomial: t^-3q^8 + t^-2q^4 + 2t^-2q^6 + 2t^-1q^2 + 2t^0q^-2 + 5t^0q^0 + 2t^0q^2 + 2t^1q^-2 + 2t^2q^-6 + t^2q^-4 + t^3q^-8

$ krhom 2,3 --link "L6a4{0,0}" --reduced --mark 1
source: link L6a4{0,0} from data/links_upto11.jsonl (LinkInfo braid; LinkInfo orientation as named)
conjugated braid: -2 1 -2 1 -2 1  (component 1 moved from bottom position 1 to 0; same closure)
  Poincare polynomial: t^-3q^6 + 3t^-2q^4 + 2t^-1q^2 + 4t^0q^0 + 2t^1q^-2 + 3t^2q^-4 + t^3q^-6
  Poincare polynomial: t^-3q^8 + t^-2q^4 + 2t^-2q^6 + 2t^-1q^2 + 2t^0q^-2 + 5t^0q^0 + 2t^0q^2 + 2t^1q^-2 + 2t^2q^-6 + t^2q^-4 + t^3q^-8

$ krhom 2,3 --link "L6a4{0,0}" --reduced --mark 2
source: link L6a4{0,0} from data/links_upto11.jsonl (LinkInfo braid; LinkInfo orientation as named)
conjugated braid: -2 1 -2 1 -2 1  (component 2 moved from bottom position 2 to 0; same closure)
  Poincare polynomial: t^-3q^6 + 3t^-2q^4 + 2t^-1q^2 + 4t^0q^0 + 2t^1q^-2 + 3t^2q^-4 + t^3q^-6
  Poincare polynomial: t^-3q^8 + t^-2q^4 + 2t^-2q^6 + 2t^-1q^2 + 2t^0q^-2 + 5t^0q^0 + 2t^0q^2 + 2t^1q^-2 + 2t^2q^-6 + t^2q^-4 + t^3q^-8

```

### `--latex` example: `krhom 2,3 --knot 3_1 --colors 2 --reduced --latex t31.tex`

```latex
% written by krhom --latex; compile with pdflatex (needs only graphicx)
\documentclass{article}
\usepackage{graphicx}
\begin{document}
% ---- krhom tables begin
\section*{knot 3\_1}
\begin{table}[ht]
\centering
\resizebox{\ifdim\width>\linewidth\linewidth\else\width\fi}{!}{%
\begin{tabular}{r|r}
$t \backslash q$ & $0$ \\
\hline
$0$ & 1 \\
\end{tabular}}
\caption{knot 3\_1: braid $1\,1\,1$ on 2 strands, colors [2], $N = 2$, reduced (strand at bottom position 0 marked); coefficients Q; rows $t$, columns $q$ (krhom's grading), total rank 1.}
\end{table}
\begin{table}[ht]
\centering
\resizebox{\ifdim\width>\linewidth\linewidth\else\width\fi}{!}{%
\begin{tabular}{r|rrrrr}
$t \backslash q$ & $-12$ & $-10$ & $-8$ & $-6$ & $-4$ \\
\hline
$0$ &  &  &  &  & 1 \\
$1$ &  &  &  &  &  \\
$2$ &  &  & 1 &  &  \\
$3$ & 1 &  &  &  &  \\
\end{tabular}}
\caption{knot 3\_1: braid $1\,1\,1$ on 2 strands, colors [2], $N = 3$, reduced (strand at bottom position 0 marked); coefficients Q; rows $t$, columns $q$ (krhom's grading), total rank 3.}
\end{table}
\clearpage
% ---- krhom tables end
\end{document}
```
