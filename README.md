# KRcalculator (krhom): colored Khovanov–Rozansky gl_N homology

`krhom` computes the Poincaré polynomial of undeformed Khovanov–Rozansky gl_N homology of an oriented link. The link is given as a braid word or a PD code. Each component is either uncolored or colored by an exterior power ∧^a C^N. At N = 2 the result is Khovanov homology with q ↦ q⁻¹. Symmetric colorings (Sym^a) and colored HOMFLY-PT homology are not computed. One command picks the engine, runs it modulo two large primes, and prints the Poincaré polynomial, total rank and Euler characteristic for every N in a list.

In practice it can compute:
- **Uncolored:** any N ≤ 32 on braids with up to about 5 strands. Table knots take seconds to minutes.
- **∧^a with a ≤ 3 on 2-strand braids** (the torus knots and links T(2,n)). The two components of a link may have different colors. Tested up to N = 10 for ∧² and N = 8 for ∧³; memory limits larger N.
- **∧² on 3-strand braids, N ≤ 4:** minutes to hours, and up to about 10 GB.
- **∧² on 4-strand braids, N ≤ 4 (experimental):** only words the library covers. These take hours per prime, so use `--reduced`.

For a knot colored by ∧², N = 2 gives `t^0q^0` (∧²C² is one-dimensional), and N = 3 equals the uncolored gl_3 homology of the knot with reversed orientation (a ↔ N−a duality). On 3 and 4 strands, N = 4 is therefore the only new value. A result is printed only if the two primes agree. [Correctness](#correctness) explains what this guarantees and what it does not.

**Authorship.** All code in this repository was written by Claude (Anthropic's AI model). Hongjian Yang (Stanford University) provided the mathematical direction.

## Quick start

```sh
git clone https://github.com/Exista-Yang/KRcalculator.git && cd KRcalculator
cargo build --release                     # edition 2024, so Rust 1.85 or newer; binaries go to target/release/
tools/fetch_library.sh l222_n4 l2222_n4   # only for colored braids on 3-4 strands; run once
tools/smoke.sh                            # checks known values in a few seconds (more with the libraries)
```

**What can I compute, and how long does it take?**

| Input | Example | Typical time |
|---|---|---|
| uncolored, ≤ 5 strands, any N ≤ 32 | `krhom 2,3,4 --knot 8_19` | seconds (12-crossing 4-strand knots: up to minutes) |
| colored, 2 strands (T(2,n)), colors ≤ 3 | `krhom 2,3,4 --braid aaa --colors 2` | ∧²: < 1 s; ∧³: seconds to an hour |
| color 2, 3 strands, N ≤ 4 (library `l222_n4`) | `krhom 2,3,4 --knot 4_1 --colors 2` | minutes; large 12-crossing knots up to hours |
| color 2, 4 strands, N ≤ 4 (library `l2222_n4`) | `krhom 2,3,4 --knot 6_1 --colors 2 --reduced --checkpoint run61` | hours per prime: check first with `--estimate` |
| anything else colored (mixed colors on ≥ 3 strands, colors ≠ 2 on 3–4 strands, ≥ 5 strands) | — | refused with a message (exit 1) |

Use `--estimate` before a 3–4-strand run to see the number of evaluation points and an extrapolated time, `--checkpoint DIR` so that an interrupted run continues where it stopped, and `--reduced` to cut the cost of 4-strand runs about sixfold.

The crate depends only on `num-bigint`, `num-integer`, `num-rational`, `num-traits` and `rand`.

**Fetching libraries.** Run `tools/fetch_library.sh` from inside the clone. It needs bash, git, tar with gzip, `sha256sum` or `shasum`, and the usual awk, sed, grep, find and mktemp.
- It fetches data commit `ee548d2` of the orphan branch `data-gl6` from the repository you cloned (`--repo URL` selects another source).
  - It downloads only the library chunks of that commit (a partial fetch of about 40 MB); if the server refuses partial fetches it falls back to a full fetch.
  - `--url BASE` instead downloads release assets over HTTPS (prepared by `tools/dev/make_release_assets.sh`) and verifies them against `data/`.
  - The branch is named after the parked (3,3,3) gl_6 work it also stores. The two libraries krhom uses are for N ≤ 4.
- It installs one directory per library and prime: 21 MB per prime for `l222_n4` and 68–75 MB per prime for `l2222_n4`.
  - The install root is `--root DIR`, else `KRHOM_LIB_DIR`, else `${XDG_CACHE_HOME:-$HOME/.cache}/krhom/lib`.
- A library built for N = 4 serves every N ≤ 4. The summands it drops vanish at every N′ ≤ 4 (see [Pruning](#how-it-works)).
- `tools/fetch_library.sh --verify-only l222_n4 l2222_n4` re-checks every installed file against its sha256.

**Examples.** The outputs below are verbatim from the recorded runs (directory `runs_calc/` on branch `archive-pre-restructure`), made on a 4-core Linux cloud machine. The 3-strand example used `KRHOM_THREADS=4`. Here `krhom` means `target/release/krhom`.

Uncolored figure-eight knot, gl_2 and gl_3:
```
$ krhom 2,3 --braid aBaB
gl_N for N in [2, 3]  strands 3  crossings 4  components 1  colors [1]
braid: 1 -2 1 -2
closure: tower (Right end first)
coefficients: Q
time 0.19s
N = 2
  Poincare polynomial: t^-2q^5 + t^-1q^1 + t^0q^-1 + t^0q^1 + t^1q^-1 + t^2q^-5
  total rank: 6
  Euler characteristic: q^5 + q^-5
N = 3
  Poincare polynomial: t^-2q^6 + t^-2q^8 + t^-1q^0 + t^-1q^2 + t^0q^-2 + t^0q^0 + t^0q^2 + t^1q^-2 + t^1q^0 + t^2q^-8 + t^2q^-6
  total rank: 11
  Euler characteristic: q^8 + q^6 - 1 + q^-6 + q^-8
```
The `closure:` line names the closure algorithm of the uncolored engine: the restriction tower, starting from the right end. `--closure exact` gives the same result.

∧²-colored right-handed trefoil on 2 strands:
```
$ krhom 2,3,4 --braid aaa --colors 2
gl_N for N in [2, 3, 4]  strands 2  crossings 3  components 1  colors [2]
braid: 1 1 1
method: singular (exact 2-strand singular Soergel engine)
coefficients: Q
time 0.09s
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
"Singular" in the `method:` line refers to singular Soergel bimodules, not to singular links. N = 2 is `t^0q^0`, as expected. N = 3 has the same ranks as the uncolored gl_3 homology of the trefoil (duality).

The same knot, reduced (excerpt). The route line says which strand is marked:
```
$ krhom 2,3,4 --braid "1 1 1" --strands 2 --colors 2 --reduced
...
method: singular (exact 2-strand singular Soergel engine), reduced (strand at bottom position 0 marked)
...
N = 4
  Poincare polynomial: t^0q^-8 + t^2q^-14 + t^2q^-12 + t^3q^-18 + t^3q^-16 + t^4q^-20 + t^5q^-24 + t^5q^-22 + t^6q^-26
  total rank: 9
```
The reduced ranks are 1, 3 and 9 at N = 2, 3 and 4.

∧²-colored figure-eight on 3 strands, using library `l222_n4` (`--braid aBaB --colors 2` is the same input):
```
$ krhom 2,3,4 --braid "1 -2 1 -2" --strands 3 --colors 2
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
  Poincare polynomial: t^-4q^16 + t^-3q^8 + t^-3q^10 + t^-2q^2 + t^-2q^6 + 2t^-2q^8 + 2t^-2q^10 + t^-2q^12 + t^-1q^-2 + 2t^-1q^0 + 2t^-1q^2 + t^-1q^4 + t^-1q^6 + t^-1q^8 + t^0q^-4 + 2t^0q^-2 + 4t^0q^0 + 2t^0q^2 + t^0q^4 + t^1q^-8 + t^1q^-6 + t^1q^-4 + 2t^1q^-2 + 2t^1q^0 + t^1q^2 + t^2q^-12 + 2t^2q^-10 + 2t^2q^-8 + t^2q^-6 + t^2q^-2 + t^3q^-10 + t^3q^-8 + t^4q^-16
  total rank: 46
  Euler characteristic: q^16 + q^12 + q^10 + q^-10 + q^-12 + q^-16
```
In the `method:` line, `prune 4` and `box N 4` mean the computation is done once at N = 4 and serves every N ≤ 4. `seed` is the seed of the random evaluation points. Per-prime check summaries go to stderr and end with:
```
krhom: F_2147483647 and F_4294967291 agree: certified over both primes (equal over Q unless both are torsion primes)
```
With `--reduced` the ranks are 1, 5 and 25. The N = 3 block equals the uncolored gl_3 block of 4_1 above, as duality predicts.

**Computing a table knot.** Every record in `data/knots_upto12.jsonl` has a braid in KnotInfo's chirality (see [Conventions](#conventions)). With `jq` installed:
```sh
w=$(jq -r 'select(.name == "5_2") | .braid.word | map(tostring) | join(" ")' data/knots_upto12.jsonl)
target/release/krhom 2,3,4 --braid "$w" --colors 2
```

## Conventions

- **Braid words.** A word is read from bottom to top, and all strands are oriented upward, so every component of the closure runs the same way. For other orientations of a link, use a PD code.
- **Positions.** Strand positions at the bottom are numbered 0, 1, …, n−1 from the left. σ_k (written `k`, or the k-th letter) crosses positions k−1 and k. In a positive crossing (`+k`) the left strand passes over, so `aaa` and `"1 1 1"` are the right-handed trefoil.
- **Chirality.**
  - The braids and PD codes in `data/*.jsonl` have KnotInfo/LinkInfo's chirality.
  - From those braids, krhom reproduces the tables' Khovanov data at N = 2 after q ↦ q⁻¹. This was checked for every record with ≤ 4 strands.
  - Knot Atlas / KnotTheory PD codes are often the mirror image. For 139 knots up to 10 crossings, among them 3_1, 5_2 and 6_1, the Knot Atlas knot is the mirror of KnotInfo's (`data/README.md`).
- **q-direction.** Positive knots sit in negative q-degree, as in Khovanov–Rozansky [1] and Lewark–Lobb [11].
  - The right-handed trefoil at N = 2 is `t^0q^-3 + t^0q^-1 + t^2q^-5 + t^3q^-9`.
  - Khovanov tables (Bar-Natan, KnotInfo/LinkInfo, KnotJob) use the opposite direction: krhom's bidegree (h, q) is the table's (h, −q). Bar-Natan's trefoil is q + q^3 + t^2q^5 + t^3q^9.
- **Crossing shifts.** Take σ_k^{±1} with colors c (left) and d (right) at the bottom, and set m = min(c, d) and x = m(N − m).
  - The colored Rickard complex [3, Def. 3.3] has terms W_0, …, W_m, where W_k = F^{(c−d+k)}E^{(k)}1_{(c,d)} for c ≥ d, mirrored for c < d.
  - For a positive crossing, W_k sits in homological degree k with q-shift q^{−(x+k)}.
  - A negative crossing is the vertical flip, in degree −k with q^{x+k}.
  - Shifts multiply graded dimensions by q^s, and webs are MOY-normalized.
  - For c = d this choice is invariant under Reidemeister I. It differs from ETW's q^{k−x} by a framing factor per crossing.
  - Mixed colors use the same formula.
- **Reduced homology.** `--reduced` marks the strand at bottom position 0.
  - In the closure, that strand's component carries a one-dimensional algebra in place of H*(Gr(a, N)), where a is its color, and its shift a(N − a) is dropped. The reduced ∧^a unknot is therefore `t^0q^0`.
  - For knots the result does not depend on the marked point (checked).
  - For PD input of a link, the marked component is whichever component sits at position 0 of the braid on the `braid:` line, and it cannot be chosen directly. To mark a particular component, pass a braid in which that component starts at position 0.
- **Unknot.** The unreduced ∧^a unknot is [N choose a]_q (uncolored: [N]). The reduced unknot is `t^0q^0`.
- **Euler characteristic.** χ = Σ (−1)^h k q^j.
  - Uncolored: χ = [N]·P(v = q^−N, z = q − q^−1), where P is HOMFLY-PT in the KnotInfo/Morton convention. Reduced χ = P(v = q^−N, z = q − q^−1).
  - Colored: unreduced χ is the ∧-colored MOY (quantum gl_N) evaluation of the closure with the crossing shifts above. For knots, χ(unreduced) = [N choose a]_q · χ(reduced).
- **gl_N.** This is the gl_N version of the theory. At N = 2, uncolored results agree with the Khovanov homology tables for both knots and links.
- **Mirror.** Mirroring maps (t, q, k) ↦ (−t, −q, k).
- **Reverse and duality.**
  - rev(w) is the word w read backwards. Its closure is the closure of w with the orientation reversed, which is the same knot unless the knot is non-invertible (such as 8_17).
  - Duality, expected [5] and used only as a check: ∧^a at N = a+1 equals uncolored gl_{a+1} of rev(w).
- **Wang's formulas.** ours(t, q) = Wang(1/t, 1/q) (`tools/wang.py`).
- **Reference values.**
  - The right-handed trefoil (above).
  - The positive Hopf link `"1 1"` at N = 2 is `t^0q^-2 + t^0q^0 + t^2q^-6 + t^2q^-4`. This is LinkInfo L2a1{1} with q ↦ q⁻¹ and is covered by the table check.
  - `"1 2 3" --colors 2 --no-simplify` gives [N choose 2]_q, with ranks 1, 3 and 6 at N = 2, 3 and 4.

## Usage

```
krhom N[,N...] (--braid WORD [--strands n] | --knot NAME | --link NAME | --pd PD | --file BRAIDS)
               [--colors c | c1,c2,...] [--reduced [--mark K]] [--estimate] [--checkpoint DIR]
               [--table] [--latex OUT] [--json OUT] [--csv OUT] [--no-simplify] [--check-library] [--lib-dir DIR]
               [--seed S] [--seed-retries K] [--mem-mb M] [--mem-cap-mb M]
               [--colored-method cable|cube|singular] [--cube] [--closure tower|exact] [--end right|left|alt]
```

- **Knot and link names.** `--knot 6_1` and `--link L2a1{1}` take the braid from `data/knots_upto12.jsonl` / `data/links_upto11.jsonl` (KnotInfo/LinkInfo chirality and orientations). Unknown names list close matches.
- **`--mark K`.** With `--reduced`, mark link component K (0-based, in the order of `--colors`) instead of the component at bottom position 0; krhom conjugates the braid, which does not change the link.
- **`--estimate`.** On the 3–4-strand route, time one evaluation point per prime and print the number of points and an extrapolated time and memory; nothing is computed.
- **`--checkpoint DIR`.** On the 3–4-strand route, save progress in DIR. Rerunning the same command after an interruption continues where it stopped and prints the same result as an uninterrupted run; a different input with the same DIR is refused.
- **`--latex OUT`.** Write the bigraded tables as a standalone LaTeX document (one tabular per N).
- **N list.** The first argument, for example `2,3,4`. Each value is an integer ≥ 1. The list is sorted and duplicates are removed.
- **Braid words.**
  - Signed integers, for example `"1 -2 1 -2"`, separated by spaces or commas.
  - KnotJob letters, for example `aBaB`: a..i are σ1..σ9 and A..I their inverses. Write the letters without spaces.
  - A word that mixes the two forms is refused.
  - `--strands n` defaults to max |generator| + 1. Extra strands become split unknots.
- **PD codes.**
  - Accepted forms: `--pd "[[1,4,2,5],[3,6,4,1],[5,2,6,3]]"`, `X[1,4,2,5] X[3,6,4,1] ...` and `PD[X[...], ...]`.
  - The convention is KnotTheory's: `X[i,j,k,l]` lists edge labels counterclockwise, starting from the incoming under-strand.
  - The code is validated, converted to a braid by Vogel's algorithm [24] and simplified. The resulting braid is printed on the `braid:` line.
- **Colors.**
  - `--colors 2` colors every component with 2. A list such as `--colors 2,1` gives one color per component. The default is all 1 (uncolored).
  - For a braid, components are numbered by the smallest bottom position they pass through. For PD input, the component that contains the smallest edge label is component 0.
  - If N is below a color, the homology is 0, and krhom prints the zero polynomial with a note.
- **`--reduced`.** See [Conventions](#conventions). Each run prints one mode, reduced or unreduced. The flag works on the uncolored and singular colored routes, but not with `--colored-method cube|cable`.
- **Simplification of colored braids.** krhom tries cyclic reduction, Markov destabilization and a deterministic random search.
  - It uses the result only if it has fewer strands, and, with `--reduced` on a link, only if the marked component is unchanged.
  - For example, `abc --colors 2 --reduced` becomes the 1-strand unknot.
  - `--no-simplify` keeps the given word.
- **`--check-library`.** For a colored word on 3–4 strands, it loads and re-certifies the library maps the word needs and reports whether the library covers the word. It does not run the box scan. For a 6_1-sized word this takes about 7 minutes and 1 GB.
- **`--table`.** After each N block, also print the bigraded table: rows t ascending, columns q ascending, entries the ranks (blank for 0). Example (∧² trefoil, N = 3):
  ```
  table (rows t, columns q):
  t\q  -14  -12  -10   -8   -6   -4   -2
    0                        1    1    1
    1
    2                   1    1
    3    1    1
  ```
- **Batch, JSON and CSV.**
  - `--file BRAIDS` takes one braid word per line (PD codes are not accepted), with an optional `NAME;` prefix. `#` starts a comment.
  - N, colors, `--strands` and `--reduced` apply to every line.
  - Each block starts with `name: NAME`. A failing line prints `name: NAME` and `error: exit K (kind)`, and the batch continues. The exit code is that of the first failing line.
  - With `--file`, the JSON and CSV files are rewritten after every finished line, so an aborted batch keeps its finished results.
  - `--json OUT` writes one object per input (an array with `--file`). It holds the input, route, library, seeds, per-prime checks, the agreement flag and the results, with terms in `results[].P.terms` as `[t,q,rank]`.
  - `--csv OUT` writes rows `name,N,reduced,t,q,rank` (names are quoted as in RFC 4180).
  - [docs/CALCULATOR.md](docs/CALCULATOR.md) lists every field.
- **Environment.**
  - `KRHOM_THREADS` sets the thread count; the default is all cores.
  - The library root is `--lib-dir DIR`, else `KRHOM_LIB_DIR`, else `${XDG_CACHE_HOME:-$HOME/.cache}/krhom/lib`.
  - `KRHOM_MEM_CAP_MB` sets the closure memory cap.
  - On Linux, run long jobs under `ulimit -v 14000000`; this is the real memory guard.
- **Advanced options.**
  - `--colored-method cube`, or its alias `--cube`, computes the full cube: one exact run over Q, exponential in the number of crossings, for tiny inputs only.
  - `--colored-method cable` cables the link through the uncolored engine.
  - `--closure tower|exact` and `--end right|left|alt` choose the uncolored closure.
  - `--seed`, `--seed-retries`, `--mem-mb` and `--mem-cap-mb` tune the 3–4-strand route.

**Output layout.** For each input, stdout prints:
1. `name: NAME` (with `--file`);
2. a header (`gl_N for N in [...]`, strands, crossings, components, colors);
3. the `braid:` line, followed by `simplified braid: WORD  (strands k)` if simplification reduced the number of strands;
4. the route line (`closure:` or `method:`; with `--reduced` it ends in `, reduced (strand at bottom position 0 marked)`), then `coefficients: Q` and the wall time;
5. one block per N with the Poincaré polynomial, total rank and Euler characteristic. With `--check-library`, a `check-library:` report is printed instead.

Progress and the per-prime checks go to stderr.

**Polynomial format.** A term `k t^h q^j` means rank k in homological degree h and quantum degree j. A coefficient of 1 is omitted. Terms are sorted by t and then by q, and the zero polynomial prints as `0`. The Euler characteristic is printed in descending powers of q.

| Exit | Meaning |
|---|---|
| 0 | success |
| 1 | usage error or unsupported input (the message lists what is supported) |
| 2 | engine or certification failure, including a caught panic (please report it, with the `--json` output) |
| 3 | the two primes disagree; no result is printed |
| 4 | library not installed, mismatched, damaged, or not covering the word (the message gives the command to run, or what is missing) |
| 5 | resource refusal: per-point memory above `--mem-mb`, or a closure estimate above `--mem-cap-mb` |

## What is supported, and what it costs

| Input (after simplification) | Engine (chosen automatically) | Supported | Relies on |
|---|---|---|---|
| all colors 1 (knots and links) | localized Soergel engine, exact restriction-tower closure | N = 1..32 and ≤ 10 strands (code limits); about 5 strands in practice | two primes; random evaluation points |
| colored, ≤ 2 strands | exact 2-strand singular Soergel engine | colors 1..3, mixed colors allowed; unreduced N ≤ 119 for color 2 and N ≤ 35 for color 3 (exact-oracle size), and within the closure memory cap (`--mem-cap-mb`, default 768 MB; exit 5 names the flag to raise) | two primes; exact arithmetic over each |
| colored, uniform color 2, 3 strands | localized box engine, library `l222_n4` | N ≤ 4; knots and non-split links | two primes; random points; Foam_N kill rule (checked, not proven) |
| colored, uniform color 2, 4 strands (experimental) | localized box engine, library `l2222_n4` | N ≤ 4; only words the library covers | as above |
| any other colored input | refused, exit 1 | mixed colors or colors ≠ 2 on 3–4 strands, including (3,3,3); ≥ 5 strands; N > 4 on 3–4 strands; colors > 3 on 2 strands; split braids on 3–4 strands | fallback: `--colored-method cube` |

**4-strand coverage.**
- The library covers these 18 knots in the braids of `tools/dev/ntrunc_t1/data/knots4_cloud.txt`: 6_1, 7_2, 7_4, 7_6, 7_7, 8_4, 8_6, 8_8, 8_11, 8_13, 8_14, 9_10, 9_11, 9_13, 9_17, 9_18, 9_20 and 9_23.
- It does not cover 8_15, 9_22, 11n_42 or the Conway knot; these exit with code 4.
- `--check-library` answers the question for any word.

**Typical cost.** "a + b" gives the two primes, F_{2^31−1} + F_4294967291. Unless noted, the runs used krhom on a 4-core Linux cloud machine.

| Input | Time | Peak memory |
|---|---|---|
| uncolored table knots, N = 2..6, unreduced and reduced (one prime, 1 thread) | 3 strands, 12-letter braids: median 1.3 s, max 13 s. 4 strands, 13-letter braids: median 27 s, max 670 s. Conway knot, N = 2..10: 1.24 s unreduced, 1.35 s reduced (laptop, 2 threads). | |
| ∧² trefoil, 2 strands, N = 2, 3, 4 / ∧³ T(2,9), N = 3..8, both modes | 0.09 s in total / 4032–5925 s per prime | T(2,9): 2.8–3.0 GB |
| ∧² 4_1, 3 strands, unreduced / reduced (4 threads) | 43.5 + 38.2 s / 16.5 + 15.8 s | 166 / 160 MB |
| ∧², other 3-strand knots, per prime, both modes (`sbox`, 1–3 threads) | 8 letters: 100–1432 s. 10 letters: 247–8568 s. Longest: 12n_114, 16764 s. | mostly ≤ 0.3 GB; some 10-letter knots 1.4–3.3 GB; 10_104 8.4 GB; 10_99 9.8 GB |
| ∧² 6_1, 4 strands, `--check-library` | 222 + 217 s | 1046 MB |
| ∧² 6_1, 7_6, 7_7, 4 strands, reduced (sharded `sbox`, 4 shards × 1 thread, 4-thread finish), per prime | 216 box points at 26–40 s each; 0.70–1.01 h wall; 2.0–3.1 core-hours of shards (6_1: 2.2, 7_6: 3.0–3.1, 7_7: 2.0–2.9) | 2.7–2.9 GB per shard |
| ∧², 4 strands, unreduced | Not run. Estimate: 1296 box points instead of 216, so about 26.4 s × 1296 ≈ 9.5 core-hours per prime, if unreduced points cost the same as reduced ones. | |

For multi-hour 4-strand words, run `--estimate` first, then use `--reduced` and `--checkpoint DIR` (an interrupted run resumes). For several machines, the sharded `sbox box` / `sbox finish` driver splits the evaluation points.
- The 6_1, 7_6 and 7_7 runs used `tools/dev/box.sh` (with `tools/dev/finish.sh`). It takes the library directory and cost table from the environment (`CACHE`, `COSTS`); see its header.
- [docs/CALCULATOR.md](docs/CALCULATOR.md) §8 describes `sbox`.

Large 3-strand knots can exceed the default memory estimates and exit with code 5. The estimate is usually several times too high, so raise `--mem-cap-mb` and guard the real memory with `ulimit -v`. On macOS, the VmHWM figures on stderr read 0 because they come from `/proc`. More timings are in docs/CALCULATOR.md §6.

## Correctness

### Coefficients and reliability

All ranks are meant as dimensions over Q.

**Two primes.** Every route except the full cube computes modulo p = 2^31 − 1 and p = 4294967291. A result is printed only if the two runs agree exactly; otherwise krhom exits with code 3. In krhom's output, "certified" means this agreement together with the checks listed below. It is not a proof:
- Over F_p a rank is never smaller than over Q. Agreement still leaves the case that both primes are bad for the input: a prime with torsion in the homology, or one at which the precomputed canonical-basis data differ from characteristic 0.
- The uncolored and 3–4-strand routes evaluate maps at seeded random points. The checks at those points hold with high probability, not with certainty.
- Results on 3–4 strands also assume the Foam_N kill rule, which has been checked but not proven.

**Full cube.** `--colored-method cube` runs once, exactly over Q.

No torsion is computed, and small characteristic is not supported.

**Checks on each route.**
- **Uncolored.** Seed 7 for the first prime and seed 11 for the second. The restriction-tower closure is exact.
- **Colored, 2 strands.** All arithmetic is exact over each F_p.
  - Decompositions are certified (Π∘Ι = id).
  - d∘d = 0 is checked after every crossing.
  - The closure is exact.
- **Colored, 3–4 strands.**
  - The library is loaded in replay-only mode. Every loaded map is re-certified (ends, degrees, scaling, Π∘Ι, forms and closure certificates), and a map that fails is rejected.
  - A generator-determinant check runs at every box point. If it fails, the points are redrawn with seed+1, at most twice.
  - A Freivalds test checks D² = 0 on the closed complex for each N′ and mode. The closure refuses pruned data without this certificate.
  - χ is compared with `decat_chi`, an independent decategorified q-Schur computation [10], for every closure. A homogeneity (scaling) check also runs. A failure of either is fatal.
- **Library integrity.**
  - At install, every chunk's sha256 must match both the registry `data/libraries.tsv` and the data commit's MANIFEST.
  - On every run, krhom checks `LIBRARY.txt`, the presence of every file and every file's build stamp, and exits with code 4 on failure.
  - Integrity does not imply soundness, which is why maps are re-certified on load.

### What the computation relies on

- **Theorems:**
  - Williamson's classification of singular Soergel bimodules and his Hom formula [14] (Thm 7.9, with erratum; restated in [15]);
  - singular light leaves [15];
  - the Khovanov–Robert–Wagner equivalence with braid-like foams [20];
  - invariance of the colored Rickard complexes under braid relations and Markov moves ([3], Thm 3.5, in Foam_N; Webster–Williamson [4] for singular Soergel bimodules).
- **Checked computationally, not proven:**
  - **The Foam_N kill rule** (B_p vanishes at gl_N if ν(p)_1 > N). It was validated on 9833 decategorified cases with 0 mismatches and by pruned = unpruned runs. krhom uses it only on 3–4 strands.
  - **Closure dimensions equal MOY evaluations** [8] on hundreds of random ladders, N ≤ 6. The identification of the closure with foam state spaces rests on this check.
  - **The degree-1 Rickard differentials are unique up to scalar** for labels ≤ 5, so they are computed, not hand-coded.
  - **The implementation's braid and Markov invariance** (on Markov and stabilization words) and its base-point independence for knots.
- **Expected, used only as a check:** the a ↔ N−a duality [5].

For uncolored links, the results agree with the Khovanov tables at N = 2 and with HOMFLY-PT at N = 3..6 (below). We have not compared the colored results with the matrix-factorization constructions of Wu [6] and Yonezawa [7].

### Validation

- **Uncolored acceptance tests** (`tests/acceptance.rs`):
  - the Khovanov trefoil and figure-eight;
  - the KR trefoil formula q^{2−2N}[N] + t²q^{−2N−1}[N−1] + t³q^{−4N−1}[N−1] at N = 2, 3, 4, also after Markov stabilization;
  - the reduced trefoil q^{−2(N−1)} + t²q^{−2(N+1)} + t³q^{−4N} for N = 2..10;
  - agreement with an independent full cube on 4 braids at N = 2, 3, and of the tower with the exact closure on 8 braids at N = 2..6;
  - the Conway knot: unreduced ranks 34, 83, 132, …, 426 at N = 2..10. A 2048-vertex full cube at N = 2 gave the identical polynomial.
- **Uncolored tables.** `tools/batch_check.py` ran the engine, over one prime, on every KnotInfo/LinkInfo record with ≤ 4 strands: 1212 knots and 2138 oriented links, with 0 failures (`results/uncolored/SUMMARY.md`).
  - At N = 2 it compares the full unreduced Khovanov polynomial, plus the reduced one for knots.
  - At N = 3..6 for knots and N = 3..5 for links it compares χ with HOMFLY-PT.
  - Records with 5 or more strands were not checked.
- **Colored, 2 strands.**
  - Full bigraded agreement with Wang's T(2,2m+1) formulas [12, 13]: the ∧² trefoil at N = 3..8 and 10, ∧² T(2,5) at N = 3..6, ∧² T(2,7) at N = 3, 4, and ∧³ T(2,3), T(2,5), T(2,7) and T(2,9) at N = 3..8 wherever a Wang block was generated.
  - The colored Hopf links (2,1), (1,2) and (2,2) agree with a closed formula, the full cube and cabling.
  - For ∧³ σ^9, all 96 closures agree over both primes.
- **∧² gl_4.** All 101 computed three-strand knots, and 6_1, 7_6 and 7_7 on 4 strands (reduced), passed:
  - agreement of both primes;
  - the χ and scaling checks;
  - N′ = 3 equal to the uncolored gl_3 homology of the reversed braid, unreduced and reduced;
  - rank 1 at N′ = 2.
- **The same box engine at colors (3,3,3)**, N = 3, 4: all 185 knots of braid index 3 pass over both primes, with χ checked on 740 closures. N = 4 equals the uncolored gl_4 homology of the reversed braid, and N = 3 gives `t^0q^0`. krhom does not expose these colors.
- **4-strand Markov words.** 11 Markov and stabilization words reproduce the 2- and 3-strand values over both primes. Among them are the unknot, the Hopf link, 4_1·σ_3 and 5_2·σ_3, computed unreduced and reduced.
- **krhom runs.**
  - The 4_1 N = 4 output is byte-identical to the earlier `sbox` result.
  - `"1 1 1 2"` and `"1 1 1 2 3"` (with `--no-simplify`) equal the 2-strand trefoil.
  - `"1 2 3"` gives [N choose 2]_q.
- **Random checker.** `tools/randcheck.py` draws random color-2 braids on 2–4 strands and checks fixed relations: the two primes, seed independence, χ, isotopy, Markov moves, mirror, the N = 3 duality, N = 2 = `t^0q^0`, Wang's closed forms, unknots and cross-engine agreement. Its recorded runs passed every case, and an injected fault is caught.
- **Test suite.** `cargo test --release` runs the test suite; `cargo test --release --features dev` adds the tests that run developer binaries. The tests that need the libraries are `#[ignore]`d by default and pass with the libraries installed.

## Results

Representative total ranks follow. Uncolored values come from the table runs; ∧² values come from the box engine and the exact 2-strand engine.

| knot | strands | uncolored, unreduced, N = 2..6 | uncolored, reduced, N = 2..6 | ∧², unreduced, N = 2/3/4 | ∧², reduced, N = 2/3/4 |
|---|---|---|---|---|---|
| 3_1 | 2 | 4 / 7 / 10 / 13 / 16 | 3 | 1 / 7 / 22 (N = 5..8, 10: 46, 79, 121, 172, 301) | 1 / 3 / 9 |
| 4_1 | 3 | 6 / 11 / 16 / 21 / 26 | 5 | 1 / 11 / 46 | 1 / 5 / 25 |
| 5_2 | 3 | 8 / 15 / 22 / 29 / 36 | 7 | 1 / 15 / 78 | 1 / 7 / 49 |
| 8_19 = T(3,4) | 3 | 8 / 19 / 30 / 41 / 52 | 5 at N = 2, then 11 | 1 / 19 / 90 | 1 / 11 / 29 |
| 10_124 = T(3,5) | 3 | 10 / 27 / 44 / 61 / 78 | 7 at N = 2, then 17 | 1 / 27 / 154 | 1 / 17 / 61 |
| 10_99 | 3 | 82 / 163 / 244 / 325 / 406 | 81 | 1 / 163 / 6886 | 1 / 81 / 6561 |
| 6_1 | 4 | 10 / 19 / 28 / 37 / 46 | 9 | not run | 1 / 9 / 81 |
| 7_6 | 4 | 20 / 39 / 58 / 77 / 96 | 19 | not run | 1 / 19 / 361 |
| 7_7 | 4 | 22 / 43 / 64 / 85 / 106 | 21 | not run | 1 / 21 / 441 |
| 11n_34 (Conway) | 4 | 34 / 83 / 132 / 181 / 230 (426 at N = 10) | 33 at N = 2, then 49 | not covered by the library | not covered by the library |

**Example: ∧² 6_1 on 4 strands.** The braid is `1 1 2 -1 -3 2 -3`. The reduced homology at N = 4 has rank 81, with both primes identical and every check passing:
```
t^-4q^12 + t^-3q^8 + t^-3q^10 + t^-2q^4 + 3t^-2q^6 + t^-2q^8 + t^-1q^0 + 4t^-1q^2 + 3t^-1q^4 + 3t^0q^-2 + 6t^0q^0 + t^0q^2 + 2t^1q^-6 + 6t^1q^-4 + 4t^1q^-2 + t^2q^-10 + 4t^2q^-8 + 6t^2q^-6 + t^2q^-4 + 2t^3q^-12 + 5t^3q^-10 + 3t^3q^-8 + t^4q^-16 + 3t^4q^-14 + 4t^4q^-12 + t^4q^-10 + t^5q^-18 + 3t^5q^-16 + 2t^5q^-14 + t^6q^-20 + 2t^6q^-18 + t^7q^-22 + t^7q^-20 + t^8q^-24
```

**∧³-colored T(2,n), unreduced** (exact 2-strand engine, both primes):

| link | N = 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|
| T(2,3) (agrees with Wang) | 1 | 10 | 46 | 136 | 307 | 586 |
| T(2,4) | 1 | 22 | 184 | 864 | 2817 | 7326 |
| T(2,6) | 1 | 28 | 292 | 1568 | 5321 | 13916 |
| T(2,8) | 1 | 34 | 424 | 2560 | 8929 | 23386 |

For ∧³, the reduced ranks of the knots T(2,3), T(2,5), T(2,7) and T(2,9) stabilize at N ≥ 6 to 27, 125, 343 and 729.

**Where the data are.**
- **∧² gl_4 of the three-strand knots** (N′ = 2, 3, 4, unreduced and reduced): `results/colored_gl4/gl4_222.json` and `gl4_222_s{0,1,2}.json`.
- **6_1, 7_6 and 7_7:** `results/colored_gl4/gl4_2222_c27_compare.json`.
- **∧³ σ^9 closures:** `docs/dev/colored_design/gl6/m2s/sigma9_33_F31.json` and `sigma9_33_F32b.json`.
- **Uncolored table runs:** `results/uncolored/`.
- [results/README.md](results/README.md) describes every results file.

**Values not found in the references below.** The following values are not found there and do not follow from a ↔ N−a duality:
- ∧² gl_4 (N = 4) of the 101 three-strand knots, and of 6_1, 7_6 and 7_7;
- ∧³ T(2,2m) at N = 5..8. At N = 3 the rank is 1, and N = 4 follows from duality.

Every other value listed here reproduces a known formula, a table or a duality.

## How it works

**Tangles.** Each crossing becomes its colored Rickard complex ([3], Def. 3.3; Rouquier complexes [2] in the uncolored case).
- Its terms are singular Bott–Samelson bimodules of ladder webs [9] over R^λ = ⊗ Sym(λ_i). These are untruncated singular Soergel bimodules [14], equivalent to braid-like foams [20].
- Each differential is computed as the unique degree-1 map between adjacent terms, up to scalar.
- Nothing at this level depends on N.

**Scan.** As in Bar-Natan's algorithm [23], crossings are added one at a time.
- Terms split into indecomposables B_p, indexed by contingency tables (double cosets).
- Degree-0 isomorphisms are cancelled, so the complex stays minimal.
- Maps are stored localized, one block per coset. They are evaluated at random points over F_p, except on 2 colored strands, where they are exact.

**Libraries (3–4 strands).** The decompositions B_p·rung ≅ ⊕ B_q⟨s⟩ and the closure generators are precomputed per (colors, N, prime) with singular Soergel calculus and light leaves [15–19], and they are re-certified on load.
- A Hom basis is accepted only if its graded dimension matches Williamson's Hom formula [14] in every degree, it keeps full rank at stacked and at fresh generic points, and it has the declared ends and degrees.
- Candidate maps come from small exact solves, composites of elementary maps and singular light leaves.
- Idempotents e_p cut out B_p where needed.

**Pruning.** At fixed N, B_p is dropped if ν(p)_1 > N, where ν(p) is the RSK shape of the column-reversed table. This is the Foam_N kill rule. It has been checked but not proven, and krhom applies it only on 3–4 strands. A summand dropped at N is also dropped at every smaller N, so a library pruned at N = 4 serves N = 2, 3 and 4.

**Closure.** The gl_N closure is the universal construction, as in Robert–Wagner [21]: Hom(R, B) modulo the radical of the pairing (f, g) ↦ τ_N(g∘f). Here τ_N is the top coefficient in ⊗ H*(Gr(λ_i, N)).
- Uncolored: an exact restriction tower gives every N from one scan.
- 3–4 strands: closure entries are sampled at points indexed by the λ_i × (N − λ_i) boxes and decoded by factorial Schur vanishing [22], with a check in every degree.

**Cost.** The cost is governed by the number of indecomposable summands kept after each crossing.
- On 2 strands it grew polynomially in every measured case. For ∧³ on σ^m it was 4, 15, 40, 85, 156, …, that is, Σ_{j≤3} m^j.
- On 3 or more strands it can grow exponentially in the number of crossings. Already at N = 2, reduced Khovanov homology has rank at least det(K). The closures of (σ1⁻¹σ2)^m have determinant L_{2m} − 2, where L is the Lucas sequence; for example, the determinant is 15125 at 20 crossings.
- The scan does not depend on N; the closure does.
- For Khovanov homology (N = 2) of much larger braids, use KnotJob [25].

## Limitations

- **Colored on 3 or more strands.** Only uniform color 2 at N ≤ 4 is supported, through the two libraries. On 4 strands this is experimental and limited to the words the library covers; 8_15, 9_22, 11n_42 and the Conway knot exit with code 4.
- **∧³ on 3 strands** (colors (3,3,3)) is parked. krhom refuses these colors.
  - Its libraries are incomplete at N = 5 (174 of 177 library build jobs) and N = 6 (345 of 377 for F_{2^31−1}, 344 of 377 for F_4294967291).
  - At N = 3 and 4 the values follow from duality.
- **Other theories.** There are no symmetric colorings, no colored HOMFLY-PT homology, and no deformed or equivariant versions.
- **No integral coefficients.** There are no Z coefficients and no torsion. Results are over Q, cross-checked by two large primes. Small characteristic is not supported.
- **Uncolored size.** About 5 strands is the practical limit, because idempotent memory grows quickly with the number of strands; a 6-strand scan reached 77 GB. The code stops at 10 strands and N = 32.
- **Colored on 2 strands.** Colors are limited to 1..3, and color 3 is slow: T(2,9) takes over an hour per prime and about 3 GB. For ∧³ at larger N, raise the closure cap with `--mem-cap-mb`, for example `--mem-cap-mb 4096`.
- **Long runs.** 4-strand words take hours (use `--checkpoint`). No knot of braid index 4 has been computed unreduced.
- **Unproven inputs.** The Foam_N kill rule is validated but not proven, and the random-point checks are probabilistic.
- **Size limits** are checked before any computation: the uncolored tower closure supports ≤ 10 strands and N ≤ 32 (`--closure exact` has no N cap), and the exact 2-strand route needs ∏ C(N, a_i) ≤ 5·10^7 (unreduced color 3: N ≤ 35; color 2: N ≤ 119). Larger inputs exit with code 1 and a message.

## Repository map

```
src/                  library crate `krhom` (module map in src/lib.rs)
  lib.rs, calc.rs     crate root; calculator support: library registry and checks, braid parser, prime comparison, JSON
  calc_ckpt.rs        resumable --checkpoint runs of the 3-4-strand route
  input/              crossing (braid generators, colored Rickard complexes), pd (PD codes -> braids)
  algebra/            field, poly, mpoly, laurent, pmat, sparse, dense, linsys, symfun (Grassmannians)
  combinat/           web, moy (MOY oracle), hecke, schur (q-Schur oracle: Hom formula, decat_chi, kill rule), coset, dcoset
  soergel/            sbim (singular Bott-Samelson bimodules), local (localization), soergel, indec
  uncolored/          engine, engine2 (uncolored route), tower (restriction-tower closure)
  colored/            xengine, cclose (exact 2-strand route); sengine, swire, smap, sexpr, tmpl (localized
                      multi-strand engine); boxclose, boxrun (fixed-N box closure, sharded driver)
  library/            homprov, sindec, smapio, leaves, leafgen, rexplan (Hom provider, certification, library, light leaves)
  reference/          cube (full cube), cabling (reference methods)
  legacy/             algebra, hom, untrunc (obsolete models, not used by the calculator)
  bin/krhom.rs        the calculator (bin/krhom/cli_ext.rs: --knot/--link, --mark, --latex, --estimate, limit checks);
                      the other files in src/bin are developer tools (feature "dev")
tools/                fetch_library.sh, smoke.sh (quick known-value check), randcheck.py, batch_check.py (table checks),
                      wang.py (Wang's formulas)
tools/dev/            developer scripts: sharded box driver, library persistence and analysis (README.md)
data/                 knot and link tables, libraries.tsv (library registry), libstatus/
results/              uncolored/ table checks, colored_gl4/ colored gl_4 results (README.md)
docs/CALCULATOR.md    full calculator reference
docs/dev/             PLAN.md (goals, model, conventions, milestone history), colored_design/ (design,
                      specifications, run reports and reviews; start with SYNTHESIS.md), README.md (index)
tests/, examples/     integration tests (with Wang references and fixtures) and small developer programs
.github/workflows/    CI: build, tests, dev build, smoke test
CHANGELOG.md          release notes
```

Old module paths (`krhom::sbim`, `krhom::xengine`, `krhom::colored::cable`, ...) still resolve through re-exports in `src/lib.rs`.

Developer tools are built with `cargo build --release --features dev` (`cargo build --release` builds only `krhom`):
- **`sbox`** is the box engine behind krhom's 3–4-strand route, with checkpointed, resumable and sharded runs: `sbox box --dir D --shard k/K`, then `sbox finish` (driver script `tools/dev/box.sh`).
- **`xscan`** is the exact 2-strand colored engine, one prime per call: `xscan LABELS WORD NS [--reduced] [--both] [--prime 31|32b]`.
- **`libbuild`** builds, replays, calibrates and checks the coverage of the per-(colors, N, prime) libraries.
- **`engine2`** is the uncolored engine (used by `tools/batch_check.py`); further tools: `twoprime`, `pd2braid`, `g1` (cost calibration), `gp0`, `m4c`, `sprof`, `schur_scan`, `g4_free`, `cube`, `engine`, `colored_cable`, `grwarm`, `b4wstamp`.
- **Development history.** The development logs (the former `runs_*` directories: scripts, queues and logs of the acceptance and result runs), one-off exploration binaries and archived scripts live on branch `archive-pre-restructure`.
- **Other branches.** `data-gl6` holds the libraries (an orphan branch with no code; `tools/fetch_library.sh` reads it).

## Data

The reference tables in `data/` come from KnotInfo and LinkInfo by C. Livingston and A. H. Moore, accessed 2026-11-09:
- `knots_upto12.jsonl`: 2977 prime knots with 3–12 crossings;
- `links_upto11.jsonl`: 4188 oriented records of 1424 prime links with up to 11 crossings.

These two files are **not covered by the code license**. KnotInfo and LinkInfo state no license and ask to be cited:

- C. Livingston, A. H. Moore, *KnotInfo: Table of Knot Invariants*, knotinfo.org (accessed 2026-11-09).
- C. Livingston, A. H. Moore, *LinkInfo: Table of Link Invariants*, knotinfo.org (accessed 2026-11-09).

In those tables, the Khovanov data come from KnotJob (D. Schütz) for knots and from A. Shumakovitch for links. Sources, chirality and conventions are described in `data/README.md`.

## License

The code is licensed under either of
- the MIT license ([LICENSE-MIT](LICENSE-MIT)), or
- the Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE)),

at your option (`MIT OR Apache-2.0`). Unless you explicitly state otherwise, any contribution intentionally submitted for inclusion in this work, as defined in the Apache-2.0 license, is dual licensed as above, without any additional terms or conditions. The files `data/knots_upto12.jsonl` and `data/links_upto11.jsonl` are excluded; see [Data](#data).

## Citation

GitHub shows a "Cite this repository" button generated from [CITATION.cff](CITATION.cff). In BibTeX:

```bibtex
@software{yang_krcalculator_2026,
  author  = {Yang, Hongjian},
  title   = {{KRcalculator} (krhom): a calculator for exterior-colored {Khovanov--Rozansky} $\mathfrak{gl}_N$ link homology},
  year    = {2026},
  version = {0.1.0},
  url     = {https://github.com/Exista-Yang/KRcalculator}
}
```

If you use the tables, please also cite KnotInfo/LinkInfo. Please cite the papers below for the methods.

## References

**Theory and conventions**
1. M. Khovanov, L. Rozansky, *Matrix factorizations and link homology*, arXiv:math/0401268.
2. R. Rouquier, *Categorification of the braid groups*, arXiv:math/0409593.
3. M. Ehrig, D. Tubbenhauer, P. Wedrich, *Functoriality of colored link homologies*, Proc. Lond. Math. Soc. 117 (2018), 996–1040, arXiv:1703.06691.
4. B. Webster, G. Williamson, *A geometric construction of colored HOMFLYPT homology*, arXiv:0905.0486.
5. H. Queffelec, D. E. V. Rose, *The sl_n foam 2-category: a combinatorial formulation of Khovanov–Rozansky homology via categorical skew Howe duality*, Adv. Math. 302 (2016), 1251–1339, arXiv:1405.5920.
6. H. Wu, *A colored sl(N) homology for links in S^3*, Dissertationes Math. 499 (2014), 217 pp., arXiv:0907.0695.
7. Y. Yonezawa, *Quantum (sl_n, ∧V_n) link invariant and matrix factorizations*, Nagoya Math. J. 204 (2011), 69–123, arXiv:0906.0220.
8. H. Murakami, T. Ohtsuki, S. Yamada, *Homfly polynomial via an invariant of colored plane graphs*, Enseign. Math. 44 (1998), 325–360.
9. S. Cautis, J. Kamnitzer, S. Morrison, *Webs and quantum skew Howe duality*, Math. Ann. 360 (2014), 351–390, arXiv:1210.6437.
10. M. Mackaay, M. Stošić, P. Vaz, *A diagrammatic categorification of the q-Schur algebra*, Quantum Topol. 4 (2013), 1–75, arXiv:1008.1348.
11. L. Lewark, A. Lobb, *New quantum obstructions to sliceness*, Proc. Lond. Math. Soc. (3) 112 (2016), 81–114, arXiv:1501.07138.
12. J. Wang, *Colored sl(N) homology and SU(N) representations*, Adv. Math. 474 (2025), 110314, arXiv:2211.08409 (Prop. 1.5: the ∧^a-colored trefoil).
13. J. Wang, *Link homology and loop homology*, arXiv:2604.25174 (Thm. 1.2: ∧^k-colored homology of T(2,2m+1)).

**Methods**

14. G. Williamson, *Singular Soergel bimodules*, Int. Math. Res. Not. IMRN 2011, no. 20, 4555–4632, arXiv:1010.1283.
15. B. Elias, H. Ko, N. Libedinsky, L. Patimo, *Singular light leaves*, J. Eur. Math. Soc. (online first 2025), doi:10.4171/JEMS/1715, arXiv:2401.03053.
16. B. Elias, H. Ko, N. Libedinsky, L. Patimo, *On reduced expressions for core double cosets*, Math. Res. Lett. 32 (2025), 1845–1876, arXiv:2402.08673; *Demazure operators for double cosets*, J. Pure Appl. Algebra 230 (2026), 108207, arXiv:2307.15021; *Subexpressions and the Bruhat order for double cosets*, arXiv:2307.15726; *The atomic Leibniz rule*, Rev. Mat. Iberoam. 42 (2026), 187–228, arXiv:2407.13128.
17. B. Elias, H. Ko, *A singular Coxeter presentation*, Proc. Lond. Math. Soc. (3) 126 (2023), 923–996, arXiv:2105.08563.
18. B. Elias, *Thick Soergel calculus in type A*, Proc. Lond. Math. Soc. 112 (2016), 924–978, arXiv:1009.2120.
19. B. Elias, G. Williamson, *Soergel calculus*, Represent. Theory 20 (2016), 295–374, arXiv:1309.0865.
20. M. Khovanov, L.-H. Robert, E. Wagner, *A foamy approach to Soergel bimodules*, arXiv:2603.22957.
21. L.-H. Robert, E. Wagner, *A closed formula for the evaluation of foams*, Quantum Topol. 11 (2020), 411–487, arXiv:1702.04140.
22. I. G. Macdonald, *Schur functions: theme and variations*, Sém. Lothar. Combin. 28 (1992), Art. B28a (factorial Schur functions).
23. D. Bar-Natan, *Fast Khovanov homology computations*, J. Knot Theory Ramifications 16 (2007), 243–255, arXiv:math/0606318.
24. P. Vogel, *Representation of links by braids: a new algorithm*, Comment. Math. Helv. 65 (1990), 104–113.

**Data and software**

25. D. Schütz, *KnotJob* (software), https://www.maths.dur.ac.uk/users/dirk.schuetz/knotjob.html (accessed 2026-11-14).
26. C. Livingston, A. H. Moore, *KnotInfo* and *LinkInfo* (see [Data](#data)).
