# CALC-1 report: the krhom calculator

Session CALC-1, 2026-11-14, branch `colored-calc` (from origin/colored-2c4-unblock a2f14fe, with
origin/dev 4a48c5f merged; origin/colored-2c4-unblock merged again at aba6192). Spec:
CALC_DECISIONS.md (binding) > CALC_SPEC_review.md required corrections > CALC_SPEC_v1.md.
Reference documentation for users: docs/CALCULATOR.md.

## 1. What was built

| Piece | Where | Priority |
|---|---|---|
| Library registry (single source of truth, 4 rows, data ee548d2, per-layer sha256, stamp layers) | data/libraries.tsv | P0 |
| Status files (C2-7 union job status, diagnostics only) | data/libstatus/l2222_n4_{f31,f32b}.txt | P0 |
| Fetch / verify tool | tools/fetch_library.sh | P0 |
| Calculator library (registry, `locate`, stamp preflight, `run_prime`, `run_prime_live`, classifier, `diagnose_missing`, `compare_primes`, braid parser, JSON) | src/calc.rs (+ 8 unit tests) | P0 |
| One-command CLI with automatic routing | src/bin/krhom.rs | P0 |
| CLI tests (default suite) + `#[ignore]` library tests | tests/calc.rs (8 + 4 ignored) | P0 / P1 |
| Updated CLI test | tests/xengine.rs `cli_singular_method_matches_cube_output` | P0 |
| sbox uses `calc::poly_json` (pure move) | src/bin/sbox.rs | P0 |
| Random-example checker | tools/randcheck.py | P1 |
| `--check-library`, batch `--file`, `--csv` | src/bin/krhom.rs | P1 |
| User reference | docs/CALCULATOR.md | P1 |
| Acceptance runs | runs_calc/, runs_random/ | — |

Not done (P2): checker Mode C (table knots), sbox reusing the calc helpers beyond `poly_json`, a
tiny in-repo library fixture, resumable `--ckpt`.

### The command

```
krhom N[,N...] (--braid WORD [--strands n] | --pd PD | --file BRAIDS) [--colors c | c1,...] [--reduced] [--json OUT] [--csv OUT] ...
```

- WORD: signed integers or KnotJob letters (a..i, A..I); mixtures refused. `--strands` defaults to
  max |generator| + 1 and must be at least that if given.
- Routes (D2): all colors 1 → Soergel engine2, two primes (seeds 7 / 11); colored ≤ 2 strands
  after simplification → exact xengine path, two primes; colored 3–4 strands, uniform color 2,
  N ≤ 4 → box route with the installed library, two primes; anything else colored → exit 1 with
  the support matrix and the `--colored-method cube` hint. `--colored-method cable|cube|singular`
  stays as an override; `cube` is now exact over Q (`khr_full_cube::<Q>`), `cable` now runs both
  primes and compares.
- stdout: `coefficients: Q` replaces the `fields:` line (D1); the per-prime certification lines go
  to stderr and `--json`.
- Exit codes 0 / 1 / 2 / 3 / 4 / 5 as in spec B; engine panics are caught and reported as exit 2.

## 2. Acceptance (spec H as amended)

All runs: release build, `KRHOM_THREADS=4`, 4-core cloud machine, under `ulimit -v 14000000`.
Outputs, stderr and JSON are in runs_calc/.

### H1. Test suite

`KRHOM_THREADS=4 cargo test --release` (runs_calc/full_suite.log): **306 passed, 0 failed, 41
ignored** (the m5i base was 287 / 37; new: 8 calc unit tests, 8 tests/calc.rs, 4 ignored library
tests, plus the C2-7 tests of the base). The 4 ignored library tests of tests/calc.rs (H3-H6 via
cargo) were also run with the libraries installed: see section 6.

### H2. Fetch tool

- `tools/fetch_library.sh l222_n4 l2222_n4`: 4 directories installed in 98 s (l222_n4: 233
  layer-A + 471 template files per prime; l2222_n4: 2873 layer-A + 3048 (F31) / 3029 (F32b)
  template files). A second run: `already installed` for all four (0.7 s, full sha256 check).
- On a scratch copy: deleting one layer-A file → `--verify-only` exits 2 (`1 listed file(s)
  missing or changed`, `run tools/fetch_library.sh --force l222_n4/f31`), and krhom exits 4:
  `damaged install: 1 of the 704 files listed in .../LIBRARY.txt are missing (e.g.
  layer_a/a-0072264e8c90a791.txt): run tools/fetch_library.sh --force l222_n4`. `--force`
  reinstalled it and `--verify-only` passed again.
- Stamp preflight (K4) on a scratch copy: one layer-A header edited → krhom exits 4 with `1 of the
  233 layer-A files of the stamped layers of l222_n4/f31 do not match this binary's route stamp
  ...`.

### H3. 4_1 on 3 strands

`krhom 2,3,4 --braid "1 -2 1 -2" --strands 3 --colors 2` (runs_calc/h3_u.out): exit 0, ranks 1, 11,
46; N=3 = `t^-2q^6 + t^-2q^8 + t^-1q^0 + t^-1q^2 + t^0q^-2 + t^0q^0 + t^0q^2 + t^1q^-2 + t^1q^0 +
t^2q^-8 + t^2q^-6`; N=4 byte-identical to the unreduced N=4 string of results/colored_gl4/gl4_222_4_1_f31.json.
`--reduced` (h3_r.out): ranks 1, 5, 25; N=3 = `t^-2q^6 + t^-1q^2 + t^0q^0 + t^1q^-2 + t^2q^-6`;
N=4 identical to the JSON's reduced string. N=2 = `t^0q^0` in both modes.

### H4. Stabilized trefoil

`"1 1 1 2" --strands 3 --no-simplify`: N blocks identical to `"1 1 1" --strands 2` in both modes
(ranks 1, 7, 22; reduced 1, 3, 9; N=3 strings as in the spec).

### H5. 4 strands

- `"1 2 3" --strands 4 --no-simplify`: ranks 1, 3, 6; reduced 1, 1, 1 (= [N 2]_q and t^0q^0).
- `"1 1 1 2 3" --strands 4 --no-simplify`: N blocks identical to the 2-strand trefoil, both modes.

### H6. `--check-library` at 4 strands

- 6_1 `"1 1 2 -1 -3 2 -3"`: exit 0, `check-library: library l2222_n4 covers this word (2 primes;
  no box scan run)`. 1119 layer-A loads per prime, 0 rejected, 94 types, closure data for every
  closable type (so the record scan was not needed, review K10); 221 s + 216 s; VmHWM 958 /
  1046 MB. This was the first 4-strand library build through e_p (2 idempotent types, M5i).
- 8_15 `"1 1 -2 1 3 2 2 2 3"`: exit 4 in 101 s, `library l2222_n4 does not cover this word: 11 jobs
  not certified for both primes: pair 926 (p=..., rung=E1.1), 949, 986, 997, 1054, 1412, 1435,
  1471, 1472, 1483, 1540` — exactly the cover list.

### H7. Refusal and exit-code matrix

From tests/calc.rs (all pass) and the manual runs:

| Input | Exit | Message (stderr, abridged) |
|---|---|---|
| 4_1 colors 2, empty library dir | 4 | `library l222_n4 (F_2147483647) is not installed under <root>: run tools/fetch_library.sh l222_n4 (about 3.5 MB per prime)` |
| 6_1 colors 2 `--check-library`, empty library dir | 4 | `... run tools/fetch_library.sh l2222_n4 ...` |
| one listed library file deleted | 4 | `damaged install: 1 of the 704 files listed ... --force l222_n4` |
| a stamped layer-A header changed | 4 | `1 of the 233 layer-A files of the stamped layers ... route stamp` |
| 8_15 colors 2 | 4 | `does not cover this word: 11 jobs ...` |
| `"1 1 2 2" --strands 3 --colors 2,2,1` | 1 | `labels [2, 2, 1] on 3 strands are not uniform; supported on 3-4 strands: ...` |
| 4_1 `--colors 3` | 1 | `no library for labels [3, 3, 3] on 3 strands; ...` |
| `"1 1 2 2 3 3 4 4" --strands 5 --colors 2` | 1 | `no library for labels [2, 2, 2, 2, 2] on 5 strands; ...` |
| `krhom 5` 4_1 colors 2 | 1 | `N = 5 is above the library's N <= 4 ...` |
| 4_1 uncolored `--colored-method singular` | 1 | `uncolored input: drop --colored-method singular; ...` |
| `"1 1" --strands 2 --colors 2,2,2` | 1 | `3 colors given, but the link has 2 components ...` |
| `"1 1" --strands 3 --colors 2 --no-simplify` | 1 | `split link (sigma_2 never occurs): compute the pieces separately` |
| `2,x` / `0` as N | 1 | `N list "2,x": ...` / `every N must be at least 1` |
| `"1 0 1"`, `a1`, `"a b"`, `"1 1.5"` | 1 | generator 0 / mixes letters and integers / no spaces / not a signed integer |
| `"1 2" --strands 2` | 1 | `--strands 2 is too small for the braid ...` |
| `--colors 0`, `--colors 2,x` | 1 | `colors are integers >= 1` |
| bad PD (3 labels; a label once) | 1 | `invalid PD code: ...` |
| unknown option; `--colored-method magic`; both `--braid` and `--pd` | 1 | usage / `must be cable, cube or singular` / `exactly one of` |
| colors 2 `--colored-method cube --reduced` | 1 | `--reduced is only supported ...` |
| `"1 1 1" --colors 4` (2 strands) | 1 | `the exact 2-strand engine handles colors 1..=3; --colored-method cube ...` |
| `KRHOM_NO_M5I=1`, box route | 1 | `validation-only environment switches ...: unset KRHOM_NO_M5I` |
| engine / certification failure | 2 | `localized engine failed on F_...: ... certification failure ...` (not triggered) |
| primes disagree | 3 | `the two primes disagree (torsion or a bug); no result printed` (not triggered) |
| resource refusal | 5 | per-point memory / memory cap messages (not triggered) |

### H8. Library untouched

After H3–H6 (including the 8_15 failure path) and H10: `find <dir> -newer <dir>/LIBRARY.txt -type
f` is empty for all four directories, the counts are unchanged (233/471, 233/471, 2873/3048,
2873/3029), and `tools/fetch_library.sh --verify-only l222_n4 l2222_n4` passes. Answer to Q4: a
replay run writes nothing into the library (layer A is never written by a replay; no template was
recomputed in any of these runs, so TmplLib stored nothing).

### H9. Not applicable (base is colored-2c4-unblock, which has the route_header filter).

### H10. Checker quick run

`tools/randcheck.py --seed 20261115 --quick` (runs_calc/h10_quick.log):

```
seed=20261115 n=10 PASS 10 FAIL 0 WARN 0 SKIP 0 TIMEOUT 0 ERROR 0; checks: primes 21, seeds 13, chi 42, isotopy 6, markov 5, mirror 8, duality-N3 18, reversal 0, N2 21, lewark 10, unknot 4, closed-form 35, cross-engine 8; 484s
```

`--inject-fault --quick --only T3-01` (runs_calc/h10_fault.log): T3-01 FAIL (O5 and O12: `N=4 ...
t^0q^-4: 2 vs 1`), exit 1.

### H11. Checker default run

`tools/randcheck.py --seed 20261115` (runs_calc/h11_default.log; runs_random/20261115/), `--jobs 2`
(2 threads per krhom), 31 min:

```
seed=20261115 n=24 PASS 24 FAIL 0 WARN 0 SKIP 0 TIMEOUT 0 ERROR 0; checks: primes 47, seeds 31, chi 94, isotopy 10, markov 15, mirror 17, duality-N3 41, reversal 4, N2 47, lewark 22, unknot 6, closed-form 72, cross-engine 16; 1871s
```

No WARN. T2 8, T3 12 (6 Mode A including 2 unknot cores, 6 Mode B: partners reverse x4,
mirror, rotate), T4 4 (all 120 T4 candidate words were covered by the library; the 4
words `1 3 2`, `-1 2 -1 -1 3`, `-1 3 -1 -1 -2`, `1 1 -3 1 -2` took 160-510 s per example). Every
box run used different node seeds for the two primes (201 and 202 + subseed mod 97) and agreed.

## 3. Timings (per prime, this branch, 4 threads)

| Run | F_2147483647 | F_4294967291 | library load (per prime) | box points | VmHWM |
|---|---|---|---|---|---|
| H3 4_1 unreduced | 43.5 s | 38.2 s | 6.4 / 5.6 s | 216 | 163 / 166 MB |
| H3 4_1 reduced | 16.5 s | 15.8 s | 5.5 / 4.9 s | 36 | 160 MB |
| H4 `1 1 1 2` unreduced | 0.9 s | 0.9 s | 0.06 s | 216 | 16 MB |
| H4 reduced | 0.2 s | 0.2 s | 0.06 s | 36 | 15 MB |
| H5 `1 2 3` unreduced | 83.0 s | 126.0 s | 26.2 / 28.1 s | 1277 | 106 / 119 MB |
| H5 `1 2 3` reduced | 46.0 s | 55.6 s | 27.1 / 30.0 s | 216 | 104 / 116 MB |
| H5 `1 1 1 2 3` unreduced | 95.9 s | 91.3 s | 27.2 / 25.4 s | 1296 | 116 / 128 MB |
| H5 `1 1 1 2 3` reduced | 41.3 s | 54.7 s | 25.9 / 27.7 s | 216 | 109 / 121 MB |
| H6 6_1 `--check-library` | 222 s | 217 s | 221 / 216 s | — | 958 / 1046 MB |
| H6 8_15 `--check-library` (exit 4) | 101 s total, including the diagnosis | | | | |

VmHWM is process-wide (the F32b number includes the F31 run in the same process).

## 4. Answers to the open questions

- **Q2** (old `c2/l2222_n4` base layer): restored first, as the MANIFEST order says. It contributes
  no final file: every one of its layer-A and template files is overwritten by a later layer.
  Final layer-A files by last writer: F31 2612 (c25) + 249 (c24) + 12 (c27) = 2873; F32b 2612 + 247
  + 14 = 2873. Templates: F31 2669 + 209 + 170 = 3048; F32b 2649 + 218 + 162 = 3029. The stamp
  preflight requires every layer-A file of layers c25, c24, c27 to match the route stamp: 2873 of
  2873 match, 0 legacy files.
- **Q3** (`set_replay_only(true)`): on by default. No run in H3–H6, H10 or H11 needed a template
  recompute. If a run fails with "not in layer A" and the diagnosis finds no missing job, krhom
  retries once with replay-only off and records `"template_recompute": true` in the JSON. On the
  8_15 path the diagnosis now runs first, so a genuinely uncovered word never triggers the retry
  (the first H6 8_15 run, before that change, did retry: 166 s; after: 101 s).
- **Q4**: see H8. Nothing is written into the library by replays.
- **Q13** (single-mode box grids): validated. Unreduced-only grids (216 / 1277 / 1296 points) and
  reduced-only grids (36 / 216 points) reproduce the production u,r values: H3 N=4 strings equal
  results/colored_gl4/gl4_222_4_1_f31.json in both modes; H4 and H5 equal the 2-strand exact path in both modes.
  Note (review K1): an unreduced-only run also closes the reduced closures from the unreduced
  extraction (χ is checked for them: `chi = decat_chi for 6 closures`); krhom does not print or
  export them.
- **Q5 / K14**: N below the color prints the zero polynomial on every route, with a stderr note,
  without running an engine (on the box route also without touching the library).
- **Q8**: deferred (documented cost; `--reduced` recommended; sbox for resumable runs).
- **Q9**: `ulimit -v` documented; no RLIMIT from krhom.
- **Q10**: seed 201 for both primes by default; `--seed-f32b` is a hidden developer flag used by
  the checker.
- **Q11**: split braids on 3–4 strands refused.

## 5. Deviations from the spec, and conservative choices

1. **`--mem-cap-mb` precedence**: an explicit `--mem-cap-mb` wins over `KRHOM_MEM_CAP_MB`; the
   environment variable wins over the default. (The spec let the environment variable win over
   the flag.) Default (K7): 4500 × max(1, threads/2) at 3 strands, 200000 at 4 strands.
2. **Library file list**: LIBRARY.txt keeps the spec's header keys plus `stamp_layers` and a
   `files:` section of `LAYER SHA256 PATH` lines (D4/K3). krhom checks presence of every listed
   file on each run (not hashes; hashes are `--verify-only`). Extra files are tolerated.
3. **Stamp preflight (K4)**: per-layer, via the registry's `stamp_layers` column (l222_n4: layer 0;
   l2222_n4: layers c25, c24, c27). Any mismatching file of those layers is exit 4.
4. **Template-recompute retry (Q3)**: automatic, once, only after the diagnosis finds no missing
   job (or for a `complete` library), recorded in the JSON.
5. **Environment hygiene (K12)**: refused on the box route: `KRHOM_FORCE_LADDER`, `KRHOM_NO_M5I`,
   `KRHOM_DCID_REF`, `KRHOM_M5I_EXACT_UNK`, `KRHOM_M5I_NO_PREFIX`, `KRHOM_M5I_NO_REROUTE`,
   `KRHOM_NO_REROUTE_FALLBACK`, `KRHOM_LEAVES_C27` (non-empty). `KRHOM_LOC_FULL` (the reference
   localization, identical values) is allowed. All `KRHOM_*` variables are recorded in `--json`.
6. **Simplification** applies on the singular route (default for colored input), also on 2
   strands; the 2-strand results are unchanged (same link). Same-strand rewrites are not used.
7. **Cable**: runs both primes (seeds 7 and 11) and compares, instead of becoming a hidden method.
8. **Cube** is always over Q now (`--field q` no longer needed; `--field` is accepted and ignored
   with a note).
9. **Unknown options are refused** (exit 1) instead of being ignored, so a typo cannot change the
   route. A repeated value option is refused too.
10. **Panics** in engines are caught and reported as exit 2 (the batch continues).
11. **`--box-live`** requires a nonempty braid on exactly 2 strands; it runs both primes (seeds
    `--seed` / `--seed-f32b`).
12. **Batch output**: each block starts with `name: NAME`; a failed line prints `error: exit K
    (kind)` on stdout and the message on stderr; the exit code is the first failing line's code.
13. **JSON**: one object per input (array for `--file`) with `input.kind/text`, `braid`,
    `used_braid`, `route`, `coefficients`, `library`, `primes[]` (with stamp preflight counts and
    `template_recompute`), `agree`, `zero_N`, `results[]`, `env`, `error`. The error object has no
    `route`.
14. **randcheck**: T3 partner length capped at 6 and per-call timeouts T3/T4 2400 s (K6); T4
    timeouts 2400 s rather than 1500 s; T2 runs N' = 2..6 on the exact path; Mode-B partners are run
    unreduced only; reduced comparisons only for knots; O11 labels B for m = 4 and for m = 3 at N ≥ 5
    and uses the |n| = 1 and n = 0 closed forms (no Hopf |n| = 2 form); `--no-anchor` skips the
    4_1 anchor (also skipped with `--only`); the quick and default outputs of the same seed would
    share `runs_random/<seed>/`, so the quick run is kept under `runs_random/20261115-quick/`.
15. The data commit is fetched by SHA with a fallback to the `data-gl6` branch (K11).

## 6. Ignored library tests (F3)

`KRHOM_THREADS=4 cargo test --release --test calc -- --ignored` with both libraries installed
(runs_calc/f3_ignored.log): 4 passed (lib_41_three_strands, lib_stabilized_trefoil,
lib_four_strands_markov, lib_check_library_4_strands), 559 s with the four tests in parallel.
After them, H8 still holds: no file newer than LIBRARY.txt in any library directory, and
`--verify-only` passes for all four.

## 7. Example commands in the new syntax (runs_calc/examples_new_syntax.out)

- `krhom 2,3,4 --braid aaa --colors 2`: letters, default strands 2, exact path, ranks 1, 7, 22.
- `krhom 2,3 --braid aBaB`: uncolored 4_1, default strands 3, Soergel engine on both primes,
  ranks 6 (Khovanov) and 11.
- `krhom 2,3,4 --braid abc --colors 2 --reduced`: default strands 4, simplified to the 1-strand
  unknot (`simplified braid: (empty)  (strands 1)`), t^0q^0 for every N.
No `--colored-method` is needed for any of them.
