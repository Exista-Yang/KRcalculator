# CALC-3 report: `--convention` removed, CALCULATOR.md refreshed, C2-7 and dev merged

Branch `colored-calc3`, created from origin/colored-calc2 566de7c; merged origin/colored-2c4-unblock
3394e28 (C2-7 results) and origin/dev ed1e117. Machine: 4 cores, `KRHOM_THREADS=4`, heavy runs under
`ulimit -v 14000000`. No certification check changed, src/leaves.rs untouched, no README.md.

## 1. `--convention` removed

The user decided not to ship it (dev ed1e117). Removed from src/bin/krhom.rs: the usage synopsis and
help line, `--convention` in `VALUE_OPTS` (now 18 value options), the `Opts::knotjob` field and its
parsing, `to_knotjob`, the stdout line `convention: knotjob (q -> q^-1)` and the JSON field
`convention`. src/calc.rs had no helper for it. `--convention` is now an unknown argument (exit 1,
`krhom: unknown argument "--convention"` + usage). `--table` is unchanged.

tests/calc.rs: `table_and_knotjob_convention` replaced by `table_and_no_convention_option`
(krhom-graded table of `aaa` at N = 2, no `convention` in stdout or JSON, `--table` in `--file`
mode incl. the zero table, `--convention knotjob` and `--convention krhom` refused with exit 1).

docs/CALCULATOR.md: synopsis, the flag-table row and the JSON field list no longer mention it; the
output section keeps one sentence: krhom's q-grading is the inverse of Bar-Natan / KnotJob /
KnotInfo (q → q⁻¹, t unchanged; `aaa` at N = 2 is `t^0q^-3 + t^0q^-1 + t^2q^-5 + t^3q^-9` here and
`q + q^3 + t^2q^5 + t^3q^9` there).

## 2. CALCULATOR.md refresh

(a) Support matrix, uncolored row: "any N ≥ 1 / any" → "1 ≤ N ≤ 32 (tower closure) / ≤ 10 strands
(tower closure); about 5 in practice" (src/tower.rs: `MAX_VARS = 10`, truncation cap assert
"N ≤ 32"; observed: `krhom 33 --braid aaa` panics in the engine with that assert, reported as
"engine panicked ... please report this" — not a clean refusal, see §5).

(b) The 4-strand projection paragraph ("0.2–1.3 core-hours reduced, 1.2–8 unreduced, never run to
completion") is replaced by the measured C2-7 reduced runs. Sources: docs/dev/colored_design/c2/C2_7_report.md
§6, results/colored_gl4/gl4_2222_c27_compare.json (status PASS, ranks 81/361/441, F31 = F32b) and the per-shard and
finish logs runs_c2/c27/k4/out/{knot}_r_{prime}_s{0..3}/box.out and {knot}_r_{prime}.err (s per
point, shard seconds, VmHWM, finish total). 4 shards × 1 thread × 54 points + finish (4 threads):

| knot | prime | s/point | wall | shard core-hours | VmHWM/shard |
|---|---|---|---|---|---|
| 6_1 | F31 / F32b | 26.3–26.5 / 26.8–27.1 | 0.75 h / 0.75 h | 2.21 / 2.15 | 2681–2682 / 2758–2768 MB |
| 7_6 | F31 / F32b | 36.5–36.8 / 39.5–39.9 | 0.97 h / 1.01 h | 2.95 / 3.14 | 2811 / 2898 MB |
| 7_7 | F31 / F32b | 35.9–36.3 / 25.7–25.9 | 0.93 h / 0.70 h | 2.90 / 2.04 | 2764 / 2850–2854 MB |

(wall = slowest shard + finish; core-hours exclude the 4-thread finish.) Unreduced: 1 296 points
instead of 216, at least 1 296 × 26.3 s ≈ 9.5 core-hours of scanning per prime for 6_1 (lower bound;
the C2-7 report's estimate is ≈ 4 h wall per (knot, prime) with 4 shards).

## 3. Merges

- origin/colored-2c4-unblock 3394e28: one conflict, src/bin/gp0.rs (`--steps` census max step).
  CALC-2 item 13 counted a replaced *final* step with its replacement's score sw; C2-7 set *detoured*
  lower-level steps to 0. Resolution keeps both intents: every replaced step (C2-7 detour at a lower
  level or CALC-2 per-predecessor final word) is counted with its replacement word's own score
  (map `repl` keyed by (level, to, from)); "steps > 1e5" uses the same effective score. gp0 is a
  diagnostic binary; no certification check is involved.
- origin/dev ed1e117: clean (LICENSE-MIT, LICENSE-APACHE, CITATION.cff, Cargo metadata, docs).

## 4. Acceptance

- `KRHOM_THREADS=4 cargo test --release --no-fail-fast` (libraries installed): **311 passed, 0 failed,
  41 ignored** (61 test binaries; runs_calc3/full_suite.log). A first run had 310/1/41: my new test
  asserted the wrong term order; fixed in d5a7505, `--test calc` 13 passed, then the full re-run above.
- `krhom 2 --braid aaa --convention knotjob` → exit 1, `krhom: unknown argument "--convention"`.
- `krhom 2,3,4 --braid aaa --colors 2 --table`: identical to colored-calc2 apart from the timing
  lines (diff with `[F_...]`/`time` lines removed is empty; output runs_calc3/table_aaa_c2.txt).
- `tools/fetch_library.sh l222_n4 l2222_n4`: installed (l2222_n4/f32b: 2873 layer_a, 3029 tmpl).
- `tools/randcheck.py --seed 20261115 --quick --out runs_calc3/random` (runs_calc3/randcheck_quick.log):
  self-tests passed in 162 s; exit 0;
  `seed=20261115 n=10 PASS 10 FAIL 0 WARN 0 SKIP 0 TIMEOUT 0 ERROR 0; checks: primes 21, seeds 13, chi 42, isotopy 6, markov 5, mirror 8, duality-N3 18, reversal 0, N2 21, lewark 10, unknot 4, closed-form 45, cross-engine 8; 454s`
  (same check counts as CALC-2).

## 5. Not done / observations

- N > 32 on the uncolored route is not refused cleanly: `krhom 33 --braid aaa` hits the tower
  assert and is reported as an engine panic. The docs now state the limit; a clean exit-1 refusal
  before the engine would be a small follow-up (not in this session's scope).
- Unreduced 4-strand runs were not attempted here (≥ 9.5 core-hours per prime).
