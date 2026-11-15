# CALC-2 report: fixes after the review of colored-calc

Branch `colored-calc2`, created from origin/colored-calc 1a9c8de with origin/dev 1914bad merged
(CALC_2_FIXES.md, the review files review_calc1/*.md and CALC_DECISIONS.md, all binding). Machine:
4 cores, 15 GB, `KRHOM_THREADS=4`, heavy runs under `ulimit -v 14000000`.

Hard limits respected: no certification check changed, src/leaves.rs untouched, no table sweep,
no library build, no new dependency, no README.md. Libraries installed with
`tools/fetch_library.sh l222_n4 l2222_n4` (the CALC-2 script; 4 directories, 60 s).

## Summary per item

| # | Item | Status |
|---|---|---|
| 1 | exact 2-strand route: cap refusal exit 5, `--mem-cap-mb` on every route, shape exit 1, docs | done |
| 2 | N below the color: zero before any refusal and without an engine | done |
| 3 | batch: `name:` before any early return; name on the error line | done |
| 4 | CSV: RFC 4180 quoting of the name | done |
| 5 | batch: JSON/CSV written per finished item | done |
| 6 | JSON `input` without a duplicate key | done (already so in code; test added) |
| 7 | randcheck: krhom exit 4 is FAIL | done |
| 8 | randcheck: harness failures exit 2 | done |
| 9 | fetch_library.sh: `|| die` on every install step, empty staging refused, tar check without SIGPIPE | done |
| 10 | fetch_library.sh / LIBRARY.txt: costs and status sources recorded and re-checked | done |
| 11 | CALCULATOR.md: real randcheck timeouts | done |
| 12 | randcheck O11: Hopf \|n\| = 2 closed form | done |
| 13 | gp0 census / c27_joint `plan_max`: replacement score counted | done |
| 14 | sindec reroute_family fail-closed for a decomposable `a` | **tried, reverted** (changes 4 certified jobs); comments fixed |
| 15 | `--table` | done |
| 16 | `--convention knotjob` (verified map q → q⁻¹) | done |

## Must fix

### 1. Exact 2-strand route: resource and shape refusals

- `xengine::khr_singular_exact` now closes with `XScan::try_close` → `try_close_guarded` →
  `cclose::try_close_exact_colored` (Err instead of a panic for the memory cap and A_λ shape
  refusals; the pruned-scan d∘d guard, a certification check, is unchanged and still panics).
  `XScan::close` and `close_guarded` are unchanged for the other callers.
- krhom maps the errors (`exact2_err`): "above the memory cap" → exit 5 with `with krhom use
  --mem-cap-mb M to raise it`; "top weight" / "labels above 16" → exit 1; anything else exit 2.
- Before any engine work krhom checks the shape (`exact2_shape`, the rule of
  `cclose::check_alg_shape`): every unmarked block of color a needs a(N − a) ≤ 255, i.e. N ≤ 256
  (color 1), 129 (color 2), 88 (color 3). Exit 1 with that text.
- `--mem-cap-mb` sets `KRHOM_MEM_CAP_MB` in `main` before the first item (every route); the box
  route keeps its own default when the flag is absent. The exact route records the effective cap
  in JSON `mem_cap_mb` (env or 768).
- docs/CALCULATOR.md §2 (formerly "any N") now gives the real limits; §5 and §6 updated.

### 2. N below the color first

In `run_item`, right after the zero-N computation: if every requested N is below the largest
color (and not `--check-library`), krhom prints `method: none (every N is below the color a: the
homology is 0, nothing to compute)` and the zero blocks, before any route refusal (no library,
colors > 3, 5 strands, mixed labels) and before any engine (the exact route no longer runs xscan
for an empty N list). JSON route `"zero"`. An N at or above the color still goes through the
refusals (`krhom 2,3,4 --braid "1 1 1" --colors 4` → exit 1). The CALC-1 test line `krhom 2 ...
--colors 4 → exit 1` contradicted this item and was changed to N = 4.

### 3. Batch names

`name: NAME` is the first line of `run_item` (before parsing). The stdout error line is
`error: NAME: exit C (kind)` (stderr: `krhom: NAME: message`, as before).

### 4. CSV quoting

`csv_field`: a name containing `,` `"` CR or LF is quoted with doubled quotes (RFC 4180):
`"L6a5{0,1}"`, `"""q"""`.

### 5. Incremental batch outputs

- CSV: the file is created (and the header written) at start, so a bad path fails before any
  work; each item's rows are appended and flushed when it finishes.
- JSON: in batch mode `[\n]` is written at start and the whole array of finished records is
  rewritten after every item atomically (temporary file in the same directory + rename; direct
  write as fallback), so the file is always a valid JSON array. Single input: the path is checked
  for writability at start, the object is written when the item finishes.
- stdout is flushed after each item.

### 6. JSON `input`

At 1a9c8de the code already wrote `"input":{"kind":"braid","text":"...","N":...,"braid":[...],...}`;
the duplicate `"braid"` key exists only in runs_calc/*.json, written by the earlier build 95ec4a9.
No rename was needed; tests/calc.rs now asserts that every `input` object has `"text"` once and
`"braid"` at most once (the array).

### 7. randcheck exit 4 = FAIL

`run_krhom` maps krhom exit 4 to FAIL (was SKIP(library)). T3 runs only with the complete
l222_n4 and T4 only on words that passed the m5i_cover filter, so every exit 4 there is a broken
install, a stamp/registry mismatch or a map rejected on load. (T4 examples without a covered
candidate are still SKIP before any krhom call.) Checked with a stub krhom that exits 4: status
FAIL.

### 8. randcheck harness failures exit 2

`HarnessError`: tools/wang.py failures (also with a 600 s timeout), the `cargo build --example
m5i_cover` and the m5i_cover run raise it; `run_example` re-raises it (it is not an example
ERROR); the convention self-test call is wrapped (an error becomes a self-test failure → exit
2); `__main__` maps HarnessError and any other uncaught harness exception to exit 2. Checked:
`randcheck.py --seed 1 --quick --tiers 2 --krhom /bin/false` → `self-test FAILED: convention
call failed: exit 1 ...`, exit 2.

### 9. fetch_library.sh guards

- Every install step has `|| die`: mktemp (the result is also checked non-empty and a
  directory), every mkdir, cp, the owners/final/LIBRARY.txt writes, git init of the temporary
  repository, and the moves. If the new copy cannot be moved into place, the old copy is moved
  back before dying; the old copy is removed only after the new one is in place.
- The tar checks read here-strings/variables (`grep -qv ... <<<"$types"`, `bad_entry="$(grep -Ev
  ... <<<"$names" || true)"`), so they cannot fail open under pipefail.
- Checked: `--root FILE` (not a directory) → `fetch_library: cannot create the install root ...`,
  exit 2.

### 10. Sources of costs.txt and status.txt

LIBRARY.txt gains `costs_source: <path in the data commit>` and `status_source: complete` or
`status_source: data/libstatus/<file> <sha256>`. `verify_installed` (used by "already
installed" and `--verify-only`) requires `costs_source` = the registry's costs path (same
data_commit + same path = same bytes; the installed costs.txt is still checked against its
recorded hash) and, for an incomplete library, the installed status.txt and the recorded hash to
equal the current checkout file. Checked on a scratch copy: an appended status line → `status.txt
missing or differs from data/libstatus/l2222_n4_f31.txt in this checkout`, exit 2; a changed
`costs_source` → `costs_source differs from the registry`, exit 2. Installs made by the CALC-1
script lack these lines: they fail verification once and are reinstalled by a plain
`tools/fetch_library.sh ID` ("does not verify; reinstalling"); `--verify-only` tells the user to
run `--force`. krhom itself ignores the new lines.

## Should fix

### 11. Timeouts

docs/CALCULATOR.md §7 now lists the per-call timeouts of the code: 120 s (2-strand exact-path
calls and the convention self-test), 300 s (`--box-live`), 600 s (uncolored O8), 2400 s (box
calls and the 4_1 anchor); and the exit-code rules of items 7 and 8.

### 12. O11 Hopf

`hopf_closed(n, a, b)` is tests/common/mod.rs `hopf_formula` (Wang Thm 3.5) in Python (q-multinomial
as a product of symmetric Gaussian binomials). O11 uses it for unreduced 2-strand words with
writhe ±2 (mirrored for −2), label C. A new self-test compares it with the exact path for the
(2,2) Hopf link at N = 2..5 (the uncolored value at N = 2, `1 + q^-2 + t^2q^-4 + t^2q^-6`, is
Khovanov's Hopf link in krhom's convention).

### 13. gp0 census

src/bin/gp0.rs: a final step replaced by a per-predecessor word counts with the replacement's
own W-aware score `sw` in `mx` and in `over` (instead of 0 and an unconditional decrement).
tests/c27_joint.rs `plan_max` does the same (computes `sw` with `step_score_w` as gp0 does), so
`per_predecessor_final_pair_1444` now really asserts the replacement step ≤ 1e5. src/leaves.rs is
unchanged. `cargo test --release --test c27_joint`: 3 passed (pair 1444's replacement step is
≤ 1e5, so the published C2-7 figure holds).

### 14. sindec reroute_family — tried and reverted

Implemented as specified (return None when the prefix type `a` is decomposable,
`!is_indecomposable_of(w[..j], a)`), then replayed the decomposable-source reroute jobs of
data/libstatus on scratch copies of the installed l2222_n4 for both primes
(`libbuild one --lam 2,2,2,2 --n 4 --field F --cache COPY --costs COPY/costs.txt --provider
m4c-lf --stamp-route --replay -- pair J`, script runs_calc2/replay14.sh):

| job | fail-closed change, F31 / F32b | reverted (this branch), F31 / F32b |
|---|---|---|
| pair 109 | ok / ok | ok / ok |
| pair 1789 | ok / ok | ok / ok |
| pair 1812 | **Refused** / **Refused** | ok / ok |
| pair 2052 | **Refused** / **Refused** | ok / ok |
| pair 2086 | **Refused** / **Refused** | ok / ok |
| pair 2275 | ok / ok | ok / ok |
| pair 2298 | **Refused** / **Refused** | ok / ok |

(The Refused text: `target ... at v^0 of ... ⊗ E2.1: replay: pairs not in layer A`: without the
reroute the job falls through to the provider route, whose pairs were never stored.) Four
certified jobs change, so per CALC_2_FIXES the behaviour change is reverted. These four jobs
reroute through a decomposable prefix type, their family sizes equal the char-X multiplicity,
and the downstream checks certify them. Kept (comments only): the reroute_family doc comment now
states what the family is when L_a is decomposable and that a mismatch fails closed, and the
M5I_STAMP comment states that the stamp does not fingerprint how e_q was built and that a stale
key fails closed at `check_pi_iota`. Evidence: runs_calc2/replay14_failclosed.{log,jsonl} (the
change) and runs_calc2/replay14.{log,jsonl} (this branch).

The alternative in the item ("compare against cnt") was not attempted: accepting a family of cnt
< mx pairs would feed fewer pairs than the target multiplicity into the downstream selection,
which is new behaviour that the replay set does not exercise.

## Opt-in output features

### 15. `--table`

After each N block (also in `--file` mode):

```
  table (rows t, columns q):
  t\q  -14  -12  -10   -8   -6   -4   -2
    0                        1    1    1
    1
    2                   1    1
    3    1    1
```

Rows t from the minimum to the maximum (empty rows kept), columns q ascending in steps of 2 when
all q have one parity (else 1), entries right-aligned, blank for 0; `  table: 0` for the zero
polynomial.

### 16. `--convention knotjob`

Verified map, before implementing:
- data/README.md §"Relation to krhom's gradings" and tools/batch_check.py's `kh` check: krhom's
  N = 2 homology in bidegree (h, q) equals the KnotInfo table entry (h, −q) (`table_poly`:
  `(h, q) → (h, −q)`); KnotInfo's Khovanov data are computed by KnotJob.
- data/knots_upto12.jsonl 3_1 (braid `[1,1,1]`): `kh` = [[0,1,1],[0,3,1],[2,5,1],[3,9,1]],
  `kh_reduced` = [[0,2,1],[2,6,1],[3,8,1]].
- krhom default: `krhom 2 --braid aaa` → `t^0q^-3 + t^0q^-1 + t^2q^-5 + t^3q^-9`, reduced
  `t^0q^-2 + t^2q^-6 + t^3q^-8`.
- With the option: `t^0q^1 + t^0q^3 + t^2q^5 + t^3q^9` (= Bar-Natan's right-handed trefoil q + q³ +
  t²q⁵ + t³q⁹ = KnotInfo `kh`) and reduced `t^0q^2 + t^2q^6 + t^3q^8` (= `kh_reduced`).

So the map is exactly q → q⁻¹ with t unchanged. It is applied to the printed polynomial, the
table, the Euler characteristic, and the `--json`/`--csv` results; stdout has
`convention: knotjob (q -> q^-1)` after `coefficients:`, and JSON has `"convention":"knotjob
(q -> q^-1)"` (`"krhom"` by default). `--convention krhom` is accepted as the default; anything
else is exit 1.

Sample (`krhom 2 --braid aaa --convention knotjob --table`):

```
coefficients: Q
convention: knotjob (q -> q^-1)
time 0.01s
N = 2
  Poincare polynomial: t^0q^1 + t^0q^3 + t^2q^5 + t^3q^9
  total rank: 4
  Euler characteristic: -q^9 + q^5 + q^3 + q^1
  table (rows t, columns q):
  t\q  1  3  5  7  9
    0  1  1
    1
    2        1
    3              1
```

## Tests

New tests in tests/calc.rs: `exact2_memory_cap_and_shape_refusals` (1),
`n_below_the_color_short_circuits_refusals` (2), `batch_names_csv_quoting_and_incremental_outputs`
(3, 4, 6), `batch_outputs_survive_a_kill` (5: SIGKILL during the second line; the JSON array and the
CSV keep the first), `table_and_knotjob_convention` (15, 16). tests/c27_joint.rs `plan_max`
(13).

TESTS_PLACEHOLDER

## Acceptance runs

DEMO_PLACEHOLDER

H_PLACEHOLDER

RANDCHECK_PLACEHOLDER

## Commits

COMMITS_PLACEHOLDER
