Verification of the 'cli' findings (worktree calc @ 1a9c8de, read-only; nothing was run)

Result: 5 of the 6 findings are REAL and 1 is NOT-REAL. One is major, the rest are minor. None is a blocker: no route prints an uncertified or wrong polynomial for a valid input.

**1. REAL, major. On the exact 2-strand route a cap or shape refusal panics and exits 2, and `--mem-cap-mb` is ignored there.**
- How it happens:
  - krhom.rs:399 builds `em`, but `em` only wraps the `Err` from `khr_singular_exact`, which only carries `xscan` errors (xengine.rs:2013).
  - The per-N work is `scan.close` (xengine.rs:1936-1937) → `close_guarded` → `close_exact_colored` (xengine.rs:1986).
  - `close_exact_colored` calls `try_close_exact_colored(...).unwrap_or_else(|e| panic!(..))` (cclose.rs:1280-1282).
  - So the cap error from `exact_precheck` (cclose.rs:1243-1251) becomes a panic. The same happens to the `check_alg_shape` "top weight > 255" error (cclose.rs:847-848), which fires for N ≥ 89 at color 3, N ≥ 130 at color 2 and N ≥ 257 at color 1.
  - `guard` (krhom.rs:241-248) turns the panic into `CalcErr::Engine("... panicked: ... (please report this)")`, which is exit 2 (calc.rs:45), not `Resource`/5 (calc.rs:48).
- The flag: `KRHOM_MEM_CAP_MB` is set from `o.mem_cap_mb` only in `box_route` (krhom.rs:480-491). On exact2 the 768 MB default applies (cclose.rs:739, 749-751), so `--mem-cap-mb` does nothing.
- Example: `krhom 30 --braid "1 1 1" --colors 3`. Labels are [3,3] and dim A_λ = C(30,3)² ≈ 1.65e7. The estimate is at least 2·dim·150 B ≈ 4.9 GB, above 768 MB, so the run panics and exits 2.
- Docs it breaks: CALCULATOR.md:89 ("any N", colors 1..3) and the exit table at CALCULATOR.md:174 (a cap refusal is exit 5; raise `--mem-cap-mb`). CALCULATOR.md:205-206 does say krhom sets the cap only on the box route.
- Not real: the knock-on that randcheck would log a false FAIL. randcheck uses only color 2 with N ≤ 6 on 2 strands (tools/randcheck.py:339, 375, 470), where dim A_λ ≤ 225 and the cap is never hit.

**2. REAL, minor. A process abort loses all `--json`/`--csv` output of a batch.**
- `jsons` and `csv_rows` are kept in memory (krhom.rs:657-658, 698) and written only after the loop (krhom.rs:700-710).
- `guard` uses `catch_unwind`, which does not catch these:
  - an allocation failure, which aborts by default (there is no custom hook, and the Cargo.toml release profile has no panic setting);
  - SIGKILL or the OOM killer.
- CALCULATOR.md:208-209 recommends `ulimit -v`, which makes this abort an expected way to fail.
- After an abort the finished items survive only on stdout, which is line-buffered.
- Overstated: "breaks D3". No implementation can continue a batch after an abort, and D3 (CALC_DECISIONS.md:47-49) does not require incremental writes. The real cost is that the JSON with the certification details is lost.

**3. REAL, minor. The Q5 rule "zero on every route, no engine" is not fully honoured.**
- (a) On exact2 nothing short-circuits an empty `run_ns`. krhom.rs:388-400 still calls `two_primes`, so `xscan` runs in full for both primes just to print 0. The box route (krhom.rs:495-498) and box-live (krhom.rs:409-410) do short-circuit. The waste grows with word length, and an `xscan` error here would turn a trivial 0 into exit 2.
- (b) Refusals fire before the zero shortcut:
  - `krhom 1,2 --braid aBaB --colors 3`: labels [3,3,3] and 4_1 stays on 3 strands. `reg_for` is None because the registry holds only l222/l2222 (calc.rs:171-172, 673-674), so it exits 1 at krhom.rs:460-462.
  - `krhom 2,3 --braid "1 1 1" --colors 4` exits 1 at krhom.rs:393-395.
  - In both cases the "homology is 0" note at krhom.rs:337 has already been printed. This contradicts D4 Q5 (CALC_DECISIONS.md:64-65) and CALCULATOR.md:94-95.

**4. REAL, minor. The CSV name field is not escaped.**
- `nm` is written raw at krhom.rs:677 and 680. Names come from the text before the first `;` (krhom.rs:217-218), which may contain commas and quotes.
- LinkInfo names with commas ship in the repo, e.g. `"L6a4{0,1}"` in data/links_upto11.jsonl.
- A line `L6a4{0,1}; ...` produces 7-column rows, so N, t, q and rank are misparsed. JSON is unaffected because `jstr` escapes.

**5. NOT-REAL. The empty braid is accepted on purpose.**
- `parse_braid` has an explicit branch that returns the empty word, the 1-strand unknot (calc.rs:550-551). `check_input` accepts an empty word (xengine.rs:1713-1732).
- The output labels it: `crossings 0` and `braid: (empty)` (krhom.rs:330-331, `sigstr` at 234-235), and the JSON records input text "". The homology printed is correct for that input.
- At most this is an input-validation suggestion, for example warning on an empty word after `NAME;` in `--file`.

**6. REAL, minor. A batch line that fails early has no name on stdout.**
- `name:` is printed at krhom.rs:327-329, after the early returns at krhom.rs:284, 287, 294-295, 304, 322 and 325.
- The batch error line at krhom.rs:688 prints only `error: exit N (kind)`. Only stderr (krhom.rs:686) and the JSON carry the name.
- tests/calc.rs:190-193 confirms it: the test checks stdout names only for the good lines, and the bad line ("bad") only on stderr.
- This contradicts D3 (CALC_DECISIONS.md:48), which asks for one block per line headed by the name.

**Spot check of the "checked and found correct" notes.** I agree that a box item's leftover `KRHOM_MEM_CAP_MB` and `KRHOM_COSTS_ANY_HOST` have no practical effect. A successful box item implies all N ≤ 4 (krhom.rs:464) and the same colors on all 2s for the whole batch. So later items are tiny color-2 exact2 or box runs, and engine2 cannot be reached.

Severity summary: major: #1. Minor: #2, #3, #4, #6. Not-real: #5, plus the randcheck sub-claim of #1 and the D3 sub-claim of #2.

Files cited: /Users/existair/Documents/KRcalculator/.claude/worktrees/calc/src/bin/krhom.rs, /Users/existair/Documents/KRcalculator/.claude/worktrees/calc/src/xengine.rs, /Users/existair/Documents/KRcalculator/.claude/worktrees/calc/src/cclose.rs, /Users/existair/Documents/KRcalculator/.claude/worktrees/calc/src/calc.rs, /Users/existair/Documents/KRcalculator/.claude/worktrees/calc/docs/CALCULATOR.md, /Users/existair/Documents/KRcalculator/.claude/worktrees/calc/docs/colored_design/calc/CALC_DECISIONS.md, /Users/existair/Documents/KRcalculator/.claude/worktrees/calc/tests/calc.rs, /Users/existair/Documents/KRcalculator/.claude/worktrees/calc/tools/randcheck.py, /Users/existair/Documents/KRcalculator/.claude/worktrees/calc/data/links_upto11.jsonl