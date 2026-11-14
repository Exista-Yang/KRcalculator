Review of src/calc.rs, src/bin/krhom.rs, tests/calc.rs and tests/xengine.rs (worktree calc @ 1a9c8de; read-only, nothing was run). I found six defects, listed from most to least severe. None of them makes the calculator print an uncertified result.

**1. MEDIUM: on the exact 2-strand route, a memory-cap refusal exits 2 ("please report this") instead of 5, and `--mem-cap-mb` has no effect there.**
- Where: krhom.rs:399-400 → xengine.rs:2012-2015 → `close` → `close_guarded` (xengine.rs:1977-1989) → `close_exact_colored` (cclose.rs:1280-1282).
- Why:
  - `close_exact_colored` panics on the cap error from cclose.rs:1243-1251 (and on the "top weight" shape error).
  - The `em` closure (krhom.rs:399) only maps `Err` strings that come out of `xscan`. The cap refusal arrives as a panic, so `guard` (krhom.rs:241-248) turns it into `CalcErr::Engine("engine (F_2147483647) panicked: close_exact_colored: ... above the memory cap of 768 MB ... (please report this)")`, exit 2.
  - `KRHOM_MEM_CAP_MB` is set only inside `box_route` (krhom.rs:480-491). On this route `--mem-cap-mb` is ignored and the 768 MB default applies.
- Docs vs behaviour: CALCULATOR.md:88 promises "any N" for colored 2-strand inputs. CALCULATOR.md:174 says a cap refusal is exit 5 and the fix is to raise `--mem-cap-mb`.
- Example: `krhom 30 --braid "1 1 1" --colors 3`. Here dim A_λ = C(30,3)² ≈ 1.6e7, the same size tests/cclose.rs:794-798 shows is above the default cap. Result: exit 2 with a "certification failure" style message.
- Knock-on: tools/randcheck.py:237 counts exit 2 as FAIL, so the checker would log a false failure.

**2. MEDIUM-LOW: a batch loses all its `--json`/`--csv` output if the process aborts.**
- Where: krhom.rs:657-710 collects JSON and CSV in memory and writes them only after the last item.
- The problem: `guard` catches unwinding panics only. Allocation failure aborts the process by default, which is exactly what happens under the `ulimit -v 14000000` that CALCULATOR.md:207-209 recommends. SIGKILL, a timeout or Ctrl-C also bypass it.
- Example: `ulimit -v 14000000; krhom 2,3,4 --file knots4.txt --colors 2 --json out.json`. If line 10 runs out of memory, out.json is never written although lines 1-9 (hours of work) finished.
- This also breaks D3 ("the batch continues") for this failure class. Only stdout keeps the finished blocks.

**3. LOW: "N below the color gives zero on every route, without running an engine" (D4 Q5, CALCULATOR.md:94-95) is not fully honoured.**
- (a) Exact 2-strand route: when every requested N is below the color, `run_ns` is empty but krhom.rs:388-400 still calls `two_primes`. That runs a full `xscan` for both primes just to print 0. The box route (krhom.rs:495-498) and box-live (409-411) do short-circuit.
- (b) Refusals fire before the zero shortcut. Both of these exit 1 even though the answer is trivially 0 and stderr has already printed the "is below the color ... 0" note at krhom.rs:337:
  - `krhom 1,2 --braid aBaB --colors 3` → "no library" (krhom.rs:460-462);
  - `krhom 2,3 --braid "1 1 1" --colors 4` → "colors 1..=3" (krhom.rs:393-395).

**4. LOW: the CSV name field is not escaped.**
- Where: krhom.rs:675-681 writes `nm` raw. Names come from the `NAME;` prefix (krhom.rs:217-219).
- Example: LinkInfo-style names such as `L6a5{0,1}`, or any name containing a comma or quote, produce rows with extra columns, so N/t/q/rank are misparsed. The JSON side is fine: `jstr` escapes everything.

**5. LOW: an empty braid is accepted silently.**
- Where: `parse_braid("")` returns `Ok(vec![])` (calc.rs:550-552), which becomes a 1-strand unknot (krhom.rs:285).
- Example: a `--file` line `6_1;` with the word missing prints unknot homology under the name 6_1 and exits 0.

**6. LOW: failed batch lines have no name on stdout.**
- Where: `name:` is printed at krhom.rs:327-329, after the early returns at 284-326 (parse, strands, color count, `--check-library`/`--reduced` refusals). The error line at krhom.rs:687-689 has no name either.
- Effect: such a line's stdout block is just `error: exit 1 (unsupported)`, so stdout cannot tell which line failed. D3 asks for one block per line headed by the name; only stderr names it.

**Checked and found correct:**
- **Routing vs D1-D4.**
  - Uncolored inputs go to engine2 on F31 and F32b, seeds 7 and 11.
  - Colored inputs on ≤ 2 strands after simplification go to xengine, two primes.
  - Colored 3-4 strands, uniform color 2: N > 4 is checked against all N before anything loads; split braids and the refused env switches give exit 1.
  - Everything else colored exits 1 with the support matrix.
  - Cube runs exactly over Q; cable now runs both primes and compares.
- **The Q claim.**
  - Every F_p route goes through `compare_primes` over the computed N, and the printed map is F31's map for exactly those N plus the zero N.
  - `--single-prime` relabels the coefficients line on all routes.
  - `--check-library` prints no polynomial.
  - A disagreement prints nothing and exits 3.
- **Panics.** No partial result is printed: N blocks come only from an `Ok` outcome, and stdout before a failure holds only header lines. No `panic=abort` profile exists, the thread pools are scoped, and library code never calls `process::exit`.
- **Simplification guard.** `from_signed` always numbers the component at bottom position 0 as 0, and `simplify` keeps component ids, so the reduced-link guard and the per-component color mapping after simplification are correct.
- **Template-recompute retry.** It cannot hide a library gap or weaken certification:
  - `replay_only` has a single consumer (smapio.rs:350). With it off, a missing hom basis is recomputed deterministically and still checked against the stored template hash (smapio.rs:353-357).
  - The provider is `NoProvider`, so no layer-A map can be computed at run time.
  - For an incomplete library the missing-job diagnosis runs first. A complete library that fails again is reported as a damaged install (exit 4).
- **Other checks.** The K1 χ post-check, the K4 per-layer stamp preflight and JSON escaping are correct. The `KRHOM_MEM_CAP_MB` value set by a box item stays set for later batch items, but this has no practical effect: a box item forces color 2 and N ≤ 4 on the whole batch.
