# CALC-2: fixes after the review of colored-calc (orchestrator, 2026-11-14)

Source: review_calc1/{cli,engine,leaves,tools}_{findings,verdict}.md (read-only review of origin/dev...origin/colored-calc
1a9c8de, every finding adversarially verified). Verdict: NO blocker; no path prints an uncertified or wrong polynomial.
Fix the items below on branch colored-calc2 (from origin/colored-calc). Keep CALC_DECISIONS.md binding. Do not change
certification checks; do not edit src/leaves.rs except where stated.

## Must fix
1. (major, cli #1) Exact 2-strand route: a closure memory-cap or shape refusal panics inside close_exact_colored
   (cclose.rs:1280-1282) and is reported as exit 2 "please report this". Make it a Resource error (exit 5) with the
   flag to raise, e.g. by using try_close_exact_colored on this route; and honour --mem-cap-mb on EVERY route (set
   KRHOM_MEM_CAP_MB before any engine work, not only in box_route). Same for the "top weight > 255" shape error:
   exit 1 (unsupported N) with a clear message. Update docs/CALCULATOR.md:88 ("any N") to the real limits.
2. (cli #3) N below the color: short-circuit to the zero polynomial BEFORE any refusal (no-library, colors > 3) and
   before running an engine (exact2 with an empty run list must not run xscan).
3. (cli #6) Batch (--file): print `name: NAME` before any early return, and include the name on the error line.
4. (cli #4) CSV: quote/escape the name field (RFC 4180).
5. (cli #2) Batch: write --json/--csv incrementally (one record appended and flushed per finished item, e.g. JSON Lines
   to OUT while running and the final JSON array at the end, or rewrite the array file after each item) so an abort
   (OOM under ulimit -v, SIGKILL) keeps the finished results.
6. JSON: `input` contains the key "braid" twice (string and array) — rename the string one (e.g. "braid_text").
7. (tools #1) randcheck: exit 4 from krhom on T3, and on T4 words that passed the m5i_cover filter, is FAIL (not SKIP).
8. (tools #5) randcheck: harness/self-test failures exit 2 as documented (wrap the convention call, wang(), the cargo
   build and m5i_cover calls).
9. (tools #2, #3) fetch_library.sh: errexit is disabled inside `install_row ... || rc=$?`; add `|| die` to every
   install step (mktemp, mkdir, cp, mv, writes) and refuse an empty staging path; fix the tar-entry check so it cannot
   fail open under pipefail (write the listing to a file and grep the file, or use a here-string).
10. (tools #4) fetch_library.sh / LIBRARY.txt: record the sha256 of status.txt and costs.txt sources and re-check
    them against the registry on "already installed" and --verify-only.

## Should fix
11. (tools #6) docs/CALCULATOR.md randcheck timeouts: state the real per-call timeouts.
12. (tools #7) randcheck O11: add the Hopf |n| = 2 closed form (tests/xengine.rs hopf_formula).
13. (leaves #1) gp0 census (`--steps`): count the replacement final-word score `sw` in `mx`/`over` instead of 0 and a
    decrement; tests/c27_joint.rs plan_max: take the max with the replacement word's score. Do NOT change leaves.rs
    planning behaviour (certified libraries depend on it).
14. (engine #1, fail-closed) sindec reroute_family: if the chosen source type `a` is decomposable, return None (or
    compare against cnt) so the multiplicity check cannot falsely reject; fix the doc comment at sindec.rs:1169 and the
    M5I_STAMP comment (sindec.rs:51-53). This must not change any certified job: re-run the M5i replay jobs of
    data/libstatus for both primes (libbuild one --replay) for the decomposable-source reroutes (109, 1789, 1812, 2052,
    2086, 2275, 2298) and confirm "ok" — if any changes, revert this item and report.

## Opt-in output features (user-facing; default output unchanged)
15. `--table`: after each N block print the bigraded table (rows t ascending, columns q ascending, entries = dims,
    blank for 0), KnotJob-style. Also accepted in --file mode.
16. `--convention knotjob`: print the Poincare polynomial, table and Euler characteristic with q -> q^{-1} (t unchanged),
    so that the right-handed trefoil (`aaa`) at N = 2 reads q + q^3 + t^2q^5 + t^3q^9 (Bar-Natan, KnotJob). VERIFY the
    map first against Bar-Natan's trefoil and against tools/batch_check.py's KnotInfo comparison (its `kh` check
    documents the table convention); if the map is not exactly q -> q^{-1}, implement the verified map or skip this item
    and report why. Print a header line `convention: knotjob (q -> q^-1)` when active; JSON records it.

## Acceptance
- KRHOM_THREADS=4 cargo test --release: 0 failed (new tests for 1-7, 15, 16).
- Re-run H3-H8 of CALC_1_report.md quickly (outputs identical to runs_calc/), and randcheck --quick (seed 20261115) PASS.
- Demonstrate items 1, 2, 3, 5 with small commands (record stdout/stderr/exit in runs_calc2/).
- Write docs/colored_design/calc/CALC_2_report.md.
