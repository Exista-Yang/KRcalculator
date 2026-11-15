Verification of the 'tools' findings at origin/colored-calc 1a9c8de. This was read-only: I read files and used read-only git (`git show ee548d2:...MANIFEST.txt`). Nothing was built or run.

## 1. randcheck turns krhom exit 4 into SKIP: REAL, minor
- **Code path:**
  - randcheck.py:237 maps exit 4 to "SKIP(library)".
  - run_example (505-507) keeps that status.
  - main returns 1 only when FAIL or ERROR is non-zero (681).
  - On the box route, exit 4 covers more than "word not covered". It also covers a registry mismatch ("does not match this krhom", calc.rs:275), a stamp mismatch (calc.rs:320-334) and "damaged install / rejected on load" (krhom.rs:520-537). The library-class strings come from calc.rs:77-88 and 382-387.
- **Scenario A (confirmed):**
  - Run with `--no-anchor` or `--only` (randcheck.py:620). The self-tests are then only the 2-strand convention check (519) and `fetch_library.sh --verify-only` (526). That check uses the checkout's data/libraries.tsv, not the registry compiled into the krhom binary (calc.rs:23).
  - Use a stale or foreign `--krhom` whose compiled registry or stamp differs from the install. Every T3/T4 box call then exits 4 and becomes SKIP. Only T2 is validated, and the run still exits 0.
  - With the anchor on, the l222_n4 mismatch is caught, because the anchor's RunFail goes to errs and the run exits 2 (539-541). An l2222_n4 mismatch (T4) is never caught. A T4 SKIP also looks exactly like the normal "no covered candidate" SKIP (464-467).
- **Scenario B (reachable):** l222_n4 is `complete`, so `missing()` returns None (krhom.rs:542-543). A krhom regression that falsely rejects a layer-A map that 4_1 does not use then ends as "damaged install", exit 4, and SKIP.
- **Severity:** minor. krhom itself fails closed, so no wrong polynomial is printed or certified, and the SKIP count and error text appear in the summary. The defect is a validator that can exit 0 without having validated any colored output.
- **Fix:** treat exit 4 as FAIL for T3, and for T4 words that passed the m5i_cover filter.
- **Sub-claim "TIMEOUT excluded from the exit code": NOT-REAL.** This is the specified behavior: CALC_SPEC_v1.md:439 says "0 if no FAIL or ERROR", and :427 says TIMEOUT is never PASS.

## 2. fetch_library.sh: errexit is off for the whole install path: REAL, minor
- **Mechanism (confirmed):** `install_row "$row" || rc=$?` (fetch_library.sh:211). Bash ignores `set -e` inside a function called on the left of `||`. The steps at 129-131, 159, 168, 172 and 197-200 have no `|| die`.
- **Scenario (confirmed):** run as root on Linux (cloud containers) with `--root`/`$KRHOM_LIB_DIR` pointing at a read-only mount or a regular file.
  - mktemp at 130 fails, so `staging=""`.
  - 131 creates /chunks and /lib/{layer_a,tmpl,MANIFESTS}.
  - 136-158 write /MANIFEST.i and /chunks/i and extract several thousand files into /lib (/usr/lib on merged-usr). 157 writes /owners and 172 writes /final.
  - 163 and 187-195 write /lib/costs.txt and /lib/LIBRARY.txt; 168 writes /lib/status.txt for l2222_n4.
  - 199 `mv /lib $dest` fails only because the destination cannot be created. 202 then dies with exit 2.
  - The EXIT trap (95) cleans nothing, because staging is empty. System files are not overwritten, since archive names are limited to layer_a/ and tmpl/.
- **Non-root:** the script dies at 137 with a misleading "no MANIFEST" message and does no harm.
- **mv sub-scenarios at 197-200:** these need contrived permission or mount layouts: a non-writable $root/$id with a writable $dest, or $root/$id on another mount. The only loss is an old copy that either failed verification or was being replaced on purpose with --force.
- **No path yields a verified bad install:** 202 still re-verifies, and tar (158) and costs (163) have `|| die`.

## 3. Tar-entry checks fail open under pipefail: REAL, minor (defense in depth)
- **Mechanism (confirmed):** fetch_library.sh:19 sets `set -o pipefail`.
  - At 146-147, `grep -q` exits at the first bad line. The upstream echo/awk then gets SIGPIPE (141), or EPIPE with status 1 if SIGPIPE is ignored. With pipefail the pipeline is non-zero, so the `if` takes the no-die branch.
- **Reachability:** listing sizes are confirmed from the ee548d2 MANIFESTs.
  - c2/l2222_n4/f31: layer_a 2575 and tmpl 2607 entries. c2/l2222_n4_c25/f31: 2612 and 2777. That is about 160 KB of names, above the 64 KB pipe buffer.
  - Even l222_n4 (704 entries, about 22 KB) is above macOS's default 16 KB pipe buffer.
  - So a disallowed entry near the start of the listing passes the explicit check.
- **Impact:** low. Chunks are sha256-pinned against the in-repo registry and the MANIFEST (141-142), so exploiting this needs a compromised repo. A leftover symlink is still caught by the `find -type f` count check (175-176), and bsdtar and GNU tar refuse to extract through symlinks.

## 4. Installed status.txt and costs.txt are not re-checked against the registry: REAL (latent), minor
- **Confirmed:**
  - verify_installed checks only the 8 keys at fetch_library.sh:71, plus the costs hash against LIBRARY.txt's own record (78-79).
  - status.txt is copied from the checkout (168) and its source is not recorded.
  - krhom's check_library_kv (calc.rs:249-263) does not compare costs or status either.
- **Trigger:** none today. It needs a future commit that edits data/libstatus/*.txt, or repoints the registry's costs column, while keeping data_commit and layers. Existing installs would then pass "already installed" (125) and `--verify-only`.
- **Effects are fail-closed only:**
  - The route stamp hashes the costs file content (sindec.rs:493, 502), so stale costs give an exit-4 stamp mismatch whose message advises --force.
  - A stale status file affects only diagnose_missing (calc.rs:497-522) and the `missing()` short-circuit (krhom.rs:541-560). It can refuse with exit 4 instead of trying the template-recompute retry.
  - No wrong results.

## 5. Self-test harness failures exit 1, not the documented 2: REAL, minor
- **Confirmed:** the convention call at randcheck.py:519 is outside any try. A RunFail, or a TypeError when no JSON is written (239), propagates out of main with a traceback and Python's exit 1.
- The same happens with `wang()` `check=True` (308), and in t4_coverage with the cargo build (554) and the `m5i_cover` call (564, `check=True`).
- **Docs and spec say 2:** docs/CALCULATOR.md:227 and :263, randcheck.py:35, CALC_SPEC_v1.md:367 and :439.
- **Impact:** the exit is still non-zero, so this is never a false pass. A wrapper that separates FAIL (1) from harness failure (2) misclassifies it.

## 6. Doc timeouts wrong: REAL (documentation only), minor
- docs/CALCULATOR.md:263-264 says "T2 120 s, T3 and T4 2400 s".
- The code uses:
  - 600 s for O8 uncolored calls (randcheck.py:300);
  - 300 s for T2 `--box-live` (366);
  - 120 s for T3/T4 2-strand core calls (393, 474).
- The 2400 s figure applies only to box calls (399, 427, 432, 475). This has no runtime effect.

## 7. O11 has no Hopf |n| = 2 closed form: REAL, minor (already disclosed)
- **Confirmed:** randcheck.py:312-327 returns None when |e| = 2. CALC_SPEC_review.md:226 (correction 11) requires the Hopf form, and CALC_DECISIONS.md:71-74 makes all corrections mandatory.
- **Disclosed:** CALC_1_report.md:246 states the omission.
- **What remains:** |writhe| = 2 T2 words still get O4, O6, O8 and O13a. The loss is coverage only, with no wrong output.

## Summary
- **Blockers:** none.
- **REAL, minor:** 1, 2, 3, 4 (latent), 5, 6 (documentation only) and 7 (disclosed spec deviation).
- **NOT-REAL:** the TIMEOUT-exit-code sub-claim of finding 1, which matches the spec.
- **Worth fixing first:** 1, since it is cheap, and 5, then 2 and 3, which are one-line `|| die` / here-string fixes.

## Files read
- /Users/existair/Documents/KRcalculator/.claude/worktrees/calc/tools/fetch_library.sh
- /Users/existair/Documents/KRcalculator/.claude/worktrees/calc/tools/randcheck.py
- /Users/existair/Documents/KRcalculator/.claude/worktrees/calc/docs/CALCULATOR.md
- /Users/existair/Documents/KRcalculator/.claude/worktrees/calc/src/calc.rs
- /Users/existair/Documents/KRcalculator/.claude/worktrees/calc/src/bin/krhom.rs
- /Users/existair/Documents/KRcalculator/.claude/worktrees/calc/src/sindec.rs
- /Users/existair/Documents/KRcalculator/.claude/worktrees/calc/data/libraries.tsv
- /Users/existair/Documents/KRcalculator/.claude/worktrees/calc/docs/colored_design/calc/CALC_SPEC_v1.md
- /Users/existair/Documents/KRcalculator/.claude/worktrees/calc/docs/colored_design/calc/CALC_SPEC_review.md
- /Users/existair/Documents/KRcalculator/.claude/worktrees/calc/docs/colored_design/calc/CALC_DECISIONS.md
- /Users/existair/Documents/KRcalculator/.claude/worktrees/calc/docs/colored_design/calc/CALC_1_report.md