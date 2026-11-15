Review of tools/fetch_library.sh, data/libraries.tsv, data/libstatus/*, tools/randcheck.py, tools/dev/k4m_status.py and docs/CALCULATOR.md at origin/colored-calc 1a9c8de. This was a read-only review: I ran nothing from the repo, used only read-only git, and parsed one JSON file once with `python3 -I`.

## Findings (most severe first)

### 1. [Medium] randcheck turns krhom exit 4 into SKIP, so certification and install failures can still exit 0
- **Where:** tools/randcheck.py:237 maps exit code 4 to "SKIP(library)". run_example (505-507) keeps it as the example status. main (681) exits 1 only on FAIL or ERROR.
- **What exit 4 means on the box route:**
  - a stamp mismatch;
  - a registry mismatch ("does not match this krhom");
  - a damaged install;
  - a map rejected on load: "loaded target rejected" and "not in layer A" are mapped to Library errors at src/calc.rs:77-88 and 382-387, and the diag messages are at src/bin/krhom.rs:520-536.
- **Why this is not a coverage gap:**
  - In T3, the l222_n4 registry status is `complete`. Any exit 4 there is a broken install or a rejected map. krhom.rs:536 says so itself ("damaged install: ... complete but a map is missing or rejected").
  - In T4, words are pre-filtered to 0 missing jobs by m5i_cover (randcheck.py:550-570, 463).
- **Failure scenario A:** `randcheck.py --seed S --no-anchor` or `--only T3-xx`. The anchor is skipped (randcheck.py:530, 620). Use a krhom whose route stamp or compiled-in registry differs from the installed library. The `--verify-only` self-test (526) uses the checkout's registry and does not check stamps, so it still passes. Every T3 and T4 box call exits 4, and the run reports `PASS <T2 count> FAIL 0 ... SKIP n` with exit 0.
- **Failure scenario B:** even with the anchor, a load-time re-certification rejection of a map that 4_1 does not use shows up only as SKIP for the affected random words, with exit 0.
- **Related:** TIMEOUT is also excluded from the exit code (681), so a run where every box call times out exits 0.
- **Fix:** treat exit 4 as FAIL or ERROR for T3, and for T4 words that passed the coverage filter. Or exit non-zero whenever SKIP or TIMEOUT is greater than 0.

### 2. [Low-medium] fetch_library.sh: `set -e` is inactive for the whole install path
- **Cause:** line 211 calls `install_row "$row" || rc=$?`. Bash ignores errexit inside a function called on the left of `||`, including ensure_commit and verify_installed.
- **Unguarded steps:**
  - 129 `mkdir -p "$root"`
  - 130 `staging="$(mktemp -d ...)"`
  - 131 `mkdir`
  - 159 and 168 `cp`
  - 172 `awk | sort`
  - 197-200 `mv` / `mv` / `rm -rf`
- **Scenario: `--root` or `$KRHOM_LIB_DIR` not writable.** mktemp fails, so `staging=""`. Then:
  - Line 131 runs `mkdir -p /chunks /lib/layer_a /lib/tmpl /lib/MANIFESTS`.
  - As a non-root user, the script dies at 137 with a misleading "no MANIFEST for c2/... in <commit>".
  - As root (typical in cloud containers) with a read-only lib mount, it writes /MANIFEST.N, /chunks/N, /owners and /final. It extracts about 5k files into /lib/layer_a and /lib/tmpl (/usr/lib on merged-usr systems), plus /lib/costs.txt and /lib/LIBRARY.txt. It then runs `mv "/lib" "$dest"` (199).
  - That mv fails only because the destination is unwritable. The EXIT trap (95) cannot clean up because `staging` is empty.
- **Scenario: mv failures.**
  - If 197 `mv "$dest" "$aside"` fails, the new copy ends up nested at `$dest/lib`, and line 202 verifies the old copy.
  - If 197 succeeds and 199 fails, line 200 still runs `rm -rf "$aside"` and deletes the old install.
- **Fix:** add `|| die` to mktemp, mkdir and mv, and assert `[ -n "$staging" ] && [ -d "$staging" ]`.

### 3. [Low] fetch_library.sh:146-147: the tar-entry checks fail open under `set -o pipefail` (line 19)
- **Line 147:** `if echo "$names" | grep -Evq RE; then die` — grep -q exits at the first bad line. If the name list is larger than the pipe buffer, `echo` gets SIGPIPE, the pipeline returns 141, and the `if` takes the no-die branch.
- **Sizes in data ee548d2:** c2/l2222_n4/f31 lists layer_a 2575 + tmpl 2607 entries, and c2/l2222_n4_c25/f31 lists 2612 + 2777. At about 31 bytes per name that is roughly 160 KB.
- **Line 146:** the same happens once awk has flushed its first output buffer.
- **Effect:** a disallowed entry in the first part of a large chunk passes the explicit check.
  - A leftover symlink is still caught by the count check at 175 (`find -type f`). A hardlink entry is not.
  - The extraction then relies only on tar's defaults (no `..`, absolute paths stripped).
- **TOCTOU:** none. The check and the extraction read the same sha-verified file `$staging/chunks/$i` inside a 0700 mktemp directory.
- **Impact is defense-in-depth only:** the chunks are pinned by sha256 in the registry.
- **Fix:** run `grep -Ev RE <<<"$names"` into a variable and test that it is non-empty, or use `! grep -Evq RE <<<"$names"`. A here-string has no separate writer process.

### 4. [Low] Installed status.txt and costs.txt are never re-checked against the registry
- **What happens:**
  - verify_installed (fetch_library.sh:71) compares only format, id, lambda, n, field, provider, data_commit and layers.
  - status.txt is copied from the checkout (168), and its source is not recorded.
  - costs.txt is checked only against the hash LIBRARY.txt recorded for it (78-79).
- **Scenario:** a commit corrects data/libstatus/l2222_n4_*.txt, or changes the registry `costs` column, but keeps the same data_commit and layers. Existing installs still print "already installed" (125) and pass `--verify-only`.
  - krhom's diagnosis then reads the stale copy (calc.rs:285-288; krhom.rs:541-553). `missing()` can refuse with exit 4 instead of trying the template-recompute retry (krhom.rs:555-560).
  - A stale costs.txt changes the route stamp (calc.rs:396-398), which gives exit 4 with a stamp-mismatch message.
- **Fix:** record the status and costs sources (and the status hash) in LIBRARY.txt and compare them in verify_installed.

### 5. [Low] Some randcheck self-test failures exit 1, not the documented 2
- **Where:** the convention self-test call (randcheck.py:519) is not inside a try. `wang()` uses `check=True` (308).
- **Effect:** a krhom error or a wang.py failure ends with a traceback and exit 1. docs/CALCULATOR.md:227 and 263, randcheck.py:35 and spec E all say exit 2 for a self-test or harness failure. A wrapper that separates 1 (FAIL) from 2 (harness failure) gets the wrong class.

### 6. [Low, doc vs code] Per-call timeouts in docs/CALCULATOR.md:263-264 are wrong
- **Docs say:** "T2 120 s, T3 and T4 2400 s".
- **Code uses:**
  - 600 s for every O8 uncolored call (randcheck.py:300);
  - 300 s for the T2 `--box-live` call (366);
  - 120 s for the T3 and T4 2-strand core calls (393, 474).

### 7. [Low, already disclosed] O11 has no Hopf |n| = 2 closed form
- **Where:** randcheck.py:312-327 returns None for even |n| ≥ 2.
- **Why it matters:** CALC_SPEC_review.md:226 (correction 11) requires the Hopf form, and CALC_DECISIONS.md:71-74 makes all corrections mandatory. CALC_1_report.md:246 discloses the omission.
- **Effect:** T2 words with |writhe| = 2 get no closed-form check. O6, O8 and O13a still apply to them.

## Checked and found correct
- **Oracles:**
  - Markov with no shift: interior σ2/σ3 insertion is a conjugate of an end stabilization.
  - Isotopy moves: braid_relation (155-163) is valid for both signs.
  - Mirror (t,q)→(−t,−q): correct over a field.
  - O8 duality: Λ²C³ ≅ (C³)* ⊗ det, and the closure of rev(w) is the reversed link.
  - qbinom2 = Gaussian [N 2] in q².
  - O11: Hopf |n| = 2 is not covered (finding 7). Unknot [N 2]; unlink [N 2]² unreduced and [N 2] reduced; wang.py already converts to ours = Wang(1/t,1/q) (wang.py:171); the m = 1 self-test anchors the convention; mirroring for writhe < 0 is right.
  - Reduced comparisons are made only for knots, or for links with the same marked component.
  - ncomp is correct.
- **Vacuous passes:** cmp_check (267) skips an N missing on either side. That is not reachable today, because krhom always emits every requested N: zero_N insertion at krhom.rs:442 and the box post-check at calc.rs:440-444. O3 requires `closures > 0`. The documented `--inject-fault --only T3-01` run is caught by O5 and O12 at N = 4.
- **fetch_library:**
  - Every chunk is checked against both the registry and the MANIFEST (141-142).
  - All 10 registry hashes and file counts match the MANIFESTs in ee548d2. Every layer is a single `part-000`.
  - The fetch-by-SHA fallback cannot install a different commit: every read uses `$commit:`, and `cat-file -e "$commit^{commit}"` is checked at 104 and 112.
  - `--force` builds the full staging copy before swapping it in.
- **Portability (macOS):** /sbin/sha256sum is the BSD tool in GNU mode and supports `-c`, `--quiet` and stdin (man page checked). bsdtar 3.5.3 and bash 3.2 constructs are fine.
- **data/libstatus:** f31 and f32b have identical keys, 1254 lines each, 1236 ok, and identical 18-job fail sets, matching CALCULATOR.md:102.
- **k4m_status.py:** no defect found that affects results; it feeds diagnostics only.
- **Other CALCULATOR.md claims** checked against the code (exit codes, seeds, memory caps, simplify guard, Q wording, stamp preflight) are consistent.

Files reviewed:
- /Users/existair/Documents/KRcalculator/.claude/worktrees/calc/tools/fetch_library.sh
- /Users/existair/Documents/KRcalculator/.claude/worktrees/calc/tools/randcheck.py
- /Users/existair/Documents/KRcalculator/.claude/worktrees/calc/docs/CALCULATOR.md
- /Users/existair/Documents/KRcalculator/.claude/worktrees/calc/data/libraries.tsv
- /Users/existair/Documents/KRcalculator/.claude/worktrees/calc/data/libstatus/
- /Users/existair/Documents/KRcalculator/.claude/worktrees/calc/tools/dev/k4m_status.py