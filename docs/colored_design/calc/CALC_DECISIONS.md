# Calculator: binding decisions (orchestrator, 2026-11-14)

Read in this order: this file, then CALC_SPEC_review.md, then CALC_SPEC_v1.md. This file wins over both where they
conflict, and the review's "Required corrections" win over the spec. Background: FIELDS_review.md (which base fields
the code can support) and KNOTJOB_reference.md (the UX reference the user named).

The user's requirement, in their words: a calculator — "给一个 strand + braid word + color，就能输出结果，可以参考
knotjob 算 Khovanov homology 的能力"; validate on some random examples, never by sweeping knot tables; "不要写这两个
大素数接口，现阶段不写这个参数了，就默认Q".

## D1. No base-field parameter; every result is over Q

- Remove `--field` and `--single-prime` from the usage text, the docs and the README. If existing tests or scripts need
  them, keep them as hidden developer flags (not printed by usage, listed only under "developer options" in
  docs/CALCULATOR.md). `--seed-f32b` (needed by the checker's oracle O2) is likewise a hidden developer flag.
- stdout replaces the `fields: ...` line by exactly `coefficients: Q`. How Q is certified (two primes agree; the
  per-prime seeds, timings and checks) goes to stderr and to `--json`, not to the main output. Update the tests that
  parse the old line.
- Every route the calculator picks must be certified over Q:
  - colored singular routes (2-strand exact, 3/4-strand box): two primes must agree (as today); exit 3 otherwise.
  - the uncolored Soergel route (engine2): NEW — run F31 and F32b (precedent: src/bin/twoprime.rs, seeds 7 and 11) and
    compare with the shared `compare_primes`; exit 3 on disagreement.
  - `--colored-method cube`: exact over Q (`khr_full_cube::<Q>`), no primes involved.
  - `--colored-method cable` runs over a single prime today: either make it run both primes and compare, or make it a
    hidden developer method that prints a clear "F_p only, not certified over Q" warning. Never suggest it to users.
- README/docs wording for the Q claim: "certified over F_p for two large primes; equal over Q unless both are torsion
  primes" (C_sound wording, FIELDS_review.md §4). Small characteristic is not supported (FIELDS_review.md §5).

## D2. Automatic routing — `--colored-method` is not needed

The user should type strands + braid + colors (+ N list, optional `--reduced`) and get an answer:
- all colors 1 → uncolored Soergel engine (any number of strands), two primes;
- colored, ≤ 2 strands (after simplification) → the existing exact singular path;
- colored, 3–4 strands, uniform color 2, N ≤ 4 → the box route with the installed library;
- any other colored input → exit 1 with the support matrix; mention `--colored-method cube` (exact over Q, tiny inputs
  only) as the advanced fallback. This overrides spec Q6 (default unchanged) — the default for colored input is now the
  singular route whenever it is supported.
- `--colored-method cable|cube|singular` remains as an advanced override.

## D3. KnotJob-style input and output (P0 unless marked)

- `--braid` accepts signed integers ("1 -2 1 -2") OR KnotJob letters ("aBaB": a..i = sigma_1..sigma_9, A..I their
  inverses, letters only, no spaces). Detect the format; reject mixtures with a clear message.
- `--strands` becomes optional: default = max |generator| + 1 (KnotJob convention). If given, it must be ≥ that.
- `--colors c` broadcasts a single value to all components (spec); per-component lists stay.
- N list as today (`krhom 2,3,4 ...`).
- Batch (P1): `--file F`, one braid per line, optional `NAME;` prefix, blank lines and `#` comments ignored; the same
  colors / N / `--reduced` apply to every line; output one block per line headed by the name; `--json` then writes a JSON
  array. Errors on one line are reported and the batch continues (non-zero exit at the end).
- Output: keep the current N blocks (Poincare polynomial sorted by t then q, total rank, Euler characteristic).
  `--csv OUT` (P1): rows `name,N,reduced,t,q,rank`.
- PD input stays (`--pd`).

## D4. Decisions on the spec's open questions

- Q1 base: branch `colored-calc` from the newest origin/colored-2c4-unblock at session start (≥ 194a731); C2-7 is still
  pushing there, merge origin/colored-2c4-unblock again before the final push if it moved (the orchestrator reviews
  both). Do not edit src/leaves.rs or src/bin/libbuild.rs.
- Q2: follow the MANIFEST restore order (l2222_n4, c25, c24, c27); report final counts.
- Q3: layer A stays replay-only (never compute layer-A maps at run time). If a run fails only because a template is
  missing, template recomputation may be enabled (deterministic, does not affect soundness); report it.
- Q4 / review K3: installs are verified by a per-file sha256 list written into LIBRARY.txt; extra files (templates
  written later) are tolerated; `--verify-only` checks only the listed files.
- Q5 / review K14: N below the color → print the zero polynomial on every route (wedge^a C^N = 0 for a > N), with a note
  on stderr; do not run the engine.
- Q8 (resumable 4-strand runs): deferred; document the expected cost (6_1 unreduced 1.2–8 core-hours) and recommend
  `--reduced` for speed.
- Q9: document `ulimit -v`; no new dependencies.
- Q10: same seed (201) for both primes by default.
- Q11: split braids on ≥ 3 strands are refused (product formula is future work).
- All 15 "Required corrections" of CALC_SPEC_review.md are mandatory (K1 chi_checked, K2 classifier, K3 file lists,
  K4 stamp preflight per layer, K5/K6 cost figures and checker timeouts, K7 thread-aware memory cap, K8 explicit method
  in tests — now moot for the default route but keep tests explicit, K9 panics, K10 check-library, K11 fetch by SHA,
  K12 env hygiene, K13 JSON escaping, K14, K15).

## D5. Scope limits (unchanged from the spec)

No table sweeps; no new library builds; no (3,3,3); certification checks unchanged; no new crate dependencies. The
random checker (spec §E) with the review's fixes is the validation: default run of about 24 examples.

## D6. Docs

The session writes docs/CALCULATOR.md (reference). The orchestrator writes the top-level README.md (English) from it
afterwards; do not create README.md in this session.
