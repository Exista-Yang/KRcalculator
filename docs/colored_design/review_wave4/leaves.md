# Review wave 4: M4c A2/A3 Leaves closure (branch m4c-closure, b266a1a; diff 595068c..HEAD)

Lens: can Leaves accept a closure family that is not a basis of Hom(R, B_p) / Hom(B_p, R)? Is acceptance only
through `certify_closure`? Are Elem constants validated, and is Elem off in production? Does `ClosureLeaves`
change non-closure behaviour? Do budget and admission rules only affect what is tried?

Read-only review. No cargo, no binaries. Paths are relative to /Users/existair/Documents/KRcalculator/.claude/worktrees/revclo.

## Verdict

No path to a false certification or a wrong closure basis was found. Acceptance is gated by `certify_closure`
everywhere: inside the provider, in sindec after every provider call, and after every layer-A load. The Elem gate
holds twice: `validate_elem` runs before use, and `assert_genuine` rejects any Elem node whose key was not
certified in this arena. Elem is off in every production constructor.

There is one engineering defect in the recommended A3 wiring (F1, major). There are also two test/coverage
mismatches (F2, F3, minor).

## Lens checks (verified)

1. **Acceptance only through certify_closure.**
   - `closure_leaves` runs `closure_precheck` first (src/leaves.rs:955), which refuses a decomposable L_p and a
     wrong top. It builds the family, commits, then runs `certify_closure(ar, lam, r0, lp, None, kept, want=c_diag(b_p), outs, ..)`
     (src/leaves.rs:971-977) and returns only on Ok.
   - `certify_closure` (src/homprov.rs:945-986) checks four things:
     - ends equal to (lp, R) or (R, lp) as ExprIds;
     - degree multiset equal to c_diag;
     - an n x n identity-coset matrix of full rank at 2 points;
     - `assert_genuine` on every map.
   - With L_p = B_p, Hom(B_p, R) is free of graded rank c_diag. The transition matrix to a true basis then has a
     determinant of degree 0, so independence plus the exact degree multiset gives a basis. Wrong ends, wrong
     degrees, dependent families and non-genuine maps are all rejected.
   - The internal rank tests in `step_family`/`build_leaves` (src/leaves.rs:753-779, 924-936) are diagnostics
     only: they can make a call fail, never accept one.
   - SIndec re-certifies every closure result, both after a fresh provider call and after a layer-A load
     (src/sindec.rs:1060-1080).
2. **Elem gate.**
   - The Elem path runs only when `n == 1 && self.use_elem` (src/leaves.rs:690). A node is added only if the key
     is already certified in the arena, or if `validate_elem` succeeds and certifies it (src/leaves.rs:692-705).
   - The key from `elem_pattern` (src/leaves.rs:499-526) equals the key that `smap::elem_spec` computes for the
     node (src/smap.rs:336-366):
     - Eps/EpsL use (la, lb);
     - Trace uses (la, b) and TraceL uses (lb, b).
     - (kind, key) fixes the window uniquely, so a validated key covers every use.
   - `validate_elem` (src/leaves.rs:532-577) requires Hom^d to be one-dimensional. It compares every block at
     every coset of the window support at 2 points, with one common nonzero ratio, and rejects a zero entry of the
     exact map paired with a nonzero closed-form entry.
   - Even without validation, `assert_genuine` rejects an Elem whose (kind, key) is not in `consts`
     (src/smap.rs:863-875).
3. **Elem off in production.**
   - The default is `use_elem: false` (src/leaves.rs:454).
   - Every production constructor uses `Leaves::new()`: `provider_for_stamped_m4c` (src/leaves.rs:1029),
     `libbuild --provider leaves` (src/bin/libbuild.rs:321) and `m4c leaves` (src/bin/m4c.rs). None has a flag to
     turn Elem on.
   - Only tests set `use_elem = true` (tests/leaves.rs:158, 328).
   - In src, `certify_elem` is called only by the validated path (src/leaves.rs:694).
4. **Layer A with Elem constants.**
   - Constants are not serialized: smapio imports only the Elem node (src/smapio.rs:396), and `consts` is filled
     only by `certify_elem`. A loaded Elem node with no validated key in the session fails `assert_genuine`.
   - Sindec then returns `DecErr::Check` (src/sindec.rs:1080). That is a hard failure, not an acceptance.
   - The node's scalar comes from the file but is checked to be nonzero (src/smap.rs:692). A nonzero multiple of a
     validated closed form is genuine.
   - The stamp changed when Elem was turned off ("elem-window-whole" became "window-whole" in 8e70c0c), so layer A
     built with Elem on lives under another header.
5. **ClosureLeaves wiring.**
   - `full`, `pairs` and `end0` go straight to the inner provider, with the same `failed_stats` propagation
     (src/leaves.rs:1049-1063).
   - A closure with e_p goes to the inner provider only (src/leaves.rs:1066-1070).
   - The only indirect effect is shared state: Leaves' exact templates enter the shared TmplLib and become cache
     hits that Search always admits (src/homprov.rs:1752). Pair routes can therefore differ, but they stay certified.
   - On a `certify_closure` failure, `closure_leaves` does not drop the committed nodes (src/leaves.rs:973-977).
     The arena grows, and nothing else changes.
6. **Budget and admission.**
   - `admit`, `step_score` and `plan_words` only decide which exact calls are made and which junction words are
     used (src/leaves.rs:589-618, 798-872).
   - Every planned word is rechecked with `reduced_for` before use (src/leaves.rs:898).
   - A refused call makes the step fail (Refused). Nothing is accepted on the strength of a budget decision.

## Findings

### F1 (major): Leaves ignores the call Budget; in the recommended `m4c-lf` order it front-loads ~1000 s per side of doomed exact work before the chain

**Where:**
- `impl HomProvider for Leaves::closure` takes `_b: &Budget` and ignores it (src/leaves.rs:1000).
- `ClosureLeaves::closure` calls `closure_leaves` without `b` (src/leaves.rs:1075).
- `Leaves::admit` admits any call in (1e5, 4e5] that has a *finite* exact-shape record (src/leaves.rs:606-613). It
  never compares the predicted time with `b.t_max_s` or with the chain's per-call ceiling (`EXACT_L3_MAX_CALL_S` =
  600 s, the libbuild watchdog; src/sindec.rs:432, 471).
- In `step_family`, admitted degrees run in ascending order before a later degree is refused; on the Whole kind the
  error is returned only at that later degree (src/leaves.rs:722-748).

**Scenario:**
- Use `--provider m4c-lf`, which the report recommends (A2_report §5), with the a2 cost table, which the report
  says should be merged. The records are docs/colored_design/m4c/a2/costs_333_m4c_a2.txt:363-366: 139 371 u at
  474-494 s and 213 148 u at 440-500 s, up to 1.96 GB.
- Take type 0,1,2/0,2,1/3,0,0, which M5c's chain certifies. Its whole-word step Hom^d(E0.3 E1.3, F0.2 E1.3 E0.5)
  has no common prefix or suffix, so only the Whole kind runs.
- Leaves computes d=13 and d=15 (about 950-1000 s), then refuses d=17 (319 734 u, no record) and fails. This
  repeats for the other side (`closure_gens` does ins and outs in one job).
- Under the libbuild default `--wall 600 --rss 2048` (src/bin/libbuild.rs:6, 447), the job is killed inside Leaves.
  The retry at 2x reruns the same deterministic Leaves work first, so it is likely killed again.
- Result: a closure job that the production chain alone certifies becomes Time/Memory in the m4c-lf rebuild (new
  stamp, so no M5c layer A is reused). A single admitted record of more than 600 s would be admitted the same way.
- No wrong result is produced; the job fails.

**Fix:**
- Thread `b` into `closure_leaves`/`admit`. Admit a record only when `t * slowdown <= min(b.t_max_s, EXACT_L3_MAX_CALL_S)`.
- Keep a running wall-time total per closure call against `b.t_max_s`.
- Alternatively, pre-check every degree's admission for a step before making any exact call, so that a step known
  to be refused at d=17 never spends time on d=13/15.
- Or keep `m4c` (Leaves last) until this is in place.

### F2 (minor): `m4c_chain_closure_222` does not test the wiring paths its doc claims

**Where:** tests/leaves.rs:405-428.
- The doc says "non-closure calls go to the inner chain", but the test makes no `pairs`, `full` or `end0` call.
- For lf=false, λ=(2,2,2) has max label 2, so the inner provider is `Chain::for_labels(2, Exact)` (Exact first,
  src/sindec.rs:462-464). That chain certifies the (2,2,2) closure types, so the Leaves-last fallback is likely
  never reached. This covers the second loop iteration and the error merging and `failed_stats` accumulation
  (src/leaves.rs:1073-1096).

**Scenario:** a regression that sends `pairs` to Leaves (always Refused), or one that breaks the fallback after
an inner failure, would pass this test.

**Fix:**
- Add one pair call through the m4c provider and assert it matches the bare chain.
- Add a fallback case with an inner provider that refuses, such as NoProvider, wrapped in `ClosureLeaves` with
  `leaves_first=false`.

### F3 (minor): the V0 span gate and G1(a) closure gate run only with `use_elem = true`, not the production configuration

**Where:** tests/leaves.rs:157-158 (V0) and 327-328 (G1(a)). Production is `use_elem=false` (src/leaves.rs:454).
- Per the report's G1(a) row, 26 of the (2,2,2) steps are Elem. In production those steps take the Window/Whole
  path, which these gates never exercise for those steps.
- The Elem-off configuration is covered only by Ex 7.39/7.40 (no Hom-dimension comparison), and by the (2,2,2)
  wiring test in leaves-first mode. (3,3) is not covered.
- Likewise, the F31 §4 regression ran types 0-53 with the Elem-on binary 77329f3 (A2_report §4).

**Scenario:** a defect in the Window/Whole step selection that appears only on those steps would show up only as
production closure failures. They would be refused or exhausted, never accepted, thanks to `certify_closure`.
These gates would not show it.

**Fix:** run `v0` and `g1a_closure` for both `use_elem` values. At least run G1(a) (3,3) and (2,2,2) with Elem off.

## Notes (not defects)

- The `validate_elem_small` "wrong window" case (tests/leaves.rs:397) is refused by `Arena::add`/`check_window`.
  It never reaches the proportionality comparison, so that comparison has no negative test. The logic was checked
  by reading it and is correct. Elem is off in production.
- `ClosureLeaves` merges Leaves' Exhausted misses with the inner chain's misses. A chain Time/Refused failure
  followed by a Leaves Exhausted failure is therefore classified Exhausted. libbuild retries Exhausted too
  (src/bin/libbuild.rs:618), so only reporting changes.
- Leaves uses the unscaled cost table, while the chain uses `scaled_costs(.., slowdown)`. This only matters for F1.
