# colored-m5 (M5a) soundness review: sindec / smapio

Branch colored-m5, worktree .claude/worktrees/revm5, HEAD bf326b7, diff ce91f0e..HEAD. Read-only; nothing
compiled or run (one 1-second Python check of the order of SCALE_LAMBDA). Lens: can a decomposition be
reported, cached or reloaded wrong without an Err or a panic?

## Verdict

No critical or major defect. The decomposition path is sound as written: every route to a stored
summand either goes through the full re-certification or is an exact composite of certified pieces, and
every failure is an Err or a release assert. There are 5 minor findings: one D14 coverage gap, one
unchecked API argument, one assertion that can never fail, one weak test, and one cache failure mode
that cannot recover by itself.

## Checked and correct (file:line)

- **Kill predicate (D3).**
  - `SIndec::killed`, `ty`, `killed_rank` and `split_mult` use `Schur::killed_at` only:
    sindec.rs:504-506, 519-523, 543-553, 556-569.
  - `DecData::killed_at` uses `kill_level` (sindec.rs:186-192, 1230-1232). This equals killed_at
    exactly: killed ⇔ max(margins, ν₁) > N (schur.rs:1436-1442).
  - Rickard truncation uses the label rule max(c,d)+k > N (sindec.rs:957-961, 1277). The cell-rule
    assert checks it against killed_at in every pruned builder (948-956).
- **Kept multiset and dropped summands.**
  - dec_one builds its targets from the kept part of `term_mult` (sindec.rs:615, 637-667), so the
    multiset is right by construction.
  - dec_composite has the release `assert_eq!(got, kept)` (sindec.rs:908-913) and rank bookkeeping
    (914-917). dec_one has rank bookkeeping at 669-673.
  - All-killed: release assert killed_rank = flows(X), no provider call (sindec.rs:618-622). The
    composite all-killed case is covered by the same rank assert with total = 0 (916-920).
- **Per-target certification, on fresh and loaded maps alike.**
  - Ends and degrees are read from `ar.info` (sindec.rs:334-345).
  - Form at rep(q) at ≥ 2 points, equal at both and invertible (645-654); `assert_genuine` (655-657).
  - π'_a = Σ_b (G⁻¹)_ab π_b gives π'_a ι_c = δ_ac (659-664). This is exact because L_q is checked
    indecomposable of type q (`ty`, 525-527), so End⁰(L_q) = k.
- **D14 law.** v(λz) = λ^u v(z) with 2u = d − sh(tgt) + sh(src).
  - This is the same unnormalized degree that `certify_added` uses (tmpl.rs:398).
  - Localized entries F = L_t Φ L_s⁻¹ are homogeneous of uniform degree.
  - For Local nodes, the window shift and the full-expression shift differ by the same prefix and
    suffix contribution on both ends, so the law holds there too.
  - Kron and Dot preserve the law.
  - Odd 2u with a nonzero entry is rejected (sindec.rs:242-251). Negative exponents are handled
    (218-225).
  - λ = 1234567 has no order ≤ 5000 mod 2³¹−1 or mod 4294967291 (checked).
  - The check is applied to every pairs target: fresh (378, 701) and layer-A load (686-692). On a
    load it runs before the form check, which then guarantees the entries at rep(q) are nonzero.
- **Neumann.**
  - Unitriangularity: any nonzero π_t∘ι_s with s_s ≤ s_t is an Err (sindec.rs:744-750). This also
    catches same-target off-diagonal leftovers.
  - The series is M = Σ(−N)^k, pw ← N·pw, with signs −,+,… (764-808). Π' = MΠ gives Π'Ι = M(I+N) = I.
    Degrees are consistent: −s_t for every Lin term (821-824).
  - Termination is bounded by n (783).
  - `commit` keeps only nodes reachable from the new π. No expression is created inside the mark
    region, so ι/X ids stay valid (smap.rs:894-996).
  - Π∘Ι = id is re-checked afterwards on fresh points over all cosets of X (838-859, called at 674-676).
  - Σιπ = e_X is checked whenever nothing is killed (860-875).
- **Composites.** Kron/Comp ends are checked against X (sindec.rs:889-903). The identity cases
  (None,None) and (Some,Some,true) give the right source by the literal word equality asserted at 628.
  Π∘Ι is re-checked (922).
- **Rickard.**
  - (ks,kt) = (j,j+1) / (j+1,j) matches the template direction (tmpl.rs:503-505, 523, 531-533).
  - The window sits at at = |L_p|, first = i (sindec.rs:966-967), and `check_window` validates the
    labels.
  - Err on a killed p or term (936, 957-961).
- **closure_gens.** `certify_closure` runs on fresh and loaded generators alike (sindec.rs:1044): ends,
  declared-degree multiset = c_diag, rank n at 2 points, genuine (homprov.rs:945-1009).
- **restrict / freeze.**
  - `restrict` returns Err for n > prune_n and allows any n on unpruned data (sindec.rs:1236-1241).
  - A filtered subfamily keeps Π∘Ι = id. The Rickard remap `(ms[a]?, mt[b]?)` is correct
    (1281-1283).
  - Every Tid used by `restrict` has a `level` entry (1150-1161).
  - Frozen lookups return Err(NotPrecomputed) and never return an empty result (1195-1209).
- **DecHom guards.**
  - `dechom_precheck` (sindec.rs:298-326) re-implements the F1/F2 checks plus the canonical-ladder
    check, using the same ladder function as `idem_needed`.
  - It does not call `homprov::pairs_precheck`, but nothing is bypassed: every backend calls
    pairs_precheck / closure_precheck itself (homprov.rs:1316, 1349, 2605, 2734), and Chain only
    delegates.
- **Layer A and smapio.**
  - The key holds λ, field (F_p), table strings, the canonical-ladder words of L_p and L_q, R, s and
    mult (sindec.rs:1077-1096). The header must match exactly, including the key (1110).
  - Templates are re-obtained by content and checked against their hash; added templates are
    re-certified (smapio.rs:321-383).
  - Nodes are rebuilt with `ar.add`, so degrees are recomputed and never read from the file.
  - The post-order export is correct for DAGs (smapio.rs:174-198).
  - Writes are atomic: temp file, then rename (sindec.rs:1134-1136).

## Findings (all minor)

### m1. Closure generators get no D14 scaling check (fresh or layer-A)
- **Where.** sindec.rs:1023-1046 (`closure_gens`) relies only on `certify_closure`, which checks the
  *declared* degree multiset (homprov.rs:967-974) and the identity-coset rank. `dechom_postcheck` /
  `check_scaling` run for pairs only (sindec.rs:330-360, 686-702). The DecHom closure adapter has no
  postcheck at all (384-397).
- **Failure.** A generator whose `MapInfo.deg` differs from its true homogeneous degree passes. This
  is the review-F3 class: a leaf template with a wrong `deg`. Such a generator is saved to layer A and
  reloaded through the same check. Its weight then enters the closure pairing with the wrong q-shift.
  Only the MBOX generator-entry scaling check, which is not built yet, would see it. This is the C_sound
  M1 "generator entry" case.
- **Fix.** In `closure_gens`, after `certify_closure`, call `check_scaling(ar, m, identity_coset,
  info.deg, z)` for every generator, on both the load and the fresh path, and count the zero entries.

### m2. Frozen `DecData::rickard_comp` ignores (c, d)
- **Where.** sindec.rs:1204-1206 discards `_c, _d` and returns position-0 components. The builder
  returns Err when (c,d) ≠ top(p)[0..2] (998-1004).
- **Failure.** On 2 strands with c ≠ d, a caller whose label bookkeeping is swapped gets Err from the
  builder but silently gets data from the frozen path. The data are the components for p's real top,
  while the caller keeps using its own (c,d), for example in `crossing_terms` shifts. This is
  harmless for the 3-colored target, where c = d.
- **Fix.** Repeat the check in `DecData`. The top labels are available from
  `types[p].lad` → `Expr::top`, as `restrict` already does at 1271-1274.

### m3. The Rickard component degree assertion can never fail
- **Where.** sindec.rs:984-990 asserts `info(m).deg == 1 + s_σ − s_τ`.
- **Why it is vacuous.** m = Comp(π_τ, Comp(φ, ι_σ)), and `Arena::add` defines the Comp degree as the
  sum of its parts (smap.rs:736-743). So the equality holds for any inputs. The M5a report lists it as
  a check.
- **Effect.** No soundness gap today: ι and π are D14-checked, and d_j's degree is asserted in
  tmpl.rs:503-506. The assertion should not be counted as Rickard degree coverage.
- **Fix.** Either drop it or replace it with `check_scaling` of each kept component at a coset where
  it is nonzero.

### m4. d7_restrict does not test the re-indexing it claims to test
- **Where.**
  - tests/sindec.rs:614-622 compares the Rickard lists of `restrict(n)` with a fresh prune-n build by
    `len()` only.
  - The summand lists are compared by raw Tid numbers taken from two independent `Schur` instances
    (518-527, `format!("{}", p)`).
  - No scan (D5) runs on restricted data.
- **Failure.** A bug in the (si, ti) remap at sindec.rs:1281-1283 that keeps counts unchanged would
  pass. Examples: swapped `ms`/`mt`, or indexing into the unfiltered list. The bug would first show up
  in M6a. The report's claim "restrict(3/4) of a prune-5 build = fresh builds" is checked only one way
  (fresh ⊆ restricted) and only by counts for Rickard.
- **Fix.**
  - Compare (table(q_σ), s_σ, table(q_τ), s_τ) per component.
  - Key by table strings, not Tids.
  - Run `mini_scan` on `restrict(3)` of the prune-5 data against xengine at N = 3.

### m5. A bad layer-A file is a permanent hard error, and files are saved before M5's own checks
- **Where.**
  - `target_pairs` saves a provider result (sindec.rs:707) *before* `dec_one` runs M5's form and
    genuineness re-certification (644-657).
  - A loaded file that fails the postcheck returns `DecErr::Cache` (688). It is not rejected and
    recomputed.
  - Only a header mismatch triggers a recompute (1110-1113). The module doc at 22-26 says otherwise.
- **Failure.** A result that passes the provider's certify_pairs but fails M5's check (a different
  point set) is persisted. So is a file written by a build whose checks were later tightened without a
  SINDEC_VERSION bump. Every later run, including M5b's resumable one-process-per-pair driver, then
  fails on that target until someone deletes the file by hand. This fails closed, so it cannot produce
  a wrong result.
- **Fix.** Save only after `dec_one`'s checks pass. On a load-side check failure, count it as
  rejected, delete or ignore the file, and recompute.

## Not defects, recorded for M5b/M6a
- Builder and frozen data report a killed item differently. The builder returns `Killed`;
  `DecData`/`restrict` return `NotPrecomputed` for a killed type, decomposition or Rickard term
  (1195-1203, 1277-1278). M6a must decide truncation from Schur and must treat NotPrecomputed as
  fatal. A blanket "Err ⇒ truncate" would drop missing data silently (C_sound m7).
- `DecData::reduced` panics ("index not built") for a type whose ladder never served as the end of a
  map, for example an intermediate type that is only an all-killed source (1217). This cannot happen
  for scan objects, because they are summands whose ι source is L_q.
- Rickard nonzero tests use a fresh `Eval` per (σ,τ) pair (sindec.rs:975). This is the cost hotspot the
  report flags for (3,3,3).
