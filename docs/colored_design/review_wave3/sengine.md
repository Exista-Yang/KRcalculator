# Review wave 3: localized engine (revm6a: src/sengine.rs, src/swire.rs, src/bin/sbox.rs)

Branch colored-m6a c95517f, read-only. Nothing was built or run. One check took under 1 s: the
multiplicative order of the scaling λ = 0x12345677 is above 5000 for both primes.

## Verdict

I found no defect in the engine core that can produce a wrong Poincaré polynomial. Gaussian
elimination, streaming eligibility, the Markowitz heap, frozen-trace replay, the Kron transport of
d_C ⊗ id, the Rickard sign, and the closure entries (values, weights, ordering) are all correct as
written; the checks are listed in "Verified correct" below.

The 6 findings are all minor:
- 3 are gaps in the driver or the guards (F1-F3).
- 2 are latent scale or merge problems for M6a-2a (F4-F5).
- 1 is a test claim with no test behind it (F6).

The highest-priority one is F1: the D14 scaling check is the only guard against the C_sound M1
failure (silent box zeros that both primes and both node sequences agree on), and it is opt-in.

## Findings

### F1 (minor). The D14 scaling check is opt-in, and box_word does not enforce it
- **Where.** src/bin/sbox.rs:80 (`if flag(a, "--scaling")`); src/swire.rs:730-829 (`box_word` never
  calls `scaling_check`, src/swire.rs:831).
- **Detail.**
  - SPEC D14 and C_sound M1 make the scaling law a required check per word and prime. It is the
    only check that sees an overestimated declared weight: the box extraction turns that into
    zeros, the residual at w + 1 is 0, and χ cannot see it. Both primes and both node sequences
    agree on that error.
  - sbox runs the check only with `--scaling`. Without the flag it still prints the closures and
    writes JSON with `"scaling": null`, Freivalds block counts and `chi_ok: true`, so the output
    looks fully certified.
  - Every recorded run passed `--scaling`. But any new caller of `box_word` (M6a-2a, or a
    production script that drops the flag) gets uncertified results with no signal.
  - Side note: `check_scaling_all` counts entries declared `None` (negative or odd weight, which
    must be 0) in `zero_both`. The reported "zero at both points" figure therefore mixes expected
    zeros with entries that carry no information.
- **Failure scenario.** A library ι/π is labelled with degree s but has true degree s − 2 (the
  M4b review F3 class), and the run omits `--scaling`. The affected P entries are declared one
  weight too high and come out of the box as 0. The closed differential loses those components,
  the homology is too large, and both primes and both node sequences report the same wrong
  polynomial.
- **Fix.** Run `scaling_check` inside `box_word`, or refuse to emit results unless it passed, and
  record its outcome in BoxWordOut / JSON. Report `None`-declared entries separately from
  zero-at-both.

### F2 (minor). A χ ≠ decat_chi mismatch is not fatal, and the compare tools ignore it
- **Where.** src/bin/sbox.rs:94-104 (prints `CHI FAIL`, stores `ok = false`); results are still
  printed (:106-110) and the process exits 0. `chi_ok` goes into the JSON (:124), but
  tools/g2b_compare.py and tools/stretch_compare.py never read it.
- **Failure scenario.** A (3,3,3) word with no m2s/engine2 reference has a wrong closure whose
  Euler characteristic is off. sbox prints the polynomial and exits 0; the only signal is one
  stderr line among many. A batch pipeline that collects the JSON results certifies it.
- **Fix.** Exit nonzero (or panic) on any χ mismatch, and make the compare tools assert `chi_ok`.

### F3 (minor, hardening). Replay checks that each recorded pivot is nonzero, but not that it equals the recorded value
- **Where.**
  - src/sengine.rs:1119-1135: replay applies the recorded `(a, b)` and requires only `scalar_of`
    = Some(c ≠ 0).
  - `Round` (src/sengine.rs:566) stores only the pivot pairs, not the values.
- **Detail.**
  - A degree-0 component between equal types is c·id with c ∈ k, the same constant at every
    point, because End⁰(B_p) = k and every component is built from certified bimodule maps.
  - A replay pivot with c(z) ≠ c(z₀) is therefore proof that the construction at z is not a
    genuine bimodule map (Kron transport, π K ι, a library block). Today it goes undetected
    until (and unless) the residual, scaling or Freivalds checks happen to see it downstream.
  - Storing one field element per pivot (≈ 8 300 pivots at σ^7) and comparing on replay costs
    nothing. It would turn this case into an exact `PointDependent` refusal at every box point,
    and at λz.
- **Failure scenario.** A point-dependent defect makes a pivot scalar a nonconstant degree-0
  rational function, for example x₁/x₄. Elimination divides by it, and the entries become
  rational in z. Detection then depends on the box residual seeing it. A degree-0 ratio also
  passes the scaling law, so the scaling check cannot catch this one.

### F4 (minor, latent for M6a-2a). The BoxAcc header does not bind the library content or the frozen trace
- **Where.**
  - src/swire.rs:238-247: `WLib.hash` covers only SWIRE_STAMP, the prime, λ, prune, type names,
    closure-generator degree lists and the words.
  - src/swire.rs:745: that hash becomes `BoxAcc::new(.., lib.hash, ..)`, i.e. the merge guard
    `BoxHeader.library_hash` (src/boxclose.rs:602).
  - BoxHeader has no hash of the record trace (pivot sequence), and `entries_hash` covers only the
    weights.
- **Detail.**
  - The closure entries are coordinates in a basis fixed by the library maps (ι/π, Rickard
    components, closure generators) and by the pivot sequence chosen at the record point.
  - Two libraries built with different provider stamps or cost tables (the host-stamp problem in
    the M6a1b report §1, §7.2) can have identical names and degrees but different bases.
  - Likewise, `--lazy-markowitz` versus the default heap gives a different pivot sequence that can
    still produce the same skeleton.
  - Shards from either pair pass `BoxAcc::merge` and get summed. This is C_sound M2's "shards from
    a different … library", and it is silent.
  - `box_word` does not save or merge today, so this is latent until box_run/box_finish.
- **Fix.** Fold the provider stamp and a digest of the layer-A keys and maps used (or DecData's own
  content hash) into `WLib.hash`, and add the frozen-trace hash (e.g. `FrozenTrace.stages`) to
  BoxHeader.

### F5 (minor, scale). box_word keeps every point's full entry vector in memory before accumulating
- **Where.** src/swire.rs:764-778: `par_map` returns `Vec<Result<(id, Vec<Vec<F>>, ..)>>` for all
  points, and only then does the loop at :782 feed BoxAcc.
- **Detail.**
  - Peak memory is (#points) × (#entries over all closed stages) field elements, plus the
    per-thread scans. That is outside `mem_point_bytes` and outside every budget counter.
  - G2B (7.5k entries × 399 points) and σ^7 reduced (34k × 20) are small.
  - (3,3,3) at N = 6 has 20³ = 8 000 unreduced points. P entries are ins_s × outs_t per edge with
    large generator counts, so the vector reaches tens of GB and the process dies with OOM instead
    of being refused.
- **Fix.** Accumulate each point as it finishes (a Mutex<BoxAcc> per stage and mode, or a
  channel), or batch the points.

### F6 (minor, test claim). tests/swire.rs advertises a negative scaling test that does not exist
- **Where.** tests/swire.rs:8: "the scaling check (D14) and its negative test". The only scaling
  use in the file is the positive call at :103, which asserts `checked > 0`.
- **Detail.**
  - The negative test that does exist (tests/boxclose.rs:480) runs `check_scaling_all` on exact
    R^λ data.
  - Nothing tests that the wired path catches a mis-declared weight end to end: scaled node
    sequence → `box_x(&sn, ..)` → `scan_point` replay → `check_scaling_all`. This path is the D14
    gate for M6a-1.
- **Fix.** Add the negative test: run `scaling_check` with one declared weight bumped to w + 1
  (or one generator degree off by 2) on (2,2) σ^3, and require an Err.

## Verified correct (file:line, for the lens questions)

**Gaussian elimination at a point.**
- `apply_pivot_inner` (src/sengine.rs:1313-1373) computes d_xy −= d_ay·c⁻¹·d_xb for
  x ∈ in(b)∖a and y ∈ out(a)∖b, blockwise per coset. The call is
  `compose_acc(dst, g = d_ay, f = d_xb, α = −c⁻¹)` with
  `gemm_acc(c(m×n) += α·a(m×k)·b(k×n))` (src/smap.rs:65). Dimensions: m = dim y, k = dim a = dim b,
  n = dim x.
- in(a) and out(b) are dropped. This is the standard lemma, with the right sign.
- The pivot is checked to be c·I at every coset with c ≠ 0 (`scalar_of`, :508-533); a non-scalar
  candidate is a hard `SErr::Pivot` error, not a skip.
- `add_map` refuses h ≠ h + 1 and foreign layouts (:922-956). Zero sums remove the edge, so the
  pattern reflects values at the point.

**Streaming eligibility (streaming = full elimination).**
- swire inserts the term-k objects, the d_C ⊗ id arrows inside term k, and the id ⊗ d arrows from
  the previous step, all within step k (src/swire.rs:417-558). So in(b) is complete for every b in
  a closed step.
- out(a) is complete iff step(a) < the current step (`eligible`, src/sengine.rs:1086-1094).
- Fill-in never reaches step cur + 1, and the arrows inserted later (out of step cur, or inside
  step cur + 1) are never in the update set of an earlier pivot. So the streaming order equals full
  elimination in the same pivot order.
- Arrows out of dead targets and into dead sources are dropped. The other two kinds raise
  `SErr::Stream`, which is right.

**Markowitz heap vs rescan.**
- `touched` = in(a) ∪ in(b) ∪ out(a) ∪ out(b) covers every object whose |in| or |out| changes,
  including cancellation and fill-in endpoints. Versions are bumped before the re-push. The
  same-key candidate index is pruned by `retain` on touched objects.
- The heap therefore pops the same min (cost, a, b) as the rescan. Duplicate current entries for
  the same (a, b) are harmless.
- `Lazy` terminates (each entry can be re-pushed at most once per pivot) and is deterministic.

**Replay (D14).**
- These are all compared against point 0: the pattern hash after every step and every round
  (objects with ids, p, h, t, and the nonzero-component set), the per-crossing determinism hash,
  the round count and stage, and `finish` requires the trace to be consumed exactly
  (src/sengine.rs:829-848, 1010-1016, 1053-1057, 1119-1144).
- A recorded pivot that becomes 0 or non-scalar is refused as `PointDependent`. So is a component
  that vanishes or appears anywhere.
- Equal final patterns imply that the replayed complex is minimal: point 0 leaves no equal-(type, t)
  edge.
- Block-level zero patterns inside a component are value-level only and correctly do not matter.

**d_C ⊗ id via the Kron rule (src/swire.rs:340-366, 432-529).**
- `find` returns an index relative to the coset and `range` returns global ids
  (src/sexpr.rs:188-213). They are used consistently: colinfo is indexed by position in
  `sfi.range(c)`, which is the row order of ι's block, and colp and rowmap are relative.
- cp is the rank of the prefix's top configuration in the top coset space shared by pa and pb.
- K ι accumulates f[rp, colp] · ι[g, :] into row(rp, g). π · KI has π's cols = kr. The block lands
  in layout(σ.q, τ.q) with its size checked.
- ι/π ends and degrees are asserted in `build_wlib` (:213-225).

**id ⊗ d.**
- The sign is (−1)^{h of the C object}, Koszul.
- j = min(kp, k) matches `build_wlib`'s (ks, kt), which is (j, j+1) for a positive crossing and
  (j+1, j) for a negative one.
- Terms are inserted ascending for a positive crossing (hdeg = k) and descending for a negative one
  (hdeg = −k), so every term arrow runs from step s to s + 1. `kept_terms` matches the web label
  max(c, d) + k for all four rung patterns (src/schur.rs:1645-1665).

**K_0 / LB / d∘d actually run.**
- At the record point, `box_word` sets `k0_check: true` (src/swire.rs:734). `Decat` is compared per
  crossing (:568-582).
- sbox turns `check_d2` on for unpruned runs (src/bin/sbox.rs:63), and `end_crossing` executes it
  (src/sengine.rs:988-993).
- Replays share the object multiset through `same_structure`.
- The LB assert is implied by K_0 equality, so it adds no information, but it is not wrong.

**Closure entries.**
- G[b][a] = h_b · g_a and P[b][a] = h^t_b · d_cid · g^s_a at the identity coset
  (src/swire.rs:596-660). These match cclose's `Pairing` contract (src/cclose.rs:2365).
- Weights are (ia + ob)/2 and (t_s − t_t + ia + ob)/2 (src/boxclose.rs:1000-1025), consistent
  with the engine's e = t_x − t_y.
- Reduced coordinates equal flows, guaranteed because `lib_shape` refuses e_p ≠ None for every
  type in the walk.

**Same skeleton at every point.**
- `same_structure` and a check of `skel.objs` and edges run per point (src/swire.rs:768-773).
- `kinds` and `weights` are pure functions of the skeleton.
- `add_point` refuses duplicates, wrong lengths, and nonzero values for `None` entries.

**Box points and scaling.**
- `box_x` uses the same roots as `block_e_values`. Nodes are distinct across blocks
  (src/boxclose.rs:178-198).
- `scaling_check` scans z_μ and λz_μ under the same frozen trace. λ has order > 5000 for both
  primes.

**NoProvider.**
- Every call returns `MissKind::Refused`. sindec turns that into an engineering `DecErr`, and
  `build_wlib` either errors or records `no_closure`, which `extract_stage` turns into a hard
  error for final types.
- K_0 at the record point would also catch an incomplete decomposition.

**Two primes and two node sequences agreeing on a wrong answer.**
- The guards are: the K_0 = decategorified scan, d∘d (unpruned), Freivalds D² on the closed
  complex, the down-set and residual checks, χ, and scaling.
- The remaining gaps are F1 (scaling is opt-in), F2 (χ is non-fatal), F3 (pivot values are not
  compared), and F4 (once shards exist).
