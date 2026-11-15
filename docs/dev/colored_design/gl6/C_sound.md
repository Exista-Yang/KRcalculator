# C_sound: soundness critique of SPEC_3C_GL6 rev 1 (and R1-R4)

Date 2026-11-11. Read-only review. Nothing was compiled or run. Code checked on dev f2d8be9 (which now contains
the M4b merge d49ee31; same code as colored-m4b 8ad4ac6 plus the review docs). Lens: can a step of the plan
report a wrong Poincare polynomial as certified?

Context from dev since the SPEC was written (PLAN.md, last section): run H gave the exact 3-colored trefoil at
gl_3..gl_6 (ranks 1, 10, 46, 136, bigraded equal to Wang, both primes). T(2,5) and larger ran out of memory at
crossing 5 (13.6 GB). This changes the schedule, not the soundness analysis below, except where noted (m5).

## 0. Verdict

**No critical flaw.** The mathematical pillars of the plan hold. I checked each against the code:

- **2-strand path (xengine).** It is exact end to end.
  - The decompositions are certified exactly: Π∘Ι = id is checked exactly (xengine.rs:654-661), Σ ranks = rank X
    (xengine.rs:635-637), and every map is a genuine bimodule map from `hom_space`, whose `kernel_lazy` stops only
    after verifying all equations.
  - d∘d = 0 is checked exactly on the unpruned scan.
  - The closure is exact.
  - The rung-by-rung transport identity of M2S is valid, because `decompose_word` composes without
    re-orthogonalizing: ι = (ι_1 ⊗ id)ι_2 and π = π_2(π_1 ⊗ id) (xengine.rs:624-630). The orthogonality check runs
    on the composites.
- **Pruning (R2 §3, D4).**
  - The "kept family suffices" argument is correct.
  - For (λ,λ) types the categorified kill rule is not an extra assumption. It follows from three facts:
    - tr_N(b_q·y) = 0 for all y; this is the unit test `kernel_of_gl_n_quotient`, which covers the (3,3,3) end
      block for N = 3..6 (schur.rs:2217-2227);
    - grdim Φ_N(B) = τ_N([B]), the CL1 closure contract;
    - positivity.
  - Kills in the intermediate blocks (4,2,3), (5,1,3), ... are not unit-tested. They are protected by the
    composite kept-multiset assertion in `dec_word`. If an intermediate summand is wrongly killed and its R_2
    product contains a kept type, the assertion fails.
  - So the weak point is the implementation, not assumption A1 (see M4).
- **P = G_y·F (D5).** It is exactly the product that cclose's rank stage already forms (`pab`, cclose.rs:1336-1362).
  P can replace F.
- **Box extraction (R3 §3).** It is exact for a homogeneous polynomial of the declared weight:
  - the vanishing of factorial Schur functions at z_μ, with nodes a_{μ_j+4−j}, is correct;
  - S is triangular for containment, so S^{-1} is too;
  - top coefficients are intrinsic;
  - nesting N' < 6 and reduced mode (ρ^0 = ∅, points with μ^0 = ∅) are correct;
  - D6 (Gram weight 40 > 36) is correctly settled: `validate` uses 2·weight = ins + outs (cclose.rs:452-456).
- **Freivalds (D7).** D_{h+1}Π^{-1}D_h is the right composite in rep/dual coordinates. It is equivalent to
  `closure_kills_map` on every component, and the p^{-3} bound is right.

What follows are gaps where an implementation that follows the SPEC literally can return a wrong polynomial
that passes the planned checks. Several of them pass two-prime agreement as well, because they are systematic.

## 1. Major

### M1. Nothing certifies that the declared degrees equal the actual degrees in the localized path. Box extraction turns an overestimated weight into silent zeros
- **Where.** R3 §3 and §5.2; SPEC D5 and MBOX `BoxAcc::new(lam, n, w, ..)`; the risk row "polynomiality of
  high-weight entries" (SPEC §5).
  - The weight of a P entry is computed from bookkeeping: w = (outs_y[b] + ins_x[a] + t_x − t_y)/2.
  - The object shifts t come from the summand records s, which are copied from `term_mult` and checked only in
    K_0.
  - The map degrees come from `ar.info`. The M4b review (`docs/dev/colored_design/review_m4b/certify.md`, F3) shows
    that `certify_pairs` does not check ends or degrees. It also says: "the form check cannot catch a degree
    error", and a cross-Arena TmplId reuse can return ι/π of degree s' ≠ s labelled s.
  - In xengine, `check_degrees` catches this exactly. In sengine no value is ever compared with the bookkeeping.
- **Failure.** Suppose an entry's true weight w' is smaller than the declared weight w.
  - Then c'_ρ = 0 for |ρ| = w and also for |ρ| = w + 1, so the residual check passes.
  - A second node sequence also gives 0, so the agreement check passes.
  - Both primes give the same zero.
  - Freivalds may still pass.
  - Result: the closed differential loses that component and the homology comes out too large. χ cannot see
    this, since it does not depend on D.
  - The same bookkeeping error also shifts objects in q.
  - SPEC §5's mitigation (spot re-extraction on a second node sequence) does not detect this failure.
- **Fix (cheap and exact in practice).** Use homogeneity under scaling: z_μ(λa) = λ·z_μ(a), so every identity-coset
  value satisfies P(λz) = λ^w P(z) and G(λz) = λ^{w_G} G(z).
  - Per word and prime, scan two extra points λ·z for a random λ of large multiplicative order. Assert the scaling
    law for every final edge entry, Gram entry and generator entry against its declared weight. An entry that is 0
    at both points gives no information; report how many such entries there are.
  - In M5, apply the same test per target: the entries of ι and π on the top-stratum coset must scale with the
    degree of `ar.info`. This closes review F3 inside the DecHom adapter.
  - Make it a gate in MBOX test (b) and in M6a-1. Add a negative test: declare one entry at w + 1, and the check
    must fail.

### M2. Streaming box accumulators are neither idempotent nor self-checking across checkpoints, shards and node redraws
- **Where.** SPEC MBOX `BoxAcc` and M6a-2 `box_run`:
  - "P values stream into BoxAcc checkpoints (accumulators + done-points), so runs are resumable and shardable";
  - "asserts the generator determinant at every point (redraws the nodes on failure)";
  - R3 §3: "Partial accumulators are summed".
- **Failure.** Any of these silently corrupts every c_ρ with ρ ⊇ μ:
  - a crash between the add and the done-set write, which re-adds point μ on resume;
  - overlapping shard ranges that are summed;
  - redrawn nodes after point 5 000 of 8 000, which mixes two node sequences;
  - shards from a different prime, node seed, library or prune_n that are summed.

  Detection is weak:
  - Residuals detect this only for entries with w < top. Unreduced P entries of weight ≥ 27, and reduced ones of
    weight ≥ 18, have no residual.
  - Two-run agreement is not a safe detector. Systematic corruption, such as an off-by-one shard boundary in both
    runs, pushes both ranks toward the same generic (too large) value, so the runs can agree on a wrong
    polynomial.
- **Fix.**
  - A checkpoint is written atomically: accumulators, done-point bitmap, and a header (prime, node-seed hash,
    library stamp, word list hash, prune_n, mode, N). Merging refuses mismatched headers and overlapping
    bitmaps.
  - `finish()` asserts that the received points contain the down-set {μ ⊆ ρ} of every output ρ and of every
    residual ρ. A run must never output ρ from an incomplete down-set, for example L(w_max) computed from Gram
    weights only, when P needs L(top).
  - Before any scan, run the generator determinant at all 8 000 points in a generator-only pass, which is cheap.
    A redraw after scanning has started invalidates the whole node sequence.
  - Add MBOX tests: double-add one point, shuffle the order, and split into shards.

### M3. M2S disk caches lack prune_n, ladder and algorithm stamps, and `dloc` is indexed by summand position
- **Where.** SPEC M2S part A item 2:
  - `wdecomp` is keyed by (table of p, rung word);
  - `dloc` is keyed by (table of p, sign, j);
  - the header is `krhom-xcache v1 field=.. labels=..`.

  `dloc` stores (si, ti, PMat) as positions in `decompose(p, ±, k)` (xengine.rs:833-860). `tensor_crossing`
  places them through `place[&(o, ks)][si]` (xengine.rs:908).
- **Failure.** In Gate A, pruned scans (`--prune 3/4/5`) and unpruned scans share `--cache DIR`.
  - A pruned summand list is the unpruned list with the killed summands removed. Its π's are normalized among
    the kept summands only.
  - An unpruned dloc loaded by a pruned run therefore attaches components to the wrong summands whenever the
    types at the shifted positions coincide.
  - In pruned mode d∘d is not checked exactly, only that the final components close to 0. Rickard and wdecomp
    collisions fail closed through the term and multiplicity asserts (xengine.rs:534, 543, 596, 653); dloc
    collisions do not.
  - R1 §3.2 also proposes switching `ladder()` to `ladder_minflow`. That changes the fibre basis under a cache
    keyed by table, so every cached PMat becomes stale.
  - R3 §6(a) caches closure state spaces and Π^{-1} per (type, N', mode), with no prime or library stamp.
- **Fix.**
  - Key and stamp every cache by: field, prune_n (None or N), the canonical-ladder rung word of p and of every
    summand q, `rungwise`, and a code stamp `XCACHE_ALG`.
  - Store the summand list (q, s) next to each dloc and assert it on load.
  - On load, recheck Π∘Ι = id exactly (cheap at ≤ 900).
  - Stamp state-space caches with the Gram-library hash and the prime. On load, recheck that Π is invertible.

### M4. The pruned 3-strand implementation has no pruned = unpruned categorified test, and the oracles planned for it cover few types
- **Where.**
  - R4 §9: "pruned vs unpruned | 3, 4 (2 strands) | 3 strands deferred".
  - The M6a-1 tests are 2-strand except "(1,1,1) words: object counts against engine2".
  - The duality check at N = 4 exercises only the 6 prune-4 types.
  - The Markov words visit only the p⊗id_3 corner (R4 R9).

  All of the new code (killed drops in dec_word and rickard_comp, `restrict`, transport, Freivalds, box closure
  of pruned data) is first run on its own 3-strand target at prune 6.
- **Failure.** A coset-block indexing or sign bug in the 3-strand transport that occurs only on general (3,3,3)
  types gives a complex that satisfies K_0, LB and D² = 0 (on a consistent but wrong complex). For 4_1 at N = 6,
  the only further checks are amphichirality and nesting.
- **Fix (cheap gates before M6a-3).**
  - (a) sengine on (1,1,1) at prune N = 2 and 3, for all 3- and 4-strand table knots, against engine2's gl_2/gl_3
    homology (validated on tables). At N = 2, B_{w0} is killed (cell (3)), so this is a real categorified test of
    dropping, Freivalds and box closure (a = 1 blocks).
  - (b) (2,2,2) pruned at N = 3 and 4 against unpruned (2,2,2) on 4_1, 5_2 and 6_2. This is feasible: G1(a)
    certified all 400 (2,2,2) targets, and every (2,2,2) canonical ladder has character b_p (review F1). (2,2,2)
    is the only λ where a 3-strand pruned = unpruned comparison is possible.
  - Report, per oracle, which kept (3,3,3) types and library pairs it touched.

### M5. The pruned closure API does not enforce valid_n = min(prune_n, box N) or D² = 0 before ranks are read
- **Where.** MBOX `ClosureDataC.valid_n` ("close refuses N' > valid_n") and M6a-2 `box_finish` ("valid_n, then
  `close_exact_colored` staged per type ..., then `closed_d2_freivalds`").
  - The homology formula dim − r_out − r_in (cclose.rs:1395-1399) is meaningful only if D² = 0.
  - xengine enforces both conditions inside `close` (xengine.rs:1276-1290) and refuses pruned scans without
    check_d2 (xengine.rs:1140-1144).
- **Failure.** M6a-2 (c) runs prune-4 scans "on nested points". If the driver reuses the N = 6 node set, it gets
  3×3-box Schur coordinates from a prune-4 scan, and box_finish sets valid_n = 6. A later `close(.., 5 or 6)` is
  then silently wrong: types killed at 4 but alive at 5 and 6 were dropped. Separately, a reporting path that
  prints the polynomial before or without the Freivalds verdict publishes an uncertified number.
- **Fix.**
  - `ClosureDataC` carries `prune_n` and `valid_n = min(prune_n, N_box)`.
  - `close_*` refuses N' > valid_n.
  - For pruned data, `close_*` also refuses unless a D² certificate (Freivalds, or P2) is recorded for that
    (N', mode).
  - Add a test: a prune-4 scan closed at N' = 5 must return Err.

## 2. Minor

- **m1. Two primes do not certify Q.** The computed complex is the reduction mod p of a complex defined over Z,
  so dim_{F_p} H ≥ dim_Q H in every bidegree. Agreement at 2^31−1 and 4294967291 makes torsion at both primes very
  unlikely, but it is not a proof. Wang (over Q) confirms T(2,n). The SPEC (§4) and the tests/data headers should
  state "over F_p for p = 2^31−1 and 4294967291; equal to Q unless both primes are torsion primes", and not "over Q".
- **m2. The rank cap from T0 (R3 §6(b), M6a-2 `box_finish`) is circular unless T0 itself runs uncapped.** If a rank
  is wrongly too large (corrupted Gram data), stopping rep selection at r_d hides it. Run T0 uncapped once per
  (type, N', mode) on the exact library file that production uses (hash in the header). After that, the cap is an
  assertion only.
- **m3. Reduced mode has no independent reference at N = 6 (only χ).** Add base-point independence for knots:
  - in xengine, `close_exact_colored(data, N, Some(1))` must equal the result with Some(0) (it costs seconds);
  - on 3 strands, mark block 1 and block 2. For the box this uses a disjoint point subset (μ^1 = ∅), so it is
    also an independent check of the extraction.
- **m4. Zero-dropping and Markowitz costs at box points.** SPEC_M4_M6a §6.3 says "zero results (every block zero)
  are dropped".
  - A small positive-degree component, for example on the rank-1 identity type, vanishes at roughly 1e6-1e8
    (component, point) pairs times deg/p. That is about 0.01-1 events per unreduced run.
  - The pivot sequence then changes, and the trace check fails. This fails closed, but it kills a multi-session
    run, and with the "redraw" policy of M2 it compounds the problem.
  - Fix: freeze the sparsity pattern from point 0, and never decide structure from values at box points.
- **m5. Guards to preserve.**
  - `close_stage(m)` must reuse `close`'s pruned guards (N' ≤ prune_n; the stage-m final_d2 closes to 0) and
    `crossings[..m]` for global_q.
  - M2S production must keep `check_d2` on. It is the only exact guard against transport and cache bugs. A fix
    for the 13.6 GB out-of-memory error must not switch it off; if it has to go, replace it by Freivalds on the
    closed complex.
- **m6. M5a adapter guards (M4b review F1/F2/F3).** D11 freezes homprov.rs at a version without these fixes, so
  `DecHom` must:
  - assert expand(ladder_char(word(lq))) == [(q,1)] with the same ladder function that `idem_needed` uses;
  - assert that lq is the ladder of q;
  - assert mult == [v^s] expand(char X)[q];
  - check the ends and degree of every map.

  Also, the `replay` of F32b must re-run M5's per-target certification, not trust F31's tier choices.
- **m7. Coverage and failing closed.**
  - The rotations of the Markov words, for example [2,1,1,1], need σ_1 pairs outside the Markov corner (D8).
  - A frozen DecData must return Err on a missing pair, never skip it.
- **m8. Conjectural oracles.** "Duality at N = 4 (expected)" and the Lewark-type (1+t)Q ≥ 0 relation can block a
  result (they fail closed) but cannot certify one. Label results checked only by them accordingly. Optionally,
  compare 4_1 with the conjectured Λ^3 superpolynomial (Gukov-Stošić) as a further cross-check.
- **m9. Risk weighting.** SPEC §5 lists "kill rule on 3 strands (assumption A1)" as the pruning risk. As argued in
  §0, the rule for (λ,λ) types is implied by tested decategorified data and the CL1 contract, and intermediate-block
  kills are protected by the composite multiset assertion. Move the weight of this risk to M4 (implementation)
  and keep the composite kept-multiset assertion mandatory in M5 (never debug-only).

## 3. Checked and fine (no action)
- Neumann among kept summands. Any subfamily of a unitriangular family stays unitriangular.
- `restrict(n')` keeps Π∘Ι = id.
- Stupid truncation of the Rickard terms equals the cell rule (asserted, xengine.rs:543).
- Reduced and box kill compatibility: setting e_0 = 0 is a ring homomorphism, so killed types stay killed.
- Identity-coset composition for P and P2: bimodule maps are block-diagonal per coset, and R^λ is supported on
  the identity coset.
- The closure-generator determinant criterion (review E). Each box point is uniformly random in the nodes, so the
  determinant vanishes with probability ≤ deg/p.
- Koszul sign (−1)^h, as in xengine.rs:906.
- The residual check catches an *underestimated* weight whenever w < top.
