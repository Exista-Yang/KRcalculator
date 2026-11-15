# M4b review: certification soundness of `src/homprov.rs`

Branch `colored-m4b` at 8ad4ac6, diff `3c3323b..HEAD`, worktree
`/Users/existair/Documents/KRcalculator/.claude/worktrees/m4b`. Read-only review: nothing was compiled or
run. Every claim below was checked against the code (file:line).

Question asked: can the provider report a result as certified when it is not a basis of the true Hom space
(wrong dimension, dependent maps, non-genuine maps, wrong degree, wrong stratum)?

## Verdict

| certificate | verdict |
|---|---|
| `certify_full_plan` / `certify_full` (Exact full, Search full, end0) | **Sound**, including the new orbit fallback (section A). |
| Genuineness of Search-built maps (window splits, extremes plus dots, commutations, rex moves, leaf isos, S4, S6) | **Sound** (section B). |
| Cost-table and window admission | **Cannot lead to a partial result being accepted** (section C). |
| Template cache reuse | **Keys are complete** for `TmplLib`. Search's own TmplId caches are not bound to an Arena (F3). |
| `certify_closure` | **Sound** when L_p is indecomposable. **Fails closed** otherwise (F4). |
| `certify_pairs` | **Sound only under End^0(L_q) = k** (or a correct e_q). The provider never checks this, and with `eq = None` and a decomposable canonical ladder it certifies wrong summand pairs (**F1, major**). |

Findings: 1 major, 5 minor. None affects the G1(a) result, because every canonical ladder of (3,3) and
(2,2,2) has character exactly b_p (SPEC §0.2 F1). F1 and F2 do affect G1(b), the (3,3,3) blocks of the
3-colored target, and any later caller that leaves out e_q.

---------------------------------------------------------------------------------------------------

## Findings

### F1 (major): pairs are certified without checking that L_q is indecomposable of type q

**Where.**
- `certify_pairs` (homprov.rs:750-770) checks only three things: the rank of the top-coset form at
  c = rep(q), that the form is the same at a second point, and genuineness.
- `form_at` (697-725) requires only that L_q has exactly one flow at c.
- Both Exact (1065-1077) and Search (`finish`, 2282-2285) end in this check.
- The g1 harness passes `eq = None` for every type: src/bin/g1.rs:342 (pairs) and :375 (closure).

**Why the check is insufficient.** The form proves a splitting only if π_j∘ι_i = G_ji·id_{L_q}, which needs
End^0(L_q) = k (or e_q End^0 e_q = k·e_q). Now take L_q = the minimal-flow canonical ladder of q with
char(L_q) = b_q + Σ m_r b_r, where every r lies below q:
- The extra B_r are not supported at rep(q), so L_q still has exactly one flow there. The `p.rows != 1`
  test (707) passes.
- Exact: the bases of Hom^s(L_q, X) and Hom^{-s}(X, L_q) give a form of rank = #B_q⟨s⟩ in X = mult, so the
  pairs are certified.
- Search S6 (`window_maps`) behaves the same way.
- Search S1/S3: leaves have character exactly b_q, so `lf.w != qw` and `leaf_iso` (2018-2071) runs. Its
  premise "both are indecomposable of the same type: any such map is an isomorphism" (2017, M4b_design.md:52)
  is false here. `pick` (2053-2063) still finds degree-0 maps with nonzero top scalar, so the form is
  nonsingular.

**Result.** Certified `Pairs` with π_j∘ι_i = G_ji·e, where e is the idempotent of B_q inside L_q. A consumer
that splits L_q⟨s⟩^mult off X with Σ ι G^{-1} π gets the wrong object: rank and character of L_q instead
of B_q. This is the "wrong stratum/object" failure.

**Same gap for a mismatched type.** Nothing checks that `lq` is the ladder of `q`. If `lq` has type q' ≠ q
and one flow at rep(q), pairs for L_{q'} are reported as type q.

**Realistic trigger.** SPEC §0.2 F1 says that for the minimal-flow ladder L_p, char(L_p) ≠ b_p for:
- (3,3,3): 1 of 55 types, 3 of 219 over the blocks a scan visits;
- (2,2,2,2): 6 of 282 types, 25 of 1 110 over the visited blocks.

G1(b) draws exactly these populations and calls with `None`. The 3-colored 3-strand gl₆ target lives in
(3,3,3).

**Fix (cheap).** In `pairs`, when `eq.is_none()`:
- require `s.expand(&s.ladder_char(&word(lq)))` to equal exactly `[(q, 1)]`. This is equivalent to
  `hom_gd(char lq, char lq, 0, 0) == [1]` together with the leading table q, i.e. End^0(L_q) = k by the
  Hom formula;
- otherwise return `Missing { kind: Refused, what: "L_q decomposable: e_q required" }`.

When `eq` is present, check at least that rank(e_q block at each coset) equals c_t(b_q)(1). Apply the same
check in `closure`: today it fails closed there, see F4. g1 should skip or flag such types instead of
counting them as certified.

### F2 (minor): `mult` is never checked against the character of X, so Search can certify a partial set

**Where.** `pairs_inner` computes `mq = expand(char X)[q]` (2279) but uses it only for s_lo/s_hi (2280).
`certify_pairs` accepts any family whose form rank equals the caller's `mult` (select_form, 728-733).

**How the backends differ.**
- Exact fails closed: its full bases give form rank = the true multiplicity, so `select_form` returns None
  whenever the two differ.
- Search accepts when its candidate set happens to reach exactly `mult`. Ways this happens:
  - pruned children: `flows_max`, `continue` in `tree` (1910);
  - stuck branches;
  - S2 uses only `io_lo.first()` / `pi_hi.first()` (2312).

**Realistic trigger.** G1(b) takes `mult` from `term_mult(p, [r])`, which is expand(b_p·[r]) (homprov.rs:69),
while X is the word L_p R. For a decomposable minimal-flow L_p (F1 populations), char(X) contains strictly
more B_q⟨s⟩, so a "certified" result is a strict subset of the summands. SPEC §4.4 says: "Partial bases
are never accepted".

**Fix.** Require `mult == [v^sdeg] mq` in `pairs_inner` and in `Exact::pairs`, or return Refused with the
two numbers.

### F3 (minor): pairs/closure do not check map ends or degrees; Search caches TmplIds without binding them to an Arena

**What is checked where.**
- `certify_full_plan` checks (src, tgt, deg) of every map (622-628).
- `certify_pairs` (750-770) checks none of them.
- `certify_closure` (772-797) checks the degree multiset but not src/tgt.
- The form check cannot catch a degree error: if ι has degree s' and π has degree −s', then π∘ι still has
  degree 0 and the form rank is unchanged.

**Degrees are correct by construction today:**
- tree shift `sh + pc.s` equals the template degrees: the split templates are certified homogeneous in
  `certify_added`, tmpl.rs:371-402;
- S2 dots add 2k;
- S6 is exact at sdeg.

**Why this is fragile.** That correctness depends on Search's per-instance caches `splits` / `comms`
(1370-1371). They store TmplIds of one `TmplLib`, and "one Search per Arena" exists only as a comment
(1351). If a `Search` or `Chain` is reused with a second Arena, which is plausible in M5/M6a with one Arena
per block or scan:
- `Local::add`'s `check_window` catches most stale ids;
- it does not catch a same-window template of another degree (Hom^s versus Hom^{s+2} on the same window
  pair);
- the leaf shift then differs from `MapInfo.deg`, and `certify_pairs` would return ι/π of degree s' ≠ sdeg
  as degree-sdeg pairs. Downstream this means a wrong q-shift.

**Two more gaps.**
- M4b_design.md:82 says `assert_genuine` "plus a check that every node's degree bookkeeping is consistent".
  `homprov::assert_genuine` (537-539) only forwards to `Arena::assert_genuine` (smap.rs:855-889), which
  checks leaf certification only.
- `detect_plan` (578) does not assert `s.lam == ar.lam`, unlike `detecting_cosets` (568). A mismatched
  `Schur` silently changes the dimension that `certify_full` compares against.

**Fix.**
- In `certify_pairs`, check `info(ι).src == lq`, `tgt == x`, `deg == sdeg`, and π: x → lq at −sdeg.
- In `certify_closure`, check src/tgt (R^λ ↔ lp).
- Assert the lambda in `detect_plan`.
- Tag `splits` / `comms` with an identity of the library (for example a nonce stored in `TmplLib`), or
  clear them when it changes.

### F4 (minor): closure with an idempotent e_p can never certify; S2 closure candidates leave the image of e_p

**The e_p path always fails.** `certify_closure` requires every identity-coset block to have exactly
n = c_diag(b_p)(1) entries (785-787), and `Selector::offer_many` skips any other candidate (1163-1165).
When e_p is needed, L_p is decomposable, and every extra summand B_r is supported at the identity coset
(c_diag(b_r)(1) ≥ 1, since the identity is the minimal double coset). So every candidate has more than n
entries. SPEC §4.3 item 7 asks for the identity-coset flows "composed with e_p", i.e. reduced coordinates
of im e_p; this is not implemented.

**S2 candidates are not re-wrapped.**
- `wrap(ep)` is applied to the S1/S6 and S4 candidates (2400-2405 and their uses).
- It is not applied to the S2 candidates Dot∘g (2433): for ins, dt∘(e_p∘g0) is not in e_p·Hom(R^λ, L_p),
  because a dot need not commute with e_p.
- Today this fails closed. It would become a false certification if reduced coordinates were added
  without re-wrapping, because Q·block projects the map and hides that it is not in Hom(R^λ, B_p).

**Impact.**
- The one decomposable (3,3,3) type of (λ, λ) cannot get closure data through this API.
- In G1(b) such a closure job (called with `None`) ends as `Exhausted` (2498). G1(b) counts that as a
  completeness failure and an M4c trigger, although the real cause is the missing idempotent.

### F5 (minor): Search reports admission refusals as `Exhausted`, which mislabels the G1(b) completeness signal

**Where.**
- A branch on which `win_ok` refuses every window becomes a "stuck branch": `choose_window` returns None
  (1876-1880).
- A refused S6 is only noted (2386).
- The final verdict is `Exhausted` (2389, 2498, 2737).
- SPEC §4.2/§4.7 defines `Exhausted` as "candidate stream ended with budget left (completeness signal)" and
  `Refused` as "an admission rule declined the call".

**Effect.** The G1(b) rerun raises l_max and c_max, neither of which changes window admission. So a pair
blocked only by the cost caps is still `Exhausted` after the rerun and is counted as a completeness failure
(red flag for U1, input to the M4c decision). The gate then does not measure what it claims to measure.

**Fix.** Count stuck branches whose cause was admission separately, and return `Refused`, or a separate
kind, when no other cause remains.

### F6 (minor): Search full ranks candidates at the planned detecting set, which is known not to be faithful, and has no fallback

**Measured counterexample.** HANDOFF_M4b.md ("Cloud session", item 3) records an exact basis on (1,3),
F0(1)E0(3)F0(2) → E0(3)F0(2)F0(1), d = 6, with witness rank 141 < 145 at M + 2 = 18 points. The witness set
with M + 2 points is therefore not a faithful detecting set in general. The SPEC §4.3 item 3 claim holds
only for the listed λ and degrees.

**Inside homprov this is not a soundness problem.** Evaluation is only ever used to prove independence or
nonvanishing.

**But `Search::full_inner` has no fallback.** It inserts candidates into a `DenseEch` at `plan.cosets`
with `plan.npoints` points (2624-2650) and stops at `kept.len() == dim`. On such a space even a complete
candidate stream can never reach `dim`, so the call ends `Exhausted` (2737), a false completeness failure.
`certify_full_plan` got the orbit fallback (640-668); the Search echelon did not.

**Downstream caution.** Do not use "witness cosets + M + 2 points" as a zero test or as faithful detection
(for example for pivot decisions in a point-based M6a scan).

**Fix.** Make the Search echelon retry at the orbit cosets (the same plan as the fallback) before reporting
`Exhausted`.

---------------------------------------------------------------------------------------------------

## A. certify_full and the orbit fallback: checked, sound

**Claim (homprov.rs:644-647):** "stacked rank = |S| at any coset set and points proves independence".

**The argument.** Suppose the k-linear map S ↦ (Eval_{p,c}(m))_{(p,c)} is well defined on Hom^d(src, tgt),
i.e. it depends on the map and not on its DAG. Then a dependency Σ a_i m_i = 0 forces the stacked rows to be
dependent, so rank = |S| implies that S is independent. Together with |S| = dim Hom^d, S is a basis. The
coset set and the number of points affect only false negatives.

**Preconditions, each verified in the code:**
1. **Same ends and degree for all maps.** Checked first (622-626), so all rows have the same column layout
   (`Eval::stacked`, smap.rs:1278-1291).
2. **Genuine leaves before ranking.** Checked at 627. Eval is then the localization functor applied to each
   node:
   - `Comp` is a block product per top coset (smap.rs:1114-1121);
   - `Lin` is linear (1122-1135);
   - `Local` is id ⊗ t ⊗ id via `window_block` (1197-1255);
   - `Dot` is diagonal in the flow basis (1105-1113).

   Bimodule maps preserve the top-coset decomposition. `poly_certificate` asserts block diagonality for added
   templates (tmpl.rs diff, the `tcos[i] != scos[j]` check). Hence the value depends only on the map. Without
   genuineness the argument would fail, because different DAGs of the "same" map could evaluate differently.
   B5 documents this.
3. **No poles.** Points have pairwise distinct coordinates in 1..2·10^9 < p = 2^31 − 1 (smap.rs:200-216).
   Every localization denominator is a product of root differences (`e_pair(..).inv()`, smap.rs:473;
   `localize`), so none vanishes.
4. **The dimension is right.** It is `[v^d] hom_gd(char src, char tgt)` (585, compared at 630). Every Exact
   call compares it with `sbim::hom_space`, and B2 compares it with hom60 (tests/homprov.rs:89-120). The
   formula is the self-dual-y EKLP form (schur.rs:1288); ladder characters are self-dual.
5. **The fallback adds only stronger evidence.** It uses `CosetSpace(lam, top(src))` orbits of the common
   tables, deduplicated (649-660), with M+2 points and then 2(M+2) points. The residual call (662, 674)
   appends columns and cannot lower the rank, so it is vacuous. That is documented and harmless.

**Search full** (2765-2771) runs the same `certify_full` after its own echelon, and end0 is full at d = 0.
Exact full certifies the whole-ladder hom_space basis (1034-1061).

## B. Genuineness of Search compositions: checked, sound

**Leaves that can occur, and why each is genuine:**
- `Id` and `Dot`. Dot is multiplication by s_ν of an intermediate alphabet, i.e. by an element of the middle
  ring, which is a bimodule endomorphism.
- `Local` over a HomSpace or Rickard template.
- `Local` over `Added` templates certified by `certify_added`:
  - entries homogeneous of the declared degree;
  - exact commutation with every top operator, which is the hom_space condition itself;
  - `poly_certificate` with block diagonality (tmpl.rs:371-402).

  `add_template` forces `certified: false` (tmpl.rs:354-357). Added templates never enter the `hom` map
  (tmpl.rs:456 inserts only hom_space ids), so they cannot pose as exact bases.

**Window splits.**
- π templates are always Added and certified (1729-1733).
- ι templates are either exact templates or Added and certified (1735-1741).
- The "extremes plus dots" PMats (1754-1817) are only raw material for those templates.

**Commutation and rex isomorphisms.**
- Exact degree-0 window maps; Hom^0 dimension 1 is required, plus g∘f = c·id (1470-1515, 1966-1984).
- If `rungs_commute` were wrong, this would fail closed or yield non-inverse pairs. Neither matters for
  certification, because every final certificate recomputes the form or rank from the final maps.
- `rungs_commute` itself (140-155) is correct: two splits from, or two merges into, the shared upright are
  associativity.

**Validity of the decomposition tree is not needed for soundness.** The final certificate (rank or form)
is computed from the final maps.

## C. Admission and partial results: checked, no leak

- **Exact.** Refusal (`admit`, 1005-1018) or a `try_hom_basis_cap` error leads to `Err(Refused)` before any
  map is produced (1023-1024). Success is returned only after `certify_*` (1042-1053, 1075-1091, 1094-1134).
- **Search.** A refused window is skipped (1461-1465, `win_ok`, `iso_ok`). Results are returned only through
  `certify_full` (2769), `certify_pairs` via `finish` (2282-2285), or `certify_closure` (2413, 2501).
  `run_call` commits only on Ok (2741-2758). `Chain` returns the first Ok and merges Missing (2832-2857).
- **Closure.** The degree loops return Err when `got < need` (1117-1119, 2498).
- The comm-window cap `win_ok || unknowns <= 20 000` (1479) ignores a measured cost-table veto. That is a cost
  risk only, not a soundness risk.

## D. Template cache keys: checked, complete

- `TmplLib.hom` is keyed by (Window(labels, rungs), Window, normalized degree) (tmpl.rs:420-423).
- The disk file name is `fnv64(field|kind|key)` (tmpl.rs:305-307). The header carries v2, alg stamp, field
  and the full key, and a mismatch is rejected (tmpl.rs:233-240), so hash collisions are caught.
- Window templates do not depend on the global λ. `Local::add` re-checks the window against the words
  (`check_window`, smap.rs:719-734).
- The `kernel_lazy` change keeps TMPL_ALG. Its output is the reduced kernel basis, which the row space and
  column order determine (linsys.rs diff; tests/hom_lazy.rs). Even a different basis would still be a valid
  basis.
- Search's `splits` / `comms` (TmplIds) are the one cache without a complete key: see F3.

## E. Closure: checked, sound when L_p is indecomposable

**The determinant argument holds.** Suppose:
- the n maps lie in Hom(R^λ, L_p), which is free of rank n with generator degrees c_diag(b_p) (B1 checks
  `hom_gens` against Nakayama);
- their degree multiset equals c_diag (775-781);
- their identity-coset n×n determinant is nonzero at a point (782-791).

Write m = A·b with A over R^λ and homogeneous entries of degree deg m_i − deg b_j. Then det A has degree 0,
so it is a nonzero constant and m is a basis. Ends are correct by construction but are not checked (F3).

## F. Multiplicities and G1 targets

- `g1_targets` (60-93) takes the summand targets from `term_mult(p, [r])`. Their convention (ι: L_q → X of
  degree s for v^s b_q) matches `xengine::decompose_direct` (xengine.rs:719-721).
- Because X and b_q are self-dual, the multiplicity polynomials are symmetric. Counts alone would therefore
  not reveal a sign error in this convention; it is consistent with xengine by inspection.
- `mult` itself is not validated (F2).
