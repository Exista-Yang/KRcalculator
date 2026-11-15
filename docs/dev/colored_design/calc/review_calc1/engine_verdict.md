I verified both engine findings by reading code and docs at 1a9c8de. Nothing was run. Neither finding is confirmed: #1 is UNCERTAIN and #2 is NOT-REAL. Neither one can produce a wrong certified result.

## Finding 1: reroute family checked against the wrong multiplicity — UNCERTAIN (minor if it ever triggers; fails closed)

**Confirmed by reading:**
- `reroute_family` picks `a` as the first type in the expansion of char(w[..j]) whose canonical ladder is w[..j] (src/sindec.rs:1187). It never checks that `a` is indecomposable, so a decomposable `a` is allowed.
- It then calls `dec_word(a, R·w[j..])` (sindec.rs:1195). For a decomposable `a` this decomposes only B_a ⊗ R·w[j..]:
  - `dec_composite` asserts that the kept multiset equals `term_mult(a, rs)` (sindec.rs:1495-1500);
  - `term_mult` expands b_a·rs only (src/schur.rs:1389-1397).
- A decomposable `a` is reachable. Prefix-propagated types have a decomposable prefix by construction (`idem_by_prefix`, sindec.rs:725-744, which also requires b_{p'}·R' = b_p). C2_4_report.md:158-159 shows 1812 rerouting through T5 and 132 rerouting through the peeled type.
- In the prefix case the family therefore has exactly cnt pairs. The M5i branch requires `a.len() == mx` (sindec.rs:978), where mx is the multiplicity in all of char X (sindec.rs:960). The doc claim at sindec.rs:1169, repeated in C2_4_report.md:75, is false when `a` is decomposable.

**Not confirmed: the trigger.** The job fails only if some kept b_q⟨s⟩ of b_p·R also occurs in extras_p·R.
- I found no instance. Small Hecke-algebra attempts did not produce one either.
- In the shipped (2,2,2,2) library this case already ran without failing:
  - 1812 (T2 ⊗ E2.1, rerouted through the decomposable T5) is certified for both primes after C2-7 (C2_7_report.md:137);
  - the other decomposable-source reroutes (109, 1789, 2052, 2086, 2275, 2298) are certified (C2_4_report.md:146-148).
- So there is a latent false rejection for future libraries, not a confirmed failure.

**The knock-on claim about :999 is NOT-REAL.** The fallback reroute (sindec.rs:993-1001) runs only when `ex` is None:
- With `ex` Some and `xq` None, `fam_tried` is true (:950), which disables the fallback.
- With `ex` Some and `xq` Some, the `xq` branch at :982 is taken first.
- With `ex` None, mx = cnt (:959). char(L_a)·w[j..] = b_p forces b_a·w[j..] = b_p, so the family has exactly cnt pairs and the :998 message cannot fire.
- When `reroute_family` itself errors or returns None, the original "not in layer A" text is kept (:995-996), so `calc::classify` (src/calc.rs:76-77) still reports Library.

## Finding 2: layer-A keys do not identify which e_q / e_p a map was built against — NOT-REAL

The code facts are correct. M5I_STAMP is a constant (sindec.rs:54) used in the keys at :1714, :1732 and :1737. Loaded pairs are not re-wrapped with the current e_q (`target_pairs_keyed`, :1241-1250). `check_pi_iota` compares Π∘Ι numerically with the builder's e_q (:1431-1432). The failure scenario still cannot be reached:
- **Every decomposable type in the shipped libraries has distinct extras in degree 0.** C2_2_report.md:100-110 says "All extra summands sit in degree 0"; T7's two extras are distinct (b_r1 + b_r2).
- **That makes e_p unique.** Hom^0 between non-isomorphic indecomposables is 0, so End^0(L_p) ≅ k^{1+#r}, and the idempotent with top scalar 1 and image ≅ B_p is unique. C2_4_report.md:55-57 states this, and the End^0 tests check it for T5 and (3,1,3).
- **So every construction gives the same map.** The prefix Kron, peeling and a `KRHOM_M5I_NO_PREFIX` run all give the same e_p, a stale file's e_q' equals the current e_q, and `check_pi_iota` passes.
- **The finding's stated condition is wrong.** It says the idempotent is not unique when extras "sit in degrees ±1 or have multiplicity". Degree-0 multiplicity still gives a semisimple k × ∏ M_m(k) with a unique such idempotent. Only extras at nonzero degree shifts break uniqueness, and no current type has one.
- **The calc route refuses `KRHOM_M5I_NO_PREFIX`** (calc.rs:650).
- **Even for a future λ with shifted extras, it fails closed** with DecErr::Check. No wrong output is possible.
- **Remaining non-defect:** the doc comment on M5I_STAMP (sindec.rs:51-53) still describes only the peel construction, although prefix propagation also builds e_p.

## Summary
- No blocker and no major issue in the engine area.
- #1: UNCERTAIN. The multiplicity mismatch is real in the code; the trigger is unconfirmed and was not hit in the shipped data. Minor and fail-closed if it ever triggers. The cheap fix is to compare against cnt and use `a.len()` as the selection width, or to return None from `reroute_family` when `a` is decomposable.
- #2: NOT-REAL.