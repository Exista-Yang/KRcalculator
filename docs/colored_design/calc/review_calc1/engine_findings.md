I reviewed the engine part of the diff by reading code only: src/sindec.rs, swire.rs, sengine.rs, local.rs, smap.rs, tmpl.rs and boxrun.rs (worktree /Users/existair/Documents/KRcalculator/.claude/worktrees/calc). Nothing was run.

**Bottom line:** I found no path that silently accepts an uncertified map, and no wrong-output path in reduced coordinates or the column-restricted localization. There are two real defects. Both make a job fail when it should not; neither produces a wrong answer.

## Findings

### 1. [PLAUSIBLE, medium: correct jobs rejected] Reroute family checked against the wrong multiplicity when the prefix type is decomposable
- **Where:** src/sindec.rs:975-977 (`a.len() != mx` in the M5i source branch) and :997-999 (same test in the fallback). `reroute_family` is at :1172-1190.
- **What happens:** `reroute_family` takes the type `a` whose canonical ladder is w[..j] (:1180) and never checks that L_a is indecomposable.
  - If L_a is decomposable, `dec_word(a, R·w[j..])` (:1189) only decomposes B_a ⊗ R·w[j..] ≅ B_p ⊗ R. The family then has cnt pairs per kept (q,s).
  - The code instead requires mx pairs, the multiplicity in char X = char(L_p)·R, which also counts the extra summands of L_p.
  - The doc comment's claim ("a basis of the multiplicity space of b_q⟨s⟩ in char X") is false in this case.
- **When it happens:** for prefix-propagated types the reroute usually lands on a decomposable prefix. C2_4_report shows exactly this: T2 = Kron(e_T5, F1.2); job 1812 reroutes through `dec_word(T5, E2.1·F1.2)`; job 132 through the peeled (decomposable) 1,1,0,0/0,1,1,0/… type.
- **Failure scenario:** p is prefix-propagated, R commutes only past the last rung of L_p, and some kept b_q⟨s⟩ of b_p·R also occurs in extras(L_p')·R'·R. Then mx > cnt and the job fails with DecErr::Check "M5i reroute of … : cnt pairs …, char X has mx". The family was in fact usable: after the E-projection, the rank-cnt selection plus Σιπ = E (:1012-1023, :1453-1459) would certify it.
- **Knock-on in replay:** a missing layer-A target that reaches the fallback (:999) gets a Check error without the original "not in layer A" text. calc::classify (calc.rs:76) then reports Engine instead of Library.
- **Not a soundness hole:** when the sizes happen to match, the downstream checks still certify the result.
- **Fix:** accept `a.len() >= cnt` and use `a.len()` as the selection width, or return None for a decomposable `a` so the other routes are tried. Keep the original error text in the :999 message.

### 2. [CONFIRMED by reading; low; fails closed] Layer-A keys do not identify which e_q / e_p a map was built against
- **Where:** M5I_STAMP is a constant (sindec.rs:54). It is the only e_q / e_p token in the pairs, idem and closure keys (:1714, :1732, :1737). It was not bumped when prefix propagation changed how e_p is built for 5 of the 7 types (3cfd150 → 4225203). KRHOM_M5I_NO_PREFIX and the choice of peel route are not stamped either.
- **Why it matters:** when the extra summands sit in degrees ±1 or have multiplicity, the idempotent with image B_q is not unique.
- **Failure scenario:**
  - A cache is filled under one e_q, for example by a KRHOM_M5I_NO_PREFIX validation run, or by a peel route that differed through cost-table or c27 differences accepted under `layer_a_route_stamp`.
  - A later run loads the `eq=m5i-peel1` pair files. They pass `dechom_postcheck` (:1241).
  - Π∘Ι then equals the stale e_q', not the builder's e_q, so `check_pi_iota` (:1431/:1437) returns DecErr::Check. Nothing is recomputed, so the job stays failed until the file is deleted.
  - Closure files are not affected: they are re-wrapped and re-certified (:1617-1622, homprov.rs:945).
- **Fix:** put a fingerprint of e_q's construction into `eq=` / `ep=`, or treat a Π∘Ι failure on loaded maps as a cache miss.

## Checked and sound

**Idempotent checks.** e² = e is checked at several points:
- `check_idempotent` (homprov.rs:837): at 2 points, on the rep(t) coset of every table of char L, with rank = c_t(b_p)(1) and top scalar 1.
- `idem_checks` (sindec.rs:751-776): on every coset at a third point, with Σ rank = rank B_p.
- `lib_shape` (sengine.rs:224-249): `reduced_of` asserts e² = e, QP = I and PQ = E on every coset at 2 fixed points (smap.rs:1392).
- Again at every scan point and coset through `RedCache` → `reduced_of`.

The rank argument (rank 1 at rep(p) and 0 for the extra summands) forces Im e ≅ B_p.

**Every new map source is re-certified downstream:**
- the prefix Kron, the peel construction and the Exact fallback;
- the X = L_q shortcut and the commutation shortcut;
- the reroute family and the xq kind=idem pairs;
- closure re-wrapping.

The certification is DecHom pre/post checks, `form_twice` (2 points plus genuineness), Neumann orthogonalization, and `check_pi_iota` at 2 points: Π∘Ι = δ with e_q on the diagonal, and Σιπ = E = Kron(e_p, R). `swap_iso` normalizes g∘f at one point only, but a wrong inverse would fail `check_pi_iota`. Every ι/π is e_q-wrapped on the side reduced coordinates need: the provider's `with_idem`, `pairs_leaves_eq`, the first/comm branches, and `orthogonalize`, which keeps π = e_q π.

**Reduced coordinates.** One `RedCache` per scan point serves lib_map, assemble, to_flows and extract_stage, so the basis is consistent within a scan.
- `to_flows` gives P_b f Q_a = e_b F e_a.
- Closure entries come out as h·e·D·e·g.
- The frozen trace is component-level (sengine.rs:471) and pivots are c·I, so both are basis-invariant.
- Rank of e at a coset is constant wherever e is defined (trace argument, n < p for F31/F32b), so the shapes match at box points.
- The reference path KRHOM_DCID_REF refuses libraries that contain idempotent types.

**Column-restricted localization** (tmpl.rs `localize_cols`, local.rs `inverse_cols`, smap.rs window_block):
- It is exactly equal to the full localization: L_t·(Φ·L_s⁻¹[:,js]) = (L_t·Φ·L_s⁻¹)[:,js] by associativity.
- `inverse_cols` is a correct pivoted solve of L_s X = E_cols. It returns None on singular input, so it panics exactly where `inverse` does.
- Columns are gathered by a loop identical to the consuming loop, and the partial results are local to the call (never stored in `Eval.wins`), so no unfilled column can be read.
- No fallback or runtime verification is needed; the "identical point hashes" test is only extra confirmation.

**c27 `route_header` filter.** It only affects acceptance under `layer_a_route_stamp` (loaded maps are re-certified) and the WLib hash. Box checkpoints made with differing libraries are still separated by the record/trace bindings (boxrun.rs:138-146). The calc route refuses KRHOM_LEAVES_C27 and uses NoProvider, so the Exact fallback cap is 0 there.

**Numeric edge cases.** Zero pivots are handled in `inverse_cols`, `swap_iso` (c = 0 is an error), and in `build_idem` / `dec_one` (singular form is an error; projected rank must equal cnt). Multiplicities are guarded as non-negative before the u128 casts. No new overflow found.
