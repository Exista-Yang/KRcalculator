# D_sinister: B5 through localized sinister leaves (design)

Date 2026-11-12. Read-only design on branch `m4c-rexplan` (4e886f6, `.claude/worktrees/revrex`); no cargo, no binaries.
Inputs: B5_census_analysis.md (the 42 worst steps, the chosen 3-upright joint words, groups A-E), M4C_PLAN §3.1/§3.3,
D1 §1.3/§1.7, L1 §2.1/§2.3, B1 §1/§4, B2 §3/§6, A1/A2, C_m4c C1/M2, `src/{leafgen,leaves,dcoset,rexplan,smap}.rs`.
New Python checks (all runs < 3 s) are in `b5/tools_s/`: `shapes.py`, `schub2.py`, `tcost.py`, `msearch.py`, `flagsub.py`,
`cover.py`; outputs `*.json` there. Paper keys as in L1 (SLL = EKLP 2401.03053).
Tags: [P] proved in the paper; [V] computed here (Python replica of krhom combinatorics, one random point over F_p for
ranks); [K] inference, decided by `certify_added` / the rank diagnostic / `certify_pairs`; [E] estimate.

## 0. Verdict

1. **Literal route (compile the sinister diagram in SLL's own charts, no exact call at all): no-go for B5.**
   - Compiling the diagram is trivial: localized, the elementary light leaf is *one* Push, its flip one Pull, the
     sprinkle one Dot (§1.1). Nothing is gained by building it from crossings and caps.
   - The cost sits in the charts. The descending half of a rung needs the giver column's cells on their own
     uprights: up to 4 blocks inside the K-window and 5 globally. That is `Pad` twice (λ·0·0), plus rex chains
     between charts on 4-5 uprights. The planner has never run there, and C_m4c M2 measured 4-block words at up to
     1.3·10⁷ flows. This is 40 of the 42 steps (B5 census §2, class a-pad).
   - **New obstruction for group C [V].** A `Local` step family is zero on every top-to-top column whose prefix or
     suffix flow differs from the target's top flow. So the number of window flows at c_w must be at least G(1).
     Pair 383's window F·E → ∅ on (3,3) has 3 such flows against G(1) = 6, so *no* window family can exist there,
     closed-form or exact.
2. **Go for B5-S, the "window sinister" route.** It puts the localized leaves directly on the 3-upright joint
   windows of the census analysis. There is no `Pad`, no new `MapNode`, and certification is unchanged. Each step
   uses one of four closed forms:
   - S2: `Pull`/`Push` of B1;
   - S1: `Push∘Dot` / `Dot∘Pull`;
   - **one new closed-form kind, `PullTw`/`PushTw`**: a Pull/Push twisted by one incidence Euler factor
     (§1.3, the window image of eq. (133)), with sprinkles as Dot/Lin nodes;
   - `Id` plus a rex chain for group E.
3. **What it buys [V].**

   | outcome | targets (of 101) |
   |---|---|
   | closed form, every template predicted ≤ 300 s | **67** |
   | exact window ≤ 10⁵ with G_loc = G | 9 |
   | Out template over the cap (458) | 1 |
   | recorded exact call of 1-4·10⁵: 137/145 (318 975), 385/467 (389 029) | 14 |
   | rank test on tiny windows: 410 (3 583), 283/293 (58) | 6 |
   | non-local: 383 | 4 |

   No closed-form template reaches a 9-root label: every window label is ≤ 6, so GrTable(r, ≤ 6) suffices.
4. **Per pair [V, worst step only].**
   - **N = 5:** 313 and 326 are closed (both move up relative to the census analysis); 410 depends on a rank test.
   - **N = 6, 32 pairs:**
     - 18 closed: 24, 78, 127, 212, 219, 227, 273, 278, 281, 285, 313, 319, 326, 332, 336, 376, 381, 411;
     - 5 with cheap exact windows: 45, 82, 292, 335, 418;
     - 458: one Out template over the cap, or a 210 002 record;
     - 137, 145, 385, 467: records;
     - 283, 293, 410: rank tests;
     - 383: `Pad` K-chart or a gap.
5. **Honest outlook [E].**
   - N = 5 complete: 70-80 %. It hinges on 410's rank test and on 326's unanalysed hidden final step (147 091).
   - N = 6: the 23 clean pairs 70-75 %; all 32 pairs 45-55 %.
   - Effort: 3-4 sessions, +1-2 if `Pad` is wanted for 383 or 410.
   - **Stopping point:** certify the closed and cheap-exact pairs, record 137/145 and 385/467 only if time allows,
     and leave 383 (plus whatever fails) as documented gaps.

## 1. Mathematics

### 1.1 The sinister ELL is a single correspondence [P pieces, K composite]
- **Ascending leaf.** ELL([m,n]): X_m∘[J+s] → X_n (Def 7.3, Lemma 8.2) uses only left-facing crossings and ccw caps.
  Both are `Push` correspondences (B1 §1: left crossing entry e(C, A∖T); Trace = Push).
  - The target X_n = [I + v₁ + … + v_k] (Lemma 7.2) consists of merges only, so it is flag-type.
  - Push-forward composes, and Euler weights telescope along matching flows. Hence
    **ELL([m,n]) = Push(X_m[J+s] → X_n)**, with entry e_{X_n}(c)/e_Y(p).
  - Its degree is Σ(shift(src) − shift(tgt)) = def[m,n], matching Lemma 7.4.
  - ELLP([m,n], f) = Push ∘ Dot(f) (Def 7.14).
- **Descending leaf.** ELL([n,m]) = (id ⊗ cw cap) ∘ (D ELL([m,n]) ⊗ id) (eq. 133). D swaps Push and Pull (Def 5.10),
  and the cw cap is E2 `Eps` = Pull to []. So ELL([n,m]) is one Pull from the flag-type X_n[−s].
  This is consistent with 1⊗ ↦ 1⊗ (Lemmas 8.9-8.10).
- **Light leaf.** SSLL = id_Z ⊗ (that) is one `Local` node. The light leaf of Def 7.28 is therefore a chain
  Local(Push∘Dot) ∘ rex ∘ Local(Pull) through the chart words Z·X_m and Z·X_n of Notation 7.22.
- **So the sinister diagram itself costs nothing.** Every difficulty is in the charts and in the rex moves between
  them (§0.1).

### 1.2 Window form: the flag-target lemma [K, standard; certification decides]
**Lemma.** Let X be a flag-type window word (one flow per coset of its support) and Y any word with the same ends. Then
φ ↦ φ(1⊗) is an isomorphism

  Hom(BS X, BS Y) ≅ { m ∈ BS(Y) : m vanishes at every flow of Y over a coset ∉ supp X }.

**Proof sketch.**
1. BS(X) ≅ H_T(X_t). Flag-type forces finite fibres over the fixed points, the map is birational, and Schubert
   varieties are normal, so this follows by Zariski's main theorem.
2. H_T(X_t) is a cyclic bimodule, generated by 1, because restriction from G/P is surjective.
3. BS(X) is torsion-free with one line per coset, so the relations of 1⊗ are detected coset by coset at supp X.
   Any m vanishing off supp X satisfies them.
4. The Out side is the dual map, with entry m(p)·e_X(c)/e_Y(p) (SLL Def 5.10, Atiyah-Bott pairing).

**Elements of BS(Y) we can write down.**
- For each rung j, let M_j be its moved set. M_j is a block of the rung's intermediate parabolic, so e_k(M_j) ∈ BS(Y).
- So is every e_k(U_u^ℓ) of an upright content at a level.
- Hence the **incidence Euler factor** ε(M_j, R) = Π_{γ∈M_j, s∈R}(x_γ − x_s) is an element of BS(Y). Here R is
  another moved set or an upright content. ε(M_j, R) vanishes exactly on the flows where M_j meets R: a cell
  condition "rung j may not take roots from R".
- This is the localized value of the cw cap/cup in eq. (133). Compare E2 `Eta`, whose value is e(B, A)
  (smap `ElemKind` doc).

**The twisted kinds.**
- **PullTw_ε(X → Y)** has entry ε(p) at (p, c), and **PushTw_ε(Y → X)** has entry ε(p)·e_X(c)/e_Y(p). Both are
  genuine when ε kills every bad flow.
- Degree: shift(Y) − shift(X) + 2·deg ε.
- With ε = 1 and supp Y ⊂ supp X this is B1's Pull/Push: S1, and S2 with the roles swapped.

**Step family.** For every d ∈ supp G, with multiplicity:
- In: ι = Dot(β) ∘ Local(PullTw_ε);
- Out: π = Local(PushTw_ε) ∘ Dot(β).

The sprinkle β = e_k(M_j) or e_k(U_u^ℓ) plays the role of SLL Def 7.33. Each family has **one template per window and
direction**; the sprinkles are DAG nodes.

### 1.3 Correctness and completeness
- **Genuineness.** Templates pass `TmplLib::certify_added` (via `leafgen::build` / `build_in`); this is unchanged
  trust. Dot and Lin nodes are genuine by construction (smap.rs). No theorem certificate is requested.
- **Invertible top-to-top block.** This is checked at the found family (search, one point) and again by
  `build_leaves`' rank diagnostic. Lemma A (D1 §1.3) then composes the steps as today. Acceptance remains only
  `certify_pairs` / `certify_closure`, plus π∘ι = id blockwise (B2 §5).
- **Locality.** Windows are used only where G_loc = G (census analysis), so global columns = window columns.
  - The search requires |T| = G(1). This held on all 13 successes [V].
  - Where G_loc ≠ G (410, 283, 293), a cheap rank test decides.
  - Where the window has fewer flows than G(1) (383), no window family can exist, by the bound in §0.1.
- **Failure is visible.** A wrong ε fails `certify_added` (`interpolate` → Err), and a deficient family fails the
  rank diagnostic. Neither can produce an accepted wrong map.

## 2. Census of the 42 worst steps under B5-S [V]
Windows are the census analysis's chosen joint words (the G_loc = G variant when it exists), after `win_split`.
- Template cost uses B2's fit `predict_template_secs` = 5e-7·M·n² + 7e-9·M³ (±2×), with n the larger flow count.
  M = `monos_of_deg` at the base entry degree: In d_min − 2·shift(b); Out d_min − 2·shift(a·ρ).
- Search: `msearch.py` tries ≤ 3 incidence factors with exact degree that cover all bad flows. It then picks
  sprinkles greedily, by ascending degree, for rank on the columns over F_{2³¹−1}.

| group | steps (targets) | B5-S realization | max template flows | predicted template s (In / Out) | max exact unknowns left |
|---|---|---|---|---|---|
| A | 2 (2) | S2 Pull/Push | 20×30 | ≤ 0.1 / 0.3 | 0 |
| B | 25 (57) | S2 on 8 (18), S1 on 3 (6), PullTw on 8 (24); exact on 6 (9): 45 v⁰, 335, 285, 292 v⁻²+v⁰, 82 v¹+v³, 418 (no flag-type target or substitute exists) | 240×16 (137) | ≤ 5 / ≤ 61 | 51 816 |
| C | 7 (19) | 313: PullTw on its G_loc = G window (283 842 as exact). 137/145: none (cover degree 4 > 2), record. 410, 283, 293: rank test. 383: non-local | 240×40 (313) | 17.5 / 235 | 318 975 (rec.), 3 583 |
| D1 | 3 (6) | S2 on all three | 100×120 | ≤ 25 / ≤ 65; 458 Out ≈ 1 900 | 210 002 (458 Out, optional) |
| D2 | 2 (9) | 326: PullTw, ε = ε(M₁, M₂) (hand-checked, §4). 385/467: none (cover degree 7 > 5), record | 120×30 | 2.9 / 25 | 389 029 (rec.) |
| E | 3 (8) | Id, then a rex chain W_{t'}·ρ → W_t (B2 planner) | – | – | detour chains (B2 budget) |

**Observations [V].**
- Every chosen window is small: 12-720 flows on the larger side, median 50.
- supp(b_w) ⊂ supp(a_w·ρ) holds in all 39 non-Id windows. The source support is strictly larger in all but the 3
  S1 windows. That is why plain Push fails there and the ε-twist is needed.
- All 13 PullTw successes need a **single** incidence factor. Examples:
  - 127 (n = 6): ε(M₂, M₃) "ρ may not re-move the root just moved", with sprinkles e₁(M₀), e₁(M₁), e₂(M₀), e₂(M₃),
    e₃(U₀²);
  - 411: ε(M₁, U₀⁰);
  - 137 (B): ε(M₀, U₁⁴).
- Exact degree matched d_min in every success. For n = 1, the relative codimension ℓ(t̄_a) − ℓ(t̄_w) predicts the
  degree in 15 of 16 cases; 410 is the exception (`schub2.py`).
- The two failures (D2 385/467, C 137/145) need a minimal leaf that is not a product of incidence factors. They would
  need a factorial-Schur / double-Schubert ε (AJS-Billey localization) [K, future].

## 3. Construction, step by step (per step with `step_score` > budget)
1. **Words.** Use the census analysis's joint pair (W'_{t'}, W'_t), with the projected Grassmannian prefix and
   G_loc = G (`leaves::joint_pair`, `dcoset::window_gloc`, as proposed there). Then `homprov::win_split` gives
   (Λ, a_w·ρ, b_w).
2. **Shape test (pure combinatorics).** Compute:
   - window flows and configurations;
   - `is_flag_type` on each side, support inclusion;
   - d_push = shift(aρ) − shift(b) and d_pull.

   The classification is:
   - S2 if aρ is flag-type, G = v^{d_pull} and n = 1;
   - S1 if both supports are nested, b is flag-type and d_min = d_push;
   - Tw candidate if b is flag-type.
3. **Twist search.**
   - Find bad flows: flows of aρ over cosets ∉ supp b.
   - Choose ε from the factor list (M_j × {M_i, U_u^ℓ}), ≤ 3 factors, with deg ε = (d_min − d_push)/2, covering the
     bad flows and not vanishing on the columns T.
   - Choose sprinkles for each d ∈ supp G, by greedy rank over T at one point.
   - Output: ε, a list of (d, β), and the predicted template seconds.
4. **Templates.** Build only the base template for each direction with `leafgen::build`/`build_in`:
   - Pull/Push (S2), Push/Pull (S1), PullTw/PushTw.
   - Cap: 300 s predicted. Above the cap, that direction uses the exact window under `Leaves::admit` (≤ 10⁵, or
     ≤ 4·10⁵ with a record). Today that applies only to 458 Out.
5. **Sprinkle DAG.**
   - e_k(U_u^ℓ) = `Dot{level ℓ, upright u, ν = 1^k}`.
   - e_k(M_j) = Σ_i (−1)^i e_{k−i}(U_g^{j−1}) h_i(U_g^j), written as `Lin` of `Comp(Dot, Dot)`, where g is the giver
     of rung j and U^{j−1} = U^j ⊔ M_j.
   - Shift both by the window offsets (`at`, `first`) to the global word.
6. **Lemma A.**
   - `build_leaves` is unchanged; only `step_family` gains the closed kinds, tried before the exact windows.
   - Detours:
     - adopt the joint words in `plan_words` when the neighbouring steps stay ≤ 10⁵;
     - otherwise, for non-final steps, a rex chain W_{t'} → W'_{t'} and W'_t → W_t via `rexplan::Planner::plan`
       (B2 budget: templates ≤ 1 200 flows, Exact ≤ 10⁵);
     - for final steps, the exit chain W'_q → L_q (as `alt_final`).
7. **Pairs.** `pairs_leaves` builds π = exit ∘ LL and ι = LL^∨ ∘ entry, then runs `certify_pairs` unchanged and checks
   π∘ι = id at every coset of L_q at two points (B4 procedure).

## 4. Code changes (no frozen file touched; no `Pad`)
| file | change |
|---|---|
| src/wsl.rs (new, "window sinister leaves") | `Region {Moved(j), Upright{u, level}}`, `Incidence {j, region}`; `shape(lab, a, b) -> Shape` (flags, supports, shifts, kind); `bad_flows`; `find_family(lab, a, b, cols, g, d_push, cap) -> Option<Family {eps, sprinkles: Vec<(i32, Sprinkle)>}>`; `sprinkle_dag(ar, e, at, first, &Sprinkle) -> MapId`; `predict_secs`. Combinatorics plus one F_p point; `confs`/`flows_full` re-implemented on `FlowIndex` masks |
| src/leafgen.rs | `CfKind::PullTw(Vec<Incidence>)`, `CfKind::PushTw(Vec<Incidence>)`. `ClosedForm::new` checks a flag-type side, that ε kills every bad flow (asserted combinatorially) and that ε ≠ 0 on some flow, and computes the degree. `winloc` gives entries ε(p) and ε(p)·e_X(c)/e_Y(p), using `euler` and the masks. `interpolate`/`build` are unchanged |
| src/leaves.rs | `StepKind::{Closed(ClosedKind), Rex}`; `step_family`: Id → Closed (wsl) → Window exact → Whole; `closed_family()` (per direction, cap fallback); `plan_words` calls `joint_pair`; `StepRec` records the kind and template seconds; the stamp becomes `LEAVES_STAMP` v2 (`…,wsl=S1S2Tw`) |
| src/dcoset.rs | `window_gloc` (local types on the window labels, then `step_gf`), as in the census-analysis plan |
| src/bin/gp0.rs | B5a: every step over budget (not only the worst), the joint search, `wsl::shape` + `find_family`, predicted template seconds, and a per-target verdict {closed, exact ≤ 10⁵, record, rank, non-local} |
| src/smapio.rs | provenance tag `[wsl Tw ε=…]` on Added templates (they are re-certified on import, B2 §6.3) |
| tests/wsl.rs (new) | gates V-S0 to V-S2 below |

Untouched: smap.rs (only Id/Dot/Lin/Comp/Kron/Local), homprov.rs, tmpl.rs, sbim.rs, linsys.rs, sindec.rs. The
theorem-level Assoc is used only where B2's planner already uses it, which today is nowhere.

## 5. Validation plan
- **V-S0, exact oracle on small λ.** λ ∈ {(2,1), (1,2,1), (2,1,2), (2,2,1), (2,2,2)}, every step of every G1(a) Lemma A
  plan, both directions. Whenever `wsl` returns a family:
  - each base template lies in the exact `try_hom_basis_cap` window (stacked rank);
  - the family has degree multiset G and an invertible top-to-top block at 2 points;
  - the composed leaves have stacked rank equal to dim Hom^d (A2 V0 style).

  Also count how often `find_family` fails while the exact window succeeds; this measures incompleteness.
- **V-S1, SLL examples.** Ex 7.39 (Trace∘Dot, which is S1), Ex 7.40/7.41 (I = (2,1), J = (1,2)), and the D2 326
  window, hand-derived: ε = Π_{s∈S}(x_γ − x_s), β ∈ {1, x_γ}.
- **V-S2, regressions with closed kinds preferred.**
  - A2 closure 12/12, both primes, with equal degree multisets and replay from layer A;
  - B2 rex-class 7 pairs (159, 71, 323, 322, 334, 419, 313-rex), with form [[1]] and π∘ι = id;
  - G1(a) closure and pairs (3,3) and (2,2,2), 480/480.
- **V-S3, B5a gate (Rust).** It must reproduce §2: 13 PullTw successes, the 2 failures, and the S1/S2/E classes. It
  lists all hidden over-budget steps (12 known, B5 census §5.3) with their classes. Combinatorics only, minutes.
- **V-S4, (3,3,3).**
  - Order: N = 5 (326, 313, then 410's rank test); N = 6, the 23 clean pairs; then the records (`m4c measure` for
    318 975 and 389 029); then 283/293.
  - Both primes. Downstream: sindec wiring, d∘d, K₀, χ, replay (B2 §6.3). G-L5 counters (DAG nodes, index bytes).

## 6. Effort (cloud sessions of about 4-5 h)
| session | content | exit |
|---|---|---|
| S1 | `wsl.rs` combinatorics, gp0 B5a (all steps, joint words, shapes, ε search, verdicts) | V-S3 |
| S2 | leafgen PullTw/PushTw, `step_family` closed kinds, sprinkle DAGs, `joint_pair` in `plan_words` | V-S0, V-S1, V-S2 |
| S3 | (3,3,3) runs: N = 5, then N = 6 clean; detours; both primes; downstream | V-S4 for the clean set |
| S4 (optional) | records 137/145 and 385/467; 410/283/293 rank tests; report | final go/no-go list |
| S5-6 (only if wanted) | `Pad` K-chart for 383 (local window 690) and as fallback for 410/C/D2 | – |

**Total: 3-4 sessions, +1-2 for `Pad`.**

## 7. Risks
| # | risk | likelihood [E] | detection | consequence |
|---|---|---|---|---|
| 1 | hidden second steps (12 known, more possible) are not closed-form | medium | V-S3 | exact ≤ 10⁵ / record / gap |
| 2 | detour chains to the joint words are not found within budget (32 of 42 exceed one 10⁵ Exact) | low-medium | B5a planner dry run | adopt the joint words in `plan_words`; otherwise the next candidate |
| 3 | the flag-target lemma is wrong for some window (it is [K]): `certify_added` fails | low | V-S0, `interpolate` Err | exact window for that step |
| 4 | template cost model off by > 2× (monomial-bound); Added templates are re-certified on every import | medium | B1/B2 logs, timing in S2 | lower the cap; exact windows; replay slower by minutes |
| 5 | rank tests fail for 410/283/293 (G_loc ≠ G) | medium | S4 | `Pad` K-chart (local 58-3 583) or a gap; N = 5 then misses 410 |
| 6 | records for 318 975 / 389 029 are censored (> 600 s or > 3 GB) | medium | `m4c measure` | gap, or `Pad` (local 838 / 1 397) |
| 7 | downstream M6a cost of sprinkle `Lin` DAGs | low | G-L5 | merge the sprinkle into the template (n templates) |
| 8 | 326's hidden final step 147 091 fails | medium | V-S3 | N = 5 incomplete without a record |
