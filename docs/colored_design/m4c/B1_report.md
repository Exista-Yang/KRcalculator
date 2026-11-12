# M4c B1 report: closed-form generators in localized form (R1 test)

Branch `m4c-gen`, from `origin/m4c-census` 595068c. New files only:
- `src/leafgen.rs`;
- `tests/leafgen.rs`;
- this report;
- one `pub mod leafgen;` line in `src/lib.rs`.

No frozen file and no `src/leaves.rs` was touched. Field: F31 (p = 2³¹−1).

## 0. Verdict on R1

**R1 holds for the generators and for the switchbacks.** Rex-type moves can be written as closed forms in krhom's
localized representation, and every closed form tried here validated:
- every one was certified by `certify_added`;
- every one lies in the exact `hom_space` window wherever that window was computed (≤ 10⁵ unknowns). All those
  windows have dimension 1, so each closed form is *the* map, up to a scalar.

The closed forms are not hand-tuned per case. They follow from one principle (§1): every Frob generator and every rex
move through a maximal coset is a **pull-back or push-forward along a Bott-Samelson correspondence**. Localized, that
is either 1, or a ratio of equivariant Euler weights, on matching flows.

**One residual risk sits outside the generators: the up-up/down-down step in the middle of the atomic braid.**
- The braid route of CORE Lemma 3.12 is ψ, then (up-up, down-down), then φ, through the fully merged expression. The
  middle step M → M′ is a degree-0 isomorphism between two flag-type words over the whole 3-cell window.
- At cells (3,3,3) that window has **1 680 flows**. The Assoc validator refuses it under the plan's ≤ 500-flow rule.
  As a Pull template it is also beyond what interpolation and `certify_added` can afford, because the certificate's
  cost grows with the monomial count, not the flow count (§4).
- Below the cap it is already slow: the 140-flow Assoc at cells (3,1,3) took 245 s to certify, and the 210-flow one
  at (2,3,2) took 1 509 s (§4). At (3,2,3) the middle step has 560 flows and is already refused.
- So the braid composite is **validated end to end at cells (1,1,1), (1,2,1), (2,1,2) and (2,2,2)**. At cells (3,1,3)
  and (2,3,2) every piece is certified, including the middle Assoc, but the composite was not run this session.
- At (3,2,3) the switchback pieces are certified and the middle Assoc (560 flows) is refused.
- At (3,3,3) the middle step is refused. The switchback pieces (400 × 84 flows) could not finish this session: they wait
  on the one-time GrTable(3, 9) build that the exact SBim needs (§3.2).
- This is a B2/B4 decision, not a mathematical obstruction. §6 gives the options.

## 1. The localized closed forms (src/leafgen.rs)

**Euler weight.** For a window word Y and a flow p (moved masks), the equivariant Euler weight is

  e_Y(p) = Π_rungs e(moved, giver ∖ moved),   e(S,T) = Π_{s∈S,t∈T} (y_s − y_t).

This uses `smap::elem_winloc`'s sign convention: Trace = 1/e(S, A∖S).

**Flag-type words.** A word is *flag-type* if it has at most one flow per top coset. Example: `canonical_max(bottom,
top)`, which merges every upright into upright 0 and then splits off uprights 1, 2, ... in order. It is the reduced
word of the window's maximal double coset.

**The kinds.** The localized maps are block-diagonal over top cosets. Each kind is a correspondence:

| kind | map | entry (q, p at the same coset c) | normalized degree |
|---|---|---|---|
| `Pull` | flag-type X → Y | 1 | shift(Y) − shift(X) |
| `Push` | Y → flag-type X | e_X(c)/e_Y(p) | shift(Y) − shift(X) |
| `Rex` | Y₁ → Y₂, both reduced for the window maximum | (1/e_{Y₁}(p)) / Σ_{p′∈c} 1/e_{Y₁}(p′) | 0 |
| `Cap` | [X₀.r, X₀′.r] → [] | 1/e(S, G∖S) where both masks are S | ra − r(g−r) |
| `Cup` | [] → [X₀.r, X₀′.r] | e(S, R) where both masks are S | ra − r(g−r) |
| `Swap` | A·B → B·A, disjoint uprights | 1 at the permuted flow | 0 |

In general, degree = 2·(entry degree in y) + shift(tgt) − shift(src).

**What each generator is:**

| generator | realization | checked by |
|---|---|---|
| E2 Eps / EpsL | Pull to [] | `elem_kinds_are_instances`: `winloc` equal to `elem_winloc` |
| E2 Eta / EtaL | Push from [] | same |
| E3 Iota | Pull | same |
| E3 Trace | Push | same |
| right-facing sideways crossing [Is+t−s] → [Is−s+t] (multiplication) | Pull | tests |
| left-facing sideways crossing (SLL eqs. (91)-(92), ∂Δ) | Push | tests |
| switchback φ (SLL Thm 6.7) | Pull: all ones, f⊗g ↦ f⊗1⊗g | tests |
| switchback ψ | Rex into X_max, the normalized push-forward | tests |
| general cups/caps E^(r)F^(r) ⇄ 1 | Cup / Cap (the Frob composites cw-cap∘ccw-cap and ccw-cup∘cw-cup) | tests |
| far commutation (CORE Lemma 3.13) | Swap | tests |

For blocks (a,b,c), Y = F0.b over (a+b, c) and X = E0.c F0.(b+c). The left-facing crossing's entry is e(C, A∖T),
and right ∘ left is multiplication by e(C, A∖T) on R^I.

**Hand derivation of ψ (the first R1 test).** Take the EKo (108) switchback on two uprights (a, m), with N = a + m:
X = E0.m F0.(N−b) and Y = E0.b F0.a.
- Write V for the top set of upright 0, S for the roots moved by E0.b, and Q = M ∖ S.
- The push-forward weight is ψ_S(V) = e(Q, A∩V)/e(Q, (M∖V)∖Q), up to sign.
- Σ_S ψ_S(V) = 1 is the Lagrange identity: the trace of a product of linear factors of top degree.
- The same entry also comes out as the composite (split-trace) ∘ Assoc ∘ (left crossing) of Frob generators.
- The Rex kind computes it as normalized 1/e_Y, which removes any sign question.

**Trust.** Nothing in leafgen.rs is trusted.
- `interpolate` turns a closed form into the exact PMat of the bimodule map. It interpolates every entry of
  Φ(y) = L_tgt⁻¹ f(y) L_src at its forced degree, then checks it at 2 extra points. It returns Err if f is not
  polynomial.
- `build` adds the PMat as an `Added` template and runs `TmplLib::certify_added` with an independent seed. That checks
  exact degree, exact commutation with every top operator, and the polynomiality certificate.
- `build` then re-checks that the certified template localizes back to the closed form at a fresh point.
- Composites are `Comp`/`Local` DAGs of certified templates and of Assoc nodes that passed `validate_assoc`.
  `Arena::assert_genuine` accepts them.

**Assoc validator** (`validate_assoc`, C_m4c M5). It accepts an Assoc only if all three checks pass:
1. both words have equal EKLP multistep expressions (`multistep`: wall-set sequence, repeats removed, interior points
   of monotone runs dropped);
2. the window flows are in bijection with equal top configurations, with equal shifts;
3. `poly_certificate` passes on `elem_winloc(Assoc)` when both windows have ≤ 500 flows.

It refuses otherwise, and only then records the constant 1 (`certify_assoc`). A switchback pair is correctly refused:
its multistep expressions differ.

## 2. Tests (tests/leafgen.rs)

"Exact" means: in the span of the `sbim::hom_space` basis (via `TmplLib::try_hom_basis_cap`), localized at 2 random
points, whenever `hom_unknowns` ≤ 10⁵.

**Default suite:**

| test | content | result |
|---|---|---|
| `r1_sideways_crossing_111_212` | crossings on blocks (1,1,1), (2,1,2), (1,2,1), (2,2,2), both directions: degree ac; certified; exact (dim 1); right∘left = e(C, A∖T) diagonal | pass |
| `r1_switchback` | EKo switchbacks (a,m,b) ∈ {(1,2,1), (1,3,1), (1,3,2), (2,2,1), (2,3,1), (2,3,2), (2,2,2)}: φ and ψ of degree 0, certified, exact (dim 1); ψ∘φ = id; φ∘ψ idempotent | pass |
| `elem_kinds_are_instances` | Pull/Push equal Eps, Eta, EpsL, EtaL, Iota and Trace (`elem_winloc`), with their degrees; Cap/Cup on (c,0) equal Trace/Iota | pass |
| `cups_and_caps_small` | general cups/caps on (a,b) ≤ (2,3), r = 1..b (E-type) and 1..a (F-type): certified, exact (dim 1); cap∘cup vanishes exactly when its degree is negative | pass |
| `mirror_crossings_and_switchbacks_small` | E-type crossings; mirrored switchbacks F0.b E0.a; a non-reduced word is refused | pass |
| `assoc_validator` | up-up through a zero upright, down-down on (4,0,0): accepted; a switchback pair: refused, nothing recorded | pass |
| `braid_small_cells` | braid composite at (1,1,1), (1,2,1), (2,1,2) (details below) | pass |
| `braid_222` | the same at (2,2,2), composite checks only | pass |
| `far_commutation` | cells (1,1,1,1), (1,2,1,1), (2,1,1,2): certified, exact (dim 1) | pass |

The braid composite checks are:
- the DAG is genuine and has degree 0;
- **1⊗ ↦ 1⊗ at every coset**;
- the top coset (antidiagonal table) has a **1×1 block equal to 1**;
- it equals the direct `Rex` closed form on the whole window at every coset;
- it lies in the exact Hom⁰(A, B) window (dim 1);
- the B→A ∘ A→B top scalar is 1, i.e. the 1×1 form has rank 1.

**Ignored runs** (results in §3 and §4): `certify_all_to_333`, `braid_pieces_large`, `braid_large`,
`far_commutation_3321`.

## 3. Results

### 3.1 Exact comparisons
- Every closed form compared with an exact window was in it.
- Every exact window was 1-dimensional. So each generator is the unique map up to a scalar, and the closed form is
  nonzero (checked).
- Degrees agree with SLL:
  - crossings: ac;
  - switchbacks and braids: 0;
  - caps/cups: ra − r(g−r);
  - Eps/Eta: ab;
  - Iota/Trace: −ab.

### 3.2 certify_added up to blocks (3,3,3)
Logs: `m4c/b1/cert333.log` (sideways crossings) and `m4c/b1/switchbacks_cups.log` (switchbacks, cups and caps).
Each line gives degree, flows, interpolation and certificate time, and monomial count. Every template attempted was
certified; **no closed form failed.**

| family | range | certified | not reached |
|---|---|---|---|
| sideways crossings, F-type and E-type, both directions | blocks (a,b,c) ∈ {1,2,3}³ | **104 of 108** | 4: blocks (3,3,3), F-type 20 × 84 and its E-type mirror (§4) |
| EKo switchbacks φ and ψ on (a,m), a ≤ 3, plus cups/caps on (3,3) | windows up to K = 8, and K = 9 for a = 1, b = 1 | **100** (94 φ/ψ, 6 cup/cap) | K = 9 windows with b ≥ 2, and cups/caps on (3,6) and (6,3) (§4) |
| braid pieces ψ_A and φ_B | cells (3,1,3), (2,3,2), (3,2,3) | 6 | cells (3,3,3), 400 × 84 (§4) |
| far commutation | cells (3,3,2,1) | 1 | — |

**Every template that is not reached is waiting on the same thing.** Each window contains a split of all 9 roots
(F0.r with giver label 9). Building the exact window SBim (`Expr::sbim`), which `certify_added` needs, then calls
`symfun::GrTable::compute(r, 9)`:
- GrTable(6, 9) was still running after 75 min;
- GrTable(3, 9) was still running after 40 min.

These are one-time, in-process builds in frozen code. The closed forms themselves cost microseconds.

Excluding first-use GrTable builds, the measured template costs are:
- certificate ≤ 20 s at ≤ 1 590 monomials, for 2-upright windows up to 200 × 56 flows;
- interpolation of the same order;
- with first-use builds for b = 8 included, one template took up to 450 s.

### 3.3 Braid composites (CORE Lemma 3.12)

| cells | A flows | ψ_A window | φ_B window | M → M′ | composite |
|---|---|---|---|---|---|
| (1,1,1) | 8 | | | Assoc validated | all checks pass, exact (dim 1) |
| (1,2,1) | 18 | 6 → 4 | 4 → 6 | Assoc validated (12 flows) | all checks pass, exact (dim 1) |
| (2,1,2) | 54 | 18 → 10 | 10 → 18 | Assoc validated (30 flows) | all checks pass, exact (dim 1) |
| (2,2,2) | 216 | 36 → 15 | 15 → 36 | Assoc validated (90 flows) | all checks pass (260 s, mostly first-use GrTable builds) |
| (3,1,3) | 320 | 80 → 35, certified 6 s | 35 → 80, certified 9 s | Assoc validated, 140 flows, **245 s** | pieces certified; composite run (`braid_large`) not run this session |
| (2,3,2) | 600 | 60 → 21, certified 0.5 s | 21 → 60, certified 1.3 s | Assoc validated, 210 flows, **1 509 s** | pieces certified |
| (3,2,3) | 2 000 | 200 → 56, certified 243 s | 56 → 200, certified 405 s | 560 flows: **refused** (> 500) | not certifiable under the current rule |
| (3,3,3) | 8 000 | 400 → 84: waiting on GrTable(3, 9) (> 40 min) | 84 → 400: not reached | **1 680 flows: refused** | not certifiable (§0) |

For comparison, the exact braid windows (C_m4c C1) have 111, 946 and 18 019 unknowns at the first three cells, 711 580 at
(2,2,2), and 1.6·10¹⁰ at (3,3,3).

### 3.4 Far commutation
Cells (3,3,2,1): a 60 × 60-flow template, certified. Interpolation took 420 s and the certificate 399 s, with 4 644
monomials (total 819 s). Its exact window is 192 484 unknowns (C1).

## 4. Cost: the bottleneck is monomials, not flows

The cost of `interpolate` and `certify_added` ≈ (#points = max monomials of one entry degree) × (n³ localization
inverse + nm² Vandermonde).
- The monomial count grows with the window's label sum and the fibre-degree spread.
- Examples:
  - the (3,2,3) braid pieces (200 × 56 flows): 1 590 monomials, certificate 32 s;
  - the far commutation on (3,3,2,1) has 4 644 monomials for only 60 flows (820 s);
  - the 140-flow Assoc on (3,1,3) took 245 s.

**One-time cost (the real (3,3,3) blocker).** The first use of a Grassmannian table (`symfun::GrTable::compute`)
happens inside `Expr::sbim` and is cached per process only:
- (r, 8) tables take 3-7 min each;
- (6, 9) took > 75 min and (3, 9) > 40 min, both unfinished when stopped.

Every `Added` template on a 9-root window needs its exact SBim for `certify_added`, so these tables must exist before
any (3,3,3) template is certified.

**Consequence for B2/B4.** The plan's rule "Assoc only when ≤ 500 flows" is not the binding constraint; time is.
Budget validations by `max_monos × n³`, not by flows.

## 5. What B1 provides to later milestones
- `ClosedForm` with kinds Pull, Push, Rex, Cap, Cup and Swap; `winloc(y)` is pure combinatorics plus field operations.
- `interpolate`, `build`, `build_in`: closed form → certified template.
- Constructors:
  - `sideways`, `sideways_abc` (both directions, any 1-rung word on 2 uprights);
  - `switchback`, for any word reduced for the window maximum (EKo (108) via `switchback_eko`; mirrors included);
  - `cup_cap`, general r;
  - `far_commute`;
  - `atomic`;
  - `braid_words` and `braid_dag` (CORE 3.12 route).
- Combinatorics: `canonical_max`, `is_flag_type`, `euler`, `multistep`, `validate_assoc`.
- General rex principle (stronger than the per-move library): **for any two reduced words Y₁, Y₂ of a window's maximal
  coset, Rex(Y₁ → Y₂) is a closed form of degree 0 with 1⊗ ↦ 1⊗.** It is genuine because it equals
  Pull(X→Y₂)∘Push(Y₁→X), and it agrees with the braid composite exactly at (1,1,1), (1,2,1), (2,1,2) and (2,2,2).
  Every EKo braid relation (switchback, atomic braid) lives on such a window. So a rex move is one template whenever its
  window is small enough to certify, and otherwise a composite.

## 6. Open problems

**For B2 (rex planner):**
1. **Large flag-type Assoc (the braid middle step).** M → M′ at (3,3,3) has 1 680 flows on bottom (3,3,3), and with the
   prefix stripped it is the down-down F0.3 F1.3 F0.3 ⇄ F0.6 F1.3 on (9,0,0). Neither validation route is affordable.
   Options:
   - (a) **a structural Assoc certificate for flag-type words**: equal multistep expressions [[I ⊂ W ⊃ J]] imply that
     both BS bimodules are R^I ⊗_{R^W} R^J, generated by 1⊗1, with the identity matching (SLL eq. (88)). This is
     theorem-backed and costs nothing, but it changes the trust rule of M4C_PLAN §3.4 and needs the user's approval;
   - (b) a planner that never needs a full 9-root down-down. Lemma 3.12 always passes through the merged expression,
     so this means cell-aligned chains that avoid the atomic braid on three 3-blocks, if CORE Thm 3.14 allows it per
     target;
   - (c) the plan's fallback R1(b): word-agnostic types.
2. **Budget model.** G-P0 should predict certification cost from max_monos and window flows (§4), not from flows
   alone.
3. Rex windows must be contiguous rung ranges on adjacent uprights (`Local`). Zero-upright placement (C m7) is free
   here: Rex/Pull/Push accept any word on the window, including words through zero uprights (as in the braid windows
   (x+y, 0, z)).

**For B4/B5:**
4. Exit/entry maps W_q ⇄ L_q can use Rex directly when the differing part of the two words sits in one certifiable
   window. Otherwise they need a chain. Each step is degree 0 with 1⊗ ↦ 1⊗, so any composite is admissible
   (Prop 4.19).
5. **GrTable(r, 9)** (> 40-75 min each, in-process cache only) gates every certified template on a 9-root window,
   including the (3,3,3) crossing (20 × 84) and the (3,3,3) braid pieces. Options:
   - persist the tables to disk (a symfun.rs change, outside B1's scope; ask first);
   - or run one long warm-up job per prime that certifies all 9-root generators in one process and stores the PMats
     (an `Added` template cache, keyed by window and closed-form kind).
6. The Cap/Cup generalization to whole-block windows with three or more uprights was not needed and is not
   implemented.

## 7. Reproduce
```
cargo test --release --test leafgen                          # default suite
cargo test --release --test leafgen -- --ignored --nocapture certify_all_to_333
cargo test --release --test leafgen -- --ignored --nocapture braid_pieces_large
cargo test --release --test leafgen -- --ignored --nocapture braid_large
cargo test --release --test leafgen -- --ignored --nocapture far_commutation_3321
```
