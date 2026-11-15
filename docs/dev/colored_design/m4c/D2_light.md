# D2 (light): M4c as targeted cell-local light leaves for closure and pairs

Date 2026-11-11. Read-only design study: nothing was compiled or run in the repository.

**Inputs:**
- L1_eklp.md, L2_alt.md, L3_code.md (this folder);
- origin/colored-m5 docs/dev/colored_design/gl6/M5b_report.md;
- HANDOFF_M4b.md, SPEC_M4_M6a.md §4, SPEC_3C_GL6_rev1.md (D11, D14), SYNTHESIS.md §2.4 and §3;
- code on dev 277a21c.

**Paper keys:**
- SLL = Elias-Ko-Libedinsky-Patimo, arXiv 2401.03053 (JEMS 2025). Its theorem numbers were re-checked against L1's label tables (tools/m4c_tools/eklp_numbers.txt).
- CORE = arXiv 2402.08673.
- EKo = arXiv 2105.08563.
- QR = arXiv 1405.5920.
- W = Williamson, arXiv 1010.1283.

**Tags:** [P] proved in the cited paper; [K] my inference; [E] estimate; [?] open.

## 0. Decision in brief

1. **Route.** Build EKLP light leaves (SLL Def 7.28), but only those the two failing classes need:
   - **closure:** every path of L_p that ends at the identity coset (Prop 8.32 at p = id: no lower terms, no Y_p, no g, no final rex move);
   - **pairs:** the paths of X that end at q with defect ∓s, composed with one top-path leaf of L_q (it replaces the rex move ρ_q);
   - **before both:** a product-rule stage (Stage 0) with no new mathematics. It certifies 5 of the 12 failing closure types (L3 §3.2).

   Double leaves are added only for the B2 validation (and later end0/e_p). No general rex planner is needed (§2.6, §2.9).
2. **Structural finding: padding is unavoidable.** No complete construction stays local inside 3 uprights at (3,3,3):
   - reduced expressions are webs that split blocks into up to 9 cells (EKo §7, Def "SCWebDiag"; CORE Thm 3.11);
   - foam fullness holds only with extra zero uprights, "there exists m′ ≥ m" (QR Prop 3.15, checked in the paper text);
   - the only 3-upright alternative builds rung-step maps L_t·r → L_u between canonical ladders. Those are as large as the failing pairs themselves (§1).

   D2 therefore builds every leaf in a **zero-padded copy** of the word (m² uprights, so 9 for (3,3,3)) and returns it through one new wrapper node `Repad`. This needs a **D11 exception**: an additive change to smap.rs and smapio.rs (§5.2).
3. **Generators.**
   - The existing E2/E3/E4 closed forms.
   - Tiny exact windows: ≤ 3 cell uprights with labels equal to cell sizes, ≤ 10⁴ unknowns, normalized by 1⊗ ↦ 1⊗.
   - Sideways crossings as Added templates from SLL eqs. (91)-(92).
   - Sprinkles taken verbatim from Def 7.33, as `Lin` of `Dot` products on cell uprights.

   No window spans 9 roots, so the 38 Refused and 7 Memory causes disappear.
4. **Certification is unchanged:** `certify_closure`, `certify_pairs`, `certify_full_plan`, genuineness, and the D14 checks in `dechom_postcheck`.
5. **Effort and main risk.**
   - About 9-11 cloud sessions. Stage 0 alone takes 1 session.
   - The main risk is memory. Cell-form intermediate words have 14-26x the flows of canonical ladders: up to 110 592 for a (λ,λ) coset, and about 2·10⁶ after one more rung.
   - A census gate (G-L0, combinatorics only) runs before any generator code (§7, §9).

## 1. Why this route and not the other candidates

| candidate (source) | complete for the failing classes? | local, no big exact window? | verdict |
|---|---|---|---|
| Elias light ladders 1510.06840 (L2 §2.4) | no: decategorified webs, not bimodule maps | n/a | rejected |
| Foam_N / state-space bases: RW Prop 4.13, Hu-Mathas, Tubbenhauer (L2 §2.6) | no: finite N, flag webs only, cyclotomic | n/a | rejected (unit-test oracles only) |
| Regular light leaves after induction, Elias thick idempotents (L2 §2.2) | yes (EW Thm 6.11) | flows ×216 per side; idempotents of w₀(S₆), w₀(S₉) | rejected (size) |
| KL/foam thick spanning set, L2's "S7" | only if U̇(gl_{m′}) → SSBim is full, which by QR Prop 3.15 needs m′ ≥ m, i.e. the same padding; spanning only, not a basis | yes, after padding | worse than D2: same infrastructure, redundant candidate streams (10⁴-10⁸ for large X [L2 E]), and no published thick minimal-diagram set for gl_m (KLMS covers sl₂ only) |
| 3-upright rung-step leaves: X_t := canonical L_t, steps L_t·r → L_u [K] | yes (triangularity, §3) | no: each step is a whole-ladder Hom of the same size as the failing pair | rejected |
| **D2: EKLP leaves, targeted, cell-local after padding** | **yes**: SLL Prop 8.32 and Cor 8.33 (closure), Thm 7.52 (pairs), Thm 7.49 (validation) | **yes**: every generator is on ≤ 3 cell uprights | **chosen** |

**Williamson's reduction** (W Lemma 7.2.1, Thm 7.2.2) remains a valid alternative to the top-path normalization of §2.9: bend X → L_q into closure leaves of X·rev(L_q). It gives the same leaves, plus bending units, so D2 does not use it.

## 2. Mathematical construction, step by step

### 2.0 Dictionary (L1 §1.2; L3 §1.4)

- **Cosets and paths.**
  - A (λ,μ)-double coset is a `schur::Table` t (rows = bottom λ, columns = top μ); ℓ(t_) = `lmin` and ℓ(t̄) = `lmax`.
  - LR(t) = rows refined by cells, RR(t) = columns refined by cells (CORE Props 3.1, 3.3); R_t = R^{LR(t)}.
  - A krhom rung E_i^(r) / F_i^(r) is the EKLP singlestep pair [−s, +s] (shift r(b−r), SLL Def 4.6). At the −s step, column j splits as (r, b−r); at the +s step, one part merges into its neighbour.
  - A path chooses, at each −s step, how every cell of the split column divides (α_i, t_ij − α_i). The +s step is forced (SLL Def 2.3).
- **Path weights** (SLL Def 3.5): defect from `lmin`/`lmax`; poly = π(LR q)/π(LR p) at merges. Σ_paths = `Schur::ladder_char` (Cor 3.7).
- **Y-form.** Y_t := [[λ ⊃ LR(t)]] ∘ atomic(ψ(t)) ∘ [[RR(t) ⊂ μ]], i.e. split the rows into cells, cross, merge into columns. It is a reduced expression of t (CORE Thm 3.11; ψ from Def 3.7 and Thm 3.8).
  - flows(Y_t) = Π_rows multinomial(λ_i; t_i·) · Π_{inversions of ψ} C(a+b, a).
- **Generators** (SLL §5.4): polynomials (Dot), cw/ccw cups and caps (E2 Eps/Eta, E3 Trace/Iota), up-up/down-down crossings (E4 Assoc), sideways crossings (eqs. 84, 91-92).
- **Duality.** D (Def 5.10) flips a diagram.

### 2.1 Stage 0: product rule (closure, no new mathematics; L3 §3.2)

1. Look for a split L_p ≃ A·B (after at most one far commutation E_i F_j, i ≠ j, `homprov::rungs_commute`). A and B are (λ,λ) words of indecomposable types a, b with `c_diag(p) = c_diag(a)·c_diag(b)`; check with `Schur::c_diag` on `ladder_char` of the factors.
2. Get the closure families of A and B from the inner chain (their own small calls: 2-upright (3,3) ladders, already certified at level 4).
3. Compose:
   - outs: `Comp(b_j, Kron(a_i, B))`;
   - ins: `Comp(Kron(a′_i, B), b′_j)`;
   - then the degree-0 commutation isomorphism, an exact 2-rung window: 20 x 20 flows for E0.3/F1.3 at labels (3,6,0).
4. Correctness. Every identity flow of A·B factors through the identity junction exactly when the counts multiply. The identity-coset matrix is then the Kronecker product of two nonsingular matrices, and `certify_closure` verifies this anyway.
5. Targets: 0,0,3/3,0,0/0,3,0 (Eta∘Eta, degree 18), 0,3,0/0,0,3/3,0,0, 1,0,2/2,0,1/0,3,0, 1,2,0/0,0,3/2,1,0, 0,3,0/1,0,2/2,0,1.

### 2.2 Padding (the enabling device) [K]

- **The padded words.**
  - λ_pad puts the real upright i at position m·i and fills the rest with zeros: (3,0,0,3,0,0,3,0,0); (3,0,3,0) for (3,3).
  - Pad(W) replaces a rung between real uprights i and i+1 by a chain of m rungs. The chain moves the same root set across the m−1 zero uprights, and only its first rung has a choice.
  - Flows biject, the shifts agree (r(b−r) + 0 + 0), and the EKLP singlestep word is the same, since zero labels are invisible (EKo §7, "the object 0 is the monoidal identity").
  - So BS(Pad W) = BS(W) with the same localization.
- **What the padding hosts.** Row i's cells live on uprights mi..mi+m−1 and column j's cells on mj..mj+m−1. So Y-forms and every Grassmannian window (SLL Thm 2.6) fit on contiguous uprights.
- **The wrapper.** All M4c maps are built inside the Arena on padded expressions. Each output is then wrapped as `Repad{inner}`: W → W′ (unpadded), where inner: Pad(W) → Pad(W′).

### 2.3 Paths, defects and sprinkles (SLL §3, §7.6) [P]

- **Enumeration.** Enumerate the sprinkled paths (t, y) of a padded word; it is `ladder_char` with the terms kept apart (`split_col` :913, `merge_cols` :931).
  - Keep only the paths whose terminus is the identity (closure) or q with the target degree (pairs). Prune any prefix that cannot reach the terminus; the reachability table comes from `ladder_char` of the suffix.
- **Sprinkles** at merge steps: f_j = ∂_{L_j} ∂_{y_j}(P_{R_{j+1}}) (Def 7.33), with P_R a staircase monomial (Notation 7.31).
  - Compute them symbolically (divided differences in ≤ 9 variables, mpoly.rs).
  - Expand them in products of Schur polynomials of the L_j blocks. In padded form every block is one upright, so this is a `Lin` of `Comp` of `Dot{level, upright, nu}`.
  - Descending steps carry f = 1. The count is Σ v^deg = c_terminus (Lemma 8.7).

### 2.4 Generator library (all ≤ 3 adjacent padded uprights)

| generator | realization | certification |
|---|---|---|
| polynomial | `Dot` (E1) | none needed |
| cw cap/cup, full block (merge-split of whole blocks) | `Elem` Eps/Eta/EpsL/EtaL | `Arena::certify_elem` per padded label pair (A7-style check against the 2-upright `hom_space`) |
| ccw cap/cup (digon) | `Elem` Trace/Iota (+L) | same |
| up-up/down-down crossing, block slides across zero uprights | `Elem` Assoc (E4) | `certify_assoc` per window pair |
| sideways crossing [Is+t−s] ⇄ [Is−s+t] | `Local` on an `add_template` built from eqs. (91)-(92): right-facing = multiplication, left-facing = f ↦ ∂Δ^I_{Ist}(f) | `TmplLib::certify_added` (degree, exact commutation, `poly_certificate`) |
| atomic braid move a_i a_{i+1} a_i ↔ a_{i+1} a_i a_{i+1} on 3 cells, far commutation on 4 cells (CORE Lemmas 3.12, 3.13) | degree-0 exact window (`try_hom_basis_cap`), normalized by top scalar 1 | genuine (HomSpace) |

A leaf may use any degree-0 map that sends 1⊗ ↦ 1⊗ wherever SLL uses a rex move: Prop 8.26 uses only this property (Thm 6.11; L1 §5.2). By Prop 4.19 such a map equals the rex move modulo Hom_{<p}.

### 2.5 Single-step leaves (SLL §7.1-7.4)

For the step [t_k, t_{k+1}] of a path:

1. **Grassmannian reduction** (Thm 2.6, eq. 11): the coset pair splits as z, m, n. Here z is the (λ, L)-coset with L = RR(q), and the pair m ⊂ n lies in W_{Js}, i.e. the two active columns.
2. **Expressions.** X_m (Notation 7.1: splits into LR(m) cells, the core of m, merges) and X_n (merges, ordered by Lemma 7.2). Z = Y-form of z.
3. **Ascending step (+s):** SSLLP = id_Z ⊗ ELLP([m,n], f) (Def 7.14, Def 7.23), placed as a `Local`/`Elem` composite at `at = |Z|`, `first` = the first cell upright of the active columns.
   - ELL is the **sinister diagram**, compiled by the right-to-left algorithm of Lemma 8.2 (Def 7.3): ccw caps (E3 Trace) and left-facing crossings (2.4), without polynomials. Its degree is def[m,n] (Lemma 7.4).
4. **Descending step (−s):** ELL([n,m]) := (id ⊗ cw cap) ∘ (D ELL([m,n]) ⊗ id_[−s]) (Def 7.16, eq. 133).
5. **Step fusion (optimization).** When the active window has ≤ 2 000 flows, compose the step's generators exactly with `tmpl::tensor_exact` or a PMat product. Register the result as one Added template (`certify_added`). Each step then interns 2 whole expressions instead of about 10.

### 2.6 Rex moves inside the recursion (cell-local only)

SLL Def 7.28 needs a rex move from the previous target Z′∘X_{n′} to the next source Z∘X_m. Both are split-cross-merge webs on cells of one coset t_k [K].

- **Normalizer.**
  - First apply up-up/down-down: E4 Assoc merges and splits in any order.
  - Then bring the cell order to ψ(t_k) by bubble sort on reduced words of ψ (CORE Thm 3.14, atomic Matsumoto). Each move is an atomic braid template (2.4).
  - The low road EKo Prop 4.28 gives the [[λ ⊃ LR]] / [[RR ⊂ J]] parts.
- **Fallback.** If a pair of forms is not related by these moves (a switchback on non-cell-aligned blocks, EKo eq. 108), use one exact degree-0 window on the differing cells (≤ 6 cell uprights). Take any element with top scalar 1 (Prop 4.19).

### 2.7 Light leaves LL_W(t): Pad(W) → Y_{term t} (SLL Def 7.28) and flips

- **Recursion.** LL_{≤k+1} = SSLLP_{k+1} ∘ rex_k ∘ `Kron`(LL_{≤k}, step_{k+1}).
  - This is bottom-up: the prefix grows through `Kron` and windows act on top, which matches L3 §1.2 workaround 1.
  - Memoize per path prefix, so the leaves of one word share their DAG.
- **Flip D** [K]:
  - reverse `Comp`; swap Eta↔Eps, Iota↔Trace, Assoc↔Assoc⁻¹; a template becomes the exact window in the reverse direction at the same degree, normalized by top scalar 1; `Dot` is unchanged.
  - Every leaf is a pure composite, with no `Lin` except sprinkle sums. A flip is therefore D(leaf) up to one nonzero scalar per leaf, which changes no rank that certification tests.
- **Top-path leaf.** LL^top_{L_q} is the leaf of the top path of Pad(L_q). It is degree 0 and sends 1⊗ ↦ 1⊗ (Lemmas 8.9, 8.10, 8.14, Thm 6.11), because L_q is reduced: char L_q = b_q gives ℓ(L_q) = ℓ(q) by Thm 4.13 (L1 §1.2).

### 2.8 Closure output (`HomProvider::closure`)

- **outs:** for every sprinkled path t of Pad(L_p) with terminus diag(λ), output `Repad{ LL(t) }`: L_p → R (Y_id = [λ], the empty word).
- **ins:** `Repad{ D LL(t) }`.
- **ep:** with e_p, wrap as in the M4b contract (`Comp(e_p, g)` / `Comp(g, e_p)`).
- **Count and degrees:** n = `Schur::c_diag(p)(1)` maps whose degree multiset is c_diag(b_p) (Lemma 8.7). This is exactly the `want_degs` of `certify_closure`.

### 2.9 Pair output (`HomProvider::pairs`)

- **π candidates:** π_t = `Repad{ D(LL^top_{L_q}) ∘ LL_X(t) }` for sprinkled paths t of X with terminus q and deg = −s. Only paths with g = 1 can contribute (§3).
- **ι candidates:** ι_{t′} = `Repad{ D(LL_X(t′)) ∘ LL^top_{L_q} }` for deg = +s.
- The number of candidates is [v^{∓s}] c_q(char X): 1 for the 7 rex-class pairs, and a handful for the others.
- No rex planner between X and L_q is needed: the top-path leaf plays ρ_q/σ_q.
- `certify_pairs` picks a rank-mult minor of the form at rep(q).

### 2.10 Double leaves (validation and later `full` / `end0`; SLL Def 7.44)

- **Family:** DLL = D(LL_tgt(t′)) ∘ g ∘ LL_src(t) for each stratum c, where:
  - t and t′ are coterminal at c;
  - g runs over the Schur-product basis of R^{LR(c)} over R^λ, as `Dot`s on the cell uprights just above [[λ ⊃ LR(c)]] in Y_c;
  - each DLL is multiplied by left R^λ monomials (`Dot` at level 0).
- Selection by `certify_full_plan`.

## 3. Why the result is complete for the failing classes

- **Assumptions** (SLL §4.1) [P]:
  - F_p is a field, with p = 2³¹−1 and 4294967291, both > 9!;
  - the gl₉ realization is faithful, reflection faithful, balanced and Demazure surjective (∂_s x_i = 1).
  - `hom_gd` already assumes Thm 4.16 / Cor 4.17.
- **Closure.**
  - diag(λ) is Bruhat-minimal, so Hom_{≮id} = Hom, Y_id = [λ] and R_id = R^λ.
  - Prop 8.32 (Thm 8.31 independence, Lemma 8.7 count, graded Nakayama) says that {LL(t,y) : term = id} is an R^λ-basis of Hom(BS L_p, R^λ). Cor 8.33 gives the flipped basis of Hom(R^λ, BS L_p).
  - SLL allows any choice of reduced expressions, rex moves, P_L and enumerations (§7.10, Thm 7.51). Prop 8.26 uses only degree 0 and 1⊗ ↦ 1⊗ (Thm 6.11), so the substitutes of 2.4 and 2.6 are admissible [K].
  - BS(L_p) ≅ B_p because `closure_precheck` asserts char = b_p.
  - A map to R is determined by its identity-coset block, since R has one flow and localization is faithful. So a basis gives n Frac-independent columns, which is exactly what `certify_closure` tests.
- **Pairs.**
  - Prop 8.32 and Cor 8.33 at q make the leaves an R^{LR(q)}-basis of Hom_{≮q}(X, Y_q) and of Hom_{≮q}(Y_q, X).
  - SLL is fibered cellular with cell algebras R_q (Thm 7.52; cellular theory of Elias-Lauda Def 2.17).
  - Therefore the multiplicity of B_q⟨s⟩ in X equals the rank of the degree-0 local intersection form on the leaves of degree ∓s [K, standard consequence; SLL intro §1.1].
  - Leaves with positive-degree R_q factors pair into the maximal ideal and drop out.
  - D(LL^top) and LL^top are 1x1 nonzero scalars u, u′ at rep(q), since both L_q and Y_q have a single flow there. So form(π_t, ι_{t′}) = u·u′·(local form), with the same rank = mult = [v^s] m_q (the mult checked by `pairs_precheck`).
- **Full** (validation only): Thm 7.49.
- **Why the failure modes cannot recur.**
  - Every Refused and Memory failure came from 3-upright windows spanning all 9 roots (L3 §2.4).
  - In D2 the only exact calls are generator windows on ≤ 3 cell uprights, plus the §2.6 fallback.
  - If an inference [K] above were wrong, certification fails visibly. It never accepts a wrong family.

## 4. Certification (unchanged; all `pub` in homprov.rs)

- **Closure:** `certify_closure` (:945) checks the degree multiset, the n x n identity-coset rank at two points, and ends; then `assert_genuine` (:550). sindec `closure_gens` (:1007) re-certifies after every layer-A load.
- **Pairs:** `pairs_precheck` (:880) and `certify_pairs` (:790): ends and degrees, `form_at` (:726) at rep(q) at two points, a rank-mult minor, point independence, genuineness. Then the sindec checks: `dec_one` (:610), its own form check, `orthogonalize`, `check_pi_iota`, and `dechom_postcheck` with the D14 degree scaling.
- **Full:** `certify_full_plan` (:635): count = `hom_gd`, stacked rank, fresh residual, orbit fallback.
- **New leaves of genuineness:**
  - `Repad` is genuine iff its inner map is.
  - Elem constants for padded label pairs are recorded through `certify_elem` after an A7-style comparison.
  - Added sideways-crossing templates pass `certify_added`.
  - Atomic braid templates are HomSpace, hence genuine.

  Nothing else is trusted.

## 5. Code changes

### 5.1 New files (D11-compatible)

- **src/dcoset.rs** (combinatorics only, no field):
  - `Perm`; `table_min_max(&Table) -> (Perm, Perm)`; `lr_cells`, `rr_cells`; `psi(&Table) -> Perm` (CORE Def 3.7);
  - `grass_reduce(p, q, step) -> (Table /*z*/, Grass{m, n, L})` (Thm 2.6);
  - `yform(&Table, &Pad) -> Vec<Rung>`, `pad_word(&[u8], &[Rung]) -> Expr`, `xm_xn(&Grass, &Pad)`;
  - `paths(&mut Schur, &Expr, Terminus, deg_filter) -> Vec<SPath{cosets, choices, defect, sprinkle_idx}>`;
  - `census(...)` (gate G-L0);
  - tests in tests/dcoset.rs.
- **src/leafgen.rs:** `GenLib<F>` with `elem(kind, labels)`, `sideways(a,b,c,dir)` (Added + `certify_added`), `atomic_braid(x,y,z)`, `far_comm`, `rex_window_fallback`, and `flip_tmpl`. Templates are cached in the shared `TmplLib` (`try_hom_basis_cap`, `add_template`).
- **src/leaves.rs:** `pub struct Leaves<F> { pad, gens: GenLib<F>, inner: Option<Box<dyn HomProvider<F>>>, memo }` implementing `homprov::HomProvider<F>` with `name() = "Leaves"`:
  - `closure`: Stage 0 `product_closure`, else §2.8;
  - `pairs`: §2.9;
  - `full` and `end0`: §2.10;
  - helpers `ll`, `ssll`, `ell_sinister` (Lemma 8.2), `rex_cells`, `sprinkle`, `flip`, `ll_top`.
  - Re-implemented private helpers (L3 §1.1): `cdiag_degrees`, `identity_coset`, `rep_rank`, the closure `Selector`.
  - It follows the `ar.mark()` / `ar.commit(mark, keep)` discipline of `run_call`.
- **src/bin/m4c.rs:**
  - `m4c census|g1a|b2|list`. g1.rs is frozen, so the G1(a) and B2 job lists are rebuilt through `schur` / `sindec::canonical_ladder` (:174).
  - The gates run here with a Leaves-only provider.

### 5.2 Edited files

- **D11 exception, additive** (ask the user first; a review of the diff is recommended):
  - **smap.rs:**
    - `Arena::with_pad(lam, pad)`; `Arena::expr` (:581) also accepts bottom == pad;
    - `MapNode::Repad { inner }`;
    - arms in `Arena::add` (pad(src) = src(inner), pad(tgt) = tgt(inner), same degree), `children`, `assert_genuine` (:855, recurse);
    - `Eval::block` (:1092): translate c with `CosetSpace(λ_pad, top_pad).rank(pad(conf))`, evaluate inner at c′, then reorder rows and columns through the flow bijection (`FlowIndex::find` on padded masks). About 150 lines.
  - **smapio.rs:** a record for `Repad` and for padded expressions; a "pad=" field in the stamp.
- **Not frozen:**
  - **sindec.rs:**
    - `provider_for_stamped` (:421) appends `Leaves` for max label ≥ 2, with stamp `chain=S0-2,Exact(l3,600s),S3-6,Leaves(v1,pad)`;
    - a new `provider_m4c(lam)` puts Leaves first, for the acceptance list.
  - **libbuild.rs:** a `--leaves-first` flag.
- **Untouched:** homprov.rs, tmpl.rs, sbim.rs, linsys.rs, g1.rs.

## 6. Validation plan (cloud; the laptop runs nothing)

1. **G-L0 census** (combinatorics, seconds). For every M5b failure, list the needed paths, the visited cosets, and the flows and FlowIndex bytes of every expression the construction interns. Decision rule:
   - if the worst job exceeds 3 GB, do §9 R2(c)-(d) before Stage 2.
2. **G-L1:**
   - Σ paths = `ladder_char` (Cor 3.7) for all canonical ladders of (3,3), (2,1,2), (2,2,2), (3,3,3) and their one-rung extensions, padded and unpadded;
   - sprinkled counts = `c_diag` and `term_mult` (Lemma 8.7);
   - ℓ(L_q) = ℓ(q) for all kept types.
3. **G-L2:** SLL examples against `sbim::hom_space`:
   - Ex 7.39: [s,∅,s] → [s], two leaves;
   - Ex 7.40 and 7.41: S₃ with I = (2,1), J = (1,2), 4 and 5 leaves;
   - Ex 7.46: 10 double leaves;
   - ELL_sinister = the exact Grassmannian window up to scalar, on all windows with ≤ 2 000 flows.
4. **G-R:** for random words of (3,3) and (2,2,2), Repad(Pad f) = f blockwise at random cosets and points (f from layer A); `assert_genuine`; smapio round trip.
5. **Backend alone (Leaves only):**
   - G1(a): (3,3) 72 pair targets + 8 closure; (2,2,2) 358 + 42; **480/480**;
   - B2 through `full`: (2,2) 45/45, **(2,1,2) 605/605** (Search alone reaches 581), then (2,2,2) 2205 (the Chain reaches 2064).
6. **(3,3,3) regression:**
   - level 4 (47/47) and level 5 (177 jobs, including the 6 M5b level-5 failures) through `provider_m4c`, both primes;
   - compare against the layer-A maps: same kept multiset, π∘ι = id.
7. **M5b level-6 list**, both primes:
   - 12 closure types, 11 Exhausted pairs, 7 Memory pairs, then the 38 Refused;
   - `libbuild replay` with a refusing provider;
   - the D14 scaling check is part of `dechom_postcheck`.

## 7. Cost and complexity for the (3,3,3) failing targets

flows(Y_t) is computed exactly by d2tools/cellflows.py, which uses the formula of §2.0. Everything else is [E].

| class (jobs) | X or L_p flows (canonical) | leaves per call | largest Y-form flows | exact calls | estimated time per job |
|---|---|---|---|---|---|
| closure, product (5) | 400-600 | n = 1 or 3 per side | none | 1 commutation window (20x20) | < 5 s |
| closure 2 700 (2) | 2 700 | 18 | Y_p 26 244 | generator windows only (≤ 10⁴ unknowns) | 10-60 s |
| closure 2 880 (2) | 2 880 | 18 | Y_p 46 656 | same | 10-60 s |
| closure 7 200 (2) | 7 200 | 27 | Y_p 104 976 | same | 30-200 s |
| closure 1,1,1/1,1,1/1,1,1 (1) | 4 320 | 54 | Y_p 110 592 | same | 60-300 s |
| pairs, rex class (7) | 240-7 200 | 1 π + 1 ι, plus LL^top of L_q | Y_q ≤ 104 976 | same | 10-120 s |
| pairs, split class (4) | 1 080-14 400 | Σ mult 2-6 per target | ≤ 104 976 | same | 30-300 s |
| Memory (7) | 240-7 200 | per target | ≤ 104 976 | same | 10-120 s |
| Refused (38) | 480-21 600 | per target | ≤ 104 976; ≤ 944 784 if the decomposable type 0,1,2/1,1,1/2,1,0 is visited | same | 30-600 s |

- **Comparison with today.** The same targets now need whole-ladder or 9-root windows of 0.3-7.5 M unknowns (closure) and 0.5-919 M unknowns (pairs; L3 §3).
- **Generator windows.** Labels are cell sizes ≤ 3, or block sizes ≤ 9 on ≤ 3 uprights and ≤ 2 rungs. That means ≤ 1 680 flows and ≤ 10⁴ unknowns, in total O(10²) distinct templates per λ.
- **Memory.**
  - FlowIndex = flows · rungs · 8 B (Mask = u64, coset.rs:12).
  - A Y-form word plus one padded rung has about 35 rungs: ≈ 31 MB at 10⁵ flows. After a 3-root rung (×20), about 2·10⁶ flows ≈ 0.6 GB for one expression.
  - The number of distinct interned expressions is about 2 per (level, coset) with step fusion, about 10 without. That gives 0.3-3 GB per job for the 7 200-flow and 4 320-flow targets. This must be measured by G-L0.

## 8. Effort (cloud sessions of about 4-5 h)

| # | content | exit criterion |
|---|---|---|
| 1 | Stage 0 product closure in leaves.rs (inner chain); G-L0 census tool | 5 closure types certified (both primes); census table |
| 2 | dcoset.rs: tables ↔ permutations, cells, ψ, Grassmannian triple, Y-form, padding, path and sprinkle enumerator | G-L1 |
| 3 | smap.rs/smapio.rs Repad + padded Arena (D11 exception) | G-R |
| 4 | leafgen.rs: Elem constants at padded labels, sideways Added templates, atomic braid and far-commutation templates, flips | unit tests against `hom_space` |
| 5 | ELL compiler (Lemma 8.2), SSLL, sprinkles, LL recursion, cell rex normalizer and fallback | G-L2 |
| 6 | closure through Leaves: (3,3), (2,2,2), then the 12 (3,3,3) types | G1(a) closure 50/50; 12/12 |
| 7 | pairs: LL^top, flips, π/ι | G1(a) 480/480 alone; 11 Exhausted + 7 Memory |
| 8 | double leaves (`full`); 38 Refused; production wiring and stamps; level-6 rerun of both primes and replay | B2 605/605; M5b list closed |
| 9-11 | buffer: memory engineering (§9 R2), fallbacks, report | |

**Decision points:**
- after session 1: does the census fit? Is the D11 exception granted?
- after session 5: if G-L2 fails, stop and re-derive the type-A ELL before going on.

Stage 0 is useful even if everything after it is deferred.

## 9. Main risks and fallbacks

- **R1. D11 governance.** smap.rs and smapio.rs are frozen.
  - Fallback A: build the (3,3,3) production arena padded from the start. This is invasive: all keys and Tids change, and layer A must be rebuilt.
  - Fallback B: Stage 0 only, and leave the remaining 63 jobs to the parallel engineering route.
- **R2. Memory of Y-form intermediates** (14-26x canonical; up to about 2·10⁶ flows per expression). Mitigations, in order:
  - (a) choose Z with the fewest flows (Z is a free choice, §7.10);
  - (b) step fusion into Added templates (§2.5);
  - (c) `Arena::drop_index`, and evaluate leaf by leaf before the family certification;
  - (d) a lazy per-coset FlowIndex in sexpr.rs (not frozen), using append-only storage behind the same methods. Evaluation at the identity coset or rep(q) touches only the cosets reachable from that terminus.
  - Run the worst jobs alone at 6-8 GB.
- **R3. Cell rex normalizer incomplete.** It fails if a needed move is a switchback on non-cell-aligned blocks (EKo eq. 108).
  - Fallback: one exact degree-0 window on the differing cells, with top scalar 1 (Prop 4.19). This is measured by the census; worst case 6 singleton cells, a regular S₆ piece of about 10⁴-10⁵ unknowns.
- **R4. ELL compiler or type-A specialization bugs** (orientation: krhom bottom = SLL left; E vs F sign; sinister order).
  - Caught by G-L2 against exact windows. A wrong leaf fails certification visibly.
  - Fallback: ELL by exact Grassmannian window when it has ≤ 2 000 flows.
- **R5. Scalars and normalization of flipped templates.**
  - Ranks are invariant under nonzero per-leaf scalars.
  - Where a generator space has dimension > 1, use the closed forms (91)-(92), never "any element".
- **R6. The form-rank statement for pairs** rests on SLL Thm 7.52 plus cellular theory [K].
  - If it failed, `certify_pairs` would fail with a rank shortfall.
  - Fallback: add D LL^top ∘ g ∘ LL candidates with degree-0 g from R_q.
- **R7. Effort overrun** of sessions 4-5, the two new algorithmic pieces. The SYNTHESIS estimate for full M4c was 10-20 days. D2 drops double leaves from the critical path and drops any general rex planner, but keeps ELL.

## 10. Open checks before coding

1. Read SLL §8.1 (Lemma 8.2) and EKo §6.1 in full, and write the type-A cell-form ELL for one Grassmannian pair from Ex 2.5/2.10 by hand.
2. Confirm that `tmpl::tensor_exact` (:618) can compose templates exactly for step fusion; otherwise drop the fusion.
3. Check that `Window::new` accepts zero labels in every position (padding chains), and that `elem_spec` handles them.
4. Run G-L0 for the 68 failures first (session 1).
