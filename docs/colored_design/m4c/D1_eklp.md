# D1: M4c "Leaves" backend, singular light leaves after EKLP (design)

Date 2026-11-11. Read-only design; nothing was compiled or run. Inputs: L1_eklp.md (paper extraction), L2_alt.md
(alternatives), L3_code.md (code map, failure analysis), origin/colored-m5 M5b_report.md / M5b_measure.md, dev 277a21c.
Paper keys as in L1: SLL = EKLP arXiv 2401.03053 (JEMS 2025), CORE = arXiv 2402.08673v2, EKo = arXiv 2105.08563,
W = Williamson IMRN 2011. Numbers re-checked against m4c/tools/m4c_tools/{eklp,core,eko}_numbers.txt.
Tags: [P] proved in the cited paper; [K] my inference, checked by a gate below; [E] estimate.

## 0. Decisions in brief

1. **Leaves are built by EKLP's recursion (SLL Def 7.28), but in localized form.** Junction expressions are *chosen*
   reduced words. Single-step leaves are *selected* from small exact `hom_space` bases by a rank test at one coset;
   they are not compiled from sinister diagrams. Lemma A (§1.3) turns EKLP's triangularity (Prop 8.26) into a numerical
   criterion that krhom can check. Prop 8.32 and Cor 8.33 guarantee that the selection never runs dry (§1.4).
   **No generators-and-relations presentation is used.** Every map is an evaluated bimodule map: an exact template, an
   E2/E3/E4 closed form, or a Dot.
2. **Closure first (stage L1).** It needs no rex move, no padding and no edit of a frozen file. It covers all 12
   Exhausted closure types. The 5 product types of L3 §3.2 are the special case where the path passes the identity at
   the product junction.
3. **Pairs (stage L2) need locality.** A pair X = L_p·ρ → L_q has an irreducible "top-prefix" step from p via ρ to q,
   which is the original Hom problem again. EKLP makes this step local by the Grassmannian reduction (SLL Thm 2.6).
   In krhom that needs junction words that pass through finer parabolics, i.e. **more uprights than λ has**.
   - Design: a right-padded working chart (λ⁺ = λ·0^z, z ≤ 6 for (3,3,3)) and one additive DAG node `Pad`.
   - Rex moves are chains of local braid-move windows: EKo eq. (108) switchbacks, CORE Lemmas 3.12/3.13, E4 Assoc.
4. **Certification is unchanged:** `certify_closure`, `certify_pairs`, `certify_full_plan`, `assert_genuine`. Leaves
   adds only leaf-level items: hom_space templates are genuine by construction, plus Elem/Assoc constants and the
   `Pad` identity.
5. **Failure semantics.** Leaves returns `Refused` when a window is not admitted and `Time`/`Memory` on budget.
   `Exhausted` occurs only if a rank test fails although every window was computed. By §1.4 that cannot happen with
   correct code, so it is a red-flag bug report, not a completeness signal.

## 1. Mathematics

### 1.1 Dictionary (details in L1 §1.2)
- A (λ,μ)-coset is a `schur::Table` t. ℓ(t_) = `lmin` and ℓ(t̄) = `lmax`. LR(t) and RR(t) are the cells read by rows
  and by columns (CORE Props 3.1, 3.3). R_t = R^{LR(t)}.
- A krhom rung on uprights (i,i+1) is two EKLP single steps: a split of the giver column (descending, any child) and
  a merge into the receiver (ascending, forced) (L1 §1.2).
- A *configuration* c is a top coset (`coset::Conf`, ranked by `CosetSpace::rank`). Its type is the table of root
  counts (λ-block × upright). rep(t) = `CosetSpace::rep`.
- **Reduced junction word** W_t of a (λ,μ)-type t: a rung word on λ with c_t(char W_t)(1) = 1 and every table of
  char W_t ≤ t (SUB Thm 2.16: Bruhat order = support). Then W_t has exactly one flow at every configuration of type t
  ("top flow") and none at configurations of type ≰ t. Any word with char = b_t qualifies. Checked by `ladder_char`.
- Normalization: krhom degrees are the normalized degrees of `MapInfo::deg`. The EKLP defect equals the krhom degree
  [K]. Gate G-L1 checks this.

### 1.2 Step generating function
For a junction t' at level μ', the rung ρ (giver column g, split (r, b−r)) and a target t at level μ:

  G_{t'→t}(v) = Σ_{u: split_col(h_{t'}) ∋ h_u, merge_cols(h_u) ∋ h_t} v^{lmax t' − lmax u} · [π(LR t)/π(LR u)](v)

This comes from schur.rs `split_col` (:913, SLL eq. 44) and `merge_cols` (:931, eq. 45) applied to single h-terms;
SLL Def 3.5 gives the defect and the poly factor. G(1) is the number of top-to-top flows from t' to t.
Deodhar (SLL Cor 3.7): Σ_{t'} c_{t'}(char X') · G_{t'→t} = c_t(char X'ρ). Gate G-L1 checks this against
`ladder_char` (:997) and `hom_gd` (:1289).

### 1.3 Lemma A (localized light leaves)
**Setting.** X = X'·ρ is a word. Fix reduced junction words W_{t'} for the types at the top level of X', and W_t at
the top level of X.

**Hypotheses.**
- (i) For every t', a family LL_{t'} ⊂ Hom(X', W_{t'}) whose block at rep(t') is square and invertible. The block
  is `Eval::block`: rows = maps, columns = flows of X' at rep(t').
- (ii) For every pair (t', t), a family S_{t'→t} ⊂ Hom(W_{t'}·ρ, W_t) whose *top-to-top matrix* at c = rep(t) is
  square and invertible. Its rows are the maps. Its columns are the configurations c' of type t' that step to c. The
  entry is the coefficient of block(S, c) on the flow (top flow of W_{t'} at c') ⊗ (the step c' → c).

**Claim.** LL_t := { S ∘ (ℓ ⊗ id_ρ) : t', ℓ ∈ LL_{t'}, S ∈ S_{t'→t} } has a square invertible block at rep(t). The
degrees add.

**Proof sketch.**
- `Kron` evaluates ℓ ⊗ id_ρ at the prefix configuration of each flow (smap.rs `kron_block` :1145). So the entry of
  row (ℓ, S) at a flow f with prefix configuration c' is Σ_g S_{(g,step)} · ℓ_{g,f'}(c'), summed over the flows g of
  W_{t'} at c'.
- This vanishes unless type(c') ≤ t', because W_{t'} has no flows elsewhere. Ordering rows by t' and columns by
  type(c') along a linear extension of the Bruhat order makes the matrix block-triangular.
- On the diagonal, W_{t'} has the single flow g₀(c'), so the entry is T_S(c')·ℓ_{f'}(c').
- Each block ℓ(c') is a W_λ-translate of the block at rep(t'), hence invertible, because genuine maps are
  R^λ-bilinear. After multiplying by blockdiag ℓ(c')^{-1}, the diagonal block becomes T ⊗ Id, which is invertible.
- Squareness is the Deodhar count of §1.2. □

**Base case.** X' = ∅ and LL_id = {Id}. The **dual statement** (columns, for maps W_t → X) holds verbatim. It builds
ι and "ins" from S^∨ ⊂ Hom(W_t, W_{t'}·ρ).

### 1.4 Existence of the step families (this is the completeness argument)
Apply SLL Prop 8.32 [P] to the expression W_{t'}·ρ with target t: rex∘LL is an R^{LR(t)}-basis of Hom_{≮t}(W_{t'}ρ, W_t).
- **Localization at rep(t) is injective on Hom_{≮t}.** A degree-0 map killing the top is lower (SLL Prop 4.19 [P]).
  In general degrees this is W Thm 7.3.5 (stalk at t) [P per L2 §2.3]. Hence the R-basis has an invertible
  rep(t)-block.
- **The top-path subfamily is the one we need.** Apply the triangularity of §1.3 to X' = W_{t'}, whose own leaves
  have the identity as their t'-diagonal. The leaves whose prefix path is the top path of W_{t'} then have an
  invertible top-to-top block. Their degrees are exactly G_{t'→t} (Lemma 8.7 [P]). Leaves through lower prefixes
  have zero top-to-top coordinates.
- **Greedy selection succeeds.** The top-to-top image of ⊕_{d ∈ supp G} Hom^d(W_{t'}ρ, W_t) at a point is spanned by
  that subfamily: polynomial multiples only rescale at a point. So greedy selection by ascending degree on a k-basis
  of each Hom^d (the `Selector` logic, homprov.rs :1396) finds a family with degree multiset G and an invertible
  top-to-top block.
- **Dual side.** Cor 8.33 [P] gives the flipped families.
- **Consequence.** For *any* choice of reduced junction words, the candidate stream contains what Lemma A needs.
  Assumptions (SLL §4.1): a field, and the gl_K realization, which is faithful, reflection-faithful and Demazure
  surjective. krhom works over F_p with p > 9! (L1 §4).

### 1.5 Closure (outs: L_p → R; ins: R → L_p)
- The target is the identity (λ,λ)-coset id. Its junction word is the empty word, and there are no lower terms
  (SLL §7.7 and its footnote: id is Bruhat-minimal).
- LL_id(L_p) has n = c_id(b_p)(1) = `c_diag`(p)(1) maps. Their degree multiset is c_diag(b_p) (Lemma 8.7 + §1.2), and
  their identity-coset block is invertible (Lemma A).
- That is exactly `certify_closure` (homprov.rs :945). The determinant argument there makes it an R^λ-basis.
- **Short junction words [K].** A junction type t_k on an identity path satisfies t_k ≤ top(ρ_1…ρ_k), and it must
  reach id through ρ_{k+1}…ρ_m.
  - So for m ≤ 6, the junction words needed have ≤ ⌈m/2⌉+1 rungs, and step calls are between words of ≤ 4 rungs.
  - Gate G-L0 verifies this per type. The fallback is the canonical ladder of t.
- **No top-prefix problem in practice:** a top prefix type can lie on an identity path (e.g. E^(b) then F^(b)),
  but then its junction word is itself short. So every closure step is a call between short words.

### 1.6 Pairs (π: X → L_q of degree −s, ι: L_q → X of degree s)
- Lemma A gives leaves of X ending at q with an invertible rep(q)-block, i.e. an R_q-basis of Hom_{≮q}(X, ·)
  (Prop 8.32). It keeps only paths of total degree −s; the prefix degree windows come from G.
- The degree-(−s) members modulo R_q⁺·(lower) are the generators. Their local intersection form with the degree-s ι
  has rank [v^s]m_q: SLL Thm 7.52 (fibered cellular) plus graded Nakayama [K, singular analogue of EW].
- `certify_pairs` (:790) checks the rank directly, so a shortfall is a visible failure, never a wrong answer.
- **Obstruction.** The step through the top prefix p of X (S_{p→q} ⊂ Hom(L_p·ρ, W_q)) is the original whole-ladder
  problem. Every final step also targets a long word for q. §1.7-1.8 make these local.

### 1.7 Rung-local Grassmannian reduction [P for single steps, K for the rung]
- Let K be the merged block of the two columns of ρ, and Q the (λ, μ'_K)-coset containing t'. The rung keeps every
  row total in K, so t ⊂ Q as well.
- SLL Thm 2.6 [P] gives z with z_ = Q_ and Grassmannian pieces in W_K: t' = z.m', t = z.m'', with the **same z**.
  The right parabolic of z is L_K⁺ = (row totals c_0, c_1, c_2 inside K; the other column). Undoing the other-column
  refinement is a commuting merge.
- Lemma 3.9 [P]: the defect and the poly of each single step equal their local values, so G_{t'→t} = G_{m'→m''}
  (gate G-L1).
- **Step factorization.** With junction words in factored form W = Z·T_K, every step is
  **S = id_Z ⊗ s with s ∈ Hom(T'_K·ρ, T''_K)**. This is a window on the K uprights only: bottom (c_0, c_1, c_2),
  N = μ'_i + μ'_{i+1} ≤ 9 roots, ≤ 4 rungs.
- The top-to-top matrix of S equals the local one up to nonzero diagonal scalings, so invertibility transfers.
- Z·T_K is a word through a level of ≤ 4 blocks (3 row parts + the other column) and needs ≥ 4 uprights. The chart
  words of §1.8 need up to #cells ≤ 9.

### 1.8 Charts and rex moves
- **Chart word** C_K(t) = Z(z)·T_K(m) is a junction word factored for column pair K:
  - Z(z) is CORE's atomic word of z (Thm 3.11 [P]): split rows into (K part, other cells), cross by a reduced word of
    the squashed permutation ψ(z), merge the other column;
  - T_K(m) is the local canonical ladder of m (3 uprights inside K) or its atomic word.
- **Rex moves.**
  - When consecutive rungs use the same K, z is unchanged and no rex move is needed.
  - When K changes, a rex C_K(t) → C_{K'}(t) is needed. Any degree-0 map with nonzero top scalar between reduced
    words of t will do, because the proof uses only "degree 0, 1⊗ ↦ 1⊗" (SLL Thm 6.11, Cor 6.12, Prop 4.19 [P];
    L1 §5.2). In Lemma A it is one more step with G = 1.
  - Construction: go through the atomic word Y_t (all cells) by local moves, each a degree-0 `hom_space` window on
    ≤ 4 adjacent blocks of cell size ≤ 3 (≤ 9 roots):
    - (a) split associativity (down-down, E4 `Assoc`);
    - (b) fork slides / type-A switchbacks with d = 2 (EKo eq. (108) [P]);
    - (c) atomic braid moves = Coxeter moves on ψ (CORE Lemmas 3.12, 3.13, Thm 3.14 [P]).
  - Existence of a chain: EKo Thm 5.30 [P]. The planner (§2.3) is new algorithmic work.
- **Entry and exit.**
  - For X no rex is needed: the recursion starts from the empty word and consumes X rung by rung.
  - For L_q (pairs), use the top-path leaf of the word L_q into its chart, L_q → C(q), and the dual leaf C(q) → L_q.
    Both have degree 0 and a nonzero top scalar, and both come from the same recursion with one leaf.
  - So π = exit ∘ LL and ι = LL^∨ ∘ entry. One exit pair per q is shared by all π and ι of that q.

### 1.9 Double leaves (full / end0, needed only for the B2 gates)
- Candidates for Hom^d(A, B): LL^∨_{r}(C_r → B) ∘ Dot_g ∘ LL_r(A → C_r), over middle types r. Here g runs over the
  Schur-monomial basis of R^{LR(r)} over R^λ, placed on the cell segments of the atomic chart word C_r (Dot on
  those uprights; SLL Def 7.44, Rem 7.48).
- Thm 7.49 [P] says that EKLP's version is a basis; our family spans the same space by §1.4 applied on both sides.
- `certify_full_plan` (:635) selects the basis at the Hom-formula count.

## 2. Code

### 2.1 Files and the D11 freeze
| file | stage | change |
|---|---|---|
| src/dcoset.rs (new) | L1 | combinatorics only: configuration → table, G_{t'→t}, reachability, reduced-word search, cells/ψ/z/m (CORE §3), chart words |
| src/leaves.rs (new) | L1, L2, L3 | `Leaves` implements `homprov::HomProvider` through public API only |
| src/sindec.rs | L1 | `provider_for_stamped` (:421): insert Leaves, extend the stamp |
| src/smap.rs (frozen by D11) | L2 | additive: padded bottoms + `MapNode::Pad` (§2.4); needs a recorded D11 exception |
| src/smapio.rs (not frozen) | L2 | `P` record for `Pad`; padded Expr bottoms |
| src/bin/libbuild.rs | L1 | `--provider leaves` (Leaves alone) for acceptance runs |
| tests/leaves.rs (new) | all | V0-V4 below |

L1 touches no frozen file. L2 needs one additive smap.rs commit, scheduled when no open branch edits smap.rs (after
M6a-1b), and reviewed alone.

### 2.2 dcoset.rs (pure combinatorics, seconds)
- `fn conf_type(lam, conf: &Conf) -> Table`, `fn rep_conf(s, t) -> u64` (re-implements private `rep_rank` :529).
- `fn step_gf(s: &mut Schur, t: Tid, r: &Rung) -> Vec<(Tid, Lp)>`: G_{t→·} via `split_col`/`merge_cols` on h_t.
- `fn reach_back(s, w: &[Rung], targets: &[Tid]) -> Vec<BTreeSet<Tid>>`: types at each level that reach a target
  (backward pass over `step_gf`), plus degree windows [dmin, dmax] per (level, type) for pair pruning.
- `fn junction_word(s, t, chart: Chart) -> Vec<Rung>`:
  - `Chart::ThreeUp`: BFS over rung words of length ≤ 5 that pass the §1.1 test (memoized `ladder_char`), with
    minimal flows; fallback `canonical_ladder`.
  - `Chart::Factored(K)` / `Chart::Atomic`: the explicit words of §1.8, with a "spread" prefix that moves whole blocks
    onto padded zero uprights (trivial rungs, one flow each).
- `fn grass(t_from: Tid, r: &Rung) -> Grass { z, m_from, m_to: Vec<Tid>, local_lam: [u8;3], k_uprights }` (SLL eq. 11).
- `fn psi(t) -> Vec<u8>` (CORE Thm 3.8), `fn rex_plan(w_from, w_to, t) -> Vec<Move>` (§2.3).

### 2.3 leaves.rs
```rust
pub enum Chart { ThreeUp, Factored }                 // L1 uses ThreeUp; L2 uses Factored (padded)
pub struct Leaves { pub chart: Chart, pub costs: Option<CostTable>, pub win_cap: (u64, u64) /* fibre, unknowns */,
                    pub max_win_s: f64, pub pad: u8, stats: SearchStats, last_fail: Option<SearchStats> }
struct Junction { level: usize, t: Tid, w: ExprId, fam: Vec<(MapId, i32)> } // rows invertible at rep(t)
struct StepFam { from: Tid, to: Tid, maps: Vec<(MapId, i32)>, gf: Lp }       // cached by (W_from, rung, W_to)
impl<F: Field> HomProvider<F> for Leaves { /* name "Leaves"; closure; pairs; full; end0 */ }
```
- **`build_leaves(ar, s, x: ExprId, targets, dir: Out|In, deg_win) -> Vec<Junction>`.** This is Lemma A level by
  level. At level k+1, for each kept target t:
  - fam_t = ⋃_{t'} { Comp(S, Kron(ℓ, [ρ])) }, using `Kron` (prefix grows on top, so no left whiskering is needed).
  - In the `In` direction: Comp(Kron(ℓ, [ρ]), S^∨).
  - Assert full row rank at rep(t) at one point (diagnostic).
- **`step_family(ar, s, w_from, ρ, w_to, gf, dir)`.** Candidates, in order:
  1. `Id` when w_to = w_from·ρ;
  2. E2/E3 `Elem` when the step is a cup or cap pattern;
  3. exact windows: for d ∈ supp(gf) ascending, `TmplLib::try_hom_basis_cap` on the minimal differing window. That is
     the whole short words (ThreeUp), or T_K·ρ → T''_K on the K uprights (Factored). Each becomes a `Local` node.

  Then greedy rank selection on the top-to-top columns: flows f of w_from·ρ at rep(to) with
  `conf_type(FlowIndex::conf_at(f, |w_from|)) == from` (sexpr.rs :229). Admission is by `CostTable::predict` or by
  the caps; otherwise `Refused`.
- **`closure(lp, p, None, outs)`.**
  - `build_leaves(lp, [id], outs ? Out : In)`, then `certify_closure` (:945).
  - Private helpers re-implemented from homprov: `cdiag_degrees` :989, `identity_coset` :538.
  - `ep = Some` → `Refused` (no kept (3,3,3) type needs it at N ≤ 6, HANDOFF_M4b).
- **`pairs(x, q, lq, None, s, mult)`.**
  - Run `pairs_precheck` (:880).
  - Build π = Comp(Pad, Comp(exit_q, leaf)) of degree −s and ι likewise; each leaf is padded at entry (`Pad: X → X⁺`).
  - Then `certify_pairs` (:790). Offer ≥ mult candidates; it selects a rank-mult minor.
- **`full`/`end0`** (L3): the double leaves of §1.9, then `certify_full_plan` / `orbit_fallback`.
- **`rex_plan`** (L2):
  1. A deterministic rewrite: split associativity, then fork slides, then bubble-sort Coxeter moves from ψ-words.
     Every intermediate word is checked as reduced for t (§1.1).
  2. If a rewrite rule does not apply: a bounded BFS over ≤ 3-rung window rewrites on padded words. Every
     intermediate word is checked reduced, at most 2 000 words.
  3. Last resort: one exact call on the differing window when admitted.

  Each move is a degree-0 `Local` template, taking the basis element with nonzero top scalar.
- **Arena discipline.** Every call runs between `ar.mark()` and `ar.commit(mark, keep)` (smap.rs :891/:894). Junction
  and step DAGs are shared across the targets of one call; the StepFam cache holds MapIds only inside one call, and
  TmplIds through `TmplLib`.

### 2.4 smap.rs change for L2 (additive, ~120 lines + tests)
- `Arena::expr` (:580) accepts `bottom == lam ++ [0; z]` for z ≤ `MAX_PAD`, instead of the strict assert.
- `MapNode::Pad { src, tgt }`: the same rungs, bottoms differ by trailing zeros, degree 0, `certified = true`.
- **Eval.** The block is the identity. This holds because `CosetSpace::rank` is invariant under trailing zero uprights
  (radix C(0,0) = 1, coset.rs :120) and `FlowIndex::build` enumerates the same moved masks. A unit test asserts the
  identical order; otherwise the block is a stored permutation.
- `assert_genuine`, `children`, `commit` and smapio handle `Pad`.
- Nothing else changes. `Local`, `Elem` and `Kron` already work on any number of uprights.

### 2.5 Wiring and stamps
- **Production chain** (`provider_for_stamped`, max label 3): `[Search S0-S2, Leaves, Exact(l3), Search S3-S6]`.
  Leaves comes before Exact (600 s calls) because its windows are small.
- **Stamp:** `chain=S0-2,Leaves(v1,chart=..,pad=z),Exact(l3),S3-6`. A new stamp invalidates layer A. A warm rebuild of
  level 6 is about 26 min per prime (M5b_report), and every load is re-certified anyway.
- **Acceptance runs:** `libbuild run --provider leaves` (Leaves alone), as G1(a) and B2 require below.

## 3. Why the failing classes are covered
| class (N = 6, both primes) | what Leaves builds | why complete | largest exact call |
|---|---|---|---|
| closure, 12 types (Exhausted) | identity-path leaves of L_p, ins and outs | §1.4-1.5 (Prop 8.32 at the Bruhat-minimal coset; Lemma 8.7 count = c_diag) | ThreeUp: ≤ 4-rung whole words at step degrees; Factored: K windows |
| pairs, rex class 71, 159, 322, 323, 325, 334, 419 | the top-path leaf of X into the chart of q, then exit_q (degree 0, top scalar ≠ 0) | §1.6-1.8; Thm 6.11, Cor 6.12; EKo Thm 5.30 | K windows and braid-move windows (≤ 9 roots) |
| pairs, split class 268, 336, 361, 368 | leaves of degree −s to q | §1.6 + Thm 7.52 [K]; `certify_pairs` checks the rank | same |
| Refused 38, Memory 7 (all pairs) | same as above; X is consumed rung by rung, never as one window | the call size is bounded by local windows, independent of X (≤ 21 600 flows) | same |

The Search window of ≥ 495 k unknowns (M5b_measure §3) never occurs. The only whole-word calls are ThreeUp closure steps
between ≤ 4-rung words. Gate G-L0 lists their sizes before any build.

## 4. Certification
- **Family level: unchanged and the only acceptance criterion.** It consists of:
  - `certify_closure`: ends, degree multiset = c_diag(b_p), identity-coset rank n at 2 points, genuineness;
  - `certify_pairs`: ends, degrees, form rank = mult at rep(q) at 2 points, genuineness;
  - `certify_full_plan`: count = Hom formula, stacked rank, fresh residual;
  - downstream `dechom_postcheck` (sindec.rs :330, D14 scaling), `dec_one`, `check_pi_iota`, and the replay of
    layer A.
- **Leaf level (genuineness).**
  - hom_space templates are genuine by construction.
  - `Elem` constants are recorded once per (kind, labels) through `certify_elem` or `certify_assoc` (A7/A8 style).
    Labels go up to 9.
  - `Pad` is an identity. A test checks it against `Id` on equal words, and a negative test uses a non-trailing zero.
- **Internal rank tests** (junction invertibility, top-to-top selection, nonzero top scalars of rex moves) are
  diagnostics. They speed up failure; they are never used to accept anything.
- **Wording (D13) is unchanged.** Results are over F_p for both primes.

## 5. Validation plan (each gate with Leaves alone, then in the chain)
- **G-L0 census** (S1, combinatorics, seconds). For the 12 closure types and the 56 failing pairs, list:
  - junction types per level, G, junction word lengths;
  - `hom_unknowns` (homprov.rs :344) of every ThreeUp step call and every K window;
  - predicted cost from costs_333.txt buckets (≤ 10³ unknowns < 1 s, 10⁴ ≤ 24 s, 10⁵ ≤ 420 s [measured]).

  Decision rule: closure stays ThreeUp if max ≤ 10⁴ unknowns; otherwise closure also uses the Factored chart.
- **G-L1** (combinatorics): Σ_{t'} c_{t'}·G = `ladder_char` for every canonical ladder of (3,3), (2,2,2), (2,1,2) and
  (3,3,3), and for each one-rung extension. Also check `sexpr::deodhar_check` (:267), G_local = G_global (§1.7), and
  the Lemma 8.7 sums = `c_diag` and `term_mult`.
- **V0** (exact oracle, tiny λ ∈ {(2), (1,1), (2,1), (1,2,1)}): leaves vs `sbim::hom_space` bases (stacked rank).
  Also SLL Ex 7.39-7.41 and 7.46 (L1 §5.4 G-L2), and the degree normalization of §1.1.
- **V1 = G1(a) with Leaves alone:** 480/480 targets ((3,3) 72 + closure 8; (2,2,2) 358 + closure 42), both primes.
  The closure part after L1, the pair part after L2.
- **V2 = B2 through Leaves alone (after L3):** (2,1,2) 605/605 (Search alone: 581) and (2,2,2) 2 205/2 205 in
  degrees 0..4 (the Chain fails at 2 064/2 205, HANDOFF_M4b). Also (2,2) 45/45 and (3,3) 80/80.
- **V3 regression:** level-n4 47/47 and a random 30-job level-5 sample with Leaves alone. Compare the certified
  multisets with the existing layer A (same targets, equal rank).
- **V4 = M5b failure list:**
  - the 12 closure types (after L1);
  - the 11 Exhausted, 7 Memory and 38 Refused pairs (after L2);

  each at both primes. The aim is 68/68 and then a full level-6 run of 377/377 with the production chain. Then replay
  with the refusing provider, and persist to data-gl6 (D12).

## 6. Cost estimates [E unless marked]
- **Closure, ThreeUp.**
  - ≤ 6 levels × ≤ 10 junction types gives ≤ 150 step families per type. Each is 1-4 exact calls between ≤ 4-rung
    words at step degrees, typically 10³-10⁴ unknowns, i.e. 0.5-20 s [costs_333 buckets].
  - With template reuse across the 12 types: **1-3 h per prime cold, minutes warm**.
  - Certification evaluates only identity-coset blocks: n ≤ 54 rows, so it is negligible.
- **Closure, Factored.**
  - K windows have N ≤ 9 roots split into ≤ 3 row parts with ≤ 4 rungs. This is the measured class of (2,2,2) and
    (3,3) windows, 10²-10⁴ unknowns, < 1-5 s.
  - **< 1 h per prime cold.**
- **Pairs, Factored.**
  - 229 uncertified targets (Σ mult 393). Each takes a pruned leaf tree of ≤ 8 levels × ≤ 20 types.
  - Distinct K and braid-move windows across the library: a few hundred to 2 000, each ≤ 5 s and cached.
  - Rex chains: ≤ 40 moves per change of K.
  - FlowIndex builds of chart words: ≤ 10⁵ flows × ≤ 50 rungs × 8 B ≤ 40 MB each. 110 592 flows is the worst atomic
    word (1,1,1/1,1,1/1,1,1, L1 §5.5). They live in an LRU of 256 MiB-1 GiB.
  - **2-6 h per prime cold** for the 56 pairs; dominated by index builds and windows, not by certification (one coset).
- **DAG size and downstream cost.**
  - A π from Leaves has 50-300 nodes, against 5-20 for Search, shared by prefix.
  - M6a evaluates maps globally per point, so per-point cost for these maps could rise 5-20×. Gate G-L5 measures it on
    the level-6 library.

## 7. Effort (cloud sessions of about 4-5 h)
| session | content | exit gate |
|---|---|---|
| S1 | dcoset.rs: conf_type, step_gf, reach_back, junction_word(ThreeUp), grass, psi; census tool | G-L0, G-L1 |
| S2 | leaves.rs core: Lemma A both directions, step_family (Id/Elem/exact), closure; libbuild `--provider leaves` | V0, V1 closure part (25 types) |
| S3 | (3,3,3): 12 closure types × 2 primes; layer-A stamp; replay | V4 closure 12/12 |
| S4 | smap.rs `Pad` + padded exprs (D11 exception), smapio, tests; chart words (spread, Z, T_K, atomic Y) | Pad and round-trip tests |
| S5 | factored steps, exit/entry leaves, rex_plan (rewrite + BFS), pairs() | V1 pairs part: 480/480 |
| S6 | (3,3,3) pairs: the 11 + 7 + 38 list, both primes | V4 68/68 |
| S7 | double leaves `full`/`end0` | V2 B2 |
| S8 | full level-6 production run, G-L5 downstream cost, report | 377/377, replay |

**Total: 8 sessions (35-40 h), plus 2 in reserve for planner debugging.**
- Closure-only M4c (S1-S3, about 3 sessions) is self-contained. It does not need the D11 exception.
- SYNTHESIS estimated 10-20 days. The saving comes from selecting step leaves by rank instead of compiling the
  sinister diagram, and from avoiding Y_p in closure.

## 8. Risks and fallbacks
1. **Rex planner (L2) is the only piece without a closed form.** The rewrite rules (§1.8 a-c) might not cover a
   configuration.
   - Fallback 1: bounded BFS.
   - Fallback 2: one exact difference window (admission-gated).
   - Fallback 3: increase the padding to the full cell chart Y_t, where CORE Thm 3.14 reduces rex moves to Coxeter
     moves on ψ.
   - Fallback 4 (outside M4c): "word-agnostic types" in sindec (L3 §5), which removes the exit rex.
2. **ThreeUp closure calls too large** (G-L0 says > 10⁴ unknowns). Run closure in the Factored chart, i.e. S4 before
   S3.
3. **Chart words blow up flows** (atomic words up to 25× the minflow ladder). Mitigations:
   - prefer the 4-upright factored chart, and refine only the columns the next rung touches;
   - per-coset evaluation only during the build; LRU `drop_index`.
4. **D11 exception refused or delayed.** Ship closure-only M4c (S1-S3). Pairs then stay with the parallel engineering
   route (calibration, S4 fix, admission).
5. **Downstream evaluation cost** (G-L5). Options:
   - memoized per-point LocMaps;
   - an optional compaction pass that replaces a certified Leaves map by a cheaper certified Search map spanning the
     same form, when one exists.
6. **Theory gaps [K].** These are:
   - the rung-level combination of SLL Thm 2.6;
   - the singular local-intersection-form rank (Thm 7.52);
   - the degree normalization.

   Each only affects whether candidates exist, never correctness: certification rejects anything wrong. G-L1 and V0
   test them first.
7. **CORE is a preprint (v2).** It is used only to *construct* chart words and to plan moves. A mistake there shows up
   as a failed reducedness or rank test, never as an accepted wrong map.

## 9. First actions
1. S1 census G-L0 on the 12 closure types: it decides between ThreeUp and Factored.
2. One `KRHOM_TRACE` run per rex-class pair on the cloud, to confirm X ≅ B_q (L3 §3.3).
3. File the D11 exception request for the `Pad` node together with S4's diff (about 120 lines), for review before
   S5 starts.
