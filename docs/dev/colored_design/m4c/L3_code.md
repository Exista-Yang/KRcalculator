# L3: M4c code map and failure analysis (read-only study, 2026-11-11)

Sources: dev 277a21c (src identical to colored-m5 except libbuild/xscan/xengine and a small sindec delta);
origin/colored-m5 docs/dev/colored_design/gl6/M5b_report.md and M5b_measure.md; origin/data-gl6 m5b/f31 run logs
(run_n6.jsonl, .targets, costs_333.txt) and the 1 469 layer-A key lines (canonical ladders L_p); tests/data/schur/cb_ref.txt
(canonical basis, c_diag). Nothing was compiled or run. Tags: [C] read in the code, [L] from the run logs, [D] derived here
by arithmetic on cb_ref / logs, [?] not verified.

## 0. Findings in brief

1. **Plug-in point.** A new `src/leaves.rs` implementing the public trait `homprov::HomProvider`, inserted into the chain built by
   `sindec::provider_for_stamped`. D11 freezes homprov.rs, smap.rs, tmpl.rs, sbim.rs, linsys.rs and g1.rs, so everything must be
   built on their public API. All certification functions are `pub` and apply unchanged (§4).
2. **Node kinds.** The existing DAG nodes are sufficient: Dot, the closed forms Elem (E2/E3/E4), Local, Comp, Lin and Kron.
   There is no prefix-whiskering node (id_P ⊗ g for a DAG g); Kron only appends identity rungs on the top side.
3. **Closed forms are never used.** E2 units and counits, E3 digon inclusion and trace, and E4 associativity exist with
   localized formulas (`smap::ElemKind`, `elem_winloc`). No provider uses them. Every map Search builds bottoms out in exact
   `hom_space` window calls. The closed-form constants are recorded only in tests (tests/m4a.rs A7/A8).
4. **Closure failures (12 types at N = 6).**
   - 11 of the 12 fail at their **lowest** c_diag degree. The S1 whole-ladder exact call has 0.3-7.5 M unknowns (or no record)
     and is refused or censored.
   - S4, the only constructive closure tier, is **unreachable** in the production chain: `Search::tiers(3,6)` has S0 and S1
     off, so its recursion never reaches a base case [C].
   - 1,1,1/1,1,1/1,1,1 fails at degree 6: the S2 dot multiples give 3 of the 4 generators.
5. **5 of the 12 failing closure types are products [D].** L_p is commutation-equivalent to A·B, where A and B are canonical
   ladders of *certified* (λ,λ) types and c_diag(p) = c_diag(A)·c_diag(B). Their generators are Comp(Kron(a_i, B), b_j),
   which is buildable today. For example, 0,0,3/3,0,0/0,3,0 = Eta∘Eta, degree 9 + 9 = 18.
   The other 7 (2 700, 2 880, 4 320 and 7 200 flows) need genuinely 3-upright light leaves.
   - For the 2 880, 4 320 and 7 200 types, c_diag(p) = (1 + v⁻² [+ v⁻⁴])·c_diag(A)c_diag(B): their identity paths pass
     through non-identity cosets at the junction.
   - The 2 700 words do not split at λ.
6. **Pair failures (11 Exhausted at N = 6).**
   - 7 of them (71, 159, 322, 323, 325, 334, 419) are rex moves [D]. X has exactly r(q) flows and one target (q, s = 0,
     mult 1), so X ≅ B_q while the word of X differs from L_q. The missing pair is a **degree-0 singular rex move X ⇄ L_q**.
   - The other 4 (268, 336, 361, 368) are genuine splits whose window tree is stuck.
7. **Refused (38) and Memory (7) failures** have one root cause. λ has 3 uprights, so Search's "local" 3-upright windows span all
   9 roots, and their exact calls reach 10⁵-10⁷ unknowns. A leaves provider makes no such call except for 3-upright rex moves (§5).
8. **What a leaves provider outputs.**
   - Pairs: the candidates ι (flipped leaves L_q → X) and π (leaves X → L_q) of defect s that end at q. Their form at rep(q) has
     rank mult, and `certify_pairs` selects from them.
   - Closure: the n flipped (or plain) leaves of the identity-coset paths of L_p. Their degree multiset is c_diag(b_p), and
     `certify_closure` checks rank n.
   - New certification is needed only for the leaf building blocks: Elem constants and rex templates.
9. **Cheap engineering findings for the parallel session** (§2.3): the S4 reachability bug; `leaf_iso` labels an admission
   failure as Exhausted; the 2x retry does not raise l_max; Exact bucket admission vetoes measured 237-394 s calls.

## 1. Where a complete provider plugs in

### 1.1 Interface and chain [C]
- `trait HomProvider<F>` (homprov.rs:461-481) has `name`, `full`, `pairs`, `closure`, `end0` and the default `failed_stats`.
  - Results: `Pairs<F>{iotas, pis, form}` (:455), `Vec<MapId>`, `SearchStats` (:388).
  - Errors: `Vec<Missing>` with `MissKind {Exhausted, Time, Memory, Refused}` (:437).
- Each call runs between `ar.mark()` and `ar.commit(mark, keep)` (smap.rs:891, 894; Search's `run_call` homprov.rs:3126). Only
  the certified results and their sub-DAGs survive. Leaves must follow the same discipline.
- `Chain<F>` (:3211): the first Ok wins, the Missing lists are merged, and the stats of failed providers are absorbed.
  - (3,3,3) does not use `Chain::for_labels` (:3216). `sindec::provider_for_stamped` (sindec.rs:421) builds
    `[Search S0-S2, Exact{max_label 3, costs, 600 s}, Search S3-S6]`.
  - The stamp `chain=S0-2,Exact(l3,600s),S3-6` goes into every layer-A header.
- Adapter `sindec::DecHom` (sindec.rs:283), a blanket impl for every HomProvider:
  - `dechom_precheck` (:298): L_q is the canonical ladder, indecomposable without e_q, and mult = [v^s] of b_q in char X;
  - `dechom_postcheck` (:330): counts, ends, degrees, and the D14 degree-scaling check at rep(q);
  - `prov_fail` (:362): any Exhausted in the chain makes the failure a completeness failure; otherwise it is engineering.
- **Wiring M4c:**
  - Push `Box::new(Leaves::..)` into that vector: after S0-S2 for production, or first for the M4c list.
  - Extend the stamp, for example `chain=S0-2,Leaves,Exact(..),S3-6`.
  - The new stamp invalidates layer A. Rebuilding level 6 takes about 26 min per prime warm (M5b_report), and every loaded map is
    re-certified anyway.
- **D11 constraint** (SPEC_3C_GL6_rev1 §D11). The planned file name `leaves.rs` is from SYNTHESIS §3 (file table).
  - Private helpers to re-implement (a few lines each):
    - `cdiag_degrees` (:989), via `Schur::c_diag` (schur.rs:1361);
    - `identity_coset` (:538) = `CosetSpace::new(λ,λ)` ranked at `roots.bottom()`;
    - `rep_rank` (:529);
    - the closure `Selector` (:1396).
  - Public and reusable:
    - certification: `certify_pairs` (:790), `certify_closure` (:945), `certify_full_plan`/`_at` (:635/:643), `detect_plan`
      (:591), `form_at` (:726), `check_idempotent` (:837), `pairs_precheck` (:880), `closure_precheck` (:912), `hom_gens` (:99);
    - word helpers: `hom_unknowns` (:344), `win_split` (:197), `commutation_class` (:158), `rungs_commute` (:143);
    - dots: `segments` (:1470), `dot_pmat` (:2545).

### 1.2 The map DAG (smap.rs) [C]

| node (smap.rs:245) | meaning | degree | genuine if | evaluation (Eval::block :1092) |
|---|---|---|---|---|
| `Id(e)` | identity | 0 | always | identity block |
| `Dot{e, level, upright, nu}` | multiply by s_ν of the roots on one segment; level 0 = left R^λ action, level = len: right action | 2\|ν\| | always | diagonal per flow |
| `Elem{kind, src, tgt, at, upright, scalar}` | closed form on the 2-upright window at rung `at` | per kind | constant recorded (`certify_elem` :821, `certify_assoc` :831) | `elem_winloc` (:434) per flow |
| `Local{src, tgt, at, first, tmpl}` | id_P ⊗ template ⊗ id_Q (any width) | tmpl.deg | template HomSpace / Rickard, or Added after `certify_added` | `window_block` (:1189): localized template |
| `Comp(g, f)`, `Lin(consts)` | composition; constant linear combination (equal degree) | sum / equal | children genuine | gemm / sum |
| `Kron(f, rungs)` | f ⊗ id on rungs appended **on top** | deg f | f genuine | `kron_block` (:1145): f at the prefix coset of each flow |

- `Arena::add` (:667) validates windows, ends and degrees.
- `assert_genuine` (:855) walks the leaves of the DAG.
- smapio.rs serializes every kind, including Elem as `M` records, so layer A can store leaves. Before an imported Elem passes
  `assert_genuine`, its constant must be recorded in the Arena.
- **Missing: left whiskering of a DAG**, i.e. g ↦ id_P ⊗ g. Two workarounds:
  1. Build maps bottom-up, so that prefixes grow through `Kron` and windows act on top.
  2. Rebuild g with every `at`/`level` shifted by |P| and with re-interned expressions P·src, P·tgt (MapNode is pub).

### 1.3 Elementary maps that already exist [C]

| map | where | ladder / foam meaning | degree | certification | used by |
|---|---|---|---|---|---|
| Dot (E1) | MapNode::Dot | polynomial on a facet segment | 2\|ν\| | none needed | S2, local_cands |
| Eps / Eta (E2) | ElemKind (:224) | [E_u^(b) F_u^(b)] over (a,b), all of upright u+1 merged and split back ↔ empty: counit / unit. Eta has value e(B,A) | +ab | A7: compared with hom_space on 2 uprights (tests/m4a.rs:760-820) | **none in src** |
| EpsL / EtaL | mirrored [F^(a) E^(a)] | same | +ab | A7 | none |
| Iota / Trace (+L) (E3) | digon [F^(b) E^(b)] over (c,0): inclusion / Frobenius trace Σ_S m(S)/e(S, A∖S) | −ab | A7 | none |
| Assoc (E4) | identity in value coordinates. Valid only when the window flows are in bijection with top configurations | 0 | `certify_assoc` per window pair (A8) | none |
| hom_space templates | `TmplLib::try_hom_basis_cap` (tmpl.rs:441) | any exact window map | any | genuine by construction | Exact; Search splits (:1909), commutations `comm_tmpls` (:1792), rex `iso_move` (:2249), `leaf_iso` (:2342), S6 `window_maps` (:2454), closure units in `closure_seeds` |
| Rickard | `TmplLib::rickard` (tmpl.rs:486) | crossing differentials d_j | | genuine | sindec `rickard_comp_at` |
| Added | `add_template` (:375) + `certify_added` (:392: degree, exact bimodule check, `poly_certificate` :677) | Search's window projections π' | | after certify | `split_compute` (:1926) |

SYNTHESIS §2.4 specified Search paths built from the E2/E3 closed forms, fusion/fission and square switches. The implemented
Search replaced all of these by exact windows. This is the design difference that matters at (3,3,3), where any 3-upright window
spans all 9 roots.

### 1.4 Combinatorics a leaves enumerator can reuse [C]
- schur.rs is the path-counting engine:
  - `split_col` (:913, EKLP (44)): h_q ↦ Σ_p v^{lmax q − lmax p} h_p;
  - `merge_cols` (:931, EKLP (45)): coefficient π(LR q)/π(LR p);
  - `ladder_char` (:997): a product of rung actions.
- A path enumerator is `ladder_char` with the terms kept apart: one path per step choice, and one per monomial of the merge
  coefficient. Gate (SYNTHESIS §3, M4c): Σ paths = ladder_char (EKLP Cor 3.7) and = `hom_gd` (Cor 4.17).
- Other schur.rs data: `hom_gd` (:1289, Thm 4.16), `c_diag` (:1361), `term_mult` (:1389), `killed_at` (:1436), and canonical
  ladders `canonical_basis(μ).ladder_minflow` (:636).
- coset.rs `CosetSpace`: `rank` (:116), `rep` (:195), `cosets_of` (:211).
- sexpr.rs `FlowIndex` (:119): a flow is one root subset per rung, i.e. a fixed point; it has `conf_at` (:229) and
  `coset_of_flow` (:250).
- #flows of X at rep(q) = c_q(char X)(1) = #paths of X ending at q. This holds as counts only; a canonical flow ↔ path bijection
  is [?].

## 2. Why Search fails

### 2.1 Closure (closure_inner :2729, closure_seeds :2859) [C]
- Degrees come from c_diag in ascending order.
- **Lowest degree d0:** only S1 (:2792, `tier = if d == d0 {1} else {6}`), via `window_maps([] → L_p)`. `win_split` against the
  empty word gives the whole ladder, so the fibre product is flows(L_p) and the unknowns are those of a whole-ladder call.
- **d > d0:**
  - S2: `Dot∘g` for chosen lower generators g, with dots on *interior* segments only (cap 200, ≤ 3 parts);
  - then S6, the same whole-ladder window at d;
  - then S4.
- **S4 is dead code in the (3,3,3) chain.**
  - S4 runs only in the `Search::tiers(3,6)` instance.
  - There `on(0)` and `on(1)` are false, so the recursive `closure_inner(p2, depth+1)` for a smaller type cannot find p2's
    lowest generator (S1 off), and the identity base case (S0, :2747) is off too.
  - Every recursion therefore ends at `depth >= 4` (:2860) with empty seeds, and the result is "0 candidates", Exhausted.
  - The level-6 targets log confirms it: closure calls hit only S1 + S2 (`tiers [0,1,k,0,0,0,0]`), never S4 [L].
- Even if reachable, S4 covers only one map shape: π_leaf ∘ (id ⊗ the lowest unit of [R, R^∨] at the top) ∘ g'. It considers
  only the first 3 smaller types by rank (`take(3)`, :2897), down to depth 4. In light-leaf terms these are only the paths that
  sit at the identity coset just before the last two rungs.
- S2's interior dots act on an identity-coset vector by flow-dependent scalars. They cannot produce generators that are not
  R_interior-multiples of lower ones.

### 2.2 Pairs (pairs_inner :2604) [C]
- S1 (extremes) and S3 (target degree) take the leaves of `tree` (:2173).
- `tree` uses `choose_window` (:2124), which looks only at the **commutation class** of the word (cap 64, :2129), never at rex
  moves. Windows have width ≤ 3, which is the full width for 3 uprights, and length ≤ 4 (`win_w`/`win_len` :1693-1694).
- `split` (:1909) needs the exact window, or the extreme times dots (`local_cands` :2076). The latter is incomplete and fails as
  "local form rank != cnt", which counts as Exhausted.
- Depth `l_max` = 6 (:2195) counts as a stuck branch without admission.
- `leaf_iso` (:2342):
  1. commutations;
  2. meet in the middle over `iso_reach` (:2276, cap 400 words). It only rewrites an indecomposable window of ≤ 4 rungs to its
     canonical local word, with `max_iso_unknowns` 6 000;
  3. one exact window on the non-common part. When that window is not admissible, it returns **`MissKind::Exhausted`** (:2371),
     although the cause is admission.

### 2.3 Engineering points (cheap, for the parallel session) [C]/[L]
- **(a) S4 is unreachable** (§2.1). Fix: run closure S4 in a Search with tiers 0-6, or let the S3-S6 recursion use S0/S1.
- **(b) leaf_iso misreports admission as Exhausted** (:2371). This is why the rex-class pairs count as "Exhausted at 2x".
- **(c) The retry does not raise l_max.** libbuild's `scaled_budget` doubles t_max, m_max, index_max and c_max. SPEC §4.7 asked
  for l_max + 2 and 4·c_max. Exhausted at 2x is therefore the same deterministic stream.
- **(d) Exact bucket admission** (M5b_report open problem 3). The 600-flow closure calls at d = 13 are measured uncensored at
  237-394 s, but a censored record in the same bucket vetoes them ("predicts inf").
  - Fixing (d) alone does not certify these types: d = 15 and d = 17 are censored (0.52-0.94 M unknowns).

### 2.4 Failure classes at N = 6 (F31 and F32b identical) [L]

| class | jobs | root cause |
|---|---|---|
| closure Exhausted | 12 types | lowest generator out of reach (11 types), or S2 incomplete (1,1,1/1,1,1/1,1,1). Exact whole-ladder calls 0.3-7.5 M unknowns; S4 unreachable |
| pairs, rex class | 71, 159, 322, 323, 325, 334, 419 | X ≅ B_q, word ≠ L_q: degree-0 rex move; `leaf_iso` fails |
| pairs, split class | 268, 336, 361, 368 | tree stuck (no decomposable ≤ 4-rung window in the commutation class, split incomplete, or depth) [?: the logs do not say which] |
| Refused | 38 pairs (X = 480-21 600) | every window of the tree ≥ 40 k unknowns without a record; large-X windows ≥ 495 308 unknowns measured infeasible |
| Memory | 195, 196, 199, 444, 445, 460, 475 (X = 240-7 200) | a window or whole-ladder exact call above 4-6.8 GB (M5b_report) |

## 3. The failing targets: sizes and what the missing maps are

### 3.1 Closure types (Exact-call cost records from costs_333.txt [L]; c_diag from cb_ref [D])

| p | L_p (canonical, layer-A keys) | flows | n | c_diag (degree:count) | Search stops at | whole-ladder exact call at that degree |
|---|---|---|---|---|---|---|
| 0,0,3/3,0,0/0,3,0 | E0.3 F0.3 E1.3 F1.3 | 400 | 1 | 18:1 | d=18 (lowest) | 852 011 unk, censored (600 s / 3.3 GB) |
| 0,3,0/0,0,3/3,0,0 | E1.3 E0.3 F1.3 F0.3 | 400 | 1 | 18:1 | 18 | same |
| 1,0,2/2,0,1/0,3,0 | E0.2 F0.2 E1.3 F1.3 | 600 | 3 | 13 15 17 | 13 | d13: 327-393 k unk, 237-394 s (vetoed by bucket); d15/17: 0.52-0.94 M censored |
| 1,2,0/0,0,3/2,1,0 | E1.3 E0.2 F1.3 F0.2 | 600 | 3 | 13 15 17 | 13 | same family |
| 0,3,0/1,0,2/2,0,1 | E1.2 E0.3 F1.2 F0.3 | 600 | 3 | 13 15 17 | 13 | same family |
| 1,0,2/1,2,0/1,1,1 | F0.2 E0.1 E1.2 E0.1 F1.2 | 2 700 | 18 | 9:1 11:3 13:5 15:5 17:3 19:1 | 9 | d9: 701 771 unk censored (> 3 GB) |
| 1,1,1/0,2,1/2,0,1 | F0.1 E1.2 F1.2 F0.1 E0.2 | 2 700 | 18 | same | 9 | (bucket inf) |
| 0,2,1/2,0,1/1,1,1 | E0.2 E1.3 E0.1 F0.1 F1.3 F0.2 | 2 880 | 18 | same | 9 | (bucket inf) |
| 1,1,1/1,0,2/1,2,0 | E0.1 F0.1 E1.3 E0.1 F1.3 F0.1 | 2 880 | 18 | same | 9 | (bucket inf) |
| 0,1,2/2,0,1/1,2,0 | E0.2 F0.2 E1.3 E0.1 F1.3 F0.1 | 7 200 | 27 | 10:1 12:3 14:6 16:7 18:6 20:3 22:1 | 10 | (no record) |
| 0,2,1/1,0,2/2,1,0 | E0.1 F0.1 E1.3 E0.2 F1.3 F0.2 | 7 200 | 27 | same | 10 | (no record) |
| 1,1,1/1,1,1/1,1,1 | E0.1 F0.1 E1.2 E0.1 F1.2 F0.1 | 4 320 | 54 | 4:1 6:4 8:9 10:13 12:13 14:9 16:4 18:1 | d=6, 3 of 4 | d4: 224 249 unk, 290-316 s, 2.2 GB (ok); d6 406 634 … d18 7 492 529, all censored |

A leaves construction evaluates each generator only at the identity coset. The identity-coset block has n entries (1-54), and
each leaf is a DAG of a few nodes per rung over ≤ 6 rungs. That is negligible next to 10⁵-10⁷ unknowns.

### 3.2 What the missing generators are

**Product types (5) [D].** Directly, or after commuting E0 past F1 (`rungs_commute`), L_p is A·B, where A and B are canonical ladders of
(λ,λ) types already certified at level 4/5 (closure keys in layer A):

| p | A | B | product |
|---|---|---|---|
| 0,0,3/3,0,0/0,3,0 | [E0.3 F0.3] (0,3,0/3,0,0/0,0,3, c_diag v⁹) | [E1.3 F1.3] (v⁹) | v¹⁸ |
| 0,3,0/0,0,3/3,0,0 | [E1.3 F1.3] (v⁹) | [E0.3 F0.3] (v⁹) | v¹⁸ |
| 1,0,2/2,0,1/0,3,0 | [E0.2 F0.2] (1,2,0/2,1,0/0,0,3: v⁴+v⁶+v⁸) | [E1.3 F1.3] (v⁹) | v¹³+v¹⁵+v¹⁷ |
| 1,2,0/0,0,3/2,1,0 | [E1.3 F1.3] (v⁹) | [E0.2 F0.2] (v⁴+v⁶+v⁸) | same |
| 0,3,0/1,0,2/2,0,1 | [E1.2 F1.2] (3,0,0/0,1,2/0,2,1: v⁴+v⁶+v⁸) | [E0.3 F0.3] (v⁹) | same |

Each row matches c_diag(p) exactly, including n = n_A·n_B.
- The generators are:
  - ins: Comp(Kron(a_i, B), b_j) : R → B → A·B;
  - outs: Comp(b'_j, Kron(a'_i, B));
  - each followed by the commutation isomorphism to L_p (`comm_tmpls`: 2-rung windows, measured cheap).
- The identity-coset matrix is the Kronecker product of two nonsingular matrices (every identity path factors when
  n_p = n_A·n_B), so `certify_closure` passes. This is a theorem-free check, since `certify_closure` verifies it.
- For 0,0,3/3,0,0/0,3,0, a and b are E2 `Eta` closed forms of degree ab = 9 on (3,3). The generator is literally two digon units.
- Search misses these maps for three reasons:
  - S1 tries the whole 3-upright ladder at once;
  - S2 needs a lower generator, and there is none;
  - S4, which has the right shape (smaller type, then a unit on top), is unreachable and capped at 3 options.

**Non-product types (7) [D].**
- **1,1,1/1,1,1/1,1,1**: A = [E0.1 F0.1] (2,1,0/1,2,0/0,0,3, v+v³+v⁵) and B = [E1.2 E0.1 F1.2 F0.1] (2,1,0/0,1,2/1,1,1,
  v⁵+2v⁷+3v⁹+2v¹¹+v¹³).
  - c_diag(p) = (1 + v⁻²)·c_diag(A)c_diag(B): 27 product generators (degrees 6-18) and 27 that are not products (degrees
    4-16, the same polynomial shifted by v⁻²).
  - The degree-4 generator Search found and the 4th degree-6 generator it missed belong to the second family.
- **2 880 types**, e.g. 1,1,1/1,0,2/1,2,0 = [E0.1 F0.1]·[E1.3 E0.1 F1.3 F0.1] (B = 2,1,0/0,0,3/1,2,0, v¹⁰+v¹²+v¹⁴):
  c_diag = (1 + v⁻²)·product, with lowest degree 9 < 11.
- **7 200 types**: (1 + v⁻² + v⁻⁴)·[(v⁴+v⁶+v⁸)(v¹⁰+v¹²+v¹⁴)], with lowest degree 10 < 14.
- **2 700 types**: the L_p words do not split at λ.
- **Interpretation** (light-leaf recursion cut at the junction): identity paths of A·B are (α: path of A ending at a coset c,
  β: path of B from c to the identity). For c ≠ id, the leaf is
  - outs: LL_β(L_c·B → R) ∘ Kron(π_α: A → L_c, B);
  - ins: the flip of the above.

  π_α is a 2-upright decomposition map of A (cheap, (3,3) labels). LL_β is a closure-type map of the *word* L_c·B, which is
  shorter than L_p. Search has neither the junction at a non-identity c nor closure maps for non-canonical words.
- The extra factors v⁻² and v⁻⁴ say how many such junction cosets contribute: 1 for 2 880 and 4 320, 2 for 7 200 [D].

### 3.3 Pairs (Exhausted at 2x) [L]; class from X flows = r(q) [D]

| pair | p ⊗ R (L_p) | X flows | targets / Σ mult | fails at s | whole-ladder Exact: fibre product, unknowns | class |
|---|---|---|---|---|---|---|
| 71 | 1,0,2/3,0,0/0,2,1 ⊗ F0.1 (E0.3 F0.2 E1.3 F1.3) | 600 | 1/1 | 0 | 360 000, 24.7 M | rex (q one of the 6 types with r = 600) |
| 159 | 2,0,1/3,0,0/0,1,2 ⊗ F0.2 (E0.3 F0.1 E1.3 F1.3) | 240 | 1/1 | 0 | 57 600, 2.51 M | rex (r = 240) |
| 322 | 1,2,0/1,1,1/1,1,1 ⊗ F1.1 (E0.1 E1.2 E0.1 F1.1 F0.2) | 4 320 | 1/1 | 0 | 18.7 M, 919 M | rex onto 1,1,1/1,1,1/1,1,1 (the only r = 4 320 type) |
| 323 | 0,3,0/2,0,1/1,1,1 ⊗ F1.1 (E0.2 E1.3 E0.1 F1.2 F0.3) | 2 880 | 1/1 | 0 | 8.29 M, 598 M | rex (r = 2 880) |
| 325 | 1,2,0/0,2,1/2,0,1 ⊗ F1.1 | 2 700 | 1/1 | 0 | 7.29 M, 545 M | rex (r = 2 700) |
| 334 | 0,3,0/1,0,2/2,1,0 ⊗ F1.1 (E0.1 E1.3 E0.2 F1.2 F0.3) | 7 200 | 1/1 | 0 | 51.8 M, 5.50 G | rex (r = 7 200) |
| 419 | 0,3,0/2,0,1/1,2,0 ⊗ F1.2 (E0.2 E1.3 E0.1 F1.1 F0.3) | 7 200 | 1/1 | 0 | 51.8 M, 5.50 G | rex |
| 336 | 0,2,1/2,0,1/1,2,0 ⊗ F1.1 | 14 400 | 2/2 | −1 | 103.7 M, 8.85 G | split: very likely B_q⟨−1⟩ ⊕ B_q⟨1⟩, r(q) = 7 200 |
| 268 | 1,1,1/1,2,0/1,0,2 ⊗ E1.1 (F0.2 E1.1 E0.2 F1.1) | 1 620 | 3/3 | 0 | 2.19 M, 123 M | split |
| 361 | 1,2,0/1,1,1/1,0,2 ⊗ E1.2 (E1.1 E0.2 F1.1 F0.2) | 1 080 | 6/6 | 0 | 291 600, 9.86 M | split |
| 368 | 1,1,1/1,2,0/1,0,2 ⊗ E1.2 | 1 620 | 4/4 | −2 | 437 400, 10.4 M | split |

- **Rex class.** The ranks of the 55 types come from cb_ref (Σ_t N(t) h_tp(1)) [D]. X flows equal the rank of a (λ,λ) type, and
  there is one kept target (q, 0, 1).
  - Unless killed summands make up the difference [?], char X = b_q. X then has a single flow at rep(q), the form is 1×1, and
    the pair is any degree-0 ι, π with nonzero top scalar (C2 M6: such a map is an isomorphism on B_q).
  - Example, pair 322: X = E0.1 E1.2 E0.1 F1.1 F0.2 F1.1 against L_q = E0.1 F0.1 E1.2 E0.1 F1.2 F0.1. After the common E0.1 the
    words differ on all 5 remaining rungs, on all 3 uprights.
  - `iso_reach` (local rewrites to canonical local words, 400 words) found no common word. The whole-difference window is not
    admissible (fibre product of order 10⁶-10⁷). The light-leaf view does not avoid this: the top path of X ends at the reduced
    expression X, so a singular rex move to L_q is unavoidable here.
- **Split class.** These need leaves of X at the extreme defect (−1, −2) or at s = 0 with mult 3-6. The rex part comes in only
  through the final iso to L_q.

## 4. What a leaves provider must output, and how certification applies unchanged

**4.1 `pairs(x, q, lq, eq=None, s, mult)`** (all kept (λ,λ) types have e = None; the only decomposable (λ,λ) ladder,
0,1,2/1,1,1/2,1,0, is killed at N ≤ 6, R4 §1).
- Output any ≥ mult candidates:
  - ι_f : L_q → X of degree s, the flipped light leaf of a path f of X with defect s ending at q, composed with the rex move
    L_q → (end expression of f);
  - π_g : X → L_q of degree −s, likewise.
- `certify_pairs` (:790) is unchanged. It checks ends and degrees, computes form[j][i] = block(π_j∘ι_i, rep(q)) at two points
  (`form_at` :726: L_q has one flow there), selects a rank-mult minor (`select_form`), requires point-independence, and runs
  `assert_genuine`.
- The rank of the local intersection form on the defect-s leaves ending at q equals mult = [v^s] m_q. This is the singular
  analogue of the Elias–Williamson local-intersection-form statement [? exact EKLP reference]. **Certification does not rely
  on it**: a rank shortfall is a visible failure, never a wrong answer.
- Downstream nothing changes:
  - `dechom_postcheck`, D14 scaling;
  - sindec `dec_one` (sindec.rs:610): its own form check at 2 points, π' = G⁻¹π;
  - `orthogonalize` (:717): unitriangularity is asserted, by degree, for any family;
  - `check_pi_iota` (:838);
  - `dec_composite` (:882), Kron-based.

**4.2 `closure(lp, p, None, outs)`.**
- Output exactly n = c_diag(b_p)(1) maps R → L_p (ins) or L_p → R (outs): the leaves of the n paths of L_p that end at the
  identity coset (the unique coset of the identity stratum).
- Their degree multiset is c_id(char L_p) = c_diag(b_p), since char L_p = b_p.
- `certify_closure` (:945) is unchanged: degree multiset, n×n identity-coset matrix of rank n at two points, genuineness.
  - The leaves are expected to be triangular in path order with nonzero diagonals (true for regular light leaves; singular case
    [?]). Certification checks the rank directly.
- No final rex move is needed for closure, because the endpoint expression is the empty ladder R.
- Rex moves are needed only inside the recursion (D-type steps).
- The stated sindec contract is met as is: `closure_gens` (sindec.rs:1007) re-runs `certify_closure` after any layer-A load.

**4.3 `full` / `end0`** (later; e_p, the B2 (2,2,2) gaps in HANDOFF_M4b).
- Double leaves LL^∨_{tgt,r} ∘ LL_{src,r} over the strata r, times R^λ monomials (Dot at level 0 on the bottom segments).
- `certify_full_plan` is unchanged: count = hom_gd, stacked rank at the detecting set, fresh residual.

**4.4 New certification, for leaf building blocks only.**
- Elem constants for every (kind, labels) used. Labels go up to 9 at (3,3,3): E2 (a,b) with a+b ≤ 9 and E3 (c;b) with c ≤ 9.
  Two routes:
  - A7-style comparison with 2-upright `hom_space` (`certify_elem`);
  - a polynomiality certificate of the closed form (`poly_certificate`, block diagonality included) plus nonvanishing, because
    these Hom spaces are 1-dimensional.
- Rex-move and fusion templates as hom_space or Added templates (`certify_added`).
- Nothing else changes: the family-level checks are the same functions and remain the only acceptance criteria. The M5b
  replay (a provider that refuses every call) works unchanged for Leaves output stored in layer A.

**4.5 Staging suggested by the data.**
- **L0, product rule** (closure): certifies the 5 product types with existing certified generators, existing nodes and
  commutation isos. No new mathematics, about 1 day.
- **L1, junction recursion** (closure, then pairs): leaves of A·B through a non-identity junction coset c, made of π_α of the
  2-upright factor and the closure of the word L_c·B. Targets: the 2 880, 4 320 and 7 200 types.
- **L2, full EKLP singular light leaves** with 3-upright rex moves: the rex-class pairs, the 2 700 types and the split class.
- **Acceptance set**: the M5b list of 12 closure types, 11 Exhausted pairs and 7 Memory pairs. The Refused pairs are optional,
  since some may clear with calibration on the parallel route.

## 5. Risks and open points
- **3-upright singular rex moves** are the only exact calls a leaves provider still needs. Their cost at (3,3,3) is unmeasured:
  the minimal degree-0 window between two reduced 5-6-rung words may itself be 10⁵⁺ unknowns. Options:
  - closed-form rex diagrams (EKLP (102)-(106) as cited by C2, not re-checked; a singular Matsumoto theorem would bound the moves);
  - chains of 2-upright fusion/switch maps;
  - by C2 M6, any degree-0 map with nonzero top scalar.
- **Word-agnostic types** (sindec redesign, not frozen): if each object carried its own reduced word, the rex class would
  disappear, because leaves give closure and decompositions for any word. The cost is that `dechom_precheck`, layer-A keys and
  DecData would be keyed by (q, word).
- **References.** EKLP numbers (Def 7.44 double-leaf shape, Thm 7.49 double leaves, Cor 3.7, Cor 4.17, Thm 4.16, eqs
  (44)/(45)) are taken from SYNTHESIS §5. The paper study must re-check them, together with the singular local-intersection-form
  and triangularity statements used in §4.
- **Pair classes.** The rex/split classification of the pairs rests on X flows = r(q) and on the target counts. One
  `KRHOM_TRACE` run per pair on the cloud would confirm it (`leaf_iso` messages vs `choose_window` / depth notes).
