# M4c A1 report: double-coset combinatorics, census G-L0, Stage 0

Date 2026-11-12. Branch `m4c-census` (from origin/dev 370f5b1). Scope: M4C_PLAN §4 milestone A1, §5 item 1, §9.

**Merge note.** origin/colored-m5 (f21e6e6, the newest library code: libbuild, M5c fixes, cost table
`costs_333_m5c.txt`) was merged into m4c-census at the start (commit e473db2). The merge was clean. All runs
below use this code and `tests/data/homprov/costs_333_m5c.txt`.

**Frozen files.** No edit to smap.rs, homprov.rs, sbim.rs, linsys.rs, tmpl.rs or sindec.rs. A1 adds new files
plus two `pub mod` lines in lib.rs. libbuild.rs is unchanged: the Stage 0 run uses its own subcommand
(`m4c stage0`).

## 1. What was built

| file | content |
|---|---|
| src/dcoset.rs (new) | pure combinatorics, described below |
| src/bin/m4c.rs (new) | `m4c census` (gate G-L0) and `m4c stage0` (Stage 0 run and replay) |
| src/leaves.rs (new) | `product_factorization`, `ProductClosure<F>`: a `HomProvider` wrapper around the production chain |
| tests/dcoset.rs (new) | gate G-L1: 6 tests |
| tests/leaves.rs (new) | Stage 0: factorizations, and certified product closure of the 5 types, F31 |
| docs/dev/colored_design/m4c/a1/census_closure.txt | full census output: summary table plus every level, junction word and step call |
| docs/dev/colored_design/m4c/a1/stage0_data/{f31,f32b} | Stage 0 layer A and templates (300 KB in all), with the run and replay logs |

Contents of src/dcoset.rs:
- `conf_type` and `rep_rank`.
- `step_terms`, `step_gf_table` and `step_gf`, with a `StepCache`:
  - G_{t'→t}(v) of one rung, written directly from SLL Def 3.5;
  - descending defect ℓ(t̄) − ℓ(ū);
  - ascending defect ℓ(t_) − ℓ(u_) − ℓ(LR t) + ℓ(LR u);
  - polynomial factor π⁺(LR t)/π⁺(LR u) as a product of Gaussian binomials;
  - this code does not share anything with `Schur::rung`.
- Deodhar recursion and junction analysis: `forward`, `reach_back`, `junction_types` (per level: forward/backward
  generating functions and degree windows).
- Permutation and Bruhat layer: `min_perm` (p_ in S_K), `perm_len`, `perm_comp`, `perm_inv`, `perm_bruhat_le`,
  `bruhat_le` (p ≤ q iff p_ ≤ q_), `table_of_perm`.
- Reducedness and junction words:
  - `reduced_for`: c_t(char W) is exactly 1, and every table of char W is ≤ t;
  - `junction_word`: prefix extension W_{t'}·ρ first (C_m4c m3), then the minimal-flow canonical ladder, then a
    bounded BFS;
  - `bfs_reduced` and `reduced_words`.
- Cells and Grassmannian reduction:
  - `cells`, `lr`, `rr`, `psi` (CORE ψ), `psi_weighted_len`;
  - `grass`: SLL Thm 2.6 at rung level. It returns Q, L = RR(Q), z, m', m'' and the local bottom.

**Stage 0 (`leaves::ProductClosure`).**
- **Factorization search.** For a closure call, it looks for a factorization L_p ~ A·B with these properties:
  - the commutation class is searched with the swap path recorded;
  - A and B are indecomposable non-identity (λ,λ) types a and b;
  - c_diag(p) = c_diag(a)·c_diag(b).
- **Generators.**
  - ins: iso ∘ Comp(Kron(α_i, B), β_j);
  - outs: Comp(β'_j, Kron(α'_i, B)) ∘ iso;
  - α and β come from the inner provider (the production chain) on the canonical ladders of a and b;
  - iso is a chain of `Local` nodes on degree-0 `hom_space` templates of the 2-rung commutation windows. Each
    window must be 1-dimensional, otherwise the call fails.
- **Acceptance** happens only through `certify_closure`: inside the provider, again in sindec after the call, and
  on every layer-A load.
- **Order.** By default the inner provider goes first (Leaves last, C_m4c M4); with `product_first` the product
  rule goes first.
- **Arena.** mark/commit is respected: a failed product call commits nothing.
- **Other calls.** `full`, `pairs` and `end0` are delegated to the inner provider.
- **Stamp suffix:** `stage0=ProductClosure(v1,inner-first|product-first)`.

## 2. Gate G-L1: passed

`KRHOM_THREADS=4 cargo test --release --test dcoset`: 6/6 in 3.2 s. The tests cover every block (λ, μ), with every
μ a composition of K with |λ| parts, for λ = (3,3), (2,1,2), (2,2,2) and (3,3,3).

| check | (3,3) | (2,1,2) | (2,2,2) | (3,3,3) |
|---|---|---|---|---|
| Deodhar sums (independent G) = `ladder_char`: canonical ladders (minflow and BFS) / words incl. every one-rung extension | 19 / 133 | 130 / 986 | 269 / 2 416 | 1 446 / 18 764 |
| flows per coset = c_t(1) (`sexpr::deodhar_check`), ladders | 19 | 130 | 269 | 1 446 |
| Lemma 8.7: Σ_t fwd·back on identity paths = c_diag(b_p) at every level, closure types | 4 | 11 | 21 | 54 (all but the decomposable 0,1,2/1,1,1/2,1,0) |
| Lemma 8.7: one-rung extensions of indecomposable L_p, Deodhar sum = Σ_q term_mult_q·b_q | 96 | 720 | 1 728 | 11 874 |
| tables: ℓ(min_perm) = lmin, table_of_perm round trip, weighted ℓ(ψ) = lmin, Bruhat-below ⇒ smaller lmin/lmax | 16 | 108 | 216 | 1 000 |
| every canonical ladder is a reduced word of its type (c_t = 1 and support Bruhat-below t) | all | all | all | all |
| Grassmannian reduction (rung level): t_ = Q_·m_ with lengths adding, same z, **G_local = G_global** | 168 steps | 1 168 | 3 136 | 28 000 |

**Degree sign.**
- E2 `Eta` on (a,b) ∈ {(1,1), (2,1), (1,2), (2,2), (3,3), (3,6)} has krhom degree +ab.
- The identity path of E0.b F0.b has defect +ab.
- So the EKLP defect equals krhom's normalized degree, with no sign flip. This agrees with `certify_closure`, which
  accepts closure maps of degree multiset c_diag = Σ_paths v^def.

**Two [K] items of D1 are now verified combinatorially** for every step of these four λ:
- the rung-level combination of SLL Thm 2.6 with Lemma 3.9 (G_local = G_global; M4C_PLAN R7);
- the normalization of degrees.

## 3. Gate G-L0: census of the closure failure list

**Command.** `m4c census` takes 0.8 s and writes docs/dev/colored_design/m4c/a1/census_closure.txt.

**Failure list.** It is the M5b list: N = 6, 12 types; N = 5, 1,1,1/1,1,1/1,1,1. The M5b runs fail identically
for both primes. colored-m5 has no "after M5c" refresh of the list yet, so it is unchanged.

**Replica check (hu.py) reproduced exactly** by `homprov::hom_unknowns` (the run asserts it):

| call | unknowns |
|---|---|
| Hom^18(∅, E0.3 F0.3 E1.3 F1.3) | 852 011 |
| Hom^13(E0.2 F0.2 E1.3 F1.3, ∅) | 392 995 |
| L_p of 1,1,1/1,1,1/1,1,1 at d = 4 | 224 249 |
| L_p of 1,1,1/1,1,1/1,1,1 at d = 6 | 406 634 |

All four agree in both directions.

**What the census plans.** Lemma A on 3 uprights (D1 §1.3), restricted to the identity paths of L_p:
- junction types per level are `junction_types(L_p, {id})`;
- junction words come from `junction_word`, with prefix extension first;
- one step call is made per junction edge t' → t with G_{t'→t} ≠ 0:
  - `Id` when W_t = W_{t'}·ρ;
  - otherwise the exact call Hom^d(W_{t'}·ρ, W_t) for every d in supp G, reported as the maximum over d and over
    both directions.
- Two cheaper realizations are listed for comparison:
  - the minimal differing window (`win_split`; complete only if the common prefix is the chart's Z part);
  - the K window of the Grassmannian reduction: local canonical words on the K cells of Q. "needs pad" means the K
    column of Q has 3 cells, so a 4th upright would be needed (D1 §1.7, C_m4c C1).

| N | type p | L_p flows | n | product (Stage 0) | M5b whole-ladder call at d0 | junctions/level (max) | max W_t flows (rungs) | steps Id / exact | max step unknowns | > 10⁴ | > 10⁵ | max window | max K-window (needs pad) | words prefix/canonical/BFS | not reduced |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 5 | 1,1,1/1,1,1/1,1,1 | 4 320 | 54 | no | 224 249 | 2 | 36 (3) | 4 / 5 | 88 301 | 1 | 0 | 7 225 | 214 (2) | 4/4/0 | - |
| 6 | 0,0,3/3,0,0/0,3,0 | 400 | 1 | 0,3,0/3,0,0/0,0,3 · 3,0,0/0,0,3/0,3,0 | 852 011 | 1 | 1 (1) | 2 / 2 | 3 184 | 0 | 0 | 750 | 750 (0) | 2/2/0 | - |
| 6 | 0,3,0/0,0,3/3,0,0 | 400 | 1 | 3,0,0/0,0,3/0,3,0 · 0,3,0/3,0,0/0,0,3 (1 swap) | 852 011 | 1 | 1 (3) | 2 / 2 | 3 184 | 0 | 0 | 3 184 | 750 (1) | 2/2/0 | - |
| 6 | 1,0,2/2,0,1/0,3,0 | 600 | 3 | 1,2,0/2,1,0/0,0,3 · 3,0,0/0,0,3/0,3,0 | 392 995 | 1 | 3 (1) | 2 / 2 | 3 184 | 0 | 0 | 787 | 787 (0) | 2/2/0 | - |
| 6 | 1,2,0/0,0,3/2,1,0 | 600 | 3 | 3,0,0/0,0,3/0,3,0 · 1,2,0/2,1,0/0,0,3 (1 swap) | 392 995 | 1 | 3 (2) | 1 / 3 | 2 983 | 0 | 0 | 2 983 | 787 (1) | 1/3/0 | - |
| 6 | 0,3,0/1,0,2/2,0,1 | 600 | 3 | 3,0,0/0,1,2/0,2,1 · 0,3,0/3,0,0/0,0,3 (1 swap) | 392 995 | 1 | 3 (2) | 1 / 3 | 5 104 | 0 | 0 | 5 104 | 750 (1) | 1/3/0 | - |
| 6 | 1,0,2/1,2,0/1,1,1 | 2 700 | 18 | no | 701 771 | 1 | 9 (2) | 2 / 3 | 13 765 | 1 | 0 | 13 765 | 787 (1) | 2/3/0 | - |
| 6 | 1,1,1/0,2,1/2,0,1 | 2 700 | 18 | no | 701 771 | 1 | 9 (2) | 2 / 3 | **101 507** | 1 | **1** | 2 401 | 787 (1) | 2/3/0 | - |
| 6 | 0,2,1/2,0,1/1,1,1 | 2 880 | 18 | no | 735 960 | 2 | 12 (3) | 3 / 5 | 14 815 | 1 | 0 | 787 | 787 (3) | 3/4/0 | - |
| 6 | 1,1,1/1,0,2/1,2,0 | 2 880 | 18 | no | 735 960 | 2 | 12 (3) | 4 / 5 | 14 815 | 1 | 0 | 1 220 | 127 (2) | 4/4/0 | - |
| 6 | 0,1,2/2,0,1/1,2,0 | 7 200 | 27 | no | 2 815 899 | 2 | 12 (3) | 3 / 6 | 14 815 | 1 | 0 | 3 191 | 3 191 (2) | 3/5/0 | - |
| 6 | 0,2,1/1,0,2/2,1,0 | 7 200 | 27 | no | 2 815 899 | 2 | 12 (3) | 4 / 5 | **111 497** | 1 | **1** | 111 497 | 787 (2) | 4/4/0 | - |
| 6 | 1,1,1/1,1,1/1,1,1 | 4 320 | 54 | no | 224 249 | 2 | 36 (3) | 4 / 5 | 88 301 | 1 | 0 | 7 225 | 214 (2) | 4/4/0 | - |

**Highlights.**
- **The Lemma A plan is small.**
  - At most 2 junction types per level.
  - Junction words have at most 3 rungs and at most 36 flows. The largest word in the plan is L_p's own prefix.
  - Over all 13 jobs there are **27 distinct exact step calls**:

    | unknowns | calls |
    |---|---|
    | ≤ 10³ | 7 |
    | 10³-10⁴ | 14 |
    | 10⁴-10⁵ | 4: 13 765, 14 815 (twice), 88 301 |
    | > 10⁵ | 2: 101 507 and 111 497 |

  - This replaces whole-ladder calls of 0.22-2.8 M unknowns.
- **Steps above the 10⁵ budget:**
  - **1,1,1/0,2,1/2,0,1**, level 3: Hom^{6,8,10}(F0.1 E1.2 F1.2, F0.1), 101 507 unknowns.
    - Its differing window is the closure of E1.2 F1.2 behind the common prefix F0.1 (2 401 unknowns).
    - Its K window needs a pad.
    - No other reduced word for 2,1,0/0,3,0/0,0,3 exists with ≤ 5 rungs: F0.1 is the only one.
  - **0,2,1/1,0,2/2,1,0**, level 4: Hom^{3,5,7,9}(E0.1 F0.1 E1.3 E0.2, E0.2 E1.3), 111 497 unknowns.
    - It has no common prefix or suffix.
    - Its K column holds all 9 roots: C_m4c C1's "third column empty" case, which needs a pad.
    - All 9 reduced words of 3,0,0/2,1,0/0,3,0 with ≤ 5 rungs were tried for the target; none is cheaper.
- **Steps between 10⁴ and 10⁵:**
  - 88 301 for both 1,1,1/1,1,1/1,1,1 jobs; its window is 7 225.
  - 13 765 and 14 815, in 4 jobs.
- **Every type has a one-flow (reduced) junction word** (C_m4c m2 / M4C_PLAN "types without a one-flow word"):
  - every junction word chosen is a prefix extension (34) or the canonical ladder (43, counting the empty word at identity levels);
  - the BFS fallback was never needed;
  - by G-L1, every canonical ladder of every (3,3,3) block is reduced in the strict sense, including the
    decomposable 0,1,2/1,1,1/2,1,0. **The m2 concern does not materialize.**
- **Product types.** The census finds exactly the 5 product types of L3 §3.2. Three need one commutation swap and two
  need none. The 7 others are not products, which tests/leaves.rs asserts.

## 4. Stage 0: 5/5 for F31 and F32b, replay 5/5

**Runs.**
- Command: `m4c stage0 --field 31|32b --product-first --cache DIR --costs tests/data/homprov/costs_333_m5c.txt`, then
  the same with `--replay` (a provider that refuses every call).
- The SIndec instances are the production ones (prune_n 6, layer A in DIR, stamp = chain stamp + Stage 0 suffix).
- Every certification below is `certify_closure` in `SIndec::closure_gens`, after the provider's own check.

| type p | factors a · b | swaps | n | F31 ins/outs | F32b ins/outs | wall F31 / F32b | replay F31 / F32b |
|---|---|---|---|---|---|---|---|
| 0,0,3/3,0,0/0,3,0 | 0,3,0/3,0,0/0,0,3 · 3,0,0/0,0,3/0,3,0 | 0 | 1 | 1/1 | 1/1 | 0.07 / 0.10 s | ok / ok |
| 0,3,0/0,0,3/3,0,0 | 3,0,0/0,0,3/0,3,0 · 0,3,0/3,0,0/0,0,3 | 1 | 1 | 1/1 | 1/1 | 0.29 / 0.26 s | ok / ok |
| 1,0,2/2,0,1/0,3,0 | 1,2,0/2,1,0/0,0,3 · 3,0,0/0,0,3/0,3,0 | 0 | 3 | 3/3 | 3/3 | 0.01 / 0.01 s | ok / ok |
| 1,2,0/0,0,3/2,1,0 | 3,0,0/0,0,3/0,3,0 · 1,2,0/2,1,0/0,0,3 | 1 | 3 | 3/3 | 3/3 | 4.1 / 3.6 s | ok / ok |
| 0,3,0/1,0,2/2,0,1 | 3,0,0/0,1,2/0,2,1 · 0,3,0/3,0,0/0,0,3 | 1 | 3 | 3/3 | 3/3 | 0.05 / 0.05 s | ok / ok |

**Result: 5/5 certified for each prime; replay 5/5 for each prime, every map loaded from layer A** (20 files per prime,
in a/stage0_data).

The degree multisets are v¹⁸ (n = 1) and v¹³ + v¹⁵ + v¹⁷ (n = 3). These are 5 of the 12 N = 6 closure gaps (all 5
are universal blockers in M5c_coverage.md). The commutation windows are 1-dimensional Hom⁰ templates.

**Inner-first mode** (the default order, "Leaves last"), F31: the 5 types are again certified through the product
rule. Their layer-A files contain the `K` (Kron) records. The production chain, however, first spends:
- 281 s on 1,0,2/2,0,1/0,3,0;
- 190 s on 1,2,0/0,0,3/2,1,0;
- 181 s on 0,3,0/1,0,2/2,0,1.

Presumably this is the Exact d = 13 whole-ladder call, which the M5c records admit, before the chain fails at d = 15/17. The log does not separate the inner tiers, so this is not verified. For
the A3 wiring, put product closure before Exact for product types, or skip the chain when a factorization exists.
The maps are certified either way.

**Persistence.** data-gl6 is not written to in A1. The Stage 0 layer A plus templates (both primes, with the run and
replay logs) are committed in docs/dev/colored_design/m4c/a1/stage0_data. They must be merged into the data-gl6 layer A
with the A3 stamp change.

## 5. Implications for A2/A3

1. **Closure through Lemma A on 3 uprights (ThreeUp) is viable** for the 7 non-product types and the N = 5 job.
   - Junction words are tiny and step calls are 10²-10⁴ unknowns, except:
     - 88 301 (a 7 225 window);
     - the two over-budget steps above.
   - The M4C_PLAN §5 decision rule ("ThreeUp if max ≤ 10⁴, else Factored") formally says Factored. In fact only 6
     distinct calls exceed 10⁴:
     - 3 at 1.4-1.5·10⁴, inside the measured cheap band (≤ 48 s, 0.5 GB);
     - 1 at 88 301;
     - 2 at 1.0-1.1·10⁵.
   - **Recommendation:** keep ThreeUp in A2/A3, with three conditions:
     - (a) use the differing window id_P ⊗ s when W_{t'}·ρ and W_t share a prefix. This applies to the 101 507 step
       (window 2 401) and the 88 301 step (window 7 225). Whether these window maps have an invertible top-to-top
       matrix is checked by Lemma A's rank diagnostic and, finally, by `certify_closure`.
     - (b) measure the single remaining 111 497 call (0,2,1/1,0,2/2,1,0, d = 3..9) on the cloud. It needs an
       uncensored exact-shape cost record to fall into the 10⁵-4·10⁵ band of the budget rule. If it is censored,
       that one type needs the B-phase closed-form steps or a padded K window.
     - (c) greedy selection over all of supp G (C m4).
2. **Expected A3 cost.** About 30 distinct exact calls per prime, mostly < 10 s. The worst are 1-2 calls of 10⁵
   (measured 10⁵-band calls: up to 420 s). The cost is minutes per prime, not hours (D1 §6 estimated 1-3 h).
3. **Stage 0 is done.** Production wiring (A3) needs `ProductClosure` around the chain in
   `sindec::provider_for_stamped` (a frozen file; scheduled in A3), with product-first for factorizable closure
   calls. That gives one stamp change, coordinated with M5c (M4C_PLAN §6).
4. **K windows.** 18 of the 49 exact steps (counted per job) need a pad (the K column of Q has 3 cells). The remaining K windows are
   ≤ 3 191 unknowns. This confirms C_m4c C1 and M6: the Factored chart cannot run without `Pad` for these steps. For
   closure, ThreeUp avoids it.
5. **Theory.** G-L1 verified, for every step of (3,3,3), the two [K] items that matter for A2-B5:
   - G_local = G_global under the rung-level Grassmannian reduction;
   - the degree sign.

## 6. Open problems

1. Two step calls at 1.0-1.1·10⁵ unknowns (§3). One has a 2 401-unknown window. For the other (111 497), measurement
   or a B-phase closed form is needed. Neither can be avoided by choosing a different reduced target word with
   ≤ 5 rungs. Varying the source junction word W_{t'} jointly with the target word was not explored.
2. Window maps (id_P ⊗ s) are not a priori complete as step families. Lemma A's top-to-top rank test (A2) decides
   per step, and any shortfall is a visible `Exhausted`, never a wrong map.
3. The census covers closure jobs only. A pairs census (G-P0: chains, pads, braid windows) is milestone B2.
4. The M5c failure list has not been refreshed on colored-m5. If M5c certifies some of the 13 closure jobs, the A3
   list shrinks. Stage 0 stays an independent route.
5. Inner-first order wastes about 3-5 min per product type in the Exact chain (§4); A3 wiring should avoid it.

## 7. Reproduce
```
cargo build --release --bin m4c
./target/release/m4c census --out docs/dev/colored_design/m4c/a1/census_closure.txt
./target/release/m4c stage0 --field 31  --product-first --cache C31 --costs tests/data/homprov/costs_333_m5c.txt
./target/release/m4c stage0 --field 31  --product-first --cache C31 --costs tests/data/homprov/costs_333_m5c.txt --replay
(same for --field 32b)
KRHOM_THREADS=4 cargo test --release --test dcoset --test leaves
```
