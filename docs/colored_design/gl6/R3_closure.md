# R3: closure Φ_N at fixed N = 6 for 3-colored 2-, 3- and 4-strand braids

Design note, 2026-11-10. Read-only study (no builds, no runs of compiled code).
Inputs: dev `src/cclose.rs` (contract v2, exact oracle, reduced mode, Δ-tower adapter, G4), `src/xengine.rs` (closure_data, pruning, `closure_kills_map`), `src/schur.rs` (GramRow, c_diag), `tools/ntrunc_t1/src/main.rs`, SPEC_M4_M6a rev 2 §6.5-6.6, SYNTHESIS §1.4 / CL2 / CL3, NTRUNC_STUDY §3-4, sibling notes gl6/R4_m6a.md (§1.2 box points, §2-3 cost model) and gl6/R1_2strand.md.
Tags: [code] read in the source; [C] computed for this note (`gl6/py/r3_*.py`, < 1 s each, from `scripts/cb.pkl` and the ntrunc scan files); [M] measured earlier; [E] estimate.
Notation: k = number of strands (blocks, all labelled 3); box_k[w] = coefficient of t^w in [6 choose 3]_t^k; L_k(w) = Σ_{w' ≤ w} box_k[w']; "weight" = q-degree / 2.

## 0. Recommendation

1. **Open issue resolved: neither a convention mismatch nor a bug.** The "Gram weight" w(p) = l(p_max) − l(w_λ) is the *top* weight of B_p's Gram matrix, in the same unit as t in [6 3]_t. B_p's Gram entries have weights (d_a + d_b)/2, where d_a and d_b run over the exponents of c_diag(b_p), so they spread down from w. Entries of weight > 36 are 0 in A_λ = H*(Gr(3,6))^⊗4 by degree, so they are never computed or stored. The 5 surviving (3,3,3,3) types with w = 37, 39, 40 behave like the validated cases (3,3,3) at N = 4 and (2,2,2,2) at N = 3, 4, where surviving types also exceed the box top (§1).
2. **The closure at fixed N needs only G and the differential data modulo I_N.** An entry of weight w then has box_k[w] unknowns: at most 777 on 3 strands (48 reduced) and 13 656 on 4 strands (777 reduced). The R^λ e-monomial count is up to 88 233 on 3 strands and 1.3e8 on 4.
3. **Extraction at Grassmannian box points** (R4 §1.2, made concrete in §3), with three refinements:
   - streaming accumulators: memory equals the output size, and points can be processed in any order, in parallel, and resumed across sessions;
   - output in Schur coordinates: these are intrinsic, so restricting them gives every N' ≤ 6 and the reduced theory (ρ^0 = ∅);
   - residual coefficients serve as a polynomiality check.
4. **Ship the closed pairing P = G_y·F instead of F.** P = H^y_id · d_id · G^x_id uses identity-coset blocks only. It needs no solve with G'_id and does not depend on G'_id being invertible at structured points. cclose's rank stage already forms exactly this product (`pab`, cclose.rs:1336-1362).
5. **Certify the pruned d∘d by Freivalds on the closed complex**: check D_{h+1} Π^{-1} D_h u = 0 for random u, per (h, q), per N' and per mode. This is equivalent to `closure_kills_map` on every component. It costs two sparse mat-vecs per block, so the scan no longer has to collect `final_d2`. `closure_kills_map` stays as the test oracle (§5.3).
6. **Rank stage:** per-type state spaces, cached and independent of the knot. They are capped by the known MOY ranks (T0) and built with sparse echelon. The Δ-tower adapter does not apply to box data (§6).
7. **Volumes:**
   - 3 strands: the Gram library is 6 800 entries, 4.1e6 box coefficients unreduced (33 MB) and 1.8e5 reduced (1.4 MB). The P data are 2.6e5 (4_1) to 2.1e7 (median knot) entries [E], i.e. 54 MB to 4.4 GB reduced. Closure data are not the 3-strand bottleneck; scan memory per point is (R4 §3).
   - 4 strands (Conway): the per-entry bounds are fine, but the library is [E] 2e6-2e7 Gram entries (≈ 16 GB reduced, ≈ 0.8 TB unreduced) and every point is a full Conway scan. Out of reach with pointwise closure data (§7).
8. **Order:** reduced first (400 points on 3 strands), then unreduced (8 000 points, which contain the 400), or M6c for unreduced.

## 1. The open issue: Gram weight 40 vs box degree 36

### 1.1 What the code measures [code]
- `schur.rs:1928-1930` (`decat_scan`):
  - w = max exponent of c_diag(b_p), the coefficient of the identity table in b_p;
  - asserts w == l(p_max) − l(w_λ);
  - `unknowns = emono_count(λ, w)`, where e_j of block i has weight j (`schur.rs:1684`).
- `xengine.rs:1044-1058`: the closure generators' degrees (ins with the ladder shift, outs likewise) must equal the exponent multiset of c_diag(b_p).
- `cclose.rs:454`: every Gram entry G[b][a] is homogeneous, with 2·weight = ins[a] + outs[b]. `ERing` gives e_{i,j} weight j and q-degree 2j (`cclose.rs:146-192`).
- `cclose.rs:754` (`gauss_binom_coeffs`) and `ntrunc_t1/src/main.rs:171-196` (`gauss`, documented "in units of degree 2") index [N choose a]_t by weight, with top a(N − a) = 9 per Gr(3,6) block.
- Consequences:
  - one unit everywhere;
  - box top for k blocks = 9k (36 for k = 4);
  - an entry's weight is (d_a + d_b)/2, and w(p) is its maximum (d_a = d_b = top exponent).
- `ntrunc_t1/src/main.rs:740-755` ("gramBox"): takes box[w(p)] at the top weight only (`bs.get(w).unwrap_or(&0)`). So it reports 0 for w > top and never sees the lower entries; this is where the apparent conflict came from. The NTRUNC figure "at most 13 656 per entry" is the box maximum (w = 18) and is a correct upper bound.

### 1.2 Why entries above the box are harmless
- Φ_N(B) = Im(G ⊗ A_λ) uses G only through its image in A_λ (cclose module doc; R3 A3). Since (A_λ)_w = 0 for w > 9k, entries of weight 37-40 are exactly 0 at N = 6.
- Survival is not in tension with this: Φ_6(B_p) ≠ 0 iff G_p mod I_6 ≠ 0. A type kept by the kill rule (μ_1 ≤ 6) must therefore have nonzero entries of weight ≤ 36.
- The spread of entry weights:
  - in every pickled basis, all exponents of c_diag(b_p) are ≡ w (mod 2) and the maximum is w [C];
  - in S_4, c_diag(b_x) = v^{l(x)} P_{e,x}(v^{-2}): 3412 → v^2 + v^4 and 4231 → v^3 + v^5 [C];
  - both match c_diag(b_p) = v^w P_{w_λ, p_+}(v^{-2}), the KL polynomial of the longest double-coset element, whose lower terms produce the lower-weight entries.
- Validated cases where surviving types exceed the box top [C] (r3_weights.py, r3_cdiag.py; exponent: multiplicity):

| λ, N | top (red. top) | surviving with w > top | c_diag(b_p) of those types |
|---|---|---|---|
| (3,3,3), 4 | 9 (6) | 2, w = 10 | {2:1, 4:2, 6:3, 8:2, 10:1} (both) |
| (2,2,2,2), 3 | 8 (6) | 7, w = 9 and 12 | {3:1,5:3,7:3,9:1} ×4; {5:1,7:2,9:1} ×2; w = 12: {2:1,4:4,6:10,8:10,10:5,12:1} |
| (2,2,2,2), 4 | 16 (12) | 2, w = 17 | {7:2, 9:8, 11:13, 13:11, 15:5, 17:1} |
| (1^6), 2 | 6 (5) | 14, w ≤ 9 (e.g. 456123) | Temperley-Lieb survivors; engine2's Khovanov closures are validated on tables |
| (3,3,3,3), 6 | 36 (27) | 5: w = 40 (J − I), 39 ×2, 37 ×2 | not computed (needs the cloud basis) |

  The (2,2,2,2) and (3,3,3) rows are inside the 9 833-pair MOY validation of the kill rule. There, only the top entry or entries vanish.
- The 5 tables:
  - w = 40: J − I (all off-diagonal entries 1);
  - w = 39: 0,0,2,1/1,0,0,2/2,1,0,0/0,2,1,0 and 0,1,2,0/0,0,1,2/2,0,0,1/1,2,0,0;
  - w = 37: 0,1,1,1/1,0,1,1/2,0,0,1/0,2,1,0 and 0,1,2,0/1,0,0,2/1,1,0,1/1,1,1,0.
  - All have μ_1 = 6.
- A sharper fact [C] (r3_minexp.py): over every pickled basis and every N, max over surviving p of min exp c_diag(b_p) ≤ the *reduced* top. The bound is attained, for example (3,3,3) N = 6: 18 = 18 and (2,2,2,2) N = 4: 12 = 12.
  - This fits the theory: a (1,1)-tangle closure is a scalar times id on the irreducible Λ^3 C^6, so the reduced closure of a type vanishes iff the unreduced one does.
  - Prediction for (3,3,3,3) at N = 6: min exponent ≤ 27 for all 790 kept types.

### 1.3 Verdict and actions
- **Verdict:** entries above the box degree are forced to vanish. "Gram weight" is the top entry weight, not the weight of every entry. Units agree everywhere (no factor 2, no shift). The object q-shifts and global_q are separate and do not change which entries vanish.
- **Actions:**
  - (a) In the cloud (3,3,3,3) basis run (800 s / 3.6 GB once), print c_diag(b_p) of every kept type and assert min exponent ≤ 27.
  - (b) Change the ntrunc / `schur_scan` Gram statistic to the max and the sum over *entries* of box[(d_a + d_b)/2].
  - (c) Run T0 (Φ_6(B_p) = MOY) on the 5 high types once their generators exist.

## 2. What the closure needs at fixed N

- The closure of contract v2 reads G and F only through their images in A_λ (`nf`, `nf_prod`), and the rank stage reads only P = G_t·F.
- At N = 6, an entry of weight w therefore carries box_k[w] numbers. Every entry with w > 9·(#unmarked blocks) is 0. So is a generator row or column whose entries all exceed the top; it can be dropped.

| unmarked blocks | top | dim A | max box[w] (at w) | L(w) at w = 9 / 13 / 18 / 22 / 27 / 36 | R^λ e-monomials at w (for comparison) |
|---|---|---|---|---|---|
| 1 | 9 | 20 | 3 (3-6) | 20 | - |
| 2 | 18 | 400 | 48 (9) | 224 / 366 / 400 | (3,3): 158 at 9; 3 866 at 22 |
| 3 | 27 | 8 000 | 777 (13, 14) | 1 285 / 4 000 / 7 150 / 7 920 / 8 000 | (3,3,3): 28 237 at 18; 88 233 at 22 |
| 4 | 36 | 160 000 | 13 656 (18) | 4 982 / 26 551 / 86 828 / 123 877 / 157 122 / 160 000 | (3,3,3,3): 5.2e6 at 27; 1.3e8 at 40; 1.9e9 at 54 |

- SPEC §6.5's dense e-monomial interpolation ("needs CL2" above 2 500 unknowns) is therefore superseded for 3 and 4 strands at fixed N. It stays for 2 strands (§4.3).

## 3. Box extraction (CL2 restricted to one N)

**Nodes and points.**
- Block i gets random distinct nodes a^(i)_1..a^(i)_6 in F_p, and all 6k nodes are distinct, so every point is regular. Prime 2 uses a fresh sequence.
- Point z_μ, for μ = (μ^1..μ^k) with μ^i ⊆ 3×3: the roots of block i are (a^(i)_{μ^i_1+3}, a^(i)_{μ^i_2+2}, a^(i)_{μ^i_3+1}).
- Points do not depend on the knot, so one scan process handles a batch of knots at each point (R4 §2).

**Vanishing.**
- Definition: s_ρ(x | a) = det[(x_j | a)^{ρ_l+3−l}] / Δ(x), with (x | a)^m = Π_{t ≤ m} (x − a_t), multiplied over the blocks.
- Then s_ρ(z_μ | a) = 0 unless ρ ⊆ μ, and s_μ(z_μ | a) = Π (node differences) ≠ 0.
- So S_{μρ} = s_ρ(z_μ) = ⊗_i S^(i) (20 × 20 per block) is triangular for containment, and so is S^{-1}.

**Top coefficients are intrinsic.**
- Write f (homogeneous, weight w) as Σ_ρ c'_ρ s_ρ(x | a). Then c'_ρ = 0 for |ρ| > w, and for |ρ| = w, c'_ρ equals the Schur coefficient c_ρ of f.
- At box points only the terms with ρ ⊆ box survive. Hence f mod I_6 = Σ_{ρ ⊆ box, |ρ| = w} c_ρ s_ρ, with c_ρ = Σ_{μ ⊆ ρ} (S^{-1})_{ρμ} f(z_μ).
- Points needed: L(w), i.e. all μ with |μ| ≤ w (each lies below some ρ of size w when w ≤ top).

**Streaming accumulators (default).**
- When a value v arrives at z_μ: acc_ρ += (S^{-1})_{ρμ} v for every ρ ⊇ μ with |ρ| = w.
- Memory per entry = box[w] numbers.
- Linear in the values, so points can come in any order, from several threads or sessions. Partial accumulators are summed, and a checkpoint is (accumulators, done-points).
- Raw samples need not be kept. R4's per-point sample files become optional.

**Per-block sweep (alternative).**
- Store the values on L(w) and apply S^(i)^{-1} block by block (175 containment pairs per block).
- Use it where streaming is too slow, i.e. 4 strands unreduced.

| unmarked blocks | streaming: mult-adds per entry (max over w) | sweep: mult-adds per entry | memory per entry: streaming / sweep |
|---|---|---|---|
| 2 | ≤ 4.2e3 | ≤ 7.0e3 | ≤ 384 B / ≤ 3.2 KB |
| 3 | ≤ 6.1e5 | ≤ 2.1e5 | ≤ 6.2 KB / ≤ 63 KB |
| 4 | ≤ 9.3e7 | ≤ 5.6e6 | ≤ 109 KB / ≤ 1.28 MB |

**Nesting, reduced theory, output.**
- N' < 6: restrict to ρ^i ⊆ 3 × (N' − 3). These are the same intrinsic coefficients, so no new points are needed.
- Reduced (block 0 marked, as in cclose and xengine): keep ρ^0 = ∅. These need only the points with μ^0 = ∅ (400 on 3 strands, 8 000 on 4).
- Output: Schur coordinates per entry. They are converted per N' into cclose's e-monomial box basis (the per-block Schur-to-e matrix of `BlockAlg`). An alternative is a Schur-coordinate entry point in `GrAlg`, which already multiplies by Pieri in Schur coordinates.
- `ClosureDataC` gains `valid_n = 6`; `close` refuses N' > 6 (as xengine does for `prune_n`).

## 4. Gram and differential data

### 4.1 Gram (no scans; independent of the knot)
- At each point and for each type: G = H_id · G_id, where G_id and H_id are the n_p × n_p identity-coset blocks of the certified ins and outs. The identity-coset block of B_p has dimension c_diag(b_p)(1) = n_p.
- Cost per point: Σ n_p³ ≤ 3.7e5 for (3,3,3), plus evaluating the generator templates at the identity coset.
- The library is one-time per prime and shared by all knots: for (3,3,3) at N = 6 it is 40 types, 344 ins, 6 800 entries [C]:
  - entry weights 0-22, mean 12.4;
  - unreduced: Σ box = 4 094 072 coefficients, max 777 and mean 602 per entry; no entry exceeds 27;
  - reduced: Σ box = 178 897, max 48 and mean 26; 160 entries (weight > 18) vanish;
  - at prune 5 / 4: 3 808 / 217 entries.
- 4_1 uses 34 of the 40 types; 8_19 and the median knots (12a_859, 12n_674) use all 40 [C].

### 4.2 Differentials: the closed pairing P
- For a final edge x → y with component d: P[b][a] = (h^y_b ∘ d ∘ g^x_a)_id = (H^y_id · d_id · G^x_id)[b][a], of weight (d^y_b + e + d^x_a)/2. It is extracted like a Gram entry.
- Why P instead of F = G'^{-1}_id d_id G_id (SPEC §6.5 step 4):
  - the closed map is u ↦ G_y F u = P u, so the rank stage needs P only;
  - P needs no solve and does not need G'_id to be invertible at structured points;
  - P's higher weights make more of its entries vanish above the top.
- Costs:
  - P needs all L(top) points (R4's L(24) = 7 965 becomes 8 000): negligible;
  - contract v2 gains an edge variant carrying P;
  - F (or exact `solve_f`) stays on xengine-size cases, with the test P = G'·F.
- One scan per point yields all P values at once. Entries of weight > top are skipped.
- Points per prime: 400 reduced and 8 000 unreduced on 3 strands, 8 000 and 160 000 on 4 strands.

### 4.3 Two strands (3,3)
- Nothing is pruned at N = 6. Keep the exact path, which is independent of N: `rlambda_bases` (3.3 s for all 4 types) plus `solve_f` (R1). The same scan then also gives N = 7 (Wang 307).
- Box mode would need 309 / 20 points, at most 48 / 3 unknowns, and 20 Gram entries (435 / 48 coefficients) [C]. It is not needed.

## 5. Certification

### 5.1 Genuineness
- Every value at a box point is a value of a genuine bimodule map:
  - frozen library: M4a templates with polynomiality certificates, M4b-certified Hom bases;
  - eliminations pivot on degree-0 scalars, which do not depend on the point;
  - closure generators are certified by their c_diag degree multiset and by the determinant criterion at every point (assert, and redraw the nodes on failure).
- All final object lists and traces must be identical at all points (SPEC §6.5.1).

### 5.2 Polynomiality and homogeneity of the extracted entries
- (i) Entries of negative or non-integer weight are asserted 0 at every point (SPEC).
- (ii) **Residuals for free.** For an entry of weight w < top, also accumulate c'_ρ for |ρ| = w + 1 (box[w+1] more numbers). They must be 0. The shared point set L(top) always contains these points. Only entries of weight ≥ top lack residuals, and those exist only in the unreduced data.
- (iii) **Node independence.** The c_ρ are intrinsic:
  - re-extract a few random (entry, ρ) on a second node sequence, using only the down-set {μ ⊆ ρ} (40-200 extra scans per spot check);
  - equality certifies polynomiality where (ii) cannot;
  - in addition, two primes on two node sequences must give the same Poincaré polynomial (R4).
- (iv) The Freivalds D² test (§5.3) checks the P data at every weight. On unpruned runs D² = 0 holds exactly, so it is a strong consistency test.
- (v) Oracles:
  - T0: Φ_{N'}(B_p) = MOY per kept type at N' = 4, 5, 6, which checks the Gram data;
  - χ per knot = `decat_chi`; K_0 per stage;
  - nested pruning: prune-6 data closed at N' = 4, 5 must equal the prune-4 / prune-5 runs;
  - duality at N = 4;
  - Wang's T(2,n) through Markov words (136 at N = 6);
  - the Lewark-type relation [6 3]·P_red − P_unred = (1+t)Q with Q ≥ 0, where both closures exist.

### 5.3 Pruned d∘d closes to 0
- **Fact:** Φ is a functor, so D² = Φ(d∘d). D² = 0 iff Φ((d∘d)_{xw}) = 0 for every pair (x, w), and that is exactly what `closure_kills_map` tests (xengine.rs:1296: Poincaré polynomial of the two-object complex with and without the edge).
- **Production test (Freivalds):**
  - Coordinates: D_{h,q} goes from reps (object, generator, box element) to duals. For each (type, degree), the state space caches the invertible r × r matrix Π (duals of reps).
  - For random u, compute D_{h+1,q} Π^{-1} D_{h,q} u and require 0. Two vectors per (h, q), per N' and per mode give a false pass with probability ≤ 2/p per (h, q) block (≈ 1e-9 over F31). Each extra vector multiplies that bound by 1/p; use 3-4 vectors where a bound like 2^{-60} is wanted.
  - Cost: two sparse mat-vecs per block, negligible next to the ranks.
  - The scan stops collecting `final_d2`.
- **Deterministic variant, for M6a acceptance on small knots:**
  - at each point, P2_{xw} = H^w_id (Σ_y d^{yw}_id d^{xy}_id) G^x_id (identity-coset blocks compose);
  - extract it like P and require all top coefficients to be 0;
  - additionally, check that it agrees with `closure_kills_map` on xengine-size cases.
- **Negative tests:** perturb one P entry, or keep a killed composite; both checks must fail.

## 6. Rank stage at N = 6

- Keep `close_exact_colored`'s structure (state spaces, groups by (h, q), rank per group, parallel), with three changes:
  - (a) State spaces per (type, N', mode) are independent of the knot. They are cached on disk together with Π^{-1}.
  - (b) Rep selection stops at the known rank r_d of Φ_{N'}(B_p)_d, which the decategorified closure gives (T0). Duals are chosen among r_d rows only, so memory is r_d × cols instead of rows × cols.
  - (c) Sparse echelon throughout. `exact_mem_estimate` must be re-fitted per staged type (R4 risk R3).
- Sizes [C/E]:
  - rows per type and degree: ≤ n_p × max box, i.e. 54 × 777 ≈ 42 000 unreduced and 54 × 48 = 2 592 reduced on 3 strands; [E] ≈ 2.7e7 and 1.6e6 on 4 strands (n_p ≈ 2e3);
  - knot level (R4 §3): 4_1 reduced 1e5-1e6 chain dimensions (minutes); unreduced 2e6-2e7 (hours).
- **Not applicable:**
  - the (Δ_λ G, F) tower (`close_delta_tower`): it needs R^λ polynomials, and box data are taken mod I_N;
  - CL3: not needed for 3 strands at fixed N; possibly relevant for 4-strand state spaces (research).

## 7. Sizes per target

Per entry:

| case | points unred / red | unknowns per entry, max (mean) unred / red | accumulators per entry unred / red |
|---|---|---|---|
| (3,3) T(2,n) | exact path (box: 309 / 20) | 48 (22) / 3 (2.4) | 384 B / 24 B |
| 3 strands, all knots | 8 000 / 400 | 777 (602) / 48 (26) | ≤ 6.2 KB (4.8) / ≤ 384 B (210) |
| 4 strands (Conway) | 160 000 / 8 000 | 13 656 (~1.2e4 [E]) / 777 (~250 [E]) | ≤ 109 KB / ≤ 6.2 KB |

Per knot (N = 6, prune 6):

| knot | final objects / types [M] | Gram library | P entries [E] | P box data unred / red [E] | binding constraint |
|---|---|---|---|---|---|
| 4_1 | 276 / 34 | ≤ 6 800 entries; 33 MB / 1.4 MB | 2.6e5 | 1.2 GB / 54 MB | 8 000 / 400 scans of 0.5 GB each (R4) |
| 8_19 | 776 / 40 | same library | 7.2e5 | 3.5 GB / 150 MB | per-point scan 4.5 GB (R4) |
| median (12a_859) | 22 452 / 40 | same library | 2.1e7 | 100 GB / 4.4 GB | per-point scan ≈ 156 GB (R4): no |
| Conway (3,3,3,3) | [E] 1e5-2e5 / ≤ 790 | [E] 2e6-2e7 entries; ≈ 0.8 TB / 16 GB | [E] ≥ 1e9 | no | everything |

- The P-entry estimate is ≈ 930 entries per final object, uncertain by a factor 3 either way. It is scaled from SPEC's (2,2,2) 8_19 (1.2e4 entries for 287 objects) by the mean n_x·n_y per type (Σ n²/types: 7.7 → 170).
- The (3,3,3,3) library estimate (Σ n²):
  - growth factors 27 and 160 for a = 2 (k = 2 → 3 → 4) and 371 for a = 3 (2 → 3) suggest about 2 200 for a = 3 (3 → 4), i.e. 1.6e7 unpruned;
  - at N = 2a, (2,2,2,2) keeps 56%, giving ~8e6;
  - max n_p ≈ 1e3-3e3 [E].
- **4-strand verdict.** The box bounds the per-entry size, and extraction and certification carry over unchanged. But:
  - each of the 8 000 reduced points is a full Conway scan (untruncated peak 1.86e7 objects [M]), far beyond 15 GB per point;
  - Gram evaluation alone is Σ n_p³ ≈ 1e10 per point [E].
  - So 4-strand 3-colored gl_6 needs closed-form per-type closure data (D1 §9.3) and a non-pointwise scan representation (M6c/M7a). It is research, not an M6a extension.

## 8. Risks and open points

- **O1.** (3,3,3,3) c_diag of the 5 high types, and the prediction "min exponent ≤ 27": unverified until the cloud basis run (§1.3a).
- **O2.** No residuals exist for unreduced entries of weight ≥ 27 on 3 strands. They rest on (iii)-(v). The reduced data always have residuals below weight 18.
- **O3.** The P-contract needs a cclose change. Until it lands, F via the identity-coset solve also works at box points, with the determinant asserted.
- **O4.** Unreduced state spaces on 3 strands (42 000 rows per degree) need the rank-capped sparse selection of §6 before any unreduced N = 6 run.
- **O5.** Reduced closures of links depend on the marked component. Knots come first; for links, mark each component in turn (the points with μ^j = ∅ are inside the unreduced set).
- **O6.** The deformed-homology shortcut of SYNTHESIS note [D] does not apply: per-block disjoint nodes are an interpolation device, not a consistent deformation along a component.

## 9. Provenance (all in `gl6/py/`, `python3 -I`, < 1 s, no compiled code run)
- `r3_weights.py LAM N`: tables, Greene μ_1, w = l(p_max) − l(w_λ), weight histogram of the surviving types, and the types with w > box top (790 / 2 008 surviving, max 40; agrees with NTRUNC §3).
- `r3_cdiag.py`, `r3_minexp.py`, `r3_kl.py`: c_diag(b_p) from `scripts/cb.pkl` for (3,3), (2,2,2), (3,3,3), (2,2,2,2), (1^4), (2,1,2): exponent pattern, min exponent vs (reduced) top, KL check in S_4.
- `r3_counts.py`, `r3_scale.py`, `r3_cost.py`: box series, L(w), e-monomials, the Gram entry histogram (6 800 entries), Σ n and Σ n² scaling, extraction costs.
- `final333_N6.txt`: final LB and LB localized dimension at N = 6 for the 185 three-strand knots (median 22 452 objects), from `ntrunc/t1/out/scan_knots3*.txt`.
