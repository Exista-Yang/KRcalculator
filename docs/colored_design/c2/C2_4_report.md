# C2-4 report: M5i (idempotents of decomposable canonical ladders) for (2,2,2,2) at N = 4 (branch colored-2c4-m5i)

Date 2026-11-14. Machine: 4 cores, 15 GB. Base: origin/dev 710c102, merged with origin/colored-2c4-fix 4553f23
(C2-5: Leaves routes, `sbox --provider`, calibrated table costs_2222_c25.txt, repaired layer A). Tags: [M] measured,
[E] estimate. Outputs: runs_c2/k4m/, tests/m5i.rs, examples/m5i_*.rs, tools/c2/k4m_status.py; library delta on
data-gl6 under c2/l2222_n4_c24/.

## Summary

- **M5i is implemented end to end** (sindec, sengine, swire; no other source file changed). Every decomposable
  canonical ladder L_p gets an idempotent e_p with image ≅ B_p; decompositions accept decomposable sources and
  targets; closure generators go through e_p; the engine works in reduced coordinates of Im e_p.
- **e_p for all 7 kept decomposable types of (2,2,2,2) at N = 4** [M], checked (e² = e and KL ranks at every table
  stratum at 2 points, e² = e on every coset and Σ ranks = rank B_p at a third, top scalar 1). Constructions: 2
  by peeling the extra summand with certified provider pairs (one through the Leaves Lemma A route, one through an
  Exact fallback of 112 k unknowns), 5 by **prefix propagation** e_p = Kron(e_{p'}, R) (the canonical ladder is the
  decomposable canonical ladder of another type p' followed by one rung with b_{p'}·R = b_p). Exact End^0
  reference where affordable: (3,1,3) type 2,0,1/0,1,0/1,0,2 and the (2,2,2,2) type of block (2,4,0,2) (End^0 =
  k², 13 k and 196 k unknowns): e_p is the unique End^0 idempotent (in the span, dim 2) — PASS.
- **NeedsIdempotent jobs** of C2-2 (43: 40 pairs + 3 closure types): **31 / 43 certified for F31 and for F32b**
  (identical sets), replay (provider refusing every call) reproduces them (section 3). The 12 left fail on a
  provider call that is not an M5i question any more (a summand that splits non-locally inside a 1 300-7 800-flow
  X, Exact 5-60 M unknowns; section 4).
- **Engine validation** [M]: (a) forced non-canonical decomposable ladders (`KRHOM_FORCE_LADDER`, 2 sets of 2-3 types
  of (2,2,2), m = 18-216 reduced ranks, one with extra summands in degrees ±1) on 4_1, 5_2 and 3_1 at (2,2,2):
  **all 6 closures identical** to the canonical-ladder results (C2-3), both primes; (b) the 11 C2-2 4-strand words
  with the M5i binary on the C2-5 library: **identical, closure by closure, both primes**; (c) χ = decat_chi and
  the D14 scaling check in every run.
- **4-strand knots: none computed.** With M5i, every class-A knot of knots4_cloud.txt (6_1, 7_x, 8_x, 9_x; all
  share the same 450 + 93 jobs) misses exactly **6 jobs**: the 4 non-M5i pairs 1235, 1444, 1470, 1894 that C2-5 could
  not repair, and the 2 M5i pairs 1769, 1812 (decomposable targets reached only through the same kind of
  non-local splitting). Every 4-strand word that visits an idempotent type also needs 1235, 1444, 1470 (checked on
  268 stabilized 3-strand knots), so no 4-strand run through e_p at (2,2,2,2) was possible this session; the engine
  path through e_p is validated at (2,2,2) instead (a).
- Full suite `KRHOM_THREADS=4 cargo test --release`: @@SUITE@@

## 1. Design as built

### 1.1 e_p (sindec `ty` → `ty_full` → `build_idem`)

`ty(p)` (kept types) and `ty_full(p)` (any type, also killed ones needed as extra summands) intern L_p = the
minimal-flow canonical ladder. If char L_p ≠ b_p:

1. **Prefix propagation** (`idem_by_prefix`): if L_p = L_{p'}·R with L_{p'} the canonical ladder of a decomposable type
   p' and term_mult(p', R) = b_p, then e_p = Kron(e_{p'}, R) (image B_{p'} ⊗ R ≅ B_p).
2. Otherwise **peeling** (C2_2_report §2 step 1): for every extra summand (r, s, cnt) of char L_p, certified pairs
   (ι_r: L_r → L_p of degree s, π_r) from the provider (layer A key `kind=idem`, e_r passed if L_r is itself
   decomposable), form at rep(r) at 2 points (equal), normalized by its inverse, orthogonalized by the Neumann series
   (degrees ±1 allowed), e_p = id − Σ ι_r π_r.
3. Checks (`idem_checks`, spec §5.2 step 6): `check_idempotent` at 2 points (e² = e and rank = c_t(b_p)(1) at
   rep(t) of every table t of char L_p, top scalar 1); e² = e on every coset at a third point and Σ_c rank e_c =
   rank B_p. Π∘Ι = δ for the peeled summands (e_r on the diagonal).

The spec's End^0 + Krylov route is not used (End^0 of the 7 types is 0.2 M-2.9 G Exact unknowns); the End^0
reference is a test (`examples/m5i_check --end0`, tests/m5i.rs): when char L_p = b_p + Σ b_r with distinct r in
degree 0, End^0(L_p) ≅ k^{1+#r}, so the idempotent with top scalar 1 and the ranks of B_p is unique and e_p must lie
in the span of the End^0 basis.

**Exact fallback.** An extra-summand splitting is global (the singular analogue of B_{s2}B_{s1}B_{s3}B_{s2} =
B_{s2s1s3s2} ⊕ B_{s2}): Search's windows never see it (stuck branches). For M5i calls only (kind=idem, an e_q target,
a decomposable source, closure through e_p) a failed chain call is retried by Exact with raised caps
(`KRHOM_M5I_EXACT_UNK`, default 1.3 M unknowns; never under the replay provider). Its output passes the same
DecHom output checks and M5's own re-certification.

### 1.2 Decompositions (sindec `dec_one`, `dec_composite`, `reroute_family`)

- **Source idempotent** E = Kron(e_p, R): the targets are the kept (q, s, cnt) of term_mult(p, R); the provider sees
  X = L_p·R with the multiplicity mx of b_q⟨s⟩ in char X (≥ cnt); ι ↦ E∘ι, π ↦ π∘E; cnt independent columns I and
  rows J of the projected form (rank must be cnt), π' = G[J][I]⁻¹ π_J. Checks: Π∘Ι = δ, Σ ιπ = E (not id) on every
  coset at 2 points; rank bookkeeping uses term_mult ranks.
- **Commutation reroute** (new; the provider cannot decompose L_p·R for the same reason as above): R commutes with
  the tail w[j..] of w = L_p (homprov's `rungs_commute`: distant rungs, or adjacent rungs of opposite direction), so
  X ≅ Y = w[..j]·R·w[j..] through degree-0 window isomorphisms (`swap_iso`: one-dimensional Hom^0 templates, g∘f
  normalized to id); w[..j] is the canonical ladder of a type a, and Y is decomposed rung by rung from a by the
  ordinary `dec_word(a, R·w[j..])`, which never meets L_p. The family (c⁻¹∘ι_Y, π_Y∘c) is a basis of each
  multiplicity space of X with Π∘Ι = δ; the source-idempotent selection above then applies. Used for every
  decomposable source when R commutes past the last rung, and as a fallback for any target the provider cannot reach
  (cycle guard). Killed summands are dropped consistently (the kill set is a tensor ideal).
- **Decomposable targets**: e_q is passed to the provider (it wraps ι∘e_q, e_q∘π); shortcuts: X = L_q (ι = π = e_q),
  X commutation-equivalent to L_q (ι = c⁻¹∘e_q), X = L_{xq} for a decomposable xq (its extra-summand pairs reused
  from layer A kind=idem). Π∘Ι checks compare with e_q on the diagonal.
- dec_word(p, []) of an idempotent type is (e_p, e_p).

### 1.3 Closure, layer A, frozen data

- `closure_gens` passes e_p (the provider wraps every generator as e_p∘g / g∘e_p and certifies rank n on the
  identity coset); generators loaded from layer A are re-wrapped with this builder's e_p (an imported copy is a
  different node) and re-certified. Exact fallback as for pairs.
- Layer-A keys: unchanged for every map without an idempotent; `eq=m5i-peel1` / `ep=m5i-peel1` (M5I_STAMP) for maps
  that contain one; `kind=idem` for the extra-summand pairs. The C2-2 / C2-5 layer A is reused as it is.
- `DecData::reduced(ev, p, c)`: identity for e_p = None; else (P, Q) = smap::reduced_of(block(e_p, c)) (pivot
  columns of e_c: P = E[:, J], Q = rref rows; QP = id_m, PQ = E_c asserted).
- `KRHOM_FORCE_LADDER="TABLE=WORD;..."` replaces canonical ladders (validation only); `KRHOM_NO_M5I` restores the
  M5a refusal; `KRHOM_M5I_NO_PREFIX`, `KRHOM_M5I_NO_REROUTE`, `KRHOM_NO_REROUTE_FALLBACK` switch the new routes off.

### 1.4 Engine (sengine `lib_shape`, swire)

- `lib_shape`: for an idempotent type the shape is the rank of e_p per coset (cosets of rank 0 dropped), read at
  two fixed generic points that must agree; rep(p) has m = 1 (asserted by `TypeShape::new`).
- swire `RedCache` (per scan point): (P, Q) per (type, coset). Every library block that ends at an idempotent type
  is taken to reduced coordinates: ι of a summand q → ι·P_q, π → Q_q·π (`assemble`), Rickard components Q_t·b·P_s
  (`lib_map`), closure generators Q_p·g, h·P_p (`extract_stage`, `generator_pass`; `gram_at_point` already did).
- d_C ⊗ id: a component f: a → b of the previous complex touching an idempotent type is taken back to flow
  coordinates P_b·f·Q_a (`to_flows`) before the Kron step; since ι = E∘ι and π = π∘E, Kron(P_b f Q_a) gives the same
  blocks as Kron(e_b f' e_a). The reference path `KRHOM_DCID_REF` refuses idempotent libraries.
- Point dependence: P, Q come from an rref of e_c at each point, i.e. a point-dependent basis of Im e_c. Every
  engine quantity is invariant under it: the sparsity pattern is per component (a block is zero iff its transform
  is), pivots are c·I between equal types (unchanged by a common change of basis), and the closure entries
  h·D·g = h e D e g do not depend on the basis. The frozen trace (record at one point, replay at all box points,
  pivot values included) passed in every run.
- The WLib hash includes the idempotent types and M5I_STAMP; `[swire] M5i: k idempotent types` is logged.

## 2. e_p of the 7 kept decomposable types of (2,2,2,2), N = 4 [M]

(examples/m5i_check; F31; ladders and characters in runs_c2/k4/idem_2222.txt and runs_c2/k4m/... logs)

| type (T) | block | flows | rank B_p | construction | time |
|---|---|---|---|---|---|
| T1 1,0,1,0/0,1,0,1/1,0,1,0/0,1,0,1 | (2,2,2,2) | 1 296 | 1 290 | peel b_{2,0,0,0/0,0,2,0/0,2,0,0/0,0,0,2}: Exact fallback 112 k unknowns (default chain) | 52 s |
| T2 1,0,1,0/0,0,1,1/1,1,0,0/0,1,0,1 | (2,2,2,2) | 576 | 516 | Kron(e_T5, F1.2) | 10 s |
| T3 0,1,0,1/1,0,1,0/0,1,0,1/1,0,1,0 | (2,2,2,2) | 7 776 | – | peel through the Leaves Lemma A pair route (m4c-lf provider, 18 s); default chain: Search exhausted, Exact 9.1 M not admitted | – |
| T4 1,0,1,0/0,1,0,1/1,1,0,0/0,1,0,1 | (2,3,1,2) | 384 | 344 | Kron(e_T5, F1.1) | < 1 s |
| T5 1,1,0,0/0,1,0,1/1,1,0,0/0,1,0,1 | (2,4,0,2) | 96 | 86 | Kron(e of 1,1,0,0/0,1,1,0/1,1,0,0/0,1,0,1, F2.1) (that type peeled) | 3 s |
| T6 1,0,1,0/1,0,0,1/0,1,1,0/0,1,1,0 | (2,2,3,1) | 576 | 456 | Kron(e of 1,0,1,0/1,0,0,1/0,1,1,0/0,0,2,0, E1.1) | 3 s |
| T7 1,0,1,0/0,0,1,1/1,1,0,0/0,1,1,0 | (2,2,3,1) | 1 152 | 832 | Kron(e of 1,0,1,0/0,0,1,1/1,1,0,0/0,0,2,0 (two extras, Exact 137 k), E1.1) | 86 s |

End^0 references (Exact): (3,1,3) 2,0,1/0,1,0/1,0,2 (48 flows, 13 k unknowns): dim 2, e_p in span — PASS; T5
(196 k unknowns, 73 s): dim 2, e_p in span — PASS. The others are 5-63 M unknowns (not computed).

Closure generators through e_p: T1 31 + 31 (2.4 s from layer A, 26 s fresh), T2 6 + 6 (Exact fallback, 390-540 s,
up to 0.5 M unknowns per degree), T3: **not certified** (Search closure exhausted at degree 4: 0 of 1 generators;
the Leaves closure route is leaves-first only without e_p).

## 3. The 43 NeedsIdempotent jobs (runs_c2/k4m/lib) [M]

Setup: `libbuild one --lam 2,2,2,2 --n 4 --field F --cache C --costs costs_2222_c25.txt --provider m4c-lf
--stamp-route -- JOB` per job (runs_c2/k4m/lib/runjobs.sh: resumable, 2 in parallel, ulimit 7 GB, 1 thread),
C = C2-5's cache (data-gl6 c2/l2222_n4_c25 restored), then the same in replay mode (`--replay`: provider refusing
every call, templates must be cached).

| prime | run | replay | job wall mean / median / max (run) |
|---|---|---|---|
| F31 | **31 / 43** | @@R31@@ | 61 / 43 / 539 s |
| F32b | **31 / 43** | @@R32@@ | 48 / 45 / 389 s |

Certified (both primes): pairs 109, 202, 403, 595, 618, 688, 972, 1042, 1175, 1190, 1209, 1371, 1390, 1404, 1417,
1418, 1458, 1528, 1629, 1789, 1791, 1837, 1882, 2052, 2082, 2086, 2275, 2298, 2368, closure T1, closure T2.
Routes: decomposable sources through the commutation reroute (109, 1789, 2052, 2086, 2275, 2298, ...), T3-source
pairs through the Leaves-peeled e_T3 (202), shortcut X = L_q (1175, 1629), decomposable targets through Search with
e_q (1190, 1371-1418, 1791, 1837, ...).

Not certified (identical for both primes):

| job | p ⊗ R | where it stops |
|---|---|---|
| 926, 949, 986, 997 (E1.1), 1435, 1472, 1483 (E1.2) | various ⊗ E1.x | target T4 (decomposable, via Kron(e_T5, F1.1)) not reachable: Search stuck, Exact 7.7-31 M; no reroute (E1.x does not commute past the last rung) |
| 1769 | 1,0,1,0/1,0,0,1/0,1,1,0/0,1,0,1 ⊗ E2.1 | reroute → dec_one(1,1,0,0/1,0,0,1/0,2,0,0/0,1,1,0, F1.2): target T6 in X = E0.1,E1.2,E2.2,E1.1,F1.1,F0.1,F2.1,F1.2 (576 flows = the flows of L_T6, X = B_T6 ⊕ B_{1,0,1,0/1,0,1,0/0,1,1,0/0,1,0,1}, not commutation-equivalent to L_T6); Exact 10.8 M |
| 1812 | T2 ⊗ E2.1 | reroute → dec_word(T5, E2.1·F1.2) → target T7 in X = E1.2,E0.1,E2.2,E1.1,F1.1,F0.1,F2.1,F1.2 (1 152 flows = L_T7's); Exact 63 M |
| 132, 1412 | T2 ⊗ E0.1, ... | Exact fallback killed at 7 GB (1.4-1.7 M unknowns; cap now 1.3 M) |
| closure T3 | | Search closure through e_T3 exhausted (degree 4) |

## 4. 4-strand knots: cover after M5i [M]

`examples/m5i_cover` (sindec::knot_jobs, upper bound = K_0 lower bound here) against the union of C2-5's F31 replay
and the M5i run (runs_c2/k4m/lib/status_f31.txt, 1 145 / 1 169 jobs ok):

| words | jobs | missing |
|---|---|---|
| class A (18 knots: 6_1, 7_2, 7_4, 7_6, 7_7, 8_4, 8_6, 8_8, 8_11, 8_13, 8_14, 9_10, 9_11, 9_13, 9_17, 9_18, 9_20, 9_23) | 450 + 93 each (the same set) | **6**: pairs 1235, 1444, 1470, 1894 (non-M5i, C2-5) + 1769, 1812 (M5i) |
| class B (8_15, 9_22, 11n_42) | 555 + 93 | 17 |
| Conway | 846 + 125 | 22 |
| T(4,2) = (σ1σ2σ3)², p6_1b | 508 + 93, 450 + 93 | 6 (as class A) |
| T(4,3), T(4,5), x2, x3 | 776-1 032 | 10-24 |
| 268 stabilized 3-strand knots (σ3^{±1} inserted anywhere into the C2-3 words) | | 162 without idempotent types: 0 missing; **every word that visits an idempotent type misses 1235, 1444, 1470** |

So no 4-strand knot (and no 4-strand word through an idempotent type) is computable at (2,2,2,2), N = 4 yet. The 6
class-A blockers end in the same kind of call: a small target (36-1 152 flows) inside an X of 576-7 776 flows whose
summand splits non-locally (Search: stuck branches; Leaves Lemma A: steps of 38-84 M unknowns; Exact 5-60 M).

## 5. Engine validation [M]

(a) **Forced decomposable ladders** (runs_c2/k4m/force; `sbox` at (2,2,2), prune 4, box N = 4, closures N' = 2, 3,
4, u and r; provider on for the new keys; references runs_c2/k3/out = C2-3):

| word | forced set | idempotent types (reduced rank m) | F31 | F32b |
|---|---|---|---|---|
| 4_1 (1 −2 1 −2) | A: 1,1,0/1,0,1/0,1,1 = E0.1,E1.1,F1.1,F0.1; 1,0,1/0,2,0/1,0,1 = F0.1,F1.1,E1.1,E0.1; 0,1,1/1,1,0/1,0,1 = F0.2,F1.1,E1.1,E0.2 (extras in degrees ±1) | 2 visited (m = 18, 48) | identical (6/6) | identical |
| 4_1 | B: 0,0,2/0,2,0/2,0,0 = F0.2,F1.2,E1.2,E0.2 (3 extras, degrees 0, ±1); 1,0,1/0,1,1/1,1,0 = F0.1,F1.2,E1.2,E0.1 | | identical | identical |
| 5_2 (1 1 1 2 −1 2) | A | | identical | identical |
| 5_2 | B | | identical | identical |
| 3_1·σ2 (1 1 1 2) | B | | = C2-2's trefoil (rank 22 / reduced 9) | same |

All runs: χ = decat_chi (6/6), D14 scaling, the frozen-trace replay at every box point.

(b) **The 11 C2-2 4-strand words** (Markov words σ1^m σ2^{±1} σ3^{±1}, σ1σ2σ1σ3, 4_1·σ3, 5_2·σ3) recomputed with
the C2-4 binary, library = C2-5 cache + M5i additions, `--replay --provider m4c-lf` (runs_c2/k4m/markov): F31
**11 / 11 identical** to C2-2 closure by closure; F32b @@VB32@@. (These words visit no idempotent type: this checks
that the engine changes leave the e = None path unchanged.)

(c) tests/m5i.rs (default: (3,1,3) e_p vs End^0, dec_rung with the source idempotent, closure through e_p, reduced
coordinates QP = id and PQ = e_p on every coset from the frozen data; (2,2,2,2) T5 by prefix propagation; ignored:
T5 vs End^0, 75 s).

## 6. Open problems

1. **Non-local splittings at large X** block the last 6 class-A jobs (4 non-M5i, 2 M5i) and 10 more M5i jobs (T4
   targets, closure T3). Concrete next steps: (i) for 1769 / 1812 the failing X has exactly the flows of L_T6 / L_T7
   and contains B_T6 / B_T7 plus small summands that Search may reach: choosing that word as the ladder of T6 / T7
   (a different minimal-flow tie-break, e.g. via `KRHOM_FORCE_LADDER` first) turns both into the X = L_q shortcut,
   with e from peeling the small summands; (ii) a Leaves pair route for decomposable targets (the m4c-lf provider
   did not try LeafPairs with e_q) and for the 36-flow targets of 1235 / 1444 / 1470 (Lemma A steps 38-84 M: a
   smaller-window plan); (iii) a reroute that also moves R past non-commuting rungs (sl2 commutation E_iF_i ≅ F_iE_i
   ⊕ [k] id) would cover the T4-target jobs.
2. 132 / 1412 need an Exact fallback of 1.4-1.7 M unknowns (> 7 GB here): raise ulimit or find a reroute.
3. The Exact fallback admits calls outside the cost table (raised caps for M5i calls); it is bounded by
   `KRHOM_M5I_EXACT_UNK` and by memory. A cost-table record per shape would be the principled admission.
4. No 4-strand engine run through an idempotent type exists yet at (2,2,2,2); the reduced-coordinate path is
   validated at (2,2,2) with forced ladders (5a). The first (2,2,2,2) one becomes possible with pairs 1235, 1444,
   1470 (all stabilized words through idempotent types need exactly those).
