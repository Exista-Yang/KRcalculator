# R2: M5 (singular decompositions) for 3-colored gl_6 with Foam_N pruning built in

Date 2026-11-10. Inputs: SPEC_M4_M6a.md rev 2 (§4, §5; F1-F7), M4b_design.md, HANDOFF_M4b.md,
NTRUNC_STUDY.md, `src/homprov.rs` / `src/smap.rs` (worktree m4b, 59eaea2; origin/colored-m4b 9ca07ad
= that commit + merge of dev, homprov.rs unchanged), `src/xengine.rs` pruning code and `src/schur.rs`
(`cell_weight`, `killed_at`) on dev, G1(a) JSON (`colored_design/g1a.json`, `g1a_33_try.json`). New
counts: `gl6/py/blocks.py` (tables + Greene DP, 0.1 s). Nothing compiled or run, no repository edits.
Tags: [V] verified (here or in the cited run), [D] derived, [E] estimate, [?] open.

## 0. Recommendation in brief

1. **Pruning lives only in M5, never in the provider.** `HomProvider` stays N-agnostic: every Hom
   problem it solves (Hom spaces in SSBim, certified against the EKLP formula) is independent of N.
   M5 decides which (q, s) targets it asks for. Killed summands get no ι, π, no idempotent, no
   flow index and no further rungs.
2. **Two cache layers.** Layer A, expensive and N-independent: the certified raw pairs per target
   (X, q, s, mult), closure generators and End^0 sets. They are keyed by table strings and word,
   never by Tid. Layer B, cheap and per prune_n: selecting the kept targets, normalizing the form,
   Neumann orthogonalization, composites, Rickard components and checks. A DecData built at
   prune_n = 6 serves N' ≤ 6 through `restrict(N')`, with no provider calls.
3. **Certification.** Σ ιπ = e_X is no longer available. It is replaced by per-target form
   certification (unchanged), kept multiset = kept part of `term_mult`, Π∘Ι = id among the kept
   summands and rank bookkeeping, plus Krull-Schmidt and K_0. Exactness then rests on the kill
   rule, so the downstream checks are mandatory: d∘d components in I_N close to 0, Wang, duality
   at N = 4, pruned = unpruned where feasible.
4. **N = 6 prunes nothing on (3,3).** The 2-strand part of M5 is plain SPEC §5 with its full
   checks. Its G1 part is already done: 80/80 targets certified, 68 s cold [V]. Build this first;
   it unblocks M6a T(2,n) at gl_6 (Wang rank 136). It must not wait for the (2,2,2) G1(a)
   failures, which are 2-colored and off the target path.
5. **(3,3,3) at N = 6 is a fixed, knot-independent library.** It has 338 (type, rung) pairs with a
   kept source, of 494 (68 %). Build it once per prime on the cloud, in nested order N = 4 (42
   pairs), then 5 (158), then 6 (338), so that an interrupted session still leaves complete
   lower levels. Gate it with a new G1(c): every kept target at N = 6, not a 20-pair sample.
6. **Idempotents are probably not needed up to 3 strands at N ≤ 6.** (3,3) has none (F1). The
   only decomposable end-block (3,3,3) ladder, 0,1,2/1,1,1/2,1,0, has μ_1 = 7, so it is killed at
   N = 6 [V, Greene DP]. The other 2 of the 3 in the visited blocks still have to be classified
   (§2.4). Keep e = P(a) as a small code path, tested on forced ladders only.

## 1. What remains at N = 6 (and N = 4, 5 for the checks)

The blocks a scan visits are λ together with the intermediate blocks of one crossing term
W_k = [E_i(k), F_i(k)] (c = d = 3: the same words for both signs, `schur::crossing_terms`).
"Pairs" are the `g1_work` (block type, rung) pairs: rung E_i(k) on the λ-block types and rung
F_i(k) on the intermediate-block types. A source p is kept iff `!killed_at(p, N)` (Greene μ_1 ≤ N
and all margins ≤ N).

| λ | blocks | types: all / kept at N = 4, 5, 6 | pairs: all / with kept source at N = 4, 5, 6 |
|---|---|---|---|
| (3,3) | (3,3) 4, (4,2) 3, (5,1) 2, (6,0) 1 | 10 / 3, 6, 10 | 18 (6 identity) / 7, 12, 18 |
| (3,3,3) | (3,3,3) 55; (4,2,3), (3,4,2) 45 each; (5,1,3), (3,5,1) 27 each; (6,0,3), (3,6,0) 10 each | 219 / 12, 58, 138 | 494 (429 non-identity) / 42, 158, 338 |
| (3,3,3,3) | 10 blocks | 10 636 / 60, 824, 3 385 | 26 700 / 252, 2 640, 9 705 |

Kept types per (3,3,3) block at N = 6: 40/55, 30/45, 15/27, 4/10 (both mirror blocks equal) [V].
This matches NTRUNC's 40 of 55 for (λ, λ) and 790 of 2 008 for (3,3,3,3).

- **(3,3).** Unchanged at N = 6: 12 non-identity pairs, 72 targets, Σ mult 121, closure 4 types
  with 8 ins and 8 outs.
  - G1(a), cold, laptop: 80/80 certified in 68 s.
  - Tiers: S1 34, S2 38, S1+S2 4, S0 2, Exact 2. The slowest pair took 22.6 s; max X flows 90.
  - So the 5 401 s exact Hom of the N = 5 xengine run is gone once Search is used.
- **(3,3,3) at N = 6.** 240 composite crossing-term decompositions (40 kept types × 6 words;
  330 unpruned), built from:
  - 240 first-rung pairs;
  - 98 second-rung pairs (of 164) on the kept intermediate types.

  Further work items:
  - Rickard component lists: 40 × 2 positions × 2 signs × 3 = 480;
  - closure generators for 40 types (55 unpruned).

  Targets: 3 151 over the unpruned non-identity pairs [V-P11]. The kept-target count needs
  `term_mult` and is the first output of `pruned_plan` (§6), which uses schur only and takes
  seconds. [E] 55-70 % of them, because killed q concentrate in the larger intermediate blocks.
- **The library is independent of the knot.** A 3-strand scan only ever holds (λ, λ) objects
  between crossings. So the library above is complete for every (3,3,3) braid at N ≤ 6. The
  knot size enters only M6a: 4_1 LB 276, 8_19 LB 776, 12a_1248 LB 7.4e5 objects.
- **(3,3,3,3) at N = 6.** 9 705 kept-source pairs: the library approach becomes a lazy, per-knot
  cache (§4.3). The full canonical basis (≈ 800 s / 3.6 GB once, cloud) must be loaded from a
  `save_basis` file.

## 2. Pruned decompositions (changes to SPEC §5.2)

`SIndecOpts` gains `prune_n: Option<u8>` (None = unpruned), fixed for the life of the builder.
`killed(q) := prune_n.is_some_and(|n| schur.killed_at(q, n))`. Use only this predicate: no ladder
labels (they fail for (3,3,3,1)).

### 2.1 `ty(p)`
- p killed: `Err(DecErr::Killed { p, n })`; no ladder, no flow index. p kept: as SPEC.
- On ≤ 3 strands a kept p has a minimal-flow ladder with labels ≤ N. Assert this as a debug
  check (true for every 2- and 3-strand block tested; NTRUNC §2).

### 2.2 `dec_rung(p, R)` for kept p
1. X = L_p·R and mult = `term_mult(p, [R])`. Split mult into kept K and killed D.
   `killed_rank = Σ_{(q, m) ∈ D} m(1)·rank(q)`.
2. **All killed** (K empty). Return `Decomp { summands: [], all_killed: true }` after asserting
   killed_rank = flows(X). No provider call. M6a drops the object.
3. **Identity shortcut.** As SPEC: a single (q, 1) with q kept, L_q = X and e = None. If that
   single q is killed, case 2 applies.
4. **Targets.** Loop over K only, in SPEC order (l_max(q) descending, then s ascending).
   - Each (q, s, cnt) is a layer-A lookup (§4) or a provider call.
   - M5 re-verifies the returned family itself: `homprov::form_at` at rep(q) at 2 points, rank
     cnt, `assert_genuine`, and `ar.info(m).deg` equal to s and −s.
   - π' = G[J][I]^{-1} π_J as in SPEC.
5. **Orthogonalize among kept summands only.** N_ts = π_t∘ι_s for kept t ≠ s is unitriangular by
   the degree argument, which holds for every subfamily; then Π' = Σ (−N)^k Π. Checks: §3.

### 2.3 `dec_word(p, rs)` (rung by rung)
- A killed summand B_r of L_p·R_1 is never tensored with R_2: its summands are all killed, since
  the killed span is a two-sided ideal with positive structure constants.
- Assertion: the multiset of (q, total s) equals the kept part of `term_mult(p, rs)`. This is a
  real, non-trivial check of the ideal property at the K_0 level, with the same code as xengine.

### 2.4 Idempotents
- Compute e_q = P(a) only for kept q with `ladder_char(L_q) != b_q`. At start-up, M5 lists
  `idem_needed(λ, N)` over all visited blocks (schur only).
- Expected:
  - (3,3): empty (F1);
  - (3,3,3) end block: empty at N ≤ 6 (0,1,2/1,1,1/2,1,0 is killed, μ_1 = 7);
  - the other 2 of the 3 in the visited blocks: [?], to be classified by the M5 report.
- If a kept q needs e_q even though L_q's excess is entirely killed, the excess still has to be
  split off (End^0(L_q) ≠ k, so the form test needs it). So do not try "L_q ≅ B_q mod I_N".
- Apply Kron(e_p, R) on the X side and pass e_q to the provider, exactly as SPEC.

### 2.5 `rickard_comp(p, c, d, sign, j)`
- Valid only if p, the term W_j and the term W_{j±1} are kept. Otherwise `Err`. M6a builds the
  stupidly truncated Rickard complex.
- Assert, as `xengine` does, `killed_at(p_k, N) == (3 + k > N)` for the term types p_k (cell
  (3+k, 3−k)).
- Components only between kept σ and τ. Zero components are dropped as before (single-map test).

### 2.6 Closure generators
`provider.closure` only for kept types of (λ, λ). Layer A, N-independent.

## 3. Certification when summands are dropped

| statement | status in pruned mode | how |
|---|---|---|
| per-target pairs: form rank = mult at rep(q), same matrix at a 2nd point, genuine, homogeneous | **required, unchanged** | M4b contract, re-verified by M5 (§5) |
| kept multiset = kept part of term_mult (one rung and composite) | **required** | assertion (§2.2-2.3) |
| Π∘Ι = id among kept summands | **required** | 2 points, all cosets of the support of X; composites re-checked at 2 points |
| Σ ι_t π_t = e_X | **not available** | replaced by the argument below |
| Σ_kept m·rank(q) + Σ_killed m·rank(q) = flows(X) | bookkeeping only (consistency of schur, not of maps) | assertion; needs the killed b_q (full basis) |
| per-stratum rank of e_kept at rep(t) = Σ_kept m·c_t(b_q)(1) | implied by Π∘Ι = id | debug only |
| all-killed X: killed_rank = flows(X) | bookkeeping | assertion |
| Rickard term pruning = cell rule | **required** | assertion per crossing |

**Why the kept family suffices** [D, as in xengine]:
1. Π∘Ι = id makes e_kept = Σ ιπ an idempotent with Im e_kept ≅ ⊕_kept B_q{s}.
2. Its complement has class [X] − Σ_kept = Σ_killed m_q b_q in K_0.
3. Classes of indecomposables form a basis (Soergel-Williamson, assumption A2 of NTRUNC), so
   with Krull-Schmidt the complement is ≅ ⊕_killed B_q{s}.
4. By the kill rule the complement is 0 in Foam_N. So (ι, π) are mutually inverse modulo I_N,
   and the transported Rickard components are correct modulo I_N.

The only non-computational input is the kill rule. It is validated decategorified on 9 833
(type, N) pairs plus the (3,3,3) N = 9 row (cloud run F3, 0 mismatches), and categorified on 2 strands.

**Downstream checks (M6a, required because of the dropped summands):**
1. K_0 after every crossing = the untruncated decategorified class with the killed types removed
   (same as xengine; equal to NTRUNC's truncated scan).
2. d∘d ≠ 0 on the nose. The final components (∈ I_N) must close to 0 at every requested
   N' ≤ prune_n (`closure_kills_map`), with check_d2 forced on, as `xengine::scan` refuses
   otherwise.
3. Pruned = unpruned homology where both are feasible: (3,3) at N = 3, 4 against xengine pruned
   and unpruned and Wang.
4. 3-strand (3,3,3) at N = 4 against uncolored gl_4 of the reversed braid (duality N = a+1,
   "expected").
5. Self-consistency: the homology from `restrict(4)` of the N = 6 library equals that of a fresh
   prune_n = 4 build.
6. Wang k = 3, N = 6 for T(2,m).

## 4. Caches, laziness and a fixed prune_n

### 4.1 Layer A (N-independent, expensive, on disk)
- Unit: one certified provider result:
  - `pairs` for a target (X word, q's table, s, mult, e flags);
  - `closure` (p, outs);
  - `end0` (L).
- Format: the Arena sub-DAG, templates by content hash, as SPEC §5.2 "Disk cache".
- Key: λ, field, table strings and rung words. Never Tids: they depend on interning order.
- Header stamps: `krhom-sindec v2`, `TMPL_ALG` (now "hs1.lcx1"), the field, and a provider stamp.
  A bump in any of them invalidates the cache.
- Loading re-runs the per-target checks at 1 point.
- Reused by every prune level and by the G1(c) harness, which should save its certified maps
  rather than only statistics.
- Cloud sessions last 5 h and lose local disk. Push layer A and the template cache as an archive
  on the work branch, or regenerate them. Report the size; [E] tens of MB per λ and prime.

### 4.2 Layer B (per prune_n, cheap, in memory)
- Contents: target selection, π' normalization, Neumann series, composites, Rickard components and
  the §3 checks.
- Rebuilt for each prune_n from layer A. It contains no provider calls once A is warm.

### 4.3 Laziness and memory
- SIndec is lazy, as SPEC: `dec_word` is called on demand by the first M6a scan. For 2 and
  3 strands, add `SIndec::prebuild(&plan)`, which visits every kept pair of `pruned_plan` in
  nested order (§6), so the library is complete before scans start.
- 4 strands stay lazy. Never intern a ladder, flow index or Arena node for a killed type. Killed
  Tids appear only inside `term_mult` output; their rank is needed for bookkeeping only.
- With NTRUNC's restricted 4-strand basis (X_q = char of the minimal-flow ladder for killed q),
  the rank assertion becomes vacuous. Keep the kept-multiset assertion.

### 4.4 Using one library at several N
- Kill sets are nested: killed at 6 implies killed at N' < 6, since I_6 ⊆ I_N'. So a family
  orthogonalized at level 6 and filtered to the summands kept at N' still has Π∘Ι = id.
- `DecData::restrict(n)` (n ≤ prune_n) therefore returns filtered Decomps and Rickard lists over
  the same Arc'd Arena. A request with n > prune_n is an Err. `DecData::prune_n()` is exposed,
  and M6a asserts its own N against it (a pruned scan is valid for N' ≤ prune_n only).
- The in-memory key (p, rung word) needs no N because prune_n is fixed per instance. Disk keys
  carry no N because layer A is N-free.

## 5. Interface with M4b

### 5.1 Points that may still change
The cloud session is fixing the stuck (2,2,2) branches and finishing G1(a). Risk is given for
each point.

| item | risk | why |
|---|---|---|
| `SearchStats` fields, `Search` tuning fields (`max_window_unknowns`, `win_w`, `win_len`, `costs`), `Search::tiers` | high | the stuck-branch fixes will add counters and knobs |
| `Budget` fields and defaults | medium | (2,2,2) pairs hit the 600 s and candidate caps |
| `MissKind` variants (e.g. a separate `Stuck`) | medium | the "2 stuck branches" are reported as `Exhausted` today |
| `Chain::for_labels` order, `Exact` admission, `CostTable` | medium | laptop-calibrated; the cloud is ≈ 1.5x slower and needs `g1 calib` |
| `HomProvider` method signatures, `Pairs` fields | low-medium | might gain a request struct or a trace or seed argument |
| `with_idem` (where e_q is applied) | low | M5 relies on it only if a kept e_q exists (§2.4) |
| `Arena::commit` remapping, `assert_genuine` leaf kinds, new `MapNode` variants (e.g. commutation isomorphisms as a node kind) | low-medium | breaks the layer-A serializer if it lives in sindec |
| `g1_targets`, `G1Pair` | low | could gain filters |
| `form_at`, `certify_one`, `stacked_rank`, `detect_plan`, `hom_gens` | low | pure certification helpers |

### 5.2 Proposed stable interface
1. **One adapter in `sindec.rs`.** It is the only code that touches `HomProvider`:
   ```rust
   pub struct PairTarget { pub x: ExprId, pub q: Tid, pub lq: ExprId, pub eq: Option<MapId>, pub s: i32, pub mult: usize }
   pub struct ProvOut { pub maps_in: Vec<MapId>, pub maps_out: Vec<MapId>, pub secs: f64, pub summary: String }
   pub struct ProvFail { pub engineering: bool /* Time | Memory | Refused */, pub text: String }
   pub trait DecHom<F: Field> {
       fn pairs(&mut self, ar: &mut Arena<F>, s: &mut Schur, t: &PairTarget, b: &Budget) -> Result<ProvOut, ProvFail>;
       fn closure(&mut self, ar: &mut Arena<F>, s: &mut Schur, lp: ExprId, p: Tid, ep: Option<MapId>, outs: bool, b: &Budget) -> Result<ProvOut, ProvFail>;
       fn end0(&mut self, ar: &mut Arena<F>, s: &mut Schur, l: ExprId, b: &Budget) -> Result<ProvOut, ProvFail>;
   }
   impl<F: Field, P: HomProvider<F> + ?Sized> DecHom<F> for P { /* the only place that destructures Pairs, Missing, SearchStats */ }
   ```
   - `Missing.kind` is matched with a wildcard. Anything not `Exhausted` counts as engineering.
   - `SearchStats` reaches M5 only as `secs` and a `Debug` summary string.
2. **M5 re-certifies what it stores.** It calls `form_at`, `assert_genuine` and `ar.info` itself
   and does not trust `Pairs.form`. So changes to provider internals cannot weaken M5's
   guarantees, and a provider returning more than `mult` maps is fine, because M5 re-selects
   with rref.
3. **The provider stays N-agnostic.** No `prune_n` argument. Search's window tree already recurses
   only into children whose character contains the requested kept q, and a killed child contains
   only killed types. So no pruning hook is needed in homprov now. An optional
   `Search::prune_hint(n)` can come after M4b merges.
4. **One provider construction point,**
   `sindec::provider_for(lam, costs: Option<&Path>, slowdown: f64) -> Box<dyn HomProvider<F>>`.
   It applies the cloud factor to the CostTable admission. The invariant "one Search per Arena"
   (TmplIds are Arena-local) is enforced by SIndec owning both.
5. **Node (de)serialization goes into `smap`** (`Arena::export(&[MapId])` / `import`), next to
   `commit`. Then a new `MapNode` variant breaks one match in one file, not the cache format in
   sindec.
6. **Branching.** Start M5 from origin/colored-m4b once the cloud session's final commit lands; it
   has `killed_at` through 9ca07ad. Until then, write sindec against `DecHom` only. Rebasing then
   touches one blanket impl. Do not add `--prune` to `src/bin/g1.rs` before M4b merges: put
   `pruned_plan` in sindec, and let g1 call it afterwards.

## 6. M5 tests and gates (changes to SPEC §5.3-5.4 and §4.7)

New planning function, `pruned_plan(s, λ, n) -> Plan`. For each pair it lists the kept and killed
targets, all-killed and identity flags, and `min_level` = max(μ_1(p), μ_1(q)) per target. It is
the input of `prebuild` and of G1(c), and its counts go into the report (fills the [E] of §1).

- **D1.** Unchanged for (2,2), (3,3), (2,1,2). For the 7 (3,3,3) blocks, add:
  - killed_at(·, N) for N = 3..9 against NTRUNC's counts;
  - the debug check "kept p has a ladder with labels ≤ N";
  - `idem_needed(λ, N)` recorded (blocking only if non-empty and End^0 cannot be certified).
- **D2.** Keep the forced ladders BS(s,s) over (1,1) and [E,F,E,F] over (2,2). Mark (3,3,3)
  0,1,2/1,1,1/2,1,0 (killed at the target N) and (2,2,2,2) (2-colored) ignored. Add a D2' case
  for every kept type in `idem_needed(3,3,3; 6)`, if any.
- **D3.**
  - (3,3), all 12 non-identity pairs: unpruned with the full SPEC checks, including Σ ιπ = e_X.
  - The same pairs pruned at N = 3, 4, 5 with the §3 checks; there the pruned multisets equal
    the filtered unpruned ones.
  - (2,2,2): becomes ignored (off the 3-colored path; G1(a) pairs 18, 20, 25, 32 fail there).
  - (3,3,3) pruned at N = 4: default if it takes under 60 s, else cloud. N = 5, 6: cloud.
- **D4.** W_k words of (3,3) at N = 3..6 and of (3,3,3) blocks at N = 4 (default or cloud) and
  N = 6 (cloud): kept multiset, Π∘Ι = id, all-killed counts.
- **D5.** Against M6.0 on (2,2), (2,1): unchanged. New: against xengine with `XOpts::prune_n` on
  (3,3) at N = 3, 4 (about 0 s there [V]): equal summand lists (q, s), and transported components
  equal up to basis change (rank per coset).
- **D6, D8, D9.** Unchanged.
- **D7.** Layer-A round trip; `restrict(4)` of a prune-6 DecData gives the summand lists of a
  fresh prune-4 build; restrict(7) is an Err; a stale `TMPL_ALG` stamp is rejected.
- **D10 (new).** `ty` and `dec_rung` on a killed p, and `rickard_comp` with a killed term, give
  Err. All-killed decompositions make no provider call (call counter). The Rickard-term pruning
  assertion holds.
- **G1(a).** Split for the 3-colored target: G1(a)-(3,3) has passed (80/80). G1(a)-(2,2,2) stays
  a 2-colored gate and does not block M5, M6a or 3-colored work.
- **G1(b)** (3,3,3) sample is replaced by **G1(c)**:
  - Scope: every kept target of `pruned_plan(3,3,3; 6)`, run on the cloud in order of
    `min_level` (4, then 5, then 6), one process per pair, watchdog 600 s and 1 GB.
  - Passes iff 100 % of each level's targets are certified. Per-level pass/fail is reported, so
    level 4 can open the M6a duality check while level 6 is still running.
  - Failures: rerun as in §4.7. An `Exhausted` result at the larger budget goes to the M4c set.
  - Certified maps are saved as layer A.
- **Done criteria.** Add to SPEC §5.4: DecStats with `killed_summands`, `all_killed`,
  `skipped_targets` and the provider seconds saved (estimated from G1); layer-A size;
  `idem_needed`; the restrict-consistency result.

## 7. What SPEC §5 can simplify or skip, given that only N = 6 (and N < 6 for checks) matters

1. **Idempotent machinery.** Implement it minimally (Krylov and Horner as specified). Test it
   only on forced ladders and on D2' cases. The (3,3,3) killed case and (2,2,2,2) leave the
   M5 done criteria.
2. **Direct (whole-word) mode and half-step adjunction.** Skip; rung by rung only (F3). S5 is
   already unimplemented.
3. **Padded (n+1)-upright model.** Skip.
4. **Unpruned (3,3,3).** Not built: 156 pairs and 81 killed types (15 in the end block) cannot
   matter for N ≤ 6. An unpruned build is needed only for (3,3), where it coincides with N = 6.
5. **`reduced()`.** It is the identity for all types when `idem_needed` is empty, so it can be a
   stub with an assertion.
6. **Closure side (M6a, recorded here).** The Δ-tower path is not needed at N = 6 for a = 3: the
   exact oracle in the Gr(3,6) box has 777 unknowns per entry (48 reduced) on 3 strands. The
   4-strand weight-to-degree question (weights up to 40, box degree 36) does not affect M5.
7. **One library per prime.** F31 and F32b each get their own layer A, since templates are
   computed per prime. The structure (targets, which tiers hit) is prime-independent. The
   F32b build can replay the F31 tier choices to save Search time ([E], optional).

## 8. Order of work and estimates

1. **M5-core on (3,3)**, unpruned = N = 6, full checks (D1, D3, D4, D5, D7-D9). [E] about 2 min
   per prime cold (G1 68 s plus orthogonalization and Rickard components). Then M6a 2-strand: T(2,m)
   at gl_6 against Wang (trefoil 136).
2. **Pruned mode plus layer A/B, restrict, pruned_plan, D10**, tested on (3,3) at N = 3, 4, 5
   against xengine (D5).
3. **G1(c) and prebuild of (3,3,3).**
   - N = 4 (42 pairs): [E] under 30 min. It enables 3-strand duality against uncolored gl_4.
   - N = 5 (158 pairs) and N = 6 (338 pairs): [E] 2-10 h per prime at 20-100 s per pair, i.e.
     2-4 cloud sessions per prime. The Search cost of (3,3,3) pairs is unmeasured (G1(b) never
     ran), which is the main schedule risk.
4. **4 strands.** Lazy per knot (§4.3), after the (3,3,3,3) basis file exists and its
   `idem_needed` is known.

## 9. Open points

- [?] The other 2 decomposable canonical ladders over the visited (3,3,3) blocks: are they killed
  at N = 6? `idem_needed` decides it in seconds.
- [?] The Search completeness rate on (3,3,3) pairs. The (2,2,2) stuck branches (pairs 18, 25,
  32) point at leaf-isomorphism search as the weak tier. G1(c) level 4 is the cheap probe.
- [?] The kept-target count and Σ mult for (3,3,3) at N = 4, 5, 6 (`pruned_plan`).
- [?] The categorified kill rule beyond 2 strands rests on assumption A1. The downstream checks
  of §3 (d∘d closing to 0, duality at N = 4) are its first 3-strand tests.
