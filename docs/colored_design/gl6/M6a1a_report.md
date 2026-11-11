# M6a-1a report: sengine core without a library (branch colored-m6a)

Date 2026-11-11. Machine: cloud container, 4 cores, 15 GB. Base: origin/colored-m5 bf326b7 merged with
origin/dev (ce91f0e) at the start (merge 1c5532b). No existing file was edited except `src/lib.rs`
(`pub mod sengine;`). New files: `src/sengine.rs`, `tests/sengine.rs`, this report.

Spec followed: SPEC_3C_GL6 rev1 §2 "M6a-1a" and brief; SPEC_M4_M6a §6.3-6.4 (state, pivots,
Markowitz, streaming rule, T5); D14 (frozen structure); R4 §7-9 (scope at 3 strands).

## 1. API (src/sengine.rs)

**Library interface.** `trait SLib<F>: Send + Sync` is the read-only part of SPEC_M4_M6a §5.1 as built in
M5a: `ar, lam, prune_n, ty, dec_word, rickard_comp_at, closure_gens, reduced`. It is implemented for
`sindec::DecData<F>` (pure delegation), so M6a-1b passes `&*Arc<DecData<F>>` as `&dyn SLib<F>`.
Helpers for the wiring: `rep_coset(&schur, lam, p)` (call it on the builder before `freeze`: Tids
are per Schur instance), `lib_shape(lib, p, rep)` (flows of the canonical ladder per coset; Err for
e_p ≠ None, which is M5i).

**Shapes and layouts.**
- `TypeShape { p, cosets, dims, rep }`: support cosets (ascending), block sizes m_c (reduced
  coordinates), rep(p) with m_rep = 1 (checked). `rank() = Σ m_c`, `s() = Σ m_c²`.
- `ShapeReg`: registered shapes and the per-scan layout cache keyed by (src type, tgt type)
  (a type determines its top block). `layout_of(src, tgt)` = joint support, rows = tgt blocks.
- Components are `smap::LocMap<F>` over these `Arc<Layout>`s (one contiguous vector per map, F10).

**Complex.** `SObj { p, h, t }`; `SCx<F>` with `objs`, `out: BTreeMap<x, BTreeMap<y, LocMap>>`,
`inn`, an incremental byte counter, and `k0`, `object_multiset`, `pattern_hash`, `stored_entries`,
`zero_block_frac`, `sum_s`, `d2_nonzero`. `compose_acc(dst, g, f, α)` = dst += α g∘f blockwise per coset
via `gemm_acc` (no allocation per block). `scalar_of(m, shape)` = Some(c) iff every block is c·I.

**Engine.** `PEngine<F>::new(mode, EngOpts)`; per crossing:
```
let prev = eng.begin_crossing()?;          // previous complex, read-only (Arc), counted in memory
for each term in the order the term differential runs:
    eng.begin_step()?;
    eng.add_obj(SObj)?            // objects of C ⊗ W_k (summands of dec_word)
    eng.add_map(x, y, LocMap)?    // d_C ⊗ id inside the term; id ⊗ d_W from the previous term
    eng.end_step()?;              // checkpoint, memory check, streaming elimination
eng.end_crossing()?;              // eliminate all, optional d∘d, determinism hash, StageStat
...
let (point_scan, final_cx) = eng.finish(point_id)?;
```
`load(cx)` sets a complex without elimination (the start object). `set_extern_bytes(b)` declares the
Eval memo / evaluated library blocks for the memory budget. `EngOpts { streaming, mem_point_bytes,
order: Heap | Rescan, check_d2, check_dead, ignore_eligibility (test only), id_coset }`.

**Outputs.** `StageStat { inserted, objects, components, pivots, rounds, stored_entries, sum_s, c_e,
zero_block_frac, bytes, peak_bytes, secs, hash }` per crossing; `MemStats { peak_bytes, peak_at,
peak_cx_bytes, peak_entries }`; `PointScan { point_id, objs, id_blocks (final components at
id_coset), trace, frozen: Arc<FrozenTrace>, stages, mem, k0, secs }` with `same_structure`;
`GramAtPoint { p, point_id, id_coset, ins, outs, g }` from `gram_at_point(lib, ev, p)`
(G[b][a] = h_b ∘ g_a at the identity coset, through `reduced`); `par_map` over KRHOM_THREADS.

**Errors.** `SErr::{Memory, PointDependent, Pivot, Stream, Shape, D2, Lib, Usage}`; `Display` gives
e.g. "per-point memory 1.8 MB above the budget 1.8 MB at stage 0 insertion: run with streaming or on
the cloud" and "point-dependent scan at stage s round r: ...".

**Synthetic complexes** (`sengine::synth`, used by tests and measurements): see §3.

## 2. Invariants (all enforced at run time, not debug-only)

1. **Pivots** (SPEC §6.4): a → b with type(a) = type(b), t(a) = t(b), h(b) = h(a)+1; the component must be
   c·I at every coset of the type (cost Σ m_c²) with c ≠ 0, else `SErr::Pivot`. Update
   d_xy −= d_ay c⁻¹ d_xb blockwise per coset, zero sums dropped.
2. **Markowitz order**: cost (|in(b)|−1)(|out(a)|−1), ties by (a, b). Lazy min-heap with per-object
   versions; after a pivot every object whose in- or out-set changed (in(a) ∪ in(b) ∪ out(a) ∪ out(b))
   gets a new version and its candidate pairs are re-pushed. Tested identical to the rescanning order.
3. **Streaming eligibility**: steps are inserted so that every term-differential arrow runs from step
   s to s+1 (ascending terms for a positive crossing, descending for a negative one: then in both
   cases the source term is inserted first). in(b) is always complete; out(a) is complete iff step(a)
   < current step. So a pivot is eligible iff step(a) < current step, all pivots after the last step.
   Run-time form: an arrow at an eliminated object is discarded only if it is an out-arrow of an
   object that died as a pivot *target* or an in-arrow of one that died as a pivot *source*
   (Gaussian elimination never reads these); any other arrow at a dead object is `SErr::Stream`.
4. **Frozen structure (D14)**: `Mode::Record` (point 0) decides pivots from values and records a
   `FrozenTrace { checkpoints, rounds: [Round { stage, pivots, pattern }], stages }`: the pattern
   hash (objects + nonzero components) after every insertion step, every elimination round's pivot
   list and pattern hash after it, and a determinism hash per crossing (sorted objects, pivot
   sequence, pattern). `Mode::Replay(Arc<FrozenTrace>)` never decides structure: it applies the
   recorded pivots (each must still be an eligible component and a nonzero scalar), and compares
   every checkpoint, round pattern and stage hash; `finish` asserts the trace was consumed exactly.
   Any difference is `SErr::PointDependent`. Soundness of the check: if every recorded pivot is a
   nonzero scalar at the new point, the result is homotopy equivalent to the input there; equal
   patterns imply it is minimal (no degree-0 component left); so a point where the minimal model
   changes cannot pass silently (it would need fewer objects than its minimal model).
5. **Memory**: `SCx` keeps an incremental byte counter (LocMap capacity × size_of F + 160 B per object +
   96 B per edge; `recount_bytes` checks it in the tests). Counted: current complex + previous complex
   (during a crossing) + shapes and layouts + declared external bytes. Checked after every
   `add_map`, every pivot with updates, every step and crossing; above `mem_point_bytes`:
   `SErr::Memory` (the refusal fires before more than one LocMap past the budget is held).
6. **K_0 / degrees**: `k0()` and `object_multiset()` are provided; the K_0 = decategorified-scan and
   LB asserts need schur's decategorified scan and are wired in M6a-1b. `check_d2` runs d∘d = 0 after
   every crossing (meant for unpruned scans only, R4 §8).

## 3. Synthetic model and tests (tests/sengine.rs; 11 default tests, 2.2 s release)

**Model.** Random type shapes (1-6 cosets out of ≤ 6, m_c ≤ 3, m_rep = 1). A component x → y has
degree e = t_x − t_y ≥ 0, and degree 0 only between equal types as c·I (Soergel: Hom^0(B_p, B_q) =
δ k, no negative degrees). A complex is P(D_min ⊕ cones)P⁻¹: D_min has arrows of degree > 0 only from
even h to h+1 (so d² = 0, minimal), cones c(x)·I, P = I + N per h with N of degree ≥ 0 strictly
increasing in (−t, index). Values are α + β·x with (α, β) from the structure seed, so the structure is
identical at generic points and degenerates exactly at the cone roots −α/β. The tensor C ⊗ W uses
product shapes (cosets (c1, c2), m multiplies), f ⊗ id and (−1)^h id ⊗ g as Kronecker blocks;
optionally every term is twisted by its own unitriangular automorphism Q_k (Q d Q⁻¹ blockwise over
terms: an isomorphic complex whose within-term pivots differ between terms, like the summand bases
of real decompositions). Known minimal model of C ⊗ W: C_min ⊗ W_min. Iso invariants compared:
object multiset, K_0, and the filtered rank profile (per coset, h and threshold s: rank of d_h on the
subcomplex t ≤ s, preserved by d and by every isomorphism of the model).

| test | what | result |
|---|---|---|
| known_homology_single_complexes | 150 complexes (988 objects, 47 contractible): eliminate to D_min (multiset, K_0, rank profile), d² = 0, no candidate left, byte counter exact | pass, 365 pivots |
| markowitz_heap_equals_rescan_order | 200 seeds × {full, streaming}: heap and rescan give identical round-by-round pivot sequences on C ⊗ W | pass, 4 780 pivots |
| streaming_equals_full_and_known_model | 200 seeds (160 with nonempty homology; 3/4 twisted): streaming = full = C_min ⊗ W_min; d² = 0; streaming peak ≤ full | pass; Σ peak 0.89× |
| streaming_sensitivity_ignored_eligibility_is_caught | eligibility ignored: refused (`Stream`) on 136/200 seeds; with the run-time check off, wrong result on 9/200 | pass (both ≥ 1) |
| identical_traces_at_two_points_and_replay | 60 seeds: independent records at 2 points give identical pivot traces; replay of point 0 at point 1 identical; 6 parallel replays (`par_map`) identical | pass |
| replay_refuses_a_degenerate_point | 210 degenerate points (cone roots; single complexes and streaming C ⊗ W): every replay refused with `PointDependent`; a record there has more objects | pass |
| replay_refuses_a_foreign_trace | trace of another complex, truncated trace: refused | pass |
| memory_accounting_and_refusal | streaming peak < full peak; a budget between them refuses full with the spec message and admits streaming; a tiny budget refuses both | pass |
| pivot_assertion_and_api_misuse | non-scalar degree-0 component → `Pivot`; h not rising, unknown type, wrong layout, double registration, steps misuse → Err | pass |
| library_trait_and_gram_at_point_on_11 | real (1,1) DecData through `&dyn SLib`: `lib_shape`, `gram_at_point` at 2 points (nonsingular; point-dependent when degrees > 0) | pass |
| measured_costs_small / _large (ignored) | cost table (§4) | pass |

Finding (recorded because it shaped the generator): with an untwisted product tensor, ignoring the
eligibility rule never changed the result on 200 seeds. Every within-term pivot then comes from a
cone of C, which is mirrored in the next term; the corrections that early elimination misses all land
on in(a′) of the mirrored pivot a′ → b′, which elimination discards. Real crossings do not have this
property (summand bases differ per term), hence the per-term twist; with it the rule is necessary
(9 wrong results without the check). The run-time check is stricter than necessary on the untwisted
model (it refuses harmless cases) — acceptable: it never refuses a correct streaming order.

Full suite `KRHOM_THREADS=4 cargo test --release`: see §6.

## 4. Measured costs on synthetic complexes (release, one thread)

Seconds are insertion + elimination of one crossing (the synthetic construction of C ⊗ W excluded);
"replay" is the same crossing replayed at a second point from the frozen trace (streaming).

| C objects | W objects | twisted | inserted | final | pivots | entries after | Σ S final | c_e final | zero-block frac | peak bytes full | peak bytes streaming | s full | s streaming | s replay |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 18 | 9 | yes | 162 | 18 | 72 | 108 | 770 | 0.14 | 0.00 | 293 684 | 225 148 | 0.001 | 0.001 | 0.001 |
| 34 | 14 | yes | 476 | 40 | 218 | 675 | 1 644 | 0.41 | 0.00 | 2 993 124 | 1 552 028 | 0.014 | 0.015 | 0.007 |
| 50 | 22 | yes | 1 100 | 60 | 520 | 1 666 | 3 562 | 0.47 | 0.02 | 20 791 288 | 13 879 832 | 0.367 | 0.308 | 0.104 |
| 100 | 32 | no | 3 200 | 160 | 1 520 | 2 634 | 8 416 | 0.31 | 0.00 | 9 455 416 | 5 189 192 | 0.062 | 0.055 | 0.026 |
| 150 | 42 | no | 6 300 | 300 | 3 000 | 17 686 | 20 124 | 0.88 | 0.17 | 34 087 216 | 21 634 344 | 0.288 | 0.310 | 0.134 |
| 200 | 52 | no | 10 400 | 480 | 4 960 | 43 560 | 35 607 | 1.22 | 0.03 | 73 667 264 | 47 159 336 | 0.663 | 0.914 | 0.291 |

(`measured_costs_small` = rows 1-2, default; `measured_costs_large` = rows 3-6, ignored, 44 s including
the synthetic construction. Twisting makes the terms dense, hence row 3's higher peak.)

Reading:
- Engine overhead is small: ≈ 60-90 µs per inserted object including elimination at 10 400 inserted
  objects; the per-entry block arithmetic dominates as blocks grow (real (3,3,3) blocks are far larger).
- **Replay is 2-3× cheaper than record** (0.29 s vs 0.66-0.91 s on the largest row): no heap, no
  candidate scans. D14's frozen structure therefore saves time at every box point after point 0.
- Streaming lowers the peak to 0.53-0.67× of full elimination on these products (≈ 1 − the fraction
  of the complex eliminable before the last term). On the largest row streaming costs 1.4× the time
  of full elimination (more rounds, each rebuilding the heap); at real sizes memory is binding, so
  streaming stays on for (3,3,3) (R4 §8).
- c_e on these synthetic complexes is 0.1-1.2, far below the 5-30 assumed in R4 §4: synthetic
  components are sparse in the block pattern. The real c_e comes from S3.

## 5. What M6a-1b must wire

1. **Shapes**: after building the DecData, register `lib_shape(lib, q, rep_coset(&builder.schur, lam, q))`
   for every type that occurs (summands of every `dec_word`, closure types). Keep a Tid → rep table
   next to the `Arc<DecData>` (DecData does not export tables).
2. **Crossing step** (SPEC §6.3 steps 1-3) into `begin_crossing` / `begin_step` / `add_obj` /
   `add_map` / `end_step` / `end_crossing`:
   - terms ascending for a positive crossing, descending for a negative one (one step per term);
   - objects (q_j, h + hdeg_k, t + dq_k + s_j) from `dec_word(p, W_k.rungs)`;
   - d_C ⊗ id: π_τ ∘ (f ⊗ id_{W_k}) ∘ ι_σ blockwise per coset into `LocMap::zeros(eng.layout(σ, τ))`
     (Kron rule of §3.4.2 on the f blocks; ι/π blocks from the per-thread `Eval`, through `reduced`);
   - id ⊗ d: `rickard_comp_at(p, i, ε, j)` evaluated by `Eval::all`, re-laid onto `eng.layout`,
     times (−1)^h, from place[(o, ks)] to place[(o, kt)];
   - `set_extern_bytes(ev.memo_bytes())` before each step.
3. **Asserts**: degree rule e(π) + e(f) + e(ι) = t_src − t_tgt from `MapInfo` degrees; K_0 = the
   decategorified scan (xengine's `decat_step`/`decat_k0` equivalents) and #objects ≥ LB per stage
   (use `SCx::k0`, `StageStat`); `check_d2` on unpruned scans only.
4. **Points**: record at point 0 with the builder's data, `freeze`, then replay at all other points with
   `Mode::Replay(point0.frozen)` via `par_map`, one `Eval` per thread; `PointScan::same_structure`.
   λz scans for the D14 scaling check are ordinary replays at the scaled point.
5. **Closure data per point**: `EngOpts.id_coset` = the identity coset (`identity_coset(ar, gens)`);
   `PointScan.id_blocks` gives f_id per final component; `gram_at_point` per final type.
6. Spec items left to M6a-1b by design: `SOpts` / `sscan` / `StageInfo` / `Measured` (the §6.2
   interface proper), the K_0/LB oracle, time budgets.

## 6. Open problems

- c_e, κ_S and seconds per unit of peak Σ S are unmeasured on real data (S3). The synthetic numbers
  above are an upper bound on overheads, not a model of real colored complexes.
- `id_blocks` stores one block per final component at one coset; the box driver (M6a-2a) may want
  only the identity-coset scalars of (λ, λ) types — trivial to narrow.
- The memory model counts BTreeMap overheads with fixed constants (160 B / object, 96 B / edge);
  validate against RSS on the first real scan.
- The heap is rebuilt from scratch at every round (once per term); fine at ≤ 4 terms per crossing.
- u32 LocMap storage (R4 S4) not done (conditional on S3).
