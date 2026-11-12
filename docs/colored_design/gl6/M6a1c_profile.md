# M6a-1c profile: the per-point scan of the localized path (before optimization)

Date 2026-11-12. Machine: cloud container, 4 cores (Xeon @ 2.10 GHz), 15 GB. Code: colored-m6a at
c9d97b9 (M6a-1b + merge of colored-m5), with profiling timers only (no change of the algorithm).
Tool: `src/bin/sprof.rs` (new): builds the wired library for the word, scans the record point
(Markowitz heap, K_0 check) and K replay points of the N = 6 box, each single-threaded, and splits
the wall time of every scan into

| part | what |
|---|---|
| dec+obj | `dec_word` lookups and `add_obj` of the summands |
| dC⊗id | d_C ⊗ id per term: Kron index, π K ι per component, coset and summand pair, `add_map` of the results |
| rickard | (−1)^h id ⊗ d: `rickard_comp_at`, `lib_map` evaluation, `add_map` |
| extract | closure entries at the identity coset (Gram, pairing) |
| k0 | decategorified scan, K_0 / LB asserts (record point only) |
| pivots | `apply_pivot` (Gaussian elimination updates) |
| pattern hash | sparsity-pattern hashes (checkpoints after every step and round, D14) |
| rest | everything else: at the record point the Markowitz heap bookkeeping |

plus counters (`KRHOM_SPROF=1`): pivot updates and multiply-adds, and for d_C ⊗ id the
(component, coset) visits, the multiply-adds of K ι and π (K ι), the summand pairs visited.
λ = (3,3), word σ^M, F31, node seed 101, replay points 133 and 266 of the 400-point box (μ-grid of N = 6).
"Library evaluation" (ι, π, Rickard and Dot/Comp blocks at a point) is memoized per point by
`Eval` (11.8 MB at every M): it is not a cost.

## 1. Time per point (seconds, one thread)

| M | point | total | dec+obj | dC⊗id | rickard | extract | pivots | pattern hash | rest |
|---|---|---|---|---|---|---|---|---|---|
| 4 | record | 4.45 | 0.002 | 2.02 | 0.06 | 0.001 | 0.45 | 0.02 | **1.90** (heap) |
| 4 | replay | 2.63 | 0.002 | **2.06** (78 %) | 0.06 | 0.001 | 0.47 | 0.02 | 0.02 |
| 5 | record | 24.4 | 0.005 | 7.99 | 0.12 | 0.002 | 2.52 | 0.07 | **13.7** (heap) |
| 5 | replay | 10.5 | 0.004 | **7.8-8.1** (75 %) | 0.11 | 0.002 | 2.36 | 0.07 | 0.03 |
| 6 | record | 104.6 | 0.012 | 25.9 | 0.27 | 0.003 | 10.3 | 0.23 | **67.8** (heap) |
| 6 | replay | 32.3-35.0 | 0.009 | **23.2-25.4** (72 %) | 0.22 | 0.004 | 8.6-9.0 | 0.21-0.23 | 0.05 |

Growth per crossing (replay): ×4.0 (4→5), ×3.2 (5→6). Pivot updates 3.0e5 / 1.5e6 / 6.0e6 at M = 4/5/6
(1.3e8 multiply-adds at M = 6: 0.7 µs per update for ≈ 21 multiply-adds, i.e. overhead-bound).

## 2. Where the time goes

**d_C ⊗ id (≈ 75 % of every replay point).** Counters at M = 4 / 6:

| M | (component, coset) visits | K ι multiply-adds | π (K ι) multiply-adds (bound) | summand pairs visited |
|---|---|---|---|---|
| 4 | 14 270 | 6.7e6 | 1.6e8 | 6.1e6 |
| 6 | 143 193 | 6.2e7 | 1.6e9 | 5.9e7 |

The per-pair construction visits ≈ 400 (ι_s, π_t) pairs per (component, coset) (≈ 20 summands on
each side), and per pair does one `Eval::block` hash lookup for π, one allocation for the
product, a scan of the nonzero rows of K ι_s, and a BTreeMap insert. The arithmetic (dense
π · K ι over the cosets of X) is 1.6e9 multiply-adds at M = 6; the overhead per pair is of the
same order. ι and π are looked up again for every component, although they depend only on
(type, term, coset).

**Elimination (pivots, ≈ 22-27 %).** 0.7-1.5 µs per update: per update a registry layout lookup,
`compose_acc` with two binary searches per coset, a BTreeMap insert of a freshly allocated zero
LocMap for every fill-in edge, an `is_zero` scan; the multiply-adds are ≈ 21 per update.

**Record point: the Markowitz heap (56-65 % of the record point).** At M = 5 the round of the
last step pushes 5.6 M pairs (2.9 M out-side + 2.6 M in-side) and pops 5.6 M; the heap reaches
5.4 M entries (130 MB, uncounted by the budget); every push and pop does 4-6 BTreeMap/HashMap
lookups (alive, version, degrees, eligibility). Cost: 10.2 s of the 13.3 s of that crossing.
This is a one-off per run, but it grows like the replay (≈ ×4-5 per crossing: M6a-1b measured
232 s at M = 7, > 20 min at M = 8).

**Not costs.** dec_word / add_obj (< 0.01 s), Rickard maps (< 1 %), closure extraction (< 0.01 s),
pattern hashes (< 1 %), the K_0 check (< 0.01 s), library evaluation (memoized per point).

## 3. Memory per point

Engine peak (the per-point budget counter; = VmHWM per thread up to allocator slack): 48 MB (M = 4),
130 MB (M = 5), 322 MB (M = 6), 756 MB (M = 7, M6a-1b), > 2 GB (M = 8, record point refused).
At M = 5 the peak is during the insertion of the last step of crossing 5: **4 332 objects and
404 370 components** (93 per object); complex 130 MB = 292 bytes per component, of which ≈ 164
bytes LocMap data (≈ 20 u64 entries) and ≈ 128 bytes fixed overhead (BTreeMap entries of `out`
and `inn`, LocMap header, allocation); the previous complex is 0.2 MB. The dense summand-pair
structure of d_C ⊗ id (every ι_s of X_a against every π_t of X_b) makes the unsimplified step
dense; the eliminations of a step can only start at the next step (streaming eligibility).

Levers by size: u32 entries (both primes < 2^32) save ≈ 28 % of the component bytes; the fixed
overhead per component is ≈ 44 %.

## 4. Plan (result-preserving only)

1. d_C ⊗ id: assemble ι = [ι_1 | ι_2 | …] and π = [π_1; π_2; …] once per (type, term, coset) at
   the point (cache per step), compute π (f ⊗ id) ι as one product per (component, coset), cut
   into summand pairs; exact arithmetic, so the blocks are bit-identical.
2. Elimination: per (x type, pivot type, y type) composition plans (block offsets), no binary
   searches, one BTreeMap entry per update.
3. Record point: the same heap order with dense per-object arrays, plus compaction of stale
   entries.
4. Then memory: u32 component storage; per-point release.

Verification of every change: value digests of all closure entries at replay points, the record
trace digest (objects, pivots, pattern hashes), the sengine/swire test suites, and the G2B
σ^1..σ^5 comparisons.
