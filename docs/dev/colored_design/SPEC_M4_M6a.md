# SPEC M4a / M4b / M5 / M6a: localized substrate, Hom provider, singular decompositions, R-full engine on <= 3 strands

Date: 2026-11-10. Author: spec agent. Repository state: `/Users/existair/Documents/KRcalculator`, branch `dev`, commit 7070814 (M3a `schur.rs`, CL1 `cclose.rs`, M6.0 `xengine.rs` merged).

Inputs, read in full: SYNTHESIS.md (all sections; binding: §1.3, §1.4, §2, §2.4, §3.1, §4, §5); R1 §2-5; D1; C1; C2 (H1-H3, M1-M6); C3 (H1, H3, E1-E5); R2 (§2.2, §2.4, §2.5, §3, §5); `m6_0_report.json`; the sources `xengine.rs`, `local.rs`, `sbim.rs`, `indec.rs`, `engine2.rs`, `cclose.rs`, `schur.rs` (API), `cube.rs` (`LocalCx`, `local_cx`, `apply_local_poly`), `crossing.rs`, `web.rs`.

New measurements in this note: probe crate `colored_design/proto_spec/` (12 binaries, path dependency on krhom, release build in `colored_design/target_spec`). Every number tagged [V-P*] comes from probe P* (Appendix M gives the command and raw output file). No repository file was edited.

**Revision 2 (2026-11-10), after two adversarial critiques** (theory/correctness and engineering/scale). Every blocker and major issue was verified, re-running the evidence where it was cheap, and all of them were adopted (one with a modification). §9 (Critique log) lists each issue with its decision, reason and the sections changed. Probes of this revision: crate `colored_design/rev_probe` (`r1_e2all`, `r2_t0`), plus re-runs of the critique probes `e2const`, `witness`, `c2_store`, `c3_interp`, `c4_e2` and `c6_flowidx`. Outputs are in `colored_design/rev/`, and numbers from them are tagged [V-R1], [V-R2] or [V-crit]. The main changes:
- **Laptop envelope** (§0.2 F6, §6.6). Under M6a the laptop covers 2 strands, (3,3) up to T(2,5), uncolored 3 strands, the 3-strand Markov family and 4_1. 5_2, 6_2, 6_3 and 8_19 are cloud jobs. 2-colored 3-strand words with ≥ 9 crossings additionally need M6c, and stay cloud jobs even then. SYNTHESIS §1.3 is superseded on this point.
- **Closure policy** (§6.5, §6.6). On the laptop, reduced closures run for N ≤ 6. Unreduced closures run wherever the exact estimate is under the cap: N ≤ 4 for the (2,2,2) knots and N = 5 for 4_1. The Auto method reports both memory estimates when it refuses.
- **Binding storage layouts** (§3.3). Localized maps are one contiguous `Vec<F>` over a shared `Layout`, with dense row-major blocks. `FlowIndex` is a packed mask array. Search runs in an Arena scratch region (`mark`/`commit`) and builds flow indices lazily.
- **Thread safety** (§3.3.4, §5.1). `Arena`, `TmplLib` and `DecData` are `Send + Sync`, which is asserted at compile time. `Rc` is used only inside the per-thread `Eval`.
- **Closed forms** (§3.4.3, A7). The forms are renormalized so that every constant is 1. This is verified on all 60 E2/E3 cases with a + b ≤ 6 and both orientations [V-R1], and A7 covers all of them by default.
- **New oracles** (§6.7):
  - T0: single-type graded closure = MOY [V-R2];
  - T9: the 3-strand Markov family against xengine, now ignored-but-required;
  - T10: cabling as a cloud reference;
  - Appendix R.2 states that the N = 3 duality sees weights ≤ 6 only.
- **Certification and providers** (§4.3, §4.5). Detecting-set counting bound with fallback. Exact admission by label and a measured cost table; for labels ≥ 3, Search S0-S2 runs before Exact.
- **Interpolation and memory** (§6.5). F samples are stored per entry at the entry's own weight. Interpolation is batched per weight with 8 residual points. Odd and negative weights are asserted zero. Memory is asserted before rescans and per point.
- **Smaller fixes**:
  - the component degree identity (§6.3);
  - the streaming test, now term by term (§6.4);
  - the polynomiality certificate shares one inverse per degree (§3.4.5);
  - the G1(b) rerun rule (§4.7).

Tags as in SYNTHESIS: [L] literature, [V] verified numerically, [D] derived and not independently checked, [?] open, [E] estimate.

Section map:
- §0 summary, measured findings, decisions;
- §1 interfaces between the milestones and changes to existing modules;
- §2 conventions;
- §3 **M4a, the first task (ready to implement)**;
- §4 M4b and gate G1;
- §5 M5;
- §6 M6a;
- §7 where each SYNTHESIS §1.4 decision is implemented and tested;
- §8 risks;
- §9 critique log (revision 2);
- Appendix M: measurements; Appendix R: reference values.

---------------------------------------------------------------------------------------------------

## 0. Summary

### 0.1 What the four milestones deliver

| milestone | delivers | new modules | gate |
|---|---|---|---|
| M4a | localized substrate: cosets; ladder expressions with zero labels and flows at a coset; a map arena of genuine global maps, evaluated one coset at a time at stacked points; elementary maps E1-E4 in closed form with constant scalars; exact window templates (hom_space, Rickard) with disk cache; polynomiality certificates; reduced coordinates of idempotent images | `coset.rs`, `sexpr.rs`, `smap.rs`, `tmpl.rs` | acceptance tests only |
| M4b | `HomProvider` with backends `Exact` and `Search`; stacked-point certification with M(p,q,d) from the oracle; genuineness bookkeeping; the G1 target enumerator and harness | `homprov.rs` | **G1** |
| M5 | rung-by-rung decompositions B_p ⊗ R = ⊕ B_q{s} as global maps (ι, π) with scalar intersection forms and orthogonalization; idempotents e_p = P(a); per-coset reduced coordinates; disk cache | `sindec.rs` | acceptance tests only |
| M6a | R-full localized scan per point (typed complex, rungwise colored Rickard insertion, Markowitz scalar elimination with an incremental heap and the streaming rule); multi-point closure data (contract v2) in engine2 style; closure by CL1. Laptop envelope: 2 strands, (3,3) up to T(2,5), uncolored 3 strands, the 3-strand Markov family and 4_1. Other (2,2,2) knots are cloud jobs (§6.6) | `sengine.rs`, `src/bin/sscan.rs` | G2 is evaluated at the end of M6a |

### 0.2 Measured findings that shape this spec

**F1. Idempotents are almost never needed in M6a's scope** [V-P2].

The canonical (minimal-flow) ladder L_p of `schur` has character exactly b_p for every type of:
- (2,2), (3,3), (2,1,2), (1,1,1), (1,1,1,1), (2,2,1,1) and (2,2,2);
- all 63 types of the five blocks a (2,2,2) scan visits.

So e_p = id throughout M6a. The exceptions:
- (3,3,3): 1 of 55 types (table 0,1,2/1,1,1/2,1,0, flows 10 500 against rank 8 820, dim End^0(L_p) = 2); 3 of 219 over the 7 blocks a scan visits;
- (2,2,2,2): 6 of 282 (dim End^0(L_p) <= 11, flows/rank <= 1.39); 25 of 1 110 over the visited blocks.

M5 still implements and tests e_p = P(a). R1's figure (sum of flows 1 171 against sum of ranks 1 045 for (2,2,2)) was for BFS ladders; the minimal-flow ladders of M3a remove the overhead.

**F2. The exact backend cannot produce the (2,2,2) decompositions, nor the (3,3) ones** [V-P3, P12; M6.0 report].

`sbim::hom_space` timings on decomposition Hom spaces Hom^s(L_q, L_p ⊗ R):

| case | fibres L_q × X | unknowns | time |
|---|---|---|---|
| (2,2,2), 6 × 36, degree 4 | 6 × 36 | 4 084 | 0.105 s |
| (2,2,2), 20 × 60, degree 2 | 20 × 60 | 21 694 | > 150 s (killed) |
| (2,2,2), 36 × 72, degree −1 (Hom is 1-dimensional) | 36 × 72 | 12 420 | > 150 s (killed) |
| (2,2,2), 24 × 432, degree −1 | 24 × 432 | 55 945 | > 120 s (killed, 463 MB) |
| (3,3), 15 × 36 (M6.0) | 15 × 36 | 2 524 | 48 s for the whole decomposition |
| (3,3), 30 × 90 (M6.0) | 30 × 90 | - | > 7 min |

Of the 358 (2,2,2) decomposition targets, 173 have fibre product <= 1e3 and 185 lie above it. Exact closure generators (`cclose::rlambda_bases`) on (2,2,2) canonical ladders: 0.25-40 s at fibre 36, 2.5-145 s at fibre 60; the four types of fibre 72-216 did not finish within a 400 s budget [V-P9].

**Consequence.** The 3-strand part of M6a and every (3,3) braid beyond the kink require the Search provider (or M4c). This confirms SYNTHESIS §1.4 and gate G1.

**F3. Rung-by-rung decomposition equals half-step adjunction in size** [V-P4].

For every summand B_q{s} of B_p ⊗ R, the graded dimensions of Hom(B_q, B_p ⊗ R) equal those of Hom(B_{q'}, B_{p+}) up to one shift. Here B_{p+} = B_p ⊗ split and B_{q'} = B_q ⊗ reverse split, in the block with one more upright. This holds on all 352 summand cases of (2,2), (3,3), (1,1,1), (2,2,2) and (2,1,2).

So the half-step pipeline does not shrink the Hom problems. It only re-expresses them, and it costs a zero-label upright and a padded `Schur`.

**Decision:** M5 decomposes rung by rung, as M6.0 did and validated. Adjunction becomes one candidate family of the Search provider (§4.4, tier S5).

**F4. Extreme decomposition maps are unique** [V-P1].

A summand family q of B_p ⊗ R occupies degrees s_lo..s_hi. The lowest-degree inclusion ι_{s_lo} ∈ Hom^{s_lo}(L_q, X) and the highest-degree projection π_{s_hi} ∈ Hom^{−s_hi}(X, L_q) are 1-dimensional:
- in 156 of the 158 (2,2,2) families;
- in all 16 (3,3) families;
- in all 6 (2,2) families.

The multiplicity at s_lo is always 1. A 1-dimensional space is certified by one nonzero genuine map, which is a single-map statement (one point suffices). This yields the provider strategy "extremes plus dots" (§4.4, tiers S1-S2).

**F5. Dense localization is cheap at (2,2,2), infeasible at (3,3,3)** [V-P5].

`local.rs` (dense L and its inverse on the whole ladder, one point):
- 4 ms for all 21 (2,2,2) canonical ladders (max fibre 216, at most 8 flows per coset);
- 22 ms for all 21 X = L_p ⊗ E_0^(1) (max fibre 432, at most 20 flows per coset);
- (3,3,3): 10.6 s at fibre 4 320, and above 1.2 GB (killed) at fibres 8 820 and 10 500.

Coset-local evaluation (M4a) is mandatory for scale-out. In M6a's scope `local.rs` remains an affordable oracle.

**F6. Closure data and per-point cost decide where 3-strand knots run. Of the (2,2,2) knots, only 4_1 runs on the laptop** [V-P6, P7; V-crit c1, c3, c4].

- **Weights.** Gram weights of the target types reach 12 for (2,2,2), i.e. 1 218 unknowns. F weights reach or exceed the Gram weight on exact 2-strand data: (2,2) T(2,7) has F weight 5 against Gram weight 4, and a component of degree 10. For (2,2,2) the e-monomial counts at weights 10..16 are 630, 882, 1 218, 1 638, 2 178, 2 838, 3 663. F by rescans needs emono(w_F) + 8 scans (§6.5).
- **Per-point scan cost.** Proxy from engine2 (`prepare_closure_data`, one thread): 0.8-2.3·10^-5 s per unit of decategorified peak Σ S. My re-run on 12a_146 gave 0.015 s per point at peak Σ S 1 173, i.e. 1.3·10^-5. This is a lower bound for the colored engine, whose per-unit cost includes the ι/π sandwich.
- **Where each case runs** (§6.6, two threads, one prime, one point set):
  - 4_1 needs 1-6 min of rescans: laptop.
  - 5_2, 6_2 and 6_3 need 7-42 min, and 8_19 needs 20-103 min: cloud under M6a's 600 s rule.
  - A 10-crossing (2,2,2) word (peak Σ S 2.9·10^6) needs 4-20 h: cloud, and only after M6c.
- **Closure.** Unreduced closure of (2,2,2) data exceeds the exact oracle's 768 MB cap at N = 6 for every target (4_1: about 2.1 GB at the measured Gram density 0.20), and at N = 5 for 8_19 and beyond. The Δ-tower is worse (6-55 GB). Reduced closures stay below 0.3 GB up to N = 6.

**F7. G1 target sizes** [V-P1, P11]:
- (2,2,2): 126 (type, rung) pairs, of which 28 are identity shortcuts; 158 summand families; 358 (q, s) targets with Σ mult = 404; closure: 21 types with Σ c_diag(1) = 48 ins and 48 outs (max 8 per type, max degree 12);
- (3,3): 18 pairs (6 identity), 16 families, 72 targets, Σ mult = 121; closure: 4 types, 8 ins and 8 outs;
- sampling populations for G1(b): (2,2,2,2) 2 191 non-identity pairs (9 255 targets); (3,3,3) 429 (3 151 targets).

**F8. One normalization makes every closed-form constant 1** [V-R1, 0.3 s release].

Cases checked: every E2 pair (a, b) with a, b ≥ 1 and a + b ≤ 6, every E3 digon (c; b) with 1 ≤ b < c ≤ 6, in both orientations (60 cases). In each case the exact `hom_space` map is 1-dimensional, and the following hold at 3 points:
- Eps, EpsL, Iota and IotaL equal the closed forms with value 1;
- Eta = e(B, A) at S = B, which is (−1)^{ab} e(A, B);
- EtaL = e(A, B) at S = A;
- Trace = Σ_S m(S)/e(S, A∖S), and TraceL is its mirror.

§3.4.3 defines the forms this way. This settles the earlier open question on the signs of E2/E3.

**F9. The per-type closure oracle is exact and q-graded** [V-R2].

For the exact xengine closure data of (2,2) σ^3, (2,1) σ^2, (3,1) σ^2 and (1,1) σ^3, the closure of the one-object complex (p, 0, 0) (no edges, no crossings, same TypeBasisC) equals `moy::closure_value` of L_p as a q-polynomial. This holds for every type and N = max label..5 (36/36, including the zeros when a label exceeds N). It becomes test T0 (§6.7).

**F10. Storage layout is a first-order cost** [V-crit c2, c6 re-run].
- Real (2,2,2) coset blocks hold 2.26 entries on average, and 112 of 175 (type, table) shapes have one flow per coset. `BTreeMap<u64, Vec<Vec<F>>>` costs 59.6 B per entry. A flat vector with shared u32 offsets costs 11.2 B.
- A FlowIndex with hash maps costs 225-285 B per flow. The packed layout of §3.3.2 costs 8 B per flow and rung.
- Both packed layouts are now binding.

**F11. Witness-only detection has a counting limit** [V-crit witness re-run].

With one witness coset per table, the bound (M+2)·E_w ≥ dim Hom^d holds for d ≤ 4 in (2,2), (3,3), (2,1,2) and (2,2,2), so M6a's scope is safe. It fails from:
- d = 11 in (2,2) (dim 69, M 20, E_w 3);
- d = 12 in (2,1,2);
- d = 15 in (3,3) and in (2,2,2) (542 of 4 923 (pair, d) with d ≤ 24).

§4.3 asserts the bound and falls back to full orbits when it fails.

### 0.3 Decisions (binding for M4a-M6a)

1. **Objects.** B_p is the `schur` minimal-flow ladder L_p plus an optional idempotent e_p, with e_p = None when ladder_char(L_p) = b_p. Ladders are expressions on n uprights with zero labels allowed and no label bound other than K. The (n+1)-upright padded model is supported by `sexpr` but not used by the M6a pipeline.
2. **Maps.** Every map stored by M5 and M6a's decomposition data is a node of a global, point-independent expression DAG (`smap::Arena`). Its leaves are certified templates, closed forms with constant scalars, and dots. Coefficients are constants of F_p. Numbers exist only inside a per-point `Eval`.
3. **Decompositions.**
   - Rung by rung (F3).
   - Certified by the scalar intersection form at a top-stratum coset. This is a single-map statement: one point decides, a second point re-verifies the identical matrix.
   - Completed by orthogonalization, global checks Π∘Ι = id and ΣιΠ = id at two points on all cosets, and multiplicities equal to `schur::term_mult`.
4. **Hom provider order.** The order depends on the labels (§4.5). `Leaves` (M4c) is used only after a G1/G2 failure.
   - Max label ≤ 2: `Exact` (admitted by caps or by a measured cost table, with a per-call ceiling of 10 s), then `Search` (tiers S0-S6, §4.4).
   - Max label ≥ 3: Search S0-S2, then Exact (cost table only), then Search S3-S6.

   Every decision about linear independence over Hom^d uses stacked evaluation:
   - at the witness cosets with M(p,q,d)+2 points when (M+2)·E_w ≥ dim Hom^d;
   - otherwise at the full W_λ-orbits of the witnesses, or with ⌈dim/E_w⌉+2 points (§4.3, F11).
5. **Engine.**
   - R-full at one point per scan.
   - Scans at many points (engine2 style), run in parallel with an immutable shared `DecData` built during the first scan.
   - Pivots and object lists must be identical at every point (asserted).
   - `Arena`, `TmplLib` and `DecData` are `Send + Sync` (`Arc` inside, no interior mutability), asserted at compile time. `Rc` lives only in the per-thread `Eval` memo.
6. **Closure data.**
   - Closure generators: the one-point determinant criterion (C2 M2), re-verified at every scan point.
   - Gram: from generator-only points by dense interpolation over R^λ e-monomials, with 8 residual points, up to 2 500 unknowns.
   - F: by rescans at one point set fixed in advance (emono(w_F) + 8 points). Each entry keeps only the samples it needs at its own weight. Equal-weight entries are interpolated in one blocked product. A time rule and a memory rule turn large cases into cloud jobs.
   - Closure: `cclose::close_exact_colored`, or `cclose::close_delta_tower` (G4 certified for all a ≤ 5, N ≤ 10) when the exact oracle refuses. Auto reports both estimates when it refuses.
   - Laptop policy: reduced N ≤ 6. Unreduced at every N where `close_plan`'s exact estimate is under the 768 MB cap: N ≤ 4 for the (2,2,2) knots, N = 5 for 4_1, and N ≤ 6 for the small Markov-family words.
7. **Laptop envelope under M6a** (F6, §6.6):
   - laptop: 2 strands, (3,3) up to T(2,5), uncolored 3 strands, the 3-strand Markov family (T9) and 4_1;
   - cloud: every other (2,2,2) knot;
   - ≥ 9 crossings at (2,2,2) also requires M6c, or another F path without rescans.
8. **Storage** (F10).
   - A localized map is one contiguous `Vec<F>` over a shared `Layout` (u32 offsets, dense row-major blocks).
   - `FlowIndex` is a packed mask array.
   - Every provider call runs in an Arena scratch region. Only certified results survive `commit`.
   - Flow indices are built only for words that are evaluated.
9. **Closed forms** are normalized so that every constant is 1 (F8). `assert_genuine` rejects any (kind, labels) without a recorded check.

---------------------------------------------------------------------------------------------------

## 1. Interfaces and changes to existing code

### 1.1 Data flow between the milestones

```
schur (M3a): types, b_p, ladder_minflow, term_mult, hom_gd, c_diag, rank, crossing_terms, decat
   │
M4a  coset::CosetSpace ── sexpr::{Expr, FlowIndex} ── tmpl::TmplLib (exact window maps, Rickard, cache)
        └──────────────── smap::{Arena (immutable DAG), Eval (per point memo), ElemKind E1-E4, Reduced}
   │
M4b  homprov::{HomProvider: Exact | Search | Chain}, certify_full, m_points, detecting_cosets, g1_targets
   │       (adds certified MapIds to the Arena; never numbers)
M5   sindec::{SIndec (builder), DecData (frozen, Arc), Decomp, idempotent, reduced}
   │       (calls HomProvider::{pairs, end0}; decomposition data are Arena nodes)
M6a  sengine::{sscan, PEngine (one point), SCx<LocMap>, closure data} ── cclose (CL1) ── Poincaré
```

### 1.2 The binding interface types (defined in M4a/M4b; used unchanged downstream)

| type | defined in | consumed by | invariant |
|---|---|---|---|
| `Expr`, `ExprId` | sexpr / smap | all | rung word with bottom labels; zero labels allowed; interned (equal words give equal ids) |
| `FlowIndex` | sexpr | smap, M5, M6a | packed: flows sorted by (top coset rank, lexicographic moved masks), all masks in one `Vec<u64>`; this is the row and column order of every block |
| `MapId`, `MapInfo` | smap | all | global map with (src, tgt, normalized degree e, certified flag, provenance); e = u + shift(tgt) − shift(src) as in `xengine::hom_norm` |
| `Block<F>`, `Layout`, `LocMap<F>` | smap | M5, M6a | `Block`: one coset, dense row-major rows × cols (rows = tgt), used in the Eval memo. `Layout`: cosets of the joint support with (rows, cols, u32 offset), shared through `Arc`. `LocMap`: one contiguous `Vec<F>` over a `Layout`, with zero blocks stored as zeros |
| `Point<F>` | smap | all | K pairwise distinct roots, deterministic from (seed, index); `id` unique within a run |
| `Eval<'a, F>` | smap | M4b, M5, M6a | per-point memo over an immutable `&Arena` |
| `HomProvider<F>` | homprov | M5, M6a | returns certified MapIds or `Vec<Missing>`; never partial |
| `DecData<F>` | sindec | M6a | frozen after the first scan; `Send + Sync`; shared read-only as `Arc<DecData>` by parallel scans; lookups return `Arc` or references |

### 1.3 Changes to existing modules (complete list; only pub-ification or additive)

| file | change | milestone |
|---|---|---|
| `src/lib.rs` | add `pub mod coset; pub mod sexpr; pub mod smap; pub mod tmpl;` | M4a |
| `src/lib.rs` | add `pub mod homprov;` | M4b |
| `src/lib.rs` | add `pub mod sindec;` | M5 |
| `src/indec.rs` | make `ptrim`, `pmul`, `psub`, `pdivrem`, `pxgcd` `pub` (polynomial helpers for e = P(a)); no other edit | M5 |
| `src/lib.rs` | add `pub mod sengine;` | M6a |
| `src/engine2.rs` | make `par_map` `pub` (thread pool honoring KRHOM_THREADS); no other edit | M6a |

Nothing else in existing modules may change. New binaries:
- `src/bin/g1.rs` (M4b, the G1 harness);
- `src/bin/sscan.rs` (M6a).

`src/bin/krhom.rs` integration is M6b. New test files: `tests/m4a.rs`, `tests/homprov.rs`, `tests/sindec.rs`, `tests/sengine.rs`, plus data under `tests/data/sengine/`.

Reused as they are:
- `local.rs` (`fixed_points`, `loc_matrix`, `inverse`, `localize_map`, `esym`, `root_offsets`: oracle and template localization);
- `sbim.rs` (`SBim`, `hom_space`);
- `cube.rs` (`LocalCx`, `local_cx`, `apply_local_poly`, already pub);
- `schur.rs`, `cclose.rs`, `xengine.rs` (oracle for M6a), `engine2.rs` (oracle), `dense.rs`, `linsys.rs`, `symfun.rs`, `mpoly.rs`, `pmat.rs`, `field.rs`, `laurent.rs`, `crossing.rs`, `web.rs`.

---------------------------------------------------------------------------------------------------

## 2. Conventions (all four milestones)

1. **Roots.** K = Σ λ_i <= 64 (masks are u64; exact bridges need K <= 16 because of `mpoly::MAXVARS`). Bottom block i owns roots off[i]..off[i]+λ_i (`local::root_offsets`).
   - A point x ∈ F^K has pairwise distinct coordinates.
   - `points(k, seed, count)` draws x_j uniformly from 1..2·10^9 with `StdRng::seed_from_u64(seed)`, redrawing on a collision.
2. **Configurations and cosets.** A configuration is `Conf(Vec<u64>)`, one root mask per upright (zero labels give mask 0). A top coset of the block (λ, μ) is a configuration with |mask_j| = μ_j. `CosetSpace::rank` is the multinomial number system.
   - The table of a coset: M_ij = popcount(B_i & T_j).
   - The witness of a table t: `rep(t)` (§3.3.1).
3. **Flows.** A flow of an expression is the vector of moved masks, one per rung, in global root numbering. Flows are grouped by top coset and sorted lexicographically inside a coset. This order is the row and column order of every `Block`.
4. **Degrees.**
   - Unnormalized u (as in `sbim::hom_space`); normalized e = u + shift(tgt) − shift(src), with shift(L) = −Σ_rungs r(b − r).
   - A component B_p{t} → B_q{t'} has e = t − t' (SYNTHESIS §1.4). Asserted on every stored component: total q-degree 0.
5. **Signs and order of the crossing complex.** Exactly as `xengine` (validated in M6.0):
   - terms from `schur::crossing_terms(c, d, positive, i)`, equal to `crossing_cx(.., n = c+d)` rungs;
   - object (p, h, t) means q^t L_p;
   - d_C ⊗ id carries no sign; id ⊗ d_j carries the Koszul sign (−1)^h of the source object;
   - positive crossings map term j to j+1, negative ones term j+1 to j;
   - closure via `cclose::global_q`.
6. **Arithmetic.** F_p with p = 2^31−1 (`F31`). Final results are also computed over F_4294967291 (`F32b`). Exact templates are computed per prime.
7. **Constants only.** A `Lin` node takes constant coefficients only. No map is ever rescaled "so that an entry is 1 at the current point" (C1 I7). Scalars that fix closed-form normalizations are found once, by comparison with exact maps, and stored as constants.
8. **Single-map versus space statements** (C1 I1; SYNTHESIS §1.4).
   - One point (plus a re-check at a second point) may decide: the zero test of one map; top-stratum scalars; intersection forms (matrices of scalars); pivots; the minimal polynomial of one element; the nonzero test certifying a 1-dimensional space.
   - Everything else (the rank of a family of maps in Hom^d, coordinates, structure maps, certification of a basis with dim > 1) uses stacked evaluation at M(p,q,d) + 2 points, plus a residual at 2 fresh points.

---------------------------------------------------------------------------------------------------

## 3. M4a, the localized substrate (FIRST TASK, ready to implement)

### 3.1 Goal, scope, permissions

**Goal.** Provide the substrate on which every later milestone builds:
- genuine bimodule maps between ladder expressions, represented as an immutable expression DAG;
- evaluation one top coset at a time at numeric points;
- closed-form elementary maps with constant scalars;
- exact window templates with a disk cache;
- certificates that make genuineness checkable.

All of it is validated against the exact oracles `sbim::hom_space`, `local::localize_map` and `schur`.

**Non-goals.** No Hom-basis search (M4b), no decompositions or idempotent construction (M5), no scan or closure (M6a), no backward flow search (forward enumeration suffices up to (2,2,2,2): 51 840 flows), no padded half-step pipeline.

**Permissions.**
- New files of §3.2, plus 4 lines in `src/lib.rs`. Change nothing else.
- Work on a new branch `colored-m4a` from `dev`; commit only as the orchestrator instructs.
- No stash, reset or checkout of other branches.

**Resource rules.**
- `export PATH="$HOME/.cargo/bin:$PATH"`.
- Every cargo command is prefixed with `KRHOM_THREADS=2 CARGO_BUILD_JOBS=2 RUST_TEST_THREADS=2 CARGO_TARGET_DIR=<scratchpad>/colored_design/target_m4a nice -n 19`.
- One heavy process at a time, under 1 GB RSS, watched by a watchdog (poll `ps -o rss=`, kill above 1.2 GB).
- Anything predicted above 10 min or 1 GB is `#[ignore]` and listed as a cloud job.

### 3.2 Files

| file | content |
|---|---|
| `src/coset.rs` | `Roots`, `Conf`, `CosetSpace` (§3.3.1) |
| `src/sexpr.rs` | `Expr`, `FlowIndex`, bridges to `SBim` and `schur` (§3.3.2) |
| `src/tmpl.rs` | `Window`, `Template`, `WinLoc`, `TmplLib` (hom_space and Rickard templates, disk cache; no interior mutability, the localization memo lives in `Eval`), `poly_certificate`, `tensor_exact` (§3.3.3, §3.4.4, §3.4.5) |
| `src/smap.rs` | `Point`, `Block`, `gemm_acc`, `Layout`, `LocMap`, `ElemKind`, `MapNode`, `MapInfo`, `Arena` (with `mark`/`commit`), `Eval`, `Reduced` (§3.3.4) |
| `src/lib.rs` | 4 `pub mod` lines |
| `tests/m4a.rs` | acceptance tests A1-A16 (§3.5) |

### 3.3 Types and functions (signatures are binding; internals free)

#### 3.3.1 `coset.rs`

```rust
//! Cosets W/W_mu as top configurations of roots (D1 §2.2), tables, witnesses, W_lambda action.
use crate::schur::Table;
pub type Mask = u64;

#[derive(Clone, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub struct Conf(pub Vec<Mask>);                       // one root mask per upright

#[derive(Clone, Debug)]
pub struct Roots { pub lam: Vec<u8>, pub off: Vec<usize>, pub k: usize }
impl Roots {
    pub fn new(lam: &[u8]) -> Roots;                  // asserts K <= 64
    pub fn bottom(&self) -> Conf;                     // block i on upright i
    pub fn block_of(&self, root: usize) -> usize;
}

#[derive(Clone, Debug)]
pub struct CosetSpace { pub roots: Roots, pub mu: Vec<u8> }   // mu.len() == lam.len(), zeros allowed
impl CosetSpace {
    pub fn new(lam: &[u8], mu: &[u8]) -> CosetSpace;  // asserts sum mu = K and equal lengths
    pub fn count(&self) -> u64;                       // K! / prod_j mu_j!
    pub fn rank(&self, c: &Conf) -> u64;              // bijection onto 0..count (multinomial number system)
    pub fn unrank(&self, r: u64) -> Conf;
    pub fn table(&self, c: &Conf) -> Table;           // M_ij = popcount(B_i & T_j)
    pub fn rep(&self, t: &Table) -> Conf;             // witness: for j = 0.., for i = 0..: the first M_ij unused roots of B_i
    pub fn cosets_of(&self, t: &Table) -> Vec<u64>;   // all N(t) cosets of table t, ascending rank
    pub fn act(&self, sigma: &[usize], c: &Conf) -> Conf;   // sigma: permutation of 0..K that preserves every bottom block
}
```

Rank algorithm: for j = 0..n−1, the rank of mask_j among the μ_j-subsets of the roots not used by masks 0..j−1 (colexicographic combinatorial number system), combined in mixed radix with radices C(remaining, μ_j). Use u128 internally and assert that the count fits u64.

#### 3.3.2 `sexpr.rs`

```rust
//! Ladder expressions with zero labels (no truncation) and their flows, grouped by top coset.
use crate::coset::{Conf, CosetSpace, Mask, Roots};
use crate::web::Rung;

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Expr { pub bottom: Vec<u8>, pub rungs: Vec<Rung> }
impl Expr {
    pub fn new(bottom: &[u8], rungs: &[Rung]) -> Result<Expr, String>;   // giver >= r at every rung; r > 0
    pub fn levels(&self) -> Vec<Vec<u8>>;          // labels after 0..=len rungs
    pub fn top(&self) -> Vec<u8>;
    pub fn shift(&self) -> i32;                    // -sum_rungs r(b - r), b = giver label (= SBim::shift)
    pub fn nflows(&self) -> u128;                  // = schur::flows(bottom, rungs)
    pub fn sbim<F: crate::field::Field>(&self) -> crate::sbim::SBim<F>; // exact bridge, K <= 16
    pub fn concat(&self, more: &[Rung]) -> Result<Expr, String>;
}

#[derive(Clone, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub struct Flow(pub Vec<Mask>);                    // owned copy of one flow (templates, tests, oracles)

/// Packed flow index (binding layout, F10): flows sorted by (top coset rank, lexicographic moved
/// masks); flow f (global id) occupies masks[f*nr .. (f+1)*nr]. No hash maps: heap memory is
/// 8·nr bytes per flow + 12 bytes per coset of the support.
pub struct FlowIndex {
    pub expr: Expr,
    pub top: CosetSpace,                           // block (bottom, top labels)
    pub nr: usize,                                 // rungs = masks per flow
    pub cosets: Vec<u64>,                          // support, ascending rank
    pub start: Vec<u32>,                           // len = cosets.len() + 1; flows at cosets[i] are start[i]..start[i+1]
    pub masks: Vec<Mask>,                          // len = nr * nflows
}
impl FlowIndex {
    pub fn build(expr: &Expr) -> FlowIndex;        // forward enumeration as local::fixed_points, then one sort
    pub fn range(&self, c: u64) -> std::ops::Range<u32>;   // global ids of the flows at c; empty outside the support
    pub fn flow(&self, f: u32) -> &[Mask];
    pub fn find(&self, c: u64, moved: &[Mask]) -> Option<u32>;  // row index inside coset c (binary search)
    pub fn support(&self) -> &[u64];               // ascending
    pub fn conf_at(&self, f: u32, level: usize) -> Conf;   // configuration after `level` rungs
    pub fn max_block(&self) -> usize;
    pub fn bytes(&self) -> usize;                  // heap bytes (A16)
}
/// Flows per coset against the singular Deodhar coefficient c_{table(c)}(ladder_char)(1),
/// and every coset of every table with nonzero coefficient present.
pub fn deodhar_check(fi: &FlowIndex, s: &mut crate::schur::Schur) -> Result<(), String>;
```

#### 3.3.3 `tmpl.rs`

```rust
//! Exact window maps (sbim::hom_space, cube::local_cx), localized per window root tuple,
//! applied locally inside larger ladders; polynomiality certificates; disk cache.
use crate::field::Field;
use crate::pmat::PMat;
use crate::web::Rung;

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Window { pub labels: Vec<u8>, pub rungs: Vec<Rung> }   // rungs relative to the window's first upright

pub type TmplId = u32;
pub enum TmplSource { HomSpace { index: usize }, Rickard { c: u8, d: u8, positive: bool, j: usize } }
pub struct Template<F: Field> {
    pub src: Window, pub tgt: Window,
    pub deg: i32,                                  // normalized
    pub phi: PMat<F>,                              // exact over R^{src.labels}, unnormalized degree
    pub source: TmplSource,
    pub src_flows: Vec<crate::sexpr::Flow>,        // window-local moved masks, sorted as in FlowIndex
    pub tgt_flows: Vec<crate::sexpr::Flow>,        // (point-independent; computed once)
}
/// Localized window map at one root tuple, sparse by column (CSC over the sorted window flows):
/// the nonzeros of column j (source flow j) are (rows[k], vals[k]) for k in col_start[j]..col_start[j+1].
/// Window maps are block diagonal over window top configurations, so columns hold few entries.
pub struct WinLoc<F> { pub nrows: u32, pub col_start: Vec<u32>, pub rows: Vec<u32>, pub vals: Vec<F> }
impl<F: Field> WinLoc<F> {
    pub fn dense(&self) -> Vec<Vec<F>>;            // rows = tgt flows, cols = src flows; tests and oracles only
}

/// Send + Sync (A16): no Rc, RefCell or Cell inside.
pub struct TmplLib<F: Field> { /* templates, (src,tgt,deg) -> ids, optional cache dir */ }
impl<F: Field> TmplLib<F> {
    pub fn new(cache_dir: Option<std::path::PathBuf>) -> Self;     // cache file names carry F::name()
    pub fn get(&self, t: TmplId) -> &Template<F>;
    /// exact basis of Hom^deg(src, tgt) (normalized degree), computed once, cached on disk
    pub fn hom_basis(&mut self, src: &Window, tgt: &Window, deg: i32) -> Vec<TmplId>;
    /// untruncated Rickard complex of a crossing on labels (c, d): term windows and d_j templates
    pub fn rickard(&mut self, c: u8, d: u8, positive: bool) -> (Vec<Window>, Vec<TmplId>);
    /// localization at window roots y (upright by upright, ascending global root id inside an upright);
    /// not memoized here (the caller, `Eval`, memoizes per point)
    pub fn localize(&self, t: TmplId, y: &[F]) -> WinLoc<F>;
    pub fn save(&self) -> std::io::Result<()>;                       // header "krhom-tmpl v1 field=.. "
}
/// Exact id_P (x) phi (x) id_Q as a PMat (for tests and certificates): phi acts on the uprights
/// first..first+w of the top of `prefix` (substitution of the prefix's top operators, generalizing
/// cube::apply_local_poly to w >= 2 uprights), then kron_id(suffix_fsize).
pub fn tensor_exact<F: Field>(prefix: &crate::sbim::SBim<F>, first: usize, t: &Template<F>, suffix_fsize: usize) -> PMat<F>;
/// Polynomiality certificate (C3 H3, §3.4.5).
pub struct CertReport { pub points: usize, pub max_monos: usize, pub distinct_degrees: usize, pub entries: usize, pub redraws: usize, pub secs: f64 }
pub fn poly_certificate<F: Field>(src: &Window, tgt: &Window, deg: i32,
        f: &mut dyn FnMut(&[F]) -> WinLoc<F>, seed: u64) -> Result<CertReport, String>;
```

#### 3.3.4 `smap.rs`

```rust
//! Global (point-independent) maps between ladder expressions as an immutable DAG; evaluation one
//! top coset at a time at a numeric point, memoized per point.
use crate::coset::Conf;
use crate::field::Field;
use crate::sexpr::{Expr, FlowIndex};
use crate::tmpl::{TmplId, TmplLib};
use std::rc::Rc;
use std::sync::Arc;

pub type ExprId = u32;
pub type MapId = u32;

/// One coset block, dense row-major: rows = tgt flows at c, cols = src flows at c.
#[derive(Clone, Debug, PartialEq)]
pub struct Block<F> { pub rows: u32, pub cols: u32, pub a: Vec<F> }           // a.len() == rows * cols
impl<F: Field> Block<F> {
    pub fn zeros(rows: usize, cols: usize) -> Self;
    pub fn identity(n: usize) -> Self;
    pub fn get(&self, i: usize, j: usize) -> F;
    pub fn mul(&self, o: &Block<F>) -> Block<F>;                              // self · o, dimension-checked
    pub fn to_rows(&self) -> Vec<Vec<F>>;                                     // tests and oracles only
}
/// c[m×n] += alpha · a[m×k] · b[k×n] on row-major slices, without allocation. Every block product
/// of M4a-M6a goes through it.
pub fn gemm_acc<F: Field>(c: &mut [F], a: &[F], b: &[F], m: usize, k: usize, n: usize, alpha: F);

/// Block shapes of all maps between two flow spaces over the cosets of their joint support. Computed
/// once per (src, tgt): expressions in M4a/M5, (src type, tgt type, top block) in M6a. Shared by every
/// map with these ends.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Layout { pub cosets: Vec<u64>, pub rows: Vec<u32>, pub cols: Vec<u32>, pub off: Vec<u32>, pub len: u32 }
impl Layout {
    pub fn new(cosets: Vec<u64>, rows: Vec<u32>, cols: Vec<u32>) -> Layout;   // offsets = prefix sums; asserts len < 2^32
    pub fn of(src: &FlowIndex, tgt: &FlowIndex) -> Layout;                 // cosets in both supports; flow counts
    pub fn index(&self, c: u64) -> Option<usize>;                          // binary search in `cosets`
}
/// A map localized at one point (binding layout, F10): one contiguous vector over a shared layout.
/// Zero blocks are stored as zeros; their fraction is measured (§6.8).
#[derive(Clone, Debug)]
pub struct LocMap<F> { pub lay: Arc<Layout>, pub data: Vec<F> }
impl<F: Field> LocMap<F> {
    pub fn zeros(lay: Arc<Layout>) -> Self;
    pub fn block(&self, i: usize) -> &[F];                                  // i indexes lay.cosets
    pub fn block_mut(&mut self, i: usize) -> &mut [F];
    pub fn is_zero(&self) -> bool;
    pub fn nonzero_blocks(&self) -> usize;
    pub fn bytes(&self) -> usize;                                           // 8·len + O(1); the shared layout is not counted
}

#[derive(Clone, Debug)]
pub struct Point<F> { pub id: u32, pub x: Vec<F> }
pub fn points<F: Field>(k: usize, seed: u64, count: usize, first_id: u32) -> Vec<Point<F>>;

/// Closed-form elementary maps (R1 §2.4), all on a 2-upright window (u, u+1) starting at rung `at`.
/// Normalized so that each equals the unique exact map (constant 1, F8, §3.4.3).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ElemKind {
    /// E2 on window [E_u^(b), F_u^(b)] over labels (a, b) (all of upright u+1 merged and split
    /// back): Eps (window -> empty, degree +ab) and Eta (empty -> window, degree +ab, value e(B,A));
    /// mirrored window [F_u^(a), E_u^(a)] (all of upright u): EpsL, EtaL (value e(A,B))
    Eps, Eta, EpsL, EtaL,
    /// E3 digon [F_u^(b), E_u^(b)] over (c, 0), and mirrored [E_u^(b), F_u^(b)] over (0, c):
    /// Iota (empty -> digon, degree -ab), Trace (digon -> empty, degree -ab, Σ_S m(S)/e(S, A∖S))
    Iota, Trace, IotaL, TraceL,
    /// E4: identity by matching flows with equal intermediate configurations outside the
    /// window and equal top configuration; only for window pairs whose flows are in bijection
    /// with their top configurations (asserted)
    Assoc,
}

pub enum MapNode<F> {
    Id(ExprId),
    /// E1: multiplication by s_nu (Schur polynomial; nu = partition, empty = 1) of the roots on
    /// upright `upright` after `level` rungs; degree 2|nu|
    Dot { e: ExprId, level: usize, upright: usize, nu: Vec<u8> },
    Elem { kind: ElemKind, src: ExprId, tgt: ExprId, at: usize, upright: usize, scalar: F },
    /// id_prefix (x) template (x) id_suffix: src = P·A·Q, tgt = P·A'·Q, A at rungs at..at+|A|
    Local { src: ExprId, tgt: ExprId, at: usize, first: usize, tmpl: TmplId },
    Comp(MapId, MapId),                                         // Comp(g, f) = g ∘ f
    Lin(Vec<(F, MapId)>),                                       // constant coefficients only
    Kron(MapId, Vec<Rung>),                                     // f (x) id_rungs
}
#[derive(Clone, Debug)]
pub struct MapInfo { pub src: ExprId, pub tgt: ExprId, pub deg: i32, pub certified: bool, pub provenance: String }

/// Send + Sync (A16): no Rc, RefCell or Cell inside; immutable once frozen into DecData (M5).
pub struct Arena<F: Field> { /* exprs (interned), Option<FlowIndex> per expr, nodes, infos, TmplLib, elem constants */ }
#[derive(Clone, Copy, Debug)]
pub struct ArenaMark { /* node, expression and index counts at the mark */ }
#[derive(Clone, Debug, Default)]
pub struct ArenaBytes { pub exprs: usize, pub indices: usize, pub indexed_exprs: usize, pub nodes: usize, pub templates: usize }
impl<F: Field> Arena<F> {
    pub fn new(lam: &[u8], lib: TmplLib<F>) -> Self;
    pub fn expr(&mut self, e: Expr) -> ExprId;                   // interns only; no FlowIndex is built
    pub fn index(&mut self, id: ExprId) -> &FlowIndex;           // builds the FlowIndex on first call
    pub fn ensure_indices(&mut self, m: MapId);                  // indices of every expression reachable from m
    pub fn drop_index(&mut self, id: ExprId);                    // frees it; rebuilt identically on demand (used by the
                                                                 // scratch-region LRU; SIndec::freeze asserts that every
                                                                 // expression of a retained node has its index built)
    pub fn expr_of(&self, id: ExprId) -> &Expr;
    pub fn flows(&self, id: ExprId) -> &FlowIndex;               // panics "index not built" (Eval relies on ensure_indices)
    pub fn lib(&mut self) -> &mut TmplLib<F>;
    /// checks compatibility (equal bottoms/tops, window match, Lin terms with equal src/tgt/deg),
    /// computes the degree, builds the indices of src and tgt; `certified` = all leaves certified
    /// (Templates are; Elem iff its (kind, labels) is recorded; Dot always; Comp/Lin/Kron: all children)
    pub fn add(&mut self, n: MapNode<F>, provenance: &str) -> Result<MapId, String>;
    pub fn info(&self, m: MapId) -> &MapInfo;
    pub fn certify_elem(&mut self, kind: ElemKind, labels: (u8, u8), scalar: F);  // records a validated constant (A7, A8)
    pub fn elem_certified(&self, kind: ElemKind, labels: (u8, u8)) -> Option<F>;
    /// Scratch regions (used by every provider call, §4.4): `commit` drops everything created after
    /// `mark` (nodes, expressions, indices) except the nodes reachable from `keep`, with their
    /// expressions and indices. It re-interns them and returns their new ids in the order of `keep`.
    /// Ids below the mark are unchanged; templates are always kept.
    pub fn mark(&self) -> ArenaMark;
    pub fn commit(&mut self, mark: ArenaMark, keep: &[MapId]) -> Vec<MapId>;
    pub fn bytes(&self) -> ArenaBytes;
}

/// Per point and per thread (not Send): memo (MapId, coset) -> Rc<Block>, template/closed-form memo
/// (TmplId or Elem, root ids) -> Rc<WinLoc>, layout memo (src, tgt) -> Arc<Layout>.
pub struct Eval<'a, F: Field> { /* &Arena, point, memos */ }
impl<'a, F: Field> Eval<'a, F> {
    pub fn new(ar: &'a Arena<F>, pt: Point<F>) -> Self;
    pub fn point(&self) -> &Point<F>;
    pub fn block(&mut self, m: MapId, c: u64) -> Rc<Block<F>>;          // zero block of the right shape if outside support
    pub fn layout(&mut self, src: ExprId, tgt: ExprId) -> Arc<Layout>;  // Layout::of, memoized
    pub fn all(&mut self, m: MapId) -> LocMap<F>;                        // over layout(src, tgt)
    pub fn is_zero(&mut self, m: MapId) -> bool;                         // single-map zero test (all cosets)
    pub fn stacked(ar: &Arena<F>, ms: &[MapId], pts: &[Point<F>], cosets: &[u64]) -> Vec<Vec<F>>;  // one flattened row per map
    pub fn memo_bytes(&self) -> usize;
    pub fn clear(&mut self);
}

/// Reduced coordinates of an idempotent's image at one coset: P (flows × m), Q (m × flows),
/// Q P = I_m, P Q = E_c (D2's U_c, V_c).
pub struct Reduced<F> { pub p: Block<F>, pub q: Block<F>, pub m: usize }
pub fn reduced_of<F: Field>(e_c: &Block<F>) -> Result<Reduced<F>, String>;     // asserts E^2 = E
```

### 3.4 Algorithms

#### 3.4.1 Flows (`FlowIndex::build`)

1. Enumerate flows forward exactly as `local::fixed_points` (successive r-subsets of the giver mask), recording moved masks.
2. Compute the top configuration of each flow, then its coset rank.
3. Sort the flows by (coset rank, moved masks). Fill `cosets`, `start` and `masks`.
4. The temporary list is dropped. No per-flow allocation survives.

Cost O(fsize · rungs · log fsize). Memory 8·nr B per flow. Examples:
- (2,2,2) X with 432 flows and 7 rungs: 24 KB;
- (3,3,3) X with 31 500 flows and 7 rungs: 1.8 MB;
- (2,2,2,2) X with 51 840 flows and 11 rungs: 4.6 MB.

The hash-map layout cost 225-285 B per flow [V-crit c6].

#### 3.4.2 Evaluation `Eval::block(m, c)` (memoized per (m, c))

Blocks are dense row-major (`Block`). Products and sums go through `gemm_acc` on slices.
- `Id`: identity of size |range(c)|.
- `Dot`: diagonal; the entry of flow f is s_ν evaluated at the roots of `conf_at(f, level)[upright]` (via `local::esym` and Jacobi-Trudi).
- `Comp(g, f)`: block(g, c) · block(f, c), with a dimension check.
- `Lin`: Σ c_i block(m_i, c).
- `Kron(f, rs)`: the source of the node is src(f)·rs.
  - For a source flow (σ, w), where σ is a flow of src(f) at coset c' (the top of src(f)) and w the suffix moved masks: entry ((τ, w), (σ, w)) = block(f, c')[τ][σ].
  - Target flows (τ, w) are looked up with `find` in the target FlowIndex at c.
- `Local { at, first, tmpl }`, src = P·A·Q, tgt = P·A'·Q. For each source flow σ = (σ_P, σ_A, σ_Q) at c:
  1. Let C = conf_at(σ, |P|). The window roots are C[first..first+w], each upright's mask as an ascending root list; y = their values.
  2. Map σ_A to window-local masks, then to its column j in `Template::src_flows` (binary search). Take `localize(tmpl, y)` from the Eval memo, keyed by (tmpl, root id tuple).
  3. For every nonzero (i, v) in column j of the CSC `WinLoc`:
     - translate window target flow i (`tgt_flows[i]`) back to global masks; the target flow is τ = (σ_P, τ_A, σ_Q);
     - add v at (row(τ), col(σ)), where row(τ) = `find(c, τ)`.

     τ lies at coset c because the window top configurations agree (block diagonality).
  4. Assert that every such τ exists in the target index.
- `Elem`: the same pattern, with closed-form window blocks (§3.4.3) produced in the same CSC format instead of templates.
- `Eval::all(m)`: allocates one `LocMap` over `layout(src, tgt)` and copies each block into its slice.

#### 3.4.3 Closed forms (R1 §2.4; normalizations fixed here so that every constant is 1 [V-R1], re-validated in A7/A8)

Window roots: A (upright u) and B (upright u+1); e(S, T) = Π_{s ∈ S, t ∈ T} (x_s − x_t), so e(B, A) = (−1)^{ab} e(A, B).
- **Eps** on [E_u^(b), F_u^(b)] over (a, b), degree +ab: a row vector with value 1 at the flow whose second rung moves S = B, 0 elsewhere.
- **Eta**, degree +ab: a column with value **e(B, A)** at S = B, 0 elsewhere. Eps∘Eta = e(B, A).
- **EpsL/EtaL**: the mirror, on [F_u^(a), E_u^(a)] over (a, b). Value 1 at S = A, and **e(A, B)** at S = A.
- **Iota** on the digon [F_u^(b), E_u^(b)] over (c, 0), degree −ab with a = c − b. Flows are indexed by the moved b-subset S ⊂ A, and Iota has value 1 at every flow.
- **Trace**, degree −ab: d(m) = Σ_S m(S) / **e(S, A∖S)**.
- **IotaL/TraceL**: the mirror digon [E_u^(b), F_u^(b)] over (0, c), with C = the roots of upright u+1. Value 1 at every flow, and Σ_S m(S)/e(S, C∖S).
- **Assoc** (E4): identity by matching flows. Defined only for window pairs (A, A') with exactly one flow per top configuration on both sides and equal supports.

**Constants.**
- Each `Elem` node carries a constant `scalar` (any nonzero constant; default 1).
- M4a measures, per (kind, labels), the constant c with closed form = c · localize(exact unique map). It asserts c is equal at 3 points and records c with `certify_elem`.
- With the normalizations above, c = 1 for all 60 E2/E3 (kind-pair, labels) cases with a + b ≤ 6 or c ≤ 6, in both orientations [V-R1: `r1_e2all`, 3 points each, 0.3 s]. E4 constants are recorded by A8 (expected 1, [D]).
- `assert_genuine` rejects an `Elem` whose (kind, labels) has no recorded constant.

**Coverage.** M6a's scope ((2,2,2), (3,3), and their blocks) needs 2-upright windows with a + b ≤ 6. Among them:
- the (2,2,2) blocks (4,0,2) and (2,4,0) need E2 over (4,1), (4,2), (2,4) and the digon (4;3);
- (3,3) needs (5,1), (4,2), (3,3), (5;b) and (6;3).

A7 certifies every pair with a + b ≤ 6 by default. Pairs with a + b ≤ 9 ((3,3,3)) are certified as an ignored test before G1(b).

#### 3.4.4 Templates (`TmplLib`)

- `hom_basis(src, tgt, deg)`: build window SBims (`SBim::identity(2·K_w, labels).add_rungs(rungs)`) and call `sbim::hom_space(src, tgt, u)` with u = deg − shift(tgt) + shift(src). Store the PMats.
  - Disk cache, one text file per call, named by a hash of (field, src, tgt, deg); the header is checked on load.
  - A window is admitted only if fsize(src)·fsize(tgt) <= 1 600 (otherwise Err; M4b decides). [V-P3: 6·36 takes 0.1 s; 15·36 for (3,3) takes up to 48 s.]
- `rickard(c, d, positive)`: `cube::local_cx(Gen{i:0, positive}, c, d, c+d)`. Assert that the term rungs equal `schur::crossing_terms` (as `xengine::local`) and that each map is the unique degree-1 map. Cached on disk.
  - (3,3) takes about 120 s per sign per prime (M6.0), hence the cache.
- `localize(t, y)`: `local::localize_map(src_w, tgt_w, phi, y)`, then reorder rows and columns from `local::fixed_points` order to sorted window-flow order. Memoized by the caller (`Eval`).

#### 3.4.5 Polynomiality certificate (C3 H3)

Input: windows (src, tgt), predicted normalized degree, and a function y ↦ WinLoc (the closed form or composite under test).

1. Let T = 2 + max over entries (β', β) of |Mono_D|, where D = u + fdeg_src[β] − fdeg_tgt[β'] and Mono_D = `PRing::monos_of_deg(D)` of the window bottom ring. D < 0 means the entry must vanish, and this is checked at every point.
2. Draw **one** point set y_1..y_{T+2} (seeded). The same set serves every entry.
3. At each point: Φ(y) = L_tgt(y)^{-1} · f(y) · L_src(y). Here L is `local::loc_matrix` in `fixed_points` order, so the WinLoc is mapped back to that order.
4. For each distinct D:
   - build the Vandermonde V_D of Mono_D at the first |Mono_D| points and invert it **once**;
   - apply V_D^{-1} to the value matrix of all entries with that D (one blocked product);
   - if some V_D is singular, redraw the **whole** point set (seed + 1) and restart. After 3 redraws, return Err.
5. Residual: every entry's polynomial must reproduce its values at the remaining points (at least 2).
6. Pass iff every entry passes. The report gives T, max |Mono_D|, the number of distinct D, entries, redraws and seconds.

Cost.
- One inverse per distinct D: |Mono_D| ≤ about 160 at window fibre 400, under 0.05 s each [V-crit c3: m = 294 takes 0.02 s].
- A certificate is dominated by the T+2 evaluations of Φ: about 1-60 s, once per template kind.
- A per-entry factorization would cost about 10^3 s at fibre 400.

### 3.5 Tests (`tests/m4a.rs`; expected values or the oracle that produces them)

Default suite budget: < 60 s debug, < 300 MB. `#[ignore]` tests are run once each and reported.

**A1. Cosets.** `count` for (λ → μ):

| λ → μ | count |
|---|---|
| (2,2,2) → (2,2,2) | 90 |
| → (3,1,2) | 60 |
| → (4,0,2) | 15 |
| → (2,3,1) | 60 |
| → (2,4,0) | 15 |
| (3,3) → (3,3) | 20 |
| → (4,2) | 15 |
| → (5,1) | 6 |
| → (6,0) | 1 |
| (2,2) → (2,2) | 6 |
| (2,1,2) → (2,1,2) | 30 |
| (1,1,1) → (1,1,1) | 6 |
| (3,3,3) → (3,3,3) | 1 680 |
| (2,2,2,2) → (2,2,2,2) | 2 520 |

Further checks:
- `unrank(rank(c)) = c` for all cosets of the rows up to 1 680;
- Σ_{t ∈ schur::tables(λ, μ)} N(t) = count, and `cosets_of(t).len() = N(t)` (e.g. (2,2,2): 21 tables, 90 cosets);
- `table(rep(t)) = t` for every t;
- `act` preserves tables.

**A2. Flows = singular Deodhar** (R1 t2.py, 280 cases). 280 random ladders (seed s = 1..280): 2-4 uprights, bottom labels 0..3 with K <= 8, 1-6 rungs with r = 1..giver (drawn as `web::random_closed_ladder` does, without closing). `deodhar_check` passes for all 280.

In addition, all canonical ladders and all X = L_p ⊗ R of the five (2,2,2) blocks (63 + 126) and of (3,3) (10 + 18) pass.

Recorded maxima of `max_block`, which must match [V-P5]:
- (2,2,2) canonical: 8;
- (2,2,2) X = L_p ⊗ E_0^(1): 20;
- (3,3) X = L_p ⊗ E_0^(2): 18.

**A3. FlowIndex vs local.rs.** For the A2 ladders, the multiset of moved-mask vectors equals that of `local::fixed_points`, and the total equals `schur::flows`.

**A4. Whole-ladder exact bridge.**
- 20 seeded pairs (A, B) of random ladders with equal bottoms and tops, K <= 6, 2-3 uprights, fsize <= 40.
- For each pair, the lowest normalized degree d with hom_gd > 0 (`schur::hom_gd` on ladder characters); `hom_basis(A, B, d)` as whole-ladder templates (first = 0, empty prefix and suffix).
- `Eval::block` at every coset equals the corresponding block of `local::localize_map(A, B, Φ, x)` at 2 points.

**A5. Local application inside a ladder.** Prefix P (1-3 random rungs on 3 uprights, K <= 6), window A → A' on uprights (0,1) or (1,2), suffix Q (0-2 rungs). Windows:
- (1,1): [E(1), F(1)] → [] degree +1 and its reverse;
- (2,2): [E(1)] → [E(1)] dots;
- (1,1) and (2,1) Rickard d_0 at the top (Q empty) and in the middle (Q nonempty);
- one 3-upright window [E_0(1), E_1(1)] → itself, degree 2.

Requirements:
- `Local` eval = localize_map(tensor_exact(P, first, t, fsize(Q))) on all cosets at 2 points;
- `Kron(f, rs)` = localize_map(f_exact.kron_id(fsize(rs))).

**A6. Dots.** `Dot` at the top level equals the localized exact top operator e_{u,j} (`SBim::top[u][j-1]`) for s_ν = e_j, and the Jacobi-Trudi product for s_ν with |ν| <= 3. At a middle level it equals Local of the window identity times a dot (via `tensor_exact`).

**A7. E2 and E3, every label pair of M6a's scope** (default; [V-R1] `r1_e2all` runs all 60 cases in 0.3 s release, under 1 MB).
- **Cases:**
  - Eps, Eta, EpsL and EtaL for every (a, b) with a, b ≥ 1 and a + b ≤ 6 (15 pairs);
  - Iota, Trace, IotaL and TraceL for every digon (c; b) with 1 ≤ b < c ≤ 6 (15 digons).
- **(i) Constants.** For each form: dim Hom = 1 in its degree (`hom_basis`), and closed form = c · localize(exact) at 3 points, with the same c at all 3. **Expected: c = 1 in every case** (F8, [V-R1]). `certify_elem` records c. Any other value fails the test, naming kind and labels.
- **(ii) Products.** Exactly at all cosets:
  - Eps∘Eta = e(B, A) and EpsL∘EtaL = e(A, B);
  - Eps and EpsL have exactly one nonzero entry, and so do Eta and EtaL [V-R1].
- **(iii) Digons.**
  - Trace∘Iota = 0.
  - Trace∘Dot(s_ν of the moved subset S)∘Iota = 0 for every ν strictly inside the b × a box (|ν| < ab).
  - For ν = (a^b) it is a nonzero constant, the same at 3 points (Atiyah-Bott on Gr(b, c)). This constant is recorded (expected ±1).
- **(iv) Ignored, run before G1(b):** the same checks for a + b ≤ 9 and c ≤ 9 ((3,3,3) windows); time and RSS reported.

**A8. E4.** For (a, b, c) ∈ {(1,1,1), (2,1,1), (1,2,1), (1,1,2), (2,2,2)} and the pairs [E_0(b), E_1(c), E_0(c)] ↔ [E_1(c), E_0(b+c)] (and the F-mirrors): dim Hom^0 = 1 (hom_gd), `Assoc` = c · exact (constant over 3 points), and Assoc composed with its reverse = id.

**A9. Rickard templates.**
- (c, d) ∈ {(1,1), (2,1), (1,2), (2,2)} both signs: default.
- (3,1), (1,3), (3,2), (2,3): default if each `local_cx` call takes < 10 s; otherwise ignored.
- (3,3): ignored, then cached.

Checks:
- d_{j+1}∘d_j = 0 localized at all cosets at 2 points;
- block diagonality (no nonzero entry between different top cosets; checked on the dense localize_map output);
- uniqueness: dim Hom^1(W_j, W_{j+1}) = 1 by `schur::hom_gd` on the term characters, and dim Hom^{<1} = 0;
- equivariance: for σ a transposition of two roots of one bottom block, block(d, σc) at σx equals block(d, c) at x after transporting flows by σ.

**A10. Polynomiality certificate.**
- Passes for every template of A4, A5, A7, A8 and A9 (default ones) and for the closed forms Eps, Eta, Iota, Trace with dots (as WinLoc functions).
- Fails for:
  - (a) Eta multiplied, at one coset only, by 1/(x_0 − x_1) (non-polynomial);
  - (b) Rickard d_0 of (2,1) rescaled at every point so that its first nonzero entry is 1 (point-dependent normalization, C1 I7);
  - (c) Eps with a sign flipped on one non-identity flow.

  The error names the entry.
- Report T and max |Mono_D| for (2,2) and (3,2).

**A11. Reduced coordinates.** Exact idempotent on BS(s,s) = [E_0(1), F_0(1), E_0(1), F_0(1)] over (1,1): ι spans Hom^{−1}(B_s, BS) (1-dimensional by hom_gd); π is the first element of the exact basis of Hom^{+1}(BS, B_s) with G = π∘ι a nonzero scalar (one exists because B_s{−1} is a summand); e = ι π / G.
- `reduced_of(e_c)` gives Q P = I and P Q = E_c at every coset.
- rank E_c = m_c(B_s) from `schur` (1 at both strata).

Repeat on (2,2) [E(1), F(1), E(1), F(1)] with the summand B_{p1}{±1} (an ignored test if hom_space there exceeds 30 s; End^0 took 19 s in D1 §8.2).

**A12. Point independence.** Two `Eval`s at different points: every template block equals `localize_map` of the same exact PMat at that point. A Lin of templates with constant coefficients commutes with evaluation at 3 points.

**A13. Cache round trip.**
- `TmplLib::save`, then a fresh lib loads identical PMats (subtract-and-is_zero) for every template of A4-A9.
- The header contains version and field.
- A corrupted header is rejected.

**A14. Determinism.** Same seeds give byte-identical blocks (hash) across two runs in one process; FlowIndex orders are identical.

**A15. Arena checks.**
- `add` rejects, each with an Err and a message:
  - Comp with mismatched expressions;
  - Lin with mixed degrees or non-constant use;
  - Local whose window does not match the src/tgt words.
- `assert_genuine` (M4b) and the `certified` flag reject an `Elem` whose (kind, labels) has no recorded constant. The Err names kind and labels; test it with (2,2) Eps before `certify_elem`.
- Scratch regions: mark, add 50 random nodes over 10 new expressions, then `commit` keeping 3 of them. Required:
  - exactly the nodes and expressions reachable from the 3 survive (counts checked);
  - ids below the mark are unchanged;
  - the 3 renumbered maps evaluate to the same blocks at 2 points as before the commit;
  - `bytes()` decreases.
- `drop_index` followed by `index` rebuilds an identical FlowIndex (hash of `masks`, `start` and `cosets`).

**A16. Storage and thread safety.**
- Compile time: `fn _ss<T: Send + Sync>() {}` called on `Arena<F31>`, `TmplLib<F31>`, `FlowIndex`, `Layout` and `LocMap<F31>`.
- `FlowIndex::bytes()` ≤ 8·nr·nflows + 12·|support| + 1 KiB for every A2 ladder. For the (2,2,2) X with 432 flows and 7 rungs it is ≤ 26 KB (the hash-map layout took 97 KB [V-crit c6]).
- `LocMap::bytes()` = 8·len + O(1). Two maps with equal ends from one `Eval` share their `Layout` (`Arc::ptr_eq`).
- `gemm_acc` equals the naive product on 100 random shapes up to 20 × 20 × 20, including alpha ≠ 1 and accumulation into a nonzero c.
- `WinLoc::dense` of every A5/A9 template equals the reordered `local::localize_map` output.

### 3.6 Performance targets (release, one thread)

- FlowIndex of a (2,2,2) X with 432 flows: < 2 ms and ≤ 26 KB.
- FlowIndex of a (2,2,2,2) X with 51 840 flows (ignored test): < 0.3 s and ≤ 4.6 MB plus the support.
- `Eval::all` of a Local map on a ladder with 432 flows: < 20 ms.
- Template localization at window fibre ≤ 36: < 1 ms.
- The default suite: < 60 s debug.

### 3.7 Done criteria and report

- `cargo test --test m4a` passes. Ignored tests are run one at a time and pass, with time and peak RSS reported.
- `cargo test` for the whole crate still passes, with no change to existing tests.
- The report contains:
  - the constants c of A7/A8 per (kind, labels) (A7: all 1 expected);
  - the certificate statistics of A10;
  - the FlowIndex and LocMap byte counts of A16 and §3.6;
  - the A2 maxima;
  - the timings of §3.6;
  - every mismatch with a minimal example;
  - the list of cloud jobs, if any.

---------------------------------------------------------------------------------------------------

## 4. M4b: Hom provider, certification, gate G1

### 4.1 Goal

Produce certified, genuine maps for four kinds of targets, without ever accepting a partial answer:
- full Hom^d bases;
- decomposition pairs;
- closure generators;
- End^0 spanning sets.

Also measure the Search provider against gate G1.

### 4.2 Interface (`src/homprov.rs`; binding)

```rust
use crate::smap::{Arena, ExprId, MapId, Point};
use crate::schur::{Schur, Tid};

#[derive(Clone, Debug)]
pub struct Budget { pub l_max: usize, pub c_max: usize, pub t_max_s: f64, pub m_max_bytes: u64, pub flows_max: u64, pub index_max_bytes: u64 }
impl Default for Budget { /* l_max 6, c_max 4096, t_max_s 600, m_max 1 GiB, flows_max 4 * max(flows(src), flows(tgt)) at call time, index_max 256 MiB */ }

#[derive(Clone, Debug)]
pub struct SearchStats { pub candidates: usize, pub evaluated: usize, pub secs: f64, pub tier_hits: [usize; 7], pub max_flows: u64, pub points: usize,
                         pub visited_words: usize, pub indexed_words: usize, pub peak_index_bytes: u64, pub peak_arena_bytes: u64,
                         pub exact_calls: usize, pub exact_secs: f64, pub detect: &'static str /* "witness" | "orbits" | "more-points" */ }
/// Why a provider call stopped without a certified answer. Exhausted = candidate stream ended with
/// budget left (completeness signal); Time/Memory = budget overrun (engineering); Refused = an
/// admission rule declined the call (Exact).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MissKind { Exhausted, Time, Memory, Refused }
#[derive(Clone, Debug)]
pub struct Missing { pub what: String, pub degree: i32, pub stratum: Option<String>, pub got: usize, pub want: usize, pub kind: MissKind }

pub struct Pairs<F> { pub iotas: Vec<MapId>, pub pis: Vec<MapId>, pub form: Vec<Vec<F>> }   // form[j][i] = top scalar of pis[j]∘iotas[i], rank = mult

pub trait HomProvider<F: crate::field::Field> {
    fn name(&self) -> &'static str;
    /// certified k-basis of Hom^d(src, tgt) (|S| = Hom formula, stacked rank, fresh residual)
    fn full(&mut self, ar: &mut Arena<F>, s: &mut Schur, src: ExprId, tgt: ExprId, d: i32, b: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>>;
    /// `mult` pairs (iota: L_q -> X degree s, pi: X -> L_q degree -s) whose intersection form at a
    /// top-stratum coset of q has rank `mult`; e_q applied on both sides when present
    fn pairs(&mut self, ar: &mut Arena<F>, s: &mut Schur, x: ExprId, q: Tid, lq: ExprId, eq: Option<MapId>, sdeg: i32, mult: usize, b: &Budget) -> Result<(Pairs<F>, SearchStats), Vec<Missing>>;
    /// n = c_diag(b_p)(1) maps R^lambda -> L_p (ins) or L_p -> R^lambda (outs) with degree multiset
    /// c_diag(b_p) and nonsingular identity-coset matrix at two points (C2 M2)
    fn closure(&mut self, ar: &mut Arena<F>, s: &mut Schur, lp: ExprId, p: Tid, ep: Option<MapId>, outs: bool, b: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>>;
    /// a spanning set of End^0(L) (certified by stacked rank = Hom formula)
    fn end0(&mut self, ar: &mut Arena<F>, s: &mut Schur, l: ExprId, b: &Budget) -> Result<(Vec<MapId>, SearchStats), Vec<Missing>>;
}
/// defaults: 1 600, 5 000, max_label 2, costs None, max_call_s 10.0 (§4.5)
pub struct Exact { pub max_fibre_product: u64, pub max_unknowns: u64, pub max_label: u8, pub costs: Option<CostTable>, pub max_call_s: f64 }
/// measured hom_space / rlambda_bases costs keyed by (sorted labels, fibre product rounded up to a
/// power of 2, degree) -> max seconds observed; built and saved (text) in M4b
pub struct CostTable { /* ... */ }
impl CostTable { pub fn load(path: &std::path::Path) -> std::io::Result<CostTable>; pub fn predict(&self, labels: &[u8], fibre_product: u64, deg: i32) -> Option<f64>; }
pub struct Search { /* tiers S0..S6 (a sub-range can be selected), TmplLib windows, path cache */ }
impl Search { pub fn tiers(lo: usize, hi: usize) -> Search; }
pub struct Chain<F>(pub Vec<Box<dyn HomProvider<F>>>);                 // first Ok wins; Missing lists merged
impl<F: crate::field::Field> Chain<F> {
    /// §4.5 order: max label <= 2: [Exact, Search S0-S6]; max label >= 3: [Search S0-S2, Exact (cost table only), Search S3-S6]
    pub fn for_labels(max_label: u8, exact: Exact) -> Chain<F>;
}

// certification primitives
pub fn hom_gens(s: &mut Schur, src_char: &crate::schur::Elt, tgt_char: &crate::schur::Elt) -> std::collections::BTreeMap<i32, i128>;
pub fn m_points(lam: &[u8], gens: &std::collections::BTreeMap<i32, i128>, d: i32) -> usize;
pub fn detecting_cosets<F: crate::field::Field>(ar: &Arena<F>, s: &mut Schur, src: ExprId, tgt: ExprId) -> Vec<u64>;
pub fn stacked_rank<F: crate::field::Field>(ar: &Arena<F>, ms: &[MapId], pts: &[Point<F>], cosets: &[u64]) -> usize;
pub fn certify_full<F: crate::field::Field>(ar: &Arena<F>, s: &mut Schur, src: ExprId, tgt: ExprId, d: i32, ms: &[MapId], seed: u64) -> Result<(), String>;
pub fn assert_genuine<F: crate::field::Field>(ar: &Arena<F>, m: MapId) -> Result<(), String>;
pub struct G1Pair { pub mu: Vec<u8>, pub p: Tid, pub rung: crate::web::Rung, pub targets: Vec<(Tid, i32, usize)>, pub identity: bool }
pub fn g1_targets(s: &mut Schur) -> (Vec<G1Pair>, Vec<Tid> /* closure types (lam, lam) */);
```

### 4.3 Certification algorithms (SYNTHESIS §1.4)

1. **Generator degrees** (`hom_gens`). gd Hom(src, tgt) = `schur::hom_gd(char src, char tgt, lo, hi)`, multiplied by Π_{i, j<=λ_i} (1 − v^{2j}). The result must be a finite polynomial with nonnegative coefficients: these are the R^λ-generator degrees d_j (asserted; lo = −4K², hi = 4K²).
2. **M** (C1 I1). M(src, tgt, d) = max over d_j <= d with d − d_j even of `emono_count(λ, (d − d_j)/2)`.
3. **Detecting set** (F11).
   - Let T = the tables t with c_t(char src)(1) · c_t(char tgt)(1) > 0. Let E_w = Σ_{t ∈ T} c_t(char src)(1) · c_t(char tgt)(1), the number of entries of the witness blocks per point.
   - **Witness set**, if (M+2)·E_w ≥ dim Hom^d: cosets `rep(t)` for t ∈ T (one witness per stratum; D3's smaller S_e only after U9 is settled), with M + 2 seeded points. This holds for every d ≤ 4 in (2,2), (3,3), (2,1,2), (2,2,2) and every d ≤ 10 in (3,3,3) [V-crit witness].
   - **Fallback A**, otherwise, if (M+2)·Σ_t N(t)·c_t c'_t ≥ dim: the full W_λ-orbits `cosets_of(t)`, with M + 2 points.
   - **Fallback B**, otherwise: the witnesses with ⌈dim/E_w⌉ + 2 points.
   - The choice is recorded in `SearchStats::detect`.
   - Stacked rank = |S| certifies independence at any detecting set. So a too-small set could only cause false negatives, and the bound check turns them into a deterministic fallback. `certify_full` asserts (rows per point)·points ≥ dim before ranking.
4. **certify_full(S).** Valid iff:
   - |S| = [v^d] hom_gd;
   - `stacked_rank` = |S| at the detecting set of item 3;
   - at 2 fresh points the rank of S together with the stacked rows is still |S| (residual: no rank drop);
   - `assert_genuine` holds for every map.
5. **1-dimensional spaces.** If [v^d] hom_gd = 1, a single map m certifies iff `Eval::is_zero(m)` is false at one point and again at a second point.
6. **Pairs.** form[j][i] = block(π_j ∘ ι_i, w_q)[0][0], where w_q = rep(q) in the block of q (B_q has one flow there).
   - The rank must equal `mult`.
   - The same matrix is required at a second point (entries are constants: a disagreement is an Err, "point-dependent map").
   - Each ι_i, π_j must be genuine.
   - Stacking is not needed: form entries are top-stratum scalars of single maps.
7. **Closure.**
   - The degree multiset of the chosen maps equals c_diag(b_p) (normalized; shift handling as in `xengine::closure_data`).
   - The n × n identity-coset value matrix (rows: identity-coset flows of L_p, composed with e_p when present) is nonsingular at point 1 and at point 2.
8. **End0.** certify_full on Hom^0(L, L).
9. **Genuineness** (C3 H3). `assert_genuine(m)` walks the DAG. Every leaf must be one of:
   - an exact template;
   - an `Elem` whose (kind, labels) has a recorded constant (`elem_certified`, A7/A8); otherwise Err naming kind and labels;
   - a `Dot`.

   Template-free composites used as moves (S3) pass the polynomiality certificate once per window kind before use.

### 4.4 The Search provider (SYNTHESIS §2.4, made concrete)

**Move library** (all genuine, cached):
- M1 dots `Dot` on any segment, s_ν with |ν| <= 4;
- M2 E2/E3 closed forms (insert or remove a merge-split pair or a digon);
- M3 E4 associativity;
- M4 rung biadjunction units and counits: the lowest-degree exact maps 1 ↔ [R, R^∨] on 2-upright windows (`hom_basis`, dim 1 asserted);
- M5 square-switch maps: exact `hom_basis` between 2-upright window words of length 2 ([E^(a), F^(b)] ↔ [F^(b'), E^(a')]) in degrees −4..4, and 3-upright windows of length <= 2 with window fibre product <= 1 600;
- M6 Rickard templates;
- M7 certified maps already in the Arena (earlier decompositions, rex isomorphisms found earlier), applied to matching windows.

**A move** rewrites a contiguous window of the current rung word (insertion = an empty window). Its degree is that of the template.

**Tiers**, tried in order (the first that reaches the target stops the search):

| tier | candidates | applies to |
|---|---|---|
| S0 | identity (X = L_q word, mult 1, s = 0); the `Id` of equal words | pairs, full |
| S1 | **extremes**: one nonzero map in a 1-dimensional space (F4). Found by path search (S3) restricted to that degree; certified by rule 5 | pairs (ι at s_lo, π at s_hi), closure (lowest in- and out-generator) |
| S2 | **dot multiples** of certified maps: ι_s = Dot_f ∘ ι_{s_lo} (dots on segments of X) and ι_{s_lo} ∘ Dot_g (segments of L_q); π_s likewise from π_{s_hi}; closure generators g_a = Dot_f ∘ g_0 on interior segments of L_p. Order: segment index, then partitions by size, then lexicographic | pairs, closure, full |
| S3 | **paths**: best-first search over rung words from src to tgt with moves M1-M7. Cost = (number of moves, |degree − target|). Intermediate words: length <= max(len src, len tgt) + 2, labels <= K, flows <= flows_max, depth <= l_max. Each arrival with total degree d' <= d, d − d' even, yields path ∘ dots of degree d − d' (S2 order) | all |
| S4 | **double-leaf shape**: for each table r with c_r(src) c_r(tgt) ≠ 0, in descending l_max: (path L_r → tgt) ∘ Dot ∘ (path src → L_r) with paths from S3 (cached per word pair) | full, end0 |
| S5 | **adjunction transfer** (F3): candidates for Hom(L_q, L_p R) from maps L_q ⊗ R^∨ → L_p composed with the unit 1 → R^∨ R (M4); certified maps of the reverse rung's decompositions (M7) are reused | pairs |
| S6 | **Exact window fallback**: id_P ⊗ φ ⊗ id_Q with φ ∈ `hom_basis` of the maximal non-common window of src and tgt, when its fibre product is <= 1 600 (P8: rarely applicable, kept for completeness) | all |

**Main loop.** Candidates are streamed lazily in tier order, deterministically (fixed seeds). Each candidate is checked by `assert_genuine` and appended:
- full / end0: to the stacked echelon at the detecting set; stop at rank = target;
- pairs: to the intersection-form matrix; stop at rank = mult;
- closure: to the greedy degree-by-degree determinant selection; stop at n.

Exhausted tiers produce `Missing { kind: Exhausted, .. }`; a budget overrun produces `Missing { kind: Time | Memory, .. }`. Either way the Err carries degree, stratum, got, want and the SearchStats. Partial bases are never accepted.

**Scratch region and memory.**
- Each provider call runs between `ar.mark()` and `ar.commit(mark, certified)`. The words visited by S3/S4 and the rejected candidates are dropped when the call returns. Only the certified results and their sub-DAGs stay in the Arena that M5 freezes.
- A FlowIndex is built (`ensure_indices`) only for words that a candidate is evaluated on, never for words merely visited.
- Inside a call, the indices of the scratch region form an LRU cache bounded by `Budget::index_max_bytes` (256 MiB). An evicted index is rebuilt identically on demand.
- With packed indices (8 B per flow and rung), 256 MiB holds about 55 (2,2,2,2) X-sized words (51 840 flows, 11 rungs, 4.6 MB) or about 140 (3,3,3) X-sized words (31 500 flows, 7 rungs, 1.8 MB).
- SearchStats records visited words, indexed words, peak index bytes and peak Arena bytes.
- Overrunning `m_max_bytes` ends the call with `MissKind::Memory`. This is an engineering failure and is kept apart from `Exhausted` in G1.

Completeness is not guaranteed [?] (U1). G1 measures it.

### 4.5 The Exact backend

- `full`/`pairs`/`closure`/`end0` through `sbim::hom_space` on the whole ladders, wrapped as whole-ladder templates.
- **Admission.** The fibre caps alone admit calls of 50-300 s once labels are ≥ 3. Example: the (3,2) one-step call 144 → 8 has fibre product 1 152 and took 275 s (M6.0 report). A call is admitted iff one of:
  - max label ≤ `max_label` (2), fibre product ≤ `max_fibre_product` (1 600), and unknowns = Σ_{β,β'} #monos(D) ≤ `max_unknowns` (5 000);
  - the `CostTable` predicts ≤ `max_call_s` (10 s) for (labels, fibre product, degree). The prediction is the maximum over the key's bucket, because cost is not a function of fibre: mirror-twin (2,2,2) types of fibre 36 took 0.30 s and 43.8 s in `rlambda_bases` [V-crit c5].

  Otherwise `Missing { kind: Refused }` is returned at once, with no call.
- **No in-process interruption.** `hom_space` cannot be interrupted inside the process, so admission is the only in-process guard.
  - The G1 harness runs one process per pair under the external watchdog. A kill is recorded as `MissKind::Time` for the running target.
  - `sscan` runs under the same watchdog in tests and scripts.
- **Chain order** (`Chain::for_labels`):
  - max label ≤ 2: Exact, then Search S0-S6;
  - max label ≥ 3: Search S0-S2, then Exact (cost table only), then Search S3-S6. The extremes and dot multiples (F4) are tried before any exact call.
- **Calibration** [V-P3, P12; M6.0 report]:
  - 6·36 with 4 084 unknowns: 0.1 s;
  - 20·60 with 21 694 unknowns: > 150 s;
  - (3,3) 15·36: up to 48 s per decomposition;
  - (3,2) rungwise: at most 9.7 s at fibre 36;
  - (3,2) one-step 144·8: 275 s, and 8·144: 68 s.

  M4b measures at least 20 calls per label class and saves the `CostTable` (unknowns, fibres, labels, degree, time, RSS) next to the template cache.
- **Closure.** `cclose::rlambda_bases` is admitted only for fibres ≤ 36, and only through the same rule (max label ≤ 2 or cost table ≤ 10 s) [V-P9: up to 40 s there].

### 4.6 Tests (`tests/homprov.rs`)

- **B1.** `hom_gens` equals the generator degrees found by exact Nakayama (`cclose::rlambda_bases` degrees) on the closed canonical ladders of (2,2), (3,3) and on the (2,2,2) types of fibre <= 18.
- **B2.** Certified dimension = Hom formula: on the 60 R1 cases (`tests/data/schur/hom60.txt`) through the Exact backend, and on all type pairs of (2,2), (3,3), (2,1,2), (2,2,2) in degrees 0..4 through the Chain (Exact, Search).
  - For type pairs between indecomposables, M <= dim R^λ in weight 2: 9 for (2,2,2) and 5 for (3,3).
  - Expected: all certified; the counts per composition are reported.
- **B3.** C1 faith table: for the 9 cases of C1 I1, the number of points at which `stacked_rank` reaches the full dimension equals the predicted M: 2, 2, 5, 2, 2, 5, 3, 9, 3 (the (2,2) E1F1E1F1 End^0 case takes about 20 s: an ignored test).
- **B4.** False-negative guard. For each B3 case, `stacked_rank` at a single point is below the dimension exactly as C1 I1 measured. B4 runs both coset sets:
  - full cosets (`cosets_of` orbits): 4/5, 4/5, 7/14, 2/3, 3/4, 36/48, 3/5, 7/21, 7/9;
  - witness cosets (the detecting set of §4.3 item 3): 4/5, 4/5, **5/14**, 2/3, 3/4, **35/48**, 3/5, **5/21**, 7/9.

  The points needed for full rank are the same for both sets: 2, 2, 5, 2, 2, 5, 3, 9, 3 (C1 I1). This documents why one point is forbidden.
- **B4b.** Detecting-set fallback. On (2,2) types (1,2) at d = 11 (dim 69, M 20, E_w 3, bound violated [V-crit]), `certify_full` with an exact basis selects Fallback A or B (recorded) and certifies. Forcing the witness set with M + 2 points fails the counting assertion, not the rank test. Ignored if hom_space at d = 11 exceeds 60 s.
- **B5.** Genuineness. An Arena map with a template replaced by its point-rescaled copy (A10 b) is rejected by `assert_genuine` (it is not in the certified library), and its polynomiality certificate fails.
- **B6.** Pairs on (2,2) and (2,1) decompositions of M6.0: the intersection-form ranks equal `term_mult`. ι, π from Search equal the exact M6.0 maps up to the Π∘Ι-normalized basis change: `localize` of the exact xengine-style maps lies in the span (stacked rank test).
- **B7.** Determinism. Two runs give identical candidate sequences, MapIds and SearchStats counters.
- **B8.** Admission and order.
  - The (3,2) one-step request 144 → 8 in degree 2 is refused without a call (`MissKind::Refused`, `exact_calls = 0`) when no cost-table entry admits it.
  - With `Chain::for_labels(3, ..)` on a (3,3) pair whose extremes S1/S2 certify, `SearchStats::exact_calls = 0`.
  - After a provider call, `Arena::bytes()` equals its value before the call plus the bytes of the certified sub-DAG (scratch region released).

### 4.7 Gate G1 (concrete)

**Harness:** `src/bin/g1.rs` (`g1 [a|b] [--pairs N] [--seed S] [--json]`). One process per pair under a watchdog (1 GB), sequential, nice 19.

**G1(a), needed for M6a.** Targets from `g1_targets`:
- every non-identity (type, rung) pair of the blocks reached by one crossing term from λ ∈ {(3,3), (2,2,2)}, each with all its summand targets;
- `closure` (ins and outs) for every type of (λ, λ).

| λ | pairs (non-identity) | (q,s) targets | Σ mult | closure types / Σ n (ins = outs) |
|---|---|---|---|---|
| (3,3) | 12 | 72 | 121 | 4 / 8 |
| (2,2,2) | 98 | 358 | 404 | 21 / 48 |

**G1(b), scale-out.** 50 non-identity pairs of (2,2,2,2) and 20 of (3,3,3), drawn uniformly with `StdRng::seed_from_u64(2026)` from the `g1_targets` populations (2 191 and 429 pairs). Plus closure for 10 random (2,2,2,2) and 5 random (3,3,3) types.

**Budget per pair** (all its targets): T_max = 10 min wall, 1 GB RSS, `Budget::default()`.

**Pass / fail:**
- **G1(a) passes** iff 100 % of its targets are certified within budget. Total wall time is reported. Predicted <= 30 min; run on the laptop when the measured mean per pair is < 15 s, otherwise as a cloud job.
- **G1(b) passes** iff ≥ 95 % of the sampled pairs certify within budget, and no pair ends with `MissKind::Exhausted` after the rerun below. G1(b) always runs on the cloud (70 pairs, up to 11.7 h, plus reruns).
- **Rerun rule.**
  - Whether the candidate stream ends depends on the Budget constants. So a pair that stops with `Exhausted` while budget remains is rerun automatically with l_max + 2 and 4·c_max (same seeds), and both runs are reported.
  - Only an `Exhausted` outcome of the larger run counts as a completeness failure (red flag for U1, and input to the M4c decision).
  - Pairs that end with `Time` or `Memory` (including watchdog kills) count as engineering failures. They count against the 95 %, but they do not trigger M4c by themselves.
  - The report separates the two classes and gives peak index bytes, peak Arena bytes and visited words per pair.

**Recorded per target:** tiers used, candidates tried against the target rank, time, max intermediate flows, points. Also the fraction of extremes found by S1 and of non-extreme degrees covered by S2 (tests the "extremes plus dots" hypothesis).

**On failure:**
- **G1(a) fails:**
  1. M6a still runs its 2-strand part with `Chain(Exact)`: labels (1,1), (2,1), (1,2), (2,2) by default, and (3,1) through the cost table (its decompositions take about 0.01 s); (3,2) as an ignored test.
  2. Only the (3,3) braids that need no Hom space (kinks, by the identity shortcut) stay in scope.
  3. M4c (EKLP singular light leaves, SYNTHESIS §3.1) starts before the 3-strand part of M6a. The 3-strand part of M6a re-runs G1(a) with `Chain(Exact, Search, Leaves)`.
  4. The failing targets are listed in the report as the M4c acceptance set.
- **G1(b) fails, G1(a) passes:** M6a proceeds. M4c starts in parallel and must pass on the G1(b) set before any 4-strand or 3-colored 3-strand work. G2 at the end of M6a re-decides.

### 4.8 Done criteria

- `tests/homprov.rs` passes (B3 and the E1F1E1F1 case as ignored tests, run once).
- G1(a) has run and its JSON is attached. G1(b) is queued as a cloud job with exact commands, commit and seeds.
- The report contains:
  - the G1 tables;
  - the Exact calibration;
  - the tier statistics;
  - a recommendation: M6a 3-strand go, or M4c first.

---------------------------------------------------------------------------------------------------

## 5. M5: singular idempotents and merge decompositions (`src/sindec.rs`)

### 5.1 Interface (binding)

```rust
use crate::homprov::{Budget, HomProvider};
use crate::schur::{Schur, Tid};
use crate::smap::{Arena, Eval, ExprId, MapId, Point, Reduced};
use crate::web::Rung;
use std::sync::Arc;      // no Rc anywhere in sindec: DecData must be Send + Sync (engine2::par_map requires Sync)

pub struct SType { pub p: Tid, pub lad: ExprId, pub e: Option<MapId>, pub rank: u128 }
#[derive(Clone)] pub struct SSummand { pub q: Tid, pub s: i32, pub iota: MapId, pub pi: MapId }
pub struct Decomp { pub p: Tid, pub rungs: Vec<Rung>, pub x: ExprId, pub summands: Vec<SSummand>, pub identity: bool }
pub struct SIndecOpts { pub cache_dir: Option<std::path::PathBuf>, pub budget: Budget, pub check_points: usize /* 2 */, pub seed: u64 }

/// Builder: owns the Schur instance, the Arena and the provider; decompositions are computed
/// lazily and cached (memory and disk).
pub struct SIndec<F: crate::field::Field> { pub schur: Schur, pub ar: Arena<F>, /* provider, caches, stats */ }
impl<F: crate::field::Field> SIndec<F> {
    pub fn new(lam: &[u8], provider: Box<dyn HomProvider<F>>, opts: SIndecOpts) -> Self;
    pub fn ty(&mut self, p: Tid) -> Result<&SType, String>;                      // canonical ladder; e_p if needed
    pub fn dec_rung(&mut self, p: Tid, r: Rung) -> Result<Arc<Decomp>, String>;
    pub fn dec_word(&mut self, p: Tid, rs: &[Rung]) -> Result<Arc<Decomp>, String>; // rung by rung; multiplicities = term_mult(p, rs)
    pub fn rickard_comp(&mut self, p: Tid, c: u8, d: u8, positive: bool, j: usize) -> Result<Arc<Vec<(usize, usize, MapId)>>, String>; // π_τ ∘ (id ⊗ d_j) ∘ ι_σ
    pub fn freeze(self) -> Arc<DecData<F>>;                                        // immutable snapshot for parallel scans
    pub fn stats(&self) -> DecStats;
}
/// Frozen data: Arena (immutable), types, decompositions, Rickard components; read-only lookups
/// that return Err("not precomputed") instead of computing. Send + Sync (D9); no interior mutability.
pub struct DecData<F: crate::field::Field> { /* ... (Arc<Decomp>, Arc<Vec<..>> inside) */ }
impl<F: crate::field::Field> DecData<F> {
    pub fn ar(&self) -> &Arena<F>;
    pub fn ty(&self, p: Tid) -> Result<&SType, String>;
    pub fn dec_word(&self, p: Tid, rs: &[Rung]) -> Result<Arc<Decomp>, String>;
    pub fn rickard_comp(&self, p: Tid, c: u8, d: u8, positive: bool, j: usize) -> Result<Arc<Vec<(usize, usize, MapId)>>, String>;
    /// computed on each call from the per-thread Eval (identity when e_p = None); callers memoize per point
    pub fn reduced(&self, ev: &mut Eval<F>, p: Tid, c: u64) -> Result<Reduced<F>, String>;
}
pub struct DecStats { pub decomps: usize, pub identity: usize, pub provider_calls: usize, pub provider_secs: f64, pub max_x_flows: u64, pub idempotents: usize }
```

### 5.2 Algorithms

**Types.** L_p = `canonical_basis(μ).ladder_minflow[p]`, interned as an Expr. If `ladder_char(L_p) == b_p`, then e = None [F1]. Otherwise the idempotent algorithm below applies.

**Idempotent e_p = P(a)** (C1 I1(a); SYNTHESIS §1.4):
1. E = provider.end0(L_p).
2. a = Σ_i c_i E_i with c_i = seeded constants in 1..10^6, scaled so that the top-stratum scalar of a at rep(p) is 1 (a constant: the top scalar of each E_i is constant).
3. At point 1, compute the minimal polynomial of the stacked block-diagonal matrix of a over the cosets rep(t) of all tables t in the support (Krylov as `indec::primitive_idempotent`).
4. Split off the (t−1)-primary component with `indec::pdivrem`/`pxgcd`; this gives the polynomial P.
5. e = P(a) as a Horner chain of `Comp`/`Lin` nodes.
6. Checks:
   - e² = e at points 2 and 3, all cosets (Err otherwise);
   - rank block(e, rep(t)) = c_t(b_p)(1) for every table t;
   - top scalar 1.
   - With ranks matching b_p at every stratum and top-stratum rank 1, Krull-Schmidt gives Im e ≅ B_p (proof in §7 item 2).
7. Never re-solve e in a basis.

**One-rung decomposition `dec_rung(p, R)`**:
1. X = L_p · R (Expr); mult = `term_mult(p, [R])`.
2. **Identity shortcut** (as xengine, with the same guard): if mult is a single (q, 1) and L_q's word = X's word and e_p = e_q = None, the result is ι = π = Id.
3. For each (q, s, cnt) in mult, ordered by (l_max(q) descending, s ascending): (I, P, G) = provider.pairs(X, q, L_q, e_q, s, cnt) with source idempotent (e_p ⊗ id_R) applied: ι ↦ Kron(e_p, R)∘ι, π ↦ π∘Kron(e_p, R).
   - Choose cnt independent rows J and columns I of G (rref).
   - π' = G[J][I]^{-1} π_J as a `Lin` with constant coefficients.
   - Re-check G at a second point: equal.
4. **Orthogonalize** (xengine Neumann series, now with global maps):
   - N_ts = π_t ∘ ι_s for t ≠ s; zero-tested at 2 points;
   - unitriangularity asserted: N_ts : B_{q_s} → B_{q_t} has degree s_s − s_t, so it vanishes unless s_s > s_t (degree 0 between different summands is excluded by Hom^0(B_q, B_q') = δ k and the normalization of each (q, s) block). This is the assertion of `xengine::decompose_direct`;
   - Π' = Σ_k (−N)^k Π as `Lin`/`Comp` nodes; the iteration count is at most the number of summands.
5. **Checks** at check_points points, all cosets of X's support:
   - π_t∘ι_s = δ_ts · (e_q or id);
   - Σ_t ι_t π_t = e_X;
   - Σ_t rank(B_{q_t}) = rank(B_p ⊗ R) (= flows(X) when all e = id);
   - every component homogeneous: degree(ι) = s, degree(π) = −s.

**Composite `dec_word(p, rs)`** (rung by rung, xengine): ι = Kron(ι_1, rs[1..]) ∘ ι_2 and π = π_2 ∘ Kron(π_1, rs[1..]), recursively. The assertion: the multiset of (q, total s) equals `term_mult(p, rs)`. Π∘Ι = id is re-checked at 2 points on all cosets (cheap).

**Rickard components `rickard_comp(p, c, d, sign, j)`.**
- Source term ks = j and target term kt = j+1 (positive), reversed for negative.
- φ = Local(template d_j at the top window, rungs of L_p as prefix).
- The result: the list of (σ, τ, Comp(π_τ, Comp(φ, ι_σ))) whose evaluation at 2 points is nonzero (zero components dropped; the zero test is a single-map statement).

**Disk cache.** One text file per (field, λ, p's table, rung word):
- the Decomp as the list of Arena nodes reachable from its ι/π (templates by content hash, constants as integers);
- version header "krhom-sindec v1".

Loading re-runs the step-5 checks at 1 point.

### 5.3 Tests (`tests/sindec.rs`)

- **D1.** For every type of (2,2), (3,3), (2,1,2) and (2,2,2) and every block they reach (63 types for (2,2,2)):
  - e_p = None (F1);
  - rank of L_p at one coset per stratum = c_t(b_p)(1);
  - End^0 dimension 1 (hom_gd).
- **D2.** Idempotents:
  - (3,3,3) type 0,1,2/1,1,1/2,1,0 (End^0(L_p) of dim 2): an ignored test if `end0` exceeds 10 min;
  - the 6 decomposable (2,2,2,2) types: cloud;
  - BS(s,s) over (1,1) and [E(1),F(1),E(1),F(1)] over (2,2) as forced non-canonical ladders.

  Requirements: e² = e at 3 fresh points; block ranks = KL values at 1; top scalar 1; rank(Im e) = rank(B_p).
- **D3.** `dec_rung` for all 98 non-identity (2,2,2) pairs, the 12 (3,3) pairs, (2,2), (2,1)/(1,2) and (2,1,2): every check of §5.2 step 5 holds. Multiplicities = `term_mult`. Default or ignored according to G1(a)'s measured time (default if the total is < 60 s).
- **D4.** `dec_word` for all crossing terms W_k of (2,2), (2,1), (1,2), (3,1), (3,2), (3,3) and (2,2,2): multiset = `term_mult(p, W_k)`, Π∘Ι = id.
- **D5.** Comparison with M6.0 on (2,2) and (2,1). For every (p, k), the summand list (q, s) equals xengine's, and the transported Rickard components of `rickard_comp` match xengine's `dloc` components up to the basis change: for every (σ, τ), the localized exact xengine component and ours have the same rank at every coset, and the transported complex (§6) is identical up to isomorphism (checked in M6a).
- **D6.** (1,1,1,1): for all 72 (x, s) with x ∈ S_4 and s ∈ {1,2,3} (SYNTHESIS wrote 84; 24 × 3 = 72), the multiplicities of `dec_word(x, [E_s(1), F_s(1)])` equal `indec::kl_mult` (`hecke.rs` ↔ permutation tables as in M3a T-e). The decompositions themselves are certified.
- **D7.** Disk cache round trip: a reloaded Decomp evaluates to identical blocks at 2 points.
- **D8.** Determinism and freeze: `freeze()` followed by `DecData::dec_word` returns the same MapIds; a request not precomputed gives Err.
- **D9.** Thread safety, as a compile-time assertion: `fn _ss<T: Send + Sync>() {}` called on `DecData<F31>` and `Arc<DecData<F31>>`. Also a runtime check: `std::thread::scope` with 2 threads, each with its own `Eval` over a shared `Arc<DecData<F31>>` and its own point, gives the same blocks as sequential evaluation. `engine2::par_map` is made pub only in M6a.

Budgets: default suite < 60 s; (2,2,2) all-pairs and (3,3,3) are ignored or cloud depending on G1(a) timings.

### 5.4 Done criteria

- Tests pass.
- The report contains:
  - DecStats per composition;
  - provider time per pair (distribution);
  - the cache size on disk;
  - the idempotent results (End^0 dims, Krylov degrees);
  - any multiplicity mismatch (blocking).

---------------------------------------------------------------------------------------------------

## 6. M6a: R-full localized engine on <= 3 strands (`src/sengine.rs`, `src/bin/sscan.rs`)

### 6.1 Scope and phases

- **M6a.1 (2 strands).** Any provider, `Chain(Exact)` suffices for labels <= 2 and (3,1). It must reproduce M6.0 exactly. It does not depend on G1.
- **M6a.2 (3 strands and (3,3)).** Needs G1(a) passed, or M4c, and T0 passed for every (2,2,2) and (3,3) type. Targets and where they run (F6, §6.6):
  - laptop:
    - 3-colored 2-strand: Hopf, T(2,3), T(2,5);
    - uncolored 3-strand regression against engine2;
    - the 2-colored 3-strand Markov family (T9);
    - 2-colored 4_1;
  - cloud jobs under M6a (scan plus 7-103 min of rescans per prime and point set): 2-colored 5_2, 6_2, 6_3 and 8_19;
  - not in M6a: 2-colored 3-strand words with ≥ 9 crossings. They need M6c (coordinates of the final components) or another F path without rescans, and remain cloud jobs even then: M6c needs emono((e − d_j)/2) points per component, 119-630 for component degrees 12-20.

### 6.2 Interface (binding)

```rust
use crate::cclose::ClosureDataC;
use crate::crossing::Gen;
use crate::field::Field;
use crate::homprov::HomProvider;
use crate::laurent::Poincare;
use std::collections::BTreeMap;

#[derive(Clone, Debug)]
pub struct SOpts {
    pub seed: u64,             // points and generator choices
    pub check_d2: bool,        // d∘d = 0 at the point after every crossing
    pub streaming: Option<bool>, // streaming elimination (§6.4); None (default): true iff the decategorified peak Σ S > 1e5
    pub points0: usize,        // initial scan points, default 24 (engine2)
    pub time_budget_s: f64,    // default 600 (env KRHOM_TIME_BUDGET_S); refuse to extend beyond it
    pub interp_cap: usize,     // dense interpolation cap, default 2500 unknowns
    pub residual: usize,       // residual points for Gram and F interpolation, default 8
    pub mem_point_bytes: u64,  // per-point budget: Eval memo + complex, default 512 MiB (asserted after every crossing)
    pub f_mem_bytes: u64,      // stored F samples, default 512 MiB (asserted before the rescans start)
    pub verbose: bool,
}
pub struct StageInfo { pub step: usize, pub labels: Vec<usize>, pub unsimplified: usize, pub objects: usize, pub components: usize,
                       pub lower_bound: u128, pub stored_entries: u64, pub sum_s: u128, pub secs: f64 }
pub struct Measured { pub c_e: f64, pub lb_ratio: f64, pub final_degree_hist: BTreeMap<i32, usize>, pub peak_rss_mb: f64,
                      pub scan_secs_per_point: f64, pub points: (usize, usize), pub weights: (u32, u32) /* max Gram, max F */,
                      pub zero_block_frac: f64, pub peak_point_bytes: u64, pub f_sample_bytes: u64, pub interp_secs: f64, pub streaming: bool }
pub struct SScan<F: Field> {
    pub labels: Vec<usize>, pub word: Vec<Gen>,
    pub data: ClosureDataC<F>,                 // contract v2, object q = t, types = schur Tids
    pub stages: Vec<StageInfo>,
    pub final_objects: Vec<(String, i32, i32)>, // (table, h, t), sorted, as XScan::final_objects
    pub k0: Vec<(String, i32, i128)>,
    pub measured: Measured,
    pub timings: (f64, f64, f64),              // (scans, closure data, total)
}
pub fn sscan<F: Field>(labels: &[usize], word: &[Gen], provider: Box<dyn HomProvider<F>>, opts: &SOpts) -> Result<SScan<F>, String>;
#[derive(Clone, Copy, Debug)] pub enum CloseMethod { Exact, DeltaTower, Auto }
pub struct ClosePlan { pub n: usize, pub reduced: bool, pub exact_bytes: u128, pub tower_bytes: u128, pub cap_bytes: u128, pub method: Option<CloseMethod> /* None = refused */ }
impl<F: Field> SScan<F> {
    /// estimates from cclose::exact_mem_estimate and delta_tower_estimate on this data, before closing
    pub fn close_plan(&self, ns: &[usize], reduced: bool) -> Vec<ClosePlan>;
    /// Auto: exact oracle if its estimate <= cap; else the (Delta G) tower if its estimate <= cap and (a, N)
    /// is G4-certified; else Err naming N, reduced/unreduced, both estimates and the cap (no bare Err)
    pub fn close(&self, ns: &[usize], reduced: bool, m: CloseMethod) -> Result<BTreeMap<usize, Poincare>, String>;
}
pub fn khr_singular<F: Field>(labels: &[usize], word: &[Gen], ns: &[usize], reduced: bool) -> Result<BTreeMap<usize, Poincare>, String>;
```

CLI: `sscan LABELS WORD NS [--reduced] [--prime 31|32b] [--streaming] [--budget SECONDS] [--provider exact|search|chain] [--json]`. It prints the stage statistics, Measured and the Poincaré polynomials, as `xscan` does.

### 6.3 One scan at one point (`PEngine`)

**State.** `SCx<F> { objs: BTreeMap<usize, (Tid, h, t)>, out: BTreeMap<usize, BTreeMap<usize, LocMap<F>>>, inn: BTreeMap<usize, BTreeSet<usize>> }`.
- LocMap blocks are in reduced coordinates (= flow coordinates when e = None) at the cosets of the current top block.
- Each component is one `LocMap` over a `Layout` shared through `Arc`. Layouts are keyed by (src type, tgt type, top block) in a per-scan `HashMap`. Rows and columns are m_c (reduced coordinates), and the cosets are the joint support.
- All block arithmetic (sandwiches, elimination updates) runs on slices with `gemm_acc`, with no allocation per block.
- One `Eval` over the DecData's Arena at the point, per thread. Reduced coordinates are memoized per point in the PEngine.

**Start.** One object (identity type of (λ, λ), 0, 0).

**Crossing g at (i, i+1)** with current labels (c, d) and sign ε: terms W_k = `crossing_terms(c, d, ε, i)`.
1. **Objects.** For every object o = (p, h, t) and term k: D = dec_word(p, W_k.rungs). Its summands become objects (q_j, h + hdeg_k, t + dq_k + s_j); place[(o, k)] = their ids.
2. **d_C ⊗ id.** For every edge o → o' with LocMap f and every k:
   - f ⊗ id_{W_k} at a new coset c'' is block-diagonal over the source flows (σ, w): entry ((τ, w), (σ, w)) = f_{c'}[τ][σ] (§3.4.2 `Kron` rule applied to numbers);
   - the component for summands (σ_j, τ_l) = π_{τ_l} ∘ (f ⊗ id) ∘ ι_{σ_j}, computed blockwise per coset into a fresh `LocMap` over the layout of (type σ_j, type τ_l). The ι/π blocks come from Eval, sandwiched by `reduced` when e ≠ None;
   - zero results (every block zero) are dropped.
3. **id ⊗ d.** For every object o and j: the components of `rickard_comp(p, c, d, ε, j)` evaluated at the point, times (−1)^h; added between place[(o, ks)] and place[(o, kt)].
4. **Checks.**
   - Degrees: every component has total q-degree 0. This is asserted through the MapInfo degrees, which add under composition: **e(π_τ) + e(f) + e(ι_σ) = t_src − t_tgt**, with e(ι_σ) = s_σ and e(π_τ) = −s_τ. The placement of step 1 gives t_src − t_tgt = e(f) + s_σ − s_τ.
   - Components raise h by 1.
   - Optional d∘d = 0 at all cosets.
   - Memory: `Eval::memo_bytes()` + Σ `LocMap::bytes()` + the layouts ≤ `mem_point_bytes` after every crossing. Otherwise Err("per-point memory ...: run with streaming or on the cloud").
5. **Elimination** (§6.4).
6. **K_0 and LB.**
   - The K_0 class of the minimal complex equals the decategorified scan (`schur` `mul_rungs`/`expand`, exactly `xengine::decat_step`/`decat_k0`);
   - #objects >= LB.
   - A violation returns Err.

**Determinism record.** After each crossing, a hash of (sorted object list, pivot sequence) is appended to `trace`. Scans at other points must reproduce the trace exactly; otherwise Err("point-dependent scan").

### 6.4 Elimination

- **Pivot candidates.** Components a → b with type(a) = type(b), t(a) = t(b) (h(b) = h(a) + 1 holds automatically).
- **Pivot scalar.** The 1 × 1 block of the component at w = rep(type) (top stratum, one reduced coordinate). Assert:
  - the block at every coset is c·I (End^0 = k; one check per pivot, cost Σ_c m_c²);
  - c ≠ 0.
- **Order: Markowitz.** Cost (|in(b)| − 1)(|out(a)| − 1), ties broken by (cost, a, b). Kept in a lazy min-heap with versioned entries: after an update, the affected pairs are re-pushed; stale entries are skipped. This replaces xengine's quadratic rescans (M6.0 review). Equivalence test against the rescanning order on random complexes.
- **Update.** For x ∈ in(b) ∖ {a}, y ∈ out(a) ∖ {b}: d_xy −= d_ay ∘ c^{-1} ∘ d_xb, blockwise per coset. Remove a and b.
- **Streaming rule** (C1 I5; SYNTHESIS §1.4). With streaming on, terms are tensored one at a time. Streaming is on by default when the decategorified peak Σ S exceeds 10^5.
  - Positive crossing: with terms 0..k+1 present, a pivot a → b is eligible only if a lies in a term ≤ k.
  - Negative crossing (the differential runs from term k+1 to k): a pivot is eligible only if b lies in a term ≤ k.
  - After the last term, all pivots are eligible.
  - Why this suffices (re-derived by critique 1): Gaussian elimination reads only in(b) and out(a), since d'_xy = d_xy − d_ay c^{-1} d_xb. For a positive crossing in(b) is always complete, and out(a) is complete iff a's term is ≤ k. The negative case is the mirror.
- **Unit test (part of T5)**, term by term. A broken eligibility rule only shows when later terms add arrows at objects that were already eliminated. So the test must build C ⊗ W incrementally.
  - Generator (200 seeds):
    - C: a random complex over (2,2) types with 3-8 objects and scalar and block components. d² = 0 by construction: a direct sum of cones of identity maps and a minimal part, conjugated by a random invertible unitriangular change of basis.
    - W: a random term complex W_0 → … → W_m (m = 1..3) with a genuine differential d_W: either the localized (2,2) Rickard components of xengine at one point, or random unitriangular cones.
    - C ⊗ W gets arrows from both d_C ⊗ id and id ⊗ d_W (Koszul sign).
  - Streaming inserts C ⊗ W_k one term at a time (ascending k for a positive crossing, descending for a negative one). After each insertion it eliminates every eligible pivot. Full-then-eliminate builds C ⊗ W first.
  - Required: identical K_0 classes, identical object multisets, and isomorphic complexes (equal ranks of every (h, type, t) block of the residual differential, and d² = 0).
  - Sensitivity: the variant with eligibility ignored must differ from full-then-eliminate (or violate d² = 0) on at least one seed. If it never differs on 200 seeds, the generator is too weak and the test fails.
- **Streaming equality on real scans** (T1, T2): streaming on and off give identical final objects, K_0 and Poincaré polynomials on every T1/T2 word.

### 6.5 Multi-point closure data (engine2 style; contract v2)

1. **Points.**
   - Scan at P[0] sequentially with the builder `SIndec`, which computes all decompositions, then `freeze()`.
   - Scan the remaining points of P_0 = P[0..points0] in parallel (`engine2::par_map`, KRHOM_THREADS) against the frozen `Arc<DecData>`.
   - All final object lists and traces must be identical.
   - The full scan point sequence P = `points(K, seed, n)` is one deterministic sequence; step 4 fixes its length before any rescan.
2. **Closure generators** per final type p (final tops = λ; otherwise Err "not closable"):
   - ins = provider.closure(L_p, ins) and outs likewise, called during the builder phase;
   - their degree multisets must equal c_diag(b_p) (as in xengine);
   - the determinant criterion is re-verified at every scan point (C2 M2).
3. **Gram.**
   - G[b][a](x) = block(h_b ∘ g_a, identity coset) (1 × 1), evaluated at generator-only points (no scan).
   - Weight w = (d_a + d_b)/2. An entry with w < 0 or w not an integer must be 0, and this is asserted at every point.
   - Interpolate in the e-monomials of R^λ of weight w (`cclose::ERing::monos_of_weight`) with a dense Vandermonde over the first m = |monos| points. The residual check uses the next `residual` (8) points.
   - One inverse per weight, applied to all entries of that weight as one blocked product (as `engine2::Interp`).
   - Refuse above `interp_cap` unknowns (Err "needs CL2").
4. **F** (critique 2: sample handling is binding).
   - w_F = max over final edges of (e + d_a − d_c)/2 (degree bound, as engine2's `maxk`). Entries whose weight is negative or not an integer must be 0 at every scan point (asserted).
   - If emono(λ, w_F) > `interp_cap`: Err "needs CL2". For (2,2,2) this happens at w_F ≥ 15 (2 838 unknowns).
   - **Point set.** n_F = emono(λ, w_F) + `residual` scan points: P[0..n_F], fixed before any rescan.
   - **Time rule.** Measure t_scan on P_0. If (n_F − points0) · t_scan / threads > time_budget_s, return Err("cloud job: predicted … s for … rescans; command: sscan …").
   - **Memory rule.** Entry e with weight w_e keeps its values at P[0..emono(w_e) + residual] only. The store therefore needs 8 · Σ_e (emono(w_e) + residual) bytes. This is computed and asserted ≤ `f_mem_bytes` before the rescans start; otherwise Err("cloud job: F samples need … MB"). Example upper bounds: 8_19, about 1.2·10^4 entries, ≤ 0.21 GB; a 10-crossing word, ≥ 2.7·10^5 entries, ≤ 4.7 GB.
   - **Rescans** run in parallel. For each entry still needing that point, a per-point scan returns F(x) = solve(G'_id(x), f_id(x) G_id(x)), where f_id is the identity-coset block of the final component (exactly `engine2` lines 500-520, with identity-coset blocks of L_p). The scan's complex and Eval are dropped as soon as these values are extracted.
   - **Interpolation.** Per weight w:
     - one inverse V_w^{-1} on P[0..emono(w)];
     - all weight-w entries form one matrix (entries × emono(w)), solved by one blocked product with V_w^{-1};
     - the residual uses P[emono(w)..emono(w) + residual], points that this weight's solve did not use.

     Cost [V-crit c3]. The inverse at m = 1 218 takes 2.4 s (about 14 s at m = 2 178, scaling as m^3). The products do the same arithmetic as per-entry applies (entries · m² multiply-adds): for 8_19 that is 1.5-5 min single-threaded at the measured 6.8 ms per entry at m = 1 218. Batching improves cache use and splits over threads; it is not an asymptotic saving.
5. **Output and closure.**
   - `ClosureDataC { lam, objs: (p, h, t), types: TypeBasisC { ins, outs, gram }, edges: F, crossings: (positive, min label) }`, followed by `validate()`.
   - `close_plan` evaluates `cclose::exact_mem_estimate` and `delta_tower_estimate` for every requested (N, reduced) before closing.
   - `close` with Auto:
     - exact oracle if its estimate ≤ `mem_cap_bytes`;
     - else the Δ-tower if its estimate ≤ cap and (a, N) is G4-certified;
     - else Err, for example "closure refused at N=6 unreduced: exact 2.1 GB, tower 53 GB, cap 768 MB; cloud job".
   - **Laptop policy** (F6, §6.6): reduced N ≤ 6. Unreduced at every N where the exact estimate is under the cap: N ≤ 4 for the (2,2,2) knots, N = 5 for 4_1, and N ≤ 6 for the small T9 words. Anything else is listed as a cloud job, with the `close_plan` numbers in the report.
6. **Two primes.** The whole pipeline runs over F31 and F32b. The Poincaré polynomials must agree for every final result.

### 6.6 Expected costs and laptop envelope (from Appendix M and the critique probes; colored rows [E] until M6a measures them)

**Models.**
- **Per-point scan time** t_scan = (0.8-2.3)·10^-5 s × decategorified peak Σ S. This is the engine2 proxy [V-crit c4; re-run on 12a_146: 1.3·10^-5]. It is a lower bound for the colored engine.
- **Rescans:** emono(w_F) + 8, with w_F = Gram weight + 0..2. Wall time is for 2 threads, one prime and one point set.
- **Per-point memory:** peak Σ S × c_e (5-30) × 11.2 B per entry (flat layout, F10).
- **Closure memory:** `cclose::exact_mem_estimate` with the measured Gram density 0.20 [V-crit c5]. The Δ-tower estimate is not density-adjusted. The source is c1_scale (`colored_design/rev/c1_*.txt`, `crit_eng/c1_out.txt`). Cap: 768 MB.

| case | LB final / peak objects | peak Σ S | t_scan per point | memory per point | max Gram w (unknowns) | scans | rescan wall | exact closure, unreduced N = 4 / 5 / 6 | reduced N = 6 | where (M6a) |
|---|---|---|---|---|---|---|---|---|---|---|
| (2,2) T(2,n), n ≤ 7 | 11-57 / ≤ 81 | small | < 0.01 s | < 1 MB | 4 (14) | 24 | seconds | small | small | laptop |
| (3,3) Hopf / T(2,3) / T(2,5) | 13 / 24 / 42; peak 189 / 802 / 2 154 | 4 569 / 18 489 / 48 849 | 0.04-0.11 / 0.15-0.43 / 0.4-1.1 s | ≤ 21 MB (14.2 B per entry) | 9 (158) | 166-310 | ≤ 3 min | small (N ≤ 7) | small | laptop |
| (1,1,1) uncolored 3-strand | ≤ 216 | small | engine2-like | small | 3 (10) | 24 | seconds | small | small | laptop |
| T9 Markov family: [1,1,1,±2], [1,1,1,1,1,2], [1,1,1,1,2], [1,1,2] and their rotations | 21-57 / 21-57 | 615-1 785 | < 0.05 s | < 1 MB | 8 (294) | 302-638 | < 30 s | ≤ 0.02 / 0.07 / 0.22 GB | 0.02 GB | laptop (ignored-but-required) |
| T9 heavier words: [1,2,1,2,2] = [2,1,2,2,2] (braid relation) | 51 / 461 | 40 639 | 0.33-0.93 s | 2-14 MB | 12 (1 218) | 1 226-2 186 | 3-17 min | 0.16 / 0.74 / 2.5 GB | 0.17 GB | laptop if the time rule admits it, else cloud |
| T9 conjugate [−2,1,1,1,2,2] | 241 / 1 345 | 114 445 | 0.9-2.6 s | 6-38 MB | 12 (1 218) | 1 226-2 186 | 9-48 min | 0.17 / 0.80 / 2.7 GB | 0.18 GB | cloud |
| 4_1 | 77 / 325 | 26 733 | 0.21-0.61 s | 1.5-9 MB | 10 (630) | 638-1 226 | 1-6 min | 0.13 / 0.62 / 2.1 GB | 0.14 GB | laptop; unreduced N = 6 on the cloud |
| 5_2 | 191 / 927 | 84 893 | 0.68-1.95 s | 5-29 MB | 12 (1 218) | 1 226-2 186 | 7-36 min | 0.19 / 0.86 / 2.9 GB | 0.19 GB | cloud |
| 6_2 | 233 / 1 115 | 100 849 | 0.81-2.32 s | 6-34 MB | 12 (1 218) | 1 226-2 186 | 8-42 min | 0.19 / 0.87 / 2.9 GB | 0.20 GB | cloud |
| 6_3 | 253 / 1 797 | 128 453 | 1.03-2.95 s | 7-43 MB | 11 (882) | 890-1 646 | 8-40 min | 0.15 / 0.71 / 2.4 GB | 0.16 GB | cloud |
| 8_19 [1,1,1,2,1,1,1,2] | 287 / 2 361 | 245 557 | 1.96-5.65 s | 14-83 MB | 12 (1 218) | 1 226-2 186 | 20-103 min | 0.19 / 0.88 / 3.0 GB | 0.20 GB | cloud |
| 10-crossing (2,2,2) word of SYNTHESIS A.2 | 8 835 / 33 369 | 2 897 749 | 23-67 s | 0.16-0.97 GB (streaming required) | 12 (1 218) | 1 226-2 186 | 4-20 h | 0.41 / 1.9 / 6.4 GB | 0.43 GB | not in M6a: M6c first, then cloud |

Δ-tower estimates for unreduced closure are 0.5-0.6 GB at N = 4, 6.1-6.6 GB at N = 5 and 53-55 GB at N = 6 for 4_1 through 8_19, and 23 / 104 GB at N = 5 / 6 for the 10-crossing word. The tower is therefore no fallback for unreduced N ≥ 5. Reduced tower estimates stay ≤ 0.12 GB up to N = 6 (≤ 1.5 GB for the 10-crossing word).

**Outcomes the pipeline reports** (each an Err with the numbers):
- "cloud job: predicted … s" (time rule, §6.5 step 4);
- "cloud job: F samples need … MB" (memory rule);
- "needs CL2": emono(w_F) > 2 500, i.e. w_F ≥ 15 for (2,2,2). P6 shows w_F up to Gram + 1 on 2 strands, so 5_2, 6_2 and 8_19 reach it if F exceeds the Gram weight by 3 or more;
- "closure refused at N=…" (`close_plan`);
- "per-point memory …" (§6.3 step 4).

Notes:
- The F rescans dominate everything from 4_1 upward.
- M6c (coordinates of the final components) removes the rescans. It still needs emono((e − d_j)/2) points per component, which is 119-630 for component degrees 12-20. So for ≥ 9 crossings M6c is necessary, but it does not bring them to the laptop. This is measured input to G3.
- SYNTHESIS §1.3 lists "2-colored 3-strand, ≤ 12 crossings: laptop per knot". That row is wrong by this table and is superseded here (critique log, §9).

### 6.7 Tests (`tests/sengine.rs`; data `tests/data/sengine/`)

Oracles:
- xengine (M6.0, exact) on 2 strands;
- engine2 (Closure::Exact) uncolored and on reversed words (duality);
- `schur::decat_chi`;
- `wang_ref.txt` plus `wang_k3_extra.txt` (copied from this spec's probe output, Appendix R);
- the Hopf formula (copy of `hopf_formula` from `tests/xengine.rs`);
- the MOY value of each type (`moy::closure_value`, T0);
- two primes;
- braid moves, in particular Markov stabilization against xengine (T9).

**Default suite** (< 60 s debug; anything slower is ignored):

- **T0. Per-type closure data = MOY** (SYNTHESIS V15; critique 1; [V-R2]). This checks the closure generators and the Gram matrix of every final type, independently of any knot.
  - For each type p whose closure data the M6a pipeline produces (generators from the provider, Gram by generator-only interpolation), build a `ClosureDataC` with the single object (p, 0, 0), no edges, no crossings and p's `TypeBasisC`.
  - Required: the t^0 part of `close_exact_colored(.., N, None)`, q-graded, equals `moy::closure_value(&Ladder { n: N, bottom: λ, rungs: L_p })` exactly, and is 0 when the ladder passes a label > N.
  - The harness was validated on exact xengine data: 36/36 equal for (2,2), (2,1), (3,1) and (1,1), N = max label..5 (`rev_probe/r2_t0`).
  - Default:
    - every type of (2,2), (2,1), (1,2), (3,1) and (1,1,1) at N = max label..6;
    - the 4 types of (3,3) at N = 3..8;
    - the (2,2,2) types of fibre ≤ 36 at N = 2..5.

    Each part is moved to the ignored set if, by the G1(a) timings, it would push the suite over 60 s.
  - Ignored-but-required (must pass before any M6a.2 target runs):
    - all 21 (2,2,2) types at N = 2..6, and N = 7..10 where `try_close_exact_colored` stays under the cap (refusals listed as cloud jobs);
    - (3,3) at N = 9, 10;
    - both primes.
  - Cross-check on the (2,2,2) types of fibre ≤ 36: the generator degree multisets of the provider equal those of `cclose::rlambda_bases` (ins and outs). The single-type closures from both bases equal MOY.
- **T1. Equality with M6.0** (`xengine::xscan`, the oracle), for:
  - (1,1) σ^{±3};
  - (2,1) and (1,2) σ^{±2}, σ^4;
  - (2,2) σ^{±1}, σ^{±3}, σ^5;
  - (3,1) and (1,3) σ^2;
  - 8 random words from the M6.0 review (labels (2,2), (2,1), (1,2), lengths 2-5, seed 11).

  Required:
  - `final_objects` identical (multiset of (table, h, t));
  - K_0 identical;
  - Poincaré for N = 1..6 unreduced and reduced (block 0) identical;
  - F31 = F32b;
  - streaming on and off identical (final objects, K_0, Poincaré).
- **T2. Uncolored 3-strand = engine2**, for:
  - 4_1 [1,-2,1,-2];
  - the trefoil stabilized [1,1,1,2];
  - 8_19 [1,1,1,2,1,1,1,2] uncolored;
  - 6 random 3-strand words (seed 12, lengths 3-7).

  Required:
  - Poincaré N = 2..6, unreduced and reduced, identical to `khr_engine2_closure(.., Closure::Exact)`;
  - final object multiset identical to engine2's (type x ∈ S_3 ↔ permutation table, the M3a T-e convention; (h, q) ↔ (h, t));
  - streaming on and off identical.

  If the q-convention of engine2's objects differs by a global constant, the test documents it once on the trefoil; Poincaré equality stays binding.
- **T3. Hopf formula (V4)** for the 13 R3 (N; a, b) cases (2;1,1), (3;1,1), (3;2,1), (3;2,2), (4;2,1), (4;2,2), (4;3,1), (4;3,2), (5;2,1), (5;2,2), (5;3,2), (6;3,2), (6;3,3), bigraded. The (3,3) and (3,2) cases are ignored when slower than the budget.
- **T4. Per stage:** K_0 = decategorified scan and #objects >= LB at every stage of every scan above (built into `sscan`; the test asserts that `stages` is populated).
- **T5. Markowitz heap = rescanning order** (unit), and the term-by-term streaming unit test of §6.4.
- **T6. Degree histogram and Measured** are populated; c_e > 0; `zero_block_frac`, `peak_point_bytes` and `f_sample_bytes` are set.
- **T6b. Refusals.** Unit tests on synthetic inputs, without scans:
  - the time rule refuses 8_19's parameters (t_scan = 2 s, w_F = 12, points0 = 24, 2 threads → predicted > 600 s) and names the command;
  - the memory rule refuses an F-sample estimate above `f_mem_bytes`;
  - `close_plan` on 4_1-sized synthetic closure data reports exact > cap at unreduced N = 6, and Auto's Err names both estimates and the cap;
  - emono(w_F) > 2 500 gives "needs CL2".

**Ignored tests** (run once each on the laptop when predicted < 10 min and < 1 GB, otherwise cloud jobs):

- **T7. 3-colored 2-strand.**
  - (3,3) Hopf: bigraded = hopf_formula, N = 3..6.
  - (3,3) T(2,3): bigraded = Wang k = 3, N = 3..7 (ranks 1, 10, 46, 136, 307).
  - (3,3) T(2,5): Wang k = 3, N = 4, 5, 6 (ranks 16, 106, 396).
  - (3,3) T(2,3) at N = 4: equal to engine2's uncolored gl_4 of the reversed word (via duality (expected)).
  - χ = decat_chi.
  - Two primes.
- **T8. 2-colored 3-strand knots** 4_1, 5_2, 6_2, 6_3, 8_19 (words from `data/knots_upto12.jsonl`). Checks:
  - χ = decat_chi for N = 2..6 (Appendix R.1 lists N = 3);
  - N = 2: t^0 q^0 (wedge^N rule, V2);
  - N = 3, unreduced and reduced: bigraded equal to engine2's uncolored gl_3 of the reversed word (Appendix R.2), labelled "via duality (expected)". This checks closure data only up to weight 6 (R.2);
  - two primes, and two point sets (seeds 7 and 8), giving identical final objects and Poincaré polynomials;
  - Lewark-type relation [N choose 2] P_red − P_unred = (1+t)Q with Q ≥ 0, asserted where both closures exist (as M6.0 after its fix);
  - T0 has passed for every final type.

  Where each runs (§6.6):
  - **Laptop: 4_1.** Three separate ignored tests, each predicted at 1-6 min of rescans: (F31, seed 7), (F32b, seed 7), (F31, seed 8). Closures: reduced N ≤ 6; unreduced N ≤ 4, and N = 5 when `close_plan` is under the cap.
  - **Cloud: 5_2, 6_2, 6_3, 8_19**, and 4_1 unreduced at N = 6. On the laptop `sscan` refuses them by the time rule, which T6b covers. The report gives the exact commands, the commit and the seeds.
- **T9. V10 and the Markov family** (ignored-but-required for M6a.2). These words test 3-strand closure data at all weights, with no N = 3 truncation, against an independent exact engine. Bigraded equality at N = 2..6, unreduced and reduced (block 0 marked on both sides):
  - Laptop (< 30 s of rescans each, §6.6):
    - [1,1,1,2] and [1,1,1,−2] equal xengine (2,2) σ^3;
    - [1,1,1,1,1,2] equals (2,2) σ^5;
    - [1,1,1,1,2] equals (2,2) σ^4. This is the T(2,4) link: components {strand 1} and {strands 2, 3}; block 0 is strand 1 in both engines;
    - [1,1,2] equals (2,2) σ^2 (Hopf link), and also `hopf_formula`;
    - the cyclic rotations [2,1,1,1] and [1,2,1,1] equal [1,1,1,2];
    - mirror: [−1,2,−1,2] equals 4_1 mirrored (P(t,q) ↦ P(1/t,1/q)), with 4_1 from T8.
  - Laptop if the time rule admits it, else cloud: braid relation [1,2,1,2,2] = [2,1,2,2,2] (3-17 min of rescans each).
  - Cloud: the conjugate [−2,1,1,1,2,2] = (2,2) σ^3; all cyclic rotations of 5_2.
  - Unreduced N = 6 is closed where `close_plan` admits it (the stabilized words need about 0.22 GB).
- **T10. Cloud references.**
  - Full cube for the 2-colored 4_1 at N = 3, 4 (SYNTHESIS §6), once it exists.
  - Cabling, validated against xengine: `colored_cable 3,4,5 --braid '1 -2 1 -2' --strands 3 --colors 2 --closure exact` (4_1) and `colored_cable 3,4,5 --braid '1 1 1 2 -1 2' --strands 3 --colors 2 --closure exact` (5_2). These are 6-strand cables: 4_1 passed 4.4 GB within 22 s on the laptop [critique 1], so they run on the cloud only, recording peak RSS.
  - Required: bigraded equality with T8 at N = 3..5, unreduced and reduced.

### 6.8 Measured outputs (recorded in the report and in `Measured`)

- Colored c_e (stored block entries / Σ S(q,q)) at peak, per case; the zero-block fraction of the LocMaps; peak per-point bytes (Eval memo + complex + layouts); whether streaming was on.
- F sample bytes, interpolation seconds, and the `close_plan` numbers per (N, reduced).
- **Calibration run** (critique 2). After G1(a) passes: one single-point scan, no closure, of the SYNTHESIS A.2 10-crossing (2,2,2) word with streaming on. It runs under the watchdog at 1.0 GB and 10 min and records time, peak RSS, peak stored entries and the stage reached. If it is killed, the record is "> 1 GB or > 10 min at stage k", and the word is confirmed as a cloud job. The measured seconds per unit of peak Σ S replace the engine2 proxy in §6.6.
- The LB/actual ratio per stage.
- The degree histogram of final components (U10).
- Peak RSS per point and scan seconds per point.
- The rescans used; max Gram and F weights.
- DecStats (decompositions, identity shortcuts, provider time).
- At the end of M6a, **G2** (SYNTHESIS §3.1): whether every Hom the M6a runs needed came from Search (no Exact fallbacks beyond fibre about 40), and whether Search certifies the G1(b) (2,2,2,2) set within budget.

### 6.9 Done criteria

- The default suite passes.
- T0 passed for every (2,2,2) and (3,3) type, including the ignored-but-required part. T9's laptop words passed.
- T7 to T10 ran: locally when within the envelope of §6.6, otherwise queued as cloud jobs with commands, commit and seeds.
- Every 3-strand result is labelled "computed" and also passed χ, duality (N = 3, weight ≤ 6), two primes and two point sets. Its types passed T0.
- The report contains:
  - the tables of §6.6 with measured values, including the calibration run;
  - mismatches with minimal examples;
  - the G2 verdict;
  - the recommendation for M6b/M6c.

---------------------------------------------------------------------------------------------------

## 7. Where each SYNTHESIS §1.4 decision lives

| # | decision | algorithm (section) | test |
|---|---|---|---|
| 1 | stacked evaluation at M(p,q,d)+2 points for every independence decision on Hom^d | `hom_gens`, `m_points`, `detecting_cosets` (witnesses with the counting bound, fallbacks A/B), `stacked_rank`, `certify_full` (§4.3 items 1-4) | B2, B3 (exact M), B4 (one point fails, full and witness cosets), B4b (fallback) |
| 2 | idempotents e = P(a), never re-solved; e² = e at every new point | §5.2 "Idempotent". Proof that ranks certify Im e ≅ B_p: Im e is a summand of L_p = B_p ⊕ (B_r's with r < p); rank 1 at the top stratum gives exactly one B_p; the stratum ranks then force the B_r multiplicities to 0 | D2; e² = e at 3 fresh points |
| 3 | polynomiality certificates for templates and closed forms | `poly_certificate` (§3.4.5), `assert_genuine` (§4.3 item 9) | A10 (passes; fails on 3 perturbations), B5 |
| 4 | constant scalars only | `Lin` constants; closed forms normalized so that every recorded constant is 1 (§3.4.3, F8); no per-point normalization anywhere | A7 (60 cases, c = 1), A8, A10(b), A12 |
| 5 | streaming-elimination rule | §6.4 eligibility (on by default above peak Σ S 10^5) | T5 (term-by-term arrival, 200 seeds; broken rule detected), T1/T2 streaming on = off |
| 6 | genuineness (rank does not certify it) | every leaf certified; composites inherit; Rickard d² = 0; pairs certified by the intersection form plus genuine leaves | A9, A10, B5, D3 checks |
| 7 | closure generators by the one-point determinant criterion | §4.3 item 7; §6.5 step 2 (re-verified at every scan point) | T0 (per type, q-graded MOY; degree multisets = `rlambda_bases`), T1, T2 (regression against xengine's exact `rlambda_bases` degrees and against engine2's `select_basis` results through identical Poincaré polynomials) |
| 8 | pivot at any coset of the top stratum (no identity coset mid-scan) | §6.4, w = rep(type) | T1 with mixed labels |
| 9 | component degree e = t_src − t_tgt; total q-degree 0 | §2 item 4, §6.3 step 4 | asserted per crossing |
| 10 | object model: ladders with zero labels, any number of uprights | `Expr` (§3.3.2) | A2 (labels 0..3) |
| 11 | point-dependence detection | determinism trace (§6.3); intersection forms re-checked at a second point (§4.3 item 6) | T8 (two point sets) |
| 12 | closure data per type: dim Φ_N(B_p) = MOY (V15) | single-object closure of each type's TypeBasisC (§6.7 T0) | T0 (q-graded, N = 2..6 or more, two primes) |
| 13 | closure data at high weight (beyond the N = 3 duality's weight 6) | Markov stabilization and conjugation against xengine | T9 (laptop words binding), T10 (cable, cloud) |

---------------------------------------------------------------------------------------------------

## 8. Risks specific to M4a-M6a

| # | risk | signal | mitigation |
|---|---|---|---|
| Q1 | Search does not span (U1), in particular the 2 (2,2,2) families whose extremes are 2-dimensional, and high-degree non-extreme summands (φ degrees 3-7, dims up to 199 [V-P4]) | G1(a) Missing lists | S3-S5 tiers; then M4c (EKLP) with the failing targets as its acceptance set |
| Q2 | Window templates for (3,3) are slow (hom_space at fibre 30: 86 s; Rickard 120 s per sign) | M4a/M4b timings | disk caches per prime; build once (cloud if needed); the closed forms E2/E3/E4 avoid most windows |
| Q3 | Closure F rescans: 4_1 1-6 min; 5_2-8_19 7-103 min; 10-crossing words 4-20 h (per prime and point set) | t_scan measured in M6a; calibration run (§6.8) | time and memory rules make them cloud jobs; M6c is a prerequisite for ≥ 9 crossings; batched interpolation |
| Q4 | Dense interpolation above about 2.5e3 unknowns (w_F ≥ 15 for (2,2,2)) | Err "needs CL2" | the CL2 spike (SYNTHESIS) before targets with w_F ≥ 15 |
| Q5 | Point-dependent behaviour (a template or Search candidate silently non-genuine) | form mismatch at the second point; trace mismatch; residuals; d² ≠ 0 | hard errors; certificates |
| Q6 | Multiplicity disagreement between Search pairs and `term_mult` | Err in dec_rung | blocking; investigate (p-canonical K8 or a provider bug); second prime |
| Q7 | Per-point memory (Eval memo plus complex) | `memo_bytes`, `LocMap::bytes`, assertion after each crossing | flat layout (11.2 B per entry); streaming on above peak Σ S 10^5; `Eval::clear` between crossings; per-point budget 0.5 GB |
| Q8 | Unreduced closure at N ≥ 5 above the 768 MB oracle cap; the Δ-tower needs 6-55 GB | `close_plan` | laptop policy (reduced N ≤ 6, unreduced N ≤ 4 or 5); cloud; CL3 / structured closure before unreduced N ≥ 6 claims |
| Q9 | `hom_space` calls of minutes that cannot be interrupted | `exact_calls`, `exact_secs`, watchdog kills | admission by label and cost table (10 s ceiling); Search S0-S2 first for labels ≥ 3; one process per G1 pair |
| Q10 | Search memory from visited words | `peak_index_bytes`, `peak_arena_bytes`, `MissKind::Memory` | scratch regions (`mark`/`commit`), lazy packed indices with an LRU bound; memory failures counted as engineering in G1(b) |
| Q11 | Gram density of Search-produced generators differs from the 0.20 measured on exact bases, which moves the closure estimates | `close_plan` actual numbers | the estimates are recomputed on the real data before closing; refusals are explicit |

---------------------------------------------------------------------------------------------------

## 9. Critique log (revision 2)

Two adversarial critiques reviewed revision 1: critique 1 (theory/correctness) and critique 2 (engineering/scale). Both returned "adopt with changes". Each blocker and major issue below was verified before deciding. Minor issues were adopted after a check of the reasoning.

**Verification runs** (all under the watchdog `colored_design/watch.sh`, nice 19, at most one process at a time, each under 60 MB and 20 s):
- `rev_probe/r1_e2all` (new): E2/E3 constants for all 60 cases with a + b ≤ 6, both orientations. Output `rev/r1_e2all.txt`: 60 OK, 0 mismatches, 0.3 s.
- `rev_probe/r2_t0` (new): single-type graded closure against MOY on exact xengine data. Output `rev/r2_t0.txt`: 36/36 equal.
- Re-runs of the critique probes (outputs in `rev/`):
  - `e2const`: identical output;
  - `witness` on (2,2), (3,3), (2,2,2) and (2,1,2): identical counts. New: (2,1,2) first violates at d = 12;
  - `c2_store`: 59.6 against 11.2 B per entry;
  - `c3_interp`: m = 630: 0.26 s; m = 1 218: 2.44 s;
  - `c4_e2` on 12a_146: 0.015 s per point, 1.3·10^-5 per unit of peak Σ S;
  - `c6_flowidx`: 225 / 243 / 285 B per flow;
  - `c1_scale` on 4_1, 5_2, 6_2, 6_3 and nine Markov-family words (`rev/c1_*.txt`).
- Source checks:
  - `src/engine2.rs:73` (`par_map` requires `Fun: Sync`);
  - `tests/xengine.rs:386-399` (the existing MOY type-dimension oracle);
  - `cclose.rs` (`exact_type_dims_colored` returns ungraded dimensions; `DEFAULT_MEM_CAP_MB` = 768);
  - the M6.0 report (the (3,2) one-step hom_space call took 275 s at fibre product 1 152).

| # | source, severity | issue (short) | verified | decision | where |
|---|---|---|---|---|---|
| E1 | crit. 2, blocker | F rescans cannot reach 5_2-8_19 within 600 s, nor 10-12 crossings on the laptop; SYNTHESIS §1.3 claims a laptop envelope up to 12 crossings | yes: c4 proxy re-run (1.3·10^-5 s per unit Σ S); c1 peak Σ S re-run for 4_1, 5_2, 6_2, 6_3 | **adopted.** Laptop envelope stated: 2 strands, (3,3), the Markov family, 4_1. 5_2-8_19 are cloud jobs. ≥ 9 crossings need M6c first and stay on the cloud. A calibration single-point scan of the A.2 10-crossing word is added. SYNTHESIS §1.3 is superseded here (this spec does not edit SYNTHESIS; the orchestrator should update it) | §0 (revision note, F6, decision 7), §6.1, §6.6, §6.8, Q3 |
| E2 | crit. 2, blocker | unreduced closure at N = 5, 6 exceeds the oracle cap; the tower is worse; Auto fails with a bare Err; no cost column | yes: c1 estimates re-run; density 0.20 taken from c5 (not re-run: 44 s and not needed for the decision) | **adopted.** `close_plan` with both estimates; Auto's Err names them; closure-cost columns in §6.6; laptop policy reduced N ≤ 6, unreduced N ≤ 4 (5 under the cap); unreduced N = 5, 6 for larger knots on the cloud | §0.3 item 6, §6.2, §6.5 step 5, §6.6, T6b, T8, Q8 |
| E3 | crit. 2, major | `BTreeMap<u64, Vec<Vec<F>>>` layout costs 7.5× the payload | yes: c2 re-run, identical | **adopted.** Binding `Block` (dense row-major), `Layout` (shared u32 offsets), `LocMap` (one `Vec<F>`); `gemm_acc` on slices; streaming on by default above peak Σ S 10^5; per-point memory assertion | §1.2, §3.3.4, §3.4.2, §6.2, §6.3, A16, Q7 |
| E4 | crit. 2, major | F sample storage and interpolation unspecified (memory and time blow-up) | yes: c3 re-run | **adopted.** Point set fixed in advance (emono(w_F) + 8); per-entry samples at the entry's own weight only; scans dropped after extraction; one inverse per weight, applied as a blocked product; residual on points unused by that weight's solve; memory formula asserted before the rescans. The critique's implied speed-up of batching is corrected: the arithmetic is the same, the gain is cache use and threads | §6.5 step 4 |
| E5 | crit. 2, major (= crit. 1 minor 6) | `DecData` returns `Rc`, so it is not `Sync`; `par_map` needs `Sync` | yes: `engine2.rs:73` | **adopted.** `Arc` everywhere in `DecData`; `Arena`/`TmplLib` without interior mutability (template memo moved to `Eval`); `reduced` returns an owned value; compile-time `Send + Sync` asserts (A16, D9) | §0.3 item 5, §1.2, §3.3.3, §3.3.4, §5.1, A16, D9 |
| E6 | crit. 2, major | Search interns every visited word with a full FlowIndex, never freed; G1(b) could fail on memory and be misread as incompleteness | yes: c6 re-run | **adopted with a modification.** Scratch regions (`mark`/`commit`); packed FlowIndex (8 B per flow and rung, about 3-4× smaller); indices built only for evaluated words, with an LRU bound of 256 MiB; `MissKind::Memory` counted as engineering in G1(b). **Not adopted:** per-coset lazy indices. They need backward flow enumeration, which is a non-goal of M4a. The per-expression laziness plus packing reaches the same order of saving | §3.3.2, §3.3.4, §3.4.1, §4.2, §4.4, §4.7, A15, Q10 |
| E7 | crit. 2, major | Exact caps admit 50-300 s hom_space calls at labels ≥ 3; Exact always runs first; calls cannot be interrupted | yes: M6.0 report (275 s and 68 s at fibre product 1 152) | **adopted.** Exact admission: max label ≤ 2 under the caps, or a measured `CostTable` with a 10 s ceiling (bucket maximum, since cost does not follow fibre); otherwise `Refused` without a call. `Chain::for_labels` puts Search S0-S2 first for labels ≥ 3. Interruption: one process per G1 pair under the watchdog (no child-process machinery in the library) | §0.3 item 4, §4.2, §4.5, B8, G1(a) fallback text, Q9 |
| T1 | crit. 1, major | no per-type closure-data oracle (dim Φ_N(B_p) = MOY) | yes: the oracle exists for M6.0 at q = 1; I validated a q-graded single-object form, 36/36 | **adopted** as T0. It is q-graded via the single-object `ClosureDataC` (no repository change; `exact_type_dims_colored` is ungraded). Default for 2-strand, (3,3) and small (2,2,2) types; ignored-but-required for all 21 (2,2,2) types before M6a.2; cross-check of degree multisets with `rlambda_bases` | §0.2 F9, §6.1, §6.7 T0, §6.9, §7 rows 7 and 12 |
| T2 | crit. 1, major | N = 3 duality sees only weight ≤ 6; no binding high-weight 3-strand oracle; cabling not mentioned | yes: weight bound re-derived; Markov-family costs measured with c1 (peak Σ S 615-1 785, Gram weight 8) | **adopted.** T9 Markov family (stabilizations, T(2,4), Hopf, rotations, mirror) is ignored-but-required and runs on the laptop in < 30 s each. The braid-relation pair runs if the time rule admits it. The conjugate [−2,1,1,1,2,2] is a cloud job (peak Σ S 114 445). T10 lists the cable commands for 4_1 and 5_2 as cloud jobs. R.2 states the weight-6 limit | §6.6 rows, §6.7 T8-T10, Appendix R.2, §7 row 13 |
| T3 | crit. 1, major | closed-form constants certified only for few label pairs; G1(a) windows need (4,1), (4,2), (2,4), (3,3), (5,1), (4;3), (5;b), (6;3) | yes: re-ran `e2const`; extended to all 60 cases (`r1_e2all`), all matching the rule | **adopted.** Forms redefined with e(B,A) and e(S, A∖S) so that every constant is 1. A7 covers every a + b ≤ 6 by default (0.3 s release). a + b ≤ 9 is ignored, before G1(b). `assert_genuine` rejects unrecorded (kind, labels). This closes the earlier open question on E2/E3 signs | §0.2 F8, §3.3.4, §3.4.3, A7, A15, §4.3 item 9 |
| T4 | crit. 1, minor | witness detecting set fails the counting bound at high degree; B4 values were full-coset values | yes: witness re-run; C1 I1 table read | **adopted.** Counting assertion with fallbacks A/B; B4 lists both value sets (witness: 5/14, 35/48, 5/21); new B4b | §0.2 F11, §0.3 item 4, §4.3 item 3, B4, B4b, §7 row 1 |
| T5 | crit. 1, minor | the degree identity had the wrong sign for π | yes: re-derived | **adopted.** e(π_τ) + e(f) + e(ι_σ) = t_src − t_tgt | §6.3 step 4 |
| T6 | crit. 1, minor | the T5 streaming test did not model term-by-term arrival | yes: re-derived | **adopted.** Incremental generator with genuine d_W; sensitivity check that fails if the broken rule is never caught; streaming on = off on every T1/T2 word | §6.4, T1, T2, T5, §7 row 5 |
| T7 | crit. 1, minor | interpolation details: zero entries of odd or negative weight; F residual count; w_F ≥ 15 above the cap | yes | **adopted.** Asserted zeros; residual fixed at 8 points; "needs CL2" outcome listed | §6.5 steps 3-4, §6.6 outcomes, Q4 |
| S1 | crit. 2, minor | polynomiality certificate factorizes per entry and redraws single points | yes: c3 timings | **adopted.** One point set; one inverse per distinct D; redraw the whole set | §3.3.3 (`CertReport`), §3.4.5 |
| S2 | crit. 2, minor | G1(b) "exhausted below budget" depends on the Budget constants | yes | **adopted.** Automatic rerun with l_max + 2 and 4·c_max; only a second `Exhausted` counts against completeness | §4.7 |

**Open items handed on:**
- SYNTHESIS §1.3 (laptop envelope row) should be corrected by the orchestrator.
- The Gram density of Search-produced generators is unmeasured until M6a (Q11).
- The colored per-unit scan cost replaces the engine2 proxy after the calibration run (§6.8).

---------------------------------------------------------------------------------------------------

## Appendix M. Measurements (probe crate `colored_design/proto_spec`, release, nice 19, one process at a time, watchdog 1.1-1.2 GB)

Commands:
- Build: `cd .../colored_design/proto_spec && CARGO_TARGET_DIR=../target_spec CARGO_BUILD_JOBS=2 nice -n 19 cargo build --release`.
- Binaries in `../target_spec/release/`.
- Raw outputs: `p1_*.txt`, `p4_*.txt`, `p10_refs.txt`, `wang_k3_extra.txt` in `proto_spec/`.

**P1 (`p1_decomp LAMBDA [v]`): rung-by-rung decomposition targets.**

Every (type, rung) pair that a crossing term on λ can produce; all types of the reached blocks; Hom dims on ladder characters by `schur::hom_gd`; M = stacked points of §4.3. Run times < 7 s, at most 81 MB.

| λ | blocks | types | pairs | identity | (q,s) targets | Σ mult | Σ dim Hom^s(L_q, X) | max X flows | max dim | max M |
|---|---|---|---|---|---|---|---|---|---|---|
| (2,2) | 3 | 6 | 9 | 4 | 18 | 21 | 139 | 12 | 35 | 14 |
| (3,3) | 4 | 10 | 18 | 6 | 72 | 121 | 5 382 | 90 | 787 | 158 |
| (1,1,1) | 3 | 12 | 18 | 8 | 18 | 18 | 40 | 6 | 5 | 3 |
| (2,2,2) | 5 | 63 | 126 | 28 | 358 | 404 | 7 330 | 432 | 290 | 39 |
| (3,3,3) | 7 | 219 | 494 | 65 | 3 151 | 5 146 | 2.65e6 | 31 500 | 51 909 | 915 |
| (2,2,2,2) | 7 | 1 110 | 2 520 | 329 | 9 255 | 10 322 | 533 795 | 51 840 | 4 810 | 176 |

Extreme maps (F4): the lowest-degree ι and highest-degree π are 1-dimensional, with multiplicity 1 at s_lo, in (2,2) 6/6, (3,3) 16/16 and (2,2,2) 156/158 families.

**P2 (`p2_indec LAMBDA [MU]`): decomposable canonical ladders.**

| λ | types | decomposable | Σ flows / Σ rank |
|---|---|---|---|
| (2,2,2) | 21 | 0 | 1 045 / 1 045 |
| (2,1,2) | 11 | 0 | 199 / 199 |
| (1,1,1,1) | 24 | 0 | 219 / 219 |
| (2,2,1,1) | 58 | 0 | 3 969 / 3 969 |
| (3,3,3) | 55 | 1 (0,1,2/1,1,1/2,1,0: 10 500 / 8 820, End^0(L_p) = 2) | 78 079 / 76 399 |
| (2,2,2,2) | 282 | 6 (max End^0 11, max ratio 1.39) | 402 109 / 397 687 |

Basis times: (2,2,2,2) 1.7 s and 81 MB.

**P3 / P12 (`p3_homtime`, `p12_unknowns`): exact hom_space.** See F2.

**P4 (`p4_adj LAMBDA [v]`): adjunction.**
- Graded Hom(B_q, B_p R) = graded Hom(B_{q'}, B_{p+}) up to one shift in all cases: (2,2) 10, (3,3) 22, (1,1,1) 20, (2,2,2) 186, (2,1,2) 114.
- q' = p+ for (2,2,2) in 132 of 386 (q,s) cases.
- Otherwise the needed φ ∈ Hom^k(B_{q'}, B_{p+}) has (k, dim) = (1, 1) × 74, (2, 1) × 42, (3, ...) × 56, (4, ...) × 30, (5, ...) × 22, (6, ...) × 14, (7, ...) × 16 (total 254), with dims up to 199.

**P5 (`p5_local LAMBDA MU [RUNG]`): local.rs** (dense, one point).

| ladders | max fsize | max flows per coset | loc_matrix + inverse |
|---|---|---|---|
| (2,2,2) canonical | 216 | 8 | 4 ms (all 21) |
| (2,2,2) X = L_p E_0^(1) | 432 | 20 | 22 ms (all 21) |
| (2,2,2) → (3,1,2) | 80 | 4 | 1 ms |
| (3,3) canonical | 30 | 3 | < 1 ms |
| (3,3) X = L_p E_0^(2) | 90 | 18 | < 1 ms |
| (3,3,3) | 4 320 | 54 | 10.6 s |
| (3,3,3) | 8 820 / 10 500 | - | killed above 1.2 GB |

**P6 (`p6_fweights`): closure-data weights** (exact 2-strand data from xengine; engine2 uncolored data).

| case | objects | max Gram w | max F w (actual / degree bound) | component degrees |
|---|---|---|---|---|
| (2,2) σ^3 | 13 | 4 | 2 / 3 | 1:4, 2:9, 3:4, 4:6, 6:1 |
| (2,2) σ^5 | 31 | 4 | 4 / 4 | 1..8 |
| (2,2) σ^7 | 57 | 4 | 5 / 5 | 1..10 |
| (2,1) σ^4 | 5 | 2 | 2 / 2 | 2:3, 4:1 |
| (3,1) σ^2 | 3 | 3 | 0 / 1 | 2:1, 3:1 |
| (1,1,1) 4_1 (engine2, 24 points, 0.1 s) | 12 | 2 | 1 / 1 | 1:12, 2:10 |
| (1,1,1) 10-crossing 3-braid | 216 | 3 | 2 / 2 | 1:274, 2:956, 3:104 |
| (1,1,1,1) Conway (2.0 s) | 202 | 4 | 3 / 3 | 1:476, 2:1032, 3:170 |

**P7 (`p7_targets`): decategorified sizes of the targets** (Gram rows with block 0 marked for "reduced").

| case | LB | types | LB ldim | LB Σ S | peak objects / ldim / Σ S | Gram weights | max unknowns (reduced) |
|---|---|---|---|---|---|---|---|
| 4_1 (2,2,2) | 77 | 16 | 3 049 | 6 073 | 325 / 17 641 / 26 733 | 0..10 | 630 (91) |
| 5_2 | 191 | 21 | 10 357 | 18 289 | 927 / 58 513 / 84 893 | 0..12 | 1 218 (140) |
| 6_2 | 233 | 21 | 11 077 | 20 025 | 1 115 / 69 433 / 100 849 | 0..12 | 1 218 (140) |
| 6_3 | 253 | 17 | 11 401 | 21 217 | 1 797 / 87 841 / 128 453 | 0..11 | 882 (112) |
| 8_19 [1,1,1,2,1,1,1,2] | 287 | 21 | 20 965 | 39 097 | 2 361 / 175 669 / 245 557 | 0..12 | 1 218 (140) |
| 8_19 [1,2]^4 | 187 | 21 | 13 933 | 27 753 | 1 847 / 130 117 / 187 221 | 0..12 | 1 218 (140) |
| (3,3) T(2,3) | 24 | 4 | 507 | 693 | 802 / 16 695 / 18 489 | 0, 5, 8, 9 | 158 (12) |
| (3,3) T(2,5) | 42 | 4 | 911 | 1 253 | 2 154 / 44 667 / 48 849 | same | 158 (12) |
| (3,3) Hopf | 13 | 4 | 265 | 373 | 189 / 3 969 / 4 569 | same | 158 (12) |
| (3,2) Hopf | 7 | 3 | 57 | 61 | 65 / 625 / 641 | 0, 4, 6 | 39 (4) |

e-monomials of R^(2,2,2) by weight 10..16: 630, 882, 1 218, 1 638, 2 178, 2 838, 3 663.

**P8 (`p8_prefix`): common bottom prefix of L_q and L_p R** ((2,2,2), 158 families). Length histogram {0: 60, 1: 41, 2: 22, 3: 28, 4: 6, 5: 1}. The remaining-window fibre products reach 34 560. Prefix windows alone (tier S6) cannot be the main source.

**P9 (`p9_rlb`): exact closure generators on (2,2,2) canonical ladders.**
- Fibre <= 18: <= 0.02 s.
- Fibre 36: 0.09, 0.16, 0.16, 0.23, 0.25, 0.80, 6.28 and 40.3 s.
- Fibre 48: 0.35 s.
- Fibre 60: 2.5 and 145.5 s.
- The 4 types of fibre 72-216 were not reached within 400 s.

**P11 (`p11_cdiag`): closure-generator targets.** Σ_p c_diag(b_p)(1) (max per type; max degree):

| λ | Σ | max per type | max degree |
|---|---|---|---|
| (2,2) | 4 | 2 | 4 |
| (3,3) | 8 | 3 | 9 |
| (2,1,2) | 18 | 4 | 8 |
| (1,1,1) | 6 | 1 | 3 |
| (2,2,2) | 48 | 8 | 12 |
| (3,3,3) | 420 | 54 | 27 |
| (2,2,2,2) | 1 748 | 63 | 24 |

---------------------------------------------------------------------------------------------------

## Appendix R. Reference values for M6a (from validated oracles; P10 and `r3scripts/wang.py`)

### R.1 Decategorified χ (`schur::decat_chi`, λ = (2,2,2), convention A)

At N = 2 all five knots give 1. At N = 3:

| knot (word) | χ at N = 3 |
|---|---|
| 4_1 [1,-2,1,-2] | q^8 + q^6 − 1 + q^-6 + q^-8 |
| 5_2 [1,1,1,2,-1,2] | q^-2 + q^-6 + q^-8 + q^-10 + q^-12 − q^-18 − q^-20 |
| 6_2 [1,1,1,-2,1,-2] | q^4 + q^2 + 1 + q^-2 − q^-4 − q^-8 + q^-16 |
| 6_3 [1,1,-2,1,-2,-2] | −q^10 + 2q^2 + 1 + 2q^-2 − q^-10 |
| 8_19 [1,1,1,2,1,1,1,2] | q^-10 + q^-12 + 2q^-14 + 2q^-16 + 2q^-18 − q^-22 − 2q^-24 − 2q^-26 − q^-28 + q^-32 |

N = 4, 5, 6: `proto_spec/p10_refs.txt`; the tests call `decat_chi` directly.

### R.2 Duality at N = 3 (engine2 `Closure::Exact`, uncolored gl_3 of the reversed word; "via duality (expected)")

The Euler characteristics agree with R.1 (checked by hand for 4_1 and 6_2).

The duality itself is sound bigraded on 2 strands [critique 1]: xscan (2,2) on '1 1 1', '-1 -1 -1' and '1 1 1 1 1' at N = 3 equals engine2 2-strand gl_3 exactly, unreduced and reduced (ranks 7/3, 7/3, 11/5).

**Limitation.** At N = 3 each label-2 block closes through H*(Gr(2,3)), whose basis has weights 0, 1 and 2. So every e-monomial of total weight > 6 vanishes in the closure. R.2 checks the 3-strand closure data only up to weight 6, while Gram weights reach 12 and F weights reach or exceed them. T0 (per type, N up to 6-10) and T9 (Markov family against xengine, N = 2..6) cover the higher weights.

| knot | unreduced (rank) | reduced (rank) |
|---|---|---|
| 4_1 | t^-2(q^6+q^8) + t^-1(q^0+q^2) + t^0(q^-2+q^0+q^2) + t^1(q^-2+q^0) + t^2(q^-8+q^-6) (11) | t^-2q^6 + t^-1q^2 + t^0q^0 + t^1q^-2 + t^2q^-6 (5) |
| 5_2 | t^0(q^-6+q^-4+q^-2) + t^1(q^-6+q^-4) + t^2(q^-12+q^-10+q^-8+q^-6) + t^3(q^-14+q^-12) + t^4(q^-14+q^-12) + t^5(q^-20+q^-18) (15) | t^0q^-4 + t^1q^-6 + t^2(q^-10+q^-8) + t^3q^-12 + t^4q^-14 + t^5q^-18 (7) |
| 6_2 | t^-2(q^2+q^4) + t^-1(q^-4+q^-2) + t^0(q^-6+q^-4+2q^-2+q^0) + t^1(q^-8+2q^-6+q^-4) + t^2(q^-12+q^-10+q^-8+q^-6) + t^3(q^-14+q^-12+q^-10+q^-8) + t^4(q^-16+q^-14) (23) | t^-2q^2 + t^-1q^-2 + t^0(q^-4+q^-2) + 2t^1q^-6 + t^2(q^-10+q^-8) + t^3(q^-12+q^-10) + t^4q^-14 (11) |
| 6_3 | t^-3(q^8+q^10) + t^-2(q^2+q^4+q^6+q^8) + t^-1(q^0+q^2+q^4+q^6) + t^0(2q^-2+3q^0+2q^2) + t^1(q^-6+q^-4+q^-2+q^0) + t^2(q^-8+q^-6+q^-4+q^-2) + t^3(q^-10+q^-8) (27) | t^-3q^8 + t^-2(q^4+q^6) + t^-1(q^2+q^4) + 3t^0q^0 + t^1(q^-4+q^-2) + t^2(q^-6+q^-4) + t^3q^-8 (13) |
| 8_19 | t^0(q^-14+q^-12+q^-10) + t^2(q^-16+q^-14) + t^3(q^-22+q^-20) + t^4(q^-20+2q^-18+q^-16) + t^5(q^-26+2q^-24+q^-22) + t^6q^-22 + t^7(q^-28+q^-26) + t^8q^-32 (19) | t^0q^-12 + t^2q^-16 + t^3q^-20 + t^4(q^-20+q^-18) + t^5(q^-24+q^-22) + t^6q^-24 + t^7(q^-28+q^-26) + t^8q^-30 (11) |

### R.3 Wang k = 3 (in addition to `tests/data/wang_ref.txt`)

Full polynomials are in `proto_spec/wang_k3_extra.txt`:
- T(2,3) at N = 7: rank 307;
- T(2,5) at N = 5: rank 106;
- T(2,5) at N = 6: rank 396.

Already in `wang_ref.txt`: T(2,3) at N = 3..6 (ranks 1, 10, 46, 136); T(2,5) at N = 4 (16).

### R.4 Hopf formula

`hopf_formula(N, a, b)` of `tests/xengine.rs` (R3 V4), for the 13 cases of T3.
