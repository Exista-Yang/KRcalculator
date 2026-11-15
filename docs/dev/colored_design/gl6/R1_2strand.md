# R1: fastest certified route to 3-colored gl_6 homology of 2-strand braids

Scope: wedge^3 on both strands, labels (3,3), words sigma^n (T(2,n) knots for odd n, T(2,2m) links,
Hopf = sigma^2, mirrors), N = 6. Sources: dev 3c3323b (xengine, sbim, schur, cclose, cube), colored-m4b
59eaea2 (homprov, smap, tmpl), the laptop G1(a) runs (`colored_design/g1a_33_try.json` warm cache and
`g1a.json` cold cache), the M6.0 logs (`m6_run33k.log`, `m6_run32.log`), and SPEC_M4_M6a section 6.6.
Unknown counts in section 4 come from `gl6/unk.py` (sbim::hom_space conventions, about 1 s CPU).
No cargo was run.

## 0. Verdict

Recommended route: option (i), plus one change in xengine. Use `homprov` as the source of the
decomposition maps, convert its DAG maps to exact `PMat`s, and keep xengine's exact engine and
exact closure. On top of that, add a **rung-by-rung transport** of d_C (x) id in
`XEngine::tensor_crossing` (option iv-a, which is needed in any case at (3,3) fibre sizes).

- The result is exact end to end. Every map xengine uses is re-verified exactly, so correctness
  does not depend on homprov's point-based certification. That certification is still WIP and
  unreviewed.
- The decomposition data that took hours to days with `sbim::hom_space` takes about 1-2 minutes
  with homprov.
- Expected cost on the cloud:
  - first run of T(2,3): 15-60 min, dominated by one-time cached data;
  - T(2,9): 0.5-1.5 h cold, 15-45 min warm;
  - one scan of sigma^10 yields Hopf, T(2,3) through T(2,10) and every N at once.

M6a (option iii) remains the 3-strand engine. The exact 2-strand results become its oracle.

## 1. The 2-strand (3,3) problem is a finite, fixed set of data

- **Closable types.** The block ((3,3),(3,3)) has 4 types, t_k = E0(k)F0(k), k = 0..3, with ranks
  1, 12, 30, 20 (`m6_run33k.log`, the `types` line). The Rickard terms are W_k = E0(k)F0(k) = L_{t_k}
  for **both** signs: `schur::crossing_terms` with c = d gives [E(k), F(k)] in each case. Only the
  differentials differ between the signs. There are 3 per sign, and `Rickard term k ... indecomposable`
  is asserted in `XEngine::local`.
- **Decompositions.** Rung by rung, `XEngine::decompose_word` (xengine.rs:606) decomposes L_p (x) E0(k)
  into canonical ladders of the blocks (4,2), (5,1), (6,0), and then L_r (x) F0(k) back to (3,3).
  - These are exactly the G1 (3,3) pairs: 18 pairs, 6 of them identity shortcuts, 72 (q, s) targets,
    sum of multiplicities 121, largest X = L_p R of fibre 90 (`homprov::g1_targets`).
  - They are cached in `wdecomp` independently of the braid and of N. Every 2-strand (3,3) braid needs
    the same 12 non-identity decompositions.
  - Composite example: t2 (x) W_2 = 6 t2 (+) 36 t3. The first level t2 (x) E0(2) has 18 summands
    (pair 9), and each of its 12 copies of E0(3)F0(1) contributes 3 copies of t3 (pair 12):
    6*30 + 36*20 = 900.
- **N.** At N = 6 nothing is pruned on 2 strands. The unpruned xengine scan is N-independent
  (`XScan::close` takes any list of N), so one scan gives N = 3..8 simultaneously. This includes
  N = 5, the target of cloud run G, and N = 4, already validated by tests/ntrunc.rs.
- **Sizes** (SPEC section 6.6; T(2,9) extrapolated about linearly):

  | | Hopf | T(2,3) | T(2,5) | T(2,9) |
  |---|---|---|---|---|
  | LB final | 13 | 24 | 42 | about 78 |
  | peak objects | 189 | 802 | 2 154 | about 5 000 |

  Gram weight is at most 9 (158 unknowns), and the closure at N = 6 is small.

## 2. What xengine needs from hom_space (exactly)

`hom_space` is called only from `XEngine::hom_norm` (xengine.rs:553). `hom_norm` is called only from
`decompose_direct` (xengine.rs:670), and only for X = L_p (x) R with a single rung R in the default
`rungwise` mode:

1. `iotas = hom_norm(L_q, X, s)` and `pis = hom_norm(X, L_q, -s)`: k-bases of the full graded pieces,
   returned as `PMat<F>`. These are exact polynomial matrices over R_lambda (bottom e-variables) in the
   SBim fibre basis built by `SBim::identity(2K, lam).add_rungs(rungs)`, unnormalized degree.
2. Intersection form G[j][i] = `scalar_of(pi_j . iota_i)`, which must be a scalar
   (End^0(L_q) = k). Then `dense::rref` selects `cnt = [v^s] m_q` columns and rows, and `invert`
   rescales pi.
3. Across all summands: Pi.Iota = I + nilpotent, inverted by a Neumann series. Asserted: unipotent,
   then Pi.Iota = id exactly. Also asserted: Σ summand ranks (+ killed ranks) = rank X, and
   multiplicities = `schur::term_mult`.
4. Downstream, only the `Summand { q, s, iota: PMat, pi: PMat }` values are used:
   - composites in `decompose_word` (kron_id, matmul);
   - `dloc`, i.e. pi (id_{L_p} (x) d_j) iota with `id_tensor_local`;
   - `tensor_crossing` for d_C (x) id: `f.kron_id(fs_k)`, sandwiched;
   - `simplify`: scalar pivots, `h.matmul(g)` updates;
   - `d2_components`;
   - `closure_data`: `cclose::rlambda_bases` per final type and `cclose::solve_f` per edge.

   None of these needs a full Hom basis. **xengine only needs `mult` maps iota and `mult` maps pi per
   (q, s), with an invertible form, as exact PMats.** This is precisely the contract of
   `HomProvider::pairs` (homprov.rs:454).
5. Exactness of the decomposition no longer depends on the provider. Suppose the PMats are bimodule
   maps, and steps 2-3 hold exactly (Pi.Iota = id and Σ ranks = rank X).
   - Then X = (+) summands (+) C, where C is the image of id - Σ iota pi.
   - C is a graded direct summand of a free R_lambda-module and has rank 0, so C = 0.
   - The EKLP formula enters only as the multiplicity cross-check.

Other costs in xengine at (3,3), measured:

- **Rickard maps.** `cube::local_cx(3,3,N=6)`, mostly the Hom^1 W_2 -> W_3 (30 -> 20) inside it,
  takes 125-141 s per sign, with 0 decomposition Hom calls (`m6_run33k.log`, `m6_test33k.log`).
  `tmpl::TmplLib::rickard` (tmpl.rs:462) already wraps `local_cx`, caches it on disk, and asserts
  uniqueness of d_j.
- **Closure data.** All 4 types took 3.3 s.
- **Transport and elimination.** Small at M6.0 sizes: (3,2) sigma^2 took 45.64 s in total, of which
  45.61 s were hom_space, for 65 objects with X <= 36. At (3,3) the composite X reaches fibre 900.
  See section 3.3.

## 3. Option (i): homprov pairs -> exact PMats -> xengine

### 3.1 Conversion is mechanical; xengine does not need to go localized

homprov's maps are point-independent DAGs in `smap::Arena`.

- **Node kinds.** Every node homprov creates goes through `ar.add` with one of
  `MapNode::{Comp, Dot, Id, Lin, Local}` (grep of homprov.rs: 41 Comp, 11 Local, 8 Id, 2 Lin, 1 Dot).
  There are **no `Elem` closed forms and no `Kron` nodes** in its outputs.
- **Templates.** Every `Local` leaf references a `tmpl::Template` that carries an exact `phi: PMat<F>`.
  - It is either a `sbim::hom_space` window kernel element (`TmplSource::HomSpace`),
  - or an `Added` window projection or inclusion (homprov.rs:1668-1676), certified by
    `certify_template` (exact commutation with every top operator plus polynomiality certificate).

Exact evaluation `exact_pmat(ar, m) -> PMat<F>` (memoized per MapId, about 80 lines, in homprov.rs or
smap.rs):

- `Id(e)`: `PMat::identity(flows(e))`.
- `Local{src, at, first, tmpl}`: `tmpl::tensor_exact(prefix_sbim(src.rungs[..at]), first, template,
  suffix_fsize)` (tmpl.rs:594). Here suffix_fsize = `schur::flows(level after the window, suffix)`.
  Test A5 in tests/m4a.rs:631 already checks `Local` against `tensor_exact` at points.
- `Dot{e, level, upright, nu}`: `homprov::dot_pmat` (homprov.rs:2149).
- `Comp(g, f)`: `exact_pmat(g).matmul(exact_pmat(f))`.
- `Lin`: Σ c_i M_i.
- `Kron(f, rs)`: `exact_pmat(f).kron_id(flows(rs))` (for completeness).
- `Elem`: refuse. If ever needed, replace it by the 1-dimensional `hom_basis` window template with the
  F8 constant 1.

Fibre bases agree:

- `Expr::sbim` (sexpr.rs:95) is `SBim::identity(2K, lab).add_rungs(rungs)`, the same construction as
  `XEngine::base.add_rungs`.
- homprov words are `schur` `ladder_minflow`. xengine's `ladder()` picks the minimal-flow, fewest-rung,
  lexicographic ladder among those with at most 2 rungs. They agree on the 7 (3,3) types that were
  printed (`m6_run33k.log`, G1 words).
- Assert `same_sbim(expr.sbim(), lad.b)` per type, or switch `ladder()` to `ladder_minflow`. Nothing
  forces xengine into the localized representation.

Shapes: the largest converted maps are 90 x 30 (X of one rung step against L_q). The 242 maps (121
summands, iota and pi each) are sparse. A cached (2,2,2) 10 x 36 template has 444 terms, about 1.2 per
entry.

### 3.2 What to build

1. `exact_pmat` (above), with tests:
   - localization of `exact_pmat` equals `Eval::all` at 2 points for every certified G1 (2,2) and
     (3,3) target;
   - exact top-operator commutation `T' Phi = Phi T` for all 6 top operators.
2. `XOpts::hom: HomSource { HomSpace, Provider { cache_dir } }`. In `decompose_direct`, for each kept
   (q, s, cnt), call `Search` (or `Chain::for_labels(3, Exact::default())`) `.pairs(ar, schur, x, q,
   lq, None, s, cnt, &Budget::default())`, then `exact_pmat` on every iota and pi, then the existing
   form/rref/Neumann code unchanged.
   - Fall back to `hom_norm` on `Missing`. This never happens for (3,3).
   - Also assert exactly: the form equals homprov's `Pairs::form` at the same scalars (cheap
     cross-check), and every converted map is `homogeneous` and commutes with the top operators.
   - Arena and TmplLib are per lambda, owned by XEngine.
3. Rickard maps from `TmplLib::rickard` (disk cache) instead of `cube::local_cx` in `XEngine::local`
   when `prune_n` is None or N >= c + d. The truncated path is unchanged.
4. Rung-by-rung transport (section 3.3) and disk caching of `dloc` (per type p, j, sign, field).
5. `XOpts::close_stages`: closure data after every crossing. The scan of sigma^10 then closes
   sigma^2..sigma^10, giving Hopf, T(2,3), the T(2,4) link, ... in one process. Closure data costs
   `rlambda_bases` once per type plus `solve_f` per edge.

Integration is a branch from dev 3c3323b merged with colored-m4b. The branches touch disjoint files:
trunc-n changed xengine, schur and tests/ntrunc; M4b changed smap, tmpl, homprov and tests/m4a. About
350-450 new lines, one cloud session including the gates of section 8.

### 3.3 The real second bottleneck: exact transport through fibre-900 composites

`tensor_crossing` forms `phi = f.kron_id(fs_k)` on all of X = L_p (x) W_k (up to 900 x 900). It then
computes `tu.pi.matmul(phi.matmul(sg.iota))` with the composite 900 x 20-30 iota and pi of every
summand pair. For a t2 -> t2 edge and k = 2 that is 42 x 42 sandwiches at dimension 900.

Estimate: 10-100 ms each with `Poly::mul`, which goes through a HashMap. That is 20-180 s per
(edge, k), and about 1-7 h per crossing at T(2,9) scale (around 150 edges, 3 non-trivial terms).

The composites in `decompose_word` are exactly iota = (iota_1 (x) id) . iota_2 and
pi = pi_2 . (pi_1 (x) id), with no re-orthogonalization. Hence

    pi_tau (f (x) id_{R1 R2}) iota_sigma = pi_{2,tau} ( g_{r r'} (x) id_{R2} ) iota_{2,sigma},
    g_{r r'} = pi_{1,r'} (f (x) id_{R1}) iota_{1,r}.

Every product then has dimension at most 90:

- 18 maps a_r = (f (x) id) iota_{1,r};
- at most 18^2 maps pi_{1,r'} a_r of size 6 x 90 x 6;
- then, for each nonzero g, sandwiches of dimension at most 60.

Estimate: 0.1-1 s per (edge, k), i.e. about 1-8 min per crossing at T(2,9) scale.

The id (x) d_j part (`dloc`) needs the 900-dimensional composites, but only once per (type, j, sign).
That is 12 combinations for positive braids, estimated 10-60 min in total on the cloud, disk-cached.

`simplify` rescans quadratically, which is still fine here: about 2.5e3 pivots x 1e4 components. If
transport is still too slow, the fallback is R-linear transport tables: T(c h) = c T(h) for c in
R_bottom. Precompute T on R-generators of Hom(B_p, B_p') (`homprov::hom_gens`, coordinates by
`solve_f`-type solves). About 3 000 transports one-time, then polynomial combinations.

### 3.4 Certification, risk, cost

- **Certification.** Fully exact:
  - converted maps checked as bimodule maps;
  - decompositions by section 2, item 5;
  - K_0 equals the decategorified scan, objects >= LB, d.d = 0 exactly, degrees checked;
  - exact closure data (determinant criterion plus c_diag degrees, `solve_f`) and
    `close_exact_colored`;
  - both primes F31 and F32b;
  - independent references: Wang (`tests/data/wang_ref.txt`, k N m = 3 6 m, and `r3scripts/wang.py`
    for more m), the Hopf formula, `decat_chi`, the N = 3, 4 pruned results, and duality at N = 4.

  homprov's stacked-point certificates are used only to find the maps.
- **Risks.**
  - Medium: transport speed. Section 3.3 is an estimate; measure it on Hopf and T(2,3) first.
  - Low: Search non-determinism or exhaustion on a cold cache. All 72 (3,3) targets certified cold;
    seed and cache are fixed.
  - Low: branch integration. colored-m4b has not compiled since the WIP commit; the cloud session
    finishing M4b fixes that first.
- **Measured one-time costs (laptop).**
  - homprov (3,3) decompositions: all 12 pairs plus the 4 closure types, cold cache: job walls
    0.26-22.6 s, sum about 69 s, peak RSS 124 MB. Warm cache: 8.5 s for all 16 jobs. Tiers
    S1/S2/S3, one `Exact` window (0.011 s).
  - Same pair under hom_space: (2,1/1,2) (x) E0(1) took 47-48 s in xengine, 0.26 s with homprov.
    The 11 769-unknown pair 5 took 5 401 s under hom_space, 0.26 s with homprov.
- **Cloud estimates (about 1.5x slower).**

  | step | cold | warm |
  |---|---|---|
  | homprov | about 2 min | under 15 s |
  | Rickard, positive sign | about 3 min | 0 (disk) |
  | `exact_pmat` and composites | about 1-3 min | |
  | `dloc` | 10-60 min | 0 (disk) |
  | each crossing of T(2,3) | under 1 min | |
  | each late crossing of T(2,9) | 1-8 min | |
  | closure data and closure, N = 3..8 | under 2 min | |

  This gives:

  | target | cold | warm |
  |---|---|---|
  | T(2,3) at N = 6 | 15-60 min | under 5 min |
  | T(2,9) | 0.5-1.5 h | 15-45 min |
  | sigma^10 with `close_stages` | about the same as T(2,9) | |

  Memory stays under 2 GB.
- **Helps 3 strands?**
  - `exact_pmat` is the exact bridge M5 and M6a need as an oracle for localized maps.
  - The scan exercises homprov's (3,3) pairs inside a real engine (early G2 evidence).
  - The exact (3,3) results are M6a.2's reference values (SPEC 6.1: M6a must reproduce M6.0).
  - The rungwise transport mirrors M6a's "rungwise colored Rickard insertion".
  - xengine itself stays 2-strand only (`check_input`, two-rung canonical ladders; exact PMats are
    infeasible on 3 strands, SPEC F2/F5).

## 4. Option (ii): speed up sbim::hom_space

**Anatomy** (sbim.rs:202-310):
- The unknowns are all coefficients of all fibre entries: Σ_{beta, beta'} #mono(d + fdeg[beta] -
  fdeg'[beta']) over R_lambda (6 generators, degrees 2, 4, 6 twice).
- There is one equation per (top operator: 6, beta', beta, monomial) of T' Phi - Phi T = 0.
- `linsys::Echelon` (sparse rows in a HashMap, restart sweeps in `reduce`, full `rref` for the
  kernel) behaves cubically in the unknowns.

**Calibration.** 2 524 unknowns (Hom 15 <-> 36, degree +-1) take about 20 s. The 5 401 s call is
pair 5, Hom between X = E0(2)F0(1)F0(1) (60) and t2 (30) at s = +-1, with 6 521 / 11 769 unknowns.
The cubic model gives about 2 000 s on the laptop, and the measured cloud value is consistent with
it.

**Totals over all 72 targets (iota and pi).** 431 538 unknowns. The largest calls are 13 538-15 033
unknowns (pairs 2, 3, 6, 12: X of fibre 60-90 against types of rank 15-20). The cubic model gives
about 51 000 s (14 h) on the laptop. Scaled by the calibration point, that is about 1-2.5 days on the
cloud, one-time but braid-independent.

| idea | gain | comment |
|---|---|---|
| packed dense or blocked modular elimination | about 2-5x at best | the system is already sparse; a dense solve of 15 000 unknowns is about 1e13 operations |
| only the needed degrees | already done | xengine asks only for s and -s |
| EKLP dimension as stopping certificate | up to about 2x | randomized kernel sampling still needs the elimination |
| known multiplicities: stop at `mult` maps with an invertible form | small | the system must still be solved once |
| restrict through known idempotents | none | the idempotents are the output |
| ladder-flip duality: pi from iota | 2x | |
| extremes plus dots, exact: hom_space only at s_lo / s_hi (1-dim, F4), others as `dot_pmat` multiples, rank by the form | 5-50x per family | extremes still reach 13 538 unknowns (pair 3, iota at s = 3); dots need not span (pair 15 needed S3 and one window) |
| strip common prefix and suffix (windows) | decisive | this is homprov's local decomposition tree, i.e. reimplementing M4b |
| generator images (Hom from a bimodule with g generators: unknowns / fibre(src) x g) | up to (fibre/g)^3 | needs presentations of the canonical ladders; research risk |

**Verdict.**
- Generic speed-ups give 2-10x, leaving 2-24 h one-time.
- The ideas that really help amount to rebuilding homprov.
- As a **fallback only**: persist `decompose_direct` results to disk (PMat text format of
  `tmpl::read_cache_file`) and let the brute-force scan accumulate them over several cloud sessions.
  This is about 50 lines, has zero design risk, and costs 1-2.5 days of cloud compute. Downstream it
  is identical to (i), including section 3.3. Exact.

## 5. Option (iii): the M6a localized engine (sengine) on 2 strands

- **Build.** All of the following, about 3-5 cloud sessions plus reviews (M4a and M4b each took about
  1-2):
  - M5 `sindec.rs`: rungwise decompositions as global maps; e_p = None on 2 strands (F1); Pi.Iota and
    Σ iota pi = id at two points;
  - M6a `sengine.rs` and `sscan`: typed complex per point; `Eval` transport over shared Layouts;
    heap Markowitz with streaming; determinism traces; scans at many points (`engine2::par_map`);
    Gram by interpolation; F by rescans;
  - tests T0-T9.

  (3,3) needs only the passing (3,3) part of G1(a).
- **Run cost** (SPEC 6.6, laptop, two threads, one prime):
  - (3,3) Hopf / T(2,3) / T(2,5): 0.04-1.1 s per point, 166-310 scans, at most 3 min.
  - T(2,9): about 0.9-2.6 s per point (peak Σ S extrapolated to about 1.1e5), 170-300 scans,
    i.e. 3-13 min.
  - Cloud x1.5, both primes x2.
  - No 900-dimensional exact products: point blocks are small.
- **Certification.**
  - Decompositions: form rank at a top-stratum coset, re-verified at a second point.
  - Closure generators: determinant criterion at every scan point.
  - Gram and F: interpolation with 8 residual points. Probabilistic in the Schwartz-Zippel sense;
    the standard of engine2, not exact.
- **Risk.** Medium-high: a new engine, interpolation caps, point-dependence traces.
- **Helps 3 strands.** Yes. It is the 3-strand engine, so it must be built anyway, but first results
  on 2 strands arrive several sessions later than (i), and they would need (i) to be cross-checked
  exactly.

## 6. Option (iv): other observations

- **(a) Rungwise transport.** Section 3.3: required by (i) and (ii), free in (iii).
- **(b) One scan, all N and all n.** The unpruned scan is N-independent. With `close_stages`, sigma^m
  gives every T(2,k), k <= m. So the whole target list (Hopf, T(2,3..9), T(2,4..10), N = 3..8,
  reduced and unreduced) is **one** cloud job per prime.
  - Run G (N = 5, pruned) is superseded by its N = 5 rows. It stays useful only as an independent
    check of the pruned path.
  - F2b is cancelled, as the user asked.
- **(c) Both signs share decompositions** (same words E0(k)F0(k)). Mirrors and negative Hopf cost only
  the negative Rickard maps (about 3 min, cached) and their `dloc`.
- **(d) Closure side.** xengine's `rlambda_bases` already takes 3.3 s for all (3,3) closure types.
  homprov's `closure()` (8 ins, 8 outs, certified in milliseconds) is not needed on 2 strands.
- **(e) Rigid or precomputed 2-strand structure.** Every needed decomposition is an instance of the
  categorified sl_2 relations E^(a)E^(b) = [a+b, a] E^(a+b) and the E^(a)F^(b) <-> F E square switch
  (thick calculus). homprov's window tree finds exactly these as small window templates.
  - Coding them by hand (thick calculus, or explicit colored 2-strand complexes from the literature)
    adds a second source of truth with sign and normalization risk.
  - It gains nothing over (i) for certification; Wang's closed formula is already the independent
    reference.
  - Precomputed rigid data alone (term_mult, Rickard uniqueness, hom_gd) does not fix the components:
    Hom^1 between the rank-30 types has dimension 6.
  - Not recommended.
- **(f) R-linear transport tables.** The section 3.3 fallback; a 2-strand special case of D3's
  Hom-basis coordinates.

## 7. Comparison

| | (i) homprov + xengine (+ rungwise transport) | (ii) faster hom_space / brute force with disk cache | (iii) M6a sengine | (iv-e) explicit 2-strand calculus |
|---|---|---|---|---|
| to build | about 400 lines, 1 session | 50 lines (cache) + 200-600 (solver) | M5 + M6a, 3-5 sessions | research plus code |
| one-time (3,3) data, cloud | about 15-70 min (homprov 2 min; Rickard 3 min; `dloc` 10-60 min) | 1-2.5 days (10x solver: 2-6 h) plus the same `dloc` | inside the scans | ? |
| T(2,3), N = 6 | 15-60 min cold, under 5 min warm | as (i) once cached | 1-3 min per prime | ? |
| T(2,9), N = 6 | 0.5-1.5 h cold, 15-45 min warm | as (i) once cached | 5-20 min per prime | ? |
| certificate | exact | exact | probabilistic interpolation | exact if derived right |
| main risk | transport speed (measure on Hopf) | compute time, cloud sessions | build time, new engine | derivation errors |
| helps 3 strands | oracle and bridge (`exact_pmat`), reference values | no | it is the 3-strand engine | no |

## 8. Recommended plan

1. After the M4b cloud session pushes colored-m4b (G1(a) (3,3) is already green), create branch
   `colored-2s` from dev and merge colored-m4b.
2. Implement `exact_pmat` with its tests (section 3.2, item 1).
3. Add `XOpts::hom = Provider`, the Rickard maps from `TmplLib::rickard`, `ladder()` agreeing with
   `ladder_minflow`, and the exact bimodule checks.
4. Add rungwise transport and the disk cache for `dloc`.
   - Test: rungwise transport gives bit-identical components to direct transport on every existing
     xengine and ntrunc word.
   - Test: `HomSource::Provider` gives identical Poincare polynomials to `HomSpace` on all (2,1),
     (2,2), (3,1), (3,2) tests.
5. Gate A (cloud, 30 min budget): (3,3) sigma^{+-1} and sigma^2.
   - Checks: Hopf formula; N = 3, 4 equal to the pruned scans; N = 5, 6 against the Hopf formula and
     chi; both primes.
   - Print per stage the seconds spent in decomposition, `dloc`, transport, elimination, d^2 and
     closure data.
   - If transport for sigma^2 exceeds 2 min, implement section 3.3's tables before going on.
6. Gate B: sigma^5 with `close_stages`. Wang m = 1, 2 at N = 6 and also N = 3, 4, 5, 7.
7. Production: sigma^10 (and sigma^9), `close_stages`, N = 3..8, reduced and unreduced, F31 and F32b.
   Compare with `r3scripts/wang.py` for k = 3, N = 6, m = 1..4 (rank 136 for the trefoil). Record
   the timings as input to M6a's cost model.
8. Then M5 and M6a for 3 strands, with these exact (3,3) results and `exact_pmat` as oracles.

**Open measurements** (none can be taken locally):
- the transport time per (edge, k) at fibre 900 versus rungwise;
- `dloc` one-time cost;
- whether Search reproduces the 72 (3,3) certifications on the cloud's cold cache within its budget
  (the laptop needed at most 23 s per pair).
