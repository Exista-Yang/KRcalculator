# M4c A2 (+ A3) report: Lemma A leaves backend, closure list

Date 2026-11-12. Branch `m4c-closure` (from origin/m4c-census 595068c). Scope: M4C_PLAN §4 milestones A2 and A3,
§3.1 (Lemma A), §3.2 (closure), §3.4, §3.6, §5 items 2-4.

**Frozen files.** No edit to smap.rs, homprov.rs, sbim.rs, linsys.rs, tmpl.rs or sindec.rs (including
`provider_for_stamped`: the A3 wiring is a separate constructor in leaves.rs, see §5). src/leafgen.rs was not
created. Changed: src/leaves.rs, src/bin/m4c.rs, src/bin/libbuild.rs, tests/leaves.rs; new data in
docs/colored_design/m4c/a2/.

## 1. What was built (src/leaves.rs)

- **`Leaves::build_leaves(x, targets, Out|In)`**: Lemma A (D1 §1.3) level by level on 3 uprights.
  - Junction types per level are `dcoset::junction_types(L_p, targets)`.
  - **Junction words are planned first** (`plan_words`, pure combinatorics): the default rule is
    `dcoset::junction_word` (prefix extension, then canonical ladder, then BFS; ∅ for the final identity type).
    If some step exceeds the 10⁵ budget, a coordinate descent runs over the reduced words (≤ 5 rungs) of each
    junction type. Its objective is the estimated unknowns (`step_score`) of the steps into and out of the type
    (minimize the max, then flows; at most 4 sweeps). This is new relative to A1: A1 varied only the target word.
    Varying the *source* word solves 3 over-budget steps (§4).
  - Out: LL_t = { S ∘ Kron(ℓ, [ρ]) }. In: { Kron(ℓ, [ρ]) ∘ S }, with S ⊂ Hom(W_t, W_{t'}·ρ).
  - Diagnostic only: every family has exactly c_t(X)(1) maps and an invertible block at rep(t) at one point.
    Otherwise the result is Exhausted.
- **`step_family`**: S_{t'→t} with degree multiset G_{t'→t} and an invertible top-to-top matrix at rep(t).
  - Columns: the flows of W_{t'}·ρ at rep(t) whose configuration after |W_{t'}| rungs has type t'.
  - Order:
    - `Id` when W_t = W_{t'}·ρ;
    - (opt-in) an E2/E3 `Elem` for a one-dimensional cup/cap step;
    - exact `hom_space` bases on the differing window (common prefix/suffix stripped);
    - exact bases on the whole words.
  - Selection: greedy by ascending degree over every degree of supp G (C_m4c m4), capped at the coefficient of
    each degree. A family is accepted only with multiset exactly G and full top-to-top rank. The window step
    falls back to the whole words.
- **Budget rule** (M4C_PLAN §0.3, in `Leaves::admit`):
  - ≤ 10⁵ unknowns: admitted;
  - 10⁵ to 4·10⁵: only with an uncensored exact-shape cost record (`CostTable::predict_exact` finite);
  - above 4·10⁵: Refused;
  - new: one-sided calls (an empty end) on windows with a label > 3 need a record above 2·10⁴ unknowns (§4).
- **`closure_leaves`** (the `HomProvider::closure` of `Leaves`):
  - `closure_precheck`, then Lemma A to the identity type with final word ∅, then `certify_closure`, inside a
    mark/commit;
  - `ep = Some` → Refused;
  - `full`, `pairs`, `end0` → Refused (A4/A5).
- **`validate_elem(kind, labels, window)`**:
  - the closed form with scalar 1 is compared with the exact `hom_space` map, which must be the unique one up
    to a scalar (Hom^deg one-dimensional);
  - the comparison is blockwise proportional at every coset at 2 points, with a nonzero ratio;
  - refused above 1.2·10⁵ unknowns;
  - `elem_pattern` maps a 2-upright step window to Eps/EpsL/Trace/TraceL (Out) or Eta/EtaL/Iota/IotaL (In).
- **Re-implemented private helpers:** `identity_coset`, `rep_rank_of`, and `cdiag_degrees` (A1).
- **Wiring (A3):** `provider_for_stamped_m4c(lam, costs, slowdown, leaves_first)` = `ClosureLeaves`, which wraps:
  - Stage 0 `ProductClosure` (product-first) around the unchanged production chain (`sindec::provider_for_stamped`);
  - `Leaves` for closure calls without e_p, after the chain (`m4c`, M4C_PLAN "Leaves last") or before it
    (`m4c-lf`);
  - every other call goes to the inner chain alone, so pair failure reports are unchanged.

  Stamp: chain stamp + `stage0=ProductClosure(v1,product-first)` + `LEAVES_STAMP(leaves-first|last)`.
- **Tools:**
  - `m4c leaves --field 31|32b --types fails|all|product|T1;T2 [--direct] [--cache DIR] [--costs PATH] [--replay] [--from I]`:
    `--direct` prints per family the steps, exact calls, unknowns and G-L5 counters; otherwise it runs SIndec
    with layer A and the Leaves stamp.
  - `m4c measure --base IN --out OUT CALL...`: in-process exact calls written as uncensored exact-shape cost
    records.
  - `libbuild ... --provider leaves|m4c|m4c-lf`.

### Finding: Elem constants are not persisted
`Elem` constants live in the arena (`certify_elem`), not in layer A. A layer-A file containing Elem nodes is
rejected on load in a fresh session as "not genuine" (replay: 4/12, seen in this session). That is correct
behaviour: nothing wrong was accepted. Elem steps are therefore **opt-in** (`use_elem`, off by default), and
production layer A holds exact templates only. The unit tests run with Elem on, so V0 and G1(a) also exercise
the Elem path. Persisting validated constants (or re-validating them at SIndec start) is an open item (§7).

**Full suite:** `KRHOM_THREADS=4 cargo test --release` passes on the final code: 224 passed, 0 failed (about
5 min).

## 2. Gates (tests/leaves.rs, both primes where stated)

| gate | result |
|---|---|
| **V0** vs `sbim::hom_space`, λ ∈ {(2) as (2,0), (1,1), (2,1), (1,2,1)} (every word with top λ of ≤ 4/4/3/3 rungs, Out and In) | **102 (word, direction) cases pass.** Leaves have degree multiset c_id(char X). In each of 5 (resp. 3) consecutive degrees, the stacked rank of leaves·(R^λ Schur basis) = dim Hom^d (exact basis) = rank of their union, so the span is equal (F31) |
| **SLL Ex 7.39** (λ = (2), [s,∅,s] = F0.1 E0.1 → [s]) | 2 leaves of degrees −1, +1, Out and In, `certify_closure` |
| **SLL Ex 7.40/7.41** (I = (2,1), J = (1,2)) | Every word (2,1) → (1,2) of ≤ 5 rungs, plus (2,1,0) → (1,2,0) of ≤ 4 rungs, to every coset, Out and In: count c_t(1), degree multiset c_t, invertible block at rep(t) at 2 points, genuine. 576 + 182 families of sizes 1-11, including words with 4 and 5 leaves in all. The paper's exact expressions are not in the repository; the claim is the family property on all words with these I, J |
| `validate_elem` | 36 closed forms (labels ≤ 3, all 8 kinds) validated for F31 and F32b; a mismatched window is rejected |
| **G1(a) closure part, Leaves alone** | **(3,3): 8/8 F31, 8/8 F32b; (2,2,2): 42/42 F31, 42/42 F32b** (< 1 s; step families (2,2,2): 72 Id, 26 Elem, 6 Window, 54 Whole) |
| A3 wiring | `m4c_chain_closure_222`: all (2,2,2) closure types through `provider_for_stamped_m4c`, leaves first and last |

## 3. A3: the closure list with the Leaves backend alone

**List.** M5c_report (colored-m5) still lists 9 closure types Exhausted at 2x at N = 6, plus 1,1,1/1,1,1/1,1,1 at
N = 5. The 3 product types it certified in a resume pass remain in the A1/M5b list of 12 + 1. All 12 distinct types
were run (1,1,1/1,1,1/1,1,1 serves both N = 5 and N = 6: same type, same canonical ladder).

**Result: 12/12 types certified for F31 and F32b, both sides, with Leaves alone; replay 12/12 for each prime from
layer A** (SIndec with prune 6, stamp `homprov[hs-kl1] maxl=3 leaves=...`, provider refusing every call in the
replay).

| type | n | degrees | wall F31 (ins+outs, direct) | steps (Id/Window/Whole) | max step unknowns | DAG nodes | index bytes | eval (all cosets, 1 point) |
|---|---|---|---|---|---|---|---|---|
| 0,0,3/3,0,0/0,3,0 | 1 | 18 | < 0.1 s | 2/0/2 | 3 184 | 9 | 20 063 | < 0.1 ms |
| 0,3,0/0,0,3/3,0,0 | 1 | 18 | 0.7 s | 2/0/2 | 3 184 | 9 | 20 608 | 0.2 ms |
| 1,0,2/2,0,1/0,3,0 | 3 | 13-17 | < 0.1 s | 2/0/2 | 2 852 | 19 | 26 639 | 0.2 ms |
| 1,2,0/0,0,3/2,1,0 | 3 | 13-17 | 0.4 s | 1/1/2 | 2 983 | 15 | 28 421 | 0.3 ms |
| 0,3,0/1,0,2/2,0,1 | 3 | 13-17 | 0.9 s | 1/1/2 | 5 104 | 19 | 27 133 | 0.1 ms |
| 1,0,2/1,2,0/1,1,1 | 18 | 9-19 | 0.3 s (2 s cold) | 2/0/3 | 13 765 | 47 | 135 335 | 1.3 ms |
| 1,1,1/0,2,1/2,0,1 | 18 | 9-19 | 0.9 s | 2/1/2 | 2 852 (window 2 401; whole 101 507) | 48 | 135 825 | 3.0 ms |
| 0,2,1/2,0,1/1,1,1 | 18 | 9-19 | 0.1 s | 3/2/3 | 2 852 | 65 | 169 155 | 1.4 ms |
| 1,1,1/1,0,2/1,2,0 | 18 | 9-19 | 0.2 s | 4/2/3 | 1 220 | 68 | 189 607 | 0.7 ms |
| 0,1,2/2,0,1/1,2,0 | 27 | 10-22 | 1.0 s | 3/2/4 | 9 504 | 95 | 456 503 | 1.3 ms |
| 0,2,1/1,0,2/2,1,0 | 27 | 10-22 | 167 s cold (1 exact call 41-60 s per side + templates) | 4/0/4 (+1 Elem in the direct run) | **111 497** (measured record) | 92 | 403 725 | 23 ms |
| 1,1,1/1,1,1/1,1,1 | 54 | 4-18 | 1.7 s | 4/2/3 | 7 225 (window; whole 88 301) | 130 | 277 523 | 5.7 ms |

Notes:
- The degree multisets are c_diag(b_p), enforced by `certify_closure`. They therefore equal those of any other
  certified closure family, including M5c's layer A for the types it has.
- Counts are per side. Wall times are the direct runs (warm template cache except where marked). SIndec runs of
  all 12 types: 167 s (F31) and 131 s (F32b) cold, dominated by the 111 497 call; replays take 0.8-0.9 s.
- The 111 497-unknown step (0,2,1/1,0,2/2,1,0, level 4, Hom^{3..9}(E0.1 F0.1 E1.3 E0.2, E0.2 E1.3), the A1
  open problem) was measured on this machine:
  - In direction: 46.7 s, 666 MB;
  - Out direction: 65.7 s, 739 MB;
  - Hom^9 has dimension 74.

  These are uncensored exact-shape records, so the call falls into the 10⁵-4·10⁵ band of the budget rule. The
  101 507 step of 1,1,1/0,2,1/2,0,1 was avoided by the differing window (2 401 unknowns), as A1 §5 (a) proposed.
- **Persistence.** docs/colored_design/m4c/a2/data/f31 and f32b hold layer A (100 KB per prime) plus the template
  cache (13 MB per prime, 9.8 MB of it the 74-dimensional Hom^9 basis). They are stored as `data_f31.tar.gz` /
  `data_f32b.tar.gz` (2.4 MB each), with the run and replay logs as text. The N = 5 job (prune 5, the same
  type) was run and replayed in the same cache: 1/1 for each prime. The cost records are in
  docs/colored_design/m4c/a2/costs_333_m4c_a2.txt (costs_333_m5c.txt + 14 records). This is not merged into
  data-gl6: that needs the coordinated stamp change (§5).

## 4. A3: (3,3,3) closure regression, every type of (λ, λ), Leaves alone

`m4c leaves --types all --direct` covers all 55 types of the (3,3,3) block, i.e. every level-4/5 closure job and
more. Cost table: a2.

| outcome | F31 | F32b |
|---|---|---|
| **certified, both sides** | **51 / 55 types (102 / 110 sides)** | **51 / 55 types (102 / 110 sides)** |
| refused: decomposable 0,1,2/1,1,1/2,1,0 (needs e_p; Leaves refuses e_p by design) | 1 | 1 |
| refused: 0,1,2/2,1,0/1,1,1, a step of **434 787** unknowns (Hom^13(E1.2, E0.3 E1.2 E0.1 F0.4)), above the hard cap of 4·10⁵, although it was measured uncensored here (179 s / 272 s, 2.3 / 3.2 GB). No word choice of ≤ 5 rungs avoids it; the descent removed the d = 11 step of 273 477 (measured: 86 / 168 s) but not this one | 1 | 1 |
| refused: 0,1,2/0,2,1/3,0,0, a one-sided step Hom^d(E0.3 E1.3, F0.2 E1.3 E0.5). d = 13 (139 371) and d = 15 (213 148) were measured in this session (440-500 s and 1.0-2.0 GB each); d = 17 (319 734) is not measured yet | 1 | 1 |
| refused: 0,0,3/0,3,0/3,0,0 (w_0, n = 1). Only one-sided calls on (3,6) remain: Hom^18(∅, F0.3 E0.3) at 42 872 unknowns ran > 25 min and 3.3 GB and was killed; the whole words give 233 282 | 1 | 1 |

Run details:
- F31 ran types 0-53 with the binary of commit 77329f3 (Elem on, no one-sided rule) and type 54 with the final
  binary.
- F32b ran all types with the commit 8e70c0c binary. Its 3 types refused by that binary's one-sided cap of 2·10⁴
  (0,2,1/0,1,2/3,0,0, 1,0,2/0,2,1/2,1,0 and 0,0,3/1,2,0/2,1,0, all at labels ≤ 3) were rerun with the final rule
  and certified.
- Every certified family passed `certify_closure`. Total F31 wall time for the certified types: 547 s.

These types are all at level ≤ 5, and the production chain certifies them (M5c). In the A3 wiring Leaves only
backs the chain up, so they are not regressions of the production library.

**Cost observation.** One-sided `hom_space` calls (an empty end) on labels (3,6) are far slower per unknown than
on labels ≤ 3. Other one-sided calls measured in this session:
- 72 117 / 83 641 unknowns: ≤ 66 s;
- 139 371 unknowns (Hom^13(E0.3 E1.3, F0.2 E1.3 E0.5)): 494 s, 958 MB;
- 273 477 unknowns: 86 s / 168 s, 1.4 / 2.1 GB.

Unknown counts alone are therefore not a reliable time budget; the per-call cost records are.

## 5. Wiring into the production chain (A3)

- **Not done in sindec.rs.** `provider_for_stamped` is unchanged, which keeps M5c's layer A and stamps valid until
  the single coordinated rebuild (M4C_PLAN §6, R10).
- **Ready instead:** `leaves::provider_for_stamped_m4c` and `libbuild --provider m4c | m4c-lf`. The one-line
  change in sindec (return `provider_for_stamped_m4c(...)` for max label 3) is left for that merge. The plan
  allows it in A3 only "if required"; it is not required to certify the list.
- **Order.** Stage 0 runs product-first, avoiding the 3-5 min Exact detour per product type (A1 §4). Leaves runs
  last by default.
- **G-L5 argues for Leaves first on closure calls.** Leaves closure maps have ≤ 130 DAG nodes and ≤ 0.46 MB of
  resident FlowIndex per family. A full localization of a whole family at one point takes ≤ 23 ms. The 24 new
  families together hold about 3.9 MB of index and take < 0.05 s per point, against the R4_m6a budget of
  t_lib ≈ 10 s per point and 2 GB for library DAG and templates at prune 6.

  With Leaves first, the chain's 262-600 s Exact closure calls would never run for these types. Leaves last (the
  plan's C M4 rule) costs only wall time in the library build. Recommendation: `m4c-lf` for closure.
  `ClosureLeaves` sends only closure calls to Leaves, so pair maps are unaffected either way.

## 6. Commits (branch m4c-closure)

| commit | content |
|---|---|
| 996aa57 | Leaves backend: `build_leaves` (Out/In), `step_family`, closure; V0 and SLL example tests |
| 3507914 | `m4c leaves` / `m4c measure`, `libbuild --provider leaves` |
| 6fcff9b | `validate_elem`, Elem step families, tests |
| 77329f3 | planned junction words with cost-aware descent; G-L5 counters; cost records of the 111 497 step |
| c4688fa | `ClosureLeaves` wiring (`provider_for_stamped_m4c`, `--provider m4c / m4c-lf`); records of the 273 477 step |
| 8e70c0c | Elem opt-in (layer A self-contained); one-sided admission; `--from` and progress |
| 2e69309 | one-sided admission only for labels > 3 (measured) |
| a1602d0 | layer A + templates of the 12 types, both primes, run/replay logs (N = 6 and N = 5) |
| 1479aa0 | this report; cost records of the 139 371 / 213 148 / 434 787 steps |

## 7. Open problems

1. **0,1,2/2,1,0/1,1,1 (level ≤ 5) with Leaves alone** needs a 434 787-unknown step, which is above the hard cap.
   The step is measured (≤ 272 s, 3.2 GB), so raising the hard cap to 4.5·10⁵ for recorded shapes would admit
   it; that is a plan decision. Other options:
   - descent over words of 6 rungs;
   - a 4-upright (padded) K window;
   - a B1 closed-form step.

   The production chain certifies this type.
2. **w_0 = 0,0,3/0,3,0/3,0,0 with Leaves alone:** one-sided (3,6) windows are pathological for `hom_space`. The
   E2 `Eta` closed form is exactly this step, but validating it also needs the 42 872-unknown call. Options: use
   `poly_certificate` for validation (M4C_PLAN §3.4), or a stored validated constant (item 3). The chain and
   Stage 0 certify this type.
3. **Elem constants in layer A:** persist validated constants (smapio is not frozen), or re-validate them when
   SIndec starts. Until then `use_elem` stays off in production.
4. **e_p closure** (decomposable L_p, 0,1,2/1,1,1/2,1,0): Leaves refuses it. Lemma A applies to the image of
   e_p; this is not needed for the failure list.
5. **Costs and admission:** step admission uses unknowns plus exact-shape records. A one-sided/two-sided split
   of the budget (as M5c did for buckets) would be more faithful (§4). New records are in a2/costs_333_m4c_a2.txt
   and should be merged into the shared cost table with the coordinated stamp change.
6. **0,1,2/0,2,1/3,0,0 with Leaves alone** needs the d = 17 one-sided record (319 734 unknowns; about 10 min
   per direction, < 3 GB expected from the d = 13/15 records).
7. **Next milestones:** A4 (`pairs` via Lemma A with W_q = L_q) and A5 (double leaves). The G1(a) pair part and B2
   gates are untouched here.

## 8. Reproduce
```
cargo build --release --bin m4c
KRHOM_THREADS=4 cargo test --release --test leaves         # V0, Ex 7.39-7.41, validate_elem, G1(a) closure, wiring
C=docs/colored_design/m4c/a2/costs_333_m4c_a2.txt
./target/release/m4c leaves --field 31 --types fails --cache D31 --costs $C            # layer A, 12/12
./target/release/m4c leaves --field 31 --types fails --cache D31 --replay              # 12/12 from layer A
./target/release/m4c leaves --field 31 --types all --direct --cache R31 --costs $C     # regression + G-L5
(same for --field 32b)
./target/release/m4c measure --base tests/data/homprov/costs_333_m5c.txt --out $C "3,3,3;E0.2,E1.3;E0.1,F0.1,E1.3,E0.2;9"
```
