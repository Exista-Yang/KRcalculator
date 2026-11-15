ADVERSARIAL VERIFICATION of "BASE FIELDS IN krhom" (worktree /Users/existair/Documents/KRcalculator/.claude/worktrees/m5i, HEAD 73cbd7c, detached, as the report says). This was read-only: grep/sed/git plus 8 web fetches. Nothing was built or run. Unless a path is given, file:line refers to this worktree.

Bottom line: the code facts are almost all CONFIRMED, apart from a few loose line cites. There are 4 substantive problems:
- (a) The "p > K!" lower bound is a misreading.
- (b) "The repo says nothing about the Hom formula at small p" is wrong.
- (c) The report misses that krhom also silently ignores `--field` values other than q.
- (d) The report mixes one-sided checks (false failures only) with checks that can falsely pass. It also does not say that the kill-rule caveat applies only to pruned runs.

## (1) Modulus representation

- Fp<const P> (field.rs:46-48), aliases (134-141), the two Field impls (Fp at 84; Q struct at 145, impl at 170-205), and F2/F3/F5/F7 used only in the unit test (field.rs:211-216; no other use in src/ or tests/): CONFIRMED.
- PrimeField only for Fp: CONFIRMED (boxclose.rs:59-71; the report says 58).
- krhom F31 + F32b for singular (krhom.rs:126-141): CONFIRMED.
- Omission: krhom only checks whether `--field` equals "q" (krhom.rs:67). So `--field 7` or `--field 32b` silently runs F31 on the cube, Soergel and cable paths, the same defect as sbox.
- sbox silent F31 fallback: CONFIRMED (sbox.rs:206-210). The JSON does record field=F_2147483647 (sbox.rs:197), but there is no error.
- libbuild assert (libbuild.rs:556-557), xscan error (xscan.rs:169-176), m4c and gp0 panic (m4c.rs:596-605, gp0.rs:84-88): CONFIRMED.
- "P < 2^32 required": CONFIRMED for mul (108), fma (130) and pow (60, 62), and for the xengine u32 packing (xengine.rs:301-311). WRONG as supporting evidence:
  - from_i64 (`P as i64`, field.rs:120) only needs P < 2^63.
  - Layer-A signed-integer storage (smapio.rs:25-26, 166-170) is also not a 2^32 limit.
- The xengine Plain fallback at 32 B vs 8 B: CONFIRMED (doc comment xengine.rs:313-319).
- Field's associated functions take no self (field.rs:8-22). So a runtime modulus must be type-level or process-global: CONFIRMED.
- A thread-local fails: CONFIRMED. par_map uses std::thread::scope with fresh threads (sengine.rs:1768-1790).
- Barrett/Montgomery, u128, and the cost of a menu of primes: [reasoning], fine. Nothing is measured.
- Libraries are per prime: CONFIRMED.
  - field= appears in the Layer-A keys (sindec.rs:1701, 1720, 1735), the smapio header (smapio.rs:54) and the template header (tmpl.rs:217).
  - "Cloud-hour" is roughly right: C2_1_report.md:144 projects about 1 h per prime for (2,2,2,2).
- xengine needs no library: CONFIRMED. It is generic and has no file I/O.
- Lower bounds:
  - smap::points: CONFIRMED (smap.rs:200-216, draw at 208). For p < K it loops forever.
    - Extra fact: for small p it can return 0 as a coordinate. Whether that matters is unknown.
    - The doc's no-poles argument assumes "1..2·10^9 < p" (review_m4b/certify.md:219). The code itself checks distinctness after reduction, so it stays sound for any p ≥ K.
  - e_pair (smap.rs:417-431) and its inverse (473), the local.rs:229 panic, and the samplers at xengine.rs:1637, engine2.rs:368 and engine.rs:343: CONFIRMED.
  - BoxNodes assert p > 8·N·#blocks (boxclose.rs:181, draws at 189), and the examples 128 and 144: CONFIRMED.
    - This is a margin assert. The real need is p−1 ≥ N·#blocks distinct nonzero nodes.
  - SCALING_LAMBDA = 0x12345677 = 305419895, divisible by 5 (swire.rs:1132-1133). SCALE_LAMBDA = 1234567 = 127·9721 (sindec.rs:275-276). The review note "no order ≤200" (review_wave3/box.md:52): CONFIRMED.
  - "Weights about 17-24": understated. The docs cite weight 22 (SYNTHESIS.md:344, 744) and entries of weight ≥27 on 3 strands (gl6/R3_closure.md:239).
  - Freivalds (cclose.rs:2913-2920, draw at 2994; k=3 at swire.rs:1128): CONFIRMED.
  - Interp (engine2.rs:90-134, panic at 130): CONFIRMED. The "8 residual points" cite is wrong: they are at engine2.rs:384, 416 and 434, not 401-404.
  - x^p·y − x·y^p: CONFIRMED (standard). Over F_p, homogeneous degree ≤ p is still determined by values on F_p^n.
  - GrTable (symfun.rs:160-162): CONFIRMED.

## (2) Exact vs randomized

- A. Cube: CONFIRMED deterministic. cube.rs imports only crossing, linsys, mpoly, pmat, sbim and sparse (cube.rs:5-13). None of them use rand.
  - The cube also uses no canonical-basis (schur/hecke) data, which the report does not note.
  - local_cx asserts the unique degree-1 map (cube.rs:48-51).
- B. xengine: CONFIRMED (header xengine.rs:1-66; det_criterion at 1636-1641 and cclose.rs:2220-2253; c_diag at 1645-1660).
  - Only Field bounds (xengine.rs:1736, 2012). "Q would type-check": UNVERIFIABLE (not compiled).
  - Exit 3 on disagreement: CONFIRMED. It is not "prints nothing": the differing values go to stderr (krhom.rs:132-137).
  - Omitted: krhom uses XOpts::default() with prune_n None (xengine.rs:126, 2013). So the singular krhom path uses no kill rule.
- C. engine2: CONFIRMED.
  - Points: ≥24 (engine2.rs:369), seed 7 (krhom.rs:156), asserts at 375 and 487, select_basis at 167.
  - R1_2strand.md:312-313 says "Probabilistic in the Schwartz-Zippel sense; the standard of engine2". That sentence describes a proposed engine, attributed to engine2's standard.
  - Tower exact: CONFIRMED. rand in tower.rs is test-only (1602ff).
  - twoprime (twoprime.rs:10-11) uses khr_engine2 with seeds 7 and 11.
  - indec (indec.rs:233-299): CONFIRMED.
- D. Box path: the citations check out:
  - homprov.rs:128-136, 643-687, 788-818, 834-873, 944-988
  - sindec.rs:645-652, 774-778
  - sengine.rs:365-389, 1248-1276, 1583, scalar_of 538-563
  - swire.rs:1052, 1180-1200 (chi_check), 1206-1230, 1438-1465
  - boxclose.rs:4-33, 1223-1252
  - SPEC_M4_M6a.md:918
- Correction on classification. Some checks are one-sided: failure gives an Err and they never give a false pass. These are:
  - the stacked-rank independence check;
  - the generator determinant pass and det_criterion;
  - pattern replay.
  - The C_sound "≤ deg/p" quote (gl6/C_sound.md:224-225) is about one of these. At small p they abort runs; they do not corrupt results.
- The following checks can falsely pass, so soundness there is probabilistic:
  - form constancy at a second point;
  - e²=e and Π∘Ι=δ at 2-3 points;
  - the scaling check at 2 box points;
  - Freivalds (≤ p^-k);
  - engine2 interpolation residuals.
  - The smap samplers draw from 1..2e9, so their effective sample size is about 2e9, not p.
- "What a result certifies": condition (c), the kill rule (schur.rs:1405), applies only to pruned runs (sbox `--prune`, sbox.rs:60; xscan `--prune`). It does not apply to krhom singular, cube, cable or engine2.

## (3) Extension field

- No extension arithmetic exists: CONFIRMED.
- Rank and dimension invariance under F_p ⊂ F_{p^k}: CONFIRMED (standard). The complex is defined over F_p, so dim_{F_{p^k}} H = dim_{F_p} H.
- It cures the field-size issues (distinctness, Schwartz-Zippel, ord λ, interpolation): CONFIRMED.
- Characteristic is unchanged: CONFIRMED.
- Obstacles (prime-subfield samplers, BoxNodes using characteristic() as size, to_u64 hashing, integer-only Layer A): CONFIRMED. Add these too:
  - the xengine u32 packing;
  - indec.rs:242.

## (4) Q

- Today: CONFIRMED (krhom.rs:67, 119-121, 144-146, 153-155, 162).
- Box over Q: the binding obstacle is compile-time (`impl<F: PrimeField> BoxNodes`, boxclose.rs:177). The assert at :181 would never be reached.
- Doc quotes (C_sound.md:174-177; SPEC_3C_GL6_rev1.md:182-184; M6a1b_report.md:3-5; review_m4b/solver.md:153; sbim.rs:297-299): CONFIRMED.
- "More primes only give an upper bound": CONFIRMED. Given C_sound's premise, dim_Q ≤ min_p dim_{F_p}.
- Multimodular reconstruction and "infeasible for 3-4 strands": [reasoning], UNVERIFIABLE (no data).
- Missed option [my reasoning, conditional on C_sound's integrality premise]:
  - In each q-degree j, Σ_i dim_{F_p} H^{i,j} = Σ r_i + 2Σ t_i ≥ |χ_j|, where χ_j is exact from decat_chi.
  - So if the F_p total in q-degree j equals |χ_j|, the Q value is proven there and there is no p-torsion. Thin and concentrated cases certify this way for free.

## (5) Characteristic caveats

- All quotes are CONFIRMED: PLAN.md:100-101; docs/dev/colored_design/SYNTHESIS.md:412 and 456; NTRUNC_STUDY.md:78; m4c/L1_eklp.md:156 and 305-310; schur.rs:1405-1411; hecke.rs:3-5; data/README.md:111 ("Only Z/2 occurs").
- The list of KL-hard-wired data (term_mult, pairs mult, check_idempotent ranks, M5i, c_diag, kill rule): CONFIRMED.
- WRONG: "By the repo's own stated assumption (p > K!) … K=8 gives 40320". L1_eklp.md:309 is a passing remark that the two primes exceed 9!; it is not a derived threshold.
  - The standard condition for invariant theory to behave as in characteristic 0 for S_K is p ∤ K!, i.e. p > K.
  - K! is also not a known sufficient bound for p-canonical = KL. Only "p ≫ 0" is known, with finitely many bad primes for each W, and bad primes in SL_n grow at least exponentially in n.
  - The validated set is exactly {2^31−1, 4294967291}. Drop "K!" from the proposed p checks.
- WRONG: "The repo says nothing about the Hom formula at small p".
  - L1_eklp.md:243-249 records the SLL standing assumptions: k is a field; the realization is faithful, balanced, generalized-Demazure-surjective and reflection faithful. There is no characteristic-0 hypothesis.
  - Thm 4.16 / Cor 4.17 (the Hom formula = hom_gd) and Thm 4.13 (indecomposables are the B_p) are stated under those assumptions.
  - L1_eklp.md:156 says Demazure surjectivity holds in every characteristic for gl_K.
  - So by the repo's reading, the Hom formula is characteristic-free. The characteristic-sensitive input is ch(B_p) = b_p, i.e. the p-canonical basis.
  - [my understanding, not verified here] The permutation realization is faithful and reflection faithful in every characteristic.
- Partial: "nothing about other primes". SPEC_M4_M6a.md:1550 (Q6) names p-canonical (K8) as a cause of a multiplicity Err.
- The report's argument that silent errors are possible only through pruning: agreed.
  - [my reasoning] Where a canonical ladder has ch = b_q, as with all 2-strand types (U13), ^p b_q ≤ b_q forces ^p b_q = b_q.
  - Elsewhere, a p-canonical difference makes pair forms or idempotent ranks fail loudly.
- Foam_N / universal-construction state spaces in small characteristic: the repo does not discuss them. Unknown.
- Web:
  - The Jensen–Williamson and Jensen RS descriptions: CONFIRMED at abstract level.
  - Added: Braden found 2-torsion in integral IC stalks in type A7 (SL_8) (Williamson, intro).
    - [my understanding] That means the 2-canonical basis differs from KL somewhere in S_8, and K = 8 is exactly (2,2,2,2).
    - Whether that element is relevant to the singular types: unknown.
  - A claim that A_n for n < 7 has no bad primes appeared only in a search snippet. Unverified.

## Corrected summary

- **Supported today:** F_{2^31−1} on every path. F_4294967291 in krhom only on the singular path (both primes compared in code), and in sbox/libbuild/xscan by flag. Exact Q only on the cube path (`--field q`). Any other `--field` value is silently F31 in krhom and sbox.
- **An arbitrary large prime p < 2^32 needs:**
  - a runtime or menu Fp (global modulus or const dispatch; u128 for p ≥ 2^32);
  - rejection of unsupported values;
  - p checks: p ≥ K for the samplers, p > 8·N·#blocks for the nodes, a fresh λ with large ord_p for each D14 constant, and p ≫ weights for interpolation;
  - new per-prime libraries for 3-4 strands;
  - the same KL assumption as today (fine for p ≫ 0, never proven for any specific p).
  - Paths with no library (cube, xengine) only need the type plumbing.
- **Small primes are not feasible on the Soergel/box/xengine paths:**
  - the samplers loop or the asserts fail;
  - one-sided checks abort often and false-pass bounds weaken;
  - most importantly, the library encodes KL data, and at bad primes the p-canonical basis differs (2-torsion is known in S_8).
  - F_{p^k} fixes only the field-size problems.
  - The cube path (deterministic, no canonical-basis data) could run small p mechanically on tiny inputs. Its validity in small characteristic is untested.
- **An honest Q option:**
  - "F_p at two (or k) primes; equals Q unless every prime is a torsion prime" (the C_sound wording);
  - plus a per-q-degree exact certificate wherever the F_p total equals |χ_j|;
  - true exact Q only via the cube today (xengine-over-Q would need wiring and measurement).

Sources:
- [Jensen–Williamson, The p-Canonical Basis for Hecke Algebras](https://arxiv.org/abs/1510.01556v2)
- [Jensen, The ABC of p-Cells](https://arxiv.org/abs/1901.02323)
- [Jensen, Cellularity of the p-Canonical Basis for Symmetric Groups, Prop 2.4(iv)](https://ar5iv.labs.arxiv.org/html/2009.11715)
- [Williamson, Schubert calculus and torsion explosion (Braden A7 2-torsion; exponential growth)](https://ar5iv.arxiv.org/html/1309.5055)
- [Elias–Qi, §1.7 citing Wil17](https://ar5iv.labs.arxiv.org/html/2005.03128)