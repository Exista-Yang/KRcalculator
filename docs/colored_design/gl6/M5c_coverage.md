# M5c coverage: which 3-strand knots a pruned (3,3,3) scan can already do at N = 5, 6

Date 2026-11-12, branch colored-m5. Tool: `libbuild cover` (src/bin/libbuild.rs) over `sindec::knot_jobs`
(src/sindec.rs), test tests/knotcov.rs. Status source: the M5b run files of both primes (data-gl6/m5b/{f31,f32b}/runs/
run_n5.jsonl, run_n6.jsonl, state of M5b: level 5 171/177, level 6 309/377). A job counts as certified only if it is
certified for **both** primes (the lists are identical anyway). No provider call is made.

## Method

For a braid word on 3 strands (all labels 3), starting from the identity type of (λ, λ), every crossing at strands
(i, i+1) is applied to every type p present: for every Rickard term W_k that the cell rule keeps at N
(max(c, d) + k ≤ N; at N = 5 the term k = 3 of the (3,3) crossing is dropped), the rungs of W_k are applied one at a
time; each (type, rung) step is a pair of `pruned_plan(N)` (i.e. a `pair K` library job), and its kept summand types
(`term_mult`, types killed at N dropped) feed the next rung. After the last crossing every non-identity type present
needs its closure generators (a `closure P` job). This is exactly how `dec_word` / `dec_composite` / `rickard_comp_at`
decompose (rung by rung, killed terms refused), so the job set is that of a localized scan of the word.

Two bounds:
- **UB** (upper bound): types united over all homological degrees, no cancellation.
- **LB** (lower bound): only types whose signed, q-graded class Σ_h (−1)^h [object] in K_0 is nonzero after the
  crossing; these are present in every complex homotopy equivalent to the stage (Gaussian elimination cannot remove
  them), so their steps through the next crossing, and the closure of the final ones, are needed by *any* scan.
  Checked against `schur::decat_scan` (unpruned K_0 class, tests/knotcov.rs).

## Result

**No 3-strand knot is fully covered at N = 5 or at N = 6, not even with the lower bound.** The pruned scan of the
smallest knot (4_1 = σ1 σ2⁻¹ σ1 σ2⁻¹) already reaches almost the whole level set: its 4 stages hold 3, 9, 19, 20 kept types
(UB) at N = 5, i.e. every kept type of (λ, λ) after 4 crossings, and the job set of every knot contains the jobs that fail.

| N | knots (3-strand list, 185) | jobs per knot UB (min-max) | failures blocking every knot (LB) | blocking failures per knot (LB) | fully covered (UB / LB) |
|---|---|---|---|---|---|
| 5 | 185 | 108-137 | 6 of 6 | 6-6 | 0 / 0 |
| 6 | 185 | 283-377 | 57 of 68 | 57-68 | 0 / 0 |

So knot coverage does not come in steps: at N = 5 the six level-5 failures (pairs 313, 322, 326, 406, 410 and closure
1,1,1/1,1,1/1,1,1) are each needed by every knot (LB), and at N = 6 57 of the 68 failures are needed by every knot.
The engineering target of M5c is therefore the failure list itself, in full: **N = 5 becomes computable for all
3-strand knots exactly when those 6 jobs are certified**; N = 6 needs at least the 57 universal blockers, and all 68
for most knots.

### Failures at N = 6 that not every knot needs (LB)

| failure | knots (of 185) whose LB needs it |
|---|---|
| p24 | 104 |
| p78 | 104 |
| p137 | 161 |
| p145 | 161 |
| p196 | 161 |
| p219 | 161 |
| p227 | 161 |
| p37 | 161 |
| p45 | 161 |
| cl 1,1,1/0,2,1/2,0,1 | 170 |
| cl 1,0,2/1,2,0/1,1,1 | 181 |

Universal N = 6 blockers (needed by every knot, LB): cl 0,0,3/3,0,0/0,3,0, cl 0,1,2/2,0,1/1,2,0, cl 0,2,1/1,0,2/2,1,0, cl 0,2,1/2,0,1/1,1,1, cl 0,3,0/0,0,3/3,0,0, cl 0,3,0/1,0,2/2,0,1, cl 1,0,2/2,0,1/0,3,0, cl 1,1,1/1,0,2/1,2,0, cl 1,1,1/1,1,1/1,1,1, cl 1,2,0/0,0,3/2,1,0, p127, p159, p163, p195, p199, p212, p268, p27, p273, p278, p281, p283, p285, p292, p293, p313, p319, p322, p323, p325, p326, p332, p334, p335, p336, p361, p368, p376, p381, p383, p385, p392, p406, p410, p411, p418, p419, p444, p445, p455, p458, p460, p467, p475, p66, p71, p82.

## Per knot, N = 5 (knots with ≤ 8 crossings; all 185 in the summary above)

Columns: jobs = pairs + closure types, UB (LB); certified UB (LB); failures in the LB / only in the UB.

| knot | crossings | jobs UB (LB) | certified UB (LB) | blocking (LB) | UB-only failures |
|---|---|---|---|---|---|
| 4_1 | 4 | 89+19 (89+15) | 102/108 (98/104) | p313, p322, p326, p406, p410, cl 1,1,1/1,1,1/1,1,1 | - |
| 5_2 | 6 | 89+19 (89+19) | 102/108 (102/108) | p313, p322, p326, p406, p410, cl 1,1,1/1,1,1/1,1,1 | - |
| 6_2 | 6 | 89+19 (89+19) | 102/108 (102/108) | p313, p322, p326, p406, p410, cl 1,1,1/1,1,1/1,1,1 | - |
| 6_3 | 6 | 91+19 (91+16) | 104/110 (101/107) | p313, p322, p326, p406, p410, cl 1,1,1/1,1,1/1,1,1 | - |
| 7_3 | 8 | 89+19 (89+19) | 102/108 (102/108) | p313, p322, p326, p406, p410, cl 1,1,1/1,1,1/1,1,1 | - |
| 7_5 | 8 | 91+19 (91+19) | 104/110 (104/110) | p313, p322, p326, p406, p410, cl 1,1,1/1,1,1/1,1,1 | - |
| 8_2 | 8 | 89+19 (89+19) | 102/108 (102/108) | p313, p322, p326, p406, p410, cl 1,1,1/1,1,1/1,1,1 | - |
| 8_5 | 8 | 113+19 (113+19) | 126/132 (126/132) | p313, p322, p326, p406, p410, cl 1,1,1/1,1,1/1,1,1 | - |
| 8_7 | 8 | 91+19 (91+19) | 104/110 (104/110) | p313, p322, p326, p406, p410, cl 1,1,1/1,1,1/1,1,1 | - |
| 8_9 | 8 | 91+19 (91+16) | 104/110 (101/107) | p313, p322, p326, p406, p410, cl 1,1,1/1,1,1/1,1,1 | - |
| 8_10 | 8 | 115+19 (115+19) | 128/134 (128/134) | p313, p322, p326, p406, p410, cl 1,1,1/1,1,1/1,1,1 | - |
| 8_16 | 8 | 118+19 (118+19) | 131/137 (131/137) | p313, p322, p326, p406, p410, cl 1,1,1/1,1,1/1,1,1 | - |
| 8_17 | 8 | 118+19 (118+16) | 131/137 (128/134) | p313, p322, p326, p406, p410, cl 1,1,1/1,1,1/1,1,1 | - |
| 8_18 | 8 | 118+19 (112+17) | 131/137 (123/129) | p313, p322, p326, p406, p410, cl 1,1,1/1,1,1/1,1,1 | - |
| 8_19 | 8 | 113+19 (113+19) | 126/132 (126/132) | p313, p322, p326, p406, p410, cl 1,1,1/1,1,1/1,1,1 | - |
| 8_20 | 8 | 113+19 (113+19) | 126/132 (126/132) | p313, p322, p326, p406, p410, cl 1,1,1/1,1,1/1,1,1 | - |
| 8_21 | 8 | 115+19 (115+19) | 128/134 (128/134) | p313, p322, p326, p406, p410, cl 1,1,1/1,1,1/1,1,1 | - |

## Per knot, N = 6 (knots with ≤ 8 crossings; all 185 in the summary above)

Columns: jobs = pairs + closure types, UB (LB); certified UB (LB); failures in the LB / only in the UB.

| knot | crossings | jobs UB (LB) | certified UB (LB) | blocking (LB) | UB-only failures |
|---|---|---|---|---|---|
| 4_1 | 4 | 244+39 (244+33) | 224/283 (220/277) | 57 (universal 57) | cl 1,0,2/1,2,0/1,1,1, cl 1,1,1/0,2,1/2,0,1 |
| 5_2 | 6 | 244+39 (244+39) | 224/283 (224/283) | 59 (universal 57 + cl 1,0,2/1,2,0/1,1,1, cl 1,1,1/0,2,1/2,0,1) | - |
| 6_2 | 6 | 244+39 (244+39) | 224/283 (224/283) | 59 (universal 57 + cl 1,0,2/1,2,0/1,1,1, cl 1,1,1/0,2,1/2,0,1) | - |
| 6_3 | 6 | 250+39 (250+35) | 230/289 (227/285) | 58 (universal 57 + cl 1,0,2/1,2,0/1,1,1) | cl 1,1,1/0,2,1/2,0,1 |
| 7_3 | 8 | 244+39 (244+39) | 224/283 (224/283) | 59 (universal 57 + cl 1,0,2/1,2,0/1,1,1, cl 1,1,1/0,2,1/2,0,1) | - |
| 7_5 | 8 | 250+39 (250+39) | 230/289 (230/289) | 59 (universal 57 + cl 1,0,2/1,2,0/1,1,1, cl 1,1,1/0,2,1/2,0,1) | - |
| 8_2 | 8 | 244+39 (244+39) | 224/283 (224/283) | 59 (universal 57 + cl 1,0,2/1,2,0/1,1,1, cl 1,1,1/0,2,1/2,0,1) | - |
| 8_5 | 8 | 324+39 (324+39) | 297/363 (297/363) | 66 (universal 57 + cl 1,0,2/1,2,0/1,1,1, cl 1,1,1/0,2,1/2,0,1, p137, p145, p196, p219, p227, p37, p45) | - |
| 8_7 | 8 | 250+39 (250+39) | 230/289 (230/289) | 59 (universal 57 + cl 1,0,2/1,2,0/1,1,1, cl 1,1,1/0,2,1/2,0,1) | - |
| 8_9 | 8 | 250+39 (250+35) | 230/289 (227/285) | 58 (universal 57 + cl 1,0,2/1,2,0/1,1,1) | cl 1,1,1/0,2,1/2,0,1 |
| 8_10 | 8 | 330+39 (330+39) | 303/369 (303/369) | 66 (universal 57 + cl 1,0,2/1,2,0/1,1,1, cl 1,1,1/0,2,1/2,0,1, p137, p145, p196, p219, p227, p37, p45) | - |
| 8_16 | 8 | 338+39 (338+39) | 309/377 (309/377) | 68 (universal 57 + cl 1,0,2/1,2,0/1,1,1, cl 1,1,1/0,2,1/2,0,1, p137, p145, p196, p219, p227, p24, p37, p45, p78) | - |
| 8_17 | 8 | 338+39 (338+35) | 309/377 (306/373) | 67 (universal 57 + cl 1,0,2/1,2,0/1,1,1, p137, p145, p196, p219, p227, p24, p37, p45, p78) | cl 1,1,1/0,2,1/2,0,1 |
| 8_18 | 8 | 338+39 (326+37) | 309/377 (295/363) | 68 (universal 57 + cl 1,0,2/1,2,0/1,1,1, cl 1,1,1/0,2,1/2,0,1, p137, p145, p196, p219, p227, p24, p37, p45, p78) | - |
| 8_19 | 8 | 324+39 (324+39) | 297/363 (297/363) | 66 (universal 57 + cl 1,0,2/1,2,0/1,1,1, cl 1,1,1/0,2,1/2,0,1, p137, p145, p196, p219, p227, p37, p45) | - |
| 8_20 | 8 | 324+39 (324+39) | 297/363 (297/363) | 66 (universal 57 + cl 1,0,2/1,2,0/1,1,1, cl 1,1,1/0,2,1/2,0,1, p137, p145, p196, p219, p227, p37, p45) | - |
| 8_21 | 8 | 330+39 (330+39) | 303/369 (303/369) | 66 (universal 57 + cl 1,0,2/1,2,0/1,1,1, cl 1,1,1/0,2,1/2,0,1, p137, p145, p196, p219, p227, p37, p45) | - |

## Reproduce

```
R=data-gl6/m5b   # worktree of origin/data-gl6
libbuild cover --lam 3,3,3 --n 5 --knots tools/dev/ntrunc_t1/data/knots3_all.txt \
    --runs $R/f31/runs/run_n5.jsonl,$R/f32b/runs/run_n5.jsonl [--max-cross 8] [--names 4_1,5_2] [--list] [--md OUT]
```
(0.4 s per N for all 185 knots.)

## After M5c

With the M5c run files of both primes (level 5: 172/177; level 6: 330/377 F31, 329/377 F32b; a job counts if
certified for both primes):

| N | fully covered (UB / LB) | blocking failures per knot (LB) | failures needed by every knot (LB) |
|---|---|---|---|
| 5 | 0 / 0 | 5 | 5 of 5: pairs 313, 322, 326, 410, closure 1,1,1/1,1,1/1,1,1 |
| 6 | 0 / 0 | 38-48 | 38 of 48 |

4_1 at N = 6 now has 243 of its 283 jobs certified (M5b: 224). Still no knot is computable; N = 5 for every 3-strand
knot is exactly the 5 jobs above away.
