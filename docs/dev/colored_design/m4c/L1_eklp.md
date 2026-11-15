# L1: EKLP singular light leaves, extracted for M4c (krhom)

Date 2026-11-11. This is a read-only study: nothing was compiled or run in the repository.
Tags:
- [P] stated and proved in the cited paper.
- [C] conjectured or only sketched there.
- [D] deferred by the authors to future work.
- [K] my inference for krhom; to be checked when implementing.

## 0. Sources and numbering

| key | paper | what I read | local copy (under m4c/papers/) |
|---|---|---|---|
| SLL | Elias, Ko, Libedinsky, Patimo, *Singular light leaves*, arXiv 2401.03053v1 (5 Jan 2024, the only version); J. Eur. Math. Soc. 2025, doi 10.4171/JEMS/1715 | all of it, from the TeX source | eklp_2401.03053/ (src, pdf, 8clean.tex, 9clean.tex) |
| CORE | the same authors, *On reduced expressions for core double cosets*, arXiv 2402.08673 (v2 11 Nov 2024) = SLL's "[EKLP3] in preparation" | §1 and §3 (type A) | eklp_2402.08673/ |
| EKo | Elias, Ko, *A singular Coxeter presentation*, arXiv 2105.08563, Proc. LMS 2023 | §4.8-4.9, §5, §6.1, §7.3 | eko_2105.08563/ |
| DEM | EKLP, *Demazure operators for double cosets*, arXiv 2307.15021 | only through SLL (DEM Thm 4.14 = SLL Thm 8.15) | - |
| SUB | EKLP, *Subexpressions and the Bruhat order for double cosets*, arXiv 2307.15726 | only through SLL (Thm 2.16, Props 3.3, 3.33, 3.34, Lemma 3.32) | - |
| LEIB | EKLP, *The atomic Leibniz rule*, arXiv 2407.13128 (polynomial forcing) | abstract only | - |
| W, ESW | Williamson, *Singular Soergel bimodules*, IMRN 2011 (with an erratum on duality); Elias, Snyder, Williamson, *On cubes of Frobenius extensions*, 2017 | only through SLL §3-§5 | - |

**Numbering.**
- I recomputed theorem numbers from the TeX: all theorem-like environments share one counter per section, and equations are numbered globally.
- I spot-checked equation numbers against the arXiv HTML and confirmed (44), (45), (76a), (77a/b), (85)-(93), (100)-(106), (111) and (119).
- The numbers agree with those already cited in SYNTHESIS.md (Thm 4.16, Cor 4.17, Def 7.44, Thm 7.49, ...).
- My helper scripts are in m4c/tools/m4c_tools/: texnum.py, and eklp_numbers.txt, eko_numbers.txt and core_numbers.txt, which list label → number.

**Status in November 2026.**
- A web search (not exhaustive) found no paper on arXiv that gives a presentation of singular Soergel bimodules by generators and relations.
- H. Ko gave a talk at ICERM in October 2025, "Diagrammatic singular Soergel bimodules". The slides are saved in ko_icerm_2025_slides/, but their text could not be extracted.

## 1. Combinatorics

### 1.1 EKLP objects (SLL §2-§3) [P]

**Double cosets and length.**
- (W,S) is a Coxeter system. For a finitary subset I ⊂ S, ℓ(I) := ℓ(w_I).
- An (I,J)-coset p means the triple (p,I,J). Its maximal element is p̄ and its minimal element is p_.
- Length: ℓ(p) := 2ℓ(p̄) − ℓ(I) − ℓ(J) (eq. 3).

**Expressions.**
- A multistep expression is L = [[I_0 ⊂ K_1 ⊃ I_1 ⊂ ... ⊂ K_m ⊃ I_m]] (eq. 4); equalities are allowed.
- Its length is ℓ(L) = −ℓ(I_0) + 2ℓ(K_1) − 2ℓ(I_1) + ... + 2ℓ(K_m) − ℓ(I_m) (eq. 5).
- L is a *reduced expression* of p if p̄ = w_{K1} w_{I1}^{-1} ⋯ w_{Km} w_{Im}^{-1} and ℓ(p) = ℓ(L) (eq. 6).
- A singlestep expression [I_0, ..., I_d] adds (+s) or removes (−s) one simple reflection at each step.
- Multistep and singlestep versions of an expression give the same Bott-Samelson bimodule BS (§4.2).
- Notation p = q.r means that concatenating reduced expressions of q and r gives a reduced expression of p (Notation 2.1).

**Coset pairs and paths.**
- A coset pair for (I,J,s) is p ⊂ q, where p is an (I,J)-coset and q is the (I,Js)-coset containing it (Notation 2.2).
- A subordinate path t ⊂ I_• (Def 2.3) is a sequence in which:
  - t_k is an (I_0, I_k)-coset, and t_0 contains the identity;
  - each consecutive pair is a coset pair.
- At a +s step, t_{k+1} is forced: it is the coset containing t_k.
- At a −s step, t_{k+1} is any (I_0, I_{k+1})-coset contained in t_k.
- The terminus is term(t) = t_d.
- Paths play the role of subexpressions (Bruhat strolls).

**Bruhat order.** p ≤ q iff p is the terminus of a path in some (equivalently, any) expression of q. Equivalently, p_ ≤ q_ (Def 2.4; SUB Thm 2.16).

**Redundancy and core.**
- LR(p) = I ∩ p_ J p_^{-1} and RR(p) = J ∩ p_^{-1} I p_ (eq. 7).
- p̄ = w_I . p_ . (w_{RR}^{-1} w_J) (eq. 9).
- The core p^core is the (LR,RR)-coset with the same minimal element.
- p = [[I ⊃ LR]] . p^core . [[RR ⊂ J]] (§2.5; EKo Prop 4.28, the "low road").

**Defect and polynomial factor of one step (Def 3.5).**
- Ascending step [p,q]:
  - def = ℓ(q_) − ℓ(p_) − ℓ(LR q) + ℓ(LR p);
  - poly = π+(LR q)/π+(LR p), the graded rank of one invariant ring over the other (Lemma 4.3).
- Descending step [q,p]: def = ℓ(q̄) − ℓ(p̄) and poly = 1.
- The two are related by def[q,p] = def[p,q] + ℓ(Js) − ℓ(J) (Lemma 3.8).

**Singular Deodhar formula (Cor 3.7).**
- H(I_•) = Σ_{t⊂I_•} poly(t) v^{def(t)} h_{term(t)}.
- One step of it is Prop 3.4 (= W Prop 2.8): eq. (44) for a split and eq. (45) for a merge.

**Grassmannian reduction (Thm 2.6, Prop 2.9, eq. 11).** Let p ⊂ q be a coset pair for (I,J,s), and put L := RR(q).
- z is the (I,L)-coset with z_ = q_.
- n is the (L,Js)-coset whose underlying set is W_{Js}; it has the reduced expression [[L ⊂ Js]].
- m ⊂ n is the (L,J)-coset with m_ = q_^{-1} p_.
- Then p = z.m and q = z.n.
- The defect and polynomial factor of [p,q] equal those of [m,n] (Lemma 3.9). Also LR(m) = q_^{-1} LR(p) q_ (eq. 15).
- [m,n] is a *Grassmannian pair* (Notation 2.7): I ⊂ Js and n contains the identity, so only W_{Js} matters.

### 1.2 Type A dictionary for krhom

Entries marked [K] are my inference; the others are cited.

**Parabolic subsets.**
- W = S_K, with K = Σλ = 9 for λ = (3,3,3).
- A finitary subset I is a composition (its blocks). Its walls are the c with s_c ∉ I, and ℓ(I) = Σ C(block, 2).
- krhom's zero labels are invisible to EKLP: two walls at the same position.
- So one EKLP expression has many krhom realizations. They are related by E4 `Assoc`, which has degree 0 and sends 1⊗ ↦ 1⊗ [K].

**Double cosets as contingency tables.**
- An (I,J)-coset is a contingency table M: rows are I-blocks (bottom λ), columns are J-blocks (top μ). In krhom this is `schur::Table`, with the same orientation.
- ℓ(p_) = `Table::lmin` and ℓ(p̄) = `Table::lmax`.
- LR(p) is I refined by the cells M_ij, read in row order; RR(p) is J refined by the cells, read in column order (CORE Props 3.1, 3.3).
- R_p = R^{LR(p)}: polynomials symmetric within each cell.
- rank of R_p over R^λ = N(p) = `Table::ncosets`; the graded dimension is encoded by `pi_lr` = Π[M_ij]!.

**Core, squashing and atomic expressions (CORE §3.3).**
- p_ permutes the cells as blocks.
- ψ(p) ∈ S_{#cells} is the permutation taking the row-major order of the nonzero cells to their column-major order (Def 3.7, Thm 3.8).
- An atomic coset swaps two adjacent cell blocks (x,y). Its unique reduced expression is [I + s_a − s_b]: merge (x,y) into x+y, then split as (y,x). This is the thick crossing, an "H"-web (Prop 3.10).
- Atomic reduced expressions of p^core correspond to reduced words of ψ(p) (Thm 3.11).
- Atomic braid relations are composites of switchbacks with up-up/down-down relations:
  - Lemma 3.12 for a_i a_{i+1} a_i;
  - Lemma 3.13 for far commutation;
  - a Matsumoto theorem for them (Thm 3.14).
- Consequence: every p has the explicit reduced expression **Y_p := [[λ ⊃ LR]] ∘ atomic(ψ(p)) ∘ [[RR ⊂ μ]]**: split into cells, cross, merge.

**krhom rungs as EKLP expressions [K].**
- schur.rs `rung` defines E_i^(r) as "split col i+1 as (r,b−r), merge cols i,i+1". This is the multistep [[μ ⊃ ν ⊂ μ']], where ν is the common refinement.
- As a singlestep expression, a rung is [−s_{c±r} + s_c] with c = μ_1 + ... + μ_i: E uses −s_{c+r}, F uses −s_{c−r}.
- Steps that do not change the subset are dropped (zero labels, or r equal to the whole upright).
- BS(rung) = R^ν as an (R^μ, R^μ')-bimodule with shift σ = r(b−r), where b is the giver's label. This follows from SLL Def 4.6, σ = (ℓ(L) + ℓ(I_0) − ℓ(I_m))/2.
- So a krhom `Expr` is a singlestep expression with at most 2·#rungs steps, and `flows` is the rank of its localization.
- SYNTHESIS §2.3 already states that krhom ladders are EKLP singular BS bimodules with the same grading. Re-check this on SLL Ex 7.39 (§5.4, gate G-L2).
- EKLP's merge-then-split [+s_c − s_c'] (the atomic crossing) is not a krhom rung. Realize it as E^(b) followed by F through a zero upright.

**Paths in krhom terms [K].**
- A path is a sequence of tables t_k of shape λ × I_k.
- +s merges two adjacent columns; the result is forced.
- −s splits one column as (x,y): choose α_i + β_i = M_ij for each row i, with Σ α_i = x.
- Defect and poly come from lmin, lmax and the cells. The sum over paths must equal `Schur::ladder_char` (Cor 3.7); this is gate G-L1.

**Is krhom's canonical ladder reduced? [K]**
- If char L_p = b_p, then BS(L_p) ≅ B_p. The lowest degree of B_p is the degree of 1⊗ in any reduced expression of p (Thm 4.13).
- Therefore ℓ(L_p) = ℓ(p) = 2·lmax − ℓ(λ) − ℓ(μ). Check this per type (it is cheap) before using L_p as a reduced expression.

## 2. Construction (SLL §5-§7; proofs in §8)

### 2.1 Generators (the 2-category Frob of ESW, SLL §5.4-5.5) and krhom elementary maps

| Frob generator | bimodule map (evaluation functor F, §5.5) | degree | krhom (src/smap.rs) |
|---|---|---|---|
| polynomial f ∈ R^I in a region labeled I | multiplication | deg f (deg V = 2) | `Dot` (s_ν on one upright); products, `Lin` |
| cw cap m: [I+s−s] → [I] | f⊗g ↦ fg | +(ℓ(Is)−ℓ(I)) = +ab for blocks a,b | E2 `Eps`/`EpsL` |
| cw cup Δ: [I] → [I+s−s] | f ↦ f·Δ^I_{Is} (coproduct, eqs. 71-73) | +ab | E2 `Eta`/`EtaL` |
| ccw cap ∂: [Is−s+s] → [Is] | Frobenius trace ∂_{w_{Is} w_I^{-1}} | −ab | E3 `Trace` (Σ_S m(S)/e(S,A∖S)) |
| ccw cup ι: [Is] → [Is−s+s] | inclusion | −ab | E3 `Iota` |
| up/down crossing [I+t+s] ↔ [I+s+t] | 1⊗1 ↦ 1⊗1 (natural isomorphism) | 0 | E4 `Assoc` |
| sideways crossings [Is+t−s] ⇄ [Is−s+t] | right-facing: multiplication; left-facing: f ↦ ∂Δ^I_{Ist}(f), eqs. (91)-(92) | ℓ(I)+ℓ(Ist)−ℓ(Is)−ℓ(It): equals ac for adjacent walls (blocks a,b,c), 0 for far walls | composite of the rows above (eq. 84), or a closed-form template |

**What krhom does not need.** The relations (85)-(93) are not needed: krhom works with evaluated maps. Every Frob diagram is a composite of krhom's E1-E4 [K].

**Duality D (Def 5.10).**
- D flips a diagram upside down and reverses its orientation; it is contravariant and preserves degree.
- On krhom DAGs [K]: reverse the order of composition; swap Eps↔Eta, EpsL↔EtaL, Trace↔Iota, TraceL↔IotaL; replace Assoc by its inverse; leave Dot unchanged; track krhom's E2/E3 normalization scalars.

**Type A closed forms.**
- Merging blocks (a,b): ∂(f) = Σ_{σ ∈ S_{a+b}/(S_a × S_b)} σ(f/μ), with μ = Π_{i∈A, j∈B}(x_i − x_j), the product-coproduct element (§5.2). The dual bases are the standard Schur-polynomial bases.
- P_L with ∂_L P_L = 1 (Notation 7.31) is a product over blocks of staircase monomials.
- In the gl_K realization ∂_s x_i = 1, so generalized Demazure surjectivity holds in every characteristic.

### 2.2 Rex moves (SLL §6; EKo §5) [P]

**Braid relations.**
- Up-up: [I+s+t] ↔ [I+t+s]. Down-down: the same with removals.
- Switchback (eq. 100): [I + u_0 − u_d] ↔ [I − u_1 + u_0 − u_2 + u_1 − ... − u_d + u_{d−1}].
- The u_i form the rotation sequence (EKo Def 5.10, Prop 5.14): u_0 = s, u_{−1} = w_{Js} t w_{Js}, u_{i+1} = w_{I_i} u_{i−1} w_{I_i}, with I_i = Js ∖ u_i.

**Switchback in type A (EKo §6.1, eq. 108).**
- In S_{n+1}, let J = S ∖ s_a and choose b with a + b ≠ n + 1. Then [J + s_a − s_b] ≅ [J − s_c + s_a − s_b + s_c], where:
  - c = a + b if a + b ≤ n, and c = a + b − n − 1 if a + b > n + 1. Always d = 2.
- Ladder reading: "merge (a, n+1−a), then split as (b, n+1−b)" ≅ split–merge–split–merge on three uprights (SLL Ex 6.1: S_10, a=3, b=6, c=9).

**Singular Matsumoto (EKo Thm 5.30).**
- Any two reduced expressions of p are related by braid relations.
- The proof is a constructive induction on ℓ(p), using EKo Prop 4.21 (the high road) and Prop 4.28 (the low road).

**Maps attached to braid relations (§6.1).**
- Up-up and down-down give crossings, which are isomorphisms (eq. 88).
- A switchback gives two maps whose composite ψ∘φ is the identity (Thm 6.7, via the circle-evaluation Lemma 6.10):
  - φ, built from right-facing crossings and ccw cups;
  - ψ, built from left-facing crossings and ccw caps.
- So a switchback gives the inclusion and projection of a direct summand, not an isomorphism.

**Properties of rex moves.**
- Every rex move has degree 0 and sends 1⊗ ↦ 1⊗ (Thm 6.11).
- A rex endomorphism ≡ id modulo End_{<p} (Cor 6.12). This uses Prop 4.19 (= SUB Prop 3.33): a degree-0 map between reduced expressions of p lies in Hom_{<p} iff it kills 1⊗.

### 2.3 The five stages [P]

**Stage 1: elementary light leaves (§7.1-§7.3).** Input: a Grassmannian pair m ⊂ n for (I,J,s) with I ⊂ Js.

*The expressions.* Fix a reduced expression K_• of m^core and enumerations of I∖LR(m) and J∖RR(m) (Notation 7.1). Then:
- X_m = [I − t_1 − ... − t_l] ∘ K_• ∘ [RR(m) + u_1 + ... + u_r];
- X_n = [I + v_1 + ... + v_k], with the order of the v_i given by Lemma 7.2.

*ELL([m,n]): X_m∘[J+s] → X_n* is the "sinister" diagram (Def 7.3, Def 8.1):
- it uses only left-facing crossings and ccw caps; each color t ∈ I is a cap, each t ∈ Js∖I a through-strand plus caps;
- no two strands cross twice, and there are no polynomials.
- It is built by the right-to-left algorithm of Lemma 8.2, which keeps the dangling strands in end-order.

*Degree.* deg ELL([m,n]) = def[m,n] = −ℓ(m_) + ℓ(LR m) − ℓ(I) (Lemma 7.4, §8.2). In type A:
- inside the merged block, the I-blocks c_i split as (k_i, c_i − k_i), with Σ k_i = a, so deg = −ℓ(m_) − Σ k_i(c_i − k_i) [K].

*Variants.*
- ELLP([m,n],f) puts f ∈ R^{LR(m)} in the leftmost LR(m)-region of the domain (Def 7.14).
- ELL([n,m]) := (id_m ⊗ cw cap) ∘ (D ELL([m,n]) ⊗ id_[−s]) (Def 7.16, eq. 133). Its degree is def[n,m] (Lemma 7.17) and it carries no polynomial.
- At a forward step, ELL([n,m]) followed by a rex move is the identity (Ex 7.18, Lemma 8.14).

*Evaluations (§8.4).*
- D ELL([m,n]) and ELL([n,m]) both send 1⊗ ↦ 1⊗ (Lemmas 8.9, 8.10).
- For b ∈ R^J, **ELLP([m,n],f)(1⊗⊗b) = ∂^{LR(m)}_I(f·∂_{m_}(b))·1⊗_n** (Algorithm 8.12, inside the proof of Lemma 8.11; here BS(X_n) ≅ R^I).

**Stage 2: singlestep light leaves (§7.4).** For any coset pair p ⊂ q:
- take z, m, n from Thm 2.6 and any reduced expression Z of z; set X_p := Z∘X_m and X_q := Z∘X_n (Notation 7.22);
- SSLL([p,q]) = id_Z ⊗ ELL([m,n]) and SSLL([q,p]) = id_Z ⊗ ELL([n,m]) (Def 7.23);
- with polynomials: id_Z ⊗ ELLP([m,n],f) for f ∈ R^{LR(m)}; degrees equal the defects (Lemma 7.24).

**Stage 3: light leaves (§7.5, Def 7.28).**
- LL(t_{≤k+1}) = SSLLP([t_k,t_{k+1}], f_{k+1}) ∘ (rex ⊗ id) ∘ (LL(t_{≤k}) ⊗ id).
- The rex move goes from the previous target X'_{t_k} to the source X_{t_k} that the next step requires.

**Stage 4: sprinkles (§7.6).** At each ascending index j:
- L_j := t_{j+1}_^{-1} LR(t_j) t_{j+1}_ ⊂ R_{j+1} := RR(t_{j+1}) (Notation 7.29, Lemma 7.30);
- y_j runs over the minimal representatives of W_{L_j}\W_{R_{j+1}} (a sprinkled path, Def 7.32);
- the sprinkle is **f_j = ∂_{L_j} ∂_{y_j}(P_{R_{j+1}})** (Def 7.33).
- The f_j form a basis of R^{L_j} over R^{R_{j+1}}, of graded rank poly[t_j,t_{j+1}] (Lemma 7.34, via DEM Thm 4.14 = Thm 8.15).
- Descending steps carry f = 1; with I_0 = ∅ there are no sprinkles at all (Rem 7.35).
- Count: Σ_{(t,y): term = p} v^{deg LL} = c_p(I_•) (Lemma 8.7).

**Stage 5: double leaves (§7.8, Def 7.44).**
- Y_p := [[I ⊃ LR(p)]] ∘ (a reduced expression of p^core) ∘ [[RR(p) ⊂ J]].
- Take a coterminal sprinkled quadruple (p, (t,y), (t',y'), g ∈ R_p) (Def 7.42).
- DLL = D LL(t',y') ∘ rex ∘ g ∘ rex ∘ LL(t,y), with g placed in the leftmost LR(p)-region of Y_p.

**Free choices.** All of the following may be chosen freely, and every choice gives a basis (§7.10):
- reduced expressions of cores and of z;
- the enumerations;
- the rex moves;
- the polynomials P_L;
- the basis B_p of R_p.

Rem 7.11-7.13 add that one may insert rex isomorphisms (crossings) at will.

## 3. Basis theorems (precise)

**Standing assumptions (§4.1).**
- k is a field. The realization is faithful, balanced, satisfies generalized Demazure surjectivity, and is reflection faithful. Such a realization exists over ℝ for every W. R = Sym(V) with deg V = 2.
- Hom^i(M,N) = Hom^0(M, N(i)). The convention p·M := ⊕ M(−i)^{a_i} is opposite to Williamson's.
- BS(L) = R^{I_0} ⊗_{R^{K_1}} R^{I_1} ⊗ ⋯ ⊗_{R^{K_m}} R^{I_m}(σ) (Def 4.6). The element 1⊗ spans its lowest degree, −σ (eq. 53).

**Results.**
- **Thm 4.13** (W Thm 7.10 and SUB Prop 3.3): the indecomposables are the B_p(m). For a reduced expression of p, the summand containing 1⊗ is B_p.
- **Thm 4.16 and Cor 4.17** (Soergel-Williamson Hom formula):
  - gdim Hom(BS I_•, BS I'_•) = Σ_p Σ_{t,t'→p} poly(t) poly(t') v^{def t + def t'} gdim R_p;
  - R_p is the standard bimodule: R^{LR(p)} with the right action twisted by p_ (Def 4.15);
  - this is what schur.rs `hom_gd` computes.
- **Prop 4.21**: Hom_{≮p}(B, BS(rex of p)) ≅ c_p·R_p.
- **Prop 8.32**: the set {rex∘LL(t,y) : term(t) = p} is a basis of Hom_{≮p}(BS I_•, BS Y_p) as a left R^{LR(p)}-module. Its proof combines independence (Thm 8.31), the count (Lemma 8.7) and graded Nakayama.
- **Cor 8.33**: the flipped light leaves form a basis of Hom_{≮p}(BS Y_p, BS I_•).
- **Thm 7.49 (main theorem)**: {DLL(T,b)} is a k-basis of Hom(BS I_•, BS I'_•) in all degrees, where T runs over coterminal sprinkled triples and b over B_p.
  - Rem 7.48/7.50: it is also a basis as a left R^I-module (take b in a basis of R_p over R^I), and similarly as a right R^J-module.
- **Lemma 8.35**: double leaves through q ≤ p span Hom_{≤p}. This holds for any downward-closed set (Rem 8.36).
- **Thm 7.51** (proof in §8.10): other choices give a unitriangular change of basis with respect to ≺_3 (Def 8.37). For light leaves alone, see Thm 8.34.
- **Thm 7.52**: SBSBim(J,I) is a fibered cellular category in the sense of [ELauda, Def 2.17], with cell algebras A_p = R_p.
  - So local intersection forms control decompositions (intro §1.1).
  - For krhom, this means the multiplicity of B_q(s) equals the rank of the degree-0 pairing between light leaves of degree −s and flipped light leaves of degree s. Factors from R_q of positive degree vanish modulo m_q [K].

**The proof engine** (reusable as a structural certificate in krhom [K]).
- Dual sequence: y_j° = w_{I_j} t_j_^{-1} t_{j+1}_ y_j w_{I_{j+1}} (Def 8.19, eq. 200).
- Test elements: b(u,z) = 1 ⊗ b_1 ⊗ ... ⊗ b_d ∈ BS(I_•), with b_k = ∂_{I_k} ∂_{z_k°}(P_{I_{k+1}}) at ascending k and b_k = 1 otherwise (Notation 8.20).
- Total order ≺: lexicographic order on y°, using a length-compatible order on W (Def 8.23, Lemma 8.25).
- **Prop 8.26**: LL(t,y)(b(u,z)) = 1⊗ if (t,y) = (u,z), and = 0 if (t,y) ≺ (u,z).
- The argument uses only three facts: rex moves have degree 0 and preserve 1⊗; Lemmas 8.9-8.11; degree counting.

## 4. What is proved, conjectured, or deferred

**[P] in SLL** (for any Coxeter system, under the §4.1 assumptions):
- Thm 2.6 (Grassmannian reduction), Cor 3.7 (singular Deodhar formula), Lemma 7.4 (sinister diagrams exist, with degree = defect);
- Thm 6.7, Thm 6.11 and Cor 6.12 (rex moves);
- Prop 8.26, Thm 8.31, Prop 8.32 and Cor 8.33 (light leaves);
- Thms 7.49, 7.51, 7.52, 8.34 (double leaves basis, change of basis, fibered cellular).
- These rest on proved results elsewhere:
  - W Thm 7.9 (Hom formula; W's duality error is fixed by an erratum);
  - DEM Thm 4.14;
  - SUB Thm 2.16, Prop 3.3, Lemma 3.32 and Props 3.33-3.34;
  - EKo Prop 4.28, Prop 5.14 and Thm 5.30;
  - the 2-category Frob and its evaluation functor F from ESW.

**[P] in CORE (type A):**
- squashing (Thm 3.8); atomic reduced expressions correspond to reduced words of ψ(p) (Thm 3.11); atomic Matsumoto (Thm 3.14);
- the atomic nil-category is equivalent to a nilHecke algebroid (Prop 3.19).

**[C] conjectured or sketched:**
- Ex 7.18: at a forward step, ELL([n,m]) is itself a rex move. Not needed: Lemma 8.14 suffices.
- Rem 8.3: the sinister cap diagram is unique up to isotopy (only sketched).
- Rem 8.4: any "reasonable" diagram gives the same ELL.
- EKo §7.3: Elias's double ladders (type A webs) are reduced expressions.

**[D] deferred:**
- A complete presentation of SSBim by generators and relations. These are ingredients 5-6 of §1.1: the additional relations, and a diagrammatic proof that double leaves span. F: Frob → SBSBim is full (by Thm 7.49) but not faithful, and the missing relations are unknown in general (§5.3).
- Rem 6.13: proving in Frob itself that rex endomorphisms ≡ id modulo lower terms requires relations specific to the Soergel setting.
- Rem 4.12 and Rem 5.12: an intersection form on singular BS bimodules, and the natural isomorphism D∘F ≅ F∘D.
- Rem 7.25: the two actions of R^{LR(p)} on SSLL agree modulo lower terms. This is "future work"; polynomial forcing is treated in LEIB.

**Consequence for krhom [K].**
- SPEC_M4_M6a named one main risk for M4c: the companion paper on singular Soergel calculus has not appeared. That risk does not apply.
- SLL proves the basis statement for evaluated bimodule maps (Thm 7.49). krhom computes exactly these, in localized form, and never uses relations.
- What krhom needs is explicit evaluation of the Frob generators, and it already has E1-E4.
- On the assumptions:
  - krhom works over F_p with p = 2^31−1 (and 4294967291) and the gl_K realization;
  - this realization is faithful, reflection faithful, balanced and Demazure surjective;
  - p is larger than K! = 9!, so the invariant theory behaves as in characteristic 0;
  - the Hom formula is already assumed by `hom_gd`.

## 5. Implementing the leaves in krhom

### 5.1 What M4c must deliver (the `HomProvider` trait, src/homprov.rs)

**`closure(lp, p, outs)`: the 12 Exhausted closure types [K].**
- The identity (λ,λ)-coset p0 = diag(λ) is Bruhat-minimal. Its only reduced expression is [λ], so Y_{p0} = [λ], and R_{p0} = R^λ.
- Hence {LL(t,y) : t ⊂ L_p, term(t) = diag(λ)} is an R^λ-basis of Hom(BS L_p, R^λ). This is Prop 8.32 with no lower terms; see also the footnote in §7.7.
- outs = false uses D(LL). The count is c_diag(b_p)(1) = `Schur::c_diag`, and the degree multiset is c_diag(b_p) (Lemma 8.7). This matches the trait's contract exactly.
- Neither Y_p, nor g, nor a final rex move is needed.

**A plausible explanation of why Search misses the top-degree generators [K].** At top degree the sprinkles f_j have y_j = id, so polynomials of maximal degree sit at several intermediate levels at once ("polynomials within the diagram are unavoidable", §1.3, Ex 7.39). Search's seeds (S4: unit/counit plus recursion of depth ≤ 4) and its single dot segment rarely produce such maps.

**`pairs(x, q, lq, sdeg, mult)`: the 11 Exhausted pairs, and also the Refused and Memory jobs [K].**
- π-candidates: ρ_q ∘ LL(t,y) with term(t) = q and deg LL = −s.
- ι-candidates: D(LL) ∘ σ_q with deg LL = +s.
- ρ_q: X_q → L_q and σ_q: L_q → X_q are any degree-0 maps that send 1⊗ ↦ 1⊗. They are unique modulo Hom_{<q} by Prop 4.19, and they induce Hom_{≮q}(X, Y_q) ≅ Hom_{≮q}(X, L_q), because BS(L_q) ≅ B_q.
- mult is the rank of the local intersection form, the coefficient of 1⊗ in π_i ι_j(1⊗), which is the existing `form_at`. No Search windows are needed.
- Types that need an idempotent (`NeedsIdempotent`) would go through Y_q instead (later).

**`full` / `end0`: double leaves (Def 7.44).** These need Y_q and a basis of R_q. They are not needed for the acceptance set.

### 5.2 Ingredients: explicit formula, or algorithm only

| ingredient | reference | explicit? | where in krhom |
|---|---|---|---|
| tables, lmin/lmax, LR/RR cells, N(p) | SLL §2.2; CORE Props 3.1/3.3 | yes | `schur::Table` |
| paths, def, poly; Σ = character | Defs 2.3, 3.5; Cor 3.7 | yes (enumeration) | new; oracle `ladder_char`, `split_col`, `merge_cols` |
| z, m, n (Grassmannian reduction) | Thm 2.6, eq. (11) | yes (permutations in S_K) | new |
| reduced expression of a core (atomic) | CORE Thm 3.11 | yes (reduced word of ψ(p)) | new |
| X_m, X_n, Z | Notation 7.1/7.22, Lemma 7.2 | yes | new |
| ELL (sinister diagram) | Def 7.3, Lemma 8.2, eq. (133) | yes (planar algorithm) | new compiler: diagram → `Elem`/`Local` nodes |
| sideways crossings | eqs. (84), (91)-(92) | yes | `tmpl.rs` template per (a,b,c) |
| sprinkles f_j, P_L, ∂_y | Notation 7.31, Def 7.33 | yes (divided differences, ≤ 9 variables) | new "polynomial in region" node, evaluated per flow |
| test elements b(u,z) | Notation 8.20, Prop 8.26 | yes | optional structural certificate |
| elementary rex moves | §6.1, eqs. (100)-(106); EKo eq. (108) | yes | templates |
| **rex paths between two given reduced expressions** | EKo Thm 5.30 (constructive induction); type A cores: CORE Thm 3.14 | algorithm only, no closed form | planner, with fallback below |
| ρ_q, σ_q | Prop 4.19, Thm 6.11 | via a rex path, or any degree-0 map preserving 1⊗ | planner or one `Exact` call on the differing window |

**Fallback for rex moves [K].**
- In §8, rex moves enter only through "degree 0 and 1⊗ ↦ 1⊗": Thm 6.11 inside Prop 8.26, and Lemma 8.7.
- So any degree-0 bimodule map between reduced expressions of one coset that preserves 1⊗ can replace a rex move. By Prop 4.19 it equals the rex move modulo Hom_{<p}. The same holds for ELL([n,m]) at forward steps (Lemma 8.14).
- On the flipped side (Cor 8.33) there is no diagram to flip, so build the reverse map separately.
- krhom's numerical certification (count, stacked rank, genuineness) covers any remaining gap.

### 5.3 Module sketch [K]

**src/dcoset.rs (combinatorics only):**
- permutations of S_K; Table → p_, p̄, LR, RR, core, ψ;
- Grassmannian reduction;
- conversion between `Expr` and singlestep words (Vec<Add(c)|Rem(c)>), with zero-upright placement;
- builders for the reduced expressions above (through-core atomic, X_m, X_n, Z);
- the braid-relation planner (EKo Thm 5.30, memoized).

**src/leaves.rs:**
- enumeration of sprinkled paths with their duals y°;
- the ELL compiler;
- SSLL, LL, D, DLL;
- the `Leaves` provider, plugged in as Chain(Exact, Search, Leaves) (SPEC_M4_M6a lines 164 and 1060; SYNTHESIS §3 lists `leaves.rs`).

**src/smap.rs:**
- a node for a polynomial in a region (product of per-cell Schur polynomials, or direct evaluation at the flow's root assignment);
- templates for sideways crossings and for switchback φ/ψ (`Local` nodes in tmpl.rs).

**Evaluation of π.** Propagate a row vector from right to left (the target R^λ has one flow). This avoids materializing intermediate matrices.

### 5.4 Gates and tests

- **G-L1** (pure combinatorics, seconds): path sums equal `ladder_char` for every (3,3,3) canonical ladder and each one-rung extension. Lemma 8.7 sums equal `c_diag` and `term_mult`.
- **G-L2** (SLL examples at small λ, with `sbim::hom_space` as the exact oracle):
  - Ex 7.39, λ = (2), [s,∅,s] → [s]: two leaves, `Trace` with the sprinkle P_s = x_1 and `Trace` with 1, of degrees ±1.
  - Ex 7.40/7.41, S_3 with I = (2,1), J = (1,2): 4 light leaves and 5 light leaves.
  - Ex 7.46: 10 double leaves from [I−s+t−t+s−s+t] (the rung word F_1^(1) E_1^(1) F_1^(1) on (2,1)) to [I−s+s+t−s].
  - Ex 2.5/2.10/7.6 (S_11): the core reduced expression (34) and ELL (124).
- **G-L3**: the M4c acceptance set from M5b_report.md, at level 6, for both primes:
  - closure types 0,3,0/1,0,2/2,0,1, 1,0,2/2,0,1/0,3,0, 1,2,0/0,0,3/2,1,0, 0,0,3/3,0,0/0,3,0, 0,3,0/0,0,3/3,0,0, 1,1,1/1,1,1/1,1,1, 0,2,1/2,0,1/1,1,1, 1,0,2/1,2,0/1,1,1, 1,1,1/0,2,1/2,0,1, 1,1,1/1,0,2/1,2,0, 0,1,2/2,0,1/1,2,0, 0,2,1/1,0,2/2,1,0;
  - pairs 71, 159, 268, 322, 323, 325, 334, 336, 361, 368, 419;
  - then the 38 Refused and 7 Memory pairs.

### 5.5 Risks [K]

1. **Size of intermediate expressions.** Through-core reduced expressions pass through the cell parabolic.
   - For p = 1,1,1/1,1,1/1,1,1, BS(Y_p) localizes to 6³·2⁹ = 110 592 flows, against 4 320 for L_p.
   - Light-leaf prefixes are tensored with the rest of the ladder.
   - Mitigations: choose Z with few flows (it is free); for closures, avoid Y_p entirely; propagate vectors instead of materializing matrices.
2. **Planning rex paths.** The number of reduced expressions grows quickly. Use type A normal forms (CORE §3.4, EKo §6.1) and the fallback in §5.2.
3. **Conventions.** krhom's E2/E3 scalars, D on DAGs, and the degree normalization (D14). These are harmless for spanning but must be consistent for ι/π.
4. **Effort.** SYNTHESIS estimates 10-20 days. The ELL compiler and the rex planner are the two new pieces of real algorithmic work; everything else is bookkeeping on tables.
