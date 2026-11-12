# L2_alt: constructive Hom bases other than EKLP, judged for krhom M4c

Date 2026-11-11. Read-only study: no cargo, no binaries. Papers were downloaded as data under
`m4c/papers/<paper>/` (PDF plus extracted text `paper.flat`); arXiv HTML was used for KRW and EKLP.

Tags:
- [L] the theorem was read in the paper text (number checked);
- [L?] the reference is from memory or from an abstract only (number not checked);
- [D] derived here, not checked;
- [E] estimate.

## 0. Bottom line

1. **No published alternative gives a basis of Hom between singular ladders that is more explicit than EKLP.**
   Every other source gives one of three things:
   - a basis of a different space: regular Soergel bimodules, sl_n webs (decategorified), or cyclotomic Foam_N state spaces;
   - an existence or fullness theorem (KRW, QR, Williamson);
   - a spanning set that is complete but redundant (Khovanov-Lauda diagrams).
2. **Two facts from the alternatives are directly useful.**
   - **(a) Completeness of krhom's local generators.** By KRW Thm 4.6 (SSBim ≅ braid-like foams), QR Prop 3.15
     (every foam between ladder webs is a ladder foam, possibly after adding zero uprights) and EKLP Thm 7.49 (double
     leaves are built from Frobenius structure maps and polynomials), every Hom between ladders is a linear combination
     of composites of krhom's elementary pieces (dots, zip/unzip, digon cup/cap) plus small templates (fusion/fission
     of consecutive rungs, thick crossings). The intermediate expressions may need extra zero uprights [D]. Search is therefore missing
     *enumeration completeness*, not generator types.
   - **(b) Reduction to the identity coset.** Williamson's translation adjunction (thesis Lemma 7.2.1, proof of
     Thm 7.2.2) turns every Hom between ladders into a Hom(R, ladder), that is, elements supported on the identity
     coset. Closure generators are the primitive case.
3. **Most implementable complete alternative: KL/foam matching enumeration.** This is the Khovanov-Lauda spanning set
   B_{i,j,λ} (KL3 Prop 3.11, "for any intermediate choices") in thick form (QR Def 2.2, Thm 3.9, Cor 3.10; KLMS
   thick calculus).
   - It is pushed into SSBim through the foamation and skew-Howe 2-functor.
   - Candidates are perfect thick matchings of E- and F-strands, with Schur decorations, built from a fixed library of
     templates on at most 3 uprights.
   - It is complete if the 2-functor is full ([L] QR Prop 3.15 for Foam_N; [D] for SSBim via KRW Thm 4.6).
   - It needs no rex moves, no cores and no exact windows larger than about 10^4 unknowns [E].
   - It covers the **closure** failures well: one ladder, evaluation only at the identity coset.
   - For **pairs** with large X its candidate count may explode. EKLP, with one light leaf per flow at rep(q), stays
     the economical basis there. Measure the matching counts first; this is combinatorics only, see §5.
4. **Not practical here:**
   - Regular reduction (B_p as an idempotent image in a regular Bott-Samelson bimodule, Elias thick calculus): rank
     blow-up by |W_λ| = 216 per side, and the idempotent comes from the rex graph of w_0(S_6) / w_0(S_9).
   - Elias light ladders: decategorified, webs rather than foams.
   - Foam_N cellular bases (Hu-Mathas, Tubbenhauer web algebras): finite N, cyclotomic, one-sided boundary. This
     does not match krhom's N-independent R^λ-linear model.
   - Moment-graph / GKM descriptions: Bott-Samelson fibres are not GKM, and BMP sheaves presuppose B_p.

## 1. What M4c must produce (from the M5b report and the code)

Failing level-6 targets of the (3,3,3) library: 38 Refused, 23 Exhausted, 7 Memory, identical for both primes.
Source: origin/colored-m5 docs/colored_design/gl6/M5b_report.md, "Level 6".

| kind | Hom space | provider entry | what certification needs |
|---|---|---|---|
| C: closure, 12 types Exhausted (e.g. 0,0,3/3,0,0/0,3,0, X = 400, degree 18; 1,1,1/1,1,1/1,1,1, X = 4 320) | Hom^d(R^λ, L_p), Hom^d(L_p, R^λ), d = top degrees of c_diag(b_p) | `homprov.rs` `Search::closure_inner` (tiers S1/S6: whole-ladder exact window; S2: interior dots × lower generators; S4: `closure_seeds`) | n = c_diag(b_p)(1) maps with that degree multiset and a nonsingular n×n identity-coset matrix at 2 points (`certify_closure`). **Only identity-coset blocks are evaluated.** |
| P: pairs, 11 Exhausted (322, 323 at degree 0; 336 at degree -1; ...) | ι ∈ Hom^s(L_q, X), π ∈ Hom^{-s}(X, L_q), X = L_p·F_1^(1), X = 2 880-14 400 flows | `Search::pairs_inner` (S0-S3, S6) | form π_j∘ι_k at rep(q) of rank `mult` (`certify_pairs`), genuineness (`assert_genuine`). **Only the rep(q) blocks are evaluated.** |
| R/M: 38 Refused, 7 Memory | the same spaces | Search window admission; exact window calls of 0.5-4 GB | the same; the cause is the size of exact windows (>= 495 k unknowns) |

Generator types that genuineness admits (`homprov.rs` `assert_genuine`):
- exact templates (`sbim::hom_space` on small windows, `tmpl.rs`);
- Rickard templates;
- `Elem` closed forms with recorded constants (`smap.rs` `ElemKind`: Eps/Eta/EpsL/EtaL = E2, Iota/Trace/IotaL/TraceL = E3, Assoc = E4);
- `Dot` (E1, Schur polynomial s_ν of the roots on an upright).

Any new method must express its maps as `Comp`/`Lin`/`Kron` DAGs over these leaves. New leaves must be small
templates.

## 2. The sources, one by one

### 2.1 Regular light leaves and double leaves (Libedinsky; Elias-Williamson)

- **What they give.**
  - Libedinsky, J. Algebra 320 (2008), "Sur la catégorie des bimodules de Soergel": light leaves, and the theorem
    that dots, trivalent vertices and 2m-valent vertices generate all maps between Bott-Samelson bimodules [L?]
    (EW cite it for fullness in Thm 6.28 [L]).
  - Libedinsky arXiv:1304.1448 (Adv. Math. 2015): double leaves, a cellular basis [L?, abstract].
  - Elias-Williamson, Soergel calculus, arXiv:1309.0865: **Thm 6.11** LL_{x,y} is a free R-basis of Hom(BS(x), BS(y))
    in D; **Prop 6.12** LL_{x,e} is a free R-basis of Hom(BS(x), 1); **Thm 6.28** D ≅ BSBim under Demazure
    surjectivity [L]. §4: rex moves of two different paths agree modulo lower terms (localized, they agree on the top
    summand) [L].
- **Ring.** R = Sym(h*) (regular, W = S_K), left module, free.
- **What is explicit.** A basis indexed by 01-subexpressions, i.e. paths in the Bruhat graph (Deodhar). Each leaf is
  a composite of U0/U1/D0/D1 steps (dot, trivalent, identity) and rex moves (2m-valent vertices).
- **What krhom would need.** Regular objects (labels 1 only). For the colored target, singular objects must first be
  induced to regular ones (§2.2).
- **Covers the failures?** Only through §2.2, which is impractical. Two parts carry over:
  - the template of the argument (EKLP generalizes it);
  - the fact that rex moves need only agree modulo lower terms (EW §4, EKLP Cor 6.12). This is the theoretical basis
    for krhom's "any degree-0 map with nonzero top scalar" leaf isomorphisms.

### 2.2 Singular objects as summands of regular ones (Elias "Thick Soergel calculus in type A")

- **What it gives.** Elias, arXiv:1009.2120 (PLMS 112 (2016) 924-978):
  - an **explicit idempotent** in End(BS(rex of w_0(S_n))) whose image is B_{w_0} = R ⊗_{R^{S_n}} R; it is built
    from the reduced-expression graph of w_0 with the Manin-Schechtman semi-orientation (§3);
  - a diagrammatic presentation of generalized Bott-Samelson bimodules R ⊗_{R^{J_1}} R ⊗ ... (§4, thick calculus);
  - induced trivial representations, i.e. one-sided singular modules (§5) [L].
- **Ring.** R = Q[x_1..x_n]. The objects are *regular* R-bimodules. The two-sided singular 2-category R^I-R^J is not
  covered.
- **What krhom would need.**
  - Replace L_p by Ind L_p = R ⊗_{R^λ} L_p ⊗_{R^λ} R.
  - Compute a regular double-leaves basis of Hom(Ind L_q, Ind X).
  - Recover Hom(L_q, X) as its image under the Demazure trace ∂_{w_λ} ⊗ ∂_{w_λ}. This is a spanning set, by the
    Frobenius adjunction (Williamson Lemma 7.2.1 [L]).
- **Practicality for λ = (3,3,3)** [D/E]:
  - |W_λ| = 216, so the localized rank grows by 216 per induced side: X = 21 600 flows becomes >= 4.7e6 per side.
  - Intermediate rungs pass through labels 6 and 9. The Elias idempotents then live on BS(rex w_0(S_6)) (length 15,
    localized rank 2^15) and BS(rex w_0(S_9)) (length 36). Their construction walks rex graphs with 292 864 (S_6) and
    far more (S_9) vertices.
  - The regular basis is redundant by a factor |W_λ|^2 = 46 656 over the singular space.
  - **Verdict: not practical.** krhom already handles singular objects natively. Its decomposable-type idempotents
    e_p (M5i, `check_idempotent`) are the only idempotents it needs.

### 2.3 Williamson, "Singular Soergel bimodules" (arXiv:1010.1283v3, IMRN 2011)

- **What it gives.** No basis, but the exact structure behind krhom's certification [L]:
  - Thm 1: indecomposables ↔ double cosets.
  - **Thm 7.2.2**: Hom(M, N) is graded free over R^I, with Hom(M, N)[-ℓ(w_J)] ≅ ⟨ch_Δ M, ch_∇ N⟩ R^I, for M a
    Bott-Samelson bimodule (in IB^J_BS) and N with a ∇-flag.
  - **Thm 7.4.1**: the same for all Soergel bimodules.
  - **Thm 7.3.5**: Hom(R_p, B) ≅ Hom(R_p, (support truncation of B at p)) and dually, up to shifts [-2ℓ(p_-)].
- **Constructive content of the proof of Thm 7.2.2.**
  - Lemma 7.2.1 (translation adjunction J ⊂ K) moves the factors of M onto N.
  - This reduces everything to M = R^I, where **Hom(R^I, N) = Γ_{W_I}(N)**: the sections of N supported on the
    identity double coset, of rank = coefficient of H_id in ch_∇ N.
- **For krhom** [D]:
  - Every target is Γ_id of a bent ladder. Hom(L_q, X) ≅ Hom(R, X·L_q^∨) up to shift, using units and counits of
    rung adjunctions, which are `ElemKind` Eta/Eps/Iota/Trace and templates.
  - The closure problem is the primitive one. A single complete generator for "elements of a ladder supported at the
    identity coset", plus bending maps, serves closure, pairs and `full`.
  - Thm 7.3.5 says pairs modulo lower terms live in the stalk at q. This is why only rep(q) blocks are needed.
  - The identity-coset rank test in `certify_closure` is exactly the statement "rank Γ_id = c_diag(b_p)(1)".
- **Covers the failures?** It is the reduction, not the generator. It does not by itself produce the top-degree
  closure generators.

### 2.4 Elias, "Light ladders and clasp conjectures" (arXiv:1510.06840, v1 only)

- **What it gives** [L, text read]:
  - A basis of **Hom spaces in Fund+(sl_n)** (sl_n webs of Cautis-Kamnitzer-Morrison) over Z[q, q^-1]: the double
    ladders. Light ladders are in §2.5.
  - **Thm 2.39** (double ladders span), Lemmas 2.40/2.41 (strong and weak light ladder lemma), §2.8 (linear
    independence).
  - Remark 2.38: "neutral ladders are to sl_n-webs what rex moves are to Soergel diagrams". The weak lemma 2.41
    allows *any* neutral ladder.
  - Intro: the double ladders basis "is part of a larger double singular light leaves basis" for singular Soergel
    bimodules of *affine* type A (via geometric Satake), announced with Williamson. That basis became EKLP.
- **Ring and space.** Webs are 1-morphisms. Their Hom spaces are decategorified: elements of Hom_{U_q(sl_n)}, i.e. of
  the q-Schur algebra under skew Howe. They are not spaces of bimodule maps.
- **For krhom.** Nothing directly. Two uses only:
  - its combinatorics (paths of dominant weights) is the decategorified shadow of krhom's flows;
  - Lemma 2.41 is an analogue of "any rex move works" (SYNTHESIS U1 / D1 Q6).
- **Covers the failures?** No.

### 2.5 Categorified skew Howe: Khovanov-Lauda, Queffelec-Rose, Lauda-Queffelec-Rose, Mackaay-Stošić-Vaz

- **Khovanov-Lauda III**, "A categorification of quantum sl(n)", arXiv:0807.3250, QT 1 (2010) [L]:
  - §3.2 **Prop 3.11**: the set B_{i,j,λ} (minimal diagrams for each strand matching, dots at fixed positions, bubble
    monomials in the rightmost region) spans HOM_U(E_i 1_λ, E_j 1_λ) "for any intermediate choices made".
  - Prop 3.12: the graded count.
  - Thm 1.3: nondegeneracy for sl_n, proved through the 2-representation on cohomology of iterated flag varieties
    (§5-6, Flag_N). That 2-representation is the type-A flag model underlying SSBim and Foam_N.
- **Queffelec-Rose**, arXiv:1405.5920, Adv. Math. 302 (2016) [L]:
  - Def 2.2: thick (divided power) U̇_Q(sl_m).
  - Def 3.1: foam 2-category nFoam_m(N).
  - Thm 3.9: 2-representation Φ_n: U̇_Q(sl_m) → nFoam_m(N). Cor 3.10: images of thick-calculus generators.
  - Lemma 3.13: web isotopies are realized by images.
  - **Prop 3.15 (fullness)**: every foam between ladder webs is Φ_{m',N}(D) for some m' >= m. Extra uprights with
    label 0 are allowed.
  - Prop 3.21 and Prop 3.22: nFoam(N) is the direct limit of N-bounded quotients of U̇_Q(gl_m).
  - Intro: Tubbenhauer adapted the Hu-Mathas basis to get explicit bases of sl_n web algebras (§2.6 below).
- **Lauda-Queffelec-Rose**, arXiv:1212.6076 (AGT 2015) [L?]: the same picture for N = 2 (Khovanov homology).
  N = 2 only.
- **Mackaay-Stošić-Vaz**, "A diagrammatic categorification of the q-Schur algebra", arXiv:1008.1348, QT 4 (2013)
  [L]:
  - Ŝ(n,d) is a quotient of U(gl_n), with a 2-representation F_Bim to bimodules (Prop 4.3, from KL Thm 6.13).
  - For n >= d, a *fully faithful* 2-functor from regular Soergel bimodules into Ŝ(n,d).
  - Freeness and End(1_λ) are only **conjectured**: Conj 5.5, Conj 5.6.
  - So MSV does **not** prove fullness onto singular Soergel bimodules. This settles SYNTHESIS §2.4's open "fullness
    result of MSV type": it is not in MSV. The relevant statement is QR Prop 3.15 for Foam_N, and KRW (§2.7) for
    SSBim.
- **What they give krhom: a complete, explicit candidate stream.** It is a spanning set, not a basis:
  - **Hom(L_1, L_2)** is the span over R^λ of the images of the thick B_{i,j,λ}, where i, j are the rung words of
    L_1, L_2. Bubbles map to R^λ. Colors beyond the words appear only in bubbles, so extra zero uprights contribute
    only scalars [D].
  - **Hom(R^λ, L)** is spanned by perfect thick matchings of the E- and F-rungs of L (cups and caps), made adjacent
    by crossings, with decorations on the strands.
  - **The local pieces are a fixed list:**
    - thick dots: Schur polynomials of the moved roots of a rung, i.e. `Dot` differences;
    - fusion and fission E^(a)E^(b) ↔ E^(a+b) on 2 uprights;
    - same-colour thick crossing (= fission∘fusion with Schur decorations, degree -ab);
    - adjacent-colour thick crossing E_i^(a)E_{i+1}^(b) → E_{i+1}^(b)E_i^(a) on 3 uprights;
    - general thick cups and caps E^(r)F^(r) ↔ 1 on 2 uprights. krhom's E2/E3 closed forms cover only "all of an
      upright".
  - Each piece is a `sbim::hom_space` window on <= 3 uprights with <= 3 rungs and labels summing to <= 9. This is
    small: for example, two rungs of size 3 have <= C(9,3)·C(6,3) = 1 680 flows [E], against >= 495 k unknowns in the
    Refused windows.
- **Completeness argument** [D]:
  - KL3 Prop 3.11 gives spanning in U̇.
  - Fullness of U̇(gl_{m'}) → SSBim between ladder words of colors 0,1: QR Prop 3.15 proves it for Foam_N. For SSBim
    (N = ∞) it follows from KRW Thm 4.6 plus the same argument applied to braid-like foams. Alternatively, from EKLP
    Thm 7.49: double leaves are composites of Frobenius structure maps of parabolic inclusions (EKLP §5) and
    polynomials, and those are QR's foam generators (3.2)-(3.4).
  - **Caveat:** this has not been checked line by line. The certification makes it moot in practice, since a
    candidate family is accepted only at the Hom-formula count.
- **Cost risk.**
  - The spanning set is redundant. Thin strands are not an option: thinning E^(r) multiplies flows by r! per rung.
    So the thick form is mandatory.
  - The number of thick matchings is a product of contingency-table counts per colour. It is small for one ladder
    (closure), but can be 10^4-10^8 for Hom(L_q, X) with X of 8-10 thick rungs [E]. It must be measured before
    building (§5).

### 2.6 Foam evaluation and state-space bases

- **Robert-Wagner**, "A closed formula for the evaluation of foams", arXiv:1702.04140, QT 11 (2020) [L]:
  - §2: evaluation of closed gl_N foams as a sum over colourings (facets coloured by subsets of {1..N}).
  - Thm 3.13: categorified MOY square relation, with LR coefficients.
  - **Thm 3.30**: universal-construction functor F: Foam_N → Pol_N, with graded rank of F(Γ) = ⟨Γ⟩_MOY. Rem 3.31:
    free over Z[X_1..X_N], with dual bases.
  - **Prop 4.13**: explicit basis of F(Γ) only for flag-type webs (Schur polynomials s_λ, λ ∈ T(a_i, Σ_{j<i} a_j),
    with an explicit dual basis). Prop 4.14: structure constants.
  - For general webs there is no explicit basis. One applies the categorified MOY isomorphisms recursively, which is
    Lewark-style web simplification.
- **Robert-Wagner**, "Symmetric Khovanov-Rozansky link homologies", arXiv:1801.02244, JEP 7 (2020) [L]:
  - ∞-equivalence of foams.
  - Lemma 3.30: every rooted foam is ∞-equivalent to a combination of tree-like foams. Lemma 3.32: tree-like foams
    with equal decorations are equivalent. Lemma 3.40: annular version.
  - Conj 4.16: disk-like foams ≅ singular Soergel bimodules (now KRW Thm 4.6).
  - Lemma 4.17: HH_* cyclicity.
  - These describe the *elements* of a singular Bott-Samelson bimodule (decorated tree-like foams = polynomials on
    edges). They do not give Hom spaces.
- **Ehrig-Tubbenhauer-Wedrich**, arXiv:1703.06691, PLMS 2018 [L]:
  - GL_N-equivariant foams (§2.3-2.4); Lemmas 2.29/2.30: colouring relations from RW.
  - Ex 2.6: Sym(A|B) is free over Sym(A ∪ B) with Schur basis s_α(A), α ∈ P(a,b). This is already the krhom digon
    convention.
  - Thm 3.5: invariance in K(SFoam).
  - No Hom basis.
- **Hu-Mathas**, arXiv:0907.2985 (Adv. Math. 225 (2010)) [L?]: graded cellular basis ψ_{st} of cyclotomic KLR algebras
  of type A.
- **Tubbenhauer**, arXiv:1404.5752 [L?, abstract and QR intro]: an explicit basis of the sl_n web algebra from it.
  These are bases of Foam_N spaces with all boundary on one side, at **finite N**, in a cyclotomic (highest-weight)
  quotient.
- **For krhom.**
  - krhom works N-independently in SSBim over R^λ and prunes with Foam_N at the end (`Schur::killed_at`).
  - Switching the Hom layer to Foam_6 would change the localization (RW colourings with N colours, not krhom's flows
    with K = 9 roots), the module structure (finite rank over Sym(X_1..X_N)) and the closure.
  - **Verdict:** useful as an oracle (RW evaluation of closures, Prop 4.13 bases for flag webs as unit tests). Not
    as M4c.
- **Covers the failures?** No.

### 2.7 Khovanov-Robert-Wagner, "A foamy approach to Soergel bimodules" (arXiv:2603.22957, 2026)

- **What it gives** [L, arXiv HTML read]:
  - **Thm 4.6**: Rep: blFoam_k → SSoe is a 2-equivalence between braid-like foams modulo ∞-equivalence and singular
    Soergel bimodules of type A.
  - Lemma 4.4: Rep(Γ) ≅ F(Γ) as graded spaces. The bimodule of a web = element-type foams. Lemmas 4.2/4.3 are
    RW20b 3.30/3.32.
  - **Lemma 4.7**: Rep is a graded isomorphism on 2HOM. **Lemma 4.8** reduces it to identity 1-morphisms using the
    adjunction isomorphisms (§3.3, (17)-(21)) and the Queffelec-Rose algorithm (**Prop 2.11** = QR18 Lemma 5.2:
    Γ ⊕ ⊕ Q·S ≃ ⊕ P·S, annular webs are stably sums of concentric circles).
  - The identity case: decorated disks span (surjectivity). Injectivity uses Prop 4.9 = QRS18 Thm 3.2 (S¹-symmetry)
    and bubble and diabolo removal (Lemmas 2.8, 2.9).
  - Graded ranks: eq. (14) (End over Q) and (15) (over A_k).
- **What it does not give.** Any explicit basis or spanning set of 2HOM(Γ_0, Γ_1).
- **Constructive reading** [D]:
  - The proof is an algorithm for **closure-type spaces**: bend to Hom(id, Γ), close annularly, simplify by
    QR-algorithm isomorphisms to concentric circles, and take decorated disks.
  - Each QR step is a local 2-upright decomposition: E^(a)F^(b) ≅ ⊕ F^(b-j)E^(a-j) (categorified commutator) and
    E^(a)E^(b) ≅ [a+b, a] E^(a+b). It needs explicit inclusions and projections, i.e. small pairs that krhom can
    template.
  - The extra Q·S on the left makes the result a spanning set after projection. The certification selects the
    basis.
  - This is an alternative closure generator. It is more delicate than §2.5: annular moves, the Q·S bookkeeping and
    long chains of isomorphisms. Its value is the **completeness of foam generators for SSBim**, the missing half of
    §2.5's argument.
- **Covers the failures?** Closure, in principle. Pairs only through bending.

### 2.8 Other relevant items

- **Elias-Ko**, "A singular Coxeter presentation", arXiv:2105.08563:
  - Thm 5.30: singular Matsumoto. Any two reduced expressions of a double coset are related by singular braid
    relations.
  - A type-A web presentation [L?, abstract].
  - These are the moves of EKLP §6. In krhom they are degree-0 maps between reduced ladders.
- **Ko**, arXiv:2312.16666 (atomic presentation). EKLP Rem 2.11 cites it (§3) for type-A core rexes.
- **EKLP follow-ups:**
  - arXiv:2307.15021 (Demazure operators for double cosets);
  - arXiv:2307.15726 (subexpressions and Bruhat order for double cosets);
  - arXiv:2407.13128 (atomic Leibniz rule, RMI 42 (2026)).
  - No published "type-A cores" companion was found (search 2026-11-11). The EKLP risk stands.
- **Moment graphs / GKM** (Fiebig; Williamson §4.4, §7.3 localization). krhom's coset evaluation is this
  localization. But the image of a ladder in ⊕_c L_c is not cut out by edge conditions: Bott-Samelson-type fibres
  have repeated tangent weights and are not GKM. The BMP-sheaf description needs B_p, not L_p. Not constructive here.
- **Hogancamp-Rose-Wedrich**, arXiv:2107.08117 (skein relation for SSBim) [L?, abstract]: explicit SSBim maps for
  Rickard complexes. Possibly a source of closed forms for thick crossings. Not read.

## 3. Comparison

| method | basis of | ring / module | explicit? | krhom work | closure (C) | pairs (P) | refused / memory (R/M) | main risk |
|---|---|---|---|---|---|---|---|---|
| EKLP singular light leaves (baseline) | Hom(BS_sing, BS_sing) (Thm 7.49) | R^I-free | yes, modulo choices (core rexes, rex moves, sprinkles) | Grassmannian reduction (Thm 2.6), cores (§2.5), elementary leaves (§7.1-7.3), sprinkles (§7.2, 7.5-7.6), rex moves (§6); expressions with up to 9 uprights | yes: leaves to the identity coset, one per identity flow | yes: one leaf per flow of X at rep(q) | yes: no big exact windows | type-A core rexes unpublished (EKLP24+); size of the combinatorics; flows of intermediate expressions with more uprights |
| Regular leaves after induction (EW Thm 6.11, Elias thick idempotents) | regular Hom(BS, BS) | R-free | yes (regular) | induction, Demazure trace, w_0(S_6)/w_0(S_9) idempotents | no (×216 per side) | no | no | infeasible size |
| Williamson adjunction | none (reduction) | R^I-free (Thm 7.2.2, 7.4.1) | the reduction is | bending maps (units/counits) | primitive case | reduces to C on X·L_q^∨ | n/a | none: glue only |
| Elias light ladders | Hom in Fund(sl_n) (webs) | Z[q^±] | yes | n/a | no | no | no | wrong level (decategorified) |
| KL / QR thick spanning set | spans HOM_U (KL3 Prop 3.11); images span SSBim Hom if full (QR Prop 3.15; KRW Thm 4.6) | R^λ (bubbles) | yes | matching enumerator, about 5 template kinds on <= 3 uprights | yes: matchings of one ladder, identity-coset evaluation only | yes, but candidate count possibly large | yes: templates <= 10^4 | redundancy for large X; fullness for SSBim not checked line by line |
| KRW / QR annular algorithm | none (equivalence) | A_k | the algorithm is | QR moves as templated pairs, annular bookkeeping | yes, in principle | via bending | yes | Q·S bookkeeping; long isomorphism chains |
| RW / ETW evaluation, Prop 4.13 | F(Γ) for flag webs | Sym(X_1..X_N) | for special webs | n/a | no | no | no | finite N; wrong model |
| Hu-Mathas / Tubbenhauer | cyclotomic KLR, sl_n web algebras | field, finite N | yes (tableaux) | re-architecture | no | no | no | wrong model |

## 4. Judgment

1. **The only basis-level competitor is EKLP itself.** All alternatives either change the setting (regular, webs,
   finite-N cyclotomic) or prove existence and fullness without a basis (Williamson, MSV, QR, KRW).
2. **The most implementable complete alternative is the thick KL/foam spanning set** (§2.5). The reasons:
   - The local pieces are 5 small template kinds, all expressible as genuine `Local`/`Elem`/`Dot` leaves.
   - Enumeration is pure combinatorics (thick matchings plus Schur decorations, filtered by degree). There are no rex
     moves, cores or Grassmannian pairs, and no window larger than about 10^4 unknowns.
   - Certification (`certify_closure`, `certify_pairs`, `certify_full`) is unchanged. Redundancy costs evaluations,
     not correctness.
   - Completeness rests on KL3 Prop 3.11 plus fullness. Fullness is QR Prop 3.15 for Foam_N and KRW Thm 4.6 for
     SSBim [D]. This is a stronger guarantee than Search's (none).
3. **Where it fits best: the closure failures (C).**
   - Candidates are matchings of one ladder's own rungs: L_p has <= ~10 thick rungs and few matchings.
   - Evaluation is needed only at the identity coset, giving n-vectors for an n×n test.
   - Greedy selection by ascending degree (the existing `Selector`) is sound: n Frac-independent homogeneous elements
     with degree multiset c_diag(b_p) form an R^λ-basis, because the determinant is a degree-0 polynomial.
4. **Where EKLP stays better: pairs (P) with large X.**
   - EKLP gives exactly one candidate per flow of X at rep(q) with the right defect, a handful.
   - KL matchings of Hom(L_q, X) may number 10^4-10^8 before the degree filter [E].
   - The degree-0 and degree -1 Exhausted pairs (322, 323, 336) are KL-μ-type lower summands. Every method needs
     *global* information on the whole ladder L_p there: full light leaves, or matchings across all of X.
5. **Recommended M4c shape** [D]:
   - **(a) A KL-span tier "S7"** in `homprov.rs`, after S6, for `closure_inner` first and then `pairs_inner` /
     `full_inner`. Its enumerator lives in a new `src/klspan.rs`: thick matchings of E/F rungs per colour,
     contingency tables between thick strands, Schur decorations by degree.
   - **(b) Bending maps** (Williamson Lemma 7.2.1), so that pairs can reuse the closure enumerator on X·L_q^∨ when
     that is cheaper.
   - **(c) EKLP singlestep and full light leaves** (`leaves.rs`, SYNTHESIS §3.1 M4c) only for the targets whose S7
     matching count, measured in §5, exceeds a budget.

   This orders the work by risk:
   - S7 needs no unpublished input.
   - EKLP's type-A core rexes can be replaced by krhom's numerically found degree-0 maps (EKLP Cor 6.12, EW §4:
     correct modulo lower terms).

## 5. Concrete next steps (cheap, no heavy compute)

1. **Matching census.** Combinatorics only, laptop-safe, seconds. For each of the 68 failing level-6 targets, use
   `schur.rs` words: canonical ladders from `sindec::canonical_ladder`, X = L_p·F_1^(1). Count:
   - (i) thick perfect matchings of L_p (closure) and of X·L_q^∨ (pairs), per degree;
   - (ii) the number of distinct local templates the matchings need.

   Decision rule: S7 for every target whose degree-filtered count is <= ~10^4; EKLP leaves for the rest.
2. **Template inventory.** Fusion/fission, thick crossings (same and adjacent colour) and general cups/caps for
   labels summing to <= 9 on <= 3 uprights.
   - Check that each is a `hom_space` window under the existing Exact caps: unknowns <= 40 000, fibre product
     <= 1 600 / 6 400.
   - Prefer closed forms where EKLP §5 gives Frobenius maps, e.g. general cups = Frobenius unit and counit of
     R^{(a+r, b-r)} ⊃ R^{(a, r, b-r)}, with constants as in `smap.rs` E2/E3.
3. **Gate for S7** (the cloud runs it):
   - B2 (2,1,2) through S7 alone must reach 605/605 (Search alone: 581/605);
   - G1(a) 480/480;
   - then the 23 Exhausted, 38 Refused and 7 Memory level-6 targets.
4. **Literature check still open:**
   - fullness of U̇(gl_m) → SSBim (N = ∞) between ladder words, stated explicitly (KRW §4.2 proof of Lemma 4.7, or
     QR16 Prop 3.15 adapted to blFoam);
   - KLMS thick spanning sets (arXiv:1006.2866 [L?]) for thick same-colour crossings;
   - HRW arXiv:2107.08117 for closed-form thick crossings in SSBim.

## 6. Files

- Papers (data): `m4c/papers/{elias_light_ladders_1510.06840, elias_thick_1009.2120, williamson_ssbim_1010.1283,
  ew_soergel_calculus_1309.0865, rw_closed_formula_1702.04140, rw_symmetric_1801.02244, etw_functoriality_1703.06691,
  qr_foam_1405.5920, msv_qschur_1008.1348, kl3_sln_0807.3250}/paper.{pdf,flat}`.
- Text extraction helper (not a paper): `m4c/tools/pdftxt.py`, `m4c/tools/kw.py`.
- KRW (2603.22957) and EKLP (2401.03053) were read through arXiv HTML. The EKLP text was truncated after §7.1, so
  theorem numbers 7.49+ are from the introduction and SYNTHESIS.
