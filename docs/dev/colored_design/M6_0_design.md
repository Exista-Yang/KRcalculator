# M6.0 design note: exact small colored engine `src/xengine.rs`

Worktree `.claude/worktrees/m6`, branch colored-m6 (base 99c2d13). Checked against SYNTHESIS
§9.2 (M6.0), §1.4 (disputed points), §2.3 (conventions), §3.1 (acceptance), §4.1 L3 / §4.2.

## Scope

2 strands (1 strand accepted as the trivial case), labels per bottom position in {1,2,3},
mixed labels allowed when the braid closes up. Exact arithmetic over a prime field F
(generic `F: Field`; tests use F31 and F32b). N-independent scan, closure per N.

## Data structures

```
Lad<F>       { rungs: Vec<Rung>, b: SBim<F> }           // canonical 2-rung ladder L_p of type p, e_p = id
Summand<F>   { q: Tid, s: i32, iota: PMat<F> /*L_q -> X*/, pi: PMat<F> /*X -> L_q*/ }
XObj         { p: Tid, h: i32, t: i32 }                  // object B_p{t} in homological degree h
XCx<F>       { objs: BTreeMap<id, XObj>,
               out:  BTreeMap<id, BTreeMap<id, PMat<F>>>, // component L_src -> L_tgt (exact, over R^lambda)
               inn:  BTreeMap<id, BTreeSet<id>> }
XEngine<F>   { lam, schur: Schur (one per lambda), ladders: Tid -> Lad,
               locals: (c,d,positive) -> cube::LocalCx (untruncated, n = c+d),
               decomp: (p, positive, k) -> Rc<Vec<Summand>>,
               dloc:   (p, positive, j) -> transported id(x)d_j components,
               decat:  BTreeMap<h, schur::Elt> (decategorified scan, for the K_0 assert) }
XScan<F>     { data: ClosureDataC<F>, stages: Vec<StageInfo>, timings, max hom_space time }
```

Type ids are schur `Tid`s of the single `Schur(lambda)` instance (bottom = lambda, top = the
current labels). A component B_p{t} -> B_q{t'} is stored as a PMat L_p -> L_q of normalized
degree e = t - t' (unnormalized u = e - shift(L_q) + shift(L_p)); every stored component is
checked to be homogeneous of that degree (total q-degree 0).

Grading: object (p, h, t) means q^t L_p (K_0 class v^t b_p). Term k of a crossing contributes
h += hdeg_k (= +-k) and t += dq_k (= -k positive, +k negative), the N-independent part of
convention A; the N-dependent -eps m(N-m) goes into `cclose::global_q`.

## Canonical ladders (T-h2)

For each type p: among all ladders with at most 2 rungs from lambda to mu (on 2 uprights, labels
unbounded) whose Schur character equals b_p exactly, take the one with fewest flows, then fewest
rungs, then lexicographic (deterministic). Assert one exists (T-h2) and fsize(L_p) = rank(B_p).
Then L_p = B_p and e_p = id. End^0(L_p) = k is a consequence (checked through the scalar asserts
below and by an explicit test on all types used).

## Order of operations per crossing (generator 1, sign eps, labels mu = (c, d))

1. `LocalCx` from `cube::local_cx(g, c, d, n = c + d)` (made `pub`): term ladders W_k (rungs,
   hdeg), unique degree-1 differentials as exact PMat over R^{(c,d)}.
2. For every old object o = (p, h, t) and term k: decomposition of X = L_p (x) W_k
   (`SBim::add_rungs`), cached per (p, eps, k):
   - multiplicities m_q(v) = `schur::term_mult(p, W_k.rungs)`;
   - for each (q, s) with c = [v^s] m_q: iota-candidates = `hom_space(L_q, X)` in normalized
     degree s, pi-candidates = `hom_space(X, L_q)` in normalized degree -s;
     intersection form G_ji = pi_j o iota_i, asserted to be a scalar multiple of id (End^0 = k),
     read exactly; assert rank G = c; choose c rows/cols, normalize pi' = G_JI^{-1} pi_J;
   - assert sum over summands of fsize(L_q) = fsize(X) (ranks add up);
   - orthogonalize across summands: M = Pi o Iota is unitriangular (pi_tau o iota_sigma has degree
     s_sigma - s_tau, zero unless s_sigma > s_tau or sigma = tau); Pi <- M^{-1} Pi by the Neumann
     series (nilpotent part, terminates); assert Pi o Iota = id exactly.
   New objects (q, h + hdeg_k, t + dq_k + s), one per summand.
   Implementation update (after measuring): the decomposition is done rung by rung (L_p (x) R_1 =
   (+) L_r{s_1}, then L_r (x) R_2 = (+) L_q{s_2}; iota = (iota_1 (x) id) o iota_2, pi = pi_2 o (pi_1 (x) id)),
   each step by the procedure above, and the composite multiplicities are asserted to equal
   term_mult(p, W_k). Reason: hom_space on L_p (x) W_k (fibre 144 for (3,2)) took up to 275 s per call,
   rung by rung the largest fibre is 36 and the (2,2) trefoil scan drops from 8 s to 0.05 s. A
   tensor product that is literally the canonical ladder of its only summand (multiplicity 1) is
   taken with iota = pi = id (no Hom space). The one-step mode stays as an option and is compared in
   a test.
3. Differential of C (x) X, transported: component sigma -> tau equals pi_tau o Phi o iota_sigma with
   - Phi = f.kron_id(fsize W_k) for every old component f: o -> o' (d_C (x) id);
   - Phi = (-1)^{h_o} id_{L_p} (x) d_j, built column by column with `cube::apply_local_poly`
     (made `pub`; local variables substituted by the top operators of L_p), the Koszul sign of
     cube.rs; the transported components are cached per (p, eps, j), the sign applied per object.
4. Gaussian elimination: pivots are components between equal type and equal t (degree 0,
   asserted to be c * id with c a scalar); Markowitz order (min (|in(B)|-1)(|out(A)|-1), ties by
   ids); update x -> y += -c^{-1} h o g. The complex of step 3 is complete before elimination
   starts, so the streaming rule of §1.4 (C1 I5) is not invoked: every pivot is eligible.
5. Asserts after the crossing: d o d = 0 on the reduced complex (option, on in tests); K_0 class
   (signed multiset of (type, t)) = the decategorified scan of the same stage (schur
   `mul_rungs` + `expand`, as `decat_scan`); #objects >= LB.

## Closure data (exact)

Final labels = lambda (else error: not closable). Per final type p (closed ladder L_p):
- `cclose::rlambda_bases` (exact graded Nakayama with freeness certificate) and the one-point
  determinant criterion `det_criterion` at a seeded random point (both must hold);
- assert: ins and outs degree multisets (normalized) = exponents of c_diag(b_p);
- Gram = `gram_of` (composition); F per final component = `cclose::solve_f` (exact over R^lambda).
`ClosureDataC { lam, objs (Tid, h, t), types (normalized degrees), edges, crossings (eps, m) }`,
`validate()`. Closure per N: `close_exact_colored(data, N, None)` and reduced with marked block 0;
global shift `cclose::global_q`. Labels > N give 0 (one scan for all N).

## Public API

- `xscan::<F>(labels, word: &[Gen], opts) -> Result<XScan<F>, String>`
- `khr_singular_exact::<F>(labels, word, ns, reduced) -> Result<BTreeMap<usize, Poincare>, String>`
- CLI: `krhom N --braid W --colors .. --colored-method singular [--reduced]` (strands <= 2, else
  a clear error).

## Checks against the synthesis

- §9.2: objects (type, h, t), 2-rung canonical ladders, e_p = id; maps exact PMat; f (x) id =
  kron_id; id (x) d via apply_local_poly with untruncated local_cx; decomposition by hom_space and
  scalar intersection forms with rank = [v^s] m_q; ranks add up; scalar pivots, Markowitz, Koszul
  (-1)^h; closure generators by exact Nakayama + one-point determinant; Gram by composition; F by
  exact solve; CL1 exact oracle. All implemented as stated. Deviations: (1) "top-down by l_max" is
  replaced by computing every (q, s) block independently and orthogonalizing with M^{-1} (equivalent:
  the l_max/shift order only matters for the unitriangularity, which is asserted); (2) rung-by-rung
  decomposition of L_p (x) W_k (see step 2; same objects, multiplicities asserted against term_mult).
- §1.4: streaming rule respected trivially (full tensor complex before elimination); idempotents
  never re-solved (e_p = id); scalars read exactly, no point-dependent normalization; component
  degree e = t - t' asserted; Rickard complexes untruncated (n = c + d).
- §2.3: convention A shifts, global_q, labels > N give 0, two primes.
- §4.2: K_0 class asserted per stage; Euler characteristic compared with `decat_chi` in tests (not
  used as evidence about differentials); rank-level oracles (Hopf formula, Wang, cube, cabling,
  duality) in tests/xengine.rs.
