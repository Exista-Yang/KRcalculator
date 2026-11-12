# Review wave 4: M4c B1 closed forms (src/leafgen.rs, tests/leafgen.rs)

Branch m4c-gen f5fe90c, worktree .claude/worktrees/revgen. Diff reviewed: 595068c..HEAD (git diff, read-only).
No cargo runs. Monomial counts below come from a 1-second generating-function script. The script was cross-checked
against the logged counts: (3,2,3) psi_A gives 1590 and far (3,3,2,1) gives 4644, both matching the logs.

## Verdict on the lens questions

**Trust path interpolate -> Added -> certify_added: sound for genuineness.**
- certify_added's exact checks (tmpl.rs:417-433) are enough by themselves: homogeneity of every term at the forced
  degree, plus exact `phi*M_src = M_tgt*phi` for every top operator.
- The PMat is a matrix over R_bottom, so left linearity is automatic. A closed form that is not a bimodule map is
  therefore refused exactly, not probabilistically.
- Caveat, not a defect: the "independent-seed polynomiality certificate" inside certify_added (tmpl.rs:435-438) runs
  on `localize(t)`, which equals `L_t * Phi * L_s^-1` for a polynomial Phi (tmpl.rs:563-567).
  - poly_certificate recomputes `L_t^-1 * f * L_s = Phi` exactly (tmpl.rs:762), so it cannot fail for an Added
    template.
  - It adds no independent evidence, and it costs about as much as the interpolation (far 3321: 399 s of 819 s).
- The only check that ties the certified PMat to the closed form is build()'s single fresh-point comparison.

**build() fresh-point re-check: present.**
- leafgen.rs:437-439 compares at `points(kw, seed^0x9e37, 1, 7)`. `point_seed` is injective in (seed, id), and
  `seed^0x9e37` cannot equal `seed+0x1eaf0000+a` or `seed^0xce570000b1b1+a`. So the point is fresh.
- `build_in` (leafgen.rs:653-661), which braid_dag uses, has no such re-check. This is harmless today:
  - interpolate is exact whenever f is polynomial;
  - u is always even (Pull 0, Push 2(ss-ts), Cap -2r(g-r), Cup 2ra, Rex/Swap 0) and fdeg is even, so the
    `nm == 0 -> continue` branch (leafgen.rs:378), which does not check that the entries are zero, is unreachable
    for these kinds.
  - Not reported as a finding.

**validate_assoc: cannot accept a wrong constant.**
- The recorded scalar is never read in evaluation. `certify_assoc` (smap.rs:831-834) only records presence, through
  `elem_rec_certified` (smap.rs:799-805).
- Eval uses the node's own scalar times `elem_winloc(Assoc)` (smap.rs:1183, 1196, 1245).
- validate_assoc runs poly_certificate on exactly that constant-1 closed form, at degree 0 with equal shifts
  (leafgen.rs:639-646). A pass means the map is genuine, with 1(x) -> 1(x).
- Multistep equality is only an extra filter.
- Key consistency: braid_dag's (me, mpe) equals the window that `elem_spec` derives. There is no common suffix
  (F0.y vs F1.x), and the window spans uprights 0..2. A mismatch would only make the node uncertified, which is
  conservative.

**Euler sign convention: consistent with smap::elem_winloc.**
- `euler` = prod over rungs of e(moved, giver \ moved) (leafgen.rs:92-102).
- Push into [] gives Trace = 1/e(S, A\S). Push from [] gives Eta = e(B, A) and EtaL = e(A, B).
- Asserted entrywise in `elem_kinds_are_instances` (tests:143-180).

**Rex normalization: correct when Y1 is reduced.**
- Entry = e_X(c) / (kappa * e_Y1(p)), which equals Pull(X->Y2) o Push(Y1->X) / kappa.
- 1(x) -> 1(x) by construction.
- psi o phi = id is checked in r1_switchback.
- For non-reduced Y1, see finding 2.

**Exact hom_space comparison: done where claimed, but weaker than the report states (finding 3).**

## Findings

### 1. (major) Flow-count gates do not bound interpolation/certificate memory; braid_dag's Pull fallback OOMs at cells (3,2,3)
- **Gates are on flows, cost is in monomials.**
  - `interpolate` admits any window with at most MAX_INTERP_FLOWS = 1200 flows (leafgen.rs:51, 311).
  - It then stores `npts = max_monos + 2` dense nt x ns matrices (`phis`, leafgen.rs:339-358).
  - It also inverts an nm x nm Vandermonde, built as an nm x 2nm augmented rref (leafgen.rs:381-382; local.rs:205-216).
  - poly_certificate inside certify_added repeats the same with `max_monos + 4` points (tmpl.rs:730-736, 787).
- **braid_dag reaches the fallback at (3,2,3).**
  - When validate_assoc refuses (more than 500 flows), braid_dag falls back to a Pull template M -> M'
    (leafgen.rs:691-695).
  - At (3,2,3), M = [E0.2, E1.3, E0.3, F0.3, F1.3, F0.2] is flag-type with 560 flows, which is at most 1200, so it is
    admitted.
  - Its max fdeg is 2*(15+6) = 42. Ring weights for labels (3,2,3) are (2,4,6,2,4,2,4,6), which gives 39 148
    monomials at degree 42.
- **Cost of that fallback:**
  - phis: 39 150 x 560 x 560 x 8 B, about 98 GB;
  - Vandermonde: about 24.5 GB;
  - roughly 2e13 field operations for the 560^3 inversions alone;
  - then certify_added repeats all of it.
- (3,3,2) and (2,3,3) are the same size (560 flows).
- **Over-cap windows still pay a heavy build before the refusal.** interpolate builds both exact SBims
  (leafgen.rs:307-308) before the flow gate at line 311. At (3,3,3), with 1680 flows, the fallback first runs
  `Expr::sbim` on the 9-root window, which waits on GrTable(3,9) (more than 40 min per the report), and only then
  returns Err.
- **Failure:** braid_dag(3,2,3,...) is a call B2 will make for 3-colored gl6. It runs about 11 min of certified
  pieces, then allocates more than 100 GB and crashes or swaps the laptop, instead of returning Err.
- The report itself says (sec. 4) that cost is driven by monomials and that (3,2,3) is "not certifiable". The code
  does not enforce either.
- **Fix:**
  - compute max_monos from fdeg ranges and ring weights before any SBim build;
  - gate interpolate and poly_certificate on npts*ns*nt and on nm^2;
  - make braid_dag return Err instead of attempting an over-budget Pull fallback;
  - move the flow gate above the sbim() calls.

### 2. (minor) Rex accepts non-reduced equal-shift pairs, then panics inside the Result-returning build path
- ClosedForm::new for Rex only checks equal shifts and equal supports (leafgen.rs:167-175). Its error text assumes
  "different shifts implies not both reduced", but the converse does not hold. Example: Y1 = Y2 = a non-reduced word,
  or two non-reduced words of equal degree.
- For non-reduced Y1, kappa(c) = sum over p in c of 1/e_Y1(p) vanishes identically: it is the pushforward of 1 along
  a map with positive-dimensional fibres.
- `winloc` then hits `assert!(!s.is_zero(), ...)` (leafgen.rs:240) inside the interpolate closure. That panic goes
  through `build` / `build_in`, which are documented as returning Err for maps that are not genuine.
- **Failure:** a B2/B4 planner probes Rex(W -> L) on a window whose differing parts are non-reduced with equal
  shift (report sec. 6.4 proposes this use). The process aborts instead of refusing.
- **Fix:** return Err from winloc's caller, or pre-check with a single point in ClosedForm::new (sum is zero means
  refuse).

### 3. (minor) Exact-window tests do not assert the "dim 1, so the closed form is *the* map" claim
- `exact_check` returns (dim, in_span) (tests:37-54). `check` asserts only in_span and prints dim (tests:68-73). A
  window above EXACT_CAP is skipped with an eprintln only.
- In-span adds little: a template that passed certify_added is genuine, so it lies in Hom^deg anyway. The check tests
  hom_space completeness more than the closed form.
- R1's substantive claim is "every exact window has dimension 1, so each closed form is the unique map up to scalar"
  (B1_report sec. 0 and sec. 3.1, and the "exact (dim 1)" column in sec. 2). It rests on log output, not on
  assertions.
- **Failure:** a regression where a window has dim 2, or where a case silently goes above the cap, still passes the
  default suite.
- **Fix:** `assert_eq!(dim, 1)`, and assert that exact_check returned Some for the cases the report lists as exact.

## Checked, no defect
- **certify_added degree formula** (tmpl.rs:416-420) is consistent with interpolate (leafgen.rs:318-324) and
  hom_space.
- **Interpolation layout:** `fd[q][pp] = d[tpos[q]][spos[pp]]` and `Phi = L_t^-1 * fd * L_s` invert localize's
  `L_t * Phi * L_s^-1` (tmpl.rs:567).
- **Interpolation check points:** every degree class is checked at at least 2 extra points, and negative-degree
  entries are checked to be zero at every point (leafgen.rs:353-357, 388-399).
- **Window sides:** the SBim and fibre basis in interpolate are the same as TmplLib::side's (`w.expr().sbim()`,
  tmpl.rs:235).
- **Pull/Push coset ranks:** these are shared between sfi and tfi (same bottom/top CosetSpace). `range(c)` is empty
  outside the support.
- **Cap/Cup:** m[0] == m[1] holds exactly on the identity coset, and the degrees ra - r(g-r) are correct. The entries
  are certified exactly in any case.
- **Swap:** with disjoint uprights the permuted flow always exists (leafgen.rs:200-202, 272).
- **canonical_max:** flag-type, and handles zero labels (braid windows (x+y, 0, z)).
- **Arena flags:** certify_template refreshes the flags. smapio re-certifies Added templates on load (smapio.rs:373-379)
  and refuses to write uncertified ones.
- **braid_case:** composite == direct Rex at every coset (one point). The direct Rex was built with build(), so it has
  the fresh-point check.
- **Ignored runs:** the (2,2,2) direct comparison now runs only in the ignored braid_large; braid_222 correctly says
  "composite checks only".
