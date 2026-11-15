# Review: colored-m5 (bf326b7): cache, serialization and tests

Scope: src/sindec.rs (layer A, DecHom, restrict), src/smapio.rs, tests/sindec.rs, diffed against
dev ce91f0e. Read-only review: no cargo, no binaries. All line numbers refer to the revm5 worktree
unless marked otherwise.

## Verdict

**Layer A and smapio are sound.** I found no way for a stale or foreign file to be loaded as valid:

- A file is accepted only on an exact match of its full header line (sindec.rs:1110). That line holds
  every stamp and the full key (1074-1076, 1102-1104). The file name is only fnv64(key) (1099-1101),
  so a collision makes two keys overwrite each other's file and recompute. It never causes a wrong load.
- The key contains the λ, the field (F::name() = "F_<P>", field.rs:125), both table strings, the
  rung words of L_p and L_q, R, s and mult (1077-1096).
- The body has its own header, holding the field, λ and TMPL_ALG (smapio.rs:40-42, 300-304).
- Templates are re-obtained on import and compared against their content hashes (smapio.rs:322-344).
  Added templates are re-certified through certify_added (smapio.rs:373-378; tmpl.rs:392-421: degree,
  exact commutation, polynomiality).
- Loaded maps then go through the same checks as fresh ones:
  - dechom_postcheck: count, ends, info-degree and D14 scaling (688);
  - the form at rep(q) at 2 points, its invertibility, and assert_genuine (645-657);
  - Neumann orthogonalization with its unitriangularity check, then Π∘Ι = id at 2 points on every
    coset (674-676);
  - for closures, certify_closure (1044).
- These checks depend on content, not on where the file came from. A file with a correct header but
  altered genuine content still gives a valid decomposition, or else it fails closed.

Other points checked and found correct:
- Writes are atomic (tmp file, then rename: 1134-1136; tmpl.rs:220-246).
- Signed coefficients round-trip exactly for P < 2^32 (field.rs:67-82, 119-121).
- λ = 1 234 567 has multiplicative order 2^30−1 mod 2^31−1, and p−1 (primitive) mod 4 294 967 291
  (computed). The scaling check is therefore not weak.

**Negative DecHom tests fail for the intended reason.** Each test asserts both `guard` and a text that
only the intended check can produce:
- "wrong mult": sindec.rs:323 only;
- "wrong lq": 306 only;
- "decomposable": 315;
- "wrong degree": 337, from the info-degree comparison in the postcheck.

**No existing test was weakened.** `git diff ce91f0e..HEAD -- tests/` on colored-m5 only adds
tests/sindec.rs. On colored-2s, tests/xengine.rs only moves parse_poincare, wang_ref and hopf_formula
into tests/common/mod.rs, unchanged.

The findings below are test gaps or fail-closed and performance issues. None produces a wrong result
today.

## Findings

### 1. [minor] D5 compares only the objects of the final minimal complex, not the differentials
- **Where:** tests/sindec.rs:176-199 and 370-418.
- **What the SPEC asks for.** M5a imports D5 from the old SPEC (SPEC_3C_GL6_rev1 §2 M5a). The old SPEC
  D5 (SPEC_M4_M6a.md:1179) requires:
  - the per-(p, k) summand lists to equal xengine's;
  - every Rickard component (σ, τ) to have the same rank as xengine's dloc at every coset.
- **What the test does.**
  - It runs a test-only mini scan, eliminating isomorphisms at pts[0], and compares the sorted
    (table, h, t) list with `xscan(..).final_objects`.
  - The 3-strand test compares only object lists, between σ1σ2σ1 and σ2σ1σ2.
  - There is no d∘d = 0 check anywhere in the mini scan.
- **What this misses.** Components between non-isomorphic summands are never compared, and they are
  what determine homology after closure.
  - Example: a wrong scalar or sign in one Rickard component (σ, τ), say a wrong normalization of one
    d_j template.
  - Ends and degrees still check out, so Arena::add and the degree assert at sindec.rs:985-990 accept
    it.
  - The minimal complex can keep the same objects with a different differential. D5 passes, and M6a's
    ranks come out wrong in the same way at both primes.
- **Fix.**
  - In the unpruned mini scan, assert d∘d = 0 at 2 points after every crossing.
  - Compare per-(p, k) summand lists with xengine, as the SPEC requires.
  - Compare the remaining components' coset ranks with xengine at each stage, at least for (2,1),
    (2,2) and (3,2).

### 2. [minor] d7_restrict checks only the lengths of Rickard lists; the positional remap is untested
- **Where:**
  - tests/sindec.rs:614-622 asserts `r.rickard_comp_at(..).len() == x.len()`;
  - DecData::restrict re-indexes (si, ti) through `remap` at sindec.rs:1281-1283;
  - no test ever uses restricted data in a scan.
- **Why it matters.** A remap bug that attaches components to the wrong kept summand would pass:
  - the summand lists are compared as sorted multisets (tests/sindec.rs:518-527, 609-613);
  - the lengths of the Rickard lists are unchanged.
- **Context.** This is exactly the positional-cache failure class of C_sound M3. The current code
  looks correct; the test would not catch a regression.
- **Fix.**
  - Compare the (si, ti) pairs, mapped to (q, s) of the source and target summands, together with the
    blocks at 2 points, against the fresh prune-n build.
  - Or run mini_scan on the restricted data and compare it with a fresh build.

### 3. [minor; outside this review's scope, latent] `rickard_comp` does not identify the strand position on 3 strands
- **Where:**
  - `SIndec::rickard_comp` (sindec.rs:998-1004) checks only that (mu[0], mu[1]) == (c, d) and always
    uses position 0;
  - `DecData::rickard_comp` (sindec.rs:1204-1206) ignores c and d completely.
- **Failure.** On (3,3,3), or on (1,1,1), the labels at positions 0 and 1 are equal.
  - A crossing at position 1 requested through the §5.1 form returns the position-0 components with no
    error.
  - SPEC M6a-1b wires crossings "via dec_word and rickard_comp".
  - The (si, ti) indices refer to the decompositions of the position-0 terms. Placed into the
    position-1 decompositions, they either panic or put wrong differentials in place.
- **Fix.** Return Err(Invalid) from both functions when lam.len() > 2, so callers must use
  rickard_comp_at. Also check (c, d) in the frozen version.

### 4. [minor; possible hang at (3,3,3) scale, not measured] Every layer-A save rewrites the whole template cache
- **Where:** `layer_a_save` calls `self.ar.lib_ref().save()` (sindec.rs:1138).
- **Cost.** `TmplLib::save` rewrites one file for every hom and Rickard call the library has seen,
  whether loaded or computed (tmpl.rs:481, 537, 578-593).
  - The cost is O(saved targets × library calls) file rewrites per process.
  - The (3,3) run is small: 80 layer-A files.
  - `SIndec::prebuild` (1055-1066) over a (3,3,3) N=6 plan writes 1 592 targets plus closures in one
    process. If Search touches thousands of hom windows, the rewrites could dominate.
- **Fix.** Save only the calls made since the last save, or save once at the end and at checkpoints.

### 5. [minor, fails closed] A cache write failure aborts a decomposition that has already been certified
- **Where:** `self.layer_a_save(&key, &all)?` (707) and `self.layer_a_save(&key, &ms)?` (1036). Both
  propagate DecErr::Cache.
- **Causes:**
  - an I/O error (a full disk in a 5-h cloud session);
  - an export error (smapio.rs:152-158 refuses non-integer coefficients, so F = Q fails as soon as any
    Lin coefficient is a fraction).
- **Effect.** The provider work and its certification are thrown away, and the run stops.
- **Fix.** Log the error and continue without caching.

## Notes (not defects)

- **Files are saved before M5's own checks, and failures on load are final.**
  - Pairs are saved at 707, before M5's own form, orthogonalization and Π∘Ι checks (645-676). Closures
    are saved at 1036, before certify_closure (1044).
  - Both were already certified by the provider (homprov.rs:1322, 1371, 2615, 2754, 2851), so
    this only matters after a bug.
  - Such a file is reloaded on every later run and fails as Check/Cache, for example at line 688. It is
    never rejected and recomputed.
  - The provider stamp (436) has no code version, so fixing the provider does not invalidate the file.
    It has to be deleted by hand.
- **The layer-A stamp breaks reuse across sessions.**
  - The stamp includes the machine-dependent `slowdown` and a hash of the cost-table file (427, 436).
  - Layer A built on one machine, or before a recalibration, is rejected in full elsewhere, although the
    content is re-certified on load and provenance does not affect soundness.
  - This works against D12's archive reuse. Consider leaving the stamp only in a comment, or dropping
    slowdown from it.
- **D14 asks for prune_n in every cache key; layer A omits it.** This is correct: layer A is N-free
  (R2 §4.1), and the targets are the same at every N.
- **Untested paths and stamps.**
  - No test covers a body with a valid header but tampered content, a template-hash mismatch, or a
    foreign field or provider stamp. Only the stale TMPL_ALG header is tested (tests/sindec.rs:581-593).
  - d7 does not assert that the "A" and "R" template paths of smapio were exercised. Whether they were
    depends on which tier served (3,2).
  - `SIndec::new` stamps "unstamped" (475). Production drivers should use `with_stamp`.
- **D5's reference engine changes when colored-2s merges.**
  - colored-2s changes `XOpts::default()` to `stream`, `compact` and `rung_transport` = true
    (rev2s src/xengine.rs:122-126).
  - tests/sindec.rs:197 (`..Default::default()`) then compares against the M2S′ streaming engine,
    not the full-build engine that was validated on the tables.
  - Pin `stream: false, compact: false, rung_transport: false` for the oracle, or run both.
