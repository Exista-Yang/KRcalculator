# Knot and link reference tables

Reference data for testing krhom: every prime knot with 3 to 12 crossings, and every
prime link with 2 or more components and up to 11 crossings, in every orientation that
LinkInfo lists. Each record has a PD code, a braid word, and invariants from KnotInfo /
LinkInfo: Jones, HOMFLY-PT, Khovanov homology, signature, Rasmussen s, and linking
numbers.

| file | records | content |
|---|---|---|
| `knots_upto12.jsonl` | 2977 | prime knots `3_1` … `12n_888` (KnotInfo); the unknot is not included |
| `links_upto11.jsonl` | 4188 | 1424 unoriented prime links `L2a1` … `L11n459`, each in all 2^(c-1) orientations listed by LinkInfo |

The format is JSON Lines: one compact JSON object per line, in the source's table order.
`tools/` holds the scripts that built and checked the tables (see *Reproducing*).

**12-crossing links are not included.** LinkInfo stops at 11 crossings, and no source we
found gives 12-crossing links together with invariants. The Knot Atlas `KnotTheory`
package does ship DT codes for them (files `L12A` and `L12N`), but without invariants.

## Counts

Knots per crossing number:

| crossings | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | total |
|---|---|---|---|---|---|---|---|---|---|---|---|
| knots | 1 | 1 | 2 | 3 | 7 | 21 | 49 | 165 | 552 (367a + 185n) | 2176 (1288a + 888n) | 2977 |

Links per crossing number. The *unoriented* row counts links; the *oriented* row counts
records (one record per orientation):

| crossings | 2 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | total |
|---|---|---|---|---|---|---|---|---|---|---|
| unoriented links | 1 | 1 | 1 | 6 | 9 | 29 | 83 | 287 | 1007 | 1424 |
| oriented records | 2 | 2 | 2 | 18 | 20 | 96 | 214 | 914 | 2920 | 4188 |

Links by number of components (unoriented / oriented): 2 components 914 / 1828,
3 components 438 / 1752, 4 components 68 / 544, 5 components 4 / 64.

Every field listed below is present in every record. KnotInfo lists s for all knots up
to 12 crossings, and LinkInfo lists it for all of its links.

## Sources

The files were downloaded on 2026-11-09 and are kept only in scratch, not in this repo.

* **KnotInfo.** C. Livingston and A. H. Moore, *KnotInfo: Table of Knot Invariants*,
  knotinfo.org (accessed 2026-11-09).
  * Download page: https://knotinfo.org/homelinks/database_download.php
  * File: https://raw.githubusercontent.com/knotinfo-org/assets/main/knotinfo_data_complete.xls.zip
    (sha256 `f6abd30a…0ae8`; the spreadsheet inside was last saved 2026-11-07).
* **LinkInfo.** C. Livingston and A. H. Moore, *LinkInfo: Table of Link Invariants*,
  knotinfo.org (accessed 2026-11-09).
  * LinkInfo is now served at https://knotinfo.org/linkinfo/. The old host
    linkinfo.sitehost.iu.edu presents a TLS certificate for a different host name.
  * File: https://raw.githubusercontent.com/knotinfo-org/assets/main/linkinfo_data_complete.xls.zip
    (sha256 `fe4c2972…be34a`; spreadsheet last saved 2026-11-04).
* **Knot Atlas data.** Used only for the `knot_atlas` comparison field and for checks.
  * Package: the `KnotTheory` Mathematica package,
    https://drorbn.net/AcademicPensieve/Projects/KnotTheory/KnotTheory.zip, linked from
    https://katlas.org/wiki/Setup (sha256 `24d8bfe9…a7a102`).
  * Files used: `PD4Knots.m`, `Jones4Knots.m`, `Jones4Knots11.m`, `Kh4Knots11.m`,
    `PD4Links.m`, `Jones4Links.m` and `Kh4Links.m`.

Terms of use: KnotInfo and LinkInfo state no license; they ask to be cited as above.
The Knot Atlas is a public wiki. If this data is used in a publication, cite KnotInfo /
LinkInfo (and, for the Khovanov data, KnotJob by D. Schütz for knots and A. Shumakovitch's
computations for links).

Sources inside the databases, as the KnotInfo / LinkInfo description pages give them:

| data | source |
|---|---|
| knot PD codes | SnapPy, from the Knot Atlas DT codes |
| knot braids | Gittings' minimal braids; Stoimenow for 11-12 crossings; Jones for `10_136` |
| knot braid index | Jones (≤ 10 crossings); Stoimenow (11-12 crossings) |
| knot Khovanov homology | computed with KnotJob (D. Schütz) |
| link braids and braid index | S. Orevkov |
| link Khovanov polynomial | A. Shumakovitch |

## Fields

### Common fields

| field | meaning |
|---|---|
| `name` | Knots: KnotInfo name (`3_1`, `10_132`, `11a_1`, `12n_242`). Links: LinkInfo oriented name such as `L6a4{1,0}`. |
| `crossings` | crossing number |
| `components` | number of components (1 for knots) |
| `alternating` | `true` or `false` |
| `pd` | PD code `[[i,j,k,l], ...]` exactly as the source gives it (KnotInfo `pd_notation`, LinkInfo `pd_notation_vector`) |
| `writhe` | writhe of `pd`, computed here |
| `braid` | `{"strands": n, "word": [...]}`: a braid whose closure is this oriented knot or link |
| `braid_index` | braid index as listed. `braid.strands == braid_index` for every record. |
| `jones` | Jones polynomial, `[[k, c], ...]` meaning Σ c·t^(k/2) |
| `homfly` | HOMFLY-PT polynomial, `[[a, b, c], ...]` meaning Σ c·v^a·z^b |
| `kh` | Khovanov homology, free part (= rational Poincaré polynomial): `[[h, q, rank], ...]` |
| `signature` | signature |
| `s` | Rasmussen s-invariant |
| `determinant` | determinant |
| `knot_atlas` | relation to the Knot Atlas entry of the same name (see *Chirality*) |
| `pd_vs_invariants` | whether `pd` has the listed invariants or those of the mirror (see *Validation*) |
| `braid_vs_invariants` | the same check for the closure of `braid` |
| `issues` | list of flags; empty for almost all records (see *Known issues*) |

### Knot-only fields

| field | meaning |
|---|---|
| `symmetry` | KnotInfo symmetry type: `reversible`, `chiral`, `negative amphicheiral`, `positive amphicheiral` or `fully amphicheiral` |
| `kh_torsion` | torsion of integral unreduced Khovanov homology: `[[h, q, n, m], ...]` means m copies of Z/n in bidegree (h, q). Only Z/2 occurs up to 12 crossings. |
| `kh_reduced`, `kh_reduced_torsion` | the same for reduced integral Khovanov homology. Reduced homology has no torsion up to 12 crossings. |
| `braid_alt` | present for 17 knots that KnotInfo lists with two braids (`10_136`, `11n_8`, `11n_121`, `11n_131` and 13 knots `12n_*`). `braid` is the one with the fewest strands (= braid index); `braid_alt` is the other, which has fewer crossings and one more strand. Both are verified. |

### Link-only fields

| field | meaning |
|---|---|
| `name_unoriented` | e.g. `L6a4`; same as the Knot Atlas link name |
| `orientation` | e.g. `[1,0]`. Entry i refers to component i+2 of the `{0,…,0}` PD: 1 means that component is reversed relative to the Knot Atlas PD. Component 1 is never reversed. Checked for all records: each linking matrix equals the `{0,…,0}` matrix with the corresponding sign changes. |
| `linking_matrix` | LinkInfo linking matrix (diagonal 0). Components are in PD order: component 0 contains the smallest edge label, the next component contains the smallest remaining label, and so on (the same numbering as krhom's `Pd::label_components`). Checked against `pd` for every record. |
| `nullity` | nullity, as listed |
| `kh` | LinkInfo Khovanov polynomial: ranks only, no torsion information. For all 1424 `{0,…,0}` records it equals the Knot Atlas rational Khovanov polynomial (`Kh4Links.m`). |

## Conventions

Every record describes **one oriented knot or link**. The PD code, the braid and all the
invariants refer to the same object; *Validation* explains how this was checked.

### PD codes and crossing signs

PD codes use the KnotTheory / Knot Atlas convention, which `src/input/pd.rs` also uses.
`X[i,j,k,l]` lists the edges counterclockwise, starting from the incoming under-strand,
so the under-strand runs i → k.

Edges are numbered consecutively along each component in the direction of its
orientation. The orientation can also be recovered by walking the under-passages.

A crossing is positive iff the over-strand runs l → j. For example, KnotInfo's
`3_1 = [[1,5,2,4],[3,1,4,6],[5,3,6,2]]` has writhe +3: it is the right-handed trefoil.

For links, `pd` is LinkInfo's code for that orientation. For the `{0,…,0}` orientation it
is *identical* to the Knot Atlas `PD[Link[...]]` code, which was checked for all 1424
links.

### Braids

`word` is a list of signed 1-based generators. +i means the strand at position i passes
OVER the strand at position i+1 (a positive, right-handed crossing); −i is the inverse.
This is krhom's convention: `BraidWord::from_signed(strands, &word)` reads it directly.

All checks closed braids with a port of `krhom::pd::braid_to_pd`. Example: `3_1` has
braid `[1,1,1]`, the right-handed trefoil.

### Jones polynomial

The Jones polynomial is stored as `[[k, c]]` meaning Σ c·t^(k/2), with V(unknot) = 1 and
the skein relation t^(-1)·V(L+) − t·V(L−) = (t^(1/2) − t^(-1/2))·V(L0).

The right trefoil is `[[2,1],[6,1],[8,-1]]` = t + t^3 − t^4. This is KnotInfo's and
KnotTheory's convention, and the same `(2·exponent, coefficient)` encoding as the
`Poly` in `tests/pd.rs`.

LinkInfo writes the Jones polynomial in x = t^(1/2), so k is the exponent of x.

### HOMFLY-PT polynomial

The HOMFLY-PT polynomial is stored as `[[a, b, c]]` meaning Σ c·v^a·z^b, with
P(unknot) = 1 and the skein relation v^(-1)·P(L+) − v·P(L−) = z·P(L0) (KnotInfo /
LinkInfo / Morton). Related conventions:

* The Jones polynomial is V = P(v = t, z = t^(1/2) − t^(-1/2)).
* KnotTheory's `HOMFLYPT[K][a,z]` uses a = v^(-1).
* The right trefoil is 2v^2 − v^4 + v^2·z^2.

### Khovanov homology

Gradings are Bar-Natan / KnotTheory / KnotInfo: h is the homological degree (the
exponent of t) and q is the quantum degree. The right trefoil is
q + q^3 + t^2·q^5 + t^3·q^9, with Z/2 at (3, 7). In this convention:

* unreduced: Σ (−1)^h·rank·q^j = (q + q^(-1))·V(t^(1/2) = −q);
* reduced (knots): Σ (−1)^h·rank·q^j = V(t^(1/2) = −q).

Over F_p, the rank in bidegree (h, q) is

    rank_free(h, q) + #{Z/p^e at (h, q)} + #{Z/p^e at (h+1, q)},

because the differential raises h. For example, the trefoil over F_2 has F_2 in both
(2, 7) and (3, 7).

### Signature and s

* Signature: the KnotInfo / LinkInfo convention, in which positive knots and links have
  **negative** signature: σ(right trefoil) = −2 and σ(positive Hopf link) = −1.
* Rasmussen s: s(right trefoil) = +2, so s = −σ for alternating knots.
* For links, LinkInfo's s is normalised so that the positive Hopf link has s = +1. In
  every record the listed s satisfies Kh^{0, s±1} ≠ 0.

### Relation to krhom's gl_N gradings (docs/dev/PLAN.md)

docs/dev/PLAN.md puts positive knots in **negative** q-degrees: krhom's N = 2 result is Khovanov
homology with q ↦ q^(-1), and homological degrees are unchanged. So:

* **Comparing with `kh`.** krhom's N = 2 homology in bidegree (h, q) corresponds to the
  table entry (h, −q).
* **Euler characteristics in krhom's convention.**
  * Reduced: χ = P(v = q^(-N), z = q − q^(-1)).
  * Unreduced: χ = [N]·P(v = q^(-N), z = q − q^(-1)), where [N] = (q^N − q^(-N))/(q − q^(-1)).
* **Euler characteristics in the Bar-Natan direction.** The corresponding formula is
  (−1)^(c−1)·P(v = q^N, z = q − q^(-1)).
* **Checks on the right trefoil.** In krhom's convention the reduced value is
  q^(-2(N-1)) + q^(-2(N+1)) − q^(-4N). At N = 2 the unreduced value is
  q^(-1) + q^(-3) + q^(-5) − q^(-9), which matches docs/dev/PLAN.md.
* **Rasmussen s.** In krhom's convention, s is read off at q-degrees −s ± 1.

## Chirality and the Knot Atlas

KnotInfo and the Knot Atlas sometimes choose opposite mirror images for the same knot
name. The `knot_atlas` field records, for knots up to 11 crossings, how this record (in
KnotInfo's chirality) relates to the Knot Atlas / KnotTheory knot of the same name
(`Knot[n,k]`, `Knot[11,Alternating,k]` = `K11a k`):

| value | meaning | count |
|---|---|---|
| `same` | the Knot Atlas Jones polynomial is the same (or, if the Jones polynomial is symmetric, the signature or Khovanov polynomial is) | 638 |
| `mirror` | this record is the mirror image of the Knot Atlas knot | 139 |
| `amphichiral` | KnotInfo lists the knot as amphicheiral | 20 |
| `undetermined` | the knot is chiral, but every invariant stored here is mirror-symmetric | 4 (`10_48`, `10_71`, `10_91`, `10_104`) |
| `null` | 12-crossing knots: the Knot Atlas has no invariant data for them | 2176 |

All 552 knots with 11 crossings are `same`. Among knots up to 10 crossings, 139 are
`mirror`; for example `3_1`, `5_1`, `5_2` and `6_1`. The Knot Atlas PDs in krhom's
`tests/pd.rs` are therefore often the mirror images of the records here.

For links, `knot_atlas` is `identical_pd` for the `{0,…,0}` orientation: the PD code
equals the Knot Atlas code character for character, so the chirality is the same. For the
other orientations it is `reoriented`.

## Validation

`tools/build.py` recomputed the following from the downloaded data, using independent
Python code: a Kauffman bracket by planar contraction, Gordon–Litherland signatures from
both checkerboard colourings, a HOMFLY-PT skein algorithm on descending diagrams, linking
numbers, and krhom's braid-closure convention.

1. **Source formats agree.** The polynomial strings and the vector encodings agree for
   the Jones, HOMFLY and Khovanov data. KnotInfo's Khovanov-over-Q data equals the free
   part of its integral data.
2. **The listed invariants agree with each other.** For every record:
   * HOMFLY → Jones (via the formula above);
   * Khovanov Euler characteristic → Jones, for both unreduced and reduced;
   * |V(−1)| = determinant;
   * Kh^{0, s±1} ≠ 0.
3. **PD versus listed invariants (`pd_vs_invariants`).** The checks are the Jones
   polynomial, the HOMFLY-PT polynomial, the signature and, for links, the linking matrix,
   each computed from `pd`. The value is:
   * `same` if they match the listed values and are not compatible with the mirror image,
     or if the knot is amphicheiral;
   * `same_or_mirror` if they match the listed values but cannot tell the mirror apart;
   * `mirror` or `conflict` otherwise. Neither occurs.

   Results: **2949 knots `same`, 28 `same_or_mirror`; 4162 links `same`, 26
   `same_or_mirror`.** For every `same_or_mirror` record, every chirality-sensitive
   invariant stored here is mirror-symmetric: Jones, HOMFLY, Khovanov homology with
   torsion, reduced Khovanov homology, σ = 0 and, for knots, s = 0. The PD code is
   therefore a valid source for all the stored values whichever chirality it is, with one
   exception: for 14 link records the stored `s` is −1, not 0, and no stored or
   recomputed invariant fixes whether that value belongs to the PD's chirality or to its
   mirror. These are `L10n36{0}`, `L10n36{1}`, `L10n59{0}`, `L10n59{1}`, all 8
   orientations of `L10n107`, `L11n247{0}` and `L11n247{1}`; treat their `s` as
   unverified. (Found by the independent validation pass.) The records are:
   * knots `10_48`, `10_71`, `10_91`, `10_104`, `12a_126`, `12a_132`, `12a_214`,
     `12a_222`, `12a_453`, `12a_674`, `12a_697`, `12a_943`, `12a_1011`, `12a_1081`,
     `12a_1088`, `12a_1199`, `12a_1222`, `12a_1231`, `12a_1235`, `12a_1250`,
     `12a_1258`, `12a_1283`, `12n_250`, `12n_278`, `12n_490`, `12n_586`, `12n_610`,
     `12n_695`;
   * links, in all their orientations: `L6a4` (the Borromean rings), `L10a140`,
     `L10a151`, `L10n36`, `L10n59`, `L10n107`, `L11n247`.
4. **Braid closure versus listed invariants (`braid_vs_invariants`).** The same checks
   are run on the braid closure; linking numbers are compared up to a permutation of the
   components. Results:
   * knots: 2949 `same`, 28 `same_or_mirror`; all 17 `braid_alt` braids are `same`;
   * links: 4145 `same`, 42 `same_or_mirror`, and 1 `conflict` (see *Known issues*).
5. **Knot Atlas self-checks.**
   * Jones polynomials computed from `PD4Knots.m` and `PD4Links.m` equal `Jones4Knots.m`
     and `Jones4Links.m` for all 249 knots and 1424 links. This confirms the Jones code,
     the sign rule and the table ordering.
6. **Spot checks against literature values.** `tools/verify_output.py` checks:
   * the table counts;
   * that PD labels are well formed;
   * the right trefoil: Jones, HOMFLY, integral Kh and reduced Kh;
   * 4_1, 5_1, and `12n_242` (the (−2,3,7) pretzel: positive braid, s = 10, σ = −8);
   * the Conway and Kinoshita–Terasaka mutants (equal Jones and HOMFLY, s = 0);
   * the positive and negative Hopf links, the Borromean rings, and T(2,4).

Records can therefore be filtered with `pd_vs_invariants == "same"`, or
`!= "conflict"`, which keeps all records.

## Known issues (`issues` field)

* `L10n99{0,0,0}`: **LinkInfo's braid is for a different orientation.** Its closure has
  the linking matrix of `L10n99{1,0,0}` / `L10n99{1,1,1}`, not that of `{0,0,0}`; the PD
  code is fine. Only the linking numbers detect this: `{0,0,0}`, `{1,0,0}` and `{1,1,1}`
  share the same Jones, HOMFLY and σ, so a Jones-only check passes. The record carries
  `"braid_closure_matches": ["L10n99{1,0,0}", "L10n99{1,1,1}"]` and
  `braid_vs_invariants: "conflict"`. Use the PD code for this record.
* `L11n459{1,1,1}`: the LinkInfo HOMFLY string has a stray trailing `}`. It was parsed
  without it, and the result agrees with the vector form and with the PD
  (`source_typo_homfly_string_trailing_brace`).
* `pd_chirality_undetermined` / `braid_chirality_undetermined`: see *Validation* item 3.
* 17 knots list two braids; see `braid_alt`.

## Reproducing

The builder is pure Python 3 (standard library only). The `.xls` reader needs `xlrd`.
The full build takes about 5 minutes and about 350 MB of memory.

```
python -m venv venv && venv/bin/pip install xlrd==2.0.1
# download and unzip the three archives listed under Sources into separate directories
venv/bin/python -I tools/extract.py knotinfo_data_complete.xls work/knotinfo_sel.json knot
venv/bin/python -I tools/extract.py linkinfo_data_complete.xls work/linkinfo_sel.json link
venv/bin/python -I tools/build.py work <unzipped KnotTheory dir>/KnotTheory out   # writes the two .jsonl files
venv/bin/python -I tools/verify_output.py out
venv/bin/python -I tools/compare_kh_links.py out/links_upto11.jsonl <unzipped KnotTheory dir>/KnotTheory
```

`tools/krlib.py` contains the shared routines: PD orientation, Kauffman bracket,
Gordon–Litherland signature, HOMFLY-PT skein algorithm, braid closure and a polynomial
parser.
