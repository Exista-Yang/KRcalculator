ADVERSARIAL VERIFICATION OF THE KNOTJOB REPORT

Overall: the report is accurate. I found no fabricated quotes or numbers. There are 2 substantive errors (S1, S2 below) and about 12 omissions or nuances a developer needs. The live Durham site still returns ECONNRESET, both to WebFetch and to curl, and WebFetch cannot reach web.archive.org. So I checked the page through the Wayback CDX index (via curl) and through the saved snapshots.

**Evidence base**
- **Zip:** the local KnotJob.zip has SHA-1 (base32) CDSZXL25MZE6YVIWRJ7I5SQR6S5N2GOQ. This is identical to the CDX digest of the only archived copy, https://web.archive.org/web/20250319065211id_/https://www.maths.dur.ac.uk/users/dirk.schuetz/KnotJob.zip. Every entry in the zip is dated 12-14-2023. Below, "[Z:path]" means a file inside that zip.
- **Page snapshots:**
  - [P23] https://web.archive.org/web/20240224091119/https://www.maths.dur.ac.uk/users/dirk.schuetz/knotjob.html
  - [P25a] same URL at /20250530040959/
  - [P25b] same URL at /20250913200424/
  - [P26] same URL at /20261002121751/ (a revisit with the same digest as the 200 captures of 20250819 and 20260904)
- **Downloaded file:** /private/tmp/claude-501/-Users-existair-Documents-KRcalculator/bf6f8870-08ad-436b-9ff7-0414655cc59a/scratchpad/verify_dl/ki_desc.html

== PROVENANCE / VERSIONS ==
- **Live site resets connections:** CONFIRMED. My own WebFetch and curl attempts failed the same way.
- **Only one archived KnotJob.zip:** CONFIRMED. A CDX prefix search shows only KnotJob.zip at 20250319065211. The old path ~dma0ds/KnotJob.zip returns a 301 to it.
- **Zip is the Dec 2023 build, Java 17, last modified 14.12.2023, source says "gray":** CONFIRMED. Sources: [P23]; [Z:knotjob/KnotJob.java] ("This version is gray"); [Z:knotjob/Options.java] (versionColor = Color.GRAY).
- **sl(3) first appears on the page modified 02.04.2025:** CONFIRMED [P25a], with two nuances:
  - No snapshot exists between 2024-02-24 and 2025-05-30, so the exact date sl(3) was added is unknown.
  - Missing from the report: the April 2025 page required **Java 21**, and its version colour image was "aqua". The Aug 2025 page requires Java 23 and uses "chart" (chartreuse). So there were at least two builds after the inspected one.
  - The Dec 2023 source already contains commented-out sl(3) hooks: calculateSlThreeHom and SlThreePolCalculator in [Z:knotjob/KnotButtons.java] and [Z:knotjob/Knobster.java].
- **Current page modified 15.08.2025, Java 23, unchanged through 2026-10-02:** CONFIRMED [P25b][P26].
- **Download links:** CONFIRMED [P26]. The page also links KnotJob_j8_mag.zip (a Java 8 build, "magenta"), KnotJob.tar and XKnotJob.tar.gz.
- **No public git repository; GPLv3 gpl.txt in the zip:** CONFIRMED for gpl.txt [Z:gpl.txt]. Web searches found no repository, though absence cannot be fully proven.
- **Old page ~dma0ds/knotjob.html:** CONFIRMED that it existed. It has returned a 301 to the users/dirk.schuetz page since 2022 (CDX).

== (1) WHAT IT COMPUTES ==
- **Home-page feature list:** CONFIRMED verbatim [P25b/P26]. That covers even/odd Kh; s "for fields up to characteristic 211 and an integral version, both knots and links"; LS Sq^1/Sq^2 even and odd; DLS refinements; Sq^2 even and odd; BLT spectral sequence up to char 211; sl(3) homology and spectral sequences; sl(3) s-invariants; Alexander, Jones, signature, determinant. "both knots and links" first appears in the 15.08.2025 version.
- **Vintage KnotJob, Lewin Strauss mod-4 methods, XKnotJob = earlier GUI version:** CONFIRMED [P23].
- **README -kb0 / 1 / q, including "q^2 is less than 2^31":** CONFIRMED [Z:README.TXT].
- **47 primes from 2 to 211:** CONFIRMED [Z:knotjob/Options.java]. Each prime's maximum power is limited by an int-overflow test divided by 2.
- **GUI coefficient choices:** CONFIRMED with a nuance [Z:knotjob/dialogs/KhovDialog.java, Knobster.java].
  - They are two separate menu items, not one dialog with four choices.
  - "Char 0 Khovanov Cohomology" offers radio buttons Integral / Rational / Local ("Choose primes to be inverted").
  - "Modular Khovanov Cohomology" offers a prime spinner and a power spinner.
  - Local is GUI-only: the command line has no way to express localization.
- **-s0 / -s<p> with p < 212:** CONFIRMED [Z:README.TXT].
- **Odd homology is reduced-only on the command line; GUI has "odd" and "odd reduced":** CONFIRMED [Z:README.TXT; Knobster.java]. The GUI buttons are "even", "even reduced", "odd", "odd reduced".
- **Graded s printed as "s.-s" / "s(2).-s":** CONFIRMED [Z:Documentation/Graded s-invariants.rtf].
- **Cohomology quote:** CONFIRMED [Z:Documentation/Khovanov cohomology.rtf].
- **Steenrod squares and Bocksteins:** CONFIRMED [Z:Documentation/Stable Homotopy Types.rtf]. The Bocksteins are for 0→Z/2→Z/2^{k+1}→Z/2^k→0 and the dual sequence.
- **-kx only in the source:** CONFIRMED [Z:knotjob/KnotJobCommand.java]. Its output label is "Extortion Order mod c :".
- **sl(3) is built on Mackaay–Vaz:** PARTLY CONFIRMED. The abstract of arXiv:2504.07806 says the paper uses Mackaay–Vaz universal sl3 homology (https://arxiv.org/abs/2504.07806). The paper does not say how knotjob implements it.
- **s_sl3 computed for p = 0, 2, 3, 5, 7; spectral sequences X^3−X−w and X^3−X^2 in char 0, and char 3; the "knotjob can calculate wE_k^{i,j} …" quote:** CONFIRMED, arXiv:2504.07806 §8.1–8.3. The char-3 polynomial is X^3−X^2−w with w = ±1.
- **No colored homology and no gl_N / KR for N ≥ 4:** CONFIRMED as an absence. Nothing on the page, in the docs or in the source mentions them; the only sl_n item is the vintage matched-diagram homotopy type [P26].

== (2) INPUT FORMATS ==
- **Torus links:** CONFIRMED [Z:Documentation/New Links.rtf]. p is the strand count, p ≤ 10, q ≤ 60, both positive; use a mirror for negatives.
- **Pretzel:** CONFIRMED. The doc's wording is "Up to 12 pretzels", not "tassels".
- **Braid quote, DT, Gauss, PD and ADT examples, the operations, unlink ≤ 10, Change Orientation:** CONFIRMED [New Links.rtf].
- **Whitehead Doubles button:** CONFIRMED [Z:knotjob/Knobster.java].
- **Strand count = index of the largest letter + 1; fewer than 2 or more than 10 strands rejected:** CONFIRMED [Z:knotjob/links/LinkCreator.java, obtainBraid/largestGen].
  - Missing: the word must be letters only. Any space, digit or comma makes Character.getNumericValue fall outside the range, and the whole word is rejected (the method returns null).
- **S1 — WRONG in part:** "you cannot declare unused extra strands; a split unknot needs Split Union."
  - obtainBraid records every strand position at or below the largest letter that no crossing touches as an unlink component (`loops` → Link.unlinkComp).
  - Example: "bbb" is a 3-strand word whose strand 1 is untouched, so it gives trefoil ⊔ unknot.
  - Only strands above the largest letter cannot be expressed.
  - Source: [Z:knotjob/links/LinkCreator.java; Z:knotjob/links/Link.java, components() = orientation.size() + unlinkComp].
- **Lowercase = +1 = positive, uppercase = −1:** CONFIRMED by the doc's "(4,-5)" example and by the code (crossins.add(±1) and crossingSigns).
- **The 3-braid paper's convention ("a and b leading to positive crossings …"):** CONFIRMED (arXiv:2501.11547 §1). That paper does not mention KnotJob.
- **Orientation for multi-component braid closures:** the report said unverified. My finding is PLAUSIBLE that the braid orientation is kept, from reading the code without running it.
  - Link.firstZero scans crossing positions in the order 0, 3, 2, 1.
  - obtainBraid stores the incoming arcs at positions 0 and 3.
  - So every component starts on an incoming arc, and all strands run in the braid direction.
  - girthMinimize copies the orientation records.
  - The doc's "somewhat randomly assigned" refers to links in general.
- **.kjb / .kts / .tkj / .gld / .dtc / .txt / .adc:** CONFIRMED [Z:Documentation/Open Links.rtf].
  - .kjb tags: CONFIRMED [Z:knotjob/LinkLoader.java]. 0: name, 1: comment, 2: crossings, 3: paths, 4: unreduced Kh, 5: reduced Kh, 9: orientations, and 6–33 for the other invariants.
  - .dtc auto-naming (e.g. 13n0001): CONFIRMED (getDTName).
  - .txt naming: the name is taken from the text before " = ". The doc's example has no space: "Trefoil = PD[[...]]".
- **.brd, one word per line, the word becomes the name, missing from the docs:** CONFIRMED [LinkLoader.addBraids; Knobster file filter]. It is also absent from the README.
- **.adt vs .adc mismatch:** CONFIRMED, and it is worse than reported [KnotJobCommand.acceptableFileName; LinkLoader.enterOtherData].
  - The command line accepts only .txt, .dtc, .adt, .brd, .kts and .kjb.
  - The loader has no "adt" branch, so a .adt file loads nothing, and a .adc file is not accepted on the command line at all.
  - So alphabetical DT, .tkj and .gld are effectively GUI-only.
  - Unknown extensions are silently ignored.
- **stdin mode:** CONFIRMED [KnotJobCommand.calculate/readOneLink].
  - It needs at least one "-" argument; with no arguments the GUI starts.
  - Each link is named "Knot". Output goes to the screen only (printFile is set to false).
  - The README contradicts this ("at least one [file] is needed").
  - The read loop uses a non-short-circuit `&` (`(readLine()!=null) & (length()!=0)`), so EOF without a blank line likely throws a NullPointerException. This is an inference from the code.
- **Bundled tables:** CONFIRMED [Z:Knots/]. Includes FourTor.kjb, KnotsSq1.kjb and StrGrade.kjb. The zipped .txt files must be unzipped before loading.
- **No lookup by name:** CONFIRMED. The menus offer only Sort (by name) and Filters; there is also "Open selected Link(s)".
- **Girth-minimized by default:** CONFIRMED. The "Minimize Girth" checkbox in EnterData defaults to on, and every file and stdin loader passes true. Stored .kjb diagrams are kept as they are.
- **Girth doc quotes and the T(8,9) advice:** CONFIRMED [Z:Documentation/Minimize Girth.rtf].
- **Export format and Szabó compatibility:** CONFIRMED [Z:Documentation/Export Links.rtf].

== (3) OUTPUT ==
- **Labels "Integral/Rational … Khovanov Homology :", "Unreduced Khovanov Homology mod N :", "Odd integral …", "S-Invariant mod p : ", "Even Sq^1 Invariant : (a,b,c,d)", "Graded s-invariant : ", "Signature : ", "Torsion of order m :":** CONFIRMED [KnotJobCommand.printTheStuff; Z:knotjob/links/LinkData.java:427].
- **Sorted by t then q; "t^2 q^5", "t^-1"; coefficients like "2 t^2 q^5":** CONFIRMED [Z:knotjob/polynomial/Coefficient.java, Polynomial.java].
- **Grading start h = −n₋, q = n₊ − 2n₋:** CONFIRMED [Z:knotjob/homology/evenkhov/EvenKhovHomology.java:81-82].
- **doubleTor at (h,q) and (h−1,q) for mod p^k:** CONFIRMED [Z:knotjob/homology/HomologyInfo.java, standardTorPoly].
- **Rebuilt trefoil example:** CONFIRMED as consistent. It matches the code's string format exactly, including the leading space after the colon.
- **KnotInfo's "computed using the program KnotJob" and vector [[0,1,0,1],[0,1,0,3],[0,1,2,5],[0,1,3,9],[2,1,3,7]]:** CONFIRMED (https://knotinfo.org/descriptions/khovanov_odd_mod2_vector.html). The page's own prose example says "[2,1,2,7]", which is a typo on KnotInfo's side.
- **Output file name such as "Rolfsen.kjb_s0_s2_sqo_ku1":** CONFIRMED [KnotJobCommand.newFileName], with corrections:
  - The file is a plain-text copy of the screen lines, not a .kjb. The command line never writes results back into .kjb.
  - Only s, sqe, sqo, sgn, k* and ko* add suffixes.
  - `-csv` turns off the normal file and writes `<name up to its first '.'>.csv`.
  - The CSV header lists only s, Sq^1, signature and the integral Kh/torsion columns, so the columns for rational or mod-p results do not line up.
- **-nf / -ns:** CONFIRMED. `-ns` does not suppress the "… finished." messages.
- **GUI viewer:** CONFIRMED [Z:knotjob/frames/ViewCohomology.java].
  - Radio buttons: integral / rational / modular.
  - "Save as LaTeX file" writes a standalone document using diagbox; the axis can be "q \ h" or "h \ q", and the q-range is selectable.
  - When integral is missing it shows "Approximation" plus "Improvable Primes : …".
  - The mod groups are derived from integral via modHomFrom.
- **sl(3) Table 2: T(4,5) E3-page with h from −12 to 0 and q from 22 to 52:** CONFIRMED (arXiv:2504.07806 Table 2).
- **sl(3) h-sign convention in the program:** UNVERIFIABLE. The sl(3) code is not available.

== (4) UI ==
- **GUI with no arguments, headless otherwise:** CONFIRMED [KnotJob.main].
- **README commands and example:** CONFIRMED.
- **Source-only flags:** CONFIRMED. Also missing from the report: the aliases -sq1o, -sq1e and -sq1s.
- **Invalid prime silently ignored:** CONFIRMED. The same applies to unknown flags and to `-kb` with no number.
- **Flag combinations are not independent:** reduced/unreduced are global switches across all coefficients, so `-kr0 -ku1` computes both reduced and unreduced. If 0 is requested, other coefficients are dropped, and for each prime only the highest power is kept (removeRedundancies).
- **Threads:** CONFIRMED, more precisely as one thread per (invariant × file), all started together.
- **"finished" messages only; time info off:** CONFIRMED [Z:knotjob/dialogs/DialogWrap.java; KnotJobCommand constructor].
- **GUI progress labels:**
  - "Crossing : ", "Girth : ", "Objects : ", "h-Level : " (h-Level only when the girth-info option = 2): CONFIRMED.
  - Abort / Cancel: CONFIRMED [Z:knotjob/AbortInfo.java]. cancel() also sets abort.
  - Calculation Time dialog in Days/Hours/Minutes/Seconds: CONFIRMED [Z:knotjob/dialogs/TimerDialog.java].
- **Caching:** CONFIRMED [EvenKhovCalculator.calculationRequired]. A calculation is skipped if the integral result exists, or if p^k' with k' ≥ k exists. GUI "Save Links" stores results in .kjb, at most 1000 links per file by default.
- **Menus File / Edit / Calculations / Filters / Diagram / Help:** CONFIRMED [Knobster.java]. The Hayden notes describe the same GUI flow: "Char 0 Khovanov Cohomology → Integral → OK".
- **Colour version names:** CONFIRMED (gray → aqua → chart).

== (5) PERFORMANCE ==
- **arXiv:2008.07410:** CONFIRMED (https://arxiv.org/abs/2008.07410). It says T(8,9) reduced odd Kh took "about half an hour on a standard PC with 4GB RAM", with the "slightly slower and less memory efficient" and "can also be calculated" quotes. It is AGT 22 (2022) 1287–1324.
- **arXiv:1811.06432 Table 2:** CONFIRMED (https://arxiv.org/abs/1811.06432). Values 0:57 / 1:17 / 1:31, maximum 4:18, "standard laptop", caption "Calculation times in minutes". Whether entries are m:ss or h:mm is genuinely ambiguous.
- **"20 crossings … 60 crossings" quote:** CONFIRMED, but it is in the paper's introduction about the Khovanov s-invariant in general, not a KnotJob benchmark.
- **sl(3) memory quote ("(6,7)-torus knot … 32GB desktop PC"):** CONFIRMED.
- **sl(3) bulk runs:** CONFIRMED (arXiv:2504.07806 §8.1). The 120,890 seventeen-crossing knots are a filtered subset: hyperbolic knots with the Fox–Milnor Alexander condition, signature 0 and s0 = 0.
- **arXiv:2412.05156:** CONFIRMED (https://arxiv.org/abs/2412.05156). Covers knots up to 16 crossings and the Manolescu–Marengon variant with a full twist on 8 strands (xo = 4). The 10-strand variant is "outside of the range for calculations".
- **arXiv:2507.15305 (Hayden):** CONFIRMED (https://arxiv.org/abs/2507.15305). Quote: Kh "of links with >50 crossings relatively quickly".
- **arXiv metadata:** CONFIRMED. 2504.07806 has v1 10 Apr 2025 and v2 15 Aug 2025; 2501.11547 has v2 3 Oct 2025.

== TAKEAWAYS ==
- **S2:** Takeaway 1 is mostly right. Correct it to say that KnotJob keeps untouched strands at or below the largest letter as split unknots (see S1); only trailing strands cannot be expressed.
- **Takeaway 2:** CONFIRMED. KnotJob ignores bad coefficients silently.
- Takeaways 3–5 are reasonable.
- The one-line answer to the user's "strands + braid word + field" question: KnotJob's command line has no braid argument. You put the word into a .brd file and run `java -jar KnotJob.jar f.brd -ku0`. The strand count is never given.

==== CORRECTED DEVELOPER REFERENCE (Dec 2023 "gray" build; the sl(3) builds after Apr 2025 are unseen) ====

**Inputs**
- **GUI New Link:**
  - Torus (p strands ≤ 10, q ≤ 60, positive only).
  - Pretzel (≤ 12).
  - Braid: letters a–i are σ1–σ9 as positive crossings; A–I are their inverses (negative). Use the standard closure.
  - DT, Gauss, PD and alphabetical-DT codes.
  - Operations: concatenate, split union, mirror, unlink (≤ 10), change orientation or drop components, Whitehead doubles.
- **Braid strand count:** inferred as the index of the largest letter + 1, between 2 and 10. Untouched strands below the largest letter become unknots. Strands above it cannot be expressed. The word must be letters only, with no spaces.
- **Girth minimization:** on by default, via a checkbox.
- **Files:**
  - GUI: .kjb (diagram plus results), .kts, .tkj, .gld, .dtc (Knotscape; names derived automatically), .txt (one PD per line, "Name = PD[...]"), .adc, .brd (one word per line; the word is the name).
  - Command line: .kjb, .kts, .txt, .dtc and .brd only. .adt is accepted but loads nothing; anything else is ignored. With no file, PD codes are read from stdin until a blank line.

**Base ring**
- **Command line:** chosen per invariant by a numeric suffix.
  - -s0 is Q; -s<p> is F_p for one of the 47 primes p ≤ 211.
  - -kb/-kr/-ku/-ko take 0 = Z, 1 = Q, or q = p^k (p ≤ 211, q^2 < 2^31).
  - Invalid values and flags are silently dropped. Z subsumes the other coefficients, and the highest power of each prime wins.
- **GUI:** Integral / Rational / Local (Z with chosen primes inverted), or a prime and power.
- **Grading:** cohomological, with h starting at −n₋ and q at n₊ − 2n₋.

**Outputs**
- **Per link:** the name, then label lines with a Poincaré polynomial in t (h) and q, sorted by t then q.
- **Torsion:** one "Torsion of order m :" line per order. Mod-p^k torsion is printed at both h and h−1.
- **Other invariants:** s as an integer, Sq^1 as (r+,s+,r−,s−), graded s as "s.-s".
- **Files:** a text file named <input>_<suffixes> (-nf turns it off), or -csv for <base>.csv. -ns turns off screen output.
- **GUI:** a bigraded table viewer (integral / rational / modular), a multi-ring "Approximation" view with "Improvable Primes", and a standalone LaTeX table export. Results are saved inside .kjb and not recomputed.

**UI**
- One jar: Swing GUI with no arguments, headless command line otherwise.
- Command line: one thread per (invariant × file); prints only "<invariant> of file <f> finished."; no per-crossing progress or timing.
- GUI progress: Crossing u/n, Girth, Objects, h-Level, with Abort (current link) and Cancel (whole batch), and an optional total-time dialog.
- Menus: File, Edit, Calculations, Filters, Diagram, Help (RTF docs).
- No colored or gl_N (N ≥ 4) homology. sl(3) exists only in builds from 2025 on, whose command-line flags and output format I could not check.
