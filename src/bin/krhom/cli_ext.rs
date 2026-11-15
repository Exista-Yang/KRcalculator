//! F-CLI extensions of krhom (kept in their own module to ease merging with other branches):
//! clean limit refusals before any engine work, `--mark K`, `--knot` / `--link NAME`,
//! `--latex OUT` and `--estimate`. See docs/CALCULATOR.md.
use super::*;

// ------------------------------------------------------------------ limits (refusals, exit 1)

/// Largest N of the tower closure: the truncation cap N - 1 must fit the packed monomials
/// (src/uncolored/tower.rs `Caps::new`: cap <= 31).
pub const TOWER_MAX_N: usize = 32;
/// The exact oracle's dimension limit (src/colored/cclose.rs `GrAlg::try_new_capped`).
pub const EXACT_MAX_DIM: u128 = 50_000_000;

/// The tower closure's limits (Soergel and cable routes, `--closure tower`): at most
/// `tower::MAX_VARS` strands and N <= 32. Checked before any engine work; exit 1.
pub fn tower_limits(o: &Opts, strands: usize, run_ns: &[usize]) -> Result<(), CalcErr> {
    if !matches!(o.closure, Closure::Tower(_)) {
        return Ok(());
    }
    let max_vars = krhom::uncolored::tower::MAX_VARS;
    if strands > max_vars {
        return Err(CalcErr::Unsupported(format!(
            "{} strands: the tower closure of the Soergel engine supports at most {} strands (one packed variable per strand, src/uncolored/tower.rs MAX_VARS); about 5 strands are practical",
            strands, max_vars
        )));
    }
    if let Some(&n) = run_ns.iter().find(|&&n| n > TOWER_MAX_N) {
        return Err(CalcErr::Unsupported(format!(
            "N = {}: the tower closure supports N <= {} (its truncation cap N - 1 must be <= 31 in the packed monomials); use --closure exact for larger N",
            n, TOWER_MAX_N
        )));
    }
    Ok(())
}

/// C(n, k) (saturating), the dimension of H*(Gr(k, n)).
pub fn binom(n: usize, k: usize) -> u128 {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut r: u128 = 1;
    for i in 0..k {
        r = r.saturating_mul((n - i) as u128) / (i as u128 + 1);
    }
    r
}

/// dim A_lambda of the exact 2-strand closure as cclose.rs computes it: the product of
/// C(N, a_i) over the unmarked blocks (the block at position 0 is marked when reduced).
pub fn exact2_alg_dim(labels: &[usize], n: usize, reduced: bool) -> u128 {
    labels.iter().enumerate().filter(|&(i, _)| !(reduced && i == 0)).map(|(_, &a)| binom(n, a)).fold(1u128, |x, y| x.saturating_mul(y))
}

/// The exact oracle's dimension limit dim A_lambda <= 5e7 (e.g. color 3 at N >= 36 and color 2
/// at N >= 120, unreduced). Checked before any engine work; exit 1.
pub fn exact2_dim(labels: &[usize], run_ns: &[usize], reduced: bool) -> Result<(), String> {
    for &n in run_ns {
        let d = exact2_alg_dim(labels, n, reduced);
        if d > EXACT_MAX_DIM {
            return Err(format!(
                "N = {} with labels {:?}{}: the exact 2-strand engine's algebra A_lambda has dimension {} (product of C(N, a) over the unmarked strands), above its limit 5e7 (the exact oracle stores dense tables of that size); e.g. color 3 needs N <= 35 and color 2 N <= 119 unreduced",
                n,
                labels,
                if reduced { ", reduced" } else { "" },
                d
            ));
        }
    }
    Ok(())
}

// ------------------------------------------------------------------ --mark K

pub fn parse_mark(v: Option<String>, reduced: bool) -> Result<Option<usize>, CalcErr> {
    let Some(v) = v else { return Ok(None) };
    let k = v.trim().parse::<usize>().map_err(|_| CalcErr::Unsupported(format!("--mark {:?}: expected a component number 0, 1, ...", v)))?;
    if !reduced {
        return Err(CalcErr::Unsupported("--mark applies to reduced homology: add --reduced".into()));
    }
    Ok(Some(k))
}

/// The method-line note of reduced homology.
pub fn red_note(o: &Opts) -> String {
    match (o.reduced, o.mark) {
        (false, _) => String::new(),
        (true, None) => ", reduced (strand at bottom position 0 marked)".into(),
        (true, Some(k)) => format!(", reduced (component {} marked: --mark {}, at bottom position 0)", k, k),
    }
}

/// Free reduction of a braid word (cancels adjacent inverse letters; the bottom component
/// labels do not change). No cyclic reduction: that would undo the conjugation below.
pub fn free_reduce(w: &[i64]) -> Vec<i64> {
    let mut st: Vec<i64> = vec![];
    for &x in w {
        if st.last() == Some(&-x) {
            st.pop();
        } else {
            st.push(x);
        }
    }
    st
}

/// `bw` conjugated so that link component `k` passes through bottom position 0. First a cyclic
/// rotation of the word (x g -> g x, same length) if one brings component k to position 0;
/// otherwise (s_1 ... s_p) bw (s_p^-1 ... s_1^-1), p the smallest bottom position of
/// component k, freely reduced. A conjugate braid has the same closure; the components keep
/// their ids (colors are per component id). Returns (braid, p).
pub fn conjugate_to_mark(bw: &BraidWord, k: usize) -> (BraidWord, usize) {
    let p = bw.components.iter().position(|&c| c == k).expect("component exists");
    if p == 0 {
        return (bw.clone(), 0);
    }
    let sw = bw.signed();
    // rotations: the bottom labels of w[r..] w[..r] are the labels after the first r letters
    let mut comps = bw.components.clone();
    let mut found = None;
    for r in 1..sw.len() {
        let i = sw[r - 1].unsigned_abs() as usize - 1;
        comps.swap(i, i + 1);
        if comps[0] == k {
            found = Some((r, comps.clone()));
            break;
        }
    }
    let (w, comps) = match found {
        Some((r, c)) => ([&sw[r..], &sw[..r]].concat(), c),
        None => {
            let pre: Vec<i64> = (1..=p as i64).collect();
            let mut w = pre.clone();
            w.extend(&sw);
            w.extend(pre.iter().rev().map(|x| -x));
            // the labels at the top of the prefix are bw's; pull them back through it
            let mut c = bw.components.clone();
            for &g in pre.iter().rev() {
                let i = g as usize - 1;
                c.swap(i, i + 1);
            }
            (free_reduce(&w), c)
        }
    };
    let mut out = BraidWord::from_signed(bw.strands, &w);
    out.components = comps;
    assert_eq!(out.components[0], k, "conjugation moved component {} to position 0", k);
    assert_eq!(out.top_components(), out.components, "conjugated braid is not a closure");
    (out, p)
}

/// --mark K: conjugate the braid (in place) so component K is at bottom position 0.
/// Returns the original position when the braid changed.
pub fn apply_mark(bw: &mut BraidWord, o: &Opts) -> Result<Option<usize>, CalcErr> {
    let Some(k) = o.mark else { return Ok(None) };
    let nc = bw.num_components();
    if nc == 1 && k != 0 {
        return Err(CalcErr::Unsupported(format!("--mark {}: a knot has one component; only --mark 0 is accepted", k)));
    }
    if k >= nc {
        return Err(CalcErr::Unsupported(format!("--mark {}: the link has {} components (numbered 0..{} by their smallest bottom strand position)", k, nc, nc - 1)));
    }
    let (b, p) = conjugate_to_mark(bw, k);
    if p == 0 {
        return Ok(None);
    }
    *bw = b;
    Ok(Some(p))
}

// ------------------------------------------------------------------ --knot / --link NAME

#[derive(Clone, Debug)]
pub struct TableBraid {
    /// "knot" | "link"
    pub kind: &'static str,
    pub name: String,
    pub file: String,
    pub strands: usize,
    pub word: Vec<i64>,
}

/// Directory of the knot and link tables: `$KRHOM_DATA_DIR`, else the source tree's data/.
pub fn data_dir() -> PathBuf {
    match std::env::var_os("KRHOM_DATA_DIR").filter(|d| !d.is_empty()) {
        Some(d) => PathBuf::from(d),
        None => PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("data"),
    }
}

pub fn table_file(kind: &str) -> &'static str {
    if kind == "knot" { "knots_upto12.jsonl" } else { "links_upto11.jsonl" }
}

/// (name, strands, word) of one JSON line of the tables (fields "name" and "braid":
/// {"strands": n, "word": [...]}); None if the line has another shape.
pub fn parse_table_line(l: &str) -> Option<(String, usize, Vec<i64>)> {
    let s = l.find("\"name\":\"")? + 8;
    let name = l[s..s + l[s..].find('"')?].to_string();
    let b = l.find("\"braid\":{\"strands\":")? + 19;
    let rest = &l[b..];
    let strands: usize = rest[..rest.find(',')?].trim().parse().ok()?;
    let ws = rest.find("\"word\":[")? + 8;
    let we = ws + rest[ws..].find(']')?;
    let word: Vec<i64> = rest[ws..we].split(',').filter(|x| !x.trim().is_empty()).map(|x| x.trim().parse().ok()).collect::<Option<_>>()?;
    Some((name, strands, word))
}

fn table_names(kind: &str) -> Result<Vec<(String, usize, Vec<i64>)>, CalcErr> {
    let p = data_dir().join(table_file(kind));
    let text = std::fs::read_to_string(&p).map_err(|e| CalcErr::Unsupported(format!("--{} needs the table {} ({}); set KRHOM_DATA_DIR to the repository's data/ directory", kind, p.display(), e)))?;
    Ok(text.lines().filter_map(parse_table_line).collect())
}

/// Levenshtein distance (for close matches).
pub fn edit_distance(a: &str, b: &str) -> usize {
    let (a, b): (Vec<char>, Vec<char>) = (a.chars().collect(), b.chars().collect());
    let mut prev: Vec<usize> = (0..=b.len()).collect();
    for i in 1..=a.len() {
        let mut cur = vec![i; b.len() + 1];
        for j in 1..=b.len() {
            cur[j] = (prev[j] + 1).min(cur[j - 1] + 1).min(prev[j - 1] + usize::from(a[i - 1] != b[j - 1]));
        }
        prev = cur;
    }
    prev[b.len()]
}

/// Up to 5 names closest to `want` (case-insensitive edit distance, then table order).
pub fn close_matches(want: &str, names: &[String]) -> Vec<String> {
    let w = want.to_lowercase();
    let mut v: Vec<(usize, usize, &String)> = names.iter().enumerate().map(|(i, n)| (edit_distance(&w, &n.to_lowercase()), i, n)).collect();
    v.sort();
    v.into_iter().take(5).map(|(_, _, n)| n.clone()).collect()
}

/// The braid of a knot (KnotInfo) or an oriented link (LinkInfo, names with the orientation
/// suffix, e.g. L2a1{1}). Exact name first, then a unique case-insensitive match; an unoriented
/// link name lists its orientations; otherwise exit 1 with up to 5 close matches.
pub fn lookup_table(kind: &'static str, name: &str) -> Result<TableBraid, CalcErr> {
    let recs = table_names(kind)?;
    let want = name.trim();
    let mk = |r: &(String, usize, Vec<i64>)| TableBraid { kind, name: r.0.clone(), file: format!("data/{}", table_file(kind)), strands: r.1, word: r.2.clone() };
    if let Some(r) = recs.iter().find(|r| r.0 == want) {
        return Ok(mk(r));
    }
    let ci: Vec<&(String, usize, Vec<i64>)> = recs.iter().filter(|r| r.0.eq_ignore_ascii_case(want)).collect();
    if ci.len() == 1 {
        return Ok(mk(ci[0]));
    }
    if kind == "link" {
        let pre = format!("{}{{", want.to_lowercase());
        let variants: Vec<&String> = recs.iter().map(|r| &r.0).filter(|n| n.to_lowercase().starts_with(&pre)).collect();
        if !variants.is_empty() {
            return Err(CalcErr::Unsupported(format!(
                "--link {}: give an oriented name; LinkInfo lists {} orientations: {}",
                want,
                variants.len(),
                variants.iter().map(|s| s.as_str()).collect::<Vec<_>>().join(", ")
            )));
        }
    }
    let other = if kind == "knot" { "link" } else { "knot" };
    if let Ok(o) = table_names(other) {
        if o.iter().any(|r| r.0.eq_ignore_ascii_case(want)) {
            return Err(CalcErr::Unsupported(format!("--{} {}: that is a {} name; use --{} {}", kind, want, other, other, want)));
        }
    }
    let names: Vec<String> = recs.iter().map(|r| r.0.clone()).collect();
    Err(CalcErr::Unsupported(format!("--{} {}: unknown {} name (data/{}); close matches: {}", kind, want, kind, table_file(kind), close_matches(want, &names).join(", "))))
}

/// The braid of a table entry (prints the source line).
pub fn table_braid(tb: &TableBraid, strands: Option<usize>) -> Result<BraidWord, CalcErr> {
    if let Some(s) = strands {
        if s != tb.strands {
            return Err(CalcErr::Unsupported(format!("--strands {} differs from the table's {} strands for {}: drop --strands", s, tb.strands, tb.name)));
        }
    }
    println!(
        "source: {} {} from {} ({})",
        tb.kind,
        tb.name,
        tb.file,
        if tb.kind == "knot" { "KnotInfo braid; KnotInfo chirality" } else { "LinkInfo braid; LinkInfo orientation as named" }
    );
    Ok(BraidWord::from_signed(tb.strands, &tb.word))
}

// ------------------------------------------------------------------ --latex OUT

pub fn tex_escape(s: &str) -> String {
    let mut o = String::new();
    for c in s.chars() {
        match c {
            '_' | '{' | '}' | '#' | '$' | '%' | '&' => {
                o.push('\\');
                o.push(c);
            }
            '\\' => o.push_str("\\textbackslash{}"),
            '^' => o.push_str("\\^{}"),
            '~' => o.push_str("\\~{}"),
            c => o.push(c),
        }
    }
    o
}

/// One tabular of a Poincare polynomial: rows t ascending, columns q ascending (step 2 when
/// every q has the same parity), blank for 0 (as `--table`).
pub fn latex_tabular(p: &Poincare) -> String {
    if p.c.is_empty() {
        return "\\emph{(the homology is 0)}".into();
    }
    let (tmin, tmax) = (p.c.keys().map(|k| k.0).min().unwrap(), p.c.keys().map(|k| k.0).max().unwrap());
    let (qmin, qmax) = (p.c.keys().map(|k| k.1).min().unwrap(), p.c.keys().map(|k| k.1).max().unwrap());
    let step = if p.c.keys().all(|k| (k.1 - qmin).rem_euclid(2) == 0) { 2 } else { 1 };
    let qs: Vec<i32> = (qmin..=qmax).step_by(step).collect();
    let mut s = format!("\\begin{{tabular}}{{r|{}}}\n", "r".repeat(qs.len()));
    s.push_str(&format!("$t \\backslash q$ & {} \\\\\n\\hline\n", qs.iter().map(|q| format!("${}$", q)).collect::<Vec<_>>().join(" & ")));
    for t in tmin..=tmax {
        let row: Vec<String> = qs.iter().map(|q| p.c.get(&(t, *q)).map_or(String::new(), |k| k.to_string())).collect();
        s.push_str(&format!("${}$ & {} \\\\\n", t, row.join(" & ")));
    }
    s.push_str("\\end{tabular}");
    s
}

/// One section (one input) of the --latex file.
pub fn latex_section(item: &Item, o: &Opts, r: &Result<Outcome, CalcErr>) -> String {
    let what = match &item.input {
        Input::Table(tb) => format!("{} {}", tb.kind, tb.name),
        Input::Braid(b) => item.name.clone().map_or(format!("braid {}", b), |n| format!("{} (braid {})", n, b)),
        Input::Pd(p) => item.name.clone().map_or(format!("PD {}", p), |n| format!("{} (PD {})", n, p)),
    };
    let mut s = format!("\\section*{{{}}}\n", tex_escape(&what));
    match r {
        Err(e) => s.push_str(&format!("Not computed: exit {} ({}): \\verb|{}|\n", e.code(), e.kind(), e.msg().lines().next().unwrap_or("").replace('|', "/"))),
        Ok(out) if out.estimate.is_some() || out.check_library.is_some() => s.push_str("No homology computed (--estimate / --check-library).\n"),
        Ok(out) => {
            let used = out.used.as_ref().or(out.orig.as_ref());
            let braid = used.map_or(String::new(), |b| format!("braid ${}$ on {} strands", b.signed().iter().map(|x| x.to_string()).collect::<Vec<_>>().join("\\,"), b.strands));
            let mode = match (o.reduced, o.mark) {
                (false, _) => "unreduced".to_string(),
                (true, None) => "reduced (strand at bottom position 0 marked)".to_string(),
                (true, Some(k)) => format!("reduced (component {} marked)", k),
            };
            for (n, p) in &out.results {
                s.push_str("\\begin{table}[ht]\n\\centering\n\\resizebox{\\ifdim\\width>\\linewidth\\linewidth\\else\\width\\fi}{!}{%\n");
                s.push_str(&latex_tabular(p));
                s.push_str("}\n");
                s.push_str(&format!(
                    "\\caption{{{}: {}, colors {}, $N = {}$, {}; coefficients {}; rows $t$, columns $q$ (krhom's grading), total rank {}.}}\n\\end{{table}}\n",
                    tex_escape(&what),
                    braid,
                    tex_escape(&format!("{:?}", out.colors)),
                    n,
                    mode,
                    tex_escape(&out.coefficients),
                    p.total_rank()
                ));
            }
            s.push_str("\\clearpage\n");
        }
    }
    s
}

/// The whole --latex file: a standalone article (graphicx only); the sections between the
/// marker comments can be pasted into another document that loads graphicx.
pub fn latex_doc(sections: &[String]) -> String {
    format!(
        "% written by krhom --latex; compile with pdflatex (needs only graphicx)\n\\documentclass{{article}}\n\\usepackage{{graphicx}}\n\\begin{{document}}\n% ---- krhom tables begin\n{}% ---- krhom tables end\n\\end{{document}}\n",
        sections.concat()
    )
}

// ------------------------------------------------------------------ --estimate (box route)

/// --estimate on the 3-4 strand box route: per prime, load the library and time the record
/// scan, the determinant pass and the replay of ONE box point; extrapolate to the full run at
/// the current thread count. Nothing is computed.
pub fn box_estimate(o: &Opts, out: &mut Outcome, a: (&RegEntry, &InstalledLib), b: (&RegEntry, Option<&InstalledLib>), wi: &[i32], run_ns: &[usize], mem_mb: u64) -> Result<BTreeMap<usize, Poincare>, CalcErr> {
    let mk = |seed: u64| BoxParams { word: wi.to_vec(), ns: run_ns.to_vec(), reduced: o.reduced, seed, retries: o.retries, mem_mb, check_only: false, replay_only: true };
    let threads = krhom::num_threads();
    let mut ests = vec![guard("estimate (F_2147483647)", || estimate_prime::<F31>(a.0, a.1, &mk(o.seed)))?];
    if let Some(inst) = b.1 {
        ests.push(guard("estimate (F_4294967291)", || estimate_prime::<F32b>(b.0, inst, &mk(o.seed2)))?);
    }
    let mut lines = vec![];
    let (mut tot, mut mem) = (0.0f64, 0u64);
    for e in &ests {
        let (wall, mb) = e.extrapolate(threads);
        tot += wall;
        mem = mem.max(mb);
        eprintln!("[{}] estimate: library {:.1}s, record scan {:.1}s, one box point {:.1}s (determinant pass {:.2}s, peak {} MB), VmHWM {} MB", e.field, e.lib_build_secs, e.record_secs, e.point_secs, e.det_secs, e.point_peak_bytes >> 20, e.vmhwm_kb / 1024);
        lines.push(format!("{}: {} box points ({}), library {:.0}s + record scan {:.0}s + {} x {:.1}s per point (replay + determinant pass) / {} threads + scaling check -> about {} wall, peak memory about {} MB", e.field, e.points, if o.reduced { "reduced grid" } else { "unreduced grid" }, e.lib_build_secs, e.record_secs, e.points, e.point_secs + e.det_secs, threads, hms(wall), mb));
    }
    lines.push(format!(
        "total: about {} wall for {} prime{} at {} threads, peak memory about {} MB (an estimate from one timed point per prime; the closure phase is not included and adds minutes on 4 strands; nothing was computed)",
        hms(tot),
        ests.len(),
        if ests.len() == 1 { "" } else { "s" },
        threads,
        mem
    ));
    out.estimate = Some(lines.join("\n"));
    Ok(BTreeMap::new())
}

pub fn hms(s: f64) -> String {
    let s = s.max(0.0).round() as u64;
    if s >= 3600 { format!("{}h{:02}m", s / 3600, (s % 3600) / 60) } else if s >= 60 { format!("{}m{:02}s", s / 60, s % 60) } else { format!("{}s", s) }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dims_as_cclose() {
        assert_eq!(binom(36, 3), 7140);
        assert_eq!(binom(120, 2), 7140);
        assert_eq!(exact2_alg_dim(&[3, 3], 35, false), 6545 * 6545);
        assert!(exact2_dim(&[3, 3], &[35], false).is_ok());
        assert!(exact2_dim(&[3, 3], &[36], false).is_err());
        assert!(exact2_dim(&[3, 3], &[36], true).is_ok());
        assert!(exact2_dim(&[2, 2], &[119], false).is_ok());
        assert!(exact2_dim(&[2, 2], &[120], false).is_err());
        assert!(exact2_dim(&[2, 1], &[200], false).is_ok());
    }

    #[test]
    fn conjugation_keeps_the_closure() {
        // every component of a few links is brought to bottom position 0
        for (n, w) in [(2, vec![1i64, 1]), (3, vec![1, 1, 2, 2]), (3, vec![1, -2, 1, -2, 1, -2]), (4, vec![1, 1, 3, 3, 2, 2]), (4, vec![2, 2, 1, 3, -2, 1, 3])] {
            let bw = BraidWord::from_signed(n, &w);
            for k in 0..bw.num_components() {
                let (c, _) = conjugate_to_mark(&bw, k);
                assert_eq!(c.components[0], k);
                assert_eq!(c.num_components(), bw.num_components());
                // per-component strand counts are kept
                for j in 0..bw.num_components() {
                    assert_eq!(c.components.iter().filter(|&&x| x == j).count(), bw.components.iter().filter(|&&x| x == j).count());
                }
                assert!(c.word.len() <= bw.word.len() + 2 * (n - 1));
            }
        }
        assert_eq!(free_reduce(&[1, 2, -2, -1, 3]), vec![3]);
    }

    #[test]
    fn table_line_and_matches() {
        let l = r#"{"name":"L2a1{1}","name_unoriented":"L2a1","braid":{"strands":2,"word":[1,1]},"x":[]}"#;
        assert_eq!(parse_table_line(l), Some(("L2a1{1}".into(), 2, vec![1, 1])));
        assert_eq!(edit_distance("4_2", "4_1"), 1);
        let names: Vec<String> = ["3_1", "4_1", "5_1", "5_2", "6_1", "10_1"].iter().map(|s| s.to_string()).collect();
        assert_eq!(close_matches("4_2", &names), ["4_1", "5_2", "3_1", "5_1", "6_1"]);
        assert_eq!(tex_escape("L2a1{0}_x"), "L2a1\\{0\\}\\_x");
    }
}
