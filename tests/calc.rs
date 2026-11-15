//! CALC-1: the krhom calculator CLI (spec F1, F3; CALC_DECISIONS.md D1-D4).
//!
//! Default tests run `krhom` with `KRHOM_LIB_DIR` set to a fresh empty directory and
//! `KRHOM_THREADS=1`: refusals and exit codes, braid input forms, color broadcast, N below a
//! color, batch / CSV / JSON plumbing, and the 2-strand cross-engine parity (`--box-live`).
//! The `#[ignore]` tests need both libraries installed (`tools/fetch_library.sh l222_n4 l2222_n4`)
//! under `KRHOM_LIB_DIR` (or the default root) and take minutes (H3-H6).
use std::path::PathBuf;
use std::process::Command;

struct Run {
    code: i32,
    out: String,
    err: String,
}

fn empty_lib() -> PathBuf {
    use std::sync::atomic::{AtomicUsize, Ordering};
    static K: AtomicUsize = AtomicUsize::new(0);
    let d = std::env::temp_dir().join(format!("krhom_calc_emptylib_{}_{}", std::process::id(), K.fetch_add(1, Ordering::SeqCst)));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}
fn krhom_env(args: &[&str], env: &[(&str, &str)]) -> Run {
    let mut c = Command::new(env!("CARGO_BIN_EXE_krhom"));
    c.args(args).env("KRHOM_THREADS", "1").env("KRHOM_LIB_DIR", empty_lib());
    for v in ["KRHOM_NO_M5I", "KRHOM_FORCE_LADDER", "KRHOM_MEM_CAP_MB"] {
        c.env_remove(v);
    }
    for (k, v) in env {
        c.env(k, v);
    }
    let o = c.output().expect("cannot run krhom");
    Run { code: o.status.code().unwrap_or(-1), out: String::from_utf8_lossy(&o.stdout).into(), err: String::from_utf8_lossy(&o.stderr).into() }
}
fn krhom(args: &[&str]) -> Run {
    krhom_env(args, &[])
}
/// The N blocks of the output (what tests compare).
fn blocks(s: &str) -> String {
    s.lines().filter(|l| l.starts_with("N =") || l.trim_start().starts_with("Poincare") || l.trim_start().starts_with("total rank") || l.trim_start().starts_with("Euler")).collect::<Vec<_>>().join("\n")
}
fn ok(r: &Run, what: &str) {
    assert_eq!(r.code, 0, "{}: exit {}\nstdout:\n{}\nstderr:\n{}", what, r.code, r.out, r.err);
}
fn refused(args: &[&str], code: i32, needle: &str) {
    let r = krhom(args);
    assert_eq!(r.code, code, "krhom {:?}: exit {} (want {})\nstdout:\n{}\nstderr:\n{}", args, r.code, code, r.out, r.err);
    assert!(r.err.contains(needle), "krhom {:?}: stderr lacks {:?}:\n{}", args, needle, r.err);
}

#[test]
fn box_route_without_library_exits_4_with_the_fetch_command() {
    refused(&["2,3,4", "--braid", "1 -2 1 -2", "--strands", "3", "--colors", "2", "--colored-method", "singular"], 4, "tools/fetch_library.sh l222_n4");
    // the default route for colored input is the same (D2), and letters work
    refused(&["2,3,4", "--braid", "aBaB", "--colors", "2"], 4, "tools/fetch_library.sh l222_n4");
    refused(&["2,3,4", "--braid", "1 1 2 -1 -3 2 -3", "--colors", "2", "--check-library"], 4, "tools/fetch_library.sh l2222_n4");
}

#[test]
fn refusals_exit_1() {
    // mixed labels on a 3-component chain (simplify cannot drop strands: strands <= components)
    refused(&["2,3", "--braid", "1 1 2 2", "--strands", "3", "--colors", "2,2,1"], 1, "not uniform");
    // (3,3,3), 5 strands, N above the library
    refused(&["3,4", "--braid", "1 -2 1 -2", "--strands", "3", "--colors", "3"], 1, "no library for labels");
    refused(&["2,3", "--braid", "1 1 2 2 3 3 4 4", "--strands", "5", "--colors", "2"], 1, "no library for labels");
    refused(&["5", "--braid", "1 -2 1 -2", "--strands", "3", "--colors", "2"], 1, "above the library");
    // uncolored input with the singular override on 3 strands
    refused(&["2,3", "--braid", "1 -2 1 -2", "--strands", "3", "--colored-method", "singular"], 1, "uncolored input");
    // wrong number of colors, split braid
    refused(&["2", "--braid", "1 1", "--strands", "2", "--colors", "2,2,2"], 1, "3 colors given");
    refused(&["2,3", "--braid", "1 1", "--strands", "3", "--colors", "2", "--no-simplify"], 1, "split link");
    // parse errors (review K9: exit 1, never a panic)
    refused(&["2,x", "--braid", "1 1 1"], 1, "N list");
    refused(&["0", "--braid", "1 1 1"], 1, "at least 1");
    refused(&["2", "--braid", "1 0 1"], 1, "generator 0");
    refused(&["2", "--braid", "a1"], 1, "mixes letters and integers");
    refused(&["2", "--braid", "a b"], 1, "no spaces");
    refused(&["2", "--braid", "1 2", "--strands", "2"], 1, "too small");
    refused(&["2", "--braid", "1 1.5"], 1, "not a signed integer");
    refused(&["2", "--braid", "1 1", "--colors", "0"], 1, "colors are integers >= 1");
    refused(&["2", "--braid", "1 1", "--colors", "2,x"], 1, "colors are integers >= 1");
    refused(&["2", "--pd", "[[1,2,3]]"], 1, "invalid PD code");
    refused(&["2", "--pd", "[[1,4,2,5],[3,6,4,1],[5,2,6,9]]"], 1, "invalid PD code");
    refused(&["2", "--braid", "1 1", "--frobnicate"], 1, "unknown argument");
    refused(&["2", "--braid", "1 1", "--colored-method", "magic"], 1, "cable, cube or singular");
    refused(&["2", "--braid", "1 1", "--pd", "[[1,4,2,5],[3,6,4,1],[5,2,6,3]]"], 1, "exactly one of");
    refused(&["2", "--braid", "1 1 1", "--colors", "2", "--colored-method", "cube", "--reduced"], 1, "--reduced is only supported");
    refused(&["4", "--braid", "1 1 1", "--colors", "4"], 1, "colors 1..=3");
    // validation-only environment switches are refused on the box route (review K12)
    let r = krhom_env(&["2,3", "--braid", "1 -2 1 -2", "--strands", "3", "--colors", "2"], &[("KRHOM_NO_M5I", "1")]);
    assert_eq!(r.code, 1, "{}", r.err);
    assert!(r.err.contains("KRHOM_NO_M5I"), "{}", r.err);
}

#[test]
fn simplification_and_broadcast() {
    // "1 1 1 2" on 3 strands destabilizes to the 2-strand trefoil
    let a = krhom(&["2,3", "--braid", "1 1 1 2", "--strands", "3", "--colors", "2"]);
    ok(&a, "1 1 1 2");
    assert!(a.out.contains("simplified braid: 1 1 1  (strands 2)"), "{}", a.out);
    let b = krhom(&["2,3", "--braid", "1 1 1", "--strands", "2", "--colors", "2"]);
    ok(&b, "1 1 1");
    assert_eq!(blocks(&a.out), blocks(&b.out));
    assert!(b.out.contains("coefficients: Q") && !b.out.contains("fields:"), "{}", b.out);
    // KnotJob letters with the default strand count
    let c = krhom(&["2,3", "--braid", "aaa", "--colors", "2"]);
    ok(&c, "aaa");
    assert_eq!(blocks(&c.out), blocks(&b.out));
    // a single color applies to every component
    let d = krhom(&["2,3", "--braid", "1 1", "--strands", "2", "--colors", "2"]);
    let e = krhom(&["2,3", "--braid", "1 1", "--strands", "2", "--colors", "2,2"]);
    ok(&d, "Hopf 2");
    ok(&e, "Hopf 2,2");
    assert_eq!(blocks(&d.out), blocks(&e.out));
    // the 1-strand unknot after destabilization: [3 2]_q has rank 3
    let f = krhom(&["3", "--braid", "1 2", "--strands", "3", "--colors", "2", "--colored-method", "singular"]);
    ok(&f, "1 2");
    assert!(f.out.contains("simplified braid:") && f.out.contains("total rank: 3"), "{}", f.out);
}

#[test]
fn uncolored_route_two_primes_and_letters() {
    let a = krhom(&["2,3", "--braid", "aBaB"]);
    ok(&a, "aBaB");
    assert!(a.out.contains("strands 3") && a.out.contains("coefficients: Q"), "{}", a.out);
    assert!(a.err.contains("F_2147483647 and F_4294967291 agree"), "{}", a.err);
    let b = krhom(&["2,3", "--braid", "1 -2 1 -2", "--strands", "3"]);
    ok(&b, "1 -2 1 -2");
    assert_eq!(blocks(&a.out), blocks(&b.out));
    // the Khovanov (N=2) figure-eight has total rank 6 unreduced
    assert!(a.out.contains("total rank: 6"), "{}", a.out);
    // PD input
    let c = krhom(&["2", "--pd", "[[1,4,2,5],[3,6,4,1],[5,2,6,3]]"]);
    ok(&c, "PD trefoil");
    assert!(c.out.contains("coefficients: Q"), "{}", c.out);
}

#[test]
fn n_below_the_color_is_zero_on_every_route() {
    let a = krhom(&["1,2,3", "--braid", "1 1 1", "--strands", "2", "--colors", "2"]);
    ok(&a, "N=1 color 2");
    assert!(a.out.contains("N = 1\n  Poincare polynomial: 0\n  total rank: 0"), "{}", a.out);
    assert!(a.err.contains("below the color"), "{}", a.err);
    // box route: no library needed when every N is below the color
    let b = krhom(&["1", "--braid", "1 -2 1 -2", "--strands", "3", "--colors", "2"]);
    ok(&b, "box route, N below the color");
    assert!(b.out.contains("N = 1\n  Poincare polynomial: 0"), "{}", b.out);
}

#[test]
fn box_live_matches_the_exact_path() {
    for red in [false, true] {
        let mut ex = vec!["2,3,4", "--braid", "1 1 1", "--strands", "2", "--colors", "2", "--colored-method", "singular"];
        if red {
            ex.push("--reduced");
        }
        let a = krhom(&ex);
        ok(&a, "exact");
        let mut lv = ex.clone();
        lv.push("--box-live");
        let b = krhom(&lv);
        ok(&b, "box-live");
        assert_eq!(blocks(&a.out), blocks(&b.out), "reduced {}", red);
        assert!(b.out.contains("coefficients: Q") && b.err.contains("agree"), "{}\n{}", b.out, b.err);
    }
    refused(&["2,3", "--braid", "1 -2 1 -2", "--strands", "3", "--colors", "2", "--box-live", "--no-simplify"], 1, "--box-live");
}

#[test]
fn json_csv_and_batch() {
    let d = empty_lib();
    let j = d.join("o.json");
    let a = krhom(&["2,3", "--braid", "1 1 1", "--strands", "2", "--colors", "2", "--json", j.to_str().unwrap()]);
    ok(&a, "json");
    let s = std::fs::read_to_string(&j).unwrap();
    for k in ["\"tool\":\"krhom\"", "\"route\":\"exact2\"", "\"coefficients\":\"Q\"", "\"agree\":true", "\"results\":[", "\"used_braid\":[1, 1, 1]", "\"error\":null"] {
        assert!(s.contains(k), "json lacks {}: {}", k, s);
    }
    assert!(s.contains("\"N\":3,\"P\":{\"str\":\"t^0q^-6 + t^0q^-4 + t^0q^-2 + t^2q^-8 + t^2q^-6 + t^3q^-14 + t^3q^-12\",\"rank\":7"), "{}", s);
    assert!(a.out.contains("total rank: 7"));
    // batch: names, a bad line (reported, the batch continues, nonzero exit), JSON array, CSV
    let f = d.join("b.txt");
    std::fs::write(&f, "# a comment\ntrefoil; aaa\n\nhopf;1 1\nbad; 1 0\nfig8; aBaB  # 4_1\n").unwrap();
    let (jb, cb) = (d.join("b.json"), d.join("b.csv"));
    let b = krhom(&["2", "--file", f.to_str().unwrap(), "--json", jb.to_str().unwrap(), "--csv", cb.to_str().unwrap()]);
    assert_eq!(b.code, 1, "{}\n{}", b.out, b.err);
    for n in ["name: trefoil", "name: hopf", "name: fig8"] {
        assert!(b.out.contains(n), "{}", b.out);
    }
    assert!(b.err.contains("bad: braid") && b.err.contains("generator 0"), "{}", b.err);
    let sj = std::fs::read_to_string(&jb).unwrap();
    assert!(sj.starts_with("[\n{") && sj.matches("\"tool\":\"krhom\"").count() == 4 && sj.contains("\"kind\":\"unsupported\""), "{}", sj);
    let sc = std::fs::read_to_string(&cb).unwrap();
    assert!(sc.starts_with("name,N,reduced,t,q,rank\n") && sc.contains("trefoil,2,false,"), "{}", sc);
    // the figure-eight Khovanov homology: rank 6 = sum of the CSV ranks
    let r8: i64 = sc.lines().filter(|l| l.starts_with("fig8,")).map(|l| l.rsplit(',').next().unwrap().parse::<i64>().unwrap()).sum();
    assert_eq!(r8, 6, "{}", sc);
}

#[test]
fn developer_single_prime_is_labelled() {
    let r = krhom(&["2,3", "--braid", "1 1 1", "--strands", "2", "--colors", "2", "--single-prime"]);
    ok(&r, "single prime");
    assert!(r.out.contains("coefficients: F_2147483647 (developer option --single-prime: not certified over Q)"), "{}", r.out);
    let u = krhom(&["--help"]);
    assert_eq!(u.code, 0);
    for hidden in ["--field", "--single-prime", "--seed-f32b"] {
        assert!(!u.out.contains(hidden), "usage shows {}", hidden);
    }
}

// ------------------------------------------------------------------ CALC-2 fixes

/// CALC-2 item 1: on the exact 2-strand route a memory-cap refusal is exit 5 naming
/// --mem-cap-mb (not a panic, exit 2), --mem-cap-mb is honoured there, and the A_λ shape limit
/// (top weight a(N - a) > 255) is exit 1 before any engine work.
#[test]
fn exact2_memory_cap_and_shape_refusals() {
    // the reviewer's example: dim A_λ = C(30,3)^2, far above the default cap
    refused(&["30", "--braid", "1 1 1", "--colors", "3"], 5, "--mem-cap-mb");
    // --mem-cap-mb applies on this route: a tiny cap refuses an input that runs by default
    let r = krhom(&["2,3", "--braid", "1 1 1", "--strands", "2", "--colors", "2", "--mem-cap-mb", "1"]);
    assert_eq!(r.code, 5, "{}\n{}", r.out, r.err);
    assert!(r.err.contains("above the memory cap of 1 MB") && !r.err.contains("please report"), "{}", r.err);
    let j = empty_lib().join("cap.json");
    let r = krhom(&["2", "--braid", "1 1 1", "--strands", "2", "--colors", "2", "--mem-cap-mb", "4000", "--json", j.to_str().unwrap()]);
    ok(&r, "cap 4000");
    assert!(std::fs::read_to_string(&j).unwrap().contains("\"mem_cap_mb\":4000"));
    // top weight: 2 * (130 - 2) = 256 > 255 (exit 1, no engine run, no panic)
    refused(&["130", "--braid", "1 1 1", "--colors", "2"], 1, "a(N - a) <= 255");
    refused(&["89", "--braid", "1 1 1", "--colors", "3"], 1, "N <= 88 for color 3");
}

/// CALC-2 item 2: N below the color is answered (zero) before any route refusal and without an
/// engine run.
#[test]
fn n_below_the_color_short_circuits_refusals() {
    for args in [
        &["1,2", "--braid", "aBaB", "--colors", "3"][..], // no library for (3,3,3)
        &["2,3", "--braid", "1 1 1", "--colors", "4"][..], // exact 2-strand colors > 3
        &["1", "--braid", "1 1 1", "--strands", "2", "--colors", "2"][..], // exact 2-strand, no xscan
        &["2,3", "--braid", "1 1 2 2 3 3 4 4", "--strands", "5", "--colors", "4"][..], // 5 strands
    ] {
        let r = krhom(args);
        ok(&r, &format!("{:?}", args));
        assert!(r.out.contains("method: none (every N is below the color"), "{:?}: {}", args, r.out);
        assert!(r.out.contains("Poincare polynomial: 0") && !r.out.contains("Poincare polynomial: t"), "{}", r.out);
        assert!(!r.err.contains("[F_2147483647]"), "{:?}: an engine ran:\n{}", args, r.err);
    }
    // an N at or above the color still goes through the refusal
    refused(&["2,3,4", "--braid", "1 1 1", "--colors", "4"], 1, "colors 1..=3");
}

/// CALC-2 items 3, 4, 5, 6: batch blocks are named before any early return and on the error line;
/// CSV names are RFC 4180-quoted; JSON/CSV are written per finished item; `input` has no duplicate key.
#[test]
fn batch_names_csv_quoting_and_incremental_outputs() {
    let d = empty_lib();
    let f = d.join("b.txt");
    std::fs::write(&f, "L6a5{0,1}; aaa\nbad; 1 0\nthree, strands; 1 2 1 -2 1\nsay \"hi\"; aBaB\n").unwrap();
    let (jb, cb) = (d.join("b.json"), d.join("b.csv"));
    let b = krhom(&["2", "--file", f.to_str().unwrap(), "--json", jb.to_str().unwrap(), "--csv", cb.to_str().unwrap(), "--colors", "1"]);
    assert_eq!(b.code, 1, "{}\n{}", b.out, b.err);
    // a parse error still gets its name block and a named error line
    assert!(b.out.contains("name: bad\nerror: bad: exit 1 (unsupported)"), "{}", b.out);
    let sc = std::fs::read_to_string(&cb).unwrap();
    assert!(sc.contains("\"L6a5{0,1}\",2,false,0,-3,1"), "{}", sc);
    assert!(sc.contains("\"say \"\"hi\"\"\",2,false,"), "{}", sc);
    for l in sc.lines().skip(1) {
        // 6 fields after RFC 4180 parsing
        let mut n = 1;
        let mut q = false;
        for c in l.chars() {
            match c {
                '"' => q = !q,
                ',' if !q => n += 1,
                _ => {}
            }
        }
        assert_eq!(n, 6, "{}", l);
    }
    let sj = std::fs::read_to_string(&jb).unwrap();
    assert!(sj.starts_with("[\n{") && sj.trim_end().ends_with(']') && sj.matches("\"tool\":\"krhom\"").count() == 4, "{}", sj);
    // `input` has each key once (the string is "text", the parsed array "braid")
    for rec in sj.split("\"tool\":\"krhom\"").skip(1) {
        let inp = &rec[rec.find("\"input\":{").unwrap()..];
        let inp = &inp[..inp.find('}').unwrap()];
        assert_eq!(inp.matches("\"braid\":").count(), inp.contains("\"braid\":[").then_some(1).unwrap_or(0), "{}", inp);
        assert!(inp.matches("\"text\":").count() == 1, "{}", inp);
    }
}

/// CALC-2 item 5: an aborted batch keeps the finished records (JSON array and CSV rows written
/// after each item).
#[test]
fn batch_outputs_survive_a_kill() {
    let d = empty_lib();
    let f = d.join("k.txt");
    // the first line is instant; the second (80 crossings, color 2: minutes) is killed
    std::fs::write(&f, format!("first; aaa\nslow; {}\n", "a".repeat(80))).unwrap();
    let (jb, cb) = (d.join("k.json"), d.join("k.csv"));
    let mut c = Command::new(env!("CARGO_BIN_EXE_krhom"));
    c.args(["5", "--file", f.to_str().unwrap(), "--colors", "2", "--json", jb.to_str().unwrap(), "--csv", cb.to_str().unwrap()]).env("KRHOM_THREADS", "1").env("KRHOM_LIB_DIR", empty_lib()).stdout(std::process::Stdio::null()).stderr(std::process::Stdio::null());
    let mut ch = c.spawn().unwrap();
    let t0 = std::time::Instant::now();
    // wait until the first record is on disk, then SIGKILL
    loop {
        let s = std::fs::read_to_string(&jb).unwrap_or_default();
        if s.contains("\"name\":\"first\"") {
            break;
        }
        assert!(t0.elapsed().as_secs() < 600, "first record not written");
        if let Some(st) = ch.try_wait().unwrap() {
            panic!("krhom ended early: {:?}", st);
        }
        std::thread::sleep(std::time::Duration::from_millis(50));
    }
    ch.kill().unwrap();
    let _ = ch.wait();
    let sj = std::fs::read_to_string(&jb).unwrap();
    assert!(sj.starts_with("[\n{") && sj.trim_end().ends_with(']') && sj.matches("\"tool\":\"krhom\"").count() == 1, "{}", sj);
    let sc = std::fs::read_to_string(&cb).unwrap();
    assert!(sc.starts_with("name,N,reduced,t,q,rank\nfirst,5,false,"), "{}", sc);
}

/// CALC-2 item 15: --table. CALC-3: --convention was removed (the user decided not to ship it);
/// krhom's q-grading stays the inverse of Bar-Natan / KnotJob / KnotInfo (aaa at N = 2 is
/// t^0q^-3 + t^0q^-1 + t^2q^-5 + t^3q^-9 here, q + q^3 + t^2q^5 + t^3q^9 there).
#[test]
fn table_and_no_convention_option() {
    let r = krhom(&["2", "--braid", "aaa", "--table"]);
    ok(&r, "table");
    assert!(!r.out.contains("convention"), "{}", r.out);
    assert!(r.out.contains("Poincare polynomial: t^0q^-3 + t^0q^-1 + t^2q^-5 + t^3q^-9"), "{}", r.out);
    assert!(r.out.contains("  table (rows t, columns q):\n  t\\q  -9  -7  -5  -3  -1\n    0               1   1\n    1\n    2           1\n    3   1\n"), "{}", r.out);
    // JSON has no convention field; default output has no table
    let d = empty_lib();
    let j1 = d.join("a.json");
    let a = krhom(&["2,3", "--braid", "1 1 1", "--strands", "2", "--colors", "2", "--json", j1.to_str().unwrap()]);
    ok(&a, "default");
    assert!(!a.out.contains("convention") && !a.out.contains("table"), "{}", a.out);
    assert!(!std::fs::read_to_string(&j1).unwrap().contains("\"convention\""));
    // --table in --file mode, zero table
    let f = d.join("t.txt");
    std::fs::write(&f, "u; 1\nz; aaa\n").unwrap();
    let r = krhom(&["1,2", "--file", f.to_str().unwrap(), "--colors", "2", "--table"]);
    ok(&r, "file table");
    assert!(r.out.contains("name: z") && r.out.contains("  table: 0") && r.out.contains("  table (rows t, columns q):"), "{}", r.out);
    // --convention is now an unknown option
    refused(&["2", "--braid", "aaa", "--convention", "knotjob"], 1, "unknown argument \"--convention\"");
    refused(&["2", "--braid", "aaa", "--convention", "krhom"], 1, "unknown argument \"--convention\"");
}

// ------------------------------------------------------------------ with the libraries (F3)

fn lib_run(args: &[&str]) -> Run {
    let mut c = Command::new(env!("CARGO_BIN_EXE_krhom"));
    c.args(args).env("KRHOM_THREADS", std::env::var("KRHOM_THREADS").unwrap_or("4".into()));
    let o = c.output().unwrap();
    Run { code: o.status.code().unwrap_or(-1), out: String::from_utf8_lossy(&o.stdout).into(), err: String::from_utf8_lossy(&o.stderr).into() }
}
fn ranks(s: &str) -> Vec<u64> {
    s.lines().filter_map(|l| l.trim().strip_prefix("total rank: ")).map(|x| x.parse().unwrap()).collect()
}

#[test]
#[ignore = "needs the l222_n4 library (tools/fetch_library.sh l222_n4); about 2 x 1.5 min"]
fn lib_41_three_strands() {
    let a = lib_run(&["2,3,4", "--braid", "1 -2 1 -2", "--strands", "3", "--colors", "2"]);
    ok(&a, "4_1");
    assert_eq!(ranks(&a.out), vec![1, 11, 46]);
    assert!(a.out.contains("t^-2q^6 + t^-2q^8 + t^-1q^0 + t^-1q^2 + t^0q^-2 + t^0q^0 + t^0q^2 + t^1q^-2 + t^1q^0 + t^2q^-8 + t^2q^-6"), "{}", a.out);
    let r = lib_run(&["2,3,4", "--braid", "aBaB", "--colors", "2", "--reduced"]);
    ok(&r, "4_1 reduced");
    assert_eq!(ranks(&r.out), vec![1, 5, 25]);
    assert!(r.out.contains("t^-2q^6 + t^-1q^2 + t^0q^0 + t^1q^-2 + t^2q^-6"), "{}", r.out);
}

#[test]
#[ignore = "needs the l222_n4 library; about 2 x 1 min"]
fn lib_stabilized_trefoil() {
    for red in [false, true] {
        let mut a = vec!["2,3,4", "--braid", "1 1 1 2", "--strands", "3", "--colors", "2", "--no-simplify"];
        let mut b = vec!["2,3,4", "--braid", "1 1 1", "--strands", "2", "--colors", "2"];
        if red {
            a.push("--reduced");
            b.push("--reduced");
        }
        let (x, y) = (lib_run(&a), lib_run(&b));
        ok(&x, "stabilized");
        ok(&y, "2-strand");
        assert_eq!(blocks(&x.out), blocks(&y.out));
        assert_eq!(ranks(&x.out), if red { vec![1, 3, 9] } else { vec![1, 7, 22] });
    }
}

#[test]
#[ignore = "needs the l2222_n4 library (tools/fetch_library.sh l2222_n4); 4 strands, several minutes"]
fn lib_four_strands_markov() {
    let a = lib_run(&["2,3,4", "--braid", "1 2 3", "--strands", "4", "--colors", "2", "--no-simplify"]);
    ok(&a, "1 2 3");
    assert_eq!(ranks(&a.out), vec![1, 3, 6]);
    let r = lib_run(&["2,3,4", "--braid", "1 1 1 2 3", "--strands", "4", "--colors", "2", "--no-simplify", "--reduced"]);
    let y = lib_run(&["2,3,4", "--braid", "1 1 1", "--strands", "2", "--colors", "2", "--reduced"]);
    ok(&r, "1 1 1 2 3 reduced");
    assert_eq!(blocks(&r.out), blocks(&y.out));
}

#[test]
#[ignore = "needs the l2222_n4 library; --check-library on 6_1 and 8_15 (minutes)"]
fn lib_check_library_4_strands() {
    let a = lib_run(&["2,3,4", "--braid", "1 1 2 -1 -3 2 -3", "--strands", "4", "--colors", "2", "--check-library"]);
    ok(&a, "6_1 check-library");
    assert!(a.out.contains("check-library: library l2222_n4 covers this word"), "{}", a.out);
    let b = lib_run(&["2,3,4", "--braid", "1 1 -2 1 3 2 2 2 3", "--strands", "4", "--colors", "2", "--check-library"]);
    assert_eq!(b.code, 4, "{}\n{}", b.out, b.err);
    for p in [926, 949, 986, 997, 1054, 1412, 1435, 1471, 1472, 1483, 1540] {
        assert!(b.err.contains(&format!("pair {} ", p)), "8_15: pair {} not listed:\n{}", p, b.err);
    }
}

// ------------------------------------------------------------------ F-CLI

/// Terms (t, q) -> k of every "Poincare polynomial:" line, in order.
fn polys(s: &str) -> Vec<std::collections::BTreeMap<(i32, i32), i64>> {
    s.lines()
        .filter_map(|l| l.trim().strip_prefix("Poincare polynomial: "))
        .map(|p| {
            let mut m = std::collections::BTreeMap::new();
            if p.trim() == "0" {
                return m;
            }
            for term in p.split(" + ") {
                let ti = term.find('t').unwrap();
                let k: i64 = if ti == 0 { 1 } else { term[..ti].parse().unwrap() };
                let qi = term.find('q').unwrap();
                let t: i32 = term[ti + 2..qi].parse().unwrap();
                let q: i32 = term[qi + 2..].parse().unwrap();
                m.insert((t, q), k);
            }
            m
        })
        .collect()
}
/// Graded Euler characteristic q -> sum (-1)^t k.
fn euler(p: &std::collections::BTreeMap<(i32, i32), i64>) -> std::collections::BTreeMap<i32, i64> {
    let mut e = std::collections::BTreeMap::new();
    for (&(t, q), &k) in p {
        *e.entry(q).or_insert(0) += if t.rem_euclid(2) == 0 { k } else { -k };
    }
    e.retain(|_, v| *v != 0);
    e
}
/// The quantum dimension [N choose a]_q (symmetric, step 2).
fn qbinom(n: usize, a: usize) -> std::collections::BTreeMap<i32, i64> {
    // subsets of size a of {0..n-1}: q-degree 2 * (sum - a(a-1)/2) - a(n-a)
    let mut m = std::collections::BTreeMap::new();
    for mask in 0u32..(1 << n) {
        if mask.count_ones() as usize != a {
            continue;
        }
        let s: i32 = (0..n).filter(|i| mask >> i & 1 == 1).map(|i| i as i32).sum();
        *m.entry(2 * (s - (a * (a - 1) / 2) as i32) - (a * (n - a)) as i32).or_insert(0) += 1;
    }
    m
}
fn lmul(a: &std::collections::BTreeMap<i32, i64>, b: &std::collections::BTreeMap<i32, i64>) -> std::collections::BTreeMap<i32, i64> {
    let mut m = std::collections::BTreeMap::new();
    for (x, u) in a {
        for (y, v) in b {
            *m.entry(x + y).or_insert(0) += u * v;
        }
    }
    m.retain(|_, v| *v != 0);
    m
}

/// F-CLI item 1: limits are refused with exit 1 before any engine work.
#[test]
fn fcli_limits_refused_before_engine_work() {
    for (args, needle) in [
        (&["33", "--braid", "aaa"][..], "N <= 32"),
        (&["2", "--braid", "1 2 3 4 5 6 7 8 9 10"][..], "at most 10 strands"),
        (&["33", "--braid", "aaa", "--colors", "2", "--colored-method", "cable"][..], "N <= 32"),
        (&["36", "--braid", "aaa", "--colors", "3"][..], "above its limit 5e7"),
        (&["120", "--braid", "aaa", "--colors", "2"][..], "above its limit 5e7"),
        (&["2,3,40", "--braid", "1 1", "--colors", "3,3"][..], "N = 40 with labels [3, 3]"),
    ] {
        let r = krhom(args);
        assert_eq!(r.code, 1, "{:?}: exit {}\n{}", args, r.code, r.err);
        assert!(r.err.contains(needle), "{:?}: {}", args, r.err);
        assert!(!r.err.contains("[F_2147483647]") && !r.err.contains("panicked"), "{:?}: an engine ran:\n{}", args, r.err);
    }
    // just inside the limits: the exact closure has no N cap; reduced color 3 at N = 36 has dim C(36,3)
    let r = krhom(&["33", "--braid", "aaa", "--closure", "exact"]);
    ok(&r, "N = 33 exact closure");
    let r = krhom(&["36", "--braid", "aaa", "--colors", "3", "--reduced", "--mem-cap-mb", "1"]);
    assert_eq!(r.code, 5, "reduced passes the dimension check and meets the memory cap: {}", r.err);
}

/// F-CLI item 2: --mark K (Hopf link with colors (2,1); a 3-component chain; Borromean rings).
#[test]
fn fcli_mark_component() {
    let run = |extra: &[&str]| {
        let mut a = vec!["2,3,4", "--braid", "1 1"];
        a.extend_from_slice(extra);
        let r = krhom(&a);
        ok(&r, &format!("{:?}", a));
        r
    };
    let unred = run(&["--colors", "2,1"]);
    let red = run(&["--colors", "2,1", "--reduced"]);
    let m0 = run(&["--colors", "2,1", "--reduced", "--mark", "0"]);
    let m1 = run(&["--colors", "2,1", "--reduced", "--mark", "1"]);
    let swapped = run(&["--colors", "1,2", "--reduced"]);
    assert_eq!(blocks(&m0.out), blocks(&red.out));
    assert!(m0.out.contains("reduced (component 0 marked"), "{}", m0.out);
    assert!(m1.out.contains("reduced (component 1 marked") && m1.out.contains("conjugated braid: 1 1  (component 1 moved from bottom position 1 to 0"), "{}", m1.out);
    // marking the color-1 component = putting it first
    assert_eq!(blocks(&m1.out), blocks(&swapped.out));
    assert_ne!(blocks(&m0.out), blocks(&m1.out));
    // hand-checkable: N = 2 marked color 2 (wedge^2 C^2 is a line) = unreduced; marked color 1 = q^-2
    assert_eq!(polys(&m0.out)[0], polys(&unred.out)[0]);
    assert!(m1.out.contains("N = 2\n  Poincare polynomial: t^0q^-2\n"), "{}", m1.out);
    // chi(unreduced) = chi(reduced, K marked) * qdim(wedge^{c_K} C^N) for every N
    let (pu, p0, p1) = (polys(&unred.out), polys(&m0.out), polys(&m1.out));
    for (i, n) in [2usize, 3, 4].iter().enumerate() {
        assert_eq!(euler(&pu[i]), lmul(&euler(&p0[i]), &qbinom(*n, 2)), "N = {} mark 0", n);
        assert_eq!(euler(&pu[i]), lmul(&euler(&p1[i]), &qbinom(*n, 1)), "N = {} mark 1", n);
    }
    // ranks too (the Hopf homology is in even t): 2 * 1, 9 = 3 * 3, 24 = 4 * 6 = 6 * 4
    assert_eq!(ranks(&unred.out), vec![2, 9, 24]);
    assert_eq!(ranks(&m0.out), vec![2, 3, 4]);
    assert_eq!(ranks(&m1.out), vec![1, 3, 6]);
    // 3-component chain (components 0, 1, 2; 1 in the middle), uncolored: every mark satisfies
    // the Euler identity; the ends are symmetric
    let u = krhom(&["2,3", "--braid", "1 1 2 2"]);
    ok(&u, "chain");
    let mut outs = vec![];
    for k in ["0", "1", "2"] {
        let r = krhom(&["2,3", "--braid", "1 1 2 2", "--reduced", "--mark", k]);
        ok(&r, &format!("chain mark {}", k));
        let (pu, pr) = (polys(&u.out), polys(&r.out));
        for (i, n) in [2usize, 3].iter().enumerate() {
            assert_eq!(euler(&pu[i]), lmul(&euler(&pr[i]), &qbinom(*n, 1)), "chain mark {} N = {}", k, n);
        }
        outs.push(blocks(&r.out));
    }
    assert_eq!(outs[0], outs[2]);
    // Borromean rings from the link table: all marks agree (symmetry)
    let b: Vec<String> = ["0", "1", "2"].iter().map(|k| blocks(&krhom(&["2,3", "--link", "L6a4{0,0}", "--reduced", "--mark", k]).out)).collect();
    assert!(!b[0].is_empty() && b[0] == b[1] && b[1] == b[2], "{:?}", b);
    // refusals
    refused(&["2", "--braid", "1 1", "--mark", "1"], 1, "add --reduced");
    refused(&["2", "--braid", "aaa", "--reduced", "--mark", "1"], 1, "only --mark 0");
    refused(&["2", "--braid", "1 1", "--reduced", "--mark", "2"], 1, "the link has 2 components");
    refused(&["2", "--braid", "1 1", "--reduced", "--mark", "x"], 1, "expected a component number");
    ok(&krhom(&["2", "--braid", "aaa", "--reduced", "--mark", "0"]), "knot mark 0");
}

/// F-CLI item 3: --knot / --link NAME from the shipped tables.
#[test]
fn fcli_knot_and_link_names() {
    let a = krhom(&["2,3", "--knot", "4_1"]);
    let b = krhom(&["2,3", "--braid", "1 -2 1 -2"]);
    ok(&a, "4_1");
    ok(&b, "1 -2 1 -2");
    assert_eq!(blocks(&a.out), blocks(&b.out));
    assert!(a.out.starts_with("name: 4_1\nsource: knot 4_1 from data/knots_upto12.jsonl (KnotInfo braid; KnotInfo chirality)\n"), "{}", a.out);
    // KnotInfo chirality: 3_1 is the braid [1,1,1]
    assert_eq!(blocks(&krhom(&["2", "--knot", "3_1"]).out), blocks(&krhom(&["2", "--braid", "aaa"]).out));
    // oriented links: L2a1{1} = sigma_1^2, L2a1{0} = sigma_1^-2
    let l = krhom(&["2", "--link", "L2a1{1}"]);
    ok(&l, "L2a1{1}");
    assert!(l.out.contains("source: link L2a1{1} from data/links_upto11.jsonl"), "{}", l.out);
    assert_eq!(blocks(&l.out), blocks(&krhom(&["2", "--braid", "1 1"]).out));
    assert_eq!(blocks(&krhom(&["2", "--link", "L2a1{0}"]).out), blocks(&krhom(&["2", "--braid", "-1 -1"]).out));
    // case-insensitive unique match
    ok(&krhom(&["2", "--link", "l2a1{1}"]), "lower case");
    refused(&["2", "--link", "L2a1"], 1, "LinkInfo lists 2 orientations: L2a1{0}, L2a1{1}");
    refused(&["2", "--knot", "4_2"], 1, "close matches: 4_1,");
    refused(&["2", "--knot", "L2a1{1}"], 1, "use --link L2a1{1}");
    refused(&["2", "--knot", "4_1", "--strands", "4"], 1, "drop --strands");
    refused(&["2", "--knot", "4_1", "--braid", "aaa"], 1, "exactly one of");
    // the JSON records the name as the input
    let j = empty_lib().join("k.json");
    ok(&krhom(&["2", "--knot", "4_1", "--json", j.to_str().unwrap()]), "json");
    let s = std::fs::read_to_string(&j).unwrap();
    assert!(s.contains("\"input\":{\"kind\":\"knot\",\"text\":\"4_1\"") && s.contains("\"braid\":[1, -2, 1, -2]"), "{}", s);
}

/// The table parser reads every record of both tables.
#[test]
fn fcli_tables_parse_completely() {
    let d = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("data");
    for (f, n) in [("knots_upto12.jsonl", 2977), ("links_upto11.jsonl", 4188)] {
        let t = std::fs::read_to_string(d.join(f)).unwrap();
        let mut k = 0;
        for l in t.lines() {
            // the shape the CLI parser relies on
            assert!(l.starts_with("{\"name\":\"") && l.contains("\"braid\":{\"strands\":") && l.contains(",\"word\":["), "{}", &l[..80]);
            k += 1;
        }
        assert_eq!(k, n, "{}", f);
    }
    // the last knot runs through the lookup
    let r = krhom(&["1", "--knot", "12n_888"]);
    ok(&r, "12n_888");
    assert!(r.out.contains("braid: 1 1 1 2 2 2 1 1 1 2 2 2"), "{}", r.out);
}

/// F-CLI item 4: --latex.
#[test]
fn fcli_latex() {
    let d = empty_lib();
    let tex = d.join("a.tex");
    let r = krhom(&["2,3", "--knot", "3_1", "--colors", "2", "--reduced", "--latex", tex.to_str().unwrap()]);
    ok(&r, "latex");
    let s = std::fs::read_to_string(&tex).unwrap();
    assert!(s.contains("\\documentclass{article}") && s.contains("\\usepackage{graphicx}") && s.trim_end().ends_with("\\end{document}"), "{}", s);
    assert_eq!(s.matches("\\begin{tabular}").count(), 2);
    assert_eq!(s.matches("\\end{tabular}").count(), 2);
    assert_eq!(s.matches("\\begin{table}").count(), s.matches("\\end{table}").count());
    assert!(s.contains("\\section*{knot 3\\_1}"), "{}", s);
    assert!(s.contains("\\caption{knot 3\\_1: braid $1\\,1\\,1$ on 2 strands, colors [2], $N = 3$, reduced (strand at bottom position 0 marked); coefficients Q;"), "{}", s);
    // N = 3 table: rows t = 0..3, columns q = -12..-4 step 2
    assert!(s.contains("$t \\backslash q$ & $-12$ & $-10$ & $-8$ & $-6$ & $-4$ \\\\\n\\hline\n$0$ &  &  &  &  & 1 \\\\\n$1$ &  &  &  &  &  \\\\\n$2$ &  &  & 1 &  &  \\\\\n$3$ & 1 &  &  &  &  \\\\\n"), "{}", s);
    // --file: one section per input, errors and zero results included; braces escaped
    let f = d.join("f.txt");
    std::fs::write(&f, "tref; aaa\nhopf_{x}; 1 1\nbad; 1 0\n").unwrap();
    let tex2 = d.join("f.tex");
    let r = krhom(&["1,2", "--file", f.to_str().unwrap(), "--colors", "2", "--latex", tex2.to_str().unwrap()]);
    assert_eq!(r.code, 1);
    let s = std::fs::read_to_string(&tex2).unwrap();
    assert_eq!(s.matches("\\section*{").count(), 3, "{}", s);
    assert!(s.contains("\\section*{hopf\\_\\{x\\} (braid 1 1)}") && s.contains("Not computed: exit 1 (unsupported)") && s.contains("(the homology is 0)"), "{}", s);
    refused(&["2", "--braid", "aaa", "--latex", "/nonexistent/dir/x.tex"], 1, "--latex /nonexistent/dir/x.tex");
}

/// F-CLI item 5: --estimate outside the box route is a note, nothing is computed.
#[test]
fn fcli_estimate_fast_routes() {
    for args in [&["2,3", "--braid", "aaa", "--colors", "2", "--estimate"][..], &["2,3", "--braid", "aBaB", "--estimate"][..]] {
        let r = krhom(args);
        ok(&r, &format!("{:?}", args));
        assert!(r.out.contains("estimate: fast route (") && r.out.contains("no estimate needed"), "{}", r.out);
        assert!(!r.out.contains("Poincare") && !r.err.contains("[F_2147483647]"), "{}\n{}", r.out, r.err);
    }
    // the box route still needs the library (exit 4 before any estimate)
    refused(&["2,3,4", "--braid", "aBaB", "--colors", "2", "--estimate"], 4, "tools/fetch_library.sh l222_n4");
    refused(&["2", "--braid", "aBaB", "--colors", "2", "--estimate", "--check-library"], 1, "exclusive");
}

#[test]
#[ignore = "needs the l222_n4 library; about 2 x 1 min"]
fn lib_fcli_knot_name_and_estimate() {
    let a = lib_run(&["2,3,4", "--knot", "4_1", "--colors", "2", "--reduced"]);
    let b = lib_run(&["2,3,4", "--braid", "1 -2 1 -2", "--colors", "2", "--reduced"]);
    ok(&a, "--knot 4_1");
    ok(&b, "--braid 4_1");
    assert_eq!(blocks(&a.out), blocks(&b.out));
    let e = lib_run(&["2,3,4", "--knot", "4_1", "--colors", "2", "--reduced", "--estimate"]);
    ok(&e, "estimate");
    assert!(e.out.contains("estimate: F_2147483647: ") && e.out.contains("box points (reduced grid)") && e.out.contains("estimate: total: about"), "{}", e.out);
    assert!(!e.out.contains("Poincare"), "{}", e.out);
}
