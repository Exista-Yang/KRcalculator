//! libbuild (M5b): plan counts of (3,3,3) and a small run → replay → resume cycle.
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_libbuild")
}

fn run(args: &[&str]) -> String {
    let o = Command::new(bin()).args(args).env("KRHOM_THREADS", "1").output().expect("libbuild");
    assert!(o.status.success(), "libbuild {:?} failed: {}", args, String::from_utf8_lossy(&o.stderr));
    String::from_utf8_lossy(&o.stdout).to_string()
}

#[test]
fn plan_333_sets() {
    let o = run(&["plan", "--lam", "3,3,3", "--n", "4"]);
    assert!(o.contains("set level-n4: 47 jobs {\"closure\": 5, \"pair\": 42}; known targets 29 (Σ mult 19)"), "{}", o);
    let o = run(&["plan", "--lam", "3,3,3", "--n", "6", "--corner", "markov-ext"]);
    assert!(o.contains("32 jobs {\"closure\": 16, \"cross\": 16}"), "{}", o);
    let o = run(&["plan", "--lam", "3,3,3", "--n", "6", "--sample", "largex"]);
    assert!(o.contains("20 jobs {\"pair\": 20}; known targets 188 (Σ mult 255); max X 21600 flows"), "{}", o);
}

#[test]
fn run_replay_resume_33() {
    let dir = std::env::temp_dir().join(format!("libbuild_test_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let cache = dir.join("cache");
    let cache = cache.to_str().unwrap();
    let out = dir.join("run.jsonl");
    let out = out.to_str().unwrap();
    let common = ["--lam", "3,3", "--n", "5", "--field", "31", "--cache", cache, "--jobs", "2", "--wall", "300", "--rss", "2048"];
    let mut a = vec!["run"];
    a.extend(common);
    a.extend(["--json", out]);
    let o = run(&a);
    assert!(o.contains("jobs certified") && o.contains("(100.0%)"), "{}", o);
    // resume: everything certified in OUT is skipped
    let o = run(&a);
    assert!(o.contains("(100.0%)"), "{}", o);
    let t = std::fs::read_to_string(out).unwrap();
    let last_header = t.lines().filter(|l| l.starts_with("{\"header\"")).last().unwrap();
    assert!(!last_header.contains("\"jobs_skipped_done\":0"), "{}", last_header);
    // replay: a provider that refuses every call; everything loads from layer A and is re-checked
    let rout = dir.join("replay.jsonl");
    let mut r = vec!["replay"];
    r.extend(common);
    r.extend(["--json", rout.to_str().unwrap()]);
    let o = run(&r);
    assert!(o.contains("(100.0%)"), "{}", o);
    let rt = std::fs::read_to_string(&rout).unwrap();
    assert!(rt.contains("\"layer_a_loaded\":") && !rt.contains("\"layer_a_saved\":1"), "replay saved new files");
    // replay with an empty layer A must fail (the provider refuses)
    std::fs::remove_dir_all(dir.join("cache").join("layer_a")).unwrap();
    std::fs::create_dir_all(dir.join("cache").join("layer_a")).unwrap();
    let rout2 = dir.join("replay2.jsonl");
    let mut r2 = vec!["replay"];
    r2.extend(common);
    r2.extend(["--json", rout2.to_str().unwrap()]);
    let o = run(&r2);
    assert!(!o.contains("(100.0%)") && o.contains("FAIL"), "{}", o);
    let _ = std::fs::remove_dir_all(&dir);
}

/// M6a-2a hardening (review wave 3 m5_lib F1-F4): a replay writing into the run's own output
/// file does not resume from the run's lines; certified counts are restricted to the current set;
/// replay never runs the exact solver (a missing template is refused, not recomputed); every
/// attempt records the commit of the binary that ran it; the default output name carries the
/// mode.
#[test]
fn resume_governing_header_and_replay_only() {
    let dir = std::env::temp_dir().join(format!("libbuild_m6a2a_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let cache = dir.join("cache");
    let cache = cache.to_str().unwrap();
    let out = dir.join("shared.jsonl");
    let out = out.to_str().unwrap();
    let common = ["--lam", "3,3", "--n", "5", "--field", "31", "--cache", cache, "--jobs", "2", "--wall", "300", "--rss", "2048"];
    let mut a = vec!["run"];
    a.extend(common);
    a.extend(["--json", out]);
    let o = run(&a);
    assert!(o.contains("(100.0%)"), "{}", o);
    let t = std::fs::read_to_string(out).unwrap();
    let header = t.lines().find(|l| l.starts_with("{\"header\"")).unwrap();
    assert!(header.contains("\"governing\":\"mode=run "), "{}", header);
    let commit = header.split("\"commit\":\"").nth(1).unwrap().split('"').next().unwrap().to_string();
    assert!(commit != "unknown" && !commit.is_empty());
    for l in t.lines().filter(|l| l.starts_with("{\"attempt\"")) {
        assert!(l.contains(&format!("\"commit\":\"{}\"", commit)), "attempt without the binary's commit: {}", l);
        assert!(l.contains("\"tmpl_computed\":"), "{}", l);
    }
    // an ok attempt line of a job outside the set, under the same header: not counted
    let fake = "{\"attempt\":1,\"job\":\"pair 99999\",\"scale\":1,\"outcome\":\"ok\",\"result\":null}\n";
    std::fs::write(out, format!("{}{}", t, fake)).unwrap();
    let o = run(&a);
    assert!(o.contains("(100.0%)"), "certified count must stay over the current set: {}", o);
    // replay into the SAME file: it must not resume from the run's lines
    let mut r = vec!["replay"];
    r.extend(common);
    r.extend(["--json", out]);
    let o = run(&r);
    assert!(o.contains("(100.0%)"), "{}", o);
    let t = std::fs::read_to_string(out).unwrap();
    let rh = t.lines().filter(|l| l.starts_with("{\"header\"")).last().unwrap();
    assert!(rh.contains("\"mode\":\"replay\"") && rh.contains("\"jobs_skipped_done\":0"), "replay resumed from run lines: {}", rh);
    // replay with the template cache gone: refused (never recomputed by the exact solver)
    let tdir = std::fs::read_dir(cache).unwrap().map(|e| e.unwrap().path()).filter(|p| p.is_dir() && !p.ends_with("layer_a")).collect::<Vec<_>>();
    assert!(!tdir.is_empty(), "no template cache directory in {}", cache);
    for d in &tdir {
        for f in std::fs::read_dir(d).unwrap() {
            let f = f.unwrap().path();
            if f.file_name().unwrap().to_string_lossy().starts_with("hom") {
                std::fs::remove_file(f).unwrap();
            }
        }
    }
    let rout = dir.join("replay_notmpl.jsonl");
    let mut r2 = vec!["replay"];
    r2.extend(common);
    r2.extend(["--json", rout.to_str().unwrap()]);
    let o = run(&r2);
    assert!(!o.contains("(100.0%)") && o.contains("FAIL"), "replay recomputed missing templates: {}", o);
    let rt = std::fs::read_to_string(&rout).unwrap();
    assert!(!rt.contains("\"tmpl_computed\":1") && !rt.lines().filter(|l| l.contains("\"tmpl_computed\":")).any(|l| !l.contains("\"tmpl_computed\":0")), "the exact solver ran in replay");
    // default output name: carries the mode
    let o2 = std::process::Command::new(bin()).args(["run", "--lam", "3,3", "--n", "5", "--field", "31", "--cache", cache, "--jobs", "2", "--max-jobs", "1"]).current_dir(&dir).env("KRHOM_THREADS", "1").output().unwrap();
    assert!(o2.status.success());
    let names: Vec<String> = std::fs::read_dir(&dir).unwrap().map(|e| e.unwrap().file_name().to_string_lossy().to_string()).collect();
    assert!(names.iter().any(|n| n.starts_with("libbuild_33_") && n.ends_with("_f31_run.jsonl")), "{:?}", names);
    let _ = std::fs::remove_dir_all(&dir);
}
