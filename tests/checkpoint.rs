//! F-CKPT: `krhom --checkpoint DIR` (src/calc_ckpt.rs). The default tests check the command-line
//! refusals (no library needed); the `#[ignore]` tests need l222_n4 installed
//! (`tools/fetch_library.sh l222_n4`) and interrupt a 3-strand run with the test hook
//! `KRHOM_CKPT_STOP_AFTER`, resume it, compare with the plain run and refuse a changed input.
use std::path::PathBuf;
use std::process::Command;

struct Run {
    code: i32,
    out: String,
    err: String,
}

fn tmp(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("krhom_ckpt_{}_{}", tag, std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    d
}

fn krhom(args: &[&str], env: &[(&str, &str)]) -> Run {
    let mut c = Command::new(env!("CARGO_BIN_EXE_krhom"));
    c.args(args).env("KRHOM_THREADS", std::env::var("KRHOM_THREADS").unwrap_or("4".into()));
    for v in ["KRHOM_CKPT_SECS", "KRHOM_CKPT_STOP_AFTER", "KRHOM_NO_M5I", "KRHOM_FORCE_LADDER"] {
        c.env_remove(v);
    }
    for (k, v) in env {
        c.env(k, v);
    }
    let o = c.output().expect("cannot run krhom");
    Run { code: o.status.code().unwrap_or(-1), out: String::from_utf8_lossy(&o.stdout).into(), err: String::from_utf8_lossy(&o.stderr).into() }
}

/// stdout without the `time` line (the only line that may differ between runs)
fn untimed(s: &str) -> String {
    s.lines().filter(|l| !l.starts_with("time ")).collect::<Vec<_>>().join("\n")
}

#[test]
fn checkpoint_refusals_exit_1() {
    let d = tmp("refuse");
    std::fs::create_dir_all(&d).unwrap();
    let f = d.join("b.txt");
    std::fs::write(&f, "a; aaa\nb; aBaB\n").unwrap();
    let dir = d.join("ck");
    let r = krhom(&["2", "--file", f.to_str().unwrap(), "--colors", "2", "--checkpoint", dir.to_str().unwrap()], &[("KRHOM_LIB_DIR", d.to_str().unwrap())]);
    assert_eq!(r.code, 1, "{}\n{}", r.out, r.err);
    assert!(r.err.contains("--checkpoint takes one braid"), "{}", r.err);
    let r = krhom(&["2,3,4", "--braid", "1 -2 1 -2", "--colors", "2", "--check-library", "--checkpoint", dir.to_str().unwrap()], &[("KRHOM_LIB_DIR", d.to_str().unwrap())]);
    assert_eq!(r.code, 1, "{}\n{}", r.out, r.err);
    assert!(r.err.contains("--checkpoint with --check-library"), "{}", r.err);
    assert!(!dir.exists());
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
#[ignore = "needs the l222_n4 library (tools/fetch_library.sh l222_n4); about 1 min"]
fn lib_checkpoint_resume_equals_plain_run() {
    let dir = tmp("resume");
    let ds = dir.to_str().unwrap();
    let args = ["2,3,4", "--braid", "1 -2 1 -2", "--colors", "2", "--reduced"];
    let mut ck: Vec<&str> = args.to_vec();
    ck.extend(["--checkpoint", ds]);
    let plain = krhom(&args, &[]);
    assert_eq!(plain.code, 0, "{}\n{}", plain.out, plain.err);
    // interrupted after 20 of the 36 points of F31 (checkpoint after every batch of 16)
    let a = krhom(&ck, &[("KRHOM_CKPT_SECS", "0"), ("KRHOM_CKPT_STOP_AFTER", "20")]);
    assert_eq!(a.code, 75, "{}\n{}", a.out, a.err);
    assert!(a.err.contains("36 box points, 0 already done, 20 scanned now"), "{}", a.err);
    assert!(dir.join("f31/manifest.txt").is_file() && dir.join("f31/seed201/run.hdr").is_file() && dir.join("f31/seed201/w0_s4_r0.acc").is_file());
    assert!(!dir.join("f32b").exists());
    // resumed: the 16 checkpointed points are not scanned again
    let b = krhom(&ck, &[]);
    assert_eq!(b.code, 0, "{}\n{}", b.out, b.err);
    assert!(b.err.contains("16 already done (resumed), 20 scanned now"), "{}", b.err);
    assert!(b.err.contains("checkpoint") && b.err.contains("f32b/seed201: 36 box points, 0 already done, 36 scanned now"), "{}", b.err);
    assert_eq!(untimed(&b.out), untimed(&plain.out));
    // a finished directory: nothing scanned, same output
    let c = krhom(&ck, &[]);
    assert_eq!(c.code, 0, "{}\n{}", c.out, c.err);
    assert_eq!(c.err.matches("36 already done (resumed), 0 scanned now").count(), 2, "{}", c.err);
    assert_eq!(untimed(&c.out), untimed(&plain.out));
    // another input on the same directory is refused before any work (exit 1)
    for (other, key) in [(vec!["2,3,4", "--braid", "1 -2 1 -2", "--colors", "2", "--checkpoint", ds], "mode"), (vec!["2,3,4", "--braid", "1 -2 1 -2 1 -2", "--colors", "2", "--reduced", "--checkpoint", ds], "braid"), (vec!["2,3", "--braid", "1 -2 1 -2", "--colors", "2", "--reduced", "--checkpoint", ds], "N"), (vec!["2,3,4", "--braid", "1 -2 1 -2", "--colors", "2", "--reduced", "--seed", "7", "--checkpoint", ds], "seed")] {
        let r = krhom(&other, &[]);
        assert_eq!(r.code, 1, "{:?}\n{}\n{}", other, r.out, r.err);
        assert!(r.err.contains("resume refused") && r.err.contains(&format!("\n  {}: directory has", key)), "{:?}: {}", other, r.err);
    }
    let _ = std::fs::remove_dir_all(&dir);
}
