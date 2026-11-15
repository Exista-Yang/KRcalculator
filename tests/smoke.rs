//! Smoke test of the calculator binary (the library-free part of tools/smoke.sh): known total ranks
//! of the uncolored figure-eight and the color-2 trefoil. Runs in well under a second (release).

use std::process::Command;

fn total_ranks(args: &[&str]) -> Vec<u64> {
    let out = Command::new(env!("CARGO_BIN_EXE_krhom"))
        .args(args)
        .output()
        .expect("run krhom");
    assert!(
        out.status.success(),
        "krhom {:?} exited with {:?}: {}",
        args,
        out.status.code(),
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout)
        .expect("utf-8 stdout")
        .lines()
        .filter_map(|l| l.strip_prefix("  total rank: "))
        .map(|r| r.trim().parse().expect("rank"))
        .collect()
}

#[test]
fn uncolored_figure_eight() {
    assert_eq!(total_ranks(&["2,3", "--braid", "aBaB"]), vec![6, 11]);
}

#[test]
fn trefoil_color_2() {
    assert_eq!(total_ranks(&["2,3,4", "--braid", "aaa", "--colors", "2"]), vec![1, 7, 22]);
}

#[test]
fn trefoil_color_2_reduced() {
    assert_eq!(
        total_ranks(&["2,3,4", "--braid", "aaa", "--colors", "2", "--reduced"]),
        vec![1, 3, 9]
    );
}
