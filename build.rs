//! Embeds the commit of the source tree into the binaries (KRHOM_BUILD_COMMIT), so that every
//! record names the code that actually ran (libbuild review F4, xscan review m2).
use std::process::Command;

fn main() {
    let run = |a: &[&str]| Command::new("git").args(a).output().ok().filter(|o| o.status.success()).map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string());
    let commit = match run(&["rev-parse", "HEAD"]) {
        Some(h) => {
            let dirty = run(&["status", "--porcelain", "--untracked-files=no"]).map(|x| !x.is_empty()).unwrap_or(false);
            if dirty { format!("{}-dirty", h) } else { h }
        }
        None => "unknown".into(),
    };
    println!("cargo:rustc-env=KRHOM_BUILD_COMMIT={}", commit);
    println!("cargo:rerun-if-changed=.git/HEAD");
    println!("cargo:rerun-if-changed=.git/refs/heads");
    println!("cargo:rerun-if-changed=src");
    println!("cargo:rerun-if-changed=build.rs");
}
