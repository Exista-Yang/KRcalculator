//! M4c B4w: carry a layer-A cache built with the production chain (M5c) over to the M4c
//! provider stamp (`leaves::provider_for_stamped_m4c`), without editing sindec.rs or
//! libbuild.rs.
//!
//! The M4c provider is the production chain plus extra routes (Stage 0 product closure, Leaves
//! closure, Leaves pairs); its stamp is the chain stamp plus a suffix. A layer-A file written
//! under the chain stamp is accepted by the M4c run only if its header carries the M4c stamp.
//! This tool rewrites the header of every layer-A file whose provider token equals the chain
//! stamp up to the cost-table token (`sindec::route_header`, as `--stamp-route` compares) by
//! appending the M4c suffix; the key and the body are unchanged. Files with any other header are
//! left alone (and reported).
//!
//! Soundness is the route-stamp argument of M5c §3: the provider stamp only steers which calls
//! a provider attempts; every loaded map is re-certified on load exactly as a fresh one (count,
//! ends, degrees, D14 scaling; Π∘Ι = id, form and genuineness; `certify_closure`), and a replay
//! with a provider that refuses every call checks that. Use with `libbuild --stamp-route`.
//!
//! usage: b4wstamp --cache DIR --lam 3,3,3 --costs PATH --provider m4c|m4c-lf [--dry]
use krhom::field::F31;
use std::path::PathBuf;

fn opt(args: &[String], k: &str) -> Option<String> {
    args.iter().position(|a| a == k).and_then(|i| args.get(i + 1).cloned())
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let cache = PathBuf::from(opt(&args, "--cache").expect("--cache"));
    let lam: Vec<u8> = opt(&args, "--lam").unwrap_or("3,3,3".into()).split(',').map(|x| x.parse().unwrap()).collect();
    let costs = opt(&args, "--costs").map(PathBuf::from);
    let pv = opt(&args, "--provider").unwrap_or("m4c".into());
    assert!(pv == "m4c" || pv == "m4c-lf", "--provider m4c|m4c-lf");
    let dry = args.iter().any(|a| a == "--dry");
    let (_, chain) = krhom::sindec::provider_for_stamped::<F31>(&lam, costs.as_deref(), 1.0);
    let (_, m4c) = krhom::leaves::provider_for_stamped_m4c::<F31>(&lam, costs.as_deref(), 1.0, pv == "m4c-lf");
    let suffix = m4c.strip_prefix(chain.as_str()).expect("the M4c stamp extends the chain stamp").to_string();
    eprintln!("chain stamp: {}\nm4c suffix: {}", chain, suffix);
    let dir = cache.join("layer_a");
    let (mut done, mut already, mut other) = (0usize, 0usize, 0usize);
    let mut names: Vec<PathBuf> = std::fs::read_dir(&dir).expect("layer_a dir").filter_map(|e| e.ok()).map(|e| e.path()).filter(|p| p.extension().is_some_and(|x| x == "txt")).collect();
    names.sort();
    for p in names {
        let text = std::fs::read_to_string(&p).unwrap();
        let (h, body) = text.split_once('\n').unwrap_or((text.as_str(), ""));
        let (Some(i), Some(j)) = (h.find(" prov="), h.find(" key=")) else {
            other += 1;
            eprintln!("skip (no prov/key): {}", p.display());
            continue;
        };
        let prov = &h[i + 6..j];
        if prov.ends_with(suffix.as_str()) {
            already += 1;
            continue;
        }
        if krhom::sindec::route_header(prov) != krhom::sindec::route_header(&chain) {
            other += 1;
            eprintln!("skip (other provider stamp {:?}): {}", prov, p.display());
            continue;
        }
        let nh = format!("{} prov={}{}{}", &h[..i], prov, suffix, &h[j..]);
        if !dry {
            let tmp = p.with_extension("restamp.tmp");
            std::fs::write(&tmp, format!("{}\n{}", nh, body)).unwrap();
            std::fs::rename(&tmp, &p).unwrap();
        }
        done += 1;
    }
    println!("b4wstamp {}: restamped {}, already m4c {}, other {} ({})", dir.display(), done, already, other, if dry { "dry run" } else { "written" });
}
