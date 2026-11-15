// PD code -> braid word (Vogel's algorithm + simplification).
// usage: pd2braid [--no-simplify] ["[[1,4,2,5],[3,6,4,1],[5,2,6,3]]"]   (reads stdin if no PD given)
// PD convention: KnotTheory / Knot Atlas (X[i,j,k,l] counterclockwise from the incoming under-strand).
// Output: number of strands, the word as signed 1-based generators (k: σ_k, -k: σ_k^{-1},
// σ_k positive = left strand over), and the link component of each bottom strand position.
use krhom::input::pd::{parse_pd, pd_to_braid, pd_to_braid_vogel};
use std::io::Read;

fn main() {
    let mut args: Vec<String> = std::env::args().skip(1).collect();
    let no_simplify = args.iter().any(|a| a == "--no-simplify");
    args.retain(|a| a != "--no-simplify");
    let input = if args.is_empty() {
        let mut s = String::new();
        std::io::stdin().read_to_string(&mut s).expect("reading stdin");
        s
    } else {
        args.join(" ")
    };
    let pd = match parse_pd(&input) {
        Ok(pd) => pd,
        Err(e) => {
            eprintln!("error: {e}");
            std::process::exit(1);
        }
    };
    let b = if no_simplify { pd_to_braid_vogel(&pd) } else { pd_to_braid(&pd) };
    let word: Vec<String> = b.signed().iter().map(|k| k.to_string()).collect();
    let comps: Vec<String> = b.components.iter().map(|c| c.to_string()).collect();
    println!("strands: {}", b.strands);
    println!("word: {}", word.join(" "));
    println!("components: {}", comps.join(" "));
}
