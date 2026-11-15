//! C2-4: per braid word (KNOTS lines "name;labels;word;N"), the library jobs of a pruned scan
//! (sindec::knot_jobs: upper bound and K_0 lower bound) missing from a status file (lines "JOB ok",
//! JOB = "pair K" | "closure TABLE"); idempotent types (M5i) marked.
//! usage: m5i_cover LAM N KNOTS STATUS
use krhom::combinat::schur::Schur;
use krhom::library::sindec::{canonical_ladder, is_indecomposable_of, knot_jobs, pruned_plan};
use std::collections::BTreeSet;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let lam: Vec<u8> = a[1].split(',').map(|x| x.parse().unwrap()).collect();
    let n: u8 = a[2].parse().unwrap();
    let ok: BTreeSet<String> = std::fs::read_to_string(&a[4]).unwrap().lines().filter(|l| l.ends_with(" ok")).map(|l| l[..l.len() - 3].to_string()).collect();
    let mut s = Schur::new(&lam);
    let plan = pruned_plan(&mut s, Some(n));
    for line in std::fs::read_to_string(&a[3]).unwrap().lines() {
        let f: Vec<&str> = line.split(';').collect();
        if f.len() < 3 {
            continue;
        }
        let word: Vec<i32> = f[2].split_whitespace().map(|x| x.parse().unwrap()).collect();
        let kj = match knot_jobs(&mut s, &plan, &word) {
            Ok(k) => k,
            Err(e) => {
                println!("{}: {}", f[0], e);
                continue;
            }
        };
        let mut miss = vec![];
        let mut miss_lb = 0;
        let mut idem = BTreeSet::new();
        for &k in &kj.pairs {
            let pp = &plan.pairs[k];
            let w = canonical_ladder(&mut s, pp.p);
            if !is_indecomposable_of(&mut s, &w, pp.p) {
                idem.insert(s.table(pp.p).fmt());
            }
            for t in &pp.kept {
                let wq = canonical_ladder(&mut s, t.q);
                if !is_indecomposable_of(&mut s, &wq, t.q) {
                    idem.insert(s.table(t.q).fmt());
                }
            }
            let key = format!("pair {}", k);
            if !ok.contains(&key) {
                miss_lb += kj.pairs_lb.contains(&k) as usize;
                miss.push(key);
            }
        }
        for &q in &kj.closure {
            let key = format!("closure {}", s.table(q).fmt());
            let w = canonical_ladder(&mut s, q);
            if w.is_empty() {
                continue;
            }
            if !is_indecomposable_of(&mut s, &w, q) {
                idem.insert(s.table(q).fmt());
            }
            if !ok.contains(&key) {
                miss_lb += kj.closure_lb.contains(&q) as usize;
                miss.push(key);
            }
        }
        println!("{} [{}]: jobs {}+{} (LB {}+{}), missing {} (LB {}), idempotent types visited {} | {}", f[0], f[2], kj.pairs.len(), kj.closure.len(), kj.pairs_lb.len(), kj.closure_lb.len(), miss.len(), miss_lb, idem.len(), miss.join(", "));
    }
}
