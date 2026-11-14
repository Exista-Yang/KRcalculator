//! C2-2: search for an indecomposable ladder of the types whose minimal-flow canonical ladder is
//! decomposable (M5i types): BFS over rung words from 1_λ, states deduplicated by character, up to
//! a depth and state cap. usage: ladder_search 2,2,2,2 N DEPTH CAP
use krhom::schur::{Elt, Schur, Tid};
use krhom::sindec::{idem_needed, is_indecomposable_of};
use krhom::web::{Dir, Rung};
use std::collections::{HashMap, HashSet, VecDeque};
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let lam: Vec<u8> = a[1].split(',').map(|x| x.parse().unwrap()).collect();
    let n: u8 = a[2].parse().unwrap();
    let depth: usize = a[3].parse().unwrap();
    let cap: usize = a[4].parse().unwrap();
    let mut s = Schur::new(&lam);
    let cases = idem_needed(&mut s, Some(n));
    // targets: Tid by table string
    let mut targets: HashMap<Tid, String> = HashMap::new();
    for c in &cases {
        let order = s.canonical_basis(&c.mu.clone()).order.clone();
        // c.mu is the block (column sums)
        for q in order {
            if s.table(q).fmt() == c.table {
                targets.insert(q, c.table.clone());
            }
        }
    }
    println!("{} targets", targets.len());
    let mut found: HashMap<Tid, Vec<Rung>> = HashMap::new();
    let mut seen: HashSet<String> = HashSet::new();
    let mut q: VecDeque<(Vec<Rung>, Elt)> = VecDeque::new();
    let u = s.unit();
    seen.insert(format!("{:?}", u));
    q.push_back((vec![], u));
    let nn = lam.len();
    let mut states = 0usize;
    while let Some((w, x)) = q.pop_front() {
        states += 1;
        if states > cap || found.len() == targets.len() {
            break;
        }
        if w.len() >= depth {
            continue;
        }
        let cur = x.mu.clone();
        for i in 0..nn - 1 {
            for dir in [Dir::E, Dir::F] {
                let giver = if dir == Dir::E { cur[i + 1] } else { cur[i] } as usize;
                for r in 1..=giver {
                    let rg = Rung { i, dir, r };
                    let y = s.rung(&x, &rg);
                    if y.t.is_empty() {
                        continue;
                    }
                    if !seen.insert(format!("{:?}", y)) {
                        continue;
                    }
                    let mut w2 = w.clone();
                    w2.push(rg);
                    // candidate: the character contains a target with coefficient 1 at its top
                    for (tq, _) in y.t.iter() {
                        if targets.contains_key(tq) && !found.contains_key(tq) && is_indecomposable_of(&mut s, &w2, *tq) {
                            println!("FOUND {} : {} (len {})", targets[tq], krhom::smapio::rungs_str(&w2), w2.len());
                            found.insert(*tq, w2.clone());
                        }
                    }
                    if q.len() < cap { q.push_back((w2, y)); }
                }
            }
        }
    }
    println!("states {}, found {} of {}", states, found.len(), targets.len());
    for (tq, t) in &targets {
        if !found.contains_key(tq) {
            println!("NOT FOUND {}", t);
        }
    }
}
