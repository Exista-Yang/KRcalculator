//! C2-4 (M5i validation): non-canonical *decomposable* ladders for engine tests. BFS over rung
//! words from 1_λ (states deduplicated by character), at the block λ (closable types) or any
//! block (`any`): a word W is a forced-ladder candidate for p when b_p occurs in char W with
//! coefficient exactly 1, char W has coefficient 1 at the table p (one flow at rep(p)), and
//! expand(char W) has more than one term. Prints the shortest per type, with flows.
//! usage: m5i_force_search LAM DEPTH CAP [any] [maxlabel]
use krhom::schur::{flows, Elt, Lp, Schur, Tid};
use krhom::web::{Dir, Rung};
use std::collections::{BTreeMap, HashSet, VecDeque};
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let lam: Vec<u8> = a[1].split(',').map(|x| x.parse().unwrap()).collect();
    let depth: usize = a[2].parse().unwrap();
    let cap: usize = a[3].parse().unwrap();
    let any = a.get(4).is_some_and(|x| x == "any");
    let maxl: u8 = a.get(5).map(|x| x.parse().unwrap()).unwrap_or(99);
    let mut s = Schur::new(&lam);
    let mut found: BTreeMap<Tid, (Vec<Rung>, String, u128)> = BTreeMap::new();
    let mut seen: HashSet<String> = HashSet::new();
    let mut q: VecDeque<(Vec<Rung>, Elt)> = VecDeque::new();
    let u = s.unit();
    seen.insert(format!("{:?}", u));
    q.push_back((vec![], u));
    let nn = lam.len();
    let mut states = 0usize;
    while let Some((w, x)) = q.pop_front() {
        states += 1;
        if states > cap {
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
                    if y.t.is_empty() || y.mu.iter().any(|&l| l > maxl) {
                        continue;
                    }
                    if !seen.insert(format!("{:?}", y)) {
                        continue;
                    }
                    let mut w2 = w.clone();
                    w2.push(rg);
                    if any || y.mu == lam {
                        let ex = s.expand(&y);
                        if ex.len() > 1 {
                            for (p, m) in &ex {
                                if *m == Lp::one() && y.coef(*p) == Lp::one() && !found.contains_key(p) {
                                    let parts: Vec<String> = ex.iter().map(|(t, c)| format!("{}{{{}}}", s.table(*t).fmt(), c.fmt_ec())).collect();
                                    let fl = flows(&lam, &w2);
                                    println!("FORCE {} = {} | len {} flows {} | char {}", s.table(*p).fmt(), krhom::smapio::rungs_str(&w2), w2.len(), fl, parts.join(" + "));
                                    found.insert(*p, (w2.clone(), parts.join(" + "), fl));
                                }
                            }
                        }
                    }
                    if q.len() < cap {
                        q.push_back((w2, y));
                    }
                }
            }
        }
    }
    println!("states {}, types with a forced decomposable ladder: {}", states, found.len());
}
