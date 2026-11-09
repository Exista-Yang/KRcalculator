// Local Rickard differentials: dimension of degree-1 Hom between consecutive terms, and d^2 = 0.
use krhom::crossing::{crossing_cx, Gen};
use krhom::field::F31;
use krhom::mpoly::MonoCache;
use krhom::sbim::{hom_space, SBim};

fn main() {
    for n in 2..=5usize {
        for c in 1..=n {
            for d in 1..=n {
                for pos in [true, false] {
                    let g = Gen { i: 0, positive: pos };
                    let cx = crossing_cx(g, c, d, n);
                    let id = SBim::<F31>::identity(n, &[c, d]);
                    let webs: Vec<SBim<F31>> = cx.terms.iter().map(|t| id.add_rungs(&t.rungs)).collect();
                    let mut mc = MonoCache::new(id.ring.clone());
                    let mut dims = vec![];
                    let mut maps = vec![];
                    for j in 0..webs.len().saturating_sub(1) {
                        // positive: W_j -> W_{j+1}; negative: W_{j+1} -> W_j
                        let (s, t) = if pos { (&webs[j], &webs[j + 1]) } else { (&webs[j + 1], &webs[j]) };
                        let u = 1 + s.shift - t.shift;
                        let h = hom_space(s, t, u, &mut mc);
                        dims.push(h.len());
                        maps.push(h.into_iter().next());
                    }
                    // d^2
                    let mut d2ok = true;
                    for j in 0..maps.len().saturating_sub(1) {
                        if let (Some(a), Some(b)) = (&maps[j], &maps[j + 1]) {
                            let comp = if pos { b.matmul(a) } else { a.matmul(b) };
                            if !comp.is_zero() {
                                d2ok = false;
                            }
                        }
                    }
                    if dims.iter().any(|&x| x != 1) || !d2ok {
                        println!("N={} c={} d={} pos={} terms={} dims={:?} d2ok={}", n, c, d, pos, webs.len(), dims, d2ok);
                    }
                }
            }
        }
        println!("N={} done", n);
    }
}
