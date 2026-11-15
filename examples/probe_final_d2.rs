use krhom::input::crossing::Gen;
use krhom::algebra::field::F31;
use krhom::colored::xengine::*;
fn gens(w: &[i32]) -> Vec<Gen> { w.iter().map(|&g| Gen { i: (g.unsigned_abs() - 1) as usize, positive: g > 0 }).collect() }
fn main() {
    let t = std::time::Instant::now();
    let mut found = 0;
    for labels in [vec![3usize, 3], vec![2,2]] {
    for len in 3..=8usize {
        for mask in 0..(1u32 << len) {
            let w: Vec<i32> = (0..len).map(|i| if mask >> i & 1 == 1 { -1 } else { 1 }).collect();
            if w[0] < 0 { continue; }
            for p in [3usize, 4] {
                if p < labels[0] { continue; }
                let t1 = std::time::Instant::now();
                let Ok(sc) = xscan::<F31>(&labels, &gens(&w), &XOpts { prune_n: Some(p), ..Default::default() }) else { continue };
                if sc.final_d2.is_empty() { continue; }
                let mut non = vec![];
                for (i, (x, y, f)) in sc.final_d2.iter().enumerate() {
                    for n in p+1..=p+2 { if !closure_kills_map(&sc.data, *x, *y, f, n, None) { non.push((i, n)); break; } }
                }
                println!("{:?} {:?} p={} objs={} final_d2={} nonclosing(i,N)={:?} scan {:.1}s total {:.1}", labels, w, p, sc.data.objs.len(), sc.final_d2.len(), non, t1.elapsed().as_secs_f64(), t.elapsed().as_secs_f64());
                if !non.is_empty() { found += 1; if found > 3 { return; } }
            }
            if t.elapsed().as_secs() > 500 { return; }
        }
    }
    }
}
