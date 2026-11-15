// Time End^0 of uncolored Bott–Samelson bimodules BS(s_i1 ... s_ik) on 4 strands.
use krhom::algebra::field::F31;
use krhom::algebra::mpoly::MonoCache;
use krhom::soergel::sbim::{hom_space, SBim};
use krhom::combinat::web::Rung;

fn bs(n: usize, strands: usize, word: &[usize]) -> SBim<F31> {
    let mut b = SBim::<F31>::identity(n, &vec![1; strands]);
    for &i in word {
        // B_{s_i}: (1,1) -> E: (2,0) -> F: (1,1)
        b = b.add_rung(&Rung::e(i, 1)).add_rung(&Rung::f(i, 1));
    }
    b
}

fn main() {
    let words: Vec<Vec<usize>> = vec![vec![0], vec![0, 1], vec![0, 1, 0], vec![1, 0, 2, 1], vec![0, 1, 0, 2, 1], vec![0, 1, 0, 2, 1, 0]];
    for w in words {
        let b = bs(4, 4, &w);
        let mut mc = MonoCache::new(b.ring.clone());
        let t = std::time::Instant::now();
        let h = hom_space(&b, &b, 0, &mut mc);
        println!("word {:?}: fsize {} End^0 dim {} time {:.2}s", w, b.fsize(), h.len(), t.elapsed().as_secs_f64());
    }
}
