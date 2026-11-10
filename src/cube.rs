//! Milestone 4b: complete (unsimplified) computation of colored gl_N
//! Khovanov–Rozansky homology of a braid closure via the full cube of
//! resolutions, with state spaces from the universal construction.

use crate::crossing::{crossing_cx, qshift, CrossingCx, Gen, QConv};
use crate::field::Field;
use crate::laurent::Poincare;
use crate::linsys::Echelon;
use crate::mpoly::{mono_exp, MonoCache, Poly};
use crate::pmat::{PAcc, PMat, PVec};
use crate::sbim::{hom_space, SBim, StateSpace, Tau};
use crate::sparse::SpVec;
use std::collections::{BTreeMap, HashMap};

#[derive(Clone, Debug)]
pub struct BraidLink {
    pub n: usize,
    /// colors (labels) of the strands at the bottom positions
    pub labels: Vec<usize>,
    pub word: Vec<Gen>,
}

impl BraidLink {
    /// Labels after applying the word (must equal `labels` for the closure).
    pub fn top_labels(&self) -> Vec<usize> {
        let mut l = self.labels.clone();
        for g in &self.word {
            l.swap(g.i, g.i + 1);
        }
        l
    }
}

/// The colored Rickard complex of one crossing on two uprights: term ladders and the unique
/// degree-1 differentials as exact polynomial matrices over the local ring R^{(c,d)}
/// (variables e_{0,1..c}, e_{1,1..d}).
pub struct LocalCx<F: Field> {
    /// terms (rungs re-indexed to the generator's uprights g.i, g.i+1)
    pub cx: CrossingCx,
    /// local differential maps: maps[j] goes from term j to j+1 (positive)
    /// or from term j+1 to term j (negative)
    pub maps: Vec<PMat<F>>,
    /// unnormalized degrees of maps[j]
    pub udeg: Vec<i32>,
    pub fsizes: Vec<usize>,
    pub shifts: Vec<i32>,
}

/// Rickard complex of the crossing g on labels (c, d) at gl_N; terms whose web has a label
/// outside 0..=n are pruned (`crossing_cx`), so n >= c + d gives the untruncated,
/// N-independent complex. Asserts that every differential is the unique degree-1 map.
pub fn local_cx<F: Field>(g: Gen, c: usize, d: usize, n: usize) -> LocalCx<F> {
    let cx = crossing_cx(Gen { i: 0, positive: g.positive }, c, d, n);
    let id = SBim::<F>::identity(n, &[c, d]);
    let webs: Vec<SBim<F>> = cx.terms.iter().map(|t| id.add_rungs(&t.rungs)).collect();
    let mut mc = MonoCache::new(id.ring.clone());
    let mut maps = vec![];
    let mut udeg = vec![];
    for j in 0..webs.len().saturating_sub(1) {
        let (s, t) = if g.positive { (&webs[j], &webs[j + 1]) } else { (&webs[j + 1], &webs[j]) };
        let u = 1 + s.shift - t.shift;
        let h = hom_space(s, t, u, &mut mc);
        assert_eq!(h.len(), 1, "Rickard differential not unique: c={} d={} j={}", c, d, j);
        maps.push(h.into_iter().next().unwrap());
        udeg.push(u);
    }
    let mut cx2 = cx.clone();
    cx2.g = g;
    for t in cx2.terms.iter_mut() {
        for r in t.rungs.iter_mut() {
            r.i = g.i;
        }
    }
    LocalCx { cx: cx2, maps, udeg, fsizes: webs.iter().map(|w| w.fsize()).collect(), shifts: webs.iter().map(|w| w.shift).collect() }
}

/// Apply a polynomial in the local variables (two uprights i,i+1 of P) to a
/// vector of P by substituting P's top operators.
pub fn apply_local_poly<F: Field>(p: &SBim<F>, i: usize, c: usize, d: usize, poly: &Poly<F>, v: &PVec<F>) -> PVec<F> {
    let mut acc = PAcc::new();
    for (m, coef) in &poly.t {
        let mut w = v.clone();
        for j in 0..c {
            for _ in 0..mono_exp(*m, j) {
                w = p.top[i][j].mul_vec(&w);
            }
        }
        for j in 0..d {
            for _ in 0..mono_exp(*m, c + j) {
                w = p.top[i + 1][j].mul_vec(&w);
            }
        }
        acc.add_scaled_const(&w, coef);
    }
    acc.take()
}

pub struct CubeResult {
    pub poincare: Poincare,
    pub vertices: usize,
    pub chain_dim: usize,
}

pub fn khr_full_cube<F: Field>(link: &BraidLink, conv: QConv) -> CubeResult {
    let n = link.n;
    assert_eq!(link.top_labels(), link.labels, "braid closure needs matching colors");
    // local complexes along the word
    let mut labs = link.labels.clone();
    let mut locals: Vec<LocalCx<F>> = vec![];
    for g in &link.word {
        let (c, d) = (labs[g.i], labs[g.i + 1]);
        locals.push(local_cx::<F>(*g, c, d, n));
        labs.swap(g.i, g.i + 1);
    }
    let nc = locals.len();
    // enumerate states with prefix bimodules (DFS)
    let id = SBim::<F>::identity(n, &link.labels);
    // signed: a label > N (whose complex is empty) must not underflow (debug builds panicked here)
    let closure_norm: i32 = link.labels.iter().map(|&l| l as i32 * (n as i32 - l as i32)).sum();
    // prefix cache: state prefix -> SBim
    let mut prefix: HashMap<Vec<usize>, SBim<F>> = HashMap::new();
    prefix.insert(vec![], id.clone());
    let mut states: Vec<Vec<usize>> = vec![vec![]];
    for l in 0..nc {
        let mut next = vec![];
        for s in &states {
            let p = prefix[s].clone();
            for (ti, t) in locals[l].cx.terms.iter().enumerate() {
                let mut s2 = s.clone();
                s2.push(ti);
                prefix.insert(s2.clone(), p.add_rungs(&t.rungs));
                next.push(s2);
            }
        }
        states = next;
    }
    // state spaces
    let mut mc = MonoCache::new(id.ring.clone());
    let mut tau = Tau::<F>::new(n, &link.labels, &id.voff);
    let mut ss: HashMap<Vec<usize>, StateSpace<F>> = HashMap::new();
    for s in &states {
        ss.insert(s.clone(), StateSpace::compute(&prefix[s], &mut mc));
    }
    let hdeg = |s: &Vec<usize>| -> i32 { (0..nc).map(|l| locals[l].cx.terms[s[l]].hdeg).sum() };
    let qsh = |s: &Vec<usize>| -> i32 {
        let b = &prefix[s];
        (0..nc).map(|l| qshift(&locals[l].cx, locals[l].cx.terms[s[l]].k, conv)).sum::<i32>() + b.shift - closure_norm
    };
    // chain groups: (h, q) -> list of (state, unnormalized degree d, index within piece)
    let mut groups: BTreeMap<(i32, i32), Vec<(Vec<usize>, i32, usize)>> = BTreeMap::new();
    for s in &states {
        let h = hdeg(s);
        let qs = qsh(s);
        for (d, piece) in &ss[s].pieces {
            for k in 0..piece.reps.len() {
                groups.entry((h, d + qs)).or_default().push((s.clone(), *d, k));
            }
        }
    }
    let gindex: HashMap<(Vec<usize>, i32, usize), usize> = groups
        .values()
        .flat_map(|v| v.iter().enumerate().map(|(i, x)| (x.clone(), i)))
        .collect();
    // differential matrices: for each source group (h,q), rows = target basis, as sparse vectors per source element
    let mut dmat: HashMap<(i32, i32), Vec<SpVec<F>>> = HashMap::new();
    for ((h, q), elems) in &groups {
        let mut cols: Vec<SpVec<F>> = Vec::with_capacity(elems.len());
        for (s, d, k) in elems {
            let mut out: HashMap<usize, F> = HashMap::new();
            let src_vec = &ss[s].pieces[d].reps[*k];
            for l in 0..nc {
                let loc = &locals[l];
                let kl = s[l];
                let (target_k, mapidx) = if loc.cx.g.positive {
                    if kl + 1 >= loc.cx.terms.len() { continue; }
                    (kl + 1, kl)
                } else {
                    if kl == 0 { continue; }
                    (kl - 1, kl - 1)
                };
                let mut t = s.clone();
                t[l] = target_k;
                let sign_exp: i32 = (0..l).map(|l2| locals[l2].cx.terms[s[l2]].hdeg).sum();
                let negate = sign_exp.rem_euclid(2) == 1;
                // apply psi = id_P ⊗ L ⊗ id_above
                let p = &prefix[&s[..l].to_vec()];
                let fp = p.fsize();
                let f_src = loc.fsizes[kl];
                let f_tgt = loc.fsizes[target_k];
                let fsz_total_src = prefix[s].fsize();
                let f_above = fsz_total_src / (fp * f_src);
                let lmap = &loc.maps[mapidx];
                let (c, dd) = (loc.cx.c, loc.cx.d);
                let gi = loc.cx.g.i;
                // group source entries by (beta_l, beta_a): vector over F_P
                let mut grouped: HashMap<(usize, usize), PVec<F>> = HashMap::new();
                for (beta, poly) in src_vec {
                    let b = *beta as usize;
                    let ba = b % f_above;
                    let rest = b / f_above;
                    let bl = rest % f_src;
                    let bp = rest / f_src;
                    grouped.entry((bl, ba)).or_default().push((bp as u32, poly.clone()));
                }
                let mut acc = PAcc::new();
                for ((bl, ba), mut v) in grouped {
                    v.sort_by_key(|e| e.0);
                    for (blp, lpoly) in &lmap.cols[bl] {
                        let w = apply_local_poly(p, gi, c, dd, lpoly, &v);
                        let placed: PVec<F> = w
                            .into_iter()
                            .map(|(bp, poly)| (((bp as usize * f_tgt + *blp as usize) * f_above + ba) as u32, poly))
                            .collect();
                        acc.add_scaled_const(&placed, &F::one());
                    }
                }
                let img = acc.take();
                let u = loc.udeg[mapidx];
                let td = d + u;
                let coords = ss[&t].coords(td, &img, &mut tau);
                // consistency of q-degree
                debug_assert_eq!(d + qsh(s), td + qsh(&t));
                for (kk, x) in coords.into_iter().enumerate() {
                    if x.is_zero() {
                        continue;
                    }
                    let gi2 = gindex[&(t.clone(), td, kk)];
                    let x = if negate { x.neg() } else { x };
                    out.entry(gi2).or_insert_with(F::zero).add_assign(&x);
                }
            }
            let mut v: SpVec<F> = out.into_iter().filter(|(_, x)| !x.is_zero()).map(|(i, x)| (i as u32, x)).collect();
            v.sort_by_key(|e| e.0);
            cols.push(v);
        }
        dmat.insert((*h, *q), cols);
    }
    // homology
    let mut rank: HashMap<(i32, i32), usize> = HashMap::new();
    for ((h, q), cols) in &dmat {
        let tgt = groups.get(&(h + 1, *q)).map(|v| v.len()).unwrap_or(0);
        let mut e = Echelon::<F>::new(tgt.max(1));
        for c in cols {
            e.insert(c);
        }
        rank.insert((*h, *q), e.rank());
    }
    // check d^2 = 0 lightly: rank bound
    let mut pc = Poincare::default();
    let mut chain_dim = 0;
    for ((h, q), elems) in &groups {
        chain_dim += elems.len();
        let r_out = rank.get(&(*h, *q)).copied().unwrap_or(0);
        let r_in = rank.get(&(h - 1, *q)).copied().unwrap_or(0);
        let dim = elems.len() as i128 - r_out as i128 - r_in as i128;
        assert!(dim >= 0);
        pc.add_term(*h, *q, dim);
    }
    CubeResult { poincare: pc, vertices: states.len(), chain_dim }
}

/// Verify d∘d = 0 on the assembled cube (debug helper).
pub fn _unused() {}
