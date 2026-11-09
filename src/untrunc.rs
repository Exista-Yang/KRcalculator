//! EXPERIMENTAL (superseded by `sbim`): untruncated singular Bott–Samelson bimodules (edge rings without
//! facet truncations), cut off above a q-degree `cap`, together with the
//! operators generating the facet-truncation ideal J. Used to test the model
//!   F_N(closure of W) = Hom_SSBim(R, B_W) / J·Hom_SSBim(R, B_W).

use crate::algebra::{tensor_id, truncation_relations, Alg};
use crate::field::Field;
use crate::hom::joint_kernel;
use crate::laurent::LPoly;
use crate::sparse::{Acc, SpMat, SpVec};
use crate::web::{Ladder, Rung};

pub struct UAlg<F: Field> {
    pub alg: Alg<F>,
    pub cap: i32,
    /// operators generating J (facet truncation relations of all edges)
    pub rels: Vec<SpMat<F>>,
    /// for each basis index: (index of R_λ monomial, composite fibre index)
    pub pb: Vec<(u32, u32)>,
    /// total fibre size (mixed radix product)
    pub ftot: usize,
    /// exponent vectors of the R_λ monomials (concatenated e-variables of all uprights)
    pub rexp: Vec<Vec<u16>>,
}

/// Restrict an algebra to basis elements of degree <= cap (quotient by the
/// ideal of high-degree elements). Returns kept indices.
fn cap_filter<F: Field>(deg: &[i32], cap: i32) -> (Vec<usize>, Vec<u32>) {
    let keep: Vec<usize> = (0..deg.len()).filter(|&i| deg[i] <= cap).collect();
    let mut newidx = vec![u32::MAX; deg.len()];
    for (k, &i) in keep.iter().enumerate() {
        newidx[i] = k as u32;
    }
    (keep, newidx)
}

fn restrict<F: Field>(m: &SpMat<F>, keep: &[usize], newidx: &[u32]) -> SpMat<F> {
    let cols = keep
        .iter()
        .map(|&i| m.cols[i].iter().filter(|(r, _)| newidx[*r as usize] != u32::MAX).map(|(r, x)| (newidx[*r as usize], x.clone())).collect())
        .collect();
    SpMat { nrows: keep.len(), ncols: keep.len(), cols }
}

impl<F: Field> UAlg<F> {
    fn apply_cap(&mut self) {
        let (keep, newidx) = cap_filter::<F>(&self.alg.deg, self.cap);
        if keep.len() == self.alg.dim {
            return;
        }
        let a = &mut self.alg;
        a.deg = keep.iter().map(|&i| a.deg[i]).collect();
        a.dim = keep.len();
        let r = |m: &SpMat<F>| restrict(m, &keep, &newidx);
        a.top = a.top.iter().map(|v| v.iter().map(&r).collect()).collect();
        a.bot = a.bot.iter().map(|v| v.iter().map(&r).collect()).collect();
        a.eq = a.eq.iter().map(&r).collect();
        self.rels = self.rels.iter().map(&r).collect();
        self.pb = keep.iter().map(|&i| self.pb[i]).collect();
    }

    /// Polynomial ring ⊗_i Sym_{λ_i} capped at degree cap.
    pub fn identity(n: usize, labels: &[usize], cap: i32) -> Self {
        let mut u = UAlg { alg: Alg::base(n), cap, rels: vec![], pb: vec![(0, 0)], ftot: 1, rexp: vec![vec![]] };
        for &a in labels {
            u = u.add_upright(a);
        }
        u
    }

    fn add_upright(self, a: usize) -> Self {
        // Sym_a = k[e_1..e_a]: build monomial basis of weight <= cap/2 for this factor
        let maxw = (self.cap / 2).max(0) as usize;
        let mut monos: Vec<Vec<usize>> = vec![];
        fn rec(i: usize, a: usize, left: usize, cur: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
            if i == a {
                out.push(cur.clone());
                return;
            }
            let w = i + 1;
            let mut k = 0;
            while k * w <= left {
                cur.push(k);
                rec(i + 1, a, left - k * w, cur, out);
                cur.pop();
                k += 1;
            }
        }
        rec(0, a, maxw, &mut vec![], &mut monos);
        let fsize = monos.len();
        let idx: std::collections::HashMap<Vec<usize>, usize> = monos.iter().cloned().enumerate().map(|(i, m)| (m, i)).collect();
        let fdeg: Vec<i32> = monos.iter().map(|m| 2 * m.iter().enumerate().map(|(i, &k)| (i + 1) * k).sum::<usize>() as i32).collect();
        let p = &self.alg;
        let dim = p.dim * fsize;
        let mut deg = Vec::with_capacity(dim);
        for pi in 0..p.dim {
            for b in 0..fsize {
                deg.push(p.deg[pi] + fdeg[b]);
            }
        }
        let mut eops = vec![];
        for j in 0..a {
            let mut cols: Vec<SpVec<F>> = Vec::with_capacity(dim);
            for pi in 0..p.dim {
                for b in 0..fsize {
                    let mut m = monos[b].clone();
                    m[j] += 1;
                    match idx.get(&m) {
                        Some(&t) => cols.push(vec![((pi * fsize + t) as u32, F::one())]),
                        None => cols.push(vec![]),
                    }
                }
            }
            eops.push(SpMat { nrows: dim, ncols: dim, cols });
        }
        let tr = |m: &SpMat<F>| tensor_id(m, fsize);
        let mut top: Vec<Vec<SpMat<F>>> = p.top.iter().map(|v| v.iter().map(&tr).collect()).collect();
        let mut bot: Vec<Vec<SpMat<F>>> = p.bot.iter().map(|v| v.iter().map(&tr).collect()).collect();
        top.push(eops.clone());
        bot.push(eops.clone());
        let eq: Vec<SpMat<F>> = p.eq.iter().map(&tr).collect();
        let mut rels: Vec<SpMat<F>> = self.rels.iter().map(&tr).collect();
        rels.extend(truncation_relations(&eops, &eq, p.n, dim));
        let mut tl = p.top_labels.clone();
        tl.push(a);
        let mut bl = p.bot_labels.clone();
        bl.push(a);
        let alg = Alg { n: p.n, dim, deg, top_labels: tl, bot_labels: bl, top, bot, eq, shift: p.shift };
        let mut rexp = Vec::with_capacity(dim);
        for pi in 0..p.dim {
            for b in 0..fsize {
                let mut e = self.rexp[self.pb[pi].0 as usize].clone();
                e.extend(monos[b].iter().map(|&k| k as u16));
                rexp.push(e);
            }
        }
        let pb: Vec<(u32, u32)> = (0..dim).map(|i| (i as u32, 0)).collect();
        let mut u = UAlg { alg, cap: self.cap, rels, pb, ftot: 1, rexp };
        u.apply_cap();
        // renumber R_λ monomials to match kept indices
        u.rexp = u.pb.iter().map(|&(pi, _)| u.rexp[pi as usize].clone()).collect();
        u.pb = (0..u.alg.dim).map(|i| (i as u32, 0)).collect();
        u
    }

    /// Free (untruncated) rung extension; adds the merged edge's truncation
    /// relations to J.
    pub fn add_rung(&self, rung: &Rung) -> Self {
        let p = &self.alg;
        let n = p.n;
        let (rec, giv) = rung.ends();
        let a = p.top_labels[rec];
        let b = p.top_labels[giv];
        let r = rung.r;
        // reuse Alg::add_rung on a copy whose eq is zero and which never truncates:
        // we emulate by calling the free part directly.
        let ext = p.add_rung_free(rung);
        let fsize = ext.1;
        let mut alg = ext.0;
        let tr = |m: &SpMat<F>| tensor_id(m, fsize);
        let mut rels: Vec<SpMat<F>> = self.rels.iter().map(&tr).collect();
        if a > 0 && r > 0 {
            rels.extend(truncation_relations(&alg.top[rec], &alg.eq, n, alg.dim));
        }
        let _ = b;
        alg.shift = p.shift - (r * (b - r)) as i32;
        let mut pb = Vec::with_capacity(alg.dim);
        for &(pi, bi) in &self.pb {
            for b in 0..fsize {
                pb.push((pi, bi * fsize as u32 + b as u32));
            }
        }
        let mut u = UAlg { alg, cap: self.cap, rels, pb, ftot: self.ftot * fsize, rexp: self.rexp.clone() };
        u.apply_cap();
        u
    }

    pub fn ladder(lad: &Ladder, cap: i32) -> Self {
        let mut u = Self::identity(lad.n, &lad.bottom, cap);
        for rg in &lad.rungs {
            u = u.add_rung(rg);
        }
        u
    }

    /// Graded dim of Hom(R,B)/J·Hom(R,B), reliable in degrees <= reliable_max.
    pub fn candidate_state_space(&self) -> (LPoly, i32) {
        let a = &self.alg;
        let maxlab = a.top_labels.iter().chain(a.bot_labels.iter()).copied().max().unwrap_or(0) as i32;
        let reliable_inv = self.cap - 2 * maxlab;
        let ops = crate::hom::closure_ops(a);
        let inv = joint_kernel(&a.deg, &ops);
        // all invariant vectors with degree <= reliable_inv
        let mut invvecs: Vec<SpVec<F>> = vec![];
        let mut by_deg: std::collections::BTreeMap<i32, Vec<SpVec<F>>> = std::collections::BTreeMap::new();
        for (d, vs) in &inv {
            if *d <= reliable_inv {
                for v in vs {
                    invvecs.push(v.clone());
                }
                by_deg.insert(*d, vs.clone());
            }
        }
        // J·Hom: images rho * v
        let mut acc = Acc::new(a.dim.max(1));
        let mut jgens: std::collections::BTreeMap<i32, Vec<SpVec<F>>> = std::collections::BTreeMap::new();
        for rho in &self.rels {
            for v in &invvecs {
                let w = rho.mul_vec(v, &mut acc);
                if !w.is_empty() {
                    let d = a.deg[w[0].0 as usize];
                    jgens.entry(d).or_default().push(w);
                }
            }
        }
        let mut out = LPoly::zero();
        for (d, vs) in &by_deg {
            // dim of span(vs) / span(jgens[d]) — jgens lie inside the invariants
            let js = jgens.get(d).cloned().unwrap_or_default();
            let rank_j = span_rank(&js, &a.deg, *d);
            let k = vs.len() as i128 - rank_j as i128;
            out.add_term(*d, k);
        }
        (out, reliable_inv)
    }
}

fn span_rank<F: Field>(vs: &[SpVec<F>], deg: &[i32], d: i32) -> usize {
    if vs.is_empty() {
        return 0;
    }
    let coords: Vec<usize> = (0..deg.len()).filter(|&i| deg[i] == d).collect();
    let loc: std::collections::HashMap<usize, usize> = coords.iter().enumerate().map(|(k, &i)| (i, k)).collect();
    let mut rows: Vec<Vec<F>> = vs
        .iter()
        .map(|v| {
            let mut r = vec![F::zero(); coords.len()];
            for (i, x) in v {
                r[loc[&(*i as usize)]] = x.clone();
            }
            r
        })
        .collect();
    crate::dense::rref(&mut rows, coords.len()).len()
}

/// Universal-construction state space of the closure of W:
/// Hom_SSBim(R,B_W) modulo the radical of (f, g) -> tau_N(g∘f), g in Hom_SSBim(B_W,R).
/// Returns graded dimension in unshifted degrees of Hom(R,B_W).
pub fn universal_state_space<F: Field>(u: &UAlg<F>) -> LPoly {
    use std::collections::HashMap;
    let a = &u.alg;
    let n = a.n;
    let lam = a.bot_labels.clone();
    assert_eq!(lam, a.top_labels);
    let top: i32 = 2 * lam.iter().map(|&l| (l * (n - l)) as i32).sum::<i32>();
    // R_λ polynomial arithmetic on exponent vectors; variable weights
    let wts: Vec<i32> = lam.iter().flat_map(|&l| 1..=l as i32).collect();
    let mdeg = |e: &Vec<u16>| -> i32 { 2 * e.iter().zip(wts.iter()).map(|(&k, &w)| k as i32 * w).sum::<i32>() };
    // all monomials of R_λ up to degree D
    fn monos_upto(wts: &[i32], maxw: i32) -> Vec<Vec<u16>> {
        let mut out = vec![];
        fn rec(w: &[i32], i: usize, left: i32, cur: &mut Vec<u16>, out: &mut Vec<Vec<u16>>) {
            if i == w.len() {
                out.push(cur.clone());
                return;
            }
            let mut k = 0;
            while k * w[i] <= left {
                cur.push(k as u16);
                rec(w, i + 1, left - k * w[i], cur, out);
                cur.pop();
                k += 1;
            }
        }
        rec(wts, 0, maxw, &mut vec![], &mut out);
        out
    }
    // tau_N on monomials of degree `top`: via truncated identity algebra
    let tr = Alg::<F>::identity(n, &lam);
    // locate top basis element of truncated algebra
    let topidx: Vec<usize> = (0..tr.dim).filter(|&i| tr.deg[i] == top).collect();
    assert_eq!(topidx.len(), 1);
    let topidx = topidx[0];
    let mut tau_cache: HashMap<Vec<u16>, F> = HashMap::new();
    let mut acc = Acc::new(tr.dim);
    let mut tau = |e: &Vec<u16>| -> F {
        if let Some(x) = tau_cache.get(e) {
            return x.clone();
        }
        let mut v: SpVec<F> = vec![(0, F::one())];
        let mut k = 0;
        for (ui, &l) in lam.iter().enumerate() {
            for j in 0..l {
                for _ in 0..e[k] {
                    v = tr.top[ui][j].mul_vec(&v, &mut acc);
                }
                k += 1;
            }
        }
        let x = v.iter().find(|(i, _)| *i as usize == topidx).map(|(_, x)| x.clone()).unwrap_or(F::zero());
        tau_cache.insert(e.clone(), x.clone());
        x
    };
    // invariants Hom(R, B_W) by degree
    let inv = joint_kernel(&a.deg, &crate::hom::closure_ops(a));
    // Hom(B_W, R) of degree delta: unknowns phi_beta in R_λ of degree deg(beta)+delta
    // fibre degrees and the columns T_j(1⊗beta)
    let mut one_idx: HashMap<u32, usize> = HashMap::new(); // beta -> index of 1⊗beta
    let zero_mono = u.rexp.iter().position(|e| e.iter().all(|&k| k == 0)).unwrap() as u32;
    for (i, &(p, b)) in u.pb.iter().enumerate() {
        if p == zero_mono {
            one_idx.insert(b, i);
        }
    }
    let betas: Vec<u32> = {
        let mut v: Vec<u32> = one_idx.keys().copied().collect();
        v.sort();
        v
    };
    let bdeg: HashMap<u32, i32> = betas.iter().map(|&b| (b, a.deg[one_idx[&b]])).collect();
    let mut result = LPoly::zero();
    for (d, fs) in &inv {
        let delta = top - d;
        // unknown layout: for each beta, monomials of degree deg(beta)+delta
        let mut unk: Vec<(u32, Vec<u16>)> = vec![];
        let mut unk_idx: HashMap<(u32, Vec<u16>), usize> = HashMap::new();
        for &b in &betas {
            let dd = bdeg[&b] + delta;
            if dd < 0 {
                continue;
            }
            for m in monos_upto(&wts, dd / 2) {
                if mdeg(&m) == dd {
                    unk_idx.insert((b, m.clone()), unk.len());
                    unk.push((b, m));
                }
            }
        }
        if unk.is_empty() {
            continue;
        }
        // constraints: for each j (top var), beta: sum_{(p',g)} c * p' * phi_g - e_j * phi_b = 0
        let mut rows: Vec<Vec<F>> = vec![];
        let mut var_offset = 0usize;
        for (ui, &l) in lam.iter().enumerate() {
            for j in 0..l {
                let evar = var_offset + j;
                for &b in &betas {
                    let col = &a.top[ui][j].cols[one_idx[&b]];
                    let bot_e = &a.bot[ui][j];
                    let _ = bot_e;
                    // equation polynomial coefficients keyed by monomial
                    let mut eqs: HashMap<Vec<u16>, Vec<F>> = HashMap::new();
                    for (ri, c) in col {
                        let (pp, gg) = u.pb[*ri as usize];
                        let pe = &u.rexp[pp as usize];
                        // p' * phi_g : for each unknown (gg, m) add c at monomial pe+m
                        for (k, (ub, um)) in unk.iter().enumerate() {
                            if *ub != gg {
                                continue;
                            }
                            let mm: Vec<u16> = pe.iter().zip(um.iter()).map(|(x, y)| x + y).collect();
                            eqs.entry(mm).or_insert_with(|| vec![F::zero(); unk.len()])[k].add_assign(c);
                        }
                    }
                    for (k, (ub, um)) in unk.iter().enumerate() {
                        if *ub != b {
                            continue;
                        }
                        let mut mm = um.clone();
                        mm[evar] += 1;
                        eqs.entry(mm).or_insert_with(|| vec![F::zero(); unk.len()])[k].sub_assign(&F::one());
                    }
                    rows.extend(eqs.into_values());
                    if rows.len() > 4 * unk.len() {
                        crate::dense::rref(&mut rows, unk.len());
                    }
                }
            }
            var_offset += l;
        }
        let homs = crate::dense::kernel(&rows, unk.len());
        if homs.is_empty() {
            continue;
        }
        // pairing matrix: f (vector over basis indices) with phi: sum f_(p,b) * tau(p * phi_b)
        let mut mat: Vec<Vec<F>> = vec![];
        for f in fs {
            let mut row = vec![F::zero(); homs.len()];
            for (i, x) in f {
                let (pp, bb) = u.pb[*i as usize];
                let pe = &u.rexp[pp as usize];
                for (h, phi) in homs.iter().enumerate() {
                    for (k, (ub, um)) in unk.iter().enumerate() {
                        if *ub != bb || phi[k].is_zero() {
                            continue;
                        }
                        let mm: Vec<u16> = pe.iter().zip(um.iter()).map(|(x, y)| x + y).collect();
                        let t = tau(&mm);
                        if !t.is_zero() {
                            row[h].fma(&x.mul(&phi[k]), &t);
                        }
                    }
                }
            }
            mat.push(row);
        }
        let rk = crate::dense::rank(&mat, homs.len());
        result.add_term(*d, rk as i128);
    }
    result
}
