//! M4c B1: closed-form generators in krhom's localized form (M4C_PLAN §3.4, milestone B1).
//!
//! # The localized picture
//! A window word Y (a ladder on w adjacent uprights) localizes to one line per flow p (a fixed
//! point: the moved root mask of every rung). Its equivariant Euler weight is
//!
//!   e_Y(p) = Π_rungs e(moved, giver ∖ moved),   e(S, T) = Π_{s∈S, t∈T} (y_s − y_t),
//!
//! evaluated at the window roots y (the convention of `smap::elem_winloc`: Trace = 1/e(S, A∖S)).
//! A word is *flag-type* when it has at most one flow per top coset (e.g. "merge everything, then
//! split": the canonical reduced word X_max of the window's maximal double coset, `canonical_max`).
//! Every generator of SLL §5 (Frob), and every rex move through a maximal coset, is one of the
//! following correspondences (geometrically: pull-back and push-forward along Bott-Samelson maps):
//!
//! | kind | src → tgt | localized entry (q at coset c, p at c) | normalized degree |
//! |---|---|---|---|
//! | `Pull` | flag-type X → Y | 1 | shift(Y) − shift(X) |
//! | `Push` | Y → flag-type X | e_X(c) / e_Y(p) | shift(Y) − shift(X) |
//! | `Rex` | Y1 → Y2, both reduced for the window maximum | (1/e_{Y1}(p)) / Σ_{p'∈c} 1/e_{Y1}(p') | 0 |
//! | `Cap` | [X0.r, X0'.r] (there and back) → [] | 1/e(S, G∖S) at the flows with equal masks | ra − r(g−r) |
//! | `Cup` | [] → [X0.r, X0'.r] | e(S, R) at the flows with equal masks | ra − r(g−r) |
//! | `Swap` | A·B → B·A (disjoint uprights) | 1 at the flow with permuted masks | 0 |
//!
//! Instances: E2 Eps/EpsL = Pull, Eta/EtaL = Push, E3 Iota = Pull, Trace = Push; the right-facing
//! sideways crossing [Is+t−s] → [Is−s+t] (multiplication) = Pull, the left-facing one (∂Δ, SLL
//! eqs. (91)-(92)) = Push, entry e(C, A∖T) of degree ac; the switchback φ (SLL Thm 6.7, EKo eq.
//! (108)) = Pull (all ones: f⊗g ↦ f⊗1⊗g), ψ = Rex into X_max; the atomic braid (CORE Lemma 3.12) =
//! Pull ∘ (up-up, down-down) ∘ Rex through the merged expression; far commutation (CORE Lemma
//! 3.13) = Swap.
//!
//! # Trust
//! Nothing here is trusted. A closed form becomes a `tmpl::Template` (`TmplSource::Added`) whose
//! exact PMat is obtained by interpolation of L_tgt⁻¹ f(y) L_src (`interpolate`); the template is a
//! genuine leaf only after `TmplLib::certify_added` (exact degree, exact commutation with every
//! top operator, independent polynomiality certificate at another seed). Composites (switchback
//! ψ∘φ, braid moves) are `Comp`/`Local` DAGs of certified templates and Elem Assoc nodes that
//! passed `validate_assoc`.
use crate::coset::{Mask, Roots};
use crate::field::Field;
use crate::local::{esym, fixed_points, inverse, loc_matrix, matmul};
use crate::mpoly::{mono_exp, Mono, Poly};
use crate::pmat::{PMat, PVec};
use crate::sexpr::{rung_ends, Expr, FlowIndex};
use crate::smap::{e_pair, elem_winloc, points, Arena, ElemKind, MapId, MapNode};
use crate::tmpl::{poly_certificate, CertReport, Template, TmplId, TmplLib, TmplSource, WinLoc, Window};
use crate::web::{Dir, Rung};
use std::collections::BTreeMap;

/// Largest window (flows on either side) for which `interpolate` builds dense localization
/// matrices (L is n×n and is inverted once per interpolation point).
pub const MAX_INTERP_FLOWS: usize = 1200;

/// Largest window for which `validate_assoc` runs the polynomiality certificate (M4C_PLAN §3.4).
pub const ASSOC_CERT_FLOWS: usize = 500;

// ------------------------------------------------------------------ words

/// labels after applying the rungs (None if a rung is invalid)
pub fn word_top(labels: &[u8], rungs: &[Rung]) -> Option<Vec<u8>> {
    Expr::new(labels, rungs).ok().map(|e| e.top())
}

/// The canonical reduced word of the maximal double coset of a window: merge every upright into
/// upright 0 (E rungs from the right), then split off uprights 1.. in order (F rungs). Zero moves
/// are skipped. Flag-type: one flow per top coset, all cosets.
pub fn canonical_max(bottom: &[u8], top: &[u8]) -> Vec<Rung> {
    assert_eq!(bottom.len(), top.len(), "canonical_max: lengths differ");
    let w = bottom.len();
    let mut out = vec![];
    let mut carry = 0usize;
    for j in (1..w).rev() {
        carry += bottom[j] as usize;
        if carry > 0 {
            out.push(Rung::e(j - 1, carry));
        }
    }
    for j in 0..w.saturating_sub(1) {
        let rest: usize = top[j + 1..].iter().map(|&x| x as usize).sum();
        if rest > 0 {
            out.push(Rung::f(j, rest));
        }
    }
    out
}

/// at most one flow per top coset
pub fn is_flag_type(fi: &FlowIndex) -> bool {
    fi.start.windows(2).all(|w| w[1] - w[0] <= 1)
}

/// e_Y(p) = Π_rungs e(moved, giver ∖ moved) at window roots y
pub fn euler<F: Field>(fi: &FlowIndex, f: u32, y: &[F]) -> F {
    let mut conf = Roots::new(&fi.expr.bottom).bottom().0;
    let mut v = F::one();
    for (r, &m) in fi.expr.rungs.iter().zip(fi.flow(f)) {
        let (rec, giv) = rung_ends(r);
        v = v.mul(&e_pair(m, conf[giv] & !m, y));
        conf[giv] &= !m;
        conf[rec] |= m;
    }
    v
}

// ------------------------------------------------------------------ closed forms

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CfKind {
    Pull,
    Push,
    Rex,
    Cap,
    Cup,
    /// src = A·B with |A| = na rungs, tgt = B·A
    Swap { na: usize },
}

/// A closed-form window map (see the module doc for the formulas).
#[derive(Clone, Debug)]
pub struct ClosedForm {
    pub name: String,
    pub kind: CfKind,
    pub src: Window,
    pub tgt: Window,
    pub deg: i32,
    pub sfi: FlowIndex,
    pub tfi: FlowIndex,
}

fn shift(w: &Window) -> i32 {
    w.expr().shift()
}

impl ClosedForm {
    /// checks the shape conditions of `kind` and computes the normalized degree
    pub fn new(name: &str, kind: CfKind, src: Window, tgt: Window) -> Result<ClosedForm, String> {
        let (se, te) = (Expr::new(&src.labels, &src.rungs)?, Expr::new(&tgt.labels, &tgt.rungs)?);
        if src.labels != tgt.labels {
            return Err(format!("{}: window bottoms differ ({:?} vs {:?})", name, src.labels, tgt.labels));
        }
        if se.top() != te.top() {
            return Err(format!("{}: window tops differ ({:?} vs {:?})", name, se.top(), te.top()));
        }
        let sfi = FlowIndex::build(&se);
        let tfi = FlowIndex::build(&te);
        let (ss, ts) = (se.shift(), te.shift());
        let there_and_back = |w: &Window| -> Result<(usize, usize, usize), String> {
            // [X_u.r, X'_u.r] on two uprights returning to the bottom labels: (r, giver label g, receiver label a)
            if w.labels.len() != 2 || w.rungs.len() != 2 || w.rungs[0].r != w.rungs[1].r || w.rungs[0].dir == w.rungs[1].dir || w.rungs[0].i != 0 || w.rungs[1].i != 0 {
                return Err(format!("{}: {:?} is not a two-rung there-and-back window", name, w.rungs));
            }
            let (rec, giv) = rung_ends(&w.rungs[0]);
            Ok((w.rungs[0].r, w.labels[giv] as usize, w.labels[rec] as usize))
        };
        let deg = match &kind {
            CfKind::Pull => {
                if !is_flag_type(&sfi) {
                    return Err(format!("{}: Pull needs a flag-type source, got {:?}", name, src.rungs));
                }
                ts - ss
            }
            CfKind::Push => {
                if !is_flag_type(&tfi) {
                    return Err(format!("{}: Push needs a flag-type target, got {:?}", name, tgt.rungs));
                }
                ss - ts
            }
            CfKind::Rex => {
                if ss != ts {
                    return Err(format!("{}: Rex between words of different shifts ({} vs {}): not both reduced", name, ss, ts));
                }
                if sfi.support() != tfi.support() {
                    return Err(format!("{}: Rex between words of different supports", name));
                }
                0
            }
            CfKind::Cap => {
                if !tgt.rungs.is_empty() {
                    return Err(format!("{}: Cap target must be empty", name));
                }
                let (r, g, _) = there_and_back(&src)?;
                -2 * (r * (g - r)) as i32 - ss
            }
            CfKind::Cup => {
                if !src.rungs.is_empty() {
                    return Err(format!("{}: Cup source must be empty", name));
                }
                let (r, _, a) = there_and_back(&tgt)?;
                2 * (r * a) as i32 + ts
            }
            CfKind::Swap { na } => {
                let na = *na;
                if na == 0 || na >= src.rungs.len() {
                    return Err(format!("{}: Swap split {} out of range", name, na));
                }
                let (a, b) = src.rungs.split_at(na);
                let want: Vec<Rung> = b.iter().chain(a).cloned().collect();
                if tgt.rungs != want {
                    return Err(format!("{}: Swap target {:?} != B·A {:?}", name, tgt.rungs, want));
                }
                let ups = |rs: &[Rung]| -> u64 { rs.iter().fold(0u64, |m, r| m | (3u64 << r.i)) };
                if ups(a) & ups(b) != 0 {
                    return Err(format!("{}: Swap blocks share an upright ({:?} | {:?})", name, a, b));
                }
                0
            }
        };
        Ok(ClosedForm { name: name.to_string(), kind, src, tgt, deg, sfi, tfi })
    }

    /// the localized map at window roots y (rows tgt flows, cols src flows, sorted flow order)
    pub fn winloc<F: Field>(&self, y: &[F]) -> WinLoc<F> {
        let (ns, nt) = (self.sfi.nflows(), self.tfi.nflows());
        let mut cols: Vec<Vec<(u32, F)>> = vec![vec![]; ns];
        let (sfi, tfi) = (&self.sfi, &self.tfi);
        match &self.kind {
            CfKind::Pull => {
                for j in 0..ns as u32 {
                    let c = sfi.coset_of_flow(j);
                    for i in tfi.range(c) {
                        cols[j as usize].push((i, F::one()));
                    }
                }
            }
            CfKind::Push => {
                for j in 0..ns as u32 {
                    let c = sfi.coset_of_flow(j);
                    for i in tfi.range(c) {
                        cols[j as usize].push((i, euler(tfi, i, y).mul(&euler(sfi, j, y).inv())));
                    }
                }
            }
            CfKind::Rex => {
                for (k, &c) in sfi.support().iter().enumerate() {
                    let r = sfi.start[k]..sfi.start[k + 1];
                    let inv: Vec<F> = r.clone().map(|j| euler(sfi, j, y).inv()).collect();
                    let mut s = F::zero();
                    for v in &inv {
                        s.add_assign(v);
                    }
                    assert!(!s.is_zero(), "{}: Σ 1/e_Y vanishes at coset {} (word not reduced for the window maximum?)", self.name, c);
                    let si = s.inv();
                    for (j, v) in r.zip(inv) {
                        let w = v.mul(&si);
                        for i in tfi.range(c) {
                            cols[j as usize].push((i, w.clone()));
                        }
                    }
                }
            }
            CfKind::Cap | CfKind::Cup => {
                let (fi, cap) = if self.kind == CfKind::Cap { (sfi, true) } else { (tfi, false) };
                let r0 = &fi.expr.rungs[0];
                let (rec, giv) = rung_ends(r0);
                let bot = Roots::new(&fi.expr.bottom).bottom().0;
                for f in 0..fi.nflows() as u32 {
                    let m = fi.flow(f);
                    if m[0] != m[1] {
                        continue;
                    }
                    if cap {
                        cols[f as usize].push((0, e_pair(m[0], bot[giv] & !m[0], y).inv()));
                    } else {
                        cols[0].push((f, e_pair(m[0], bot[rec], y)));
                    }
                }
            }
            CfKind::Swap { na } => {
                let na = *na;
                for j in 0..ns as u32 {
                    let m = sfi.flow(j);
                    let moved: Vec<Mask> = m[na..].iter().chain(&m[..na]).cloned().collect();
                    let i = tfi.locate(&moved).expect("Swap: permuted flow is not a target flow");
                    cols[j as usize].push((i, F::one()));
                }
            }
        }
        let mut col_start = vec![0u32];
        let (mut rows, mut vals) = (vec![], vec![]);
        for mut c in cols {
            c.sort_unstable_by_key(|e| e.0);
            for (i, v) in c {
                rows.push(i);
                vals.push(v);
            }
            col_start.push(rows.len() as u32);
        }
        WinLoc { nrows: nt as u32, col_start, rows, vals }
    }

    /// flows of (src, tgt)
    pub fn flows(&self) -> (usize, usize) {
        (self.sfi.nflows(), self.tfi.nflows())
    }
}

// ------------------------------------------------------------------ exact PMat by interpolation

fn positions(fps: &[crate::local::FixedPoint], fi: &FlowIndex) -> Vec<usize> {
    fps.iter().map(|p| fi.locate(&p.moved).expect("fixed point is not a flow") as usize).collect()
}

/// The exact PMat (over the window bottom ring, fibre bases of `sbim`) of the bimodule map whose
/// localization is `f`: Φ(y) = L_tgt(y)⁻¹ f(y) L_src(y), every entry interpolated as a polynomial
/// of its forced degree from (#monomials) points and checked at 2 more. Err if an entry is not
/// polynomial (the closed form is not a bimodule map of degree `deg`) or the window is too large.
pub fn interpolate<F: Field>(src: &Window, tgt: &Window, deg: i32, f: &mut dyn FnMut(&[F]) -> WinLoc<F>, seed: u64) -> Result<PMat<F>, String> {
    let sb = src.expr().sbim::<F>();
    let tb = tgt.expr().sbim::<F>();
    let (sfi, tfi) = (FlowIndex::build(&src.expr()), FlowIndex::build(&tgt.expr()));
    let (ns, nt) = (sfi.nflows(), tfi.nflows());
    if ns.max(nt) > MAX_INTERP_FLOWS {
        return Err(format!("interpolate: window of {} x {} flows exceeds {}", nt, ns, MAX_INTERP_FLOWS));
    }
    let sfp = fixed_points(&sb);
    let tfp = fixed_points(&tb);
    let spos = positions(&sfp, &sfi);
    let tpos = positions(&tfp, &tfi);
    let u = deg - tb.shift + sb.shift;
    let ring = sb.ring.clone();
    let mut by_d: BTreeMap<i32, Vec<(usize, usize)>> = BTreeMap::new();
    let mut neg: Vec<(usize, usize)> = vec![];
    for bp in 0..nt {
        for b in 0..ns {
            let d = u + sb.fdeg[b] - tb.fdeg[bp];
            if d < 0 {
                neg.push((bp, b));
            } else {
                by_d.entry(d).or_default().push((bp, b));
            }
        }
    }
    let monos: BTreeMap<i32, Vec<Mono>> = by_d.keys().map(|&d| (d, ring.monos_of_deg(d))).collect();
    let max_monos = monos.values().map(|v| v.len()).max().unwrap_or(0);
    let npts = max_monos + 2;
    let kw: usize = src.labels.iter().map(|&a| a as usize).sum();
    let off = crate::local::root_offsets(&src.labels.iter().map(|&a| a as usize).collect::<Vec<_>>());
    'redraw: for attempt in 0..3u64 {
        let pts = points::<F>(kw, seed.wrapping_add(0x1eaf_0000 + attempt), npts, 0);
        let mut phis: Vec<Vec<Vec<F>>> = Vec::with_capacity(npts);
        let mut evals: Vec<Vec<Vec<F>>> = Vec::with_capacity(npts);
        for p in &pts {
            let y = &p.x;
            let wl = f(y);
            if wl.nrows as usize != nt || wl.ncols() != ns {
                return Err(format!("interpolate: WinLoc shape {}x{} != {}x{}", wl.nrows, wl.ncols(), nt, ns));
            }
            let d = wl.dense();
            let fd: Vec<Vec<F>> = (0..nt).map(|q| (0..ns).map(|pp| d[tpos[q]][spos[pp]].clone()).collect()).collect();
            let ls = loc_matrix(&sb, &sfp, y);
            let lt = loc_matrix(&tb, &tfp, y);
            let ltinv = inverse(&lt).ok_or("interpolate: L_tgt singular")?;
            let phi = matmul(&matmul(&ltinv, &fd), &ls);
            for &(bp, b) in &neg {
                if !phi[bp][b].is_zero() {
                    return Err(format!("interpolate: entry ({}, {}) of negative degree is nonzero: not a map of degree {}", bp, b, deg));
                }
            }
            phis.push(phi);
            evals.push(src.labels.iter().enumerate().map(|(i, &l)| esym(crate::coset::low_mask(l as usize) << off[i], y)).collect());
        }
        let mono_val = |k: usize, m: Mono| -> F {
            let mut v = F::one();
            let mut var = 0usize;
            for (i, &l) in src.labels.iter().enumerate() {
                for j in 0..l as usize {
                    for _ in 0..mono_exp(m, var) {
                        v = v.mul(&evals[k][i][j + 1]);
                    }
                    var += 1;
                }
            }
            v
        };
        let mut cols: Vec<PVec<F>> = vec![vec![]; ns];
        for (&d, ents) in &by_d {
            let ms = &monos[&d];
            let nm = ms.len();
            if nm == 0 {
                continue;
            }
            let vd: Vec<Vec<F>> = (0..nm).map(|k| ms.iter().map(|&m| mono_val(k, m)).collect()).collect();
            let vinv = match inverse(&vd) {
                Some(v) => v,
                None => continue 'redraw,
            };
            let vals: Vec<Vec<F>> = (0..nm).map(|k| ents.iter().map(|&(bp, b)| phis[k][bp][b].clone()).collect()).collect();
            let coef = matmul(&vinv, &vals);
            for k in nm..npts {
                let row: Vec<F> = ms.iter().map(|&m| mono_val(k, m)).collect();
                for (e, &(bp, b)) in ents.iter().enumerate() {
                    let mut s = F::zero();
                    for (mi, r) in row.iter().enumerate() {
                        s.fma(r, &coef[mi][e]);
                    }
                    if s != phis[k][bp][b] {
                        return Err(format!("interpolate: entry ({}, {}) of degree {} is not a polynomial (closed form is not a bimodule map)", bp, b, d));
                    }
                }
            }
            for (e, &(bp, b)) in ents.iter().enumerate() {
                let mut t: Vec<(Mono, F)> = (0..nm).filter(|&mi| !coef[mi][e].is_zero()).map(|mi| (ms[mi], coef[mi][e].clone())).collect();
                if t.is_empty() {
                    continue;
                }
                t.sort_unstable_by_key(|x| x.0);
                cols[b].push((bp as u32, Poly { t }));
            }
        }
        for c in cols.iter_mut() {
            c.sort_unstable_by_key(|e| e.0);
        }
        return Ok(PMat { nrows: nt, ncols: ns, cols });
    }
    Err("interpolate: Vandermonde singular after 3 redraws".into())
}

/// Result of turning a closed form into a certified template.
#[derive(Clone, Debug)]
pub struct Built {
    pub tmpl: TmplId,
    pub deg: i32,
    pub flows: (usize, usize),
    pub cert: CertReport,
    pub interp_secs: f64,
}

/// Interpolate the closed form, add it as an `Added` template and certify it
/// (`TmplLib::certify_added` with an independent seed). The template is a genuine leaf iff Ok.
pub fn build<F: Field>(lib: &mut TmplLib<F>, cf: &ClosedForm, seed: u64) -> Result<Built, String> {
    let t0 = std::time::Instant::now();
    let phi = interpolate::<F>(&cf.src, &cf.tgt, cf.deg, &mut |y: &[F]| cf.winloc(y), seed).map_err(|e| format!("{}: {}", cf.name, e))?;
    let interp_secs = t0.elapsed().as_secs_f64();
    let tm = Template::new(cf.src.clone(), cf.tgt.clone(), cf.deg, phi, TmplSource::Added { certified: false });
    let id = lib.add_template(tm);
    let cert = lib.certify_added(id, seed ^ 0xce57_0000_b1b1).map_err(|e| format!("{}: certify_added: {}", cf.name, e))?;
    // the certified template localizes to the closed form (fresh point)
    let y = points::<F>(cf.src.labels.iter().map(|&a| a as usize).sum(), seed ^ 0x9e37, 1, 7)[0].x.clone();
    if lib.localize(id, &y).dense() != cf.winloc(&y).dense() {
        return Err(format!("{}: certified template does not localize to the closed form", cf.name));
    }
    Ok(Built { tmpl: id, deg: cf.deg, flows: cf.flows(), cert, interp_secs })
}

// ------------------------------------------------------------------ generator constructors

fn win(labels: &[u8], rungs: &[Rung]) -> Window {
    Window::new(labels, rungs)
}

/// Sideways crossings for a one-rung word Y on two uprights (SLL eqs. (84), (91)-(92)):
/// (right-facing: X_max → Y, multiplication = Pull; left-facing: Y → X_max, ∂Δ = Push), where
/// X_max = `canonical_max`. For Y = F0.b over (a+b, c) this is [Is+t−s] ⇄ [Is−s+t] on blocks
/// (a, b, c), degree ac; for Y = E0.b over (a, b+c) the mirror image.
pub fn sideways(labels: &[u8], y: Rung) -> Result<(ClosedForm, ClosedForm), String> {
    let top = word_top(labels, &[y]).ok_or("sideways: invalid rung")?;
    let x = canonical_max(labels, &top);
    let r = ClosedForm::new(&format!("cross-right {:?} {:?}", labels, y), CfKind::Pull, win(labels, &x), win(labels, &[y]))?;
    let l = ClosedForm::new(&format!("cross-left {:?} {:?}", labels, y), CfKind::Push, win(labels, &[y]), win(labels, &x))?;
    Ok((r, l))
}

/// The F-type sideways crossing on blocks (a, b, c): window (a+b, c), Y = [F0.b], X = [E0.c, F0.(b+c)].
pub fn sideways_abc(a: u8, b: u8, c: u8) -> Result<(ClosedForm, ClosedForm), String> {
    sideways(&[a + b, c], Rung::f(0, b as usize))
}

/// Switchback φ: X_max → Y (Pull) and ψ: Y → X_max (Rex, normalized push) for any word Y that is
/// a reduced expression of the window's maximal coset (SLL Thm 6.7; ψ∘φ = id).
pub fn switchback(labels: &[u8], y: &[Rung]) -> Result<(ClosedForm, ClosedForm), String> {
    let top = word_top(labels, y).ok_or("switchback: invalid word")?;
    let x = canonical_max(labels, &top);
    let phi = ClosedForm::new(&format!("phi {:?} {:?}", labels, y), CfKind::Pull, win(labels, &x), win(labels, y))?;
    let psi = ClosedForm::new(&format!("psi {:?} {:?}", labels, y), CfKind::Rex, win(labels, y), win(labels, &x))?;
    if phi.deg != 0 || psi.deg != 0 {
        return Err(format!("switchback: {:?} is not a reduced expression of the window maximum (degrees {}, {})", y, phi.deg, psi.deg));
    }
    Ok((phi, psi))
}

/// The type A switchback of EKo eq. (108) on two uprights (a, m), N = a+m, 1 ≤ b < m:
/// [J + s_a − s_b] = E0.m F0.(N−b) (merge, split as (b, N−b)) ⇄ [J − s_c + s_a − s_b + s_c] =
/// E0.b F0.a (needs b ≤ m and a ≥ 1).
pub fn switchback_eko(a: u8, m: u8, b: u8) -> Result<(ClosedForm, ClosedForm), String> {
    switchback(&[a, m], &[Rung::e(0, b as usize), Rung::f(0, a as usize)])
}

/// General cw/ccw cap and cup E^(r)F^(r) ⇄ 1 on two uprights (a, b): the there-and-back window
/// [first, back] with first = E0.r (giver upright 1) or F0.r (giver upright 0).
pub fn cup_cap(labels: &[u8], first: Rung) -> Result<(ClosedForm, ClosedForm), String> {
    let back = Rung { i: first.i, dir: if first.dir == Dir::E { Dir::F } else { Dir::E }, r: first.r };
    let w = win(labels, &[first, back]);
    let e = win(labels, &[]);
    let cap = ClosedForm::new(&format!("cap {:?} {:?}", labels, first), CfKind::Cap, w.clone(), e.clone())?;
    let cup = ClosedForm::new(&format!("cup {:?} {:?}", labels, first), CfKind::Cup, e, w)?;
    Ok((cap, cup))
}

/// Far commutation A·B → B·A for rung blocks on disjoint uprights (CORE Lemma 3.13).
pub fn far_commute(labels: &[u8], a: &[Rung], b: &[Rung]) -> Result<ClosedForm, String> {
    let src: Vec<Rung> = a.iter().chain(b).cloned().collect();
    let tgt: Vec<Rung> = b.iter().chain(a).cloned().collect();
    ClosedForm::new(&format!("far {:?} {:?}|{:?}", labels, a, b), CfKind::Swap { na: a.len() }, win(labels, &src), win(labels, &tgt))
}

/// Atomic crossing of the blocks on uprights (u, u+1) with labels (x, y): merge, split as (y, x).
pub fn atomic(u: usize, x: u8, y: u8) -> Vec<Rung> {
    let mut v = vec![];
    if y > 0 {
        v.push(Rung::e(u, y as usize));
    }
    if x > 0 {
        v.push(Rung::f(u, x as usize));
    }
    v
}

/// The atomic braid H0 H1 H0 ⇄ H1 H0 H1 on cells (x, y, z) (CORE Lemma 3.12) in krhom words,
/// with the switchback windows of the route through the merged expression:
/// A = a1 a2 a1 = E0.y F0.x | E1.z F1.x E0.z | F0.y, B = a2 a1 a2 = E1.z | F1.y E0.z F0.x E1.y | F1.x.
#[derive(Clone, Debug)]
pub struct BraidWords {
    pub cells: (u8, u8, u8),
    pub a: Vec<Rung>,
    pub b: Vec<Rung>,
    /// merged words: M = A with its switchback window collapsed, M' = B likewise
    pub m: Vec<Rung>,
    pub mp: Vec<Rung>,
    /// (at, window labels, long word, short word) for A and for B
    pub win_a: (usize, Vec<u8>, Vec<Rung>, Vec<Rung>),
    pub win_b: (usize, Vec<u8>, Vec<Rung>, Vec<Rung>),
}

pub fn braid_words(x: u8, y: u8, z: u8) -> Result<BraidWords, String> {
    if x == 0 || y == 0 || z == 0 {
        return Err("braid_words: cells must be positive".into());
    }
    let (xu, yu, zu) = (x as usize, y as usize, z as usize);
    let a = vec![Rung::e(0, yu), Rung::f(0, xu), Rung::e(1, zu), Rung::f(1, xu), Rung::e(0, zu), Rung::f(0, yu)];
    let b = vec![Rung::e(1, zu), Rung::f(1, yu), Rung::e(0, zu), Rung::f(0, xu), Rung::e(1, yu), Rung::f(1, xu)];
    let la = vec![x + y, 0, z];
    let long_a = a[1..5].to_vec();
    let short_a = canonical_max(&la, &word_top(&la, &long_a).ok_or("braid: bad A window")?);
    let lb = vec![x, y + z, 0];
    let long_b = b[1..5].to_vec();
    let short_b = canonical_max(&lb, &word_top(&lb, &long_b).ok_or("braid: bad B window")?);
    let m: Vec<Rung> = a[..1].iter().chain(&short_a).chain(&a[5..]).cloned().collect();
    let mp: Vec<Rung> = b[..1].iter().chain(&short_b).chain(&b[5..]).cloned().collect();
    Ok(BraidWords { cells: (x, y, z), a, b, m, mp, win_a: (1, la, long_a, short_a), win_b: (1, lb, long_b, short_b) })
}

// ------------------------------------------------------------------ Assoc validator

/// wall set (cut positions 1..K−1) of a composition; zero labels give repeated cuts (a set)
fn walls(l: &[u8]) -> Vec<usize> {
    let k: usize = l.iter().map(|&a| a as usize).sum();
    let mut s = 0usize;
    let mut out = vec![];
    for &a in &l[..l.len().saturating_sub(1)] {
        s += a as usize;
        if s > 0 && s < k && out.last() != Some(&s) {
            out.push(s);
        }
    }
    out
}

fn subset(a: &[usize], b: &[usize]) -> bool {
    a.iter().all(|x| b.contains(x))
}

/// The EKLP multistep expression of a krhom word (M4C_PLAN §3.4 "up-up/down-down"): the sequence
/// of parabolic subsets (as wall sets) visited by the rungs (each rung: split the giver, merge),
/// with repeats removed and the interior points of monotone runs dropped. Two words with equal
/// multistep expressions have canonically isomorphic BS bimodules (SLL §4.2, eq. (88)).
pub fn multistep(labels: &[u8], rungs: &[Rung]) -> Result<Vec<Vec<usize>>, String> {
    let e = Expr::new(labels, rungs)?;
    let lv = e.levels();
    let mut seq: Vec<Vec<usize>> = vec![walls(&lv[0])];
    for (t, r) in rungs.iter().enumerate() {
        let (rec, giv) = rung_ends(r);
        // refinement: the moved part sits next to the receiver
        let mut nu: Vec<u8> = vec![];
        for (j, &a) in lv[t].iter().enumerate() {
            if j == giv {
                let (mv, st) = (r.r as u8, a - r.r as u8);
                if rec < giv {
                    nu.push(mv);
                    nu.push(st);
                } else {
                    nu.push(st);
                    nu.push(mv);
                }
            } else {
                nu.push(a);
            }
        }
        seq.push(walls(&nu));
        seq.push(walls(&lv[t + 1]));
    }
    loop {
        let n0 = seq.len();
        seq.dedup();
        let mut out: Vec<Vec<usize>> = vec![];
        for i in 0..seq.len() {
            if i > 0 && i + 1 < seq.len() {
                let (p, c, n) = (&seq[i - 1], &seq[i], &seq[i + 1]);
                let mono = (subset(p, c) && subset(c, n)) || (subset(c, p) && subset(n, c));
                if mono {
                    continue;
                }
            }
            out.push(seq[i].clone());
        }
        seq = out;
        if seq.len() == n0 {
            break;
        }
    }
    Ok(seq)
}

/// Validate an E4 Assoc window pair (ws → wt) and record its constant 1 in the arena
/// (M4C_PLAN §3.4, C_m4c M5): (1) the multistep expressions are equal; (2) both windows' flows are
/// in bijection with equal top configurations and the shifts agree (the `Arena::add` conditions);
/// (3) the polynomiality certificate of the closed form when both windows have ≤ 500 flows.
/// Refuses otherwise; never records an unvalidated constant.
pub fn validate_assoc<F: Field>(ar: &mut Arena<F>, ws: &Expr, wt: &Expr, seed: u64) -> Result<CertReport, String> {
    if ws.bottom != wt.bottom {
        return Err("validate_assoc: window bottoms differ".into());
    }
    let (ma, mb) = (multistep(&ws.bottom, &ws.rungs)?, multistep(&wt.bottom, &wt.rungs)?);
    if ma != mb {
        return Err(format!("validate_assoc: multistep expressions differ: {:?} vs {:?}", ma, mb));
    }
    let (a, b) = (FlowIndex::build(ws), FlowIndex::build(wt));
    if a.support() != b.support() || a.max_block() != 1 || b.max_block() != 1 {
        return Err("validate_assoc: window flows are not in bijection with equal top configurations".into());
    }
    if ws.shift() != wt.shift() {
        return Err("validate_assoc: shifts differ".into());
    }
    if a.nflows().max(b.nflows()) > ASSOC_CERT_FLOWS {
        return Err(format!("validate_assoc: {} flows > {} (refused: no certificate)", a.nflows(), ASSOC_CERT_FLOWS));
    }
    let rep = poly_certificate::<F>(&Window::new(&ws.bottom, &ws.rungs), &Window::new(&wt.bottom, &wt.rungs), 0, &mut |y: &[F]| elem_winloc(ElemKind::Assoc, &a, &b, y), seed)?;
    ar.certify_assoc(ws, wt, F::one());
    Ok(rep)
}

// ------------------------------------------------------------------ composites

/// Add a closed form to the arena's library and certify it through the arena (flags refreshed).
pub fn build_in<F: Field>(ar: &mut Arena<F>, cf: &ClosedForm, seed: u64) -> Result<Built, String> {
    let t0 = std::time::Instant::now();
    let phi = interpolate::<F>(&cf.src, &cf.tgt, cf.deg, &mut |y: &[F]| cf.winloc(y), seed).map_err(|e| format!("{}: {}", cf.name, e))?;
    let interp_secs = t0.elapsed().as_secs_f64();
    let tm = Template::new(cf.src.clone(), cf.tgt.clone(), cf.deg, phi, TmplSource::Added { certified: false });
    let id = ar.lib().add_template(tm);
    let cert = ar.certify_template(id, seed ^ 0xce57_0000_b1b1).map_err(|e| format!("{}: certify_added: {}", cf.name, e))?;
    Ok(Built { tmpl: id, deg: cf.deg, flows: cf.flows(), cert, interp_secs })
}

/// The braid move A → B (CORE Lemma 3.12) as a DAG: Local(φ_B: short_b → long_b) ∘ (M → M') ∘
/// Local(ψ_A: long_a → short_a). The middle step M → M' (up-up and down-down) is an Elem Assoc
/// when `validate_assoc` accepts it, otherwise a certified Pull template between the two
/// flag-type merged words. Returns (map, report lines).
pub fn braid_dag<F: Field>(ar: &mut Arena<F>, bw: &BraidWords, seed: u64) -> Result<(MapId, Vec<String>), String> {
    let lam = vec![bw.cells.0, bw.cells.1, bw.cells.2];
    let mut log = vec![];
    let (at_a, la, long_a, short_a) = bw.win_a.clone();
    let (at_b, lb, long_b, short_b) = bw.win_b.clone();
    let psi_a = ClosedForm::new("braid psi_A", CfKind::Rex, win(&la, &long_a), win(&la, &short_a))?;
    let phi_b = ClosedForm::new("braid phi_B", CfKind::Pull, win(&lb, &short_b), win(&lb, &long_b))?;
    let ta = build_in(ar, &psi_a, seed)?;
    log.push(format!("psi_A {:?}->{:?} flows {:?} certified", long_a, short_a, ta.flows));
    let tb = build_in(ar, &phi_b, seed + 1)?;
    log.push(format!("phi_B {:?}->{:?} flows {:?} certified", short_b, long_b, tb.flows));
    let ea = ar.expr(Expr::new(&lam, &bw.a)?);
    let em = ar.expr(Expr::new(&lam, &bw.m)?);
    let emp = ar.expr(Expr::new(&lam, &bw.mp)?);
    let eb = ar.expr(Expr::new(&lam, &bw.b)?);
    let s1 = ar.add(MapNode::Local { src: ea, tgt: em, at: at_a, first: 0, tmpl: ta.tmpl }, "braid psi_A")?;
    let s3 = ar.add(MapNode::Local { src: emp, tgt: eb, at: at_b, first: 0, tmpl: tb.tmpl }, "braid phi_B")?;
    // middle: M -> M'
    let (me, mpe) = (Expr::new(&lam, &bw.m)?, Expr::new(&lam, &bw.mp)?);
    let s2 = match validate_assoc(ar, &me, &mpe, seed + 2) {
        Ok(_) => {
            log.push(format!("M -> M' Elem Assoc validated ({} flows)", FlowIndex::build(&me).nflows()));
            ar.add(MapNode::Elem { kind: ElemKind::Assoc, src: em, tgt: emp, at: 0, upright: 0, scalar: F::one() }, "braid up-up/down-down")?
        }
        Err(e) => {
            log.push(format!("M -> M' Assoc refused ({}); using a Pull template", e));
            let mid = ClosedForm::new("braid M->M'", CfKind::Pull, win(&lam, &bw.m), win(&lam, &bw.mp))?;
            let tm = build_in(ar, &mid, seed + 3)?;
            ar.add(MapNode::Local { src: em, tgt: emp, at: 0, first: 0, tmpl: tm.tmpl }, "braid M->M'")?
        }
    };
    let c1 = ar.add(MapNode::Comp(s2, s1), "braid")?;
    let c2 = ar.add(MapNode::Comp(s3, c1), "braid")?;
    Ok((c2, log))
}
