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
use crate::combinat::coset::{Mask, Roots};
use crate::algebra::field::Field;
use crate::soergel::local::{esym, fixed_points, inverse, loc_matrix, matmul};
use crate::algebra::mpoly::{mono_exp, Mono, Poly};
use crate::algebra::pmat::{PMat, PVec};
use crate::colored::sexpr::{rung_ends, Expr, FlowIndex};
use crate::colored::smap::{e_pair, elem_winloc, points, Arena, ElemKind, MapId, MapNode};
use crate::colored::tmpl::{poly_certificate, CertReport, Template, TmplId, TmplLib, TmplSource, WinLoc, Window};
use crate::combinat::web::{Dir, Rung};
use std::collections::BTreeMap;

/// Largest window (flows on either side) for which `interpolate` builds dense localization
/// matrices (L is n×n and is inverted once per interpolation point).
pub const MAX_INTERP_FLOWS: usize = 1200;

/// Largest window for which `validate_assoc` runs the polynomiality certificate (M4C_PLAN §3.4).
pub const ASSOC_CERT_FLOWS: usize = 500;

/// Memory/time gate of `interpolate` and of the polynomiality certificate (review wave 4,
/// finding 1): at most this many monomials in one entry degree (the Vandermonde is nm × nm,
/// inverted as an nm × 2nm rref; nm³ field operations).
pub const MAX_INTERP_MONOS: u64 = 6_000;
/// ... and at most this many stored localization entries, (#points) · ns · nt, where
/// #points = max monomials + 4 (the certificate's count; interpolate uses + 2).
pub const MAX_INTERP_CELLS: u128 = 250_000_000;

/// Predicted size of an interpolation / polynomiality certificate of a degree-`deg` map
/// src → tgt, computed from the words alone (no SBim, no GrTable): fibre degrees of a window
/// word range over [0, −2·shift], so the entry degrees range over
/// [u − maxf_tgt, u + maxf_src] with u = deg − shift_tgt + shift_src.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct InterpCost {
    pub ns: u128,
    pub nt: u128,
    pub max_monos: u64,
    pub cells: u128,
}

/// The gate itself: Err (without building anything) when the window exceeds MAX_INTERP_FLOWS
/// flows on either side, MAX_INTERP_MONOS monomials in some entry degree, or MAX_INTERP_CELLS
/// stored entries. Called before any SBim is built by `interpolate`, `build`, `build_in`,
/// `validate_assoc` and `theorem_vs_numeric`.
pub fn interp_gate(src: &Window, tgt: &Window, deg: i32) -> Result<InterpCost, String> {
    let se = Expr::new(&src.labels, &src.rungs)?;
    let te = Expr::new(&tgt.labels, &tgt.rungs)?;
    let (ns, nt) = (se.nflows(), te.nflows());
    if ns.max(nt) > MAX_INTERP_FLOWS as u128 {
        return Err(format!("interpolate: window of {} x {} flows exceeds {} (refused before any SBim build)", nt, ns, MAX_INTERP_FLOWS));
    }
    let (ss, ts) = (se.shift(), te.shift());
    let (maxfs, maxft) = (-2 * ss, -2 * ts);
    let u = deg - ts + ss;
    let mut max_monos = 0u64;
    let mut d = (u - maxft).max(0);
    while d <= u + maxfs {
        max_monos = max_monos.max(crate::library::rexplan::monos_of_deg(&src.labels, d));
        d += 1;
    }
    let cells = (max_monos as u128 + 4).saturating_mul(ns).saturating_mul(nt);
    let c = InterpCost { ns, nt, max_monos, cells };
    if max_monos > MAX_INTERP_MONOS {
        return Err(format!("interpolate: {} monomials in one entry degree exceed {} ({} x {} flows; refused before any SBim build)", max_monos, MAX_INTERP_MONOS, nt, ns));
    }
    if cells > MAX_INTERP_CELLS {
        return Err(format!("interpolate: {} localization entries ({} points x {} x {} flows) exceed {} (refused before any SBim build)", cells, max_monos + 4, nt, ns, MAX_INTERP_CELLS));
    }
    Ok(c)
}

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

fn positions(fps: &[crate::soergel::local::FixedPoint], fi: &FlowIndex) -> Vec<usize> {
    fps.iter().map(|p| fi.locate(&p.moved).expect("fixed point is not a flow") as usize).collect()
}

/// The exact PMat (over the window bottom ring, fibre bases of `sbim`) of the bimodule map whose
/// localization is `f`: Φ(y) = L_tgt(y)⁻¹ f(y) L_src(y), every entry interpolated as a polynomial
/// of its forced degree from (#monomials) points and checked at 2 more. Err if an entry is not
/// polynomial (the closed form is not a bimodule map of degree `deg`) or the window is too large.
pub fn interpolate<F: Field>(src: &Window, tgt: &Window, deg: i32, f: &mut dyn FnMut(&[F]) -> WinLoc<F>, seed: u64) -> Result<PMat<F>, String> {
    interp_gate(src, tgt, deg)?;
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
    let off = crate::soergel::local::root_offsets(&src.labels.iter().map(|&a| a as usize).collect::<Vec<_>>());
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
            evals.push(src.labels.iter().enumerate().map(|(i, &l)| esym(crate::combinat::coset::low_mask(l as usize) << off[i], y)).collect());
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
    interp_gate(&Window::new(&ws.bottom, &ws.rungs), &Window::new(&wt.bottom, &wt.rungs), 0).map_err(|e| format!("validate_assoc: {}", e))?;
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
    // review wave 4 finding 1: decide the middle step's route before building anything, and
    // refuse (Err) when neither the numeric Assoc certificate nor the Pull fallback is within
    // the interpolation budget (e.g. cells (3,2,3): 560 flows, 39 148 monomials).
    {
        let (wm, wmp) = (win(&lam, &bw.m), win(&lam, &bw.mp));
        let fl = Expr::new(&lam, &bw.m)?.nflows().max(Expr::new(&lam, &bw.mp)?.nflows());
        let assoc_ok = fl <= ASSOC_CERT_FLOWS as u128 && interp_gate(&wm, &wmp, 0).is_ok();
        if !assoc_ok {
            if let Err(e) = interp_gate(&wm, &wmp, 0) {
                return Err(format!("braid_dag {:?}: middle step M -> M' over budget for the Assoc certificate and the Pull fallback: {}", bw.cells, e));
            }
        }
        interp_gate(&psi_a.src, &psi_a.tgt, psi_a.deg).map_err(|e| format!("braid psi_A: {}", e))?;
        interp_gate(&phi_b.src, &phi_b.tgt, phi_b.deg).map_err(|e| format!("braid phi_B: {}", e))?;
    }
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

// =============================================================================================
// M4c B2: theorem-level Assoc certificate for flag-type words (user decision 2026-11-12)
// =============================================================================================
//
// **Statement used (SLL §4.2, eq. (88)).** If two words ws, wt on the same window have the same
// EKLP multistep expression and that expression is a single up-down [[I ⊂ W ⊃ J]] (flag type),
// then BS(ws) ≅ BS(wt) ≅ R^I ⊗_{R^W} R^J, both generated by 1⊗1, and the isomorphism with
// 1⊗1 ↦ 1⊗1 localizes to the identity matching of flows (one flow per top coset on both sides,
// with equal top configurations). So the E4 `Assoc` closed form with scalar 1 is a genuine
// bimodule map, and it is accepted **without a numeric certificate**.
//
// Conditions checked (all combinatorial, microseconds):
//   (1) equal window bottoms and tops;
//   (2) equal multistep expressions (`multistep`), of up-down shape: after the reductions of
//       `multistep` the sequence has ≤ 3 terms and, with 3 terms, the middle wall set is
//       contained in both ends (an up step followed by a down step);
//   (3) both words flag-type (≤ 1 flow per top coset), equal supports (flow bijection with equal
//       top configurations), equal normalization shifts (degree 0).
//
// Cross-checks (the conditions of the approval): `theorem_vs_numeric` compares the verdict with
// `validate_assoc`'s polynomiality certificate wherever that is affordable; `spot_check_assoc`
// tests the bimodule property of each (3,3,3) use at random points without the GrTable
// multiplication tables (light localization: only the fibre bases are needed).
//
// **Marking.** Every window pair certified this way is entered in a process-wide registry
// (`thm_assoc_registered`); Elem nodes added through `add_thm_assoc` carry `THM_ASSOC_TAG` in
// their provenance; `relies_on_theorem` walks a DAG and reports the theorem-certified Assoc
// nodes it uses; `smapio::export` writes a `thm-assoc <node>` record for each, and
// `smapio::import` re-runs the theorem check on those records (and only on those).

use crate::colored::sexpr::Expr as SExpr;
use std::collections::HashSet;
use std::sync::{Mutex, OnceLock};

/// Provenance tag (and report label) of maps that rely on the theorem-level Assoc certificate.
pub const THM_ASSOC_TAG: &str = "[thm-assoc SLL(88)]";

/// Result of the theorem-level check.
#[derive(Clone, Debug)]
pub struct ThmAssoc {
    pub multistep: Vec<Vec<usize>>,
    pub flows: usize,
}

fn thm_registry() -> &'static Mutex<HashSet<(SExpr, SExpr)>> {
    static REG: OnceLock<Mutex<HashSet<(SExpr, SExpr)>>> = OnceLock::new();
    REG.get_or_init(|| Mutex::new(HashSet::new()))
}

/// Is the (src window, tgt window) pair certified by the theorem (in this process)?
pub fn thm_assoc_registered(ws: &SExpr, wt: &SExpr) -> bool {
    thm_registry().lock().unwrap().contains(&(ws.clone(), wt.clone()))
}

/// Is a multistep expression a single up-down [[I ⊂ W ⊃ J]] (or a monotone part of one)?
pub fn is_updown(ms: &[Vec<usize>]) -> bool {
    match ms.len() {
        0..=2 => true,
        3 => subset(&ms[1], &ms[0]) && subset(&ms[1], &ms[2]),
        _ => false,
    }
}

/// The theorem-level check (conditions (1)-(3) above). Nothing is recorded.
pub fn theorem_assoc_check(ws: &SExpr, wt: &SExpr) -> Result<ThmAssoc, String> {
    if ws.bottom != wt.bottom {
        return Err("thm-assoc: window bottoms differ".into());
    }
    if ws.top() != wt.top() {
        return Err("thm-assoc: window tops differ".into());
    }
    let (ma, mb) = (multistep(&ws.bottom, &ws.rungs)?, multistep(&wt.bottom, &wt.rungs)?);
    if ma != mb {
        return Err(format!("thm-assoc: multistep expressions differ: {:?} vs {:?}", ma, mb));
    }
    if !is_updown(&ma) {
        return Err(format!("thm-assoc: multistep expression {:?} is not a single up-down [[I ⊂ W ⊃ J]]", ma));
    }
    let (a, b) = (FlowIndex::build(ws), FlowIndex::build(wt));
    if !is_flag_type(&a) || !is_flag_type(&b) {
        return Err("thm-assoc: a word is not flag-type (more than one flow at a top coset)".into());
    }
    if a.support() != b.support() {
        return Err("thm-assoc: supports differ (no flow bijection with equal top configurations)".into());
    }
    if ws.shift() != wt.shift() {
        return Err("thm-assoc: shifts differ".into());
    }
    Ok(ThmAssoc { multistep: ma, flows: a.nflows() })
}

/// Theorem-level certificate: check, record the constant 1 for the window pair, register it.
pub fn certify_assoc_theorem<F: Field>(ar: &mut Arena<F>, ws: &SExpr, wt: &SExpr) -> Result<ThmAssoc, String> {
    let r = theorem_assoc_check(ws, wt)?;
    ar.certify_assoc(ws, wt, F::one());
    thm_registry().lock().unwrap().insert((ws.clone(), wt.clone()));
    Ok(r)
}

/// Add an E4 Assoc node src -> tgt (window at rung `at`, upright `upright`) certified by the
/// theorem; its provenance carries `THM_ASSOC_TAG`. Err (and nothing certified) if the theorem's
/// conditions fail. If the window pair already has a numerically validated constant, it is used
/// (the node is then not marked).
pub fn add_thm_assoc<F: Field>(ar: &mut Arena<F>, src: crate::colored::smap::ExprId, tgt: crate::colored::smap::ExprId, at: usize, upright: usize, what: &str) -> Result<MapId, String> {
    let m = ar.add(MapNode::Elem { kind: ElemKind::Assoc, src, tgt, at, upright, scalar: F::one() }, &format!("{} {}", what, THM_ASSOC_TAG))?;
    let (ws, wt) = {
        let (a, b) = ar.elem_windows(m).ok_or("add_thm_assoc: no window record")?;
        (a.clone(), b.clone())
    };
    if ar.assoc_certified(&ws, &wt).is_some() && !thm_assoc_registered(&ws, &wt) {
        return Ok(m);
    }
    certify_assoc_theorem(ar, &ws, &wt)?;
    Ok(m)
}

/// The Assoc nodes of the DAG of m whose constant comes from the theorem (marker for reports).
pub fn relies_on_theorem<F: Field>(ar: &Arena<F>, m: MapId) -> Vec<MapId> {
    let mut out = vec![];
    let mut stack = vec![m];
    let mut seen = HashSet::new();
    while let Some(x) = stack.pop() {
        if !seen.insert(x) {
            continue;
        }
        match ar.node(x) {
            MapNode::Elem { kind: ElemKind::Assoc, .. } => {
                if let Some((ws, wt)) = ar.elem_windows(x) {
                    if thm_assoc_registered(ws, wt) {
                        out.push(x);
                    }
                }
            }
            MapNode::Comp(g, f) => {
                stack.push(*g);
                stack.push(*f);
            }
            MapNode::Lin(t) => stack.extend(t.iter().map(|x| x.1)),
            MapNode::Kron(f, _) => stack.push(*f),
            _ => {}
        }
    }
    out.sort_unstable();
    out
}

/// Cross-check (i): theorem verdict against the numeric certificate on one window pair.
#[derive(Clone, Debug)]
pub struct AssocCross {
    pub src: Vec<Rung>,
    pub tgt: Vec<Rung>,
    pub labels: Vec<u8>,
    pub flows: usize,
    /// Ok = the theorem accepts
    pub theorem: Result<(), String>,
    /// Some(true/false) = the polynomiality certificate of the identity matching passed/failed;
    /// None = not run (flows above `cap` or no flow bijection)
    pub numeric: Option<bool>,
    pub numeric_secs: f64,
}
impl AssocCross {
    /// agreement: theorem accepts ⇒ numeric passes; numeric fails ⇒ theorem refuses
    pub fn consistent(&self) -> bool {
        !(self.theorem.is_ok() && self.numeric == Some(false))
    }
}

/// Run the theorem check and (if both windows biject and have ≤ cap flows) the numeric
/// polynomiality certificate of the identity matching on the same pair.
pub fn theorem_vs_numeric<F: Field>(labels: &[u8], a: &[Rung], b: &[Rung], cap: usize, seed: u64) -> Result<AssocCross, String> {
    let (ws, wt) = (SExpr::new(labels, a)?, SExpr::new(labels, b)?);
    let theorem = theorem_assoc_check(&ws, &wt).map(|_| ());
    let (fa, fb) = (FlowIndex::build(&ws), FlowIndex::build(&wt));
    let flows = fa.nflows().max(fb.nflows());
    let biject = fa.support() == fb.support() && fa.max_block() == 1 && fb.max_block() == 1 && ws.shift() == wt.shift();
    let (mut numeric, mut secs) = (None, 0.0);
    if biject && flows <= cap && interp_gate(&Window::new(labels, a), &Window::new(labels, b), 0).is_ok() {
        let t0 = std::time::Instant::now();
        let r = poly_certificate::<F>(&Window::new(labels, a), &Window::new(labels, b), 0, &mut |y: &[F]| elem_winloc(ElemKind::Assoc, &fa, &fb, y), seed);
        secs = t0.elapsed().as_secs_f64();
        numeric = Some(r.is_ok());
    }
    Ok(AssocCross { src: a.to_vec(), tgt: b.to_vec(), labels: labels.to_vec(), flows, theorem, numeric, numeric_secs: secs })
}

// ------------------------------------------------------------------ light localization

/// Localization data of a window word without the GrTable multiplication tables: per sorted
/// flow, the row of L (values of the fibre basis monomials, mixed radix over rungs as `sbim`),
/// and the fibre degrees `fdeg` (as `SBim::fdeg`) and the shift.
pub struct LightLoc {
    pub fi: FlowIndex,
    pub fdeg: Vec<i32>,
    pub shift: i32,
    bases: Vec<Vec<Vec<u16>>>,
}

impl LightLoc {
    pub fn new(e: &SExpr) -> LightLoc {
        let fi = FlowIndex::build(e);
        let mut labels: Vec<usize> = e.bottom.iter().map(|&a| a as usize).collect();
        let mut fdeg = vec![0i32];
        let mut bases = vec![];
        let mut shift = 0i32;
        for rg in &e.rungs {
            let (rec, giv) = rung_ends(rg);
            let (r, b) = (rg.r, labels[giv]);
            let (basis, wt) = crate::algebra::symfun::fibre_basis(r, b);
            let mut nf = Vec::with_capacity(fdeg.len() * basis.len());
            for &d in &fdeg {
                for &w in &wt {
                    nf.push(d + 2 * w as i32);
                }
            }
            fdeg = nf;
            shift -= (r * (b - r)) as i32;
            bases.push(basis);
            labels[giv] -= r;
            labels[rec] += r;
        }
        LightLoc { fi, fdeg, shift, bases }
    }
    pub fn n(&self) -> usize {
        self.fdeg.len()
    }
    /// L(y) in sorted flow order (rows = flows, columns = fibre basis)
    pub fn matrix<F: Field>(&self, y: &[F]) -> Vec<Vec<F>> {
        let n = self.n();
        assert_eq!(n, self.fi.nflows(), "LightLoc: fibre rank != flows");
        let fs: Vec<usize> = self.bases.iter().map(|b| b.len()).collect();
        (0..n as u32)
            .map(|f| {
                let mv = self.fi.flow(f);
                let vals: Vec<Vec<F>> = self
                    .bases
                    .iter()
                    .zip(mv.iter())
                    .map(|(basis, &c)| {
                        let e = esym(c, y);
                        basis
                            .iter()
                            .map(|al| {
                                let mut v = F::one();
                                for (j, &k) in al.iter().enumerate() {
                                    for _ in 0..k {
                                        v = v.mul(&e[j + 1]);
                                    }
                                }
                                v
                            })
                            .collect()
                    })
                    .collect();
                (0..n)
                    .map(|beta| {
                        let mut rem = beta;
                        let mut v = F::one();
                        for k in (0..fs.len()).rev() {
                            let a = rem % fs[k];
                            rem /= fs[k];
                            v = v.mul(&vals[k][a]);
                        }
                        v
                    })
                    .collect()
            })
            .collect()
    }
}

/// Solve A X = B (A n×n, B n×k) by Gaussian elimination; None if A is singular.
pub fn solve_multi<F: Field>(mut a: Vec<Vec<F>>, mut b: Vec<Vec<F>>) -> Option<Vec<Vec<F>>> {
    let n = a.len();
    let k = if n == 0 { 0 } else { b[0].len() };
    for col in 0..n {
        let piv = (col..n).find(|&r| !a[r][col].is_zero())?;
        a.swap(col, piv);
        b.swap(col, piv);
        let inv = a[col][col].inv();
        for j in col..n {
            a[col][j] = a[col][j].mul(&inv);
        }
        for j in 0..k {
            b[col][j] = b[col][j].mul(&inv);
        }
        let (prow_a, prow_b) = (a[col].clone(), b[col].clone());
        for r in 0..n {
            if r == col || a[r][col].is_zero() {
                continue;
            }
            let f = a[r][col].neg();
            for j in col..n {
                let x = prow_a[j].clone();
                a[r][j].fma(&f, &x);
            }
            for j in 0..k {
                let x = prow_b[j].clone();
                b[r][j].fma(&f, &x);
            }
        }
    }
    Some(b)
}

/// Report of a bimodule-property spot check.
#[derive(Clone, Debug)]
pub struct SpotReport {
    pub flows: usize,
    pub cols: Vec<usize>,
    pub cosets_checked: usize,
    pub points: usize,
    pub entries: usize,
    pub max_deg: i32,
    pub secs: f64,
}

/// Bimodule-property spot check of a degree-`deg` closed form `f` (WinLoc in sorted flow order)
/// between window words a -> b, at random points/cosets (approval condition (ii)):
/// (a) commutation with the localized action of the top generators at `ncos` random top cosets
///     (f maps the flows of a coset into the flows of the same coset; at a random point
///     f·D_src = D_tgt·f for the diagonal top-variable action D);
/// (b) polynomiality of Φ = L_tgt⁻¹ f L_src on `ncols` sampled fibre columns (column 0 plus
///     random ones of fibre degree ≤ `dcap`): along a random line y(τ) = a + τ·b every sampled
///     entry of degree D is a polynomial in τ of degree ≤ D/2 (D < 0: zero), tested with
///     D/2 + 2 points. Needs only the fibre bases (no GrTable multiplication tables).
pub fn spot_check<F: Field>(a: &SExpr, b: &SExpr, deg: i32, f: &mut dyn FnMut(&[F]) -> WinLoc<F>, ncos: usize, ncols: usize, dcap: i32, seed: u64) -> Result<SpotReport, String> {
    use rand::{rngs::StdRng, Rng, SeedableRng};
    let t0 = std::time::Instant::now();
    let (ls, lt) = (LightLoc::new(a), LightLoc::new(b));
    let (ns, nt) = (ls.n(), lt.n());
    let kw: usize = a.bottom.iter().map(|&x| x as usize).sum();
    let mut rng = StdRng::seed_from_u64(seed);
    // (a) cosets and commutation with the diagonal top action
    let y0 = points::<F>(kw, seed ^ 0x5b07, 1, 0).remove(0).x;
    let wl = f(&y0);
    if wl.nrows as usize != nt || wl.ncols() != ns {
        return Err(format!("spot check: WinLoc shape {}x{} != {}x{}", wl.nrows, wl.ncols(), nt, ns));
    }
    let supp = ls.fi.support().to_vec();
    let mut cos_checked = 0;
    for _ in 0..ncos.min(supp.len()) {
        let c = supp[rng.gen_range(0..supp.len())];
        for j in ls.fi.range(c) {
            let (s, e) = (wl.col_start[j as usize] as usize, wl.col_start[j as usize + 1] as usize);
            for k in s..e {
                let i = wl.rows[k];
                if lt.fi.coset_of_flow(i) != c && !wl.vals[k].is_zero() {
                    return Err(format!("spot check: flow {} of coset {} maps into coset {} (not right-linear)", j, c, lt.fi.coset_of_flow(i)));
                }
            }
        }
        cos_checked += 1;
    }
    // (b) polynomiality along a line on sampled columns
    let u = deg - lt.shift + ls.shift;
    let mut cols: Vec<usize> = vec![0];
    let low: Vec<usize> = (1..ns).filter(|&j| ls.fdeg[j] <= dcap).collect();
    let mut tries = 0;
    while cols.len() < ncols.min(low.len() + 1) && tries < 1000 {
        let j = low[rng.gen_range(0..low.len())];
        if !cols.contains(&j) {
            cols.push(j);
        }
        tries += 1;
    }
    let mind_t = *lt.fdeg.iter().min().unwrap_or(&0);
    let max_deg = cols.iter().map(|&j| u + ls.fdeg[j] - mind_t).max().unwrap_or(0).max(0);
    let npts = (max_deg / 2) as usize + 2;
    let base = points::<F>(kw, seed ^ 0x11e0, 2, 0);
    let (pa, pb) = (&base[0].x, &base[1].x);
    let mut taus: Vec<F> = vec![];
    let mut vals: Vec<Vec<Vec<F>>> = vec![]; // [point][row][sampled col]
    let mut p = 0i64;
    while taus.len() < npts {
        p += 1;
        let tau = F::from_i64(p * 7919 + 13);
        let y: Vec<F> = pa.iter().zip(pb).map(|(x, d)| x.add(&tau.mul(d))).collect();
        let fl = f(&y).dense();
        let lsm = ls.matrix(&y);
        // rhs = f · L_src[:, cols]
        let rhs: Vec<Vec<F>> = (0..nt)
            .map(|i| {
                cols.iter()
                    .map(|&j| {
                        let mut s = F::zero();
                        for q in 0..ns {
                            if !fl[i][q].is_zero() {
                                s.fma(&fl[i][q], &lsm[q][j]);
                            }
                        }
                        s
                    })
                    .collect()
            })
            .collect();
        let Some(phi) = solve_multi(lt.matrix(&y), rhs) else { continue };
        taus.push(tau);
        vals.push(phi);
    }
    let mut entries = 0;
    for (ci, &j) in cols.iter().enumerate() {
        for bp in 0..nt {
            let d = u + ls.fdeg[j] - lt.fdeg[bp];
            entries += 1;
            if d < 0 {
                if vals.iter().any(|v| !v[bp][ci].is_zero()) {
                    return Err(format!("spot check: entry ({}, {}) of negative degree {} is nonzero", bp, j, d));
                }
                continue;
            }
            let m = (d / 2) as usize + 1; // interpolation nodes
            if m + 1 > taus.len() {
                continue;
            }
            // Newton divided differences on the first m nodes, test the rest
            let xs = &taus[..];
            let mut coef: Vec<F> = (0..m).map(|k| vals[k][bp][ci].clone()).collect();
            for lvl in 1..m {
                for k in (lvl..m).rev() {
                    let num = coef[k].sub(&coef[k - 1]);
                    let den = xs[k].sub(&xs[k - lvl]);
                    coef[k] = num.mul(&den.inv());
                }
            }
            for k in m..taus.len() {
                let x = &xs[k];
                let mut v = coef[m - 1].clone();
                for l in (0..m - 1).rev() {
                    v = v.mul(&x.sub(&xs[l])).add(&coef[l]);
                }
                if v != vals[k][bp][ci] {
                    return Err(format!("spot check: entry ({}, {}) of degree {} is not polynomial along a line", bp, j, d));
                }
            }
        }
    }
    Ok(SpotReport { flows: ns.max(nt), cols, cosets_checked: cos_checked, points: taus.len(), entries, max_deg, secs: t0.elapsed().as_secs_f64() })
}

/// Spot check (ii) of a theorem-certified Assoc pair (identity matching).
pub fn spot_check_assoc<F: Field>(ws: &SExpr, wt: &SExpr, ncos: usize, ncols: usize, dcap: i32, seed: u64) -> Result<SpotReport, String> {
    let (a, b) = (FlowIndex::build(ws), FlowIndex::build(wt));
    spot_check::<F>(ws, wt, 0, &mut |y: &[F]| elem_winloc(ElemKind::Assoc, &a, &b, y), ncos, ncols, dcap, seed)
}

/// The braid move A → B (CORE Lemma 3.12) with the middle step M → M' certified by the
/// theorem-level Assoc certificate (`add_thm_assoc`), never numerically: Local(φ_B) ∘ Assoc ∘
/// Local(ψ_A). The result relies on the theorem (`relies_on_theorem` lists the middle node).
pub fn braid_dag_thm<F: Field>(ar: &mut Arena<F>, bw: &BraidWords, seed: u64) -> Result<(MapId, Vec<String>), String> {
    let lam = vec![bw.cells.0, bw.cells.1, bw.cells.2];
    let mut log = vec![];
    let (at_a, la, long_a, short_a) = bw.win_a.clone();
    let (at_b, lb, long_b, short_b) = bw.win_b.clone();
    let psi_a = ClosedForm::new("braid psi_A", CfKind::Rex, win(&la, &long_a), win(&la, &short_a))?;
    let ta = build_in(ar, &psi_a, seed)?;
    log.push(format!("psi_A {:?}->{:?} flows {:?} certified", long_a, short_a, ta.flows));
    let phi_b = ClosedForm::new("braid phi_B", CfKind::Pull, win(&lb, &short_b), win(&lb, &long_b))?;
    let tb = build_in(ar, &phi_b, seed + 1)?;
    log.push(format!("phi_B {:?}->{:?} flows {:?} certified", short_b, long_b, tb.flows));
    let ea = ar.expr(SExpr::new(&lam, &bw.a)?);
    let em = ar.expr(SExpr::new(&lam, &bw.m)?);
    let emp = ar.expr(SExpr::new(&lam, &bw.mp)?);
    let eb = ar.expr(SExpr::new(&lam, &bw.b)?);
    let s1 = ar.add(MapNode::Local { src: ea, tgt: em, at: at_a, first: 0, tmpl: ta.tmpl }, "braid psi_A")?;
    let s3 = ar.add(MapNode::Local { src: emp, tgt: eb, at: at_b, first: 0, tmpl: tb.tmpl }, "braid phi_B")?;
    let s2 = add_thm_assoc(ar, em, emp, 0, 0, "braid up-up/down-down")?;
    log.push(format!("M -> M' Assoc by theorem {} ({} flows)", THM_ASSOC_TAG, FlowIndex::build(&SExpr::new(&lam, &bw.m)?).nflows()));
    let c1 = ar.add(MapNode::Comp(s2, s1), "braid")?;
    let c2 = ar.add(MapNode::Comp(s3, c1), &format!("braid {}", THM_ASSOC_TAG))?;
    Ok((c2, log))
}

// ------------------------------------------------------------------ rex chains (B2/B4)

/// How each move of a realized chain was certified.
#[derive(Clone, Debug, Default)]
pub struct ChainReport {
    pub templates: usize,
    pub template_secs: f64,
    pub assoc_numeric: usize,
    pub assoc_theorem: usize,
    pub lines: Vec<String>,
}

/// Realize a rex chain (src/library/rexplan.rs) as a DAG words[0] → words.last() of degree 0:
/// Collapse = Local(Rex ψ template), Expand = Local(Pull φ template), Swap = Local(Swap
/// template), each certified by `build_in` (`certify_added`); Assoc = Elem Assoc, validated
/// numerically (`validate_assoc`) when its window has ≤ `ASSOC_CERT_FLOWS` flows, otherwise by
/// the theorem-level certificate (`add_thm_assoc`, marked). Templates are cached per window pair
/// within the arena's library by the caller's `cache` map.
pub fn realize_chain<F: Field>(ar: &mut Arena<F>, ch: &crate::library::rexplan::Chain, cache: &mut std::collections::HashMap<(crate::library::rexplan::MoveKind, Vec<u8>, Vec<Rung>, Vec<Rung>), TmplId>, seed: u64) -> Result<(MapId, ChainReport), String> {
    realize_chain_opts(ar, ch, cache, seed, false)
}

/// `realize_chain` with `persistable`: numerically certified Assoc moves (whose constants live
/// in the arena only and are not part of layer A, A2 finding) are realized instead as the same
/// identity matching written as a Pull closed-form template on the move's window, certified by
/// `build_in` (`certify_added`; Added templates are re-certified by smapio on import). The
/// theorem-level Assoc (marked, re-checked on import) is unchanged. Used by the production pair
/// route (B4w), whose maps must survive a layer-A save and load.
pub fn realize_chain_opts<F: Field>(ar: &mut Arena<F>, ch: &crate::library::rexplan::Chain, cache: &mut std::collections::HashMap<(crate::library::rexplan::MoveKind, Vec<u8>, Vec<Rung>, Vec<Rung>), TmplId>, seed: u64, persistable: bool) -> Result<(MapId, ChainReport), String> {
    use crate::library::rexplan::MoveKind;
    let lam = ar.lam().to_vec();
    let mut rep = ChainReport::default();
    let e0 = ar.expr(SExpr::new(&lam, &ch.words[0])?);
    let mut cur = ar.add(MapNode::Id(e0), "rex chain start")?;
    for (k, mv) in ch.moves.iter().enumerate() {
        let src = ar.expr(SExpr::new(&lam, &ch.words[k])?);
        let tgt = ar.expr(SExpr::new(&lam, &ch.words[k + 1])?);
        let step = match mv.kind {
            MoveKind::Assoc if persistable && mv.labels.iter().map(|&a| a as u128).sum::<u128>() > 0 && SExpr::new(&mv.labels, &mv.src).map(|e| e.nflows() <= ASSOC_CERT_FLOWS as u128).unwrap_or(false) => {
                let key = (MoveKind::FlagMatch, mv.labels.clone(), mv.src.clone(), mv.tgt.clone());
                let t = match cache.get(&key) {
                    Some(&t) => t,
                    None => {
                        let cf = ClosedForm::new("rex assoc (as Pull template)", CfKind::Pull, win(&mv.labels, &mv.src), win(&mv.labels, &mv.tgt))?;
                        if cf.deg != 0 {
                            return Err(format!("realize_chain: {} has degree {}", mv.fmt(), cf.deg));
                        }
                        let b = build_in(ar, &cf, seed + 200 + k as u64)?;
                        rep.templates += 1;
                        rep.template_secs += b.interp_secs + b.cert.secs;
                        rep.lines.push(format!("{} as Pull template, certified ({:.2} s)", mv.fmt(), b.interp_secs + b.cert.secs));
                        cache.insert(key, b.tmpl);
                        b.tmpl
                    }
                };
                ar.add(MapNode::Local { src, tgt, at: mv.at, first: mv.first, tmpl: t }, "rex chain assoc (template)")?
            }
            MoveKind::Assoc => {
                let m = ar.add(MapNode::Elem { kind: ElemKind::Assoc, src, tgt, at: mv.at, upright: mv.first, scalar: F::one() }, "rex chain assoc")?;
                let (ws, wt) = {
                    let (a, b) = ar.elem_windows(m).ok_or("realize_chain: Assoc without window")?;
                    (a.clone(), b.clone())
                };
                if ar.assoc_certified(&ws, &wt).is_some() {
                    rep.lines.push(format!("{} (cached constant)", mv.fmt()));
                    m
                } else if FlowIndex::build(&ws).nflows() <= ASSOC_CERT_FLOWS {
                    validate_assoc(ar, &ws, &wt, seed + k as u64)?;
                    rep.assoc_numeric += 1;
                    rep.lines.push(format!("{} numeric", mv.fmt()));
                    m
                } else {
                    let mt = add_thm_assoc(ar, src, tgt, mv.at, mv.first, "rex chain assoc")?;
                    rep.assoc_theorem += 1;
                    rep.lines.push(format!("{} {}", mv.fmt(), THM_ASSOC_TAG));
                    mt
                }
            }
            MoveKind::Exact => {
                // exact Hom^0 window: a basis element with 1⊗ ↦ c·1⊗, c ≠ 0, scaled to c = 1
                let (ws, wt) = (win(&mv.labels, &mv.src), win(&mv.labels, &mv.tgt));
                let t0 = std::time::Instant::now();
                let ids = ar.lib().try_hom_basis_cap(&ws, &wt, 0, u128::MAX)?;
                let kw: usize = mv.labels.iter().map(|&a| a as usize).sum();
                let y = points::<F>(kw, seed ^ 0xe8ac, 1, 0).remove(0).x;
                let mut pick: Option<(TmplId, F)> = None;
                for &t in &ids {
                    let d = ar.lib_ref().localize(t, &y).dense();
                    // image of 1⊗ (all ones)
                    let img: Vec<F> = d.iter().map(|row| row.iter().fold(F::zero(), |acc, v| acc.add(v))).collect();
                    if let Some(c) = img.first().cloned() {
                        if !c.is_zero() && img.iter().all(|v| *v == c) {
                            pick = Some((t, c));
                            break;
                        }
                    }
                }
                let (t, c) = pick.ok_or_else(|| format!("realize_chain: no exact basis element of {} maps 1⊗ to a nonzero multiple", mv.fmt()))?;
                rep.lines.push(format!("{} exact ({} basis elements, {:.2} s)", mv.fmt(), ids.len(), t0.elapsed().as_secs_f64()));
                let l = ar.add(MapNode::Local { src, tgt, at: mv.at, first: mv.first, tmpl: t }, "rex chain exact")?;
                if c.is_one() { l } else { ar.add(MapNode::Lin(vec![(c.inv(), l)]), "rex chain exact (normalized)")? }
            }
            kind => {
                let key = (kind, mv.labels.clone(), mv.src.clone(), mv.tgt.clone());
                let t = match cache.get(&key) {
                    Some(&t) => t,
                    None => {
                        let (ws, wt) = (win(&mv.labels, &mv.src), win(&mv.labels, &mv.tgt));
                        let cf = match kind {
                            MoveKind::Collapse => ClosedForm::new("rex collapse", CfKind::Rex, ws, wt)?,
                            MoveKind::Expand => ClosedForm::new("rex expand", CfKind::Pull, ws, wt)?,
                            MoveKind::FlagMatch => ClosedForm::new("rex flag match", CfKind::Pull, ws, wt)?,
                            MoveKind::Swap => {
                                let na = 1;
                                ClosedForm::new("rex swap", CfKind::Swap { na }, ws, wt)?
                            }
                            MoveKind::Assoc | MoveKind::Exact => unreachable!(),
                        };
                        if cf.deg != 0 {
                            return Err(format!("realize_chain: {} has degree {}", mv.fmt(), cf.deg));
                        }
                        let b = build_in(ar, &cf, seed + 100 + k as u64)?;
                        rep.templates += 1;
                        rep.template_secs += b.interp_secs + b.cert.secs;
                        rep.lines.push(format!("{} certified ({:.2} s)", mv.fmt(), b.interp_secs + b.cert.secs));
                        cache.insert(key, b.tmpl);
                        b.tmpl
                    }
                };
                ar.add(MapNode::Local { src, tgt, at: mv.at, first: mv.first, tmpl: t }, "rex chain move")?
            }
        };
        cur = if matches!(ar.node(cur), MapNode::Id(_)) { step } else { ar.add(MapNode::Comp(step, cur), "rex chain")? };
    }
    Ok((cur, rep))
}
