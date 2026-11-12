//! M5a: text export and import of Arena sub-DAGs (layer-A cache of `sindec`), on the public
//! `smap` / `tmpl` API only (SPEC_3C_GL6 rev1 D11: smap.rs and tmpl.rs are frozen).
//!
//! Format (one item per line):
//! ```text
//! krhom-smapio v1 field=<F::name()> lam=<labels> tmpl_alg=<TMPL_ALG>
//! exprs <n>
//! e <rungs>                         (rungs "E0.1,F1.2", "-" for none; bottom = lam)
//! tmpls <n>
//! H <src window key> <tgt window key> <deg> <index> <hash>     (hom_space basis element)
//! R <c> <d> <positive 0|1> <j> <hash>                           (Rickard differential)
//! A <src key> <tgt key> <deg> <nrows> <ncols> <nterms> <hash>   (added template; terms follow)
//! t <col> <row> <mono hex> <coef>
//! nodes <n>
//! I <e> | D <e> <level> <upright> <nu> | L <src> <tgt> <at> <first> <tmpl>
//! C <g> <f> | N <k> (<coef> <m>)* | K <f> <rungs> | M <kind> <src> <tgt> <at> <upright> <coef>
//! roots <m> <ids...>
//! end
//! ```
//! Expressions, templates and nodes are numbered locally (in order of appearance). Templates are
//! referenced by content: a hom_space template by its call (windows, degree, basis index) and a
//! Rickard template by (c, d, sign, j), both re-obtained from the library on import (disk cache
//! or recomputation, which is deterministic) and checked against the stored content hash; an
//! added template (Search window projections) is stored in full and re-certified on import
//! (`Arena::certify_template`). Coefficients are written as signed integers, so only prime
//! fields are supported (`F::from_i64` parses them back).

use crate::field::Field;
use crate::mpoly::Poly;
use crate::pmat::{PMat, PVec};
use crate::sexpr::Expr;
use crate::smap::{Arena, ElemKind, ExprId, MapId, MapNode};
use crate::tmpl::{Template, TmplId, TmplSource, Window, TMPL_ALG};
use crate::web::{Dir, Rung};
use std::collections::HashMap;
use std::sync::atomic::{AtomicBool, Ordering};

/// Replay-only import (libbuild review F3): an `H` template line whose Hom basis is neither in
/// memory nor in the template disk cache is refused instead of recomputed by the exact solver.
static REPLAY_ONLY: AtomicBool = AtomicBool::new(false);

/// Set the process-wide replay-only mode of [`import`] (libbuild `replay` children).
pub fn set_replay_only(on: bool) {
    REPLAY_ONLY.store(on, Ordering::SeqCst);
}
pub fn replay_only() -> bool {
    REPLAY_ONLY.load(Ordering::SeqCst)
}

pub const SMAPIO_VERSION: &str = "krhom-smapio v1";

/// header line of an export for this field and lambda
pub fn header<F: Field>(lam: &[u8]) -> String {
    format!("{} field={} lam={} tmpl_alg={}", SMAPIO_VERSION, F::name(), labels_str(lam), TMPL_ALG)
}

pub fn labels_str(l: &[u8]) -> String {
    l.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")
}

pub fn rung_str(r: &Rung) -> String {
    format!("{}{}.{}", if r.dir == Dir::E { 'E' } else { 'F' }, r.i, r.r)
}

/// "E0.1,F1.2" ("-" for the empty word)
pub fn rungs_str(rs: &[Rung]) -> String {
    if rs.is_empty() {
        return "-".into();
    }
    rs.iter().map(rung_str).collect::<Vec<_>>().join(",")
}

pub fn parse_rung(s: &str) -> Result<Rung, String> {
    let mut ch = s.chars();
    let dir = match ch.next() {
        Some('E') => Dir::E,
        Some('F') => Dir::F,
        _ => return Err(format!("bad rung {:?}", s)),
    };
    let rest: String = ch.collect();
    let (i, r) = rest.split_once('.').ok_or_else(|| format!("bad rung {:?}", s))?;
    Ok(Rung { i: i.parse().map_err(|_| format!("bad rung {:?}", s))?, dir, r: r.parse().map_err(|_| format!("bad rung {:?}", s))? })
}

pub fn parse_rungs(s: &str) -> Result<Vec<Rung>, String> {
    if s == "-" || s.is_empty() {
        return Ok(vec![]);
    }
    s.split(',').map(parse_rung).collect()
}

fn parse_labels(s: &str) -> Result<Vec<u8>, String> {
    if s.is_empty() {
        return Ok(vec![]);
    }
    s.split(',').map(|x| x.parse::<u8>().map_err(|_| format!("bad labels {:?}", s))).collect()
}

/// window key "1,1|E0.1,F0.1" back to a Window
fn parse_window(s: &str) -> Result<Window, String> {
    let (l, r) = s.split_once('|').ok_or_else(|| format!("bad window key {:?}", s))?;
    let labels = parse_labels(l)?;
    let rungs = if r.is_empty() { vec![] } else { r.split(',').map(parse_rung).collect::<Result<Vec<_>, _>>()? };
    Ok(Window { labels, rungs })
}

fn fnv64(h: &mut u64, s: &str) {
    for b in s.bytes() {
        *h ^= b as u64;
        *h = h.wrapping_mul(0x0100_0000_01b3);
    }
}

/// content hash of an exact PMat (shape, then every (col, row, mono, coefficient))
pub fn pmat_hash<F: Field>(m: &PMat<F>) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    fnv64(&mut h, &format!("{}x{};", m.nrows, m.ncols));
    for (j, col) in m.cols.iter().enumerate() {
        for (i, p) in col {
            for (mono, c) in &p.t {
                fnv64(&mut h, &format!("{} {} {:x} {};", j, i, mono, c));
            }
        }
    }
    h
}

/// content hash of a template (windows, degree, PMat)
pub fn template_hash<F: Field>(t: &Template<F>) -> u64 {
    let mut h = pmat_hash(&t.phi);
    fnv64(&mut h, &format!("{}>{}@{}", t.src.key(), t.tgt.key(), t.deg));
    h
}

fn kind_name(k: ElemKind) -> &'static str {
    use ElemKind::*;
    match k {
        Eps => "Eps",
        Eta => "Eta",
        EpsL => "EpsL",
        EtaL => "EtaL",
        Iota => "Iota",
        Trace => "Trace",
        IotaL => "IotaL",
        TraceL => "TraceL",
        Assoc => "Assoc",
    }
}
fn kind_of(s: &str) -> Result<ElemKind, String> {
    use ElemKind::*;
    Ok(match s {
        "Eps" => Eps,
        "Eta" => Eta,
        "EpsL" => EpsL,
        "EtaL" => EtaL,
        "Iota" => Iota,
        "Trace" => Trace,
        "IotaL" => IotaL,
        "TraceL" => TraceL,
        "Assoc" => Assoc,
        _ => return Err(format!("unknown Elem kind {:?}", s)),
    })
}

fn coef_str<F: Field>(c: &F) -> Result<String, String> {
    let s = format!("{}", c);
    if s.parse::<i64>().is_err() {
        return Err(format!("smapio: coefficient {:?} is not an integer (only prime fields are supported)", s));
    }
    Ok(s)
}
fn parse_coef<F: Field>(s: &str) -> Result<F, String> {
    Ok(F::from_i64(s.parse::<i64>().map_err(|_| format!("bad coefficient {:?}", s))?))
}

fn children<F>(n: &MapNode<F>) -> Vec<MapId> {
    match n {
        MapNode::Comp(g, f) => vec![*g, *f],
        MapNode::Lin(t) => t.iter().map(|x| x.1).collect(),
        MapNode::Kron(f, _) => vec![*f],
        _ => vec![],
    }
}

/// Export the sub-DAG reachable from `roots` (children before parents).
pub fn export<F: Field>(ar: &Arena<F>, roots: &[MapId]) -> Result<String, String> {
    // post-order DFS
    let mut order: Vec<MapId> = vec![];
    let mut state: HashMap<MapId, u8> = HashMap::new();
    for &r in roots {
        let mut stack: Vec<(MapId, bool)> = vec![(r, false)];
        while let Some((m, done)) = stack.pop() {
            if done {
                if state.get(&m) != Some(&2) {
                    state.insert(m, 2);
                    order.push(m);
                }
                continue;
            }
            if state.contains_key(&m) {
                continue;
            }
            state.insert(m, 1);
            stack.push((m, true));
            for c in children(ar.node(m)).into_iter().rev() {
                if !state.contains_key(&c) {
                    stack.push((c, false));
                }
            }
        }
    }
    let mut eidx: HashMap<ExprId, usize> = HashMap::new();
    let mut elines: Vec<String> = vec![];
    let mut tidx: HashMap<TmplId, usize> = HashMap::new();
    let mut tlines: Vec<String> = vec![];
    let mut nidx: HashMap<MapId, usize> = HashMap::new();
    let mut nlines: Vec<String> = vec![];
    let mut ex = |e: ExprId, elines: &mut Vec<String>| -> usize {
        if let Some(&i) = eidx.get(&e) {
            return i;
        }
        let i = elines.len();
        elines.push(format!("e {}", rungs_str(&ar.expr_of(e).rungs)));
        eidx.insert(e, i);
        i
    };
    for &m in &order {
        let line = match ar.node(m) {
            MapNode::Id(e) => format!("I {}", ex(*e, &mut elines)),
            MapNode::Dot { e, level, upright, nu } => {
                let nus = if nu.is_empty() { "-".to_string() } else { labels_str(nu) };
                format!("D {} {} {} {}", ex(*e, &mut elines), level, upright, nus)
            }
            MapNode::Elem { kind, src, tgt, at, upright, scalar } => {
                format!("M {} {} {} {} {} {}", kind_name(*kind), ex(*src, &mut elines), ex(*tgt, &mut elines), at, upright, coef_str(scalar)?)
            }
            MapNode::Local { src, tgt, at, first, tmpl } => {
                let ti = match tidx.get(tmpl) {
                    Some(&i) => i,
                    None => {
                        let t = ar.lib_ref().get(*tmpl);
                        let h = template_hash(t);
                        let l = match &t.source {
                            TmplSource::HomSpace { index } => format!("H {} {} {} {} {:016x}", t.src.key(), t.tgt.key(), t.deg, index, h),
                            TmplSource::Rickard { c, d, positive, j } => format!("R {} {} {} {} {:016x}", c, d, *positive as u8, j, h),
                            TmplSource::Added { certified } => {
                                if !certified {
                                    return Err(format!("smapio: template {} is an uncertified added template", tmpl));
                                }
                                let nt: usize = t.phi.cols.iter().map(|c| c.iter().map(|(_, p)| p.t.len()).sum::<usize>()).sum();
                                let mut s = format!("A {} {} {} {} {} {} {:016x}", t.src.key(), t.tgt.key(), t.deg, t.phi.nrows, t.phi.ncols, nt, h);
                                for (j, col) in t.phi.cols.iter().enumerate() {
                                    for (i, p) in col {
                                        for (mono, c) in &p.t {
                                            s.push_str(&format!("\nt {} {} {:x} {}", j, i, mono, coef_str(c)?));
                                        }
                                    }
                                }
                                s
                            }
                        };
                        let i = tlines.len();
                        tlines.push(l);
                        tidx.insert(*tmpl, i);
                        i
                    }
                };
                format!("L {} {} {} {} {}", ex(*src, &mut elines), ex(*tgt, &mut elines), at, first, ti)
            }
            MapNode::Comp(g, f) => format!("C {} {}", nidx[g], nidx[f]),
            MapNode::Lin(terms) => {
                let mut s = format!("N {}", terms.len());
                for (c, x) in terms {
                    s.push_str(&format!(" {} {}", coef_str(c)?, nidx[x]));
                }
                s
            }
            MapNode::Kron(f, rs) => format!("K {} {}", nidx[f], rungs_str(rs)),
        };
        nidx.insert(m, nlines.len());
        nlines.push(line);
    }
    let mut out = String::new();
    out.push_str(&header::<F>(ar.lam()));
    out.push('\n');
    out.push_str(&format!("exprs {}\n", elines.len()));
    for l in &elines {
        out.push_str(l);
        out.push('\n');
    }
    out.push_str(&format!("tmpls {}\n", tlines.len()));
    for l in &tlines {
        out.push_str(l);
        out.push('\n');
    }
    out.push_str(&format!("nodes {}\n", nlines.len()));
    for l in &nlines {
        out.push_str(l);
        out.push('\n');
    }
    let rs: Vec<String> = roots.iter().map(|r| nidx[r].to_string()).collect();
    out.push_str(&format!("roots {} {}\n", roots.len(), rs.join(" ")));
    out.push_str("end\n");
    Ok(out)
}

/// Import an export into `ar` (same lambda and field; the header is checked). Returns the new
/// ids of the roots, in order. Templates are re-obtained from the library and checked against
/// their content hash; added templates are re-certified.
pub fn import<F: Field>(ar: &mut Arena<F>, text: &str) -> Result<Vec<MapId>, String> {
    let mut lines = text.lines();
    let mut next = || lines.next().ok_or_else(|| "smapio: unexpected end of data".to_string());
    let h = next()?;
    let want = header::<F>(ar.lam());
    if h != want {
        return Err(format!("smapio: header mismatch: {:?} != {:?}", h, want));
    }
    let lam = ar.lam().to_vec();
    let count = |l: &str, kw: &str| -> Result<usize, String> { l.strip_prefix(kw).and_then(|x| x.trim().parse().ok()).ok_or_else(|| format!("smapio: expected {:?}, got {:?}", kw, l)) };
    let ne = count(next()?, "exprs ")?;
    let mut exprs: Vec<ExprId> = Vec::with_capacity(ne);
    for _ in 0..ne {
        let l = next()?;
        let w = l.strip_prefix("e ").ok_or_else(|| format!("smapio: bad expr line {:?}", l))?;
        let e = Expr::new(&lam, &parse_rungs(w)?)?;
        exprs.push(ar.expr(e));
    }
    let nt = count(next()?, "tmpls ")?;
    let mut tmpls: Vec<TmplId> = Vec::with_capacity(nt);
    for _ in 0..nt {
        let l = next()?;
        let f: Vec<&str> = l.split_whitespace().collect();
        let hash_of = |s: &str| u64::from_str_radix(s, 16).map_err(|_| format!("smapio: bad hash in {:?}", l));
        let id = match f.first().copied() {
            Some("H") if f.len() == 6 => {
                let (s, t) = (parse_window(f[1])?, parse_window(f[2])?);
                let deg: i32 = f[3].parse().map_err(|_| format!("bad line {:?}", l))?;
                let index: usize = f[4].parse().map_err(|_| format!("bad line {:?}", l))?;
                if replay_only() && !ar.lib_ref().has_hom(&s, &t, deg) {
                    return Err(format!("smapio (replay): hom basis {} -> {} degree {} is not in the template cache; replay never runs the exact solver", s.key(), t.key(), deg));
                }
                let ids = ar.lib().try_hom_basis_cap(&s, &t, deg, u128::MAX)?;
                let id = *ids.get(index).ok_or_else(|| format!("smapio: hom basis {} -> {} degree {} has {} elements, index {}", s.key(), t.key(), deg, ids.len(), index))?;
                if template_hash(ar.lib_ref().get(id)) != hash_of(f[5])? {
                    return Err(format!("smapio: template content changed: hom {} -> {} degree {} index {}", s.key(), t.key(), deg, index));
                }
                id
            }
            Some("R") if f.len() == 6 => {
                let c: u8 = f[1].parse().map_err(|_| format!("bad line {:?}", l))?;
                let d: u8 = f[2].parse().map_err(|_| format!("bad line {:?}", l))?;
                let pos = f[3] == "1";
                let j: usize = f[4].parse().map_err(|_| format!("bad line {:?}", l))?;
                let (_, ids) = ar.lib().rickard(c, d, pos);
                let id = *ids.get(j).ok_or_else(|| format!("smapio: Rickard ({},{}) has no d_{}", c, d, j))?;
                if template_hash(ar.lib_ref().get(id)) != hash_of(f[5])? {
                    return Err(format!("smapio: Rickard template content changed: ({},{}) {} d_{}", c, d, pos, j));
                }
                id
            }
            Some("A") if f.len() == 8 => {
                let (s, t) = (parse_window(f[1])?, parse_window(f[2])?);
                let deg: i32 = f[3].parse().map_err(|_| format!("bad line {:?}", l))?;
                let nrows: usize = f[4].parse().map_err(|_| format!("bad line {:?}", l))?;
                let ncols: usize = f[5].parse().map_err(|_| format!("bad line {:?}", l))?;
                let nterms: usize = f[6].parse().map_err(|_| format!("bad line {:?}", l))?;
                let mut cols: Vec<Vec<(u32, Vec<(u128, F)>)>> = vec![vec![]; ncols];
                for _ in 0..nterms {
                    let tl = next()?;
                    let g: Vec<&str> = tl.split_whitespace().collect();
                    if g.len() != 5 || g[0] != "t" {
                        return Err(format!("smapio: bad term line {:?}", tl));
                    }
                    let j: usize = g[1].parse().map_err(|_| format!("bad term {:?}", tl))?;
                    let i: u32 = g[2].parse().map_err(|_| format!("bad term {:?}", tl))?;
                    let mono = u128::from_str_radix(g[3], 16).map_err(|_| format!("bad term {:?}", tl))?;
                    let c: F = parse_coef(g[4])?;
                    if j >= ncols || i as usize >= nrows {
                        return Err(format!("smapio: term out of range {:?}", tl));
                    }
                    let col = &mut cols[j];
                    if col.last().map(|e| e.0) != Some(i) {
                        col.push((i, vec![]));
                    }
                    col.last_mut().unwrap().1.push((mono, c));
                }
                let cols: Vec<PVec<F>> = cols.into_iter().map(|c| c.into_iter().map(|(i, t)| (i, Poly { t })).collect()).collect();
                let phi = PMat { nrows, ncols, cols };
                let tm = Template::new(s, t, deg, phi, TmplSource::Added { certified: false });
                if template_hash(&tm) != hash_of(f[7])? {
                    return Err("smapio: added template content hash mismatch".into());
                }
                let id = ar.lib().add_template(tm);
                ar.certify_template(id, 0x5eed_5a10).map_err(|e| format!("smapio: added template re-certification failed: {}", e))?;
                id
            }
            _ => return Err(format!("smapio: bad template line {:?}", l)),
        };
        tmpls.push(id);
    }
    let nn = count(next()?, "nodes ")?;
    let mut nodes: Vec<MapId> = Vec::with_capacity(nn);
    let gx = |i: &str, v: &Vec<ExprId>| -> Result<ExprId, String> { i.parse::<usize>().ok().and_then(|k| v.get(k).copied()).ok_or_else(|| format!("smapio: bad expr ref {:?}", i)) };
    for _ in 0..nn {
        let l = next()?;
        let f: Vec<&str> = l.split_whitespace().collect();
        let gn = |i: &str| -> Result<MapId, String> { i.parse::<usize>().ok().and_then(|k| nodes.get(k).copied()).ok_or_else(|| format!("smapio: bad node ref {:?}", i)) };
        let num = |s: &str| -> Result<usize, String> { s.parse::<usize>().map_err(|_| format!("smapio: bad number in {:?}", l)) };
        let node = match (f.first().copied(), f.len()) {
            (Some("I"), 2) => MapNode::Id(gx(f[1], &exprs)?),
            (Some("D"), 5) => MapNode::Dot { e: gx(f[1], &exprs)?, level: num(f[2])?, upright: num(f[3])?, nu: if f[4] == "-" { vec![] } else { parse_labels(f[4])? } },
            (Some("M"), 7) => MapNode::Elem { kind: kind_of(f[1])?, src: gx(f[2], &exprs)?, tgt: gx(f[3], &exprs)?, at: num(f[4])?, upright: num(f[5])?, scalar: parse_coef(f[6])? },
            (Some("L"), 6) => MapNode::Local { src: gx(f[1], &exprs)?, tgt: gx(f[2], &exprs)?, at: num(f[3])?, first: num(f[4])?, tmpl: *tmpls.get(num(f[5])?).ok_or_else(|| format!("smapio: bad template ref in {:?}", l))? },
            (Some("C"), 3) => MapNode::Comp(gn(f[1])?, gn(f[2])?),
            (Some("N"), _) if f.len() >= 2 => {
                let k = num(f[1])?;
                if f.len() != 2 + 2 * k {
                    return Err(format!("smapio: bad Lin line {:?}", l));
                }
                let mut t = vec![];
                for a in 0..k {
                    t.push((parse_coef::<F>(f[2 + 2 * a])?, gn(f[3 + 2 * a])?));
                }
                MapNode::Lin(t)
            }
            (Some("K"), 3) => MapNode::Kron(gn(f[1])?, parse_rungs(f[2])?),
            _ => return Err(format!("smapio: bad node line {:?}", l)),
        };
        nodes.push(ar.add(node, "smapio import")?);
    }
    let l = next()?;
    let f: Vec<&str> = l.split_whitespace().collect();
    if f.first() != Some(&"roots") || f.len() < 2 {
        return Err(format!("smapio: bad roots line {:?}", l));
    }
    let k: usize = f[1].parse().map_err(|_| format!("smapio: bad roots line {:?}", l))?;
    if f.len() != 2 + k {
        return Err(format!("smapio: bad roots line {:?}", l));
    }
    let mut out = vec![];
    for x in &f[2..] {
        out.push(x.parse::<usize>().ok().and_then(|i| nodes.get(i).copied()).ok_or_else(|| format!("smapio: bad root {:?}", x))?);
    }
    if next()? != "end" {
        return Err("smapio: missing end line".into());
    }
    Ok(out)
}
