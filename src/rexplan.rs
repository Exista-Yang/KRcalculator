//! M4c B2: singlestep forms and a constructive rex-chain planner (M4C_PLAN §3.3, milestone B2).
//!
//! Pure combinatorics (no bimodule maps). Given two krhom words A, B that are both reduced
//! expressions of the same double coset q (`dcoset::reduced_for`), the planner finds a chain
//! A = Y_0 → Y_1 → … → Y_k = B of **local rex moves**, each realizable as a certified closed form
//! of B1 (src/leafgen.rs) or as a theorem-level Assoc (B2):
//!
//! | move | window | map (forward) | trust |
//! |---|---|---|---|
//! | `Collapse` | sub reduced for the window maximum → `canonical_max` (flag-type) | Rex ψ (normalized push) | template, `certify_added` |
//! | `Expand` | `canonical_max` → sub reduced for the window maximum | Pull φ | template, `certify_added` |
//! | `Assoc` | flag-type sub → flag-type sub with equal multistep | E4 Assoc, 1⊗1 ↦ 1⊗1 | numeric certificate (≤ 500 flows) else theorem (SLL eq. 88, up-down shape only), marked |
//! | `FlagMatch` | flag-type sub → flag-type sub with equal support and shift | Pull (identity matching) | template, `certify_added` |
//! | `Swap` | adjacent rung blocks on disjoint uprights (≥ 4 uprights only) | Swap | template |
//!
//! Every move has degree 0 and maps 1⊗ ↦ 1⊗ (SLL Thm 6.11), so every composite of moves between
//! reduced words of q has top scalar 1 at q and is an admissible exit/entry map (SLL Prop 4.19,
//! Cor 6.12). The windows are contiguous rung ranges on the minimal range of adjacent uprights
//! the rungs touch (`Local` nodes); a window is "reduced for its maximum" when its word is a
//! reduced expression of the maximal (window bottom, window top) double coset.
//!
//! **Search.** Bidirectional breadth-first search over words reduced for q, from A with forward
//! moves {Collapse, Assoc, Swap} and from B with the same moves (a backward Collapse is an Expand
//! in the chain). This is the constructive content of EKo Thm 5.30 / CORE Thm 3.14 in the form
//! used here: both ends are driven towards flag-type normal forms by collapses (high road), and
//! flag-type pieces are matched by up-up/down-down moves (EKo eq. (88)); switchbacks (EKo eq.
//! (108)) are Collapse/Expand pairs on 2-upright windows, atomic braids (CORE Lemma 3.12) are
//! Collapse · Assoc · Expand on 3 uprights. Every visited word is checked reduced for q
//! (`reduced_for`: c_q(char) = 1 and support ≤ q, gate m9). A search failure is reported, never
//! papered over: it is a census verdict (G-P0), not a proof that no chain exists.
//!
//! **Singlestep forms.** `singlestep` writes a word as EKLP's sequence of ±s steps (L1 §1: a rung
//! E_i^(r)/F_i^(r) is [−s_{c±r} + s_c], trivial steps dropped), the language of the report.
use crate::dcoset::{reduced_for, table_of_perm, top_labels, valid_word, Reduced};
use crate::leafgen::{canonical_max, is_updown, multistep};
use crate::schur::{flows, Schur, Table, Tid};
use crate::sexpr::{rung_ends, Expr, FlowIndex};
use crate::web::{Dir, Rung};
use std::collections::{HashMap, VecDeque};

// ---------------------------------------------------------------------------------------------
// singlestep forms

/// One EKLP step: `Down(c)` adds the wall c (−s_c, split), `Up(c)` removes it (+s_c, merge).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SStep {
    Down(usize),
    Up(usize),
}

fn wall_set(l: &[u8]) -> Vec<usize> {
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

/// The singlestep expression of a krhom word: per rung, the walls added by the split of the
/// giver (moved part next to the receiver) and then the walls removed by the merge.
pub fn singlestep(labels: &[u8], rungs: &[Rung]) -> Result<Vec<SStep>, String> {
    let e = Expr::new(labels, rungs)?;
    let lv = e.levels();
    let mut out = vec![];
    for (t, r) in rungs.iter().enumerate() {
        let (rec, giv) = rung_ends(r);
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
        let (w0, wm, w1) = (wall_set(&lv[t]), wall_set(&nu), wall_set(&lv[t + 1]));
        for &c in &wm {
            if !w0.contains(&c) {
                out.push(SStep::Down(c));
            }
        }
        for &c in &wm {
            if !w1.contains(&c) {
                out.push(SStep::Up(c));
            }
        }
    }
    Ok(out)
}

/// "[I −s3 +s6 …]" with I given by the bottom wall set
pub fn singlestep_str(labels: &[u8], rungs: &[Rung]) -> String {
    match singlestep(labels, rungs) {
        Ok(v) => {
            let mut s = format!("[{:?}", wall_set(labels));
            for st in v {
                match st {
                    SStep::Down(c) => s.push_str(&format!(" −s{}", c)),
                    SStep::Up(c) => s.push_str(&format!(" +s{}", c)),
                }
            }
            s.push(']');
            s
        }
        Err(e) => format!("(invalid: {})", e),
    }
}

// ---------------------------------------------------------------------------------------------
// windows

/// The maximal (rows, cols) double coset: the table of the longest permutation.
pub fn max_table(rows: &[u8], cols: &[u8]) -> Table {
    let k: usize = rows.iter().map(|&a| a as usize).sum();
    let w0: Vec<u8> = (0..k).rev().map(|a| a as u8).collect();
    table_of_perm(&w0, rows, cols)
}

/// A rung window of a word: rungs [at, end) on uprights [first, first + labels.len()).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Win {
    pub at: usize,
    pub end: usize,
    pub first: usize,
    /// window bottom labels (labels of the word's level `at` on the window uprights)
    pub labels: Vec<u8>,
    /// the window's rungs, shifted to start at upright 0
    pub sub: Vec<Rung>,
}

pub fn shift_rungs(rs: &[Rung], d: isize) -> Vec<Rung> {
    rs.iter().map(|r| Rung { i: (r.i as isize + d) as usize, dir: r.dir, r: r.r }).collect()
}

/// Window of the rungs [at, end) of w (bottom λ) on the minimal range of uprights they touch.
pub fn window_of(lam: &[u8], w: &[Rung], at: usize, end: usize) -> Win {
    let lv = top_labels(lam, &w[..at]);
    let u0 = w[at..end].iter().map(|r| r.i).min().unwrap();
    let u1 = w[at..end].iter().map(|r| r.i + 1).max().unwrap();
    Win { at, end, first: u0, labels: lv[u0..=u1].to_vec(), sub: shift_rungs(&w[at..end], -(u0 as isize)) }
}

/// Replace the window of w by `new` (window coordinates).
pub fn splice(w: &[Rung], win: &Win, new: &[Rung]) -> Vec<Rung> {
    let mut out = w[..win.at].to_vec();
    out.extend(shift_rungs(new, win.first as isize));
    out.extend_from_slice(&w[win.end..]);
    out
}

/// Equivariant monomial count of the window bottom ring in degree d (variables e_{i,j} of degree
/// 2j, j = 1..labels_i): the number of interpolation points of `leafgen::interpolate` and
/// `poly_certificate` per entry degree.
pub fn monos_of_deg(labels: &[u8], d: i32) -> u64 {
    if d < 0 || d % 2 != 0 {
        return if d == 0 { 1 } else { 0 };
    }
    let n = (d / 2) as usize;
    let mut c = vec![0u64; n + 1];
    c[0] = 1;
    for &l in labels {
        for j in 1..=l as usize {
            for x in j..=n {
                c[x] = c[x].saturating_add(c[x - j]);
            }
        }
    }
    c[n]
}

/// Maximal entry degree of a degree-0 window map src -> tgt: max fibre degree of src = −2·shift.
pub fn max_entry_deg(labels: &[u8], src: &[Rung]) -> i32 {
    -2 * Expr::new(labels, src).map(|e| e.shift()).unwrap_or(0)
}

/// Predicted seconds of interpolation + `certify_added` of a closed-form template on a window
/// with n_max flows on the larger side and M monomials (fit to the B1 logs, ± factor 2:
/// 200×56 flows / 1 590 monomials 51 s; 225×70 / 1 485: 62 s; 175×70 / 2 114: 90 s;
/// 60×60 / 4 644: 819 s).
pub fn predict_template_secs(n_max: usize, monos: u64) -> f64 {
    let (n, m) = (n_max as f64, monos as f64);
    5e-7 * m * n * n + 7e-9 * m * m * m
}

// ---------------------------------------------------------------------------------------------
// moves and chains

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum MoveKind {
    /// Rex ψ: window word reduced for its maximum → canonical_max
    Collapse,
    /// Pull φ: canonical_max → window word reduced for its maximum
    Expand,
    /// E4 Assoc between flag-type window words with equal multistep expressions
    Assoc,
    /// identity matching between flag-type window words with equal supports and shifts but
    /// different multistep expressions: a Pull template (certified by `certify_added`)
    FlagMatch,
    /// far commutation (disjoint uprights)
    Swap,
    /// exact degree-0 window map between two reduced expressions of the same window coset
    /// (the basis element with top scalar ≠ 0, normalized to 1): `hom_space`, budget-checked
    Exact,
}

/// One move of a chain, in the direction of the chain.
#[derive(Clone, Debug)]
pub struct Move {
    pub kind: MoveKind,
    pub at: usize,
    pub first: usize,
    pub labels: Vec<u8>,
    pub src: Vec<Rung>,
    pub tgt: Vec<Rung>,
    pub src_flows: u128,
    pub tgt_flows: u128,
    /// monomials at the maximal entry degree (template moves)
    pub monos: u64,
    /// Exact moves: unknowns of the degree-0 window call
    pub unknowns: u128,
}

impl Move {
    pub fn max_flows(&self) -> u128 {
        self.src_flows.max(self.tgt_flows)
    }
    /// does this move need a certified template (Collapse, Expand, Swap)?
    pub fn is_template(&self) -> bool {
        matches!(self.kind, MoveKind::Collapse | MoveKind::Expand | MoveKind::Swap | MoveKind::FlagMatch)
    }
    /// Assoc above the numeric certificate's flow limit: relies on the theorem
    pub fn needs_theorem(&self) -> bool {
        self.kind == MoveKind::Assoc && self.max_flows() > crate::leafgen::ASSOC_CERT_FLOWS as u128
    }
    /// Assoc above the numeric limit whose multistep expression is not up-down: neither the
    /// numeric nor the theorem-level certificate applies
    pub fn uncertifiable(&self) -> bool {
        self.needs_theorem() && !multistep(&self.labels, &self.src).map(|m| is_updown(&m)).unwrap_or(false)
    }
    pub fn predicted_secs(&self) -> f64 {
        if self.is_template() { predict_template_secs(self.max_flows() as usize, self.monos) } else { 0.0 }
    }
    pub fn fmt(&self) -> String {
        format!(
            "{:?}@{}/u{} {:?} [{}] -> [{}] ({}x{} flows{})",
            self.kind,
            self.at,
            self.first,
            self.labels,
            crate::smapio::rungs_str(&self.src),
            crate::smapio::rungs_str(&self.tgt),
            self.src_flows,
            self.tgt_flows,
            if self.is_template() { format!(", {} monos", self.monos) } else if self.needs_theorem() { ", theorem".to_string() } else if self.kind == MoveKind::Exact { format!(", {} unknowns", self.unknowns) } else { String::new() }
        )
    }
}

/// A rex chain: words[0] = from, words.last() = to, moves[i]: words[i] -> words[i+1].
#[derive(Clone, Debug)]
pub struct Chain {
    pub words: Vec<Vec<Rung>>,
    pub moves: Vec<Move>,
}

impl Chain {
    pub fn len(&self) -> usize {
        self.moves.len()
    }
    pub fn is_empty(&self) -> bool {
        self.moves.is_empty()
    }
    pub fn max_template_flows(&self) -> u128 {
        self.moves.iter().filter(|m| m.is_template()).map(|m| m.max_flows()).max().unwrap_or(0)
    }
    pub fn max_assoc_flows(&self) -> u128 {
        self.moves.iter().filter(|m| m.kind == MoveKind::Assoc).map(|m| m.max_flows()).max().unwrap_or(0)
    }
    pub fn max_exact_unknowns(&self) -> u128 {
        self.moves.iter().filter(|m| m.kind == MoveKind::Exact).map(|m| m.unknowns).max().unwrap_or(0)
    }
    pub fn theorem_moves(&self) -> usize {
        self.moves.iter().filter(|m| m.needs_theorem()).count()
    }
    pub fn max_word_flows(&self, lam: &[u8]) -> u128 {
        self.words.iter().map(|w| flows(lam, w)).max().unwrap_or(0)
    }
    pub fn predicted_secs(&self) -> f64 {
        // distinct template windows are certified once
        let mut seen = std::collections::HashSet::new();
        self.moves.iter().filter(|m| m.is_template() && seen.insert((m.kind, m.labels.clone(), m.src.clone(), m.tgt.clone()))).map(|m| m.predicted_secs()).sum()
    }
}

/// Planner limits.
#[derive(Clone, Debug)]
pub struct PlanOpts {
    /// longest window subword for Collapse moves
    pub max_collapse: usize,
    /// longest window subword (and enumerated replacement) for Assoc moves
    pub max_assoc: usize,
    /// largest template window (flows on either side) a Collapse/Expand may use
    pub max_template_flows: u128,
    /// largest word (flows) the chain may pass through
    pub max_word_flows: u128,
    /// words longer than this are not visited
    pub max_len: usize,
    /// visited words per side
    pub max_visit: usize,
    /// Exact moves allowed up to this many unknowns (0 = no Exact moves)
    pub exact_budget: u128,
    /// longest enumerated replacement word of an Exact move (2-upright windows; 3-upright: one less)
    pub max_exact: usize,
}

impl Default for PlanOpts {
    fn default() -> Self {
        PlanOpts { max_collapse: 6, max_assoc: 4, max_template_flows: crate::leafgen::MAX_INTERP_FLOWS as u128, max_word_flows: 1_000_000, max_len: 14, max_visit: 200_000, exact_budget: 0, max_exact: 4 }
    }
}

/// Flag-type window word data.
#[derive(Clone, Debug)]
pub struct FlagInfo {
    pub multistep: Vec<Vec<usize>>,
    pub support: Vec<u64>,
    pub shift: i32,
}

/// Memo tables of the planner (window classifications are shared between targets).
#[derive(Default)]
pub struct Planner {
    schurs: HashMap<Vec<u8>, Schur>,
    win_max: HashMap<(Vec<u8>, Vec<Rung>), bool>,
    flag: HashMap<(Vec<u8>, Vec<Rung>), Option<FlagInfo>>,
    /// (labels, top) -> support -> flag-type words with ≤ max_assoc rungs
    assoc_cands: HashMap<(Vec<u8>, Vec<u8>), HashMap<Vec<u64>, Vec<Vec<Rung>>>>,
    /// (labels, top) -> window type (table) -> words reduced for it (≤ max_exact rungs)
    exact_cands: HashMap<(Vec<u8>, Vec<u8>), HashMap<Table, Vec<Vec<Rung>>>>,
    win_type: HashMap<(Vec<u8>, Vec<Rung>), Option<Table>>,
    pub stats: PlannerStats,
}

#[derive(Clone, Debug, Default)]
pub struct PlannerStats {
    pub reduced_checks: usize,
    pub window_checks: usize,
    pub visited: usize,
}

impl Planner {
    pub fn new() -> Planner {
        Planner::default()
    }

    fn schur(&mut self, labels: &[u8]) -> &mut Schur {
        self.schurs.entry(labels.to_vec()).or_insert_with(|| Schur::new(labels))
    }

    /// is the window word a reduced expression of the window's maximal double coset?
    pub fn reduced_for_max(&mut self, labels: &[u8], sub: &[Rung]) -> bool {
        let key = (labels.to_vec(), sub.to_vec());
        if let Some(&b) = self.win_max.get(&key) {
            return b;
        }
        self.stats.window_checks += 1;
        let r = if !valid_word(labels, sub) {
            false
        } else {
            let top = top_labels(labels, sub);
            let t = max_table(labels, &top);
            let s = self.schur(labels);
            let tid = s.intern(t);
            reduced_for(s, sub, tid).ok()
        };
        self.win_max.insert(key, r);
        r
    }

    /// Some(info) if the window word is flag-type (≤ 1 flow per top coset)
    pub fn flag_info(&mut self, labels: &[u8], sub: &[Rung]) -> Option<FlagInfo> {
        let key = (labels.to_vec(), sub.to_vec());
        if let Some(v) = self.flag.get(&key) {
            return v.clone();
        }
        let r = (|| {
            let e = Expr::new(labels, sub).ok()?;
            if flows(labels, sub) > 4_000_000 {
                return None;
            }
            let fi = FlowIndex::build(&e);
            if fi.max_block() != 1 {
                return None;
            }
            Some(FlagInfo { multistep: multistep(labels, sub).ok()?, support: fi.support().to_vec(), shift: e.shift() })
        })();
        self.flag.insert(key, r.clone());
        r
    }

    /// flag-type words (≤ max_assoc rungs) of a window with the given top, by support
    fn assoc_candidates(&mut self, labels: &[u8], top: &[u8], max_assoc: usize) -> &HashMap<Vec<u64>, Vec<Vec<Rung>>> {
        let key = (labels.to_vec(), top.to_vec());
        if !self.assoc_cands.contains_key(&key) {
            let mut by: HashMap<Vec<u64>, Vec<Vec<Rung>>> = HashMap::new();
            let words = all_words(labels, top, max_assoc);
            for w in words {
                if let Some(fi) = self.flag_info(labels, &w) {
                    by.entry(fi.support).or_default().push(w);
                }
            }
            self.assoc_cands.insert(key.clone(), by);
        }
        &self.assoc_cands[&key]
    }

    /// The Assoc move old -> new (words differing inside rungs ≥ at) on smap's window (common
    /// prefix and suffix stripped), if the stripped window words are flag-type with equal
    /// multistep expressions (up-down: numeric or theorem certificate; other shapes, e.g. a
    /// down-up [[I ⊃ W ⊂ J]] or identical singlestep expressions distributed differently over the
    /// rungs: numeric certificate only, `Move::needs_theorem` then reports it as uncertifiable
    /// above 500 flows).
    pub fn assoc_move(&mut self, lam: &[u8], old: &[Rung], new: &[Rung], at: usize, max_template_flows: u128) -> Option<Move> {
        let mut a0 = at;
        while a0 < old.len() && a0 < new.len() && old[a0] == new[a0] {
            a0 += 1;
        }
        let mut sfx = 0;
        while sfx < old.len() - a0 && sfx < new.len() - a0 && old[old.len() - 1 - sfx] == new[new.len() - 1 - sfx] {
            sfx += 1;
        }
        let (a, b) = (&old[a0..old.len() - sfx], &new[a0..new.len() - sfx]);
        if a.is_empty() && b.is_empty() {
            return None;
        }
        let lo = a.iter().chain(b).map(|r| r.i).min()?;
        let hi = a.iter().chain(b).map(|r| r.i + 1).max()?;
        let lv = top_labels(lam, &old[..a0]);
        let labels = lv[lo..=hi].to_vec();
        let (sa, sb) = (shift_rungs(a, -(lo as isize)), shift_rungs(b, -(lo as isize)));
        let ma = self.flag_info(&labels, &sa)?;
        let mb = self.flag_info(&labels, &sb)?;
        if ma.support != mb.support || ma.shift != mb.shift || top_labels(&labels, &sa) != top_labels(&labels, &sb) {
            return None;
        }
        let (fs, ft) = (flows(&labels, &sa), flows(&labels, &sb));
        if ma.multistep == mb.multistep {
            return Some(Move { kind: MoveKind::Assoc, at: a0, first: lo, labels, src: sa, tgt: sb, src_flows: fs, tgt_flows: ft, monos: 0, unknowns: 0 });
        }
        // different multistep expressions: identity matching as a certified Pull template
        if fs.max(ft) > max_template_flows {
            return None;
        }
        let monos = monos_of_deg(&labels, max_entry_deg(&labels, &sa));
        Some(Move { kind: MoveKind::FlagMatch, at: a0, first: lo, labels, src: sa, tgt: sb, src_flows: fs, tgt_flows: ft, monos, unknowns: 0 })
    }

    /// the window type t the window word is a reduced expression of (if any)
    pub fn window_type(&mut self, labels: &[u8], sub: &[Rung]) -> Option<Table> {
        let key = (labels.to_vec(), sub.to_vec());
        if let Some(v) = self.win_type.get(&key) {
            return v.clone();
        }
        let r = (|| {
            if !valid_word(labels, sub) {
                return None;
            }
            let s = self.schur(labels);
            let ch = s.ladder_char(sub);
            let ex = s.expand(&ch);
            for (t, _) in &ex {
                if reduced_for(s, sub, *t).ok() {
                    return Some(s.table(*t).clone());
                }
            }
            None
        })();
        self.win_type.insert(key, r.clone());
        r
    }

    fn exact_candidates(&mut self, labels: &[u8], top: &[u8], maxlen: usize) -> &HashMap<Table, Vec<Vec<Rung>>> {
        let key = (labels.to_vec(), top.to_vec());
        if !self.exact_cands.contains_key(&key) {
            let mut by: HashMap<Table, Vec<Vec<Rung>>> = HashMap::new();
            for w in all_words(labels, top, maxlen) {
                if let Some(t) = self.window_type(labels, &w) {
                    by.entry(t).or_default().push(w);
                }
            }
            self.exact_cands.insert(key.clone(), by);
        }
        &self.exact_cands[&key]
    }

    fn reduced_q(&mut self, s: &mut Schur, w: &[Rung], q: Tid, memo: &mut HashMap<Vec<Rung>, bool>) -> bool {
        if let Some(&b) = memo.get(w) {
            return b;
        }
        self.stats.reduced_checks += 1;
        let b = reduced_for(s, w, q).ok();
        memo.insert(w.to_vec(), b);
        b
    }

    /// forward moves from w (each: the move and the new word), words reduced for q only
    fn neighbours(&mut self, s: &mut Schur, w: &[Rung], q: Tid, opts: &PlanOpts, memo: &mut HashMap<Vec<Rung>, bool>) -> Vec<(Move, Vec<Rung>)> {
        let lam = s.lam.clone();
        let mut out = vec![];
        let n = w.len();
        for at in 0..n {
            for end in at + 1..=n.min(at + opts.max_collapse.max(opts.max_assoc)) {
                let win = window_of(&lam, w, at, end);
                let top = top_labels(&win.labels, &win.sub);
                let len = end - at;
                // Collapse: sub reduced for the window maximum -> canonical_max
                if len <= opts.max_collapse {
                    let cm = canonical_max(&win.labels, &top);
                    if cm != win.sub && self.reduced_for_max(&win.labels, &win.sub) {
                        let (fs, ft) = (flows(&win.labels, &win.sub), flows(&win.labels, &cm));
                        if fs.max(ft) <= opts.max_template_flows {
                            let nw = splice(w, &win, &cm);
                            if nw.len() <= opts.max_len && flows(&lam, &nw) <= opts.max_word_flows && self.reduced_q(s, &nw, q, memo) {
                                let monos = monos_of_deg(&win.labels, max_entry_deg(&win.labels, &win.sub));
                                out.push((Move { kind: MoveKind::Collapse, at, first: win.first, labels: win.labels.clone(), src: win.sub.clone(), tgt: cm, src_flows: fs, tgt_flows: ft, monos, unknowns: 0 }, nw));
                            }
                        }
                    }
                }
                // Assoc: flag-type up-down sub -> another with the same multistep
                if len <= opts.max_assoc.max(opts.max_collapse) {
                    if let Some(fi) = self.flag_info(&win.labels, &win.sub) {
                        let cands: Vec<Vec<Rung>> = self.assoc_candidates(&win.labels, &top, opts.max_assoc).get(&fi.support).cloned().unwrap_or_default();
                        for z in cands {
                            if z == win.sub {
                                continue;
                            }
                            let nw = splice(w, &win, &z);
                            if nw.len() > opts.max_len || flows(&lam, &nw) > opts.max_word_flows {
                                continue;
                            }
                            // the Assoc window as smap sees it: common prefix and suffix stripped
                            let Some(mv) = self.assoc_move(&lam, w, &nw, at, opts.max_template_flows) else { continue };
                            if self.reduced_q(s, &nw, q, memo) {
                                out.push((mv, nw));
                            }
                        }
                    }
                }
            }
        }
        // Exact moves: two reduced expressions of the same window coset, budget-checked
        if opts.exact_budget > 0 {
            for at in 0..n {
                for end in at + 1..=n.min(at + opts.max_exact) {
                    let win = window_of(&lam, w, at, end);
                    let maxlen = if win.labels.len() > 2 { opts.max_exact.saturating_sub(1).max(end - at) } else { opts.max_exact };
                    if end - at > maxlen {
                        continue;
                    }
                    let top = top_labels(&win.labels, &win.sub);
                    let Some(t) = self.window_type(&win.labels, &win.sub) else { continue };
                    let cands: Vec<Vec<Rung>> = self.exact_candidates(&win.labels, &top, maxlen).get(&t).cloned().unwrap_or_default();
                    for z in cands {
                        if z == win.sub {
                            continue;
                        }
                        let u = crate::homprov::hom_unknowns(&win.labels, &win.sub, &z, 0);
                        if u > opts.exact_budget {
                            continue;
                        }
                        let nw = splice(w, &win, &z);
                        if nw.len() > opts.max_len || flows(&lam, &nw) > opts.max_word_flows || !self.reduced_q(s, &nw, q, memo) {
                            continue;
                        }
                        let (fs, ft) = (flows(&win.labels, &win.sub), flows(&win.labels, &z));
                        out.push((Move { kind: MoveKind::Exact, at, first: win.first, labels: win.labels.clone(), src: win.sub.clone(), tgt: z, src_flows: fs, tgt_flows: ft, monos: 0, unknowns: u }, nw));
                    }
                }
            }
        }
        // far commutation of adjacent rungs on disjoint uprights
        for at in 0..n.saturating_sub(1) {
            let (a, b) = (w[at], w[at + 1]);
            let ua = (a.i, a.i + 1);
            let ub = (b.i, b.i + 1);
            if ua.1 < ub.0 || ub.1 < ua.0 {
                let mut nw = w.to_vec();
                nw.swap(at, at + 1);
                if self.reduced_q(s, &nw, q, memo) {
                    let win = window_of(&lam, w, at, at + 2);
                    let sw = shift_rungs(&[b, a], -(win.first as isize));
                    let f = flows(&win.labels, &win.sub);
                    let monos = monos_of_deg(&win.labels, max_entry_deg(&win.labels, &win.sub));
                    out.push((Move { kind: MoveKind::Swap, at, first: win.first, labels: win.labels.clone(), src: win.sub.clone(), tgt: sw, src_flows: f, tgt_flows: f, monos, unknowns: 0 }, nw));
                }
            }
        }
        out
    }

    /// A rex chain from → to between two words reduced for q (bidirectional BFS). Err with the
    /// search statistics if none is found within the limits.
    pub fn plan(&mut self, s: &mut Schur, q: Tid, from: &[Rung], to: &[Rung], opts: &PlanOpts) -> Result<Chain, String> {
        let mut memo: HashMap<Vec<Rung>, bool> = HashMap::new();
        for (nm, w) in [("from", from), ("to", to)] {
            let r: Reduced = reduced_for(s, w, q);
            if !r.ok() {
                return Err(format!("plan: the {} word [{}] is not reduced for q ({:?})", nm, crate::smapio::rungs_str(w), r));
            }
            memo.insert(w.to_vec(), true);
        }
        if from == to {
            return Ok(Chain { words: vec![from.to_vec()], moves: vec![] });
        }
        // side 0 from `from`, side 1 from `to`; parent maps word -> (parent word, move)
        let mut par: [HashMap<Vec<Rung>, Option<(Vec<Rung>, Move)>>; 2] = [HashMap::new(), HashMap::new()];
        let mut qs: [VecDeque<Vec<Rung>>; 2] = [VecDeque::new(), VecDeque::new()];
        par[0].insert(from.to_vec(), None);
        par[1].insert(to.to_vec(), None);
        qs[0].push_back(from.to_vec());
        qs[1].push_back(to.to_vec());
        let mut meet: Option<Vec<Rung>> = None;
        'outer: while !qs[0].is_empty() || !qs[1].is_empty() {
            for side in 0..2 {
                // expand one BFS layer of this side
                let layer: Vec<Vec<Rung>> = qs[side].drain(..).collect();
                for w in layer {
                    if par[side].len() > opts.max_visit {
                        break 'outer;
                    }
                    for (mv, nw) in self.neighbours(s, &w, q, opts, &mut memo) {
                        if par[side].contains_key(&nw) {
                            continue;
                        }
                        par[side].insert(nw.clone(), Some((w.clone(), mv)));
                        self.stats.visited += 1;
                        if par[1 - side].contains_key(&nw) {
                            meet = Some(nw);
                            break 'outer;
                        }
                        qs[side].push_back(nw);
                    }
                }
            }
        }
        let Some(m) = meet else {
            return Err(format!("plan: no chain within limits (visited {} + {} words; max_len {}, collapse ≤ {}, assoc ≤ {}, template ≤ {} flows)", par[0].len(), par[1].len(), opts.max_len, opts.max_collapse, opts.max_assoc, opts.max_template_flows));
        };
        // forward half: from -> m
        let mut fw: Vec<(Vec<Rung>, Move)> = vec![];
        let mut cur = m.clone();
        while let Some(Some((p, mv))) = par[0].get(&cur).cloned() {
            fw.push((p.clone(), mv));
            cur = p;
        }
        fw.reverse();
        let mut words: Vec<Vec<Rung>> = vec![from.to_vec()];
        let mut moves: Vec<Move> = vec![];
        for (_, mv) in fw {
            let nw = splice(words.last().unwrap(), &Win { at: mv.at, end: mv.at + mv.src.len(), first: mv.first, labels: mv.labels.clone(), sub: mv.src.clone() }, &mv.tgt);
            words.push(nw);
            moves.push(mv);
        }
        debug_assert_eq!(words.last().unwrap(), &m);
        // backward half: m -> to, reversing the moves found from `to`
        let mut cur = m;
        while let Some(Some((p, mv))) = par[1].get(&cur).cloned() {
            // mv: p -> cur (found from `to`); in the chain: cur -> p
            let rev = Move {
                kind: match mv.kind {
                    MoveKind::Collapse => MoveKind::Expand,
                    MoveKind::Expand => MoveKind::Collapse,
                    k => k,
                },
                at: mv.at,
                first: mv.first,
                labels: mv.labels.clone(),
                src: mv.tgt.clone(),
                tgt: mv.src.clone(),
                src_flows: mv.tgt_flows,
                tgt_flows: mv.src_flows,
                // the Pull φ template has the same window pair (monomials of the larger side)
                monos: if mv.is_template() { mv.monos.max(monos_of_deg(&mv.labels, max_entry_deg(&mv.labels, &mv.tgt))) } else { 0 },
                unknowns: if mv.kind == MoveKind::Exact { crate::homprov::hom_unknowns(&mv.labels, &mv.tgt, &mv.src, 0) } else { 0 },
            };
            let nw = splice(&cur, &Win { at: rev.at, end: rev.at + rev.src.len(), first: rev.first, labels: rev.labels.clone(), sub: rev.src.clone() }, &rev.tgt);
            debug_assert_eq!(nw, p);
            words.push(nw);
            moves.push(rev);
            cur = p;
        }
        debug_assert_eq!(words.last().unwrap(), &to.to_vec());
        Ok(Chain { words, moves })
    }
}

/// All valid words on a window of `labels` with ≤ maxlen rungs ending at `top` (no rung of size
/// 0; repeated inverse pairs allowed).
pub fn all_words(labels: &[u8], top: &[u8], maxlen: usize) -> Vec<Vec<Rung>> {
    let mut out = vec![];
    let n = labels.len();
    let mut cur: Vec<Rung> = vec![];
    fn rec(l: &mut Vec<u8>, top: &[u8], n: usize, left: usize, cur: &mut Vec<Rung>, out: &mut Vec<Vec<Rung>>) {
        if l.as_slice() == top {
            out.push(cur.clone());
        }
        if left == 0 {
            return;
        }
        // prune: each rung changes two adjacent labels; the L1 distance shrinks by ≤ 2r per rung
        for i in 0..n.saturating_sub(1) {
            for dir in [Dir::E, Dir::F] {
                let (rec_u, giv) = if dir == Dir::E { (i, i + 1) } else { (i + 1, i) };
                for r in 1..=l[giv] as usize {
                    l[giv] -= r as u8;
                    l[rec_u] += r as u8;
                    cur.push(Rung { i, dir, r });
                    rec(l, top, n, left - 1, cur, out);
                    cur.pop();
                    l[giv] += r as u8;
                    l[rec_u] -= r as u8;
                }
            }
        }
    }
    let mut l = labels.to_vec();
    rec(&mut l, top, n, maxlen, &mut cur, &mut out);
    out
}
