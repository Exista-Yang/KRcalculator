//! Universal structure constants for relative Grassmannians.
//!
//! Gr(r,b) over the base k[c_1..c_b] is the ring k[e_1..e_r, f_1..f_{b-r}] in which
//! c(t) = e(t) f(t), i.e. c_j = sum_i e_i f_{j-i}. It is free over k[c] with basis
//! { e^alpha : alpha in N^r, |alpha| <= b-r } (e-monomials, unitriangular to Schur
//! polynomials in a box). We tabulate multiplication by e_j and by f_j in this basis,
//! with coefficients polynomials in c (integers).

use crate::dense::rref;
use crate::field::{Field, Fp};
use crate::poly::{monomials_of_weight, MPoly, Mono};
use std::collections::HashMap;
use std::sync::{Mutex, OnceLock};

/// Polynomial in c_1..c_b with integer coefficients: list of (exponent vector, coeff).
pub type CPoly = Vec<(Mono, i64)>;

#[derive(Clone, Debug)]
pub struct GrTable {
    pub r: usize,
    pub b: usize,
    /// basis exponent vectors alpha (length r)
    pub basis: Vec<Mono>,
    /// weight (= half the q-degree) of each basis element
    pub wt: Vec<u32>,
    /// e_mul[j-1][a] = list of (beta index, coefficient polynomial in c)
    pub e_mul: Vec<Vec<Vec<(usize, CPoly)>>>,
    /// f_mul[j-1][a] for j = 1..b-r
    pub f_mul: Vec<Vec<Vec<(usize, CPoly)>>>,
}

fn weight(m: &[u16]) -> u32 {
    m.iter().enumerate().map(|(i, &k)| (i as u32 + 1) * k as u32).sum()
}

impl GrTable {
    pub fn fsize(&self) -> usize {
        self.basis.len()
    }

    pub fn compute(r: usize, b: usize) -> GrTable {
        assert!(r <= b);
        let s = b - r;
        // fibre basis
        let mut basis: Vec<Mono> = Vec::new();
        for tot in 0..=s {
            // all alpha with |alpha| = tot
            let mut cur = vec![0u16; r];
            fn rec(i: usize, left: usize, cur: &mut Vec<u16>, out: &mut Vec<Mono>) {
                if i == cur.len() {
                    if left == 0 {
                        out.push(cur.clone());
                    }
                    return;
                }
                for k in 0..=left {
                    cur[i] = k as u16;
                    rec(i + 1, left - k, cur, out);
                }
                cur[i] = 0;
            }
            if r == 0 {
                if tot == 0 {
                    basis.push(vec![]);
                }
            } else {
                rec(0, tot, &mut cur, &mut basis);
            }
        }
        basis.sort_by_key(|m| (weight(m), m.clone()));
        let wt: Vec<u32> = basis.iter().map(|m| weight(m)).collect();
        let bindex: HashMap<Mono, usize> = basis.iter().cloned().enumerate().map(|(i, m)| (m, i)).collect();
        let _ = bindex;

        // polynomial ring in e_1..e_r, f_1..f_s
        let nv = r + s;
        let vw: Vec<u32> = (1..=r as u32).chain(1..=s as u32).collect();
        let evar = |i: usize| -> MPoly {
            // e_i, i in 0..=r
            if i == 0 { MPoly::one(nv) } else if i <= r { MPoly::var(nv, i - 1) } else { MPoly::zero(nv) }
        };
        let fvar = |i: usize| -> MPoly {
            if i == 0 { MPoly::one(nv) } else if i <= s { MPoly::var(nv, r + i - 1) } else { MPoly::zero(nv) }
        };
        let cpoly: Vec<MPoly> = (0..=b)
            .map(|j| {
                let mut p = MPoly::zero(nv);
                for i in 0..=j {
                    p = p.add(&evar(i).mul(&fvar(j - i)));
                }
                p
            })
            .collect();
        let cw: Vec<u32> = (1..=b as u32).collect();

        // embed alpha (e-exponents) into a monomial of the big ring
        let emono = |a: &Mono| -> Mono {
            let mut m = vec![0u16; nv];
            m[..r].copy_from_slice(a);
            m
        };

        // targets: (kind, j, a) -> polynomial
        let mut targets: Vec<((u8, usize, usize), MPoly)> = Vec::new();
        for j in 1..=r {
            for (ai, a) in basis.iter().enumerate() {
                let mut m = emono(a);
                m[j - 1] += 1;
                targets.push(((0, j, ai), MPoly::monomial(m)));
            }
        }
        for j in 1..=s {
            for (ai, a) in basis.iter().enumerate() {
                let mut m = emono(a);
                m[r + j - 1] += 1;
                targets.push(((1, j, ai), MPoly::monomial(m)));
            }
        }
        // group by weight
        let mut by_w: HashMap<u32, Vec<usize>> = HashMap::new();
        for (ti, (_, p)) in targets.iter().enumerate() {
            let w = weight_big(p.t.keys().next().unwrap(), &vw);
            by_w.entry(w).or_default().push(ti);
        }
        let mut cmemo: HashMap<Mono, MPoly> = HashMap::new();
        let mut results: Vec<Vec<(usize, CPoly)>> = vec![Vec::new(); targets.len()];
        for (&d, tis) in &by_w {
            // columns: pairs (gamma, beta)
            let mut pairs: Vec<(Mono, usize)> = Vec::new();
            for (bi, &bw) in wt.iter().enumerate() {
                if bw > d {
                    continue;
                }
                for g in monomials_of_weight(&cw, d - bw) {
                    pairs.push((g, bi));
                }
            }
            let rows_m = monomials_of_weight(&vw, d);
            assert_eq!(rows_m.len(), pairs.len(), "Gr({},{}) weight {}: not square", r, b, d);
            let ridx: HashMap<Mono, usize> = rows_m.iter().cloned().enumerate().map(|(i, m)| (m, i)).collect();
            let np = pairs.len();
            // Expand the columns once (integer coefficients).
            let mut colpolys: Vec<Vec<(usize, i64)>> = Vec::with_capacity(np);
            for (g, bi) in pairs.iter() {
                let cg = cmemo
                    .entry(g.clone())
                    .or_insert_with(|| {
                        let mut p = MPoly::one(nv);
                        for (j, &k) in g.iter().enumerate() {
                            p = p.mul(&cpoly[j + 1].pow(k as u32));
                        }
                        p
                    })
                    .clone();
                let full = cg.mul(&MPoly::monomial(emono(&basis[*bi])));
                colpolys.push(full.t.iter().map(|(m, &c)| (ridx[m], c)).collect());
            }
            let tpolys: Vec<Vec<(usize, i64)>> =
                tis.iter().map(|&ti| targets[ti].1.t.iter().map(|(m, &c)| (ridx[m], c)).collect()).collect();
            let sol1 = solve_mod::<2147483647>(&colpolys, &tpolys, rows_m.len(), r, b);
            let sol2 = solve_mod::<1000000007>(&colpolys, &tpolys, rows_m.len(), r, b);
            assert_eq!(sol1, sol2, "Gr({},{}): structure constants not small integers", r, b);
            for (k, &ti) in tis.iter().enumerate() {
                let mut res: HashMap<usize, CPoly> = HashMap::new();
                for pc in 0..np {
                    let xi = sol1[k][pc];
                    if xi == 0 {
                        continue;
                    }
                    res.entry(pairs[pc].1).or_default().push((pairs[pc].0.clone(), xi));
                }
                let mut v: Vec<(usize, CPoly)> = res.into_iter().collect();
                v.sort_by_key(|x| x.0);
                for (_, cp) in v.iter_mut() {
                    cp.sort();
                }
                results[ti] = v;
            }
        }
        let mut e_mul = vec![vec![Vec::new(); basis.len()]; r];
        let mut f_mul = vec![vec![Vec::new(); basis.len()]; s];
        for (ti, ((kind, j, ai), _)) in targets.iter().enumerate() {
            if *kind == 0 {
                e_mul[j - 1][*ai] = results[ti].clone();
            } else {
                f_mul[j - 1][*ai] = results[ti].clone();
            }
        }
        GrTable { r, b, basis, wt, e_mul, f_mul }
    }

    /// Cached access. Tables with b ≥ `persist_min_b()` are also persisted on disk
    /// (see [`GrTable::load_or_compute`]); the in-process cache is checked first.
    pub fn get(r: usize, b: usize) -> &'static GrTable {
        static CACHE: OnceLock<Mutex<HashMap<(usize, usize), &'static GrTable>>> = OnceLock::new();
        let m = CACHE.get_or_init(|| Mutex::new(HashMap::new()));
        let mut g = m.lock().unwrap();
        if let Some(t) = g.get(&(r, b)) {
            return t;
        }
        let t: &'static GrTable = Box::leak(Box::new(GrTable::load_or_compute(r, b)));
        g.insert((r, b), t);
        t
    }

    // -----------------------------------------------------------------------------------------
    // disk persistence (M4c B2 engineering decision): the GrTable(r, 9) builds take 40-75 min.
    //
    // File `<dir>/gr_<r>_<b>.bin`, little endian:
    //   magic "KRGRTBL\0" | u32 format version | u32 r | u32 b | u64 payload length |
    //   u64 FNV-1a hash of the payload | payload.
    // The payload is the full table (basis, weights, e_mul, f_mul). `load` verifies the magic,
    // version, (r, b), length, hash, and that the basis/weights equal the recomputed fibre basis
    // (cheap); on any mismatch it returns None and the table is recomputed (and rewritten).
    // Directory: $KRHOM_GRCACHE, default <crate>/target/grcache; KRHOM_GRCACHE=off disables.
    // Persisted only for b ≥ $KRHOM_GRCACHE_MIN_B (default 7): smaller tables cost seconds.

    /// Directory of the on-disk cache, None if disabled.
    pub fn cache_dir() -> Option<std::path::PathBuf> {
        match std::env::var("KRHOM_GRCACHE") {
            Ok(v) if v == "off" || v == "0" => None,
            Ok(v) if !v.is_empty() => Some(std::path::PathBuf::from(v)),
            _ => Some(std::path::PathBuf::from(concat!(env!("CARGO_MANIFEST_DIR"), "/target/grcache"))),
        }
    }

    fn persist_min_b() -> usize {
        std::env::var("KRHOM_GRCACHE_MIN_B").ok().and_then(|v| v.parse().ok()).unwrap_or(7)
    }

    pub fn cache_path(dir: &std::path::Path, r: usize, b: usize) -> std::path::PathBuf {
        dir.join(format!("gr_{}_{}.bin", r, b))
    }

    /// Load from disk (verified), else compute and persist.
    pub fn load_or_compute(r: usize, b: usize) -> GrTable {
        let dir = if b >= Self::persist_min_b() { Self::cache_dir() } else { None };
        if let Some(d) = &dir {
            let p = Self::cache_path(d, r, b);
            if let Some(t) = Self::load(&p, r, b) {
                return t;
            }
        }
        let t = GrTable::compute(r, b);
        if let Some(d) = &dir {
            let _ = t.save(d); // best effort: a failed write only costs a recompute next time
        }
        t
    }

    /// Serialize the table payload (deterministic).
    pub fn payload(&self) -> Vec<u8> {
        let mut o = Vec::new();
        fn u32w(o: &mut Vec<u8>, x: u32) {
            o.extend_from_slice(&x.to_le_bytes());
        }
        fn mono(o: &mut Vec<u8>, m: &[u16]) {
            u32w(o, m.len() as u32);
            for &x in m {
                o.extend_from_slice(&x.to_le_bytes());
            }
        }
        fn tab(o: &mut Vec<u8>, t: &[Vec<Vec<(usize, CPoly)>>]) {
            u32w(o, t.len() as u32);
            for rows in t {
                u32w(o, rows.len() as u32);
                for row in rows {
                    u32w(o, row.len() as u32);
                    for (bi, cp) in row {
                        u32w(o, *bi as u32);
                        u32w(o, cp.len() as u32);
                        for (g, c) in cp {
                            mono(o, g);
                            o.extend_from_slice(&c.to_le_bytes());
                        }
                    }
                }
            }
        }
        u32w(&mut o, self.r as u32);
        u32w(&mut o, self.b as u32);
        u32w(&mut o, self.basis.len() as u32);
        for m in &self.basis {
            mono(&mut o, m);
        }
        for &w in &self.wt {
            u32w(&mut o, w);
        }
        tab(&mut o, &self.e_mul);
        tab(&mut o, &self.f_mul);
        o
    }

    fn from_payload(p: &[u8]) -> Option<GrTable> {
        struct Rd<'a> {
            p: &'a [u8],
            i: usize,
        }
        impl<'a> Rd<'a> {
            fn take(&mut self, n: usize) -> Option<&'a [u8]> {
                let s = self.p.get(self.i..self.i.checked_add(n)?)?;
                self.i += n;
                Some(s)
            }
            fn u32(&mut self) -> Option<u32> {
                Some(u32::from_le_bytes(self.take(4)?.try_into().ok()?))
            }
            fn i64(&mut self) -> Option<i64> {
                Some(i64::from_le_bytes(self.take(8)?.try_into().ok()?))
            }
            fn mono(&mut self) -> Option<Mono> {
                let n = self.u32()? as usize;
                if n > 64 {
                    return None;
                }
                let mut m = Vec::with_capacity(n);
                for _ in 0..n {
                    m.push(u16::from_le_bytes(self.take(2)?.try_into().ok()?));
                }
                Some(m)
            }
            fn tab(&mut self, nb: usize) -> Option<Vec<Vec<Vec<(usize, CPoly)>>>> {
                let nj = self.u32()? as usize;
                let mut t = Vec::with_capacity(nj.min(64));
                for _ in 0..nj {
                    let na = self.u32()? as usize;
                    if na != nb {
                        return None;
                    }
                    let mut rows = Vec::with_capacity(na);
                    for _ in 0..na {
                        let nr = self.u32()? as usize;
                        let mut row = Vec::with_capacity(nr.min(1 << 16));
                        for _ in 0..nr {
                            let bi = self.u32()? as usize;
                            if bi >= nb {
                                return None;
                            }
                            let nc = self.u32()? as usize;
                            let mut cp = Vec::with_capacity(nc.min(1 << 16));
                            for _ in 0..nc {
                                let g = self.mono()?;
                                let c = self.i64()?;
                                cp.push((g, c));
                            }
                            row.push((bi, cp));
                        }
                        rows.push(row);
                    }
                    t.push(rows);
                }
                Some(t)
            }
        }
        let mut rd = Rd { p, i: 0 };
        let r = rd.u32()? as usize;
        let b = rd.u32()? as usize;
        let nb = rd.u32()? as usize;
        if nb > 1 << 24 {
            return None;
        }
        let mut basis = Vec::with_capacity(nb);
        for _ in 0..nb {
            basis.push(rd.mono()?);
        }
        let mut wt = Vec::with_capacity(nb);
        for _ in 0..nb {
            wt.push(rd.u32()?);
        }
        let e_mul = rd.tab(nb)?;
        let f_mul = rd.tab(nb)?;
        if rd.i != p.len() || e_mul.len() != r || f_mul.len() != b.checked_sub(r)? {
            return None;
        }
        Some(GrTable { r, b, basis, wt, e_mul, f_mul })
    }

    /// Write `<dir>/gr_<r>_<b>.bin` atomically (temp file + rename).
    pub fn save(&self, dir: &std::path::Path) -> std::io::Result<std::path::PathBuf> {
        std::fs::create_dir_all(dir)?;
        let pl = self.payload();
        let mut o = Vec::with_capacity(pl.len() + 40);
        o.extend_from_slice(GR_MAGIC);
        o.extend_from_slice(&GR_FORMAT.to_le_bytes());
        o.extend_from_slice(&(self.r as u32).to_le_bytes());
        o.extend_from_slice(&(self.b as u32).to_le_bytes());
        o.extend_from_slice(&(pl.len() as u64).to_le_bytes());
        o.extend_from_slice(&fnv1a64(&pl).to_le_bytes());
        o.extend_from_slice(&pl);
        let path = Self::cache_path(dir, self.r, self.b);
        let tmp = dir.join(format!("gr_{}_{}.bin.tmp{}", self.r, self.b, std::process::id()));
        std::fs::write(&tmp, &o)?;
        std::fs::rename(&tmp, &path)?;
        Ok(path)
    }

    /// Read and verify a persisted table; None on any mismatch (header, hash, basis).
    pub fn load(path: &std::path::Path, r: usize, b: usize) -> Option<GrTable> {
        let d = std::fs::read(path).ok()?;
        if d.len() < 36 || &d[0..8] != GR_MAGIC {
            return None;
        }
        let u32at = |i: usize| u32::from_le_bytes(d[i..i + 4].try_into().unwrap());
        let u64at = |i: usize| u64::from_le_bytes(d[i..i + 8].try_into().unwrap());
        if u32at(8) != GR_FORMAT || u32at(12) as usize != r || u32at(16) as usize != b {
            return None;
        }
        let len = u64at(20) as usize;
        let h = u64at(28);
        let pl = &d[36..];
        if pl.len() != len || fnv1a64(pl) != h {
            return None;
        }
        let t = GrTable::from_payload(pl)?;
        if t.r != r || t.b != b {
            return None;
        }
        // the fibre basis and weights are cheap to recompute: they must agree
        let (basis, wt) = fibre_basis(r, b);
        if t.basis != basis || t.wt != wt {
            return None;
        }
        Some(t)
    }
}

const GR_MAGIC: &[u8; 8] = b"KRGRTBL\0";
/// Bump when `compute` or the payload layout changes.
const GR_FORMAT: u32 = 1;

/// FNV-1a 64-bit hash.
pub fn fnv1a64(d: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &x in d {
        h ^= x as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

/// The fibre basis {e^α : |α| ≤ b − r} sorted by (weight, α), and the weights (as `compute`).
pub fn fibre_basis(r: usize, b: usize) -> (Vec<Mono>, Vec<u32>) {
    let s = b - r;
    let mut basis: Vec<Mono> = Vec::new();
    fn rec(i: usize, left: usize, cur: &mut Vec<u16>, out: &mut Vec<Mono>) {
        if i == cur.len() {
            if left == 0 {
                out.push(cur.clone());
            }
            return;
        }
        for k in 0..=left {
            cur[i] = k as u16;
            rec(i + 1, left - k, cur, out);
        }
        cur[i] = 0;
    }
    for tot in 0..=s {
        if r == 0 {
            if tot == 0 {
                basis.push(vec![]);
            }
        } else {
            let mut cur = vec![0u16; r];
            rec(0, tot, &mut cur, &mut basis);
        }
    }
    basis.sort_by_key(|m| (weight(m), m.clone()));
    let wt = basis.iter().map(|m| weight(m)).collect();
    (basis, wt)
}


/// Solve sum_pc x[pc] * col[pc] = target (for each target) modulo P, returning
/// symmetric representatives.
fn solve_mod<const P: u64>(
    cols: &[Vec<(usize, i64)>],
    targets: &[Vec<(usize, i64)>],
    nrows: usize,
    r: usize,
    b: usize,
) -> Vec<Vec<i64>> {
    let np = cols.len();
    let nt = targets.len();
    let mut mat = vec![vec![Fp::<P>(0); np + nt]; nrows];
    for (pc, c) in cols.iter().enumerate() {
        for &(i, x) in c {
            mat[i][pc] = Fp::<P>::from_i64(x);
        }
    }
    for (k, c) in targets.iter().enumerate() {
        for &(i, x) in c {
            mat[i][np + k] = Fp::<P>::from_i64(x);
        }
    }
    let piv = rref(&mut mat, np + nt);
    assert!(piv.len() == np && piv.iter().enumerate().all(|(i, &p)| i == p), "Gr({},{}): singular", r, b);
    (0..nt).map(|k| (0..np).map(|pc| mat[pc][np + k].signed()).collect()).collect()
}

fn weight_big(m: &[u16], w: &[u32]) -> u32 {
    m.iter().zip(w.iter()).map(|(&k, &x)| k as u32 * x).sum()
}

/// Binomial coefficient.
pub fn binom(n: usize, k: usize) -> usize {
    if k > n {
        return 0;
    }
    let mut r = 1usize;
    for i in 0..k {
        r = r * (n - i) / (i + 1);
    }
    r
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn table_sizes() {
        for b in 0..=6 {
            for r in 0..=b {
                let t = GrTable::get(r, b);
                assert_eq!(t.fsize(), binom(b, r));
            }
        }
    }
    #[test]
    fn persist_roundtrip_and_corruption() {
        let dir = std::env::temp_dir().join(format!("krhom_grcache_test_{}", std::process::id()));
        for (r, b) in [(1usize, 2usize), (2, 4), (3, 5), (0, 3), (4, 4)] {
            let t = GrTable::compute(r, b);
            let p = t.save(&dir).unwrap();
            let u = GrTable::load(&p, r, b).expect("load");
            assert_eq!(u.payload(), t.payload());
            assert_eq!((u.basis.clone(), u.wt.clone()), fibre_basis(r, b));
            assert!(GrTable::load(&p, r + 1, b).is_none(), "wrong (r, b) accepted");
            // flip one payload byte: the hash must reject it
            let mut d = std::fs::read(&p).unwrap();
            let k = d.len() - 3;
            d[k] ^= 0x55;
            std::fs::write(&p, &d).unwrap();
            assert!(GrTable::load(&p, r, b).is_none(), "corrupted file accepted");
            // truncated
            std::fs::write(&p, &d[..d.len() / 2]).unwrap();
            assert!(GrTable::load(&p, r, b).is_none());
        }
        let _ = std::fs::remove_dir_all(&dir);
    }
    #[test]
    fn projective_line() {
        // Gr(1,2): basis 1, e1 ; e1*e1 = c1 e1 - c2
        let t = GrTable::get(1, 2);
        assert_eq!(t.basis, vec![vec![0], vec![1]]);
        let r = &t.e_mul[0][1];
        // expect (0, -c2), (1, c1)
        assert_eq!(r.len(), 2);
        assert_eq!(r[0], (0, vec![(vec![0, 1], -1)]));
        assert_eq!(r[1], (1, vec![(vec![1, 0], 1)]));
        // f1 = c1 - e1
        let f = &t.f_mul[0][0];
        assert_eq!(f, &vec![(0, vec![(vec![1, 0], 1)]), (1, vec![(vec![0, 0], -1)])]);
    }
}
