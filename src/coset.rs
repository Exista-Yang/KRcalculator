//! Cosets W/W_mu as top configurations of roots (D1 §2.2), tables, witnesses, W_lambda action.
//!
//! Roots are numbered 0..K; bottom block i owns roots off[i]..off[i]+lambda_i (as
//! `local::root_offsets`). A configuration is one root mask per upright; a top coset of the block
//! (lambda, mu) is a configuration with |mask_j| = mu_j. `CosetSpace::rank` is the multinomial
//! number system: upright 0 is the most significant digit, digit j is the colexicographic rank of
//! mask_j among the mu_j-subsets of the roots not used by masks 0..j-1, with radix
//! C(remaining, mu_j).
use crate::schur::Table;
use std::sync::OnceLock;

pub type Mask = u64;

#[derive(Clone, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub struct Conf(pub Vec<Mask>); // one root mask per upright

/// Mask of the n lowest bits (n <= 64).
#[inline]
pub fn low_mask(n: usize) -> Mask {
    if n >= 64 { u64::MAX } else { (1u64 << n) - 1 }
}

/// Binomial coefficients C(n, k) for n, k <= 64 (u128).
pub fn binom(n: usize, k: usize) -> u128 {
    static TAB: OnceLock<Vec<Vec<u128>>> = OnceLock::new();
    let t = TAB.get_or_init(|| {
        let mut t = vec![vec![0u128; 65]; 65];
        for n in 0..=64 {
            t[n][0] = 1;
            for k in 1..=n {
                t[n][k] = t[n - 1][k - 1] + if k < n { t[n - 1][k] } else { 0 };
            }
        }
        t
    });
    if k > n || n > 64 { 0 } else { t[n][k] }
}

#[derive(Clone, Debug)]
pub struct Roots {
    pub lam: Vec<u8>,
    pub off: Vec<usize>,
    pub k: usize,
}

impl Roots {
    /// asserts K <= 64
    pub fn new(lam: &[u8]) -> Roots {
        let mut off = Vec::with_capacity(lam.len());
        let mut k = 0usize;
        for &l in lam {
            off.push(k);
            k += l as usize;
        }
        assert!(k <= 64, "Roots: K = {} > 64", k);
        Roots { lam: lam.to_vec(), off, k }
    }
    /// root mask of bottom block i
    pub fn block_mask(&self, i: usize) -> Mask {
        low_mask(self.lam[i] as usize) << self.off[i]
    }
    /// block i on upright i
    pub fn bottom(&self) -> Conf {
        Conf((0..self.lam.len()).map(|i| self.block_mask(i)).collect())
    }
    pub fn block_of(&self, root: usize) -> usize {
        assert!(root < self.k, "root {} out of range (K = {})", root, self.k);
        (0..self.lam.len()).find(|&i| root >= self.off[i] && root < self.off[i] + self.lam[i] as usize).unwrap()
    }
}

#[derive(Clone, Debug)]
pub struct CosetSpace {
    pub roots: Roots,
    pub mu: Vec<u8>,
} // mu.len() == lam.len(), zeros allowed

impl CosetSpace {
    /// asserts sum mu = K and equal lengths (and that the count fits u64)
    pub fn new(lam: &[u8], mu: &[u8]) -> CosetSpace {
        assert_eq!(lam.len(), mu.len(), "CosetSpace: lambda {:?} and mu {:?} have different lengths", lam, mu);
        let roots = Roots::new(lam);
        let s: usize = mu.iter().map(|&m| m as usize).sum();
        assert_eq!(s, roots.k, "CosetSpace: sum mu {:?} != K of lambda {:?}", mu, lam);
        let cs = CosetSpace { roots, mu: mu.to_vec() };
        let c = cs.count128();
        assert!(c <= u64::MAX as u128, "CosetSpace: count {} does not fit u64", c);
        cs
    }
    fn count128(&self) -> u128 {
        let mut rem = self.roots.k;
        let mut c: u128 = 1;
        for &m in &self.mu {
            c = c.checked_mul(binom(rem, m as usize)).expect("coset count overflow");
            rem -= m as usize;
        }
        c
    }
    /// K! / prod_j mu_j!
    pub fn count(&self) -> u64 {
        self.count128() as u64
    }
    /// radices C(remaining before j, mu_j)
    fn radices(&self) -> Vec<u128> {
        let mut rem = self.roots.k;
        self.mu
            .iter()
            .map(|&m| {
                let r = binom(rem, m as usize);
                rem -= m as usize;
                r
            })
            .collect()
    }
    /// bijection onto 0..count (multinomial number system)
    pub fn rank(&self, c: &Conf) -> u64 {
        debug_assert_eq!(c.0.len(), self.mu.len());
        let full = low_mask(self.roots.k);
        let mut used: Mask = 0;
        let mut rem = self.roots.k;
        let mut r: u128 = 0;
        for (j, &m) in c.0.iter().enumerate() {
            debug_assert_eq!(m.count_ones() as u8, self.mu[j], "rank: mask {:#x} on upright {} has the wrong size", m, j);
            debug_assert_eq!(m & used, 0, "rank: overlapping masks");
            let avail = full & !used;
            let mut cr: u128 = 0;
            let mut bits = m;
            let mut i = 0usize;
            while bits != 0 {
                let b = bits.trailing_zeros();
                bits &= bits - 1;
                i += 1;
                let below = if b == 0 { 0 } else { avail & low_mask(b as usize) };
                cr += binom(below.count_ones() as usize, i);
            }
            let radix = binom(rem, self.mu[j] as usize);
            r = r * radix + cr;
            rem -= self.mu[j] as usize;
            used |= m;
        }
        r as u64
    }
    pub fn unrank(&self, r: u64) -> Conf {
        let rad = self.radices();
        let n = self.mu.len();
        let mut digits = vec![0u128; n];
        let mut x = r as u128;
        for j in (0..n).rev() {
            digits[j] = x % rad[j];
            x /= rad[j];
        }
        assert_eq!(x, 0, "unrank: {} out of range", r);
        let full = low_mask(self.roots.k);
        let mut used: Mask = 0;
        let mut out = Vec::with_capacity(n);
        for j in 0..n {
            let avail = full & !used;
            let roots: Vec<u32> = (0..64).filter(|&b| avail >> b & 1 == 1).collect();
            let k = self.mu[j] as usize;
            let mut d = digits[j];
            let mut mask: Mask = 0;
            for i in (1..=k).rev() {
                // largest p with C(p, i) <= d
                let mut p = i - 1;
                while binom(p + 1, i) <= d {
                    p += 1;
                }
                d -= binom(p, i);
                mask |= 1u64 << roots[p];
            }
            out.push(mask);
            used |= mask;
        }
        Conf(out)
    }
    /// M_ij = popcount(B_i & T_j)
    pub fn table(&self, c: &Conf) -> Table {
        let n = self.mu.len();
        let mut e = vec![0u8; n * n];
        for i in 0..n {
            let b = self.roots.block_mask(i);
            for j in 0..n {
                e[i * n + j] = (b & c.0[j]).count_ones() as u8;
            }
        }
        Table { nr: n as u8, nc: n as u8, e }
    }
    fn check_table(&self, t: &Table) {
        let n = self.mu.len();
        assert!(t.nr as usize == n && t.nc as usize == n, "table {} has the wrong shape for {:?} -> {:?}", t.fmt(), self.roots.lam, self.mu);
        assert_eq!(t.rows(), self.roots.lam, "table {}: row sums != lambda", t.fmt());
        assert_eq!(t.cols(), self.mu, "table {}: column sums != mu", t.fmt());
    }
    /// witness: for j = 0.., for i = 0..: the first M_ij unused roots of B_i
    pub fn rep(&self, t: &Table) -> Conf {
        self.check_table(t);
        let n = self.mu.len();
        let mut next: Vec<usize> = self.roots.off.clone();
        let mut out = vec![0u64; n];
        for j in 0..n {
            for i in 0..n {
                for _ in 0..t.get(i, j) {
                    out[j] |= 1u64 << next[i];
                    next[i] += 1;
                }
            }
        }
        Conf(out)
    }
    /// all N(t) cosets of table t, ascending rank
    pub fn cosets_of(&self, t: &Table) -> Vec<u64> {
        self.check_table(t);
        let n = self.mu.len();
        let mut out = vec![];
        let mut free: Vec<Mask> = (0..n).map(|i| self.roots.block_mask(i)).collect();
        let mut cur = vec![0u64; n];
        fn rec(cs: &CosetSpace, t: &Table, j: usize, i: usize, free: &mut Vec<Mask>, cur: &mut Vec<Mask>, out: &mut Vec<u64>) {
            let n = cs.mu.len();
            if j == n {
                out.push(cs.rank(&Conf(cur.clone())));
                return;
            }
            if i == n {
                rec(cs, t, j + 1, 0, free, cur, out);
                return;
            }
            let k = t.get(i, j) as usize;
            for s in subsets(free[i], k) {
                free[i] &= !s;
                cur[j] |= s;
                rec(cs, t, j, i + 1, free, cur, out);
                cur[j] &= !s;
                free[i] |= s;
            }
        }
        rec(self, t, 0, 0, &mut free, &mut cur, &mut out);
        out.sort_unstable();
        out
    }
    /// sigma: permutation of 0..K that preserves every bottom block
    pub fn act(&self, sigma: &[usize], c: &Conf) -> Conf {
        assert_eq!(sigma.len(), self.roots.k, "act: sigma has length {} != K = {}", sigma.len(), self.roots.k);
        debug_assert!((0..self.roots.k).all(|r| self.roots.block_of(r) == self.roots.block_of(sigma[r])), "act: sigma does not preserve the bottom blocks");
        Conf(c.0.iter().map(|&m| permute_mask(sigma, m)).collect())
    }
}

/// image of a root mask under a permutation of roots
pub fn permute_mask(sigma: &[usize], m: Mask) -> Mask {
    let mut out = 0u64;
    let mut bits = m;
    while bits != 0 {
        let b = bits.trailing_zeros() as usize;
        bits &= bits - 1;
        out |= 1u64 << sigma[b];
    }
    out
}

/// All k-subsets of the bits of `mask` (lexicographic in bit positions).
pub fn subsets(mask: Mask, k: usize) -> Vec<Mask> {
    let bits: Vec<u32> = (0..64).filter(|&b| mask >> b & 1 == 1).collect();
    let nb = bits.len();
    let mut out = vec![];
    if k > nb {
        return out;
    }
    if k == 0 {
        out.push(0);
        return out;
    }
    let mut idx: Vec<usize> = (0..k).collect();
    loop {
        out.push(idx.iter().fold(0u64, |m, &i| m | 1u64 << bits[i]));
        let mut i = k;
        loop {
            if i == 0 {
                return out;
            }
            i -= 1;
            if idx[i] != i + nb - k {
                break;
            }
            if i == 0 {
                return out;
            }
        }
        idx[i] += 1;
        for j in i + 1..k {
            idx[j] = idx[j - 1] + 1;
        }
    }
}
