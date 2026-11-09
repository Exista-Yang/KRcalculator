//! Laurent polynomials in one variable q with integer coefficients, and
//! two-variable Laurent polynomials in (q, t) used for Poincaré polynomials.

use std::collections::BTreeMap;
use std::fmt;

#[derive(Clone, PartialEq, Eq, Default, Hash, PartialOrd, Ord)]
pub struct LPoly {
    /// exponent -> nonzero coefficient
    pub c: BTreeMap<i32, i128>,
}

impl LPoly {
    pub fn zero() -> Self {
        LPoly { c: BTreeMap::new() }
    }
    pub fn one() -> Self {
        Self::mono(0, 1)
    }
    pub fn mono(e: i32, k: i128) -> Self {
        let mut c = BTreeMap::new();
        if k != 0 {
            c.insert(e, k);
        }
        LPoly { c }
    }
    pub fn is_zero(&self) -> bool {
        self.c.is_empty()
    }
    pub fn add_term(&mut self, e: i32, k: i128) {
        if k == 0 {
            return;
        }
        let v = self.c.entry(e).or_insert(0);
        *v = v.checked_add(k).expect("LPoly overflow");
        if *v == 0 {
            self.c.remove(&e);
        }
    }
    pub fn add(&self, o: &LPoly) -> LPoly {
        let mut r = self.clone();
        for (&e, &k) in &o.c {
            r.add_term(e, k);
        }
        r
    }
    pub fn sub(&self, o: &LPoly) -> LPoly {
        let mut r = self.clone();
        for (&e, &k) in &o.c {
            r.add_term(e, -k);
        }
        r
    }
    pub fn scale(&self, k: i128) -> LPoly {
        let mut r = LPoly::zero();
        for (&e, &v) in &self.c {
            r.add_term(e, v.checked_mul(k).expect("LPoly overflow"));
        }
        r
    }
    pub fn shift(&self, s: i32) -> LPoly {
        LPoly { c: self.c.iter().map(|(&e, &k)| (e + s, k)).collect() }
    }
    pub fn mul(&self, o: &LPoly) -> LPoly {
        let mut r = LPoly::zero();
        for (&e1, &k1) in &self.c {
            for (&e2, &k2) in &o.c {
                r.add_term(e1 + e2, k1.checked_mul(k2).expect("LPoly overflow"));
            }
        }
        r
    }
    /// q -> q^{-1}
    pub fn bar(&self) -> LPoly {
        LPoly { c: self.c.iter().map(|(&e, &k)| (-e, k)).collect() }
    }
    pub fn min_deg(&self) -> Option<i32> {
        self.c.keys().next().copied()
    }
    pub fn max_deg(&self) -> Option<i32> {
        self.c.keys().next_back().copied()
    }
    pub fn eval_at_one(&self) -> i128 {
        self.c.values().sum()
    }
    /// Exact division; panics if not divisible. Divisor must have leading coefficient +-1.
    pub fn div_exact(&self, d: &LPoly) -> LPoly {
        assert!(!d.is_zero(), "division by zero LPoly");
        if self.is_zero() {
            return LPoly::zero();
        }
        let (dmax, dlead) = d.c.iter().next_back().map(|(&e, &k)| (e, k)).unwrap();
        assert!(dlead == 1 || dlead == -1, "div_exact needs unit leading coefficient");
        let dmin = d.min_deg().unwrap();
        let mut rem = self.clone();
        let mut quo = LPoly::zero();
        while !rem.is_zero() {
            let (rmax, rlead) = rem.c.iter().next_back().map(|(&e, &k)| (e, k)).unwrap();
            let rmin = rem.min_deg().unwrap();
            if rmax - dmax < rmin - dmin {
                panic!("div_exact: not divisible");
            }
            let e = rmax - dmax;
            let k = rlead * dlead; // dlead = +-1 so 1/dlead = dlead
            quo.add_term(e, k);
            for (&de, &dk) in &d.c {
                rem.add_term(de + e, -dk * k);
            }
        }
        quo
    }
    /// Quantum integer [n] = q^{n-1} + q^{n-3} + ... + q^{1-n}.
    pub fn qint(n: i64) -> LPoly {
        let mut r = LPoly::zero();
        if n > 0 {
            let n = n as i32;
            let mut e = n - 1;
            while e >= 1 - n {
                r.add_term(e, 1);
                e -= 2;
            }
        } else if n < 0 {
            return Self::qint(-n).scale(-1);
        }
        r
    }
    pub fn qfact(n: i64) -> LPoly {
        let mut r = LPoly::one();
        for k in 1..=n {
            r = r.mul(&Self::qint(k));
        }
        r
    }
    /// Symmetric quantum binomial [n choose k].
    pub fn qbinom(n: i64, k: i64) -> LPoly {
        if k < 0 || k > n || n < 0 {
            return LPoly::zero();
        }
        Self::qfact(n).div_exact(&Self::qfact(k).mul(&Self::qfact(n - k)))
    }
}

impl fmt::Display for LPoly {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.c.is_empty() {
            return write!(f, "0");
        }
        let mut first = true;
        for (&e, &k) in self.c.iter().rev() {
            let sign = if k < 0 { "-" } else { "+" };
            let a = k.abs();
            if first {
                if k < 0 {
                    write!(f, "-")?;
                }
            } else {
                write!(f, " {} ", sign)?;
            }
            first = false;
            match (a, e) {
                (_, 0) => write!(f, "{}", a)?,
                (1, _) => write!(f, "q^{}", e)?,
                _ => write!(f, "{}q^{}", a, e)?,
            }
        }
        Ok(())
    }
}
impl fmt::Debug for LPoly {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self)
    }
}

/// Poincaré polynomial: (homological degree t, quantum degree q) -> dimension.
#[derive(Clone, PartialEq, Eq, Default)]
pub struct Poincare {
    pub c: BTreeMap<(i32, i32), i128>,
}

impl Poincare {
    pub fn add_term(&mut self, t: i32, q: i32, k: i128) {
        if k == 0 {
            return;
        }
        let v = self.c.entry((t, q)).or_insert(0);
        *v += k;
        if *v == 0 {
            self.c.remove(&(t, q));
        }
    }
    /// Graded Euler characteristic sum (-1)^t dim q^j.
    pub fn euler(&self) -> LPoly {
        let mut r = LPoly::zero();
        for (&(t, q), &k) in &self.c {
            r.add_term(q, if t.rem_euclid(2) == 0 { k } else { -k });
        }
        r
    }
    pub fn total_rank(&self) -> i128 {
        self.c.values().sum()
    }
}

impl fmt::Display for Poincare {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.c.is_empty() {
            return write!(f, "0");
        }
        let mut parts = vec![];
        for (&(t, q), &k) in &self.c {
            let coeff = if k == 1 { String::new() } else { format!("{}", k) };
            parts.push(format!("{}t^{}q^{}", coeff, t, q));
        }
        write!(f, "{}", parts.join(" + "))
    }
}
impl fmt::Debug for Poincare {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn qbinom_values() {
        let b = LPoly::qbinom(4, 2);
        // q^4 + q^2 + 2 + q^-2 + q^-4
        assert_eq!(b.c.get(&0), Some(&2));
        assert_eq!(b.c.get(&4), Some(&1));
        assert_eq!(b.eval_at_one(), 6);
        assert_eq!(b, b.bar());
        let x = LPoly::qint(3).mul(&LPoly::qint(5));
        assert_eq!(x.div_exact(&LPoly::qint(3)), LPoly::qint(5));
    }
}
