//! Coefficient fields: prime fields F_p (p < 2^32) and the rationals Q.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Signed, ToPrimitive, Zero};
use std::fmt;

pub trait Field: Clone + PartialEq + Eq + std::hash::Hash + fmt::Debug + fmt::Display + Send + Sync + 'static {
    fn zero() -> Self;
    fn one() -> Self;
    fn is_zero(&self) -> bool;
    fn add(&self, o: &Self) -> Self;
    fn sub(&self, o: &Self) -> Self;
    fn mul(&self, o: &Self) -> Self;
    fn neg(&self) -> Self;
    /// Multiplicative inverse; panics on zero.
    fn inv(&self) -> Self;
    fn from_i64(n: i64) -> Self;
    /// 0 for Q.
    fn characteristic() -> u64;
    fn name() -> String;

    #[inline]
    fn is_one(&self) -> bool {
        *self == Self::one()
    }
    #[inline]
    fn add_assign(&mut self, o: &Self) {
        *self = self.add(o);
    }
    #[inline]
    fn sub_assign(&mut self, o: &Self) {
        *self = self.sub(o);
    }
    #[inline]
    fn mul_assign(&mut self, o: &Self) {
        *self = self.mul(o);
    }
    /// self += a*b
    #[inline]
    fn fma(&mut self, a: &Self, b: &Self) {
        *self = self.add(&a.mul(b));
    }
}

/// Prime field F_P. Requires P < 2^32 so that products fit in u64.
#[derive(Copy, Clone, PartialEq, Eq, Hash, Default)]
pub struct Fp<const P: u64>(pub u64);

impl<const P: u64> Fp<P> {
    #[inline]
    pub fn new(x: u64) -> Self {
        Fp(x % P)
    }
    pub fn pow(&self, mut e: u64) -> Self {
        let mut b = self.0;
        let mut r = 1u64 % P;
        while e > 0 {
            if e & 1 == 1 {
                r = r * b % P;
            }
            b = b * b % P;
            e >>= 1;
        }
        Fp(r)
    }
    /// Symmetric representative in (-P/2, P/2].
    pub fn signed(&self) -> i64 {
        if self.0 > P / 2 { self.0 as i64 - P as i64 } else { self.0 as i64 }
    }
}

impl<const P: u64> fmt::Debug for Fp<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.signed())
    }
}
impl<const P: u64> fmt::Display for Fp<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.signed())
    }
}

impl<const P: u64> Field for Fp<P> {
    #[inline]
    fn zero() -> Self {
        Fp(0)
    }
    #[inline]
    fn one() -> Self {
        Fp(1 % P)
    }
    #[inline]
    fn is_zero(&self) -> bool {
        self.0 == 0
    }
    #[inline]
    fn add(&self, o: &Self) -> Self {
        let s = self.0 + o.0;
        Fp(if s >= P { s - P } else { s })
    }
    #[inline]
    fn sub(&self, o: &Self) -> Self {
        Fp(if self.0 >= o.0 { self.0 - o.0 } else { self.0 + P - o.0 })
    }
    #[inline]
    fn mul(&self, o: &Self) -> Self {
        Fp(self.0 * o.0 % P)
    }
    #[inline]
    fn neg(&self) -> Self {
        Fp(if self.0 == 0 { 0 } else { P - self.0 })
    }
    fn inv(&self) -> Self {
        assert!(self.0 != 0, "inverse of zero in F_{}", P);
        self.pow(P - 2)
    }
    #[inline]
    fn from_i64(n: i64) -> Self {
        Fp(n.rem_euclid(P as i64) as u64)
    }
    fn characteristic() -> u64 {
        P
    }
    fn name() -> String {
        format!("F_{}", P)
    }
    #[inline]
    fn fma(&mut self, a: &Self, b: &Self) {
        self.0 = (self.0 + a.0 * b.0) % P;
    }
}

/// The default large prime 2^31 - 1.
pub type F31 = Fp<2147483647>;
/// A second large prime (largest prime below 2^32) for two-prime checks.
pub type F32b = Fp<4294967291>;
pub type F2 = Fp<2>;
pub type F3 = Fp<3>;
pub type F5 = Fp<5>;
pub type F7 = Fp<7>;

/// Rational numbers (arbitrary precision).
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct Q(pub BigRational);

impl fmt::Debug for Q {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}
impl fmt::Display for Q {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl Q {
    pub fn to_i64(&self) -> Option<i64> {
        if self.0.is_integer() { self.0.to_integer().to_i64() } else { None }
    }
    pub fn frac(n: i64, d: i64) -> Self {
        Q(BigRational::new(BigInt::from(n), BigInt::from(d)))
    }
    pub fn is_negative(&self) -> bool {
        self.0.is_negative()
    }
}

impl Field for Q {
    fn zero() -> Self {
        Q(BigRational::zero())
    }
    fn one() -> Self {
        Q(BigRational::one())
    }
    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
    fn add(&self, o: &Self) -> Self {
        Q(&self.0 + &o.0)
    }
    fn sub(&self, o: &Self) -> Self {
        Q(&self.0 - &o.0)
    }
    fn mul(&self, o: &Self) -> Self {
        Q(&self.0 * &o.0)
    }
    fn neg(&self) -> Self {
        Q(-&self.0)
    }
    fn inv(&self) -> Self {
        assert!(!self.0.is_zero(), "inverse of zero in Q");
        Q(self.0.recip())
    }
    fn from_i64(n: i64) -> Self {
        Q(BigRational::from_integer(BigInt::from(n)))
    }
    fn characteristic() -> u64 {
        0
    }
    fn name() -> String {
        "Q".to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn fp_basic() {
        let a = F31::from_i64(-3);
        let b = a.inv();
        assert!(a.mul(&b).is_one());
        assert_eq!(F7::from_i64(10), F7::from_i64(3));
        assert_eq!(F2::from_i64(3).0, 1);
    }
    #[test]
    fn q_basic() {
        let a = Q::frac(2, 3);
        assert_eq!(a.mul(&a.inv()), Q::one());
        assert_eq!(Q::from_i64(4).add(&Q::from_i64(-4)), Q::zero());
    }
}
