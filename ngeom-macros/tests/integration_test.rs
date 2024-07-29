use core::ops::{Add, Mul, Neg, Sub};

pub trait ScalarRing:
    Clone
    + Copy
    + Neg<Output = Self>
    + Add<Self, Output = Self>
    + Mul<Self, Output = Self>
    + Sub<Self, Output = Self>
{
    fn zero() -> Self;
    fn one() -> Self;
}

pub trait Reverse {
    fn reverse(self) -> Self;
}

pub trait Dual {
    type Output;
    fn dual(self) -> Self::Output;
}

pub trait Conjugate {
    fn conjugate(self) -> Self;
}

pub trait Norm {
    type Output;
    fn norm(self) -> Self::Output;
}

pub trait NormInfinite {
    type Output;
    fn inorm(self) -> Self::Output;
}

pub trait Project<T> {
    fn project(self, r: T) -> Self;
}

pub trait Reflect<T> {
    fn reflect(self, r: T) -> Self;
}

pub trait Transform<T> {
    fn transform(self, r: T) -> Self;
}

pub trait Commutator<T> {
    type Output;
    fn cross(self, r: T) -> Self::Output;
}

mod pga4d {
    use super::{Commutator, ScalarRing, Transform};
    use ngeom_macros::gen_algebra;

    gen_algebra!(1, 1, 1, 0);
}

#[test]
fn test_pga() {}
