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
    fn two() -> Self {
        Self::one() + Self::one()
    }
}

pub trait Sqrt {
    type Output;
    fn sqrt(self) -> Self::Output;
}

pub trait Reverse {
    fn reverse(self) -> Self;
}

pub trait Dual {
    type Output;
    fn dual(self) -> Self::Output;
}

pub trait NormSquared {
    type Output;
    fn bulk_norm_squared(self) -> Self::Output;
    fn weight_norm_squared(self) -> Self::Output;
}

pub trait Norm: NormSquared<Output: Sqrt>
where
    Self: Sized,
{
    fn bulk_norm(self) -> <Self::Output as Sqrt>::Output {
        self.bulk_norm_squared().sqrt()
    }

    fn weight_norm(self) -> <Self::Output as Sqrt>::Output {
        self.weight_norm_squared().sqrt()
    }
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
    use super::{Commutator, Dual, Project, Reflect, Reverse, ScalarRing, Transform};
    use ngeom_macros::gen_algebra;

    gen_algebra!(1, 1, 1, 0);
}

#[test]
fn test_pga() {}
