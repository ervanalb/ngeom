pub mod scalar {
    use core::ops::{Add, Mul, Neg, Sub};
    pub trait Ring:
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
        fn sqrt(self) -> Self;
    }

    pub trait Trig {
        fn sin(self) -> Self;
        fn cos(self) -> Self;
        fn sinc(self) -> Self;
    }

    impl Ring for f32 {
        fn zero() -> f32 {
            0.
        }
        fn one() -> f32 {
            1.
        }
    }

    impl Sqrt for f32 {
        fn sqrt(self) -> f32 {
            self.sqrt()
        }
    }

    impl Trig for f32 {
        fn sin(self) -> f32 {
            self.sin()
        }
        fn cos(self) -> f32 {
            self.cos()
        }
        fn sinc(self) -> f32 {
            let self_adj = self.abs() + f32::EPSILON;
            self_adj.sin() / self_adj
        }
    }

    impl Ring for f64 {
        fn zero() -> f64 {
            0.
        }
        fn one() -> f64 {
            1.
        }
    }

    impl Sqrt for f64 {
        fn sqrt(self) -> f64 {
            self.sqrt()
        }
    }

    impl Trig for f64 {
        fn sin(self) -> f64 {
            self.sin()
        }
        fn cos(self) -> f64 {
            self.cos()
        }
        fn sinc(self) -> f64 {
            let self_adj = self.abs() + f64::EPSILON;
            self_adj.sin() / self_adj
        }
    }

    impl Ring for i8 {
        fn zero() -> i8 {
            0
        }
        fn one() -> i8 {
            1
        }
    }

    impl Ring for i16 {
        fn zero() -> i16 {
            0
        }
        fn one() -> i16 {
            1
        }
    }

    impl Ring for i32 {
        fn zero() -> i32 {
            0
        }
        fn one() -> i32 {
            1
        }
    }

    impl Ring for i64 {
        fn zero() -> i64 {
            0
        }
        fn one() -> i64 {
            1
        }
    }
}

pub mod blade {
    use crate::scalar::Sqrt;

    pub trait Reverse {
        fn reverse(self) -> Self;
    }

    pub trait Dual {
        type Output;
        fn dual(self) -> Self::Output;
    }

    pub trait NormSquared {
        type Output;
        fn norm_squared(self) -> Self::Output;
    }

    pub trait INormSquared {
        type Output;
        fn inorm_squared(self) -> Self::Output;
    }

    pub trait Norm: NormSquared<Output: Sqrt>
    where
        Self: Sized,
    {
        fn norm(self) -> Self::Output {
            self.norm_squared().sqrt()
        }
    }
    pub trait INorm: INormSquared<Output: Sqrt>
    where
        Self: Sized,
    {
        fn inorm(self) -> Self::Output {
            self.inorm_squared().sqrt()
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

    pub trait Exp {
        type Output;
        fn exp(self) -> Self::Output;
    }
}

use crate::blade::Dual;

pub struct Never;

pub trait Algebra
where
    Self::Scalar: Dual<Output = Self::ScalarDual>,
{
    type Scalar;
    type Point;
    type Line;
    type Plane;
    type Space;

    type ScalarDual;
    type PointDual;
    type LineDual;
    type PlaneDual;
    type SpaceDual;

    type Motor;
    type MotorDual;
}

mod pga2d {
    use crate::blade::{
        Commutator, Dual, Exp, INorm, INormSquared, Norm, NormSquared, Project, Reflect, Reverse,
        Transform,
    };
    use crate::scalar::{Ring, Sqrt, Trig};
    use ngeom_macros::gen_algebra;

    gen_algebra!(1, 1, 0);

    impl<T: Ring + Sqrt + Trig> Exp for Bivector<T> {
        type Output = Even<T>;
        fn exp(self) -> Even<T> {
            // This formula works because a normalized bivector squares to -1
            // allowing us to treat it like the imaginary unit
            let theta = self.norm();
            self * theta.sinc() + theta.cos()
        }
    }
}

mod pga3d {
    use crate::blade::{
        Commutator, Dual, Exp, INorm, INormSquared, Norm, NormSquared, Project, Reflect, Reverse,
        Transform,
    };
    use crate::scalar::{Ring, Sqrt, Trig};
    use ngeom_macros::gen_algebra;

    gen_algebra!(1, 1, 1, 0);

    impl<T: Ring + Sqrt + Trig> Exp for Bivector<T> {
        type Output = Even<T>;
        fn exp(self) -> Even<T> {
            // This formula works because a normalized bivector squares to -1
            // allowing us to treat it like the imaginary unit
            let theta = self.norm();
            self * theta.sinc() + theta.cos()
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn construction() {
        // 2D
        let _point = pga2d::Bivector {
            a01: 1.,
            a20: 4.,
            a12: 3.,
        };
        let _line = pga2d::Vector {
            a0: 3.,
            a1: 4.,
            a2: 1.,
        };

        // 3D
        let _point = pga3d::Trivector {
            a021: 1.,
            a013: 5.,
            a032: 4.,
            a123: 3.,
        };
        // TODO: line
        let _plane = pga3d::Vector {
            a0: 3.,
            a1: 4.,
            a2: 5.,
            a3: 1.,
        };
    }
}
