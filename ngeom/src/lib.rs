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

    pub trait Norm {
        type Output;
        fn norm(self) -> Self::Output;
    }

    pub trait INorm {
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

    pub trait Exp {
        type Output;
        fn exp(self) -> Self::Output;
    }

    pub trait Hat {
        // Typically, the .hat() function returns Self,
        // but since normalizing a blade is not exception-free,
        // the output type is left to the implementer
        // so they may choose Option<Self> or similar.
        type Output;
        fn hat(self) -> Self::Output;
    }
}

use crate::blade::Dual;

pub struct Never;

pub trait Algebra
where
    Self::ScalarDual: Dual<Output = Self::Scalar>,
{
    type Scalar;
    type Point;
    type Line;
    type Plane;
    type E3;
    type E4;

    type ScalarDual; // AKA pseudoscalar
    type PointDual; // AKA hyperplane
    type LineDual;
    type PlaneDual;
    type E3Dual;
    type E4Dual;

    type Motor;
    type MotorDual; // AKA flector
}

pub mod pga2d {
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

    // Construction functions
    pub fn i<T: Ring>() -> Pseudoscalar<T> {
        Pseudoscalar { a012: T::one() }
    }

    pub fn point<T: Ring>(x: T, y: T) -> Bivector<T> {
        Bivector {
            a01: T::one(),
            a20: y,
            a12: x,
        }
    }

    pub fn point_ideal<T: Ring>(x: T, y: T) -> Bivector<T> {
        Bivector {
            a01: T::zero(),
            a20: y,
            a12: x,
        }
    }

    pub fn point_homogeneous<T: Ring>(x: T, y: T, w: T) -> Bivector<T> {
        Bivector {
            a01: w,
            a20: y,
            a12: x,
        }
    }

    /// Line with equation ax + by + c = 0
    pub fn line<T: Ring>(a: T, b: T, c: T) -> Vector<T> {
        Vector {
            a0: a,
            a1: b,
            a2: c,
        }
    }
}

pub mod pga3d {
    use crate::blade::{
        Commutator, Dual, Exp, INorm, INormSquared, Norm, NormSquared, Project, Reflect, Reverse,
        Transform,
    };
    use crate::scalar::{Ring, Sqrt, Trig};
    use crate::Algebra;
    use core::marker::PhantomData;
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

    // Construction functions
    pub fn i<T: Ring>() -> Pseudoscalar<T> {
        Pseudoscalar { a0123: T::one() }
    }

    pub fn point<T: Ring>(x: T, y: T, z: T) -> Trivector<T> {
        Trivector {
            a021: T::one(),
            a013: z,
            a032: y,
            a123: x,
        }
    }

    pub fn point_ideal<T: Ring>(x: T, y: T, z: T) -> Trivector<T> {
        Trivector {
            a021: T::zero(),
            a013: z,
            a032: y,
            a123: x,
        }
    }

    pub fn point_homogeneous<T: Ring>(x: T, y: T, z: T, w: T) -> Trivector<T> {
        Trivector {
            a021: w,
            a013: z,
            a032: y,
            a123: x,
        }
    }

    /// Line with equation ax + by + c = 0
    pub fn plane<T: Ring>(a: T, b: T, c: T, d: T) -> Vector<T> {
        Vector {
            a0: a,
            a1: b,
            a2: c,
            a3: d,
        }
    }

    struct PGA3D<T: Ring + Sqrt + Trig>(PhantomData<T>);

    impl<T: Ring + Sqrt + Trig> Algebra for PGA3D<T> {
        type Scalar = T;
        type Point = Trivector<T>;
        type Line = Bivector<T>;
        type Plane = Vector<T>;
        type E3 = Pseudoscalar<T>;
        type E4 = crate::Never;

        type ScalarDual = Pseudoscalar<T>;
        type PointDual = Vector<T>;
        type LineDual = Bivector<T>;
        type PlaneDual = Trivector<T>;
        type E3Dual = T;
        type E4Dual = crate::Never;

        type Motor = Even<T>;
        type MotorDual = Odd<T>;
    }
}

#[cfg(test)]
mod test {
    use super::blade::{Hat, Norm};
    use super::scalar::Ring;
    use super::*;
    use core::ops::{Div, Mul};

    trait IsClose {
        fn is_close(self, rhs: Self) -> bool;
    }

    impl IsClose for f32 {
        fn is_close(self, rhs: f32) -> bool {
            (self - rhs).abs() < 1e-5
        }
    }

    impl<B: Copy + Norm<Output: Ring> + Mul<<B as Norm>::Output, Output = B>> Hat for B
    where
        <B as Norm>::Output: Div<<B as Norm>::Output, Output = <B as Norm>::Output>,
    {
        type Output = Self;
        fn hat(self) -> Self {
            self * (<B as Norm>::Output::one() / self.norm())
        }
    }

    #[test]
    fn construction() {
        // 2D
        let _point = pga2d::point(3., 4.);
        let _line = pga2d::line(3., 4., 1.);

        // 3D
        let _point = pga3d::point(3., 4., 5.);
        // TODO: line
        let _plane = pga3d::plane(3., 4., 5., 1.);
    }

    #[test]
    fn metric() {
        // 2D distance between points
        let p1 = pga2d::point(10., 10.);
        let p2 = pga2d::point(13., 14.);
        let dist = (p1 & p2).norm();
        assert!(dist.is_close(5.));

        // 3D distance between points
        let p1 = pga3d::point(10., 10., 10.);
        let p2 = pga3d::point(13., 14., 10.);
        let dist = (p1 & p2).norm();
        assert!(dist.is_close(5.));

        // 2D angle of intersecting lines
        let l1 = (pga2d::point::<f32>(0., 0.) & pga2d::point(5., 5.)).hat();
        let l2 = (pga2d::point(0., 0.) & pga2d::point(-10., 10.)).hat();
        let angle = (l1 ^ l2).norm().atan2(l1 | l2);
        assert!(angle.is_close(0.25 * core::f32::consts::TAU));
    }
}
