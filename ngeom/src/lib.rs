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

use crate::blade::{Dual, Norm};
use crate::scalar::{Ring, Sqrt};
use core::marker::PhantomData;

pub struct Never;

pub trait Algebra
where
    Self::ScalarDual: Dual<Output = Self::Scalar>,
    Self::Point: core::ops::BitAnd<Self::Point, Output = Self::Line>,
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

    fn origin() -> Self::Point;
    fn i() -> Self::ScalarDual;
}

pub trait AlgebraWithSqrt: Algebra<Line: Norm<Output = Self::Scalar>> {}

pub trait Blade  {
    type Algebra: Algebra;
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
}

pub struct PGA2D<T: Ring>(PhantomData<T>);

impl<T: Ring> PGA2D<T> {
    // Construction functions
    pub fn point([x, y]: [T; 2]) -> pga2d::Bivector<T> {
        pga2d::Bivector {
            a01: T::one(),
            a20: y,
            a12: x,
        }
    }

    pub fn point_ideal([x, y]: [T; 2]) -> pga2d::Bivector<T> {
        pga2d::Bivector {
            a01: T::zero(),
            a20: y,
            a12: x,
        }
    }

    pub fn point_homogeneous([x, y, w]: [T; 3]) -> pga2d::Bivector<T> {
        pga2d::Bivector {
            a01: w,
            a20: y,
            a12: x,
        }
    }

    /// Line with equation ax + by + c = 0
    pub fn line(a: T, b: T, c: T) -> pga2d::Vector<T> {
        pga2d::Vector {
            a0: a,
            a1: b,
            a2: c,
        }
    }
}

impl<T: Ring> Algebra for PGA2D<T> {
    type Scalar = T;
    type Point = pga2d::Bivector<T>;
    type Line = pga2d::Vector<T>;
    type Plane = pga2d::Pseudoscalar<T>;
    type E3 = crate::Never;
    type E4 = crate::Never;

    type ScalarDual = pga2d::Pseudoscalar<T>;
    type PointDual = pga2d::Vector<T>;
    type LineDual = pga2d::Bivector<T>;
    type PlaneDual = T;
    type E3Dual = crate::Never;
    type E4Dual = crate::Never;

    type Motor = pga2d::Even<T>;
    type MotorDual = pga2d::Odd<T>;

    fn origin() -> pga2d::Bivector<Self::Scalar> {
        pga2d::Bivector {
            a01: Self::Scalar::one(),
            a20: Self::Scalar::zero(),
            a12: Self::Scalar::zero(),
        }
    }

    fn i() -> pga2d::Pseudoscalar<T> {
        pga2d::Pseudoscalar {
            a012: Self::Scalar::one(),
        }
    }
}

impl<T: Ring + Sqrt> AlgebraWithSqrt for PGA2D<T> {}

impl<T: Ring> Blade for pga2d::Bivector<T> {
    type Algebra = PGA2D<T>;
}
//impl<T: Ring + Sqrt> BladeWithSqrt for pga2d::Bivector<T> {}

pub mod pga3d {
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

pub struct PGA3D<T: Ring>(PhantomData<T>);

impl<T: Ring> PGA3D<T> {
    pub fn point([x, y, z]: [T; 3]) -> pga3d::Trivector<T> {
        pga3d::Trivector {
            a021: T::one(),
            a013: z,
            a032: y,
            a123: x,
        }
    }

    pub fn point_ideal([x, y, z]: [T; 3]) -> pga3d::Trivector<T> {
        pga3d::Trivector {
            a021: T::zero(),
            a013: z,
            a032: y,
            a123: x,
        }
    }

    pub fn point_homogeneous([x, y, z, w]: [T; 4]) -> pga3d::Trivector<T> {
        pga3d::Trivector {
            a021: w,
            a013: z,
            a032: y,
            a123: x,
        }
    }

    /// Line with equation ax + by + c = 0
    pub fn plane(a: T, b: T, c: T, d: T) -> pga3d::Vector<T> {
        pga3d::Vector {
            a0: a,
            a1: b,
            a2: c,
            a3: d,
        }
    }
}

impl<T: Ring> Algebra for PGA3D<T> {
    type Scalar = T;
    type Point = pga3d::Trivector<T>;
    type Line = pga3d::Bivector<T>;
    type Plane = pga3d::Vector<T>;
    type E3 = pga3d::Pseudoscalar<T>;
    type E4 = crate::Never;

    type ScalarDual = pga3d::Pseudoscalar<T>;
    type PointDual = pga3d::Vector<T>;
    type LineDual = pga3d::Bivector<T>;
    type PlaneDual = pga3d::Trivector<T>;
    type E3Dual = T;
    type E4Dual = crate::Never;

    type Motor = pga3d::Even<T>;
    type MotorDual = pga3d::Odd<T>;

    fn i() -> pga3d::Pseudoscalar<Self::Scalar> {
        pga3d::Pseudoscalar {
            a0123: Self::Scalar::one(),
        }
    }

    fn origin() -> pga3d::Trivector<Self::Scalar> {
        pga3d::Trivector {
            a021: Self::Scalar::one(),
            a013: Self::Scalar::zero(),
            a032: Self::Scalar::zero(),
            a123: Self::Scalar::zero(),
        }
    }
}

impl<T: Ring + Sqrt> AlgebraWithSqrt for PGA3D<T> {}

#[cfg(test)]
mod test {
    use super::blade::{Exp, Hat, Norm, Transform};
    use super::scalar::Ring;
    use super::AlgebraWithSqrt;
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

    impl<B: Hat<Output: core::ops::BitAnd<<B as Hat>::Output, Output: Norm<Output = f32>>>> IsClose
        for B
    {
        fn is_close(self, rhs: B) -> bool {
            (self.hat() & rhs.hat()).norm() < 1e-5
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
        let _point = PGA2D::point([3., 4.]);

        // 3D
        let _point = PGA3D::point([3., 4., 5.]);
    }

    #[test]
    fn metric() {
        // 2D distance between points
        let p1 = PGA2D::point([10., 10.]);
        let p2 = PGA2D::point([13., 14.]);
        let dist = (p1 & p2).norm();
        assert!(dist.is_close(5.));

        // 3D distance between points
        let p1 = PGA3D::point([10., 10., 10.]);
        let p2 = PGA3D::point([13., 14., 10.]);
        let dist = (p1 & p2).norm();
        assert!(dist.is_close(5.));

        // 2D angle of intersecting lines
        let l1 = (PGA2D::<f32>::origin() & PGA2D::point([5., 5.])).hat();
        let l2 = (PGA2D::origin() & PGA2D::point([-10., 10.])).hat();
        let angle = (l1 ^ l2).norm().atan2(l1 | l2);
        assert!(angle.is_close(0.25 * core::f32::consts::TAU));

        // 3D angle of intersecting planes
        let pl1 =
            (PGA3D::<f32>::origin() & PGA3D::point([0., 0., 1.]) & PGA3D::point([5., 5., 0.]))
                .hat();
        let pl2 =
            (PGA3D::origin() & PGA3D::point([0., 0., 1.]) & PGA3D::point([-10., 10., 0.])).hat();
        let angle = (pl1 ^ pl2).norm().atan2(pl1 | pl2);
        assert!(angle.is_close(0.25 * core::f32::consts::TAU));
    }

    #[test]
    fn metric_generic() {
        trait TestAlgebra: Algebra {
            fn p1() -> Self::Point;
            fn p2() -> Self::Point;
            fn d() -> Self::Scalar;
        }

        impl TestAlgebra for PGA3D<f32> {
            fn p1() -> Self::Point {
                PGA3D::point([10., 10., 0.])
            }
            fn p2() -> Self::Point {
                PGA3D::point([13., 14., 0.])
            }
            fn d() -> Self::Scalar {
                5.
            }
        }

        impl TestAlgebra for PGA2D<f32> {
            fn p1() -> Self::Point {
                PGA2D::point([10., 10.])
            }
            fn p2() -> Self::Point {
                PGA2D::point([13., 14.])
            }
            fn d() -> Self::Scalar {
                5.
            }
        }

        fn metric<A: AlgebraWithSqrt<Scalar: IsClose> + TestAlgebra>() {
            let p1 = A::p1();
            let p2 = A::p2();
            let dist = (p1 & p2).norm();
            assert!(dist.is_close(A::d()));
        }

        metric::<PGA2D<f32>>();
        metric::<PGA3D<f32>>();
    }

    #[test]
    fn motor() {
        // 2D translation
        let p1 = PGA2D::point([10., 10.]);

        let center = PGA2D::point_ideal([1., 0.]) * 2.5;
        let motor = center.exp();

        let p2 = p1.transform(motor);
        assert!(p2.is_close(PGA2D::point([10., 15.])));

        // 2D rotation
        let p1 = PGA2D::point([10., 0.]);

        // This motor goes the wrong way??
        let center = PGA2D::origin() * (0.125 * core::f32::consts::TAU);
        let motor = center.exp();
        dbg!(motor);

        // This one works
        let l1 = PGA2D::origin() & PGA2D::point([5., 5.]);
        let l2 = PGA2D::origin() & PGA2D::point([-10., 10.]);
        let motor = (l2 * l1).hat() + 1.;
        dbg!(motor);

        let p2 = p1.transform(motor);
        dbg!(p2);
        assert!(p2.is_close(PGA2D::point([0., 10.])));
        assert!(false);
    }
}
