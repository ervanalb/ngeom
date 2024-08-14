pub mod scalar {
    use core::ops::{Add, Mul, Neg, Sub};
    pub trait Ring: // TODO: add RHS?? unify with Blade??
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
    use crate::scalar::Ring;

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

    pub trait Blade:
        Sized
        + core::ops::Mul<Self::Scalar, Output = Self>
        + core::ops::Add<Self, Output = Self>
        + core::ops::Sub<Self, Output = Self>
        + core::ops::Neg<Output = Self>
        + NormSquared<Output = Self::Scalar>
        + Dual
    //where
    //    <Self as Dual>::Output: Blade<Scalar = Self::Scalar> + Dual<Output = Self>,
    {
        type Scalar: Ring;
    }
}

pub mod pga2d {
    use crate::blade::{
        Blade, Commutator, Dual, Exp, INorm, INormSquared, Norm, NormSquared, Project, Reflect,
        Reverse, Transform,
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

    pub fn i<T: Ring>() -> Pseudoscalar<T> {
        Pseudoscalar { a012: T::one() }
    }

    pub fn origin<T: Ring>() -> Bivector<T> {
        Bivector {
            a01: T::one(),
            a20: T::zero(),
            a12: T::zero(),
        }
    }

    // Construction functions
    pub fn point<T: Ring>([x, y]: [T; 2]) -> Bivector<T> {
        Bivector {
            a01: T::one(),
            a20: y,
            a12: x,
        }
    }

    pub fn point_ideal<T: Ring>([x, y]: [T; 2]) -> Bivector<T> {
        Bivector {
            a01: T::zero(),
            a20: y,
            a12: x,
        }
    }

    pub fn point_homogeneous<T: Ring>([x, y, w]: [T; 3]) -> Bivector<T> {
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

    impl<T: Ring> Blade for Bivector<T> {
        type Scalar = T;
    }
    impl<T: Ring> Blade for Vector<T> {
        type Scalar = T;
    }
}

pub mod pga3d {
    use crate::blade::{
        Blade, Commutator, Dual, Exp, INorm, INormSquared, Norm, NormSquared, Project, Reflect,
        Reverse, Transform,
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

    pub fn i<T: Ring>() -> Pseudoscalar<T> {
        Pseudoscalar { a0123: T::one() }
    }

    pub fn origin<T: Ring>() -> Trivector<T> {
        Trivector {
            a021: T::one(),
            a013: T::zero(),
            a032: T::zero(),
            a123: T::zero(),
        }
    }

    pub fn point<T: Ring>([x, y, z]: [T; 3]) -> Trivector<T> {
        Trivector {
            a021: T::one(),
            a013: z,
            a032: y,
            a123: x,
        }
    }

    pub fn point_ideal<T: Ring>([x, y, z]: [T; 3]) -> Trivector<T> {
        Trivector {
            a021: T::zero(),
            a013: z,
            a032: y,
            a123: x,
        }
    }

    pub fn point_homogeneous<T: Ring>([x, y, z, w]: [T; 4]) -> Trivector<T> {
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

    impl<T: Ring> Blade for Trivector<T> {
        type Scalar = T;
    }
    impl<T: Ring> Blade for Bivector<T> {
        type Scalar = T;
    }
    impl<T: Ring> Blade for Vector<T> {
        type Scalar = T;
    }
}

#[cfg(test)]
mod test {
    use super::blade::{Blade, Exp, Hat, Norm, Transform};
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
        let _point = pga2d::point([3., 4.]);

        // 3D
        let _point = pga3d::point([3., 4., 5.]);
    }

    #[test]
    fn metric() {
        // 2D distance between points
        let p1 = pga2d::point([10., 10.]);
        let p2 = pga2d::point([13., 14.]);
        let dist = (p1 & p2).norm();
        assert!(dist.is_close(5.));

        // 3D distance between points
        let p1 = pga3d::point([10., 10., 10.]);
        let p2 = pga3d::point([13., 14., 10.]);
        let dist = (p1 & p2).norm();
        assert!(dist.is_close(5.));

        // 2D angle of intersecting lines
        let l1 = (pga2d::origin::<f32>() & pga2d::point([5., 5.])).hat();
        let l2 = (pga2d::origin() & pga2d::point([-10., 10.])).hat();
        let angle = (l1 ^ l2).norm().atan2(l1 | l2);
        assert!(angle.is_close(0.25 * core::f32::consts::TAU));

        // 3D angle of intersecting planes
        let pl1 =
            (pga3d::origin::<f32>() & pga3d::point([0., 0., 1.]) & pga3d::point([5., 5., 0.]))
                .hat();
        let pl2 =
            (pga3d::origin() & pga3d::point([0., 0., 1.]) & pga3d::point([-10., 10., 0.])).hat();
        let angle = (pl1 ^ pl2).norm().atan2(pl1 | pl2);
        assert!(angle.is_close(0.25 * core::f32::consts::TAU));
    }

    #[test]
    fn metric_generic() {
        fn metric<P: Blade<Scalar: IsClose>>(p1: P, p2: P, d: P::Scalar)
        where
            P: core::ops::BitAnd<P, Output: Blade<Scalar = P::Scalar> + Norm<Output = P::Scalar>>,
        {
            assert!((p1 & p2).norm().is_close(d));
        }

        metric(pga2d::point([10., 10.]), pga2d::point([13., 14.]), 5.);
        metric(
            pga3d::point([10., 10., 0.]),
            pga3d::point([13., 14., 0.]),
            5.,
        );
    }

    #[test]
    fn motor() {
        // 2D translation
        let p1 = pga2d::point([10., 10.]);

        let center = pga2d::point_ideal([1., 0.]) * 2.5;
        let motor = center.exp();

        let p2 = p1.transform(motor);
        assert!(p2.is_close(pga2d::point([10., 15.])));

        // 2D rotation
        let p1 = pga2d::point([10., 0.]);

        // This motor goes the wrong way??
        let center = pga2d::origin() * (0.125 * core::f32::consts::TAU);
        let motor = center.exp();
        dbg!(motor);

        // This one works
        let l1 = pga2d::origin() & pga2d::point([5., 5.]);
        let l2 = pga2d::origin() & pga2d::point([-10., 10.]);
        let motor = (l2 * l1).hat() + 1.;
        dbg!(motor);

        let p2 = p1.transform(motor);
        dbg!(p2);
        assert!(p2.is_close(pga2d::point([0., 10.])));
        assert!(false);
    }
}
