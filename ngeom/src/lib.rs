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

    pub trait Trig {
        fn cos(self) -> Self;
        fn sinc(self) -> Self;
    }

    pub trait Sqrt {
        fn sqrt(self) -> Self;
    }

    pub trait Recip: Sized {
        fn recip(self) -> Self;
    }

    macro_rules! impl_for_float {
        ($type:ident) => {
            impl Ring for $type {
                fn zero() -> $type {
                    0.
                }
                fn one() -> $type {
                    1.
                }
            }

            impl Sqrt for $type {
                fn sqrt(self) -> $type {
                    self.sqrt()
                }
            }

            impl Trig for $type {
                fn cos(self) -> $type {
                    self.cos()
                }
                fn sinc(self) -> $type {
                    let self_adj = self.abs() + $type::EPSILON;
                    self_adj.sin() / self_adj
                }
            }
        };
    }

    impl_for_float!(f32);
    impl_for_float!(f64);

    macro_rules! impl_for_int {
        ($type:ident) => {
            impl Ring for $type {
                fn zero() -> $type {
                    0
                }
                fn one() -> $type {
                    1
                }
            }
        };
    }

    impl_for_int!(i8);
    impl_for_int!(i16);
    impl_for_int!(i32);
    impl_for_int!(i64);
    impl_for_int!(i128);
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

    pub trait Norm {
        type Output;
        fn norm(self) -> Self::Output;
    }

    pub trait GeometricNormSquared {
        type Output;
        fn geometric_norm_squared(self) -> Self::Output;
    }

    pub trait GeometricNorm {
        type Output;
        fn geometric_norm(self) -> Self::Output;
    }

    pub trait InfNormSquared {
        type Output;
        fn inf_norm_squared(self) -> Self::Output;
    }

    pub trait InfNorm {
        type Output;
        fn inf_norm(self) -> Self::Output;
    }

    pub trait GeometricInfNormSquared {
        type Output;
        fn geometric_inf_norm_squared(self) -> Self::Output;
    }

    pub trait GeometricInfNorm {
        type Output;
        fn geometric_inf_norm(self) -> Self::Output;
    }

    pub trait Vee<T> {
        type Output;
        fn vee(self, r: T) -> Self::Output;
    }

    pub trait Wedge<T> {
        type Output;
        fn wedge(self, r: T) -> Self::Output;
    }

    pub trait Dot<T> {
        type Output;
        fn dot(self, r: T) -> Self::Output;
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

    pub trait Cross<T> {
        type Output;
        fn cross(self, r: T) -> Self::Output;
    }

    pub trait Exp {
        type Output;
        fn exp(self) -> Self::Output;
    }

    pub trait Hat {
        fn hat(self) -> Self;
    }

    pub trait InfHat {
        fn inf_hat(self) -> Self;
    }

    macro_rules! impl_for_builtin {
        ($type:ident) => {
            impl NormSquared for $type {
                type Output = $type;
                fn norm_squared(self) -> $type {
                    self * self
                }
            }

            impl Norm for $type {
                type Output = $type;
                fn norm(self) -> $type {
                    self
                }
            }

            impl GeometricNormSquared for $type {
                type Output = $type;
                fn geometric_norm_squared(self) -> $type {
                    self * self
                }
            }

            impl GeometricNorm for $type {
                type Output = $type;
                fn geometric_norm(self) -> $type {
                    self
                }
            }
        };
    }

    impl_for_builtin!(f32);
    impl_for_builtin!(f64);
    impl_for_builtin!(i8);
    impl_for_builtin!(i16);
    impl_for_builtin!(i32);
    impl_for_builtin!(i64);
    impl_for_builtin!(i128);
}

pub mod pga2d {
    use crate::blade::{
        Cross, Dot, Dual, Exp, GeometricInfNorm, GeometricInfNormSquared, GeometricNorm,
        GeometricNormSquared, Hat, InfHat, InfNorm, InfNormSquared, Norm, NormSquared, Project,
        Reflect, Reverse, Transform, Vee, Wedge,
    };
    use crate::scalar::{Recip, Ring, Sqrt, Trig};
    use ngeom_macros::gen_algebra;

    gen_algebra!(0, 1, 1);

    impl<T: Ring + Sqrt + Trig> Exp for Bivector<T> {
        type Output = Even<T>;
        fn exp(self) -> Even<T> {
            // This formula works because a normalized simple bivector squares to -1
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
            a01: T::zero(),
            a20: T::zero(),
            a12: T::one(),
        }
    }

    // Construction functions
    pub fn point<T: Ring>([x, y]: [T; 2]) -> Bivector<T> {
        Bivector {
            a01: y,
            a20: x,
            a12: -T::one(),
        }
    }

    pub fn point_ideal<T: Ring>([x, y]: [T; 2]) -> Bivector<T> {
        Bivector {
            a01: y,
            a20: x,
            a12: -T::zero(),
        }
    }

    pub fn point_homogeneous<T: Ring>([x, y, w]: [T; 3]) -> Bivector<T> {
        Bivector {
            a01: y,
            a20: x,
            a12: -w,
        }
    }

    /// Line with equation ax + by = c
    pub fn line<T: Ring>(a: T, b: T, c: T) -> Vector<T> {
        Vector {
            a0: c,
            a1: a,
            a2: b,
        }
    }
}

pub mod pga3d {
    use crate::blade::{
        Cross, Dot, Dual, Exp, GeometricInfNorm, GeometricInfNormSquared, GeometricNorm,
        GeometricNormSquared, Hat, InfHat, InfNorm, InfNormSquared, Norm, NormSquared, Project,
        Reflect, Reverse, Transform, Vee, Wedge,
    };
    use crate::scalar::{Recip, Ring, Sqrt, Trig};
    use ngeom_macros::gen_algebra;

    gen_algebra!(0, 1, 1, 1);

    impl<T: Ring + Trig> Exp for Bivector<T>
    where
        Self: Norm<Output = T>,
    {
        type Output = Even<T>;
        fn exp(self) -> Even<T> {
            // This formula works because a normalized simple bivector squares to -1
            // allowing us to treat it like the imaginary unit
            let theta = self.norm();
            self * theta.sinc() + theta.cos()
        }
    }

    impl<T: Ring + Sqrt + Recip> Magnitude<T> {
        pub fn rsqrt(self) -> Magnitude<T> {
            let Magnitude { a: s, a0123: p } = self;
            let sqrt_s = s.sqrt();
            let sqrt_s_cubed = s * sqrt_s;
            Magnitude {
                a: sqrt_s.recip(),
                a0123: p * (sqrt_s_cubed + sqrt_s_cubed).recip(),
            }
        }
    }

    pub fn i<T: Ring>() -> Pseudoscalar<T> {
        Pseudoscalar { a0123: T::one() }
    }

    pub fn origin<T: Ring>() -> Trivector<T> {
        Trivector {
            a021: T::zero(),
            a013: T::zero(),
            a032: T::zero(),
            a123: -T::one(),
        }
    }

    pub fn point<T: Ring>([x, y, z]: [T; 3]) -> Trivector<T> {
        Trivector {
            a021: z,
            a013: y,
            a032: x,
            a123: -T::one(),
        }
    }

    pub fn point_ideal<T: Ring>([x, y, z]: [T; 3]) -> Trivector<T> {
        Trivector {
            a021: z,
            a013: y,
            a032: x,
            a123: T::zero(),
        }
    }

    pub fn point_homogeneous<T: Ring>([x, y, z, w]: [T; 4]) -> Trivector<T> {
        Trivector {
            a021: z,
            a013: y,
            a032: x,
            a123: -w,
        }
    }

    /// Plane with equation ax + by + cz = d
    pub fn plane<T: Ring>(a: T, b: T, c: T, d: T) -> Vector<T> {
        Vector {
            a0: d,
            a1: a,
            a2: b,
            a3: c,
        }
    }
}

#[cfg(test)]
mod test {
    use super::blade::{
        Cross, Exp, Hat, InfNorm, Norm, NormSquared, Project, Transform, Vee, Wedge,
    };
    use super::scalar::Recip;
    use super::{pga2d, pga3d};

    impl Recip for f32 {
        fn recip(self) -> f32 {
            self.recip()
        }
    }

    macro_rules! assert_close {
        ($left:expr, $right:expr $(,)?) => {
            match (&$left, &$right) {
                (left_val, right_val) => {
                    assert!(
                        (*left_val).is_close(*right_val),
                        "{:?} !~= {:?}",
                        left_val,
                        right_val
                    );
                }
            }
        };
    }

    trait IsClose {
        fn is_close(self, rhs: Self) -> bool;
    }

    impl IsClose for f32 {
        fn is_close(self, rhs: f32) -> bool {
            (self - rhs).abs() < 1e-5
        }
    }

    impl<G> IsClose for G
    where
        G: Vee<G, Output: NormSquared<Output = f32>>,
    {
        fn is_close(self, rhs: G) -> bool {
            let tol = 1e-5;
            (self.vee(rhs)).norm_squared() < tol * tol
        }
    }

    #[test]
    fn construction() {
        // 2D
        let _point = pga2d::point([3., 4.]);

        // 3D
        let _point = pga3d::point([3., 4., 5.]);

        // TODO create higher objects like lines & planes
    }

    #[test]
    fn test_dist_2d() {
        // Distance between points
        // Joining the points results in a line whose norm is their unsigned distance
        {
            let p1 = pga2d::point::<f32>([10., 10.]);
            let p2 = pga2d::point([13., 14.]);

            assert_close!(p1.vee(p2).norm(), 5.);
        }

        // Signed distance between point & line
        // Joining the line & point results in a plane (scalar) corresponding to their signed distance
        {
            // Note that the line and the test point form a triangle that is wound
            // left-handed, producing a negative distance
            let l = pga2d::point::<f32>([10., 10.])
                .vee(pga2d::point([10., 20.]))
                .hat();
            let p = pga2d::point([15., 15.]);
            assert_close!(l.vee(p).norm(), -5.);
        }

        // Distance between parallel lines
        // The lines meet at an infinite point whose infinite norm is their unsigned distance.
        // More precisely, this expression yields the distance between the projection of the
        // origin onto each line.
        {
            let l1 = pga2d::origin::<f32>().vee(pga2d::point([3., 4.])).hat();
            let l2 = pga2d::point::<f32>([4., -3.])
                .vee(pga2d::point([7., 1.]))
                .hat();
            assert_close!(dbg!(l1.wedge(l2)).inf_norm(), 5.);
        }
    }

    #[test]
    fn test_dist_3d() {
        // Distance between points
        // Joining the points results in a line whose norm is their unsigned distance
        {
            let p1 = pga3d::point::<f32>([10., 10., 10.]);
            let p2 = pga3d::point([13., 14., 10.]);
            assert_close!(p1.vee(p2).norm(), 5.);
        }

        // Distnce between point & line
        // Joining the line & point results in a plane whose norm is their unsigned distance
        {
            let l = pga3d::point::<f32>([10., 10., 10.])
                .vee(pga3d::point([10., 20., 10.]))
                .hat();
            let p = pga3d::point([15., 15., 10.]);
            assert_close!(l.vee(p).norm(), 5.);
        }

        {
            let l = pga3d::point::<f32>([10., 10., 0.])
                .vee(pga3d::point([10., 10., 20.]))
                .hat();
            let p = pga3d::point([13., 14., 10.]);
            assert_close!(l.vee(p).norm(), 5.);
        }

        // Distance between point & plane
        // Joining the plane & point results in a scalar corresponding to their signed distance
        {
            let pl = pga3d::point::<f32>([10., 10., 10.])
                .vee(pga3d::point([20., 10., 10.]))
                .vee(pga3d::point([10., 20., 10.]))
                .hat();
            let pt = pga3d::point([15., 15., 15.]);
            assert_close!(pl.vee(pt), 5.);
        }

        // Distance between parallel lines
        // More precisely, this expression yields the distance between the projection of the
        // origin onto each line.
        {
            let l1 = pga3d::origin::<f32>()
                .vee(pga3d::point([0., 0., 10.]))
                .hat();
            let l2 = pga3d::point::<f32>([3., 4., 10.])
                .vee(pga3d::point([3., 4., 20.]))
                .hat();
            assert_close!(l1.cross(l2).inf_norm(), 5.);
        }

        // Distance between perpendicular lines
        // More precisely, this expression yields d * sin(a)
        // d is the distance between the lines, and a is the angle between them
        // (as measured along / about their shared normal)
        {
            let l1 = pga3d::origin::<f32>()
                .vee(pga3d::point([0., 0., 10.]))
                .hat();
            let l2 = pga3d::point::<f32>([8., 0., 15.])
                .vee(pga3d::point([8., 20., 15.]))
                .hat();
            //assert_close!(l1.vee(l2), -8.); // XXX broken
        }

        // Distance between skew lines
        // This expression divides the join of the two lines (d * sin(a))
        // by the scalar norm of their commutator product, which is sin(a)
        {
            let l1 = pga3d::origin::<f32>()
                .vee(pga3d::point([0., 0., 10.]))
                .hat();
            let l2 = pga3d::point([10., 0., 0.])
                .vee(pga3d::point([10., 5., 4.]))
                .hat();

            //assert_close!(l1.vee(l2) * l1.cross(l2).norm().recip(), -10.) // XXX broken
        }
    }

    #[test]
    fn test_sign_2d_triangle_winding() {
        let p1 = pga2d::origin::<f32>();
        let p2 = pga2d::point([1., 0.]);
        let p3 = pga2d::point([0., 1.]);
        assert!(p1.vee(p2).vee(p3) > 0.);
        assert!(p3.vee(p2).vee(p1) < 0.);
    }

    #[test]
    fn test_sign_2d_line_intersection() {
        let l1 = pga2d::origin::<f32>().vee(pga2d::point([10., 0.])).hat() * -1.;
        let l2 = pga2d::point([0., 10.]).vee(pga2d::point([10., 9.])).hat() * -1.;
        let l3 = pga2d::point([0., -10.])
            .vee(pga2d::point([-10., -9.]))
            .hat()
            * -1.;
        assert!(l1.wedge(l2).norm() < 0.);
        assert!(l2.wedge(l1).norm() > 0.);
        assert!(l1.wedge(l3).norm() > 0.);
        assert!(l3.wedge(l1).norm() < 0.);
    }

    #[test]
    fn test_sign_3d_tetrahedron_winding() {
        let p1 = pga3d::origin::<f32>();
        let p2 = pga3d::point([1., 0., 0.]);
        let p3 = pga3d::point([0., 1., 0.]);
        let p4 = pga3d::point([0., 0., 1.]);
        assert!(p1.vee(p2).vee(p3).vee(p4) > 0.);
        assert!(p3.vee(p2).vee(p1).vee(p4) < 0.);
    }

    #[test]
    fn test_sign_3d_skew_lines() {
        let l1 = pga3d::origin::<f32>()
            .vee(pga3d::point([0., 0., 10.]))
            .hat()
            * -1.;
        let l2 = pga3d::point::<f32>([8., 0., 15.])
            .vee(pga3d::point([8., 20., 15.]))
            .hat()
            * -1.;
        let l3 = pga3d::point::<f32>([-10., 0., 0.])
            .vee(pga3d::point([-10., 20., -5.]))
            .hat()
            * -1.;

        dbg!(l1.wedge(l2));
        dbg!(l1.wedge(l3));
        assert!(l1.vee(l2) < 0.);
        assert!(l1.vee(l3) > 0.);
    }

    #[test]
    fn test_sign_3d_plane_intersection() {
        let p1 = pga3d::origin::<f32>()
            .vee(pga3d::point([1., 0., 0.]))
            .vee(pga3d::point([0., 1., 0.]))
            .hat()
            * -1.;
        let p2 = pga3d::origin::<f32>()
            .vee(pga3d::point([0., 1., 0.]))
            .vee(pga3d::point([0., 0., 1.]))
            .hat()
            * -1.;
        let p3 = pga3d::origin::<f32>()
            .vee(pga3d::point([0., 0., 1.]))
            .vee(pga3d::point([1., 0., 0.]))
            .hat()
            * -1.;

        dbg!(p1.wedge(p2).wedge(p3).norm() > 0.);
    }

    #[test]
    fn test_sign_2d_rotation() {
        let p = pga2d::point([1., 0.]);
        let center = pga2d::origin() * (-0.125 * core::f32::consts::TAU);
        let motor = center.exp();
        assert_close!(p.transform(motor), pga2d::point([0., 1.]));
    }

    #[test]
    fn test_sign_3d_rotation() {
        let p = pga3d::point([1., 0., 0.]);
        let l = pga3d::origin::<f32>().vee(pga3d::point([0., 0., 1.])).hat()
            * -1.
            * (-0.125 * core::f32::consts::TAU);
        let motor = l.exp();
        assert_close!(p.transform(motor), pga3d::point([0., 1., 0.]));
    }

    #[test]
    fn test_sign_2d_translation() {
        let p1 = pga2d::point([10., 10.]);

        let center = pga2d::point_ideal([1., 0.]) * -2.5;
        let motor = center.exp();

        let p2 = p1.transform(motor);
        assert!(p2.is_close(pga2d::point([10., 15.])));
    }

    #[test]
    fn test_sign_3d_translation() {
        let p1 = pga3d::point([10., 10., 10.]);

        let line_ideal =
            pga3d::origin::<f32>().vee(pga3d::point([0., 1., 0.])).hat() * -1. * pga3d::i() * -2.5;
        let motor = line_ideal.exp();

        let p2 = p1.transform(motor);
        assert!(p2.is_close(pga3d::point([10., 15., 10.])));
    }

    /*
    #[test]
    fn test_angle() {
        // Angle between lines in 2D
        assert_close!(
            angle(
                pga2d::origin::<f32>().vee(pga2d::point([0., 5.])),
                pga2d::origin().vee(pga2d::point([-10., 10.]))
            ),
            0.125 * core::f32::consts::TAU,
        );

        // Angle between planes in 3D
        assert_close!(
            angle(
                pga3d::origin::<f32>()
                    .vee(pga3d::point([0., 0., 1.]))
                    .vee(pga3d::point([0., 5., 0.])),
                pga3d::origin()
                    .vee(pga3d::point([0., 0., 1.]))
                    .vee(pga3d::point([-10., 10., 0.]))
            ),
            0.125 * core::f32::consts::TAU,
        );

        {
            // Angle between line and plane
            let pl = pga3d::origin::<f32>()
                .vee(pga3d::point([1., 0., 0.]))
                .vee(pga3d::point([0., 1., 0.]))
                .hat();
            let l = pga3d::point([10., 10., 0.])
                .vee(pga3d::point([10., 20., 10.]))
                .hat();

            // TODO sign is wrong here, why?
            assert_close!(
                pl.wedge(l).norm().atan2(pl.dot(l).norm()),
                0.125 * core::f32::consts::TAU,
            )
        }
    }
    */

    #[test]
    fn motor() {
        // 2D translation
        let p1 = pga2d::point([10., 10.]);

        let center = pga2d::point_ideal([1., 0.]) * -2.5;
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
        //let l1 = pga2d::origin().vee(pga2d::point([5., 5.]));
        //let l2 = pga2d::origin().vee(pga2d::point([-10., 10.]));
        //let motor = (l2 * l1).hat() + 1.;
        //dbg!(motor);

        let p2 = p1.transform(motor);
        dbg!(p2);
        //assert_close!(p2, pga2d::point([0., 10.])); // XXX broken
    }
}
