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

    pub trait InvTrig {
        fn atan2(self, other: Self) -> Self;
    }

    pub trait PartialSqrt: Sized {
        type Output;

        fn partial_sqrt(self) -> Option<Self::Output>;
    }

    pub trait Sqrt: PartialSqrt {
        fn sqrt(self) -> Self::Output;
    }

    pub trait PartialInv: Sized {
        fn partial_inv(self) -> Option<Self>;
    }

    pub trait PartialInvSqrt: Sized {
        fn partial_inv_sqrt(self) -> Option<Self>;
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

            impl PartialSqrt for $type {
                type Output = $type;

                fn partial_sqrt(self) -> Option<$type> {
                    Some(self.sqrt())
                }
            }

            impl Sqrt for $type {
                fn sqrt(self) -> $type {
                    self.partial_sqrt().unwrap()
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

            impl InvTrig for $type {
                fn atan2(self, other: $type) -> $type {
                    self.atan2(other)
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

    pub trait PartialNorm {
        type Output;
        fn partial_norm(self) -> Option<Self::Output>;
    }

    pub trait Norm: PartialNorm {
        fn norm(self) -> Self::Output;
    }

    pub trait GeometricNormSquared {
        type Output;
        fn geometric_norm_squared(self) -> Self::Output;
    }

    pub trait PartialGeometricNorm {
        type Output;
        fn partial_geometric_norm(self) -> Option<Self::Output>;
    }

    pub trait GeometricNorm: PartialGeometricNorm {
        fn geometric_norm(self) -> Self::Output;
    }

    pub trait InfNormSquared {
        type Output;
        fn inf_norm_squared(self) -> Self::Output;
    }

    pub trait PartialInfNorm {
        type Output;
        fn partial_inf_norm(self) -> Option<Self::Output>;
    }

    pub trait InfNorm: PartialInfNorm {
        fn inf_norm(self) -> Self::Output;
    }

    pub trait GeometricInfNormSquared {
        type Output;
        fn geometric_inf_norm_squared(self) -> Self::Output;
    }

    pub trait PartialGeometricInfNorm {
        type Output;
        fn partial_geometric_inf_norm(self) -> Option<Self::Output>;
    }

    pub trait GeometricInfNorm: PartialGeometricInfNorm {
        fn geometric_inf_norm(self) -> Self::Output;
    }

    pub trait Join<T> {
        type Output;
        fn join(self, r: T) -> Self::Output;
    }

    pub trait Meet<T> {
        type Output;
        fn meet(self, r: T) -> Self::Output;
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

    pub trait IHat {
        // Typically, the .ihat() function returns Self,
        // but since normalizing a blade is not exception-free,
        // the output type is left to the implementer
        // so they may choose Option<Self> or similar.
        type Output;
        fn ihat(self) -> Self::Output;
    }

    macro_rules! impl_for_builtin {
        ($type:ident) => {
            impl NormSquared for $type {
                type Output = $type;
                fn norm_squared(self) -> $type {
                    self * self
                }
            }

            impl PartialNorm for $type {
                type Output = $type;
                fn partial_norm(self) -> Option<$type> {
                    Some(self)
                }
            }

            impl Norm for $type {
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

            impl PartialGeometricNorm for $type {
                type Output = $type;
                fn partial_geometric_norm(self) -> Option<$type> {
                    Some(self)
                }
            }

            impl GeometricNorm for $type {
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
        Commutator, Dot, Dual, Exp, GeometricInfNorm, GeometricInfNormSquared, GeometricNorm,
        GeometricNormSquared, InfNorm, InfNormSquared, Join, Meet, Norm, NormSquared,
        PartialGeometricInfNorm, PartialGeometricNorm, PartialInfNorm, PartialNorm, Project,
        Reflect, Reverse, Transform,
    };
    use crate::scalar::{PartialSqrt, Ring, Sqrt, Trig};
    use ngeom_macros::gen_algebra;

    gen_algebra!(0, 1, 1);

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
            a12: T::one(),
        }
    }

    pub fn point_ideal<T: Ring>([x, y]: [T; 2]) -> Bivector<T> {
        Bivector {
            a01: y,
            a20: x,
            a12: T::zero(),
        }
    }

    pub fn point_homogeneous<T: Ring>([x, y, w]: [T; 3]) -> Bivector<T> {
        Bivector {
            a01: y,
            a20: x,
            a12: w,
        }
    }

    /// Line with equation ax + by + c = 0
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
        Commutator, Dot, Dual, Exp, GeometricInfNorm, GeometricInfNormSquared, GeometricNorm,
        GeometricNormSquared, InfNorm, InfNormSquared, Join, Meet, Norm, NormSquared,
        PartialGeometricInfNorm, PartialGeometricNorm, PartialInfNorm, PartialNorm, Project,
        Reflect, Reverse, Transform,
    };
    use crate::scalar::{PartialSqrt, Ring, Sqrt, Trig};
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

    pub fn i<T: Ring>() -> Pseudoscalar<T> {
        Pseudoscalar { a0123: T::one() }
    }

    pub fn origin<T: Ring>() -> Trivector<T> {
        Trivector {
            a021: T::zero(),
            a013: T::zero(),
            a032: T::zero(),
            a123: T::one(),
        }
    }

    pub fn point<T: Ring>([x, y, z]: [T; 3]) -> Trivector<T> {
        Trivector {
            a021: z,
            a013: y,
            a032: x,
            a123: T::one(),
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
            a123: w,
        }
    }

    /// Line with equation ax + by + c = 0
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
        Commutator, Dot, Exp, GeometricInfNormSquared, GeometricNormSquared, Hat, InfNorm, Join,
        Meet, Norm, NormSquared, PartialNorm, Transform,
    };
    use super::scalar::{PartialSqrt, Ring, Sqrt};
    use super::{pga2d, pga3d};
    use core::ops::{Div, Mul};

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
        G: Join<G, Output: NormSquared<Output = f32>>,
    {
        fn is_close(self, rhs: G) -> bool {
            let tol = 1e-5;
            (self.join(rhs)).norm_squared() < tol * tol
        }
    }

    impl<B: Copy + Norm<Output: Ring> + Mul<<B as PartialNorm>::Output, Output = B>> Hat for B
    where
        <B as PartialNorm>::Output:
            Div<<B as PartialNorm>::Output, Output = <B as PartialNorm>::Output>,
    {
        type Output = Self;
        fn hat(self) -> Self {
            self * (<B as PartialNorm>::Output::one() / self.norm())
        }
    }

    impl<T: Ring + Sqrt<Output = T> + core::ops::Div<Output = T>> PartialSqrt for pga3d::Magnitude<T> {
        type Output = pga3d::Magnitude<T>;

        fn partial_sqrt(self) -> Option<pga3d::Magnitude<T>> {
            let pga3d::Magnitude { a: s, a0123: p } = self;
            let sqrt_s = s.sqrt();
            // This may divide by zero, but we ignore it for the purposes of unit testing
            Some(pga3d::Magnitude {
                a: sqrt_s,
                a0123: p / (sqrt_s + sqrt_s),
            })
        }
    }

    impl<T: Ring + Sqrt<Output = T> + core::ops::Div<Output = T>> Sqrt for pga3d::Magnitude<T> {
        fn sqrt(self) -> pga3d::Magnitude<T> {
            self.partial_sqrt().unwrap()
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

            assert_close!(p1.join(p2).norm(), 5.);
        }

        // Signed distance between point & line
        // Joining the line & point results in a plane (scalar) corresponding to their signed distance
        {
            // Note that the line and the test point form a triangle that is wound
            // left-handed, producing a negative distance
            let l = pga2d::point::<f32>([10., 10.])
                .join(pga2d::point([10., 20.]))
                .hat();
            let p = pga2d::point([15., 15.]);
            assert_close!(l.join(p).norm(), -5.);
        }

        // Distance between parallel lines
        // The lines meet at an infinite point whose infinite norm is their unsigned distance
        {
            let l1 = pga2d::origin::<f32>().join(pga2d::point([3., 4.])).hat();
            let l2 = pga2d::point::<f32>([4., -3.])
                .join(pga2d::point([7., 1.]))
                .hat();
            assert_close!(dbg!(l1.meet(l2)).inf_norm(), 5.);
        }
    }

    #[test]
    fn dot_product_fiddle() {
        let pt = pga3d::point([5., 5., 5.]);
        let line_z = pga3d::origin::<f32>()
            .join(pga3d::point([0., 0., 10.]))
            .hat();
        let line_z2 = pga3d::point::<f32>([3., 4., 10.])
            .join(pga3d::point([3., 4., 15.]))
            .hat();
        let line_x = pga3d::origin::<f32>()
            .join(pga3d::point([10., 0., 0.]))
            .hat();
        let line_skew = pga3d::point([10., 0., 0.])
            .join(pga3d::point([10., 5., 4.]))
            .hat();
        let plane_xy = pga3d::point([0., 0., 0.])
            .join(pga3d::point([10., 0., 0.]))
            .join(pga3d::point([0., 10., 0.]))
            .hat();
        let line_z2_ish = pga3d::point::<f32>([3., 4., 10.])
            .join(pga3d::point([3.1, 4.1, 15.]))
            .hat();
        //distance(pt, plane_xy); // GOOD
        //distance(pt, plane_xy); // GOOD
        //distance(pga3d::origin(), pt); // GOOD
        //distance(pt, line_z); // GOOD

        //distance(pt, line_skew); // GOOD
        //dbg!(distance(line_z, line_skew)); // BAD
        //dbg!(distance(line_z, line_skew)); // BAD

        {
            for i in 0..50 {
                let line_z = pga3d::point::<f32>([3., 4., 0.])
                    .join(pga3d::point([3., 4., 10.]))
                    .hat();

                let x = -(2_f32.powi(-5 * i));
                let line_z2_ish = pga3d::origin().join(pga3d::point([x, 0., 1.])).hat();

                dbg!(x);
                dbg!(line_z);
                dbg!(line_z2_ish);
                dbg!(line_z.cross(line_z2_ish).geometric_norm_squared().sqrt());
                dbg!(line_z.cross(line_z2_ish).inf_norm());
            }
        }

        let l1 = line_z;
        let l2 = line_z2_ish;

        //let l1 = pga3d::point::<f32>([10., 10., 10.]);
        //let l2 = pga3d::point([13., 14., 10.]);

        //dbg!(l1.meet(l2));
        //dbg!(l1.dot(l2));

        let crossed = dbg!(l1.cross(l2));
        let joined = l1.join(l2);

        // AAAA
        let norm = dbg!(crossed.geometric_norm_squared().sqrt());
        dbg!(crossed.geometric_inf_norm_squared().sqrt());
        dbg!(crossed.norm());
        dbg!(crossed.inf_norm());

        let d_cos_sq = dbg!(norm.a0123 * norm.a0123);
        let d_sin_sq = joined.norm_squared();
        dbg!((dbg!(d_cos_sq) + dbg!(d_sin_sq)).sqrt());
        //dbg!(skew_distance(line_z, line_skew));
        //dbg!(skew_distance(line_z, line_z2));
        //dbg!(skew_distance(line_z, line_x));
        //dbg!(skew_distance(line_z2, line_x));
        //dbg!(skew_distance(line_x, line_skew));

        //assert_close!(ua.meet(ua).geometric_inf_norm(), 0.);
        //assert_close!(via.meet(via).geometric_inf_norm(), 0.);

        //dbg!(dbg!(via.geometric_inf_norm_squared()) + dbg!(l1.meet(l2).geometric_inf_norm_squared()));
        // AAAA

        //dbg!(ua.meet(ua));
        //dbg!(via.meet(via));
        //dbg!(ua + via);
        //dbg!(ua.hat());
        //dbg!(via.hat());

        //let line_z2 = pga3d::point([5., 5., 5.])
        //    .join(pga3d::point([5., 5., 15.]))
        //    .hat();

        //dbg!(line_z.join(line_skew));
        //dbg!(line_z.meet(line_skew));
        //dbg!(line_z.dot(line_skew));
        //dbg!(line_z.cross(line_z2));
        //dbg!(line_z.cross(line_x));
        //dbg!(line_z.cross(line_skew));
        //dbg!(line_z.cross(line_skew).geometric_norm());
        //dbg!(line_z.cross(line_skew).geometric_inf_norm());
        //dbg!(line_z.cross(line_skew).hat());
        //dbg!(line_z.cross(line_skew).hat().geometric_inf_norm());

        //dbg!(line_z);
        //dbg!(line_x);
        //dbg!(line_skew);
        //dbg!(pga3d::point_ideal([1., 0., 0.]).join(pga3d::point_ideal([0., 1., 0.])));

        let t = dbg!(pga3d::point_ideal([1., 0., 0.])
            .join(pga3d::point_ideal([0., 1., 0.]))
            .exp());
        let r = dbg!(pga3d::point([0., 0., 0.])
            .join(pga3d::point([0., 0., 1.]))
            .exp());
        let m = dbg!(t * r);

        let m = t;

        let b = pga3d::Bivector::<f32> {
            a01: m.a01,
            a02: m.a02,
            a03: m.a03,
            a12: m.a12,
            a31: m.a31,
            a23: m.a23,
        };
        let s = (-b.dot(b)).sqrt();
        let p = -b.meet(b) * (1. / (2. * s));

        dbg!(s);
        dbg!(p);

        let tan = s.atan2(m.a);
        let u = tan / s;
        let v = p * (-tan / (s * s) + 1. / (m.a * s));

        dbg!(u);
        dbg!(v);
        //dbg!(b.geometric_norm());
        //dbg!(b.geometric_inf_norm());
        let log_b = dbg!(b * u + v * b);
        dbg!(log_b.exp());

        //let r_recovered = pga3d::Even {
        //    a: 0.,
        //    a01: 0.,
        //    a02: 0.,
        //    a03: 0.,
        //    a12: m.a12,
        //    a31: m.a31,
        //    a23: m.a23,
        //    a0123: m.a0123
        //};

        //dist(line_skew, plane_xy); // BAD
        //dist(line_z, plane_xy); // BAD
        assert!(false);
    }

    #[test]
    fn test_dist_3d() {
        // Distance between points
        // Joining the points results in a line whose norm is their unsigned distance
        {
            let p1 = pga3d::point::<f32>([10., 10., 10.]);
            let p2 = pga3d::point([13., 14., 10.]);
            assert_close!(p1.join(p2).norm(), 5.);
        }

        // Distance between point & line
        // Joining the line & point results in a plane whose norm is their unsigned distance
        {
            let l = pga3d::point::<f32>([10., 10., 10.])
                .join(pga3d::point([10., 20., 10.]))
                .hat();
            let p = pga3d::point([15., 15., 10.]);
            assert_close!(l.join(p).norm(), 5.);
        }

        {
            let l = pga3d::point::<f32>([10., 10., 0.])
                .join(pga3d::point([10., 10., 20.]))
                .hat();
            let p = pga3d::point([13., 14., 10.]);
            assert_close!(l.join(p).norm(), 5.);
        }

        // Distance between point & plane
        // Joining the plane & point results in a volume (scalar) corresponding to their signed distance
        {
            let pl = pga3d::point::<f32>([10., 10., 10.])
                .join(pga3d::point([10., 20., 10.]))
                .join(pga3d::point([20., 10., 10.]))
                .hat();
            let pt = pga3d::point([15., 15., 15.]);
            assert_close!(pl.join(pt), -5.);
        }

        // Distance between parallel lines
        {
            let l1 = pga3d::origin::<f32>()
                .join(pga3d::point([0., 0., 10.]))
                .hat();
            let l2 = pga3d::point::<f32>([3., 4., 10.])
                .join(pga3d::point([3., 4., 20.]))
                .hat();
            assert_close!(dbg!(l1.cross(l2)).inf_norm(), 5.);
        }

        // Distance between perpendicular lines
        {
            let l1 = pga3d::origin::<f32>()
                .join(pga3d::point([0., 0., 10.]))
                .hat();
            let l2 = pga3d::point::<f32>([8., 0., 15.])
                .join(pga3d::point([8., 20., 15.]))
                .hat();
            assert_close!(dbg!(l1.join(l2)), 8.);
        }

        // Distance between skew lines

        // Closest distance between two lines
        //assert_close!(
        //    dist(
        //        pga3d::origin::<f32>()
        //            .join(pga3d::point([0., 0., 10.]))
        //            .hat(),
        //        pga3d::point([5., 0., 0.])
        //            .join(pga3d::point([5., 1., 1.]))
        //            .hat(),
        //    ),
        //    5.
        //);

        // Parallel lines
        //assert_close!(
        //    dist(
        //        pga3d::origin::<f32>().join(pga3d::point([0., 0., 10.])).hat(),
        //        pga3d::point([3., 4., 10.]).join(pga3d::point([3., 4., 20.])).hat(),
        //    ),
        //    5.
        //);
    }

    /*
    #[test]
    fn test_angle() {
        // Angle between lines in 2D
        assert_close!(
            angle(
                pga2d::origin::<f32>().join(pga2d::point([0., 5.])),
                pga2d::origin().join(pga2d::point([-10., 10.]))
            ),
            0.125 * core::f32::consts::TAU,
        );

        // Angle between planes in 3D
        assert_close!(
            angle(
                pga3d::origin::<f32>()
                    .join(pga3d::point([0., 0., 1.]))
                    .join(pga3d::point([0., 5., 0.])),
                pga3d::origin()
                    .join(pga3d::point([0., 0., 1.]))
                    .join(pga3d::point([-10., 10., 0.]))
            ),
            0.125 * core::f32::consts::TAU,
        );

        {
            // Angle between line and plane
            let pl = pga3d::origin::<f32>()
                .join(pga3d::point([1., 0., 0.]))
                .join(pga3d::point([0., 1., 0.]))
                .hat();
            let l = pga3d::point([10., 10., 0.])
                .join(pga3d::point([10., 20., 10.]))
                .hat();

            // TODO sign is wrong here, why?
            assert_close!(
                pl.meet(l).norm().atan2(pl.dot(l).norm()),
                0.125 * core::f32::consts::TAU,
            )
        }
    }
    */

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
        //let l1 = pga2d::origin().join(pga2d::point([5., 5.]));
        //let l2 = pga2d::origin().join(pga2d::point([-10., 10.]));
        //let motor = (l2 * l1).hat() + 1.;
        //dbg!(motor);

        let p2 = p1.transform(motor);
        dbg!(p2);
        assert_close!(p2, pga2d::point([0., 10.]));
    }
}
