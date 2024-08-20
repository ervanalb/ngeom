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
        fn cos(self) -> Self;
        fn sinc(self) -> Self;
    }

    pub trait InvTrig {
        fn atan2(self, other: Self) -> Self;
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

    pub trait INormSquared {
        type Output;
        fn inorm_squared(self) -> Self::Output;
    }

    pub trait INorm {
        type Output;
        fn inorm(self) -> Self::Output;
    }

    pub trait GeometricINormSquared {
        type Output;
        fn geometric_inorm_squared(self) -> Self::Output;
    }

    pub trait GeometricINorm {
        type Output;
        fn geometric_inorm(self) -> Self::Output;
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
        // but since geometric_normalizing a blade is not exception-free,
        // the output type is left to the implementer
        // so they may choose Option<Self> or similar.
        type Output;
        fn hat(self) -> Self::Output;
    }

    pub trait IHat {
        // Typically, the .ihat() function returns Self,
        // but since geometric_normalizing a blade is not exception-free,
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

            impl Norm for $type {
                type Output = $type;
                fn norm(self) -> $type {
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
        Commutator, Dot, Dual, Exp, GeometricINorm, GeometricINormSquared, GeometricNorm,
        GeometricNormSquared, INorm, INormSquared, Join, Meet, Norm, NormSquared, Project, Reflect,
        Reverse, Transform,
    };
    use crate::scalar::{Ring, Sqrt, Trig};
    use ngeom_macros::gen_algebra;

    gen_algebra!(0, 1, 1);

    impl<T: Ring + Sqrt + Trig> Exp for Bivector<T> {
        type Output = Even<T>;
        fn exp(self) -> Even<T> {
            // This formula works because a geometric_normalized bivector squares to -1
            // allowing us to treat it like the imaginary unit
            let theta = self.geometric_norm();
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
        Commutator, Dot, Dual, Exp, GeometricINorm, GeometricINormSquared, GeometricNorm,
        GeometricNormSquared, INorm, INormSquared, Join, Meet, Norm, NormSquared, Project, Reflect,
        Reverse, Transform,
    };
    use crate::scalar::{Ring, Sqrt, Trig};
    use ngeom_macros::gen_algebra;

    gen_algebra!(0, 1, 1, 1);

    //impl<T: Ring + Sqrt + Trig> Exp for Bivector<T> {
    //    type Output = Even<T>;
    //    fn exp(self) -> Even<T> {
    //        // This formula works because a geometric_normalized simple bivector squares to -1
    //        // allowing us to treat it like the imaginary unit
    //        let theta = self.geometric_norm();
    //        self * theta.sinc() + theta.cos()
    //    }
    //}

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

    // Decompose a bivector into two simple bivectors that sum to it.
    // The first is a euclidean line, and the second is a perpendicular ideal line.
    // If the input bivector is called a,
    // define ahat := the geometric_normalized copy of a (a simple bivector.)
    // The decomposition of a is a = u * ahat + v * I * ahat
    //pub fn decompose(a: Bivector<f32>) -> (Bivector<f32>, Bivector<f32>) {
    //// Option 1
    ////let s = dbg!(-a.dot(a));
    //let s = dbg!(a.geometric_norm_squared());
    //let pi = dbg!(a.meet(a));
    //let v_ahat_perp = pi * (-0.5 / s) * a;
    //let u_ahat = a - pi * (-0.5 / s) * a;
    //dbg!(("Option 1", u_ahat, v_ahat_perp));

    ////// Option 2: idecompose
    //let a2 = a.dual();
    ////let s = dbg!(-a2.dot(a2));
    //let s2 = dbg!(a.geometric_inorm_squared());
    ////let s = dbg!(a2.geometric_norm());
    //let pi = dbg!(a.meet(a));
    //let v_a2hat_perp = pi * (-0.5 / s2) * a2;
    //let u_a2hat = a2 - pi * (-0.5 / s2) * a2;

    //dbg!(("Option 2", v_a2hat_perp.dual(), u_a2hat.dual()));

    //// Exception-free approach to calculating v * I * ahat
    //    let a_dual = a.dual();
    //    let s = dbg!(a.geometric_norm_squared());
    //    let s_inf = dbg!(a.geometric_inorm_squared());
    //    let pi = dbg!(a.meet(a));
    //    let n_euc = pi * a; // * (-0.5 / s)
    //    let n_inf = (a_dual * (-s_inf * 2.) - pi * a_dual).dual(); // * (-0.5 / s2);
    //    let r2 = (n_euc + n_inf) * (-0.5 / (s + s_inf));
    //    let r1 = a - r2;

    //    dbg!((r1, r2))
    //}
}

pub mod metric {
    use crate::blade::{Dot, INorm, INormSquared, Join, Meet, Norm, NormSquared};
    use crate::scalar::InvTrig;

    /// Computes the square of the distance between two geometric_normalized geometric entities.
    /// such as two points, or a point and a line in 2D.
    pub fn perpendicular_distance_squared<G1, G2, SCALAR>(g1: G1, g2: G2) -> SCALAR
    where
        G1: Join<G2, Output: NormSquared<Output = SCALAR>>,
    {
        (g1.join(g2)).norm_squared()
    }

    /// Computes the (signed) distance between two geometric_normalized perpendicular geometric entities
    /// such as two points, or a point and a line in 2D.
    ///
    /// Some kinds of entities will produce a distance that is always positive,
    /// such as the distance between points.
    /// Others will produce a signed distance,
    /// such as the distance between a point and a plane in 3D.
    pub fn perpendicular_distance<G1, G2, SCALAR>(g1: G1, g2: G2) -> SCALAR
    where
        G1: Join<G2, Output: Norm<Output = SCALAR>>,
    {
        (g1.join(g2)).norm()
    }

    //pub fn dist2<G1, G2, SCALAR>(g1: G1, g2: G2) -> SCALAR
    //where
    //    G1: Commutator<G2>,
    //    <G1 as Commutator<G2>>::Output: GeometricINorm + GeometricNorm + Copy,
    //    <<G1 as Commutator<G2>>::Output as GeometricINorm>::Output:
    //        core::ops::Div<<<G1 as Commutator<G2>>::Output as GeometricNorm>::Output, Output = SCALAR>,
    //{
    //    let c = g1.cross(g2);
    //    c.geometric_inorm() / c.geometric_norm()
    //}

    /// Computes the squared distance between two geometric_normalized parallel geometric entities
    /// such as two parallel lines in 2D, or a plane and parallel line in 3D.
    pub fn parallel_distance_squared<G1, G2, SCALAR>(g1: G1, g2: G2) -> SCALAR
    where
        G1: Meet<G2, Output: INormSquared<Output = SCALAR>>,
    {
        (g1.meet(g2)).inorm_squared()
    }

    /// Computes the (signed) distance between two geometric_normalized parallel geometric entities
    /// such as two parallel lines in 2D, or a plane and parallel line in 3D.
    pub fn parallel_distance<G1, G2, SCALAR>(g1: G1, g2: G2) -> SCALAR
    where
        G1: Meet<G2, Output: INorm<Output = SCALAR>>,
    {
        (g1.meet(g2)).inorm()
    }

    /// Computes the (4-quadrant) angle between two geometric_normalized geometric entities
    /// such as two lines in 2D, or two planes in 3D
    pub fn angle<G1, G2, SCALAR>(g1: G1, g2: G2) -> SCALAR
    where
        G1: Copy,
        G2: Copy,
        G1: Dot<G2, Output = SCALAR>,
        G1: Meet<G2, Output: Norm<Output = SCALAR>>,
        SCALAR: InvTrig,
    {
        g1.meet(g2).norm().atan2(g1.dot(g2))
    }
}

#[cfg(test)]
mod test {
    use super::blade::{Commutator, Dot, Exp, Hat, Join, Meet, Norm, NormSquared, Transform};
    use super::metric::{angle, perpendicular_distance};
    use super::scalar::Ring;
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

        // TODO create higher objects like lines & planes
    }

    #[test]
    fn test_dist_2d() {
        // Distance between points
        assert_close!(
            perpendicular_distance(pga2d::point::<f32>([10., 10.]), pga2d::point([13., 14.])),
            5.
        );
        assert_close!(
            perpendicular_distance(
                pga3d::point::<f32>([10., 10., 10.]),
                pga3d::point([13., 14., 10.])
            ),
            5.
        );

        // Signed dist2ance between point & line
        // Note that the line and the test point form a triangle that is wound
        // left-handed, producing a negative dist2ance
        assert_close!(
            perpendicular_distance(
                pga2d::point::<f32>([10., 10.])
                    .join(pga2d::point([10., 20.]))
                    .hat(),
                pga2d::point([15., 15.])
            ),
            -5.
        );
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
            .join(pga3d::point([10., 10., 10.]))
            .hat();
        let plane_xy = pga3d::point([0., 0., 0.])
            .join(pga3d::point([10., 0., 0.]))
            .join(pga3d::point([0., 10., 0.]))
            .hat();
        perpendicular_distance(pt, plane_xy); // GOOD
        perpendicular_distance(pga3d::origin(), pt); // GOOD
        perpendicular_distance(pt, line_z); // GOOD
        perpendicular_distance(pt, line_skew); // GOOD
        dbg!(perpendicular_distance(line_z, line_skew)); // BAD
        dbg!(perpendicular_distance(line_z, line_skew)); // BAD

        let l1 = line_z;
        let l2 = line_z2;

        dbg!(l1.meet(l2));
        dbg!(l1.dot(l2));

        let crossed = dbg!(l1.cross(l2));

        // AAAA
        //let (ua, via) = dbg!(pga3d::decompose(crossed));

        //assert_close!(ua.meet(ua).geometric_inorm(), 0.);
        //assert_close!(via.meet(via).geometric_inorm(), 0.);

        //dbg!(dbg!(via.geometric_inorm_squared()) + dbg!(l1.meet(l2).geometric_inorm_squared()));
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
        //dbg!(line_z.cross(line_skew).geometric_inorm());
        //dbg!(line_z.cross(line_skew).hat());
        //dbg!(line_z.cross(line_skew).hat().geometric_inorm());

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
        dbg!(b.geometric_norm());
        dbg!(b.geometric_inorm());
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
        assert_close!(
            perpendicular_distance(pga3d::point([10., 10., 10.]), pga3d::point([13., 14., 10.])),
            5.
        );

        // Distance between point & line
        assert_close!(
            perpendicular_distance(
                pga3d::point::<f32>([10., 10., 10.])
                    .join(pga3d::point([10., 20., 10.]))
                    .hat(),
                pga3d::point([15., 15., 10.])
            ),
            5.
        );
        assert_close!(
            perpendicular_distance(
                pga3d::point::<f32>([10., 10., 0.])
                    .join(pga3d::point([10., 10., 20.]))
                    .hat(),
                pga3d::point([13., 14., 10.])
            ),
            5.
        );

        // Signed distance between point & plane
        // Note how triangle is wound left-handed WRT the point, to produce a negative distance
        assert_close!(
            perpendicular_distance(
                pga3d::point::<f32>([10., 10., 10.])
                    .join(pga3d::point([10., 20., 10.]))
                    .join(pga3d::point([20., 10., 10.]))
                    .hat(),
                pga3d::point([15., 15., 15.])
            ),
            -5.
        );

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
        //let l1 = pga2d::origin().join(pga2d::point([5., 5.]));
        //let l2 = pga2d::origin().join(pga2d::point([-10., 10.]));
        //let motor = (l2 * l1).hat() + 1.;
        //dbg!(motor);

        let p2 = p1.transform(motor);
        dbg!(p2);
        assert_close!(p2, pga2d::point([0., 10.]));
    }
}
