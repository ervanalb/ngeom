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

    pub trait AntiRecip {
        fn anti_recip(self) -> Self;
    }

    pub trait AntiSqrt {
        fn anti_sqrt(self) -> Self;
    }

    pub trait AntiMul<RHS = Self> {
        type Output;
        fn anti_mul(self, rhs: RHS) -> Self::Output;
    }

    pub trait AntiTrig {
        fn anti_cos(self) -> Self;
        fn anti_sinc(self) -> Self;
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

    pub trait AntiReverse {
        fn anti_reverse(self) -> Self;
    }

    pub trait Bulk {
        fn bulk(self) -> Self;
    }

    pub trait Weight {
        fn weight(self) -> Self;
    }

    pub trait BulkDual {
        type Output;
        fn bulk_dual(self) -> Self::Output;
    }

    pub trait WeightDual {
        type Output;
        fn weight_dual(self) -> Self::Output;
    }

    pub trait RightComplement {
        type Output;
        fn right_complement(self) -> Self::Output;
    }

    pub trait LeftComplement {
        type Output;
        fn left_complement(self) -> Self::Output;
    }

    pub trait BulkNormSquared {
        type Output;
        fn bulk_norm_squared(self) -> Self::Output;
    }

    pub trait WeightNormSquared {
        type Output;
        fn weight_norm_squared(self) -> Self::Output;
    }

    pub trait BulkNorm {
        type Output;
        fn bulk_norm(self) -> Self::Output;
    }

    pub trait WeightNorm {
        type Output;
        fn weight_norm(self) -> Self::Output;
    }

    pub trait Wedge<T> {
        type Output;
        fn wedge(self, r: T) -> Self::Output;
    }

    pub trait AntiWedge<T> {
        type Output;
        fn anti_wedge(self, r: T) -> Self::Output;
    }

    pub trait Dot<T> {
        type Output;
        fn dot(self, r: T) -> Self::Output;
    }

    pub trait AntiDot<T> {
        type Output;
        fn anti_dot(self, r: T) -> Self::Output;
    }

    pub trait WedgeDot<T> {
        type Output;
        fn wedge_dot(self, r: T) -> Self::Output;
    }

    pub trait AntiWedgeDot<T> {
        type Output;
        fn anti_wedge_dot(self, r: T) -> Self::Output;
    }

    pub trait Project<T> {
        fn project(self, r: T) -> Self;
    }

    pub trait Transform<T> {
        fn transform(self, r: T) -> Self;
    }

    pub trait AntiCommutator<T> {
        type Output;
        fn anti_commutator(self, r: T) -> Self::Output;
    }

    pub trait ExpAntiWedgeDot {
        type Output;
        fn exp_anti_wedge_dot(self) -> Self::Output;
    }

    pub trait Normalized {
        fn normalized(self) -> Self;
    }

    pub trait Unitized {
        fn unitized(self) -> Self;
    }

    macro_rules! impl_for_builtin {
        ($type:ident) => {
            impl BulkNormSquared for $type {
                type Output = $type;
                fn bulk_norm_squared(self) -> $type {
                    self * self
                }
            }

            impl BulkNorm for $type {
                type Output = $type;
                fn bulk_norm(self) -> $type {
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
        AntiCommutator, AntiDot, AntiReverse, AntiWedge, AntiWedgeDot, Bulk, BulkDual, BulkNorm,
        BulkNormSquared, Dot, ExpAntiWedgeDot, LeftComplement, Normalized, Project, Reverse,
        RightComplement, Transform, Unitized, Wedge, WedgeDot, Weight, WeightDual, WeightNorm,
        WeightNormSquared,
    };
    use crate::scalar::{AntiMul, AntiRecip, AntiSqrt, AntiTrig, Recip, Ring, Sqrt, Trig};
    use ngeom_macros::gen_algebra;

    gen_algebra!(0, 1, 1);

    impl<T: Ring + Sqrt> AntiSqrt for AntiScalar<T> {
        fn anti_sqrt(self) -> AntiScalar<T> {
            AntiScalar {
                a012: self.a012.sqrt(),
            }
        }
    }

    impl<T: Ring + Recip> AntiRecip for AntiScalar<T> {
        fn anti_recip(self) -> AntiScalar<T> {
            AntiScalar {
                a012: self.a012.recip(),
            }
        }
    }

    impl<T: Ring + Trig> AntiTrig for AntiScalar<T> {
        fn anti_cos(self) -> AntiScalar<T> {
            AntiScalar {
                a012: self.a012.cos(),
            }
        }
        fn anti_sinc(self) -> AntiScalar<T> {
            AntiScalar {
                a012: self.a012.sinc(),
            }
        }
    }

    impl<T: Ring + core::cmp::PartialEq> PartialEq for AntiScalar<T> {
        fn eq(&self, other: &Self) -> bool {
            self.a012.eq(&other.a012)
        }
    }

    impl<T: Ring + core::cmp::Eq> Eq for AntiScalar<T> {}

    impl<T: Ring + core::cmp::PartialOrd> PartialOrd for AntiScalar<T> {
        fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
            self.a012.partial_cmp(&other.a012)
        }
    }

    impl<T: Ring + core::cmp::Ord> Ord for AntiScalar<T> {
        fn cmp(&self, other: &Self) -> core::cmp::Ordering {
            self.a012.cmp(&other.a012)
        }
    }

    impl<T: Ring + Sqrt + Trig> ExpAntiWedgeDot for Vector<T> {
        type Output = Odd<T>;
        fn exp_anti_wedge_dot(self) -> Odd<T> {
            // This formula works because a unitized simple vector squares to -1
            // under the anti-wedge-dot product
            // allowing us to treat it like the imaginary unit
            let theta = self.weight_norm();
            -self.anti_mul(theta.anti_sinc()) + theta.anti_cos()
        }
    }

    pub fn anti<T: Ring>(a012: T) -> AntiScalar<T> {
        AntiScalar { a012 }
    }

    pub fn origin<T: Ring>() -> Vector<T> {
        Vector {
            a0: T::one(),
            a1: T::zero(),
            a2: T::zero(),
        }
    }

    // Construction functions
    pub fn point<T: Ring>([x, y]: [T; 2]) -> Vector<T> {
        Vector {
            a0: T::one(),
            a1: x,
            a2: y,
        }
    }

    pub fn point_ideal<T: Ring>([x, y]: [T; 2]) -> Vector<T> {
        Vector {
            a0: T::zero(),
            a1: x,
            a2: y,
        }
    }

    pub fn point_homogeneous<T: Ring>([x, y, w]: [T; 3]) -> Vector<T> {
        Vector {
            a0: w,
            a1: x,
            a2: y,
        }
    }

    /// Line with equation ax + by - c = 0
    pub fn line<T: Ring>(a: T, b: T, c: T) -> Bivector<T> {
        Bivector {
            a12: c,
            a20: a,
            a01: b,
        }
    }
}

pub mod pga3d {
    use crate::blade::{
        AntiCommutator, AntiDot, AntiReverse, AntiWedge, AntiWedgeDot, Bulk, BulkDual, BulkNorm,
        BulkNormSquared, Dot, ExpAntiWedgeDot, LeftComplement, Normalized, Project, Reverse,
        RightComplement, Transform, Unitized, Wedge, WedgeDot, Weight, WeightDual, WeightNorm,
        WeightNormSquared,
    };
    use crate::scalar::{AntiMul, AntiRecip, AntiSqrt, AntiTrig, Recip, Ring, Sqrt, Trig};
    use ngeom_macros::gen_algebra;

    gen_algebra!(0, 1, 1, 1);

    impl<T: Ring + Sqrt> AntiSqrt for AntiScalar<T> {
        fn anti_sqrt(self) -> AntiScalar<T> {
            AntiScalar {
                a0123: self.a0123.sqrt(),
            }
        }
    }

    impl<T: Ring + Recip> AntiRecip for AntiScalar<T> {
        fn anti_recip(self) -> AntiScalar<T> {
            AntiScalar {
                a0123: self.a0123.recip(),
            }
        }
    }

    impl<T: Ring + Trig> AntiTrig for AntiScalar<T> {
        fn anti_cos(self) -> AntiScalar<T> {
            AntiScalar {
                a0123: self.a0123.cos(),
            }
        }
        fn anti_sinc(self) -> AntiScalar<T> {
            AntiScalar {
                a0123: self.a0123.sinc(),
            }
        }
    }

    impl<T: Ring + core::cmp::PartialEq> PartialEq for AntiScalar<T> {
        fn eq(&self, other: &Self) -> bool {
            self.a0123.eq(&other.a0123)
        }
    }

    impl<T: Ring + core::cmp::Eq> Eq for AntiScalar<T> {}

    impl<T: Ring + core::cmp::PartialOrd> PartialOrd for AntiScalar<T> {
        fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
            self.a0123.partial_cmp(&other.a0123)
        }
    }

    impl<T: Ring + core::cmp::Ord> Ord for AntiScalar<T> {
        fn cmp(&self, other: &Self) -> core::cmp::Ordering {
            self.a0123.cmp(&other.a0123)
        }
    }

    impl<T: Ring + Sqrt + Trig> ExpAntiWedgeDot for Bivector<T> {
        type Output = Even<T>;
        fn exp_anti_wedge_dot(self) -> Even<T> {
            // This formula works because a normalized simple vector squares to -1
            // under the anti-wedge-dot product
            // allowing us to treat it like the imaginary unit
            let theta = self.weight_norm();
            -self.anti_mul(theta.anti_sinc()) + theta.anti_cos()
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

    pub fn anti<T: Ring>(a0123: T) -> AntiScalar<T> {
        AntiScalar { a0123 }
    }

    pub fn origin<T: Ring>() -> Vector<T> {
        Vector {
            a0: T::one(),
            a1: T::zero(),
            a2: T::zero(),
            a3: T::zero(),
        }
    }

    pub fn point<T: Ring>([x, y, z]: [T; 3]) -> Vector<T> {
        Vector {
            a0: T::one(),
            a1: x,
            a2: y,
            a3: z,
        }
    }

    pub fn point_ideal<T: Ring>([x, y, z]: [T; 3]) -> Vector<T> {
        Vector {
            a0: T::zero(),
            a1: x,
            a2: y,
            a3: z,
        }
    }

    pub fn point_homogeneous<T: Ring>([x, y, z, w]: [T; 4]) -> Vector<T> {
        Vector {
            a0: w,
            a1: x,
            a2: y,
            a3: z,
        }
    }

    /// Plane with equation ax + by + cz - d = 0
    pub fn plane<T: Ring>(a: T, b: T, c: T, d: T) -> Trivector<T> {
        Trivector {
            a123: d,
            a032: a,
            a013: b,
            a021: c,
        }
    }
}

#[cfg(test)]
mod test {
    use super::blade::{
        AntiCommutator, AntiWedge, BulkNorm, BulkNormSquared, ExpAntiWedgeDot, Normalized,
        RightComplement, Transform, Unitized, Wedge, WeightDual, WeightNorm, WeightNormSquared,
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

    impl<T: RightComplement<Output = f32>> IsClose for T {
        fn is_close(self, rhs: T) -> bool {
            (self.right_complement() - rhs.right_complement()).abs() < 1e-5
        }
    }

    /*
    impl<
            T: core::ops::Sub<
                T,
                Output: Copy
                            + BulkNormSquared<Output = f32>
                            + WeightNormSquared<Output: RightComplement<Output = f32>>,
            >,
        > IsClose for T
    {
        fn is_close(self, rhs: T) -> bool {
            let tol = 1e-5;
            let diff = self - rhs;
            // Taking the right complement allows us to get the weight norm as a scalar
            // rather than an antiscalar
            diff.bulk_norm_squared() < tol * tol
                && diff.weight_norm_squared().right_complement() < tol * tol
        }
    }
    */

    impl IsClose for pga2d::Vector<f32> {
        fn is_close(self, rhs: Self) -> bool {
            let tol = 1e-5;
            let diff = self - rhs;
            // Taking the right complement allows us to get the weight norm as a scalar
            // rather than an antiscalar
            diff.bulk_norm_squared() < tol * tol
                && diff.weight_norm_squared().right_complement() < tol * tol
        }
    }

    impl IsClose for pga3d::Vector<f32> {
        fn is_close(self, rhs: Self) -> bool {
            let tol = 1e-5;
            let diff = self - rhs;
            // Taking the right complement allows us to get the weight norm as a scalar
            // rather than an antiscalar
            diff.bulk_norm_squared() < tol * tol
                && diff.weight_norm_squared().right_complement() < tol * tol
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

            assert_close!(p1.wedge(p2).weight_norm(), pga2d::anti(5.));
        }

        // Signed distance between point & line
        // Joining the line & point results in a plane (scalar) corresponding to their signed distance
        {
            // Note that the line and the test point form a triangle that is wound
            // left-handed, producing a negative distance
            let l = pga2d::point::<f32>([10., 10.])
                .wedge(pga2d::point([10., 20.]))
                .unitized();
            let p = pga2d::point([15., 15.]);
            assert_close!(l.wedge(p).weight_norm(), pga2d::anti(-5.));
        }

        // Distance between parallel lines
        // The lines meet at an infinite point whose infinite norm is their unsigned distance.
        // More precisely, this expression yields the distance between the projection of the
        // origin onto each line.
        {
            let l1 = pga2d::origin::<f32>()
                .wedge(pga2d::point([3., 4.]))
                .unitized();
            let l2 = pga2d::point::<f32>([4., -3.])
                .wedge(pga2d::point([7., 1.]))
                .unitized();
            assert_close!(dbg!(l1.anti_wedge(l2)).bulk_norm(), 5.);
        }
    }

    #[test]
    fn test_dist_3d() {
        // Distance between points
        // Joining the points results in a line whose norm is their unsigned distance
        {
            let p1 = pga3d::point::<f32>([10., 10., 10.]);
            let p2 = pga3d::point([13., 14., 10.]);
            assert_close!(p1.wedge(p2).weight_norm(), pga3d::anti(5.));
        }

        // Distnce between point & line
        // Joining the line & point results in a plane whose norm is their unsigned distance
        {
            let l = pga3d::point::<f32>([10., 10., 10.])
                .wedge(pga3d::point([10., 20., 10.]))
                .unitized();
            let p = pga3d::point([15., 15., 10.]);
            assert_close!(l.wedge(p).weight_norm(), pga3d::anti(5.));
        }

        {
            let l = pga3d::point::<f32>([10., 10., 0.])
                .wedge(pga3d::point([10., 10., 20.]))
                .unitized();
            let p = pga3d::point([13., 14., 10.]);
            assert_close!(l.wedge(p).weight_norm(), pga3d::anti(5.));
        }

        // Distance between point & plane
        // Joining the plane & point results in a scalar corresponding to their signed distance
        {
            let pl = pga3d::point::<f32>([10., 10., 10.])
                .wedge(pga3d::point([20., 10., 10.]))
                .wedge(pga3d::point([10., 20., 10.]))
                .unitized();
            let pt = pga3d::point([15., 15., 15.]);
            assert_close!(pl.wedge(pt), pga3d::anti(5.));
        }

        // Distance between parallel lines
        // More precisely, this expression yields the distance between the projection of the
        // origin onto each line.
        {
            let l1 = pga3d::origin::<f32>()
                .wedge(pga3d::point([0., 0., 10.]))
                .unitized();
            let l2 = pga3d::point::<f32>([3., 4., 10.])
                .wedge(pga3d::point([3., 4., 20.]))
                .unitized();
            assert_close!(l1.anti_commutator(l2).bulk_norm(), 5.);
        }

        // Distance between perpendicular lines
        // More precisely, this expression yields d * sin(a)
        // d is the distance between the lines, and a is the angle between them
        // (as measured along / about their shared normal)
        {
            let _l1 = pga3d::origin::<f32>()
                .wedge(pga3d::point([0., 0., 10.]))
                .unitized();
            let _l2 = pga3d::point::<f32>([8., 0., 15.])
                .wedge(pga3d::point([8., 20., 15.]))
                .unitized();
            //assert_close!(l1.wedge(l2), pga3d::anti(-8.)); // XXX broken
        }

        // Distance between skew lines
        // This expression divides the join of the two lines (d * sin(a))
        // by the scalar norm of their commutator product, which is sin(a)
        {
            let _l1 = pga3d::origin::<f32>()
                .wedge(pga3d::point([0., 0., 10.]))
                .unitized();
            let _l2 = pga3d::point([10., 0., 0.])
                .wedge(pga3d::point([10., 5., 4.]))
                .unitized();

            //assert_close!(l1.wedge(l2) + l1.anti_commutator(l2).weight_norm(), -10.) // XXX broken
        }
    }

    #[test]
    fn test_sign_2d_triangle_winding() {
        let p1 = pga2d::origin::<f32>();
        let p2 = pga2d::point([1., 0.]);
        let p3 = pga2d::point([0., 1.]);
        assert!(p1.wedge(p2).wedge(p3) > pga2d::anti(0.));
        assert!(p3.wedge(p2).wedge(p1) < pga2d::anti(0.));
    }

    #[test]
    fn test_sign_2d_line_intersection() {
        let l1 = pga2d::origin::<f32>()
            .wedge(pga2d::point([10., 0.]))
            .unitized();
        let l2 = pga2d::point([0., 10.])
            .wedge(pga2d::point([10., 9.]))
            .unitized();
        let l3 = pga2d::point([0., -10.])
            .wedge(pga2d::point([-10., -9.]))
            .unitized();

        assert!(l1.anti_wedge(l2).weight_norm() < pga2d::anti(0.));
        assert!(l2.anti_wedge(l1).weight_norm() > pga2d::anti(0.));
        assert!(l1.anti_wedge(l3).weight_norm() > pga2d::anti(0.));
        assert!(l3.anti_wedge(l1).weight_norm() < pga2d::anti(0.));
    }

    #[test]
    fn test_sign_3d_tetrahedron_winding() {
        let p1 = pga3d::origin::<f32>();
        let p2 = pga3d::point([1., 0., 0.]);
        let p3 = pga3d::point([0., 1., 0.]);
        let p4 = pga3d::point([0., 0., 1.]);
        assert!(p1.wedge(p2).wedge(p3).wedge(p4) > pga3d::anti(0.));
        assert!(p3.wedge(p2).wedge(p1).wedge(p4) < pga3d::anti(0.));
    }

    #[test]
    fn test_sign_3d_skew_lines() {
        // Note that this is intentionally backwards!
        // It's impossible to reconcile the sign of this
        // with the sign of the plane-point version,
        // so we add the negative sign here
        let l1 = pga3d::origin::<f32>()
            .wedge(pga3d::point([0., 0., 10.]))
            .unitized();
        let l2 = pga3d::point::<f32>([8., 0., 15.])
            .wedge(pga3d::point([8., 20., 15.]))
            .unitized();
        let l3 = pga3d::point::<f32>([-10., 0., 0.])
            .wedge(pga3d::point([-10., 20., -5.]))
            .unitized();

        assert!(l1.wedge(l2) > pga3d::anti(0.));
        assert!(l1.wedge(l3) < pga3d::anti(0.));
    }

    #[test]
    fn test_sign_3d_plane_intersection() {
        let p1 = pga3d::origin::<f32>()
            .wedge(pga3d::point([1., 0., 0.]))
            .wedge(pga3d::point([0., 1., 0.]))
            .unitized();
        let p2 = pga3d::origin::<f32>()
            .wedge(pga3d::point([0., 1., 0.]))
            .wedge(pga3d::point([0., 0., 1.]))
            .unitized();
        let p3 = pga3d::origin::<f32>()
            .wedge(pga3d::point([0., 0., 1.]))
            .wedge(pga3d::point([1., 0., 0.]))
            .unitized();

        dbg!(p1.anti_wedge(p2).anti_wedge(p3).weight_norm() > pga3d::anti(0.));
    }

    #[test]
    fn test_2d_line_normal() {
        let l = pga2d::origin::<f32>()
            .wedge(pga2d::point([1., 0.]))
            .unitized();
        let expected_dir = pga2d::point_ideal([0., 1.]);

        let d = l.weight_dual();

        assert!((d - expected_dir).weight_norm_squared() < pga2d::anti(1e-5 * 1e-5));
    }

    #[test]
    fn test_plane_normal() {
        let p = pga3d::origin::<f32>()
            .wedge(pga3d::point([1., 0., 0.]))
            .wedge(pga3d::point([0., 1., 0.]))
            .unitized();

        let expected_dir = pga3d::point_ideal([0., 0., 1.]);

        let d = p.weight_dual();

        assert!((d - expected_dir).weight_norm_squared() < pga3d::anti(1e-5 * 1e-5));
    }

    #[test]
    fn test_3d_line_normal() {
        let l = pga3d::origin::<f32>()
            .wedge(pga3d::point([0., 1., 0.]))
            .unitized();

        let inf_line = l.weight_dual();

        let expected_inf_line =
            pga3d::point_ideal([0., 0., 1.]).wedge(pga3d::point_ideal([1., 0., 0.]));

        assert!((inf_line - expected_inf_line).weight_norm_squared() < pga3d::anti(1e-5 * 1e-5));
    }

    #[test]
    fn test_2d_rotation() {
        let p = pga2d::point([1., 0.]);
        let center = pga2d::origin() * (0.125 * core::f32::consts::TAU);
        let motor = center.exp_anti_wedge_dot();
        assert_close!(p.transform(motor), pga2d::point([0., 1.]));
    }

    #[test]
    fn test_3d_rotation() {
        let p = pga3d::point([1., 0., 0.]);
        let l = pga3d::origin::<f32>()
            .wedge(pga3d::point([0., 0., 1.]))
            .unitized()
            * (0.125 * core::f32::consts::TAU);
        dbg!(l);
        let motor = l.exp_anti_wedge_dot();
        assert_close!(p.transform(motor), pga3d::point([0., 1., 0.]));
    }

    #[test]
    fn test_2d_translation() {
        let p1 = pga2d::point([10., 10.]);

        let center = pga2d::point_ideal([-1., 0.]) * 2.5;
        let motor = center.exp_anti_wedge_dot();

        let p2 = p1.transform(motor);
        assert_close!(p2, pga2d::point([10., 15.]));
    }

    #[test]
    fn test_3d_translation() {
        let p1 = pga3d::point([10., 10., 10.]);

        let line_ideal = pga3d::point_ideal([0., 0., 1.])
            .wedge(pga3d::point_ideal([1., 0., 0.]))
            .normalized()
            * 2.5;
        dbg!(line_ideal);

        let motor = line_ideal.exp_anti_wedge_dot();

        let p2 = p1.transform(motor);
        assert_close!(p2, pga3d::point([10., 15., 10.]));
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
    */
}
