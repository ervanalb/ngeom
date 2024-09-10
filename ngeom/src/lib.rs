//! ngeom is a library for doing geometry in N dimensions.
//!
//! ngeom provides primitives for points, lines, planes, etc.,
//! along with rigid transformations (rotations, translations, and reflections.)
//!
//! ngeom uses [homogeneous coordinates](https://en.wikipedia.org/wiki/Homogeneous_coordinates) to express ideal/infinite points,
//! ideal/infinite lines, etc. and to provide for exception-free [meet](ops::Meet) & [join](ops::Join) operations.
//!
//! ngeom is generic over the [scalar] datatype, and can be used with f32, f64, or custom datatypes.
//! It can even use integers, with some restrictions--
//! most functionality requires only [scalar addition and multiplication](scalar::Ring).
//!
//! ngeom ships with modules for [2D](pga2d) and [3D](pga3d) Euclidean space.
//! Modules for spaces with different dimensionality or metric can be written with the help of
//! auto-generated operations implementations from the [provided macro](ngeom_macros::gen_algebra).
//! With some effort, usage code can be written to be generic
//! over the dimensionality or even the metric of the space.
//!
//! Under the hood, ngeom is able to make these generalizations by using [geometric algebra](https://en.wikipedia.org/wiki/Geometric_algebra).
//! Understanding geometric algebra is helpful but optional,
//! as the functions are named for their geometric interpretation.
//! Lower-level functions named after their [geometric algebra expressions](algebraic_ops) are available if needed.

/// Traits that govern the scalar data type used by ngeom
///
/// Because geometric operations are generally done through various sums and products,
/// the scalar datatype needs only be a [Ring](scalar::Ring) for most functionality to work.
pub mod scalar {
    use core::ops::{Add, Mul, Neg, Sub};

    /// A scalar datatype which is closed under addition and multiplication.
    ///
    /// see <https://en.wikipedia.org/wiki/Ring_(mathematics)>
    /// `Ring` is implemented for `f32`, `f64`, and `i8` through `i128`
    ///
    /// `Ring` requires that its datatype is `Copy` to avoid the need to clone or borrow when writing
    /// mathematical expressions. If your scalar datatype is expensive to copy,
    /// consider implementing `Ring` on a reference-counted container that is `Copy`.
    pub trait Ring:
        Clone
        + Copy
        + Neg<Output = Self>
        + Add<Self, Output = Self>
        + Mul<Self, Output = Self>
        + Sub<Self, Output = Self>
    {
        /// The additive identity
        fn zero() -> Self;

        /// The multiplicative identity
        fn one() -> Self;
    }

    /// A scalar datatype which can represent fractional values such as ½.
    ///
    /// In homogeneous coordinates, rational numbers are rarely needed due to the `w` coordinate
    /// which in many cases can act as the denominator of an integer coordinate.
    /// However, since the application of a motor causes a motion
    /// that is twice that of the geometry it is built from,
    /// it is often desirable to divide an angle or a distance by 2
    /// which, in most cases, cannot be accomplished
    /// by altering the weight of the support geometry.
    ///
    /// Failing to implement `Rational` means that certain functions for constructing motors
    /// will not be available.
    ///
    /// `Rational` comes implemented for `f32` and `f64`.
    pub trait Rational {
        /// A scalar value that when added to itself, equals [one](Ring::one)
        fn one_half() -> Self;
    }

    /// A scalar datatype which is closed under the square root function.
    ///
    /// Taking the square root is used frequently in geometry,
    /// such as when taking [norms](crate::ops::WeightNorm)
    /// or halving the motion of a motor.
    ///
    /// Failing to implement `Sqrt` means that some [norms](crate::ops::WeightNorm) may not be available.
    /// Norms whose square root can be taken symbolically will still be available.
    /// [Squared norms](crate::ops::WeightNormSquared) will always be available.
    ///
    /// `Sqrt` comes implemented for `f32` and `f64`.
    ///
    /// ## `sqrt()` of negative numbers
    ///
    /// When given a negative value,
    /// this function must either return a valid scalar datatype (e.g. `f32::NaN`)
    /// or panic. There are no other provisions for exception handling at this level.
    ///
    /// That being said, all uses of `sqrt()` within the library
    /// are on values that are guaranteed by design to be non-negative,
    /// making internal use NaN-free.
    pub trait Sqrt {
        // This scalar's positive square root. Only valid for non-negative numbers.
        fn sqrt(self) -> Self;
    }

    /// A scalar datatype which implements trigonometric functions.
    ///
    /// Taking sines and cosines of angles is used frequently in geometry,
    /// such as when constructing [rotors](crate::pga3d::axis_angle).
    ///
    /// Failing to implement `Trig` means that motors may only be constructed from geometry,
    /// and things like slerp will not be possible.
    ///
    /// Unlike most of the other scalar traits here,
    /// `Trig` allows for a different input and output type.
    /// This allows the library user to draw a distinction between
    /// linear quantities (e.g. output of `sin()`)
    /// and angular quantities (e.g. input of `sin()`.)
    ///
    /// `Trig` comes implemented for `f32` → `f32` and `f64` → `f64`.
    pub trait Trig {
        // The scalar datatype for linear quantities
        // (output of `sin()` and `cos()`)
        type Output;

        // The cosine of a scalar (in radians)
        fn cos(self) -> Self::Output;

        // The sine of a scalar (in radians)
        fn sin(self) -> Self::Output;

        // Computes sin(x) / x
        // (including at `0`, where the result should be `1`)
        fn sinc(self) -> Self::Output;
    }

    /// A scalar datatype whose reciprocal can be taken.
    /// Implementing this trait enables some convenience functions,
    /// but is generally not safe given the possibility of dividing by zero.
    ///
    /// In homogeneous coordinates, division is rarely needed due to the `w` coordinate,
    /// which in many cases can act as the denominator of an integer coordinate.
    /// Division becomes necessary when it comes time to
    /// [unitize](crate::ops::Unitized) or [normalize](crate::ops::Normalized) geometry,
    /// projecting it down so that its [weight](crate::ops::WeightNorm) or [bulk norm](crate::ops::BulkNorm) becomes unity.
    /// Failing to implement `Recip` means that these convenience methods will not be available.
    ///
    /// `Recip` comes implemented for `f32` → `f32` and `f64` → `f64`.
    ///
    /// ## `recip()` of `0`
    ///
    /// When given an input of zero,
    /// this function must return a valid scalar datatype (e.g. `f32::NaN`) or panic.
    /// There are no other provisions for exception handling at this level.
    ///
    /// For scalar datatypes that may be zero, this operation is NOT NaN-free,
    /// including internal to the library.
    ///
    /// An easy but error-prone solution is to check before calling `recip()`:
    ///
    /// ```
    /// use ngeom::pga2d::*;
    /// use ngeom::ops::*;
    /// use ngeom::scalar::*;
    ///
    /// /// Return the unitized point of intersection of two unitized lines,
    /// /// or None if they are parallel
    /// fn try_meet(l1: Bivector<f32>, l2: Bivector<f32>) -> Option<Vector<f32>> {
    ///     const PARALLEL_EPSILON: f32 = 1e-5;
    ///
    ///     let result = l1.meet(l2);
    ///     let norm: f32 = result.weight_norm().into();
    ///     if norm.abs() < PARALLEL_EPSILON { return None; } // Don't forget!
    ///     Some(result * norm.recip())
    /// }
    /// ```
    ///
    /// If you desire stronger NaN-free enforcement, consider using branded float types.
    ///
    /// ```ignore
    /// // TODO this doesn't work yet!
    /// use ngeom::pga2d::*;
    /// use ngeom::ops::*;
    /// use ngeom::scalar::*;
    ///
    /// /// f32 wrapper around real numbers that does not implement Recip
    /// #[derive(Ring, Sqrt)]
    /// struct R32(pub f32);
    ///
    /// /// f32 wrapper that doesn't contain a value close to zero
    /// #[derive(Sqrt<Output=R32>, Recip<Output=R32>)]
    /// struct NonZeroR32(R32);
    /// impl NonZeroR32 {
    ///     pub fn new_checked(value: R32, min: R32) -> Option<NonZeroR32> {
    ///         if value.abs() > min {
    ///             NonZeroR32(value)
    ///         } else {
    ///             None
    ///         }
    ///     }
    /// }
    ///
    /// /// Return the unitized point of intersection of two unitized lines,
    /// /// or None if they are parallel
    /// fn try_meet(l1: Bivector<R32>, l2: Bivector<R32>) -> Option<Vector<R32>> {
    ///     const PARALLEL_EPSILON: R32 = 1e-5;
    ///
    ///     let result = l1.meet(l2);
    ///     let norm: R32 = result.weight_norm().into();
    ///     let norm = NonZeroR32.new_checked(norm, PARALLEL_EPSILON)?;
    ///     // Forgetting the above check would cause a compile-time error
    ///     // since R32 doesn't implement Recip
    ///     Some(result * norm.recip())
    /// }
    /// ```
    pub trait Recip: Sized {
        type Output;
        fn recip(self) -> Self::Output;
    }

    pub trait AntiMul<RHS = Self> {
        type Output;
        fn anti_mul(self, rhs: RHS) -> Self::Output;
    }

    pub trait AntiRecip {
        type Output;
        fn anti_recip(self) -> Self::Output;
    }

    pub trait AntiSqrt {
        fn anti_sqrt(self) -> Self;
    }

    pub trait AntiTrig {
        type Output;

        fn anti_cos(self) -> Self::Output;
        fn anti_sin(self) -> Self::Output;
        fn anti_sinc(self) -> Self::Output;
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

            impl Rational for $type {
                fn one_half() -> $type {
                    0.5
                }
            }

            impl Sqrt for $type {
                fn sqrt(self) -> $type {
                    self.sqrt()
                }
            }

            impl Trig for $type {
                type Output = $type;

                fn cos(self) -> $type {
                    self.cos()
                }
                fn sin(self) -> $type {
                    self.sin()
                }
                fn sinc(self) -> $type {
                    let self_adj = self.abs() + $type::EPSILON;
                    self_adj.sin() / self_adj
                }
            }

            impl Recip for $type {
                type Output = $type;

                // This is not NaN-free, even for internal usage!
                fn recip(self) -> $type {
                    self.recip()
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

/// Low-level geometric algebra operations
pub mod algebraic_ops {
    pub trait Reverse {
        type Output;
        fn reverse(self) -> Self;
    }

    pub trait AntiReverse {
        type Output;
        fn anti_reverse(self) -> Self;
    }

    pub trait Bulk {
        type Output;
        fn bulk(self) -> Self::Output;
    }

    pub trait Weight {
        type Output;
        fn weight(self) -> Self::Output;
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

    pub trait BulkContraction<T> {
        type Output;
        fn bulk_contraction(self, r: T) -> Self::Output;
    }

    pub trait WeightContraction<T> {
        type Output;
        fn weight_contraction(self, r: T) -> Self::Output;
    }

    pub trait BulkExpansion<T> {
        type Output;
        fn bulk_expansion(self, r: T) -> Self::Output;
    }

    pub trait WeightExpansion<T> {
        type Output;
        fn weight_expansion(self, r: T) -> Self::Output;
    }

    pub trait AntiCommutator<T> {
        type Output;
        fn anti_commutator(self, r: T) -> Self::Output;
    }

    pub trait ExpAntiWedgeDot {
        type Output;
        fn exp_anti_wedge_dot(self) -> Self::Output;
    }
}

/// Geometric operations
pub mod ops {
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

    pub trait Join<T> {
        type Output;
        fn join(self, r: T) -> Self::Output;
    }

    pub trait Meet<T> {
        type Output;
        fn meet(self, r: T) -> Self::Output;
    }

    pub trait Compose<T> {
        type Output;
        fn compose(self, r: T) -> Self::Output;
    }

    pub trait Normal {
        type Output;
        fn normal(self) -> Self::Output;
    }

    pub trait SubsetOrthogonalTo<T> {
        type Output;
        fn subset_orthogonal_to(self, r: T) -> Self::Output;
    }

    pub trait SupersetOrthogonalTo<T> {
        type Output;
        fn superset_orthogonal_to(self, r: T) -> Self::Output;
    }

    pub trait Projection<T> {
        type Output;
        fn projection(self, r: T) -> Self::Output;
    }

    pub trait AntiProjection<T> {
        type Output;
        fn anti_projection(self, r: T) -> Self::Output;
    }

    pub trait CentralProjection<T> {
        type Output;
        fn central_projection(self, r: T) -> Self::Output;
    }

    pub trait CentralAntiProjection<T> {
        type Output;
        fn central_anti_projection(self, r: T) -> Self::Output;
    }

    pub trait Transform<T> {
        type Output;
        fn transform(self, r: T) -> Self;
    }

    pub trait Normalized {
        type Output;
        fn normalized(self) -> Self::Output;
    }

    pub trait Unitized {
        type Output;
        fn unitized(self) -> Self::Output;
    }

    pub trait MotorTo<T> {
        type Output;
        fn motor_to(self, r: T) -> Self::Output;
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

/// Projective geometry for 2D Euclidean space
pub mod pga2d {
    use crate::algebraic_ops::*;
    use crate::ops::*;
    use crate::scalar::*;
    use ngeom_macros::gen_algebra;

    gen_algebra!(0, 1, 1);

    impl<T> From<T> for AntiScalar<T> {
        fn from(value: T) -> AntiScalar<T> {
            AntiScalar { a012: value }
        }
    }

    macro_rules! impl_for_builtin {
        ($type:ident) => {
            impl From<AntiScalar<$type>> for $type {
                fn from(value: AntiScalar<$type>) -> $type {
                    value.a012
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

    impl<T: Ring + Sqrt + Trig> ExpAntiWedgeDot for Vector<T> {
        type Output = AntiEven<T>;
        fn exp_anti_wedge_dot(self) -> AntiEven<T> {
            // This formula works because a unitized simple vector squares to -1
            // under the anti-wedge-dot product
            // allowing us to treat it like the imaginary unit
            //
            panic!();
            //let theta = self.weight_norm();
            //-self.anti_mul(theta.anti_sinc()) + theta.anti_cos()
        }
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

    pub fn ideal_point<T: Ring>([x, y]: [T; 2]) -> Vector<T> {
        Vector {
            a0: T::zero(),
            a1: x,
            a2: y,
        }
    }

    pub fn homogeneous_point<T>([x, y, w]: [T; 3]) -> Vector<T> {
        Vector {
            a0: w,
            a1: x,
            a2: y,
        }
    }

    /// Line with equation ax + by - c = 0
    pub fn line<T>(a: T, b: T, c: T) -> Bivector<T> {
        Bivector {
            a12: c,
            a20: a,
            a01: b,
        }
    }

    /// Rotor about unitized point p by the given angle
    pub fn axis_angle<T: Ring, A: Ring + Rational + Trig<Output = T>>(
        axis: Vector<T>,
        phi: A,
    ) -> AntiEven<T> {
        let half_phi = phi * A::one_half();
        axis * half_phi.sin() + AntiScalar::from(half_phi.cos())
    }

    /// Rotor about point p by twice its weight
    pub fn rotor<T: Ring + Trig<Output = T>>(p: Vector<T>) -> AntiEven<T> {
        let half_phi = p.weight_norm();
        p.anti_mul(half_phi.anti_sinc()) + half_phi.anti_cos()
    }

    /// Translator towards ideal point p by its magnitude
    pub fn translator<T: Ring + Rational>(p: Vector<T>) -> AntiEven<T> {
        origin().wedge(p).weight_dual() * T::one_half() + AntiScalar::from(T::one())
    }
}

/// Projective geometry for 3D Euclidean space
pub mod pga3d {
    use crate::algebraic_ops::*;
    use crate::ops::*;
    use crate::scalar::*;
    use ngeom_macros::gen_algebra;

    gen_algebra!(0, 1, 1, 1);

    impl<T: Ring> From<T> for AntiScalar<T> {
        fn from(value: T) -> AntiScalar<T> {
            AntiScalar { a0123: value }
        }
    }

    macro_rules! impl_for_builtin {
        ($type:ident) => {
            impl From<AntiScalar<$type>> for $type {
                fn from(value: AntiScalar<$type>) -> $type {
                    value.a0123
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

    impl<T: Ring + Sqrt + Trig> ExpAntiWedgeDot for Bivector<T> {
        type Output = AntiEven<T>;
        fn exp_anti_wedge_dot(self) -> AntiEven<T> {
            // This formula works because a normalized simple vector squares to -1
            // under the anti-wedge-dot product
            // allowing us to treat it like the imaginary unit
            panic!();
            //let theta = self.weight_norm();
            //-self.anti_mul(theta.anti_sinc()) + theta.anti_cos()
        }
    }

    impl<T: Ring + Sqrt + Recip<Output = T>> Magnitude<T> {
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

    pub fn ideal_point<T: Ring>([x, y, z]: [T; 3]) -> Vector<T> {
        Vector {
            a0: T::zero(),
            a1: x,
            a2: y,
            a3: z,
        }
    }

    pub fn homogeneous_point<T>([x, y, z, w]: [T; 4]) -> Vector<T> {
        Vector {
            a0: w,
            a1: x,
            a2: y,
            a3: z,
        }
    }

    /// Plane with equation ax + by + cz - d = 0
    pub fn plane<T>(a: T, b: T, c: T, d: T) -> Trivector<T> {
        Trivector {
            a123: d,
            a032: a,
            a013: b,
            a021: c,
        }
    }

    /// Rotor about unitized line l by the given angle
    pub fn axis_angle<T: Ring, A: Ring + Rational + Trig<Output = T>>(
        axis: Bivector<T>,
        phi: A,
    ) -> AntiEven<T> {
        let half_phi = phi * A::one_half();
        axis * half_phi.sin() + AntiScalar::from(half_phi.cos())
    }

    /// Rotor about line l by its weight
    pub fn rotor<T: Ring + Rational + Trig<Output = T> + Sqrt>(l: Bivector<T>) -> AntiEven<T> {
        let half_phi = l.weight_norm() * T::one_half();
        l.anti_mul(half_phi.anti_sinc()) + half_phi.anti_cos()
    }

    /// Translator towards ideal point p by its magnitude
    pub fn translator<T: Ring + Rational>(p: Vector<T>) -> AntiEven<T> {
        origin().wedge(p).weight_dual() * T::one_half() + AntiScalar::from(T::one())
    }
}

#[cfg(test)]
mod test {
    use crate::ops::*;
    use crate::{pga2d, pga3d};

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

    impl IsClose for pga2d::AntiScalar<f32> {
        fn is_close(self, rhs: pga2d::AntiScalar<f32>) -> bool {
            <f32>::from(self - rhs).abs() < 1e-5
        }
    }

    impl IsClose for pga3d::AntiScalar<f32> {
        fn is_close(self, rhs: pga3d::AntiScalar<f32>) -> bool {
            <f32>::from(self - rhs).abs() < 1e-5
        }
    }

    impl IsClose for pga2d::Vector<f32> {
        fn is_close(self, rhs: Self) -> bool {
            let tol = 1e-5;
            let diff = self - rhs;
            // Taking the right complement allows us to get the weight norm as a scalar
            // rather than an antiscalar
            diff.bulk_norm_squared() < tol * tol && diff.weight_norm_squared() < (tol * tol).into()
        }
    }

    impl IsClose for pga2d::Bivector<f32> {
        fn is_close(self, rhs: Self) -> bool {
            let tol = 1e-5;
            let diff = self - rhs;
            diff.bulk_norm_squared() < tol * tol && diff.weight_norm_squared() < (tol * tol).into()
        }
    }

    impl IsClose for pga3d::Vector<f32> {
        fn is_close(self, rhs: Self) -> bool {
            let tol = 1e-5;
            let diff = self - rhs;
            diff.bulk_norm_squared() < tol * tol && diff.weight_norm_squared() < (tol * tol).into()
        }
    }

    impl IsClose for pga3d::Bivector<f32> {
        fn is_close(self, rhs: Self) -> bool {
            let tol = 1e-5;
            let diff = self - rhs;
            diff.bulk_norm_squared() < tol * tol && diff.weight_norm_squared() < (tol * tol).into()
        }
    }

    impl IsClose for pga3d::Trivector<f32> {
        fn is_close(self, rhs: Self) -> bool {
            let tol = 1e-5;
            let diff = self - rhs;
            diff.bulk_norm_squared() < tol * tol && diff.weight_norm_squared() < (tol * tol).into()
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

            assert_close!(p1.join(p2).weight_norm(), (5.).into());
        }

        // Signed distance between point & line
        // Joining the line & point results in a plane (scalar) corresponding to their signed distance
        {
            // Note that the line and the test point form a triangle that is wound
            // left-handed, producing a negative distance
            let l = pga2d::point::<f32>([10., 10.])
                .join(pga2d::point([10., 20.]))
                .unitized();
            let p = pga2d::point([15., 15.]);
            assert_close!(l.join(p).weight_norm(), (-5.).into());
        }

        // Distance between parallel lines
        // The lines meet at an infinite point whose infinite norm is their unsigned distance.
        // More precisely, this expression yields the distance between the projection of the
        // origin onto each line.
        {
            let l1 = pga2d::origin::<f32>()
                .join(pga2d::point([3., 4.]))
                .unitized();
            let l2 = pga2d::point::<f32>([4., -3.])
                .join(pga2d::point([7., 1.]))
                .unitized();
            assert_close!(dbg!(l1.meet(l2)).bulk_norm(), 5.);
        }
    }

    #[test]
    fn test_dist_3d() {
        // Distance between points
        // Joining the points results in a line whose norm is their unsigned distance
        {
            let p1 = pga3d::point::<f32>([10., 10., 10.]);
            let p2 = pga3d::point([13., 14., 10.]);
            assert_close!(p1.join(p2).weight_norm(), (5.).into());
        }

        // Distnce between point & line
        // Joining the line & point results in a plane whose norm is their unsigned distance
        {
            let l = pga3d::point::<f32>([10., 10., 10.])
                .join(pga3d::point([10., 20., 10.]))
                .unitized();
            let p = pga3d::point([15., 15., 10.]);
            assert_close!(l.join(p).weight_norm(), (5.).into());
        }

        {
            let l = pga3d::point::<f32>([10., 10., 0.])
                .join(pga3d::point([10., 10., 20.]))
                .unitized();
            let p = pga3d::point([13., 14., 10.]);
            assert_close!(l.join(p).weight_norm(), (5.).into());
        }

        // Distance between point & plane
        // Joining the plane & point results in a scalar corresponding to their signed distance
        {
            let pl = pga3d::point::<f32>([10., 10., 10.])
                .join(pga3d::point([20., 10., 10.]))
                .join(pga3d::point([10., 20., 10.]))
                .unitized();
            let pt = pga3d::point([15., 15., 15.]);
            assert_close!(pl.join(pt), (5.).into());
        }

        // Distance between parallel lines
        // More precisely, this expression yields the distance between the projection of the
        // origin onto each line.
        {
            use super::algebraic_ops::AntiCommutator;
            let l1 = pga3d::origin::<f32>()
                .join(pga3d::point([0., 0., 10.]))
                .unitized();
            let l2 = pga3d::point::<f32>([3., 4., 10.])
                .join(pga3d::point([3., 4., 20.]))
                .unitized();
            assert_close!(l1.anti_commutator(l2).bulk_norm(), 5.);
        }

        // Distance between perpendicular lines
        // More precisely, this expression yields d * sin(a)
        // d is the distance between the lines, and a is the angle between them
        // (as measured along / about their shared normal)
        {
            let _l1 = pga3d::origin::<f32>()
                .join(pga3d::point([0., 0., 10.]))
                .unitized();
            let _l2 = pga3d::point::<f32>([8., 0., 15.])
                .join(pga3d::point([8., 20., 15.]))
                .unitized();
            //assert_close!(l1.join(l2), pga3d::anti(-8.)); // XXX broken
        }

        // Distance between skew lines
        // This expression divides the join of the two lines (d * sin(a))
        // by the scalar norm of their commutator product, which is sin(a)
        {
            let _l1 = pga3d::origin::<f32>()
                .join(pga3d::point([0., 0., 10.]))
                .unitized();
            let _l2 = pga3d::point([10., 0., 0.])
                .join(pga3d::point([10., 5., 4.]))
                .unitized();

            //assert_close!(l1.join(l2) + l1.anti_commutator(l2).weight_norm(), -10.) // XXX broken
        }
    }

    #[test]
    fn test_sign_2d_triangle_winding() {
        let p1 = pga2d::origin::<f32>();
        let p2 = pga2d::point([1., 0.]);
        let p3 = pga2d::point([0., 1.]);
        assert!(p1.join(p2).join(p3) > (0.).into());
        assert!(p3.join(p2).join(p1) < (0.).into());
    }

    #[test]
    fn test_sign_2d_line_intersection() {
        let l1 = pga2d::origin::<f32>()
            .join(pga2d::point([10., 0.]))
            .unitized();
        let l2 = pga2d::point([0., 10.])
            .join(pga2d::point([10., 9.]))
            .unitized();
        let l3 = pga2d::point([0., -10.])
            .join(pga2d::point([-10., -9.]))
            .unitized();

        assert!(l1.meet(l2).weight_norm() < (0.).into());
        assert!(l2.meet(l1).weight_norm() > (0.).into());
        assert!(l1.meet(l3).weight_norm() > (0.).into());
        assert!(l3.meet(l1).weight_norm() < (0.).into());
    }

    #[test]
    fn test_sign_3d_tetrahedron_winding() {
        let p1 = pga3d::origin::<f32>();
        let p2 = pga3d::point([1., 0., 0.]);
        let p3 = pga3d::point([0., 1., 0.]);
        let p4 = pga3d::point([0., 0., 1.]);
        assert!(p1.join(p2).join(p3).join(p4) > (0.).into());
        assert!(p3.join(p2).join(p1).join(p4) < (0.).into());
    }

    #[test]
    fn test_sign_3d_skew_lines() {
        // Note that this is intentionally backwards!
        // It's impossible to reconcile the sign of this
        // with the sign of the plane-point version,
        // so we add the negative sign here
        let l1 = pga3d::origin::<f32>()
            .join(pga3d::point([0., 0., 10.]))
            .unitized();
        let l2 = pga3d::point::<f32>([8., 0., 15.])
            .join(pga3d::point([8., 20., 15.]))
            .unitized();
        let l3 = pga3d::point::<f32>([-10., 0., 0.])
            .join(pga3d::point([-10., 20., -5.]))
            .unitized();

        assert!(l1.join(l2) > (0.).into());
        assert!(l1.join(l3) < (0.).into());
    }

    #[test]
    fn test_sign_3d_plane_intersection() {
        let p1 = pga3d::origin::<f32>()
            .join(pga3d::point([1., 0., 0.]))
            .join(pga3d::point([0., 1., 0.]))
            .unitized();
        let p2 = pga3d::origin::<f32>()
            .join(pga3d::point([0., 1., 0.]))
            .join(pga3d::point([0., 0., 1.]))
            .unitized();
        let p3 = pga3d::origin::<f32>()
            .join(pga3d::point([0., 0., 1.]))
            .join(pga3d::point([1., 0., 0.]))
            .unitized();

        dbg!(p1.meet(p2).meet(p3).weight_norm() > (0.).into());
    }

    #[test]
    fn test_2d_line_normal() {
        let l = pga2d::origin::<f32>()
            .join(pga2d::point([1., 0.]))
            .unitized();
        let expected_dir = pga2d::ideal_point([0., 1.]);

        let d = l.normal();

        assert!((d - expected_dir).weight_norm_squared() < (1e-5 * 1e-5).into());
    }

    #[test]
    fn test_plane_normal() {
        let p = pga3d::origin::<f32>()
            .join(pga3d::point([1., 0., 0.]))
            .join(pga3d::point([0., 1., 0.]))
            .unitized();

        let expected_dir = pga3d::ideal_point([0., 0., 1.]);

        let d = p.normal();

        assert!((d - expected_dir).weight_norm_squared() < (1e-5 * 1e-5).into());
    }

    #[test]
    fn test_3d_line_normal() {
        let l = pga3d::origin::<f32>()
            .join(pga3d::point([0., 1., 0.]))
            .unitized();

        let inf_line = l.normal();

        let expected_inf_line =
            pga3d::ideal_point([0., 0., 1.]).join(pga3d::ideal_point([1., 0., 0.]));

        assert!((inf_line - expected_inf_line).weight_norm_squared() < (1e-5 * 1e-5).into());
    }

    #[test]
    fn test_2d_rotation() {
        let p = pga2d::point([1., 0.]);
        let motor = pga2d::axis_angle(pga2d::origin(), 0.25 * core::f32::consts::TAU);
        assert_close!(p.transform(motor), pga2d::point([0., 1.]));
    }

    #[test]
    fn test_3d_rotation() {
        let p = pga3d::point([1., 0., 0.]);
        let l = pga3d::origin::<f32>()
            .join(pga3d::point([0., 0., 1.]))
            .unitized();
        let motor = pga3d::axis_angle(l, 0.25 * core::f32::consts::TAU);
        assert_close!(p.transform(motor), pga3d::point([0., 1., 0.]));
    }

    #[test]
    fn test_2d_translation() {
        let p1 = pga2d::point([10., 10.]);

        let motor = pga2d::translator(pga2d::ideal_point([0., 5.]));

        let p2 = p1.transform(motor);
        assert_close!(p2, pga2d::point([10., 15.]));
    }

    #[test]
    fn test_3d_translation() {
        let p1 = pga3d::point([10., 10., 10.]);

        let motor = pga3d::translator(pga3d::ideal_point([0., 5., 0.]));

        let p2 = p1.transform(motor);
        assert_close!(p2, pga3d::point([10., 15., 10.]));
    }

    #[test]
    fn test_2d_translation_from_geo() {
        let p1 = pga2d::point([10., 10.]);

        let motor = pga2d::point([50., 50.]).motor_to(pga2d::point([50., 52.5]));

        let p2 = p1.transform(motor);
        assert_close!(p2, pga2d::point([10., 15.]));
    }

    #[test]
    fn test_2d_rotation_from_geo() {
        let p = pga2d::point([1., 0.]);

        let rot_l1 = pga2d::origin().join(pga2d::point([10., 0.])).unitized();
        let rot_l2 = pga2d::origin().join(pga2d::point([10., 10.])).unitized();
        let motor = rot_l1.motor_to(rot_l2);

        assert_close!(p.transform(motor), pga2d::point([0., 1.]));
    }

    #[test]
    fn test_2d_compose_rotations() {
        let p = pga2d::point([10., 10.]);

        let translate_up_5 = pga2d::translator(pga2d::ideal_point([0., 5.]));
        let rotate_90 = pga2d::axis_angle(pga2d::origin(), 0.25 * core::f32::consts::TAU);

        let up_then_rotate = translate_up_5.compose(rotate_90);
        assert_close!(p.transform(up_then_rotate), pga2d::point([-15., 10.]));

        let rotate_then_up = rotate_90.compose(translate_up_5);
        assert_close!(p.transform(rotate_then_up), pga2d::point([-10., 15.]));
    }

    #[test]
    fn test_2d_motor_between_parallel_lines() {
        let rot_l1 = pga2d::origin().join(pga2d::point([10., 0.])).unitized();
        let rot_l2 = pga2d::point([0., 10.])
            .join(pga2d::point([10., 10.]))
            .unitized();
        let motor = rot_l1.motor_to(rot_l2);

        let p = pga2d::origin();
        assert_close!(p.transform(motor), pga2d::point([0., 20.]));
    }

    #[test]
    fn test_2d_superset_orthogonal_to() {
        let l = pga2d::origin().join(pga2d::point([1., 1.])).unitized();
        let p = pga2d::point([0., 2.]);
        let l2 = p.superset_orthogonal_to(l);

        let expected_l2 = pga2d::point([1., 1.]).join(p).unitized();

        assert_close!(l2, expected_l2);
    }

    #[test]
    fn test_3d_superset_orthogonal_to() {
        let pl = pga3d::origin()
            .join(pga3d::point([0., 0., 10.]))
            .join(pga3d::point([1., 1., 10.]))
            .unitized();
        let p = pga3d::point([0., 2., 10.]);
        let l = p.superset_orthogonal_to(pl);

        let expected_l = pga3d::point([1., 1., 10.]).join(p).unitized();

        assert_close!(l, expected_l);
    }

    #[test]
    fn test_2d_projection() {
        let l = pga2d::origin().join(pga2d::point([1., 1.])).unitized();
        let p = pga2d::point([0., 2.]);
        let p2 = p.projection(l).unitized();

        assert_close!(p2, pga2d::point([1., 1.]));
    }

    #[test]
    fn test_3d_projection_pt_onto_plane() {
        let pl = pga3d::origin()
            .join(pga3d::point([0., 0., 10.]))
            .join(pga3d::point([1., 1., 10.]))
            .unitized();
        let p = pga3d::point([0., 2., 10.]);
        let p2 = p.projection(pl).unitized();

        assert_close!(p2, pga3d::point([1., 1., 10.]));
    }

    #[test]
    fn test_2d_anti_projection() {
        let l = pga2d::origin().join(pga2d::point([1., 1.])).unitized();
        let p = pga2d::point([0., 2.]);
        let l2 = l.anti_projection(p).unitized();

        let expected_l2 = p.join(pga2d::ideal_point([1., 1.])).unitized();
        assert_close!(l2, expected_l2);
    }

    #[test]
    fn test_3d_anti_projection_plane_onto_pt() {
        let pl = pga3d::origin()
            .join(pga3d::point([0., 0., 10.]))
            .join(pga3d::point([1., 1., 10.]))
            .unitized();
        let p = pga3d::point([0., 2., 10.]);
        let pl2 = pl.anti_projection(p).unitized();

        let expected_pl2 = p
            .join(pga3d::ideal_point([0., 0., 1.]))
            .join(pga3d::ideal_point([1., 1., 10.]))
            .unitized();

        assert_close!(pl2, expected_pl2);
    }

    #[test]
    fn test_3d_projection_line_onto_plane() {
        let pl = pga3d::point([0., 0., 10.])
            .join(pga3d::point([1., 0., 10.]))
            .join(pga3d::point([0., 1., 10.]))
            .unitized();
        let l = pga3d::point([0., 20., 20.])
            .join(pga3d::point([1., 20., 20.]))
            .unitized();

        let l2 = l.projection(pl).unitized();

        let expected_l2 = pga3d::point([0., 20., 10.])
            .join(pga3d::point([1., 20., 10.]))
            .unitized();

        assert_close!(l2, expected_l2);
    }

    #[test]
    fn test_3d_central_projection_line_onto_plane() {
        let pl = pga3d::point([0., 0., 10.])
            .join(pga3d::point([1., 0., 10.]))
            .join(pga3d::point([0., 1., 10.]))
            .unitized();
        let l = pga3d::point([0., 20., 20.])
            .join(pga3d::point([1., 20., 20.]))
            .unitized();

        let l2 = l.central_projection(pl).unitized();

        let expected_l2 = pga3d::point([0., 10., 10.])
            .join(pga3d::point([1., 10., 10.]))
            .unitized();

        assert_close!(l2, expected_l2);
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
                pl.join(l).norm().atan2(pl.dot(l).norm()),
                0.125 * core::f32::consts::TAU,
            )
        }
    }

    #[test]
    fn motor() {
        // 2D translation
        let p1 = pga2d::point([10., 10.]);

        let center = pga2d::ideal_point([1., 0.]) * -2.5;
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
