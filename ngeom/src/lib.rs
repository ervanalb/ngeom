#![cfg_attr(not(feature = "std"), no_std)]

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
//! ngeom is `no_std`-compatible, with most functionality available,
//! and the option to implement the rest.
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

    /// A scalar datatype whose absolute value can be taken.
    pub trait Abs {
        type Output;

        /// Computes the absolute value of a scalar.
        fn abs(self) -> Self::Output;
    }

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
        + Abs<Output = Self>
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
    /// meaning its use within the library is NaN-free.
    pub trait Sqrt {
        // This scalar's positive square root
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
    /// For floating point datatypes, this operation is NOT NaN-free,
    /// including in its use within the library.
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
    /// ```
    /// use ngeom::{impl_ops_for_scalar, impl_pga2d_for_scalar};
    /// use ngeom::pga2d::*;
    /// use ngeom::ops::*;
    /// use ngeom::scalar::*;
    ///
    /// /// f32 wrapper around real numbers
    /// /// that does not implement `Recip`
    /// // TODO write derive macro for Abs, Ring, Sqrt, etc.
    /// #[derive(Clone, Copy, PartialEq, PartialOrd)]
    /// struct R32(pub f32);
    /// impl core::ops::Add<R32> for R32 {
    ///     type Output = R32;
    ///     fn add(self, r: R32) -> R32 { R32(self.0 + r.0) }
    /// }
    /// impl core::ops::Sub<R32> for R32 {
    ///     type Output = R32;
    ///     fn sub(self, r: R32) -> R32 { R32(self.0 - r.0) }
    /// }
    /// impl core::ops::Mul<R32> for R32 {
    ///     type Output = R32;
    ///     fn mul(self, r: R32) -> R32 { R32(self.0 * r.0) }
    /// }
    /// impl core::ops::Neg for R32 {
    ///     type Output = R32;
    ///     fn neg(self) -> R32 { R32(-self.0) }
    /// }
    /// impl Abs for R32 {
    ///     type Output = R32;
    ///     fn abs(self) -> R32 { R32(self.0.abs()) }
    /// }
    /// impl Ring for R32 {
    ///     fn zero() -> R32 { R32(0.) }
    ///     fn one() -> R32 { R32(1.) }
    /// }
    /// impl Sqrt for R32 {
    ///     fn sqrt(self) -> R32 { R32(self.0.sqrt()) }
    /// }
    /// impl_ops_for_scalar!(R32); // TODO replace with derive?
    /// impl_pga2d_for_scalar!(R32);
    ///
    /// /// f32 wrapper around real numbers that aren't close to zero
    /// /// allowing for NaN-free `Recip`
    /// struct NonZeroR32(f32);
    /// impl NonZeroR32 {
    ///     pub fn new_checked(value: R32, min: R32) -> Option<NonZeroR32> {
    ///         if value.abs() > min {
    ///             Some(NonZeroR32(value.0))
    ///         } else {
    ///             None
    ///         }
    ///     }
    /// }
    /// impl Recip for NonZeroR32 {
    ///     type Output = R32;
    ///     fn recip(self) -> R32 {
    ///         R32(self.0.recip())
    ///     }
    /// }
    ///
    /// /// Return the unitized point of intersection of two unitized lines,
    /// /// or None if they are parallel
    /// fn try_meet(l1: Bivector<R32>, l2: Bivector<R32>) -> Option<Vector<R32>> {
    ///     const PARALLEL_EPSILON: R32 = R32(1e-5);
    ///
    ///     let result = l1.meet(l2);
    ///     let norm: R32 = result.weight_norm().into();
    ///     let norm = NonZeroR32::new_checked(norm, PARALLEL_EPSILON)?;
    ///     // Forgetting the above check would cause a compile-time error
    ///     // since R32 doesn't implement Recip
    ///     Some(result * norm.recip())
    /// }
    /// ```
    pub trait Recip: Sized {
        type Output;
        fn recip(self) -> Self::Output;
    }

    /// The dual of [Abs] for anti-scalars.
    ///
    /// e.g. (-𝟜).anti_abs() = 𝟜
    ///
    /// This will automatically be implemented for anti-scalars
    /// when their backing scalar type implements `Abs`.
    pub trait AntiAbs {
        type Output;
        fn anti_abs(self) -> Self::Output;
    }

    /// Anti-scalar multiplication
    ///
    /// e.g. 𝟜.anti_mul(𝟚) = 𝟠
    ///
    /// Note the distinction between this and scalar multiplication:
    /// * 𝟜 * 2 = 𝟠
    /// * 4 * 𝟚 = 𝟠
    ///
    /// This will automatically be implemented for anti-scalars
    /// when their backing scalar type implements `Mul`.
    pub trait AntiMul<RHS = Self> {
        type Output;
        fn anti_mul(self, rhs: RHS) -> Self::Output;
    }

    /// The dual of [Recip] for anti-scalars.
    ///
    /// e.g. 𝟚.anti_recip() = ½𝟙
    ///
    /// This will automatically be implemented for anti-scalars
    /// when their backing scalar type implements `Recip`.
    pub trait AntiRecip {
        type Output;
        fn anti_recip(self) -> Self::Output;
    }

    /// The dual of [Sqrt] for anti-scalars.
    ///
    /// e.g. 𝟜.anti_sqrt() = 𝟚
    ///
    /// This will automatically be implemented for anti-scalars
    /// when their backing scalar type implements `Sqrt`.
    pub trait AntiSqrt {
        fn anti_sqrt(self) -> Self;
    }

    /// The dual of [Trig] for anti-scalars.
    ///
    /// e.g. (⅓ π𝟙).anti_cos() = ½𝟙
    ///
    /// This will automatically be implemented for anti-scalars
    /// when their backing scalar type implements `Trig`.
    pub trait AntiTrig {
        type Output;

        fn anti_cos(self) -> Self::Output;
        fn anti_sin(self) -> Self::Output;
        fn anti_sinc(self) -> Self::Output;
    }

    macro_rules! impl_for_float {
        ($type:ident) => {
            #[cfg(feature = "std")]
            impl Abs for $type {
                type Output = $type;
                fn abs(self) -> $type {
                    self.abs()
                }
            }
            #[cfg(not(feature = "std"))]
            impl Abs for $type {
                type Output = $type;
                fn abs(self) -> $type {
                    if self < 0. {
                        -self
                    } else {
                        self
                    }
                }
            }

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

            #[cfg(feature = "std")]
            impl Sqrt for $type {
                fn sqrt(self) -> $type {
                    self.sqrt()
                }
            }

            #[cfg(feature = "std")]
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
            impl Abs for $type {
                type Output = $type;
                fn abs(self) -> $type {
                    self.abs()
                }
            }

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
///
/// There are a few conflicting conventions for geometric algebra notation.
/// This library tends to use the formulation put forth by
/// [Dr. Eric Lengyel](https://projectivegeometricalgebra.org/)
pub mod algebraic_ops {
    /// The reverse operator Ã
    ///
    /// see <https://rigidgeometricalgebra.org/wiki/index.php?title=Reverses>
    pub trait Reverse {
        type Output;
        fn reverse(self) -> Self;
    }

    /// The anti-reverse operator A̰
    ///
    /// see <https://rigidgeometricalgebra.org/wiki/index.php?title=Reverses>
    pub trait AntiReverse {
        type Output;
        fn anti_reverse(self) -> Self;
    }

    /// The bulk operator A●
    ///
    /// see <https://rigidgeometricalgebra.org/wiki/index.php?title=Bulk_and_weight>
    pub trait Bulk {
        type Output;
        fn bulk(self) -> Self::Output;
    }

    /// The weight operator A○
    ///
    /// see <https://rigidgeometricalgebra.org/wiki/index.php?title=Bulk_and_weight>
    pub trait Weight {
        type Output;
        fn weight(self) -> Self::Output;
    }

    /// The bulk dual operator A★
    ///
    /// See <https://rigidgeometricalgebra.org/wiki/index.php?title=Duals>
    pub trait BulkDual {
        type Output;
        fn bulk_dual(self) -> Self::Output;
    }

    /// The weight dual operator A☆
    ///
    /// See <https://rigidgeometricalgebra.org/wiki/index.php?title=Duals>
    pub trait WeightDual {
        type Output;
        fn weight_dual(self) -> Self::Output;
    }

    /// The right complement operator A̅
    ///
    /// See <https://rigidgeometricalgebra.org/wiki/index.php?title=Complements>
    pub trait RightComplement {
        type Output;
        fn right_complement(self) -> Self::Output;
    }

    /// The left complement operator A̲
    ///
    /// See <https://rigidgeometricalgebra.org/wiki/index.php?title=Complements>
    pub trait LeftComplement {
        type Output;
        fn left_complement(self) -> Self::Output;
    }

    /// The wedge product from exterior algebra, A ∧ B
    ///
    /// See <https://rigidgeometricalgebra.org/wiki/index.php?title=Exterior_products>
    pub trait Wedge<T> {
        type Output;
        fn wedge(self, r: T) -> Self::Output;
    }

    /// The anti-wedge product (vee) from exterior algebra, A ∨ B
    ///
    /// See <https://rigidgeometricalgebra.org/wiki/index.php?title=Exterior_products>
    pub trait AntiWedge<T> {
        type Output;
        fn anti_wedge(self, r: T) -> Self::Output;
    }

    /// The dot product A • B
    ///
    /// See <https://rigidgeometricalgebra.org/wiki/index.php?title=Dot_products>
    pub trait Dot<T> {
        type Output;
        fn dot(self, r: T) -> Self::Output;
    }

    /// The anti-dot product A ∘ B
    ///
    /// See <https://rigidgeometricalgebra.org/wiki/index.php?title=Dot_products>
    pub trait AntiDot<T> {
        type Output;
        fn anti_dot(self, r: T) -> Self::Output;
    }

    /// The geometric product A ⟑ B
    ///
    /// See <https://rigidgeometricalgebra.org/wiki/index.php?title=Geometric_products>
    pub trait WedgeDot<T> {
        type Output;
        fn wedge_dot(self, r: T) -> Self::Output;
    }

    /// The geometric anti-product A ⟇ B
    ///
    /// See <https://rigidgeometricalgebra.org/wiki/index.php?title=Geometric_products>
    pub trait AntiWedgeDot<T> {
        type Output;
        fn anti_wedge_dot(self, r: T) -> Self::Output;
    }

    /// The bulk contraction operator A ∨ B★
    ///
    /// See <https://rigidgeometricalgebra.org/wiki/index.php?title=Interior_products>
    pub trait BulkContraction<T> {
        type Output;
        fn bulk_contraction(self, r: T) -> Self::Output;
    }

    /// The weight contraction operator A ∨ B☆
    ///
    /// See <https://rigidgeometricalgebra.org/wiki/index.php?title=Interior_products>
    pub trait WeightContraction<T> {
        type Output;
        fn weight_contraction(self, r: T) -> Self::Output;
    }

    /// The bulk expansion operator A ∧ B★
    ///
    /// See <https://rigidgeometricalgebra.org/wiki/index.php?title=Interior_products>
    pub trait BulkExpansion<T> {
        type Output;
        fn bulk_expansion(self, r: T) -> Self::Output;
    }

    /// The weight expansion operator A ∧ B☆
    ///
    /// See <https://rigidgeometricalgebra.org/wiki/index.php?title=Interior_products>
    pub trait WeightExpansion<T> {
        type Output;
        fn weight_expansion(self, r: T) -> Self::Output;
    }

    /// The commutator anti-product ½(A ⟇ B - B ⟇ A)
    ///
    /// See <https://rigidgeometricalgebra.org/wiki/index.php?title=Commutators>
    pub trait AntiCommutator<T> {
        type Output;
        fn anti_commutator(self, r: T) -> Self::Output;
    }

    /// The exponential operator e^A under the geometric anti-product ⟇
    ///
    /// See <https://rigidgeometricalgebra.org/wiki/index.php?title=Motor#Exponential_Form>
    pub trait ExpAntiWedgeDot {
        type Output;
        fn exp_anti_wedge_dot(self) -> Self::Output;
    }
}

/// Geometric operations
pub mod ops {
    /// The square of the bulk norm of an element
    ///
    /// See [BulkNorm] for more details.
    ///
    /// This trait avoids the square root which is sometimes necessary to calculate [BulkNorm]
    /// and is therefore always available,
    /// even on scalar types that do not implement [Sqrt](crate::scalar::Sqrt).
    pub trait BulkNormSquared {
        type Output;
        fn bulk_norm_squared(self) -> Self::Output;
    }

    /// The square of the weight norm of an element
    ///
    /// See [WeightNorm] for more details.
    ///
    /// This trait avoids the square root which is sometimes necessary to calculate [WeightNorm]
    /// and is therefore always available,
    /// even on scalar types that do not implement [Sqrt](crate::scalar::Sqrt).
    pub trait WeightNormSquared {
        type Output;
        fn weight_norm_squared(self) -> Self::Output;
    }

    /// The bulk norm of an element
    ///
    /// Thanks to homogeneous coordinates,
    /// all elements can be multiplied by a non-zero scalar
    /// without changing the geometry they represent.
    /// In other words, `ideal_point([3, 4])` represents the same direction as
    /// `ideal_point([3, 4]) * 2`.
    /// In the case of ideal points, the bulk norm gives their "length":
    ///
    /// ```
    /// use ngeom::pga2d::*;
    /// use ngeom::ops::*;
    /// use ngeom::scalar::*;
    ///
    /// let p1 = ideal_point([3., 4.]);
    /// let p2 = p1 * 2.;
    /// assert_eq!(p1.bulk_norm(), (5.).into());
    /// assert_eq!(p2.bulk_norm(), (10.).into());
    /// ```
    ///
    /// Note that ideal elements have a [weight norm](WeightNorm) of zero.
    ///
    /// Here is another example of using `bulk_norm()`
    /// to extract information from an ideal point:
    ///
    /// ```
    /// use ngeom::pga2d::*;
    /// use ngeom::ops::*;
    ///
    /// // Construct two unitized parallel lines
    /// let l1 = point([0., 0.])
    ///     .join(point([10., 0.]))
    ///     .unitized();
    /// let l2 = point([0., 10.])
    ///     .join(point([10., 10.]))
    ///     .unitized();
    ///
    /// // Intersect the two parallel lines to get a point at infinity
    /// let p = l1.meet(l2);
    /// assert_eq!(p.weight_norm(), (0.).into());
    ///
    /// // The bulk norm of the intersection point is the distance between the lines
    /// assert_eq!(p.bulk_norm(), 10.);
    /// ```
    ///
    /// Depending on the element, this trait may or may not require taking a square root.
    /// Consider [BulkNormSquared] if this is an issue.
    pub trait BulkNorm {
        type Output;
        fn bulk_norm(self) -> Self::Output;
    }

    /// The weight norm of an element
    ///
    /// Thanks to homogeneous coordinates,
    /// all elements can be multiplied by a non-zero scalar
    /// without changing the geometry they represent.
    /// In other words, `point([2, 3])` is the same location in space as
    /// `point([2, 3]) * 5`. This extra projective factor
    /// can be retrieved by the function `.weight_norm()`:
    ///
    /// ```
    /// use ngeom::pga2d::*;
    /// use ngeom::ops::*;
    ///
    /// let p1 = point([2, 3]);
    /// let p2 = p1 * 5;
    /// assert_eq!(p1.weight_norm(), (1).into());
    /// assert_eq!(p2.weight_norm(), (5).into());
    /// ```
    ///
    /// Note that ideal elements have a weight norm of zero,
    /// but similar information can be extracted using the [bulk norm](BulkNorm).
    ///
    /// The weight norm is an anti-scalar since its sign changes under reflection.
    /// It may be cast into a scalar for the purposes of arithmetic using `.into()`.
    ///
    /// Many functions require that an element's weight norm be equal to 𝟙,
    /// i.e. require the element to be [unitized](Unitized).
    ///
    /// ```
    /// use ngeom::pga2d::*;
    /// use ngeom::ops::*;
    ///
    /// // Join two points into a line
    /// let p1 = point([10., 10.]);
    /// let p2 = point([13., 14.]);
    /// let l = p1.join(p2);
    ///
    /// // The weight norm of the line is the distance between the points
    /// assert_eq!(l.weight_norm(), (5.).into());
    /// ```
    ///
    /// Depending on the element, this trait may or may not require taking a square root.
    /// Consider [WeightNormSquared] if this is an issue.
    pub trait WeightNorm {
        type Output;
        fn weight_norm(self) -> Self::Output;
    }

    /// The higher-dimensional geometry containing its two operands, similar to a union.
    ///
    /// For example, this function will join two points into a line,
    /// or a point and line into a plane.
    ///
    /// The order of the operands will generally not affect the location of the resultant geometry,
    /// but may affect the sign of its directionality,
    /// according to the winding order in which it was built.
    ///
    /// ```
    /// use ngeom::pga3d::*;
    /// use ngeom::ops::*;
    ///
    /// // Join two points into a line
    /// // The direction of the line goes from p1 to p2,
    /// // towards +X in this case
    /// let p1 = point([0., 0., 0.]);
    /// let p2 = point([10., 0., 0.]);
    /// let l = p1.join(p2);
    ///
    /// // Join a line and a point into a plane
    /// // The direction of the plane goes from l to p3,
    /// // which produces the +XY plane in this case
    /// let p3 = point([5., 5., 0.]);
    /// let pl = l.join(p3);
    ///
    /// // Join a plane and a point into a volume
    /// // The direction of the volume is oriented according to the right hand rule.
    /// let p4 = point([3., 2., 15.]);
    /// let v = pl.join(p4);
    ///
    /// // Since we built the volume right-handed,
    /// // the resulting value should be positive
    /// assert!(v > (0.).into());
    ///
    /// // Reversing the final product produces a left-handed volume instead
    /// let v_rev = p4.join(pl);
    /// assert!(v_rev < (0.).into());
    /// ```
    ///
    /// `Join` is exception-free.
    /// Joining coincident geometry will result in an ideal (infinite) element.
    ///
    /// ```
    /// use ngeom::pga3d::*;
    /// use ngeom::ops::*;
    ///
    /// // Join two coincident points
    /// let p1 = point([4., 5., 6.]);
    /// let p2 = point([4., 5., 6.]) * 2.; // Homogeneous scaling doesn't matter
    /// let l = p1.join(p2);
    ///
    /// assert_eq!(l.weight_norm(), (0.).into())
    /// ```
    pub trait Join<T> {
        type Output;
        fn join(self, r: T) -> Self::Output;
    }

    /// The lower-dimensional geometry shared between its two operands, i.e. intersection.
    ///
    /// For example, in 2D, this function will meet two lines at a point.
    /// In 3D, it will meet a plane and a line at a point,
    /// or two planes at a line.
    ///
    /// The order of the operands will generally not affect the location of the resultant geometry,
    /// but may affect the sign of its directionality,
    /// according to the winding order in which it was built.
    ///
    /// ```
    /// use ngeom::pga3d::*;
    /// use ngeom::ops::*;
    ///
    /// // Join three points to get the XY plane
    /// let xy = point([0., 0., 0.])
    ///     .join(point([1., 0., 0.]))
    ///     .join(point([0., 1., 0.]));
    ///
    /// // Join two points to get a line travelling in the -Z direction
    /// let l = point([3., 4., 10.])
    ///     .join(point([3., 4., 9.]));
    ///
    /// // Meet the plane and the line at a point
    /// let pt = xy.meet(l);
    /// assert_eq!(pt.unitized(), point([3., 4., 0.]));
    /// ```
    ///
    /// `Meet` is exception-free.
    /// Meeting parallel geometry will result in an ideal (infinite) element,
    /// which may contain useful information such as the parallel distance.
    ///
    /// ```
    /// use ngeom::pga2d::*;
    /// use ngeom::ops::*;
    ///
    /// // Construct two unitized parallel lines
    /// let l1 = point([0., 0.])
    ///     .join(point([10., 0.]))
    ///     .unitized();
    /// let l2 = point([0., 10.])
    ///     .join(point([10., 10.]))
    ///     .unitized();
    ///
    /// // Intersect the two parallel lines to get a point at infinity
    /// let p = l1.meet(l2);
    /// assert_eq!(p.weight_norm(), (0.).into());
    ///
    /// // The bulk norm of the intersection point is the distance between the lines
    /// assert_eq!(p.bulk_norm(), 10.);
    /// ```
    pub trait Meet<T> {
        type Output;
        fn meet(self, r: T) -> Self::Output;
    }

    /// Compose motors A and B into a new motor whose motion is the result of applying A then B (extrinsically)
    ///
    /// Note that motors are composed left-to-right.
    /// This is the opposite convention of quaternions or matrices, which compose right-to-left.
    pub trait Compose<T> {
        type Output;
        fn compose(self, r: T) -> Self::Output;
    }

    /// Get the ideal element orthogonal to the given element.
    ///
    /// For example:
    /// * Given a line in 2D, get the ideal point orthogonal to the line
    ///   (on the positive side according to the right-hand rule)
    /// * Given a plane in 3D, get the ideal point orthogonal to the plane
    ///   (on the positive side according to the right-hand rule)
    /// * Given a line in 3D, get the ideal line orthogonal to it
    ///   (wrapping in the positive direction according to the right hand rule)
    pub trait Normal {
        type Output;
        fn normal(self) -> Self::Output;
    }

    /// Retrieve the lower-dimensional geometry that is contained within A and is orthogonal to B.
    ///
    /// This operation generally returns ideal (infinite) elements.
    ///
    /// For example, in 3D, given a plane and a line, it can find the direction (ideal point) within that plane
    /// that is orthogonal to a given line.
    ///
    /// This operation is an intermediate step in computing the [AntiProjection].
    /// [Joining](Join) B with the result produces the anti-projection of A onto B.
    pub trait SubsetOrthogonalTo<T> {
        type Output;
        fn subset_orthogonal_to(self, r: T) -> Self::Output;
    }

    /// Retrieve the higher-dimensional geometry that contains A and is orthogonal to B.
    ///
    /// Some examples in 3D:
    /// * Find the line containing the given point A and orthogonal to the given plane B
    /// * Find the plane containing the given point A and orthogonal to the given line B
    /// * Find the plane containing the given line A and orthogonal to the given plane B
    ///
    /// This operation is an intermediate step in computing the [Projection].
    /// [Intersecting](Meet) B with the result produces the projection of A onto B.
    pub trait SupersetOrthogonalTo<T> {
        type Output;
        fn superset_orthogonal_to(self, r: T) -> Self::Output;
    }

    /// Project a lower-dimensional element A orthogonally onto a higher-dimensional element B.
    ///
    /// Some examples in 3D:
    /// * Project a point orthogonally onto a plane
    /// * Project a point orthogonally onto a line
    /// * Project a line orthogonally onto a plane
    ///
    /// If you are looking to project higher-dimensional geometry onto lower-dimensional geometry
    /// (e.g. project a plane onto a point)
    /// use [AntiProjection].
    pub trait Projection<T> {
        type Output;
        fn projection(self, r: T) -> Self::Output;
    }

    /// Project a higher-dimensional element A orthogonally onto a lower-dimensional element B.
    ///
    /// Some examples in 3D:
    /// * Project a plane orthogonally onto a point
    /// * Project a line orthogonally onto a point
    /// * Project a plane orthogonally onto a line
    ///
    /// If you are looking to project lower-dimensional geometry onto higher-dimensional geometry
    /// (e.g. project a point onto a plane)
    /// use [Projection].
    pub trait AntiProjection<T> {
        type Output;
        fn anti_projection(self, r: T) -> Self::Output;
    }

    /// Project a lower-dimensional element A onto a higher-dimensional element B
    /// with respect to the origin
    ///
    /// Some examples in 3D:
    /// * Project a point onto a plane, along the line joining the origin and the point
    /// * Project a line onto a plane, along the plane joining the origin and the line
    pub trait CentralProjection<T> {
        type Output;
        fn central_projection(self, r: T) -> Self::Output;
    }

    /// Project a higher-dimensional element A onto a lower-dimensional element B
    /// with respect to the origin
    ///
    /// This operation does not appear to be geometrically useful
    /// but is included here for completeness
    pub trait CentralAntiProjection<T> {
        type Output;
        fn central_anti_projection(self, r: T) -> Self::Output;
    }

    /// Transform element A by motor or flector B
    ///
    /// Note that transforming by a flector will result in negated output.
    /// This is not geometrically meaningful due to homogeneous coordinates,
    /// but you may wish to add a negative sign when performing reflections.
    /// TODO test this
    pub trait Transform<T> {
        type Output;
        fn transform(self, r: T) -> Self;
    }

    /// Homogeneously scale an element so that its [bulk norm](BulkNorm) is 1.
    ///
    /// Normalization is most useful for ideal elements whose [weight norm](WeightNorm) is 0.
    /// For example, it can turn an ideal point (direction) into a vector with unit length.
    ///
    /// Normalization is generally not a useful operation to apply to regular geometry--
    /// you probably want it to be [Unitized] instead.
    pub trait Normalized {
        type Output;
        fn normalized(self) -> Self::Output;
    }

    /// Homogeneously scale an element so that its [bulk norm](WeightNorm) is 𝟙.
    ///
    /// All non-infinite elements can be unitized without changing the geometry they represent,
    /// and many functions expect their inputs to be unitized.
    ///
    /// For points, this has the effect of dividing by the projective coordinate,
    /// in effect moving the point onto the projective hyperplane.
    pub trait Unitized {
        type Output;
        fn unitized(self) -> Self::Output;
    }

    /// Compute a motor that performs twice the motion needed to take element A to element B
    ///
    /// This function is the basis of creating motors from geometry.
    /// For example, in 3D:
    /// * Given two planes, create a rotor about the line where they intersect, by twice the angle
    ///   between them
    /// * Given two lines, create a motor about their shared normal,
    ///   consisting of rotation by twice the angle between them
    ///   and translation by twice the distance between them
    /// * Given two points, create a translator along the line between them,
    ///   by twice the distance between them
    ///
    /// Like with quaternions, a motor which performs 1x the desired motion can be computed
    /// by taking the square root.
    pub trait MotorTo<T> {
        type Output;
        fn motor_to(self, r: T) -> Self::Output;
    }

    #[macro_export]
    macro_rules! impl_ops_for_scalar {
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

    impl_ops_for_scalar!(f32);
    impl_ops_for_scalar!(f64);
    impl_ops_for_scalar!(i8);
    impl_ops_for_scalar!(i16);
    impl_ops_for_scalar!(i32);
    impl_ops_for_scalar!(i64);
    impl_ops_for_scalar!(i128);
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

    #[macro_export]
    macro_rules! impl_pga2d_for_scalar {
        ($type:ident) => {
            impl From<AntiScalar<$type>> for $type {
                fn from(value: AntiScalar<$type>) -> $type {
                    value.a012
                }
            }
        };
    }

    impl_pga2d_for_scalar!(f32);
    impl_pga2d_for_scalar!(f64);
    impl_pga2d_for_scalar!(i8);
    impl_pga2d_for_scalar!(i16);
    impl_pga2d_for_scalar!(i32);
    impl_pga2d_for_scalar!(i64);
    impl_pga2d_for_scalar!(i128);

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

    #[macro_export]
    macro_rules! impl_pga3d_for_scalar {
        ($type:ident) => {
            impl From<AntiScalar<$type>> for $type {
                fn from(value: AntiScalar<$type>) -> $type {
                    value.a0123
                }
            }
        };
    }

    impl_pga3d_for_scalar!(f32);
    impl_pga3d_for_scalar!(f64);
    impl_pga3d_for_scalar!(i8);
    impl_pga3d_for_scalar!(i16);
    impl_pga3d_for_scalar!(i32);
    impl_pga3d_for_scalar!(i64);
    impl_pga3d_for_scalar!(i128);

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

#[cfg(all(test, feature = "std"))]
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
