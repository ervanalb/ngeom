//! Traits that govern the scalar data type used by ngeom
//!
//! Because geometric operations are generally done through various sums and products,
//! the scalar datatype needs only be a [Ring] for most functionality to work.

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
///
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
    + Default
{
    /// The additive identity
    fn zero() -> Self {
        Self::default()
    }

    /// The multiplicative identity
    fn one() -> Self {
        Self::from_integer(1)
    }

    /// Construct an integer scalar
    fn from_integer(i: isize) -> Self;
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
pub trait Rational: Ring {
    /// Construct a rational scalar
    /// from an integer numerator and integer denominator
    fn from_fraction(numerator: isize, denominator: isize) -> Self;

    /// A scalar value that when multiplied by 2 equals [one](Ring::one)
    fn one_half() -> Self {
        Self::from_fraction(1, 2)
    }

    /// A scalar value that when multiplied by 3 equals [one](Ring::one)
    fn one_third() -> Self {
        Self::from_fraction(1, 3)
    }

    /// A scalar value that when multiplied by 4 equals [one](Ring::one)
    fn one_fourth() -> Self {
        Self::from_fraction(1, 4)
    }

    /// A scalar value that when multiplied by 6 equals [one](Ring::one)
    fn one_sixth() -> Self {
        Self::from_fraction(1, 6)
    }
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
pub trait Sqrt: Ring {
    // The scalar datatype for the square root
    type Output;

    // This scalar's positive square root
    fn sqrt(self) -> <Self as Sqrt>::Output;
}

/// A scalar datatype which implements trigonometric functions.
///
/// Taking sines and cosines of angles is used frequently in geometry,
/// such as when constructing [rotors](crate::re3::axis_angle).
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
    type Output: Ring;

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
/// use ngeom::re2::*;
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
/// use ngeom::{impl_ops_for_scalar, impl_re2_for_scalar};
/// use ngeom::re2::*;
/// use ngeom::ops::*;
/// use ngeom::scalar::*;
///
/// /// f32 wrapper around real numbers
/// /// that does not implement `Recip`
/// // TODO write derive macro for Abs, Ring, Sqrt, etc.
/// #[derive(Clone, Copy, PartialEq, PartialOrd, Default)]
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
///     fn from_integer(i: isize) -> R32 { R32(i as f32) }
/// }
/// impl Sqrt for R32 {
///     type Output = R32;
///     fn sqrt(self) -> R32 { R32(self.0.sqrt()) }
/// }
/// impl_ops_for_scalar!(R32); // TODO replace with derive?
/// impl_re2_for_scalar!(R32);
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
pub trait Recip {
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
    type Output;
    fn anti_sqrt(self) -> Self::Output;
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
            fn from_integer(i: isize) -> $type {
                i as $type
            }
        }

        impl Rational for $type {
            fn from_fraction(numerator: isize, denominator: isize) -> $type {
                numerator as $type / denominator as $type
            }
        }

        #[cfg(feature = "std")]
        impl Sqrt for $type {
            type Output = $type;
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
            fn from_integer(i: isize) -> $type {
                i.try_into().expect("Integer out of range")
            }
        }
    };
}

impl_for_int!(i8);
impl_for_int!(i16);
impl_for_int!(i32);
impl_for_int!(i64);
impl_for_int!(i128);
