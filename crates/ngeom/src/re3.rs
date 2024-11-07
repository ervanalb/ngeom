//! Rigid Euclidean 3D geometry
//!
//! This module contains primitives for doing rigid geometry in 3D Euclidean space.
//!
//! It is *rigid* because its primitives can express transformations that preserve length--
//! such as translations, rotations, and reflections.
//! It cannot inherently represent non-rigid transformations, such as dilation or shear--
//! you will have to eject into a matrix representation (or similar) if you want to do those.
//!
//! It is *3D* because there are three spatial dimensions.
//! A fourth projective dimension is added, resulting in homogeneous coordinates.
//!
//! It is *Euclidean* because the space has zero curvature,
//! i.e. the parallel postulate holds.

use crate::algebraic_ops::*;
use crate::ops::*;
use crate::scalar::*;
use ngeom_macros::geometric_algebra;

geometric_algebra! {
    basis![w, x, y, z];
    metric![0, 1, 1, 1];

    /// e.g. a point in space or an ideal (infinite) point
    ///
    /// ## As geometry
    /// Geometrically, a `Vector` can be either:
    /// * A point in space, e.g. `Vector {x, y, z, w: 1}`, which represents a location.
    ///   See the [`Vector::point([x, y, z])`](crate::ops::Point)
    ///   or [`Vector::origin()`](crate::ops::Origin) constructors.
    /// * An ideal point, or point at infinity, e.g. `Vector {x, y, z, w: 0}`, which representing a direction.
    ///   See the [`Vector::ideal_point([x, y, z])`](crate::ops::IdealPoint),
    ///   [`Vector::x_hat()`](crate::ops::XHat),
    ///   [`Vector::y_hat()`](crate::ops::YHat), or
    ///   [`Vector::z_hat()`](crate::ops::ZHat) constructors.
    ///
    /// All `Vector`s represent some kind of point.
    ///
    /// ## As a transformation
    /// When interpreted as a transformation, a point performs an inversion--
    /// equivalent to three reflections through the mutually orthogonal planes
    /// that meet at that point.
    ///
    /// ## Example Operations
    /// * Two points [join](crate::ops::Join) into a [line](Bivector) containing both.
    /// * Two unitized points [join](crate::ops::Join) into a line whose [weight norm](crate::ops::WeightNorm) is the distance between them.
    /// * Three unitized points [join](crate::ops::Join) into a [plane](Trivector).
    /// * Four unitized points [join](crate::ops::Join) into a parallelepiped's [signed volume](AntiScalar).
    /// * The sum of two unitized points is their midpoint.
    /// * [Unitizing](crate::ops::Unitized) a point divides by the projective coordinate
    ///   to make w=1 without changing its location
    /// * [Normalizing](crate::ops::Normalized) an ideal point divides by its length
    ///   to make a unit vector.
    /// * The [motor between](crate::ops::MotorTo) two unitized points produces a translation
    ///   by twice their separation
    ///
    #[multivector]
    #[derive(Clone, Copy, Default, Debug, PartialEq, Eq)]
    pub struct Vector<T> {
        pub w: T,
        pub x: T,
        pub y: T,
        pub z: T,
    }

    /// e.g. a line in space or an ideal (infinite) line
    ///
    /// ## As geometry
    /// Geometrically, a `Bivector` can be either:
    /// * A line in space
    /// * An ideal line, which can be thought to encircle the space infinitely far away,
    ///   like the horizon.
    /// Only simple bivectors (bivectors that can be built from the wedge product of vectors)
    /// represent lines. See this page on the
    /// [geometric constraint](https://rigidgeometricalgebra.org/wiki/index.php?title=Geometric_constraint).
    ///
    /// ## As a transformation
    /// When interpreted as a transformation,
    /// a line performs a 180-degree rotation about itself.
    ///
    /// ## Example Operations
    /// * A line and a [point](Vector) [join](crate::ops::Join) into a [plane](Trivector)
    /// * A line and a [plane](Trivector) [meet](crate::ops::Meet) at a [point](Vector).
    /// * Two unitized lines [join](crate::ops::Join) into the [signed volume](AntiScalar)
    ///   of the parallelepiped built from unit segments on each
    /// * The [motor between](crate::ops::MotorTo) two unitized lines produces a rotation about their
    ///   shared normal line by twice the angle between them
    ///
    #[multivector]
    #[derive(Clone, Copy, Default, Debug, PartialEq, Eq)]
    pub struct Bivector<T> {
        pub wx: T,
        pub wy: T,
        pub wz: T,
        pub xy: T,
        pub yz: T,
        pub zx: T,
    }

    /// e.g. a plane in space or an ideal (infinite) plane
    ///
    /// ## As geometry
    /// Geometrically, a `Trivector` can be either:
    /// * A plane in space
    /// * An ideal plane, which can be thought to enclose the space infinitely far away,
    ///   like the sky.
    /// All `Trivector`s represent some kind of plane.
    ///
    /// ## As a transformation
    /// When interpreted as a transformation,
    /// a plane performs a reflection through itself.
    ///
    /// ## Example Operations
    /// * Two planes lines [meet](crate::ops::Meet) at a [line](Bivector)
    /// * A plane and a [line](Bivector) [meet](crate::ops::Meet) at a [point](Vector).
    /// * The [motor between](crate::ops::MotorTo) two unitized planes produces a rotation about their
    ///   intersection line by twice the angle between them
    ///
    #[multivector]
    #[derive(Clone, Copy, Default, Debug, PartialEq, Eq)]
    pub struct Trivector<T> {
        pub wxy: T,
        pub wzx: T,
        pub wyz: T,
        pub xyz: T,
    }

    /// e.g. signed volume
    ///
    /// An antiscalar or [Pseudoscalar](https://en.wikipedia.org/wiki/Pseudoscalar)
    /// is a quantity that behaves like a scalar
    /// except in the case of improper isometry such as reflection,
    /// where it gains a sign flip.
    ///
    /// Antiscalars are a distinct type from scalars, but behave similarly.
    /// In this documentation, they are typically written in blackboard bold,
    /// e.g. 𝟙, 𝟚, 𝟛.
    ///
    /// Operations on antiscalars that are dual to those on scalars
    /// are prefixed with `anti_` e.g. [anti_sqrt()](crate::scalar::AntiSqrt),
    /// [anti_recip()](crate::scalar::AntiRecip)
    ///
    /// ## As geometry
    /// * Geometrically, an Antiscalar can be thought to represent a signed volume.
    /// All `AntiScalar`s represent some amount of signed volume.
    ///
    /// ## As a transformation
    /// * When interpreted as a transformation, the antiscalar 𝟙 is the identity transformation.
    ///
    /// ## Example Operations
    /// * Antiscalars typically result from taking the [weight norm](crate::ops::WeightNorm).
    /// * Antiscalars may be cast to scalars using `.into()` (and vice-versa)
    #[multivector]
    #[derive(Clone, Copy, Default, Debug, PartialEq, Eq, PartialOrd, Ord)]
    pub struct AntiScalar<T> {
        pub wxyz: T,
    }

    /// e.g. a homogeneous magnitude
    ///
    /// `DualNumber` holds the sum of a scalar and [antiscalar](AntiScalar).
    /// See [Wikipedia](https://en.wikipedia.org/wiki/Dual_number)
    /// for more details on the mathematical properties of dual numbers.
    ///
    /// Dual numbers appear in several places:
    /// * The [anti-commutator product](crate::algebraic_ops::AntiCommutator) between two unitized skew [lines](Bivector)
    ///   results in a non-simple bivector,
    ///   but this bivector can be factored into the ⟇ product of a unitized line and a dual number.
    ///   That line is the common normal,
    ///   the dual number's antiscalar part is the distance between them,
    ///   and dual number's scalar part is the angle between them.
    /// * A [motor](AntiEven) computed by numerical methods
    ///   may diverge from the motor manifold.
    ///   This means that the motor composed with its inverse transformation (A ⟇ A̰)
    ///   is not the identity motor (𝟙)--
    ///   Instead it will be a dual number.
    ///   It can be brought back onto the motor manifold
    ///   by taking the ⟇ product of it and the reciprocal square root (under ⟇) of this dual number.
    #[multivector]
    #[derive(Clone, Copy, Default, Debug, PartialEq, Eq)]
    pub struct DualNumber<T> {
        pub a: T,
        pub wxyz: T,
    }

    /// e.g. a motor
    ///
    /// `AntiEven` holds the sum of a scalar, [bivector](Bivector) and [antiscalar](AntiScalar).
    /// It is so named because these are the elements whose antigrade is even
    /// (Ag(scalar) = 4, Ag(bivector) = 2, Ag(antiscalar) = 0)
    /// An alternative interpretation is that it stores transformations
    /// that are composed of an even number of reflections.
    ///
    /// ## As geometry
    /// `AntiEven` should probably not be used to represent geometric objects.
    /// It is not geometrically meaningful if it is mixed-grade.
    ///
    /// ## As a transformation
    /// This struct can hold any motor, i.e. it can describe any proper isometry
    /// (any combination of rotation & translation.)
    /// It cannot represent an odd number of reflections, e.g. a single reflection,
    /// which would result in an improper isometry.
    /// For that, see [AntiOdd].
    ///
    /// Generally, a motor [composed](crate::ops::Compose) with its [inverse transformation](crate::ops::InverseTransformation)
    /// should equal the identity motor 𝟙.
    /// If it does not, it can be multiplied by the reciprocal square root (under ⟇) of this product
    /// to return it to the motor manifold
    /// (similar to normalizing a quaternion.)
    ///
    /// ## Example Operations
    /// * [Composing](crate::ops::Compose) motors A and B results in a third motor whose motion is equivalent to A followed by B
    /// * Composing a motor A and a flector B results in a [flector](AntiOdd) whose transformation is equivalent to A followed by B
    /// * A motor can be made to move the opposite direction using [InverseTransformation](crate::ops::InverseTransformation)
    /// * A motor's reference frame can be changed using [Transform](crate::ops::Transform) or [TransformInverse](crate::ops::TransformInverse)
    /// * The square root of a motor (under ⟇) moves half of the original amount
    #[multivector]
    #[derive(Clone, Copy, Default, Debug, PartialEq, Eq)]
    pub struct AntiEven<T> {
        pub a: T,
        pub wx: T,
        pub wy: T,
        pub wz: T,
        pub xy: T,
        pub yz: T,
        pub zx: T,
        pub wxyz: T,
    }

    /// e.g. a flector
    ///
    /// `AntiOdd` holds the sum of a [vector](Vector) and [trivector](Trivector).
    /// It is so named because these are the elements whose antigrade is odd
    /// (Ag(vector) = 3, Ag(trivector) = 1.)
    /// An alternative interpretation is that it stores transformations
    /// that are composed of an odd number of reflections.
    ///
    /// ## As geometry
    /// `AntiOdd` should probably not be used to represent geometric objects.
    /// It is not geometrically meaningful if it contains
    /// both a vector part and a trivector part.
    ///
    /// ## As a transformation
    /// This struct can hold any flector, i.e. it can describe any improper isometry
    /// (a reflection followed by any combination of rotation & translation)
    /// It cannot represent an even number of reflections, e.g. no reflections,
    /// which would result in a proper isometry.
    /// For that, see [AntiEven].
    ///
    /// Generally, a flector [composed](crate::ops::Compose) with its [inverse transformation](crate::ops::InverseTransformation)
    /// should equal the identity motor 𝟙.
    /// If it does not, it can be multiplied by the reciprocal square root (under ⟇) of this product
    /// to return it to the flector manifold
    /// (similar to normalizing a quaternion.)
    ///
    /// ## Example Operations
    /// * [Composing](crate::ops::Compose) flectors A and B results in a *[motor](AntiEven)* whose motion is equivalent to A followed by B
    /// * Composing a flector A and a motor B results in a flector whose transformation is equivalent to A followed by B
    /// * A flector can be made to move the opposite direction (while preserving the flip) using [InverseTransformation](crate::ops::InverseTransformation)
    /// * A flector's reference frame can be changed using [Transform](crate::ops::Transform) or [TransformInverse](crate::ops::TransformInverse)
    #[multivector]
    #[derive(Clone, Copy, Default, Debug, PartialEq, Eq)]
    pub struct AntiOdd<T> {
        pub x: T,
        pub y: T,
        pub z: T,
        pub w: T,
        pub xyz: T,
        pub wxy: T,
        pub wzx: T,
        pub wyz: T,
    }

    /// e.g. for expressing affine transformations
    ///
    /// This may be used to transform geometry in non-rigid ways,
    /// such as non-uniform scale & shear.
    ///
    /// It can also represent rigid operations like rotation, translation, and reflection,
    /// but using a [motor](AntiEven) or [flector](AntiOdd) for those
    /// may give better ergonomics and results.
    ///
    /// A linear operator is fully described by where it takes the basis vectors.
    /// It can be thought of as a matrix, where each struct field is a column.
    #[linear_operator]
    #[derive(Clone, Copy, Default, Debug, PartialEq, Eq)]
    pub struct LinearOperator<T> {
        /// The output of this operator when applied to the X basis vector
        pub x: Vector<T>,
        /// The output of this operator when applied to the Y basis vector
        pub y: Vector<T>,
        /// The output of this operator when applied to the Z basis vector
        pub z: Vector<T>,
        /// The output of this operator when applied to the W basis vector
        pub w: Vector<T>,
    }
}

impl<T: Ring> From<T> for AntiScalar<T> {
    fn from(value: T) -> AntiScalar<T> {
        AntiScalar { wxyz: value }
    }
}

#[macro_export]
macro_rules! impl_re3_for_scalar {
    ($type:ident) => {
        impl From<AntiScalar<$type>> for $type {
            fn from(value: AntiScalar<$type>) -> $type {
                value.wxyz
            }
        }
    };
}

impl_re3_for_scalar!(f32);
impl_re3_for_scalar!(f64);
impl_re3_for_scalar!(i8);
impl_re3_for_scalar!(i16);
impl_re3_for_scalar!(i32);
impl_re3_for_scalar!(i64);
impl_re3_for_scalar!(i128);

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

impl<T: Ring + Sqrt<Output = T> + Recip<Output = T>> DualNumber<T> {
    pub fn rsqrt(self) -> DualNumber<T> {
        let DualNumber { a: s, wxyz: p } = self;
        let sqrt_s = s.sqrt();
        let sqrt_s_cubed = s * sqrt_s;
        DualNumber {
            a: sqrt_s.recip(),
            wxyz: p * (sqrt_s_cubed + sqrt_s_cubed).recip(),
        }
    }
}

impl<T: Ring> Origin for Vector<T> {
    fn origin() -> Vector<T> {
        Vector {
            x: T::zero(),
            y: T::zero(),
            z: T::zero(),
            w: T::one(),
        }
    }
}

impl<T: Ring> XHat for Vector<T> {
    fn x_hat() -> Vector<T> {
        Vector {
            x: T::one(),
            y: T::zero(),
            z: T::zero(),
            w: T::zero(),
        }
    }
}

impl<T: Ring> YHat for Vector<T> {
    fn y_hat() -> Vector<T> {
        Vector {
            x: T::zero(),
            y: T::one(),
            z: T::zero(),
            w: T::zero(),
        }
    }
}

impl<T: Ring> ZHat for Vector<T> {
    fn z_hat() -> Vector<T> {
        Vector {
            x: T::zero(),
            y: T::zero(),
            z: T::one(),
            w: T::zero(),
        }
    }
}

impl<T: Ring> IdentityTransformation for AntiScalar<T> {
    fn identity_transformation() -> AntiScalar<T> {
        AntiScalar { wxyz: T::one() }
    }
}

impl<T: Ring> IdentityTransformation for AntiEven<T> {
    fn identity_transformation() -> AntiEven<T> {
        AntiScalar { wxyz: T::one() }.into()
    }
}

impl<T: Ring> Point<[T; 3]> for Vector<T> {
    fn point([x, y, z]: [T; 3]) -> Vector<T> {
        Vector {
            x,
            y,
            z,
            w: T::one(),
        }
    }
}

impl<T: Ring> IdealPoint<[T; 3]> for Vector<T> {
    fn ideal_point([x, y, z]: [T; 3]) -> Vector<T> {
        Vector {
            x,
            y,
            z,
            w: T::zero(),
        }
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
pub fn rotor<T: Ring + Rational + Trig<Output = T> + Sqrt<Output = T>>(
    l: Bivector<T>,
) -> AntiEven<T> {
    let half_phi = l.weight_norm() * T::one_half();
    l.anti_mul(half_phi.anti_sinc()) + half_phi.anti_cos()
}

/// Translator towards ideal point p by its magnitude
pub fn translator<T: Ring + Rational>(p: Vector<T>) -> AntiEven<T> {
    Vector::<T>::origin().wedge(p).weight_dual() * T::one_half() + AntiScalar::from(T::one())
}
