//! 2D Vector geometry
//!
//! This module contains primitives for working with 2D vectors and antiscalars
//! in the vanilla metric [1, 1]. There is no projective dimension.
//!
//! If you want to work with geometry in 2D such as points & lines, you probably want the (re2) module instead,
//! which adds a third projective dimension.

use crate::algebraic_ops::*;
use crate::ops::*;
use crate::scalar::*;
use ngeom_macros::geometric_algebra;

geometric_algebra! {
    basis![x, y];
    metric![1, 1];

    /// A vector with two components
    #[multivector]
    #[derive(Clone, Copy, Default, Debug, PartialEq, Eq)]
    pub struct Vector<T> {
        /// The coefficient on the X basis vector
        pub x: T,
        /// The coefficient on the Y basis vector
        pub y: T,
    }

    /// A bivector / antiscalar with one component
    #[multivector]
    #[derive(Clone, Copy, Default, Debug, PartialEq, Eq, PartialOrd, Ord)]
    pub struct AntiScalar<T> {
        pub xy: T,
    }

    /// e.g. for expressing affine transformations
    ///
    /// A linear operator is fully described by where it takes the basis vectors.
    /// It can be thought of as a 2x2 matrix, where each struct field is a column.
    #[linear_operator]
    #[derive(Clone, Copy, Default, Debug, PartialEq, Eq)]
    pub struct LinearOperator<T> {
        /// The output of this operator when applied to the X basis vector
        pub x: Vector<T>,
        /// The output of this operator when applied to the Y basis vector
        pub y: Vector<T>,
    }
}

impl<T> From<T> for AntiScalar<T> {
    fn from(value: T) -> AntiScalar<T> {
        AntiScalar { xy: value }
    }
}

#[macro_export]
macro_rules! impl_v2_for_scalar {
    ($type:ident) => {
        impl From<AntiScalar<$type>> for $type {
            fn from(value: AntiScalar<$type>) -> $type {
                value.xy
            }
        }
    };
}

impl_v2_for_scalar!(f32);
impl_v2_for_scalar!(f64);
impl_v2_for_scalar!(i8);
impl_v2_for_scalar!(i16);
impl_v2_for_scalar!(i32);
impl_v2_for_scalar!(i64);
impl_v2_for_scalar!(i128);

impl<T: Ring> Origin for Vector<T> {
    fn origin() -> Vector<T> {
        Vector {
            x: T::zero(),
            y: T::zero(),
        }
    }
}

impl<T: Ring> XHat for Vector<T> {
    fn x_hat() -> Vector<T> {
        Vector {
            x: T::one(),
            y: T::zero(),
        }
    }
}

impl<T: Ring> YHat for Vector<T> {
    fn y_hat() -> Vector<T> {
        Vector {
            x: T::zero(),
            y: T::one(),
        }
    }
}
