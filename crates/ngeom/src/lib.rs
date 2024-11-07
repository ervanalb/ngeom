#![cfg_attr(not(feature = "std"), no_std)]

//! ngeom is a library for doing geometry in N dimensions.
//!
//! ngeom ships with modules for rigid Euclidean [2D](re2) and [3D](re3) geometry.
//! The built-in modules provide primitives for points, lines, and planes,
//! along with rigid transformations (rotations, translations, and reflections.)
//!
//! If you prefer to own your geometric primitives,
//! or if you desire a space with different dimensionality, metric, handedness, etc.
//! then you can use the [provided macro](geometric_algebra)
//! to generate all the necessary [geometric operations](ops) on your own `struct`s.
//!
//! ngeom uses [homogeneous coordinates](https://en.wikipedia.org/wiki/Homogeneous_coordinates) to express ideal/infinite points,
//! ideal/infinite lines, etc. and to provide for exception-free [meet](ops::Meet) & [join](ops::Join) operations.
//!
//! ngeom is generic over the [scalar] datatype, and can be used with `f32`, `f64`, or custom datatypes.
//! It can even use integers, with some restrictions--
//! most functionality requires only [scalar addition and multiplication](scalar::Ring).
//!
//! With some effort, usage code can be written to be generic
//! over the dimensionality or even the metric of the space.
//!
//! ngeom does not use SIMD intrinsics,
//! but rather emits compact expressions and structs
//! that are easy for the compiler to optimize into SIMD instructions
//! in an architecture-independent way.
//!
//! ngeom is dependency-free.
//! ngeom is `no_std`-compatible, with most functionality available,
//! and the option to implement the rest.
//!
//! Under the hood, ngeom is able to make these generalizations by using [geometric algebra](https://en.wikipedia.org/wiki/Geometric_algebra).
//! Understanding geometric algebra is helpful but optional,
//! as the functions are named for their geometric interpretation.
//! Lower-level functions named after their [geometric algebra expressions](algebraic_ops) are available if needed.

/// Generate geometric operations on the given set of structs
///
/// This low-level macro can be used to generate different geometric algebras
/// if the pre-made modules for [2D](re2) or [3D](re3) rigid Euclidean geometry
/// are not sufficient.
///
/// Here are some reasons you may want to use this macro:
/// * You want ownership of the geometry structs,
///   e.g. so you can implement traits on them
///   within the constraints of the [orphan rule](https://doc.rust-lang.org/book/ch10-02-traits.html)
/// * You want to rename the structs or their fields
/// * You want a left-handed space
/// * You want multivector structs with different subsets of the basis elements
/// * You want a different metric, to implement hyperbolic space or spherical space
/// * You want a different number of dimensions, e.g. 1D or 5D
///
/// Careful! This macro is a bit finnicky with regards to naming and ordering.
/// See the notes in the following example.
///
/// ```
/// use ngeom::geometric_algebra;
///
/// // These traits must be brought into scope before invoking the macro:
/// use ngeom::scalar::*;
/// use ngeom::ops::*;
/// use ngeom::algebraic_ops::*;
///
/// // Define a geometric algebra with the given basis vectors, metric,
/// // and structs representing multivectors.
/// geometric_algebra! {
///     // Use the basis![] pseudo-macro to call out the names of the basis vectors.
///     // This also sets the sign of the anti-scalar,
///     // which corresponds to the handedness of the space.
///     basis![w, x, y];
///
///     // Basis vectors may be single letters,
///     // or may be a common prefix plus a single character e.g. [e_0, e_1, e_2]
///     // Higher-grade basis elements retain one copy of the prefix,
///     // e.g. e_21, e_012
///
///     // Use the metric![] pseudo-macro to call out the metric of the space--
///     // this is a list of numbers parallel to the basis,
///     // where each number is the dot product of that basis vector with itself.
///     metric![0, 1, 1];
///
///     #[multivector] // Use this pseudo-macro to call out that a struct is a multivector
///     #[derive(Clone, Copy)] // Multivectors must be Copy for use in math expressions
///     struct Vector<T> { // Multivectors must have a single generic parameter named T
///         x: T, // each field must be type T (the scalar type)
///         y: T,
///         w: T,
///         // Any basis components that are not present are assumed to be zero.
///     }
///
///     #[multivector]
///     #[derive(Clone, Copy)]
///     struct Bivector<T> {
///         wx: T, // the coefficient on the basis bivector w ∧ x
///         wy: T, // the coefficient on the basis bivector w ∧ y
///         xy: T, // the coefficient on the basis bivector x ∧ y
///
///         // The spelling of the field name sets the multiplication order,
///         // which may introduce a negative sign.
///         // (This won't fundamentally change the behavior of the algebra,
///         // which is defined fully by the metric and sign of the anti-scalar.)
///     }
///
///     #[multivector]
///     #[derive(Clone, Copy)]
///     struct AntiScalar<T> {
///         wxy: T, // the coefficient on the basis trivector w ∧ x ∧ y
///
///         // The spelling of this field should match the order set by basis![].
///         // Otherwise, anti-scalar operations like anti_mul() and anti_sqrt()
///         // won't get automatically implemented.
///     }
///
///     // Compound object for elements with even antigrade
///     // (which can represent any motor)
///     #[multivector]
///     #[derive(Clone, Copy)]
///     struct AntiEven<T> {
///         w: T,
///         x: T,
///         y: T,
///         wxy: T,
///     }
///
///     // Compound object for elements with odd antigrade
///     // (which can represent any flector)
///     #[multivector]
///     #[derive(Clone, Copy)]
///     struct AntiOdd<T> {
///         a: T, // A field whose name isn't composed of basis vectors
///               // is assumed to be the scalar component
///         wx: T,
///         wy: T,
///         xy: T,
///     }
///
///     // The order that the structs are defined is important--
///     // the return type of a geometric operation will be the
///     // *first struct* that can represent the output.
///     // For example, were `AntiEven` defined above `Vector`,
///     // all geometric operations that return a vector
///     // would instead use `AntiEven` for their output type.
///     // For this reason, more comprehensive compound structs
///     // should be defined last.
///
///
///     // You may optionally define a struct that will serve as a linear operator.
///     // Requirements here are strict:
///     // this struct must contain a field for every basis vector,
///     // and each field must itself be a vector.
///     #[linear_operator]
///     #[derive(Clone, Copy, Default, Debug, PartialEq, Eq)]
///     pub struct LinearOperator<T> {
///         pub x: Vector<T>,
///         pub y: Vector<T>,
///         pub w: Vector<T>,
///     }
/// }
///
/// // The geometric_algebra! macro will emit the contained code verbatim
/// // (with the pseudo-macros stripped out.)
/// // It will then append `impl` blocks for all applicable traits in
/// // `ngeom::ops` and `ngeom::algebraic_ops`
/// // across all of the given multivector structs.
pub use ngeom_macros::geometric_algebra;

pub mod scalar;
pub mod algebraic_ops;
pub mod ops;
pub mod re2;
pub mod re3;

mod test;
