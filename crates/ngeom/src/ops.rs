//! Geometric operations

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
/// use ngeom::re2::*;
/// use ngeom::ops::*;
/// use ngeom::scalar::*;
///
/// let p1 = Vector::ideal_point([3., 4.]);
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
/// use ngeom::re2::*;
/// use ngeom::ops::*;
///
/// // Construct two unitized parallel lines
/// let l1 = Vector::point([0., 0.])
///     .join(Vector::point([10., 0.]))
///     .unitized();
/// let l2 = Vector::point([0., 10.])
///     .join(Vector::point([10., 10.]))
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
/// use ngeom::re2::*;
/// use ngeom::ops::*;
///
/// let p1 = Vector::point([2, 3]);
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
/// use ngeom::re2::*;
/// use ngeom::ops::*;
///
/// // Join two points into a line
/// let p1 = Vector::point([10., 10.]);
/// let p2 = Vector::point([13., 14.]);
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
/// use ngeom::re3::*;
/// use ngeom::ops::*;
///
/// // Join two points into a line
/// // The direction of the line goes from p1 to p2,
/// // towards +X in this case
/// let p1 = Vector::point([0., 0., 0.]);
/// let p2 = Vector::point([10., 0., 0.]);
/// let l = p1.join(p2);
///
/// // Join a line and a point into a plane
/// // The direction of the plane goes from l to p3,
/// // which produces the +XY plane in this case
/// let p3 = Vector::point([5., 5., 0.]);
/// let pl = l.join(p3);
///
/// // Join a plane and a point into a volume
/// // The direction of the volume is oriented according to the right hand rule.
/// let p4 = Vector::point([3., 2., 15.]);
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
/// use ngeom::re3::*;
/// use ngeom::ops::*;
///
/// // Join two coincident points
/// let p1 = Vector::point([4., 5., 6.]);
/// let p2 = Vector::point([4., 5., 6.]) * 2.; // Homogeneous scaling doesn't matter
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
/// use ngeom::re3::*;
/// use ngeom::ops::*;
///
/// // Join three points to get the XY plane
/// let xy = Vector::point([0., 0., 0.])
///     .join(Vector::point([1., 0., 0.]))
///     .join(Vector::point([0., 1., 0.]));
///
/// // Join two points to get a line travelling in the -Z direction
/// let l = Vector::point([3., 4., 10.])
///     .join(Vector::point([3., 4., 9.]));
///
/// // Meet the plane and the line at a point
/// let pt = xy.meet(l);
/// assert_eq!(pt.unitized(), Vector::point([3., 4., 0.]));
/// ```
///
/// `Meet` is exception-free.
/// Meeting parallel geometry will result in an ideal (infinite) element,
/// which may contain useful information such as the parallel distance.
///
/// ```
/// use ngeom::re2::*;
/// use ngeom::ops::*;
///
/// // Construct two unitized parallel lines
/// let l1 = Vector::point([0., 0.])
///     .join(Vector::point([10., 0.]))
///     .unitized();
/// let l2 = Vector::point([0., 10.])
///     .join(Vector::point([10., 10.]))
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
///
/// If you find yourself sandwiching `compose()` operations
/// to re-interpret transformations in the context of a different reference frame,
/// consider using [Transform] or [TransformInverse] instead.
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
/// Element A can itself be a motor or flector.
/// In that case, this operation reinterprets transformation A
/// from being intrinsic to B to being extrinsic.
///
/// In other words, `a.transform_inverse(b)` is equivalent to, but cheaper than,
/// `b.inverse_transformation().compose(a).compose(b)`
///
/// Transforming an element by a flector will negate that element.
/// This is geometrically equivalent due to homogeneous coordinates,
/// but you may wish to add a negative sign when performing reflections.
pub trait Transform<T> {
    type Output;
    fn transform(self, r: T) -> Self::Output;
}

/// Apply, to element A, the inverse of the transformation described by motor or flector B
///
/// If you need to invert the motor or flector before it is applied,
/// e.g. in a chain of composed transforms,
/// consider using [InverseTransformation] instead.
///
/// But also note that element A can itself be a motor or flector.
/// In that case, this operation reinterprets transformation A
/// from being extrinsic to B to being intrinsic.
///
/// In other words, `a.transform_inverse(b)` is equivalent to, but cheaper than,
/// `b.compose(a).compose(b.inverse_transformation())`
///
/// Transforming an element by a flector will negate that element.
/// This is geometrically equivalent due to homogeneous coordinates,
/// but you may wish to add a negative sign when performing reflections.
pub trait TransformInverse<T> {
    type Output;
    fn transform_inverse(self, r: T) -> Self::Output;
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
///
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

/// The motor or flector whose transformation is the inverse of the given one
///
/// For motors, the direction of the motion will be reversed.
/// For flectors, the motion will be reversed, but the flip will remain.
///
/// If you only need to invert the motor or flector when it gets applied,
/// consider using [TransformInverse] instead.
pub trait InverseTransformation {
    type Output;
    fn inverse_transformation(self) -> Self::Output;
}

/// Constructor for a motor that performs no motion
pub trait IdentityMotor {
    /// Construct a motor that performs no motion
    fn identity_motor() -> Self;
}

/// Constructor for a unitized point at the origin
pub trait Origin {
    /// Construct a unitized point at the origin
    fn origin() -> Self;
}

/// Constructor for a normalized vector (ideal point) in the X direction
pub trait XHat {
    /// Construct a [normalized](Normalized) vector (ideal point) in the X direction
    /// (first spatial dimension)
    ///
    /// This can be used to construct points in a dimension-agnostic way
    /// for spaces that are at least 1-dimensional.
    /// ```
    /// use ngeom::ops::*;
    /// use core::ops::Add;
    ///
    /// fn get_point<VECTOR: Copy + Add<VECTOR, Output=VECTOR> + Origin + XHat>() -> VECTOR {
    ///     VECTOR::origin()
    ///       + VECTOR::x_hat()
    ///       + VECTOR::x_hat()
    /// }
    /// ```
    fn x_hat() -> Self;
}

/// Constructor for a normalized vector (ideal point) in the Y direction
pub trait YHat {
    /// Construct a [normalized](Normalized) vector (ideal point) in the Y direction
    /// (second spatial dimension)
    ///
    /// This can be used to construct points in a dimension-agnostic way
    /// for spaces that are at least 2-dimensional.
    /// ```
    /// use ngeom::ops::*;
    /// use core::ops::Add;
    ///
    /// fn get_point<VECTOR: Copy + Add<VECTOR, Output=VECTOR> + Origin + XHat + YHat>() -> VECTOR {
    ///     VECTOR::origin()
    ///       + VECTOR::x_hat()
    ///       + VECTOR::x_hat()
    ///       + VECTOR::x_hat()
    ///       + VECTOR::y_hat()
    ///       + VECTOR::y_hat()
    /// }
    /// ```
    fn y_hat() -> Self;
}

/// Constructor for a normalized vector (ideal point) in the Z direction
pub trait ZHat {
    /// Construct a [normalized](Normalized) vector (ideal point) in the Z direction
    /// (third spatial dimension)
    ///
    /// This can be used to construct points in a dimension-agnostic way
    /// for spaces that are at least 3-dimensional.
    /// ```
    /// use ngeom::ops::*;
    /// use core::ops::Add;
    ///
    /// fn get_point<VECTOR: Copy + Add<VECTOR, Output=VECTOR> + Origin + XHat + YHat + ZHat>() -> VECTOR {
    ///     VECTOR::origin()
    ///       + VECTOR::x_hat()
    ///       + VECTOR::x_hat()
    ///       + VECTOR::x_hat()
    ///       + VECTOR::y_hat()
    ///       + VECTOR::y_hat()
    ///       + VECTOR::z_hat()
    /// }
    /// ```
    fn z_hat() -> Self;
}

/// Constructor for a unitized vector corresponding to a geometric point
/// at the given coordinates
pub trait Point<C> {
    /// Construct a unitized vector corresponding to a geometric point
    /// at the given coordinates
    fn point(x: C) -> Self;
}

/// Constructor for a vector corresponding to an ideal (infinite) point
/// at the given coordinates
pub trait IdealPoint<C> {
    /// Construct a vector corresponding to an ideal (infinite) point
    /// at the given coordinates
    fn ideal_point(x: C) -> Self;
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

