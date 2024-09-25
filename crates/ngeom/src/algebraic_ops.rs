//! Low-level geometric algebra operations
//!
//! There are a few conflicting conventions for geometric algebra notation.
//! This library tends to use the formulation put forth by
//! [Dr. Eric Lengyel](https://projectivegeometricalgebra.org/)
//!
//! Consider using the aliases in the [ops] module when available,
//! for code that reflects the geometric interpretation.
//! (e.g. when intersecting two planes, prefer `l1.meet(l2)` over `l1.anti_wedge(l2)`)

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

