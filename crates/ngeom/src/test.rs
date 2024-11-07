#![cfg(all(test, feature = "std"))]

use crate::ops::*;
use crate::{re2, re3};

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

impl IsClose for re2::AntiScalar<f32> {
    fn is_close(self, rhs: re2::AntiScalar<f32>) -> bool {
        <f32>::from(self - rhs).abs() < 1e-5
    }
}

impl IsClose for re3::AntiScalar<f32> {
    fn is_close(self, rhs: re3::AntiScalar<f32>) -> bool {
        <f32>::from(self - rhs).abs() < 1e-5
    }
}

impl IsClose for re2::Vector<f32> {
    fn is_close(self, rhs: Self) -> bool {
        let tol = 1e-5;
        let diff = self - rhs;
        // Taking the right complement allows us to get the weight norm as a scalar
        // rather than an antiscalar
        diff.bulk_norm_squared() < tol * tol && diff.weight_norm_squared() < (tol * tol).into()
    }
}

impl IsClose for re2::Bivector<f32> {
    fn is_close(self, rhs: Self) -> bool {
        let tol = 1e-5;
        let diff = self - rhs;
        diff.bulk_norm_squared() < tol * tol && diff.weight_norm_squared() < (tol * tol).into()
    }
}

impl IsClose for re3::Vector<f32> {
    fn is_close(self, rhs: Self) -> bool {
        let tol = 1e-5;
        let diff = self - rhs;
        diff.bulk_norm_squared() < tol * tol && diff.weight_norm_squared() < (tol * tol).into()
    }
}

impl IsClose for re3::Bivector<f32> {
    fn is_close(self, rhs: Self) -> bool {
        let tol = 1e-5;
        let diff = self - rhs;
        diff.bulk_norm_squared() < tol * tol && diff.weight_norm_squared() < (tol * tol).into()
    }
}

impl IsClose for re3::Trivector<f32> {
    fn is_close(self, rhs: Self) -> bool {
        let tol = 1e-5;
        let diff = self - rhs;
        diff.bulk_norm_squared() < tol * tol && diff.weight_norm_squared() < (tol * tol).into()
    }
}

#[test]
fn test_dist_2d() {
    // Distance between points
    // Joining the points results in a line whose norm is their unsigned distance
    {
        let p1 = re2::Vector::<f32>::point([10., 10.]);
        let p2 = re2::Vector::point([13., 14.]);

        assert_close!(p1.join(p2).weight_norm(), (5.).into());
    }

    // Signed distance between point & line
    // Joining the line & point results in a plane (scalar) corresponding to their signed distance
    {
        // Note that the line and the test point form a triangle that is wound
        // left-handed, producing a negative distance
        let l = re2::Vector::<f32>::point([10., 10.])
            .join(re2::Vector::point([10., 20.]))
            .unitized();
        let p = re2::Vector::point([15., 15.]);
        assert_close!(l.join(p).weight_norm(), (-5.).into());
    }

    // Distance between parallel lines
    // The lines meet at an infinite point whose infinite norm is their unsigned distance.
    // More precisely, this expression yields the distance between the projection of the
    // origin onto each line.
    {
        let l1 = re2::Vector::<f32>::origin()
            .join(re2::Vector::point([3., 4.]))
            .unitized();
        let l2 = re2::Vector::point([4., -3.])
            .join(re2::Vector::point([7., 1.]))
            .unitized();
        assert_close!(dbg!(l1.meet(l2)).bulk_norm(), 5.);
    }
}

#[test]
fn test_dist_3d() {
    // Distance between points
    // Joining the points results in a line whose norm is their unsigned distance
    {
        let p1 = re3::Vector::<f32>::point([10., 10., 10.]);
        let p2 = re3::Vector::point([13., 14., 10.]);
        assert_close!(p1.join(p2).weight_norm(), (5.).into());
    }

    // Distnce between point & line
    // Joining the line & point results in a plane whose norm is their unsigned distance
    {
        let l = re3::Vector::<f32>::point([10., 10., 10.])
            .join(re3::Vector::point([10., 20., 10.]))
            .unitized();
        let p = re3::Vector::point([15., 15., 10.]);
        assert_close!(l.join(p).weight_norm(), (5.).into());
    }

    {
        let l = re3::Vector::<f32>::point([10., 10., 0.])
            .join(re3::Vector::point([10., 10., 20.]))
            .unitized();
        let p = re3::Vector::point([13., 14., 10.]);
        assert_close!(l.join(p).weight_norm(), (5.).into());
    }

    // Distance between point & plane
    // Joining the plane & point results in a scalar corresponding to their signed distance
    {
        let pl = re3::Vector::<f32>::point([10., 10., 10.])
            .join(re3::Vector::point([20., 10., 10.]))
            .join(re3::Vector::point([10., 20., 10.]))
            .unitized();
        let pt = re3::Vector::point([15., 15., 15.]);
        assert_close!(pl.join(pt), (5.).into());
    }

    // Distance between parallel lines
    // More precisely, this expression yields the distance between the projection of the
    // origin onto each line.
    {
        use super::algebraic_ops::AntiCommutator;
        let l1 = re3::Vector::<f32>::origin()
            .join(re3::Vector::point([0., 0., 10.]))
            .unitized();
        let l2 = re3::Vector::<f32>::point([3., 4., 10.])
            .join(re3::Vector::point([3., 4., 20.]))
            .unitized();
        assert_close!(l1.anti_commutator(l2).bulk_norm(), 5.);
    }

    // Distance between perpendicular lines
    // More precisely, this expression yields d * sin(a)
    // d is the distance between the lines, and a is the angle between them
    // (as measured along / about their shared normal)
    {
        let _l1 = re3::Vector::<f32>::origin()
            .join(re3::Vector::point([0., 0., 10.]))
            .unitized();
        let _l2 = re3::Vector::<f32>::point([8., 0., 15.])
            .join(re3::Vector::point([8., 20., 15.]))
            .unitized();
        //assert_close!(l1.join(l2), re3::anti(-8.)); // XXX broken
    }

    // Distance between skew lines
    // This expression divides the join of the two lines (d * sin(a))
    // by the scalar norm of their commutator product, which is sin(a)
    {
        let _l1 = re3::Vector::<f32>::origin()
            .join(re3::Vector::point([0., 0., 10.]))
            .unitized();
        let _l2 = re3::Vector::point([10., 0., 0.])
            .join(re3::Vector::point([10., 5., 4.]))
            .unitized();

        //assert_close!(l1.join(l2) + l1.anti_commutator(l2).weight_norm(), -10.) // XXX broken
    }
}

#[test]
fn test_sign_2d_triangle_winding() {
    let p1 = re2::Vector::<f32>::origin();
    let p2 = re2::Vector::point([1., 0.]);
    let p3 = re2::Vector::point([0., 1.]);
    assert!(p1.join(p2).join(p3) > (0.).into());
    assert!(p3.join(p2).join(p1) < (0.).into());
}

#[test]
fn test_sign_2d_line_intersection() {
    let l1 = re2::Vector::<f32>::origin()
        .join(re2::Vector::point([10., 0.]))
        .unitized();
    let l2 = re2::Vector::point([0., 10.])
        .join(re2::Vector::point([10., 9.]))
        .unitized();
    let l3 = re2::Vector::point([0., -10.])
        .join(re2::Vector::point([-10., -9.]))
        .unitized();

    assert!(l1.meet(l2).weight_norm() < (0.).into());
    assert!(l2.meet(l1).weight_norm() > (0.).into());
    assert!(l1.meet(l3).weight_norm() > (0.).into());
    assert!(l3.meet(l1).weight_norm() < (0.).into());
}

#[test]
fn test_sign_3d_tetrahedron_winding() {
    let p1 = re3::Vector::<f32>::origin();
    let p2 = re3::Vector::point([1., 0., 0.]);
    let p3 = re3::Vector::point([0., 1., 0.]);
    let p4 = re3::Vector::point([0., 0., 1.]);
    assert!(p1.join(p2).join(p3).join(p4) > (0.).into());
    assert!(p3.join(p2).join(p1).join(p4) < (0.).into());
}

#[test]
fn test_sign_3d_skew_lines() {
    // Note that this is intentionally backwards!
    // It's impossible to reconcile the sign of this
    // with the sign of the plane-point version,
    // so we add the negative sign here
    let l1 = re3::Vector::<f32>::origin()
        .join(re3::Vector::point([0., 0., 10.]))
        .unitized();
    let l2 = re3::Vector::<f32>::point([8., 0., 15.])
        .join(re3::Vector::point([8., 20., 15.]))
        .unitized();
    let l3 = re3::Vector::<f32>::point([-10., 0., 0.])
        .join(re3::Vector::point([-10., 20., -5.]))
        .unitized();

    assert!(l1.join(l2) > (0.).into());
    assert!(l1.join(l3) < (0.).into());
}

#[test]
fn test_sign_3d_plane_intersection() {
    let p1 = re3::Vector::<f32>::origin()
        .join(re3::Vector::point([1., 0., 0.]))
        .join(re3::Vector::point([0., 1., 0.]))
        .unitized();
    let p2 = re3::Vector::<f32>::origin()
        .join(re3::Vector::point([0., 1., 0.]))
        .join(re3::Vector::point([0., 0., 1.]))
        .unitized();
    let p3 = re3::Vector::<f32>::origin()
        .join(re3::Vector::point([0., 0., 1.]))
        .join(re3::Vector::point([1., 0., 0.]))
        .unitized();

    dbg!(p1.meet(p2).meet(p3).weight_norm() > (0.).into());
}

#[test]
fn test_2d_line_normal() {
    let l = re2::Vector::<f32>::origin()
        .join(re2::Vector::point([1., 0.]))
        .unitized();
    let expected_dir = re2::Vector::ideal_point([0., 1.]);

    let d = l.normal();

    assert!((d - expected_dir).weight_norm_squared() < (1e-5 * 1e-5).into());
}

#[test]
fn test_plane_normal() {
    let p = re3::Vector::<f32>::origin()
        .join(re3::Vector::point([1., 0., 0.]))
        .join(re3::Vector::point([0., 1., 0.]))
        .unitized();

    let expected_dir = re3::Vector::ideal_point([0., 0., 1.]);

    let d = p.normal();

    assert!((d - expected_dir).weight_norm_squared() < (1e-5 * 1e-5).into());
}

#[test]
fn test_3d_line_normal() {
    let l = re3::Vector::<f32>::origin()
        .join(re3::Vector::point([0., 1., 0.]))
        .unitized();

    let inf_line = l.normal();

    let expected_inf_line =
        re3::Vector::ideal_point([0., 0., 1.]).join(re3::Vector::ideal_point([1., 0., 0.]));

    assert!((inf_line - expected_inf_line).weight_norm_squared() < (1e-5 * 1e-5).into());
}

#[test]
fn test_2d_rotation() {
    let p = re2::Vector::point([1., 0.]);
    let motor = re2::axis_angle(re2::Vector::origin(), 0.25 * core::f32::consts::TAU);
    assert_close!(p.transform(motor), re2::Vector::point([0., 1.]));
}

#[test]
fn test_3d_rotation() {
    let p = re3::Vector::point([1., 0., 0.]);
    let l = re3::Vector::<f32>::origin()
        .join(re3::Vector::point([0., 0., 1.]))
        .unitized();
    let motor = re3::axis_angle(l, 0.25 * core::f32::consts::TAU);
    assert_close!(p.transform(motor), re3::Vector::point([0., 1., 0.]));
}

#[test]
fn test_2d_translation() {
    let p1 = re2::Vector::point([10., 10.]);

    let motor = re2::translator(re2::Vector::ideal_point([0., 5.]));

    let p2 = p1.transform(motor);
    assert_close!(p2, re2::Vector::point([10., 15.]));
}

#[test]
fn test_3d_translation() {
    let p1 = re3::Vector::point([10., 10., 10.]);

    let motor = re3::translator(re3::Vector::ideal_point([0., 5., 0.]));

    let p2 = p1.transform(motor);
    assert_close!(p2, re3::Vector::point([10., 15., 10.]));
}

#[test]
fn test_2d_translation_from_geo() {
    let p1 = re2::Vector::point([10., 10.]);

    let motor = re2::Vector::point([50., 50.]).motor_to(re2::Vector::point([50., 52.5]));

    let p2 = p1.transform(motor);
    assert_close!(p2, re2::Vector::point([10., 15.]));
}

#[test]
fn test_2d_rotation_from_geo() {
    let p = re2::Vector::point([1., 0.]);

    let rot_l1 = re2::Vector::origin()
        .join(re2::Vector::point([10., 0.]))
        .unitized();
    let rot_l2 = re2::Vector::origin()
        .join(re2::Vector::point([10., 10.]))
        .unitized();
    let motor = rot_l1.motor_to(rot_l2);

    assert_close!(p.transform(motor), re2::Vector::point([0., 1.]));
}

#[test]
fn test_2d_compose_rotations() {
    let p = re2::Vector::point([10., 10.]);

    let translate_up_5 = re2::translator(re2::Vector::ideal_point([0., 5.]));
    let rotate_90 = re2::axis_angle(re2::Vector::origin(), 0.25 * core::f32::consts::TAU);

    let up_then_rotate = translate_up_5.compose(rotate_90);
    assert_close!(p.transform(up_then_rotate), re2::Vector::point([-15., 10.]));

    let rotate_then_up = rotate_90.compose(translate_up_5);
    assert_close!(p.transform(rotate_then_up), re2::Vector::point([-10., 15.]));
}

#[test]
fn test_2d_motor_between_parallel_lines() {
    let rot_l1 = re2::Vector::origin()
        .join(re2::Vector::point([10., 0.]))
        .unitized();
    let rot_l2 = re2::Vector::point([0., 10.])
        .join(re2::Vector::point([10., 10.]))
        .unitized();
    let motor = rot_l1.motor_to(rot_l2);

    let p = re2::Vector::origin();
    assert_close!(p.transform(motor), re2::Vector::point([0., 20.]));
}

#[test]
fn test_2d_superset_orthogonal_to() {
    let l = re2::Vector::origin()
        .join(re2::Vector::point([1., 1.]))
        .unitized();
    let p = re2::Vector::point([0., 2.]);
    let l2 = p.superset_orthogonal_to(l);

    let expected_l2 = re2::Vector::point([1., 1.]).join(p).unitized();

    assert_close!(l2, expected_l2);
}

#[test]
fn test_3d_superset_orthogonal_to() {
    let pl = re3::Vector::origin()
        .join(re3::Vector::point([0., 0., 10.]))
        .join(re3::Vector::point([1., 1., 10.]))
        .unitized();
    let p = re3::Vector::point([0., 2., 10.]);
    let l = p.superset_orthogonal_to(pl);

    let expected_l = re3::Vector::point([1., 1., 10.]).join(p).unitized();

    assert_close!(l, expected_l);
}

#[test]
fn test_2d_projection() {
    let l = re2::Vector::origin()
        .join(re2::Vector::point([1., 1.]))
        .unitized();
    let p = re2::Vector::point([0., 2.]);
    let p2 = p.projection(l).unitized();

    assert_close!(p2, re2::Vector::point([1., 1.]));
}

#[test]
fn test_3d_projection_pt_onto_plane() {
    let pl = re3::Vector::origin()
        .join(re3::Vector::point([0., 0., 10.]))
        .join(re3::Vector::point([1., 1., 10.]))
        .unitized();
    let p = re3::Vector::point([0., 2., 10.]);
    let p2 = p.projection(pl).unitized();

    assert_close!(p2, re3::Vector::point([1., 1., 10.]));
}

#[test]
fn test_2d_anti_projection() {
    let l = re2::Vector::origin()
        .join(re2::Vector::point([1., 1.]))
        .unitized();
    let p = re2::Vector::point([0., 2.]);
    let l2 = l.anti_projection(p).unitized();

    let expected_l2 = p.join(re2::Vector::ideal_point([1., 1.])).unitized();
    assert_close!(l2, expected_l2);
}

#[test]
fn test_3d_anti_projection_plane_onto_pt() {
    let pl = re3::Vector::origin()
        .join(re3::Vector::point([0., 0., 10.]))
        .join(re3::Vector::point([1., 1., 10.]))
        .unitized();
    let p = re3::Vector::point([0., 2., 10.]);
    let pl2 = pl.anti_projection(p).unitized();

    let expected_pl2 = p
        .join(re3::Vector::ideal_point([0., 0., 1.]))
        .join(re3::Vector::ideal_point([1., 1., 10.]))
        .unitized();

    assert_close!(pl2, expected_pl2);
}

#[test]
fn test_3d_projection_line_onto_plane() {
    let pl = re3::Vector::point([0., 0., 10.])
        .join(re3::Vector::point([1., 0., 10.]))
        .join(re3::Vector::point([0., 1., 10.]))
        .unitized();
    let l = re3::Vector::point([0., 20., 20.])
        .join(re3::Vector::point([1., 20., 20.]))
        .unitized();

    let l2 = l.projection(pl).unitized();

    let expected_l2 = re3::Vector::point([0., 20., 10.])
        .join(re3::Vector::point([1., 20., 10.]))
        .unitized();

    assert_close!(l2, expected_l2);
}

#[test]
fn test_3d_central_projection_line_onto_plane() {
    let pl = re3::Vector::point([0., 0., 10.])
        .join(re3::Vector::point([1., 0., 10.]))
        .join(re3::Vector::point([0., 1., 10.]))
        .unitized();
    let l = re3::Vector::point([0., 20., 20.])
        .join(re3::Vector::point([1., 20., 20.]))
        .unitized();

    let l2 = l.central_projection(pl).unitized();

    let expected_l2 = re3::Vector::point([0., 10., 10.])
        .join(re3::Vector::point([1., 10., 10.]))
        .unitized();

    assert_close!(l2, expected_l2);
}

#[test]
fn test_linear_operator() {
    // Non-uniform scale + translate
    let m: re2::LinearOperator<f32> = re2::LinearOperator {
        x: re2::Vector::x_hat(),
        y: re2::Vector::y_hat() * 2.,
        w: re2::Vector::origin() + re2::Vector::x_hat() * 3.,
    };

    println!("{:?}", m);

    let l1 = re2::Vector::point([0., 0.]).join(re2::Vector::point([1., 1.]));
    let l2 = re2::Vector::point([10., 0.]).join(re2::Vector::point([9., 1.]));

    let transform_lines = l1.transform(m).meet(l2.transform(m)).unitized();
    let transform_points = l1.meet(l2).transform(m).unitized();
    assert_close!(transform_lines, transform_points);
}

#[test]
fn test_linear_operator_rotation() {
    // Rotation about Z by 90 degrees expressed as an AntiEven

    let motor = re3::axis_angle(
        re3::Vector::origin().join(re3::Vector::z_hat()),
        0.25 * core::f32::consts::TAU,
    );

    let matrix = re3::LinearOperator {
        x: re3::Vector::y_hat(),
        y: -re3::Vector::x_hat(),
        z: re3::Vector::z_hat(),
        w: re3::Vector::origin(),
    };

    // Rotate a point
    let p = re3::Vector::point([3., 4., 5.]);
    assert_close!(p.transform(motor), p.transform(matrix));

    // Rotate a line
    let l = re3::Vector::point([10., 0., 0.]).join(re3::Vector::point([9., 1., 3.]));
    assert_close!(l.transform(motor), l.transform(matrix));
}

/*
#[test]
fn test_angle() {
    // Angle between lines in 2D
    assert_close!(
        angle(
            re2::Vector::<f32>::origin().vee(re2::Vector::point([0., 5.])),
            re2::Vector::origin().vee(re2::Vector::point([-10., 10.]))
        ),
        0.125 * core::f32::consts::TAU,
    );

    // Angle between planes in 3D
    assert_close!(
        angle(
            re3::Vector::<f32>::origin()
                .vee(re3::Vector::point([0., 0., 1.]))
                .vee(re3::Vector::point([0., 5., 0.])),
            re3::Vector::origin()
                .vee(re3::Vector::point([0., 0., 1.]))
                .vee(re3::Vector::point([-10., 10., 0.]))
        ),
        0.125 * core::f32::consts::TAU,
    );

    {
        // Angle between line and plane
        let pl = re3::Vector::<f32>::origin()
            .vee(re3::Vector::point([1., 0., 0.]))
            .vee(re3::Vector::point([0., 1., 0.]))
            .hat();
        let l = re3::Vector::point([10., 10., 0.])
            .vee(re3::Vector::point([10., 20., 10.]))
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
    let p1 = re2::Vector::point([10., 10.]);

    let center = re2::Vector::ideal_point([1., 0.]) * -2.5;
    let motor = center.exp();

    let p2 = p1.transform(motor);
    assert!(p2.is_close(re2::Vector::point([10., 15.])));

    // 2D rotation
    let p1 = re2::Vector::point([10., 0.]);

    // This motor goes the wrong way??
    let center = re2::Vector::origin() * (0.125 * core::f32::consts::TAU);
    let motor = center.exp();
    dbg!(motor);

    // This one works
    //let l1 = re2::Vector::origin().vee(re2::Vector::point([5., 5.]));
    //let l2 = re2::Vector::origin().vee(re2::Vector::point([-10., 10.]));
    //let motor = (l2 * l1).hat() + 1.;
    //dbg!(motor);

    let p2 = p1.transform(motor);
    dbg!(p2);
    //assert_close!(p2, re2::Vector::point([0., 10.])); // XXX broken
}
*/

/*
use ngeom_macros::geometric_algebra;

#[geometric_algebra()]
mod pga_test {
    struct Vector<T> {}
}
*/
