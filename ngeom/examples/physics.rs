use kiss3d::camera::{ArcBall, Camera};
use kiss3d::light::Light;
use kiss3d::nalgebra::{Point3, Vector3};
use kiss3d::planar_camera::PlanarCamera;
use kiss3d::post_processing::post_processing_effect::PostProcessingEffect;
use kiss3d::renderer::Renderer;
use kiss3d::window::{State, Window};

use ngeom::algebraic_ops::*;
use ngeom::ops::*;
use ngeom::pga2d::*;
use ngeom::scalar::*;

use core::ops::{Add, Mul};

struct AppState {
    arc_ball_camera: ArcBall,
    physics: PhysicsState,
}

#[derive(Clone, Copy, Debug)]
struct PhysicsState {
    cube_g: AntiEven<f32>,   // g, The cube's orientation (& position) in space
    cube_l_b: Bivector<f32>, // L_B, the cube's angular (& linear) momentum in its local reference
}

// The inverse of the moment of inertia
// Maps L_B to ω_B
fn a_inv(l_b: Bivector<f32>) -> Vector<f32> {
    // TODO add inertia tensor
    l_b.left_complement()
}

impl PhysicsState {
    fn d_dt(&self) -> PhysicsState {
        let w_b = a_inv(self.cube_l_b);
        PhysicsState {
            cube_g: w_b.anti_wedge_dot(self.cube_g) * 0.5,
            cube_l_b: self.gravity_torque()
                + self.spring_torque()
                + self.drag_torque()
                + w_b.anti_commutator(self.cube_l_b),
        }
    }

    fn gravity_torque(&self) -> Bivector<f32> {
        let g_s = ideal_point([0., -1.]);
        let g_b = g_s.reverse_transform(self.cube_g);

        // Apply gravity to the center of mass
        origin().join(g_b)
    }

    fn spring_torque(&self) -> Bivector<f32> {
        let p_b = point([0.5, 0.5]); // Cube attachment point, in body frame
        let a_s = origin(); // World attachment point, in space frame
        let k = -2.; // Spring constant

        let a_b = a_s.reverse_transform(self.cube_g);
        let spring_line = a_b.join(p_b);
        spring_line * k
    }

    fn drag_torque(&self) -> Bivector<f32> {
        let w_b = a_inv(self.cube_l_b);

        //let linear_k = -0.1;
        //let rotational_k = -0.05;

        let drag = w_b.right_complement();

        // TODO different linear & rotational drag
        drag * -0.05
        //Bivector {
        //    a01: drag.a01 * linear_k,
        //    a20: drag.a20 * linear_k,
        //    a12: drag.a12 * rotational_k,
        //}
    }
}

impl Add for PhysicsState {
    type Output = PhysicsState;
    fn add(self, r: PhysicsState) -> PhysicsState {
        PhysicsState {
            cube_g: self.cube_g + r.cube_g,
            cube_l_b: self.cube_l_b + r.cube_l_b,
        }
    }
}

impl Mul<f32> for PhysicsState {
    type Output = PhysicsState;
    fn mul(self, r: f32) -> PhysicsState {
        PhysicsState {
            cube_g: self.cube_g * r,
            cube_l_b: self.cube_l_b * r,
        }
    }
}

impl Unitized for PhysicsState {
    type Output = PhysicsState;
    fn unitized(self) -> PhysicsState {
        PhysicsState {
            cube_g: self.cube_g.unitized(),
            cube_l_b: self.cube_l_b,
        }
    }
}

impl State for AppState {
    fn step(&mut self, window: &mut Window) {
        self.physics = rk4(|y: PhysicsState| y.d_dt(), self.physics, 0.06_f32);
        self.physics = self.physics.unitized();
        draw_axes(window);
        draw_hypercube(self, window);
    }

    fn cameras_and_effect_and_renderer(
        &mut self,
    ) -> (
        Option<&mut dyn Camera>,
        Option<&mut dyn PlanarCamera>,
        Option<&mut dyn Renderer>,
        Option<&mut dyn PostProcessingEffect>,
    ) {
        (Some(&mut self.arc_ball_camera), None, None, None)
    }
}

// Runge-kutta 4th-order integrator
fn rk4<T: Ring + Rational, Y: Copy + Mul<T, Output = Y> + Add<Y, Output = Y>, F: Fn(Y) -> Y>(
    f: F,
    y: Y,
    h: T,
) -> Y {
    let k1 = f(y);
    let k2 = f(y + k1 * (T::one_half() * h));
    let k3 = f(y + k2 * (T::one_half() * h));
    let k4 = f(y + k3 * h);

    y + (k2 + k3 + (k1 + k4) * T::one_half()) * (h * T::one_third())
}

fn draw_axes(window: &mut Window) {
    window.draw_line(
        &Point3::new(0., 0., 0.),
        &Point3::new(1., 0., 0.),
        &Point3::new(1., 0.0, 0.0),
    );
    window.draw_line(
        &Point3::new(0., 0., 0.),
        &Point3::new(0., 1., 0.),
        &Point3::new(0.0, 1., 0.0),
    );
    window.draw_line(
        &Point3::new(0., 0., 0.),
        &Point3::new(0., 0., 1.),
        &Point3::new(0., 0., 1.),
    );
}

fn draw_hypercube(state: &AppState, window: &mut Window) {
    let points: Vec<Vector<f32>> = vec![
        point([-0.5, -0.5]),
        point([-0.5, 0.5]),
        point([0.5, 0.5]),
        point([0.5, -0.5]),
    ];

    let edges: Vec<(usize, usize)> = vec![(0, 1), (1, 2), (2, 3), (3, 0)];

    let points_transformed: Vec<_> = points
        .iter()
        .map(|p| p.transform(state.physics.cube_g))
        .collect();

    let get_point3 = |i: usize| -> Point3<f32> {
        let v = points_transformed[i].unitized();
        Point3::new(v.a1, v.a2, 0.)
    };

    // Draw cube
    for (i1, i2) in edges {
        window.draw_line(&get_point3(i1), &get_point3(i2), &Point3::new(1., 1., 1.));
    }

    // Draw spring
    window.draw_line(
        &Point3::new(0., 0., 0.),
        &get_point3(2),
        &Point3::new(1., 1., 1.),
    );
}

fn main() {
    let mut window = Window::new("ngeom N-dimensional Physics Demo");

    window.set_light(Light::StickToCamera);

    let mut arc_ball = ArcBall::new(Point3::new(0., 0., 10.), Point3::origin());
    arc_ball.set_up_axis(Vector3::new(0., 1., 0.));

    let state = AppState {
        arc_ball_camera: arc_ball,
        physics: PhysicsState {
            cube_g: identity_motor().into(),
            cube_l_b: Bivector {
                a20: 1.,
                a01: 1.,
                a12: 0.1,
            } * 1.,
        },
    };

    window.render_loop(state)
}
