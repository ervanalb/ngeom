use kiss3d::camera::{ArcBall, Camera};
use kiss3d::light::Light;
use kiss3d::nalgebra::{Point3, Vector3};
use kiss3d::planar_camera::PlanarCamera;
use kiss3d::post_processing::post_processing_effect::PostProcessingEffect;
use kiss3d::renderer::Renderer;
use kiss3d::window::{State, Window};

use core::ops::{Add, Mul};
use ngeom::algebraic_ops::*;
use ngeom::ops::*;
use ngeom::scalar::*;
use ngeom::{pga2d, pga3d};
use std::fmt::Debug;

struct AppState<SPACE: Space + 'static> {
    arc_ball_camera: ArcBall,
    space: SPACE,
    physics: PhysicsState<SPACE>,
}

trait Space {
    type Vector: Clone
        + Copy
        + Debug
        + Default
        + Origin
        + Join<Self::Vector, Output = Self::Bivector>
        + YHat
        + Mul<f32, Output = Self::Vector>
        + Add<Output = Self::Vector>
        + ReverseTransform<Self::AntiEven, Output = Self::Vector>
        + Transform<Self::AntiEven, Output = Self::Vector>;
    type Bivector: Clone
        + Copy
        + Debug
        + Default
        + LeftComplement<
            Output: Copy
                        + AntiWedgeDot<Self::AntiEven, Output = Self::AntiEven>
                        + AntiCommutator<Self::Bivector, Output = Self::Bivector>
                        + RightComplement<Output = Self::Bivector>,
        > + Add<Output = Self::Bivector>
        + Mul<f32, Output = Self::Bivector>;
    type AntiEven: Clone
        + Copy
        + Debug
        + IdentityMotor
        + Mul<f32, Output = Self::AntiEven>
        + Add<Output = Self::AntiEven>
        + Unitized<Output = Self::AntiEven>;

    // Constants
    fn cube_vertices(&self) -> &[Self::Vector];
    fn cube_edges(&self) -> &[(usize, usize)];

    // Convert VECTOR to Point3
    fn into_point3(&self, x: Self::Vector) -> Point3<f32>;
}

struct Space2d {
    vertices: Vec<pga2d::Vector<f32>>,
    edges: Vec<(usize, usize)>,
}

impl Space2d {
    fn new() -> Self {
        Space2d {
            vertices: vec![
                pga2d::Vector::point([-0.5, -0.5]),
                pga2d::Vector::point([-0.5, 0.5]),
                pga2d::Vector::point([0.5, -0.5]),
                pga2d::Vector::point([0.5, 0.5]),
            ],
            edges: vec![(0, 1), (1, 3), (3, 2), (2, 0)],
        }
    }
}

impl Space for Space2d {
    type Vector = pga2d::Vector<f32>;
    type Bivector = pga2d::Bivector<f32>;
    type AntiEven = pga2d::AntiEven<f32>;

    fn cube_vertices(&self) -> &[pga2d::Vector<f32>] {
        &self.vertices
    }
    fn cube_edges(&self) -> &[(usize, usize)] {
        &self.edges
    }
    fn into_point3(&self, p: pga2d::Vector<f32>) -> Point3<f32> {
        Point3::new(p.x, p.y, 0.)
    }
}

struct Space3d {
    vertices: Vec<pga3d::Vector<f32>>,
    edges: Vec<(usize, usize)>,
}

impl Space3d {
    fn new() -> Self {
        Space3d {
            vertices: vec![
                pga3d::Vector::point([-0.5, -0.5, -0.5]),
                pga3d::Vector::point([-0.5, -0.5, 0.5]),
                pga3d::Vector::point([-0.5, 0.5, -0.5]),
                pga3d::Vector::point([-0.5, 0.5, 0.5]),
                pga3d::Vector::point([0.5, -0.5, -0.5]),
                pga3d::Vector::point([0.5, -0.5, 0.5]),
                pga3d::Vector::point([0.5, 0.5, -0.5]),
                pga3d::Vector::point([0.5, 0.5, 0.5]),
            ],
            edges: vec![
                (0, 1),
                (1, 3),
                (3, 2),
                (2, 0),
                (4, 5),
                (5, 7),
                (7, 6),
                (6, 4),
                (0, 4),
                (1, 5),
                (2, 6),
                (3, 7),
            ],
        }
    }
}

impl Space for Space3d {
    type Vector = pga3d::Vector<f32>;
    type Bivector = pga3d::Bivector<f32>;
    type AntiEven = pga3d::AntiEven<f32>;

    fn cube_vertices(&self) -> &[pga3d::Vector<f32>] {
        &self.vertices
    }
    fn cube_edges(&self) -> &[(usize, usize)] {
        &self.edges
    }
    fn into_point3(&self, p: pga3d::Vector<f32>) -> Point3<f32> {
        Point3::new(p.x, p.y, p.z)
    }
}

#[derive(Debug)]
struct PhysicsState<SPACE: Space> {
    cube_g: SPACE::AntiEven,   // g, The cube's orientation (& position) in space
    cube_l_b: SPACE::Bivector, // L_B, the cube's angular (& linear) momentum in its local reference
}

// Handwritten clone required due to a quirk in derive(Clone)
// see: https://smallcultfollowing.com/babysteps/blog/2022/04/12/implied-bounds-and-perfect-derive/
impl<SPACE: Space> Clone for PhysicsState<SPACE> {
    fn clone(&self) -> Self {
        PhysicsState {
            cube_g: self.cube_g,
            cube_l_b: self.cube_l_b,
        }
    }
}

impl<SPACE: Space> Copy for PhysicsState<SPACE> {}

impl<SPACE: Space> PhysicsState<SPACE> {
    fn new() -> Self {
        Self {
            cube_g: SPACE::AntiEven::identity_motor(),
            cube_l_b: Default::default(),
        }
    }

    // Calculate the current angular velocity ω_B
    // from the current angular momentum L_B
    // by applying the inverse of the inertia tensor
    fn w_b(&self) -> <SPACE::Bivector as LeftComplement>::Output {
        // TODO add inertia tensor
        self.cube_l_b.left_complement()
    }

    fn d_dt(&self, space: &SPACE) -> Self {
        let w_b = self.w_b();
        let t_b = self.gravity_torque(space) + self.spring_torque(space) + self.drag_torque(space);

        Self {
            cube_g: w_b.anti_wedge_dot(self.cube_g) * 0.5,
            cube_l_b: w_b.anti_commutator(self.cube_l_b) + t_b,
        }
    }

    fn gravity_torque(&self, _space: &SPACE) -> SPACE::Bivector {
        let g_s = SPACE::Vector::y_hat() * -1.;
        let p_b = SPACE::Vector::origin(); // Apply gravity to center of mass
        let g_b = g_s.reverse_transform(self.cube_g);

        p_b.join(g_b)
    }

    fn spring_torque(&self, space: &SPACE) -> SPACE::Bivector {
        let p_b = space.cube_vertices()[0];
        let a_s = SPACE::Vector::origin(); // World attachment point, in space frame
        let k = -0.8; // Spring constant

        let a_b = a_s.reverse_transform(self.cube_g);
        let spring_line = a_b.join(p_b);
        spring_line * k
    }

    fn drag_torque(&self, _space: &SPACE) -> SPACE::Bivector {
        let w_b = self.w_b();

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

    fn draw(&self, space: &SPACE, window: &mut Window) {
        let vertices_transformed: Vec<_> = space
            .cube_vertices()
            .iter()
            .map(|p| p.transform(self.cube_g))
            .collect();

        let get_point3 = |i| space.into_point3(vertices_transformed[i]);

        // Draw cube
        for &(i1, i2) in space.cube_edges() {
            window.draw_line(&get_point3(i1), &get_point3(i2), &Point3::new(1., 1., 1.));
        }

        // Draw spring
        window.draw_line(
            &Point3::new(0., 0., 0.),
            &get_point3(0),
            &Point3::new(1., 1., 1.),
        );
    }
}

impl<SPACE: Space> Add for PhysicsState<SPACE> {
    type Output = PhysicsState<SPACE>;
    fn add(self, r: PhysicsState<SPACE>) -> PhysicsState<SPACE> {
        PhysicsState {
            cube_g: self.cube_g + r.cube_g,
            cube_l_b: self.cube_l_b + r.cube_l_b,
        }
    }
}

impl<SPACE: Space> Mul<f32> for PhysicsState<SPACE> {
    type Output = PhysicsState<SPACE>;
    fn mul(self, r: f32) -> PhysicsState<SPACE> {
        PhysicsState {
            cube_g: self.cube_g * r,
            cube_l_b: self.cube_l_b * r,
        }
    }
}

impl<SPACE: Space> Unitized for PhysicsState<SPACE> {
    type Output = PhysicsState<SPACE>;
    fn unitized(self) -> PhysicsState<SPACE> {
        PhysicsState {
            cube_g: self.cube_g.unitized(),
            cube_l_b: self.cube_l_b,
        }
    }
}

impl<SPACE: Space> State for AppState<SPACE> {
    fn step(&mut self, window: &mut Window) {
        self.physics = rk4(|y| y.d_dt(&self.space), self.physics, 0.06_f32);
        self.physics = self.physics.unitized();
        draw_axes(window);
        self.physics.draw(&self.space, window);
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

fn main() {
    let mut window = Window::new("ngeom N-dimensional Physics Demo");

    window.set_light(Light::StickToCamera);

    let mut arc_ball = ArcBall::new(Point3::new(0., 0., 10.), Point3::origin());
    arc_ball.set_up_axis(Vector3::new(0., 1., 0.));

    let state = AppState {
        arc_ball_camera: arc_ball,
        space: Space2d::new(),
        physics: PhysicsState::new(),
    };

    window.render_loop(state)
}
