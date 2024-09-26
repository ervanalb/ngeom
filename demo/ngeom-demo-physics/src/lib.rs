use core::ops::{Add, Mul};
use eframe::{egui, epaint};
use std::fmt::Debug;

use ngeom::algebraic_ops::*;
use ngeom::ops::*;
use ngeom::scalar::*;
use ngeom::{re2, re3};

#[cfg(target_arch = "wasm32")]
mod web;

#[cfg(target_arch = "wasm32")]
pub use web::*;

pub struct MyApp<SPACE: Space> {
    space: SPACE,
    physics: PhysicsState<SPACE>,
}

impl<SPACE: Space> Default for MyApp<SPACE> {
    fn default() -> Self {
        MyApp {
            space: SPACE::new(),
            physics: PhysicsState::new(),
        }
    }
}

impl<SPACE: Space> MyApp<SPACE> {
    pub fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        MyApp::default()
    }
}

pub trait Space {
    type Vector: Clone
        + Copy
        + Debug
        + Default
        + Origin
        + Join<Self::Vector, Output = Self::Bivector>
        + YHat
        + Mul<f32, Output = Self::Vector>
        + Add<Output = Self::Vector>
        + TransformInverse<Self::AntiEven, Output = Self::Vector>
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

    // Convert VECTOR to Vec2 for egui
    fn into_vec2(&self, x: Self::Vector) -> egui::Vec2;

    fn new() -> Self;
}

pub struct Space2d {
    vertices: Vec<re2::Vector<f32>>,
    edges: Vec<(usize, usize)>,
}

impl Space for Space2d {
    type Vector = re2::Vector<f32>;
    type Bivector = re2::Bivector<f32>;
    type AntiEven = re2::AntiEven<f32>;

    fn cube_vertices(&self) -> &[re2::Vector<f32>] {
        &self.vertices
    }
    fn cube_edges(&self) -> &[(usize, usize)] {
        &self.edges
    }
    fn into_vec2(&self, p: re2::Vector<f32>) -> egui::Vec2 {
        egui::vec2(p.x, p.y)
    }
    fn new() -> Self {
        Space2d {
            vertices: vec![
                re2::Vector::point([-0.5, -0.5]),
                re2::Vector::point([-0.5, 0.5]),
                re2::Vector::point([0.5, -0.5]),
                re2::Vector::point([0.5, 0.5]),
            ],
            edges: vec![(0, 1), (1, 3), (3, 2), (2, 0)],
        }
    }
}

pub struct Space3d {
    vertices: Vec<re3::Vector<f32>>,
    edges: Vec<(usize, usize)>,
}

impl Space for Space3d {
    type Vector = re3::Vector<f32>;
    type Bivector = re3::Bivector<f32>;
    type AntiEven = re3::AntiEven<f32>;

    fn cube_vertices(&self) -> &[re3::Vector<f32>] {
        &self.vertices
    }
    fn cube_edges(&self) -> &[(usize, usize)] {
        &self.edges
    }
    fn into_vec2(&self, p: re3::Vector<f32>) -> egui::Vec2 {
        egui::vec2(p.x, p.y)
    }
    fn new() -> Self {
        Space3d {
            vertices: vec![
                re3::Vector::point([-0.5, -0.5, -0.5]),
                re3::Vector::point([-0.5, -0.5, 0.5]),
                re3::Vector::point([-0.5, 0.5, -0.5]),
                re3::Vector::point([-0.5, 0.5, 0.5]),
                re3::Vector::point([0.5, -0.5, -0.5]),
                re3::Vector::point([0.5, -0.5, 0.5]),
                re3::Vector::point([0.5, 0.5, -0.5]),
                re3::Vector::point([0.5, 0.5, 0.5]),
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
        let g_b = g_s.transform_inverse(self.cube_g);

        p_b.join(g_b)
    }

    fn spring_torque(&self, space: &SPACE) -> SPACE::Bivector {
        let p_b = space.cube_vertices()[0];
        let a_s = SPACE::Vector::origin(); // World attachment point, in space frame
        let k = -0.8; // Spring constant

        let a_b = a_s.transform_inverse(self.cube_g);
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

    fn draw(&self, space: &SPACE, painter: &egui::Painter, clip_rect: egui::Rect) {

        let clip_center = clip_rect.center();
        let clip_scale = clip_rect.width().min(clip_rect.height());

        let into_screen_space = |p| {
            let mut p2 = space.into_vec2(p);
            // Adjust camera
            p2.y += 1.8;
            p2 *= 0.2;
            p2.y *= -1.;
            // Convert to egui coordinates
            p2 = p2 * clip_scale;
            clip_center + p2
        };

        let vertices_transformed: Vec<_> = space
            .cube_vertices()
            .iter()
            .map(|p| {
                let p = p.transform(self.cube_g);
                into_screen_space(p)
            })
            .collect();

        // Draw cube
        for &(i1, i2) in space.cube_edges() {
            painter.line_segment(
                [vertices_transformed[i1], vertices_transformed[i2]],
                epaint::Stroke::new(2., epaint::Color32::DARK_BLUE),
            );
        }

        // Draw spring
        painter.line_segment(
            [
                vertices_transformed[0],
                into_screen_space(SPACE::Vector::origin()),
            ],
            epaint::Stroke::new(2., epaint::Color32::DARK_RED),
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

impl<SPACE: Space> MyApp<SPACE> {
    fn step(&mut self) {
        self.physics = rk4(|y| y.d_dt(&self.space), self.physics, 0.06_f32);
        self.physics = self.physics.unitized();
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

impl<SPACE: Space> eframe::App for MyApp<SPACE> {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.step();

        ctx.set_visuals(egui::Visuals::light());
        egui::CentralPanel::default().show(ctx, |ui| {
            let (_, rect) = ui.allocate_space(ui.available_size());
            let painter = ui.painter().with_clip_rect(rect);
            self.physics.draw(&self.space, &painter, rect);
        });
        ctx.request_repaint();
    }
}
