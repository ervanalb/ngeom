use kiss3d::camera::{ArcBall, Camera};
use kiss3d::light::Light;
use kiss3d::nalgebra::{Point3, Vector3};
use kiss3d::planar_camera::PlanarCamera;
use kiss3d::post_processing::post_processing_effect::PostProcessingEffect;
use kiss3d::renderer::Renderer;
use kiss3d::window::{State, Window};

use ngeom::algebraic_ops::*;
use ngeom::ops::*;
use ngeom::scalar::*;
use ngeom::{pga2d, pga3d};

use core::ops::{Add, Mul};

struct AppState {
    arc_ball_camera: ArcBall,
    physics: PhysicsState3D,
}

#[derive(Clone, Copy, Debug)]
struct PhysicsState2D {
    cube_g: pga2d::AntiEven<f32>, // g, The cube's orientation (& position) in space
    cube_l_b: pga2d::Bivector<f32>, // L_B, the cube's angular (& linear) momentum in its local reference
}

fn d_dt<
    T: Rational,
    ANTIEVEN: Copy + core::ops::Mul<T, Output = ANTIEVEN>,
    BIVECTOR: Copy
        + LeftComplement<
            Output: Copy
                        + AntiWedgeDot<ANTIEVEN, Output = ANTIEVEN>
                        + AntiCommutator<BIVECTOR, Output = BIVECTOR>,
        > + Add<BIVECTOR, Output = BIVECTOR>,
>(
    (g, l_b): (ANTIEVEN, BIVECTOR),
    w_b: <BIVECTOR as LeftComplement>::Output,
    t_b: BIVECTOR, // External torques
) -> (ANTIEVEN, BIVECTOR) {
    (
        w_b.anti_wedge_dot(g) * T::one_half(),
        w_b.anti_commutator(l_b) + t_b,
    )
}

fn spring_torque<
    T: Ring,
    VECTOR: Copy + ReverseTransform<ANTIEVEN, Output = VECTOR> + Join<VECTOR, Output = BIVECTOR>,
    BIVECTOR: Copy + core::ops::Mul<T, Output = BIVECTOR>,
    ANTIEVEN,
>(
    p_b: VECTOR,
    a_s: VECTOR,
    g: ANTIEVEN,
    k: T,
) -> BIVECTOR {
    let a_b = a_s.reverse_transform(g);
    let spring_line = a_b.join(p_b);
    spring_line * k
}

fn gravity_torque<
    VECTOR: Copy + ReverseTransform<ANTIEVEN, Output = VECTOR> + Join<VECTOR, Output = BIVECTOR>,
    BIVECTOR,
    ANTIEVEN,
>(
    g_s: VECTOR,
    p_b: VECTOR,
    g: ANTIEVEN,
) -> BIVECTOR {
    let g_b = g_s.reverse_transform(g);

    p_b.join(g_b)
}

impl PhysicsState2D {
    fn new() -> PhysicsState2D {
        PhysicsState2D {
            cube_g: pga2d::identity_motor().into(),
            cube_l_b: Default::default(),
        }
    }

    // Calculate the current angular velocity ω_B
    // from the current angular momentum L_B
    // by applying the inverse of the inertia tensor
    fn w_b(&self) -> pga2d::Vector<f32> {
        // TODO add inertia tensor
        self.cube_l_b.left_complement()
    }

    fn d_dt(&self) -> PhysicsState2D {
        let w_b = self.w_b();
        let t_b = self.gravity_torque() + self.spring_torque() + self.drag_torque();

        let (g_dot, l_b_dot) = d_dt::<f32, _, _>((self.cube_g, self.cube_l_b), w_b, t_b);

        PhysicsState2D {
            cube_g: g_dot,
            cube_l_b: l_b_dot,
        }
    }

    fn gravity_torque(&self) -> pga2d::Bivector<f32> {
        let g_s = pga2d::ideal_point([0., -1.]);
        let p_b = pga2d::origin(); // Apply gravity to center of mass
        gravity_torque(g_s, p_b, self.cube_g)
    }

    fn spring_torque(&self) -> pga2d::Bivector<f32> {
        let p_b = pga2d::point([0.5, 0.5]); // Cube attachment point, in body frame
        let a_s = pga2d::origin(); // World attachment point, in space frame
        let k = -2.; // Spring constant
        spring_torque(p_b, a_s, self.cube_g, k)
    }

    fn drag_torque(&self) -> pga2d::Bivector<f32> {
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

    fn draw(&self, window: &mut Window) {
        let points: Vec<pga2d::Vector<f32>> = vec![
            pga2d::point([-0.5, -0.5]),
            pga2d::point([-0.5, 0.5]),
            pga2d::point([0.5, 0.5]),
            pga2d::point([0.5, -0.5]),
        ];

        let edges: Vec<(usize, usize)> = vec![(0, 1), (1, 2), (2, 3), (3, 0)];

        let points_transformed: Vec<_> = points.iter().map(|p| p.transform(self.cube_g)).collect();

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
}

#[derive(Clone, Copy, Debug)]
struct PhysicsState3D {
    cube_g: pga3d::AntiEven<f32>, // g, The cube's orientation (& position) in space
    cube_l_b: pga3d::Bivector<f32>, // L_B, the cube's angular (& linear) momentum in its local reference
}

impl PhysicsState3D {
    fn new() -> PhysicsState3D {
        PhysicsState3D {
            cube_g: pga3d::identity_motor().into(),
            cube_l_b: Default::default(),
        }
    }

    // Calculate the current angular velocity ω_B
    // from the current angular momentum L_B
    // by applying the inverse of the inertia tensor
    fn w_b(&self) -> pga3d::Bivector<f32> {
        // TODO add inertia tensor
        self.cube_l_b.left_complement()
    }

    fn d_dt(&self) -> PhysicsState3D {
        let w_b = self.w_b();
        let t_b = self.gravity_torque() + self.spring_torque() + self.drag_torque();

        let (g_dot, l_b_dot) = d_dt::<f32, _, _>((self.cube_g, self.cube_l_b), w_b, t_b);

        PhysicsState3D {
            cube_g: g_dot,
            cube_l_b: l_b_dot,
        }
    }

    fn gravity_torque(&self) -> pga3d::Bivector<f32> {
        let g_s = pga3d::ideal_point([0., -1., 0.]);
        let p_b = pga3d::origin(); // Apply gravity to center of mass
        gravity_torque(g_s, p_b, self.cube_g)
    }

    fn spring_torque(&self) -> pga3d::Bivector<f32> {
        let p_b = pga3d::point([0.5, 0.5, 0.5]); // Cube attachment point, in body frame
        let a_s = pga3d::origin(); // World attachment point, in space frame
        let k = -2.; // Spring constant
        spring_torque(p_b, a_s, self.cube_g, k)
    }

    fn drag_torque(&self) -> pga3d::Bivector<f32> {
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

    fn draw(&self, window: &mut Window) {
        let points: Vec<pga3d::Vector<f32>> = vec![
            pga3d::point([-0.5, -0.5, -0.5]),
            pga3d::point([-0.5, -0.5, 0.5]),
            pga3d::point([-0.5, 0.5, -0.5]),
            pga3d::point([-0.5, 0.5, 0.5]),
            pga3d::point([0.5, -0.5, -0.5]),
            pga3d::point([0.5, -0.5, 0.5]),
            pga3d::point([0.5, 0.5, -0.5]),
            pga3d::point([0.5, 0.5, 0.5]),
        ];

        let edges: Vec<(usize, usize)> = vec![
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
        ];

        let points_transformed: Vec<_> = points.iter().map(|p| p.transform(self.cube_g)).collect();

        let get_point3 = |i: usize| -> Point3<f32> {
            let v = points_transformed[i].unitized();
            Point3::new(v.a1, v.a2, v.a3)
        };

        // Draw cube
        for (i1, i2) in edges {
            window.draw_line(&get_point3(i1), &get_point3(i2), &Point3::new(1., 1., 1.));
        }

        // Draw spring
        window.draw_line(
            &Point3::new(0., 0., 0.),
            &get_point3(7),
            &Point3::new(1., 1., 1.),
        );
    }
}

impl Add for PhysicsState2D {
    type Output = PhysicsState2D;
    fn add(self, r: PhysicsState2D) -> PhysicsState2D {
        PhysicsState2D {
            cube_g: self.cube_g + r.cube_g,
            cube_l_b: self.cube_l_b + r.cube_l_b,
        }
    }
}

impl Add for PhysicsState3D {
    type Output = PhysicsState3D;
    fn add(self, r: PhysicsState3D) -> PhysicsState3D {
        PhysicsState3D {
            cube_g: self.cube_g + r.cube_g,
            cube_l_b: self.cube_l_b + r.cube_l_b,
        }
    }
}

impl Mul<f32> for PhysicsState2D {
    type Output = PhysicsState2D;
    fn mul(self, r: f32) -> PhysicsState2D {
        PhysicsState2D {
            cube_g: self.cube_g * r,
            cube_l_b: self.cube_l_b * r,
        }
    }
}

impl Mul<f32> for PhysicsState3D {
    type Output = PhysicsState3D;
    fn mul(self, r: f32) -> PhysicsState3D {
        PhysicsState3D {
            cube_g: self.cube_g * r,
            cube_l_b: self.cube_l_b * r,
        }
    }
}

impl Unitized for PhysicsState2D {
    type Output = PhysicsState2D;
    fn unitized(self) -> PhysicsState2D {
        PhysicsState2D {
            cube_g: self.cube_g.unitized(),
            cube_l_b: self.cube_l_b,
        }
    }
}

impl Unitized for PhysicsState3D {
    type Output = PhysicsState3D;
    fn unitized(self) -> PhysicsState3D {
        PhysicsState3D {
            cube_g: self.cube_g.unitized(),
            cube_l_b: self.cube_l_b,
        }
    }
}

impl State for AppState {
    fn step(&mut self, window: &mut Window) {
        self.physics = rk4(|y| y.d_dt(), self.physics, 0.06_f32);
        self.physics = self.physics.unitized();
        draw_axes(window);
        self.physics.draw(window);
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
        physics: PhysicsState3D::new(),
    };

    window.render_loop(state)
}
