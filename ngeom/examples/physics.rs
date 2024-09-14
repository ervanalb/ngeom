use kiss3d::camera::{ArcBall, Camera};
use kiss3d::light::Light;
use kiss3d::nalgebra::{Point3, Vector3};
use kiss3d::planar_camera::PlanarCamera;
use kiss3d::post_processing::post_processing_effect::PostProcessingEffect;
use kiss3d::renderer::Renderer;
use kiss3d::window::{State, Window};

use ngeom::ops::*;
use ngeom::pga2d::*;

struct AppState {
    arc_ball_camera: ArcBall,
    cube_pose: AntiEven<f32>,
}

impl State for AppState {
    fn step(&mut self, window: &mut Window) {
        self.cube_pose = self.cube_pose.compose(axis_angle(origin(), 0.06_f32));
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
        .map(|p| p.transform(state.cube_pose))
        .collect();

    let get_point3 = |i: usize| -> Point3<f32> {
        let v = points_transformed[i].unitized();
        Point3::new(v.a1, v.a2, 0.)
    };

    for (i1, i2) in edges {
        window.draw_line(&get_point3(i1), &get_point3(i2), &Point3::new(1., 1., 1.));
    }
}

fn main() {
    let mut window = Window::new("ngeom N-dimensional Physics Demo");

    window.set_light(Light::StickToCamera);

    let mut arc_ball = ArcBall::new(Point3::new(3., -10., 3.), Point3::origin());
    arc_ball.set_up_axis(Vector3::new(0., 0., 1.));

    let state = AppState {
        arc_ball_camera: arc_ball,
        cube_pose: identity_motor().into(),
    };

    window.render_loop(state)
}
