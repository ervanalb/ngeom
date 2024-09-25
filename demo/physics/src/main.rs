use ngeom_demo_physics::*;

fn main() {
    let (window, state) = AppState::<Space2d>::new();

    window.render_loop(state)
}
