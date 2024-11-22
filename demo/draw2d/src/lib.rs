pub mod output_writer;
pub mod render;
pub mod vector_graphics;

pub enum OutputTargetDescriptor {
    File(String),
    Memory,
}

pub enum Output {
    File,
    Memory(Vec<u8>),
}

impl Output {
    pub fn evcxr_display(&self) {
        use base64::{engine::general_purpose, Engine};

        match self {
            Output::File => {}
            Output::Memory(video_bytes) => {
                let video_b64 = general_purpose::STANDARD.encode(&video_bytes);
                let html = format!(
                    r#"<video src="data:video/mp4;base64,{}" autoplay controls loop></video>"#,
                    video_b64
                );
                println!("EVCXR_BEGIN_CONTENT text/html\n{}\nEVCXR_END_CONTENT", html);
            }
        }
    }
}

pub struct AnimationDescriptor {
    pub resolution: [usize; 2],
    pub fps: f64,
    pub output_target: OutputTargetDescriptor,
}
