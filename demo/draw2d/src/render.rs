use crate::output_writer::{OutputWriter, OutputWriterError};
use crate::{AnimationDescriptor, Output};
use pollster::FutureExt;
use rayon::prelude::*;
use std::sync::{mpsc, Arc, Condvar, Mutex};

pub type LayerIndex = usize;

pub enum LayerResolution {
    MatchCanvas,
}

pub struct FrameDescriptor {
    pub frame: usize,
    pub t: f64,
}

pub struct LayerDescriptor {
    pub resolution: LayerResolution,
}

pub struct Renderer {
    pub animation: AnimationDescriptor,
    pub device: wgpu::Device,
    pub queue: wgpu::Queue,
    pub layer_descriptors: Vec<LayerDescriptor>,
}

impl Renderer {
    pub fn new_from_adapter(adapter: &wgpu::Adapter, animation: AnimationDescriptor) -> Self {
        let (device, queue) = adapter
            .request_device(&Default::default(), None)
            .block_on()
            .unwrap();

        Self {
            animation,
            device,
            queue,
            layer_descriptors: vec![],
        }
    }

    pub fn new(animation: AnimationDescriptor) -> Self {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
            #[cfg(not(target_arch = "wasm32"))]
            backends: wgpu::Backends::PRIMARY,
            #[cfg(target_arch = "wasm32")]
            backends: wgpu::Backends::GL,
            ..Default::default()
        });

        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::default(),
                compatible_surface: None,
                force_fallback_adapter: false,
            })
            .block_on()
            .unwrap();

        Self::new_from_adapter(&adapter, animation)
    }

    pub fn add_layer(&mut self, layer: LayerDescriptor) -> LayerIndex {
        let ix = self.layer_descriptors.len();
        self.layer_descriptors.push(layer);
        ix
    }

    pub fn remove_layer(&mut self, layer_index: LayerIndex) {
        self.layer_descriptors.remove(layer_index);
    }

    fn frame_renderer(&self) -> FrameRenderer {
        FrameRenderer {
            animation: &self.animation,
            device: &self.device,
            queue: &self.queue,
        }
    }

    pub fn render(
        &mut self,
        frame_function: impl Sync + Send + Fn(&FrameRenderer, FrameDescriptor) -> Option<FrameResult>,
    ) -> Result<Output, OutputWriterError> {
        // Fire up the output writer
        let mut output_writer = OutputWriter::new(&self.animation)?;

        // Channel with enough capacity to hold an item from each thread
        let (tx, rx) =
            mpsc::sync_channel(std::thread::available_parallelism().map_or(8, |n| n.get()));

        let gate = Arc::new((Mutex::new(0), Condvar::new()));

        let fps = self.animation.fps;

        rayon::scope(|s| {
            s.spawn(move |_| {
                (0..)
                    .par_bridge()
                    .map_init(
                        || self.frame_renderer(),
                        move |frame_renderer, i| {
                            let frame_descriptor = FrameDescriptor {
                                frame: i,
                                t: i as f64 / fps,
                            };

                            frame_function(&frame_renderer, frame_descriptor)
                                .map(|frame_result| (i, frame_result))
                        },
                    )
                    .while_some()
                    .for_each_with((tx, gate), |(tx, gate), (i, frame_result)| {
                        let (lock, cond) = &**gate;

                        {
                            let mut guard =
                                cond.wait_while(lock.lock().unwrap(), |v| *v < i).unwrap();
                            let _ = tx.send(frame_result);
                            *guard = i + 1;
                        }

                        cond.notify_all();
                    });
            });

            // Write output
            for FrameResult(data) in rx {
                output_writer.write(&data)?;
            }
            output_writer.close()
        })
    }
}

pub struct FrameRenderer<'a> {
    pub animation: &'a AnimationDescriptor,
    pub device: &'a wgpu::Device,
    pub queue: &'a wgpu::Queue,
}

#[derive(Debug)]
pub struct FrameResult(pub Vec<u8>);
