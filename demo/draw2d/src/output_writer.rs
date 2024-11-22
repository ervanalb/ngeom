use crate::{AnimationDescriptor, Output, OutputTargetDescriptor};
use std::fs;
use std::io;
use std::io::prelude::*;
use std::process::{Child, ChildStdin, Command, Stdio};
use tempfile::NamedTempFile;

#[derive(Debug)]
pub enum OutputWriterError {
    OsError(io::Error),
    FfmpegError(String),
}

pub struct OutputWriter {
    subprocess: Child,
    input: ChildStdin,
    stderr_tempfile: NamedTempFile,
    output_tempfile: Option<NamedTempFile>,
}

impl OutputWriter {
    pub fn new(descriptor: &AnimationDescriptor) -> Result<OutputWriter, OutputWriterError> {
        let size_option = format!("{}x{}", descriptor.resolution[0], descriptor.resolution[1]);
        let framerate_option = format!("{}", descriptor.fps);
        let input_options = vec![
            "-f",
            "rawvideo",
            "-video_size",
            &size_option,
            "-pixel_format",
            "rgb24",
            "-framerate",
            &framerate_option,
            "-i",
            "-",
        ];
        let mut output_options = vec!["-c:v", "h264", "-pixel_format", "yuv420p"];

        let stderr_tempfile =
            NamedTempFile::with_suffix(".log").map_err(|e| OutputWriterError::OsError(e))?;

        let output_tempfile = match &descriptor.output_target {
            OutputTargetDescriptor::File(_) => None,
            OutputTargetDescriptor::Memory => Some(
                NamedTempFile::with_suffix(".mp4").map_err(|e| OutputWriterError::OsError(e))?,
            ),
        };

        match &descriptor.output_target {
            OutputTargetDescriptor::File(filename) => {
                output_options.push(filename);
            }
            OutputTargetDescriptor::Memory => {
                let path = output_tempfile.as_ref().unwrap().path().to_str().unwrap();
                output_options.push(path);
            }
        }

        let mut options = vec!["-y"];
        options.extend(input_options);
        options.extend(output_options);

        let stderr = stderr_tempfile
            .reopen()
            .map_err(|e| OutputWriterError::OsError(e))?;

        let mut subprocess = Command::new("ffmpeg")
            .args(&options)
            .stdin(Stdio::piped())
            .stderr(stderr)
            .stdout(match &descriptor.output_target {
                OutputTargetDescriptor::File(_) => Stdio::inherit(),
                OutputTargetDescriptor::Memory => Stdio::piped(),
            })
            .spawn()
            .map_err(|e| OutputWriterError::OsError(e))?;

        let input = subprocess.stdin.take().unwrap();

        Ok(Self {
            subprocess,
            input,
            stderr_tempfile,
            output_tempfile,
        })
    }

    fn ffmpeg_stderr(f: &NamedTempFile) -> OutputWriterError {
        let path = f.path().to_str().unwrap();
        match fs::read_to_string(path) {
            Ok(stderr) => OutputWriterError::FfmpegError(stderr),
            Err(os_err) => OutputWriterError::OsError(os_err),
        }
    }

    pub fn write(&mut self, data: &[u8]) -> Result<(), OutputWriterError> {
        self.input
            .write_all(data)
            .map_err(|_| Self::ffmpeg_stderr(&self.stderr_tempfile))
    }

    pub fn close(mut self) -> Result<Output, OutputWriterError> {
        drop(self.input); // Close stdin

        let rc = self.subprocess.wait().expect("Failed to finish subprocess");

        if !rc.success() {
            panic!("Ffmpeg finished with error code {:?}", rc);
        }

        Ok(match self.output_tempfile {
            Some(output_tempfile) => {
                let path = output_tempfile.path().to_str().unwrap();
                let bytes =
                    fs::read(path).map_err(|_| Self::ffmpeg_stderr(&self.stderr_tempfile))?;
                Output::Memory(bytes)
            }
            None => Output::File,
        })
    }
}
