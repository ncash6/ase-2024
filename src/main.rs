use std::{fs::File, io::Write};
use std::time::Instant;

use fast_convolver::{FastConvolver, ConvolutionMode};

mod ring_buffer;
mod fast_convolver;

fn show_info() {
    eprintln!("MUSI-6106 Assignment Executable");
    eprintln!("(c) 2024 Stephen Garrett & Ian Clester");
}

fn main() {
   show_info();

    // Parse command line arguments
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 1 {
        eprintln!("Usage: {} <input wav> <impulse response wav>", args[0]);
        return
    }

    // Open the input wave file - Assume mono audio sources
    let input_wav = hound::WavReader::open(&args[1]).unwrap();
    let impulse_response_wav = hound::WavReader::open(&args[2]).unwrap();

    // Wave to samples
    let input_samp: Vec<f32> = input_wav.into_samples::<i16>().map(|s| s.unwrap() as f32).collect();
    let impulse_response_samp: Vec<f32> = impulse_response_wav.into_samples::<i16>().map(|s| s.unwrap() as f32).collect();

    // Define block_size and combined length
    let block_size = 1024;

    let ir_len = impulse_response_samp.len();
    let input_len = input_samp.len();
    let combined_len = input_len + ir_len - 1;

    // Compare runtime performance
    let mut t_convolver = FastConvolver::new(&impulse_response_samp, ConvolutionMode::TimeDomain);
    let mut output_t = vec![0.0; combined_len];

    let mut f_convolver = FastConvolver::new(&impulse_response_samp, ConvolutionMode::FrequencyDomain { block_size });
    let mut output_f = vec![0.0; combined_len];
   
   // Time domain convolution
    let t_c_start = Instant::now();
    t_convolver.process(&input_samp, &mut output_t);
    let t_c_elapsed = t_c_start.elapsed();

    // Frequency domain convolution
    let f_c_start = Instant::now();
    f_convolver.process(&input_samp, &mut output_f);
    let f_c_elapsed = f_c_start.elapsed();

    // Comparison output in markdown file
    let mut markdown_file = File::create("fast_convolver_performance_results.md").unwrap();
    writeln!(markdown_file, "# Time vs. Frequency Domain Convolution Performance Comparison").unwrap();
    writeln!(markdown_file, "## Time Domain").unwrap();
    writeln!(markdown_file, "Elapsed Time: {:?}", t_c_elapsed).unwrap();
    writeln!(markdown_file, "## Frequency Domain").unwrap();
    writeln!(markdown_file, "Elapsed Time: {:?}", f_c_elapsed).unwrap();
}