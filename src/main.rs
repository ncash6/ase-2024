use std::{fs::File, io::Write};
use hound::{WavWriter, WavSpec};
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
    if args.len() < 5 {
        eprintln!("Usage: {} <input wav> <impulse response wav> [time|frequency] <output wav>", args[0]);
        return
    }

    // Open the input wave file - Assume mono audio sources
    let input_wav = hound::WavReader::open(&args[1]).unwrap();
    let spec = input_wav.spec();
    let impulse_response_wav = hound::WavReader::open(&args[2]).unwrap();

    // Open output wave file
    let mut writer = hound::WavWriter::create(&args[4], spec).unwrap();

    // Define block_size and combined length
    let block_size = 1024;

    // Command line convolution mode argument (Time or Frequency)
    let convolution_mode = match args[3].as_str() {
        "time" => ConvolutionMode::TimeDomain,
        "frequency" => ConvolutionMode::FrequencyDomain { block_size },
        _ => {
            eprintln!("Invalid convolution mode. Acceptable options are time or frequency.");
            return;
        }
    };

    // Wave to samples
    let input_samp: Vec<f32> = input_wav.into_samples::<i16>().map(|s| s.unwrap() as f32).collect();
    let impulse_response_samp: Vec<f32> = impulse_response_wav.into_samples::<i16>().map(|s| s.unwrap() as f32).collect();

    let ir_len = impulse_response_samp.len();
    let input_len = input_samp.len();
    let combined_len = input_len + ir_len - 1;


    // Process audio with fast convolver
    let mut out_convolve = FastConvolver::new(&impulse_response_samp, convolution_mode);
    let mut out_buffer = vec![0.0; combined_len];
    out_convolve.process(&input_samp, &mut out_buffer);

    // Write processed samples to wave file
    for ob in &out_buffer {
        let out_ob_i16 = (*ob * (1 << 15) as f32) as i16;
        writer.write_sample(out_ob_i16).unwrap();

    }
    
    
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