use std::{fs::File, io::Write};

use hound::{WavWriter, WavSpec};

mod ring_buffer;
mod vibrato;
mod lfo;

fn show_info() {
    eprintln!("MUSI-6106 Assignment Executable");
    eprintln!("(c) 2024 Stephen Garrett & Ian Clester");
}

fn main() {
   show_info();

    // Parse command line arguments
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("Usage: {} <input wave filename> <output text filename>", args[0]);
        return
    }

    // Open the input wave file
    let mut reader = hound::WavReader::open(&args[1]).unwrap();
    let spec = reader.spec();
    let num_channels = spec.channels as usize;
    let sample_rate = spec.sample_rate;

    // Initialize vibrato implementation
    let vib_frequency_hz = 5.0;
    let vib_amplitude = 0.1;
    let vib_delay_sec = 0.01;
    let buffer_capacity = 1024;
    let mut vibrato = vibrato::Vibrato::new(
        vib_frequency_hz,
        vib_amplitude,
        sample_rate as f32,
        vib_delay_sec,
        buffer_capacity,
        num_channels,

    );

    // Read audio data and write it to the output text file (one column per channel)
    let spec = WavSpec {
        channels: num_channels as u16,
        sample_rate,
        bits_per_sample: 16,
        sample_format: hound::SampleFormat::Int,
    };

    let mut out = WavWriter::create(&args[2], spec).unwrap();

    // Write vibrato processed audio to output text file (one column per channel)
    for sample in reader.samples::<i32>() {
        let sample = sample.unwrap() as f32 / (1 << 15) as f32;

        // Initialize input and output buffers
        let mut in_buffer = vec![vec![0.0; 1]; num_channels];
        let mut out_buffer = vec![vec![0.0; 1]; num_channels];

        // Apply vibrato effect to input audio
        in_buffer[0][0] = sample;

        // Buffers to slices
        let in_slice: Vec<&[f32]> = in_buffer.iter().map(|channel| &channel[..]).collect();
        let mut out_slice: Vec<&mut [f32]> = out_buffer.iter_mut().map(|channel| &mut channel[..]).collect();

        vibrato.process(&in_slice, &mut out_slice);

        // Write processed audio to text file
        for ch in &out_slice {
            for  &out_sample in ch.iter() {
                let out_samp_i16 = (out_sample * (1 << 15) as f32) as i16;
                out.write_sample(out_samp_i16).unwrap(); 
            }
        }
    }
}
