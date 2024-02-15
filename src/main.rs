use std::{fs::File, io::Write};

mod comb_filter;

fn show_info() {
    eprintln!("MUSI-6106 Assignment Executable");
    eprintln!("(c) 2024 Stephen Garrett & Ian Clester");
}

fn main() {
   show_info();

    // Parse command line arguments
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 5 {
        eprintln!("Usage: {} <input wave filename> <output text filename>", args[0]);
        return
    }

    // Open the input wave file
    let mut reader = hound::WavReader::open(&args[1]).unwrap();
    let spec = reader.spec();
    let channels = spec.channels; // number of channels in input file

    // TODO: Modify this to process audio in blocks using your comb filter and write the result to an audio file.

    // Define block size
    let block_size = 1024;

    // Open output wave file
    let mut writer = hound::WavWriter::create(&args[2], spec).unwrap();

    // Assign sample rate from input
    let sample_rate = spec.sample_rate as f32; 

    // Command line comb filter argument (FIR or IIR)
    let comb_filter_type = match args[3].as_str() { 
        "FIR" => comb_filter::FilterType::FIR,
        "IIR" => comb_filter::FilterType::IIR,
        _ => {
            eprintln!("Invalid filter type. Acceptable options as FIR or IIR.");
            return;
        }
    };

    // Command line filter parameter type (Gain or Delay)
    let comb_filter_parameter = match args[4].as_str() {
        "Gain" => comb_filter::FilterParam::Gain,
        "Delay" => comb_filter::FilterParam::Delay,
        _ => {
            eprintln!("Invalid filter type. Acceptable options as Gain or Delay.");
            return;
        }
    }; 

    // Command line filter parameter argument 
    let comb_filter_value: f32 = match args[5].parse() {
        Ok(value) => value,
        Err(_) => {
            eprintln!("Invalid effect parameter value. Floating point number required.");
            return;
        }
    }; 

    // Set max_delay_sec
    let mut max_delay = 0.0;

    if let comb_filter::FilterParam::Delay = comb_filter_parameter {
        max_delay = comb_filter_value;
    }

    // Instantiate CombFilter
    let mut comb_filter = comb_filter::CombFilter::new(comb_filter_type, max_delay, sample_rate, channels as usize);
    

    // Store wave file audio samples in buffer
    let mut in_buffer =  vec![0.0; block_size * channels as usize];
    let mut out_buffer = vec![0.0; block_size * channels as usize];
    
    while let Some(samp_read) = reader.samples::<i16>().next() {
        let samp: i16 = match samp_read {
            Ok(samp) => samp,
            Err(err) => {
                eprintln!("Error reading samples: {}", err);
                return;
            }   
        };

        // Normalization
        for n in 0..block_size {
            in_buffer[n] = samp as f32 / std::i16::MAX as f32;
        }

        // Block audio and process using comb filter
        for chnl in 0..channels {
            let start_indx = chnl as usize * block_size;
            let end_indx = (chnl as usize + 1) * block_size;
            let in_chnl = &in_buffer[start_indx..end_indx];
            let out_chnl = &mut out_buffer[start_indx..end_indx];
            comb_filter.process(&[in_chnl], &mut[out_chnl]);
        }

        // Write processed block to wave file
        for ob in &out_buffer {
            let out_samp = (*ob * std::i16::MAX as f32) as i16; // Denormalization
            if let Err(err)= writer.write_sample(out_samp) {
                eprintln!("Error wrinting sample to wave file: {}", err);
                return;
            }
        }
    }

    // Read audio data and write it to the output text file (one column per channel)
    let mut out = File::create(&args[2]).expect("Unable to create file");
    for (i, sample) in reader.samples::<i16>().enumerate() {
        let sample = sample.unwrap() as f32 / (1 << 15) as f32;
        write!(out, "{}{}", sample, if i % channels as usize == (channels - 1).into() { "\n" } else { " " }).unwrap();
    }
}
