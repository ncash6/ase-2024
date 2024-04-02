use std::{fs::File, io::Write};
use hound;

fn show_info() {
    eprintln!("MUSI-6106 Assignment Executable");
    eprintln!("(c) 2024 Stephen Garrett & Ian Clester");
}

fn main() -> std::io::Result<()> {
   show_info();

    // Parse command line arguments
    // First argument is input .wav file, second argument is output text file.
    let args: Vec<String> = std::env::args().collect();
    
    // TODO: your code here
    
    // TODO: your code here; see `hound::WavReader::open`.
    let mut wav_audio = hound::WavReader::open(&args[1]).unwrap();
    let mut text_out = File::create(&args[2]).unwrap();

    // Open the input wave file and determine number of channels
    let num_channel = wav_audio.spec();
    
    // Read audio data and write it to the output text file (one column per channel)
    let mut stereo_channel = 0;
    
    // TODO: your code here; we suggest using `hound::WavReader::samples`, `File::create`, and `write!`.
    for pnt in wav_audio.samples::<i16>() {
        //       Remember to convert the samples to floating point values and respect the number of channels!
        let float_pnt = pnt.unwrap() as f32;
        if stereo_channel == num_channel.channels {
            writeln!(&mut text_out)?;
            stereo_channel = 0
        }
        write!(&mut text_out, "{} ", float_pnt/65565.0)?; // Normalization
        stereo_channel += 1;
    }
    Ok(()) 
}

//comment
