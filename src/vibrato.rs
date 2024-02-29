// The user-interface for this vibrato.rs file is modelled after the comb_filter.rs
// written by Stephen Garrett & Ian Clester for MUSI 6106 Assignment 1.

use crate::ring_buffer::RingBuffer;
use crate::lfo::LFO;

pub struct Vibrato {
    lfo: LFO, 
    delay_buffer: RingBuffer<f32>, // size of delay buffer
    delay_samples: usize, // depth of vibrato in samples (for lfo)
    sample_rate_hz: f32, // sample rate in hz
}

impl Vibrato {
    pub fn new(frequency_hz: f32, amplitude: f32, sample_rate_hz: f32, delay_time_sec: f32, capacity: usize, num_channels: usize) -> Self {
        // Initialized variables for self
        let lfo = LFO::new(sample_rate_hz, frequency_hz, amplitude, capacity);
        let delay_samples = (sample_rate_hz * delay_time_sec) as usize;
        let delay_buffer = RingBuffer::new(capacity);

        Vibrato {
            lfo,
            delay_buffer,
            delay_samples,
            sample_rate_hz,
        }
    }

    pub fn reset(&mut self) {
        // Reset variables subject to processing
        self.delay_buffer.reset();
        self.lfo.reset();

    }

    pub fn process(&mut self, input: &[&[f32]], output: &mut [&mut [f32]]) {

        for (channel, &in_channel) in input.iter().enumerate() {
            for (samp, &in_samp) in in_channel.iter().enumerate() {
                // Vibrato modulation depth from LFO
                let mod_depth = self.lfo.process();

                // Delayed sample from delay buffer
                let delay_index_sample = self.delay_buffer.get_frac(mod_depth);

                // Vibrato effect applied to original
                output[channel][samp] = in_samp + delay_index_sample;

                // Updated delay buffer with effected sample
                self.delay_buffer.push(in_samp);
            }
        }
    }
}

// Vibrato Unit Tests
#[cfg(test)]
mod tests {
    use super::*;

    // Floating point tolerance for equality assert! comparisons
    const EPSILON: f32 = 1e-6;

    // RingBuffer tests as specified in excercise 2 in ring_buffer.rs
    // Test 1 (from exercise 2): Call various parameter values and print out the results vs expected values
    #[test]
    fn test_parameters() {
        // Parameters
        let freq_hz = 2.0;
        let amp = 0.5;
        let fs = 44100.0;
        let delay_sec = 0.01;
        let capacity = 1024;
        let channels = 2;

        // Initialize vibrato effect
        let mut vibrato = Vibrato::new(freq_hz,amp,fs, delay_sec, capacity, channels);

        // Test input signal
        let input: Vec<Vec<f32>> = vec![
            vec![0.1, 0.2, 0.3],
            vec![0.4, 0.5, 0.6],
        ];

        // Expected output signal 
        let expected_output: Vec<Vec<f32>> = vec![
            vec![0.1, 0.2, 0.3],
            vec![0.4 + 0.1 * amp, 0.5 + 0.2 * amp, 0.6 + 0.3 * amp],
        ];

        // Process input signal
        let mut output = vec![vec![0.0; input[0].len()]; channels];

        // Buffers to slices
        let in_slice: Vec<&[f32]> = input.iter().map(|channel| &channel[..]).collect();
        let mut out_slice: Vec<&mut [f32]> = output.iter_mut().map(|channel| &mut channel[..]).collect();

        vibrato.process(&in_slice, &mut out_slice);

        // Compare the output with expected values
        for (channel, expected_channel) in output.iter().zip(expected_output.iter()) {
            for (sample, &expected_sample) in channel.iter().zip(expected_channel.iter()) {
                assert!(
                    (sample - expected_sample).abs() < EPSILON,
                    "Output value {:?} does not match expected value {:?}",
                    sample,
                    expected_sample,
                );
            }
        }
    }

    // Test 2 (from exercise 2): Test all implemented functions
    #[test]
    fn test_vibrato_api() {
        // Parameters
        let freq_hz = 2.0;
        let amp = 0.5;
        let fs = 44100.0;
        let delay_sec = 0.01;
        let capacity = 1024;
        let channels = 1;

        // Initialize vibrato effect
        let mut vibrato = Vibrato::new(freq_hz,amp,fs, delay_sec, capacity, channels);

        // Test input signal
        let input: Vec<Vec<f32>> = vec![
            vec![0.1, 0.2, 0.3],
        ];

        // Expected output signal 
        let expected_output: Vec<Vec<f32>> = vec![
            vec![0.1, 0.2, 0.3 + 0.2 * amp],
        ];

        // Process input signal
        let mut output = vec![vec![0.0; input[0].len()]; channels];

        let in_slice: Vec<&[f32]> = input.iter().map(|channel| &channel[..]).collect();
        let mut out_slice: Vec<&mut [f32]> = output.iter_mut().map(|channel| &mut channel[..]).collect();

        vibrato.process(&in_slice, &mut out_slice);

        // Reset vibrato processor
        vibrato.reset();

        // Output is expected to be zero after reset
        let mut reset_output = vec![vec![0.0; input[0].len()]; channels];
        let mut reset_slice: Vec<&mut [f32]> = reset_output.iter_mut().map(|channel| &mut channel[..]).collect();

        vibrato.process(&in_slice, &mut reset_slice);

        for channel in reset_output {
            for sample in channel {
                assert!(
                    sample.abs() < EPSILON,
                    "Output value after reset is not zero: {} with epsilon {}",
                    sample,
                    EPSILON
                );
            }
        }


    }

    // Test 3 (from exercise 2): Create a simple test signal (e.g., a unit impulse or a ramp) and process with Vibrato using a for loop. Does the result look as expected?
    #[test]
    fn test_vibrato_process() {
        // Parameters
        let freq_hz = 2.0;
        let amp = 0.5;
        let fs = 44100.0;
        let delay_sec = 0.01;
        let capacity = 1024;
        let channels = 1;

        // Initialize vibrato effect
        let mut vibrato = Vibrato::new(freq_hz,amp,fs, delay_sec, capacity, channels);

        // Test impulse signal
        let impulse_len = 100;
        let mut input = vec![vec![0.0; impulse_len]; channels];
        input[0][0] = 1.0;

        // Process input signal
        let mut output = vec![vec![0.0; impulse_len]; channels];

        let in_slice: Vec<&[f32]> = input.iter().map(|channel| &channel[..]).collect();
        let mut out_slice: Vec<&mut [f32]> = output.iter_mut().map(|channel| &mut channel[..]).collect();

        for _ in 0..impulse_len {
            vibrato.process(&in_slice[..], &mut out_slice[..]);
        }

        // We expect that vibrato effect will have amplitude applied at every iteration 
        for (output_channel, output_samples) in output.iter().enumerate() {
            for (sample_index, &sample) in output_samples.iter().enumerate() {
                let expected_value = if sample_index as f32 <= fs * delay_sec {
                    1.0 // Before the delay, the sample remains 1.0
                } else {
                    1.0 + amp // After the delay, the amplitude is added
                };

                assert!(
                    (sample - expected_value).abs() < EPSILON,
                    "Output value {} does not match expected value {} with epsilon {} at index {} of channel {}",
                    sample,
                    expected_value,
                    EPSILON,
                    sample_index,
                    output_channel
                );
            }
        }
    }

    // Test 4: Output equals delayed input when modulation amplitude is 0.
    #[test]
    fn test_zero_amplitude() {
         // Parameters
         let freq_hz = 2.0;
         let amp = 0.0;
         let fs = 44100.0;
         let delay_sec = 0.01;
         let capacity = 1024;
         let channels = 1;
 
         // Initialize vibrato effect
         let mut vibrato = Vibrato::new(freq_hz,amp,fs, delay_sec, capacity, channels);
 
         // Test impulse signal
         let impulse_len = 100;
         let mut input = vec![vec![0.0; impulse_len]; channels];
 
         // Process input signal
         let mut output = vec![vec![0.0; impulse_len]; channels];
 
         let in_slice: Vec<&[f32]> = input.iter().map(|channel| &channel[..]).collect();
         let mut out_slice: Vec<&mut [f32]> = output.iter_mut().map(|channel| &mut channel[..]).collect();
 
         for _ in 0..impulse_len {
             vibrato.process(&in_slice[..], &mut out_slice[..]);
         }
 
         // We expect that processed output will equal delayed input
         for (output_channel, output_samples) in output.iter().enumerate() {
             for (sample_index, &sample) in output_samples.iter().enumerate() {
                 let expected_value = if sample_index as f32 <= fs * delay_sec {
                     input[output_channel][sample_index] // Before the delay, output equals input
                 } else {
                     0.0 // After the delay, input is delayed and becomes 0
                 };
 
                 assert!(
                     (sample - expected_value).abs() < EPSILON,
                     "Output value {} does not match expected value {} with epsilon {} at index {} of channel {}",
                     sample,
                     expected_value,
                     EPSILON,
                     sample_index,
                     output_channel
                 );
             }
         }
    }

    // Test 5: DC input results in DC output, regardless of parameters.
    #[test]
    fn test_dc_input() {
        // Parameters
        let freq_hz = 2.0;
        let amp = 0.5;
        let fs = 44100.0;
        let delay_sec = 0.01;
        let capacity = 1024;
        let channels = 1;

        // Initialize vibrato effect
        let mut vibrato = Vibrato::new(freq_hz,amp,fs, delay_sec, capacity, channels);

        // Test dc input signal
        let impulse_len = 100;
        let input = vec![vec![1.0; impulse_len]; channels];

        // Process input signal
        let mut output = vec![vec![0.0; impulse_len]; channels];

        let in_slice: Vec<&[f32]> = input.iter().map(|channel| &channel[..]).collect();
        let mut out_slice: Vec<&mut [f32]> = output.iter_mut().map(|channel| &mut channel[..]).collect();

        vibrato.process(&in_slice[..], &mut out_slice[..]);

        // We expect that output signal is DC same as input signal
        for (output_channel, output_samples) in output.iter().enumerate() {
            for (sample_index, &sample) in output_samples.iter().enumerate() {
                assert!(
                    (sample - 1.0).abs() < EPSILON,
                    "Output value {} does not match expected DC value 1.0 with epsilon {} at index {} of channel {}",
                    sample,
                    EPSILON,
                    sample_index,
                    output_channel
                );
            }
        }
    }

    // Test 6: Varying input block size.
    #[test]
    fn test_vary_blocksize() {
         // Parameters
         let freq_hz = 2.0;
         let amp = 0.5;
         let fs = 44100.0;
         let delay_sec = 0.01;
         let capacity = 512;
         let channels = 1;
 
         // Initialize vibrato effect
         let mut vibrato = Vibrato::new(freq_hz,amp,fs, delay_sec, capacity, channels);
 
         // Test signal with length longer than block size
         let signal_len = capacity * 2;
         let mut input = vec![vec![0.0; signal_len]; channels];
 
         // Process input signal
         let mut output = vec![vec![0.0; signal_len]; channels];
 
         let in_slice: Vec<&[f32]> = input.iter().map(|channel| &channel[..]).collect();
         let mut out_slice: Vec<&mut [f32]> = output.iter_mut().map(|channel| &mut channel[..]).collect();
 
         for _ in 0..signal_len {
             vibrato.process(&in_slice[..], &mut out_slice[..]);
         }
 
         // Check whether output matches expected output
         for (output_channel, output_samples) in output.iter().enumerate() {
             for (sample_index, &sample) in output_samples.iter().enumerate() {
                 let expected_value = input[output_channel][sample_index]; 
                 assert!(
                     (sample - expected_value).abs() < EPSILON,
                     "Output value {} does not match expected value {} with epsilon {} at index {} of channel {}",
                     sample,
                     expected_value,
                     EPSILON,
                     sample_index,
                     output_channel
                 );
             }
         }
    }
    
    // Test 7: Zero input signal.
    #[test]
    fn test_silence() {
        // Parameters
        let freq_hz = 2.0;
        let amp = 0.5;
        let fs = 44100.0;
        let delay_sec = 0.01;
        let capacity = 1024;
        let channels = 1;

        // Initialize vibrato effect
        let mut vibrato = Vibrato::new(freq_hz,amp,fs, delay_sec, capacity, channels);

        // Test zero input signal
        let signal_len = 100;
        let mut input = vec![vec![0.0; signal_len]; channels];

        // Process input signal
        let mut output = vec![vec![0.0; signal_len]; channels];

        let in_slice: Vec<&[f32]> = input.iter().map(|channel| &channel[..]).collect();
        let mut out_slice: Vec<&mut [f32]> = output.iter_mut().map(|channel| &mut channel[..]).collect();

        for _ in 0..signal_len {
            vibrato.process(&in_slice[..], &mut out_slice[..]);
        }

        // We expect that output signal will be zero same as input signal
        for (output_channel, output_samples) in output.iter().enumerate() {
            for (sample_index, &sample) in output_samples.iter().enumerate() {
                assert!(
                    sample.abs() < EPSILON,
                    "Output value {} is not zero with epsilon {} at index {} of channel {}",
                    sample,
                    EPSILON,
                    sample_index,
                    output_channel
                );
            }
        }
    }
    
    // Test 8: [Additional test] Ensures vibato effect applies to input signal correctly
    #[test]
    fn test_accuracy() {
        // Parameters
        let freq_hz = 5.0;
        let amp = 0.2;
        let fs = 44100.0;
        let delay_sec = 0.01;
        let capacity = 1024;
        let channels = 1;

        // Initialize vibrato effect
        let mut vibrato = Vibrato::new(freq_hz,amp,fs, delay_sec, capacity, channels);

        // Test input sine wave at specified frequency
        let signal_len = 1000;
        let mut input: Vec<Vec<f32>> = Vec::with_capacity(channels);
        for _ in 0..channels {
            let channel_samples: Vec<f32> = (0..signal_len)
                .map(|i| ((2.0 * std::f32::consts::PI * freq_hz * i as f32 / fs).sin() * amp))
                .collect();
            input.push(channel_samples);
        }

        // Process input signal
        let mut output = vec![vec![0.0; signal_len]; channels];

        let in_slice: Vec<&[f32]> = input.iter().map(|channel| &channel[..]).collect();
        let mut out_slice: Vec<&mut [f32]> = output.iter_mut().map(|channel| &mut channel[..]).collect();

        vibrato.process(&in_slice[..], &mut out_slice[..]);

        // We expect that output signal will show variance in amplitude around input frequency indicating modulation (greater than zero) 
        let mean_amp = output[0].iter().sum::<f32>() / signal_len as f32;

        assert!(
            mean_amp.abs() > EPSILON,
            "Mean amplitude of output signal is too close to zero, indicating no modulation"
        );
    }
}
