// The user-interface for this vibrato.rs file is modelled after the comb_filter.rs
// written by Stephen Garrett & Ian Clester for MUSI 6106 Assignment 1.

use crate::ring_buffer::RingBuffer;
use crate::lfo::LFO;

pub struct Vibrato {
    lfo: LFO, 
    delay_buffer: RingBuffer<f32>, // size of delay buffer
    // delay_samples: usize, // depth of vibrato in samples (for lfo)
    // sample_rate_hz: f32, // sample rate in hz
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
            // delay_samples,
            // sample_rate_hz,
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

    // RingBuffer tests as specified in excercise 2 in ring_buffer.rs
    // Test 1 (from exercise 2): Call various parameter values and print out the results vs expected values.
    #[test]
    fn test_parameters() {
    }

    // Test 2 (from exercise 2): Test all implemented functions
    #[test]
    fn test_vibrato_api() {
    }

    // Test 3 (from exercise 2): Create a simple test signal (e.g., a unit impulse or a ramp) and process with Vibrato using a for loop. Does the result look as expected?
    #[test]
    fn test_vibrato_process() {
    }

    // Test 4: Output equals delayed input when modulation amplitude is 0.
    #[test]
    fn test_zero_amplitude() {
    }

    // Test 5: DC input results in DC output, regardless of parameters.
    #[test]
    fn test_dc_input() {
    }

    // Test 6: Varying input block size.
    #[test]
    fn test_vary_blocksize() {
    }
    
    // Test 7: Zero input signal.
    #[test]
    fn test_silence() {
    }
    
    // Test 8: [Additional test] Ensure ringbuffers reset upon command
    #[test]
    fn test_reset() {
    }
}
