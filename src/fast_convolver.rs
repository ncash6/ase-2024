use crate::ring_buffer::{RingBuffer};
use rustfft::{FftPlanner, num_complex::Complex, num_traits::Zero};

pub struct FastConvolver<'a> {
    n_ir: usize,
    impulse_response: &'a [f32],
    mode: ConvolutionMode,
    // Adding ring buffer as part of the state to manage internal state reset
    ring_buffer: Option<RingBuffer<f32>>,
}

#[derive(Debug, Clone, Copy)]
pub enum ConvolutionMode {
    TimeDomain,
    FrequencyDomain { block_size: usize },
}

impl<'a> FastConvolver<'a> {
    pub fn new(impulse_response: &'a [f32], mode: ConvolutionMode) -> Self {
        let n_ir = impulse_response.len(); // Length of Impulse Response

        FastConvolver {
            impulse_response,
            mode,
            n_ir,
            ring_buffer: None, // Initialize without an active ring buffer
        }
    }

    pub fn reset(&mut self) {
        // Resets the internal state of the convolver
        if let Some(rb) = &mut self.ring_buffer {
            rb.reset(); // Resetting the existing ring buffer if present
        } else {
            // If no ring buffer is present, initialize it based on the mode
            match self.mode {
                ConvolutionMode::TimeDomain => {
                    self.ring_buffer = Some(RingBuffer::new(self.n_ir));
                },
                ConvolutionMode::FrequencyDomain { block_size } => {
                    self.ring_buffer = Some(RingBuffer::new(block_size));
                },
            }
        }
    }

    pub fn process(&mut self, input: &[f32], output: &mut [f32]) {
        match self.mode {
            ConvolutionMode::TimeDomain => {
                let n_input = input.len();

                // Initialize or get the existing ring buffer
                let ring_buffer = self.ring_buffer.get_or_insert_with(|| RingBuffer::new(self.n_ir));

                // Push Impulse Response into RingBuffer
                for &val in self.impulse_response.iter() {
                    ring_buffer.push(val);
                }

                // Iterate over each sample for convolution
                for i in 0..n_input {
                    let mut current_sum = 0.0;
                    for j in 0..self.n_ir {
                        if i >= j {
                            current_sum += input[i - j] * ring_buffer.get(j);
                        }
                    }
                    output[i] = current_sum;
                }
            },
            ConvolutionMode::FrequencyDomain { block_size } => {
                let n_input = input.len();
                let mut index_output = 0;

                // Initialize or get the existing ring buffer
                let ring_buffer = self.ring_buffer.get_or_insert_with(|| RingBuffer::new(block_size));

                for start_blk in (0..n_input).step_by(block_size) {
                    let end_blk = std::cmp::min(start_blk + block_size, n_input);
                    let blk_len = end_blk - start_blk;

                    let mut signal_blk = vec![Complex::<f32>::zero(); block_size];
                    for (i, &val) in input[start_blk..end_blk].iter().enumerate() {
                        signal_blk[i] = Complex::new(val, 0.0);
                    }

                    let mut fft_plan = FftPlanner::<f32>::new();
                    let fft = fft_plan.plan_fft(block_size, rustfft::FftDirection::Forward);
                    fft.process(&mut signal_blk);

                    let mut ir_blk = vec![Complex::<f32>::zero(); block_size];
                    for (i, &val) in self.impulse_response.iter().enumerate() {
                        ir_blk[i] = Complex::new(val, 0.0);
                    }

                    fft.process(&mut ir_blk);

                    for i in 0..block_size {
                        signal_blk[i] *= ir_blk[i];
                    }

                    let mut ifft_plan = FftPlanner::<f32>::new();
                    let ifft = ifft_plan.plan_fft(block_size, rustfft::FftDirection::Inverse);
                    ifft.process(&mut signal_blk);

                    let scale = 1.0 / block_size as f32;
                    for (i, val) in signal_blk.iter().enumerate().take(blk_len) {
                        output[index_output + i] = val.re * scale;
                    }

                    index_output += blk_len;
                }
            }
        }
    }

    pub fn max_flush_output_len(&self) -> usize {
        self.n_ir - 1
    }

   pub fn flush(&mut self, output: &mut [f32]) {
    match self.mode {
        ConvolutionMode::TimeDomain => {
            // Time domain flush logic
        },
        ConvolutionMode::FrequencyDomain { block_size } => {
            // Ensure there's a ring buffer and it's handled correctly
            if let Some(ring_buffer) = &self.ring_buffer {
                let remaining_samples = ring_buffer.len();
                if remaining_samples > 0 {
                    let mut temp_block = vec![Complex::<f32>::zero(); block_size];
                    for i in 0..remaining_samples {
                        temp_block[i] = Complex::new(ring_buffer.get(i), 0.0);
                    }
                    // Perform FFT, process, and IFFT as usual
                    let mut fft_plan = FftPlanner::<f32>::new();
                    let fft = fft_plan.plan_fft(block_size, rustfft::FftDirection::Forward);
                    fft.process(&mut temp_block);

                    let mut ir_temp = vec![Complex::<f32>::zero(); block_size];
                    for (i, &val) in self.impulse_response.iter().enumerate() {
                        ir_temp[i] = Complex::new(val, 0.0);
                    }
                    fft.process(&mut ir_temp);

                    for i in 0..block_size {
                        temp_block[i] *= ir_temp[i];
                    }

                    let ifft = fft_plan.plan_fft(block_size, rustfft::FftDirection::Inverse);
                    ifft.process(&mut temp_block);

                    let scale = 1.0 / block_size as f32;
                    for i in 0..remaining_samples {
                        output[i] = temp_block[i].re * scale;
                    }
                }
            }
        }
    }
}

}



#[cfg(test)]
mod tests {
    use super::*;
    use rand::{Rng, thread_rng};

    /// Helper function to generate a random impulse response of a given length.
    fn generate_random_impulse_response(length: usize) -> Vec<f32> {
        let mut rng = thread_rng();
        (0..length).map(|_| rng.gen::<f32>()).collect()
    }

    /// Helper function to generate an impulse signal with an impulse at a specified index within a given length.
    fn generate_impulse_signal(length: usize, impulse_position: usize) -> Vec<f32> {
        let mut signal = vec![0.0; length];
        if impulse_position < length {
            signal[impulse_position] = 1.0;
        }
        signal
    }

    #[test]
    fn identity_time_domain() {
        let ir = generate_random_impulse_response(51);
        let input_signal = generate_impulse_signal(10, 3);
        let mut convolver = FastConvolver::new(&ir, ConvolutionMode::TimeDomain);
        let mut output = vec![0.0; 10];

        convolver.process(&input_signal, &mut output);

        // Checking the output values for correctness
        for i in 0..10 {
            let expected_output = if i >= 3 { ir[i - 3] } else { 0.0 };
            assert!((output[i] - expected_output).abs() < 1e-5, "Mismatch at index {}: expected {}, got {}", i, expected_output, output[i]);
        }
    }

    #[test]
    fn identity_frequency_domain() {
        let ir = generate_random_impulse_response(51);
        let input_signal = generate_impulse_signal(10, 3);
        let mut convolver = FastConvolver::new(&ir, ConvolutionMode::FrequencyDomain { block_size: 64 });
        let mut output = vec![0.0; 10];

        convolver.process(&input_signal, &mut output);

        // Checking the output values for correctness
        for i in 0..10 {
            let expected_output = if i >= 3 { ir[i - 3] } else { 0.0 };
            assert!((output[i] - expected_output).abs() < 1e-5, "Mismatch at index {}: expected {}, got {}", i, expected_output, output[i]);
        }
    }

    #[test]
    fn flush_time_domain() {
        let ir = generate_random_impulse_response(51);
        let input_signal = generate_impulse_signal(10, 3);
        let mut convolver = FastConvolver::new(&ir, ConvolutionMode::TimeDomain);
        let mut output = vec![0.0; 10];
        let mut flush_output = vec![0.0; 51]; // Capture the full tail

        convolver.process(&input_signal, &mut output);
        convolver.flush(&mut flush_output);

        // The flush output should contain the remaining part of the impulse response
        let expected_flush = &ir[7..]; // Starts from 7 because the last 7 elements of IR should appear in flush
        assert_eq!(flush_output[..44], expected_flush[..44], "Flush output does not match expected tail.");
    }

    #[test]
    fn flush_frequency_domain() {
        let ir = generate_random_impulse_response(51);
        let input_signal = generate_impulse_signal(10, 3);
        let mut convolver = FastConvolver::new(&ir, ConvolutionMode::FrequencyDomain { block_size: 64 });
        let mut output = vec![0.0; 10];
        let mut flush_output = vec![0.0; 51]; // Capture the full tail

        convolver.process(&input_signal, &mut output);
        convolver.flush(&mut flush_output);

        // The flush output should contain the remaining part of the impulse response
        let expected_flush = &ir[7..]; // Starts from 7 because the last 7 elements of IR should appear in flush
        assert_eq!(flush_output[..44], expected_flush[..44], "Flush output does not match expected tail.");
    }

   #[test]
fn blocksize_variations() {
    let block_sizes = [1, 13, 1023, 2048, 1, 17, 5000, 1897];
    let ir = generate_random_impulse_response(51);
    let input_signal = vec![1.0; 10000];

    for &block_size in &block_sizes {
        let mut convolver = FastConvolver::new(&ir, ConvolutionMode::FrequencyDomain { block_size });
        let mut output = vec![0.0; 10000];

        convolver.process(&input_signal, &mut output);
        // Simple verification could involve checking a specific result or feature of the output
        // Example: checking the initial convolution result for correctness
        let expected_value = ir[0];  // Simplistic check for first convolution result
        for i in 0..51 {  // Checking only the first few outputs for simplicity
            assert!((output[i] - expected_value).abs() < 1e-5, "Mismatch in output for block size {}, at index {}", block_size, i);
        }
    }
}
}
