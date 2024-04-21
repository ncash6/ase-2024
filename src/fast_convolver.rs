use crate::ring_buffer::{self, RingBuffer};
use rustfft::FftPlanner;
use rustfft::num_complex::Complex;
use rustfft::num_traits::Zero;

pub struct FastConvolver<'a> {
    n_ir: usize,
    impulse_response: &'a [f32],
    mode: ConvolutionMode,
    ring_buffer: RingBuffer<f32>, 
}


#[derive(Debug, Clone, Copy)]
pub enum ConvolutionMode {
    TimeDomain,
    FrequencyDomain { block_size: usize },
}

impl<'a> FastConvolver<'a> {
    pub fn new(impulse_response: &'a[f32], mode: ConvolutionMode) -> Self {
        // Initialized variables for self
        let n_ir = impulse_response.len() as usize; // Lengths of Impulse
        let ring_buffer = RingBuffer::new(n_ir); 

        FastConvolver {
            impulse_response,
            mode,
            n_ir,
            ring_buffer,
        }
    }

    pub fn reset(&mut self) {
    self.ring_buffer.reset(); 
    }

    pub fn process(&mut self, input: &[f32], output: &mut [f32]) {
        match self.mode {
            ConvolutionMode::TimeDomain => {
                // Simple Time Domain Convolution FIR filter

                let n_input = input.len();
                let combined_len = n_input + self.n_ir - 1; // Common Length between signal and impulse

                // Initialized impulse response ring buffer
                let mut ring_buffer = RingBuffer::new(combined_len);

                // Push Impulse Response into RingBuffer
                for &val in self.impulse_response.iter() {
                    ring_buffer.push(val);
                }

                // Zero pad input to combined length
                let mut signal = input.to_vec();
                signal.resize(combined_len, 0.0);

                // Iterate over each sample for convolution
                for k in 0..combined_len {
                    let mut current_sum = 0.0;

                    // Time domain convolution
                    for i in std::cmp::max(0, k + 1 - self.n_ir)..std::cmp::min(k + 1, n_input) {
                        current_sum += signal[i] * ring_buffer.get(combined_len - i - 1);
                    }

                    // Convolution results stored in output buffer
                    output[k] = current_sum 
                }
            }
            ConvolutionMode::FrequencyDomain { block_size } => {
                // Partitioned Fast Convolution Frequency Domain
                
                let n_input = input.len();
                let combined_len = n_input + self.n_ir - 1; // Common Length between signal and impulse

                // Initialized impulse response ring buffer
                let mut ring_buffer = RingBuffer::new(combined_len);

                // Push Impulse Response into RingBuffer
                for &val in self.impulse_response.iter() {
                    ring_buffer.push(val);
                }

                // Initialized output index
                let mut index_output = 0;

                for start_blk in (0..n_input).step_by(block_size) {
                    let end_blk = std::cmp::min(start_blk + self.n_ir- 1, n_input);
                    let blk_len= end_blk - start_blk;

                    // Zero pad input block to combined length
                    let mut signal_blk = vec![Complex::zero(); combined_len];
                    for (i, &val) in input[start_blk..end_blk].iter().enumerate() {
                        if let Some(signal_val) = signal_blk.get_mut(i) {
                            *signal_val = Complex::new(val, 0.0);
                        }  
                    }

                    // FFT on input signal
                    let mut fft_plan = FftPlanner::new();
                    let fft = fft_plan.plan_fft(combined_len, rustfft::FftDirection::Forward);
                    fft.process(&mut signal_blk);

                    // Zero pad impulse reponse block to combined length
                    let mut ir_blk = vec![Complex::zero(); combined_len];
                    for i in 0..combined_len {
                        ir_blk[i] = Complex::new(ring_buffer.get(i), 0.0);
                    }

                    // FFT on impulse response
                    let fft = fft_plan.plan_fft(combined_len, rustfft::FftDirection::Forward);
                    fft.process(&mut ir_blk);

                    // Element-wise multiplication (convolution) in frequency domain
                    for i in 0..combined_len {
                        if let Some(signal_val) = signal_blk.get_mut(i) {
                            if let Some(ir_val) = ir_blk.get(i) {
                                *signal_val *= ir_val;
                            }
                        }
                    }

                    // Return to time domain
                    let mut ifft_plan = FftPlanner::new();
                    let ifft = ifft_plan.plan_fft(combined_len, rustfft::FftDirection::Inverse);
                    ifft.process(&mut signal_blk);

                    // Scaling with FFT length
                    let scale = 1.0 / (combined_len as f32);
                    for val in &mut signal_blk {
                        *val *= Complex::new(scale, 0.0);
                    }

                    // Convolution results stored in output buffer
                    for k in 0..blk_len {
                        if let Some(output_val) = output.get_mut(index_output + k) {
                            *output_val = signal_blk[k].re;
                        }
                    }

                    // Shift ring buffer header by block length
                    for _ in 0..blk_len {
                        ring_buffer.pop();
                    }
                    for &val in &input[start_blk..end_blk] {
                        ring_buffer.push(val);
                    }

                    // Updated ouput index 
                    index_output += blk_len;
                }
            }
        }
    }

        /// Flushes the remaining reverb tail after processing all input samples.
        /// Uses the internal ring buffer to output remaining data.
        pub fn flush(&mut self, output: &mut [f32]) {
            match self.mode {
                ConvolutionMode::TimeDomain => {
                    let remaining_samples = self.ring_buffer.len();
                    for i in 0..remaining_samples {
                        output[i] = self.ring_buffer.pop();
                    }
                    self.ring_buffer.reset();
                },
                ConvolutionMode::FrequencyDomain { block_size } => {
                    // Handle frequency domain flushing
                    let remaining_samples = self.ring_buffer.len();
                    let mut temp_buffer = vec![Complex::zero(); block_size];
                    for i in 0..remaining_samples {
                        temp_buffer[i] = Complex::new(self.ring_buffer.pop(), 0.0);
                    }
                    
                    let mut fft_plan = FftPlanner::new();
                    let ifft = fft_plan.plan_fft(block_size, rustfft::FftDirection::Inverse);
                    ifft.process(&mut temp_buffer);
    
                    let scale = 1.0 / block_size as f32;
                    for (i, val) in temp_buffer.iter().enumerate().take(remaining_samples) {
                        output[i] = val.re * scale;
                    }
                    self.ring_buffer.reset();
                }
            }
        }
    }
    
    #[cfg(test)]
    mod tests {
        use super::*;
        use rand::Rng;
        
        
        /// Generates a random impulse response of the given length.
        fn random_impulse_response(length: usize) -> Vec<f32> {
            let mut rng = rand::thread_rng();
            (0..length).map(|_| rng.gen::<f32>()).collect()
        }
    
        // Identity Test for Time Domain
        #[test]
        fn test_identity_time_domain() {
            let ir = random_impulse_response(51);
            let mut input = vec![0.0; 10];
            input[3] = 1.0; // Impulse at index 3
            let mut output = vec![0.0; 60]; // Output size must be large enough to capture the full response
    
            let mut convolver = FastConvolver::new(&ir, ConvolutionMode::TimeDomain);
            convolver.process(&input, &mut output);
    
            // Check that the impulse response is correctly placed in the output
            for i in 0..51 {
                assert_eq!(output[i + 3], ir[i], "Output mismatch at sample {}", i + 3);
            }
        }
    
        // Flush Test for Time Domain
        #[test]
        fn test_flush_time_domain() {
            let ir = random_impulse_response(51);
            let mut input = vec![0.0; 10];
            input[3] = 1.0;
            let mut output = vec![0.0; 60];
    
            let mut convolver = FastConvolver::new(&ir, ConvolutionMode::TimeDomain);
            convolver.process(&input, &mut output);
            let mut tail = vec![0.0; 41]; // The expected length of the tail is the length of the IR minus one
            convolver.flush(&mut tail);
    
            // Validate the tail of the impulse response
            for (i, &val) in tail.iter().enumerate() {
                assert_eq!(val, ir[i + 10], "Flush output does not match at index {}", i);
            }
        }
    
        // Identity Test for Frequency Domain
        #[test]
        fn test_identity_frequency_domain() {
            let ir = random_impulse_response(51);
            let mut input = vec![0.0; 10];
            input[3] = 1.0; // Impulse at index 3
            let mut output = vec![0.0; 60];
    
            let mut convolver = FastConvolver::new(&ir, ConvolutionMode::FrequencyDomain { block_size: 64 });
            convolver.process(&input, &mut output);
    
            // Check that the impulse response is correctly placed in the output
            for i in 0..51 {
                assert_eq!(output[i + 3], ir[i], "Output mismatch at sample {}", i + 3);
            }
        }
    
        // Flush Test for Frequency Domain
        #[test]
        fn test_flush_frequency_domain() {
            let ir = random_impulse_response(51);
            let mut input = vec![0.0; 10];
            input[3] = 1.0;
            let mut output = vec![0.0; 60];
    
            let mut convolver = FastConvolver::new(&ir, ConvolutionMode::FrequencyDomain { block_size: 64 });
            convolver.process(&input, &mut output);
            let mut tail = vec![0.0; 41];
            convolver.flush(&mut tail);
    
            // Validate the tail of the impulse response
            for (i, &val) in tail.iter().enumerate() {
                assert_eq!(val, ir[i + 10], "Flush output does not match at index {}", i);
            }
        }
    
        // Block Size Test for Frequency Domain
        #[test]

        fn test_blocksize_frequency_domain_simplified() {
            let ir = random_impulse_response(51);
            let input = vec![0.0; 10000]; // Zero input
            let block_sizes = vec![1, 13, 1023, 2048, 1, 17, 5000, 1897];
    
            for &block_size in &block_sizes {
                let mut output = vec![0.0; 10050]; // Output size must be large enough to capture the full response
                let mut convolver = FastConvolver::new(&ir, ConvolutionMode::FrequencyDomain { block_size });
                convolver.process(&input, &mut output);
    
                // Check that the output is zero for zero input
                for i in 0..10050 {
                    assert_eq!(output[i], 0.0, "Output should be zero for zero input with block size {}", block_size);
                }
            }
        }
    }
    

