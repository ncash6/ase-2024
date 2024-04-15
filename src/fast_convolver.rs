use crate::ring_buffer::RingBuffer;
use rustfft::FftPlanner;
use rustfft::num_complex::Complex;
use rustfft::num_traits::Zero;

pub struct FastConvolver<'a> {
    n_ir: usize,
    impulse_response: &'a[f32],
    mode: ConvolutionMode,
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
        
        FastConvolver {
            impulse_response,
            mode,
            n_ir,
        }
    }

    pub fn reset(&mut self) {
        // Note: The way I have things setup now means none of the self variables are changed to be reset
        // Reset may be useful for the flush function.
        todo!("implement")
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
                        signal_blk[i] = Complex::new(val, 0.0);
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
                        signal_blk[i] *= ir_blk[i];
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
                        output[index_output + 1] = signal_blk[k].re;
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

    pub fn flush(&mut self, output: &mut [f32]) {
        todo!("implement")
    }

}

// TODO: add tests
