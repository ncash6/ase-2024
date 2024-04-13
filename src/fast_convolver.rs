use crate::ring_buffer::RingBuffer;

struct FastConvolver {
    n_ir: usize,
    impulse_response: &'static[f32],
    mode: ConvolutionMode,
}

#[derive(Debug, Clone, Copy)]
pub enum ConvolutionMode {
    TimeDomain,
    FrequencyDomain { block_size: usize },
}

impl FastConvolver {
    pub fn new(impulse_response: &'static[f32], mode: ConvolutionMode) -> Self {
        // Initialized variables for self
        let n_ir = impulse_response.len() as usize; // Lengths of Impulse
        
        FastConvolver {
            impulse_response,
            mode,
            n_ir,
        }
    }

    pub fn reset(&mut self) {
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

                    output[k] = current_sum // Convolution results stored in output buffer
                }
            }
            ConvolutionMode::FrequencyDomain { block_size } => {
                // Partitioned Fast Convolution Frequency Domain
                
            }
        }
    }

    pub fn flush(&mut self, output: &mut [f32]) {
        todo!("implement")
    }

    // TODO: feel free to define other functions for your own use
}

// TODO: add tests
