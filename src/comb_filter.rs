pub struct CombFilter {
    // TODO: your code here
    delay_len: usize,
    delay_line: Vec<Vec<f32>>,
    filter_coefficients: Vec<f64>,
    feedback: f32,
    filter_type: FilterType,
    sample_rate_hz: f32,
    num_channels: usize,
}

#[derive(Debug, Clone, Copy)]
pub enum FilterType {
    FIR,
    IIR,
}

#[derive(Debug, Clone, Copy)]
pub enum FilterParam {
    Gain,
    Delay,
}

#[derive(Debug, Clone)]
pub enum Error {
    InvalidValue { param: FilterParam, value: f32 }
}

impl CombFilter {
    pub fn new(filter_type: FilterType, max_delay_secs: f32, sample_rate_hz: f32, num_channels: usize) -> Self {
        // Initialized variables for self
        let delay_len = (max_delay_secs * sample_rate_hz) as usize;
        let delay_line = vec![vec![0.0; delay_len]; num_channels];
        let filter_coefficients = vec![0.0; delay_len];

        CombFilter {
            delay_len,
            delay_line,
            filter_coefficients,
            sample_rate_hz,
            num_channels,
            feedback: 0.5, // default feedback value
            filter_type,
        }
    }

    pub fn reset(&mut self) {
        // Zeroed delay line
        for channel in &mut self.delay_line {
            for sample in channel {
                *sample = 0.0;
            }
        }
        // Reset mutable states
        self.feedback = 0.5;
    }

    pub fn process(&mut self, input: &[&[f32]], output: &mut [&mut [f32]]) {
        match self.filter_type {
            FilterType::FIR => {
                // Matlab FIR Implementation
                for (channel, & in_channel) in input.iter().enumerate() {
                    let fir_delay_length = self.delay_line[channel].len();
                    let filter_coeff_length = self.filter_coefficients.len();

                    for (samp, &in_samp) in in_channel.iter().enumerate() {
                        let mut out_del_sample = 0.0;

                        // FIR difference equation
                        for n in 0..filter_coeff_length {
                            if samp >= n {
                                out_del_sample += in_channel[samp - n] * (self.filter_coefficients[n] as f32);
                            }
                        }
                        output[channel][samp] = out_del_sample;

                        // delay line update
                        let mut update_delay_line = vec![0.0; fir_delay_length];
                        update_delay_line[0] = in_samp;
                        update_delay_line[1..].copy_from_slice(&self.delay_line[channel][0..(fir_delay_length - 1)]);
                        self.delay_line[channel] = update_delay_line;
                    }
                }
            }
            FilterType::IIR => {
                // Matlab IIR Implementation
                for (channel, & in_channel) in input.iter().enumerate() {
                    let iir_delay_length = self.delay_line[channel].len();

                    for (samp, &in_samp) in in_channel.iter().enumerate() {
                        // IIR difference equation
                        let out_del_sample = self.delay_line[channel][iir_delay_length - 1] * self.feedback;
                        output[channel][samp] = in_samp + out_del_sample;

                        // delay line update
                        let mut update_delay_line = vec![0.0; iir_delay_length];
                        update_delay_line[0] = output[channel][samp];
                        update_delay_line[1..].copy_from_slice(&self.delay_line[channel][0..(iir_delay_length - 1)]);
                        self.delay_line[channel] = update_delay_line;
                    }
                }
            }   
        }
        todo!();
    }

    pub fn set_param(&mut self, param: FilterParam, value: f32) -> Result<(), Error> {
        // Optional parameters are gain and delay
        match param {
            FilterParam::Gain => {
                if value <= 0.0 {
                    return Err(Error::InvalidValue {param,value});
                }

                self.feedback = value; 
            }
            FilterParam::Delay => {
                if value <= 0.0 {
                    return Err(Error::InvalidValue {param,value});
                }
                self.delay_len = (value * self.sample_rate_hz) as usize;
            }   
        }
        Ok(())
    }

    pub fn get_param(&self, param: FilterParam) -> f32 {
        // Optional parameters are gain and delay
        match param {
            FilterParam::Gain => {
                self.feedback
            }
            FilterParam::Delay => {
                self.delay_line.len() as f32
            }
        }
    }
    
    // TODO: feel free to define other functions for your own use
}

// TODO: feel free to define other types (here or in other modules) for your own use
