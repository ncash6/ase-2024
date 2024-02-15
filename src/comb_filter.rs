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
        let feedback = 0.5;

        CombFilter {
            delay_len,
            delay_line,
            filter_coefficients,
            sample_rate_hz,
            num_channels,
            feedback, // default feedback value
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

// Define code behavior tests here

#[cfg(test)]
mod tests{
    // use core::num;

    use super::*;

    // Test 1: FIR - Output is zero if input freq matches feedforward
    #[test]
    fn test_feedforward() {
        
        // Parameters to create new FIR filter
        let sample_rate = 44100.0; // in hz
        let max_delay = 0.1; // in sec
        let num_channels = 1;
        let feedforward_hz = 1000.0;

        let mut fir_comb_filter = CombFilter::new(FilterType::FIR, max_delay,sample_rate,num_channels);

        // FIR coefficients
        let mut filter_coefficients = fir_comb_filter.filter_coefficients;
        let feedforward_index = (feedforward_hz * max_delay) as usize;
        filter_coefficients[feedforward_index] = 1.0;
        fir_comb_filter.filter_coefficients = filter_coefficients;

        // Input frequency matches feedforward 
        let signal_freq = feedforward_hz;

        // Sine wave generated as test input
        let sine_dur = 1.0; // in sec
        let num_samp = (sine_dur * sample_rate) as usize;
        let mut in_signal = vec![0.0; num_samp];

        for n in 0..num_samp {
            in_signal[n] = (2.0 * std::f32::consts::PI * signal_freq * (n as f32) / sample_rate).sin()
        }

        // Parameters to process input signal
        let mut out_signal = vec![0.0; num_samp];

        fir_comb_filter.process(&[&in_signal], &mut [&mut out_signal]);

        // Assert output is zero throughout signal
        for &sample in &out_signal {
            assert_eq!(sample, 0.0);
        }
    }

    // Test 2: IIR - amount of magnitude increase/decrease if input freq matches feedback
    #[test]
    fn test_feedback() {
        // Parameters to create new IIR filter
        let sample_rate = 44100.0; // in hz
        let max_delay = 0.1; // in sec
        let num_channels = 1;
        let feedback_hz = 1000.0;

        let mut iir_comb_filter = CombFilter::new(FilterType::IIR, max_delay,sample_rate,num_channels);

        // IIR coefficients
        let feedback_coefficient = (2.0 * std::f32::consts::PI * feedback_hz / sample_rate).cos();
        iir_comb_filter.feedback = feedback_coefficient;

        // Input frequency matches feedforward 
        let signal_freq = feedback_hz;

        // Sine wave generated as test input
        let sine_dur = 1.0; // in sec
        let num_samp = (sine_dur * sample_rate) as usize;
        let mut in_signal = vec![0.0; num_samp];

        for n in 0..num_samp {
            in_signal[n] = (2.0 * std::f32::consts::PI * signal_freq * (n as f32) / sample_rate).sin()
        }

        // Parameters to process input signal
        let mut out_signal = vec![0.0; num_samp];

        iir_comb_filter.process(&[&in_signal], &mut [&mut out_signal]);

        // Assert output magnitude increases/decreases if input matches
        let in_mag = (in_signal.iter().map(|&x| x * x).sum::<f32>() / in_signal.len() as f32).sqrt();
        let out_mag = (out_signal.iter().map(|&x| x * x).sum::<f32>() / out_signal.len() as f32).sqrt();
        
        if feedback_coefficient > 0.0 {
            assert!(out_mag > in_mag, "We expect output magnitude to increase with positive feedback.");
        } else {
            assert!(in_mag > out_mag, "We expect output magnitude to decrease with positive feedback.");
        }
    }

    // Test 3: FIR/IIR - correct result for VARYING input block size
    #[test]
    fn test_blocksize() {

        // Parameters to create new FIR/IIR filter
        let sample_rate = 44100.0; // in hz
        let max_delay = 0.1; // in sec
        let num_channels = 1;
        let signal_freq = 1000.0; // in Hz
        let block_size = 512; // Different block size other than 1024
        // let choose_filter = FilterType::FIR; // FIR filter option 
        let choose_filter= FilterType::IIR; // Uncomment IIR filter option 

        let mut choose_comb_filter = CombFilter::new(choose_filter, max_delay,sample_rate,num_channels);
       
        // IIR coefficients conditional
        if let FilterType::IIR = choose_filter {
            let feedback_coefficient = (2.0 * std::f32::consts::PI * signal_freq / sample_rate).cos();
            choose_comb_filter.feedback = feedback_coefficient;
        }
        
        // Sine wave generated as test input
        let sine_dur = 1.0; // in sec
        let num_samp = (sine_dur * sample_rate) as usize;
        let num_blocks = num_samp / block_size;
        let mut in_signal = vec![0.0; num_samp];

        for n in 0..num_samp {
            in_signal[n] = (2.0 * std::f32::consts::PI * signal_freq * (n as f32) / sample_rate).sin()
        }
        
        let mut out_signal = vec![0.0; num_samp];

        // Block audio and process using comb filter
        for block_indx in 0..num_blocks {
            let start_indx = block_indx * block_size;
            let end_indx = (block_indx as usize + 1) * block_size;
            let in_block = &in_signal[start_indx..end_indx];
            let out_block = &mut out_signal[start_indx..end_indx];
            choose_comb_filter.process(&[in_block], &mut[out_block]);
        }

        // Assert output is correct length with block size not 1024
        assert_eq!(out_signal.len(), num_samp, "Output signal length should match length of the input signal");
    }

    // Test 4: FIR/IIR - correct processing for zero input signal
    #[test]
    fn test_empty_input() {

        // Parameters to create new FIR/IIR filter
        let sample_rate = 44100.0; // in hz
        let max_delay = 0.1; // in sec
        let num_channels = 1;
        let signal_freq = 1000.0; // in Hz
        let sine_dur = 1.0; // in sec
        let num_samp = (sine_dur * sample_rate) as usize;
        let choose_filter = FilterType::FIR; // FIR filter option 
        // let choose_filter= FilterType::IIR; // Uncomment IIR filter option 

        let mut choose_comb_filter = CombFilter::new(choose_filter, max_delay,sample_rate,num_channels);
       
        // IIR coefficients conditional
        if let FilterType::IIR = choose_filter {
            let feedback_coefficient = (2.0 * std::f32::consts::PI * signal_freq / sample_rate).cos();
            choose_comb_filter.feedback = feedback_coefficient;
        }
        
        // Zero input (silence) as test input
        let in_signal = vec![0.0; num_samp];
        
        // Parameters to process input signal
        let mut out_signal = vec![0.0; num_samp];

        // Process silence test using comb filter
        choose_comb_filter.process(&[&in_signal], &mut [&mut out_signal]);

        // Assert output is correct length with block size not 1024
        assert_eq!(out_signal, in_signal, "Output signal should match zeroes of input signal");
    }

    // Test 5: (MEANINGFUL Test) Software Test - ensure tested API functions work
    #[test]
    fn test_api() {
        let sample_rate = 44100.0; // in hz
        let max_delay = 0.1; // in sec
        let num_channels = 2;
        let choose_filter = FilterType::FIR; // FIR filter option 
        // let choose_filter= FilterType::IIR; // Uncomment IIR filter option 

        let mut choose_comb_filter = CombFilter::new(choose_filter, max_delay,sample_rate,num_channels);
        
        // Test for new function
        match (choose_comb_filter.filter_type, choose_filter) {
            (FilterType::FIR, FilterType::FIR) => {
                assert_eq!(choose_comb_filter.delay_len, (max_delay * sample_rate) as usize);
                assert_eq!(choose_comb_filter.delay_line.len(), num_channels);
                assert_eq!(choose_comb_filter.filter_coefficients.len(), (max_delay * sample_rate) as usize);
                assert_eq!(choose_comb_filter.sample_rate_hz, sample_rate);
                assert_eq!(choose_comb_filter.num_channels, num_channels);
                assert_eq!(choose_comb_filter.feedback, 0.5);
            }
            (FilterType::IIR, FilterType::IIR) => {
                assert_eq!(choose_comb_filter.delay_len, (max_delay * sample_rate) as usize);
                assert_eq!(choose_comb_filter.delay_line.len(), num_channels);
                assert_eq!(choose_comb_filter.filter_coefficients.len(), (max_delay * sample_rate) as usize);
                assert_eq!(choose_comb_filter.sample_rate_hz, sample_rate);
                assert_eq!(choose_comb_filter.num_channels, num_channels);
                assert_eq!(choose_comb_filter.feedback, 0.5);
            }
            _=> {
                panic!("Nonmatching filter types: {:?} and {:?}", choose_comb_filter.filter_type, choose_filter);
            }
        }

        // Test for set parameter function
        let gain_feedback = 0.7;
        let delay_value = 0.35;

        choose_comb_filter.set_param(FilterParam::Gain, gain_feedback).unwrap();
        assert_eq!(choose_comb_filter.feedback, gain_feedback);

        choose_comb_filter.set_param(FilterParam::Delay, delay_value).unwrap();
        assert_eq!(choose_comb_filter.delay_len, (delay_value * sample_rate) as usize);

        // Test for reset function
        choose_comb_filter.feedback = 0.42;

        choose_comb_filter.reset();
        assert_eq!(choose_comb_filter.feedback, 0.5);

        // Test for get parameter function
        assert_eq!(choose_comb_filter.get_param(FilterParam::Delay), (max_delay * sample_rate));
        assert_eq!(choose_comb_filter.get_param(FilterParam::Gain), 0.5);


    }
}