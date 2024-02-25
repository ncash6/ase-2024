// The user-interface for this vibrato.rs file is modelled after the comb_filter.rs
// written by Stephen Garrett & Ian Clester for MUSI 6106 Assignment 1.

pub struct Vibrato {
    // TODO: your code here
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

impl Vibrato {
    pub fn new(filter_type: FilterType, max_delay_secs: f32, sample_rate_hz: f32, num_channels: usize) -> Self {
        todo!("implement")
    }

    pub fn reset(&mut self) {
        todo!("implement")
    }

    pub fn process(&mut self, input: &[&[f32]], output: &mut [&mut [f32]]) {
        todo!("implement");
    }

    pub fn set_param(&mut self, param: FilterParam, value: f32) -> Result<(), Error> {
        todo!("implement")
    }

    pub fn get_param(&self, param: FilterParam) -> f32 {
        todo!("implement")
    }

    // TODO: feel free to define other functions for your own use
}

// TODO: feel free to define other types (here or in other modules) for your own use