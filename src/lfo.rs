use crate::ring_buffer::RingBuffer;

pub struct LFO {
    sample_rate_hz: f32, // Sample rate in Hz
    frequency_hz: f32,
    modfreq: f32, // LFO frequency in samples
    phase: f32, 
    amplitude: f32, // LFO depth in samples
    buffer: RingBuffer<f32>,
}

#[derive(Debug, Clone, Copy)]
pub enum FilterParam {
    Frequency,
    Phase,
    Amplitude,
}

#[derive(Debug, Clone)]
pub enum Error {
    InvalidValue { param: FilterParam, value: f32 }
}

impl LFO {
    pub fn new(sample_rate_hz: f32, frequency_hz: f32, amplitude: f32, capacity: usize) -> Self {
        // Initialized variables for self
        let phase = 0.0;

        LFO {
            sample_rate_hz,
            frequency_hz,
            modfreq: frequency_hz / sample_rate_hz,
            phase,
            amplitude,
            buffer: RingBuffer::new(capacity),
        }
    }

    pub fn reset(&mut self) {
        // Reset parameters to default states
        self.buffer.reset();
        self.frequency_hz = 2.0;
        self.phase = 0.0;
        self.amplitude = 0.5;
    }

    pub fn process(&mut self) -> f32 {
        // Update phase increment based on LFO frequency 
        self.phase += self.modfreq as f32;

         // Phase limited to range of [0,1)
        if self.phase >= 1.0 {
            self.phase -= 1.0;
        }

        // Calculate modulation depth using LFO depth
        let value = self.amplitude * self.phase.sin(); 

        // Return Modulation Depth as LFO waveform
        value
    }

    pub fn set_param(&mut self, param: FilterParam, value: f32) -> Result<(), Error> {
        // Optional parameters are frequency and phase
        match param {
            FilterParam::Frequency => {
                if value <= 0.0 {
                    return Err(Error::InvalidValue {param,value});
                }

                self.frequency_hz = value; 
            }
            FilterParam::Phase => {
                if value <= 0.0 {
                    return Err(Error::InvalidValue {param,value});
                }
                self.phase = value;
            }   
            FilterParam::Amplitude => {
                if value <= 0.0 {
                    return Err(Error::InvalidValue {param,value});
                }
                self.phase = value;
            }   
        }
        Ok(())
    }

    pub fn get_param(&self, param: FilterParam) -> f32 {
        // Optional parameters are frequency and phase
        match param {
            FilterParam::Frequency => {
                self.frequency_hz
            }
            FilterParam::Phase => {
                self.phase
            }
            FilterParam::Amplitude => {
                self.amplitude
            }
        }
    }
}