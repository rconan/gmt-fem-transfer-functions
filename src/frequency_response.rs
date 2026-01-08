//! Frequency response functionalities

use indicatif::{ParallelProgressIterator, ProgressStyle};
use rayon::prelude::*;
use std::{f64::consts::PI, ops::Mul};

use crate::{
    data::{Cartesian2Polar, FrequencyResponseData, FrequencyResponseVec},
    if64,
};

const DPI: f64 = 2f64 * PI;

/// Frequency sampling options
///
/// The frequencies units is Hz
#[derive(Debug, Clone, clap::Subcommand)]
#[command(
    subcommand_help_heading = "Transfer functions sampling frequencies [Hz]",
    subcommand_value_name = "SAMPLING FREQUENCIES"
)]
pub enum Frequencies {
    /// a single frequency
    Single { value: f64 },
    /// logarithmic (log base 10) sampling of the interval `[lower,upper]` with `n` samples
    LogSpace {
        #[arg(short, long)]
        lower: f64,
        #[arg(short, long)]
        upper: f64,
        #[arg(short)]
        n: usize,
    },
    /// regular sampling of the interval `[lower,upper]` with `n` samples
    LinSpace {
        #[arg(short, long)]
        lower: f64,
        #[arg(short, long)]
        upper: f64,
        #[arg(short)]
        n: usize,
    },
    /// a given set of frequencies
    Set {
        #[arg(short, long)]
        values: Vec<f64>,
    },
}
impl From<f64> for Frequencies {
    fn from(value: f64) -> Self {
        Frequencies::Single { value }
    }
}
impl From<Vec<f64>> for Frequencies {
    fn from(values: Vec<f64>) -> Self {
        Frequencies::Set { values }
    }
}
impl From<&Vec<f64>> for Frequencies {
    fn from(values: &Vec<f64>) -> Self {
        Frequencies::Set {
            values: values.clone(),
        }
    }
}
impl From<&Self> for Frequencies {
    fn from(value: &Self) -> Self {
        value.clone()
    }
}
impl Frequencies {
    pub fn logspace(lower: f64, upper: f64, n: usize) -> Self {
        Self::LogSpace { lower, upper, n }
    }
    pub fn linspace(lower: f64, upper: f64, n: usize) -> Self {
        Self::LinSpace { lower, upper, n }
    }
}

/// Frequency response interface definition
pub trait FrequencyResponse {
    /// Transfer function type
    type Output;

    /// Returns the frequency response
    ///
    /// The argument is the imaginary frequency in radians
    fn j_omega(&self, jw: if64) -> Self::Output;
    /// Returns the frequencies and the frequency response
    ///
    /// The argument is frequencies in Hz
    fn frequency_response<T: Into<Frequencies>>(&self, nu: T) -> FrequencyResponseVec<Self::Output>
    where
        <Self as FrequencyResponse>::Output: Cartesian2Polar + Send,
        <<Self as FrequencyResponse>::Output as Cartesian2Polar>::Output: Send,
        Self: Sync,
    {
        let frequencies: Frequencies = nu.into();
        let style = ProgressStyle::with_template("|{bar} {pos}|")
            .unwrap()
            .progress_chars("-.-");
        let data = match frequencies {
            Frequencies::Single { value: nu } => {
                let jw = if64::new(0f64, DPI * nu);
                vec![FrequencyResponseData::new(nu, self.j_omega(jw))]
            }
            Frequencies::LogSpace { lower, upper, n } => {
                assert!(upper > lower);
                let log_step = (upper.log10() - lower.log10()) / (n - 1) as f64;
                (0..n)
                    .into_par_iter()
                    .progress_with_style(style)
                    .map(|i| {
                        let log_nu = lower.log10() + log_step * i as f64;
                        let nu = 10f64.powf(log_nu);
                        let jw = if64::new(0f64, DPI * nu);
                        FrequencyResponseData::new(nu, self.j_omega(jw))
                    })
                    .collect()
            }
            Frequencies::LinSpace { lower, upper, n } => {
                assert!(upper > lower);
                let step = (upper - lower) / (n - 1) as f64;
                (0..n)
                    .into_par_iter()
                    .progress_with_style(style)
                    .map(|i| {
                        let nu = lower + step * i as f64;
                        let jw = if64::new(0f64, DPI * nu);
                        FrequencyResponseData::new(nu, self.j_omega(jw))
                    })
                    .collect()
            }
            Frequencies::Set { values: nu } => nu
                .into_par_iter()
                .progress_with_style(style)
                .map(|nu| {
                    let jw = if64::new(0f64, DPI * nu);
                    FrequencyResponseData::new(nu, self.j_omega(jw))
                })
                .collect(),
        };
        FrequencyResponseVec::new(data)
    }
    /// Returns the first derivation of the frequency response
    fn j_omega_first(&self, jw: if64) -> <<Self as FrequencyResponse>::Output as Mul<if64>>::Output
    where
        <Self as FrequencyResponse>::Output: Mul<if64>,
    {
        self.j_omega(jw) * jw
    }
    /// Returns the second derivation of the frequency response
    fn j_omega_second(
        &self,
        jw: if64,
    ) -> <<<Self as FrequencyResponse>::Output as Mul<if64>>::Output as Mul<if64>>::Output
    where
        <Self as FrequencyResponse>::Output: Mul<if64>,
        <<Self as FrequencyResponse>::Output as Mul<if64>>::Output: Mul<if64>,
    {
        self.j_omega_first(jw) * jw
    }
}

/// First order low-pass
///
/// *GMT-DOC-XXXX: ASM segment modal tranfer function*, Eq.(1)
#[derive(Debug)]
pub struct FirstOrderLowPass {
    corner_frequency_hz: f64,
}
impl FirstOrderLowPass {
    pub fn new() -> Self {
        Self {
            corner_frequency_hz: 4e3,
        }
    }
}
impl FrequencyResponse for FirstOrderLowPass {
    type Output = if64;
    fn j_omega(&self, jw: if64) -> Self::Output {
        jw / (1f64 + jw / (DPI * self.corner_frequency_hz))
    }
}

/// 4th-order bessel filter
///
/// *GMT-DOC-XXXX: ASM segment modal tranfer function*, Eq.(2)
#[derive(Debug)]
pub struct BesselFilter {
    w_bf: f64,
    beta: [f64; 5],
}
impl BesselFilter {
    pub fn new() -> Self {
        Self {
            w_bf: DPI * 2.2e3,
            beta: [1f64, 3.20108587, 4.39155033, 3.12393994, 1f64],
        }
    }
}
impl FrequencyResponse for BesselFilter {
    type Output = if64;
    fn j_omega(&self, jw: if64) -> Self::Output {
        let num = self.beta[0] * self.w_bf.powi(4);
        let denom = self
            .beta
            .iter()
            .enumerate()
            .fold(if64::new(0f64, 0f64), |a, (i, b)| {
                a + b * self.w_bf.powi(4 - i as i32) * jw.powi(i as i32)
            });
        num / denom
    }
}

/// Proportional-integral compensator
///
/// *GMT-DOC-XXXX: ASM segment modal tranfer function*, Eq.(3)
#[derive(Debug)]
pub struct PICompensator {
    kp: f64,
    ki: f64,
}
impl PICompensator {
    pub fn new() -> Self {
        Self { kp: 7e4, ki: 5e5 }
    }
}
impl FrequencyResponse for PICompensator {
    type Output = if64;
    fn j_omega(&self, jw: if64) -> Self::Output {
        self.kp + self.ki / jw
    }
}

#[cfg(test)]
mod tests {
    // use std::fs::File;

    use super::*;

    #[test]
    fn folp_tf() {
        let folp = FirstOrderLowPass::new();

        let tf = folp.frequency_response(Frequencies::logspace(1., 8e3, 1000));

        // let mut file = File::create("folp_tf.pkl").unwrap();
        // serde_pickle::to_writer(&mut file, &(nu, tf), Default::default()).unwrap();
    }

    #[test]
    fn bessel_tf() {
        let bessel = BesselFilter::new();

        let tf = bessel.frequency_response(Frequencies::logspace(1., 8e3, 1000));

        // let mut file = File::create("bessel_tf.pkl").unwrap();
        // serde_pickle::to_writer(&mut file, &(nu, tf), Default::default()).unwrap();
    }

    #[test]
    fn pic_tf() {
        let pic = PICompensator::new();

        let tf = pic.frequency_response(Frequencies::logspace(1., 8e3, 1000));

        // let mut file = File::create("pic_tf.pkl").unwrap();
        // serde_pickle::to_writer(&mut file, &(nu, tf), Default::default()).unwrap();
    }
}
