use crate::{
    data::{Cartesian2Polar, FrequencyResponseData},
    if64,
};
use indicatif::{ParallelProgressIterator, ProgressStyle};
use rayon::prelude::*;

use super::DPI;

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
    /// structural model natural frequencies
    Structural {
        /// minimum natural frequency [Hz]
        #[arg(long)]
        min: Option<f64>,
        /// maximum natural frequency [Hz]
        #[arg(long)]
        max: Option<f64>,
        /// natural frequency index
        #[arg(short, long)]
        index: Option<Vec<usize>>,
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
    /// frequencies logarithmic sampling
    pub fn logspace(lower: f64, upper: f64, n: usize) -> Self {
        Self::LogSpace { lower, upper, n }
    }
    /// frequencies linear sampling
    pub fn linspace(lower: f64, upper: f64, n: usize) -> Self {
        Self::LinSpace { lower, upper, n }
    }
    /// returns the [data](FrequencyResponseData) associated to the frequencies
    pub fn data<T, F>(self, func: F) -> Vec<FrequencyResponseData<T>>
    where
        T: Cartesian2Polar + Send,
        <T as Cartesian2Polar>::Output: Send,
        F: Fn(f64, if64) -> FrequencyResponseData<T> + Sync,
    {
        let style = ProgressStyle::with_template("|{bar} {pos}|")
            .unwrap()
            .progress_chars("-.-");
        match self {
            Frequencies::Single { value: nu } => {
                let jw = if64::new(0f64, DPI * nu);
                vec![func(nu, jw)]
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
                        func(nu, jw)
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
                        func(nu, jw)
                    })
                    .collect()
            }
            Frequencies::Set { values: nu } => nu
                .into_par_iter()
                .progress_with_style(style)
                .map(|nu| {
                    let jw = if64::new(0f64, DPI * nu);
                    func(nu, jw)
                })
                .collect(),
            _ => panic!("frequencies not set, aborting"),
        }
    }
}
